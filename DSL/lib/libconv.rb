require_relative '../../boast/GenericConvolution3'
require_relative '../../boast/WaveletFilters'
require_relative '../../boast/PoissonFilters'
require_relative '../../boast/Synthesis.rb'
require_relative '../../boast/MagicFilter2.rb'
require_relative '../../boast/AnaRotPer-2.rb'
require_relative 'libconv/dimension'
require_relative 'libconv/data_space'
require_relative 'libconv/boundary_condition'


module LibConv
#default values
  @optims = {unroll_range: [1,1], mod_arr_test: false, tt_arr_test: false, openmp: true}
  @precisions=[4,8]
  @all_wavelet_families=["SYM8", "SYM4"]
  @wavelet_families=@all_wavelet_families
  @all_operations=["MF","IMF","DWT","RTOS1","IDWT","S1TOR","D1","D2","NABLA"]
  @operations=@all_operations
  @all_families=["s0s0", "s0s0_dot", "s0s1", "s1s0"]
  @families=@all_families
  @dimensions=["1d"]
  @bench=false
  @link_with_simgrid=false
  @refs=true
  #read input file if present
  if (ARGV.length == 1) then
    input = ARGV[0]
    f = YAML::load(File::open(input).read)

    f.each{|key,value|
      if key == :libconv then
        value.each{|key2, value2|
          if key2 == :precisions
            @precisions = value2
          elsif key2 == :wavelet_families
            @wavelet_families = value2
          elsif key2 == :operations
            @operations = value2
          elsif key2 == :families
            @families = value2
          elsif key2 == :dimensions
            @dimensions = value2
          elsif key2 == :simgrid
            @link_with_simgrid = value2
          elsif key2 == :bench
            @bench = value2
          elsif key2 == :refs
            @refs = value2
          else
            raise "unknown key in libconv config file : "+key2
          end
        }
      elsif key == :optims
        @optims= value
      else
        raise "unknown key in libconv config file : "+key
      end  
    }
  end

  require_relative 'libconv/lib_conv_kernel'


  puts "Libconv : generating convolutions in "+ ((BOAST::get_lang == BOAST::C) ? "C" : "Fortran")
  puts "Operations : "+ @operations.inspect
  puts "Precisions : "+ @precisions.inspect
  puts "Filters families : "+ @wavelet_families.inspect
  puts "Optimizations : "+ @optims.inspect
  if @optims.has_key? :vector_length and BOAST::get_lang != BOAST::C then puts "WARNING : vectorization asked in Fortran, will be ignored" end 
  puts "Output dir : " + @foldername + " config file " + input
  if @from_cache then puts "WARNING : using cached convolutions. Purge .cache dir if not wanted. Headers and brokers will be generated in current dir" end

if not LibConv.from_cache then
  Dir.mkdir("#{@foldername}") if not Dir.exist?(@foldername)
  FileUtils::copy(ARGV[0],@foldername)
end

  def self.const_name(operation, config)

    wavelet_family = config[:wavelet_family]
    precision = config[:precision]
    case precision
    when 4
      precision_name = "S"
    when 8
      precision_name = "D"
    end

    return "#{precision_name}_#{wavelet_family}_#{operation}"

  end


  CONFIGURATION = BOAST::BruteForceOptimizer::new( BOAST::OptimizationSpace::new( {precision: @precisions, wavelet_family: @wavelet_families } ) )

  CONFIGURATION.each { |config|
    wavelet_family = config[:wavelet_family]
    wavelet_name = wavelet_family.downcase
    precision = config[:precision]

    BOAST.push_env( default_real_size: precision ) {
      wvals = const_get(wavelet_family+"_LP")
      mfvals = const_get(wavelet_family+"_MF")
      
    if @operations != nil then
      @operations.each{ |op| 
        name = const_name(op, config)
        case op
        when "D2"  
          kf2 = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_d2', const_get(wavelet_family+"_D2"), const_get(wavelet_family+"_D2").length/2)
          const_set(name, KineticKernel1d.new( kf2, :kd2))
        when "D1"
          kf1 = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_d1', const_get(wavelet_family+"_D1"), const_get(wavelet_family+"_D1").length/2)
          const_set(name, KineticKernel1d.new( kf1, :kd1))
        when "DWT"
          wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name, wvals)
          const_set(name, WaveletKernel1d.new( wf, :dwt))
        when "IDWT"
          iwf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name, wvals)
          const_set(name, WaveletKernel1d.new( iwf, :idwt))
        when "MF"
          cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
          const_set(name, MagicFilterKernel1d.new( cf, :mf))
        when "IMF"
          icf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_imd', mfvals, mfvals.length/2 - 1)
          const_set(name, MagicFilterKernel1d.new( icf, :imf))
        when "S1TOR"
          cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
          s1s0wf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
          const_set(name, WaveletKernel1d.new( s1s0wf, :s1tor ))
        when "RTOS1"
          cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
          s0s1wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
          const_set(name, WaveletKernel1d.new( s0s1wf, :rtos1 ))
        when "NABLA"
          #valid but done later as it's a different family
        else
          raise "Unknown operation"+op
        end
      }
    end
    }
  }

  if @operations != nil and @operations.include? "NABLA"
  then
    CONFIGURATION_POISSON = BOAST::BruteForceOptimizer::new( BOAST::OptimizationSpace::new( {precision: @precisions, wavelet_family: ["NABLA"]}))
    CONFIGURATION_POISSON.each { |config|
      precision = config[:precision]
      BOAST.push_env( default_real_size: precision ) {
        nords=[2,4,6,8,16]
        nords.each{ |nord_n|
          filter=const_get("NABLA_"+nord_n.to_s)
          conv_filter = GenericConvolution::PoissonFilter::new('poisson'+nord_n.to_s,filter.each_slice(nord_n+1).to_a,nord_n)
          const_set(const_name("PS"+nord_n.to_s, config), PoissonKernel1d.new( conv_filter, :ps))
        }
      }
    }
  end

  if @refs == true
  then
    #print_reference_kernels()
    filename="refs.f90"
    Dir.mkdir("#{@foldername}") if not Dir.exist?(@foldername)
    File::open(@foldername+filename,"w") {|f|
      sref = BOAST::synthesis_per_ref
      f.puts sref
      aref = BOAST::analysis_free_ref
      f.puts aref
      mf_ref = magicfilter_ref
      f.puts mf_ref

    }
    @kernels.push(filename)
  end

end
require_relative 'libconv/brokers'
