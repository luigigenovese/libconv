require 'fileutils'
require 'yaml'

module BigDFT

  class << self
    attr_accessor :optims
    attr_reader   :from_cache
    attr_reader   :kernels
    attr_reader   :foldername
  end
  @optims = { unroll_range: [1, 1], mod_arr_test: false, tt_arr_test: false, openmp: true }
  @foldername=Time.now().strftime("%Y-%m-%d_%H-%M-%S/")
  @kernels=[]
  unless File.directory?('.cache')
    FileUtils.mkdir_p('.cache')
  end

  begin
    optims = YAML.load_file('.cache/optims.yml')
    @from_cache = if optims == @optims
                    true
                  else
                    false
                  end
  rescue
    @from_cache = false
  end

  if not @from_cache
    File.open('.cache/optims.yml', 'w') { |f|
      f.write(YAML.dump(@optims.to_h))
    }
  end

  def self.dump_to_file(kernel)
    filename=kernel.procedure.name
    case get_lang
    when C
      suffix = ".c"
    when FORTRAN
      suffix = ".f90"
    end
    filedump="#{filename}#{suffix}"
    Dir.mkdir("#{@foldername}") if not Dir.exist?(@foldername)
    File::open(@foldername+filedump,"w") {|f|
      f.puts kernel
    }
    @kernels.push(filedump)
  end
  class LibConvKernel

    attr_reader :kernel
    attr_reader :cost_procedure
    attr_reader :dims_procedure
    attr_reader :op
    attr_reader :filter
    attr_reader :default_options

    def initialize(conv_operation, op, optims)
      @op = op
      if BigDFT.from_cache and File.exist?(".cache/#{conv_operation.procedure_name}.#{RbConfig::CONFIG['DLEXT']}")
        BOAST.push_env( ffi: true ) {
          p = conv_operation.empty_procedure
          @kernel = BOAST::CKernel.new
          @kernel.procedure = p
          @kernel.cost_function = ->(*args) { conv_operation.cost(*args) }
          @kernel.build(library_path: ".cache/#{p.name}.#{RbConfig::CONFIG['DLEXT']}")
        }
      else
        conv_operation.optimize(GenericConvolution::GenericOptimization::new(optims)) if optims

        p, subops = conv_operation.procedure

        @kernel = BOAST::CKernel.new

        GenericConvolution::print_header

        subops.each_value do |sub_op|
          BOAST.pr sub_op
          puts 'chosen:' + sub_op.name
        end
        BOAST.pr p
        @kernel.procedure = p
        @kernel.cost_function = ->(*args) { conv_operation.cost(*args) }
        @kernel.build
        @kernel.dump_module( ".cache/#{p.name}.#{RbConfig::CONFIG['DLEXT']}" )

        BigDFT.dump_to_file(@kernel)
      end
      @cost_procedure = conv_operation.empty_procedure(:cost)
      @dims_procedure = conv_operation.empty_procedure(:dims)
    end

    def dims_from_in(in_dim, bc)
      case bc
      when GenericConvolution::BC::PERIODIC, GenericConvolution::BC::FREE
        [ in_dim, in_dim, in_dim ]
      when GenericConvolution::BC::GROW
        [ in_dim, in_dim, in_dim + buffer_increment ]
      when GenericConvolution::BC::SHRINK
        [ in_dim, in_dim, in_dim - buffer_increment ]
      else
        raise "Unsupported BC: #{bc}"
      end
    end

    def dims_from_out(out_dim, bc)
      case bc
      when GenericConvolution::BC::PERIODIC, GenericConvolution::BC::FREE
        [ out_dim, out_dim, out_dim ]
      when GenericConvolution::BC::GROW
        [ out_dim - buffer_increment, out_dim - buffer_increment, out_dim ]
      when GenericConvolution::BC::SHRINK
        [ out_dim + buffer_increment, out_dim + buffer_increment, out_dim ]
      else
        raise "Unsupported BC: #{bc}"
      end
    end

  end

  class MagicFilterKernel1d < LibConvKernel
    def initialize(filter, op, optims = BigDFT.optims)
      @filter = filter
      conv_operation = GenericConvolution::GenericConvolutionOperator1d.new(filter, ld: true, narr: true, a_x: true, a_y: true, a: true)
      @default_options = { narr: 1, a_x: 0.0, a_y: 0.0, a: 1.0 }
      super( conv_operation, op, optims )
    end

    def buffer_increment
      @filter.buffer_increment
    end

    def dimension_space
      :r
    end

    def run(idim, bc, source, target, **options)
      opts = @default_options.merge(options)
#      puts @kernel.procedure.name
#      puts([source.system.dimension,
#                  idim,
#                  source.dimensions,
#                  bc,
#                  source.leading_dimensions,
#                  target.leading_dimensions,
#                  opts[:narr],
#                  source.data_space.data,
#                  target.data_space.data,
#                  opts[:a],
#                  opts[:a_x],
#                  opts[:a_y]].collect(&:inspect))
      @kernel.run(source.system.dimension,
                  idim,
                  source.dimensions,
                  bc,
                  source.leading_dimensions,
                  target.leading_dimensions,
                  opts[:narr],
                  source.data_space.data,
                  target.data_space.data,
                  opts[:a],
                  opts[:a_x],
                  opts[:a_y])
    end

  end

  class WaveletKernel1d < LibConvKernel
    def initialize(wavelet_filter, op, optims = BigDFT.optims)
      @filter = wavelet_filter
      conv_operation = GenericConvolution::GenericConvolutionOperator1d.new(wavelet_filter, a: true, a_y: true, ld: true, narr: true)
      @default_options = { narr: 1, a_y: 0.0, a: 1.0 }
      super( conv_operation, op, optims )
    end

    def buffer_increment
      @filter.buffer_increment
    end

    def dimension_space
      :s1
    end

    def run(idim, bc, source, target, **options)
      opts = @default_options.merge(options)
#      puts @kernel.procedure.name
#      puts([source.system.dimension,
#                  idim,
#                  source.dimensions,
#                  bc,
#                  source.leading_dimensions,
#                  target.leading_dimensions,
#                  opts[:narr],
#                  source.data_space.data,
#                  target.data_space.data,
#                  opts[:a],
#                  opts[:a_y]].collect(&:inspect))
      @kernel.run(source.system.dimension,
                  idim,
                  source.dimensions,
                  bc,
                  source.leading_dimensions,
                  target.leading_dimensions,
                  opts[:narr],
                  source.data_space.data,
                  target.data_space.data,
                  opts[:a],
                  opts[:a_y])
    end

  end

class PoissonKernel1d < LibConvKernel
def initialize(conv_filter, op, optims = BigDFT.optims)
    @filter = conv_filter
    conv_operation = GenericConvolution::GenericConvolutionOperator1d.new(conv_filter, ld: false, narr: false, a: true, poisson: true)
    @default_options = { narr: 0, a_y: 0.0, a: 1.0 }
    super( conv_operation, op, optims )
    end


    def buffer_increment
      @filter.buffer_increment
    end

    def dimension_space
      :s0
    end

    def run(idim, bc, source, target, **options)
      opts = @default_options.merge(options)
      @kernel.run(source.system.dimension,
                  idim,
                  source.dimensions,
                  bc,
                  source.leading_dimensions,
                  target.leading_dimensions,
                  source.data_space.data,
                  target.data_space.data,
                  opts[:a])
    end
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

  CONFIGURATION = BOAST::BruteForceOptimizer::new( BOAST::OptimizationSpace::new( {precision: [4, 8], wavelet_family: ["SYM8"] } ) )

  CONFIGURATION.each { |config|
    wavelet_family = config[:wavelet_family]
    wavelet_name = wavelet_family.downcase
    precision = config[:precision]

    BOAST.push_env( default_real_size: precision ) {
      wvals = const_get(wavelet_family+"_LP")
      mfvals = const_get(wavelet_family+"_MF")

      wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name, wvals)
      const_set(const_name("DWT", config), WaveletKernel1d.new( wf, :dwt))

      iwf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name, wvals)
      const_set(const_name("IDWT", config), WaveletKernel1d.new( iwf, :idwt))

      cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
      const_set(const_name("MF", config), MagicFilterKernel1d.new( cf, :mf))

      icf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_imd', mfvals, mfvals.length/2 - 1)
      const_set(const_name("IMF", config), MagicFilterKernel1d.new( icf, :imf))

      s1s0wf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
      const_set(const_name("S1TOR", config), WaveletKernel1d.new( s1s0wf, :s1tor ))

      s0s1wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
      const_set(const_name("RTOS1", config), WaveletKernel1d.new( s0s1wf, :rtos1 ))
      
    nords=[2,4,6,8,16]
    nords.each{ |nord_n|
        filter=const_get("NABLA_"+nord_n.to_s)
        conv_filter = GenericConvolution::PoissonFilter::new('poisson'+nord_n.to_s,filter.each_slice(nord_n+1).to_a,nord_n)
        const_set(const_name("PS"+nord_n.to_s, config), PoissonKernel1d.new( conv_filter, :ps))
    }
    
    }
  }


end
