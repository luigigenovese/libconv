require 'fileutils'
require 'yaml'

module BigDFT

  class << self
    attr_accessor :optims
    attr_reader   :from_cache
  end
  @optims = { unroll_range: [1, 1], mod_arr_test: false, tt_arr_test: false }

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

  class LibConvKernel

    attr_reader :kernel
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
      end
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

#    def in_dim(n, bc)
#      case bc
#      when GenericConvolution::BC::PERIODIC, GenericConvolution::BC::GROW
#        n
#      when GenericConvolution::BC::SHRINK
#        n + buffer_increment
#      else
#        raise "Unsupported BC: #{n}"
#      end
#    end
#
#    def in_dims(idim, dimensions, boundary_conditions)
#      in_dims = dimensions.dup
#      in_dims[idim] = in_dim(dimensions[idim], boundary_conditions[idim])
#      in_dims
#    end
#
#    def out_dim(n, bc)
#      case bc
#      when GenericConvolution::BC::PERIODIC, GenericConvolution::BC::SHRINK
#        n
#      when GenericConvolution::BC::GROW
#        n + buffer_increment
#      else
#        raise "Unsupported BC: #{n}"
#      end
#    end
#
#    def out_dims(idim, dimensions, boundary_conditions)
#      out_dims = dimensions.dup
#      out_dims[idim] = out_dim(dimensions[idim], boundary_conditions[idim])
#      out_dims
#    end
#
#    # to modify for vector implementations
#    def out_ld(n, bc)
#      out_dim(n, bc)
#    end
#
#    def out_lds(idim, dimensions, boundary_conditions)
#      out_lds = dimensions.dup
#      out_lds[idim] = out_ld(dimensions[idim], boundary_conditions[idim])
#      out_lds
#    end
#
#    def in_ld(n, bc)
#      in_dim(n, bc)
#    end
#
#    def in_lds(idim, dimensions, boundary_conditions)
#      in_lds = dimensions.dup
#      in_lds[idim] = in_ld(dimensions[idim], boundary_conditions[idim])
#      in_lds
#    end
#
#    def run(dimension_index,
#            dimensions,
#            boundary_conditions,
#            source,
#            destination,
#            *coefficients)
#      @kernel.run(dimensions.length,
#                  dimension_index,
#                  in_dims(dimension_index, dimensions, boundary_conditions),
#                  boundary_conditions,
#                  in_lds(dimension_index, dimensions, boundary_conditions),
#                  out_lds(dimension_index, dimensions, boundary_conditions),
#                  1,
#                  source,
#                  destination,
#                  *coefficients)
#    end

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

  CONFIGURATION = BOAST::BruteForceOptimizer::new( BOAST::OptimizationSpace::new( {precision: [4, 8], wavelet_family: ["SYM8"] } ) )

  CONFIGURATION.each { |config|
    wavelet_family = config[:wavelet_family]
    wavelet_name = wavelet_family.downcase
    precision = config[:precision]
    case precision
    when 4
      precision_name = "S"
    when 8
      precision_name = "D"
    end

    BOAST.push_env( default_real_size: precision ) {
      wvals = const_get(wavelet_family+"_LP")
      mfvals = const_get(wavelet_family+"_MF")

      wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name, wvals)
      const_set(precision_name+"_DWT", WaveletKernel1d.new( wf, :dwt))

      wf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name, wvals)
      const_set(precision_name+"_IDWT", WaveletKernel1d.new( wf, :idwt))

      cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
      const_set(precision_name+"_MF", MagicFilterKernel1d.new( cf, :mf))

      icf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_imd', mfvals, mfvals.length/2 - 1)
      const_set(precision_name+"_IMF", MagicFilterKernel1d.new( icf, :imf))

      cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
      wf = GenericConvolution::WaveletFilterRecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
      const_set(precision_name+"_S1TOR", WaveletKernel1d.new( wf, :s1tor ))

      cf = GenericConvolution::ConvolutionFilter.new(wavelet_name+'_md', mfvals.reverse, mfvals.length/2)
      wf = GenericConvolution::WaveletFilterDecompose.new(wavelet_name+"icomb", wvals, convolution_filter: cf)
      const_set(precision_name+"_RTOS1", WaveletKernel1d.new( wf, :rtos1 ))
    }
  }


end
