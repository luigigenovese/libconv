require 'fileutils'
require 'yaml'

module LibConv

  class << self
    attr_accessor :optims
    attr_reader   :from_cache
    attr_reader   :kernels
    attr_reader   :foldername
  end
#  @optims = { unroll_range: [1, 1], mod_arr_test: false, tt_arr_test: false, openmp: true }
#  @optims = { unroll_range: [2,8,2], mod_arr_test: false, tt_arr_test: false, openmp: true, vector_length: [2,4]}
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
    attr_reader :align_procedure
    attr_reader :dims_procedure
    attr_reader :op
    attr_reader :filter
    attr_reader :default_options

    def initialize(conv_operation, op, optims)
      @op = op
      if LibConv.from_cache and File.exist?(".cache/#{conv_operation.procedure_name}.#{RbConfig::CONFIG['DLEXT']}")
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

        LibConv.dump_to_file(@kernel)
      end
      @cost_procedure = conv_operation.empty_procedure(:cost)
      @dims_procedure = conv_operation.empty_procedure(:dims)
      @align_procedure = conv_operation.empty_procedure(:align)
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
    def initialize(filter, op, optims = LibConv.optims)
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
    def initialize(wavelet_filter, op, optims = LibConv.optims)
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


  class KineticKernel1d < LibConvKernel
    def initialize(filter, op, optims = LibConv.optims)
      @filter = filter
      conv_operation = GenericConvolution::GenericConvolutionOperator1d.new(filter,kinetic: :inplace, ld: true, narr: true, a_x: true, a_y: true, a: true, dot_in: true, transpose: true)
      @default_options = { narr: 1.0, a_x: 1.0/3.0, a_y: 1.0, a: 1.0, dot_in: 0.0 }
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
                  opts[:a_y],
                  opts[:dot_in])
    end

  end


class PoissonKernel1d < LibConvKernel
def initialize(conv_filter, op, optims = LibConv.optims)
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
end
