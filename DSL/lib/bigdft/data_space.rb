module BigDFT

  class DataSpace
    attr_reader :data
    attr_reader :dimensions
    attr_reader :dimension
    attr_reader :precision
    attr_reader :alignment

    def initialize(*shape, **options)
      @precision = options.fetch(:precision, 8)
      @alignment = options.fetch(:alignment, 32)
      @dimensions = shape.dup
      @dimension = shape.length
      case @precision
      when 4
        type = NArray::SFLOAT
      when 8
        type = NArray::FLOAT
      else
        raise "Unsupported precision (#{@precision})!"
      end
      @data = ANArray.new(type, @alignment, *@dimensions)
      @data.random! if options[:random]
    end

  end

end
