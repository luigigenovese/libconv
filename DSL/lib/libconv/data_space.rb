module LibConv

  class DataSpace
    attr_reader :data
    attr_reader :dimensions
    attr_reader :dimension
    attr_reader :precision
    attr_reader :alignment
    attr_reader :m

    def initialize(*shape, **options)
      @precision = options.fetch(:precision, 8)
      @alignment = options.fetch(:alignment, 32)
      @m         = options.fetch(:m, 1)
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
      dims = @dimensions
      dims += [@m] if @m > 1
      @data = ANArray.new(type, @alignment, *dims)
      @data.random! if options[:random]
    end

  end

end
