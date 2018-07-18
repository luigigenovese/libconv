module BOAST

  class WaveFunction
    attr_accessor :data_space
    attr_reader   :system
    attr_reader   :spaces
    attr_reader   :m
    attr_reader   :dot_in
    
    def initialize(system, **options)
      @system = system
      @spaces = options.fetch(:spaces, [@system.reference_space]*@system.dimension)
      @m      = options.fetch(:m, 1)
      @dot      = options.fetch(:dot_in, 0)
      @precision = options.fetch(:precision, 8)
      raise "Invalid space dimension #{@spaces.length}!" if @spaces.length != @system.dimension
      @data_space = DataSpace::new(*@system.shapes(@spaces), random: options[:random], m: @m, precision: @precision )
    end

    def [](number)
      raise "Invalid wave function number #{number}!" if number < 0 || number >= @m
      nw = self.class::new(@system, spaces: @spaces)
      mask = [true]*@data_space.data.dimension
      mask[-1] = number if @m > 1
      nw.data_space.data[] = @data_space.data[*mask]
      return nw
    end

    def to(target_spaces)
      target_spaces = [target_spaces].flatten
      target_spaces *= @system.dimension if target_spaces.length == 1
      raise "Invalid space dimension #{target_spaces.length}!" if target_spaces.length != @system.dimension
      source_spaces = @spaces.dup
      source = self
      target_spaces.each_with_index { |s,i|
        next if s == source_spaces[i]
        target_spaces = source_spaces.dup
        target_spaces[i] = s
        target = self.class::new(@system, spaces: target_spaces, m: @m, precision: @precision)
        op_opt = { wavelet_family: @system.wavelet_family, precision: @precision }
        operator = @system.spaces[source_spaces[i]].transition_operator(s)[op_opt]
        bc = @system.bc_from_transition(i, source_spaces[i], s)
        operator.run( i, bc, source, target, narr: @m, dot_in:@dot )
        source_spaces = target_spaces
        source = target
      }
      return source
    end

    def leading_dimensions
      return NArray[*@system.leading_dimensions(@spaces)]
    end

    def shape
      return NArray[*@system.shapes(@spaces)]
    end

    def restricted_data
      return data_space.data[*restricted_slice]
    end

    def restricted_data=(other)
      return data_space.data[*restricted_slice] = other
    end

    def restricted_shape
      return @system.data_shapes(@spaces)
    end

    def restricted_slice
      data_shapes = @system.data_shapes(@spaces)
      ranges = data_shapes.collect { |s|
        0...s
      }
      ranges.push(true) if @m > 1
      return ranges
    end

    def dimensions
      return NArray[*@system.dimensions(@spaces)]
    end

  end

end
