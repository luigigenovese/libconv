module BOAST

  class WaveFunction
    attr_accessor :data_space
    attr_reader   :system
    attr_reader   :spaces

    def initialize(system, **options)
      @system = system
      @spaces = options.fetch(:spaces, [@system.reference_space]*@system.dimension)
      raise "Invalid space dimension #{@spaces.length}!" if @spaces.length != @system.dimension
      @data_space = DataSpace::new(*@system.shapes(@spaces), random: options[:random] )
    end

    def to(*target_spaces)
      target_spaces *= @system.dimension if target_spaces.length == 1
      raise "Invalid space dimension #{target_spaces.length}!" if target_spaces.length != @system.dimension
      source_spaces = @spaces.dup
      source = self
      target_spaces.each_with_index { |s,i|
        next if s == source_spaces[i]
        target_spaces = source_spaces.dup
        target_spaces[i] = s
        target = WaveFunction::new(@system, spaces: target_spaces)
        operator = @system.spaces[source_spaces[i]].transition_operator(s)
        bc = @system.bc_from_transition(i, source_spaces[i], s)
        operator.run( i, bc, source, target )
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
      data_shapes = @system.data_shapes(@spaces)
      ranges = data_shapes.collect { |s|
        0...s
      }
      return data_space.data[*ranges]
    end

    def dimensions
      return NArray[*@system.dimensions(@spaces)]
    end

  end

end
