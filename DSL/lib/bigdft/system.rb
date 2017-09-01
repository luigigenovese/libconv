module BigDFT

  class OrderedTransitions

    def self.inverse( tr )
      if tr == :shrink
        :grow
      elsif tr == :grow
        :shrink
      else
        tr
      end
    end

    def initialize(*args)
      @tr = args.last
      @states = args[0..-2]
      @transitions = Hash::new { |hash, key| hash[key] = {} }
      @states[0..-2].each_with_index { |s, i|
        @states[(i+1)..-1].each { |d|
          @transitions[s][d] = @tr
        }
      }
      reverse_tr = OrderedTransitions.inverse(@tr)
      reverse_states = @states.reverse
      reverse_states[0..-2].each_with_index { |s, i|
        reverse_states[(i+1)..-1].each { |d|
          @transitions[s][d] = reverse_tr
        }
      }
    end

    def transition(source, destination)
      return @transitions[source][destination]
    end

  end

  class System
    attr_reader :reference_space
    attr_reader :boundary_conditions
    attr_reader :grow_direction
    attr_reader :dimension
    attr_reader :spaces

    def initialize(reference_dimensions, boundary_conditions, reference_space, options = {})
      @reference_dimensions = reference_dimensions.collect { |d| Dimension::new(d, reference_space) }
      @dimension            = reference_dimensions.length
      @spaces               = options.fetch(:spaces, { s1: S1, s0: S0, r: R })
      @reference_space      = reference_space
      @boundary_conditions  = boundary_conditions
      @padding              = options.fetch(:padding, { s1: 2, s0: 2, r: 4 })
      raise "Boundary conditions and reference dimensions must have the same arity (#{boundary_conditions.length} != #{@dimension})!" if boundary_conditions.length != @dimension
      @transitions          = options.fetch(:transitions, OrderedTransitions::new(:s1, :s0, :r, :grow))
    end

    def reference_dimensions
      @reference_dimensions.collect { |d| d.get_dim(@reference_space) }
    end

    def bc_from_transition(idim, space1, space2)
      if @boundary_conditions[idim] == BC::Per then
        GenericConvolution::BC::PERIODIC
      else
        tr = @transitions.transition(space1, space2)
        if tr == :grow then
          GenericConvolution::BC::GROW
        elsif tr == :shrink
          GenericConvolution::BC::SHRINK
        else
          raise "Unknown transition: #{tr.inspect}!"
        end
      end
    end

    def dimensions(space)
      space = get_space(space)

      @reference_dimensions.each_with_index.collect do |d,i|
        if space[i] == @reference_space
          d.get_dim(@reference_space)
        else
          op = @spaces[@reference_space].transition_operator(space[i])
          Dimension::new(op.dims_from_in(d.get_dim(op.dimension_space), bc_from_transition(i, @reference_space, space[i])).last, op.dimension_space).get_dim(space[i])
        end
      end
    end

    def shapes(space)
      space = get_space(space)
      dims = dimensions(space)
      dims.each_with_index.collect { |d,i|
        Dimension::new(d, space[i]).get_shape(space[i], @padding[space[i]])
      }.reduce([], :+)
    end


    def leading_dimensions(space)
      space = get_space(space)
      dims = dimensions(space)
      dims.each_with_index.collect { |d,i|
        Dimension::new(d, space[i]).get_ld(space[i], @padding[space[i]])
      }
    end

    private

    def get_space(space)
      unless space.kind_of?(Array)
        space = [space] * @dimension
      else
        raise "Invalid space dimension #{space.length} ( != #{@dimension} )!" if space.length != @dimension
      end
      return space
    end

  end

end
