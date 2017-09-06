module BigDFT

  class Dimension
    attr_reader :dim
    attr_reader :space
    
    def initialize( dim, space )
      @dim = dim
      raise "Unknown space #{space.inspect}!" unless [:s1, :s0, :r].include?(space)
      @space = space
    end

    def r_dim
      if @space == :r then
        @dim
      else
      	@dim * 2
      end
    end

    def r_ld(padding = 0)
      pad(r_dim, padding)
    end

    def r_shape(padding = 0)
      [ r_ld(padding) ]
    end

    def s1_dim
      if @space == :r then
        raise "Invalid dimension #{@dim} for s1" if @dim%2 != 0
        @dim / 2
      else
        @dim
      end
    end

    def s1_ld(padding = 0)
      pad(s1_dim, padding)
    end

    def s1_shape(padding = 0)
      [ s1_ld(padding)/2, 2]
    end

    alias s0_dim s1_dim

    def s0_ld(padding = 0)
      pad(s0_dim, padding)
    end

    def s0_shape(padding = 0)
      [ 2, s0_ld(padding)/2 ]
    end

    def get_dim( space )
      send("#{space}_dim")
    end

    def get_ld(space, padding = 0)
      send("#{space}_ld", padding)
    end

    def get_shape(space, padding = 0)
      send("#{space}_shape", padding)
    end

    private

    def pad(d, padding)
      return d if padding == 0
      remainder = d % padding
      d += padding - remainder if remainder > 0
      d
    end

  end

end
