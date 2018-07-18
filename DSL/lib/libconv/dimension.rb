module LibConv

  class Dimension
    attr_reader :dim
    attr_reader :space
    
    def initialize( dim, space )
      @dim = dim
      raise "Unknown space #{space.inspect}!" unless [:s1, :s0, :r, :d1, :d2].include?(space)
      @space = space
    end

    def r_dim
      @dim
    end

    def r_ld(padding = 0)
      pad(r_dim, padding)
    end

    def r_shape(padding = 0)
      [ r_ld(padding) ]
    end

    def r_data_shape
      [ r_dim ]
    end

    def s1_dim
      raise "Invalid dimension #{@dim} for s1" if @dim%2 != 0
      @dim
    end

    def s1_ld(padding = 0)
      pad(s1_dim, padding)
    end

    def s1_shape(padding = 0)
      [ s1_ld(padding)/2, 2]
    end

    def s1_data_shape
      [ s1_dim/2, 2 ]
    end

    alias s0_dim s1_dim

    def s0_ld(padding = 0)
      pad(s0_dim, padding)
    end

    def s0_shape(padding = 0)
      [ 2, s0_ld(padding)/2 ]
    end

    def s0_data_shape
      [ 2, s0_dim/2 ]
    end
    
    alias d1_dim s1_dim
    alias d1_shape s0_shape
    alias d1_ld s0_ld
    alias d1_data_shape s0_data_shape
    alias d2_dim s1_dim
    alias d2_shape s0_shape
    alias d2_ld s0_ld
    alias d2_data_shape s0_data_shape
    
    def get_dim( space )
      send("#{space}_dim")
    end

    def get_ld(space, padding = 0)
      send("#{space}_ld", padding)
    end

    def get_shape(space, padding = 0)
      send("#{space}_shape", padding)
    end

    def get_data_shape(space)
      send("#{space}_data_shape")
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
