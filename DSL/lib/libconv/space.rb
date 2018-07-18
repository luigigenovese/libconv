module LibConv

  module Space
    module SpaceMethods
      attr_reader :space
      def transition_operator(space)
        raise "Unknown transition to #{space.inspect} from #{self.space.inspect}!"
      end

      def bind_ops(op)
        table = {}
        CONFIGURATION.each { |config|
          table[config] = LibConv::const_get(LibConv::const_name(op, config))
        }
        return table
      end
    end
    extend SpaceMethods

    def self.included( other )
      other.extend( SpaceMethods )
    end

  end

  module S1
    include Space

    module S1Methods
      attr_reader :s1s0
      attr_reader :s1r

      def transition_operator(space)
        case space
        when :s0
          @s1s0
        when :r
          @s1r
        else
          super
        end
      end
    end

    extend S1Methods
    @space = :s1

    def self.included( other )
      other.extend( S1Methods )
    end

    @s1s0 = bind_ops("IDWT")
    @s1s0.default = D_SYM8_IDWT
    @s1r = bind_ops("S1TOR")
    @s1r.default = D_SYM8_S1TOR

  end

  module S0
    include Space

    module S0Methods
      attr_reader :s0r
      attr_reader :s0s1
      attr_reader :s0d1
      attr_reader :s0d2

      def transition_operator(space)
        case space
        when :r
          @s0r
        when :s1
          @s0s1
        when :d1
          @s0d1
        when :d2
          @s0d2
        else
          super
        end
      end

    end

    extend S0Methods
    @space = :s0

    def self.included( other )
      other.extend( S0Methods )
    end

    @s0r = bind_ops("MF")
    @s0r.default = D_SYM8_MF
    @s0s1 = bind_ops("DWT")
    @s0s1.default = D_SYM8_DWT
    @s0d1 = bind_ops("D1")
    @s0d1.default = D_SYM8_D1
    @s0d2 = bind_ops("D2")
    @s0d2.default = D_SYM8_MF
    
  end

  module R
    include Space

    module RMethods
      attr_reader :rs0
      attr_reader :rs1

      def transition_operator(space)
        case space
        when :s0
          @rs0
        when :s1
          @rs1
        else
          super
        end
      end

    end

    extend RMethods
    @space = :r

    def self.included( other )
      other.extend( RMethods )
    end

    @rs0 = bind_ops("IMF")
    @rs0.default = D_SYM8_IMF
    @rs1 = bind_ops("RTOS1")
    @rs1.default = D_SYM8_RTOS1

  end

end
