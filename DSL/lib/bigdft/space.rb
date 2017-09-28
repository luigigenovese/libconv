module BigDFT

  module Space
    module SpaceMethods
      attr_reader :space
      def transition_operator(space)
        raise "Unknown transition to #{space.inspect} from #{self.space.inspect}!"
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

  end

  module D_S1
    include S1

    @s1s0 = D_IDWT
    @s1r = D_S1TOR

  end

  module S0
    extend Space

    module S0Methods
      attr_reader :s0r
      attr_reader :s0s1

      def transition_operator(space)
        case space
        when :r
          @s0r
        when :s1
          @s0s1
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

  end

  module D_S0
    include S0

    @s0r = D_MF
    @s0s1 = D_DWT

  end

  module R
    extend Space

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

  end

  module D_R
    include R

    @rs0 = D_IMF
    @rs1 = D_RTOS1

  end

end
