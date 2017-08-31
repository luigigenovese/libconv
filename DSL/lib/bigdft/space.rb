module BigDFT

  module Space
    class << self
      attr_reader :space
      def transition_operator(space)
        raise "Unknown transition to #{space.inspect} from #{self.space.inspect}!"
      end
    end
  end

  module S1
    extend Space

    class << self
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
    @s1s0 = IDWT
    @s1r = S1TOR

  end

  module S0
    extend Space

    class << self
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

    @s0r = MF
    @s0s1 = DWT

  end

  module R
    extend Space

    class << self
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

    @rs0 = IMF
    @rs1 = RTOS1

  end

end
