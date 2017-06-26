module BigDFT

  module Arithmetic

    def coerce(other)
      return [Constant::new(other), self] if other.kind_of?(Numeric)
      return [other, self]
    end

    def +( expr )
      expr = Constant::new(expr) if expr.kind_of?(Numeric)
      return Expressions::Addition::new(self, expr)
    end

    def *( expr )
      expr = Constant::new(expr) if expr.kind_of?(Numeric)
      return Expressions::Product::new(self, expr)
    end

    def to_real
      return Expressions::InverseWaveletTransform::new(self)
    end

    def to_wavelet
      return Expressions::WaveletTransform::new(self)
    end

  end

  module DataSpace
    attr_reader :dimension
    attr_reader :space
    attr_reader :boundary_conditions
    def assert_space(space)
      raise "Wrong space for #{self.class}, wanted #{space.inspect}, I am #{@space.inspect}!" unless space == @space
    end
  end

  module Expressions

    class BaseExpression
      include DataSpace
      include Arithmetic
    end

    class ComplexExpression < BaseExpression
    end

    class S0S0 < ComplexExpression
      attr_accessor :a, :ax, :ay
      attr_reader :x
      attr_accessor :y
      def initialize(x)
        x.assert_space(:s0)
        @x = x
        @y = nil
        @a = @ax = @ay = nil
        @space = :s0
      end
    end

    class MagicFilter < ComplexExpression
    end

    class InverseMagicFilter < ComplexExpression
    end

    class SxSy < ComplexExpression
      attr_accessor :a, :ay
      attr_reader :x
      attr_accessor :y
      def initialize(x)
        @x = x
        @y = nil
        @a = @ay = nil
      end
    end

    class WaveletTransform < SxSy
      def initialize(x)
        x.assert_space(:s0)
        super
        @space = :s1
      end
    end

    class InverseWaveletTransform < SxSy
      def initialize(x)
        x.assert_space(:s1)
        super
        @space = :s0
      end
    end

    class Affectation < BaseExpression
      attr_reader :lvalue
      attr_reader :rvalue
      def initialize( lvalue, rvalue )
        @lvalue = lvalue
        @rvalue = rvalue 
      end

    end

    class Addition < BaseExpression
      attr_reader :lmember
      attr_reader :rmember
      def initialize( lmember, rmember)
        @lmember = lmember
        @rmember = rmember
      end
    end

    class Product < BaseExpression
      attr_reader :lmember
      attr_reader :rmember
      def initialize( lmember, rmember)
        @lmember = lmember
        @rmember = rmember
      end
    end

  end

  class Constant
    include Arithmetic
    attr_reader :value
    def initialize(value)
      @value = value
    end
  end

  class WaveFunction
    include Arithmetic
    include DataSpace

    attr_reader :expression
    

    def initialize( dimensions, space )
      @dimensions = dimensions
      @space = space
      @boundary_conditions = nil
      @expression = nil
    end

    def <<(rvalue)
      @expression = Expressions::Affectation::new(self, rvalue)
    end

  end

end
