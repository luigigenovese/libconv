require 'rgl/implicit'
module BigDFT

  class OpBase

    def name
      return self.class.name.split("::").last
    end

    def to_s
      return name << ": #{object_id}"
    end

  end

  class Operator < OpBase

    def each_child
      if block_given?
        self
      else
        to_enum(:each_child)
      end
    end

    def each
      if block_given?
        yield self
        self
      else
        to_enum(:each)
      end
    end

  end

  class Wavelet < Operator
  end
  W = Wavelet

  def self.W(*args, &block)
    W::new *args, &block
  end

  def W(*args, &block)
    W::new *args, &block
  end

  class InverseWavelet < Operator
  end
  IW = InverseWavelet

  def self.IW(*args, &block)
    IW::new *args, &block
  end

  def IW(*args, &block)
    IW::new *args, &block
  end

  class MagicFilter < Operator
  end
  MF = MagicFilter

  def self.MF(*args, &block)
    MF::new *args, &block
  end

  def MF(*args, &block)
    MF::new *args, &block
  end

  class InverseMagicFilter < Operator
  end
  IMF = InverseMagicFilter

  def self.IMF(*args, &block)
    IMF::new *args, &block
  end

  def IMF(*args, &block)
    IMF::new *args, &block
  end

  class Operation < OpBase

    attr_reader :left, :right

    def initialize(left, right)
      @left = left
      @right = right
    end

    def each_child
      if block_given?
        yield @left
        yield @right
        self
      else
        to_enum(:each_child)
      end
    end

    def each(&block)
      if block
        block.call self
        @left.each(&block)
        @right.each(&block)
        self
      else
        to_enum(:each)
      end
    end

  end

  class Mul < Operation
  end

  class Add < Operation
  end

  class Exp < Operation
  end

  module Arithmetic
    def **(other)
      Exp::new(self, other)
    end

    def *(other)
      Mul::new(self, other)
    end

    def +(other)
      Add::new(self, other)
    end
  end

  class OpBase
    include Arithmetic
  end

  class Variable
    include Arithmetic
    attr_reader :name

    def initialize(name)
      @name = name
    end

    def to_s
      return name.to_s << ": #{object_id}"
    end

    def each_child
      to_enum(:each_child) unless block_given?
      self
    end

    def each
      if block_given?
        yield self
        self
      else
        to_enum(:each)
      end
    end
  end

  def self.expression_graph(exp)
    RGL::ImplicitGraph.new { |g|
      g.directed = true
      g.vertex_iterator { |b|
        exp.each(&b)
      }
      g.adjacent_iterator { |x, b|
        x.each_child(&b)
      }
    }
  end

end
