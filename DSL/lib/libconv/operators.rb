require 'rgl/implicit'
module LibConv

  class OpBase

    def name
      return self.class.name.split("::").last
    end

    def to_s
      return name
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
      if right.kind_of?(OpBase)
        @left = left
        @right = right
      else
        @right, @left = left.coerce(right)
      end
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

  class Any < OpBase
  end

  class Variable < OpBase
    attr_reader :name

    def initialize(name)
      @name = name
    end

    def to_s
      return name.to_s
    end

    def coerce(other)
      if other.kind_of?(Fixnum)
        [Constant::new(other), self]
      else
        super
      end
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

  class Constant < Variable
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
      def g.vertex_id(v)
        v.object_id
      end
    }
  end

  class TreeExp
    attr_reader :operator
    attr_reader :children

    def initialize( operator, *children )
      @operator = operator
      @children = children
    end

    def name
      return @operator.name
    end

    def to_s
      return name
    end

    def each(&block)
      if block_given?
        yield self
        @children.each { |c|
          if c.kind_of? TreeExp
            c.each(&block)
          else
            yield c
          end
        }
        self
      else
        to_enum(:each)
      end
    end

    def each_child(&block)
      if block_given?
        @children.each(&block)
        self
      else
        to_enum(:each_child)
      end
    end

  end

  module TreeArithmetic
    def **(other)
      TreeExp::new(Exp, self, other)
    end

    def *(other)
      TreeExp::new(Exp, self, other)
    end

    def +(other)
      TreeExp::new(Add, self, other)
    end
  end

  class OpBase
    class << self
      include TreeArithmetic
    end

    def match(tree_exp)
      return true if tree_exp == Any
      if tree_exp.kind_of? TreeExp
        return false unless self.kind_of? tree_exp.operator
        kids = each_child.to_a
        return false if kids.length != tree_exp.children.length
        kids.each_with_index { |c, i|
          return false unless c.match(tree_exp.children[i])
        }
      else
        return false unless (self.kind_of? tree_exp)
      end
      return true
    end

    def scan(tree_exp)
      matches = []
      matches.push self if match(tree_exp)
      each_child { |c|
        matches += c.scan(tree_exp)
      }
      return matches
    end
  end

  def self.subtree_graph(exp)
    RGL::ImplicitGraph.new { |g|
      g.directed = true
      g.vertex_iterator { |b|
        exp.each(&b)
      }
      g.adjacent_iterator { |x, b|
        x.each_child(&b) if x.kind_of? TreeExp
      }
      def g.vertex_id(v)
        v.object_id
      end
    }
  end

end
