require 'rgl/implicit'
module BigDFT

  class Operation
    attr_reader :left, :right
    def initialize(left, right)
      @left = left
      @right = right
    end

    def name
      return self.class.name.split("::").last
    end

    def to_s
      return name << ": #{object_id}"
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
        yield self
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

  module Arithmetic
    def *(other)
      Mul::new(self, other)
    end

    def +(other)
      Add::new(self, other)
    end
  end

  class Operation
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

    def each(&block)
      if block
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
        x.each_child { |y|
          b.call(y)
        }
      }
    }
  end

end
