[ '../lib', 'lib' ].each { |d| $:.unshift(d) if File::directory?(d) }
require 'bigdft/operators'
include BigDFT
require 'minitest/autorun'
require 'rgl/dot'

class TestArithmetic < Minitest::Test
  def test_expression
    a = Variable::new(:a)
    b = Variable::new(:b)
    c = Variable::new(:c)

    exp = a*MF()**W()**c + b

    g = BigDFT.expression_graph(exp)
    g.write_to_graphic_file('svg', 'toto')

    vertexes = g.each.select { |v|
      v.kind_of?(MagicFilter) && !(g.each_adjacent(v).select { |c|
        c.kind_of?(WaveletTransform)
      }.empty?)
    }

    p vertexes

    vertexes = g.each.select { |v|
      v.kind_of?(Exp) && v.left.kind_of?(MagicFilter) && v.right
    }

  end
end
