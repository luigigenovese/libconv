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

    exp = a*MF()**W()**(c*2) + 2*b
    subtree = MF**W**Any

    BigDFT.expression_graph(exp).write_to_graphic_file('svg', 'toto')

    BigDFT.subtree_graph(subtree).write_to_graphic_file('svg', 'tata')

    res = exp.scan(subtree)
    BigDFT.expression_graph(res.first).write_to_graphic_file('svg', 'tutu')
  end
end
