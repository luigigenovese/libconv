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

    exp = a*b*c + c

    g = BigDFT::expression_graph(exp)
    g.write_to_graphic_file('svg')

  end
end
