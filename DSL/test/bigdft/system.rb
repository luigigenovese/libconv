class TestSystem < Minitest::Test

  def setup
    @s1 = System::new([42, 28, 24], [BC::Free]*3, :s1)
    @s2 = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    @s3 = System::new([99, 56, 61], [BC::Free, BC::Per, BC::Free], :r)
  end

  def test_system_new
    s = System::new([42, 28, 24], [BC::Free]*3, :s1)
    assert_equal([42, 28, 24], s.reference_dimensions)
    assert_equal([BC::Free, BC::Free, BC::Free], s.boundary_conditions)
    assert_equal(3, s.dimension)

    s = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    assert_equal([42, 28, 24], s.reference_dimensions)
    assert_equal([BC::Free, BC::Per, BC::Free], s.boundary_conditions)
    assert_equal(3, s.dimension)

    s = System::new([99, 56, 61], [BC::Free, BC::Per, BC::Free], :r)
    assert_equal([99, 56, 61], s.reference_dimensions)
    assert_equal([BC::Free, BC::Per, BC::Free], s.boundary_conditions)
    assert_equal(3, s.dimension)
  end

  def test_dimensions
    assert_equal([42, 28, 24], @s1.dimensions(:s1))
    assert_equal([42 + 14, 28 + 14, 24 + 14], @s1.dimensions(:s0))
    assert_equal([42 + 30, 28 + 30, 24 + 30], @s1.dimensions(:r))
    assert_equal([42, 28 + 14, 24 + 30], @s1.dimensions([:s1, :s0, :r]))

    assert_equal([42, 28, 24], @s2.dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14], @s2.dimensions(:s1))
    assert_equal([42 + 15, 28, 24 + 15], @s2.dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15], @s2.dimensions([:s1, :s0, :r]))

    assert_equal([99, 56, 61], @s3.dimensions(:r))
    assert_equal([99-15, 56, 61-15], @s3.dimensions(:s0))
    #assert_equal([(99-15)/2, 56/2, (61-15)/2], @s3.dimensions(:s1))
  end

  def test_shapes
    assert_equal([21 + 1, 2, 14, 2, 12, 2], @s1.shapes(:s1))
    assert_equal([2, 21 + 7, 2, 14 + 7 + 1, 2, 12 + 7 + 1], @s1.shapes(:s0))
    assert_equal([42 + 30, 28 + 30 + 2, 24 + 30 + 2], @s1.shapes(:r))
    assert_equal([21 + 1, 2, 2, 14 + 7 + 1, 24 + 30 + 2], @s1.shapes([:s1, :s0, :r]))

    assert_equal([2, 22, 2, 14, 2, 12], @s2.shapes(:s0))
    assert_equal([21 - 7, 2, 14, 2, 12 - 7 + 1, 2], @s2.shapes(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], @s2.shapes(:r))
    assert_equal([21 - 7, 2, 2, 14, 24 + 15 + 1], @s2.shapes([:s1, :s0, :r]))
  end

  def test_leading_dimensions
    assert_equal([42 + 2, 28, 24], @s1.leading_dimensions(:s1))
    assert_equal([42 + 14, 28 + 14 + 2, 24 + 14 + 2], @s1.leading_dimensions(:s0))
    assert_equal([42 + 30, 28 + 30 + 2, 24 + 30 + 2], @s1.leading_dimensions(:r))
    assert_equal([42 + 2, 28 + 14 + 2, 24 + 30 +2], @s1.leading_dimensions([:s1, :s0, :r]))

    assert_equal([42 + 2, 28, 24], @s2.leading_dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14 + 2], @s2.leading_dimensions(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], @s2.leading_dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15 + 1], @s2.leading_dimensions([:s1, :s0, :r]))
  end

end
