class TestSystem < Minitest::Test

  def setup
    @s1 = System::new([42, 28, 24], [BC::Free]*3, :s1)
    @s2 = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    @s3 = System::new([99, 56, 61], [BC::Free, BC::Per, BC::Free], :r)
    @s4 = System::new([100, 28, 54], [BC::Free, BC::Per, BC::Free], :s1, transitions:  OrderedTransitions::new(:s1, :s0, :r, :shrink))
  end

  def test_system_new
    assert_equal([42, 28, 24], @s1.reference_dimensions)
    assert_equal([BC::Free, BC::Free, BC::Free], @s1.boundary_conditions)
    assert_equal(3, @s1.dimension)

    assert_equal([42, 28, 24], @s2.reference_dimensions)
    assert_equal([BC::Free, BC::Per, BC::Free], @s2.boundary_conditions)
    assert_equal(3, @s3.dimension)

    assert_equal([99, 56, 61], @s3.reference_dimensions)
    assert_equal([BC::Free, BC::Per, BC::Free], @s3.boundary_conditions)
    assert_equal(3, @s3.dimension)

    assert_equal([100, 28, 54], @s4.reference_dimensions)
    assert_equal([BC::Free, BC::Per, BC::Free], @s4.boundary_conditions)
    assert_equal(3, @s4.dimension)
  end

  def test_dimensions
    assert_equal([42, 28, 24], @s1.dimensions(:s1))
    assert_equal([42 + 14, 28 + 14, 24 + 14], @s1.dimensions(:s0))
    assert_equal([42 + 29, 28 + 29, 24 + 29], @s1.dimensions(:r))
    assert_equal([42, 28 + 14, 24 + 29], @s1.dimensions([:s1, :s0, :r]))

    assert_equal([42, 28, 24], @s2.dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14], @s2.dimensions(:s1))
    assert_equal([42 + 15, 28, 24 + 15], @s2.dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15], @s2.dimensions([:s1, :s0, :r]))

    assert_equal([99, 56, 61], @s3.dimensions(:r))
    assert_equal([99-15, 56, 61-15], @s3.dimensions(:s0))
    assert_equal([99-29, 56, 61-29], @s3.dimensions(:s1))

    assert_equal([100, 28, 54], @s4.dimensions(:s1))
    assert_equal([100- 14, 28, 54 - 14], @s4.dimensions(:s0))
    assert_equal([100- 29, 28, 54 - 29], @s4.dimensions(:r))
    assert_equal([100, 28, 54 - 29], @s4.dimensions([:s1, :s0, :r]))
  end

  def test_shapes
    assert_equal([21 + 1, 2, 14, 2, 12, 2], @s1.shapes(:s1))
    assert_equal([2, 21 + 7, 2, 14 + 7 + 1, 2, 12 + 7 + 1], @s1.shapes(:s0))
    assert_equal([42 + 29 + 1, 28 + 29 + 3, 24 + 29 + 3], @s1.shapes(:r))
    assert_equal([21 + 1, 2, 2, 14 + 7 + 1, 24 + 29 + 3], @s1.shapes([:s1, :s0, :r]))

    assert_equal([2, 22, 2, 14, 2, 12], @s2.shapes(:s0))
    assert_equal([21 - 7, 2, 14, 2, 12 - 7 + 1, 2], @s2.shapes(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], @s2.shapes(:r))
    assert_equal([21 - 7, 2, 2, 14, 24 + 15 + 1], @s2.shapes([:s1, :s0, :r]))

    assert_equal([50, 2, 14, 2, 27 + 1, 2], @s4.shapes(:s1))
    assert_equal([2, 50 - 7 + 1, 2, 14, 2, 27 - 7], @s4.shapes(:s0))
    assert_equal([100 - 29 + 1, 28, 54 - 29 + 3], @s4.shapes(:r))
    assert_equal([50, 2, 2, 14, 54 - 29 + 3], @s4.shapes([:s1, :s0, :r]))
  end

  def test_leading_dimensions
    assert_equal([42 + 2, 28, 24], @s1.leading_dimensions(:s1))
    assert_equal([42 + 14, 28 + 14 + 2, 24 + 14 + 2], @s1.leading_dimensions(:s0))
    assert_equal([42 + 29 + 1, 28 + 29 + 3, 24 + 29 + 3], @s1.leading_dimensions(:r))
    assert_equal([42 + 2, 28 + 14 + 2, 24 + 29 + 3], @s1.leading_dimensions([:s1, :s0, :r]))

    assert_equal([42 + 2, 28, 24], @s2.leading_dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14 + 2], @s2.leading_dimensions(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], @s2.leading_dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15 + 1], @s2.leading_dimensions([:s1, :s0, :r]))

    assert_equal([100, 28, 54 + 2], @s4.leading_dimensions(:s1))
    assert_equal([100 - 14 + 2, 28, 54 - 14], @s4.leading_dimensions(:s0))
    assert_equal([100 - 29 + 1, 28, 54 - 29 + 3], @s4.leading_dimensions(:r))
    assert_equal([100, 28, 54 - 29 + 3], @s4.leading_dimensions([:s1, :s0, :r]))
  end

  def test_1d
    s1 = System::new([42], [BC::Free], :s1)
    assert_equal([42], s1.reference_dimensions)
    assert_equal([BC::Free], s1.boundary_conditions)
    assert_equal(1, s1.dimension)

    assert_equal([42], s1.dimensions(:s1))
    assert_equal([42 + 14], s1.dimensions(:s0))
    assert_equal([42 + 29], s1.dimensions(:r))

    assert_equal([21 + 1, 2], s1.shapes(:s1))
    assert_equal([2, 21 + 7], s1.shapes(:s0))
    assert_equal([42 + 29 + 1], s1.shapes(:r))

    assert_equal([42 + 2], s1.leading_dimensions(:s1))
    assert_equal([42 + 14], s1.leading_dimensions(:s0))
    assert_equal([42 + 29 + 1], s1.leading_dimensions(:r))

  end

end
