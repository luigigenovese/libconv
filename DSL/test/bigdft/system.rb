class TestSystem < Minitest::Test

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
    s = System::new([42, 28, 24], [BC::Free]*3, :s1)
    assert_equal([42, 28, 24], s.dimensions(:s1))
    assert_equal([42 + 14, 28 + 14, 24 + 14], s.dimensions(:s0))
    assert_equal([42 + 30, 28 + 30, 24 + 30], s.dimensions(:r))
    assert_equal([42, 28 + 14, 24 + 30], s.dimensions([:s1, :s0, :r]))

    s = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    assert_equal([42, 28, 24], s.dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14], s.dimensions(:s1))
    assert_equal([42 + 15, 28, 24 + 15], s.dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15], s.dimensions([:s1, :s0, :r]))

    s = System::new([99, 56, 61], [BC::Free, BC::Per, BC::Free], :r)
    assert_equal([99, 56, 61], s.dimensions(:r))
    assert_equal([99-15, 56, 61-15], s.dimensions(:s0))
    #assert_equal([(99-15)/2, 56/2, (61-15)/2], s.dimensions(:s1))
  end

  def test_shapes
    s = System::new([42, 28, 24], [BC::Free]*3, :s1)
    assert_equal([21 + 1, 2, 14, 2, 12, 2], s.shapes(:s1))
    assert_equal([2, 21 + 7, 2, 14 + 7 + 1, 2, 12 + 7 + 1], s.shapes(:s0))
    assert_equal([42 + 30, 28 + 30 + 2, 24 + 30 + 2], s.shapes(:r))
    assert_equal([21 + 1, 2, 2, 14 + 7 + 1, 24 + 30 + 2], s.shapes([:s1, :s0, :r]))

    s = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    assert_equal([2, 22, 2, 14, 2, 12], s.shapes(:s0))
    assert_equal([21 - 7, 2, 14, 2, 12 - 7 + 1, 2], s.shapes(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], s.shapes(:r))
    assert_equal([21 - 7, 2, 2, 14, 24 + 15 + 1], s.shapes([:s1, :s0, :r]))
  end

  def test_leading_dimensions
    s = System::new([42, 28, 24], [BC::Free]*3, :s1)
    assert_equal([42 + 2, 28, 24], s.leading_dimensions(:s1))
    assert_equal([42 + 14, 28 + 14 + 2, 24 + 14 + 2], s.leading_dimensions(:s0))
    assert_equal([42 + 30, 28 + 30 + 2, 24 + 30 + 2], s.leading_dimensions(:r))
    assert_equal([42 + 2, 28 + 14 + 2, 24 + 30 +2], s.leading_dimensions([:s1, :s0, :r]))

    s = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s0)
    assert_equal([42 + 2, 28, 24], s.leading_dimensions(:s0))
    assert_equal([42 - 14, 28, 24 - 14 + 2], s.leading_dimensions(:s1))
    assert_equal([42 + 15 + 3, 28, 24 + 15 + 1], s.leading_dimensions(:r))
    assert_equal([42 - 14, 28, 24 + 15 + 1], s.leading_dimensions([:s1, :s0, :r]))
  end

end
