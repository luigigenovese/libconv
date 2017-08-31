class TestSystem < Minitest::Test

  def test_system_new
    s = System::new([42, 28, 23], [BC::Free]*3, :s1)
    assert_equal([42, 28, 23], s.reference_dimensions)
    assert_equal([BC::Free, BC::Free, BC::Free], s.boundary_conditions)
    assert_equal(3, s.dimension)
  end

  def test_dimensions
    s = System::new([42, 28, 23], [BC::Free]*3, :s1)
    assert_equal([42, 28, 23], s.dimensions(:s1))
    assert_equal([42 + 7, 28 + 7, 23 + 7], s.dimensions(:s0))
    assert_equal([(42 + 15)*2, (28 + 15)*2, (23 + 15)*2], s.dimensions(:r))
    assert_equal([42, 28 + 7, (23 + 15)*2], s.dimensions([:s1, :s0, :r]))
  end

  def test_shapes
    s = System::new([42, 28, 23], [BC::Free]*3, :s1)
    assert_equal([42, 2, 28, 2, 23 + 1, 2], s.shapes(:s1))
    assert_equal([2, 42 + 7 + 1, 2, 28 + 7 + 1, 2, 23 + 7], s.shapes(:s0))
    assert_equal([(42 + 15)*2 + 2, (28 + 15)*2 + 2, (23 + 15)*2], s.shapes(:r))
    assert_equal([42, 2, 2, 28 + 7 + 1, (23 + 15)*2], s.shapes([:s1, :s0, :r]))
  end

  def test_leading_dimensions
    s = System::new([42, 28, 23], [BC::Free]*3, :s1)
    assert_equal([42, 28, 23 + 1], s.leading_dimensions(:s1))
    assert_equal([42 + 7 + 1, 28 + 7 + 1, 23 + 7], s.leading_dimensions(:s0))
    assert_equal([(42 + 15)*2 + 2, (28 + 15)*2 + 2, (23 + 15)*2], s.leading_dimensions(:r))
    assert_equal([42, 28 + 7 + 1, (23 + 15)*2], s.leading_dimensions([:s1, :s0, :r]))
  end

end
