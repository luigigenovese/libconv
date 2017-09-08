class TestWaveFunction < Minitest::Test
  def setup
    @s = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s1)
  end

  def test_new
    w = WaveFunction::new(@s)
    assert_equal( [42, 28, 24], w.dimensions.to_a )
    assert_equal( [21 + 1, 2, 14, 2, 12, 2], w.data_space.data.shape )
    assert_equal( [21, 2, 14, 2, 12, 2], w.restricted_data.shape )
    assert_equal( 0.0, w.data_space.data[0] )
  end

  def test_new_random
    w = WaveFunction::new(@s, random: true)
    refute_equal( 0.0, w.data_space.data[0] )
  end

  def test_new_spaces
    w = WaveFunction::new(@s, spaces: [:s0, :r, :s1])
    assert_equal( [42 + 14, 28, 24], w.dimensions.to_a )

    assert_equal( [2, (42 + 14)/2, 28, 12, 2], w.data_space.data.shape )
    assert_equal( [2*(42 + 14)/2, 28, 12*2], w.leading_dimensions.to_a )
  end

  def test_to_s1s0
    w = WaveFunction::new(@s, random: true)
    w2 = w.to([:s0, :s1, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert( (w.restricted_data - w3.restricted_data).abs.max < 10e-16 )

    w2 = w.to([:s1, :s0, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert( (w.restricted_data - w3.restricted_data).abs.max < 10e-16 )


    w2 = w.to([:s1, :s1, :s0])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert( (w.restricted_data - w3.restricted_data).abs.max < 10e-16 )

    w2 = w.to([:s0, :s0, :s0])
    w2p = w.to([:s0, :s1, :s1])
    w3p = w2p.to([:s0, :s0, :s1])
    w4p = w3p.to([:s0, :s0, :s0])
    assert( (w2.restricted_data - w4p.restricted_data).abs.max < 10e-16 )
    w3 = w2.to([:s1, :s1, :s1])
    assert( (w.restricted_data - w3.restricted_data).abs.max < 10e-16 )
  end

end
