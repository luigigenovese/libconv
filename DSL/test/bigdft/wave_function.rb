class TestWaveFunction < Minitest::Test

  def assert_wfn_equal(w1, w2, epsilon)
    assert( (w1.restricted_data - w2.restricted_data).abs.max < epsilon )
  end

  def setup
    seed = srand
    NArray.srand(seed)
    srand(seed)
    @s1 = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s1)
    @s0 = System::new([56, 28, 38], [BC::Free, BC::Per, BC::Free], :s0)

    @c = System::new([86, 28, 40], [BC::Free, BC::Per, BC::Free], :s0, transitions:  OrderedTransitions::new(:s1, :s0, :r, :shrink))

    @d1 = System::new([42, 28, 40], [BC::Free, BC::Per, BC::Free], :s1, transitions:  OrderedTransitions::new(:s1, :s0, :r, :shrink))
    @d0 = System::new([28, 28, 26], [BC::Free, BC::Per, BC::Free], :s0, transitions:  OrderedTransitions::new(:s1, :s0, :r, :shrink))
  end

  def test_new
    w = WaveFunction::new(@s1)
    assert_equal( [42, 28, 24], w.dimensions.to_a )
    assert_equal( [21 + 1, 2, 14, 2, 12, 2], w.data_space.data.shape )
    assert_equal( [21, 2, 14, 2, 12, 2], w.restricted_data.shape )
    assert_equal( 0.0, w.data_space.data[0] )
  end

  def test_new_random
    w = WaveFunction::new(@s1, random: true)
    refute_equal( 0.0, w.data_space.data[0] )
  end

  def test_new_spaces
    w = WaveFunction::new(@s1, spaces: [:s0, :r, :s1])
    assert_equal( [42 + 14, 28, 24], w.dimensions.to_a )

    assert_equal( [2, (42 + 14)/2, 28, 12, 2], w.data_space.data.shape )
    assert_equal( [2*(42 + 14)/2, 28, 12*2], w.leading_dimensions.to_a )
  end

  def test_to_s1s0
    5.times do |i|
    w = WaveFunction::new(@s1, random: true)
    w2 = w.to([:s0, :s1, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to([:s1, :s0, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert_wfn_equal( w, w3, 10e-16 )


    w2 = w.to([:s1, :s1, :s0])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s1, :s1, :s1])
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to([:s0, :s0, :s0])
    w2p = w.to([:s0, :s1, :s1])
    w3p = w2p.to([:s0, :s0, :s1])
    w4p = w3p.to([:s0, :s0, :s0])
    assert_wfn_equal( w2, w4p, 10e-16 )
    w3 = w2.to([:s1, :s1, :s1])
    assert_wfn_equal( w, w3, 10e-15 )
    end
  end

  def test_to_s0s1
    5.times do |i|
    w = WaveFunction::new(@c, random: true)

    w2 = w.to([:s1, :s0, :s0])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s0, :s0, :s0])
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to([:s0, :s1, :s0])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s0, :s0, :s0])
    assert_wfn_equal( w, w3, 10e-16 )


    w2 = w.to([:s0, :s0, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w3 = w2.to([:s0, :s0, :s0])
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to([:s1, :s1, :s1])
    w2p = w.to([:s1, :s0, :s0])
    w3p = w2p.to([:s1, :s1, :s0])
    w4p = w3p.to([:s1, :s1, :s1])
    assert_wfn_equal( w2, w4p, 10e-16 )
    w3 = w2.to([:s0, :s0, :s0])
    assert_wfn_equal( w, w3, 10e-15 )
    end
  end

  def test_precision
    wd = WaveFunction::new(@s1, random: true)
    ws = WaveFunction::new(@s1, precision: 4)

    ws.restricted_data = wd.restricted_data

    assert_equal( wd.restricted_data.to_ptr.size, ws.restricted_data.to_ptr.size * 2)

    wd2 = wd.to(:r)
    ws2 = ws.to(:r)

    assert_equal( wd2.restricted_data.to_ptr.size, ws2.restricted_data.to_ptr.size * 2)

    assert_wfn_equal( wd2, ws2, 10e-7 )
  end

  def test_to_s1r_per
    w = WaveFunction::new(@s1, random: true)
    w2 = w.to([:s1, :r, :s1])
    w3 = w.to([:s1, :s0, :s1]).to([:s1, :r, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_rs1_per
    w = WaveFunction::new(@s1, spaces: [:s1, :r, :s1], random: true)
    w2 = w.to([:s1, :s1, :s1])
    w3 = w.to([:s1, :s0, :s1]).to([:s1, :s1, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_s1r_grow
    w = WaveFunction::new(@s1, random: true)
    w2 = w.to([:r, :s1, :s1])
    w3 = w.to([:s0, :s1, :s1]).to([:r, :s1, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_rs1_shrink
    w = WaveFunction::new(@s1, spaces: [:r, :s1, :s1], random: true)
    w2 = w.to([:s1, :s1, :s1])
    w3 = w.to([:s0, :s1, :s1]).to([:s1, :s1, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_s1r_shrink
    w = WaveFunction::new(@d1, random: true)
    w2 = w.to([:r, :s1, :s1])
    w3 = w.to([:s0, :s1, :s1]).to([:r, :s1, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_rs1_grow
    w = WaveFunction::new(@d1, spaces: [:r, :s1, :s1], random: true)
    w3 = w.to([:s0, :s1, :s1]).to([:s1, :s1, :s1])
    w2 = w.to([:s1, :s1, :s1])
    assert_wfn_equal( w2, w3, 10e-15 )
  end

  def test_to_s0s1_free
    w = WaveFunction::new(@d1, spaces: [:s0, :s1, :s1], random: true)
    d2 = System::new(w.dimensions.to_a, @d1.boundary_conditions, :s1, transitions:  OrderedTransitions::new(:s1, :s0, :r, :discard))
    wd = WaveFunction::new(d2, spaces: [:s0, :s1, :s1])
    wd.data_space.data[] = w.data_space.data[]

    w2 = w.to(:s1)
    wd2 = wd.to(:s1)

    assert( (w2.restricted_data[3..-5, true, true, true, true, true] - wd2.restricted_data).abs.max < 10e-16 )

    w2.data_space.data[0..2, true, true, true, true, true] = 0.0
    w2.data_space.data[(w2.restricted_shape[0]-4)..-1, true, true, true, true, true] = 0.0

    w3  =  w2.to([:s0, :s1, :s1])
    wd3 = wd2.to([:s0, :s1, :s1])

    assert_wfn_equal( w3, wd3, 10e-15 )
  end

  def test_to_s0r_free
    w = WaveFunction::new(@d1, spaces: [:r, :s1, :s1])
    d2 = System::new((w.dimensions + NArray[1,0,0]).to_a, @d1.boundary_conditions, :s1, transitions:  OrderedTransitions::new(:s1, :s0, :r, :discard))
    wd = WaveFunction::new(d2, spaces: [:r, :s1, :s1])
    w.restricted_data = NArray::float(*w.restricted_shape).random
    wd.data_space.data[*w.restricted_slice] = w.restricted_data

    w2 = w.to([:s0, :s1, :s1])
    wd2 = wd.to([:s0, :s1, :s1])

    assert( (w2.restricted_data[true, 4..-4, true, true, true, true] - wd2.restricted_data[true, true, true, true, true, true]).abs.max < 10e-16 )

    w2.data_space.data[true, 0..3, true, true, true, true] = 0.0
    w2.data_space.data[true, (w2.restricted_shape[1]-3)..-1, true, true, true, true] = 0.0

    w3 = w2.to([:r, :s1, :s1])
    wd3 = wd2.to([:r, :s1, :s1])

    assert( (w3.restricted_data - wd3.restricted_data[0..-2, true, true, true, true]).abs.max < 10e-16 )

  end

  def test_to_s1r
    w = WaveFunction::new(@s1, random: true)
    w2 = w.to([:r, :s1, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )
    w4 = w.to([:s0, :s1, :s1]).to([:r, :s1, :s1])
    assert_wfn_equal( w4, w2, 10e-15 )

    w4 = w.to(:r)
    w3 = w.to(:s0).to(:r)
    assert_wfn_equal( w4, w3, 10e-15 )

    w4 = w4.to(:s1)
    w3 = w3.to(:s0).to(:s1)
    assert_wfn_equal( w4, w3, 10e-15 )

    wp = WaveFunction::new(@s0, spaces: [:s1, :s1, :s1])
    wp.data_space.data[] = w.data_space.data[]
    w3 = wp.to([:s0, :s1, :s1]).to([:r, :s1, :s1])
    assert_wfn_equal( w3, w2, 10e-15 )

    w = WaveFunction::new(@d1, random: true)
    w2 = w.to([:r, :s1, :s1])
    refute_equal( w.data_space.data, w2.data_space.data )
    assert_equal( w2.shape.to_a, w2.data_space.data.shape )

    wp = WaveFunction::new(@d0, spaces: [:s1, :s1, :s1])
    wp.data_space.data[] = w.data_space.data[]
    w3 = wp.to([:s0, :s1, :s1]).to([:r, :s1, :s1])
    assert_equal( w3.shape.to_a, w3.data_space.data.shape )
    assert_wfn_equal( w3, w2, 10e-15 )
    
  end

  def test_wrong_dimensions
    w = WaveFunction::new(@s1, random: true)
    assert_raises(RuntimeError, "Invalid space dimension 2!") { w2 = w.to([:s0, :s0]) }
  end

  def test_1d
    s1 = System::new([42], [BC::Free], :s1)
    w = WaveFunction::new(s1, random: true)
    w2 = w.to([:s0])
    w3 = w2.to([:s1])
    assert( (w.restricted_data - w3.restricted_data).abs.max < 10e-16 )
  end

  def test_m
    s1 = System::new([42, 28, 24], [BC::Free, BC::Per, BC::Free], :s1)
    w = WaveFunction::new(s1, random: true, m: 3)
    assert_equal( w.data_space.data.shape, w.shape.to_a + [3] )
    w2 = w.to([:s0, :s1, :s1])
    w3 = w2.to([:s1, :s1, :s1])
    assert_equal( w.restricted_data.shape, w.restricted_shape.to_a + [3] )
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to([:s0, :s0, :s0])
    w3 = w2.to([:s1, :s1, :s1])
    assert_wfn_equal( w, w3, 10e-15 )

    ws = w[1]
    ws3 = ws.to([:s0, :s0, :s0]).to([:s1, :s1, :s1])
    assert_wfn_equal( ws3, w3[1], 10e-15 )

    ws3 = ws.to(:r).to(:s1)
    w3 = w.to(:r).to(:s1)
    assert_wfn_equal( ws3, w3[1], 10e-15 )
  end

  def test_small
    s1 = System::new([4,6,2], [BC::Free, BC::Per, BC::Free], :s1)
    w = WaveFunction::new(s1, random: true)
    w2 = w.to(:s0)
    w3 = w2.to(:s1)
    assert_wfn_equal( w, w3, 10e-16 )

    w2 = w.to(:r)
    w3 = w.to(:s0).to(:r)
    assert_wfn_equal( w2, w3, 10e-15 )

    w2 = w2.to(:s1)
    w3 = w3.to(:s0).to(:s1)
    assert_wfn_equal( w2, w3, 10e-15 )
  end

end
