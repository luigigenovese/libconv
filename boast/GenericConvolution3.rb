require 'BOAST'
require 'narray'
include BOAST

register_funccall("modulo")
register_funccall("min")
register_funccall("max")

module GenericConvolution

class Filter
  # List of the floating point values of the convolution filter
  attr_reader :fil_array
  # central point of the filter
  attr_reader :center
  attr_reader :length
  # name of the filter
  attr_reader :name
  attr_reader :base_name

  attr_reader :unroll_inner
  attr_accessor :bc

  def elem_number
    n = (@unroll + @vec_len - 1)/@vec_len
    return n <= 0 ? 1 : n
  end

  def init_optims(tt_arr, mod_arr, unroll_inner, unroll, vec_len)
    @tt_arr = tt_arr
    @unroll_inner = unroll_inner
    @unroll = unroll
    @vec_len = vec_len
    @tt_n = elem_number
    @tt = init_tt
    @filter_val = init_filter_val
    if mod_arr then
      @mods = Int("mod_arr", :allocate => true,
                      :dim => [Dim(get_mods_lower, get_mods_upper)])
    else
      @mods = nil
    end
  end

  def get_mods_lower
    return lowfil - upfil
  end

  def get_mods_upper
    return upfil - lowfil - 1
  end

  def get_mods
    return @mods
  end

  def get_tt
    return @tt
  end

  def get_filter_val
    return @filter_val
  end

  def decl_filter_val
    decl *([@filter_val].flatten)
  end

  def get_input_dim( dim, options = {} )
    if @bc.shrink? then
      if options[:ld] then
        get_input_dim_ld_shrink( dim )
      else
        get_input_dim_shrink( dim )
      end
    else
      get_input_dim_std( dim )
    end
  end

  def get_output_dim( dim, options = {} )
    if @bc.grow? then
      if options[:ld] then
        get_output_dim_ld_grow( dim )
      else
        get_output_dim_grow( dim )
      end
    else
      get_output_dim_std( dim )
    end
  end

  def get_loop_start_end( side, dim, iterator )
    register_funccall("min")
    register_funccall("max")
    loop_start = lowfil
    loop_end   = upfil
    if @bc.free? then
      if side == :begin then
        loop_start = max(-iterator, lowfil)
      elsif side == :end then
        loop_end   = min(upfil, dim - 1 - iterator)
      end
    end
    return [loop_start, loop_end]
  end

  def input_filter_index( iterator, filter_index, side, dim )
    if @bc.free? or side == :center then
      return filter_index + iterator
    else
      if @mods then
        if side == :end then
          return @mods[ filter_index + iterator - dim]
        else
          return @mods[ filter_index + iterator]
        end
      else
        return filter_index + iterator - ((filter_index + iterator + dim * 2 )/dim - 2) * dim
      end
    end
  end

end

class ConvolutionFilter < Filter
  # BOAST variables
  # Filter array (to be used on BOAST functions)
  attr_reader :fil
  # extremes of the filter, calculated via its central point (integers)
  attr_reader :lowfil_val, :upfil_val
  # extremes of the filter, calculated via its central point (BOAST object)
  attr_reader :lowfil, :upfil

  def initialize(name, filt, center)
    @fil_array = filt.dup
    arr = ConstArray::new(@fil_array,Real)
    @fil = Real("#{name}_fil",:constant => arr,:dim => [ Dim((-center),(@fil_array.length - center -1)) ])
    @lowfil_val = -center
    @upfil_val = @fil_array.length - center - 1
    @lowfil = Int("lowfil",:constant => @lowfil_val)
    @upfil = Int("upfil",:constant => @upfil_val)
    @center = center
    @name = name
    @base_name = "s0s0"
    @length = @fil_array.length
  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return Real("tt", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len)
    else
      return (0..(@tt_n-1)).collect { |index| Real("tt#{index}", :vector_length => @vec_len) }
    end
  end

  def set_zero_tt( tt_ind )
    pr @tt[tt_ind].set( 0.0 )
  end

  def init_filter_val
    return @tt[0].copy("filter_val")
  end

  def set_filter_val( index, options = {} )
    pr @filter_val === Set(@fil[index], @tt[0])
  end

  def compute_values( tt_ind, indexes, data_in)
    out = @tt[tt_ind]
    pr out === FMA(Load(data_in[*indexes].set_align(out.type.total_size), out), @filter_val, out) # + @x[*i_in]*@filter.fil[l]
  end

  def post_process_and_store_values( tt_ind, indexes, data_out, transpose, options = {} )
    i_out = indexes.rotate(transpose)
    i_in = options[:indexes_in]
    out = @tt[tt_ind]
    pr out === out * Set(options[:a], out) if options[:a]

    finish_block = lambda {
      pr options[:dot_in_tmp] === FMA(Load(options[:data_in][*i_in].set_align(out.type.total_size), out), out, options[:dot_in_tmp]) if options[:dot_in_tmp] #reduction !
      pr out === FMA(Load(options[:data_in][*i_in], out), Set(options[:a_x], out), out) if options[:a_x]
    }
    if @bc.grow? and (options[:dot_in_tmp] or options[:a_x]) and options[:side] != :center then
      if options[:side] == :begin then
        pr If(options[:position] >= 0, &finish_block)
      elsif options[:side] == :end then
        pr If(options[:position] < options[:dim], &finish_block)
      end
    else
      finish_block.call
    end

    #to be controlled in the case of non-orthorhombic cells for kinetic operations
    pr out === out + Load(options[:data_in2][*i_in], out) if options[:data_in2]

    if options[:accumulate] or (options[:kinetic] == :inplace and not options[:zero_out])  then
      pr out === out + Load(data_out[*i_out].set_align(out.type.total_size), out)
    elsif options[:a_y] then
      pr out === FMA(Set(options[:a_y], out), Load(data_out[*i_out].set_align(out.type.total_size), out), out)
    end
    pr data_out[*i_out].set_align(out.type.total_size) === out
    pr options[:data_out2][*i_out] === Load(options[:data_in][*i_in], out) if options[:data_out2]
  end

  def get_input_dim_std( dim )
    return Dim(0, dim - 1)
  end

  alias get_output_dim_std get_input_dim_std

  def get_input_dim_shrink( dim )
    return Dim( lowfil, dim + upfil  - 1)
  end

  def get_output_dim_grow( dim )
    return Dim(-upfil,  dim - lowfil - 1)
  end

  def get_input_dim_ld_shrink( dim )
    return Dim(lowfil, dim + lowfil - 1)
  end

  def get_output_dim_ld_grow( dim )
    return Dim( -upfil, dim - upfil - 1)
  end

  def decl_filters
    decl fil
  end

  def cost
    return @length * 2
  end

end #class ConvolutionFilter

class PoissonFilter < ConvolutionFilter
  attr_reader :filters_array
  attr_reader :filters

  def initialize(name, filts, nord)
    @center = nord/2
    tmp_array = []
    @filters=[]
    index=0
    filts.each_with_index { |e,i|
      @filters[i]=ConvolutionFilter::new(name+"_#{i}", e, @center)
        e.each{|val|
            tmp_array[index]=val
            index= index+1
        }
    }
    arr = ConstArray::new(tmp_array,Real)
    @filters_array = Real("#{name}_fil",:constant => arr,:dim => [ Dim(0,(2*@center+1)*(2*@center+1)-1) ])

    super(filters[@center].name, @filters[@center].fil_array, @center)
    @name = name

    @fil_array = filts.dup

  end

  def decl_filters
    decl filters_array if @bc.nper?
    decl @filters[@center].fil
  end

  def set_filter_val( index, options = {} )
    if @bc.nper? then
      if options[:side] == :center then
        super
      elsif options[:side] == :begin then
        pr @filter_val === Set(@filters_array[((options[:position] - (upfil)                       )*(2*@center+1)) + (@center)*(2*@center+1) + index + @center], @tt[0])
      elsif options[:side] == :end then
        pr @filter_val === Set(@filters_array[((options[:position] + (@center) - options[:dim] + 1 )*(2*@center+1)) + (@center)*(2*@center+1) + index + @center], @tt[0])
      end
    else
      super
    end
  end

  def input_filter_index( iterator, filter_index, side, dim )
    if @bc.nper? and side != :center then
      if side == :end then
        return filter_index + lowfil + dim - 1
      else
        return filter_index + upfil
      end
    else
      super
    end
  end

end

class WaveletFilter < Filter
  attr_reader :low, :high
  attr_reader :low_even, :low_odd
  attr_reader :high_even, :high_odd
  attr_reader :convolution
  def initialize(name, filt1, opts = {})
    @convolution = opts[:convolution_filter]
    @fil_array = filt1.dup
    if @convolution then
      #initialize low filter with the nonperiodic BC recipe for evaluating length
      @center, @low, @high = self.combine_wavelet_filters(name,@fil_array,@convolution,@convolution.center)
    else
      @center = @fil_array.length/2
      @center -= @center%2
      @low = ConvolutionFilter::new(name+"_l", @fil_array, @center)
      filt2 = []
      @fil_array.each_with_index { |e,i|
        if i % 2 == 0 then
          filt2.push(e)
        else
          if e.is_a?(String) then
            if e[0] == '-' then
              filt2.push(e[1..-1])
            else
              filt2.push("-"+e)
            end
          else
            filt2.push(-e)
          end
        end
      }
      filt2.reverse!
      @high = ConvolutionFilter::new(name+"_h", filt2, @center)

    end

    @length = @low.fil_array.length
    @name = name
  end

  def combine_wavelet_filters(name,fil_array,convolution,convcenter)
    require 'bigdecimal'
    lp = @fil_array.collect { |e| BigDecimal.new(e) }
    hp = lp.each_with_index.collect { |e,i| i % 2 == 0 ? e : -e }.reverse
    cv = @convolution.fil_array.reverse.collect { |e| BigDecimal.new(e) }

    lp_length = lp.length
    hp_length = hp.length
    cv_length = @convolution.length

    lp_center = lp_length / 2
    hp_center = hp_length / 2
    cv_center = convcenter

    cv_bounds = [ -cv_center, cv_length - cv_center - 1 ]
    lp_bounds = [ -lp_center, lp_length - lp_center - 1 ]
    hp_bounds = [ -hp_center, hp_length - hp_center - 1 ]

    lp_cv_bounds = [ cv_bounds[0] - lp_bounds[1], cv_bounds[1] - lp_bounds[0] ]
    hp_cv_bounds = [ cv_bounds[0] - hp_bounds[1], cv_bounds[1] - hp_bounds[0] ]

    lp_cv_bounds[0] -= 1 if lp_cv_bounds[0] % 2 != 0
    hp_cv_bounds[0] -= 1 if hp_cv_bounds[0] % 2 != 0
    lp_cv_bounds[1] += 1 if lp_cv_bounds[1] % 2 == 0
    hp_cv_bounds[1] += 1 if hp_cv_bounds[1] % 2 == 0

    lp_cv_center = - lp_cv_bounds[0]
    hp_cv_center = - hp_cv_bounds[0]

    lp_cv = (lp_cv_bounds[0]..lp_cv_bounds[1]).collect { |i|
      sum = BigDecimal.new(0)
      (lp_bounds[0]..lp_bounds[1]).each { |j|
        sum += lp[j+lp_center] * cv[i-j-1+cv_center] if (0...cv_length).include?(i-j-1+cv_center)
      }
      sum
    }

    hp_cv = (hp_cv_bounds[0]..hp_cv_bounds[1]).collect { |i|
      sum = BigDecimal.new(0)
      (hp_bounds[0]..hp_bounds[1]).each { |j|
        sum += hp[j+hp_center] * cv[i-j-1+cv_center] if (0...cv_length).include?(i-j-1+cv_center)
      }
      sum
    }
    center = lp_cv_center
    low = ConvolutionFilter::new(name+"_l", lp_cv.collect { |e| e.truncate(30).to_s }, lp_cv_center)
    high = ConvolutionFilter::new(name+"_h", hp_cv.collect { |e| e.truncate(30).to_s }, hp_cv_center)
    return center,low,high
  end

  def split_filters(name,lfa,hfa,reverse,center_half)

    if (reverse) then
      lowin=lfa.reverse
      higin=hfa.reverse
      r='r'
    else
      lowin=lfa
      higin=hfa
      r=''
    end

    filt_1 = lowin.values_at(*(0..(lfa.length-1)).step(2).collect)
    low_e = ConvolutionFilter::new(name+"_l"+r+"e", filt_1, center_half)

    filt_2 = lowin.values_at(*(1..(lfa.length-1)).step(2).collect)
    low_o = ConvolutionFilter::new(name+"_l"+r+"o", filt_2, center_half)

    filt_3 = higin.values_at(*(0..(hfa.length-1)).step(2).collect)
    high_e = ConvolutionFilter::new(name+"_h"+r+"e", filt_3, center_half)

    filt_4 = higin.values_at(*(1..(hfa.length-1)).step(2).collect)
    high_o = ConvolutionFilter::new(name+"_h"+r+"o", filt_4, center_half)
    return low_e,low_o,high_e,high_o
  end

  def bc=(bc)
    super
    fill_even_and_odd_filters
    return bc
  end

  def fill_even_and_odd_filters
    if @low_even then
      return #just do it the first time
    end
    if @convolution then
      if @bc and (@bc.grow? or @bc.shrink?) then
        convcntr=@convolution.center-1
      else
        convcntr=@convolution.center
      end
      #must rebuild the low and the high with shifted center
      @center, @low, @high = self.combine_wavelet_filters(@name,@fil_array,@convolution,convcntr)
    end
    if @reverse
      cntr=@low.fil_array.length - @center - 1
    else
      cntr=@center
    end
    @low_even, @low_odd, @high_even, @high_odd = split_filters(name,@low.fil_array,@high.fil_array,@reverse,cntr/2)
  end

  private :fill_even_and_odd_filters

  def lowfil(wavelet=nil)
    return @low_even.lowfil
  end

  def upfil(wavelet=nil)
    return @low_even.upfil
  end

  def init_filter_val
    return (0..3).collect { |index| @tt[0][0].copy("filter_val#{index}") }
  end

  def compute_values( tt_ind, indexes, data)
    out_even = @tt[0][tt_ind]
    out_odd  = @tt[1][tt_ind]
    i_in0 = indexes[0].flatten
    i_in1 = indexes[1].flatten
    pr out_even === FMA(Load(data[*i_in0], out_even), @filter_val[0], out_even)
    pr out_odd  === FMA(Load(data[*i_in0], out_odd ), @filter_val[1], out_odd )
    pr out_even === FMA(Load(data[*i_in1], out_even), @filter_val[2], out_even)
    pr out_odd  === FMA(Load(data[*i_in1], out_odd ), @filter_val[3], out_odd )
  end

  def post_process_and_store_values( tt_ind, indexes, data, transpose, options = {} )
    out_even = @tt[0][tt_ind]
    out_odd  = @tt[1][tt_ind]
    i_out0 = indexes[0].rotate(transpose).flatten
    i_out1 = indexes[1].rotate(transpose).flatten

    if options[:a] then
      pr out_even === out_even * Set( options[:a], out_even)
      pr out_odd  === out_odd  * Set( options[:a], out_odd )
    end

    if options[:accumulate] then
      pr out_even === Load(data[*i_out0], out_even) + out_even
      pr out_odd  === Load(data[*i_out1], out_odd ) + out_odd
    elsif options[:a_y] then
      pr out_even === FMA(Load(data[*i_out0], out_even), Set(options[:a_y], out_even), out_even)
      pr out_odd  === FMA(Load(data[*i_out1], out_odd ), Set(options[:a_y], out_odd ), out_odd )
    end

    pr data[*i_out0] === out_even
    pr data[*i_out1] === out_odd
  end

  def set_zero_tt( tt_ind )
    pr @tt[0][tt_ind].set( 0.0 )
    pr @tt[1][tt_ind].set( 0.0 )
  end

  def decl_filters
    decl low_even.fil
    decl low_odd.fil
    decl high_even.fil
    decl high_odd.fil
  end

  def cost
    return @length * 2 * 2
  end

end

class WaveletFilterDecompose < WaveletFilter

  def initialize(name, filt1, opts = {})
    super

    @base_name = "s0s1"

    @reverse=false
    #@low_even, @low_odd, @high_even, @high_odd = split_filters(name,@low.fil_array,@high.fil_array,false,@center/2)

  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return [ Real("lt", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len),
               Real("ht", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len) ]
    else
      return [ (0..(@tt_n-1)).collect { |index| Real("lt#{index}", :vector_length => @vec_len) },
               (0..(@tt_n-1)).collect { |index| Real("ht#{index}", :vector_length => @vec_len) } ]
    end
  end

  def set_filter_val( index, options = {} )
    pr @filter_val[0] === Set(@low_even.fil[index], @tt[0][0])
    pr @filter_val[1] === Set(@high_even.fil[index], @tt[0][0])
    pr @filter_val[2] === Set(@low_odd.fil[index], @tt[0][0])
    pr @filter_val[3] === Set(@high_odd.fil[index], @tt[0][0])
  end

  def get_input_dim_std( dim )
    return [ Dim(0, 1), Dim( 0, dim - 1 ) ]
  end

  def get_input_dim_shrink( dim )
    return [ Dim(0, 1), Dim( lowfil, dim + upfil - 1 ) ]
  end

  def get_input_dim_ld_shrink( dim )
    return [ Dim(0, 1), Dim( lowfil, dim + lowfil - 1 ) ]
  end

  def get_output_dim_std( dim )
    return [ Dim( 0, dim - 1 ), Dim(0, 1) ]
  end

  def get_output_dim_grow( dim )
    return [ Dim( -upfil, dim - lowfil - 1 ), Dim(0, 1) ]
  end

  def get_output_dim_ld_grow( dim )
    return [ Dim( -upfil, dim - upfil - 1 ), Dim(0, 1) ]
  end

  def convert_output_indexes(indexes, processed_index)
    tmp = [[], []]
    indexes.each_with_index{ |indx, i|
      if i == processed_index then
        tmp[0][i] = [indx, 0]
        tmp[1][i] = [indx, 1]
      else
        tmp[0][i] = indx
        tmp[1][i] = indx
      end
    }
    return tmp
  end

  def convert_input_indexes(indexes, processed_index)
    tmp = [[], []]
    indexes.each_with_index{ |indx, i|
      if i == processed_index then
        tmp[0][i] = [0, indx]
        tmp[1][i] = [1, indx]
      else
        tmp[0][i] = indx
        tmp[1][i] = indx
      end
    }
    return tmp
  end

end

class WaveletFilterRecompose < WaveletFilter

  def initialize(name, filt1, opts = {})
    super

    @base_name = "s1s0"

    @reverse = true
    #@low_even, @low_odd, @high_even, @high_odd = split_filters(name,@low.fil_array,@high.fil_array,true,(@low.fil_array.length - @center - 1)/2)

  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return [ Real("et", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len),
               Real("ot", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len) ]
    else
      return [ (0..(@tt_n-1)).collect { |index| Real("et#{index}", :vector_length => @vec_len) },
               (0..(@tt_n-1)).collect { |index| Real("ot#{index}", :vector_length => @vec_len) } ]
    end
  end

  def set_filter_val( index, options = {} )
    pr @filter_val[0] === Set(@low_odd.fil[index], @tt[0][0])
    pr @filter_val[1] === Set(@low_even.fil[index], @tt[0][0])
    pr @filter_val[2] === Set(@high_odd.fil[index], @tt[0][0])
    pr @filter_val[3] === Set(@high_even.fil[index], @tt[0][0])
  end

  def get_input_dim_std( dim )
    return [ Dim( 0, dim - 1 ), Dim(0, 1) ]
  end

  def get_input_dim_shrink( dim )
    return [ Dim( lowfil, dim + upfil - 1 ), Dim(0, 1) ]
  end

  def get_input_dim_ld_shrink( dim )
    return [ Dim( lowfil, dim + lowfil - 1 ), Dim(0, 1) ]
  end

  def get_output_dim_std( dim )
    return [ Dim(0, 1), Dim( 0, dim - 1 ) ]
  end

  def get_output_dim_grow( dim )
    return [ Dim(0, 1), Dim( -upfil, dim - lowfil - 1 ) ]
  end

  def get_output_dim_ld_grow( dim )
    return [ Dim(0, 1), Dim( -upfil, dim - upfil - 1 ) ]
  end

  def convert_output_indexes(indexes, processed_index)
    tmp = [[], []]
    indexes.each_with_index{ |indx, i|
      if i == processed_index then
        tmp[0][i] = [0, indx]
        tmp[1][i] = [1, indx]
      else
        tmp[0][i] = indx
        tmp[1][i] = indx
      end
    }
    return tmp
  end

  def convert_input_indexes(indexes, processed_index)
    tmp = [[], []]
    indexes.each_with_index{ |indx, i|
      if i == processed_index then
        tmp[0][i] = [indx, 0]
        tmp[1][i] = [indx, 1]
      else
        tmp[0][i] = indx
        tmp[1][i] = indx
      end
    }
    return tmp
  end

end

# determine the BC of the convolutions
#        Typical values are 0 : periodic BC, the size of input and output arrays are identical
#                           1 : Free BC, grow: the size of the output array is equal to the
#                                              size of input array plus the one of the filter
#                          -1 : Free BC, shrink: the size of the input array is equal to the one
#                                        of output array plus the filter
#                     Given a convolution and its inverse, in free BC -1 is the inverse of 1
#                      but not viceversa as the number of point treated is lower.
#                         10:  Free BC, the size of input and output arrays are identical
#                              (loss of information)
#                          -2 : Non periodic BC, for nabla operators
class BoundaryConditions
  # conditions names
  PERIODIC = 0
  GROW = 1
  SHRINK = -1
  NPERIODIC = -2
  FREE = 2

  CONDITIONS = [PERIODIC, GROW, SHRINK]

  # original id of the boundary conditions
  attr_reader :id
  # name of the boundary condition, used for the name of the routines
  attr_reader :name
  # determine if the boundary condition is free or not
  attr_reader :free
  alias free? free
  # determine if the convolution is a grow-like type (cannot be true if shrink is true)
  attr_reader :grow
  alias grow? grow
  # determine if the convolution is a shrink-like type (cannot be true if grow is true)
  attr_reader :shrink
  alias shrink? shrink
  # determine if the convolution is non periodic
  attr_reader :nper
  alias nper? nper
  # determine if the convolution discards data
  attr_reader :discard
  alias discard? discard

  def initialize(ibc)
    @id     = ibc
    @free   = (ibc != PERIODIC && ibc != NPERIODIC)
    @nper   = (ibc == NPERIODIC)
    @grow   = (ibc == GROW)
    @shrink = (ibc == SHRINK)
    @discard = (ibc == FREE)
    if not @free then
      if not @nper then
        @name = 'p'
      else
        @name = 'np'
      end
    else
      @name = 'f'
      if @grow then
        @name += 'g'
      elsif @shrink then
        @name += 's'
      elsif @discard then
        @name += 'd'
      end
    end
  end
end #class BoundaryConditions

BC = BoundaryConditions


#handle the choice of the best kernels in a given architecture
class GenericOptimization

  attr_reader :repeat
  attr_reader :dimensions
  attr_reader :openmp

  class DataRange
    def initialize(start,stop,step)
      @range = start..stop
      @step = step
    end
    def each(&block)
      return @range.step(@step,&block)
    end
  end

  class ParamSpace
    attr_accessor :space
    def initialize(space={})
      @space=space
    end
    def size
      return @space.size
    end
    def points
      pts=[]
      space2 = @space.dup
      dimension,value = space2.shift
      space2=ParamSpace::new(space2)
      value.each{ |val|
        pts.push({dimension.to_s.chomp("_range").to_sym => val})
      }
      if space2.size == 0 then
        return pts
      else
        pts3=[]
        pts.each{ |p1|
          space2.each { |p2|
            pts3.push(p1.dup.update(p2))
          }
        }
        return pts3
      end
    end
    def each(&block)
      return self.points.each(&block)
    end
  end

  def initialize(options = {})
    vector_length = 1
    vector_length = options[:vector_length] if options[:vector_length]
    unroll_range = 1
    unroll_range = options[:unroll_range] if options[:unroll_range]
    mod_arr_test = true
    mod_arr_test = options[:mod_arr_test] if not options[:mod_arr_test].nil?
    tt_arr_test = false
    tt_arr_test = options[:tt_arr_test] if not options[:tt_arr_test].nil?
    unrolled_dim_index_test = false
    unrolled_dim_index_test = options[:unrolled_dim_index_test] if not options[:unrolled_dim_index_test].nil?
    unroll_inner_test = false
    unroll_inner_test = options[:unroll_inner_test] if not options[:unroll_inner_test].nil?
    @repeat = 3
    @repeat = options[:repeat] if options[:repeat]
    @dimensions = [124,132,130]
    @dimensions = [options[:dimensions]].flatten if options[:dimensions]
    @openmp = true
    @openmp = false if options[:openmp] == false

    unrl_rng=[unroll_range].flatten
    if unrl_rng.length == 2 then
      unrl_rng=[*unrl_rng,1]
    elsif unrl_rng.length == 1 then
      unrl_rng=[1,*unrl_rng,1]
    end
    space={}
    space[:unroll_range] = DataRange::new(*unrl_rng[0..2])
    space[:mod_arr_range] = [true]
    space[:mod_arr_range] = [true,false] if mod_arr_test
    space[:tt_arr_range] = [false]
    space[:tt_arr_range] = [true,false] if tt_arr_test
    space[:unrolled_dim_index_range] = [0]
    space[:unrolled_dim_index_range] = [0,1] if unrolled_dim_index_test
    space[:unroll_inner_range] = [true]
    space[:unroll_inner_range] = [true, false] if unroll_inner_test
    space[:vector_length_range] = [vector_length].flatten
    @space=ParamSpace::new(space)
  end

  def each(&block)
    return @space.each(&block)
  end

end

class Convolution1dShape

  attr_reader :dim_indexes
  attr_reader :bc
  attr_reader :transpose
  attr_reader :filter
  attr_reader :ld
  attr_reader :length

  attr_reader :processed_dim_index
  attr_reader :processed_dim
  attr_reader :non_processed_dim_indexes
  attr_reader :non_processed_dims

  attr_reader :dim_n
  attr_reader :dims, :dims_in, :dims_out

  attr_reader :dimx, :dimy

  attr_reader :line_start, :line_end
  attr_reader :border_low, :border_high

  attr_reader :iterators

  attr_reader :vector_length
  attr_reader :unroll_length
  attr_reader :unrolled_dim_index

  def initialize(dim_indexes, bc, transpose, filter, ld)
    @dim_indexes = dim_indexes
    @length = @dim_indexes.length
    @bc = bc
    @transpose = transpose
    @filter = filter
    @ld = ld

    @processed_dim_index = @dim_indexes[-1]
    @non_processed_dim_indexes = @dim_indexes[0..-2]

    generate_dims
    compute_dimx_dimy
    compute_inner_loop_boundaries

    @iterators =  (1..@length).collect{ |index| Int("i#{index}")}
  end

  def to_s
    return @dim_indexes.join('')
  end

  def vars
    vars = @dims.dup
    if @ld then
      vars.push( @dims_in[processed_dim_index], @dims_out[processed_dim_index] )
    end
    return vars
  end

  def set_exploration( unrolled_dim_index, unroll_length, vector_length, dot_in )
    @unroll_length = 1
    @unroll_length = unroll_length if unroll_length

    @unrolled_dim_index = @non_processed_dim_indexes[0]
    @unrolled_dim_index = non_processed_dim_indexes[unrolled_dim_index] if @length > 2 and unrolled_dim_index

    @vector_length = 1
    @vector_length = vector_length if not dot_in and @unrolled_dim_index == 0 and vector_length and vector_length <= @filter.length and @unroll_length % vector_length == 0
  end

  def for_parameters( dim_index, in_reliq, options = {} )
    if in_reliq then
      first, last, step = for_parameters_reliq( dim_index )
    else
      first, last, step = for_parameters_main( dim_index )
    end
    options[:step] = step
    options[:openmp] = true if dim_index == @non_processed_dim_indexes.first
    return [ @iterators[dim_index], first, last, options ]
  end

  def iterator( dim_index )
    return @iterators[dim_index]
  end

  def reliq?( dim_index )
    return ( @unrolled_dim_index == dim_index and @vector_length < @unroll_length )
  end

  def output_indexes( unroll_index )
    i_out = index_unroll( unroll_index )
    if @filter.kind_of?( WaveletFilter ) then
      return @filter.convert_output_indexes(i_out, @processed_dim_index) if @filter.kind_of?( WaveletFilter )
    else
      return i_out
    end
  end

  def input_indexes( unroll_index )
    i_in = index_unroll( unroll_index )
    if @filter.kind_of?( WaveletFilter ) then
      return @filter.convert_input_indexes(i_in, @processed_dim_index) if @filter.kind_of?( WaveletFilter )
    else
      return i_in
    end
  end

  def input_filter_indexes( unroll_index, filter_index, side )
    i_in = index_unroll( unroll_index )
    i_in[@processed_dim_index] = @filter.input_filter_index( i_in[@processed_dim_index], filter_index, side, @processed_dim )
    if @filter.kind_of?( WaveletFilter ) then
      return @filter.convert_input_indexes(i_in, @processed_dim_index) if @filter.kind_of?( WaveletFilter )
    else
      return i_in
    end
  end

  private

  def index_unroll( unroll_index )
    i = @iterators.dup
    i[@unrolled_dim_index] += unroll_index
    return i
  end

  def for_parameters_main( dim_index )
    first_index = 0
    if @unrolled_dim_index == dim_index then
      last_index = @dims[dim_index] - @unroll_length
      step = @unroll_length
    else
      last_index = @dims[dim_index] - 1
      step = 1
    end
    return [ first_index, last_index, step ]
  end

  def for_parameters_reliq( dim_index )
    return nil unless @unrolled_dim_index == dim_index
    first_index = (@dims[dim_index]/@unroll_length)*@unroll_length
    last_index = @dims[dim_index] - 1
    step = @vector_length
    return [ first_index, last_index, step ]
  end

  def generate_dims
    @dims = [nil]*@length
    @dim_n = Int(:n, :dir => :in)
    @dims[processed_dim_index] = @dim_n
    non_processed_dim_indexes.each { |dim|
      @dims[dim] = Int( "ndat#{dim}", :dir => :in)
    }
    @dims_in  = @dims.dup
    @dims_out = @dims.dup
    if @ld then
      @dims_in[processed_dim_index] = Int(:nx, :dir => :in)
      @dims_out[processed_dim_index] = Int(:ny, :dir => :in)
    end

    @processed_dim = @dims[@processed_dim_index]
    @non_processed_dims = @non_processed_dim_indexes.collect { |indx|
      @dims[indx]
    }
  end

  def compute_dimx_dimy
    @dimx = @dims_in.collect{ |dim|
      if dim.name.match("ndat") then
        Dim(0, dim - 1)
      else
        @filter.get_input_dim( dim, :ld => @ld )
      end
    }
    @dimx.flatten!
    @dimy = @dims_out.collect{ |dim|
      if dim.name.match("ndat") then
        Dim(0, dim - 1)
      else
        @filter.get_output_dim( dim, :ld => @ld )
      end
    }
    @dimy.rotate!(@transpose).flatten!
  end

  def compute_inner_loop_boundaries
    if @bc.grow? then
      @line_start = -@filter.upfil
      @line_end = @dim_n - @filter.lowfil - 1
    else
      @line_start = 0
      @line_end = @dim_n - 1
    end
    @border_low = -@filter.lowfil
    @border_high = @dim_n - @filter.upfil
  end

end

class ConvolutionOperator1d
  # Convolution filter
  attr_reader :filter
  # dimensions
  attr_reader :dims, :dim_n
  # input array, unchanged on exit
  attr_reader :in
  # output array
  # wavelet:     y <- [a] wt_fil (x) in + [ a_y * y ]
  # convolution: y <- [a]    fil (x) in + [ a_y * y ] + [ a_x * in ]
  attr_reader :y
  # reduce the array or not
  attr_reader :reduce
  # constants of the convolution
  attr_reader :a_x, :a, :a_y
  # value of <in | conv (in)>
  attr_reader :dotp
  # order of the dimensions to be convolved (las one is the processed dim)
  attr_reader :dim_indexes
  # variables of the procedure
  attr_reader :vars
  # options
  attr_reader :options
  attr_reader :base_name
  # Creates new 1d convolution
  #
  # ==== Attributes
  #
  # * +filter+ - ConvolutionFilter object corresponding to the operations to be applied on data
  # * +bc+ Boundary conditions: control the way in which the convolution has to be applied.
  #
  # * +options+ - Hash table of allowed options (see options descritpion)
  #
  # ==== Options
  #
  # * +:a+ - constant of y <- a wt_fil (x) in + [ a_y * y ] or y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 1
  # * +:a_y+ - constant of y <- a wt_fil (x) in + [ a_y * y ] or y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 0
  # * +:a_x+ - constant of y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 0
  # * +:ld+ - leading dimensions enable
  # * +:wavelet+ - specify a wavelet operation, :decompose or :recompose
  def initialize(filter, bc, dim_indexes, options={})
    @bc = bc
    @filter = filter.dup
    @filter.bc = @bc
    @transpose = options[:transpose]
    @dim_indexes = dim_indexes
    @ld = options[:ld]
    @kinetic = options[:kinetic]
    @wavelet = options[:wavelet]
    @poisson = options[:poisson]

    @shape = Convolution1dShape::new(dim_indexes, bc, @transpose, @filter, @ld)

    @vars = @shape.vars


    @vars.push @x = Real("x",:dir => :in, :dim => @shape.dimx, :restrict => true)
    if @kinetic and @kinetic != :inplace and not options[:zero_out_work] then
      if @bc.grow? then
        @vars.push @x2 = Real("x2",:dir => :in, :dim => @shape.dimy, :restrict => true)
      else
        @vars.push @x2 = Real("x2",:dir => :in, :dim => @shape.dimx, :restrict => true)
      end
    end
    @vars.push @y = Real("y",:dir => options[:a_y] ? :inout : :out, :dim => @shape.dimy, :restrict => true)
    if @kinetic and @transpose != 0 then
      @vars.push @y2 =  Real("y2", :dir => :out, :dim => @shape.dimy, :restrict => true)
    end
    @vars.push @a = Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x] #and init
    if options[:a_y] then
      if options[:a_y] == 1.0 then
        @accumulate = true
      else
        @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y]
      end
    end
    @vars.push @dot_in = Real("dot_in",:dir => :out) if options[:dot_in]
    @dot_in_tmp = nil
    @cost = Int("cost", :dir => :out)
    @l = Int :l
    @options = options
    @base_name = ""
    @base_name += "s_" if default_real_size == 4
    @base_name += "d_" if default_real_size == 8
    @base_name += @filter.base_name + "_"
    @base_name += @filter.name + "_" + @bc.name + "_#{@shape}"
    @base_name += "_a" if @a
    @base_name += "_ain" if @a_x
    @base_name += "_ay" if @a_y
    @base_name += "_acc" if @accumulate
    @base_name += "_dotin" if @dot_in
    @base_name += "_ld" if @ld
    @base_name += "_x2" if @x2
  end

  def params(dim, index=@shape.processed_dim_index)
    if @wavelet then
      dim[index] /= 2
    end
    vars=[]
    varsin=[]
    varsout=[]
    nd = { @shape.processed_dim.name => dim[index] } #, "ndat" => 1, "ndat1" => 1, "ndat2" => 1 }
    if @shape.length == 2 then
      nd[@shape.non_processed_dims.first.name] = 1
      dim.each_index { |indx|
        nd[@shape.non_processed_dims.first.name] *= dim[indx] if indx != index
      }
    else
      nd[@shape.non_processed_dims[0].name] = 1
      nd[@shape.non_processed_dims[1].name] = 1
      dim.each_index { |indx|
        nd[@shape.non_processed_dims[0].name] *= dim[indx] if indx < index
        nd[@shape.non_processed_dims[1].name] *= dim[indx] if indx > index
      }
    end
    n_push = lambda { |varsi, varso|
      s_n = nd[@shape.processed_dim.name]
      if @bc.grow? then
        varsi.push(s_n)
        if @wavelet then
          varso.push(s_n + @filter.low.length - 1)
        else
          varso.push(s_n + @filter.length - 1)
        end
      elsif @bc.shrink?
        if @wavelet then
          varsi.push(s_n + @filter.low.length - 1)
        else
          varsi.push(s_n + @filter.length - 1)
        end
        varso.push(s_n)
      else
        varsi.push(s_n)
        varso.push(s_n)
      end
    }
    @shape.dims.each { |dim|
      vars.push(nd[dim.name])
      if dim.name == @shape.processed_dim.name then
        n_push.call(varsin, varsout)
      else
        varsin.push(nd[dim.name])
        varsout.push(nd[dim.name])
      end
    }
    varsout.rotate!(@transpose)
    #input and output arrays
    if @ld then
      n_push.call(vars, vars)
    end
    case default_real_size
    when 4
      type = NArray::SFLOAT
    when 8
      type = NArray::FLOAT
    else
      raise "Unsupported precision!"
    end

    align = 64
    if @wavelet then
      vars.push(ANArray::new(type, align, *varsin,2).random!)
      vars.push(ANArray::new(type, align, *varsout,2))
    else
      vars.push(ANArray::new(type, align, *varsin).random!)
      if @x2 then
        if @bc.grow? then
          vars.push(ANArray::new(type, align, *varsout).random!)
        else
          vars.push(ANArray::new(type, align, *varsin).random!)
        end
      end
      vars.push(ANArray::new(type, align, *varsout))
      vars.push(ANArray::new(type, align, *varsout)) if @kinetic and @transpose != 0
    end
    #accessory scalars
    nscal=0
    nscal+=1 if @a
    nscal+=1 if @a_x
    nscal+=1 if @a_y
    nscal.times{vars.push(0.5)}
    vars.push(0.0) if @dot_in
    return vars
  end

  def optimize(opt_space)
    opt_space=GenericOptimization::new if not opt_space
    t_best=Float::INFINITY
    p_best_optim = nil
    already_tested = {}
    opt_space.each{ |optim|
      next if optim[:unrolled_dim_index] == 1 and @shape.length < 3
      #next if optim[:mod_arr] and @bc.free?
      #puts optim
      kernel = CKernel::new
      GenericConvolution.print_header
      p = self.procedure(optim)
      pr p
      kernel.procedure = p
      next if already_tested[p.name]
      #kernel.print #if @bc.free?
      kernel.build(:openmp => opt_space.openmp)
      dimensions = opt_space.dimensions
      par = nil
      if dimensions.length < @shape.length then
        dimensions += [dimensions[0]]*(@shape.length-dimensions.length)
      end
      stats_a = []
      par = self.params(dimensions.dup)
      #puts par.inspect
      opt_space.repeat.times {
        stats_a.push kernel.run(*par)
      }
      stats_a.sort_by! { |a| a[:duration] }
      stats = stats_a.first
      #puts *par[0...@dims.length]
      if get_verbose then
        puts "#{optim} - [#{par[0...@shape.length].join(", ")}] - #{kernel.procedure.name}: #{stats[:duration]*1.0e3} ms #{self.cost(*par[0...@shape.length]) / (stats[:duration]*1.0e9)} GFlops"
        puts optim
      end
      t_min = stats[:duration]
      puts "#{kernel.procedure.name}: #{t_min*1.0e3} ms #{self.cost(*par[0...@shape.length]) / (t_min*1.0e9)} GFlops"
      already_tested[p.name] = true
      if t_best > t_min then
        t_best = t_min
        p_best_optim = optim
      end
    }
    return self.procedure(p_best_optim)
  end

  def cost(*dimens)
    n = dimens[@shape.processed_dim_index]
    ndat = 1
    @shape.non_processed_dim_indexes.each { |indx|
      ndat *= dimens[indx]
    }
    return n * @filter.cost * ndat
  end

  def get_constants
    return [@filter.lowfil, @filter.upfil]
  end

  def stos(bool)
   if bool then
    return 't'
   else
    return 'f'
   end	
  end	

  def procedure(options={})
    #(unroll, unrolled_dim, use_mod, tt_arr)
    #default values
    register_funccall("modulo")
    mod_arr = true
    tt_arr = false
    unroll_inner = true
    mod_arr = options[:mod_arr] if not options[:mod_arr].nil?
    tt_arr = options[:tt_arr] if not options[:tt_arr].nil?
    unroll_inner = options[:unroll_inner] if not options[:unroll_inner].nil?

    @shape.set_exploration( options[:unrolled_dim_index], options[:unroll], options[:vector_length], @dot_in )

    mod_arr = false if @bc.free?
    util = options[:util]

    function_name = @base_name
    function_name += "_#{@shape.unrolled_dim_index}u#{@shape.unroll_length}_v#{@shape.vector_length}"
    function_name += '_'+stos(mod_arr)+'_'+stos(tt_arr)+'_'+stos(unroll_inner)

    function_name += "_" + util.to_s if util

    if util == :cost then
      return Procedure(function_name, @shape.dims + [@cost] ){
        pr @cost === @shape.processed_dim * @filter.cost * @shape.non_processed_dims.inject(&:*)
      }
    end

    @filter.init_optims(tt_arr, mod_arr, unroll_inner, @shape.unroll_length, @shape.vector_length)

    return Procedure(function_name, vars, :constants => get_constants ){
      @filter.decl_filters
      decl *@shape.iterators
      decl @l
      decl *([@filter.get_tt].flatten)
      @filter.decl_filter_val
      if @filter.get_mods then
        decl @filter.get_mods
        pr For(@l, @filter.get_mods_lower, @filter.get_mods_upper) {
          pr @filter.get_mods[@l] === modulo(@l, @shape.dim_n)
        }
      end
      if @options[:dot_in] and @shape.vector_length > 1 then
        decl @dot_in_tmp = @dot_in.copy("dot_in_tmp", :vector_length => @shape.vector_length, :dir => nil, :direction => nil)
      else
        @dot_in_tmp = @dot_in
      end
      pr @dot_in.set(0.0) if @options[:dot_in]
      pr OpenMP::Parallel(default: :shared, reduction: (@options[:dot_in] ? {"+" => @dot_in} : nil ), private: @shape.iterators + [@l] + [@filter.get_tt] + ( @filter.get_filter_val ? [@filter.get_filter_val].flatten : [] )) {
        convolution1d
      }
    }
  end

  #here follows the internal operations for the convolution 1d
  def convolution1d
    fors = @shape.non_processed_dim_indexes.collect { |index|
      For( * @shape.for_parameters( index, false ) )
    }
    fors.each { |f| opn f }
    conv_lines( @shape.unroll_length )
    fors.each_with_index.reverse_each { |f, index|
      close f
      if @shape.reliq?( @shape.non_processed_dim_indexes[index] ) then
        pr For( * @shape.for_parameters( @shape.non_processed_dim_indexes[index], true ) ) {
          fors[(index+1)..-1].each { |f_reliq| opn f_reliq }
          conv_lines( @shape.vector_length )
          fors[(index+1)..-1].reverse.each { |f_reliq| close f_reliq }
        }
      end
    }
  end

  def conv_lines( tlen )
    # the shrink operation contains the central part only
    iter = @shape.iterators[@shape.processed_dim_index]
    if @bc.shrink? then
      pr For(iter, @shape.line_start, @shape.line_end) {
        for_conv(:center, tlen)
      }
    else
      pr For(iter, @shape.line_start, @shape.border_low - 1) {
        for_conv(:begin, tlen)
      }
      pr For(iter, @shape.border_low, @shape.border_high - 1) {
        for_conv(:center, tlen)
      }
      pr For(iter, @shape.border_high, @shape.line_end) {
        for_conv(:end, tlen)
      }
    end
  end

  def init_values(tlen)
    (0...tlen).step(@shape.vector_length).each_with_index{ |_,tt_ind|
      #WARNING: the eks conditional here can be relaxed
      @filter.set_zero_tt(tt_ind)
    }
  end

  def compute_values(side, tlen)
    (0...tlen).step(@shape.vector_length).each_with_index{ |ind,tt_ind|

      i_in = @shape.input_filter_indexes( ind, @l, side )

      @filter.compute_values(tt_ind, i_in, @x)
    }
  end


  def post_process_and_store_values(side, tlen)
    (0...tlen).step(@shape.vector_length).each_with_index{ |ind,tt_ind|
      i_out = @shape.output_indexes(ind)
      i_in = @shape.input_indexes(ind)
      @filter.post_process_and_store_values( tt_ind,
                                             i_out,
                                             @y,
                                             @transpose,
                                             :side => side,
                                             :position => @shape.iterators[@shape.processed_dim_index],
                                             :dim => @shape.processed_dim,
                                             :accumulate => @accumulate,
                                             :a => @a,
                                             :a_y => @a_y,
                                             :a_x => @a_x,
                                             :dot_in_tmp => @dot_in_tmp,
                                             :kinetic => @kinetic,
                                             :zero_out => @options[:zero_out],
                                             :data_in => @x,
                                             :data_in2 => @x2,
                                             :indexes_in => i_in,
                                             :data_out2 => @y2 )
    }
  end

  def for_conv(side, tlen)

    iters = @shape.iterators

    init_values(tlen)

    loop_start, loop_end = @filter.get_loop_start_end( side, @shape.processed_dim, @shape.iterators[@shape.processed_dim_index] )

    pr For( @l, loop_start, loop_end, :unroll => @filter.unroll_inner ) {
      @filter.set_filter_val(@l, :side => side, :position => iters[@shape.processed_dim_index], :dim => @shape.processed_dim)
      compute_values(side, tlen)
    }

    post_process_and_store_values(side, tlen)

  end

end

class GenericConvolutionOperator1d
  attr_accessor :procs
  attr_reader :needed_subops
  def initialize(filter,options={})
    @filter = filter
    @ld = options[:ld]
    @narr = options[:narr]
    @wavelet = options[:wavelet]
    @kinetic = options[:kinetic]
    @poisson = options[:poisson]

    @vars = []
    @vars.push @ndim  = Int( "d",    :dir => :in )
    @vars.push @idim  = Int( "idim", :dir => :in )
    @vars.push @dims  = Int( "n",    :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    @vars.push @bc    = Int( "bc",   :dir => :in )
    if @ld then
      @vars.push @nx  = Int( "nx", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
      @vars.push @ny = Int( "ny", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    end
    @vars.push @narr     = Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = Real("x",  :dir => :in,  :restrict => true, :dim => [ Dim() ] )
    @vars.push @y     = Real("y",  :dir => options[:a_y] ? :inout : :out, :restrict => true, :dim => [ Dim() ] )
    @vars.push @a = Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y] and options[:a_y] != 1.0
    @vars.push @dot_in = Real("dot_in",:dir => :out) if options[:dot_in]
    @cost = Int( "cost", :dir => :out )

    @transpose = 0
    @options=options.dup
    @options[:transpose] = 0
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a)
    opts.delete(:a_x)
    opts.delete(:a_y)
    dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
    opt_base = []
    opt_base.push( { :a   => @options[:a]   } ) if @options[:a]
    opt_base.push( { :a_x => @options[:a_x] } ) if @options[:a_x]
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y] and @options[:a_y] != 1.0
    opts_bases = []
    (0..opt_base.length).each { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        ch.update( { :a_y => @options[:a_y] } ) if @options[:a_y] and @options[:a_y] == 1.0
        opts_bases.push(ch)
      }
    }
    if @poisson then
      conditions = BC::CONDITIONS << BC::NPERIODIC
    else
      conditions = BC::CONDITIONS
    end
    conditions.each { |bc|
      dim_indexes_a.each { |dim_indexes|
        opts_bases.each { |opt|
          op = opt.dup
          op.update(opts)
          p = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, op)
          @needed_subops[p.base_name] = p
          @procs[p.base_name] = p.procedure
          @procs[p.base_name+"_cost"] = p.procedure( :util => :cost )
        }
      }
    }
    p = self.procedure(:cost).first
    @procs[p.name] = p
  end

  def optimize(opt_space=nil)
    @needed_subops.each { |name, subop|
      @procs[name] = subop.optimize(opt_space)
    }
  end

  def cost( idim, n, bc, m = 1 )
    ndim = n.length
    if ndim == 1 then
      dims = [ n[0], 1 ]
      dim_indexes = [1, 0]
    else
      ni = n[idim]
      if idim == 0 then
        ndat_left = nil
      else
        ndat_left = 1
        n[0...idim].each { |val| ndat_left *= val }
      end
      if idim == ndim - 1 then
        ndat_right = nil
      else
        ndat_right = 1
        n[(idim+1)..-1].each { |val| ndat_right *= val }
      end
      if idim == 0 then
        dim_indexes = [1, 0]
      elsif idim == ndim - 1 then
        dim_indexes = [0, 1]
      else
        dim_indexes = [2, 0, 1]
      end
      dims = []
      dims.push( ndat_left ) if ndat_left
      dims.push( ni )
      dims.push( ndat_right ) if ndat_right
    end
    cost = ConvolutionOperator1d::new(@filter, BC::new(bc[idim]), dim_indexes, @options).cost( *dims )
    return cost * m
  end

  def procedure_name(util = nil)
    function_name = ""
    function_name += "d_" if default_real_size == 8
    function_name += "s_" if default_real_size == 4
    function_name += @filter.base_name
    function_name += "_1d_"
    function_name += @filter.name
    function_name += "_dotin" if @dot_in
    function_name += "_" + util.to_s if util
    function_name
  end

  def empty_procedure(util = nil)
    function_name = procedure_name(util)

    vv = @vars
    vv += [ @cost ] if util == :cost

    p = Procedure( function_name, vv )
  end

  def procedure(util = nil)

    function_name = procedure_name(util)

    vv = @vars
    vv += [ @cost ] if util == :cost

    p = Procedure( function_name, vv ) {
      ndat_left = Int "ndat_left"
      ndat_right = Int "ndat_right"
      nti = Int "nti" if @narr
      nto = Int "nto" if @narr
      i = Int "i"
      j = Int "j"
      tmp_cost = Int "c"
      decl i, ndat_left, ndat_right
      decl tmp_cost if util == :cost
      decl nti, nto, j if @narr
      if @narr and @ld then
        pr nti === @nx[@idim]
        pr nto === @ny[@idim]
      elsif @narr then
        pr If(@bc[i] == BC::SHRINK => lambda {
          if @wavelet then
            pr nti === @dims[@idim] + @filter.low.length - 1
          else
            pr nti === @dims[@idim] + @filter.length - 1
          end
          pr nto === @dims[@idim]
        }, @bc[i] == BC::GROW => lambda {
          pr nti === @dims[@idim]
          if @wavelet then
            pr nto === @dims[@idim] + @filter.low.length - 1
          else
            pr nto === @dims[@idim] + @filter.length - 1
          end
        })
      end
      if @narr and @wavelet then
        pr nti === nti * 2
        pr nto === nto * 2
      end
      dims = []
      dim_indexes = []
      dats = []
      if @narr then
        f = For(j, 0, @narr-1)
        dats[0] = (@x[nti*j+1]).address
        dats[1] = (@y[nto*j+1]).address
      else
        dats[0] = @x
        dats[1] = @y
      end

      selected_bc = nil

      print_call_generic = lambda { |bc, a, a_x, a_y|
        opts = @options.dup
        opts.delete(:a) if not a
        opts.delete(:a_x) if not a_x
        opts.delete(:a_y) if not a_y
        opts.update( { :a_y => 1.0 } ) if @options[:a_y] and @options[:a_y] == 1.0
        vars = []
        vars.push( @a ) if a
        vars.push( @a_x ) if a_x
        vars.push( @a_y ) if a_y
        vars.push( @dot_in ) if @dot_in
        vars.push( tmp_cost.address ) if util == :cost
        lds = []
        lds.push( @nx[@idim] ) if @ld
        lds.push( @ny[@idim] ) if @ld
        procname = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, opts).base_name
        procname += "_" + util.to_s if util
        args = dims + lds + dats + vars
        if util == :cost then
          pr @cost === 0
          args = dims + [tmp_cost.address]
        end
        opn f if @narr
          pr @procs[procname].call( *args )
          pr @cost === @cost + tmp_cost if util == :cost
        close f if @narr
      }

      print_call_param_a_y = lambda { |bc, a, a_x|
        if @a_y then
          pr If(@a_y == 0.0 => lambda {
            print_call_generic.call(bc, a, a_x, false)
          }, else: lambda {
            print_call_generic.call(bc, a, a_x, true)
          })
        else
          print_call_generic.call(bc, a, a_x, false)
        end
      }

      print_call_param_a_x = lambda { |bc, a|
        if @a_x then
          pr If(@a_x == 0.0 => lambda {
            print_call_param_a_y.call(bc, a, false)
          }, else: lambda {
            print_call_param_a_y.call(bc, a, true )
          })
        else
          print_call_param_a_y.call(bc, a, false)
        end
      }

      print_call_param_a = lambda { |bc|
        if @a then
          pr If(@a == 1.0 => lambda {
            print_call_param_a_x.call(bc, false)
          }, else: lambda {
            print_call_param_a_x.call(bc, true )
          })
        else
          print_call_param_a_x.call(bc, false)
        end
      }

      print_call = lambda {
        case_args = {
          BC::PERIODIC => lambda {
            print_call_param_a.call( BC::PERIODIC )
          },
          BC::GROW => lambda {
            print_call_param_a.call( BC::GROW )
          },
          BC::SHRINK => lambda {
            print_call_param_a.call( BC::SHRINK )
          }
        }
        case_args[BC::NPERIODIC] = lambda { print_call_param_a.call( BC::NPERIODIC ) } if @poisson
        pr Case( @bc, case_args)
      }
      pr If( @idim == 0 => lambda {
        pr ndat_right === 1
        pr For( i, 1, @ndim - 1 ) {
          if @ld and util != :cost then
            pr ndat_right === ndat_right * @nx[i]
          else
            pr ndat_right === ndat_right * @dims[i]
          end
          pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_right
          pr nto === nto * ndat_right
        end
        dim_indexes = [1,0]
        dims = [@dims[@idim], ndat_right]
        print_call.call
      }, @idim == @ndim - 1 => lambda {
        pr ndat_left === 1
        pr For( i, 0, @ndim - 2 ) {
          if @ld and util != :cost then
            pr ndat_left === ndat_left * @nx[i]
          else
            pr ndat_left === ndat_left * @dims[i]
          end
          pr ndat_left === ndat_left * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_left
          pr nto === nto * ndat_left
        end
        dim_indexes = [0,1]
        dims = [ndat_left, @dims[@idim]]
        print_call.call
      }, else: lambda {
        pr ndat_left === 1
        pr ndat_right === 1
        pr For( i, 0, @idim - 1 ) {
          if @ld and util != :cost then
            pr ndat_left === ndat_left * @nx[i]
          else
            pr ndat_left === ndat_left * @dims[i]
          end
          pr ndat_left === ndat_left * 2 if @wavelet
        }
        pr For( i, @idim + 1, @ndim - 1 ) {
          if @ld and util != :cost then
            pr ndat_right === ndat_right * @nx[i]
          else
            pr ndat_right === ndat_right * @dims[i]
          end
          pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_left * ndat_right
          pr nto === nto * ndat_left * ndat_right
        end
        dim_indexes = [2,0,1]
        dims = [ndat_left, @dims[@idim], ndat_right]
        print_call.call
      })
    }
    return [ p, @procs ]
  end

end

class GenericConvolutionOperator
  attr_accessor :procs
  attr_reader :needed_subops
  def initialize(filter,options={})
    @filter = filter
    @ld = options[:ld]
    @narr = options[:narr]
    @wavelet = options[:wavelet]
    @kinetic = options[:kinetic]

    @vars = []
    @vars.push @ndim  = Int( "d",  :dir => :in )
    @vars.push @dims  = Int( "n",  :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    @vars.push @bc    = Int( "bc", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    if @ld then
      @vars.push @nx  = Int( "nx", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
      @vars.push @ny = Int( "ny", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    end
    @vars.push @narr     = Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = Real("x",  :dir => :in,  :restrict => true, :dim => [ Dim() ] )
    @vars.push @y     = Real("y",  :dir => options[:a_y] ? :inout : :out, :restrict => true, :dim => [ Dim() ] )
    @vars.push @w1 = Real("w1", :dir => :inout, :restrict => true, :dim => [ Dim() ] ) if options[:work]
    @vars.push @w2 = Real("w2", :dir => :inout, :restrict => true, :dim => [ Dim() ] ) if options[:work]
    @vars.push @a = Real("a",:dir => :in,:dim => [ Dim(0, @ndim - 1)]) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y]
    @vars.push @dot_in = Real("dot_in",:dir => :out,:dim =>[ Dim(0, @ndim - 1)]) if options[:dot_in]

    @transpose = 0
    @transpose = options[:transpose] if options[:transpose]
    @options=options
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a)
    opts.delete(:a_x)
    opts.delete(:a_y)
    opts.delete(:zero_out_work)
    dim_indexes_a = []
    if @transpose == 0 then
      dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
    elsif @transpose == 1
      dim_indexes_a = [ [1, 0] ]
    elsif @transpose == -1
      dim_indexes_a = [ [0, 1] ]
    end
    opt_base = []
    init_options = {}
    init_options[:a_x] = @options[:a_x] if @options[:a_x]
    init_options[:zero_out_work] = @options[:zero_out_work] if @options[:zero_out_work]
    opt_base.push( init_options ) if init_options.length > 0
    opt_base.push( { :a => @options[:a] } ) if @options[:a]
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y]
    opts_bases = []
    (0..opt_base.length).each { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        opts_bases.push(ch)
      }
    }
    BC::CONDITIONS.each { |bc|
      dim_indexes_a.each{ |dim_indexes|
        opts_bases.each { |opt|
          op = opt.dup
          op.update(opts)
          p = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, op)
          @needed_subops[p.base_name] = p
          @procs[p.base_name] = p.procedure
        }
      }
    }
  end

  def optimize(opt_space=nil)
    @needed_subops.each { |name, subop|
      @procs[name] = subop.optimize(opt_space)
    }
  end

  def cost(n, bc, m = 1)
    dims_actual = []
    compute_ni_ndat = lambda { |indx|
      indx = n.length - indx - 1 if @transpose == -1
      ndat = 1
      (0..(n.length - 1)).each { |i|
        ndat *= dims_actual[i] if i != indx
      }
      ni_ndat = [ n[indx], ndat ]
      ni_ndat.reverse! if @transpose == -1 or (@transpose == 0 and indx = n.length - 1 )
      indexes = [ 1, 0]
      indexes.reverse! if @transpose == -1 or (@transpose == 0 and indx = n.length - 1 )
      return [ni_ndat, indexes]
    }
    compute_ndat_ni_ndat2 = lambda { |indx|
      ndat = 1
      ndat2 = 1
      (0..(n.length - 1)).each { |i|
        ndat *= dims_actual[i] if i < indx
        ndat2 *= dims_actual[i] if i > indx
      }
      ni = n[indx]
      ndat_ni_ndat2 = [ndat, ni, ndat2]
      indexes = [2, 0, 1]
      return [ndat_ni_ndat2, indexes]
    }
    (0...n.length).each { |indx|
      if bc[indx] == BC::SHRINK then
        if @wavelet then
          dims_actual[indx] = n[indx] * 2 + @filter.length - 2
        else
          dims_actual[indx] = n[indx] + @filter.length - 1
        end
      else
        if @wavelet then
          dims_actual[indx] = n[indx] * 2
        else
          dims_actual[indx] = n[indx]
        end
      end
    }
    change_dims = lambda { |indx|
      if bc[indx] == BC::GROW then
        if @wavelet then
          dims_actual[indx] = n[indx] * 2 + @filter.length  - 2
        else
          dims_actual[indx] = n[indx] + @filter.length - 1
        end
      else
        if @wavelet then
          dims_actual[indx] = n[indx] * 2
        else
          dims_actual[indx] = n[indx]
        end
      end
    }
    if n.length == 1 then
      d = [ n[0], 1 ]
      d_indexes = [ 1, 0 ]
      dims.reverse! if @transpose == -1
      dim_indexes.reverse! if @transpose == -1
      cost = ConvolutionOperator1d::new(@filter, BC::new(bc[0]), d_indexes, @options).cost( *d )
    else
      cost = 0
      dims, dim_indexes = compute_ni_ndat.call(0)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[0]), dim_indexes, @options).cost( *dims )
      change_dims.call(0)
      dims_left = n.length - 1
      while dims_left > 1 do
        if @transpose == 0 then
          dims, dim_indexes = compute_ndat_ni_ndat2.call(n.length-dims_left)
        else
          dims, dim_indexes = compute_ni_ndat.call(n.length-dims_left)
        end
        cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), dim_indexes, @options).cost( *dims )
        change_dims.call(n.length-dims_left)
        dims_left -= 1
      end
      dims, dim_indexes = compute_ni_ndat.call(n.length-1)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), dim_indexes, @options).cost( *dims )
    end
    return cost * m
  end

  def procedure
    function_name = ""
    function_name += "d_" if default_real_size == 8
    function_name += "s_" if default_real_size == 4
    if @wavelet then
      if @wavelet == :decompose then
        function_name += "dwt_"
      else
        function_name += "idwt_"
      end
    end
    function_name += @filter.name
    function_name += "_ld" if @ld
    function_name += "_narr" if @narr
    p = Procedure(function_name,@vars) {
      dims_actual = Int( "dims_actual", :allocate => true, :dim => [ Dim(0, @ndim - 1) ] ) if get_lang == FORTRAN
      dims_actual = Int( "dims_actual", :allocate => true, :dim => [ Dim(0, 16) ] ) if get_lang == C
      dims_left   = Int "dims_left"
      ni = Int "ni"
      ndat = Int "ndat"
      ndat2 = Int "ndat2"
      ndat_tot_in = Int "nti" if @narr
      ndat_tot_out = Int "nto" if @narr
      i = Int "i"
      j = Int "j"
      decl i, j, dims_actual, dims_left, ni, ndat, ndat2
      decl ndat_tot_in, ndat_tot_out if @narr
      dims = []
      dim_indexes = []
      pr dims_left === @ndim
      pr For( i, 0, @ndim - 1 ) {
        if @ld then
          if @wavelet then
            pr dims_actual[i] === @nx[i] * 2
          else
            pr dims_actual[i] === @nx[i]
          end
        else
          pr If(@bc[i] == BC::SHRINK => lambda {
            if @wavelet then
              pr dims_actual[i] === @dims[i] * 2 + @filter.length - 2
            else
              pr dims_actual[i] === @dims[i] + @filter.length - 1
            end
          }, else: lambda {
            if @wavelet then
              pr dims_actual[i] === @dims[i] * 2
            else
              pr dims_actual[i] === @dims[i]
            end
          })
        end
      }
      compute_ni_ndat = lambda { |indx|
        indx = @ndim - 1 - indx if @transpose == -1
        pr ni === @dims[indx]
        pr ndat === 1
        pr For(j, 0, @ndim - 1) {
          pr If( j != indx ) {
            pr ndat === ndat * dims_actual[j]
          }
        }
        d = [ ni, ndat ]
        d_indexes = [ 1, 0]
        d.reverse! if @transpose == -1
        d_indexes.reverse! if @transpose == -1
        return [ d, d_indexes ]
      }
      compute_ndat_ni_ndat2 = lambda { |indx|
        pr ni === @dims[indx]
        pr ndat === 1
        pr ndat2 === 1
        pr For(j, 0, @ndim - 1) {
          pr If( j < indx => lambda {
            pr ndat === ndat * dims_actual[j]
          }, j > indx => lambda {
            pr ndat2 === ndat2 * dims_actual[j]
          })
        }
        return [ [ndat, ni, ndat2], [2, 0, 1] ]
      }

      opts = @options.dup
      opts.delete(:a_x)
      opts.delete(:a_y)
      opts.delete(:zero_out_work)

      print_call = lambda { |indx, init, last, datas, multi_conv|
        vars = dims
        if multi_conv then
          if dim_indexes.length == 2 then
            pr ndat_tot_in === dims[dim_indexes[0]]
          else
            pr ndat_tot_in === dims[dim_indexes[0]] * dims[dim_indexes[1]]
          end
          pr ndat_tot_out === ndat_tot_in
          pr ndat_tot_in === ndat_tot_in * dims_actual[indx]
          if @ld then
            if @wavelet then
              pr ndat_tot_out === ndat_tot_out * @ny[indx] * 2
            else
              pr ndat_tot_out === ndat_tot_out * @ny[indx]
            end
          end
          f = For(j, 0, @narr-1)
        end
        if @ld then
          vars.push @nx[indx]
          vars.push @ny[indx]
        end
        indx = @ndim - 1 - indx if @transpose == -1
        vars2 = []
        opt = {}
        vars2.push( @a[indx] ) if @options[:a]
        if init and @options[:a_x] then
          vars2.push( @a_x )
          opt[:a_x] = @options[:a_x]
        end
        if init and @options[:zero_out_work]
          opt[:zero_out_work] = @options[:zero_out_work]
        end
        if last and @options[:a_y] then
          vars2.push( @a_y )
          opt[:a_y] = @options[:a_y]
        end
        vars2.push( @dot_in[indx] ) if @options[:dot_in]
        dats = []
        opt.update( opts )
        case_args = {
          BC::PERIODIC => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::PERIODIC), dim_indexes, opt).base_name

            if multi_conv then
              pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            pr dims_actual[indx] === @ny[indx]  if @ld
          },
          BC::GROW => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::GROW), dim_indexes, opt).base_name

            if multi_conv then
              if not @ld then
                if @wavelet then
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 2 )
                else
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 1 )
                end
              end
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            if @ld then
              pr dims_actual[indx] === @ny[indx]  if @ld
            else
              if @wavelet then
                pr dims_actual[indx] === dims_actual[indx] + @filter.length - 2
              else
                pr dims_actual[indx] === dims_actual[indx] + @filter.length - 1
              end
            end
          },
          BC::SHRINK => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::SHRINK), dim_indexes, opt).base_name

            if multi_conv then
              if not @ld then
                if @wavelet then
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 2 )
                else
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 1 )
                end
              end
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            if @ld then
              pr dims_actual[indx] === @ny[indx]  if @ld
            else
              if @wavelet then
                pr dims_actual[indx] === dims_actual[indx] - @filter.length + 2
              else
                pr dims_actual[indx] === dims_actual[indx] - @filter.length + 1
              end
            end
          }
        }
        case_args[BC::NPERIODIC] = lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::NPERIODIC), dim_indexes, opt).base_name

          if multi_conv then
            pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            pr f
          else
            dats = datas.dup
          end
          pr @procs[procname].call( *vars, *dats, *vars2 )
          close f if multi_conv
          pr dims_actual[indx] === @ny[indx]  if @ld
        } if @poisson
        pr Case( @bc[indx], case_args)
      }

      pr If( @ndim == 1 => lambda {
        conv_number = 1
        conv_number = @narr if @narr
        dims = [ @dims[0], conv_number ]
        dim_indexes = [ 1, 0 ]
        dims.reverse! if @transpose == -1
        dim_indexes.reverse! if @transpose == -1
        datas = [ @x, @y ]
        print_call.call( 0, true, true, datas, false )
      }, else: lambda {
        dims, dim_indexes = compute_ni_ndat.call( 0 )
        datas = [ @x, @w1 ]
        datas = [ @x, @y ] if not @options[:work]
        print_call.call( 0, true, false, datas, @narr )
        pr dims_left === dims_left - 1
        pr i === 1
        pr While( dims_left > 2 ) {
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w1, @w2 ]
          else
            datas = [ @w1, @w2 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @narr )
          pr i === i + 1
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w2, @w1 ]
          else
            datas = [ @w2, @w1 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @narr )
          pr i === i + 1
          pr dims_left === dims_left - 2
        }
        pr If( dims_left == 2 => lambda {
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w1, @w2 ]
          else
            datas = [ @w1, @w2 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @narr )
          pr i === i + 1
          dims, dim_indexes = compute_ni_ndat.call( i )
          if @transpose == 0 then
            dims.reverse!
            dim_indexes.reverse!
          end
          if @kinetic then
            datas = [ @x, @w2, @y ]
          else
            datas = [ @w2, @y ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, true, datas, @narr )
        }, else: lambda {
          dims, dim_indexes = compute_ni_ndat.call( i )
          if @transpose == 0 then
            dims.reverse!
            dim_indexes.reverse!
          end
          if @kinetic then
            datas = [ @x, @w1, @y ]
          else
            datas = [ @w1, @y ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, true, datas, @narr )
        })
      })
    }
    return [ p, @procs ]
  end
end

class ConvolutionOptimization
  # apply transposition paradigm or not in different directions: +1,0,-1
  attr_reader :transpose
  # order of the treated dimensions
  attr_reader :dim_order
  # 2-uple containing the unrolling lengths and the unrolling dims
  attr_reader :unroll
  # use the mod_arr strategy for the convolutions
  attr_reader :use_mod
  # use the tt_arr strategy for the convolutions (temporary variables are arrays and not scalars)
  attr_reader :tt_arr
  def initialize(convolution,options)

    ndim = convolution.ndim

    @transpose = 0
    @transpose = options[:transpose] if options[:transpose]

    @use_mod =  [ false ] * ndim
    convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc.free?) } if options[:use_mod]

    @tt_arr = convolution.dims.collect { false }
    if options[:tt_arr] then
      ttopt=[options[:tt_arr]].flatten
      if ttopt.length == 1 then
        @tt_arr = [ttopt[0]] * ndim
      elsif ttopt.length == ndim then
        @tt_arr = ttopt
      else
        raise 'Incoherent dimensions specified in tt_arr options: #{ndim}, #{ttopt.length}'
      end
    end

    convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc.free?) }

    @dim_order=(0...ndim).collect{|i| i}
    @dim_order.reverse!  if @transpose == -1
    @dim_order = options[:dim_order] if options[:dim_order] and @transpose == 0

    if @transpose ==1 then
      unrolled_dim = [ 1 ] * ndim
    elsif @transpose == -1 then
      unrolled_dim = [ 0 ] * ndim
    else
      hdim = ndim / 2
      unrolled_dim = (1..ndim).collect { |i| i < hdim ? 2 : 0 }
      unrolled_dim = [1] + unrolled_dim
      if options[:unrolled_dim] then
        raise 'Incoherent dimensions specified in unrolling options: #{ndim}, #{options[:unrolled_dim].length}' if options[:unrolled_dim].length != ndim
        unrolled_dim = options[:unrolled_dim]
      end
    end
    if options[:unroll] then
      unro = [options[:unroll]].flatten
      if unro.length == 1 then
        unroll_lengths = unro * ndim
      elsif unro.length == ndim then
        unroll_lengths = unro
      else
        raise 'Incoherent dimensions specified in unrolling options: #{ndim}, #{unro.length}'
      end
    else
        unroll_lengths = [1] * ndim
    end
    @unroll = [ unroll_lengths, unrolled_dim ]
  end
end

def self.print_header(macro = false)
  if get_lang == C then
    get_output.puts "#include <immintrin.h>" if get_architecture == X86
    get_output.puts "#include <arm_neon.h>" if get_architecture == ARM
    if macro then
      get_output.print "#define modulo(a, b) ((a+b)%(b))\n"
      get_output.print "#define min( a, b) ((a) < (b) ? a : b)\n"
      get_output.print "#define max( a, b) ((a) > (b) ? a : b)\n"
    else
      get_output.print "static inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      get_output.print "static inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      get_output.print "static inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end
  end
end

end
