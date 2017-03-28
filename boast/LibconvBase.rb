#define the basic operations for libconv generation including the
#reference operations coming from BigDFT sources
require './GenericConvolution2.rb'
require './WaveletFilters.rb'
require './MagicFilter2.rb'
require './AnaRotPer-2.rb'
require './Synthesis.rb'

openmp = true
def MagicFilter1d(filter, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(filter, :ld => true, :narr => true, :a_x => true,:a_y=>true, :a => true)
  conv_operation.optimize(optims) if optims

  p, subops = conv_operation.procedure

  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op|
    BOAST::pr op
    puts "chosen:"+ op.name
  }
  BOAST::pr p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
  return kernel
end
def Wavelet1d(wavelet_filter, direction, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(wavelet_filter, :wavelet => direction, :a => true, :a_y => true, :ld => true, :narr => true)
  conv_operation.optimize(optims) if optims

  p, subops = conv_operation.procedure

  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op|
    BOAST::pr op
    puts "chosen:"+ op.name
  }
  BOAST::pr p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
  return kernel
end

def dataspace(actual_ld,random=false)
  type = NArray::FLOAT
  out = ANArray::new(type, 32, *actual_ld)
  out.random! if random
  return out
end


optims = GenericOptimization::new(:unroll_range => [1,1], :mod_arr_test => false, :tt_arr_test => false)
def get_reference_kernels()
  sref = BOAST::synthesis_per_ref
  sref.build(verbose => true)
  
  aref = BOAST::analysis_free_ref
  aref.build(verbose=>true)
  
  #now take the real space and apply the magic filters
  mf_ref = magicfilter_ref
  mf_ref.build(verbose => true)
  
  return sref,aref,mf_ref
end

def get_reference_output(kernel,ndims,bc,input)
  #case bc periodic
  type = NArray::FLOAT
  outputt = ANArray::new(type, 32, 2*ndims[1], 2*ndims[2], 2*ndims[0])
  kernel.run(ndims[0],4*ndims[1]*ndims[2],input,outputt)
  return outputt.transpose(2,0,1)
end

def allclose(ref1,ref2,epsilon=10e-12)
  diff = (ref1 - ref2).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts 'ok'
end

def shift_realspace(ndims,src,dest)
  n1=ndims[0]
  n2=ndims[1]
  n3=ndims[2]
  dest[1..(2*n1-1),0..(2*n2-1),0..(2*n3-1)]=src[0..(2*n1-2),0..(2*n2-1),0..(2*n3-1)]
  dest[0..0,0..(2*n2-1),0..(2*n3-1)]       =src[(2*n1-1)..(2*n1-1),0..(2*n2-1),0..(2*n3-1)]
end

def outdim(bc,n,op)
  case bc
  when BC::PERIODIC 
    return n
  when BC::GROW
    case op
    when :mf
      return n+15
    when :iwt
      return n+7
    when :s1tor
      return n+15
    else
      puts op.op
      raise 'op not yet defined for grow'
    end
  when BC::SHRINK
    case op
    when :imf
      return n-15
    when :dwt
      return n-7
    when :rtos1
      return n-15
    else
      puts op.op
      raise 'op not yet defined for shrink'
    end
  else
    raise 'unknown bc'
  end
end

def get_ld(idim,ndims,bc,op)
  ld_dest=ndims.dup
  ld_dest[idim]=outdim(bc,ld_dest[idim],op)
  return ld_dest
end

class LibConvOp
  attr_reader :op
  def initialize(op)
    @op=op
  end
  def inverse
    case @op
    when :dwt
      return :iwt
    when :iwt
      return :dwt
    when :mf
      return :imf
    when :imf
      return :mf
    when :s1tor
      return :rtos1
    when :rtos1
      return :s1tor
    else
      return @op
    end
  end
end

#class LibConvDatas1
#  attr_reader :ndim
#  attr_reader :bc
#  attr_reader :leading_dimensions
#  def initialize(ndim,bc=BC::PERIODIC,random=false,fill=nil)
#    @bc=bc
#    @ndim=ndim
#    @leading_dimensions=get_ld
#  end
#end

class LibConvKernel
  attr_reader :kernel
  attr_reader :op
  def initialize(kernel,op)
    @kernel=kernel
    @op=op
  end
  def run(idim,ndims,bc,src,dest,*as)
    case bc
    when BC::GROW,BC::PERIODIC
      ld_src=ndims
      ld_dest=get_ld(idim,ndims,bc,@op.op)
    when BC::SHRINK
      ld_dest=ndims
      ld_src=get_ld(idim,ndims,BC::GROW,@op.inverse)
    end
    @kernel.run(ndims.length,idim,ndims,bc,ld_src,ld_dest,1,src,dest,*as)
  end
  def in_lds(ndim,lds,shrink)
    ld_tmp=lds.dup
    if shrink then
      ld_tmp[idim]=outdim(BC::GROW,ndim,@op.inverse)
    else
      ld_tmp[idim]=ld
    end
    return ld_tmp
  end
  def out_lds(ndim,lds,shrink)
    ld_tmp=lds.dup
    if shrink then
      return ld_tmp[idim]=ld
    else
      return ld_tmp[idim]=outdim(bc,ndim,@op)
    end
    return ld_tmp
  end
  def dump_to_file(directory="src/")
    filename=@kernel.procedure.name
    case BOAST::get_lang
    when BOAST::C
      suffix = ".c"
    when BOAST::FORTRAN
      suffix = ".f90"
    end
    filedump="#{filename}#{suffix}"
    Dir.mkdir(directory) if not Dir.exist?(directory)
    File::open(directory+filedump,"w") {|f|
      f.puts @kernel
    }
    return filedump
  end
end

def generate_dwt(optims=nil)
  wave_newdwt=WaveletFilterDecompose::new("sym#{SYM8_LP.length/2}", SYM8_LP)
  dwt = Wavelet1d( wave_newdwt, :decompose, optims )
  return LibConvKernel::new(dwt,LibConvOp::new(:dwt))
end

def generate_iwt(optims=nil)
  wave_newiwt=WaveletFilterRecompose::new("sym#{SYM8_LP.length/2}", SYM8_LP)
  iwt = Wavelet1d( wave_newiwt, :recompose, optims )
  return LibConvKernel::new(iwt,LibConvOp::new(:iwt))
end

def generate_mf(optims=nil)
  conv_filter = ConvolutionFilter::new('sym8_md',SYM8_MF.reverse,8)
  mf = MagicFilter1d( conv_filter, optims )
  return LibConvKernel::new(mf,LibConvOp::new(:mf))
end

def generate_imf(optims=nil)
  conv_filteri = ConvolutionFilter::new('sym8_imd',SYM8_MF,7)
  imf = MagicFilter1d( conv_filteri, optims )
  return LibConvKernel::new(imf,LibConvOp::new(:imf))
end

def generate_s1tor(optims=nil)
  conv_filteritmp = ConvolutionFilter::new('sym8_md',SYM8_MF,8)
  icomb=WaveletFilterRecompose::new("symicomb#{SYM8_LP.length/2}", SYM8_LP,:convolution_filter=>conv_filteritmp)
  iwticomb = Wavelet1d( icomb, :recompose, optims )
  return LibConvKernel::new(iwticomb,LibConvOp::new(:s1tor))
end

def generate_rtos1(optims=nil)
  conv_filteritmp = ConvolutionFilter::new('sym8_md',SYM8_MF,8)
  icombi=WaveletFilterDecompose::new("symicomb#{SYM8_LP.length/2}", SYM8_LP,:convolution_filter=>conv_filteritmp)
  dwticomb = Wavelet1d( icombi, :decompose, optims )
  return LibConvKernel::new(dwticomb,LibConvOp::new(:rtos1))
end

