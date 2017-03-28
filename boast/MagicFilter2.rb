require "BOAST"
require 'narray'
require "./GenericConvolution2.rb" 

def magicfilter_ref( invert = false, free = false )
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  lowfil=-8
  lupfil=7
  function_name = "magicfilter"
  if not free then
    function_name += "_per"
  else
    function_name += "_free"
  end
  function_name += "_t" if invert
  function_name += "_ref"
  n = BOAST::Int("n",:dir => :in, :signed => false)
  ndat = BOAST::Int("ndat",:dir => :in, :signed => false)

  dim_in_min = 0
  dim_in_max = n-1
  dim_out_min = 0
  dim_out_max = n-1
  if free then
    if invert then
      dim_in_min = lowfil
      dim_in_max = n-1+lupfil
    else
      dim_out_min = -lupfil
      dim_out_max = n-1-lowfil
    end
  end

  x = BOAST::Real("x",:dir => :in, 
                  :dim => [ BOAST::Dim(dim_in_min, dim_in_max), BOAST::Dim(ndat) ] )
  y = BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim(ndat), BOAST::Dim(dim_out_min, dim_out_max) ] )
  p = BOAST::Procedure(function_name, [n,ndat,x,y])
  kernel.code.print(File::read("magicfilter_refs.f90"))
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end

def MagicFilter(conv_filter, optims=GenericOptimization::new)
    

  conv_operation = GenericConvolutionOperator::new(conv_filter, :transpose => 1, :work => true, :ld => true, :narr => true)

  #test of 1d kernels optimizations in view of many-d
  conv_operation.optimize(optims)

  p, subops= conv_operation.procedure()

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
