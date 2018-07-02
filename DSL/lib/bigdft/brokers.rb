#This generates high level brokers for libconv, as well as header files to use for compilation against libconv, and a Makefile.
module Brokers
  prec = BOAST::Int("prec", :dir => :in)
  op = BOAST::Int("op", :dir => :in, :reference => 1)
  d = BOAST::Int("d", :dir => :in, :reference => 1)
  idim = BOAST::Int("idim", :dir => :in, :reference => 1)
  n = BOAST::Int("n", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ])
  bc = BOAST::Int("bc", :dir => :in, :reference => 1)
  nx = BOAST::Int("nx", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ])
  ny = BOAST::Int("ny", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ])
  narr = BOAST::Int("narr", :dir => :in, :reference => 1)
  kernels=[]
  precisions=[4,8]
  wavelet_families=["SYM8", "SYM4"]
  operations=["MF","IMF","DWT","RTOS1","IDWT","S1TOR","D1","D2"]
  families=["s0s0", "s0s0_dot", "s0s1", "s1s0"]
  dimensions=["1d"]

  def self.const_name(wavelet_family, operation)
    return "#{precision_name}_#{wavelet_family}_#{operation}"

  end

  id=0
  wavelet_families.each{ |wav_fam|
    operations.each{ |op|
    const_set("#{wav_fam}_#{op}", BOAST::Int("#{wav_fam}_#{op}", :constant => id) )
    id +=1
    }
  }
#generate both Fortran and C header files every time
    foldername=""
    if not BigDFT.from_cache then
      foldername = BigDFT::foldername
    end

  ["C", "Fortran"].each { |l|
    if l == "C" then
      header = "libconv.h"
    else
      header = "libconvf.h"
    end 
    File::open(foldername+header,"w") { |f|
    set_output(f)

    id=0

    if l == "C" then
      f.puts "enum ops{"
      wavelet_families.each{ |wav_fam|
        operations.each{ |op|
          f.print "#{wav_fam}_#{op}"
          f.print "," if op != operations.last
        }
      }
      f.puts "};"
    else 
      wavelet_families.each{ |wav_fam|
#       f.puts "module brokers"
        operations.each{ |op|
        f.print "integer :: " if id%3 == 0
        f.print "#{wav_fam}_#{op}"
        if id%3 != 2 and  !(op == operations.last and wav_fam == wavelet_families.last)
          f.print ", "
        else
          f.puts ""
        end
        id+=1
        }
      }
      id=0
      wavelet_families.each{ |wav_fam|
        operations.each{ |op|
          f.puts "parameter(#{wav_fam}_#{op}=#{id})"
          id+=1
        }
      }
#       f.puts "end module brokers"
    end
    }
  }
  
  filename = "brokers.f90"
  if BOAST::get_lang == BOAST::C then
    filename = "brokers.c"
  end
  File::open(foldername+filename,"w") { |f|
  set_output(f)
  #interfaces, commented for now
  if BOAST::get_lang == BOAST::C then
    f.puts  "#include <stdint.h>"
  end

#  families.each{ |family|
#      dimensions.each{ |dimension|
#          precisions.each{ |precision|
#              case precision
#              when 4
#                  precision_name = "s_"
#              when 8
#                  precision_name = "d_"
#              end
#              broker_name = "#{family}_#{dimension}"
#              if BOAST::get_lang == BOAST::FORTRAN then
#                  f.puts "INTERFACE #{broker_name}" if precision == 4
#                  f.puts "module procedure #{precision_name}#{broker_name}"
#                  f.puts "END INTERFACE #{broker_name}" if precision == 8
#              else
#                f.puts "#define #{broker_name}(d, idim, n, bc, nx, ny, narr, x, ...) _Generic((x),\\" if precision == 4
#                f.puts "float*: #{precision_name}#{broker_name} (d, idim, n, bc, nx, ny, narr, x, __VA_ARGS__),\\" if precision == 4
#                f.puts "double*: #{precision_name}#{broker_name} (d, idim, n, bc, nx, ny, narr, x, __VA_ARGS__))" if precision == 8
#              end
#              }
#          }
#      }
      
#  f.puts "contains" if BOAST::get_lang == BOAST::FORTRAN
  
wavelet_families.each{ |wav_fam|
  families.each{ |family|
      case family
      when "s0s0"
          operations=[const_get("#{wav_fam}_MF"), const_get("#{wav_fam}_IMF")]
      when "s0s0_dot"
          operations=[const_get("#{wav_fam}_D1"), const_get("#{wav_fam}_D2")]
      when "s0s1"
          operations=[const_get("#{wav_fam}_DWT"), const_get("#{wav_fam}_RTOS1")]
      when "s1s0"
          operations=[const_get("#{wav_fam}_IDWT"), const_get("#{wav_fam}_S1TOR")]
      end

      precisions.each{ |precision|
          case precision
          when 4
              precision_name = "s"
          when 8
              precision_name = "d"
          end
          
          generate_kernel = Proc.new  { |util|
          
              BOAST.push_env( default_real_size: precision ){
                  dimensions.each{ |dimension|
                      kernels=[]
                      function_name = "#{precision_name}_#{family}_#{dimension}"
                      function_name += "_#{util}" if util
                      function_name += "_" if BOAST::get_lang == BOAST::C
                      a = BOAST::Real("a", :dir => :in, :reference => 1 )
                      a_x = BOAST::Real("a_x", :dir => :in, :reference => 1 )
                      a_y = BOAST::Real("a_y", :dir => :in, :reference => 1 )
                      dot_in = Real("dot_in",:dir => :inout, :dim => [ BOAST::Dim(1)])
                      x = BOAST::Real("x", :dir => :in, :dim => [ BOAST::Dim()] )
                      y=BOAST::Real("y", :dir => :inout, :dim => [ BOAST::Dim()] )
                      cost = BOAST::Int("cost", :dir => :out, :dim => [ BOAST::Dim(1)])
                      dims = BOAST::Int("dims", :dir => :out, :dim => [ BOAST::Dim(0, d - 1)])
                      
                      operations.each{ |operation|
                          if util then 
                              kernels.push BOAST::const_get("#{precision_name}".upcase+"_#{operation.name}")
                          else 
                              kernels.push BOAST::const_get("#{precision_name}".upcase+"_#{operation.name}").kernel
                          end
                      }
                      
                      kernel = BOAST::CKernel::new(:kernels => kernels)
                      
                      vars = [op, d, idim, n, bc]
                      vars+=[nx, ny, narr] if util != "dims"
                      vars+=[x, y] if not util
                      vars.push a
                      vars.push a_x if family == "s0s0" or family == "s0s0_dot"
                      vars.push a_y
                      vars.push dot_in if family == "s0s0_dot"
                      vars.push cost if util == "cost"
                      vars.push dims if util == "dims"
                      p = BOAST::Procedure(function_name, vars){
                          vars.shift
                          case_args={}
                          operations.each_with_index{|op, i|
                              extra=nil
                              if util == "cost" 
                                  proc = kernels[i].cost_procedure
                                  extra = cost
                              elsif  util == "dims"
                                  proc = kernels[i].dims_procedure
                              else
                                  proc = kernels[i].procedure
                              end
                              case_args[op] = lambda {
                                  BOAST::pr proc.call(*vars)
                              }
                          }
                          BOAST::pr BOAST::Case( op, case_args)
                      }
                      BOAST::pr p
                      kernel.procedure = p
                      f.puts kernel
                      }
                  }
              }
              
              generate_kernel.call()
              generate_kernel.call("cost")
              generate_kernel.call("dims")
          }
      }
  }
}


# write Makefile.am to foldername

if not BigDFT.from_cache then
    compiler_options = BOAST::get_compiler_options

    makelines="""
    lib_LIBRARIES = libconv.a

    #temporary compiling line for gfortran
    AM_FCFLAGS = -I. #{compiler_options[:FCFLAGS]} -march=#{BOAST::get_model} #{BOAST::get_openmp_flags[compiler_options[:FC]]}

    #temporary compiling line for gcc
    AM_CFLAGS = -I. #{compiler_options[:CFLAGS]} -march=#{BOAST::get_model} #{BOAST::get_openmp_flags[compiler_options[:CC]]}
    
    libconv_a_SOURCES = """

    BigDFT::kernels.each{|f| makelines+=f.to_str+"\\\n"}

    makelines+=" #{filename}\n"

#    makelines+="libconv_a_OBJECTS = $(libconv_a_SOURCES:.f90=.o)\n"

#    makelines+="%.o: %.f90
#	    gfortran $(AM_FCFLAGS) -c $< \n
#    libconv.a: $(libconv_a_OBJECTS)
#	    gfortran $(AM_FCFLAGS) -c -o $@ $^"
    File::open("#{BigDFT::foldername}Makefile.am","w") {|f|
      f.puts makelines
    }
end

#    return kernel

end 
