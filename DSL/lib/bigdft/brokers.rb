#    def run(idim, bc, source, target, **options)


#def brokers

  prec = BOAST::Int("prec", :dir => :in)
  op = BOAST::Int("op", :dir => :in)
#  wavelet_fam = BOAST::Int("wavelet_fam", :dir => :in)
    
  d = BOAST::Int("d", :dir => :in)
  idim = BOAST::Int("idim", :dir => :in)
  n = BOAST::Int("n", :dir => :in)
  bc = BOAST::Int("bc", :dir => :in)
  nx = BOAST::Int("nx", :dir => :in)
  ny = BOAST::Int("ny", :dir => :in)
  narr = BOAST::Int("narr", :dir => :in)

  kernels=[]
  #TODO : hashtable to get kernels
#  operations=["DWT", "IDWT", "MF", "IMF", "S1TOR", "RTOS1"]
  precisions=[4,8]
  wavelet_families=["SYM8"]
    
  families=["s0s0", "s0s1", "s1s0"]
  dimensions=["1d"]
  
  filename = "broker.f9O"
  if BOAST::get_lang == BOAST::C then
    filename = "broker.c"
  end
  
  
  File::open(filename,"w") { |f|
  set_output(f)


#TODO : hashmap
  SYM8_MF2 = BOAST::Int("SYM8_MF")
  SYM8_IMF = BOAST::Int("SYM8_IMF")
  SYM8_DWT = BOAST::Int("SYM8_DWT")
  SYM8_RTOS1 = BOAST::Int("SYM8_RTOS1")
  SYM8_IDWT = BOAST::Int("SYM8_IDWT")
  SYM8_S1TOR = BOAST::Int("SYM8_S1TOR")
  SYM8_D12 = BOAST::Int("SYM8_D1")
  SYM8_D22 = BOAST::Int("SYM8_D2")
  
  if BOAST::get_lang == BOAST::C then 
      f.puts "enum ops {SYM8_MF, SYM8_IMF, SYM8_DWT,SYM8_RTOS1,SYM8_IDWT,SYM8_S1TOR,SYM8_D1,SYM8_D2};"
  else 
      f.puts "parameter(SYM8_MF=0)"
      f.puts "parameter(SYM8_IMF=1)"
      f.puts "parameter(SYM8_DWT=2)"
      f.puts "parameter(SYM8_RTOS1=3)"
      f.puts "parameter(SYM8_IDWT=4)"
      f.puts "parameter(SYM8_S1TOR=5)"
      f.puts "parameter(SYM8_D1=6)"
      f.puts "parameter(SYM8_D2=7)"
  end
  
  #interfaces, only for fortran
  

  families.each{ |family|
      dimensions.each{ |dimension|
          precisions.each{ |precision|
              case precision
              when 4
                  precision_name = "s_"
              when 8
                  precision_name = "d_"
              end
              broker_name = "#{family}_#{dimension}"
              if BOAST::get_lang == BOAST::FORTRAN then
                  f.puts "INTERFACE #{broker_name}" if precision == 4
                  f.puts "module procedure #{precision_name}#{broker_name}"
                  f.puts "END INTERFACE #{broker_name}" if precision == 8
              else
                f.puts "#define #{broker_name}(d, idim, n, bc, nx, ny, narr, x, ...) _Generic((x),\\" if precision == 4
                f.puts "float*: #{precision_name}#{broker_name} (d, idim, n, bc, nx, ny, narr, x, __VA_ARGS__),\\" if precision == 4
                f.puts "double*: #{precision_name}#{broker_name} (d, idim, n, bc, nx, ny, narr, x, __VA_ARGS__))" if precision == 8
              end
              }
          }
      }
  

#s0s0 : MF, IMF, + D1, D2
  families.each{ |family|
      case family
      when "s0s0"
          operations=[SYM8_MF2, SYM8_IMF]
      when "s0s1"
          operations=[SYM8_DWT, SYM8_RTOS1]
      when "s1s0"
          operations=[SYM8_IDWT, SYM8_S1TOR]
      end

      precisions.each{ |precision|
          case precision
          when 4
              precision_name = "s"
          when 8
              precision_name = "d"
          end
          
          generate_kernel = Proc.new  { |suffix|
          
              BOAST.push_env( default_real_size: precision ){
                  dimensions.each{ |dimension|
                      kernels=[]
                      function_name = "#{precision_name}_#{family}_#{dimension}"
                      function_name += "_#{suffix}" if suffix

                      a = BOAST::Real("a", :dir => :in )
                      a_x = BOAST::Real("a_x", :dir => :in )
                      a_y = BOAST::Real("a_y", :dir => :in )
                      x = BOAST::Real("x", :dir => :in, :dim => [ BOAST::Dim()] )
                      y=BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim()] )
                      
                      operations.each{ |operation|
                          if suffix then 
                              kernels.push BOAST::const_get("#{precision_name}".upcase+"_#{operation.name}")
                          else 
                              kernels.push BOAST::const_get("#{precision_name}".upcase+"_#{operation.name}").kernel
                          end
                      }
                      
                      kernel = BOAST::CKernel::new(:kernels => kernels)
                      
                      vars = [op, d, idim, n, bc, nx, ny, narr, x, y, a]
                      vars.push a_x if family == "s0s0"
                      vars.push a_y
                      p = BOAST::Procedure(function_name, vars){
                          case_args={}
                          operations.each_with_index{|op, i|
                              if suffix == "cost" then
                                  proc = kernels[i].cost_procedure
                              else
                                  proc = kernels[i].procedure
                              end
                              case_args[op] = lambda {
                              if family == "s0s0" then
                                  BOAST::pr proc.call(d, idim, n, bc, nx, ny, narr, x, y, a,a_x, a_y)
                              else
                                  BOAST::pr proc.call(d, idim, n, bc, nx, ny, narr, x, y, a, a_y)
                              end
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
          }
      }
#  }
  }#file

#    return kernel

#end 
