#This generates high level brokers for libconv, as well as header files to use for compilation against libconv, and a Makefile.

  D_DESC = "Number of dimensions"
  OP_DESC = "Operator. see #ops"
  OP_DESC_HELPER = "Operator. see #ops. If set to negative values, query mode: ny will return the output dimensions, cost and alignment values for -op"
  TP_DESC = "Precision (4 for single, 8 for double)"
  IDIM_DESC = "Treated dimension"
  NARR_DESC = "Number of consecutive applications of the operator on x and y"
  SF_DESC = "Array of one resolution level by wavelet transform"
  SW_DESC = "Array of two resolution levels, of dimension 2**d*product(nsw)"
  X_DESC = "Input array - size nx"
  Y_DESC = "Output array - size ny"
  N_DESC = ""
  BC_DESC = "Boundary condition. See #bcs"
  BCS_DESC = "Boundary conditions. See #bcs"
  NX_DESC = "Input sizes."
  NY_DESC = "Output sizes. If query, has to be allocated to d + 2, and will return dimensions for y, cost, and alignment necessary"
  A_DESC = "Global multiplier"
  AX_DESC = "Added multiplier on x"
  AY_DESC = "Added multiplier on y - 1.0 means aggregation"
  DOTIN_DESC = "Scalar product"
  COST_DESC = "Estimated cost in flops of the operation"
  ALIGN_DESC = "Alignment needed for arrays"
  DIMS_DESC = "Dimensions for output array"
  WORK_DESC = "Work array - Same dimensions as y"
  DOTINS_DESC = "Scalar product for each dimension"
  AS_DESC = "Global multiplier for each dimension"

  def get_comment(type, util, precision)
    case precision
      when 4
        precision_fullname = "single"
      when 8
        precision_fullname = "double"
      else
        precision_fullname = "mixed"
    end
    case type
      when "1d"
        full_comment = ["normal wavelet transform in only one direction, #{precision_fullname} precision", "y := a * fil \\|X_i\\| x(x_1,...,x_i,...,x_d,j)  + a_y * y"]
      when "1ds"
        full_comment = ["sum of multi dimensional convolutions in the different direction, #{precision_fullname} precision","y = sum_i (a_i * fil \\|X_i\\| x(x_1,...,x_i,...,x_d,i))+ a_x * x + a_y * y"]
      when "md"
        full_comment = ["application of a separable multi-dimensional convolution on the input array, #{precision_fullname} precision","y := a * fil \\|X\\| x "]
    end
    if not util
      comment = full_comment
    else
      comment = []
      comment.push(util + " helper for " + full_comment[0])
    end
    return comment
  end

  def print_filters
      #print filters
      @wavelet_families.each{ |wav_fam|
        ["MF", "LP", "D1", "D2"].each{ |fil|
          filter =BOAST::const_get(wav_fam+"_"+fil)
          filt = GenericConvolution::ConvolutionFilter.new(wav_fam+'_'+fil, filter, filter.length/2-1)
          filt.decl_filters
        }
      }
  end
  
  def print_bcs(f)
    conds = {"LIBCONV_BC_PERIODIC" => GenericConvolution::BC::PERIODIC, 
             "LIBCONV_BC_ZERO" => GenericConvolution::BC::GROW, 
             "LIBCONV_BC_FREE_GROW" => GenericConvolution::BC::SHRINK, 
             "LIBCONV_BC_FREE_SHRINK" => GenericConvolution::BC::FREE, 
             "LIBCONV_BC_INTERVAL" => GenericConvolution::BC::NPERIODIC}
    if(BOAST::get_lang == BOAST::C) then
      f.puts "enum bcs{"
      conds.each_with_index{ |(name, val), i|
        f.print "#{name}=#{val}"
        f.print "," if i != (conds.length() -1)
      }
      f.puts "};"
    else
      id=0
      conds.each_with_index{ |(name, val), i|
      f.print "integer :: " if id%3 == 0
      f.print "#{name}"
      if id%3 != 2 and  i != (conds.length() -1)
        f.print ", "
      else
        f.puts ""
      end
      id+=1
      }
      id=0
      conds.each{ |name, val|
        f.puts "parameter(#{name}=#{val})"
        id+=1
      }
    end
    
  end

def print_ops(f)
    if(BOAST::get_lang == BOAST::C) then
      id=1
      f.puts "enum ops{"
      @wavelet_families.each{ |wav_fam|
        @operations.each{ |op|
          f.print "#{wav_fam}_#{op} = #{id}"
          f.print "," if op != @operations.last
          id+=1
        }
      }
      f.puts "};"
    else
      id=0
      @wavelet_families.each{ |wav_fam|
#       f.puts "module brokers"
        @operations.each{ |op|
          f.print "integer :: " if id%3 == 0
          f.print "#{wav_fam}_#{op}"
          if id%3 != 2 and  !(op == @operations.last and wav_fam == @wavelet_families.last)
            f.print ", "
          else
            f.puts ""
          end
          id+=1
        }
      }
      id=1
      @wavelet_families.each{ |wav_fam|
        @operations.each{ |op|
          f.puts "parameter(#{wav_fam}_#{op}=#{id})"
          id+=1
        }
      }
    end
end

def print_headers(fold)
  #generate both Fortran and C header files every time
  lang = BOAST::get_lang
  [BOAST::C, BOAST::FORTRAN].each { |l|
    if l == BOAST::C then
      header = "libconv.h"
    else
      header = "libconvf.h"
    end 
    File::open(fold+header,"w") { |f|
    set_output(f)
    BOAST::set_lang(l)
    print_ops(f)
    print_filters
    print_bcs(f)
    BOAST::set_lang(lang)
    }
  }

end

def print_broker_1d(f)
  d = BOAST::Int("d", :dir => :in, :reference => 1, :comment => D_DESC)
  idim = BOAST::Int("idim", :dir => :in, :reference => 1, :comment => IDIM_DESC)
  n = BOAST::Int("n", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => N_DESC)
  bc = BOAST::Int("bc", :dir => :in, :reference => 1, :comment => BC_DESC)
  nx = BOAST::Int("nx", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NX_DESC)
  ny = BOAST::Int("ny", :dir => :inout,:dim => [ BOAST::Dim(0, d + 1) ], :comment => NY_DESC)
  narr = BOAST::Int("narr", :dir => :in, :reference => 1, :comment => NARR_DESC)
  kernels=[]
  @families.each{ |family|
    ops=[]
    @wavelet_families.each{ |wav_fam|
      case family
        when "s0s0"
          ops.push(const_get("#{wav_fam}_MF")) if @operations.include? "MF"
          ops.push(const_get("#{wav_fam}_IMF")) if @operations.include? "IMF"
        when "s0s0_dot"
          ops.push(const_get("#{wav_fam}_D1")) if @operations.include? "D1"
          ops.push(const_get("#{wav_fam}_D2")) if @operations.include? "D2"
        when "s0s1"
          ops.push(const_get("#{wav_fam}_DWT")) if @operations.include? "DWT"
          ops.push(const_get("#{wav_fam}_RTOS1")) if @operations.include? "RTOS1"
        when "s1s0"
          ops.push(const_get("#{wav_fam}_IDWT")) if @operations.include? "IDWT"
          ops.push(const_get("#{wav_fam}_S1TOR")) if @operations.include? "S1TOR"
      end
    }
    next if ops.length ==0
    
    @precisions.each{ |precision|
      case precision
        when 4
          precision_name = "s"
        when 8
          precision_name = "d"
        else
          puts "unknown precision :"+precision
      end
      generate_kernel = Proc.new  { |util|
        BOAST.push_env( default_real_size: precision ){
          @dimensions.each{ |dimension|
            kernels=[]
            func_name = "#{precision_name}_#{family}_#{dimension}"
            func_name += "_#{util}" if util
            function_name = func_name
            function_name += "_" if BOAST::get_lang == BOAST::C
            op = BOAST::Int("op", :dir => :in, :reference => 1, :comment => util ? OP_DESC : OP_DESC_HELPER )
            a = BOAST::Real("a", :dir => :in, :reference => 1, :comment => A_DESC)
            a_x = BOAST::Real("a_x", :dir => :in, :reference => 1, :comment => AX_DESC)
            a_y = BOAST::Real("a_y", :dir => :in, :reference => 1, :comment => AY_DESC)
            dot_in = Real("dot_in",:dir => :inout, :dim => [ BOAST::Dim(1)], :comment => DOTIN_DESC)
            x = BOAST::Real("x", :dir => :in, :dim => [ BOAST::Dim()], :comment => X_DESC)
            y = BOAST::Real("y", :dir => :inout, :dim => [ BOAST::Dim(2)], :comment => Y_DESC )
            cost = BOAST::Int("cost", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => COST_DESC)
            alignment = BOAST::Int("alignment", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => ALIGN_DESC)
            dims = BOAST::Int("dims", :dir => :out, :dim => [ BOAST::Dim(0, d - 1)], :comment => DIMS_DESC)
            temp_util = BOAST::Int("tmp", :dim => [ BOAST::Dim(1)])
            ops.each{ |operation|
              if util then 
                kernels.push const_get("#{precision_name}".upcase+"_#{operation.name}")
              else 
                kernels.push const_get("#{precision_name}".upcase+"_#{operation.name}").kernel
              end
            }

            kernel = BOAST::CKernel::new(:kernels => kernels)
            get_args=lambda{|util|
              vars = [op, d, idim, n, bc]
              vars+=[nx, ny, narr] if not util or util == "cost"
              vars+=[x, y] if not util
              vars.push a
              vars.push a_x if family == "s0s0" or family == "s0s0_dot"
              vars.push a_y
              vars.push dot_in if family == "s0s0_dot" and not util
              vars.push cost if util == "cost"
              vars.push dims if util == "dims"
              vars.push alignment if util == "align"
              return vars
            }

            get_proc=lambda{|util|
              func = func_name
              func += "_" + util if util
              func += "_" if BOAST::get_lang == BOAST::C
              kern = const_get(func.upcase)
              return kern.procedure
            }
            vars = get_args.call(util)
            p = BOAST::Procedure(function_name, vars, :comment => get_comment("1d", util, precision)){
              decl temp_util unless util
              BOAST::pr temp_util === 0 unless util
              pr_case=lambda{
                case_args={}
                vars.shift
                ops.each_with_index{|op, i|
                  extra=nil
                  if util == "cost" 
                    proc = kernels[i].cost_procedure
                    extra = cost
                  elsif  util == "dims"
                    proc = kernels[i].dims_procedure
                  elsif  util == "align"
                    proc = kernels[i].align_procedure
                  else
                    proc = kernels[i].procedure
                  end
                  case_args[op] = lambda {
                    BOAST::pr proc.call(*vars)
                  }
                }
                BOAST::pr BOAST::Case( op, case_args)
              }

              #TODO don't generate separate utils anymore ?
              if not util
                args_cost = get_args.call("cost")
                args_dims = get_args.call("dims")
                args_align = get_args.call("align")
                BOAST::pr BOAST::If(op >= 0 => lambda{
                  pr_case.call()
                },else: lambda{
                  args_dims[0]=-op
                  args_dims[-1]=ny
                  BOAST::pr get_proc.call("dims").call(*args_dims)
                  args_cost[0]=-op
                  index = BOAST::get_lang==BOAST::FORTRAN ? 1 : 0
                  args_cost[-1]=ny[d]
                  BOAST::pr temp_util[index]===ny[d]
                  BOAST::pr get_proc.call("cost").call(*args_cost)
                  BOAST::pr ny[d]===ny[d]+temp_util[index]
                  args_align[0]=-op
                  args_align[-1]=ny[d+1]
                  BOAST::pr temp_util[index]===ny[d+1]
                  BOAST::pr get_proc.call("align").call(*args_align)
                  BOAST::pr ny[d+1]===BOAST::Max(temp_util[index], ny[d+1])
                }
              )
              else
                pr_case.call()
              end
            }
            BOAST::pr p
            kernel.procedure = p
            f.puts kernel
            LibConv.const_set(function_name.upcase, kernel)
          }
        }
      }
      generate_kernel.call("cost")
      generate_kernel.call("dims")
      generate_kernel.call("align")
      generate_kernel.call()
    }
  }
end

#multidimensionnal versions of convolutions
def print_brokers_md(f)
  d = BOAST::Int("d", :dir => :in, :reference => 1, :comment => D_DESC)
  idim = BOAST::Int("idim", :dir => :in, :reference => 1, :comment => IDIM_DESC)
  n = BOAST::Int("n", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => N_DESC)
  bcs = BOAST::Int("bcs", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => BCS_DESC)
  nx = BOAST::Int("nx", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NX_DESC)
  ny = BOAST::Int("ny", :dir => :inout,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NY_DESC)
  narr = BOAST::Int("narr", :dir => :in, :reference => 1, :comment => NARR_DESC)
  kernels=[]
  @families.each{ |family|
    next if family == "s0s0_dot"
    @precisions.each{ |precision|
      case precision
        when 4
          precision_name = "s"
        when 8
          precision_name = "d"
        else
          puts "unknown precision :"+precision
      end
      generate_kernel = Proc.new  { |util|
        BOAST.push_env( default_real_size: precision ){
            kernels=[]
            function_name = "#{precision_name}_#{family}_md"
            function_name += "_#{util}" if util
            function_name += "_" if BOAST::get_lang == BOAST::C
            op = BOAST::Int("op", :dir => :in, :reference => 1, :comment => util ? OP_DESC : OP_DESC_HELPER )
            a = BOAST::Real("a", :dir => :in, :reference => 1, :comment => A_DESC)
            x = BOAST::Real("x", :dir => :in, :dim => [ BOAST::Dim()], :comment => X_DESC)
            y = BOAST::Real("y", :dir => :inout, :dim => [ BOAST::Dim()], :comment => Y_DESC)
            work = BOAST::Real("work", :dir => :inout, :dim => [ BOAST::Dim()], :comment => WORK_DESC )
            cost = BOAST::Int("cost", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => COST_DESC)
            temp_util = BOAST::Int("tmp", :dim => [ BOAST::Dim(1)])
            alignment = BOAST::Int("alignment", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => ALIGN_DESC)
            dims = BOAST::Int("dims", :dir => :out, :dim => [ BOAST::Dim(0, d - 1)], :comment => DIMS_DESC)
            idim = BOAST::Int("idim")
            kern1d_name = "#{precision_name}_#{family}_1D"
            kern1d_name += "_#{util}" if util
            kern1d_name += "_" if BOAST::get_lang == BOAST::C
            kern1d=const_get(kern1d_name.upcase)


            kernel = BOAST::CKernel::new(:kernels => kern1d)
            vars = [op, d, n, bcs]
            idim_index=2 #will be inserted later
            bcs_index=vars.size
            vars+=[nx, ny, narr] if not util or util == "cost"
            vars+=[x, y, work] if not util
            x_index=vars.size-2
            y_index=vars.size-1
            vars.push a
            vars.push cost if util == "cost"
            vars.push dims if util == "dims"
            vars.push alignment if util == "align"
            if util then
              a_index=-2
            else
              a_index=-1
            end
            p = BOAST::Procedure(function_name, vars, :comment => get_comment("md", util, precision)){
              decl temp_util if util == "align" or util == "cost"
              vars1=vars
              vars1-=[work]
              vars1.insert(idim_index,idim)
              vars1.insert(a_index,BOAST::Real("a_x", :constant => 0.0, :reference => 1)) if family == "s0s0"
              vars1.insert(a_index,BOAST::Real("a_y", :constant => 0.0, :reference => 1))                
                printcall = lambda{|idim,x,y|
                  vars1[idim_index]=idim-1
                  vars1[bcs_index]=bcs[idim-1]
                  vars1[x_index]=x if not util
                  vars1[y_index]=y if not util
                  vars1[vars1.size-1]=temp_util if util == "align" or util == "cost"
                  BOAST::pr kern1d.procedure.call(*vars1)
                  if util == "align"
                    BOAST::pr alignment===BOAST::Max(alignment, temp_util)
                  elsif  util == "cost"
                    BOAST::pr cost===cost+temp_util
                  end
                }
                if util == "align"
                    BOAST::pr alignment===0
                    BOAST::pr temp_util===0
                  elsif  util == "cost"
                    BOAST::pr cost===0
                    BOAST::pr temp_util===0
                  end
                index = BOAST::get_lang==BOAST::FORTRAN ? 1 : 0
                BOAST::pr BOAST::If( 2*(d/2) == d => lambda {
                  printcall.call(1,x,work)
                  BOAST::pr BOAST::For(idim,1,((d/2)-1)){
                    printcall.call(2*idim,work,y)
                    printcall.call(2*idim+1,y,work)
                  }
                  printcall.call(d,work,y)
                }, else: lambda{
                  printcall.call(1,x,y)
                  BOAST::pr BOAST::For(idim,1,(d/2)){
                    printcall.call(2*idim,y,work)
                    printcall.call(2*idim+1,work,y)
                  }
                })
            }
            BOAST::pr p
            kernel.procedure = p
            f.puts kernel
            LibConv.const_set(function_name.upcase, kernel)
        }
      }
      generate_kernel.call()
      generate_kernel.call("cost")
      generate_kernel.call("dims")
      generate_kernel.call("align")
    }
  }
end

#multidimensionnal versions of convolutions without work array
def print_brokers_1ds(f)
  d = BOAST::Int("d", :dir => :in, :reference => 1, :comment => D_DESC)
  idim = BOAST::Int("idim", :dir => :in, :reference => 1, :comment => IDIM_DESC)
  n = BOAST::Int("n", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => N_DESC)
  bcs = BOAST::Int("bcs", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => BCS_DESC)
  nx = BOAST::Int("nx", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NX_DESC)
  ny = BOAST::Int("ny", :dir => :inout,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NY_DESC)
  narr = BOAST::Int("narr", :dir => :in, :reference => 1, :comment => NARR_DESC)
  kernels=[]
  @families.each{ |family|
    @precisions.each{ |precision|
      case precision
        when 4
          precision_name = "s"
        when 8
          precision_name = "d"
        else
          puts "unknown precision :"+precision
      end
      generate_kernel = Proc.new  { |util|
        BOAST.push_env( default_real_size: precision ){
#          @dimensions.each{ |dimension|
            kernels=[]
            function_name = "#{precision_name}_#{family}_1ds"
            function_name += "_#{util}" if util
            function_name += "_" if BOAST::get_lang == BOAST::C
            op = BOAST::Int("op", :dir => :in, :reference => 1, :comment => util ? OP_DESC : OP_DESC_HELPER )
            as = BOAST::Real("as", :dir => :in, :dim => [ BOAST::Dim(0, d - 1)], :comment => AS_DESC )
            a_x = BOAST::Real("a_x", :dir => :in, :reference => 1, :comment => AX_DESC )
            a_y = BOAST::Real("a_y", :dir => :in, :reference => 1, :comment => AY_DESC )
            dot_ins = Real("dot_ins",:dir => :inout, :dim => [ BOAST::Dim(0, d - 1)], :comment => DOTINS_DESC)
            x = BOAST::Real("x", :dir => :in, :dim => [ BOAST::Dim()], :comment => X_DESC)
            y = BOAST::Real("y", :dir => :inout, :dim => [ BOAST::Dim()], :comment => Y_DESC)
            cost = BOAST::Int("cost", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => COST_DESC)
            temp_util = BOAST::Int("tmp", :dim => [ BOAST::Dim(1)])
            alignment = BOAST::Int("alignment", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => ALIGN_DESC)
            dims = BOAST::Int("dims", :dir => :out, :dim => [ BOAST::Dim(0, d - 1)], :comment => DIMS_DESC)
            idim = BOAST::Int("idim")
            kern1d_name = "#{precision_name}_#{family}_1D"
            kern1d_name += "_#{util}" if util
            kern1d_name += "_" if BOAST::get_lang == BOAST::C
            kern1d=const_get(kern1d_name.upcase)


            kernel = BOAST::CKernel::new(:kernels => kern1d)
            vars = [op, d, n, bcs]
            idim_index=2 #will be inserted later
            bcs_index=vars.size
            vars+=[nx, ny, narr] if not util or util == "cost"
            vars+=[x, y] if not util
            x_index=vars.size-1
            y_index=vars.size
            vars.push as
            as_index=vars.size
            vars.push a_x if family == "s0s0" or family == "s0s0_dot"
            vars.push a_y
            vars.push dot_ins if family == "s0s0_dot" and not util
            dot_ins_index=vars.size if family == "s0s0_dot" and not util
            vars.push cost if util == "cost"
            vars.push dims if util == "dims"
            vars.push alignment if util == "align"
            
            p = BOAST::Procedure(function_name, vars, :comment => get_comment("1d", util, precision)){
              decl temp_util if util == "align" or util == "cost"
              vars1=vars
              vars1.insert(idim_index,idim)
                
                printcall = lambda{|idim,x,y|
                  vars1[idim_index]=idim
                  vars1[bcs_index]=bcs[idim]
                  vars1[as_index]=as[idim]
                  vars1[dot_ins_index]=dot_ins[idim] if family == "s0s0_dot" and not util
                  vars1[x_index]=x if not util
                  vars1[y_index]=y if not util
                  vars1[vars1.size-1]=temp_util if util == "align" or util == "cost"
                  BOAST::pr kern1d.procedure.call(*vars1)
                  if util == "align"
                    BOAST::pr alignment===BOAST::Max(alignment, temp_util)
                  elsif  util == "cost"
                    BOAST::pr cost===cost+temp_util
                  end
                }
                if util == "align"
                    BOAST::pr alignment===0
                    BOAST::pr temp_util===0
                  elsif  util == "cost"
                    BOAST::pr cost===0
                    BOAST::pr temp_util===0
                  end
                
                  BOAST::pr BOAST::For(idim,0,d-1){
                    printcall.call(idim,x,y)
                  }
            }
            BOAST::pr p
            kernel.procedure = p
            f.puts kernel
            LibConv.const_set(function_name.upcase, kernel)
#          }
        }
      }
      generate_kernel.call()
      generate_kernel.call("cost")
      generate_kernel.call("dims")
      generate_kernel.call("align")
    }
  }
end

#brokers to handle precisions
def print_entrypoints(f)
  d = BOAST::Int("d", :dir => :in, :reference => 1, :comment => D_DESC)
  idim = BOAST::Int("idim", :dir => :in, :reference => 1, :comment => IDIM_DESC)
  n = BOAST::Int("n", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => N_DESC)
  bc = BOAST::Int("bc", :dir => :in, :reference => 1, :comment => BC_DESC)
  bcs = BOAST::Int("bcs", :dir => :in, :reference => 1,:dim => [ BOAST::Dim(0, d - 1) ], :comment => BCS_DESC)
  nx = BOAST::Int("nx", :dir => :in,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NX_DESC)
  ny = BOAST::Int("ny", :dir => :inout,:dim => [ BOAST::Dim(0, d - 1) ], :comment => NY_DESC)
  narr = BOAST::Int("narr", :dir => :in, :reference => 1, :comment => NARR_DESC)
  kernels=[]
  @families.each{ |family|
    generate_kernel = Proc.new  { |util|
      ["1d", "md", "1ds"].each{ |dimension|
        next if family == "s0s0_dot" and dimension == "md"
        kernels=[]
        function_name = "#{family}_#{dimension}"
        function_name += "_#{util}" if util
        function_name += "_" if BOAST::get_lang == BOAST::C

        op = BOAST::Int("op", :dir => :in, :reference => 1, :comment => util ? OP_DESC : OP_DESC_HELPER )
        prec = BOAST::Int("prec", :dir => :in, :reference => 1, :comment => TP_DESC)
        a = BOAST::Real("a", :dir => :in, :reference => 1, :comment => A_DESC)
        a_x = BOAST::Real("a_x", :dir => :in, :reference => 1, :comment => AX_DESC )
        a_y = BOAST::Real("a_y", :dir => :in, :reference => 1, :comment => AY_DESC )
        dot_in = Real("dot_in",:dir => :inout, :dim => [ BOAST::Dim(1)], :comment => DOTIN_DESC)
        libconv_generic_kind = BOAST::Int("libconv_generic_kind", :constant => 8, :reference => 1)
        x = BOAST::Real("x", :dir => :in, :size=>libconv_generic_kind.name, :dim => [ BOAST::Dim()], :comment => X_DESC)
        y=BOAST::Real("y", :dir => :inout, :size=>libconv_generic_kind.name, :dim => [ BOAST::Dim()], :comment => Y_DESC)
        work = BOAST::Real("work", :dir => :inout, :dim => [ BOAST::Dim()], :comment => WORK_DESC)
        alignment = BOAST::Int("alignment", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => ALIGN_DESC)
        cost = BOAST::Int("cost", :dir => :out, :dim => [ BOAST::Dim(1)], :comment => COST_DESC)
        dims = BOAST::Int("dims", :dir => :out, :dim => [ BOAST::Dim(0, d - 1)], :comment => DIMS_DESC)
        @precisions.each{ |precision|
          case precision
            when 4
              precision_name = "s"
            when 8
              precision_name = "d"
            else
              puts "unknown precision :"+precision
          end
          subname = "#{precision_name}_#{family}_#{dimension}"
          subname += "_#{util}" if util
          subname += "_" if BOAST::get_lang == BOAST::C
          kernels.push const_get(subname.upcase)
        }
        kernel = BOAST::CKernel::new(:kernels => kernels)
        vars = [op, prec, d]
        vars.push idim if dimension == "1d"
        vars.push n
        if dimension == "1d" then
          vars.push bc
        else
          vars.push bcs
        end 
        vars+=[nx, ny, narr] if not util or util == "cost"
#        vars.push libconv_generic_kind
        vars+=[x, y] if not util
        vars.push work if not util and dimension == "md"
        vars.push a
        vars.push a_x if dimension != "md" and (family == "s0s0" or family == "s0s0_dot")
        vars.push a_y if dimension != "md"
        vars.push dot_in if family == "s0s0_dot" and not util
        vars.push cost if util == "cost"
        vars.push dims if util == "dims"
        vars.push alignment if util == "align"
        p = BOAST::Procedure(function_name, vars,:constants => [libconv_generic_kind], :comment => get_comment(dimension, util, "mixed")){
#          decl libconv_generic_kind
          vars.delete_at(1)
          case_args={}
          @precisions.each_with_index{|precision, i|
            extra=nil
            proc = kernels[i].procedure
            case_args[precision] = lambda {
              BOAST::pr proc.call(*vars)
            }
          }
          BOAST::pr BOAST::Case( prec, case_args)
        }
        BOAST::pr p
        kernel.procedure = p
        f.puts kernel
      }
    }
    generate_kernel.call()
    generate_kernel.call("cost")
    generate_kernel.call("dims")
    generate_kernel.call("align")
  }
end


# write Makefile.am to foldername
def print_makefile(extra_files)
    compiler_options = BOAST::get_compiler_options

    makelines="""
    lib_LIBRARIES = libconv.a

    #temporary compiling line for gfortran
    AM_FCFLAGS = -I. #{compiler_options[:FCFLAGS]} -march=#{BOAST::get_model} #{BOAST::get_openmp_flags[compiler_options[:FC]]}

    #temporary compiling line for gcc
    AM_CFLAGS = -I. #{compiler_options[:CFLAGS]} -march=#{BOAST::get_model} #{BOAST::get_openmp_flags[compiler_options[:CC]]}
    
    libconv_a_SOURCES = """

    LibConv::kernels.each{|f| makelines+=f.to_str+"\\\n"}

    makelines+=extra_files.join("\\\n")

#    makelines+="libconv_a_OBJECTS = $(libconv_a_SOURCES:.f90=.o)\n"

#    makelines+="%.o: %.f90
#	    gfortran $(AM_FCFLAGS) -c $< \n
#    libconv.a: $(libconv_a_OBJECTS)
#	    gfortran $(AM_FCFLAGS) -c -o $@ $^"
foldername=""
  if not LibConv.from_cache then
    foldername = LibConv::foldername
  end
    File::open("#{foldername}Makefile.am","w") {|f|
      f.puts makelines
    }
end

def set_ids()
  id=1
  @wavelet_families.each{ |wav_fam|
    @operations.each{ |op|
    const_set("#{wav_fam}_#{op}", BOAST::Int("#{wav_fam}_#{op}", :constant => id) )
    id +=1
    }
  }
end


module LibConv

  set_ids()
  foldername=""
  if not LibConv.from_cache then
    foldername = LibConv::foldername
  end
  print_headers(foldername)
  suffix=".f90"
  if BOAST::get_lang == BOAST::C then
    suffix=".c"
  end

  brokers_filename = "brokers#{suffix}"
  File::open(foldername+brokers_filename,"w") { |f|
    set_output(f)
    if BOAST::get_lang == BOAST::C then
      f.puts  "#include <stdint.h>"
    end
    print_broker_1d(f)
    print_brokers_md(f)
    print_brokers_1ds(f)
   }

  entrypoints_filename = "libconv#{suffix}"
  if BOAST::get_lang != BOAST::C then
    File::open(foldername+entrypoints_filename,"w") { |f|
      set_output(f)
      print_entrypoints(f)
    }
  end
  print_makefile([brokers_filename, entrypoints_filename])



#    return kernel
end


#  @families.each{ |family|
#      @dimensions.each{ |dimension|
#          @precisions.each{ |precision|
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
