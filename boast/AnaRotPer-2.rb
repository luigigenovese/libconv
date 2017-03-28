require 'BOAST'
require 'rubygems'
require 'narray'
module BOAST

  def BOAST::analysis_free_ref
    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::FORTRAN)
    kernel = CKernel::new
    kernel.lang = BOAST::FORTRAN
    function_name = "analysis_free_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = -7
    dim_in_max = n*2+6
    dim_out_min = 0
    dim_out_max = n*2-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine analysis_free_ref(n,ndat,x,y)
  implicit none
  integer, intent(in) :: n,ndat
  real(kind=8), dimension(-7:2*n+6,ndat), intent(in) :: x
  real(kind=8), dimension(ndat,0:2*n-1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(kind=8) :: ci,di
  real(kind=8), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955d0, & 
       -0.00054213233180001068935d0, 0.031695087811525991431d0, & 
       0.0076074873249766081919d0, -0.14329423835127266284d0, & 
       -0.061273359067811077843d0, 0.48135965125905339159d0,  & 
       0.77718575169962802862d0,0.36444189483617893676d0, &
       -0.051945838107881800736d0,-0.027219029917103486322d0, &
       0.049137179673730286787d0,0.0038087520138944894631d0, &
       -0.014952258337062199118d0,-0.00030292051472413308126d0, &
       0.0018899503327676891843d0 /
  data cg  / -0.0018899503327676891843d0, &
       -0.00030292051472413308126d0, 0.014952258337062199118d0, &
       0.0038087520138944894631d0, -0.049137179673730286787d0, &
       -0.027219029917103486322d0, 0.051945838107881800736d0, &
       0.36444189483617893676d0, -0.77718575169962802862d0, &
       0.48135965125905339159d0, 0.061273359067811077843d0, &
       -0.14329423835127266284d0, -0.0076074873249766081919d0, &
       0.031695087811525991431d0, 0.00054213233180001068935d0, &
       -0.0033824159510050025955d0  /

  do j=1,ndat
     do i=0,n-1
        ci=0.d0
        di=0.d0
        do l=-7,8
           ci=ci+ch(l)*x(l+2*i,j)
           di=di+cg(l)*x(l+2*i,j)
        enddo
        y(j,i)=ci
        y(j,n+i)=di
     enddo
  enddo

END SUBROUTINE analysis_free_ref
EOF
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end
  def BOAST::analysis_per_ref
    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::FORTRAN)
    kernel = CKernel::new
    kernel.lang = BOAST::FORTRAN
    function_name = "analysis_per_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n*2-1
    dim_out_min = 0
    dim_out_max = n*2-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print <<EOF
subroutine analysis_per_ref(n,ndat,x,y)
  !use module_base
  implicit none
  integer, parameter :: wp=kind(1.d0)

!dee
!  integer :: iend_test,count_rate_test,count_max_test,istart_test

  integer, intent(in) :: n,ndat
  real(wp), dimension(0:2*n+1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:2*n+1), intent(out) :: y
  !local variables
  integer :: i,j,k,l
  real(wp) :: ci,di
  real(wp), dimension(-7:8) :: ch,cg
  !       Daubechy S16
  data ch  /  -0.0033824159510050025955_wp, & 
       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
       0.77718575169962802862_wp,0.36444189483617893676_wp, &
       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
       0.0018899503327676891843_wp /
  data cg  / -0.0018899503327676891843_wp, &
       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
       0.36444189483617893676_wp, -0.77718575169962802862_wp, &
       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
       -0.0033824159510050025955_wp  /

  integer mod_arr(-7:2*n+8)   
  real(wp) :: ci1,ci2,ci3,ci4,ci5,ci6,ci7,ci8
  real(wp) :: di1,di2,di3,di4,di5,di6,di7,di8

  !write(*,*) 'ana_rot_per executed'

  call fill_mod_arr(mod_arr,-7,2*n+8,2*n+2)

!dee
!call system_clock(istart_test,count_rate_test,count_max_test)

!$omp parallel default (private) shared(x,y,cg,ch,ndat,n,mod_arr)
!$omp do 

  do j=0,ndat/8-1
     do i=0,n
        ci1=0.e0_wp
        ci2=0.e0_wp
        ci3=0.e0_wp
        ci4=0.e0_wp
        ci5=0.e0_wp
        ci6=0.e0_wp
        ci7=0.e0_wp
        ci8=0.e0_wp

        di1=0.e0_wp
        di2=0.e0_wp
        di3=0.e0_wp
        di4=0.e0_wp
        di5=0.e0_wp
        di6=0.e0_wp
        di7=0.e0_wp
        di8=0.e0_wp

        do l=-7,8
           k= mod_arr(l+2*i)

           ci1=ci1+ch(l)*x(k,j*8+1)
           ci2=ci2+ch(l)*x(k,j*8+2)
           ci3=ci3+ch(l)*x(k,j*8+3)
           ci4=ci4+ch(l)*x(k,j*8+4)
           ci5=ci5+ch(l)*x(k,j*8+5)
           ci6=ci6+ch(l)*x(k,j*8+6)
           ci7=ci7+ch(l)*x(k,j*8+7)
           ci8=ci8+ch(l)*x(k,j*8+8)

           di1=di1+cg(l)*x(k,j*8+1)
           di2=di2+cg(l)*x(k,j*8+2)
           di3=di3+cg(l)*x(k,j*8+3)
           di4=di4+cg(l)*x(k,j*8+4)
           di5=di5+cg(l)*x(k,j*8+5)
           di6=di6+cg(l)*x(k,j*8+6)
           di7=di7+cg(l)*x(k,j*8+7)
           di8=di8+cg(l)*x(k,j*8+8)
        end do
        y(j*8+1,    i)=ci1
        y(j*8+2,    i)=ci2
        y(j*8+3,    i)=ci3
        y(j*8+4,    i)=ci4
        y(j*8+5,    i)=ci5
        y(j*8+6,    i)=ci6
        y(j*8+7,    i)=ci7
        y(j*8+8,    i)=ci8

        y(j*8+1,n+1+i)=di1
        y(j*8+2,n+1+i)=di2
        y(j*8+3,n+1+i)=di3
        y(j*8+4,n+1+i)=di4
        y(j*8+5,n+1+i)=di5
        y(j*8+6,n+1+i)=di6
        y(j*8+7,n+1+i)=di7
        y(j*8+8,n+1+i)=di8
     end do
  end do

  !$omp end do
  
  !$omp do  

  do j=(ndat/8)*8+1,ndat
     do i=0,n
        ci=0.e0_wp
        di=0.e0_wp
        do l=-7,8
           k= mod_arr(l+2*i)
           ci=ci+ch(l)*x(k    ,j)
           di=di+cg(l)*x(k    ,j)
        end do
        y(j,i)=ci
        y(j,n+1+i)=di
     end do
  end do
  !$omp end do

!$omp end parallel

!dee
!call system_clock(iend_test,count_rate_test,count_max_test)
!write(*,*) 'elapsed time on ana rot per',(iend_test-istart_test)/(1.d0*count_rate_test)

END SUBROUTINE analysis_per_ref
subroutine fill_mod_arr(arr,nleft,nright,n)
  implicit none
  integer,intent(in) :: nleft,nright,n
  integer,intent(out) :: arr(nleft:nright)
  integer :: i
  
  if (nleft >= -n) then
     do i=nleft,-1
        arr(i)=n+i
     end do
  else
     do i=nleft,-1
        arr(i)=modulo(i,n)
     end do
  endif
  
  do i=max(0,nleft),min(n-1,nright)
     arr(i)=i
  end do
  
  if (nright < 2*n) then
     do i=n,nright
        arr(i)=i-n
     end do
  else
     do i=n,nright
        arr(i)=modulo(i,n)
     end do
  endif
END SUBROUTINE fill_mod_arr



!subroutine analysis_per_ref(n,ndat,x,y)
!  !use module_base
!  implicit none
!  integer, parameter :: wp=kind(1.d0)
!  integer, intent(in) :: n,ndat
!  real(wp), dimension(0:2*n-1,ndat), intent(in) :: x
!  real(wp), dimension(ndat,0:2*n-1), intent(out) :: y
!  !local variables
!  integer :: i,j,k,l
!  real(wp) :: ci,di
!  real(wp), dimension(-7:8) :: ch,cg
!  !       Daubechy S16
!  data ch  /  -0.0033824159510050025955_wp, & 
!       -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
!       0.0076074873249766081919_wp, -0.14329423835127266284_wp, & 
!       -0.061273359067811077843_wp, 0.48135965125905339159_wp,  & 
!       0.77718575169962802862_wp,0.36444189483617893676_wp, &
!       -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
!       0.049137179673730286787_wp,0.0038087520138944894631_wp, &
!       -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
!       0.0018899503327676891843_wp /
!
!!cg(l)=(-1)**l ch(1-l)
!  data cg  / -0.0018899503327676891843_wp, &
!       -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
!       0.0038087520138944894631_wp, -0.049137179673730286787_wp, &
!       -0.027219029917103486322_wp, 0.051945838107881800736_wp, &
!        0.36444189483617893676_wp, -0.77718575169962802862_wp, &
!       0.48135965125905339159_wp, 0.061273359067811077843_wp, &
!       -0.14329423835127266284_wp, -0.0076074873249766081919_wp, &
!       0.031695087811525991431_wp, 0.00054213233180001068935_wp, &
!       -0.0033824159510050025955_wp  /
!
!  do j=1,ndat
!
!     do i=0,n-1
!        ci=0.e0_wp
!        di=0.e0_wp
!        do l=-7,8
!           k=modulo(l+2*i,2*n)
!            ci=ci+ch(l)*x(k    ,j)
!            di=di+cg(l)*x(k    ,j)
!        enddo
!        y(j,i)=ci
!        y(j,n+i)=di
!     enddo
!
!  enddo
!END SUBROUTINE analysis_per_ref
EOF
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end

  def BOAST::Analysis(filt, center, unroll, free=false )
    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang
    function_name = "analysis"
    if free then
      function_name += "_free"
    else 
      function_name += "_per"
    end
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    lowfil = Variable::new("lowfil",Int,{:constant => -center}) #-7 for daub16
    upfil = Variable::new("upfil",Int,{:constant => filt.length - center - 1}) # 8 for daub16

    if free then
      dim_in_min = lowfil #-7 in free BC
      dim_in_max = n*2 -2 + upfil #2*n+6 in free BC
    else
      dim_in_min = 0
      dim_in_max = n*2-1
    end
    dim_out_min = 0
    dim_out_max = n*2-1

    #potentially not needed
    #    if free then
    #      lowlimit=lowfil-1 #(dim_out_min-1)/2 #0 in periodic, -4 in free BC
    #      uplimit=n-2+upfil #(dim_out_max-1)/2 #n-1 in periodic, n+2 in free BC
    #    else
    #      lowlimit=0
    #      uplimit=n-2
    #    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    ci = [Variable::new("ci1",Real)]
    2.upto(unroll) { |index|
      ci.push(Variable::new("ci#{index}",Real))
    }
    di = [Variable::new("di1",Real)]
    2.upto(unroll) { |index|
      di.push(Variable::new("di#{index}",Real))
    }
    arr = ConstArray::new(filt.reverse,Real)

    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(-center, -center -1 + filt.length) ]})


    if BOAST::get_lang == C then
      @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      @@output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      @@output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end
    
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil]) {
      
      i.decl
      j.decl
      k.decl
      l.decl
      fil.decl
      ci.each{ |s| s.decl }
      di.each{ |s| s.decl }

      #tentative of defining the generic analysis loop into a common block
      #define the principal convolution operation
      inner_op = lambda { |ci,di,fBC|
        ci.each{ |s| (s === 0.0).print }
        di.each{ |s| (s === 0.0).print }
        For::new(l,lowfil,upfil,step: 2) {
          (k === l + i*2).print if fBC
          ci.each_index{ |ind|
            #(k === FuncCall::new( "modulo", l+ i*2, n*2)).print if not fBC
            #alternative form for modulo, uses only integer divisions
            (k ===  l+ i*2 - (  (l+ i*2 + n*10 )/(n*2) - 5) * n* 2 ).print if not fBC
            (ci[ind] === ci[ind] + fil[l]*x[k,j+ind]).print
            (di[ind] === di[ind] - fil[-l +1 ]*x[k,j+ind]).print
            if fBC then
              (ci[ind] === ci[ind] + fil[l+1]*x[k+1,j+ind]).print
              (di[ind] === di[ind] + fil[-l]*x[k+1,j+ind]).print
            else
              #(k === FuncCall::new( "modulo", l+ i*2+1, n*2)).print
              (k ===  l+ i*2+1 - (  (l+ i*2 +1 + n*10 )/(n*2) - 5) * n* 2 ).print if not fBC
              (ci[ind] === ci[ind] + fil[l+1]*x[k,j+ind]).print
              (di[ind] === di[ind] + fil[-l]*x[k,j+ind]).print
            end
          }
        }.unroll
        ci.each_index { |ind|
          (y[j+ind,i] === ci[ind]).print
        }
        di.each_index { |ind|
          (y[j+ind,n+i] === di[ind]).print
        }
      }


      analysis_d = lambda { |ci,di,js,ntot,fBC|
        unro=ci.length
        #external loop, with unrolling. the rest is excluded
        @@output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
        if unro > 1 then
          #forJ1 = For::new(j,js,ntot-(unro-1), step: unro) #do not forget to add the rest
          forJ1 = For::new(j,js,(ntot/unro)*unro, step: unro) #do not forget to add the rest
        else
          forJ1 = For::new(j,js,ntot)
        end
        #print the do part
        forJ1.print
        #in the case of free BC there is no start or end
        if fBC then
          For::new(i,0,n-1) {
            inner_op.call(ci,di,fBC)
          }.print
        else
          #left border
          For::new(i,0,-lowfil/2) {
            inner_op.call(ci,di,false)
          }.print
          #center
          For::new(i,-lowfil/2+1,n-1-upfil/2) {
            inner_op.call(ci,di,true)
          }.print
          #right border
          For::new(i,n-upfil/2,n-1) {
            inner_op.call(ci,di,false)
          }.print
        end      
        #end do for the external loop
        forJ1.close
        @@output.print("!$omp end do\n") if BOAST::get_lang == BOAST::FORTRAN
      }

      #body of the routine
      @@output.print("!$omp parallel default (private) shared(x,y,fil,ndat,n)\n") if BOAST::get_lang == BOAST::FORTRAN
      analysis_d.call(ci,di,1,ndat,free)
      #remaining part after unrolling
      if unroll>1 then
        analysis_d.call([ci[0]],[di[0]],(ndat/unroll)*unroll+1,ndat,free)
      end
      @@output.print("!$omp end parallel\n") if BOAST::get_lang == BOAST::FORTRAN
      }
    p.print
    kernel.procedure = p
    return kernel
  end
end


##FILTER = ["0.0018899503327676891843",
##          "-0.00030292051472413308126",
##          "-0.014952258337062199118",
##          "0.0038087520138944894631",
##          "0.049137179673730286787",
##          "-0.027219029917103486322",
##          "-0.051945838107881800736",
##          "0.36444189483617893676",
##          "0.77718575169962802862",
##          "0.48135965125905339159",
##          "-0.061273359067811077843",
##          "-0.14329423835127266284",
##          "0.0076074873249766081919",
##          "0.031695087811525991431",
##          "-0.00054213233180001068935",
##          "-0.0033824159510050025955"]
##
##n1 = 124
##n2 = 132
##n3 = 130
##input = NArray.float(n1+14,n2,n3).random
##output_ref = NArray.float(n2,n3,n1)
##output = NArray.float(n2,n3,n1)
##epsilon = 10e-15
##BOAST::set_lang( BOAST::FORTRAN )
##k = BOAST::analysis_free_ref
##stats = k.run(n1/2, n2*n3, input, output_ref)
##puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
##
##(1..12).each{ |unroll|
##  k = BOAST::Analysis(FILTER,7,unroll,true)
##  #k.print
##  #k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fbounds-check",:LDFLAGS => "-lgfortran"})
##  k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -openmp",:LDFLAGS => "-openmp"})
##  stats = k.run(n1/2, n2*n3, input, output)
##  stats = k.run(n1/2, n2*n3, input, output)
##  diff = (output_ref - output).abs
##  diff.each { |elem|
##    puts "Warning: residue too big: #{elem}" if elem > epsilon
##  }
##  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
##}
##BOAST::set_lang( BOAST::C )
##(1..12).each{ |unroll|
##  k = BOAST::Analysis(FILTER,7,unroll,true)
##
##  stats = k.run(n1/2, n2*n3, input, output)
##  stats = k.run(n1/2, n2*n3, input, output)
##  diff = (output_ref - output).abs
##  diff.each { |elem|
##    puts "Warning: residue too big: #{elem}" if elem > epsilon
##  }
##  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
##}
##
##n1 = 124
##n2 = 132
##n3 = 130
##input = NArray.float(n1,n2,n3).random
##output_ref = NArray.float(n1,n2,n3)
##output = NArray.float(n1,n2,n3)
##epsilon = 10e-15
##BOAST::set_lang( BOAST::FORTRAN )
##k = BOAST::analysis_per_ref
##k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -openmp",:LDFLAGS => "-openmp",:LD => "ifort"})
##stats = k.run(n1/2-1, n2*n3, input, output_ref)
##puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
##
##(1..14).each{ |unroll|
##  #k = BOAST::analysis_per_ref
##  k = BOAST::Analysis(FILTER,7,unroll,false)
##  #k.print
##  #k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fbounds-check",:LDFLAGS => "-lgfortran"})
##  k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -openmp",:LDFLAGS => "-openmp",:LD => "ifort"})
##  #k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -g -C",:LDFLAGS => "",:LD => "ifort"})
##
##  stats = k.run(n1/2, n2*n3, input, output)
##  diff = (output_ref - output).abs
##  diff.each { |elem|
##    puts "Warning: residue too big: #{elem}" if elem > epsilon
##  }
##  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
##}
