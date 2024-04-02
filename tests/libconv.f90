!> test program for the convolution library
program libconv
  implicit none
  include 'libconvf.h'
  integer :: flops, op
  real(kind=8) :: tol
  integer, dimension(3) :: ndims, bc, ibc, nout, nout2
  real(kind=8), dimension(:,:,:), allocatable :: x, y, tmp1, tmp2, x2
  integer(kind=8) :: size
  real(kind=4) :: stol

  !to be read from command line (TODO)
  ndims = [35, 180, 80] 
  bc = [0, 0, 0] 
  ibc = -bc

  size=product(int(ndims,kind=8))

  op=SYM8_MF

  !call s0s0_3d_info(SYM8_IMF, bc, ndims, nout, flops)    
  call op_3d_info(8, op, bc, ndims, nout, flops)
  ! call s0s0_3d_info(op, bc, ndims, nout, flops)
  ! call s0s0_3d_info(SYM8_IMF, -bc, nout, nout2, flops)    
  allocate(x   (ndims(1),ndims(2),ndims(3)))
  allocate(tmp1(nout(1),ndims(2),ndims(3)))
  allocate(tmp2(nout(1),nout(2),ndims(3)))
  allocate(y   (nout(1),nout(2),nout(3)))
  allocate(x2  (ndims(1),ndims(2),ndims(3)))

  x=1.0d0
  y=0.0d0
       

  ! call s0s0_3d(op,bc,ndims,x,y,tmp1,tmp2)
  ! call s0s0_3d(SYM8_IMF,ibc,nout,y,x2,tmp2,tmp1)

  call op_3d(8,SYM8_MF , bc,ndims,1.0d0, 0.0d0, 0.0d0,x,y,tmp1,tmp2)
  call op_3d(8,SYM8_IMF,ibc,nout,1.0d0, 0.0d0, 0.0d0,y,x2,tmp2,tmp1)
  call allclose(size,x,x2,tol)

  print *,'tolerance',tol 
  
  tmp1=0.d0
  tmp2=0.0d0
  y=0.0d0
  x=0.0d0
  call to_simple(size,x2,x)
  x2=0.0d0
  call op_3d(4,SYM8_MF , bc,ndims,1.0e0, 0.0e0, 0.0e0,x,y,tmp1,tmp2)
  call op_3d(4,SYM8_IMF,ibc,nout,1.0e0, 0.0e0, 0.0e0,y,x2,tmp2,tmp1)
  call allclose_s(size,x,x2,stol)

  print *,'tolerance (single precision)',stol

  deallocate(x2)
  deallocate(y)
  deallocate(tmp2)
  deallocate(tmp1)
  deallocate(x)

end program libconv

subroutine to_simple(size, src, dest)
  implicit none
  integer(kind=8), intent(in) :: size
  real(kind=8), dimension(size), intent(in) :: src
  real(kind=4), dimension(size), intent(out) :: dest
  !local variables
  integer(kind=8) :: i

  do i=1,size
     dest(i) = real(src(i), kind=4)
  end do
end subroutine


subroutine allclose(size, b1, b2, tol)
  implicit none
  integer(kind=8), intent(in) :: size
  real(kind=8), dimension(size), intent(in) :: b1, b2
  real(kind=8), intent(out) :: tol
  !local variables
  integer(kind=8) :: i
  tol=0.0_8

  do i=1,size
    tol=max(tol, abs(b1(i)-b2(i)))
  end do
  if (tol > 1.e-8) print *,b1(1:5),b2(1:5)
end subroutine

subroutine allclose_s(size, b1, b2, tol)
  implicit none
  integer(kind=8), intent(in) :: size
  real(kind=4), dimension(size), intent(in) :: b1, b2
  real(kind=4), intent(out) :: tol
  !local variables
  integer(kind=8) :: i
  tol=0.0_4

  do i=1,size
    tol=max(tol, abs(b1(i)-b2(i)))
  end do
  if (tol > 1.e-4) print *,b1(1:5),b2(1:5)
end subroutine


subroutine s0s0_3d_info(op, bc, n, nout, flops)
  implicit none
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  integer(kind=4), intent(out), dimension(0:2) :: nout
  integer(kind=4), intent(out) :: flops
  !local variables
  integer :: i, cost
  real(kind=8) :: a, ax, ay
  integer, dimension(0:2) :: dimsin, dimsout

  !those values should not be important for the dims
  a=1.d0
  ax=1.d0
  ay=0.d0

  dimsin = n
  dimsout=dimsin
  flops=0
   do i=0,2
     call d_s0s0_1d_dims(op,3,i,dimsin,bc(i), a, ax, ay, dimsout)
     nout(i)=dimsout(i)
     print*,i,'dims in ', dimsin
     print*,i,'dims out', dimsout
     call d_s0s0_1d_cost(op,3,i,dimsin,bc(i),&
          dimsin,&
          dimsout,1,&
          a, ax, ay, cost)
     print *, "cost for call ",i,":", cost, " flops"
   dimsin=dimsout
   flops=flops+cost
 end do

end subroutine s0s0_3d_info


subroutine s0s0_3d_in_place(op,bc,n,x,y)
  implicit none
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  real(kind=8), intent(inout), dimension(*) :: x, y

  call s0s0_3d(op,bc,n,n,x,y,y,x)
end subroutine

subroutine s0s0_3d(op,bc,n,x,y,tmp1,tmp2)
  implicit none
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  real(kind=8), intent(inout), dimension(*) :: x, y ,tmp1, tmp2
  !local variables
  integer, dimension(0:2) :: dimsin, dimsout
  real(kind=8) :: a, ax, ay

  a=1.d0
  ax=0.d0
  ay=0.d0

  dimsin = n
  dimsout=dimsin
  call d_s0s0_1d_dims(op,3,0,dimsin,bc(0), a, ax, ay, dimsout)
  print *,dimsin
  print *,dimsout,'out'
  call d_s0s0_1d(op,3,0,dimsin,bc(0),&
          dimsin,&
          dimsout,1,&
          x, tmp1, &
          a, ax, ay)  
  dimsin=dimsout
  call d_s0s0_1d_dims(op,3,1,dimsin,bc(1), a, ax, ay, dimsout)
  call d_s0s0_1d(op,3,1,dimsin,bc(1),&
          dimsin,&
          dimsout,1,&
          tmp1,tmp2,&
          a, ax, ay)
  dimsin=dimsout
  call d_s0s0_1d_dims(op,3,2,dimsin,bc(2), a, ax, ay, dimsout)  
  call d_s0s0_1d(op,3,2,dimsin,bc(2),&
          dimsin,&
          dimsout,1,&
          tmp2,y,&
          a, ax, ay)

end subroutine s0s0_3d

!!$            do i=0,2
!!$              call d_s0s0_1d_dims(SYM8_IMF,3,i,dimsin,bc(i), 1.0_wp, 1.0_wp,0.0_wp, dimsout)
!!$!              print*,i,'dims in ', dimsin
!!$!              print*,i,'dims out', dimsout
!!$
!!$            call d_s0s0_1d_cost(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   1.0_wp, 1.0_wp,0.0_wp, cost)
!!$!            print *, "cost for call ",i,":", cost, " flops"
!!$            if(i/=1) then
!!$              call d_s0s0_1d(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   psir(1,idx),w%y_c(1,idx),&
!!$                   1.0_wp, 1.0_wp,0.0_wp)
!!$            else
!!$              call d_s0s0_1d(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   w%y_c(1,idx),psir(1,idx),&
!!$                   1.0_wp, 1.0_wp,0.0_wp)
!!$            endif
!!$            dimsin=dimsout
!!$          end do

!!$              do i=0,2
!!$                call d_s0s1_1d_dims(SYM8_DWT,3,i,dimsin,bc(i), 1.0_wp,0.0_wp, dimsout)
!!$                print*,i,'dims in ', dimsin
!!$                print*,i,'dims out', dimsout
!!$
!!$                call d_s0s1_1d_cost(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   1.0_wp,0.0_wp, cost)
!!$                print *, "cost for call ",i,":", cost, " flops"
!!$                if(i==1) then
!!$                  call d_s0s1_1d(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   psir(1,idx),w%y_c(1,idx),&
!!$                   1.0_wp,0.0_wp)
!!$                else
!!$                  call d_s0s1_1d(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   w%y_c(1,idx),psir(1,idx),&
!!$                   1.0_wp,0.0_wp)
!!$                endif
!!$                dimsin=dimsout
!!$             end do

!!$            do i=0,2
!!$              call d_s0s0_1d_dims(SYM8_IMF,3,i,dimsin,bc(i), 1.0_wp, 1.0_wp,0.0_wp, dimsout)
!!$!              print*,i,'dims in ', dimsin
!!$!              print*,i,'dims out', dimsout
!!$
!!$            call d_s0s0_1d_cost(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   1.0_wp, 1.0_wp,0.0_wp, cost)
!!$            print *, "cost for call ",i,":", cost, " flops"
!!$            if(i/=1) then
!!$              call d_s0s0_1d(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   psir(1,idx),w%y_c(1,idx),&
!!$                   1.0_wp, 1.0_wp,0.0_wp)
!!$            else
!!$              call d_s0s0_1d(SYM8_IMF,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   w%y_c(1,idx),psir(1,idx),&
!!$                   1.0_wp, 1.0_wp,0.0_wp)
!!$            endif
!!$            dimsin=dimsout
!!$          end do

!!$              do i=0,2
!!$                call d_s0s1_1d_dims(SYM8_DWT,3,i,dimsin,bc(i), 1.0_wp,0.0_wp, dimsout)
!!$                print*,i,'dims in ', dimsin
!!$                print*,i,'dims out', dimsout
!!$
!!$                call d_s0s1_1d_cost(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   1.0_wp,0.0_wp, cost)
!!$                print *, "cost for call ",i,":", cost, " flops"
!!$                if(i==1) then
!!$                  call d_s0s1_1d(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   psir(1,idx),w%y_c(1,idx),&
!!$                   1.0_wp,0.0_wp)
!!$                else
!!$                  call d_s0s1_1d(SYM8_DWT,3,i,dimsin,bc(i),&
!!$                   dimsin,&
!!$                   dimsout,1,&
!!$                   w%y_c(1,idx),psir(1,idx),&
!!$                   1.0_wp,0.0_wp)
!!$                endif
!!$                dimsin=dimsout
!!$             end do


