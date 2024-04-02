subroutine generic_op_3d_info(procedure_dims, procedure_cost, op, bc, n, nout, flops)
  implicit none
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  integer(kind=4), intent(out), dimension(0:2) :: nout
  integer(kind=4), intent(out) :: flops
  external :: procedure_dims, procedure_cost
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
     call procedure_dims(op,3,i,dimsin,bc(i), a, ax, ay, dimsout)
     nout(i)=dimsout(i)
     print*,i,'dims in ', dimsin
     print*,i,'dims out', dimsout
     call procedure_cost(op,3,i,dimsin,bc(i),&
          dimsin,&
          dimsout,1,&
          a, ax, ay, cost)
     print *, "cost for call ",i,":", cost, " flops"
   dimsin=dimsout
   flops=flops+cost
 end do

end subroutine generic_op_3d_info

subroutine generic_op_3d(procedure_dims, procedure_op,op,bc,n,x,y,tmp1,tmp2)
  implicit none
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  real(kind=8), intent(inout), dimension(*) :: x, y ,tmp1, tmp2
  external :: procedure_dims, procedure_op
  !local variables
  integer, dimension(0:2) :: dimsin, dimsout
  real(kind=8) :: a, ax, ay

  a=1.d0
  ax=0.d0
  ay=0.d0

  dimsin = n
  dimsout=dimsin
  call procedure_dims(op,3,0,dimsin,bc(0), a, ax, ay, dimsout)
  print *,dimsin
  print *,dimsout,'out'
  call procedure_op(op,3,0,dimsin,bc(0),&
          dimsin,&
          dimsout,1,&
          x, tmp1, &
          a, ax, ay)  
  dimsin=dimsout
  call procedure_dims(op,3,1,dimsin,bc(1), a, ax, ay, dimsout)
  call procedure_op(op,3,1,dimsin,bc(1),&
          dimsin,&
          dimsout,1,&
          tmp1,tmp2,&
          a, ax, ay)
  dimsin=dimsout
  call procedure_dims(op,3,2,dimsin,bc(2), a, ax, ay, dimsout)  
  call procedure_op(op,3,2,dimsin,bc(2),&
          dimsin,&
          dimsout,1,&
          tmp2,y,&
          a, ax, ay)

end subroutine generic_op_3d

subroutine op_3d(cp,op,bc,n, a, ax, ay,x,y,tmp1,tmp2)
  implicit none
  include 'libconvf.h'
  integer(kind=4), intent(in) :: cp
  integer(kind=4), intent(in) :: op
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  real(kind=8), intent(inout), dimension(*) :: x, y ,tmp1, tmp2
  real(kind=8) :: a, ax, ay !no intent to avoid copy in the stack
  external :: d_s0s0_1d_dims,d_s0s0_1d, s_s0s0_1d_dims,s_s0s0_1d
  external :: d_s1s0_1d_dims,d_s1s0_1d, s_s1s0_1d_dims,s_s1s0_1d
  external :: d_s0s1_1d_dims,d_s0s1_1d, s_s0s1_1d_dims,s_s0s1_1d

  select case(op)
  case(SYM8_MF,SYM8_IMF,SYM8_D1,SYM8_D2)
    call precision_switch(d_s0s0_1d_dims,d_s0s0_1d, s_s0s0_1d_dims,s_s0s0_1d)
  case(SYM8_IDWT,SYM8_S1TOR)
    call precision_switch(d_s1s0_1d_dims,d_s1s0_1d, s_s1s0_1d_dims,s_s1s0_1d)
  case(SYM8_DWT,SYM8_RTOS1)
    call precision_switch(d_s0s1_1d_dims,d_s0s1_1d, s_s0s1_1d_dims,s_s0s1_1d)
  end select

 contains

  subroutine precision_switch(procd1, procd2, procs1, procs2)
   implicit none
   external :: procs1, procs2, procd1, procd2
   select case(cp)
   case(4)
    call gen_op(procs1, procs2)
   case(8)
    call gen_op(procd1, procd2)
   end select

  end subroutine precision_switch


 subroutine gen_op(procedure_dims, procedure_op)
  implicit none
  external :: procedure_dims, procedure_op
  !local variables
  integer, dimension(0:2) :: dimsin, dimsout

  ! a=1.d0
  ! ax=0.d0
  ! ay=0.d0

  dimsin = n
  dimsout=dimsin
  call procedure_dims(op,3,0,dimsin,bc(0), a, ax, ay, dimsout)
  print *,dimsin
  print *,dimsout,'out'
  call procedure_op(op,3,0,dimsin,bc(0),&
          dimsin,&
          dimsout,1,&
          x, tmp1, &
          a, ax, ay)  
  dimsin=dimsout
  call procedure_dims(op,3,1,dimsin,bc(1), a, ax, ay, dimsout)
  call procedure_op(op,3,1,dimsin,bc(1),&
          dimsin,&
          dimsout,1,&
          tmp1,tmp2,&
          a, ax, ay)
  dimsin=dimsout
  call procedure_dims(op,3,2,dimsin,bc(2), a, ax, ay, dimsout)  
  call procedure_op(op,3,2,dimsin,bc(2),&
          dimsin,&
          dimsout,1,&
          tmp2,y,&
          a, ax, ay)

end subroutine gen_op

end subroutine op_3d



subroutine op_3d_info(cp, op, bc, n, nout, flops)
  implicit none
  include 'libconvf.h'
  integer(kind=4), intent(in) :: op, cp
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  integer(kind=4), intent(out), dimension(0:2) :: nout
  integer(kind=4), intent(out) :: flops
  external :: d_s0s0_1d_dims,d_s0s0_1d, s_s0s0_1d_dims,s_s0s0_1d
  external :: d_s1s0_1d_dims,d_s1s0_1d, s_s1s0_1d_dims,s_s1s0_1d
  external :: d_s0s1_1d_dims,d_s0s1_1d, s_s0s1_1d_dims,s_s0s1_1d
  external :: d_s0s0_1d_cost, s_s0s0_1d_cost
  external :: d_s1s0_1d_cost, s_s1s0_1d_cost
  external :: d_s0s1_1d_cost, s_s0s1_1d_cost


  select case(op)
  case(SYM8_MF,SYM8_IMF,SYM8_D1,SYM8_D2)
    call precision_switch(d_s0s0_1d_dims,d_s0s0_1d_cost, s_s0s0_1d_dims,s_s0s0_1d_cost)
  case(SYM8_IDWT,SYM8_S1TOR)
    call precision_switch(d_s1s0_1d_dims,d_s1s0_1d_cost, s_s1s0_1d_dims,s_s1s0_1d_cost)
  case(SYM8_DWT,SYM8_RTOS1)
    call precision_switch(d_s0s1_1d_dims,d_s0s1_1d_cost, s_s0s1_1d_dims,s_s0s1_1d_cost)    
  end select

 contains

  subroutine precision_switch(procd1, procd2, procs1, procs2)
   implicit none
   external :: procs1, procs2, procd1, procd2
   select case(cp)
   case(4)
    call gen_info(procs1, procs2) !, sa, sax, say)
   case(8)
    call gen_info(procd1, procd2) !, da, dax, day)
   end select

  end subroutine precision_switch


  subroutine gen_info(procedure_dims,procedure_cost)
  !local variables
  integer :: i
  real(kind=8) :: a, ax, ay
  integer, dimension(0:2) :: dimsin, dimsout
  integer, dimension(1) :: cost


  !those values should not be important for the dims
  a=1.d0
  ax=1.d0
  ay=0.d0

  dimsin = n
  dimsout=dimsin
  flops=0
   do i=0,2
     call procedure_dims(op,3,i,dimsin,bc(i), a, ax, ay, dimsout)
     nout(i)=dimsout(i)
     print*,i,'dims in ', dimsin
     print*,i,'dims out', dimsout
     call procedure_cost(op,3,i,dimsin,bc(i),&
          dimsin,&
          dimsout,1,&
          a, ax, ay, cost)
     print *, "cost for call ",i,":", cost, " flops"
   dimsin=dimsout
   flops=flops+cost(1)
 end do
 end subroutine gen_info

end subroutine op_3d_info

subroutine op_3d_in_place(cp,op,bc,n,x,y)
  implicit none
  integer(kind=4), intent(in) :: op, cp
  integer(kind=4), dimension(0:2), intent(in) :: bc
  integer(kind=4), intent(in), dimension(0:2) :: n
  real(kind=8), intent(inout), dimension(*) :: x, y

  call op_3d(cp,op,bc,n,n,x,y,y,x)
end subroutine

