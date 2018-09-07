require 'rubygems'
require 'BOAST'

def kinetic_per_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "kinetic_per_ref"
  n1 = BOAST::Int("n1", :dir => :in)
  n2 = BOAST::Int("n2", :dir => :in)
  n3 = BOAST::Int("n3", :dir => :in)
  hgrid = BOAST::Real("hgrid", :dir => :in, :dim => [BOAST::Dim(3)] )
  c = BOAST::Real("c", :dir => :in)
  x = BOAST::Real("x", :dir => :in,  :dim => [ BOAST::Dim(0, n1 - 1),  BOAST::Dim(0, n2 - 1), BOAST::Dim(0, n3 - 1)] )
  y = BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim(0, n1 - 1),  BOAST::Dim(0, n2 - 1), BOAST::Dim(0, n3 - 1)] )
  
  p = BOAST::Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,c])
  kernel.code.print <<EOF
subroutine kinetic_per_ref(n1,n2,n3,hgrid,x,y,c)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(in) :: x
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
!  scale(:)=real(-.5_gp/hgrid(:)**2,wp)
  scale(:)=1.0_gp
  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  
  do i3=0,n3-1
     ! (1/2) d^2/dx^2
     do i2=0,n2-1
        do i1=0,n1-1
           tt=x(i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1-1
        do i2=0,n2-1
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2-1
     do i1=0,n1-1
        do i3=0,n3-1
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
  enddo
  
END SUBROUTINE kinetic_per_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


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
  kernel.code.print <<EOF
!> Simple non-optimized version of the major convolution routines
subroutine magicfilter_free_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-lupfil:n1-1-lowfil), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  do j=1,ndat
     do i=-lupfil,n1-1-lowfil
        tt=0.e0_wp
        do l=max(-i,lowfil),min(lupfil,n1-1-i)
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_free_ref


!> Simple non-optimized version of the major convolution routines
subroutine magicfilter_free_t_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, parameter :: lowfil=-7,lupfil=8
  integer, intent(in) :: n1,ndat
  real(wp), dimension(lowfil:n1-1+lupfil,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1-1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /

  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_free_t_ref

subroutine magicfilter_per_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n1,ndat
  real(kind=8), dimension(0:(n1-1),ndat), intent(in) :: x
  real(kind=8), dimension(ndat,0:(n1-1)), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
  real(kind=8) :: tt

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(kind=8) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /
  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)  
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_per_ref

subroutine magicfilter_per_t_ref(n1,ndat,x,y)
  !use module_base
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1-1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: i,j,k,l
  real(wp) :: tt
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /


  do j=1,ndat
     do i=0,n1-1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)
           tt=tt+x(k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_per_t_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end



