!! a qct normal mode analysis method
!! based on unpublished manuscript by Bowman's group
! --------------------------------------------------------------------
subroutine calc_rovib(n,mass,x0,p0,AM,AJ,erot,evib,emod,molecule)
  implicit none
  integer,intent(in) :: n       ! number of atoms
  real*8,intent(in)  :: mass(n)    ! mass for each atom
  real*8,intent(in)  :: x0(3,n) ! cartetion coordinate
  real*8,intent(in)  :: p0(3,n) ! moment
  character(kind=1,len=6),intent(in) :: molecule

  real*8 :: JA(3),AM,AJ   ! AM-angular moment, AJ-quantum number
  real*8 :: IA(3,3) ! moment-of-inertia tensor
  real*8 :: II(3)   ! diagonal elements of diagonalized moment-of-inertia tensor
  real*8 :: LL(3,3) ! orthogonal transform matrix
  real*8 :: WA(3)   ! angular velocity
  real*8 :: WW(3)   ! angular velocity after orthogonal transformation
  real*8 :: F(3*n,3*n) ! force constant matrix
  real*8 :: C(3*n,3*n) ! orthogonal matrix, FC=LmC
  real*8 :: Lm(3*n)
  real*8 :: OM(3*n)

  real*8 :: et,erot,evib    ! translational, rotational and vibrational energy
  real*8 :: epot    ! potential energy 

  real*8 :: x(3,n),v(3,n),p(3,n),x_com(3),v_com(3)
  real*8 :: xe(3,n),u(3*n),b(3*n),u0(3*n)
  real*8 :: q(3*n) ! displacement
  real*8 :: QQ(3*n),PP(3*n) ! mass-scaled coordinates and velocities in normal-mode
  real*8 :: emod(3*n)

  real*8 :: tmp1,tmp3(3)
  integer :: i,j,k,ipiv(3),info
  real*8 :: work(9*n-1)

! print*,"--------------------------------------------"
  if("#"//trim(molecule).eq."#") then
    print*,"molecule not defined in calc_rovib"
    stop
  endif

!---> move the origin to center of mass (com),
!---> and eliminate the translational velocity of com
  do j=1,3
    x_com(j)=dot_product(x0(j,1:n),mass(1:n))/sum(mass(1:n))
    v_com(j)=sum(p0(j,1:n))/sum(mass(1:n))
  enddo

  do i=1,n
    x(1:3,i)=x0(1:3,i)-x_com(1:3)
    v(1:3,i)=p0(1:3,i)/mass(i)-v_com(1:3)
    p(1:3,i)=v(1:3,i)*mass(i)
  enddo

!---> to eliminate the rotation, rotate all the velocities
!---> calculate angular moment
  JA=0.d0
  do i=1,n
    call dcross(x(1:3,i),p(1:3,i),tmp3)
    JA=JA+tmp3
  enddo

  !! JA is the angular moment, |JA| = dsqrt(dot_product(JA,JA))
  AM=dsqrt(dot_product(JA,JA)) ! angular moment
  AJ=0.5d0*(-1.d0+dsqrt(1.d0+4.d0*AM*AM)) ! rotational quantum number

!---> calculate moment-of-inertia tensor
  call inertia_tensor(n,mass,x,IA)

  !! II = LL+ * IA * LL, LL is orthogonal transform matrix
  LL=IA; call dsyev('V','U',3,LL,3,II,work,8,info)

  !! WA : angular velocity in center-of-mass frame, before orthgonal transformation
  !! WA = IA^-1 * JA, solve WA in linear equations IA * WA = JA
  WA=JA; call dsysv('U',3,1,IA,3,ipiv,WA,3,work,8,info) !! caution: IA & WA changed
  if(info.ne.0) stop 'ERROR: angular velocity solve failed in dsysv'

  !! WW : angular velocity in center-of-mass frame, after orthgonal transformation
  !! WW = LL+ * WA
  call dgemv('C',3,3,1.d0,LL,3,WA,1,0.d0,WW,1)

  !! calculate total translational energy, rotational energy and vibrational energy
  erot=0.5d0*dot_product(II,WW**2)
  et=0.d0
  do i=1,n
    et=et+0.5d0*mass(i)*dot_product(v(1:3,i),v(1:3,i))
  enddo

  call geten(n,x,epot,molecule)

  evib=et-erot+epot

! print*,et*627.51d0,"kcal/mol, total translational energy"
! print*,epot*627.51d0,"kcal/mol, potential energy"
! print*,(et+epot)*627.51d0,"kcal/mol, total energy"
! print*,erot*627.51d0,"kcal/mol, rotational energy"
! print*,evib*627.51d0,"kcal/mol, vibrational energy"

!---> rotate velocities to the molecule with no net rotation
  do i=1,n
    call dcross(WA,x(1:3,i),tmp3)
    v(1:3,i)=v(1:3,i)-tmp3
  enddo
! print*,"--------------------------------------------"

!---> define the mass-scaled coordinates and velocities
  xe=x; call optimize(3*n,xe,molecule) ! optimize to the nearest local minimum configuration

  do i=1,n
  do j=1,3
    u(i*3-3+j)=sqrt(mass(i))*x(j,i) ! mass-scaled coordinates
    b(i*3-3+j)=sqrt(mass(i))*v(j,i) ! mass-scaled velocities
    u0(i*3-3+j)=sqrt(mass(i))*xe(j,i) ! mass-scaled coordinates at the local minimum
  enddo
  enddo
  q=u-u0 ! mass-scaled displacement coordinates

!---> calculate the force constant matrix
  call force_constant(n,mass,xe,u0,F,molecule)

  C=F; call dsyev('V','U',3*n,C,3*n,Lm,work,9*n-1,info) ! C^t * F * C = Lm * I

  ! print*,"harmonic frequencies:"
  do i=1,3*n
    if(Lm(i).lt.0.d0) then
      OM(i)=-dsqrt(-Lm(i))
  !   write(*,'(f12.3,a6)') -OM(i)*219474.542d0,"i cm-1"
    else
      OM(i)=dsqrt(Lm(i))
  !   write(*,'(f12.3,a6)') OM(i)*219474.542d0,"  cm-1"
    endif
  enddo
! write(*,'(<3*n>f10.5)') C

  call dgemv('C',3*n,3*n,1.d0,C,3*n,q,1,0.d0,QQ,1)
  call dgemv('C',3*n,3*n,1.d0,C,3*n,b,1,0.d0,PP,1)

  do i=1,3*n
    emod(i)=0.5d0*PP(i)*PP(i)+0.5d0*OM(i)*OM(i)*QQ(i)*QQ(i)
  enddo

! i=7;print*,emod(i)*627.51d0,"kcal/mol -- umbrella mod energy"

  return

contains

  subroutine dcross(a,b,c)
  implicit none
  real*8,intent(in) :: a(3),b(3)
  real*8,intent(out):: c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  return
  end subroutine dcross

  subroutine inertia_tensor(n,m,x,IA)
  implicit none
  integer,intent(in) :: n ! number of atoms
  real*8,intent(in)  :: m(n) ! mass
  real*8,intent(in)  :: x(3,n) ! coordinate in center-of-mass frame
  real*8,intent(out) :: IA(3,3) ! moment-of-inertia tensor
  integer :: i
  IA=0.d0
  do i=1,n
  IA(1,1)=IA(1,1)+m(i)*(x(2,i)**2+x(3,i)**2)
  IA(2,2)=IA(2,2)+m(i)*(x(3,i)**2+x(1,i)**2)
  IA(3,3)=IA(3,3)+m(i)*(x(1,i)**2+x(2,i)**2)
  IA(1,2)=IA(1,2)-m(i)*x(1,i)*x(2,i)
  IA(1,3)=IA(1,3)-m(i)*x(1,i)*x(3,i)
  IA(2,3)=IA(2,3)-m(i)*x(2,i)*x(3,i)
  enddo
  IA(2,1)=IA(1,2)
  IA(3,1)=IA(1,3)
  IA(3,2)=IA(2,3)
  return
  end subroutine inertia_tensor

  subroutine optimize(n,x0,molecule)
  implicit none
  integer,intent(in) :: n ! n is number of atoms * 3
  real*8,intent(inout) :: x0(n)
  character*6,intent(in) :: molecule
  real*8 :: x(n),v0,v1,v2,vnew,dv(n),dx(n)
  integer :: i,icyc
  real*8 :: speed=2.d0
  real*8,parameter :: dr=4.d-5
  x=x0
  icyc=0
  101 continue
  icyc=icyc+1
  x0=x

  call geten(n/3,x,v0,molecule)
  do i=1,n
    x=x0
    x(i)=x(i)-dr
    call geten(n/3,x,v1,molecule)
    x(i)=x(i)+2.d0*dr
    call geten(n/3,x,v2,molecule)
    dv(i)=(v2-v1)/2.d0/dr
  enddo
  dx=-dv
  102 continue
  x=x+dx*speed
  call geten(n/3,x,vnew,molecule)
  if(v0.lt.vnew) then
    x=x-dx*speed
  ! print*,icyc," - opt speed halfed"
    speed=speed*0.5d0
    goto 102
  endif
  if((v0-vnew)/speed.ge.1.d-8) goto 101
  !print*,icyc,vnew*27.21138505d0*1000d0," meV, -- optimization finished"
  if(vnew*27.21138505d0*1000d0.gt.5.d0) stop "ERROR: wrong local minima"
  x0=x
  return
  end subroutine optimize

  subroutine force_constant(n,m,xe,u0,F,molecule)
  implicit none
  integer,intent(in) :: n ! n is number of atoms
  real*8,intent(in)  :: m(n) ! mass
  real*8,intent(in)  :: xe(3*n) ! cartesian coordinate of local minima
  real*8,intent(in)  :: u0(3*n) ! mass-scaled cartesian coordinate of local minima
  real*8,intent(out) :: F(3*n,3*n) ! force constant matrix
  character*6,intent(in) :: molecule
  real*8 :: D(3*n) ! first order derivatives
  integer :: i,j,k1,k2
  real*8 :: x(3*n),v0,v1,v2,v(3*n,3*n,4)
  real*8,parameter :: dr=1.d-4

  !call geten(n,xe,v0,molecule)
  !do i=1,3*n
  !x=xe
  !x(i)=x(i)-dr
  !call geten(n,x,v1,molecule)
  !x(i)=x(i)+2.d0*dr
  !call geten(n,x,v2,molecule)
  !D(i)=(v2-v1)/2.d0/dr
  !print*,[v1,v0,v2]*27.21138505d0*1000.d0," meV"
  !enddo
  !print*,"first derivatives at local minimum, atomic units"
  !print*,D
  v=0.d0
  do i=1,3*n
  do j=i,3*n
    x=xe; x(i)=x(i)+dr; x(j)=x(j)+dr; call geten(n,x,v(i,j,1),molecule)
    x=xe; x(i)=x(i)-dr; x(j)=x(j)-dr; call geten(n,x,v(i,j,2),molecule)
    x=xe; x(i)=x(i)+dr; x(j)=x(j)-dr; call geten(n,x,v(i,j,3),molecule)
    x=xe; x(i)=x(i)-dr; x(j)=x(j)+dr; call geten(n,x,v(i,j,4),molecule)
  enddo
  enddo

  do i=1,3*n
  do j=i,3*n
    k1=(i-1)/3+1
    k2=(j-1)/3+1
    F(i,j)=(v(i,j,1)+v(i,j,2)-v(i,j,3)-v(i,j,4))/4.d0/dr/dr
    F(i,j)=F(i,j)/sqrt(m(k1))/sqrt(m(k2))
    F(j,i)=F(i,j)
  ! if(i.eq.j) print*,F(i,j)
  enddo
  enddo
  !print*,"force constants at local minimum, unit of VENUS96"
  !print*,F
  return
  end subroutine force_constant

end subroutine calc_rovib
! --------------------------------------------------------------------

  subroutine geten(n,x,v,molecule) ! x,v in a.u.
  implicit none
  integer,intent(in) :: n ! number of atoms
  real*8,intent(in)  :: x(3,n)
  real*8,intent(out) :: v
  character(kind=1,len=6),intent(in) :: molecule

  real*8 :: dvdx(3,n)

  if(n.eq.3 .and. trim(adjustl(molecule)).eq."OH2") then

    call pes_oh2_au(x,v,dvdx)

 !elseif(n.eq.4 .and. trim(adjustl(molecule)).eq."CH3") then
 !  call nn_potch3_cart(x,v)
 !elseif(n.eq.3 .and. trim(adjustl(molecule)).eq."CO2") then
   !call pot_co2_cart(x,v)
 !elseif(n.eq.5 .and. trim(adjustl(molecule)).eq."CH4") then
   !call surf_ch4(v,x)
 !elseif(n.eq.2 .and. trim(adjustl(molecule)).eq."HCl") then
 !  call pot1d_hcl( dsqrt(dot_product(x(1:3,1)-x(1:3,2),x(1:3,1)-x(1:3,2))) , v )
 !elseif(n.eq.2 .and. trim(adjustl(molecule)).eq."HF") then
 !  call nn_pothf( dsqrt(dot_product(x(1:3,1)-x(1:3,2),x(1:3,1)-x(1:3,2))) , v )
  else
    print*, "pes not found"
    print*,n,trim(adjustl(molecule))
    stop
  endif

  return
  end subroutine geten


 ! include "/home/chenjun/CO2/pes-nn-f12a-tz/d-pes/co2_nn.f90" ! co2 pes
 ! include "/home/chenjun/ch3/nn_potch3_cart.f90" ! ch3 pes
 ! link /home/chenjun/ch4/pes-pipnn/libch4pipnn.a
