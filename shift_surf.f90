
subroutine shift_surf(x0,xs)
  implicit none
  integer,parameter :: natom=35
  real*8 :: x0(3,natom),xs(3,natom)
  real*8 :: r0,xref04(3),xref10(3),xref35(3),xref29(3),x_old(3),x_ref(3),x_shift(3),water_shiftx,water_shifty
  integer :: i,j

  ! C-C bond length: r0=1.4246d0 Angstrom
  ! unit cell:
  !     r0*sqrt(3)      0               0
  !    -r0*sqrt(3)/2    r0*3/2          0
  ! the 4x4 cell:
  !     r0*sqrt(3)*4    0               0
  !    -r0*sqrt(3)/2*4  r0*3/2*4        0
  !     0               0               20
  !
  ! the SouthWest C (iatom= 4):  0              0       2.0
  ! the SouthEast C (iatom=10):  r0*sqrt(3)*3.5 r0*0.5  2.0
  ! the NorthEast C (iatom=35):  r0*sqrt(3)*2   r0*5    2.0
  ! the NorthWest C (iatom=29): -r0*sqrt(3)*1.5 r0*4.5  2.0
  !
  r0=1.4246d0
  xref04=[ 0.d0, 0.d0, 2.d0]
  xref10=[ r0*sqrt(3.d0)*3.5d0, r0*0.5d0, 2.d0]
  xref35=[ r0*sqrt(3.d0)*2.0d0, r0*5.0d0, 2.d0]
  xref29=[-r0*sqrt(3.d0)*1.5d0, r0*4.5d0, 2.d0]

  xs=x0

  do j=1,3
    x_old(j)=(x0(j,4)+x0(j,10)+x0(j,35)+x0(j,29))/4.d0
    x_ref(j)=(xref04(j)+xref10(j)+xref35(j)+xref29(j))/4.d0
    x_shift(j)=x_ref(j)-x_old(j)
    xs(j,:)=xs(j,:)+x_shift(j)
  enddo

  water_shifty= floor(xs(2,1)/(r0*6.d0)) * r0*6.d0 * -1.d0
  water_shiftx= floor(xs(2,1)/(r0*6.d0)) * -1.d0*r0*sqrt(3.d0)*2.d0 * -1.d0
  xs(2,1:3)=xs(2,1:3)+water_shifty
  xs(1,1:3)=xs(1,1:3)+water_shiftx

  water_shiftx= floor((xs(1,1)+xs(2,1)/sqrt(3.d0))/(r0*sqrt(3.d0)*4.d0)) * (r0*sqrt(3.d0)*4.d0) * -1.d0
  xs(1,1:3)=xs(1,1:3)+water_shiftx

  return
end subroutine

