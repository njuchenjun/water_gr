#define natom 3
#define nbond (natom*(natom-1)/2)
! ====================================================================
! original program from liulan
! q:    cartesian in bohr
! vout: energy in hartree
! dvdq: derivatives in hartree/bohr
! --------------------------------------------------------------------
subroutine energy1_au(q,vout,dvdq)
  implicit none
  real*8,intent(in) :: q(3,natom)
  real*8,intent(out) :: vout,dvdq(3,natom)

  integer i,j,k
  real*8 :: rr(nbond),v1,v2,delta=1.d-3
  real*8 :: r(nbond),dr(3,natom,nbond),dv(3,natom),dvr(nbond),dvrx(nbond)
  integer,save :: mb(nbond),nb(nbond)
  integer,save :: init=0
  real*8 :: rstore(nbond,2*nbond+1),vstore(2*nbond+1)

  if(init.eq.0) then
    k=0
    do i=1,natom
    do j=i+1,natom
      k=k+1
      mb(k)=i
      nb(k)=j
    enddo
    enddo
    init=1
  endif

  call intern(q,mb,nb,r,dr)
  rstore(1:nbond,1)=r

  do i=1,nbond
    rr=r
    rr(i)=rr(i)-delta
    rstore(1:nbond,i*2)=rr
    rr(i)=rr(i)+2.d0*delta
    rstore(1:nbond,i*2+1)=rr
  enddo

  do i=1,2*nbond+1
    call potdriv_h2o(rstore(1:nbond,i),vstore(i)) !! O-H1, O-H2, H1-H1 all in a.u
  enddo

  vout=vstore(1)
  do i=1,nbond
    v1=vstore(i*2)
    v2=vstore(i*2+1)
    dvr(i)=(v2-v1)/2.d0/delta
  enddo

  call dvdx(mb,nb,dr,dvr,dv)
  dvdq=dv

  return
end subroutine energy1_au
! --------------------------------------------------------------------
subroutine intern(q,mb,nb,r,dr)
! ---- calculates bondlengths and their derivatives
  implicit real*8 (a-h, o-z)
  real*8,intent(in) :: q(3,natom)
  integer,intent(in) :: mb(nbond),nb(nbond)
  real*8,intent(out) :: r(nbond),dr(3,natom,nbond)
  real*8 :: qx(3)
  do i=1,nbond
    qx(1:3)=q(1:3,mb(i))-q(1:3,nb(i))
! ---- bond lengths
    r(i)=dsqrt(dot_product(qx,qx))
! ---- bond length derivatives
    dr(1,mb(i),i)=qx(1)/r(i)
    dr(2,mb(i),i)=qx(2)/r(i)
    dr(3,mb(i),i)=qx(3)/r(i)

    dr(1,nb(i),i)=-dr(1,mb(i),i)
    dr(2,nb(i),i)=-dr(2,mb(i),i)
    dr(3,nb(i),i)=-dr(3,mb(i),i)
  enddo
  return
end subroutine intern
! --------------------------------------------------------------------
subroutine dvdx(mb,nb,dr,dvr,dv)
! ---- calculates dv/dx from dv/dr*dr/dx
  implicit none
  integer,intent(in) :: mb(nbond),nb(nbond)
  real*8,intent(in) :: dr(3,natom,nbond),dvr(nbond)
  real*8,intent(out) :: dv(3,natom)
  integer :: i
  dv=0.d0
  do i=1,nbond
     dv(1,mb(i))=dv(1,mb(i))+dvr(i)*dr(1,mb(i),i)
     dv(1,nb(i))=dv(1,nb(i))+dvr(i)*dr(1,nb(i),i)
     dv(2,mb(i))=dv(2,mb(i))+dvr(i)*dr(2,mb(i),i)
     dv(2,nb(i))=dv(2,nb(i))+dvr(i)*dr(2,nb(i),i)
     dv(3,mb(i))=dv(3,mb(i))+dvr(i)*dr(3,mb(i),i)
     dv(3,nb(i))=dv(3,nb(i))+dvr(i)*dr(3,nb(i),i)
  enddo
  return
end subroutine dvdx
! ====================================================================
