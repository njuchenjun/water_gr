      subroutine pes_oh2_au(x,v,dvdx)
      implicit none
      integer,parameter :: natom=3
      real*8,intent(in) :: x(3,natom)  !! cartesian of O H H in a.u.
      real*8,intent(out) :: v,dvdx(3,natom) !! energy and gradient in a.u.
      real*8 :: q(3,natom),vout,dvdq(3,natom) !! a.u
      real*8,parameter :: angstrom=0.5291772083d0
      real*8,parameter :: ev=27.21138505d0
      q=x ! now q in a.u.
      call energy1_au(q,vout,dvdq) ! all in a.u.
      v=vout ! now v in a.u.
      dvdx=-dvdq ! now dvdx in a.u.
      return
      end subroutine

      subroutine pes_oh2(x,v,dvdx)
      implicit none
      integer,parameter :: natom=3
      real*8,intent(in) :: x(3,natom)  !! cartesian of O H H in Angstrom
      real*8,intent(out) :: v,dvdx(3,natom) !! energy and gradient in eV(/Angstrom)
      real*8 :: q(3,natom),vout,dvdq(3,natom) !! a.u
      real*8,parameter :: angstrom=0.5291772083d0
      real*8,parameter :: ev=27.21138505d0
      q=x/angstrom ! now q in a.u.
      call energy1_au(q,vout,dvdq) ! all in a.u.
      v=vout*ev ! now v in eV
      dvdx=-dvdq*ev/angstrom ! now dvdx in eV/angstrom
      return
      end subroutine

! Training set : 5344 points, exchange r1 and r2 to 10688 points
! /home/chenjun/h2o/7-pes-cc-cbs/ae-ccsd/nnfit/W07.txt 3-30-30-1
!  Total rmse_07=   0.0009 meV; train=   0.0006; validation=   0.0023

! provide bond in bohr and angle in degree, return energy in Hartree
! r1=r(O-H1), r2=r(O-H2), theta=angle(H1-O-H2)
! r=1.80856326933172d0;theta=104.611249045011d0;v0=-76.4393607300896
! --------------------------------------------------------------------
      subroutine h2o_nn(r1,r2,theta,v_out)
      implicit none
      integer,parameter :: kind_use=8
      real(kind=kind_use),intent(in) :: r1,r2,theta
      real(kind=kind_use),intent(out) :: v_out
      real(kind=kind_use) :: r3
      r3=dsqrt(r1**2+r2**2-2.D0*r1*r2*dcos(theta/180.d0*dacos(-1.d0)))
      call potdriv_h2o([r1,r2,r3],v_out)
      return
      end

! --------------------------------------------------------------------
! provide bond distances of OH1 OH2 HH in bohr, return energy in Hartree
! r: O-H1, O-H2, H1-H2
      subroutine potdriv_h2o(r_in,v_out)
      implicit none
      integer,parameter :: kind_use=8
      real(kind=kind_use),intent(in) :: r_in(3)
      real(kind=kind_use),intent(out) :: v_out
      real*8 :: r1(3),r2(3),v1,v2
      real*8,parameter :: vmin=-76.4393607300896d0-1.4210854715202D-014 ! opt

      r1([1,2,3])=r_in([1,2,3])
      r2([1,2,3])=r_in([2,1,3])
      call nninter(r1,v1)
      call nninter(r2,v2)
      v_out=(v1+v2)/2.D0
      v_out=v_out-vmin

      return
      end
! --------------------------------------------------------------------
      subroutine nsim(rin,vx,ndim,neu1,neu2,nxin,nxw1,nxb1,
     &                nxw2,nxb2,nxw3,nxb3,nxout)
!     for two-hidden-layer feedforward neural network, I-J-K-1
!     two blas routines used: dgemv, ddot
      implicit none
      integer,intent(in) :: ndim,neu1,neu2
      real*8,intent(in) :: rin(ndim),nxin(2,ndim),
     &                     nxw1(ndim,neu1),nxb1(neu1),
     &                     nxw2(neu1,neu2),nxb2(neu2),
     &                     nxw3(neu2),nxb3,nxout(2)
      real*8,intent(out) :: vx
      integer :: i,j,k
      real*8 :: r(ndim),ax(neu1),bx(neu2),vrange,rrange(ndim),
     &          rtmp,rt1(neu1),rt2(neu2)
      real*8,external :: ddot

      vx=0.d0
      r=rin
      vrange=nxout(2)-nxout(1)
      ! mapminmax [-1,1]
      do i=1,ndim
        rrange(i)=nxin(2,i)-nxin(1,i)
        r(i)=2.d0*(r(i)-nxin(1,i))/rrange(i)-1.d0
      end do

      ! 1st layer
      rt1=nxb1
      call dgemv('t',ndim,neu1,1.d0,nxw1,ndim,r,1,1.d0,rt1,1)
      do j=1,neu1
        ax(j)=dtanh(rt1(j))
      end do

      ! 2nd layer
      rt2=nxb2
      call dgemv('t',neu1,neu2,1.d0,nxw2,neu1,ax,1,1.d0,rt2,1)
      do k=1,neu2
        bx(k)=dtanh(rt2(k))
      end do

      ! output layer
      vx=nxb3+ddot(neu2,nxw3,1,bx,1)

      !reverse map
      vx=vrange*(vx+1.d0)/2.d0+nxout(1)

      return
      end

! --------------------------------------------------------------------
      module nnmod
      implicit none
      save
!     character*99 :: n01fa='/home/chenjun/h2o/7-pes-cc-cbs/ae-ccsd/nnfit/W07.txt'
      character*99 :: n01fa='pes_h2ocbs.txt'
      integer,parameter :: ndim=3
      integer,parameter :: n01s1=30,n01s2=30

c     01 double hidden layers
      real*8 :: n01w1a(ndim,n01s1),n01b1a(n01s1),
     &n01w2a(n01s1,n01s2),n01b2a(n01s2),n01w3a(n01s2),
     &n01b3a,n01ina(2,ndim),n01outa(2)
      contains
      subroutine nninit
      open(158301,file=trim(n01fa),status='old')
      read(158301,*) n01w1a,n01b1a,n01w2a,n01b2a,
     &               n01w3a,n01b3a,n01ina,n01outa
      close(158301)
      return
      end subroutine nninit
      end module nnmod
! --------------------------------------------------------------------
      subroutine nninter(r0,vx)
      use nnmod
      implicit none
      real*8,intent(in) :: r0(ndim)
      real*8,intent(out) :: vx
      real*8 :: r(ndim),tmp,va,vb,vc,vd,ve,vf
      integer,save :: init_h2opes=0
      integer i,j,k
      if (init_h2opes.eq.0) then
        init_h2opes=1
        call nninit
      endif

      r(1:ndim)=r0(1:ndim)
      call nsim(r,va,ndim,n01s1,n01s2,n01ina,n01w1a,n01b1a,
     &n01w2a,n01b2a,n01w3a,n01b3a,n01outa)

      vx=va
      return
      end
! --------------------------------------------------------------------
