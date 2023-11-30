
program main
  implicit none
  real*8,parameter :: amu=1822.88848477d0
  real*8,parameter :: ang=0.5291772083d0
  real*8,parameter :: ev=27.21138505d0
  real*8,parameter :: kcal=23.0605d0
  real*8,parameter :: tenfs2au=0.0024189d0
  integer,parameter :: natom=35
  
  real*8 :: mass(natom)
  character(kind=1,len=2) :: atom(natom),atmp
  character(kind=1,len=6),parameter :: molecule="OH2   "
  character(kind=1,len=99) :: fname
  integer :: fid
  character(kind=1,len=4) :: fid_c

  real*8 :: x(3,natom),p(3,natom) ! mass: amu, v: angstrom/10fs, x: angstrom
  real*8 :: at,ae ! time: 10fs, total energy: kcal/mol
  real*8 :: erot_oh2,evib_oh2,ek_oh2,emod_oh2(9)
  real*8 :: AM_oh2,J_oh2
  integer i,j,id0

  real*8 :: xs(3,natom),e_carbon,dedx_carbon(3,32),e_oh2,dedx_oh2(3,3),x_int(3,natom),e_int,dedx_int(3,natom)

  integer :: ioh2(3)
  real*8 :: rtmp(4),vector(3),rmax

  mass(1)=15.99491461956d0
  mass(2:3)=1.00782503207d0
  mass(4:35)=12.d0
  atom(1)="O "
  atom(2:3)="H "
  atom(4:35)="C "
  ioh2=[1,2,3]

  fid=2003
  open(901,file='x3.txt',status='old',action='read')
  read(901,*)
  read(901,*)
  do i=1,35
     read(901,*) atmp,x(1:3,i)
  enddo
  close(901)

  do i=10,100,2
     xs=x
     xs(3,1:3)=xs(3,1:3)+(i-40)*0.05d0

   ! call shift_surf(x,xs)

     call gwrite(fid,natom,at,ae,atom,xs,p)
     call gwrite_poscar(fid+1000,natom,at,ae,atom,xs,p)
     
     call pes_graphene(xs(:,4:35),e_carbon,dedx_carbon)
     call pes_oh2(xs(:,1:3),e_oh2,dedx_oh2)
     x_int(:,1:32)=xs(:,4:35); x_int(:,33:35)=xs(:,1:3)
     call pes_int(x_int,e_int,dedx_int)
     
     write(*,'(5f16.8)') xs(3,1)-2.d0,e_carbon+e_oh2+e_int,e_carbon,e_oh2,e_int

  enddo

end program

