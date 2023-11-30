subroutine pes_graphene(x,e,dedx)
  use nnmod_pes1, only: init_pes1,evaluate_pes1,deallocate_pes1
  implicit none
  integer,parameter :: natom=32
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: e,dedx(3,natom)

  real*8 :: c(3,natom),v,dvdc(3,natom)
  character(kind=1,len=2),parameter :: elems(natom)=[&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C "]

  integer,save :: ifirst_graphene=0
  integer :: ifrac,iforce
  ifrac=0
  iforce=1
  if(ifirst_graphene.eq.0) then
    call init_pes1
    ifirst_graphene=1
  endif

  c=x
  call evaluate_pes1(ifrac,iforce,natom,elems,c,v,dvdc) ! PES of graphene
  e=v
  dedx=dvdc(1:3,1:natom)
  return
end subroutine

subroutine pes_int(x,e,dedx)
  use nnmod_pes2, only: init_pes2,evaluate_pes2,deallocate_pes2
  implicit none
  integer,parameter :: natom=35
  real*8,intent(in) :: x(3,natom)
  real*8,intent(out) :: e,dedx(3,natom)

  real*8 :: c(3,natom),v,dvdc(3,natom)
  character(kind=1,len=2),parameter :: elems(natom)=[&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "C ","C ","C ","C ","C ","C ","C ","C ",&
             "O ","H ","H "]

  integer,save :: ifirst_int=0
  integer :: ifrac,iforce
  ifrac=0
  iforce=1
  if(ifirst_int.eq.0) then
    call init_pes2
    ifirst_int=1
  endif

  c=x
  call evaluate_pes2(ifrac,iforce,natom,elems,c,v,dvdc) ! PES of graphene
  e=v
  dedx=dvdc(1:3,1:natom)
  return
end subroutine
