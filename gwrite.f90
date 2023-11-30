subroutine gwrite(fid,natom,at,ae,atom,x,p)
  implicit none
  integer :: fid,natom
  real*8 :: at,ae,x(3,natom),p(3,natom)
  integer :: i
  character(kind=1,len=2) :: atom(natom)
  write(fid,'(i5)') natom
  write(fid,'(f16.8," t/fs",f12.3)') ae,at
  do i=1,natom
     write(fid,'(a2,6f16.8)') atom(i),x(1:3,i),p(1:3,i)
  enddo
  return
end subroutine
