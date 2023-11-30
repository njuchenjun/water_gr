subroutine gwrite_poscar(fid,natom,at,ae,atom,x,p)
  implicit none
  integer :: fid,natom
  real*8 :: at,ae,x(3,natom),p(3,natom)
  integer :: i
  character(kind=1,len=2) :: atom(natom)

  write(fid,'(a)') "water on c32"
  write(fid,'(a)') "   1.000000000"
  write(fid,'(a)') "   9.869918322    0.000000000    0.000000000"
  write(fid,'(a)') "  -4.934959161    8.547600000    0.000000000"
  write(fid,'(a)') "   0.000000000    0.000000000   20.000000000"
  write(fid,'(a)') " O H C"
  write(fid,'(a)') " 1 2 32"
  write(fid,'(a)') "Selective dynamics"
  write(fid,'(a)') "Cartesian"
  do i=1,natom
     write(fid,'(3f16.8,a)') x(1:3,i)," T T T"
  enddo

  return
end subroutine
