program forces2text !DOC !PROG
!DOC !FILE !PROG Convert the binary frame with forces trajectory to the text format
use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use io, only : io_open, io_close
use FloatingPoint

!DOC !PROG Usage:  forces2text  forces.ftab  > forces.txt'
!DOC !PROG Arguments:
!DOC !PROG  :: forces.ftab  one frame of force trajectory ( NAtom * 12 bytes, Natom*3 floats)   '
!DOC !PROG :: forces.txt  file which contains the text representation of forces components (3 columns)

implicit none

character(SYSTEM_STRING_LENGTH) :: binary_file,  output_file
character(SYSTEM_STRING_LENGTH) :: tmp
integer :: hin,hout

integer(4), dimension(13) :: stat_data
integer(4) :: file_size 
integer :: i
integer :: natom

!real(8), parameter :: kT_kcal_mol = 1 ! does not matter, because we do not save epsilon to xyz (ergo only coordinates matter )
real(8) :: fx,fy,fz


if ( iargc() < 1 ) then

   write(*,*)  'Usage:  forces2text  forces.ftab '
   write(*,*)  '   forces.ftab - one frame of force trajectory ( NAtom * 12 bytes, Natom*3 floats)   '

   stop
end if

call getarg(1,binary_file)
!call getarg(2,output_file)

hin = io_open(binary_file,'b')

call fstat(hin, stat_data)
file_size = stat_data(8)
natom = file_size/12;

!hout = io_open(output_file,'w')

do i=1,natom

   fx =  read_float(hin)
   fy =  read_float(hin)
   fz = read_float(hin)

   write(*,*)  fx,fy,fz

end do



call io_close(hin)
!call io_close(hout)



end program forces2text
