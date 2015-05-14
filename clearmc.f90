program clearmc !DOC  
use parameters
!DOC !FILE !PROG the program is used to clear the files created during the simulation, such as parameter.prm, trajectory and frames files.
!DOC !PROG Deletes all the files listed in parameters.prm namely: energy_file,boxlength_file,frames_file,traj_file
!DOC !PROG Is used before re-start the simulation.
!
!DOC !PROG  Usage:  clearmc parameters.prm
!
!DOC !PROG Arguments:
!DOC !PROG ::  parameters.prm   the file with the parameters of the simulation

character(SYSTEM_STRING_LENGTH) :: parameters_file

if(  iargc() < 1) then

   write(*,*) 'Usage: clearmc parameters.prm'
   write(*,*) '  Deletes all the files listed in parameters.prm namely: energy_file,boxlength_file,frames_file,traj_file'
   return
end if

   
call getarg(1,parameters_file)
call Parameters_init(parameters_file,0)  ! MAxMolecule is used for the BoxLength calculation (from density, if required)

call rm_file(energy_file)
call rm_file(boxlength_file)
call rm_file(traj_file)
call rm_file(frames_file)

contains

subroutine rm_file(fname) !DOC
!DOC remove the file with a given name
!DOC Parameters:
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   character(SYSTEM_STRING_LENGTH) :: fname !DOC name of the file
   integer :: stat

   open(unit=1234, iostat=stat, file=fname, status='old')
   if( stat .eq. 0 ) close(1234,status='delete')

end subroutine



end program clearmc


