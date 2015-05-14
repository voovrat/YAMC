Module AtomicData !DOC  
!DOC !FILE The module contains the data structures and the functions to manipulate the coordinated of the atoms
!
!


   implicit none 

   Type :: TAtomicData !DOC
!DOC Fields:
       real(8), dimension(:), allocatable :: xx,yy,zz   !DOC !coordinates of atoms  in internal coordinates: BoxLength == 1
       real(8), dimension(:), allocatable :: sigma,epsilon,charge !DOC !LJ parameters of the atoms in internal coordinates: sigma/BoxLength, epsilon/kT, charge/e
       
       real(8), dimension(:), allocatable :: hard_core_angstr !DOC !in angstroems, because internal units are very inconvenient

       character(4),dimension(:),allocatable :: atomnames  !DOC !labels of the atoms 
       integer,dimension(:),allocatable :: molnum_by_atomnum !DOC !the index which can be used to detemine the molecule number of the atom with a given number
       real(8) :: BoxLength !DOC !in Angstroems
       integer :: natom  !DOC !number of atoms
       integer :: nalloc !DOC !allocated size of the arrays 
       logical :: hasAtomnames !DOC !indicates that the labels of the atoms was read from the input files (can be used for export into xyz format)
   End Type TAtomicData 
   
Contains

subroutine AtomicData_alloc(this,nalloc,BoxLength,allocAtomNames) !DOC
!DOC Allocate the arrays in the AtomicData structure
!DOC Parameters:

   Type(TAtomicData) :: this  !DOC !the AtomicData structure  (contains the arrays to be allocated)
   integer,intent(in)  :: nalloc  !DOC !size of the arrays 
   real(8),intent(in) :: BoxLength  !DOC Length of the box (in Angstroems)
   logical,intent(in) :: allocAtomNames !DOC whether or not the atom label arrays should be allocated 

   allocate( this % xx(nalloc) )
   allocate( this % yy(nalloc) )
   allocate( this % zz(nalloc) )
   allocate( this % sigma(nalloc) )
   allocate( this % epsilon(nalloc) )
   allocate( this % charge(nalloc) )
   allocate( this % molnum_by_atomnum(nalloc) )
  
   allocate( this % hard_core_angstr(nalloc) )

   if( allocAtomNames ) then
      allocate( this % atomnames(nalloc) )
      this % hasAtomnames = .TRUE.
   else
      this % hasAtomnames = .FALSE.
   end if

   this % natom = 0
   this % nalloc = nalloc
   this % BoxLength = BoxLength
end subroutine


subroutine AtomicData_dealloc(this) !DOC 
!DOC  Deallocate the AtomicData structure
!DOC Parameters:
   Type(TAtomicData) :: this !DOC !AtomicData structure 

   deallocate( this % xx )
   deallocate( this % yy )
   deallocate( this % zz )
   deallocate( this % sigma )
   deallocate( this % epsilon )
   deallocate( this % charge )
   deallocate( this % molnum_by_atomnum )

   deallocate( this % hard_core_angstr )

   if( this % hasAtomnames ) then
      deallocate( this % atomnames )
   end if 

   this % natom = 0
   this % nalloc = 0
   this % hasAtomnames = .FALSE.
end subroutine
 
subroutine AtomicData_save_to_xyz(this,filename) !DOC
!DOC writes the AtomicData structure to the file in xyz format
!DOC Parameters:
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use io, only : io_open,io_close
   use error
   implicit none
   Type(TAtomicData) :: this  !DOC ! AtomicData structure 
   character(SYSTEM_STRING_LENGTH), intent(in) :: filename  !DOC !the name of the file to save the data 
!   integer,intent(in) :: n  !number of molecules
!   real(8),dimension(:),intent(in) :: xx,yy,zz ! coordinates
!   character(4),dimension(:),intent(in) :: atomnames

! variables
   integer :: i
   integer :: hfile
   real(8) :: L ! size of the box

   if(.not.this % hasAtomnames) then
      write(error_message,*) 'AtomicData_save_to_xyz: atomnames not present!'
      call error_throw(ERROR_PARAMETER)
      return
   end if
 
   L = this % BoxLength

   hfile = io_open(filename,'w')
   if(error_code /= 0 ) return
!   open(33,file = filename, form='formatted' )

   ! first line - numer of atoms
   write(hfile,*) this % natom
   ! second line - name of the molecule
   write(hfile,*) 'Box of molecules'

   do i=1,this % natom

    
      write(hfile,*) this % atomnames(i) , this % xx(i)*L,   this % yy(i)*L,   this % zz(i)*L

   end do

   close(hfile)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module AtomicData
