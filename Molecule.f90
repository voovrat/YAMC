Module Molecule !DOC 
!DOC !FILE The module contains the datastructure to store the structure of the molecule
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   implicit none


   Type :: TMolecule !DOC
 !DOC Fields:
        integer :: Natoms = 0 !DOC number of atoms
        
        real(8), dimension(:), allocatable :: x,y,z  !DOC coordinates of atoms in angstroems
        real(8), dimension(:), allocatable :: sigma  !DOC sigmas in ansgtroems
        real(8), dimension(:), allocatable :: epsilon !DOC epsilon kcal/mol
        real(8), dimension(:), allocatable :: charge !DOC charges of atoms
        real(8), dimension(:), allocatable :: mass  !DOC mass of atoms
 
        real(8), dimension(:), allocatable :: hard_core !DOC hard_core diameter of atoms in angstroems
 
        character(4), dimension(:), allocatable :: atomnames !DOC names of the atoms

        character(SYSTEM_STRING_LENGTH) :: filename !DOC name of the file 
        logical :: hasNames !DOC whether or not the atomnames were initialized
   END TYPE TMolecule
   
 
CONTAINS

subroutine Molecule_init(this ) !DOC
!DOC set the molecule to "zero" state
!DOC Parameters:
   Type(TMolecule) :: this !DOC Molecule

   this % Natoms = 0 
   this % hasNames = .FALSE.

end subroutine

subroutine Molecule_allocate(this, Natoms) !DOC
!DOC allocate the arrays in Molecule structure
!DOC Parameters:
   Type(TMolecule) :: this !DOC Molecule 
   integer,intent(in) :: Natoms    !DOC number of atoms

   call Molecule_deallocate(this)

   allocate(this % x(Natoms) )
   allocate(this % y(Natoms) )
   allocate(this % z(Natoms) )
   allocate(this % sigma(Natoms) )
   allocate(this % epsilon(Natoms) )
   allocate(this % charge(Natoms) )

   allocate(this % atomnames(Natoms) )
   allocate(this % mass(Natoms) )

   allocate(this % hard_core(Natoms) )

   this % Natoms = Natoms

end subroutine


subroutine Molecule_deallocate(this) !DOC 
!DOC Deallocate the arrays in Molecule structure
!DOC Parameters:
   Type(TMolecule) :: this !DOC Molecule

   if( this % Natoms.gt.0) then

      deallocate(this % x)
      deallocate(this % y)
      deallocate(this % z)
      deallocate(this % sigma)
      deallocate(this % epsilon )
      deallocate(this % charge ) 
 
      deallocate(this % atomnames )
      deallocate(this % mass )

      deallocate( this % hard_core )

      this % Natoms = 0
   end if
end subroutine 

subroutine Molecule_copy(dest,src) !DOC
!DOC Make a copy dest = src. Deallocates and allocates the arrays, if necessary
!DOC Parameters:
    Type(TMolecule) :: dest,src !DOC destination and source molecules
    
     if(dest % Natoms /= src % Natoms ) then
        call Molecule_deallocate(dest)
        call Molecule_allocate( dest, src % Natoms )
     end if 
     
     dest % x(:) = src % x(:)
     dest % y(:) = src % y(:)
     dest % z(:) = src % z(:)
     dest % sigma(:) = src % sigma(:)
     dest % epsilon(:) = src % epsilon(:)
     dest % charge(:) = src % charge(:)
     dest % mass(:) = src % mass(:)

     dest % hard_core(:) = src % hard_core(:)

end subroutine

!subroutine Molecule_init_mass(this)
!   use error
!   use periodic_table, only : element_by_name,ptable,table_is_initialized,init_periodic_table
!   Type(TMolecule) :: this
!
!   integer :: i,el
!   
!   if ( .not. this % hasNames ) then
!      write(error_message,*) 'Molecule_init_mass: no atom names'
!      call error_throw(ERROR_INITIALIZATION) 
!      return
!   end if
!
!   if ( table_is_initialized == 0 ) then
!      call init_periodic_table
!   end if
!
!   do i=1,this % Natoms 
!
!      el = element_by_name( this % atomnames(i)(1:2) )
!     
!      if ( el .le. 0) then
!         write(error_message,*) 'Molecule_init_mass: unknown element', this % atomnames(i)(1:2) 
!         call error_throw(ERROR_PARAMETER)
!         return
!      end if
!
!      this % mass(i) = ptable(el) % mass
!
!   end do
!
!end subroutine


subroutine Molecule_read_from_file(this, filename) !DOC 
!DOC read the molecule from mol file
   use io, only : io_count_file_lines,io_open,io_close
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use string

   implicit none
!DOC Parameters:
   Type(TMolecule) :: this !DOC Molecule
   character(SYSTEM_STRING_LENGTH),intent(in) :: filename !DOC name of the mol file
   !  File format : rism, i.e.  x y z sigma epsilon charge
 
   integer :: nlines
   integer :: hfile
   integer :: i
   integer :: L

   L = str_get_next_pos(filename,1,' ');


   this % filename(1:L) = filename(1:L)


   nlines = io_count_file_lines(filename)

 !  write(*,*) 'Molecule_read_from_file',filename(1:20),' nlines=',nlines

   call Molecule_allocate(this, nlines)

!   write(*,*) 'this % Natoms=',this % NAtoms

   hfile =  io_open(filename,'r')

   do i=1,nlines
      read(hfile,*) this % atomnames(i), this % x(i), this % y(i), this % z(i), this % sigma(i), & 
                    this % epsilon(i), this % charge(i) ,this % mass(i),this % hard_core(i)

      if( this % hard_core(i) < 1e-9 ) this % hard_core(i) = this % sigma(i) * sqrt(0.4) 

   end do  

   call io_close(hfile)

   this % NAtoms = nlines
   this % hasNames = .TRUE. 

end subroutine

subroutine Molecule_centrate(this) !DOC
!DOC centrate the molecule to the center of mass (i.e. make the coordinates of the center of mass (0,0,0)
   use geometry, only : center_of_mass
   implicit none
!DOC Parameters: 
  Type(TMolecule) :: this !DOC Molecule
 
   real(8) :: x_center,y_center,z_center 
   integer :: i

   call center_of_mass(this % x,this % y,this % z,this % mass,this % NAtoms,&
                       x_center, y_center, z_center)
   
   do i=1,this % NAtoms
      this % x(i) = this % x(i) - x_center
      this % y(i) = this % y(i) - y_center
      this % z(i) = this % z(i) - z_center
   end do

end subroutine


!   OBSOLETE
!
!subroutine Molecule_read_atomnames(this, filename)
!   use error
!   use io, only : io_count_file_lines, io_open, io_close
!   use SystemSettings, only : SYSTEM_STRING_LENGTH
!
!   implicit none
!   Type(TMolecule) :: this
!   character(SYSTEM_STRING_LENGTH),intent(in) :: filename
!
!   ! file format: Natom lines with atom names on them
!   integer :: nlines
!   integer :: hfile
!   integer :: i
!   character(SYSTEM_STRING_LENGTH) :: str  
!
!   nlines = io_count_file_lines(filename)
!
!   if (nlines /= this % Natoms ) then
!       write(error_message,*) 'Molecule_read_names: number of lines in file ',filename(1:20),' is ',nlines, &
!                              ', should be equivalent to number of atoms ', this % Natoms
!       call error_throw(ERROR_PARAMETER)
!       return
!   end if
!
!   if( this % hasNames ) then
!      deallocate(this % atomnames )
!   end if
!
!   allocate(this % atomnames(nlines))
!
!   hfile = io_open(filename,'r')
! 
!   do i=1,nlines
!      read(hfile,*) str
!      this % atomnames(i) = str(1:4)
!   end do
!
!   this % hasNames = .TRUE.
!
!   call io_close(hfile)
!
!end subroutine

subroutine Molecule_write(this,hfile) !DOC 
!DOC write the molecule to the mol file
   Type(TMolecule),intent(in) :: this
   integer, optional :: hfile ! hfile = 0 means stdout
   
   integer :: h      
   integer :: i

   if(.not.present( hfile  )) then
       h = 0
   else
       h = hfile
   end if
  
   do i = 1,this % Natoms
     
      if (h == 0) then
         write(*,*) this % atomnames(i), this % x(i), this % y(i), this % z(i) , &
                     this % sigma(i), this % epsilon(i), this % charge(i), this % mass(i)
      else
         write(h,*) this % atomnames(i), this % x(i), this % y(i), this % z(i) , &
                    this % sigma(i), this % epsilon(i), this % charge(i), this % mass(i)
      end if

   end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE Molecule


