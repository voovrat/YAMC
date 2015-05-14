Module composition !DOC 
use Molecule
use SystemSettings, only : SYSTEM_STRING_LENGTH
use io, only : io_open, io_close 
use MoleculeHandler 
!DOC !FILE the module contains the data structure and the functions to deal with the system composition (numbers of the molecules of each type)


implicit none
Type :: TComposition  !DOC
!DOC Fields:
    integer :: n_types  !DOC number of different types of molecules
    integer,dimension(:),pointer ::  mol_types   !DOC types (indeces in MoleculeHandler array) of the molecules.  
    integer,dimension(:),pointer ::  mol_numbers !DOC numbers of the molecules of each type
    integer :: nalloc !DOC allocated size of the arrays

End Type TComposition


contains

subroutine Composition_nulify(this) !DOC
!DOC set the composition to the "zero" state. Is used before the Composition_read_from_file, where the arrays are allocated automatically.
!DOC Parameters:
   Type(TComposition) :: this  !DOC composition structure

   this % nalloc = 0
   this % n_types = 0

   nullify(this % mol_types)
   nullify(this % mol_numbers)

end subroutine

subroutine Composition_alloc(this,nalloc) !DOC
!DOC allocate the composition structure 
!DOC Parameters:
   Type(TComposition) :: this     !DOC  composition structure
   integer,intent(in) :: nalloc   !DOC  size of the arrays to be allocated

   call Composition_nulify(this)
   
   allocate( this % mol_types( nalloc ) )
   allocate( this % mol_numbers( nalloc ) )

   this % nalloc = nalloc

end subroutine


subroutine Composition_dealloc(this) !DOC
!DOC deallocate the composition structure
!DOC Parameters:
   Type(TComposition) :: this   !DOC  composition structure

   deallocate( this % mol_types )
   deallocate( this % mol_numbers ) 

   call Composition_nulify( this )

end subroutine

subroutine Composition_dealloc_molecules(this) !DOC
!DOC deallocate the molecules which were allocated when read form the file. Can be used before the Composition_dealloc
!DOC Parameters:
   Type(TComposition) :: this !DOC composition structure
   integer :: i
   Type(TMolecule),pointer :: mol_ptr
   integer :: typ

   do i=1,this % n_types

      typ = this % mol_types(i)

      call MoleculeHandler_getMolecule(typ,mol_ptr)
      call Molecule_deallocate(mol_ptr)
      call MoleculeHandler_release(typ)
   
   end do
     

end subroutine

subroutine Composition_read_from_file(this, fname, dont_load_molecules ) !DOC
!DOC read the composition from file
!DOC The format of the composition file: first line - number of species, 
!DOC each next line - name of the mol file  and the number of molecules.
!DOC Parameters:
   Type(TComposition) :: this  !DOC composition structure
   character(SYSTEM_STRING_LENGTH),intent(in) :: fname  !DOC name of the composition file
   logical,intent(in),optional :: dont_load_molecules  !DOC flag which shows that the structures of the molecules should not be read 
integer :: h
   character(SYSTEM_STRING_LENGTH) :: mol_file
   integer :: nmol
   integer :: i
   integer :: mol_id
   integer :: ntypes
   

   Type(TMolecule),pointer :: mol_ptr



   h = io_open(fname,'r')
   
   read(h,*)  ntypes

   if ( ntypes .gt. this % nalloc ) then

      if( this % nalloc .gt. 0 ) then
         call Composition_dealloc( this ) 
      end if 
      
      call Composition_alloc( this, ntypes )
      
   end if      

   do i=1,ntypes

      read(h,*)  mol_file,nmol
     
      if( ( .not. present(dont_load_molecules) ) .or. ( .not. dont_load_molecules ) ) then

         mol_id = MoleculeHandler_getFreeHandler() 
         call MoleculeHandler_occupy( mol_id, mol_ptr )  ! just gives the reference to mol_list(mol_id) 
 
         call Molecule_read_from_file( mol_ptr, mol_file )
         call Molecule_centrate( mol_ptr )

      end if

      this % mol_types(i) = mol_id
      this % mol_numbers(i) = nmol


   end do

   this % n_types = ntypes

   call io_close(h)    

end subroutine

function Composition_count_molecules( this ) !DOC
!DOC Count the molecules in the composition
   integer :: Composition_count_molecules
!DOC Parameters:
   Type(TComposition),intent(in) :: this !DOC composition structure
!DOC Return value:
!DOC   number of the molecules in the composition
   integer :: nmol

   nmol = SUM( this % mol_numbers( 1 : this % n_types ) )
   Composition_count_molecules = nmol
 
end function


function Composition_count_atoms( this ) !DOC
!DOC Count the total number of atoms in the composition
   use MoleculeHandler
   use Molecule
   integer :: Composition_count_atoms
! return the number of atoms in the molecule table
!DOC Parameters:
   Type(TComposition),intent(in) :: this !DOC the composition structure
!DOC Return value:
!DOC  the total number of atoms in all molecules in the composition 
   integer :: i, typ, natom
   Type(TMolecule),pointer :: mol_ptr  

   natom = 0
   do i=1,this % n_types

      typ = this % mol_types(i)
      call MoleculeHandler_getMolecule(typ, mol_ptr) 

      natom = natom + mol_ptr % NAtoms * this % mol_numbers(i)

   end do
  
   Composition_count_atoms = natom

end function



End Module composition


