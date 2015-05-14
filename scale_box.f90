Module ScaleBox !DOC
!DOC !FILE change the size of the box (and thus the atom coordinates)

contains

subroutine scale_box( this, mol_tab, NewBoxLength )  !DOC 
!DOC change the size of the box. The centers of mass remain the same in the relative coordinates, but the bond length change.
   use AtomicData
   use MoleculeHandler
   use MoleculeTable 
   use Molecule
   use geometry, only : center_of_mass
!DOC Parameters:
   Type(TAtomicData),target :: this !DOC AtomicData (coordinates)
   Type(TMoleculeTable),intent(in) :: mol_tab !DOC MoleculeTable
   real(8),intent(in) :: NewBoxLength !DOC new length of the box
   real(8),dimension(:),pointer :: xmol,ymol,zmol  

   integer :: i,j,first_atom,last_atom,natom,mol_type
   real(8) :: x_center,y_center,z_center
   Type(TMolecule),pointer :: mol_ptr
 
   real(8) :: scale_coeff 

   scale_coeff =  this % BoxLength / NewBoxLength

   do i = 1,mol_tab % nmol
 
      first_atom = mol_tab % first_atom(i)
      last_atom = mol_tab % last_atom(i)
      natom = last_atom - first_atom + 1

      if ( natom == 1 ) cycle
      
      mol_type = mol_tab % mol_type(i)

      xmol => this % xx(first_atom:last_atom)
      ymol => this % yy(first_atom:last_atom)
      zmol => this % zz(first_atom:last_atom)

      call MoleculeHandler_getMolecule(mol_type,mol_ptr)
      call center_of_mass(xmol,ymol,zmol,mol_ptr % mass,natom, x_center, y_center, z_center)

      ! Coordinates of the center of mass should be the same
      ! but atomic coordinates will be scaled:
      ! OA [ NewBoxLength ] = 0A [BoxLength ] * BoxLength / NewBoxLength
     
      do j=first_atom,last_atom

         this % xx(j) =  ( this % xx(j) - x_center ) * scale_coeff + x_center
         this % yy(j) =  ( this % yy(j) - y_center ) * scale_coeff + y_center
         this % zz(j) =  ( this % zz(j) - z_center ) * scale_coeff + z_center

      end do

   end do

   this % sigma(:) = this % sigma(:) * scale_coeff

   this % BoxLength = NewBoxLength

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module ScaleBox
