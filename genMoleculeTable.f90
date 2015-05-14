program genMoleculeTable !DOC !PROG
use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use MoleculeTable
use io, only : io_open, io_close
use MRandom, only : randomize

!DOC !PROG !FILE generate the positions of the atoms for the given composition
!DOC !PROG Usage:  genMoleculeTable  system.composition  output.moltab  [ nbytes_xyz nbytes_ang  [random_seed] ]
!DOC !PROG Arguments:
!DOC !PROG  ::  output.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi) 
!DOC !PROG  ::  nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) 
!DOC !PROG  ::   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) 
!DOC !PROG :: random_seed - optional argument. The seed for the random number generator



implicit none

character(SYSTEM_STRING_LENGTH) :: input_file,output_file
character(SYSTEM_STRING_LENGTH) :: tmp
integer :: nbytes_xyz,nbytes_ang
Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
integer :: i
integer :: hfile
integer :: rnd_seed

if ( iargc() < 2 ) then

   write(*,*)  'Usage:  genMoleculeTable  system.composition  output.moltab  [ nbytes_xyz nbytes_ang  [random_seed] ]'
   write(*,*)  '   output.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
   write(*,*)  '   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '

   stop
end if

call getarg(1,input_file)
call getarg(2,output_file)

if ( iargc() .ge. 4 ) then

   call getarg(3,tmp)
   read(tmp,*) nbytes_xyz
 
   call getarg(4,tmp)
   read(tmp,*) nbytes_ang

else
   
   nbytes_xyz = 2
   nbytes_ang = 2

end if

if ( iargc() .ge. 5 ) then
   call getarg(5,tmp)
   read(tmp,*) rnd_seed
   call randomize(rnd_seed)
end if

!write(*,*)  'input=',input_file(1:20), 'output=',output_file(1:20), 'nbytes_xyz=',nbytes_xyz, 'nbytes_ang=',nbytes_ang

! read the input
call Composition_nulify(comp)
call Composition_read_from_file(comp, input_file )

! create the molecule table
call MoleculeTable_nulify(mol_tab)
call MoleculeTable_placeMoleculesToGrid(mol_tab, comp % mol_types, comp % mol_numbers, comp % n_types )

! save the molecule table
hfile = io_open(output_file,'b')
call MoleculeTable_save_binary(mol_tab,hfile,nbytes_xyz,nbytes_ang)
call io_close(hfile)

! deallocate everything

call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)


end program genMoleculeTable
