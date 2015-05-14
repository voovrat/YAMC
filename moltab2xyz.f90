program moltab2xyz !DOC !PROG
!DOC !PROG !FILE convert the binary file with molecular positions to the xyz format
use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use MoleculeTable
use AtomicData
use io, only : io_open, io_close
!DOC !PROG Usage:  moltab2xyz  system.composition  input.moltab output.xyz [ BoxLength [ nbytes_xyz nbytes_ang ] ]'
!DOC !PROG Arguments:
!DOC !PROG  ::  input.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
!DOC !PROG  ::  nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
!DOC !PROG  ::  nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
!DOC !PROG  ::  BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 '



implicit none

character(SYSTEM_STRING_LENGTH) :: composition_file,binary_file,output_file
character(SYSTEM_STRING_LENGTH) :: tmp
integer :: nbytes_xyz,nbytes_ang
Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TAtomicData) :: atomic_data
integer :: i
integer :: hfile
real(8) :: V,BoxLength
real(8),parameter :: DENSITY = 33.3  ! particles/nm^3
integer,parameter :: MAX_ATOM = 10000

real(8), parameter :: kT_kcal_mol = 1 ! does not matter, because we do not save epsilon to xyz (ergo only coordinates matter )

if ( iargc() < 3 ) then

   write(*,*)  'Usage:  moltab2xyz  system.composition  input.moltab output.xyz [ BoxLength [ nbytes_xyz nbytes_ang ] ]'
   write(*,*)  '   input.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
   write(*,*)  '   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
   write(*,*)  '   BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 '


   stop
end if

call getarg(1,composition_file)
call getarg(2,binary_file)
call getarg(3,output_file)


if ( iargc() .ge. 4 ) then
 
   call getarg(4,tmp)
   read(tmp,*) BoxLength

else
  
   BoxLength = 0

end if

if ( iargc() .ge. 5 ) then

   call getarg(5,tmp)
   read(tmp,*) nbytes_xyz
 
   call getarg(6,tmp)
   read(tmp,*) nbytes_ang

else
   
   nbytes_xyz = 2
   nbytes_ang = 2

end if

!write(*,*)  'input=',input_file(1:20), 'output=',output_file(1:20), 'nbytes_xyz=',nbytes_xyz, 'nbytes_ang=',nbytes_ang

! read the system composition
call Composition_nulify(comp)
call Composition_read_from_file(comp, composition_file )

! read the molecule table
call MoleculeTable_nulify(mol_tab)

hfile = io_open(binary_file,'b')
call MoleculeTable_load_binary(mol_tab, hfile, comp, nbytes_xyz,nbytes_ang)
call io_close(hfile)


! calculate the BoxLength (if necessary) 

if ( BoxLength == 0 ) then

   ! DENSITY = N / V --->  V = N / DENSITY,  L = (V)^(1/3) 
   V = mol_tab % nmol / DENSITY
   BoxLength = V**(1.0/3)*10

end if

! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MAX_ATOM, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! save xyz 
call AtomicData_save_to_xyz(atomic_data,output_file)

! deallocate everything

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)


end program moltab2xyz
