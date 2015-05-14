program moltab_text2bin !DOC !PROG 
!DOC !PROG !FILE Convert the text file x y z theta phi psi to the binary moltab format
use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use MoleculeTable
use io, only : io_open, io_close
use string

!DOC !PROG Usage:  moltab_text2bin  system.composition  input.moltext BoxLength output.moltab [ nbytes_xyz nbytes_ang ] '
!DOC !PROG Arguments:
!DOC !PROG ::   input.moltext- in text format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
!DOC !PROG ::   BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 '
!DOC !PROG ::   output.molbin - bin  file where of the same format'
!DOC !PROG ::   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
!DOC !PROG ::   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '


implicit none

character(SYSTEM_STRING_LENGTH) :: composition_file,input_file,output_file,holes_file
character(SYSTEM_STRING_LENGTH) :: tmp
integer :: nbytes_xyz,nbytes_ang
Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
integer :: i,ii
integer :: hfile
real(8) :: V,BoxLength
integer :: MaxAtom,MaxMolecule


if ( iargc() < 4 ) then

   write(*,*)  'Usage:  moltab_text2bin  system.composition  input.moltext BoxLength output.moltab [ nbytes_xyz nbytes_ang ] '
   write(*,*)  '   input.moltext- in text format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 '
   write(*,*)  '   output.molbin - bin  file where of the same format'
   write(*,*)  '   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
   write(*,*)  '   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '

   stop
end if

call getarg(1,composition_file)
call getarg(2,input_file)
call getarg(3,tmp)
read(tmp,*) BoxLength

call getarg(4,output_file)



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

MaxAtom = Composition_count_atoms(comp)
MaxMolecule = Composition_count_molecules(comp)

! read the molecule table
call MoleculeTable_nulify(mol_tab)


hfile = io_open(input_file,'r')
call MoleculeTable_load_text(mol_tab, hfile, comp, BoxLength)
call io_close(hfile)

hfile = io_open(output_file,'b')
call MoleculeTable_save_binary(mol_tab,hfile,nbytes_xyz,nbytes_ang)
call io_close(hfile)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program moltab_text2bin
