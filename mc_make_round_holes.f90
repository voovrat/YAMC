program mc_make_round_holes !DOC !PROG
!DOC !PROG !FILE  create spherical holes in the given box with molecules
use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use MoleculeTable
use MoleculeHandler
use Molecule
use AtomicData
use io, only : io_open, io_close
use string

!DOC !PROG Usage:  mc_make_round_holes  system.composition  input.moltab BoxLength holes.txt output.moltab  [ nbytes_xyz nbytes_ang ] 

!DOC !PROG Arguments:

!DOC !PROG ::   input.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  
!DOC !PROG ::   BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 
!DOC !PROG ::   holes.txt  -  first line - number of holes, next lines:  x y z R  
!DOC !PROG ::   output.moltab - binary file where the molecules which intersect with holes
!DOC !PROG ::   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) 
!DOC !PROG ::   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) 


implicit none

character(SYSTEM_STRING_LENGTH) :: composition_file,input_file,output_file,holes_file
character(SYSTEM_STRING_LENGTH) :: tmp
integer :: nbytes_xyz,nbytes_ang
Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TMoleculeTable) :: mol_tab_out
Type(TAtomicData) :: atomic_data
integer :: i,ii
integer :: hfile
real(8) :: V,BoxLength
real(8),parameter :: DENSITY = 33.3  ! particles/nm^3
integer :: MaxAtom,MaxMolecule

real(8), parameter :: kT_kcal_mol = 1 ! does not matter, because we do not save epsilon to xyz (ergo only coordinates matter )

integer :: Nholes 
real(8), dimension(:), allocatable :: x_hole,y_hole,z_hole,r_hole
logical, dimension(:), allocatable :: remove_list

integer :: rest_count 

Type(TMolecule),pointer :: mol_ptr
integer :: mol_id, last_mol_id, mol_count

if ( iargc() < 5 ) then

   write(*,*)  'Usage:  mc_make_round_holes  system.composition  input.moltab BoxLength holes.txt output.moltab  [ nbytes_xyz nbytes_ang ] '
   write(*,*)  '   input.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   BoxLength  - in Angstroems. If not given - is recalculated to 33.3 particles/nm^3 '
   write(*,*)  '   holes.txt  -  first line - number of holes, next lines:  x y z R  '
   write(*,*)  '   output.moltab - binary file where the molecules which intersect with holes'
   write(*,*)  '   nbytes_xyz - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '
   write(*,*)  '   nbytes_ang - bytes per fixed point for coordinats ( default: 2, allowed: 1,2,3) '

   stop
end if

call getarg(1,composition_file)
call getarg(2,input_file)
call getarg(3,tmp)
read(tmp,*) BoxLength

call getarg(4,holes_file)
call getarg(5,output_file)



if ( iargc() .ge. 6 ) then

   call getarg(6,tmp)
   read(tmp,*) nbytes_xyz
 
   call getarg(7,tmp)
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


hfile = io_open(input_file,'b')
call MoleculeTable_load_binary(mol_tab, hfile, comp, nbytes_xyz,nbytes_ang)
call io_close(hfile)


! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

call read_holes_file

allocate(remove_list(MaxMolecule))

call fill_remove_list

rest_count = 0

do i=1,MaxMolecule

!   if ( remove_list(i) ) write(*,*) i,mol_tab % x(i) * BoxLength, mol_tab % y(i) * BoxLength, mol_tab % z(i) * BoxLength 
   if ( .NOT. remove_list(i) ) rest_count = rest_count + 1

end do

!write(*,*) 'rest_count=',rest_count

call MoleculeTable_alloc(mol_tab_out,rest_count)

mol_tab_out % nmol = rest_count
ii=1
last_mol_id = -1
mol_count = 0

write(*,*) comp % n_types 

do i=1,MaxMolecule

  if( .NOT. remove_list(i) ) then

     mol_id = mol_tab % mol_type(i)
    

     mol_count = mol_count + 1

     mol_tab_out % x(ii) = mol_tab % x(i)
     mol_tab_out % y(ii) = mol_tab % y(i)
     mol_tab_out % z(ii) = mol_tab % z(i)

     mol_tab_out % theta(ii) = mol_tab % theta(i)
     mol_tab_out % phi(ii) = mol_tab % phi(i)
     mol_tab_out % psi(ii) = mol_tab % psi(i)

     ii = ii + 1

     if ( (mol_id .ne. last_mol_id) .and. (last_mol_id > 0) .or. (i == MaxMolecule) ) then
    
         if ( i== MaxMolecule )  last_mol_id = mol_id
             
         call MoleculeHandler_getMolecule(last_mol_id,mol_ptr)


         write(*,*)  mol_ptr % filename ( 1:str_get_next_pos(mol_ptr % filename,1,' ') ),' ',mol_count

         mol_count = 1

         last_mol_id = mol_id

     end if ! write comp


  end if ! .NOT.remove_list


end do


hfile = io_open(output_file,'b')
call MoleculeTable_save_binary(mol_tab_out,hfile,nbytes_xyz,nbytes_ang)
call io_close(hfile)



! save xyz 
!call AtomicData_save_to_xyz(atomic_data,output_file)

! deallocate everything

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)

contains


subroutine read_holes_file !DOC
!DOC read the file with holes
!DOC first line - number of holes
!DOC next lines - x y z R  (in angstroems)
   integer :: hfile
   integer :: i  

   hfile = io_open(holes_file,'r')

   read(hfile,*)  Nholes
   
   allocate(x_hole(Nholes))
   allocate(y_hole(Nholes))
   allocate(z_hole(Nholes))
   allocate(R_hole(Nholes))

   do i=1,Nholes

       read(hfile,*)  x_hole(i),y_hole(i),z_hole(i),R_hole(i)

   end do
   
   call io_close(hfile)

end subroutine

subroutine fill_remove_list !DOC
!DOC find indeces of molecules to be removed

   integer :: i,j
   integer :: imol
   real(8) :: rc, r_HS
   real(8) :: dx,dy,dz

   remove_list(:) = .FALSE.

   do i = 1,MaxAtom

       r_HS = atomic_data % hard_core_angstr(i) / 2
       imol = atomic_data % molnum_by_atomnum(i) 

 !      write(*,*) i,Nholes

       do j = 1,Nholes

 !         write(*,*) x_hole(j),y_hole(j),z_hole(j),R_hole(j)
 
          dx = atomic_data % xx(i) * BoxLength - x_hole(j)
          dy = atomic_data % yy(i) * BoxLength - y_hole(j)
          dz = atomic_data % zz(i) * BoxLength - z_hole(j) 

          rc = SQRT( dx**2 + dy**2 + dz**2 )
 
!          write(*,*) dx,dy,dz,rc

          if  ( rc < r_hole(j) + r_HS )  remove_list(imol) = .TRUE.

       end do     

   end do

end subroutine fill_remove_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program mc_make_round_holes
