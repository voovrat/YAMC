program mccont !DOC !PROG
!DOC !PROG !FILE Continue the MC simulation
use RhoSquared
use SumSinCosKR
use parameters
use MoleculeTable
use AtomicData
use LJTypes
use MC , only : initmc, movemc,clearmc,RealSpace_TotalSum,KSpace_TotalSum,mc_calc_total_sums_dbg,mc_calc_total_sums
use MC , only : U_tail_lj6,U_tail_lj12,U_tail_coulomb, U_total_ew,U_total_lj,U_total,U_ext
use io, only : write_integer_array,write_real_array,write_xyz_array,io_open,io_close,write_real_matrix
!use MCLuc, only : initmc_Luc,movemc_Luc,MCLuc_init_input,MCLucInput,piston ,write_energy_LUC
use runmc, only : run_mc,ifirst
use MRandom

!DOC !PROG Usage:  mccont  parameters.prm system.composition  [ extra_parameters ]
!DOC !PROG Arguments:
!DOC !PROG ::   system.composition - number and types of molecules in the system
!DOC !PROG ::  extra_parameters - coma separated  prm1=val1,prm2=val2,...   NO SPACES ALLOWED! 
!DOC !PROG   first frame is the last frame from the traj_file (given in parameters)'
!DOC !PROG first boxlength is taken from the last frame in frames_file
 

implicit none

character(SYSTEM_STRING_LENGTH) :: parameters_file,composition_file,moltab_file
character(SYSTEM_STRING_LENGTH) :: tmpstr
character(SYSTEM_STRING_LENGTH) :: extra_parameters_string
logical :: extra_parameters = .FALSE.

Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TAtomicData) :: atomic_data
Type(TLJTypes) :: lj_types

integer :: i,j
integer :: hfile,nk
integer :: MaxMolecule,MaxAtom
integer :: imol

logical :: accepted

real(8) :: rnd
real(8) :: sigma2,xlj
real(8),dimension(:),pointer :: pa,pb,pc
Type(TSumSinCosKR),pointer :: Coul,LJ
real(8) :: uu_lj

real(8),dimension(:),allocatable :: xx,yy,zz ! for test movemc

call parse_command_line
call read_input_files


call randomize(rnd_seed) 
call initmc(comp, mol_tab, atomic_data, lj_types )


write(*,*) 'Starting from the cycle',ifirst
call run_mc(mol_tab,max_cycle)

!call test_mcvol
  

!   write(*,*) accepted
!end do

!call clearmc

call deallocate_all


contains 


subroutine write_energy_ME !DOC
!DOC write the energy components

write(*,*) '***************  ME **************** '
uu_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 
write(*,*) 'U_REAL ew:',RealSpace_TotalSum % uu_ew,'lj:',uu_lj,'lj6:',&
                        RealSpace_TotalSum % uu_lj6,'lj12:',RealSpace_TotalSum % uu_lj12
write(*,*) 'U_FOUR ew:',KSpace_TotalSum % energy_coulomb,'lj:',KSpace_TotalSum % energy_LJ
write(*,*) 'U_TAIL ew:',U_tail_coulomb,'lj:',U_tail_lj12 - U_tail_lj6,'lj6:',U_tail_lj6,'lj12:',U_tail_lj12
write(*,*) 'U_EXT  ew:',U_ext
write(*,*) 'U_TOTL ew:',U_total_ew,'lj:',U_total_lj
write(*,*) 'U_TOT tot:',U_total

end subroutine

subroutine parse_command_line !DOC
!DOC read the command line arguments

if ( iargc() < 2 ) then

   write(*,*)  '  mccont: continue MC simulation'
   write(*,*)  'Usage:  mccont  parameters.prm system.composition  [ extra_parameters ]'
   write(*,*)  '   system.composition - number and types of molecules in the system'
   write(*,*)  '   extra_parameters - coma separated  prm1=val1,prm2=val2,...   NO SPACES ALLOWED! '
   write(*,*)  '- first frame is the last frame from the traj_file (given in parameters)'
   write(*,*)  '-first boxlength is taken from the last frame in frames_file'

   stop
end if

call getarg(1,parameters_file)
call getarg(2,composition_file)


if ( iargc() .ge. 3 ) then
    call getarg(3,extra_parameters_string)
    extra_parameters = .TRUE.
else
    extra_parameters = .FALSE.
end if


end subroutine


subroutine read_input_files !DOC
!DOC read the input files 
   use runmc, only : read_last_frame
   use SystemSettings, only : SEEK_SET

   integer :: nframes
   integer :: block_size

call Composition_nulify(comp)
call Composition_read_from_file(comp, composition_file )

MaxAtom = Composition_count_atoms(comp)
MaxMolecule = SUM( comp % mol_numbers(1:comp % n_types) )

! read and initialize the  parameters

if ( .not.extra_parameters ) then
     call Parameters_init(parameters_file,MaxMolecule)  ! MAxMolecule is used for the BoxLength calculation (from density, if required)
else
     call Parameters_init(parameters_file,MaxMolecule,extra_parameters_string)
end if

! read the frames file 
call read_last_frame(frames_file,nframes)
! recalc parameters for the new boxlength
call Parameters_recalc(BoxLength)


! read the molecule table
call MoleculeTable_nulify(mol_tab)


hfile = io_open(traj_file,'b') ! the last frame from the traj file

! use output_nbytes because they were used for storing the trajectory
block_size = (3* output_nbytes_xyz + 3*output_nbytes_ang) * MaxMolecule
call fseek(hfile, (nframes-1) * block_size, SEEK_SET) 

call MoleculeTable_load_binary( mol_tab, hfile, comp, output_nbytes_xyz, output_nbytes_ang )
call io_close(hfile)

! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )

end subroutine


subroutine deallocate_all !DOC
!DOC deallocate everything

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)

call LJTypes_dealloc(lj_types)


end subroutine


end program mccont
