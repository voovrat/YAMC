program mc_accum_luc !DOC !PROG
!DOC !PROG !FILE calculate the projections
!DOC !PROG
!DOC !PROG  mc_accum_luc: do statistical processing of the simulation data'
!DOC !PROG Usage:  mc_accum_luc  parameters.prm system.composition output.proj first_frame[-last_frame] [step_size] [maxmn]
!DOC !PROG  :: system.composition - number and types of molecules in the system
!DOC !PROG ::   output.proj - projections file  
!DOC !PROG
!DOC !PROG  NOTE: quantities like pressure, compressibility etc are calculated incorrectly
!DOC !PROG  NOTE: density is not updates and thus is calculated incorrectly (always 0.0332891 particles/nm^3).
!DOC !PROG  to get a correct assymptote, the projections should be re-normalized to the real inverse density  mean(V^-1) which can be get from frames.dat file




use parameters
use MoleculeTable
use AtomicData
use LJTypes
use MCAccumLuc, only : accumu_init, accumu_set_frame, accumu, resul, write_projections

use runmc, only : read_traj_frame,n_mv_steps, n_mv_accepted


use EwaldSumKSpace
use EwaldSumRealSpace
use EwaldSumExternal
use RhoSquared
use SumSinCosKR

!use MC , only : initmc, movemc,clearmc,RealSpace_TotalSum,KSpace_TotalSum,mc_calc_total_sums_dbg,mc_calc_total_sums
!use MC , only : U_tail_lj6,U_tail_lj12,U_tail_coulomb, U_total_ew,U_total_lj,U_total,U_ext
!use io, only : write_integer_array,write_real_array,write_xyz_array,io_open,io_close,write_real_matrix
!use MCLuc, only : initmc_Luc,movemc_Luc,MCLuc_init_input,MCLucInput,piston ,write_energy_LUC
!use runmc, only : run_mc,ifirst
!use MRandom

implicit none

character(SYSTEM_STRING_LENGTH) :: parameters_file,composition_file,output_file
character(SYSTEM_STRING_LENGTH) :: tmpstr

character(SYSTEM_STRING_LENGTH) :: extra_parameters_string
logical :: extra_parameters = .FALSE.

Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TAtomicData),target :: atomic_data
Type(TLJTypes) :: lj_types

integer :: maxmn

integer :: i,j
integer :: hfile,nk
integer :: MaxMolecule,MaxAtom
integer :: hfram,stat,htraj
integer :: nfram
integer :: nskip,step_size,maxframe

Type(TEwaldSumKSpace),target :: kspace_sum
Type(TEwaldSumRealSpace) :: real_sum
Type(TEwaldSumExternal) :: ext_sum

!real(8),dimension(:),allocatable :: fx,fy,fz,fx_ext,fy_ext,fz_ext

Type(TSinCosKR),dimension(:),allocatable :: sincos_ew,sincos_lj



call parse_command_line
call read_input_files
call allocate_ewald

!allocate( fx(MaxAtom) )
!allocate( fy(MaxAtom) )
!allocate( fz(MaxAtom) )


call accumu_init(comp,BoxLength,maxmn)

! first frame is already read
nfram = 1
do while ( .TRUE. )
 
  
 if( (nfram .ge. nskip) .and. (mod(nfram,step_size) == 0 ) ) then
   write(*,*) 'Fram:',nfram

 
   call calc_sums
   call accumu_set_frame(comp,atomic_data,BoxLength,kspace_sum,real_sum,ext_sum,n_mv_steps,n_mv_accepted)
   call accumu

  else
     write(*,*) 'Fram',nfram,' ---> skiped'
  end if

!   call store_forces

   call load_frame
   
   if( stat /= 0) exit

   nfram = nfram + 1

   if ( nfram > maxframe) exit

end do

call resul

hfile = io_open(output_file,'w')
call write_projections(hfile)
call io_close(hfile)

call deallocate_all


contains 


!subroutine write_energy_ME
!
!write(*,*) '***************  ME **************** '
!uu_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 
!write(*,*) 'U_REAL ew:',RealSpace_TotalSum % uu_ew,'lj:',uu_lj,'lj6:',&
!                        RealSpace_TotalSum % uu_lj6,'lj12:',RealSpace_TotalSum % uu_lj12
!write(*,*) 'U_FOUR ew:',KSpace_TotalSum % energy_coulomb,'lj:',KSpace_TotalSum % energy_LJ
!write(*,*) 'U_TAIL ew:',U_tail_coulomb,'lj:',U_tail_lj12 - U_tail_lj6,'lj6:',U_tail_lj6,'lj12:',U_tail_lj12
!write(*,*) 'U_EXT  ew:',U_ext
!write(*,*) 'U_TOTL ew:',U_total_ew,'lj:',U_total_lj
!write(*,*) 'U_TOT tot:',U_total
!
!end subroutine

subroutine parse_command_line !DOC
!DOC read the command line arguments
 use string

if ( iargc() < 4 ) then

   write(*,*)  '  mc_accum_luc: do statistical processing of the simulation data'
   write(*,*)  'Usage:  mc_accum_luc  parameters.prm system.composition output.proj first_frame[-last_frame] [step_size] [maxmn]'
   write(*,*)  '   system.composition - number and types of molecules in the system'
   write(*,*)  '   output.proj - projections file  '

   stop
end if

call getarg(1,parameters_file)
call getarg(2,composition_file)
call getarg(3,output_file)
call getarg(4,tmpstr)

if ( str_get_next_pos(tmpstr,1,'-') > 0) then
  call str_subs(tmpstr,'-',' ')
  read(tmpstr,*) nskip,maxframe
else
   read(tmpstr,*) nskip
   maxframe=-1
end if

if ( iargc() >= 5 ) then
   call getarg(5,tmpstr)
   read(tmpstr,*) step_size
else
   step_size = 1
end if

if ( iargc() >=6 ) then
   call getarg(6,tmpstr)
   read(tmpstr,*) maxmn
else
   maxmn = 4
end if

extra_parameters = .FALSE.
!call getarg(3,output_file)

!if ( iargc() .ge. 4 ) then
!    call getarg(4,extra_parameters_string)
!    extra_parameters = .TRUE.
!else
!    extra_parameters = .FALSE.
!end if


end subroutine


subroutine read_input_files !DOC
!DOC read input files
   use SystemSettings, only : SEEK_SET
   use io, only : io_open,io_close


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

hfram = io_open(frames_file,'r')

! read the frames file 
call read_traj_frame(hfram,stat)
! recalc parameters for the new boxlength
call Parameters_recalc(BoxLength)

! read the molecule table
call MoleculeTable_nulify(mol_tab)
htraj = io_open(traj_file,'b') 
call MoleculeTable_load_binary( mol_tab, htraj, comp, output_nbytes_xyz, output_nbytes_ang )
!call io_close(hfile)

! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )

end subroutine

subroutine allocate_ewald !DOC 
!DOC allocate the ewald sums 
   integer :: i

   call EwaldSumRealSpace_alloc_full(real_sum, atomic_data, lj_types )
   call EwaldSumKSpace_alloc(kspace_sum,atomic_data,lj_types)
   call EwaldSumExternal_init(ext_sum, atomic_data % xx, atomic_data % yy,atomic_data % zz,atomic_data % charge,MaxAtom)

   allocate( sincos_ew(MaxAtom) )
   allocate( sincos_lj(MaxAtom) )

   do i=1,MaxAtom

      call SinCosKR_alloc( sincos_ew(i), kspace_sum % grid )
      call SinCosKR_alloc( sincos_lj(i), kspace_sum % grid )

   end do

   !allocate(fx_ext(MaxAtom))
   !allocate(fy_ext(MaxAtom)) 
   !allocate(fz_ext(MaxAtom))

end subroutine

subroutine deallocate_ewald !DOC 
!DOC deallocate the ewald sums

  call EwaldSumRealSpace_dealloc_full( real_sum )
  call EwaldSumKSpace_dealloc(kspace_sum )

  do i=1,MaxAtom
    call SinCosKR_dealloc( sincos_ew(i) )
    call SinCosKR_dealloc( sincos_lj(i) )
  end do

  deallocate( sincos_ew )
  deallocate( sincos_lj )  

 ! deallocate(fx_ext)
 ! deallocate(fy_ext)
 ! deallocate(fz_ext)

end subroutine


subroutine calc_sums !DOC 
!DOC calculate the ewald sums
  use ForceKSpace
 
  integer :: i,ityp

  real(8),dimension(:),pointer :: xx,yy,zz,charge
  Type(TSumSinCosKR),pointer :: sumsincos_ew
  Type(TSumSinCosKR),dimension(:),pointer :: sumsincos_lj

  real(8) :: fx_kspace_coulomb,fy_kspace_coulomb,fz_kspace_coulomb
  real(8) :: fx_kspace_LJ,fy_kspace_LJ,fz_kspace_LJ
  logical :: overlap


! Real Space forces
  call EwaldSumRealSpace_calc( real_sum,overlap )

  if ( overlap ) then
     write(*,*) 'OVERLAP! WTF?'
  end if

!  fx(1:MaxAtom) = real_sum % Fx(1:MaxAtom)
!  fy(1:MaxAtom) = real_sum % Fy(1:MaxAtom)
!  fz(1:MaxAtom) = real_sum % Fz(1:MaxAtom)

! KSpace forces 
  call EwaldSumKSpace_initBeta( kspace_sum )
  call EwaldSumKSpace_initBetaLJ( kspace_sum )

  call RhoSquared_nulify( kspace_sum % rho_squared_total )

  xx => atomic_data % xx
  yy => atomic_data % yy
  zz => atomic_data % zz
  charge => atomic_data % charge

!  sumsincos_coulomb_old => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
!   sumsincos_LJ_old => KSpace_TotalSum % rho_squared_total % sumsincos_LJ

  ! Calculate the total sum saving the individual summands (which are used in force calculation)
  do i=1,MaxAtom ! local atom index

      ityp = lj_types % type_by_index(i)

!!      write(*,*) 'addAtom: i',i,'ityp',ityp,'ii',ii,'imol',imol
      call RhoSquared_addAtom( kspace_sum % rho_squared_total, xx(i), yy(i), zz(i), charge(i), ityp, &
                               sincos_ew(i), sincos_lj(i) )
   end do
!
!   sumsincos_ew => kspace_sum % rho_squared_total % sumsincos_coulomb
!   sumsincos_lj => kspace_sum % rho_squared_total % sumsincos_LJ

!   ! calculating the forces 
! !  do i=1,MaxAtom
! 
! !     ityp = lj_types % type_by_index(i)
!
! !     call ForceKSpace_coulomb( kspace_sum % beta , sumsincos_ew,  sincos_ew(i),&
!                                fx_kspace_coulomb, fy_kspace_coulomb, fz_kspace_coulomb ) 
!    
!     if ( ityp > 0) then
!         call ForceKSpace_LJ(kspace_sum % beta_LJ, sumsincos_lj, lj_types % NType, ityp,  sincos_lj(i),&
!                                fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ )
!
!      else
!          fx_kspace_LJ = 0
!          fy_kspace_LJ = 0
!          fz_kspace_LJ = 0
!      end if 
!
!      fx(i) = fx(i) + fx_kspace_coulomb + fx_kspace_LJ
!      fy(i) = fy(i) + fy_kspace_coulomb + fy_kspace_LJ
!      fz(i) = fz(i) + fz_kspace_coulomb + fz_kspace_LJ
!
!   end do  
!
! External Part

  call EwaldSumExternal_calc_mu(ext_sum)
!  call external_sum_calc_forces( ext_sum % mu_x, ext_sum % mu_y,ext_sum %  mu_z, charge, MaxAtom, fx_ext,fy_ext,fz_ext )

!  fx(1:MaxAtom) = fx(1:MaxAtom) + fx_ext(1:MaxAtom)
!  fy(1:MaxAtom) = fy(1:MaxAtom) + fy_ext(1:MaxAtom)
!  fz(1:MaxAtom) = fz(1:MaxAtom) + fz_ext(1:MaxAtom) 

end subroutine

subroutine load_frame !DOC
!DOC load next trajectory frame
   use runmc, only : read_traj_frame

real(8) :: old_box_length
!hfram = io_open(frames_file,'r')

old_box_length = BoxLength

! read the frames file 
call read_traj_frame(hfram,stat)
! recalc parameters for the new boxlength

if (stat /=0 ) return


call Parameters_recalc(BoxLength)

! read the molecule table
!call MoleculeTable_nulify(mol_tab)
!htraj = io_open(traj_file,'b') 
call MoleculeTable_load_binary( mol_tab, htraj, comp, output_nbytes_xyz, output_nbytes_ang )
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
!call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )
call LJTypes_scale_sigma(lj_types, old_box_length / BoxLength )

end subroutine

!subroutine store_forces
!  use SystemSettings, only : SEEK_END
!  use FloatingPoint, only : write_float
!  use parameters, only : BoxLength 
!  integer :: hfile,i
!
!  hfile = io_open(output_file,'b')
!  call fseek(hfile,0,SEEK_END)
!
!  do i=1,MaxAtom
!
!    call write_float(hfile,fx(i)/BoxLength)
!    call write_float(hfile,fy(i)/BoxLength)
!    call write_float(hfile,fz(i)/BoxLength)
!
!  end do
!
!  call io_close(hfile)
!
!end subroutine
!

subroutine deallocate_all !DOC 
!DOC deallocate everything


call deallocate_ewald

!deallocate( fx )
!deallocate( fy )
!deallocate( fz )

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)

call LJTypes_dealloc(lj_types)


end subroutine


end program mc_accum_luc
