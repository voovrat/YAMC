Module MC !DOC
!DOC !FILE The Module which defines the functions and data structures needed to perform the monte carlo steps 
use EwaldSumRealSpace
use EwaldSumKSpace
use ForceKSpace
use SystemSettings, only : TRealArrayPointer 
use MonteCarloMove
use EwaldSumExternal
use MoleculeTable
use RhoSquared
use SumSinCosKR
use EwaldSumTails
use composition , only : TComposition
! Monte Carlo main

implicit none

Type(TComposition),pointer :: mc_system_composition

Type(TEwaldSumRealSpace),target :: RealSpace_TotalSum
Type(TEwaldSumKSpace),target :: KSpace_TotalSum 

Type(TEwaldSumRealSpace) :: RealSpace_MolSum_old
Type(TEwaldSumRealSpace) :: RealSpace_MolSum_new

Type(TEwaldSumRealSpace) :: RealSpace_MolSum_old1,RealSpace_MolSum_old2 ! for mc_xchange
Type(TEwaldSumRealSpace) :: RealSpace_MolSum_new1,RealSpace_MolSum_new2  
Type(TEwaldSumRealSpace) :: RealSpace_MolSum_old12,RealSpace_MolSum_new12


Type(TEwaldSumExternal) :: Ext_TotalSum
! for move we have only one sum
Type(TEwaldSumExternal) :: Ext_MolSum_old, Ext_MolSum_new

! for xchange we need two sums
Type(TEwaldSumExternal) :: Ext_MolSum_old1, Ext_MolSum_new1
Type(TEwaldSumExternal) :: Ext_MolSum_old2, Ext_MolSum_new2

real(8) :: U_tail, U_tail_coulomb,U_tail_lj6,U_tail_lj12
real(8) :: U_ext
real(8) :: U_total,U_total_ew, U_total_lj

!real(8) :: RealSpace_MolEnergy_old ! [U_{real}]_M^{old} = sum_{j in M} u_j

!real(8), dimension(:) :: fx_old,fy_old,fz_old
!real(8), dimension(:) :: fx_new,fy_new,fz_new
! in k space where is no difference between one molecule and two molecules: the sums will have more terms that's it
Type(TRhoSquared),target :: sumsincos_mol_old
Type(TRhoSquared),target :: sumsincos_mol_new
Type(TRhoSquared),pointer :: delta_sumsincos   ! use the same memory as sumsincos_mol_new
Type(TRhoSquared),target :: sumsincos_total_new


Type(TComposition),pointer :: mc_composition

Type(TSinCosKR),dimension(:),allocatable :: sincos_kr_coulomb, sincos_kr_LJ  ! q_i sin(sx x kx + sy y ky + sz z kz) for i in Molecule
                  ! used in the kspace - force calculations 
! for xchange we don't have k-forces --> don't need these arrays

real(8),dimension(:),allocatable :: xnew,ynew,znew  ! new molecule coordinates

real(8),dimension(:),allocatable :: xnew1,ynew1,znew1 , xnew2,ynew2,znew2 ! for xchange
real(8),dimension(:),allocatable :: xtmp1,ytmp1,ztmp1 , xtmp2,ytmp2,ztmp2 ! for temporary storage of old coordinates 
                                                                          ! before the calculation of RealSum_new

real(8),dimension(:),allocatable :: fx_old,fy_old,fz_old ! old forces (sum of real and kspace parts ( recalculated )
real(8),dimension(:),allocatable :: fx_ext_old,fy_ext_old,fz_ext_old 
real(8),dimension(:),allocatable :: fx_new,fy_new,fz_new ! old forces (sum of real and kspace parts ( recalculated )
real(8),dimension(:),allocatable :: fx_ext_new,fy_ext_new,fz_ext_new 

integer :: MaxMolAtom


Type :: TMCState !DOC
!DOC The structure is used to save the current MC state before the change of the volume
!DOC Fields:
   real(8) :: kspace_energy, kspace_energy_coulomb, kspace_energy_lj !DOC KSpace energy components
   real(8) :: rspace_uu_ew, rspace_uu_lj6, rspace_uu_lj12  !DOC Real space energy components
   real(8) :: rspace_uu_ew_intra, rspace_uu_lj6_intra, rspace_uu_lj12_intra!DOC Real space intramolecular components

   real(8) :: U_tail, U_tail_coulomb,U_tail_lj6,U_tail_lj12 !DOCtotal energies
   real(8) :: U_ext  !DOC --
   real(8) :: U_total,U_total_ew, U_total_lj !DOC --


   real(8),dimension(:),allocatable :: xx,yy,zz !DOC coordinates
   real(8),dimension(:),allocatable :: sigma !DOC sigmas
   real(8),dimension(:),allocatable :: fx,fy,fz,uu !DOC real-forces and energy per atom
   Type(TRhoSquared) :: rho_squared_total !DOC SUM of sin(kr), cos(kr)

   real(8) :: mu_x,mu_y,mu_z !DOC mu=SUM q_i*r_i

End Type

Type(TMCState) :: mc_stored_state


contains

subroutine MCState_alloc(this, natom, grid, lj_types ) !DOC
!DOC Allocate the storage for MCState structure
!DOC Parameters 
   use FourierGrid, only : TFourierGrid
   use LJTypes, only : TLJTypes
   use RhoSquared, only : RhoSquared_alloc   

   Type(TMCState) :: this !DOC MCState
   integer,intent(in) :: natom !DOC number of atoms
   Type(TFourierGrid) :: grid !DOC KSpace grid
   Type(TLJTypes) :: lj_types !DOC LJTypes

   allocate( this % xx(natom) )
   allocate( this % yy(natom) )
   allocate( this % zz(natom) )

   allocate( this % sigma(natom) )
   
   allocate( this % fx(natom) )
   allocate( this % fy(natom) )
   allocate( this % fz(natom) )
  
   allocate( this % uu(natom) )

   call RhoSquared_alloc( this % rho_squared_total, grid, lj_types )

end subroutine

subroutine MCState_dealloc(this) !DOC 
!DOC Deallocate the MCState
   use RhoSquared, only : RhoSquared_dealloc
   Type(TMCState) :: this !DOC MCState
    
   deallocate( this % xx )
   deallocate( this % yy )
   deallocate( this % zz )
 
   deallocate( this % sigma )
 
   deallocate( this % fx )
   deallocate( this % fy )
   deallocate( this % fz )
   deallocate( this % uu )
  
   call RhoSquared_dealloc( this % rho_squared_total )
   
end subroutine


subroutine mc_save_state( state ) !DOC 
!DOC Save the current state. Used in volume change
!DOC Parameters: 
   use RhoSquared, only : RhoSquared_copy
   Type(TMCState) :: state !DOC MCState structure

   Type(TAtomicData),pointer :: atomic_data


   state % U_tail         =   U_tail
   state % U_tail_coulomb =   U_tail_coulomb
   state % U_tail_lj6     =   U_tail_lj6
   state % U_tail_lj12    =   U_tail_lj12
   state % U_ext          =   U_ext
   state % U_total        =   U_total
   state % U_total_ew     =   U_total_ew
   state % U_total_lj     =   U_total_lj       



   state % kspace_energy         = KSpace_TotalSum % energy
   state % kspace_energy_coulomb = KSpace_TotalSum % energy_coulomb
   state % kspace_energy_lj      = KSpace_TotalSum % energy_lj

   state % rspace_uu_ew    = RealSpace_TotalSum % uu_ew
   state % rspace_uu_lj6   = RealSpace_TotalSum % uu_lj6
   state % rspace_uu_lj12  = RealSpace_TotalSum % uu_lj12

   state % rspace_uu_ew_intra    = RealSpace_TotalSum % uu_ew_intra
   state % rspace_uu_lj6_intra   = RealSpace_TotalSum % uu_lj6_intra
   state % rspace_uu_lj12_intra  = RealSpace_TotalSum % uu_lj12_intra


   atomic_data => RealSpace_TotalSum % atomic_data

   state % xx(:) = atomic_data % xx(:)
   state % yy(:) = atomic_data % yy(:)
   state % zz(:) = atomic_data % zz(:)

   state % sigma(:) = atomic_data % sigma(:) 

   state % fx(:) = RealSpace_TotalSum % fx(:)
   state % fy(:) = RealSpace_TotalSum % fy(:)
   state % fz(:) = RealSpace_TotalSum % fz(:)

   state % uu(:) = RealSpace_TotalSum % uu(:)

   call RhoSquared_copy( state % rho_squared_total, KSpace_TotalSum % rho_squared_total )

   state % mu_x = Ext_TotalSum % mu_x 
   state % mu_y = Ext_TotalSum % mu_y
   state % mu_z = Ext_TotalSum % mu_z 

end subroutine


subroutine mc_restore_state( state ) !DOC
!DOC restore the state from the saved ones
   use RhoSquared, only : RhoSquared_copy
!DOC Parameters:
   Type(TMCState) :: state !DOC MCState structure
   Type(TAtomicData),pointer :: atomic_data  

   U_tail         =  state % U_tail
   U_tail_coulomb =  state % U_tail_coulomb
   U_tail_lj6     =  state % U_tail_lj6
   U_tail_lj12    =  state % U_tail_lj12
   U_ext          =  state % U_ext
   U_total        =  state % U_total
   U_total_ew     =  state % U_total_ew
   U_total_lj     =  state % U_total_lj       

   KSpace_TotalSum % energy        =  state % kspace_energy          
   KSpace_TotalSum % energy_coulomb=  state % kspace_energy_coulomb  
   KSpace_TotalSum % energy_lj     =  state % kspace_energy_lj       

   RealSpace_TotalSum % uu_ew   = state % rspace_uu_ew     
   RealSpace_TotalSum % uu_lj6  = state % rspace_uu_lj6    
   RealSpace_TotalSum % uu_lj12 = state % rspace_uu_lj12   

   RealSpace_TotalSum % uu_ew_intra   = state % rspace_uu_ew_intra     
   RealSpace_TotalSum % uu_lj6_intra  = state % rspace_uu_lj6_intra    
   RealSpace_TotalSum % uu_lj12_intra = state % rspace_uu_lj12_intra   


   atomic_data => RealSpace_TotalSum % atomic_data

   atomic_data % xx(:) = state % xx(:)
   atomic_data % yy(:) = state % yy(:)
   atomic_data % zz(:) = state % zz(:)

   atomic_data % sigma(:) = state % sigma(:)

   RealSpace_TotalSum % fx(:) = state % fx(:)
   RealSpace_TotalSum % fy(:) = state % fy(:) 
   RealSpace_TotalSum % fz(:) = state % fz(:)

   RealSpace_TotalSum % uu(:) = state % uu(:)

   call RhoSquared_copy( KSpace_TotalSum % rho_squared_total, state % rho_squared_total)

   Ext_TotalSum % mu_x = state % mu_x
   Ext_TotalSum % mu_y = state % mu_y
   Ext_TotalSum % mu_z = state % mu_z 


end subroutine

subroutine mc_alloc(comp,mol_table,atomic_data,lj_types) !DOC 
!DOC allocate the arrays needed for the MC step 
!DOC Parameters:
!  use io, only : write_real_array
  implicit none
  Type(TComposition),intent(in),target :: comp !DOC composition
  Type(TMoleculeTable),intent(in) :: mol_table !DOC MoleculeTable
  Type(TAtomicData),intent(in) :: atomic_data !DOC AtomicData
  Type(TLJTypes),intent(in) :: lj_types !DOC LJTypes

  integer :: i
  integer :: nk ! number of grid points in kspace

  logical :: overlap

  mc_system_composition => comp   

! Total Sum (alloc, calculate)
  call EwaldSumRealSpace_alloc_full(RealSpace_TotalSum, atomic_data, lj_types )
  !call EwaldSumRealSpace_calc(RealSpace_TotalSum,overlap)   ! init Fx,Fy,Fz, uu


  call EwaldSumKSpace_alloc(KSpace_TotalSum,atomic_data,lj_types)
  !all EwaldSumKSpace_init(KSpace_TotalSum) ! init grid, beta, rho_squared_total (i.e. sumsincos... )

  !call EwaldSumExternal_init(Ext_TotalSum, atomic_data % xx, atomic_data % yy, atomic_data % zz, &
   !                          atomic_data % charge, atomic_data % Natom)

  !call EwaldSumExternal_calc_mu(Ext_TotalSum)

  !U_tail_coulomb = ewald_sum_coulomb_tail(atomic_data % charge, atomic_data % natom )
 
  !call  ewald_sum_lj_tails( comp, lj_types, U_tail_lj6, U_tail_lj12)

  !U_tail = U_tail_coulomb + U_tail_lj12 - U_tail_lj6
  !U_ext = EwaldSumExternal_calc_energy( Ext_TotalSum )

  !U_total_ew = RealSpace_TotalSum % uu_ew + KSpace_TotalSum % energy_coulomb + U_tail_coulomb + U_ext
  !U_total_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 + &
  !             KSpace_TotalSum % energy_lj + U_tail_lj12 - U_tail_lj6

  !U_total = U_total_ew + U_total_lj 

 ! Allocate Real Partial Sums (old,new)

  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_old, atomic_data, lj_types)
  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_new, atomic_data, lj_types)

  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_old1, atomic_data, lj_types)
  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_new1, atomic_data, lj_types)

  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_old2, atomic_data, lj_types)
  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_new2, atomic_data, lj_types)

  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_old12, atomic_data, lj_types)
  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_new12, atomic_data, lj_types)





! Allocate sincos_kr_coulomb, sincos_kr_LJ (used for the force calculations)
  MaxMolAtom = mol_table % MaxMolAtom

  allocate( sincos_kr_coulomb(MaxMolAtom * 2) )  ! *2 for exchange step
  allocate( sincos_kr_LJ(MaxMolAtom * 2) )

  nk = KSpace_TotalSum % grid % nk 
  
  do i=1,MaxMolAtom
     call SinCosKR_alloc( sincos_kr_coulomb(i) , KSpace_TotalSum % grid  ) 
     call SinCosKR_alloc( sincos_kr_LJ(i) , KSpace_TotalSum % grid )
  end do 

! Allocate x,y,znew,fx,fy,fz_old,new

  allocate(xnew( MaxMolAtom  ) )  
  allocate(ynew( MaxMolAtom  ) )  
  allocate(znew( MaxMolAtom  ) )


  allocate(xnew1( MaxMolAtom  ) )  
  allocate(ynew1( MaxMolAtom  ) )  
  allocate(znew1( MaxMolAtom  ) )

  allocate(xnew2( MaxMolAtom  ) )  
  allocate(ynew2( MaxMolAtom  ) )  
  allocate(znew2( MaxMolAtom  ) )

 
  allocate(xtmp1( MaxMolAtom  ) )  
  allocate(ytmp1( MaxMolAtom  ) )  
  allocate(ztmp1( MaxMolAtom  ) )

  allocate(xtmp2( MaxMolAtom  ) )  
  allocate(ytmp2( MaxMolAtom  ) )  
  allocate(ztmp2( MaxMolAtom  ) )


  allocate(fx_old( MaxMolAtom ) )
  allocate(fy_old( MaxMolAtom ) )  
  allocate(fz_old( MaxMolAtom ) )

  allocate(fx_ext_old( MaxMolAtom ) )
  allocate(fy_ext_old( MaxMolAtom ) )  
  allocate(fz_ext_old( MaxMolAtom ) )

  allocate(fx_ext_new( MaxMolAtom ) )
  allocate(fy_ext_new( MaxMolAtom ) )  
  allocate(fz_ext_new( MaxMolAtom ) )


  allocate(fx_new( MaxMolAtom ) )
  allocate(fy_new( MaxMolAtom ) )  
  allocate(fz_new( MaxMolAtom ) )

  call RhoSquared_alloc( sumsincos_mol_old, KSpace_TotalSum % grid, lj_types )
  call RhoSquared_alloc( sumsincos_mol_new, KSpace_TotalSum % grid, lj_types )
  call RhoSquared_alloc( sumsincos_total_new, KSpace_TotalSum % grid, lj_types )

  call MCState_alloc( mc_stored_state, atomic_data % natom, KSpace_TotalSum % grid, lj_types ) 

end subroutine





!subroutine mc_calc_total_sums_dbg(mol_table,atomic_data,lj_types)
!
!!  use io, only : write_real_array
!
!  implicit none
!!  Type(TComposition),intent(in) :: mc_system_composition
!  Type(TMoleculeTable),intent(in) :: mol_table
!  Type(TAtomicData),intent(in) :: atomic_data
!  Type(TLJTypes),intent(in) :: lj_types
!
!  integer :: i
!  integer :: nk ! number of grid points in kspace
!
!  logical :: overlap
!
!! Total Sum (alloc, calculate)
!!->!!  call EwaldSumRealSpace_calc(RealSpace_TotalSum,overlap)   ! init Fx,Fy,Fz, uu
!
!!  call EwaldSumKSpace_init(KSpace_TotalSum) ! init grid, beta, rho_squared_total (i.e. sumsincos... )
!     call EwaldSumKSpace_initBeta(KSpace_TotalSum) ! beta,beta6,beta12
!!     call EwaldSumKSpace_initBetaLJ(KSpace_TotalSum)
!!i
!! initialize coulomb and LJ energies ( and sin/cos sums )
!! 
!!      call EwaldSumKSpace_calc_total_energy( KSpace_TotalSum )
!
!
!!  call EwaldSumExternal_init(Ext_TotalSum, atomic_data % xx, atomic_data % yy, atomic_data % zz, &
!!                             atomic_data % charge, atomic_data % Natom)
!
!!  call EwaldSumExternal_calc_mu(Ext_TotalSum)
!
!!  U_tail_coulomb = ewald_sum_coulomb_tail(atomic_data % charge, atomic_data % natom )
! 
!!  call  ewald_sum_lj_tails( mc_system_composition, lj_types, U_tail_lj6, U_tail_lj12)
!
!!  U_tail = U_tail_coulomb + U_tail_lj12 - U_tail_lj6
!
!
!!  U_total = mc_total_energy()
!
!end subroutine
!

subroutine mc_calc_total_sums_dbg(mol_table,atomic_data,lj_types,scale_beta,overlap) !DOC
!DOC DEBUG ONLY

!  use io, only : write_real_array

  implicit none
!  Type(TComposition),intent(in) :: mc_system_composition
  Type(TMoleculeTable),intent(in) :: mol_table
  Type(TAtomicData),intent(in) :: atomic_data
  Type(TLJTypes),intent(in) :: lj_types
  real(8),intent(in),optional :: scale_beta
  

  integer :: i
  integer :: nk ! number of grid points in kspace

  logical,intent(out),optional :: overlap
  logical :: overlap_local

! Total Sum (alloc, calculate)
  call EwaldSumRealSpace_calc(RealSpace_TotalSum,overlap_local)   ! init Fx,Fy,Fz, uu

  if( present(overlap) ) overlap = overlap_local

  return

  if( .not.present(scale_beta) ) then
     call EwaldSumKSpace_init(KSpace_TotalSum) ! init grid, beta, rho_squared_total (i.e. sumsincos... )
  else
     call EwaldSumKSpace_init(KSpace_TotalSum,scale_beta)
  end if 

  call EwaldSumExternal_init(Ext_TotalSum, atomic_data % xx, atomic_data % yy, atomic_data % zz, &
                             atomic_data % charge, atomic_data % Natom)

  call EwaldSumExternal_calc_mu(Ext_TotalSum)

  U_tail_coulomb = ewald_sum_coulomb_tail(atomic_data % charge, atomic_data % natom )
 
  call  ewald_sum_lj_tails( mc_system_composition, lj_types, U_tail_lj6, U_tail_lj12)

  U_tail = U_tail_coulomb + U_tail_lj12 - U_tail_lj6

  U_total = mc_total_energy()

end subroutine


subroutine mc_calc_total_sums(mol_table,atomic_data,lj_types,scale_beta,overlap) !DOC
!DOC calculate the KSpace and RealSpace sums
!DOC Parameters:
!  use io, only : write_real_array

  implicit none
!  Type(TComposition),intent(in) :: mc_system_composition
  Type(TMoleculeTable),intent(in) :: mol_table !DOC Molecule Table (used to know how the atoms are grouped to the molecules)
  Type(TAtomicData),intent(in) :: atomic_data !DOC AtomicData : atomic coordinates
  Type(TLJTypes),intent(in) :: lj_types !DOC LJTypes for each pair of atom types
  real(8),intent(in),optional :: scale_beta !DOC if the volume is changed - we need to scale the beta function. Otherwise, it can be re-calculated

  integer :: i
  integer :: nk ! number of grid points in kspace

  logical,intent(out),optional :: overlap
  logical :: overlap_local

! Total Sum (alloc, calculate)
  !call EwaldSumRealSpace_alloc_full(RealSpace_TotalSum, atomic_data, lj_types )
  call EwaldSumRealSpace_calc(RealSpace_TotalSum,overlap_local)   ! init Fx,Fy,Fz, uu

  if( present(overlap) ) overlap = overlap_local

  !call EwaldSumKSpace_alloc(KSpace_TotalSum,atomic_data,lj_types)

  if( .not.present(scale_beta) ) then
     call EwaldSumKSpace_init(KSpace_TotalSum) ! init grid, beta, rho_squared_total (i.e. sumsincos... )
  else
     call EwaldSumKSpace_init(KSpace_TotalSum,scale_beta)
  end if 

  call EwaldSumExternal_init(Ext_TotalSum, atomic_data % xx, atomic_data % yy, atomic_data % zz, &
                             atomic_data % charge, atomic_data % Natom)

  call EwaldSumExternal_calc_mu(Ext_TotalSum)

  U_tail_coulomb = ewald_sum_coulomb_tail(atomic_data % charge, atomic_data % natom )
 
  call  ewald_sum_lj_tails( mc_system_composition, lj_types, U_tail_lj6, U_tail_lj12)

  U_tail = U_tail_coulomb + U_tail_lj12 - U_tail_lj6


  !U_ext = EwaldSumExternal_calc_energy( Ext_TotalSum )

  !U_total_ew = RealSpace_TotalSum % uu_ew + KSpace_TotalSum % energy_coulomb + U_tail_coulomb + U_ext
  !U_total_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 + &
  !             KSpace_TotalSum % energy_lj + U_tail_lj12 - U_tail_lj6

  !U_total = U_total_ew + U_total_lj 
  U_total = mc_total_energy()
!  U_ext = EwaldSumExternal_calc_energy( Ext_TotalSum )

!  U_total_ew = RealSpace_TotalSum % uu_ew + KSpace_TotalSum % energy_coulomb + U_tail_coulomb + U_ext
!  U_total_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 + &
!               KSpace_TotalSum % energy_lj + U_tail_lj12 - U_tail_lj6

!  U_total = U_total_ew + U_total_lj 
 
 ! Allocate Real Partial Sums (old,new)

!  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_old, atomic_data, lj_types)
!  call EwaldSumRealSpace_alloc_partial( RealSpace_MolSum_new, atomic_data, lj_types)

! Allocate sincos_kr_coulomb, sincos_kr_LJ (used for the force calculations)
!  MaxMolAtom = mol_table % MaxMolAtom

!  allocate( sincos_kr_coulomb(MaxMolAtom) )
!  allocate( sincos_kr_LJ(MaxMolAtom) )

!  nk = KSpace_TotalSum % grid % nk 
  
!  do i=1,MaxMolAtom
!     call SinCosKR_alloc( sincos_kr_coulomb(i) , KSpace_TotalSum % grid  ) 
!     call SinCosKR_alloc( sincos_kr_LJ(i) , KSpace_TotalSum % grid )
!  end do 

! Allocate x,y,znew,fx,fy,fz_old,new

!  allocate(xnew( MaxMolAtom ) )
!  allocate(ynew( MaxMolAtom ) )  
!  allocate(znew( MaxMolAtom ) )

!  allocate(fx_old( MaxMolAtom ) )
!  allocate(fy_old( MaxMolAtom ) )  
!  allocate(fz_old( MaxMolAtom ) )
!
!  allocate(fx_ext_old( MaxMolAtom ) )
!  allocate(fy_ext_old( MaxMolAtom ) )  
!  allocate(fz_ext_old( MaxMolAtom ) )
!
!  allocate(fx_ext_new( MaxMolAtom ) )
!  allocate(fy_ext_new( MaxMolAtom ) )  
!  allocate(fz_ext_new( MaxMolAtom ) )
!
!
!  allocate(fx_new( MaxMolAtom ) )
!  allocate(fy_new( MaxMolAtom ) )  
!  allocate(fz_new( MaxMolAtom ) )
!
!  call RhoSquared_alloc( sumsincos_mol_old, KSpace_TotalSum % grid, lj_types )
!  call RhoSquared_alloc( sumsincos_mol_new, KSpace_TotalSum % grid, lj_types )
!  call RhoSquared_alloc( sumsincos_total_new, KSpace_TotalSum % grid, lj_types )
!
!  call MCState_alloc( mc_stored_state, atomic_data % natom, KSpace_TotalSum % grid, lj_types ) 


end subroutine

subroutine initmc(comp,mol_table,atomic_data,lj_types) !DOC
!DOC initialize the data
!DOC Parameters:
!  use io, only : write_real_array
  implicit none
  Type(TComposition),intent(in),target :: comp !DOC composition
  Type(TMoleculeTable),intent(in) :: mol_table !DOC molecule table (how the atoms are grouped to the molecules)
  Type(TAtomicData),intent(in) :: atomic_data !DOC coordinates of atoms
  Type(TLJTypes),intent(in) :: lj_types !DOC LJ parameters for each pair of atom types

  call mc_alloc(comp,mol_table,atomic_data,lj_types)
  call mc_calc_total_sums(mol_table,atomic_data,lj_types)

  mc_composition => comp

end subroutine

subroutine clearmc

  implicit none
  integer :: i,nk

  deallocate(xnew)
  deallocate(ynew)
  deallocate(znew)

  deallocate(xnew1)
  deallocate(ynew1)
  deallocate(znew1)

  deallocate(xnew2)
  deallocate(ynew2)
  deallocate(znew2)

  deallocate(xtmp1)
  deallocate(ytmp1)
  deallocate(ztmp1)

  deallocate(xtmp2)
  deallocate(ytmp2)
  deallocate(ztmp2)


  deallocate(fx_old )
  deallocate(fy_old )  
  deallocate(fz_old )

  deallocate(fx_ext_old )
  deallocate(fy_ext_old )  
  deallocate(fz_ext_old )

  deallocate(fx_ext_new )
  deallocate(fy_ext_new )  
  deallocate(fz_ext_new )


  deallocate(fx_new )
  deallocate(fy_new )  
  deallocate(fz_new )

  nk = KSpace_TotalSum % grid % nk 
  
  do i=1,MaxMolAtom
     call SinCosKR_dealloc( sincos_kr_coulomb(i) ) 
     call SinCosKR_dealloc( sincos_kr_LJ(i) )
  end do 

  deallocate( sincos_kr_coulomb )
  deallocate( sincos_kr_LJ )

  call EwaldSumRealSpace_dealloc_full(RealSpace_TotalSum)
  call EwaldSumKSpace_dealloc(KSpace_TotalSum)  

  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_old)
  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_new)

  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_old1)
  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_new1)

  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_old2)
  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_new2)

  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_old12)
  call EwaldSumRealSpace_dealloc_partial(RealSpace_MolSum_new12)


  call RhoSquared_dealloc( sumsincos_mol_old )
  call RhoSquared_dealloc( sumsincos_mol_new )
  call RhoSquared_dealloc( sumsincos_total_new )

  call MCState_dealloc( mc_stored_state )

end subroutine

function mc_total_energy() !DOC
!DOC calculate the total energy (using the already calculated ewald sums).
!DOC Return value:
!DOC  total energy of the system
   real(8) :: mc_total_energy

  U_tail = U_tail_coulomb + U_tail_lj12 - U_tail_lj6
  U_ext = EwaldSumExternal_calc_energy( Ext_TotalSum )

  U_total_ew = RealSpace_TotalSum % uu_ew + KSpace_TotalSum % energy_coulomb + U_tail_coulomb + U_ext
  U_total_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 + &
               KSpace_TotalSum % energy_lj + U_tail_lj12 - U_tail_lj6

  U_total = U_total_ew + U_total_lj 
  
  mc_total_energy = U_total

end function

 
subroutine mcvol( mol_table, accepted , force_new_box_length, force_accepted) !DOC
!DOC Volume step: try to change the volume of the box
!DOC Parameters:
   use parameters, only : log_max_volume_scaling,BoxLength,pressure_angstr_kT,parameters_recalc,xclb
   use MRandom, only : rand
   use ScaleBox, only : scale_box

   Type(TMoleculeTable),intent(in) :: mol_table !DOC MoleculeTable
   Type(TAtomicData),pointer :: atomic_data !DOC coordinates of atoms
   Type(TLJTypes),pointer :: lj_types !DOC LJ parameters
   logical :: accepted
   real(8),intent(in),optional :: force_new_box_length  !DOC force setting the new box length. For Debug Only! (usually omitted)
   logical,intent(in),optional :: force_accepted !DOC  set accepted to force_accepted  debug only!  (usually omitted)


   real(8) :: xclb_old,xclb_new

   real(8) :: old_volume,new_volume,log_new_volume,OldBoxLength,NewBoxLength
   integer :: nmol
  
   real(8) :: U_tot_old,U_tot_new,dU,dV,U_tot_old1
   real(8) ::  arg,proba
   real(8) :: lambda

   logical :: overlap
   ! random walk in lnV
   ! see D.Frenkel & B. Smit, "Understanding Molecular Simulaions", sec. 5.4.2, Algorithm 11.
  
   atomic_data => RealSpace_TotalSum % atomic_data
   lj_types => RealSpace_TotalSum % lj_types

!   call mc_calc_total_sums(mol_table,atomic_data,lj_types)


   OldBoxLength = BoxLength  
   old_volume = BoxLength**3

   lambda = (2*rand() - 1.d0)
   log_new_volume = log(old_volume) + lambda * log_max_volume_scaling
   ! new_volume = volume * max_volume_scaling ** lambda

   new_volume = exp(log_new_volume)
   NewBoxLength = new_volume**(1.d0/3.d0)

   if( present(force_new_box_length) ) then
       NewBoxLength = force_new_box_length
   end if 

   U_tot_old = mc_total_energy()
 

!   call mc_calc_total_sums(mol_table,atomic_data,lj_types)

!   U_tot_old1 = mc_total_energy()

 !  write(*,*) 'U_tot_old:',U_tot_old,U_tot_old1


   ! store the old configuration
   call mc_save_state( mc_stored_state )

   ! update the box length
   call scale_box(  atomic_data, mol_table, NewBoxLength )

   xclb_old = xclb
   call parameters_recalc(NewBoxLength)
   xclb_new = xclb

   call LJTypes_scale_sigma( lj_types, OldBoxLength / NewBoxLength ) 
   call EwaldSumKSpace_initBetaLJ( KSpace_TotalSum )
   
   ! calculate the new sums 
   call mc_calc_total_sums( mol_table, atomic_data, lj_types, xclb_new / xclb_old, overlap)

   if( .not. overlap ) then
   
      U_tot_new = mc_total_energy()  
   
      ! acceptance rate: 
      ! acc(old->new) = min ( 1, exp(-beta dU - beta P dV  + (N+1) ln V_new/V_old ) )
      ! dU = U_new - U_old,  dV = V_new - V_old 
         
      ! beta = 1/kT. In our case the energy is in kT already
      dU = U_tot_new - U_tot_old 
      dV = new_volume - old_volume 
   
      nmol = mol_table % nmol
      arg = -dU - pressure_angstr_kT * dV + (nmol+1) * log( new_volume / old_volume )
      proba = exp(arg)

   else
      proba = -1.d0
   end if
 
!   write(*,*) 'MC_VOL V_old',old_volume,'V_new',new_volume,'U_tot_old',U_tot_old,'U_tot_new',U_tot_new
  

   if( present( force_accepted ) ) then
      if( force_accepted ) then
          proba = 1. 
      end if
   end if
   
   if( ( proba > 1 ) .or. ( rand() < proba ) ) then 
      accepted = .TRUE.  
   else
      call mc_restore_state( mc_stored_state )
      atomic_data % BoxLength = OldBoxLength
      call parameters_recalc( OldBoxLength ) 
      call LJTypes_scale_sigma( RealSpace_TotalSum % lj_types, NewBoxLength / OldBoxLength )
      call EwaldSumKSpace_initBetaLJ( KSpace_TotalSum ) 
      KSpace_TotalSum % beta(:) = KSpace_TotalSum % beta(:) * xclb_old / xclb_new
      U_tot_new = mc_total_energy()
      accepted = .FALSE.
   end if 
 
!   call initmc
!   call initmc(mc_composition,mol_table,atomic_data,lj_types)



!    call mc_calc_total_sums(mol_table,atomic_data,lj_types)


end subroutine 





subroutine movemc( imol, mol_table,  accepted ) !DOC 
!DOC Monte Carlo move step
!DOC Parameters:
  use MRandom
  use MoleculeTable
  use MoleculeHandler
  use geometry
!  use io,only : write_real_array,write_xyz_array,io_open,io_close
  use SystemSettings, only : SYSTEM_STRING_LENGTH

  implicit none
  integer :: imol  !DOC number of molecule to be tried to move (should be chosen randomly)
  Type(TMoleculeTable) :: mol_table !DOC Molecule Table
  logical :: accepted !DOC output: whether or not the move was accepted
  Type(TMonteCarloMove) :: mc_move_data

  real(8) :: fx_kspace_coulomb, fy_kspace_coulomb,fz_kspace_coulomb  ! used ad output of ForceKSpace_cuoulomb
  real(8) :: fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ   ! outputs for ForceKSpace_LJ

  integer ::  first_atom, last_atom
  integer :: imol_type

  Type(TRotation) :: rot
  Type(TRotMatrix) :: rot_matrix

  integer :: i,ii
  integer :: nk
 
!  integer :: hfile
!  character(SYSTEM_STRING_LENGTH) :: filename

 ! POINTERS
  real(8),dimension(:),pointer :: beta
  Type(TRealArrayPointer),dimension(:),pointer :: beta_LJ

  real(8),dimension(:),pointer :: xx,yy,zz,charge  ! to RealSpace_TotalSum % xx,yy,zz
  real(8),dimension(:),pointer :: xx_old,yy_old,zz_old ! to xx,yy,zz(first_atom:last_atom)
  Type(TMolecule),pointer :: mol_ptr

  real(8) :: DeltaU_fourier,DeltaU_real,dU_ext,du_total
  real(8) :: du_k_coulomb,du_k_LJ
  real(8) :: deltx,delty,deltz

  real(8) :: RealSpace_MolEnergy_new,RealSpace_MolEnergy_old

  real(8) :: dmu_x,dmu_y,dmu_z
  real(8) :: mu_x_new,mu_y_new,mu_z_new

  real(8) :: torq_x,torq_y,torq_z

  real(8) :: fx_total_old,fy_total_old,fz_total_old
  real(8) :: torq_x_old,torq_y_old,torq_z_old
  real(8) :: fx_total_new,fy_total_new,fz_total_new
  real(8) :: torq_x_new,torq_y_new,torq_z_new

  integer :: ityp
  integer,dimension(:),pointer :: type_by_index
  integer :: natom,ntotal
  integer :: ntype

  logical :: overlap

  Type(TSumSinCosKR),pointer :: sumsincos_coulomb_old,sumsincos_coulomb_new
  Type(TSumSinCosKR),dimension(:),pointer :: sumsincos_LJ_old , sumsincos_LJ_new

  real(8) :: xr,yr,zr
  real(8) :: x_center, y_center, z_center
  real(8) :: x_center_new,y_center_new,z_center_new

  character(SYSTEM_STRING_LENGTH) :: filename

  nk = KSpace_TotalSum % grid % nk  

  ntype = KSpace_TotalSum % lj_types % NType
  type_by_index => KSpace_TotalSum % lj_types % type_by_index

  beta => KSpace_TotalSum % beta
  beta_LJ => KSpace_TotalSum % beta_LJ

  xx => RealSpace_TotalSum % atomic_data % xx
  yy => RealSpace_TotalSum % atomic_data % yy
  zz => RealSpace_TotalSum % atomic_data % zz
  charge => RealSpace_TotalSum % atomic_data % charge

! FOR TEST ONLY
      call EwaldSumExternal_calc_mu(Ext_TotalSum)

! Step 1:  Select the molecule $M$.
!  imol = random(mol_table % nmol )
  first_atom = mol_table % first_atom ( imol ) 
  last_atom = mol_table % last_atom ( imol )  
  natom = last_atom - first_atom + 1 
  ntotal = RealSpace_TotalSum % atomic_data % natom

  xx_old => xx(first_atom:last_atom)
  yy_old => yy(first_atom:last_atom)
  zz_old => zz(first_atom:last_atom)

  imol_type = mol_table % mol_type(imol)
  call MoleculeHandler_getMolecule(imol_type,mol_ptr)
  
! Step 2: [U_{real}]_M^{old} = \sum_{j\in M} u_j

  RealSpace_MolEnergy_old = SUM( RealSpace_TotalSum % uu(first_atom:last_atom) )

  ! copy, not pointer, because after that we add F_kspace_old to fx,y,z_old
  fx_old(1:natom) = RealSpace_TotalSum % Fx( first_atom:last_atom )
  fy_old(1:natom) = RealSpace_TotalSum % Fy( first_atom:last_atom )
  fz_old(1:natom) = RealSpace_TotalSum % Fz( first_atom:last_atom )

! Step 3: Calculate [S_m]_M, [C_m]_M
!     [S_m^t]_M, [C_m^t]_M for $t \in T(M)
!
!  AND k-space forces
!
   call RhoSquared_nulify( sumsincos_mol_old )

   sumsincos_coulomb_old => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
   sumsincos_LJ_old => KSpace_TotalSum % rho_squared_total % sumsincos_LJ

   ii= first_atom   ! absolute atom index
   do i=1,natom ! local atom index

      ityp = type_by_index(ii)

!      write(*,*) 'addAtom: i',i,'ityp',ityp,'ii',ii,'imol',imol
      call RhoSquared_addAtom( sumsincos_mol_old, xx(ii), yy(ii), zz(ii), charge(ii), type_by_index(ii), &
                               sincos_kr_coulomb(i), sincos_kr_LJ(i) )
      ! ok, now  sincos_kr_coulomb(i) sincos_kr_LJ(i) are filled. In a sence, we don't need here the array
      ! i.e. size(sincos_kr_coulomb) can be 1
      ! but we anyway need these arrays in new! 
 
      call ForceKSpace_coulomb( beta , sumsincos_coulomb_old,  sincos_kr_coulomb(i),&
                                fx_kspace_coulomb, fy_kspace_coulomb, fz_kspace_coulomb ) 
    
     if ( ityp > 0) then
         call ForceKSpace_LJ(beta_LJ,sumsincos_LJ_old, ntype, ityp,  sincos_kr_LJ(i),&
                                fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ )

      else
          fx_kspace_LJ = 0
          fy_kspace_LJ = 0
          fz_kspace_LJ = 0
      end if 


      fx_old(i) = fx_old(i) + fx_kspace_coulomb + fx_kspace_LJ
      fy_old(i) = fy_old(i) + fy_kspace_coulomb + fy_kspace_LJ
      fz_old(i) = fz_old(i) + fz_kspace_coulomb + fz_kspace_LJ

      ii = ii + 1
   end do  

! Step 3A: Calculate ExternalMedia forces

!   U_ext_old = EwaldSumExternal_calc_energy( Ext_TotalSum ) nobody needs it

   call external_sum_calc_forces( Ext_TotalSum % mu_x, Ext_TotalSum % mu_y, &
                                  Ext_TotalSum % mu_z, charge(first_atom:last_atom), & 
                                  natom, fx_ext_old, fy_ext_old, fz_ext_old )

   fx_old(1:natom) = fx_old(1:natom) + fx_ext_old(1:natom)
   fy_old(1:natom) = fy_old(1:natom) + fy_ext_old(1:natom)
   fz_old(1:natom) = fz_old(1:natom) + fz_ext_old(1:natom)

   call EwaldSumExternal_init( Ext_MolSum_old, xx_old, yy_old, zz_old, charge(first_atom:last_atom), natom )
   call EwaldSumExternal_calc_mu( Ext_MolSum_old )


!Step 4:  Calculate total force and total torque of the molecule $M$:

      fx_total_old = SUM(fx_old(1:natom))
      fy_total_old = SUM(fy_old(1:natom))
      fz_total_old = SUM(fz_old(1:natom))


      call center_of_mass( xx_old, yy_old, zz_old, mol_ptr % mass, natom , x_center, y_center, z_center  )

      torq_x_old = 0
      torq_y_old = 0
      torq_z_old = 0

      do i=1,natom
         
         xr = xx_old(i) - x_center
         yr = yy_old(i) - y_center
         zr = zz_old(i) - z_center         

         call prod_vect(xr,yr,zr, fx_old(i) ,fy_old(i), fz_old(i), torq_x, torq_y, torq_z)
          
         torq_x_old = torq_x_old + torq_x
         torq_y_old = torq_y_old + torq_y
         torq_z_old = torq_z_old + torq_z 

      end do

       

!Step 5: Select the new position of the molecule

   call MonteCarloMove_init_old( mc_move_data, fx_total_old,fy_total_old,fz_total_old, & 
                                               torq_x_old,  torq_y_old, torq_z_old );

   call MonteCarloMove_chooseDisplacement( mc_move_data, deltx,delty,deltz)

   !write(*,*) 'deltxyz: ',deltx,delty,deltz

!   if ( natom > 1 ) then 
   call MonteCarloMove_chooseRotation( mc_move_data, rot ) ! I do this, because otherwise accept/decline conditions possibly will differ

   x_center_new = x_center + deltx - ANINT(x_center + deltx)
   y_center_new = y_center + delty - ANINT(y_center + delty)
   z_center_new = z_center + deltz - ANINT(z_center + deltz)
 
   if ( natom > 1 ) then

      call fill_rot_matrix(rot,rot_matrix)

      do i=1,natom

        xnew(i) = xx_old(i) - x_center
        ynew(i) = yy_old(i) - y_center
        znew(i) = zz_old(i) - z_center

        call rotate_vect(rot_matrix, xnew(i),ynew(i), znew(i),xr,yr,zr )

        xnew(i) = xr + x_center_new
        ynew(i) = yr + y_center_new
        znew(i) = zr + z_center_new

      end do
   
   else
        xnew(1) = x_center_new
        ynew(1) = y_center_new
        znew(1) = z_center_new
   end if

!Step 6: item For $p\in M$, $j \notin M$ calculate
!          u_p^new, (uu)
!          u_j]_M^{new} ( uu_back) 
!           [U_{real}]_M^{new}  =\sum_{p \in M} u_p^{new}
!          F_p^new(real)
!          [\mathbf{F}_j]_M^{new} ( f_p_back) (new) 
!

  call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_new, xnew, ynew, znew, first_atom, last_atom )
  call EwaldSumRealSpace_calc( RealSpace_MolSum_new, overlap)      

  if ( overlap ) then 
     accepted = .FALSE.
!     write(*,*) 'overlap'
     return
  end if 
  
  RealSpace_MolEnergy_new = SUM( RealSpace_MolSum_new % uu(1: RealSpace_MolSum_new % natom) )  
 
!Step 7: \Delta U_{real} = [U_{real}]_M^{new} - [U_{real}]_M^{old}

  DeltaU_real = RealSpace_MolEnergy_new - RealSpace_MolEnergy_old

!Step 8:  Calculate [S_{\mathbf{m}}]_M^{new}, [C_{\mathbf{m}}]_M^{new}$ and
!        [S_{\mathbf{m}}^t]_M^{new}, [C_{\mathbf{m}}^t]_M^{new} for $t \in T(M)
!     

   call RhoSquared_nulify( sumsincos_mol_new )

   ii= first_atom  ! absolute atom index
   do i=1,natom ! local atom index

      ityp = type_by_index(ii)

      call RhoSquared_addAtom( sumsincos_mol_new, xnew(i), ynew(i), znew(i), charge(ii), type_by_index(ii), & 
                               sincos_kr_coulomb(i), sincos_kr_LJ(i)  ) ! delta_sumsincos = sumsincos_mol_new
      ! save sincos_kr for each atom to use them in the force calculations
      ! ( there we need both: new C_m,S_m ( known only after this cycle) and sincos_kr ( calculated, but lost in RhoSquared)
      ! thus, we need to save sincos_kr somewhere... 
!  (OBSOLETE) --->  call SinCosKR_copy( sincos_kr_coulomb (i) , sumsincos_mol_new % sincos_kr_coulomb ) ! dst = src
!  (OBSOLETE) --->  call SinCosKR_copy( sincos_kr_LJ(i) , sumsincos_mol_new % sincos_kr_LJ )

      ii = ii + 1  
    end do      

!Step 9:

   !filename = 'sumsincos_mol_old.rho2'
   !call RhoSquared_save_to_file( sumsincos_mol_old,filename )
   !filename = 'sumsincos_mol_new.rho2'
   !call RhoSquared_save_to_file( sumsincos_mol_new, filename)
   !filename = 'sumsincos_total_old.rho2'
   !call RhoSquared_save_to_file( KSpace_TotalSum % rho_squared_total,filename)
  

   delta_sumsincos => sumsincos_mol_new 
   call RhoSquared_sub( delta_sumsincos, sumsincos_mol_old ) ! delta_sumsincos = sumsincos_mol_new - sumsincos_mol_old

   !filename = 'delta_sumsincos.rho2'
   !call RhoSquared_save_to_file( delta_sumsincos,filename)

!   sumsincos_coulomb_new => delta_sumsincos %  sumsincos_LJ(1)
!
!   sumsincos_coulomb_new => sumsincos_mol_old %  sumsincos_LJ(1)

! Step 9A : Calculate new Fourier Forces

   call RhoSquared_copy( sumsincos_total_new , KSpace_TotalSum % rho_squared_total)
   call RhoSquared_add( sumsincos_total_new , delta_sumsincos )

   sumsincos_coulomb_new => sumsincos_total_new % sumsincos_coulomb
   sumsincos_LJ_new => sumsincos_total_new % sumsincos_LJ

   ii = first_atom
   do i=1,natom

      ityp = type_by_index(ii)
      call ForceKSpace_coulomb( beta , sumsincos_coulomb_new, sincos_kr_coulomb(i) ,&
                                fx_kspace_coulomb, fy_kspace_coulomb, fz_kspace_coulomb ) 

      if( ityp > 0) then
         call ForceKSpace_LJ(beta_LJ,sumsincos_LJ_new, ntype, ityp, sincos_kr_LJ(i)  ,&
                                fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ ) 
      else
          fx_kspace_LJ = 0
          fy_kspace_LJ = 0
          fz_kspace_LJ = 0
      end if 

      fx_new(i) = RealSpace_MolSum_new % Fx(i) + fx_kspace_coulomb + fx_kspace_LJ
      fy_new(i) = RealSpace_MolSum_new % Fy(i) + fy_kspace_coulomb + fy_kspace_LJ
      fz_new(i) = RealSpace_MolSum_new % Fz(i) + fz_kspace_coulomb + fz_kspace_LJ

      ii = ii + 1
   end do  

!Step 10:  Calculate \Delta U_{fourier}
!Step 11: $\Delta U_{fourier} = \Delta U_{fourier}^C + \Delta U_{fourier}^{LJ}$

  du_k_coulomb = EwaldSumKSpace_calc_dU_coulomb( KSpace_TotalSum, delta_sumsincos % sumsincos_coulomb ) 
  du_k_LJ = EwaldSumKSpace_calc_dU_LJ(KSpace_TotalSum,delta_sumsincos % sumsincos_LJ, delta_sumsincos % type_is_present )  

  !write(*,*) 'du_k_LJ=',du_k_LJ

  DeltaU_fourier = du_k_coulomb + du_k_LJ


! Setp 11 A: new energy and forces for external field

  call EwaldSumExternal_init(Ext_MolSum_new, xnew, ynew, znew, charge(first_atom:last_atom), natom )
  call EwaldSumExternal_calc_mu(Ext_MolSum_new )

  dmu_x =  Ext_MolSum_new % mu_x - Ext_MolSum_old % mu_x 
  dmu_y =  Ext_MolSum_new % mu_y - Ext_MolSum_old % mu_y 
  dmu_z =  Ext_MolSum_new % mu_z - Ext_MolSum_old % mu_z 

  dU_ext = EwaldSumExternal_calc_dU( Ext_TotalSum, dmu_x, dmu_y, dmu_z )

! Step 11B: New forces

  mu_x_new = Ext_TotalSum % mu_x + dmu_x
  mu_y_new = Ext_TotalSum % mu_y + dmu_y 
  mu_z_new = Ext_TotalSum % mu_z + dmu_z
 
  call external_sum_calc_forces( mu_x_new, mu_y_new, mu_z_new, charge(first_atom:last_atom), natom, &
                                 fx_ext_new, fy_ext_new, fz_ext_new )

  fx_new(1:natom) = fx_new(1:natom) + fx_ext_new(1:natom)
  fy_new(1:natom) = fy_new(1:natom) + fy_ext_new(1:natom)
  fz_new(1:natom) = fz_new(1:natom) + fz_ext_new(1:natom)

!Step 12:  Calculate new forces and torques for molecule $M$: 

   fx_total_new = SUM(fx_new(1:natom))
   fy_total_new = SUM(fy_new(1:natom))
   fz_total_new = SUM(fz_new(1:natom))

   call center_of_mass( xnew, ynew, znew, mol_ptr % mass, natom , x_center, y_center, z_center  )

   torq_x_new = 0
   torq_y_new = 0
   torq_z_new = 0

   do i=1,natom
      
      xr = xnew(i) - x_center
      yr = ynew(i) - y_center
      zr = znew(i) - z_center         

      call prod_vect(xr,yr,zr, fx_new(i) ,fy_new(i), fz_new(i), torq_x, torq_y, torq_z)
       
      torq_x_new = torq_x_new + torq_x
      torq_y_new = torq_y_new + torq_y
      torq_z_new = torq_z_new + torq_z 

   end do
 
!Step 13: Accept or decline the move  

    du_total = DeltaU_real + DeltaU_fourier + dU_ext

    call MonteCarloMove_init_new( mc_move_data, fx_total_new, fy_total_new, fz_total_new, &
                                           torq_x_new,  torq_y_new, torq_z_new )

    call MonteCarloMove_accept_or_decline( mc_move_data , du_total, accepted )


!    write(*,*) 'MY: du:',du_total,'F:', fx_total_new,fy_total_new,fz_total_new,&
!                                  'T:',torq_x_new,torq_y_new,torq_z_new
   

!Step 14: if accepted then: 

   if ( accepted) then

      ! Update k-space sums:
      call RhoSquared_add( KSpace_TotalSum % rho_squared_total, delta_sumsincos  )

      !filename = 'sumsincos_total_new.rho2'
      !call RhoSquared_save_to_file( KSpace_TotalSum % rho_squared_total,filename)

      ! calculate old back forces and atom energies 
      call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_old, xx_old, yy_old, zz_old, first_atom, last_atom )
      call EwaldSumRealSpace_calc( RealSpace_MolSum_old, overlap )      

      ! Update the total real energy and forces   
      ! (back forces 

      do i=1,ntotal
 
!         if ( imol == type_by_index(i) ) cycle ! not account the atoms of the selected molecule
!  WRONG! type_by_index is type of the molecule, not ID! 

         RealSpace_TotalSum % uu(i) = RealSpace_TotalSum % uu(i) + &
                                      RealSpace_MolSum_new % uu_back(i) - RealSpace_MolSum_old % uu_back(i)

         RealSpace_TotalSum % Fx(i) = RealSpace_TotalSum % Fx(i) + &
                                      RealSpace_MolSum_new % Fx_back(i) - RealSpace_MolSum_old % Fx_back(i) 

         RealSpace_TotalSum % Fy(i) = RealSpace_TotalSum % Fy(i) + &
                                      RealSpace_MolSum_new % Fy_back(i) - RealSpace_MolSum_old % Fy_back(i)


         RealSpace_TotalSum % Fz(i) = RealSpace_TotalSum % Fz(i) +&
                                      RealSpace_MolSum_new % Fz_back(i) - RealSpace_MolSum_old % Fz_back(i)
  
      end do

          

      RealSpace_TotalSum % uu(first_atom:last_atom) = RealSpace_MolSum_new % uu(1:natom) 

      RealSpace_TotalSum % Fx(first_atom:last_atom) = RealSpace_MolSum_new % Fx(1:natom)
      RealSpace_TotalSum % Fy(first_atom:last_atom) = RealSpace_MolSum_new % Fy(1:natom)
      RealSpace_TotalSum % Fz(first_atom:last_atom) = RealSpace_MolSum_new % Fz(1:natom)

      ! Update k-space energy:

    !  write(*,*) 'energy_coulomb:',KSpace_TotalSum % energy_coulomb,'+',du_k_coulomb
    !  write(*,*) 'energy_lj:',KSpace_TotalSum % energy_lj,'+',du_k_lj

      KSpace_TotalSum % energy_coulomb = KSpace_TotalSum % energy_coulomb + du_k_coulomb
      KSpace_TotalSum % energy_LJ = KSpace_TotalSum % energy_lj + du_k_lj
      KSpace_TotalSum % energy = KSpace_TotalSum % energy_coulomb + KSpace_TotalSum % energy_lj


      ! RealSpace energies

    !  write(*,*) 'uu_ew:',RealSpace_TotalSum % uu_ew,'+',RealSpace_MolSum_new % uu_ew - RealSpace_MolSum_old % uu_ew 
    !  write(*,*) 'uu_lj6:',RealSpace_TotalSum % uu_lj6,'+',RealSpace_MolSum_new % uu_lj6 - RealSpace_MolSum_old % uu_lj6 
    !  write(*,*) 'uu_lj12:',RealSpace_TotalSum % uu_lj12,'+',RealSpace_MolSum_new % uu_lj12 - RealSpace_MolSum_old % uu_lj12 

 
      RealSpace_TotalSum % uu_ew = RealSpace_TotalSum % uu_ew + RealSpace_MolSum_new % uu_ew - RealSpace_MolSum_old % uu_ew 
      RealSpace_TotalSum % uu_lj6 = RealSpace_TotalSum % uu_lj6 + RealSpace_MolSum_new % uu_lj6 - RealSpace_MolSum_old % uu_lj6 
      RealSpace_TotalSum % uu_lj12 = RealSpace_TotalSum % uu_lj12 + RealSpace_MolSum_new % uu_lj12 - RealSpace_MolSum_old % uu_lj12 

     ! U_Ext

      Ext_TotalSum % mu_x = mu_x_new
      Ext_TotalSum % mu_y = mu_y_new
      Ext_TotalSum % mu_z = mu_z_new

!      call EwaldSumExternal_calc_mu(Ext_TotalSum)

      U_total = mc_total_energy()

      xx(first_atom:last_atom) = xnew(1:natom)
      yy(first_atom:last_atom) = ynew(1:natom)
      zz(first_atom:last_atom) = znew(1:natom)

   end if ! accepted

end subroutine 

 
subroutine mc_xchange( imol1, imol2, mol_table,  accepted ) !DOC
!DOC Monte Carlo exchange step
  use MRandom
  use MoleculeTable
  use MoleculeHandler
  use geometry, only : center_of_mass
  use io,only : write_real_array,write_xyz_array,io_open,io_close
  use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC Parameters:
  implicit none

  integer,intent(in) :: imol1,imol2 !DOC numers of the molecules to exchange
  Type(TMoleculeTable) :: mol_table !DOC MoleculeTable
  logical :: accepted !DOC output: whether or not the exchange was accepted
  Type(TMonteCarloMove) :: mc_move_data

  real(8) :: fx_kspace_coulomb, fy_kspace_coulomb,fz_kspace_coulomb  ! used ad output of ForceKSpace_cuoulomb
  real(8) :: fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ   ! outputs for ForceKSpace_LJ

  integer :: first_atom1, last_atom1, first_atom2,last_atom2
  integer :: imol_type1,imol_type2

!  Type(TRotation) :: rot
!  Type(TRotMatrix) :: rot_matrix

  integer :: i,ii
  integer :: nk
 
  integer :: hfile
  character(SYSTEM_STRING_LENGTH) :: filename

 ! POINTERS
  real(8),dimension(:),pointer :: beta
  Type(TRealArrayPointer),dimension(:),pointer :: beta_LJ

  real(8),dimension(:),pointer :: xx,yy,zz,charge  ! to RealSpace_TotalSum % xx,yy,zz
  real(8),dimension(:),pointer :: xx_old1,yy_old1,zz_old1 ! to xx,yy,zz(first_atom:last_atom)
  real(8),dimension(:),pointer :: xx_old2,yy_old2,zz_old2 ! to xx,yy,zz(first_atom:last_atom)
 
  Type(TMolecule),pointer :: mol_ptr1,mol_ptr2

  real(8) :: DeltaU_fourier,DeltaU_real,dU_ext,du_total
  real(8) :: du_k_coulomb,du_k_LJ
 ! real(8) :: deltx,delty,deltz
  real(8) :: DeltaU12

  real(8) :: RealSpace_MolEnergy_new,  RealSpace_MolEnergy_old

  real(8) :: dmu1_x,dmu1_y,dmu1_z
  real(8) :: dmu2_x,dmu2_y,dmu2_z
  real(8) :: dmu_x,dmu_y,dmu_z ! dmu = dmu1 + dmu2
  real(8) :: mu_x_new,mu_y_new,mu_z_new

  real(8) :: fx_norecalc(1000),fy_norecalc(1000),fz_norecalc(1000),uu_norecalc(1000)
  real(8) :: SFx,SFy,SFz,Suu  

!  real(8) :: torq_x,torq_y,torq_z

!  real(8) :: fx_total_old,fy_total_old,fz_total_old
!  real(8) :: torq_x_old,torq_y_old,torq_z_old
!  real(8) :: fx_total_new,fy_total_new,fz_total_new
!  real(8) :: torq_x_new,torq_y_new,torq_z_new

  integer :: ityp1, ityp2
  integer,dimension(:),pointer :: type_by_index
  integer :: natom1,natom2,ntotal
  integer :: ntype

  logical :: overlap1,overlap2,overlap12

  Type(TSumSinCosKR),pointer :: sumsincos_coulomb_old, sumsincos_coulomb_new
  Type(TSumSinCosKR),dimension(:),pointer :: sumsincos_LJ_old , sumsincos_LJ_new

  real(8) :: xr,yr,zr
  real(8) :: x_center1, y_center1, z_center1
  real(8) :: x_center2, y_center2, z_center2

  real(8) :: U_tot_1,U_tot_0
!  character(SYSTEM_STRING_LENGTH) :: filename

  nk = KSpace_TotalSum % grid % nk  

  ntype = KSpace_TotalSum % lj_types % NType
  type_by_index => KSpace_TotalSum % lj_types % type_by_index

  beta => KSpace_TotalSum % beta
  beta_LJ => KSpace_TotalSum % beta_LJ

  xx => RealSpace_TotalSum % atomic_data % xx
  yy => RealSpace_TotalSum % atomic_data % yy
  zz => RealSpace_TotalSum % atomic_data % zz
  charge => RealSpace_TotalSum % atomic_data % charge

! FOR TEST ONLY
      call EwaldSumExternal_calc_mu(Ext_TotalSum)

! Step 1:  Select the molecule $M$.
!  imol = random(mol_table % nmol )
  first_atom1 = mol_table % first_atom ( imol1 ) 
  last_atom1 = mol_table % last_atom ( imol1 )  
  natom1 = last_atom1 - first_atom1 + 1

  first_atom2 = mol_table % first_atom ( imol2 ) 
  last_atom2 = mol_table % last_atom ( imol2 )  
  natom2 = last_atom2 - first_atom2 + 1

!   write(*,*) 'first_atom1',first_atom1,'last_atom1',last_atom1,'natom1',natom1
!    write(*,*) 'first_atom2',first_atom2,'last_atom2',last_atom2,'natom2',natom2
  

  ntotal = RealSpace_TotalSum % atomic_data % natom

  xx_old1 => xx(first_atom1:last_atom1)
  yy_old1 => yy(first_atom1:last_atom1)
  zz_old1 => zz(first_atom1:last_atom1)

  xx_old2 => xx(first_atom2:last_atom2)
  yy_old2 => yy(first_atom2:last_atom2)
  zz_old2 => zz(first_atom2:last_atom2)

  imol_type1 = mol_table % mol_type(imol1)
  call MoleculeHandler_getMolecule(imol_type1,mol_ptr1)
  
  imol_type2 = mol_table % mol_type(imol2)
  call MoleculeHandler_getMolecule(imol_type2,mol_ptr2)

! Step 2: [U_{real}]_M^{old} = \sum_{j\in M} u_j

  RealSpace_MolEnergy_old = SUM( RealSpace_TotalSum % uu(first_atom1:last_atom1) )
  RealSpace_MolEnergy_old = RealSpace_MolEnergy_old + SUM( RealSpace_TotalSum % uu(first_atom2:last_atom2) ) 

 ! We don't need forces 
  ! copy, not pointer, because after that we add F_kspace_old to fx,y,z_old
!  fx_old(1:natom) = RealSpace_TotalSum % Fx( first_atom:last_atom )
!  fy_old(1:natom) = RealSpace_TotalSum % Fy( first_atom:last_atom )
!  fz_old(1:natom) = RealSpace_TotalSum % Fz( first_atom:last_atom )

! Step 3: Calculate [S_m]_M, [C_m]_M
!     [S_m^t]_M, [C_m^t]_M for $t \in T(M)
!
!  AND k-space forces
!
   call RhoSquared_nulify( sumsincos_mol_old )

   sumsincos_coulomb_old => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
   sumsincos_LJ_old => KSpace_TotalSum % rho_squared_total % sumsincos_LJ


! Note : we don't need any forces ==> everything is ok with sincos_kr_coulomb(1)
   ii= first_atom1   ! absolute atom index
   do i=1,natom1 ! local atom index

      call RhoSquared_addAtom( sumsincos_mol_old, xx(ii), yy(ii), zz(ii), charge(ii), type_by_index(ii), &
                               sincos_kr_coulomb(1), sincos_kr_LJ(1) )
      ii = ii + 1
   end do  

   ii= first_atom2   ! absolute atom index
   do i=1,natom2 ! local atom index

      call RhoSquared_addAtom( sumsincos_mol_old, xx(ii), yy(ii), zz(ii), charge(ii), type_by_index(ii), &
                               sincos_kr_coulomb(1), sincos_kr_LJ(1) )
      ii = ii + 1
   end do  



! Step 3A: Calculate ExternalMedia forces

   call EwaldSumExternal_init( Ext_MolSum_old1, xx_old1, yy_old1, zz_old1, charge(first_atom1:last_atom1), natom1 )
   call EwaldSumExternal_calc_mu( Ext_MolSum_old1 )

   call EwaldSumExternal_init( Ext_MolSum_old2, xx_old2, yy_old2, zz_old2, charge(first_atom2:last_atom2), natom2 )
   call EwaldSumExternal_calc_mu( Ext_MolSum_old2 )

!Step 4:  Calculate total force and total torque of the molecule $M$:
! Don't need this 
      

!Step 5: Select the new position of the molecule

!   call MonteCarloMove_init_old( mc_move_data, fx_total_old,fy_total_old,fz_total_old, & 
!                                               torq_x_old,  torq_y_old, torq_z_old );

!   call MonteCarloMove_chooseDisplacement( mc_move_data, deltx,delty,deltz)

   !write(*,*) 'deltxyz: ',deltx,delty,deltz

!   if ( natom > 1 ) then 
!   call MonteCarloMove_chooseRotation( mc_move_data, rot ) ! I do this, because otherwise accept/decline conditions possibly will differ

   call center_of_mass( xx_old1, yy_old1, zz_old1, mol_ptr1 % mass, natom1 , x_center1, y_center1, z_center1  )
   call center_of_mass( xx_old2, yy_old2, zz_old2, mol_ptr2 % mass, natom2 , x_center2, y_center2, z_center2  )

!   x_center1 = x_center2
!   y_center1= y_center2
!   z_center1 = z_center2

   
   ! xchanging the centers of mass of the molecules
   xnew1(1:natom1) = xx_old1(1:natom1) - x_center1 + x_center2
   ynew1(1:natom1) = yy_old1(1:natom1) - y_center1 + y_center2
   znew1(1:natom1) = zz_old1(1:natom1) - z_center1 + z_center2
  
   xnew2(1:natom2) = xx_old2(1:natom2) - x_center2 + x_center1
   ynew2(1:natom2) = yy_old2(1:natom2) - y_center2 + y_center1
   znew2(1:natom2) = zz_old2(1:natom2) - z_center2 + z_center1

!Step 6: item For $p\in M$, $j \notin M$ calculate
!          u_p^new, (uu)
!          u_j]_M^{new} ( uu_back) 
!           [U_{real}]_M^{new}  =\sum_{p \in M} u_p^{new}
!          F_p^new(real)
!          [\mathbf{F}_j]_M^{new} ( f_p_back) (new) 
!


! to compute the new energy properly, we need to have the new coordinates in the system.
! (this was not the case for 1 molecule, because we did not consider intermoleculear interactions)
! now our "molecule" consists of two, and the interactions between these two molecules are essential

! so, we store the "old" coordinates
  xtmp1(1:natom1) = xx(first_atom1:last_atom1)
  ytmp1(1:natom1) = yy(first_atom1:last_atom1)
  ztmp1(1:natom1) = zz(first_atom1:last_atom1)

  xtmp2(1:natom2) = xx(first_atom2:last_atom2)
  ytmp2(1:natom2) = yy(first_atom2:last_atom2)
  ztmp2(1:natom2) = zz(first_atom2:last_atom2)
   
! change them to "new" ones
  
  xx(first_atom1:last_atom1) = xnew1(1:natom1) 
  yy(first_atom1:last_atom1) = ynew1(1:natom1)
  zz(first_atom1:last_atom1) = znew1(1:natom1)

  xx(first_atom2:last_atom2) = xnew2(1:natom2)
  yy(first_atom2:last_atom2) = ynew2(1:natom2)
  zz(first_atom2:last_atom2) = znew2(1:natom2)

! now we can compute the energies of the molecules

  call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_new1, xnew1, ynew1, znew1, first_atom1, last_atom1 )
  call EwaldSumRealSpace_calc( RealSpace_MolSum_new1, overlap1)      

  call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_new2, xnew2, ynew2, znew2, first_atom2, last_atom2 )
  call EwaldSumRealSpace_calc( RealSpace_MolSum_new2, overlap2)      

  call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_new12, xnew1, ynew1, znew1, first_atom1, last_atom1 )
  call EwaldSumRealSpace_calc( RealSpace_MolSum_new12, overlap12, first_atom2,last_atom2)      


! and restore the initial coordinates 
 
  xx(first_atom1:last_atom1) = xtmp1(1:natom1) 
  yy(first_atom1:last_atom1) = ytmp1(1:natom1) 
  zz(first_atom1:last_atom1) = ztmp1(1:natom1) 
                              
  xx(first_atom2:last_atom2) = xtmp2(1:natom2) 
  yy(first_atom2:last_atom2) = ytmp2(1:natom2) 
  zz(first_atom2:last_atom2) = ztmp2(1:natom2) 


  if ( overlap1 .OR. overlap2 ) then 
     accepted = .FALSE.
!     write(*,*) 'OVERLAP',imol1,imol2
!     write(*,*)  'new1',xnew1(1),ynew1(1),znew1(1)
!     write(*,*)  'new2',xnew2(1),ynew2(1),znew2(1)
!     write(*,*)  'old1',xx_old1(1),yy_old1(1),zz_old1(1)
!     write(*,*)   'old2',xx_old2(1),yy_old2(1),zz_old2(1)
      write(*,'(A,$)') 'OX'
 !    stop
     return
  end if 

  
  RealSpace_MolEnergy_new = SUM( RealSpace_MolSum_new1 % uu(1:natom1) )  
  RealSpace_MolEnergy_new = RealSpace_MolEnergy_new + SUM( RealSpace_MolSum_new2 % uu(1:natom2) )
 
!Step 7: \Delta U_{real} = [U_{real}]_M^{new} - [U_{real}]_M^{old}

  DeltaU_real = RealSpace_MolEnergy_new - RealSpace_MolEnergy_old

!Step 8:  Calculate [S_{\mathbf{m}}]_M^{new}, [C_{\mathbf{m}}]_M^{new}$ and
!        [S_{\mathbf{m}}^t]_M^{new}, [C_{\mathbf{m}}^t]_M^{new} for $t \in T(M)
!     

   call RhoSquared_nulify( sumsincos_mol_new )

   ii= first_atom1  ! absolute atom index
   do i=1,natom1 ! local atom index

      call RhoSquared_addAtom( sumsincos_mol_new, xnew1(i), ynew1(i), znew1(i), charge(ii), type_by_index(ii), & 
                               sincos_kr_coulomb(1), sincos_kr_LJ(1)  ) ! delta_sumsincos = sumsincos_mol_new
      ii = ii + 1  
    end do      

   ii= first_atom2  ! absolute atom index
   do i=1,natom2 ! local atom index

      call RhoSquared_addAtom( sumsincos_mol_new, xnew2(i), ynew2(i), znew2(i), charge(ii), type_by_index(ii), & 
                               sincos_kr_coulomb(1), sincos_kr_LJ(1)  ) ! delta_sumsincos = sumsincos_mol_new
      ii = ii + 1  
    end do      


!Step 9:

   !filename = 'sumsincos_mol_old.rho2'
   !call RhoSquared_save_to_file( sumsincos_mol_old,filename )
   !filename = 'sumsincos_mol_new.rho2'
   !call RhoSquared_save_to_file( sumsincos_mol_new, filename)
   !filename = 'sumsincos_total_old.rho2'
   !call RhoSquared_save_to_file( KSpace_TotalSum % rho_squared_total,filename)
  

   delta_sumsincos => sumsincos_mol_new 
   call RhoSquared_sub( delta_sumsincos, sumsincos_mol_old ) ! delta_sumsincos = sumsincos_mol_new - sumsincos_mol_old

   !filename = 'delta_sumsincos.rho2'
   !call RhoSquared_save_to_file( delta_sumsincos,filename)

!   sumsincos_coulomb_new => delta_sumsincos %  sumsincos_LJ(1)
!
!   sumsincos_coulomb_new => sumsincos_mol_old %  sumsincos_LJ(1)

! Step 9A : Calculate new Fourier Forces

   call RhoSquared_copy( sumsincos_total_new , KSpace_TotalSum % rho_squared_total)
   call RhoSquared_add( sumsincos_total_new , delta_sumsincos )

   sumsincos_coulomb_new => sumsincos_total_new % sumsincos_coulomb
   sumsincos_LJ_new => sumsincos_total_new % sumsincos_LJ

!   ii = first_atom
!   do i=1,natom
!
!      ityp = type_by_index(ii)
!      call ForceKSpace_coulomb( beta , sumsincos_coulomb_new, sincos_kr_coulomb(i) ,&
!                                fx_kspace_coulomb, fy_kspace_coulomb, fz_kspace_coulomb ) 
!
!      if( ityp > 0) then
!         call ForceKSpace_LJ(beta_LJ,sumsincos_LJ_new, ntype, ityp, sincos_kr_LJ(i)  ,&
!                                fx_kspace_LJ, fy_kspace_LJ, fz_kspace_LJ ) 
!      else
!          fx_kspace_LJ = 0
!          fy_kspace_LJ = 0
!          fz_kspace_LJ = 0
!      end if 
!
!      fx_new(i) = RealSpace_MolSum_new % Fx(i) + fx_kspace_coulomb + fx_kspace_LJ
!      fy_new(i) = RealSpace_MolSum_new % Fy(i) + fy_kspace_coulomb + fy_kspace_LJ
!      fz_new(i) = RealSpace_MolSum_new % Fz(i) + fz_kspace_coulomb + fz_kspace_LJ
!
!      ii = ii + 1
!   end do  

!Step 10:  Calculate \Delta U_{fourier}
!Step 11: $\Delta U_{fourier} = \Delta U_{fourier}^C + \Delta U_{fourier}^{LJ}$

  du_k_coulomb = EwaldSumKSpace_calc_dU_coulomb( KSpace_TotalSum, delta_sumsincos % sumsincos_coulomb ) 
  du_k_LJ = EwaldSumKSpace_calc_dU_LJ(KSpace_TotalSum,delta_sumsincos % sumsincos_LJ, delta_sumsincos % type_is_present )  

  !write(*,*) 'du_k_LJ=',du_k_LJ

  DeltaU_fourier = du_k_coulomb + du_k_LJ


! Setp 11 A: new energy and forces for external field

  call EwaldSumExternal_init(Ext_MolSum_new1, xnew1, ynew1, znew1, charge(first_atom1:last_atom1), natom1 )
  call EwaldSumExternal_calc_mu(Ext_MolSum_new1 )

  dmu1_x =  Ext_MolSum_new1 % mu_x - Ext_MolSum_old1 % mu_x 
  dmu1_y =  Ext_MolSum_new1 % mu_y - Ext_MolSum_old1 % mu_y 
  dmu1_z =  Ext_MolSum_new1 % mu_z - Ext_MolSum_old1 % mu_z 

  call EwaldSumExternal_init(Ext_MolSum_new2, xnew2, ynew2, znew2, charge(first_atom2:last_atom2), natom2 )
  call EwaldSumExternal_calc_mu(Ext_MolSum_new2 )

  dmu2_x =  Ext_MolSum_new2 % mu_x - Ext_MolSum_old2 % mu_x 
  dmu2_y =  Ext_MolSum_new2 % mu_y - Ext_MolSum_old2 % mu_y 
  dmu2_z =  Ext_MolSum_new2 % mu_z - Ext_MolSum_old2 % mu_z 

  dmu_x = dmu1_x + dmu2_x
  dmu_y = dmu1_y + dmu2_y
  dmu_z = dmu1_z + dmu2_z

  dU_ext = EwaldSumExternal_calc_dU( Ext_TotalSum, dmu_x, dmu_y, dmu_z )
         

! Step 11B: New forces

  mu_x_new = Ext_TotalSum % mu_x + dmu_x
  mu_y_new = Ext_TotalSum % mu_y + dmu_y 
  mu_z_new = Ext_TotalSum % mu_z + dmu_z
 
!  call external_sum_calc_forces( mu_x_new, mu_y_new, mu_z_new, charge(first_atom:last_atom), natom, &
!                                 fx_ext_new, fy_ext_new, fz_ext_new )

!  fx_new(1:natom) = fx_new(1:natom) + fx_ext_new(1:natom)
!  fy_new(1:natom) = fy_new(1:natom) + fy_ext_new(1:natom)
!  fz_new(1:natom) = fz_new(1:natom) + fz_ext_new(1:natom)

!Step 12:  Calculate new forces and torques for molecule $M$: 

! Skip this 


!Step 13: Accept or decline the move  

    du_total = DeltaU_real + DeltaU_fourier + dU_ext

!    call MonteCarloMove_init_new( mc_move_data, fx_total_new, fy_total_new, fz_total_new, &
!                                           torq_x_new,  torq_y_new, torq_z_new )

!    call MonteCarloMove_accept_or_decline( mc_move_data , du_total, accepted )

    call MonteCarloMove_accept_or_decline_Simple( du_total, accepted ) 


!    write(*,*) 'MY: du:',du_total,'F:', fx_total_new,fy_total_new,fz_total_new,&
!                                  'T:',torq_x_new,torq_y_new,torq_z_new
   

!Step 14: if accepted then: 

   if ( accepted) then

      ! Update k-space sums:
      call RhoSquared_add( KSpace_TotalSum % rho_squared_total, delta_sumsincos  )

      !filename = 'sumsincos_total_new.rho2'
      !call RhoSquared_save_to_file( KSpace_TotalSum % rho_squared_total,filename)

      ! calculate old back forces and atom energies 
      call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_old1, xx_old1, yy_old1, zz_old1, first_atom1, last_atom1 )
      call EwaldSumRealSpace_calc( RealSpace_MolSum_old1, overlap1 )      

      call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_old2, xx_old2, yy_old2, zz_old2, first_atom2, last_atom2 )
      call EwaldSumRealSpace_calc( RealSpace_MolSum_old2, overlap2 )      

      call EwaldSumRealSpace_set_molecule( RealSpace_MolSum_old12, xx_old1, yy_old1, zz_old1, first_atom1, last_atom1 )
      call EwaldSumRealSpace_calc( RealSpace_MolSum_old12, overlap12, first_atom2, last_atom2)      


      ! Update the total real energy and forces   
      ! (back forces 

      do i=1,ntotal
 
!         if ( imol1 == type_by_index(i) ) cycle ! not account the atoms of the selected molecule

         RealSpace_TotalSum % uu(i) = RealSpace_TotalSum % uu(i) + &
                                      RealSpace_MolSum_new1 % uu_back(i) - RealSpace_MolSum_old1 % uu_back(i) + &
                                      RealSpace_MolSum_new2 % uu_back(i) - RealSpace_MolSum_old2 % uu_back(i)

         RealSpace_TotalSum % Fx(i) = RealSpace_TotalSum % Fx(i) + &
                                      RealSpace_MolSum_new1 % Fx_back(i) - RealSpace_MolSum_old1 % Fx_back(i) + & 
                                      RealSpace_MolSum_new2 % Fx_back(i) - RealSpace_MolSum_old2 % Fx_back(i) 

         RealSpace_TotalSum % Fy(i) = RealSpace_TotalSum % Fy(i) + &
                                      RealSpace_MolSum_new1 % Fy_back(i) - RealSpace_MolSum_old1 % Fy_back(i) + &
                                      RealSpace_MolSum_new2 % Fy_back(i) - RealSpace_MolSum_old2 % Fy_back(i)


         RealSpace_TotalSum % Fz(i) = RealSpace_TotalSum % Fz(i) +&
                                      RealSpace_MolSum_new1 % Fz_back(i) - RealSpace_MolSum_old1 % Fz_back(i) + &
                                      RealSpace_MolSum_new2 % Fz_back(i) - RealSpace_MolSum_old2 % Fz_back(i)

      end do

          

      RealSpace_TotalSum % uu(first_atom1:last_atom1) = RealSpace_MolSum_new1 % uu(1:natom1) 

      RealSpace_TotalSum % Fx(first_atom1:last_atom1) = RealSpace_MolSum_new1 % Fx(1:natom1)
      RealSpace_TotalSum % Fy(first_atom1:last_atom1) = RealSpace_MolSum_new1 % Fy(1:natom1)
      RealSpace_TotalSum % Fz(first_atom1:last_atom1) = RealSpace_MolSum_new1 % Fz(1:natom1)


      RealSpace_TotalSum % uu(first_atom2:last_atom2) = RealSpace_MolSum_new2 % uu(1:natom2) 

      RealSpace_TotalSum % Fx(first_atom2:last_atom2) = RealSpace_MolSum_new2 % Fx(1:natom2)
      RealSpace_TotalSum % Fy(first_atom2:last_atom2) = RealSpace_MolSum_new2 % Fy(1:natom2)
      RealSpace_TotalSum % Fz(first_atom2:last_atom2) = RealSpace_MolSum_new2 % Fz(1:natom2)


      ! Update k-space energy:

    !  write(*,*) 'energy_coulomb:',KSpace_TotalSum % energy_coulomb,'+',du_k_coulomb
    !  write(*,*) 'energy_lj:',KSpace_TotalSum % energy_lj,'+',du_k_lj

      KSpace_TotalSum % energy_coulomb = KSpace_TotalSum % energy_coulomb + du_k_coulomb
      KSpace_TotalSum % energy_LJ = KSpace_TotalSum % energy_lj + du_k_lj
      KSpace_TotalSum % energy = KSpace_TotalSum % energy_coulomb + KSpace_TotalSum % energy_lj


      ! RealSpace energies

    !  write(*,*) 'uu_ew:',RealSpace_TotalSum % uu_ew,'+',RealSpace_MolSum_new % uu_ew - RealSpace_MolSum_old % uu_ew 
    !  write(*,*) 'uu_lj6:',RealSpace_TotalSum % uu_lj6,'+',RealSpace_MolSum_new % uu_lj6 - RealSpace_MolSum_old % uu_lj6 
    !  write(*,*) 'uu_lj12:',RealSpace_TotalSum % uu_lj12,'+',RealSpace_MolSum_new % uu_lj12 - RealSpace_MolSum_old % uu_lj12 


    !  we need to account for the 12 energy, because now it is calculated twice
    !  indeed: we have two selected molecules, let they are 1 and 2
    !    u1^old = u1^old(others) + u12^old
    !    u2^old = u2^old(others) + u12^old
    !    u1^new = u1^new(others) + u12^new
    !    u2^new = u2^new(others) + u12^new
    ! where ui^(others) = SUM_j=3..N uij
    !  then
    !   u1^new + u2^new - u1^old - u2^old = 
    !   u1^new(others) + u2^new(others) - u1^old(others) - u2^old(others) + 2u12^new - 2u12^old 
    ! so, we count the DeltaU12 twice...
    ! we need to subtract it
    !
    ! in principle, U12 can be obtained from the uu_back
    !    indeed, it is   SUM(new1 % uu_back(new2) ) - SUM(new2 % uu_back(new1) )
    !  However, if we want to count lj and ew separately, we need something more intelegent
    ! that's why we are calculating the molecule-molecule (12) interaction separately (see above RealSpace_MolSum_new12, RealSpace_MolSum_old12 ) 
    
     

      RealSpace_TotalSum % uu_ew = RealSpace_TotalSum % uu_ew + &
                                   RealSpace_MolSum_new1 % uu_ew - RealSpace_MolSum_old1 % uu_ew + & 
                                   RealSpace_MolSum_new2 % uu_ew - RealSpace_MolSum_old2 % uu_ew - &
                                   ( RealSpace_MolSum_new12 % uu_ew - RealSpace_MolSum_old12 % uu_ew ) 


      RealSpace_TotalSum % uu_lj6 = RealSpace_TotalSum % uu_lj6 + &
                                    RealSpace_MolSum_new1 % uu_lj6 - RealSpace_MolSum_old1 % uu_lj6 + &
                                    RealSpace_MolSum_new2 % uu_lj6 - RealSpace_MolSum_old2 % uu_lj6 - &
                                    ( RealSpace_MolSum_new12 % uu_lj6 - RealSpace_MolSum_old12 % uu_lj6 )


      RealSpace_TotalSum % uu_lj12 = RealSpace_TotalSum % uu_lj12 + & 
                                     RealSpace_MolSum_new1 % uu_lj12 - RealSpace_MolSum_old1 % uu_lj12 + &
                                     RealSpace_MolSum_new2 % uu_lj12 - RealSpace_MolSum_old2 % uu_lj12 - &
                                     (RealSpace_MolSum_new12 % uu_lj12 - RealSpace_MolSum_old12 % uu_lj12 )

     ! U_Ext

      Ext_TotalSum % mu_x = mu_x_new
      Ext_TotalSum % mu_y = mu_y_new
      Ext_TotalSum % mu_z = mu_z_new

!      call EwaldSumExternal_calc_mu(Ext_TotalSum)

      U_total = mc_total_energy()
!      U_tot_0 = U_total  ! be careful : U_total is global!

      xx(first_atom1:last_atom1) = xnew1(1:natom1)
      yy(first_atom1:last_atom1) = ynew1(1:natom1)
      zz(first_atom1:last_atom1) = znew1(1:natom1)

      xx(first_atom2:last_atom2) = xnew2(1:natom2)
      yy(first_atom2:last_atom2) = ynew2(1:natom2)
      zz(first_atom2:last_atom2) = znew2(1:natom2)

!     filename = 'forces_norecalc.txt'
!     hfile = io_open(filename,'w')
!     call write_xyz_array(hfile,RealSpace_TotalSum % Fx, RealSpace_TotalSum % Fy, RealSpace_TotalSum % Fz, &
!                                RealSpace_TotalSum % atomic_data % natom ) 
!     call io_close(hfile)
!
!     filename = 'uu_norecalc.txt'
!      hfile = io_open(filename,'w')
!     call write_real_array(hfile,RealSpace_TotalSum % uu,RealSpace_TotalSum % atomic_data % natom ) 
!     call io_close(hfile)
!
!     fx_norecalc(1:ntotal) = RealSpace_TotalSum % Fx(1:ntotal)
!     fy_norecalc(1:ntotal) = RealSpace_TotalSum % Fy(1:ntotal)
!     fz_norecalc(1:ntotal) = RealSpace_TotalSum % Fz(1:ntotal)
!
!     uu_norecalc(1:ntotal) = RealSpace_TotalSum % uu(1:ntotal)
!

!      call mc_calc_total_sums_dbg(mol_table,RealSpace_TotalSum % atomic_data,RealSpace_TotalSum % lj_types)
!     call EwaldSumRealSpace_calc_dbg(RealSpace_TotalSum,overlap1)   ! init Fx,Fy,Fz, uu

     
     

      
 !    U_tot_1 = mc_total_energy()

 !    filename = 'forces_recalc.txt'
 !    hfile = io_open(filename,'w')
 !    call write_xyz_array(hfile,RealSpace_TotalSum % Fx, RealSpace_TotalSum % Fy, RealSpace_TotalSum % Fz, ntotal )
 !    call io_close(hfile)

  !   filename = 'uu_recalc.txt'
  !    hfile = io_open(filename,'w')
  !   call write_real_array(hfile,RealSpace_TotalSum % uu,ntotal ) 
  !   call io_close(hfile)
    
   !   write(*,*) 'xchg diff:',U_tot_0 - U_tot_1
   !   write(*,*) 'oldsum12',RealSpace_MolSum_old12 % uu_ew, RealSpace_MolSum_old12 % uu_lj6, RealSpace_MolSum_old12 % uu_lj12
   !   write(*,*) 'newsum12',RealSpace_MolSum_new12 % uu_ew, RealSpace_MolSum_new12 % uu_lj6, RealSpace_MolSum_new12 % uu_lj12
   

!     SFx=0d0
!      SFy = 0d0
!      SFz = 0d0
!      Suu = 0d0
!      do i=1,ntotal
!
!         SFx = SFx + dabs( RealSpace_TotalSum % Fx(i) - fx_norecalc(i) )
!         SFy = SFy + dabs( RealSpace_TotalSum % Fy(i) - fy_norecalc(i) )
!         SFz = SFz + dabs( RealSpace_TotalSum % Fz(i) - fz_norecalc(i) )
!         Suu = Suu + dabs( RealSpace_TotalSum % uu(i) - uu_norecalc(i) )
!
!      end do 
!
!      if( (SFx + SFy + SFz + Suu + U_total - U_tot_1 )  > 1e-9) stop      
!!      write(*,*)  (SFx + SFy + SFz + Suu + U_total - U_tot_1 ) 
!     write(*,*) 'new1 % uu_back(new2)',SUM(RealSpace_MolSum_new1 % uu_back(first_atom2:last_atom2))
!     write(*,*) 'new2 % uu_back(new1)',SUM(RealSpace_MolSum_new2 % uu_back(first_atom1:last_atom1))
!     write(*,*) 'old1 % uu_back(old2)',SUM(RealSpace_MolSum_old1 % uu_back(first_atom2:last_atom2))
!     write(*,*) 'old2 % uu_back(old1)',SUM(RealSpace_MolSum_old2 % uu_back(first_atom1:last_atom1))
!
!
!           
!
!      stop
   end if ! accepted

end subroutine 


subroutine check_distances(xx,yy,zz,xxn,yyn,zzn,natom) !DOC
!DOC DEBUG ONLY: checks that the distances between the atoms in two arrays are the same
!DOC Parameters:
  use error
  real(8),dimension(:),intent(in) :: xx,yy,zz,xxn,yyn,zzn !DOC coordinates of the atoms
  integer,intent(in) :: natom !DOC number of atoms

  integer :: i,j
  real(8) :: r2_old,r2_new

  do i=1,natom
    do j=i,natom

        r2_old = (xx(i) - xx(j))**2 + (yy(i) - yy(j))**2 + (zz(i) - zz(j))**2
        r2_new = (xxn(i) - xxn(j))**2 + (yyn(i) - yyn(j))**2 + (zzn(i)-zzn(j))**2

        if( abs(r2_old - r2_new) > 1d-6) then

            call error_throw(ERROR_PARAMETER)

        end if

    end do
  end do


end subroutine


subroutine check_consistency(mol_tab) !DOC
!DOC DEBUG ONLY: cheks whether the distances between the atoms in the molecules stored in mol_tab % atomic_data are the same as the distances in the input files
!DOC Parameters:
   use parameters
!   use MC, only : RealSpace_TotalSum,check_distances
   use Molecule
   use MoleculeHandler   
   use error
   use io, only : write_real_array,write_xyz_array,io_open,io_close
   use SystemSettings, only : SYSTEM_STRING_LENGTH

   Type(TMoleculeTable),intent(in) :: mol_tab !DOC MoleculeTable
   Type(TMolecule),pointer :: mol_ptr 

   integer :: i,nmol,first_atom,last_atom
   real(8),dimension(:),pointer :: xx,yy,zz
   real(8),dimension(10) :: xx1,yy1,zz1

   integer :: hfile
   character(SYSTEM_STRING_LENGTH) :: filename

   xx => RealSpace_TotalSum % atomic_data % xx
   yy => RealSpace_TotalSum % atomic_data % yy
   zz => RealSpace_TotalSum % atomic_data % zz

   call error_set_catch(ERROR_PARAMETER)

   nmol = mol_tab % nmol
   do i=1,nmol

       call MoleculeHandler_getMolecule( mol_tab % mol_type(i), mol_ptr )

       first_atom = mol_tab % first_atom(i)
       last_atom = mol_tab % last_atom(i)

       xx1 = xx(first_atom:last_atom) * BoxLength
       yy1 = yy(first_atom:last_atom) * BoxLength
       zz1 = zz(first_atom:last_atom) * BoxLength 


       call check_distances( mol_ptr % x, mol_ptr % y, mol_ptr % z, xx1, yy1, zz1, mol_ptr % NAtoms )

       if (error_code /= 0) then

           write(*,*) 'Consistency check failed at molecule #',i
           call write_xyz_array(0,xx1,yy1,zz1,mol_ptr % NAtoms )

           filename = 'xyz_cons.coors'
           hfile = io_open(filename,'w')
           call write_xyz_array(hfile,xx,yy,zz,mol_tab % last_atom(nmol))
           call io_close(hfile)

           call error_clear_catch(ERROR_PARAMETER)

           call error_throw(ERROR_PARAMETER)

           return
       end if

   end do
   
   call error_clear_catch(ERROR_PARAMETER)

end subroutine

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

End Module MC
