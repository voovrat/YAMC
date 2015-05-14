Module EwaldSumKSpace !DOC
use FourierGrid
use SumSinCosKR
use LJTypes
use AtomicData
use SystemSettings, only : TRealArrayPointer
use RhoSquared

!DOC !FILE This Module contains structures and functions which deal with the K-Space component of the Ewald Sum (both LJ and Coulomb)

! Calculates
!
!  U_fourier = U_fourier^C + U_fourier^LJ
!
!  Where
!      U_fourier^C = SUM_m beta(k_m) |rho(k_m)|^2
!      U_fourier^LJ = SUM_m (beta6(k_m) + beta12(k_m)) F_2(k_m)  
!
!
!  beta = 4pi/k_m exp(-k_m^2 / 4 alpha^2) 
!
!  beta6 = alpha^3/V B_6(b = pi m/L/alpha )
!  B6(b) = pi^1.5 / 3 ( b^3 sqrt(pi) erfc(b) + (0.5-b^2) exp(-b^2) )
!
!  beta12 = alpha^9 /V B12(b) 
!  B12(b) = pi^1.5 / (945*120) * ( -16b^9 sqrt(pi) erfc(b) +
!                                  +(16b^8 - 8b^6 + 12b^4 - 30b^2 + 105 )exp(-b^2) )
!
!
!  |rho(k_m)|^2 = (SUM_j q_j exp( 2 pi m/L r_j) ) ^2 
!               = C_m^2 + S_m^2
!
!  where C_m = SUM_j cos(2 pi m/L r_j) 
!        S_m = SUM_j sin(2 pi m/L r_j)
!
!
!  F_2(k_m) = SUM_sj A_sj exp( 2pi m/L (r_s - r_j) )
!           = SUM_t1t2 A_t1t2 F^t1(k_m) * F^t2(-km) 
!           = SUM_t1t2 A_t1t2 (C^t1_m * C^t2_m + S^t1_m S^t2_m)
!
!  where C^t_m = SUM_j=1^N(t) cos( 2 pi m/L * r_j^t )
!        S^t_m = SUM_j=1^N(t) sin( 2 pi m/L & r_j^t )
!
! To reduce number of operations the sums are also rewritten for mx,my,mz >= 0:
!
! U_fourier^C = SUM_mx,my,mz>=0 (1 + sgn(mx)) beta(k_m) 
!                               ( rho(+++)^2 + sgn(my)rho(+-+)^2 + sgn(mz)rho(++-)^2 + sgn(my*mz)rho(+--)^2)
! U_fourier^LJ = U_fourier_LJ12 - U_fourier_LJ12
! U_fourier_LJ6 = SUM_mx,my,mz>= (1 + sgn(mx)) beta6(k_m) 
!                               ( F2^LJ6(+++) + sgn(my)F2^LJ6(+-+) + sgn(mz)F2^LJ6(++-) + sgn(my*mz)F2^LJ6(+--) )
! 
! U_fourier_LJ12 = SUM_mx,my,mz>= (1 + sgn(mx)) beta12(k_m) 
!                               ( F2^LJ12(+++) + sgn(my)F2^LJ12(+-+) + sgn(mz)F2^LJ12(++-) + sgn(my*mz)F2^LJ12(+--) )
!
implicit none


Type :: TEwaldSumKSpace  !DOC
!DOC Fields:
   Type(TFourierGrid) :: grid !DOC K-Space grid

   Type(TRhoSquared) :: rho_squared_total  !DOC RhoSquared structure which contains the SUM cos(kR), SUM_i sin(k*R_i) for Coulomb and LJ (see Module RhoSquared) 
  
   Type(TAtomicData),pointer :: atomic_data !DOC Pointer to the AtomicData structure
   Type(TLJTypes),pointer :: lj_types !DOC LJTypes structure which contains the individual and pair LJ parameters for the atoms

   real(8),dimension(:),allocatable :: beta, beta6, beta12  !DOC beta(k) (for coulomb), beta6(k) and beta12(k) - for LJ Ewald 
   Type(TRealArrayPointer),dimension(:),allocatable :: beta_LJ   !DOC beta_LJ^t1t2(k) = 4 epsilon( sigma^12 beta12(k) - sigma6 beta6(k))
                                                          ! where epsilon = epsilon(t1,t2), sigma=sigma(t1,t2)
                                                          ! 2D array is implemented via 1D:
                                                          !   beta_LJ(t1 * NType + t2) is beta_LJ(t1,t2)
                                                          !   elements are physically allocated only for t2>=t1. 
                                                          !  element (t1,t2) points the same address in memory as (t2,t1)

   real(8) :: energy,energy_coulomb,energy_LJ !DOC total energy and energy components

   !DOC Temporary arrays for delta_LJ_energy calculations (EwaldSumKSpace_calculate_dU_LJ)
   ! (2 C,S(sxsysz) + DeltaC,S(sxsysz) )
   real(8),dimension(:),allocatable :: twoC_plus_dC_ppp,twoC_plus_dC_ppm,twoC_plus_dC_pmp,twoC_plus_dC_pmm  !DOC
   real(8),dimension(:),allocatable :: twoS_plus_dS_ppp,twoS_plus_dS_ppm,twoS_plus_dS_pmp,twoS_plus_dS_pmm !DOC


End Type TEwaldSumKSpace


contains

subroutine EwaldSumKSpace_alloc(this,atomic_data,lj_types) !DOC
!DOC Allocate the EwaldSumKSpace structure
   use parameters, only : kmax,alpha
!DOC Parameters:
   Type(TEwaldSumKSpace) :: this !DOC EwaldSumKSpace structure
   Type(TAtomicData),intent(in),target :: atomic_data !DOC AtomicData (coordinates + charges)
   Type(TLJTypes),intent(in),target :: lj_types!DOC LJTypes array (sigma,epsilon for each pair)
     
   integer :: nk ! number of points in k-space
   integer :: ntypes  
   integer :: t1,t2   ! type iterators
   integer :: offset_12,offset_21
 
   this % atomic_data => atomic_data
   this % lj_types => lj_types

! allocate kspace grid   
   call FourierGrid_calcNalloc(kmax,nk) 
   call FourierGrid_alloc(this % grid,nk)  

 ! initialize k-vectors  (because otherwise rho_squared_total cannot be allocated)
  call FourierGrid_init(this % grid ,kmax)
 
! allocate beta
   allocate( this % beta(nk) )
   allocate( this % beta6(nk) )
   allocate( this % beta12(nk) )

! allocate beta_LJ
! it is 2D array stored in 1D, element (t1,t2) corresponds to offset (t1-1)*ntypes + t2
! beta_LJ is symmetric, so only elements with t2>=t1 are allocated, and (t2,t1) just points to (t1,t2)
   ntypes = lj_types % NType
   allocate( this % beta_LJ( ntypes**2) )

   do t1 = 1,ntypes
   do t2 = t1,ntypes

      offset_12 =  (t1-1)*ntypes + t2
      allocate( this % beta_LJ( offset_12 ) % ptr(nk) )

      if ( t1 /= t2 ) then
         offset_21 = (t2-1)*ntypes + t1
         this % beta_LJ( offset_21 ) % ptr => this % beta_LJ( offset_12 ) % ptr
      end if

   end do
   end do

  ! allocate rho_squared_total
  call RhoSquared_alloc(this % rho_squared_total, this % grid, lj_types)

  allocate(this % twoC_plus_dC_ppp(nk) )
  allocate(this % twoC_plus_dC_pmp(nk) )
  allocate(this % twoC_plus_dC_ppm(nk) )
  allocate(this % twoC_plus_dC_pmm(nk) )

  allocate(this % twoS_plus_dS_ppp(nk) )
  allocate(this % twoS_plus_dS_pmp(nk) )
  allocate(this % twoS_plus_dS_ppm(nk) )
  allocate(this % twoS_plus_dS_pmm(nk) )

end subroutine

subroutine EwaldSumKSpace_dealloc(this) !DOC
!DOC Deallocate the EwaldSumKSpace structure
!DOC Parameters:
   Type(TEwaldSumKSpace) :: this !DOC EwaldSumKSpace structure
  
   integer :: t1,t2
   integer :: offset_12
   integer :: ntypes

! deallocate beta
   deallocate( this % beta )
   deallocate( this % beta6 )
   deallocate( this % beta12 )

! allocate beta_LJ
! it is 2D array stored in 1D, element (t1,t2) corresponds to offset (t1-1)*ntypes + t2
! beta_LJ is symmetric, so only elements with t2>=t1 are allocated, and (t2,t1) just points to (t1,t2)
   ntypes = this % lj_types % NType

   do t1 = 1,ntypes
   do t2 = t1,ntypes

      offset_12 =  (t1-1)*ntypes + t2
      deallocate( this % beta_LJ(offset_12) % ptr )

   end do
   end do

   deallocate( this % beta_LJ )

  ! allocate rho_squared_total
  call RhoSquared_dealloc(this % rho_squared_total)

  deallocate(this % twoC_plus_dC_ppp )
  deallocate(this % twoC_plus_dC_pmp )
  deallocate(this % twoC_plus_dC_ppm )
  deallocate(this % twoC_plus_dC_pmm )

  deallocate(this % twoS_plus_dS_ppp )
  deallocate(this % twoS_plus_dS_pmp )
  deallocate(this % twoS_plus_dS_ppm )
  deallocate(this % twoS_plus_dS_pmm )

end subroutine


subroutine EwaldSumKSpace_init(this,scale_beta) !DOC
   use parameters, only : kmax 
!DOC initialize all the arrays needed for the calculation
!DOC Parameters:
   Type(TEwaldSumKSpace) :: this !DOC EwaldSumKSpace structure
!   integer,intent(in) :: kmax
   real(8),intent(in),optional :: scale_beta !DOC if is used within the volume change step the beta(k) should be scaled. Then scale_beta is the coefficient. Otherwise, it can be not given (optional) or zero. Then the beta(k) will be re-calculated

   integer :: natom
   integer :: nk ! number of grid points in kspace

  
! initialize beta functions
 if( ( present(scale_beta) ).and. (scale_beta > 1e-9 ) ) then
   this % beta(:) = this % beta(:) * scale_beta
 else
   call EwaldSumKSpace_initBeta(this) ! beta,beta6,beta12
   call EwaldSumKSpace_initBetaLJ(this)
 end if
!i
! initialize coulomb and LJ energies ( and sin/cos sums )
! 
   call EwaldSumKSpace_calc_total_energy( this )

! initialize sumsincos arrays for coulomb and LJ ( C(+++), C(++-) etc arrays... see in the begining of module)

end subroutine 

subroutine EwaldSumKSpace_calc_total_energy( this ) !DOC
!DOC Calculate the Ewald KSpace component of the energy
!DOC Parameters:
   Type(TEwaldSumKSpace) :: this  !DOC EwaldSumKSpace structure
   real(8),dimension(:),pointer :: xx,yy,zz,charge 
   integer,dimension(:),pointer :: type_by_index   
   integer :: natom

   natom = this % atomic_data % natom 
 
   xx => this % atomic_data % xx
   yy => this % atomic_data % yy
   zz => this % atomic_data % zz

   charge => this % atomic_data % charge
   type_by_index => this % lj_types % type_by_index 


   call RhoSquared_nulify( this % rho_squared_total )
   call RhoSquared_addAtoms( this % rho_squared_total, xx, yy, zz, charge, type_by_index, natom )

!   call RhoSquared_calcSumSinCos_coulomb(this % rho_squared_total,xx,yy,zz,chrge,natom) ! coulomb
!   call RhoSquared_calcSumSinCos_LJ(this % rho_squared_total,xx,yy,zz,type_by_index,natom) ! LJ

! calculate rho_squared (for coulomb) 
   call RhoSquared_calc_rho_squared(this % rho_squared_total )
!   call RhoSquared_calc_F2(this) ! obsolete, se 

! now we are ready to calculate the enery

   this % energy_coulomb = EwaldSumKSpace_calc_coulomb_energy(this)
   this % energy_LJ = EwaldSumKSpace_calc_LJ_energy(this)

   this % energy = this % energy_coulomb + this % energy_LJ

end subroutine


subroutine EwaldSumKSpace_initBeta(this) !DOC 
!DOC Calculate the beta function. Is called from EwaldSumKSpace_init
!DOC Parameters:
   use constants, only : pi,two_pi
   use parameters, only : alpha, xclb
   use Functions, only : Erfc_Luc_bet6_bet12
! computes beta,beta6,beta12
   Type(TEwaldSumKSpace),target :: this !DOC EwaldSumKSpace structure
   
   real(8) :: pi_over_alpha  ! pi/alpha
   real(8) :: pi_over_alpha_squared  ! pi^2 / alpha^2
   real(8) :: b  ! b = pi m / alpha L  = pi k /alpha (our notation)
   real(8) :: b2 ! pi^2 / alpha^2 * k^2 = b^2
   real(8) :: exp_minus_b2  ! exp(-pi^2/alpha^2 * k^2) = exp(-b^2)?  b = pi h/ alpha = pi m / (alpha*L)
                            ! in our notation: k==m,   L==1
   real(8) :: erfc_b        ! erfc(b) returned by Erfc_Luc_bet6_bet12
   real(8) :: B6,B12        ! B6(b), B12(b) functions (see the definitions above) calculated by Erfc_Luc_bet6_bet12


   integer :: kk

! temporaries and pointers
   integer :: nk 
   real(8),dimension(:),pointer :: xk,xk2  !  m,m^2 (k^2 in program notation)  taken form grid object 
   integer,dimension(:),pointer :: kx   ! taken from grid
 
!   alpha = this % alpha

   xk2 => this % grid % xk2
   xk => this % grid % xk
   kx => this % grid % kx
 
   pi_over_alpha = pi / alpha
   pi_over_alpha_squared = pi_over_alpha**2

   do kk=1,this % grid % nk
      
      ! pialp2 = pi^2 / alpha^2
      b2 = pi_over_alpha_squared * xk2(kk) ! pi^2 / alpha^2 * k^2 = b^2 

      !  k_m = 2*pi*m/L 
      !  exp( - k_m^2 / 4 alpha^2 ) = exp( - 4 pi^2 m^2 / L^2 / 4 alpha^2) = 
      !    = exp( - pi^2 * m^2 / L^2 /alpha^2 ) 
      !  now, in our notation m <-> k, L=1
      !    exp ( -k_m^2 / 4 alpha^2 ) = exp( - pi^2 *k^2 / alpha^2 ) 
      !                       
      exp_minus_b2 = EXP(-b2)     ! exp(-pi^2/alpha^2 * k^2) 
      this % beta(kk) = xclb * exp_minus_b2 / xk2(kk) / two_pi          ! 1/2piL pour charge-charge
                          ! bet(kk) = 1/2pi/k^2 * exp(-pi^2 / alpha^2 * k^2 ) 
                          ! why 1/2pi: 1/2 is the sum prefactor
                          ! 4pi/k_m^2 = 4pi/ (4pi^2 L^2/m^2) = 1/ (pi * m^2) = 1 / (pi* k^2) (in our notation)
      !xk=SQRT(xk2)          

!     write(*,*) 'b2',b2,'exp(-b2)',exp_minus_b2,'beta',this % beta(kk)

      b=pi_over_alpha * xk(kk)       ! pi/alpha * k
   
   
      ! calculation of beta6 and beta12 for LJ sums in Fourier space
      ! 
      !      erfk=erfc_Luc(pialpk,pialpk2,eek)
      !      bet6(kk)=pi45/3.d0*xk2*xk*(pi12*erfk+(0.5d0/pialpk2-1.d0)/pialpk*eek)
      !      bet12(kk)=pi105/1080.d0*xk2**4*xk*(-16.d0*pi12/105.d0*erfk+           &
      !         (((((1.d0/pialpk2-2.d0/7.d0)/pialpk2+4.d0/35.d0)/pialpk2-8.d0/105.d0)/pialpk2+16.d0/105.d0))/pialpk*eek)
      !      nouveau: on peut utiliser un DL ou/et une integration numerique partielle  luc84p150
      !
      ! ( for details:  see definition of erfc_Luc_bet6_bet12) (in Module Functions)
      ! erfc_b, B6, B12 - output parameters
      call erfc_Luc_bet6_bet12(b,b2,exp_minus_b2,erfc_b,B6,B12)

      this % beta6(kk) = alpha**3*B6  ! alpha^3 *   pi^1.5/3V * b^3 (sqrt(pi) erfc(b) + (1/2b^3 -1b)exp(-b^2) ) 
      this % beta12(kk) = alpha**9*B12 ! alpha^9 *  1/(945*120) *pi^1.5 *  b^9 * 
                           !            ( -16 sqrt(pi) erfc(b) 
                           !              + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) 
                           !             ) 

      if(kx(kk)>0) then                     ! pour compter -kx
                                            ! coefficient (1+sgn(mx)) is included to the beta (see the definition of sums above) 
          this % beta(kk) = 2.d0*this % beta(kk)
          this % beta6(kk) = 2.d0*this % beta6(kk)
          this % beta12(kk) = 2.d0*this % beta12(kk)

          
      end if

   end do ! grid (kk)
 
end subroutine ! init beta

subroutine EwaldSumKSpace_initBetaLJ(this) !DOC
!DOC Calculate the betaLJ(k)
!DOC betaLJ(t1t2) = epsilon(t1,t2)* ( sigma(t1,t2)^12* beta12 - sigma(t1,t2)^6 * beta6)
!DOC Parameters:
   Type(TEwaldSumKSpace) :: this !DOC EwaldSumKSpace structure

   integer :: t1,t2
   integer :: ntype
   Type(TLJTypes),pointer :: lj_types
   real(8) :: ALJ6,ALJ12   
   integer :: base_offset, offset

   lj_types => this % lj_types
   ntype = lj_types % NType

   base_offset = 0  ! base_offset = (t1-1) * ntype
   do t1 = 1,ntype

     do t2 = t1,ntype
   
          offset = base_offset + t2 ! (t1-1) * ntype + t2

          ALJ6  = lj_types % LJ6Tab( offset )
          ALJ12 = lj_types % LJ12Tab( offset )
   
          this % beta_LJ( offset ) % ptr(:) = ALJ12 * this % beta12(:) - ALJ6 * this % beta6(:)
          ! ( beta_LJ(t2,t1) % ptr points to the same place) 
  
      end do

      base_offset = base_offset + ntype

   end do

end subroutine


function EwaldSumKSpace_calc_coulomb_energy(this) !DOC
!DOC calculate the KSpace coulomb energy component
    real(8) :: EwaldSumKSpace_calc_coulomb_energy
!DOC which is SUM beta(k) rho_squared(k)
!DOC Parameters:
    Type(TEwaldSumKSpace),intent(in),target :: this !DOC EwaldSumKSpace structure
!DOC Return value:
!DOC coulomb energy component
    real(8),dimension(:),pointer :: beta
    real(8),dimension(:),pointer :: rho_squared

    real(8) :: U
    integer :: kk,nk

    beta => this % beta
    rho_squared => this % rho_squared_total % rho_squared

    nk = this % grid % nk

   U = 0
   do kk=1,nk
      U = U + beta(kk) * rho_squared(kk)
  !    write(*,*) 'beta=',beta(kk),'rho_squared=',rho_squared(kk)
   end do

   EwaldSumKSpace_calc_coulomb_energy =  U

end function

function EwaldSumKSpace_calc_LJ_energy(this) !DOC 
   real(8) :: EwaldSumKSpace_calc_LJ_energy
!DOC calculates LJ energy using the beta_LJ:
!     U = SUM_t1 SUM t2>=t1 (2- delta_t1t2) SUM_k betaLJ_t1t2(k) * (
!                                 ( C^t1+++ C^t2+++ + S^t1+++ S^t2+++)
!                       +sgn(ky) *( C^t1+-+ C^t2+-+ + S^t1+-+ S^t2+-+) 
!                       +sgn(kz) *( C^t1++- C^t2++- + S^t1++- S^t2++-) 
!                     +sgn(kykz) *( C^t1+-- C^t2+-- + S^t1+-- S^t2+--) ) 
!
!DOC  note : C^t(sxsysz),S^t(sxsysz) should be pre-calculated ( rho_square_total )
!DOC Parameters:   
   Type(TEwaldSumKSpace),intent(in),target :: this !DOC EwaldSumKSpace strucutre
!DOC Return value:
!DOC LJ component of the energy
   integer :: t1,t2
   integer :: ntype
   integer :: base_offset ! (t1-1)*ntypes
!   integer :: offset ! (t1-1)*ntypes + t2
   integer :: kk  ! iterator in k-space
   integer :: nk ! number of points in k-sapce
   real(8) :: Sum_k ! inner sum_k
   real(8) :: sumsincos_local ! very inner sum of all combinations of C1C2 (for each t1,t2,k)
   real(8) :: U_LJ ! total sum

 
   integer,dimension(:),pointer :: kx,ky,kz 

   real(8),dimension(:),pointer :: C1_ppp,C1_pmp,C1_ppm,C1_pmm ! C^t1(sxsysz)
   real(8),dimension(:),pointer :: S1_ppp,S1_pmp,S1_ppm,S1_pmm ! S^t1(sxsysz)

   real(8),dimension(:),pointer :: C2_ppp,C2_pmp,C2_ppm,C2_pmm ! C^t2(sxsysz)
   real(8),dimension(:),pointer :: S2_ppp,S2_pmp,S2_ppm,S2_pmm ! S^t2(sxsysz)

   Type(TSumSinCosKR),dimension(:),pointer :: sumsincos_LJ
   real(8),dimension(:),pointer :: betaLJ_t1t2

   kx => this % grid % kx
   ky => this % grid % ky
   kz => this % grid % kz

   nk = this % grid % nk

   ntype = this % lj_types % NType

   sumsincos_LJ => this % rho_squared_total % sumsincos_LJ

   U_LJ = 0

   base_offset = 0
   do t1=1,ntype

      C1_ppp => sumsincos_LJ(t1) % sumcos_ppp   ! C^t1(+++)
      S1_ppp => sumsincos_LJ(t1) % sumsin_ppp   ! S^t1(+++)
 
      C1_pmp => sumsincos_LJ(t1) % sumcos_pmp   ! C^t1(+-+)
      S1_pmp => sumsincos_LJ(t1) % sumsin_pmp   ! S^t1(+-+)

      C1_ppm => sumsincos_LJ(t1) % sumcos_ppm   ! C^t1(++-)
      S1_ppm => sumsincos_LJ(t1) % sumsin_ppm   ! S^t1(++-)
  
      C1_pmm => sumsincos_LJ(t1) % sumcos_pmm   ! C^t1(+--)
      S1_pmm => sumsincos_LJ(t1) % sumsin_pmm   ! S^t1(+--)
   

      do t2=t1,ntype
  
         betaLJ_t1t2 => this % beta_LJ( base_offset + t2) % ptr 

         C2_ppp => sumsincos_LJ(t2) % sumcos_ppp   ! C^t2(+++)
         S2_ppp => sumsincos_LJ(t2) % sumsin_ppp   ! S^t2(+++)
    
         C2_pmp => sumsincos_LJ(t2) % sumcos_pmp   ! C^t2(+-+)
         S2_pmp => sumsincos_LJ(t2) % sumsin_pmp   ! S^t2(+-+)
   
         C2_ppm => sumsincos_LJ(t2) % sumcos_ppm   ! C^t2(++-)
         S2_ppm => sumsincos_LJ(t2) % sumsin_ppm   ! S^t2(++-)
     
         C2_pmm => sumsincos_LJ(t2) % sumcos_pmm   ! C^t2(+--)
         S2_pmm => sumsincos_LJ(t2) % sumsin_pmm   ! S^t2(+--)

         Sum_k = 0 

         do kk=1,nk

            sumsincos_local = C1_ppp(kk) * C2_ppp(kk) + S1_ppp(kk) * S2_ppp(kk)
            if( ky(kk) > 0 ) sumsincos_local = sumsincos_local + C1_pmp(kk) * C2_pmp(kk) + S1_pmp(kk) * S2_pmp(kk)
            if( kz(kk) > 0 ) sumsincos_local = sumsincos_local + C1_ppm(kk) * C2_ppm(kk) + S1_ppm(kk) * S2_ppm(kk)
            if( (ky(kk) > 0) .AND. (kz(kk) > 0) ) then
                sumsincos_local = sumsincos_local + C1_pmm(kk) * C2_pmm(kk) + S1_pmm(kk) * S2_pmm(kk)
            end if 

            Sum_k = Sum_k + betaLJ_t1t2(kk) * sumsincos_local

         end do ! kk

         if ( t1 /= t2 ) then ! (2 - delta_t1t2) sum prefactor  ( to account that Sum_t1t2 = Sum_t2t1 ) 
              Sum_k= Sum_k + Sum_k  
         end if 
   
         U_LJ = U_LJ + Sum_k

      end do ! t2

      base_offset = base_offset + ntype
   end do ! t1

   EwaldSumKSpace_calc_LJ_energy = U_LJ

end function

function  EwaldSumKSpace_calc_dU_coulomb( this, delta_sumsincos_coulomb ) !DOC
!DOC Calculates the coulomb energy change dU
   real(8) :: EwaldSumKSpace_calc_dU_coulomb
!DOC Parameters:
   Type(TEwaldSumKSpace),target :: this !DOC EwaldSumKSpace structure
   Type(TSumSinCosKR),target :: delta_sumsincos_coulomb !DOC changes of the sums of sins and cos (for the given molecule) (see Module SumSinCosKR)
!DOC Return value:
!DOC   coulomb energy change dU

! Delta U_C = SUM_m>=0 (1+sgn(mx)) beta(k_m)*(
!                       ( (2 C+++ + DeltaC+++)*DeltaC+++ + (2 S+++ + DeltaS+++) DeltaS+++ )                          
!             sgn(ky) * ( (2 C+-+ + DeltaC+-+)*DeltaC+-+ + (2 S+-+ + DeltaS+-+) DeltaS+-+ )
!             sgn(kz) * ( (2 C++- + DeltaC++-)*DeltaC++- + (2 S++- + DeltaS++-) DeltaS++- )
!           sgn(ky*kz)* ( (2 C+-- + DeltaC+--)*DeltaC+-- + (2 S+-- + DeltaS+--) DeltaS+-- ) )
   integer :: kk,nk
   integer,dimension(:),pointer :: ky,kz
   real(8),dimension(:),pointer :: beta

   real(8) :: local_sum,total_sum 

   real(8),dimension(:),pointer :: C_ppp,C_ppm,C_pmp,C_pmm,S_ppp,S_ppm,S_pmp,S_pmm ! S,C(sx,sy,sz)
   real(8),dimension(:),pointer :: DeltaC_ppp,DeltaC_ppm,DeltaC_pmp,DeltaC_pmm
   real(8),dimension(:),pointer :: DeltaS_ppp,DeltaS_ppm,DeltaS_pmp,DeltaS_pmm 

   real(8) :: two_C,two_S,dC,dS ! temporary variables

   nk = this % grid % nk
   ky => this % grid % ky
   kz => this % grid % kz

   beta => this % beta

   C_ppp => this % rho_squared_total % sumsincos_coulomb % sumcos_ppp
   C_ppm => this % rho_squared_total % sumsincos_coulomb % sumcos_ppm
   C_pmp => this % rho_squared_total % sumsincos_coulomb % sumcos_pmp
   C_pmm => this % rho_squared_total % sumsincos_coulomb % sumcos_pmm

   S_ppp => this % rho_squared_total % sumsincos_coulomb % sumsin_ppp
   S_ppm => this % rho_squared_total % sumsincos_coulomb % sumsin_ppm
   S_pmp => this % rho_squared_total % sumsincos_coulomb % sumsin_pmp
   S_pmm => this % rho_squared_total % sumsincos_coulomb % sumsin_pmm

   DeltaC_ppp => delta_sumsincos_coulomb % sumcos_ppp
   DeltaC_ppm => delta_sumsincos_coulomb % sumcos_ppm
   DeltaC_pmp => delta_sumsincos_coulomb % sumcos_pmp
   DeltaC_pmm => delta_sumsincos_coulomb % sumcos_pmm

   DeltaS_ppp => delta_sumsincos_coulomb % sumsin_ppp
   DeltaS_ppm => delta_sumsincos_coulomb % sumsin_ppm
   DeltaS_pmp => delta_sumsincos_coulomb % sumsin_pmp
   DeltaS_pmm => delta_sumsincos_coulomb % sumsin_pmm

   total_sum = 0
   do kk=1,nk

      two_C = C_ppp(kk) + C_ppp(kk) ! this should be faster than 2*C_ppp(kk)
      dC = DeltaC_ppp(kk)
      two_S = S_ppp(kk) + S_ppp(kk)
      dS = DeltaS_ppp(kk)
      local_sum = ( two_C + dC) * dC + (two_S + dS) * dS

      if( ky(kk) > 0 ) then
         two_C = C_pmp(kk) + C_pmp(kk) ! this should be faster than 2*C_ppp(kk)
         dC = DeltaC_pmp(kk)
         two_S = S_pmp(kk) + S_pmp(kk)
         dS = DeltaS_pmp(kk)
         local_sum = local_sum + ( two_C + dC) * dC + (two_S + dS) * dS
      end if

      if( kz(kk) > 0 ) then
         two_C = C_ppm(kk) + C_ppm(kk) ! this should be faster than 2*C_ppp(kk)
         dC = DeltaC_ppm(kk)
         two_S = S_ppm(kk) + S_ppm(kk)
         dS = DeltaS_ppm(kk)
         local_sum = local_sum + ( two_C + dC) * dC + (two_S + dS) * dS
      end if

      if( (ky(kk) > 0) .AND. (kz(kk) > 0 ) ) then
         two_C = C_pmm(kk) + C_pmm(kk) ! this should be faster than 2*C_ppp(kk)
         dC = DeltaC_pmm(kk)
         two_S = S_pmm(kk) + S_pmm(kk)
         dS = DeltaS_pmm(kk)
         local_sum = local_sum + ( two_C + dC) * dC + (two_S + dS) * dS
      end if

!      write(*,*) 'twoC=',two_C,'dC=',dC
      total_sum = total_sum + local_sum * beta(kk)

   end do
 
   EwaldSumKSpace_calc_dU_coulomb = total_sum

end function


function  EwaldSumKSpace_calc_dU_LJ( this, delta_sumsincos_LJ,type_is_present ) !DOC 
!DOC calculate the LJ energy change
   real(8) :: EwaldSumKSpace_calc_dU_LJ
!DOC Parameters:
   Type(TEwaldSumKSpace),target :: this !EwaldSumKSpace structure
   Type(TSumSinCosKR),dimension(:),intent(in),target :: delta_sumsincos_LJ !DOC changes of sums of sin and cos (for a given molecule).
   logical,dimension(:) :: type_is_present !DOC logical array which indicates which LJ diameters are present in the delta_sumsincos_LJ
!DOC Return value:
!DOC Change of the LJ energy dU_LJ

! Delta U_LJ = SUM_t1 SUM_(t2 in M) SUM_m>=0 (1+sgn(mx)) beta(k_m)*(
!                       ( (2 C^t1+++ + DeltaC^t1+++)*DeltaC^t2+++ + (2 S^t1+++ + DeltaS^t1+++) DeltaS^t2+++ )                          
!             sgn(ky) * ( (2 C^t1+-+ + DeltaC^t1+-+)*DeltaC^t2+-+ + (2 S^t1+-+ + DeltaS^t1+-+) DeltaS^t2+-+ )
!             sgn(kz) * ( (2 C^t1++- + DeltaC^t1++-)*DeltaC^t2++- + (2 S^t1++- + DeltaS^t1++-) DeltaS^t2++- )
!           sgn(ky*kz)* ( (2 C^t1+-- + DeltaC^t1+--)*DeltaC^t2+-- + (2 S^t1+-- + DeltaS^t1+--) DeltaS^t2+-- ) )
   integer :: t1,t2
   integer :: base_offset ! (t1-1) * ntype
   integer :: ntype
   integer :: kk,nk
   integer,dimension(:),pointer :: ky,kz
   real(8),dimension(:),pointer :: betaLJ_t1t2

   real(8) :: local_sum,total_sum 

   real(8),dimension(:),pointer :: C1_ppp,C1_ppm,C1_pmp,C1_pmm,S1_ppp,S1_ppm,S1_pmp,S1_pmm ! S,C(sx,sy,sz)
   real(8),dimension(:),pointer :: DeltaC1_ppp,DeltaC1_ppm,DeltaC1_pmp,DeltaC1_pmm
   real(8),dimension(:),pointer :: DeltaC2_ppp,DeltaC2_ppm,DeltaC2_pmp,DeltaC2_pmm

   real(8),dimension(:),pointer :: DeltaS1_ppp,DeltaS1_ppm,DeltaS1_pmp,DeltaS1_pmm 
   real(8),dimension(:),pointer :: DeltaS2_ppp,DeltaS2_ppm,DeltaS2_pmp,DeltaS2_pmm 

   real(8),dimension(:),pointer :: twoC_plus_dC_ppp,twoC_plus_dC_ppm,twoC_plus_dC_pmp,twoC_plus_dC_pmm
   real(8),dimension(:),pointer :: twoS_plus_dS_ppp,twoS_plus_dS_ppm,twoS_plus_dS_pmp,twoS_plus_dS_pmm

   real(8) :: two_C,two_S,dC,dS ! temporary variables

   nk = this % grid % nk
   ky => this % grid % ky
   kz => this % grid % kz

   twoC_plus_dC_ppp => this % twoC_plus_dC_ppp
   twoC_plus_dC_pmp => this % twoC_plus_dC_pmp
   twoC_plus_dC_ppm => this % twoC_plus_dC_ppm
   twoC_plus_dC_pmm => this % twoC_plus_dC_pmm

   twoS_plus_dS_ppp => this % twoS_plus_dS_ppp
   twoS_plus_dS_pmp => this % twoS_plus_dS_pmp
   twoS_plus_dS_ppm => this % twoS_plus_dS_ppm
   twoS_plus_dS_pmm => this % twoS_plus_dS_pmm

   ntype = this % lj_types % NType

   total_sum = 0
   base_offset = 0 ! (t1-1)*ntype
   do t1 = 1,ntype
      
      C1_ppp => this % rho_squared_total % sumsincos_LJ(t1) % sumcos_ppp
      C1_ppm => this % rho_squared_total % sumsincos_LJ(t1) % sumcos_ppm
      C1_pmp => this % rho_squared_total % sumsincos_LJ(t1) % sumcos_pmp
      C1_pmm => this % rho_squared_total % sumsincos_LJ(t1) % sumcos_pmm
   
      S1_ppp => this % rho_squared_total % sumsincos_LJ(t1) % sumsin_ppp
      S1_ppm => this % rho_squared_total % sumsincos_LJ(t1) % sumsin_ppm
      S1_pmp => this % rho_squared_total % sumsincos_LJ(t1) % sumsin_pmp
      S1_pmm => this % rho_squared_total % sumsincos_LJ(t1) % sumsin_pmm
   
      DeltaC1_ppp => delta_sumsincos_LJ(t1) % sumcos_ppp
      DeltaC1_ppm => delta_sumsincos_LJ(t1) % sumcos_ppm
      DeltaC1_pmp => delta_sumsincos_LJ(t1) % sumcos_pmp
      DeltaC1_pmm => delta_sumsincos_LJ(t1) % sumcos_pmm
   
      DeltaS1_ppp => delta_sumsincos_LJ(t1) % sumsin_ppp
      DeltaS1_ppm => delta_sumsincos_LJ(t1) % sumsin_ppm
      DeltaS1_pmp => delta_sumsincos_LJ(t1) % sumsin_pmp
      DeltaS1_pmm => delta_sumsincos_LJ(t1) % sumsin_pmm
 
      twoC_plus_dC_ppp(1:nk) = C1_ppp(1:nk) + C1_ppp(1:nk) + DeltaC1_ppp(1:nk)
      twoC_plus_dC_ppm(1:nk) = C1_ppm(1:nk) + C1_ppm(1:nk) + DeltaC1_ppm(1:nk)
      twoC_plus_dC_pmp(1:nk) = C1_pmp(1:nk) + C1_pmp(1:nk) + DeltaC1_pmp(1:nk)
      twoC_plus_dC_pmm(1:nk) = C1_pmm(1:nk) + C1_pmm(1:nk) + DeltaC1_pmm(1:nk)
          
      twoS_plus_dS_ppp(1:nk) = S1_ppp(1:nk) + S1_ppp(1:nk) + DeltaS1_ppp(1:nk)
      twoS_plus_dS_ppm(1:nk) = S1_ppm(1:nk) + S1_ppm(1:nk) + DeltaS1_ppm(1:nk)
      twoS_plus_dS_pmp(1:nk) = S1_pmp(1:nk) + S1_pmp(1:nk) + DeltaS1_pmp(1:nk)
      twoS_plus_dS_pmm(1:nk) = S1_pmm(1:nk) + S1_pmm(1:nk) + DeltaS1_pmm(1:nk)
 
      do t2=1,ntype

         if(.not.type_is_present(t2)) cycle

         betaLJ_t1t2 => this % beta_LJ( base_offset + t2 ) % ptr

         DeltaC2_ppp => delta_sumsincos_LJ(t2) % sumcos_ppp
         DeltaC2_ppm => delta_sumsincos_LJ(t2) % sumcos_ppm
         DeltaC2_pmp => delta_sumsincos_LJ(t2) % sumcos_pmp
         DeltaC2_pmm => delta_sumsincos_LJ(t2) % sumcos_pmm
      
         DeltaS2_ppp => delta_sumsincos_LJ(t2) % sumsin_ppp
         DeltaS2_ppm => delta_sumsincos_LJ(t2) % sumsin_ppm
         DeltaS2_pmp => delta_sumsincos_LJ(t2) % sumsin_pmp
         DeltaS2_pmm => delta_sumsincos_LJ(t2) % sumsin_pmm
         
         do kk=1,nk
      
            local_sum =  twoC_plus_dC_ppp(kk) * DeltaC2_ppp(kk) + twoS_plus_dS_ppp(kk) * DeltaS2_ppp(kk)
      
            if( ky(kk) > 0 ) then
               local_sum = local_sum +  twoC_plus_dC_pmp(kk) * DeltaC2_pmp(kk) + twoS_plus_dS_pmp(kk) * DeltaS2_pmp(kk)
            end if
      
            if( kz(kk) > 0 ) then
               local_sum = local_sum +  twoC_plus_dC_ppm(kk) * DeltaC2_ppm(kk) + twoS_plus_dS_ppm(kk) * DeltaS2_ppm(kk)
            end if
      
            if( (ky(kk) > 0) .AND. (kz(kk) > 0 ) ) then
               local_sum = local_sum +  twoC_plus_dC_pmm(kk) * DeltaC2_pmm(kk) + twoS_plus_dS_pmm(kk) * DeltaS2_pmm(kk)
            end if
      
            total_sum = total_sum + local_sum * betaLJ_t1t2(kk)
      
         end do ! kk

      end do ! t2

      base_offset = base_offset + ntype

   end do  ! t1

   EwaldSumKSpace_calc_dU_LJ = total_sum

end function



!  *************   OBSOLETE  ********************
!
!! U_c = SUM_k rho^2(k) * (1 + sgn(mx)) beta(k)  ( where (1+ sgn(mx)) is already included in beta array in our case. see initBeta )
!pure function EwaldSumKSpace_calc_coulomb_energy(this)
!! note:  this % rho_squared should be calculated before (by calc_rho_squared for example)
!! note:  beta should be pre-calculated (by initBeta normally) 
!   real(8) :: EwaldSumKSpace_calc_coulomb_energy
!
!   Type(TEwaldSumKSpace),intent(in) :: this
!
!   EwaldSumKSpace_calc_coulomb_energy = EwaldSumKSpace_calc_energy_component(this, this % beta, this % rho_squared_total % rho_squared )
!
!end function 
!
!pure function EwaldSumKSpace_calc_LJ6_energy(this)
!   real(8) :: EwaldSumKSpace_calc_LJ6_energy
!
!   Type(TEwaldSumKSpace),intent(in) :: this
!
!   EwaldSumKSpace_calc_LJ6_energy = EwaldSumKSpace_calc_energy_component(this, this % beta6, this % rho_squared_total % F2_LJ6 )
!
!end function 
!
!pure function EwaldSumKSpace_calc_LJ12_energy(this)
!   real(8) :: EwaldSumKSpace_calc_LJ12_energy
!
!   Type(TEwaldSumKSpace),intent(in) :: this
!
!   EwaldSumKSpace_calc_LJ12_energy = EwaldSumKSpace_calc_energy_component(this, this % beta12, this % rho_squared_total % F2_LJ12 )
!
!end function 
!
!pure function EwaldSumKSpace_calc_LJ_energy(this)
!   real(8) :: EwaldSumKSpace_calc_LJ_energy
!
!   Type(TEwaldSumKSpace),intent(in) :: this
!
!   EwaldSumKSpace_calc_LJ_energy = EwaldSumKSpace_calc_LJ12_energy(this) - EwaldSumKSpace_calc_LJ6_energy(this)
!
!end function 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

EndModule EwaldSumKSpace
