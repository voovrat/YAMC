Module RhoSquared !DOC 
!DOC !FILE the module contains the coordinate-dependent sums of sin and cos for the Ewald sumation in KSpace
use SumSinCosKR
use FourierGrid
use LJTypes

implicit none

Type :: TRhoSquared  !DOC
! everything which is coordinate dependent
! normally : U_fourier = SUM_k beta(k) rho^2(k) , where rho^2(k) is coordinate dependent
!           in LJ case instrad of rho^2 we have F2, but it is the same
!
!  why to introduce this class: 
!    in movemc we need to have particular sums for one molecule --> the global object which total sums is not suitable
!DOC Fields:
   Type(TLJTypes),pointer :: lj_types !DOC LJTypes 
   Type(TFourierGrid),pointer :: grid !DOC KSpace grid

   Type(TSumSinCosKR) :: sumsincos_coulomb !DOC SUM_i sin(kR_i), SUM_i cos(kR_i)
   real(8),dimension(:),allocatable :: rho_squared   !DOC rho^2 = rho(+++)^2 + sgn(ky)rho(+-+)^2 + sgn(kz)rho(++-)^2 + sgn(kykz)rho(+--)^2

   Type(TSumSinCosKR),dimension(:),allocatable :: sumsincos_LJ !DOC SUM_i sin(kR_i) for each type of atoms

    ! real(8),dimension(:),allocatable :: F2_LJ6,F2_LJ12  
    ! NOTE:  Calculation of F2_LJ6, F2_LJ12 is impractical.
    ! for calculation the energy it is better to use the formula:
    !
    ! U_LJ = sum_k>=0 sum_t1 sum_t2>t1 (2-delta_t1t2) beta_LJ^t1t2(k) 
    !                   (     ( C^t1+++ C^t2+++ + S^t1+++ S^t2+++) +
    !               + sgn(ky) ( C^t1+-+ C^t2+-+ + S^t1+-+ S^t2+-+) 
    !               + sgn(kz) ( C^t1++- C^t2++- + S^t1++- S^t2++-)
    !            + sgn(ky*kz) ( C^t1+-- C^t2+-- + S^t1+-- S^t2+--) )
    !
    !   where beta_LJ^12 = ( 4 epsilon sigma^12 beta_12(k) - 4epsilon sigma^6 beta_6(k) )
    !       sigma = sigma_12 = (sigma1+sigma2)/2, epsilon= epsilon12= sqrt(epsilon1*epsilon2)
    !

   logical,dimension(:),allocatable :: type_is_present !DOC to accelerate calc F2:  if not type_is_present --> no summation needed

   ! temporary arrays
   Type(TSinCosXYZ) :: sincos_xyz  !DOC temporary array containing the sin(kR) cos(kR)
!   Type(TSinCosKR) :: sincos_kr_coulomb, sincos_kr_LJ

End Type TRhoSquared

contains 

subroutine RhoSquared_alloc(this,grid,lj_types) !DOC
!DOC Allocate the RhoSquare arrays
!DOC Parameters:
   Type(TRhoSquared) :: this !DOC RhoSquare structure

   Type(TFourierGrid),target :: grid !DOC KSpace grid
   Type(TLJTypes),target :: lj_types !DOC LJTypes

   integer :: ntype 
   integer :: i
   integer :: nk

   


   this % lj_types => lj_types
   this % grid => grid

   nk = grid % nk

    call SinCosXYZ_alloc( this % sincos_xyz , grid % kmax )  ! kmax (sic!), not nk

    call SumSinCosKR_alloc( this % sumsincos_coulomb,  grid )
!   call SinCosKR_alloc( this % sincos_kr_coulomb, grid )   
!   call SinCosKR_alloc( this % sincos_kr_LJ, grid )   
  

   ntype =  lj_types % NType

   ! fill sumsincos_LJ table for t1t2 
   allocate( this % sumsincos_LJ( ntype ) )
 
   do i=1,ntype

      call SumSinCosKR_alloc( this % sumsincos_LJ(i),  grid )

   end do

   allocate( this % rho_squared(nk) )
   allocate( this % type_is_present(ntype) )


 
end subroutine


subroutine RhoSquared_dealloc(this) !DOC 
!DOC deallocate the RhoSquared structure
!DOC Parameters:
   use io, only : write_real_array
   Type(TRhoSquared),target :: this !DOC RhoSquared structure

   integer :: ntype
   integer :: i

   ntype = this % lj_types % NType

   nullify( this % lj_types )
   nullify( this % grid )

   call SinCosXYZ_dealloc( this % sincos_xyz )
   call SumSinCosKR_dealloc( this % sumsincos_coulomb )

   do i=1,ntype
      
      call SumSinCosKR_dealloc( this % sumsincos_LJ(i) )

   end do

   deallocate( this % sumsincos_LJ )
   deallocate( this % rho_squared )
   deallocate( this % type_is_present )

end subroutine 


subroutine RhoSquared_nulify( this ) !DOC
!DOC Set RhoSquared to the "zero" state
!DOC Parameters:
   Type(TRhoSquared) :: this !DOC RhoSquared
   
   integer :: i,ntype

    ntype = this % lj_types % NType

   call SumSinCosKR_nulify( this % sumsincos_coulomb )

   do i=1,ntype
      call SumSinCosKR_nulify( this % sumsincos_LJ(i) )
   end do

   this % type_is_present(:) = .FALSE.

end subroutine

!
! ! U_fourier^C = SUM_mx,my,mz>=0 (1 + sgn(mx)) beta(k_m) 
!                               ( rho(+++)^2 + sgn(my)rho(+-+)^2 + sgn(mz)rho(++-)^2 + sgn(my*mz)rho(+--)^2)
subroutine RhoSquared_addAtom(this,xx,yy,zz,charge,lj_type,sincos_kr_coulomb,sincos_kr_LJ) !DOC 
!DOC Add new atom to the sums
!DOC Parameters:
   Type(TRhoSquared) :: this !DOC RhoSquared structure
   real(8),intent(in) :: xx,yy,zz,charge !DOC coordinates and charge
   integer,intent(in) :: lj_type !DOC type of atom

   Type(TSinCosKR) :: sincos_kr_coulomb, sincos_kr_LJ !DOC output: sin(kRi),cos(kri) for coulomb and LJ case

!   Type(TSinCosXYZ) :: sincos_xyz
   Type(TSumSinCosKR),pointer :: sumsincos  
   integer :: i

   if( charge /= 0.0 ) then

      call SinCosXYZ_fill( this % sincos_xyz,xx,yy,zz, this % grid % kmax)
      call SinCosKR_fill( sincos_kr_coulomb, this % sincos_xyz,  charge )
      call SumSinCosKR_addAtom(this %  sumsincos_coulomb, sincos_kr_coulomb )

   else
      call SinCosKR_nulify( sincos_kr_coulomb )
   end if

   if( lj_type > 0 ) then
!      call SinCosKR_fill( this % sincos_kr_LJ, this % sincos_xyz, 1)

      if ( charge /= 0.0 ) then
         call SinCosKR_copy( sincos_kr_LJ,  sincos_kr_coulomb)! copy(dst,src) 
         call SinCosKR_mulScalar( sincos_kr_LJ, 1.0/charge)  ! sincos_kr_LJ = sincos_kr_coulomb / charge
      else
         call SinCosXYZ_fill( this % sincos_xyz,xx,yy,zz,this % grid % kmax )
         call SinCosKR_fill( sincos_kr_LJ, this % sincos_xyz, 1.0d0)
      end if

      call SumSinCosKR_addAtom( this % sumsincos_LJ(lj_type), sincos_kr_LJ )
      this % type_is_present( lj_type ) = .TRUE.
   end if 

end subroutine


subroutine RhoSquared_addAtoms(this,xx,yy,zz,charge,type_by_index,natom) !DOC
!DOC Add many atoms
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use io,only : io_open,io_close,write_xyz_array
!DOC Parameters:
   Type(TRhoSquared),target:: this !DOC RhoSquared
   real(8),dimension(:),intent(in) :: xx,yy,zz,charge !DOC coordinates and charges
   integer,dimension(:),intent(in) :: type_by_index !DOC type by atom number
   integer,intent(in) :: natom !DOC number of atoms

   Type(TSinCosKR) :: sincos_kr_coulomb,sincos_kr_LJ

   real(8),dimension(:),pointer :: pa,pb,pc

   integer :: i

   call SinCosKR_alloc( sincos_kr_coulomb, this % grid )
   call SinCosKR_alloc( sincos_kr_LJ, this % grid )


   do i=1,natom
!   do i=3,3
  
      call RhoSquared_addAtom(this,xx(i),yy(i),zz(i),charge(i),type_by_index(i),sincos_kr_coulomb,sincos_kr_LJ )


!         write(tmpstr,'(AI3.3A)') 'cosk',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!
!         pa => this % sincos_xyz % cos_xkx
!         pb => this % sincos_xyz % cos_yky
!         pc => this % sincos_xyz % cos_zkz 
!
!         call write_xyz_array(hfile,pa,pb,pc,this % grid % kmax)
!         call io_close(hfile)
!
!         write(tmpstr,'(AI3.3A)') 'sink',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!
!         pa => this % sincos_xyz % sin_xkx
!         pb => this % sincos_xyz % sin_yky
!         pc => this % sincos_xyz % sin_zkz 
!
!         call write_xyz_array(hfile,pa,pb,pc, this % grid % kmax )
!         call io_close(hfile)

   !       write(tmpstr,'(AI3.3A)') 'cos_coulomb_ppp_',i,'.txt'
    !      call write_real_array


   end do

   call SinCosKR_dealloc( sincos_kr_coulomb )
   call SinCosKR_dealloc( sincos_kr_LJ )

end subroutine

subroutine RhoSquared_calc_rho_squared(this) !DOC 
!DOC calculate rho^2
!DOC note: this % sumsincos_coulomb should be calculated before (see calc_sumsincos_coulomb)
!DOC not: this % rho_squared should be allocated before
!
!DOC rho^2 = rho(+++)^2 + sgn(ky) rho(+-+)^2 + sgn(kz) rho(++-)^2 + sgn(kykz) rho(+--)^2
!
!DOC Parameters:
   Type(TRhoSquared),target :: this !DOC RhoSquared structure

   real(8),dimension(:),pointer :: rho_squared
   real(8),dimension(:),pointer :: C_ppp,C_ppm,C_pmp,C_pmm  ! C(+++),C(++-),C(+-+),C(+--) 
   real(8),dimension(:),pointer :: S_ppp,S_ppm,S_pmp,S_pmm

   integer :: nk,kk
   integer,dimension(:),pointer :: ky,kz
   logical :: sgn_ky,sgn_kz       

   ky => this % grid % ky
   kz => this % grid % kz

   rho_squared => this % rho_squared

   C_ppp => this % sumsincos_coulomb % sumcos_ppp
   C_ppm => this % sumsincos_coulomb % sumcos_ppm
   C_pmp => this % sumsincos_coulomb % sumcos_pmp
   C_pmm => this % sumsincos_coulomb % sumcos_pmm
 
   S_ppp => this % sumsincos_coulomb % sumsin_ppp
   S_ppm => this % sumsincos_coulomb % sumsin_ppm
   S_pmp => this % sumsincos_coulomb % sumsin_pmp
   S_pmm => this % sumsincos_coulomb % sumsin_pmm


   nk = this % grid % nk

   do kk=1,nk

      rho_squared(kk) = C_ppp(kk)**2 + S_ppp(kk)**2
       
      sgn_ky = ( ky(kk) > 0)
      sgn_kz = ( kz(kk) > 0)

      if( sgn_ky ) rho_squared(kk) = rho_squared(kk) + C_pmp(kk)**2 + S_pmp(kk)**2
      if( sgn_kz ) rho_squared(kk) = rho_squared(kk) + C_ppm(kk)**2 + S_ppm(kk)**2
      if( sgn_ky .and. sgn_kz) rho_squared(kk) = rho_squared(kk) + C_pmm(kk)**2 + S_pmm(kk)**2
      
   end do
   
end subroutine

subroutine RhoSquared_copy( dst, src ) !DOC
!DOC Copy data in the RhoSquared structure 
!DOC Parameters:
   Type(TRhoSquared) :: dst,src !DOC destinatioion, source
   integer :: i

   call SumSinCosKR_copy( dst % sumsincos_coulomb , src % sumsincos_coulomb )
   do i=1,dst % lj_types % NType

      call SumSinCosKR_copy( dst % sumsincos_LJ(i) , src % sumsincos_LJ(i) )

   end do

end subroutine 

subroutine RhoSquared_add( this, rho_squared ) !DOC
!DOC Add two RhoSquared structures, i.e. this=this+rho_squared
!DOC Parameters:
    Type(TRhoSquared) :: this, rho_squared !DOC summands
    integer :: i

    call SumSinCosKR_addSum( this % sumsincos_coulomb, rho_squared % sumsincos_coulomb )
   
    do i=1,this % lj_types % NType
 
       call SumSinCosKR_addSum( this % sumsincos_LJ(i), rho_squared % sumsincos_LJ(i) )

    end do

end subroutine

subroutine RhoSquared_sub( this, rho_squared )!DOC
!DOC substract: this=this-rho_squared \
!DOC Parameters:
    Type(TRhoSquared) :: this, rho_squared !DOC operands
    integer :: i

    call SumSinCosKR_subSum( this % sumsincos_coulomb, rho_squared % sumsincos_coulomb )
   
    do i=1,this % lj_types % NType
 
       call SumSinCosKR_subSum( this % sumsincos_LJ(i), rho_squared % sumsincos_LJ(i) )

    end do

end subroutine

subroutine RhoSquared_save_to_file(this,filename) !DOC 
!DOC Save to file
   use io, only : io_open,io_close
   use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC For debug only
!DOC  format columns: 
!DOC   coulomb cos_ppp,cos_pmp,... sin_ppp,sin_pmp,... 
!DOC  LJ1 cos_ppp,cos_pmp,... sin_ppp,sin_pmp,...
!DOC  LJ2 com_ppp,cos_pmp,... 
!DOC  ...
!DOC Parameters:
   Type(TRhoSquared),intent(in),target :: this !DOC RhoSquared structure
   character(SYSTEM_STRING_LENGTH),intent(in) :: filename !DOC name of the fle
   
   Type(TSumSinCosKR),pointer :: ssc 
   integer :: hfile
   integer :: i,k
   

   hfile = io_open(filename,'w')
   
   do k=1,this % grid % nk

       ssc => this % sumsincos_coulomb

       write(hfile,'(F20.10,F20.10,F20.10,F20.10,$)') ssc % sumcos_ppp(k), ssc % sumcos_pmp(k),&
                                                      ssc % sumcos_ppm(k), ssc % sumcos_pmm(k)

       write(hfile,'(F20.10,F20.10,F20.10,F20.10,$)') ssc % sumsin_ppp(k), ssc % sumsin_pmp(k),&
                                                      ssc % sumsin_ppm(k), ssc % sumsin_pmm(k)

       do i=1,this % lj_types % NType 

          ssc => this % sumsincos_LJ(i)
          write(hfile,'(F20.10,F20.10,F20.10,F20.10,$)') ssc % sumcos_ppp(k), ssc % sumcos_pmp(k),&
                                                         ssc % sumcos_ppm(k), ssc % sumcos_pmm(k)

          write(hfile,'(F20.10,F20.10,F20.10,F20.10,$)') ssc % sumsin_ppp(k), ssc % sumsin_pmp(k),&
                                                         ssc % sumsin_ppm(k), ssc % sumsin_pmm(k)
       end do ! i 

       write(hfile,*) 

   end do ! k

   call io_close(hfile)
end subroutine

!! U_fourier^LJ = SUM_mx,my,mz>=0 (1 + sgn(mx)) (beta6(k_m) + beta12(k_m)) 
!!                               ( F2(+++) + sgn(my)F2(+-+) + sgn(mz)F2(++-) + sgn(my*mz)F2(+--) )
!! 
!! F2(+++) = SUM_t1t2 A_t1t2 ( C^t1(+++) C^t2(+++) + S^t1(+++) S^t2(+++) ) 
!! F2(++-) = SUM_t1t2 A_t1t2 ( C^t1(++-) C^t2(++-) + S^t1(++-) S^t2(++-) ) 
!! F2(+-+) = SUM_t1t2 A_t1t2 ( C^t1(+-+) C^t2(+-+) + S^t1(+-+) S^t2(+-+) ) 
!! F2(+--) = SUM_t1t2 A_t1t2 ( C^t1(+--) C^t2(+--) + S^t1(+--) S^t2(+--) ) 
!!                                
!!  C^t(sxsysz), S^t(sxsysz) calculated by SumSinCosKR
!!
!subroutine RhoSquared_calcSumSinCosLJ(this,xx,yy,zz,type_by_index,natom)
!!
!!  for the first time  first_atom = 1, atom_count = atom_data % natom 
!!  but for movemc we need to calculate partial sums for one molecule only
!!  for that purpose we use first_index and n_index
!!
!   Type(TRhoSquared) :: this
!   integer,dimension(:),intent(in) :: type_by_index ! types of the atoms 
!   real(8),dimension(:),intent(in) :: xx,yy,zz ! coordinates of the  atoms
!    integer :: natom ! number of atoms 
!
!   integer :: ityp, iatom   ! iterators for types and atoms inside one type
!   integer :: atom_index ! atom_index read from the index array
!   integer :: ntype 
!   
!
!!   Type(TLJTypes),pointer :: lj_types
!!   integer,dimension(:),pointer :: idx_array
!   logical,dimension(:),pointer :: type_is_present   
!
!   Type(TSumSinCosKR),pointer :: sumsincos
!
!
!   type_is_present => this % type_is_present
!   ntype = this % lj_types % NType
!  
!   ! fill sumsincos_LJ table for t1t2 
!   ! sumsincos_LJ should be already allocated 
!
!   ! set sums to zero
!   do ityp=1,ntype
!
!      call SumSinCosKR_nufify( this % sumsincos_LJ(ityp) )
!
!      type_is_present(ityp) = .FALSE.
!   end do
!
!  type_is_present(:) = .FALSE.
!
!
!   do iatom=1,natom
! 
!      ityp = type_by_index(iatom)
!      if ( ityp == 0 ) cycle   
!
!      type_is_present( ityp ) = .TRUE.
!
!      SumSinCorKR_addAtom( this % sumsincos_LJ(ityp), xx(iatom), yy(iatom), zz(iatom) )
!     
!   end do ! iatom
!
!end subroutine
!
    ! NOTE:  Calculation of F2_LJ6, F2_LJ12 is impractical.
    ! for calculation the energy it is better to use the formula:
    !
    ! U_LJ = sum_k>=0 sum_t1 sum_t2>t1 (2-delta_t1t2) beta_LJ^t1t2(k) 
    !                   (     ( C^t1+++ C^t2+++ + S^t1+++ S^t2+++) 
    !               + sgn(ky) ( C^t1+-+ C^t2+-+ + S^t1+-+ S^t2+-+) 
    !               + sgn(kz) ( C^t1++- C^t2++- + S^t1++- S^t2++-)
    !            + sgn(ky*kz) ( C^t1+-- C^t2+-- + S^t1+-- S^t2+--) )
    !
    !   where beta_LJ^12 = ( 4 epsilon sigma^12 beta_12(k) - 4epsilon sigma^6 beta_6(k) )
    !       sigma = sigma_12 = (sigma1+sigma2)/2, epsilon= epsilon12= sqrt(epsilon1*epsilon2)
    !


! *********************************************************************************************
! ------------> Thus use the  code below for educational purposes only (it is slow) 
!**********************************************************************************************

!
!   ! calculate F2(sxsysz) = SUM_t1t2 ( C^t1 C^t2 + S^t1
!subroutine RhoSquared_calc_F2_LJcomponent(this,lj_tab,F2)
!
!! note: this % sumsincos_LJ should be calculated before (see calcSuSinCosLJ)
!! not: this % rho_squared should be allocated before
!!
!! F2(k) =  F2(+++) + sgn(ky) F2(+-+) + sgn(kz) F2(++-) + sgn(kykz) F2(+--)
!!
!! where F2(sxsysz) = SUM_t1t2 A_t1t2 ( C^t1(sxsysz) C^t2(sxsysz) + S^t1(sxsysz) S^t2(sxsysz) )
!!
!!  Actually, we calculate
!! F2(k) = SUM_t1t2 A_t1t2 (         ( C^t1(+++) C^t2(+++) + S^t1(+++) S^t2(+++) )
!!                           +sgn(ky)( C^t1(+-+) C^t2(+-+) + S^t1(+-+) S^t2(+-+) )
!!                           +sgn(kz)( C^t1(++-) C^t2(++-) + S^t1(++-) S^t2(++-) ) 
!!                         +sgn(kykz)( C^t1(+--) C^t2(+--) + S^t1(+--) S^t2(+--) )
!!                          )
!!
!!  A_t1t2 = A_t1t2_LJ6  or  A_t1t2_LJ12 
!!  the type is defined by LJTab array which is either lj_types % LJ6Tab or lj_types % LJ12Tab
!
!   Type(TRhoSquared) :: this
!   real(8),dimension(:),intent(in) :: lj_tab ! either lj_types % LJ6Tab or lj_types % LJ12Tab
!   real(8),dimension(:),intent(out) :: F2  ! either this % F2_LJ6 or this % F2_LJ12
!
!   integer :: t1,t2
!   integer :: ntyp
!
!   real(8),dimension(:),pointer :: F2
!   real(8) :: Sum_kk  ! temporary  sum_kk = C^t1(+++) C^t2(+++) + S^t1(+++) S^t2(+++) + sgn(ky) * ...... (see above) inside brackets 
!   real(8) :: A_t1t2 ! coefficients for 1/r6 and 1/r^12 potentials
!
!   real(8),dimension(:),pointer :: C1_ppp,C1_ppm,C1_pmp,C1_pmm  ! C^t1(+++),C^t1(++-),C^t1(+-+),C^t1(+--) 
!   real(8),dimension(:),pointer :: C2_ppp,C2_ppm,C2_pmp,C2_pmm  ! C^t2(+++),C^t2(++-),C^t2(+-+),C^t2(+--) 
!
!   real(8),dimension(:),pointer :: S1_ppp,S1_ppm,S1_pmp,S1_pmm  ! 
!   real(8),dimension(:),pointer :: S2_ppp,S2_ppm,S2_pmp,S2_pmm  ! 
!
!   integer :: nk,kk
!   integer,dimension(:),pointer :: ky,kz
!   logical :: sgn_ky,sgn_kz       
!   
!   logical,dimension(:),pointer :: type_is_present
!
!   ky => this % grid % ky
!   kz => this % grid % kz
! 
!   type_is_present => this % type_is_present
!
!   ntyp = this % lj_types % NTypes
!
!   F2(:) = 0
!
!   do t1 = 1,ntyp
!
!   if( .not.type_is_present(t1) ) cycle
!
!   C1_ppp => this % sumsincos_LJ(t1) % sumcos_ppp
!   C1_ppm => this % sumsincos_LJ(t1) % sumcos_ppm
!   C1_pmp => this % sumsincos_LJ(t1) % sumcos_pmp
!   C1_pmm => this % sumsincos_LJ(t1) % sumcos_pmm
! 
!   S1_ppp => this % sumsincos_LJ(t1) % sumsin_ppp
!   S1_ppm => this % sumsincos_LJ(t1) % sumsin_ppm
!   S1_pmp => this % sumsincos_LJ(t1) % sumsin_pmp
!   S1_pmm => this % sumsincos_LJ(t1) % sumsin_pmm
! 
!
!   do t2 = t1,ntyp ! ok, sum is symmetric, because A_12 = A_21
!  
!      if( .not. type_is_present(t2) ) cycle
!
!      C2_ppp => this % sumsincos_LJ(t2) % sumcos_ppp
!      C2_ppm => this % sumsincos_LJ(t2) % sumcos_ppm
!      C2_pmp => this % sumsincos_LJ(t2) % sumcos_pmp
!      C2_pmm => this % sumsincos_LJ(t2) % sumcos_pmm
!    
!      S2_ppp => this % sumsincos_LJ(t2) % sumsin_ppp
!      S2_ppm => this % sumsincos_LJ(t2) % sumsin_ppm
!      S2_pmp => this % sumsincos_LJ(t2) % sumsin_pmp
!      S2_pmm => this % sumsincos_LJ(t2) % sumsin_pmm
!   
!      nk = this % grid % nk
!   
!
!      A_t1t2 = LJTypes_Koef(this % lj_types, t1, t2, lj_tab) 
!
!      do kk=1,nk
!
!         Sum_kk = 0   
!         Sum_kk = C1_ppp(kk) * C2_ppp(kk) + S1_ppp(kk) * S2_ppp(kk)
!          
!         sgn_ky = ( ky > 0)
!         sgn_kz = ( kz > 0)
!   
!         if( sgn_ky ) Sum_kk = Sum_kk + C1_pmp(kk) * C2_pmp(kk) + S1_pmp(kk) * S2_pmp(kk)
!         if( sgn_kz ) Sum_kk = Sum_kk + C1_ppm(kk) * C2_ppm(kk) + S1_ppm(kk) * S2_ppm(kk)
!         if( sgn_ky .and. sgn_kz ) Sum_kk = Sum_kk + C1_pmm(kk) * C2_pmm(kk) + S1_pmm(kk) * S2_pmm(kk)
! 
!         if ( t1 /= t2) Sum_kk = 2*Sum_kk  ! to account A_t2t1 = A_t1t2
! 
!         F2(kk) = F2(kk) +  A_t1t2 * Sum_kk
!    
!      end do ! kk
!
!
!   end do ! t2
!   end do ! t1
!   
!end subroutine
!
!subroutine RhoSquared_calc_F2LJ6(this)  ! osolete, do not use it. 
!   Type(TRhoSquared) :: this
!
!   call RhoSquared_calc_F2_LJcomponent(this, this % lj_types % LJ6Tab, this % F2_LJ6 )
!
!end subroutine
!
!subroutine RhoSquared_calc_F2LJ12(this) ! do not use this function. see the explanations above
!   Type(TRhoSquared) :: this
!
!   call RhoSquared_calc_F2_LJcomponent(this, this % lj_types % LJ12Tab, this % F2_LJ12 )
!
!end subroutine
!
!subroutine RhoSquared_calc_F2(this) ! do not use this function
!   Type(TRhoSquared) :: this
!
!   call RhoSquared_calc_F2LJ12(this)
!   call RhoSquared_calc_F2LJ6(this)
!  
!end subroutine






End Module
