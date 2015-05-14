Module LJTypes !DOC
!DOC !FILE Lennard Jones parameters for the pairs of atom types
! 
! In k-space ewald summantion we have a term F2 = SUM_sj A_sj exp( i k * (r_s-r_j) )
! For charges A_ij = q_i q_j, and the sum is F2 = SUM_sj A_sj exp( i k (r_s - r_j) ) = (SUM_s q_s exp(ik r_s))^2 = F^2
! where F requires only Natom summations
! This is not the case for the general A_ij, (in particular - for LJ) 
! But if we have only M atom types the sum can be rewritten
! F2 = SUM_t1 SUM_t2 A_t1t2  SUM_s=1^N(t1) SUM_j=1^(N(t2)) exp(i k (r_s^t1 - r_j^t2)) =
! = SUM_t1 SUM_t2 ( SUM_s exp(ik r_s^t1 ) ) * conj( SUM_j exp(i k r_j^t2) ) = SUM_t1 SUM_t2 F_t1 * F_t2
! where calculation of F_1... F_tM also require only N summations 
! 
! To do this, we need to build indeces for molecules which have the same A_sj
! This module is about that
! 

!const
  implicit none

  integer, parameter :: MAX_LJ_TYPE = 20  ! maximum number of types supported

! types
   Type :: LJIndex !DOC
!DOC the structure to keep the indeces of the atoms which have specific atom type
!DOC Fileds:
      real(8) :: sigma,epsilon !DOC sigma and epsilon parameters
      integer,dimension(:),allocatable :: idx !DOC array of indeces
      integer :: n,nalloc !DOC number of indeces and allocated size of array
   End Type

   Type :: TLJTypes  !DOC
  !DOC Fields:
      Type(LJIndex),dimension(MAX_LJ_TYPE) :: index_by_type !DOC get LJIndex for specific atom type
      integer,dimension(:),allocatable  :: type_by_index !DOC get type by atom index
      integer :: NType !DOC number of atom types
   
      ! usage: LJTypes_LJ6Tab(i) % array(j)  -- LJ6(i,j) coefficient
      real(8),dimension(:),allocatable :: LJ6Tab  !DOC table of LJ coefficients near  1/r^6, e.g. 4epsilon_ij sigma_ij^6
      real(8),dimension(:),allocatable :: LJ12Tab !DOC table of LJ coefficient near 1/r^12

      real(8),dimension(:),allocatable :: four_epsilon_tab  !DOC table for 4*epsilon_ij
      real(8),dimension(:),allocatable :: sigma6_tab !DOC table for sigma_ij^6
      real(8),dimension(:),allocatable :: sigma2_tab  !DOC table for sigma_ij^2. Used for checking the overlap 
      ! HOW TO USE THESE TABLES : 
      !  base_offset = 0
      !  do t1 = 1,ntype
      !    do t2 = t1,ntype ! or t1=1,ntype,  or t1=whatever... 
      !
      !       offset = base_offset + t2
      !       sigma6 = lj_types % sigma6_tab( offset ) 
      !       four_epsilon = lj_types % four_epsilon( offset )
      !       A_LJ6 = lj_types % LJ6Tab( offset )
      !       A_LJ12 = lj_types % LJ12Tab( offset )
      !
      !           <.. ( do something ) ..>
      ! 
      !    end do
      !
      !    base_offset = base_offset + ntype
      !
      !  end do

   End Type TLJTypes

contains

subroutine LJIndex_alloc(this,nalloc,sigma,epsilon) !DOC
!DOC allocate the LJIndex
!DOC Parameters:
   Type(LJIndex) :: this !DOC LJIndex
   integer :: nalloc !DOC size to be allocated
   real(8) :: sigma,epsilon !DOC values of sigma and epsilon
   

   allocate(this % idx(nalloc))
   this % n = 0
   this % sigma = sigma
   this % epsilon = epsilon
   this % nalloc = nalloc

end subroutine

subroutine LJIndex_dealloc(this) !DOC 
!DOC Dealloate the LJIndex
   Type(LJIndex) :: this !DOC LJIndex

   deallocate( this % idx )
   this % nalloc = 0

end subroutine

subroutine LJIndex_addIndex(this,i) !DOC
!DOC add index to the LJIndex array
   use error 
!DOC Parameters:
   Type(LJIndex) :: this  !DOC LJIndex
   integer :: i !DOC new index

   if( this % n + 1 .gt. this % nalloc ) then
       write(error_message,*) 'LJType_addIndex: Cannot add idx ',i,': allocated size exceeded'
       call error_throw(ERROR_LIMITS)
       return
   end if   

   this % n = this % n + 1
   this % idx(this % n) = i

end subroutine

subroutine LJTypes_dealloc(this) !DOC
!DOC deallocate the LJTypes array
!DOC Parameters:
   Type(TLJTypes) :: this !DOC LJTypes structure

   integer :: i
   integer :: ntype,ntype2

   ntype = this % NType

   do i=1,ntype
       call LJIndex_dealloc(this % index_by_type(i))
   end do

   deallocate(this % type_by_index )
   deallocate(this % LJ6Tab)
   deallocate(this % LJ12Tab)
   deallocate(this % four_epsilon_tab)
   deallocate(this % sigma6_tab)
   deallocate(this % sigma2_tab)

end subroutine


subroutine LJTypes_fill(this,sigma_array,epsilon_array,nsigma) !DOC
!DOC Allocate and fill the LJType arrays
!DOC Parameters:
   Type(TLJTypes) :: this !DOC LJTypes structure
   real(8), dimension(:) :: sigma_array, epsilon_array !DOC arrays for simga and epsilon for each atom
   integer :: nsigma !DOC size of the input arrays

   integer :: i,typ,Ntyp
   real(8),dimension(MAX_LJ_TYPE) :: sigma_types,epsilon_types
   integer,dimension(MAX_LJ_TYPE) :: sigma_count
   logical :: found

! count number of different types (to be able to allocate the arrays)
   Ntyp = 0
   do i=1,nsigma

      if ( sigma_array(i) .lt. 1e-8 ) cycle

      found = .FALSE.
      do typ=1,Ntyp
         if ( (dabs( sigma_array(i) - sigma_types(typ)) < 1d-8) .and. &
              (dabs( epsilon_array(i) - epsilon_types(typ)) < 1d-8) ) then
            sigma_count(typ) = sigma_count(typ) + 1
            found = .TRUE.
            exit
         end if 
      end do

      if(.not.found) then
         Ntyp = Ntyp + 1
         sigma_types(Ntyp) = sigma_array(i)
         epsilon_types(Ntyp) = epsilon_array(i)
         sigma_count(Ntyp) = 1
      end if 

   end do
! allocate the indices

   do typ = 1,Ntyp
      call LJIndex_alloc(this % index_by_type(typ),sigma_count(typ),sigma_types(typ),epsilon_types(typ) )
   end do
! allocate type_by_index array

  allocate(this % type_by_index(nsigma))

  this % type_by_index(:) = 0

! fill the indices
   do i=1,nsigma

      do typ = 1,Ntyp
         if( ( dabs(sigma_array(i) - sigma_types(typ)) < 1d-8) .and. &
             ( dabs(epsilon_array(i) - epsilon_types(typ)) < 1d-8 ) ) then
             call LJIndex_addIndex(this % index_by_type(typ),i)
             this % type_by_index(i) = typ
             exit
         end if
      end do
   
   end do   

   this % NType = Ntyp

   call LJTypes_allocate_LJ_tab(this)
   call LJTypes_fill_LJ_tab(this)

end subroutine

subroutine LJTypes_fill_LJ_tab_index(this,i,j,sigma12,epsilon12)   !DOC
!DOC auxilarly subroutine to fill LJ6 LJ12, 4epsilon and other arrays
!DOC Parameters:
   Type(TLJTypes) :: this !DOC LJTypes
   integer,intent(in) :: i,j !DOC pair of indeces
   real(8), intent(in) :: sigma12,epsilon12 !DOC sigma and epsilon for this pair
     
   real(8) :: sigma_power_2,sigma_power_6 ! sigma^6 (why not sigma6: sigma12 can be both: sigma_{12} and sigma^12 
   integer :: offset12,offset21

   offset12 = (i-1) * this % NType + j

   sigma_power_2 = sigma12**2
   sigma_power_6 = sigma_power_2**3

   this % LJ6Tab( offset12 ) = 4*epsilon12*sigma_power_6
   this % LJ12Tab( offset12 ) = this % LJ6Tab( offset12 ) * sigma_power_6  ! 4 epsilon12 * sigma^6 * sigma^6 = 4 epsilon12 * sigma^12   

   this % sigma2_tab( offset12 ) = sigma_power_2
   this % sigma6_tab( offset12 ) = sigma_power_6
   this % four_epsilon_tab( offset12 ) = 4 * epsilon12  


   if ( i==j ) return

   offset21 = (j-1) * this % NType + i

   this % LJ6Tab( offset21 ) = this % LJ6Tab( offset12 )
   this % LJ12Tab( offset21 ) = this % LJ12Tab( offset21 ) 
   this % sigma2_tab( offset21 ) = sigma_power_2
   this % sigma6_tab( offset21 ) = sigma_power_6
   this % four_epsilon_tab( offset21 ) = 4 * epsilon12  

end subroutine

! use for changes of the BoxLength
subroutine LJTypes_scale_sigma( this, scale_coeff ) 
!DOC scales all sigma and epsilon. Can be called to re-calculate the LJTypes without re-allocation (for example when the BoxLength changes)
!DOC Parameters:
   Type(TLJTypes) :: this !DOC LJTypes structure
   real(8),intent(in) :: scale_coeff !DOC scale coefficient

   integer :: i
      
   do i=1,this % NType
      this % index_by_type(i) % sigma = this % index_by_type(i) % sigma * scale_coeff
   end do

   call LJTypes_fill_LJ_tab(this) 

end subroutine

subroutine LJTypes_allocate_LJ_tab(this) !DOC
!DOC allocate the LJ Tabs
!DOC Parameters:
   implicit none
   Type(TLJTypes) :: this   !DOC LJTypes
   integer :: ntype,ntype2 

   ntype = this % NType
   ntype2 = ntype**2

   allocate( this % LJ6Tab( ntype2 ) )
   allocate( this % LJ12Tab( ntype2 ) ) 
   allocate( this % four_epsilon_tab( ntype2 ) )
   allocate( this % sigma6_tab( ntype2 ) )  
   allocate( this % sigma2_tab( ntype2 ) )

end subroutine


subroutine LJTypes_fill_LJ_tab(this) !DOC
!DOC use Lorentz-Berthelot rules to fill LJ6 and LJ12 tables 
!DOC Parameters:
   Type(TLJTypes) :: this   !DOC LJTypes structure
 
   integer :: i,j
   integer :: ntype,ntype2
   real(8) :: sigma1,epsilon1,sigma2,epsilon2,sigma12, epsilon12

   ntype = this % NType
   ntype2 = ntype**2

   do i=1,this % NType

      ! simga_ii = sigma_i, epsilon_ii = epsilon_i

      sigma1 = this % index_by_type(i) % sigma
      epsilon1 = this % index_by_type(i) % epsilon
  
      call LJTypes_fill_LJ_tab_index(this,i,i,sigma1,epsilon1)

      do j=i+1,this % NType
 
         sigma2 = this % index_by_type(j) % sigma
         epsilon2 = this % index_by_type(j) % epsilon 

         sigma12 = ( sigma1 + sigma2 ) / 2
         epsilon12 = dsqrt(epsilon1*epsilon2)

         call LJTypes_fill_LJ_tab_index(this,i,j,sigma12,epsilon12)
     
      end do
   end do
      

end subroutine

 ! it is cool to use the functions below, by they multiply all the time (t1-1) * ntype + t2
 ! because all the coefficients are used in the cycles
 ! it is much better to use them like that:
 ! 
 !  base_offset = 0
 !  do t1 = 1,ntype
 !    do t2 = t1,ntype
 !
 !       offset = base_offset + t2
 !       sigma6 = lj_types % sigma6_tab( offset ) 
 !       four_epsilon = lj_types % four_epsilon( offset )
 !       A_LJ6 = lj_types % LJ6Tab( offset )
 !       A_LJ12 = lj_types % LJ12Tab( offset )
 !
 !           <.. ( do something ) ..>
 ! 
 !    end do
 !
 !    base_offset = base_offset + ntype
 !
 !  end do

 ! it is actually much better to 

!pure function LJTypes_Koef(this,i,j,array)  ! works for both: LJ6 and LJ12
!!   use error
!   real(8) :: LJTypes_Koef
!
!   Type(TLJTypes),intent(in) :: this
!   integer,intent(in) :: i,j
!   real(8),dimension(:),intent(in) :: array
!   
!! say thanks to "pure" function... 
!!   if( (i>this % NType) .or. (j>this % NType) ) then
!!       write(error_message,*) 'LJTypes_Koef: element (',i,',',j,') does not exist'
!!       call error_throw(ERROR_LIMITS)
!!       return
!!   end if
!   
!   LJTypes_Koef = array( (i-1) * this % Ntype + j )
!
!end function
!
!pure function LJTypes_LJ6(this,i,j)
!   real(8) :: LJTypes_LJ6 
!
!   Type(TLJTypes),intent(in) :: this
!   integer,intent(in) :: i,j
!
!   LJTypes_LJ6 = LJTypes_Koef(this,i,j,this % LJ6Tab)
!
!end function
!
!pure function LJTypes_LJ12(this,i,j)
!   real(8) :: LJTypes_LJ12
!
!   Type(TLJTypes),intent(in) :: this
!   integer,intent(in) :: i,j
!
!   LJTypes_LJ12 = LJTypes_Koef(this,i,j,this % LJ12Tab)
!
!end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

End Module
