Module SumSinCosKR !DOC
use FourierGrid
!DOC !FILE SUM_i cos(kR_i), SUM_i sin(kR_i). Used in EwaldSums in KSpace

!
!  In ewald summation we can rewrite the fourier part sum in a form
!    U_fourier = SUM_k beta(k) F_2(k) 
! 
!  where beta(k) differs for 1/r, 1/r^6, 1/r^12 potentials
!  ( for 1/r it is 4pi/k_m exp(-k_m^2/4alpha^2) )
!  and 
!  F_2(k_m) = SUM_sj A_sj exp(i*k_m * ( r_s - r_j) )
!  (here i=sqrt(-1) )
!  
!  Now, in our notation  m==k <> k_m = 2pi m  / L
!              and L = 1
!
!  if A_sj = q_s q_j then F_2(k_m) = F(k_m) F(-k_m) 
!  where 
!         F(k_m) = SUM_j q_j exp ( 2 pi i  m / L * r_j )
!
!  For other potentials, we do separate summation over the different atom types
!
!  F_2 = SUM_t1 SUM_t2 A_t1t2 SUM_s^N1 SUM_j^N2 exp(2 pi m/L * (r_s^t1 - r_j^t2) = 
!  = SUM_t1 SUM_t2 A_t1t2  F^t1(k_m) * F^t2(-k_m)
!
!  where F_tp(k_m) = SUM_s^Np exp(2 pi i m/L * r_s^tp )
!
! which is the sum of the same kind as in coulomb ewald
!
! Our aim in this module is to calculate this sum
! To avoid complex numbers 
!
!  F^p(k_m) = C_p(k_m) + i*S_p(k_m)
!  where 
!
! C_p(k_m) = SUM_s^Np cos( 2 pi m/L * r_s ) 
! S_p(k_m) = SUM_s^Np sin( 2 pi m/L * r_s )
!
! F^1(k_m) * F^2(-k_m) = (C1 + i S1)(C2 - iS2) = C1C2 + S1S2 - i C1S2 + i S1C2
! 
! but normally beta(k) is symetric, i.e. beta(k) = beta(-k) = beta(|k|)
! 
! Thus, in the sum SUM_k_m F_2(k_m) beta(k_m) = 2 SUM_kx>0 (F_2(k_m) + F_2(-km)) beta(k_m)
!
! F_2 = SUM_t1t2 F^t1 F^t2
!
! for two types 1,2  the terms with i will cancel out
!
! So, in any case the only problem is to calculate SUM_s cos(k_m r_s) and SUM_s sin(k_m r_s)
!
! we can also use symmetry of beta again...
!
! in that case SUM_km = 2 SUM_(kx,ky,kz>0)  ( F_2(+++) + F2(++-) + F2(+-+) + F2(+--))
!
!  where F2(+-+) = SUM cos(xkx - yky + zkz) + i SUM sin( xkx - yky + zkz) 
! ( pluses and minuses in the brackets mean the signs of xkx, yky and zkz in the sum)
!
! so, we need to calculate 8 sums, 4 sin/cos pairs for +++,++-, +-+ and +--

implicit none

! auxilarly type used by SumSinCosKR
Type :: TSinCosXYZ !DOC
!DOC temporary arrays for calculation
!DOC the reason to create them: they can be once allocated and then re-used
!DOC Fields:
   integer :: nk  !DOC number of elements allocated

   real(8),dimension(:),allocatable :: cos_xkx,cos_yky,cos_zkz !DOC cos
   real(8),dimension(:),allocatable :: sin_xkx,sin_yky,sin_zkz !DOC sin

End Type TSinCosXYZ

Type :: TSinCosKR  !DOC
!DOC  sin,cos(sx*x*kx + sy*y*ky + sz*z*kz ) 
!DOC Fields:
   Type(TFourierGrid),pointer :: grid  !DOCKSpace grid 

   real(8),dimension(:),allocatable :: cos_ppp  !DOC  cos ( xkx + yky + zkz )  (+++)
   real(8),dimension(:),allocatable :: sin_ppp  !DOC  sin ( xkx + yky + zkz )  (+++)

   real(8),dimension(:),allocatable :: cos_ppm  !DOC  cos ( xkx + yky - zkz )  (++-)
   real(8),dimension(:),allocatable :: sin_ppm  !DOC  sin ( xkx + yky - zkz )  (++-)

   real(8),dimension(:),allocatable :: cos_pmp  !DOC  cos ( xkx - yky + zkz )  (+-+)
   real(8),dimension(:),allocatable :: sin_pmp  !DOC  sin ( xkx - yky + zkz )  (+-+)

   real(8),dimension(:),allocatable :: cos_pmm  !DOC  cos ( xkx - yky - zkz )  (+--)
   real(8),dimension(:),allocatable :: sin_pmm  !DOC  sin ( xkx - yky - zkz )  (+--)


End Type

! now: the procedure is the following:
!
! call SumSinCosKR_nulify(sumsincos)
!   
! do i=1,Natom
!   call SinCosXYZ_fill(sincos_xyz,x(i),y(i),z(i),kmax) ! calculate sin,cos(x*kx), sin,cos(y*ky), sin,cos(z*kz)
!   call SinCosKR_fill(sincos_kr, sincos_xyz, charge)   ! calculate sin,cos(sx x kx + sy y ky + sz z kz )
!   call SumSinCosKR_addAtom(sumsincos,sincos_kr)       ! add sin,cos... to the sum
! end do
!
! but there are alternative operations
! e.g SumSinCosKR_addSum
!     SumSinCosKR_subSum
!   

Type :: TSumSinCosKR !DOC
!DOC SUM sin(kR), SUM cos(kR)
!DOC Fields:
   Type(TFourierGrid),pointer :: grid    !DOC KSpace grid 

   real(8),dimension(:),allocatable :: sumcos_ppp  !DOC SUM_s cos ( xkx + yky + zkz )  (+++)
   real(8),dimension(:),allocatable :: sumsin_ppp  !DOC SUM_s sin ( xkx + yky + zkz )  (+++)

   real(8),dimension(:),allocatable :: sumcos_ppm  !DOC SUM_s cos ( xkx + yky - zkz )  (++-)
   real(8),dimension(:),allocatable :: sumsin_ppm  !DOC SUM_s sin ( xkx + yky - zkz )  (++-)

   real(8),dimension(:),allocatable :: sumcos_pmp  !DOC SUM_s cos ( xkx - yky + zkz )  (+-+)
   real(8),dimension(:),allocatable :: sumsin_pmp  !DOC SUM_s sin ( xkx - yky + zkz )  (+-+)

   real(8),dimension(:),allocatable :: sumcos_pmm  !DOC SUM_s cos ( xkx - yky - zkz )  (+--)
   real(8),dimension(:),allocatable :: sumsin_pmm  !DOC SUM_s sin ( xkx - yky - zkz )  (+--)

   ! temporary arrays
!   Type(TSinCosXYZ) :: sincos_xyz 

End Type TSumSinCosKR

! when calculating the sum: 
!  
!  if(kx>0)  then beta = 2*beta   (account for kx<0)
!  
!  som = sumcos_ppp(kk)^2 + sumsin_ppp(kk)^2
!  if( ky>0)  som +=  sumcos_pmp^2 + sumsin_pmp^2  
!  if( kz>0)  som +=  sumcom_ppm^2 + sumsin_ppm^2
!  if( ky>0 and kz>0 ) som += sumcos_pmm^2 + sumsin_pmm^2
!

contains

!!!!!!!!!!!!!!!!! SinCosXYZ

! *************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!! SinCosXYZ  !!!!!!!!!!!!!!!!!!!!!!!!!
! *************************************************************

subroutine SinCosXYZ_alloc(this,kmax) !DOC
!DOC allocate SinCosXYZ structure
!DOC Parameters:
   Type(TSinCosXYZ) :: this !DOC SinCosXYZ
   integer,intent(in) :: kmax !DOC max |k|

   allocate( this % cos_xkx(0:kmax) )
   allocate( this % cos_yky(0:kmax) )
   allocate( this % cos_zkz(0:kmax) )
   
   allocate( this % sin_xkx(0:kmax) )
   allocate( this % sin_yky(0:kmax) )
   allocate( this % sin_zkz(0:kmax) )

end subroutine

subroutine SinCosXYZ_dealloc(this) !DOC 
!DOC Deallocate SinCosXYZ
!DOC Parameters:
   Type(TSinCosXYZ) :: this !DOC SinCosXYZ

   deallocate( this % cos_xkx )
   deallocate( this % cos_yky )
   deallocate( this % cos_zkz )
   
   deallocate( this % sin_xkx )
   deallocate( this % sin_yky )
   deallocate( this % sin_zkz )

end subroutine

subroutine SinCosXYZ_fill(this,x,y,z,kmax) !DOC
!DOC initialize SinCosXYZ with x,y,z and given max |k|
   use constants, only : two_pi
!DOC Parameters:
   Type(TSinCosXYZ) :: this !DOC SinCosXYZ
   real(8),intent(in) :: x,y,z !DOC x,y,z
   integer,intent(in) :: kmax !DOC max |k|

   real(8) :: sx,sy,sz ! sin(2pix,y,z) 
   real(8) :: cx,cy,cz ! cos(2pix,2piy,2piz)

   integer :: k

   cx=COS(two_pi*x)  ! cos ( 2 pi *1 /L * x)  where L=1
   sx=SIN(two_pi*x)  
   cy=COS(two_pi*y)
   sy=SIN(two_pi*y)
   cz=COS(two_pi*z)
   sz=SIN(two_pi*z)
 
   ! at kx=0: cos(0) = 1, sin(0) = 0
   this % cos_xkx(0)=1.  
   this % sin_xkx(0)=0.
   this % cos_yky(0)=1.
   this % sin_yky(0)=0.
   this % cos_zkz(0)=1.
   this % sin_zkz(0)=0.
 
    ! Table of the functions
    ! cos(2pi m /L x) = cos( x*k_m )
    ! sin(k_m x), cos(k_m y) ... 
   do k=1,kmax
      this % cos_xkx(k) = this % cos_xkx(k-1)*cx - this % sin_xkx(k-1)*sx   
      this % sin_xkx(k) = this % sin_xkx(k-1)*cx + this % cos_xkx(k-1)*sx

      this % cos_yky(k) = this % cos_yky(k-1)*cy - this % sin_yky(k-1)*sy   
      this % sin_yky(k) = this % sin_yky(k-1)*cy + this % cos_yky(k-1)*sy

      this % cos_zkz(k) = this % cos_zkz(k-1)*cz - this % sin_zkz(k-1)*sz   
      this % sin_zkz(k) = this % sin_zkz(k-1)*cz + this % cos_zkz(k-1)*sz
   end do
 
end subroutine

! ******************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!! SinCosKR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *****************************************************************************

subroutine SinCosKR_alloc(this,grid) !DOC 
!DOC Allocate SinCosKR structure
!DOC Parameters:
   Type(TSinCosKR) :: this !DOC SinCosKR 
   Type(TFourierGrid),intent(in),target :: grid !DOC KSpace grid

   integer :: nk

   this % grid => grid

   nk = grid % nk
   allocate( this % cos_ppp( nk ) )
   allocate( this % sin_ppp( nk ) )

   allocate( this % cos_ppm( nk ) )
   allocate( this % sin_ppm( nk ) )
 
   allocate( this % cos_pmp( nk ) )
   allocate( this % sin_pmp( nk ) )
   
   allocate( this % cos_pmm( nk ) )
   allocate( this % sin_pmm( nk ) )

!   call SinCosKR_nulify(this)  <---  nobody needs it

end subroutine

!subroutine SinCosKR_nulify(this)
!   Type(TSinCosKR) :: this 
!
!   this % cos_ppp(:) = 0
!   this % sin_ppp(:) = 0
!   
!   this % cos_ppm(:) = 0
!   this % sin_ppm(:) = 0
!
!   this % cos_pmp(:) = 0
!   this % sin_pmp(:) = 0
!
!   this % cos_pmm(:) = 0
!   this % sin_pmm(:) = 0
!
!end subroutine 


subroutine SinCosKR_dealloc(this) !DOC
!DOC Deallocate SinCosKR structure
!DOC Parameters:
   Type(TSinCosKR) :: this  !DOC SinCosKR

   deallocate( this % cos_ppp )
   deallocate( this % sin_ppp )

   deallocate( this % cos_ppm )
   deallocate( this % sin_ppm )
 
   deallocate( this % cos_pmp )
   deallocate( this % sin_pmp )
   
   deallocate( this % cos_pmm )
   deallocate( this % sin_pmm )

   nullify(this % grid )

end subroutine



subroutine SinCosKR_fill(this,sincos_xyz,charge) !DOC 
!DOC  fill the sumcos for a given atom i.e. xkx+yky+zkz, xkx +yky - zkz etc... 
!DOC Parameters:
   Type(TSinCosKR),target :: this !DOC SinCosKR
!   real(8),intent(in) :: x,y,z
   Type(TSinCosXYZ),intent(in),target :: sincos_xyz   !DOC sin,cos(xkx,yky,zkz)
                                    !DOC to be filled before with 
                                    !DOC  SinCosXYZ_fill(sincos_xyz,x,y,z, this % grid % kmax ) 
   real(8),intent(in) :: charge  !DOC only for Coulomb sum. Use charge=1 for LJ 
   


! we do calculate the sin and cos using the reccurent formulae, i.e. 
! cos(xkx + yky + z(kz+dz) ) = cos( xkx+yky+zkz) * cos(zdz) - sin(xkx+yky+zkz)*sin(zdz)
! sin(xkx + yky + z(kz+dz) ) = sin( xkx+yky+zkz) * cos(zdz) + cos(xkx+yky+zkz)*sin(zdz)
   real(8) :: q_sin_xkx_plus_yky,  q_cos_xkx_plus_yky   ! charge*sin(xkx + yky) ,  charge*cos(xkx + yky)
   real(8) :: q_sin_xkx_minus_yky, q_cos_xkx_minus_yky  ! charge*sin(xkx - yky) ,  charge*cos(xkx - yky)


! temporary  variables
! in a sence, one can use only 4 of them. but this way it is more clear.
! the price is some 100 bytes in stack...
   real(8) :: q_cosxy_cosz   ! charge * cos(xkx + yky) * cos(zkz)
   real(8) :: q_sinxy_sinz   ! charge * sin(xkx + yky) * sin(zkz)
   real(8) :: q_cosxy_sinz   ! charge * cos(xkx + yky) * sin(zkz)
   real(8) :: q_sinxy_cosz   ! charge * sin(xkx + yky) * cos(zkz)

   real(8) :: q_cosxmy_cosz ! charge * cos(xkx - yky) * cos(zkz)
   real(8) :: q_sinxmy_sinz ! charge * sin(xkx - yky) * sin(zkz)
   real(8) :: q_cosxmy_sinz ! charge * cos(xkx - yky) * sin(zkz)
   real(8) :: q_sinxmy_cosz ! charge * sin(xkx - yky) * cos(zkz)

   real(8) :: q_cosx_cosy ! charge * cos(xkx) * cos(yky)
   real(8) :: q_sinx_siny ! charge * sin(xkx) * sin(yky)
   real(8) :: q_cosx_siny ! charge * cos(xkx) * sin(yky)
   real(8) :: q_sinx_cosy ! charge * sin(xkx) * cos(yky)

! k-space iterator
   integer :: kk,nk 
   integer :: kx,ky,kz

   real(8),dimension(:),pointer :: sin_xkx, sin_yky, sin_zkz  ! it seems that pointers will both: be clearer and work faster
   real(8),dimension(:),pointer :: cos_xkx, cos_yky, cos_zkz

   real(8),dimension(:),pointer :: cos_ppp, cos_ppm, cos_pmp, cos_pmm 
   real(8),dimension(:),pointer :: sin_ppp, sin_ppm, sin_pmp, sin_pmm 
 
 
   nk = this % grid % nk

!   sincos_xyz => this % sincos_xyz
   sin_xkx => sincos_xyz % sin_xkx
   sin_yky => sincos_xyz % sin_yky
   sin_zkz => sincos_xyz % sin_zkz 

   cos_xkx => sincos_xyz % cos_xkx
   cos_yky => sincos_xyz % cos_yky
   cos_zkz => sincos_xyz % cos_zkz 

   cos_ppp => this % cos_ppp
   sin_ppp => this % sin_ppp

   cos_ppm => this % cos_ppm
   sin_ppm => this % sin_ppm

   cos_pmp => this % cos_pmp
   sin_pmp => this % sin_pmp

   cos_pmm => this % cos_pmm
   sin_pmm => this % sin_pmm


   ! for kx,ky,kz=0: cos(0) = 1, sin(0) = 0

   q_cos_xkx_plus_yky = charge
   q_sin_xkx_plus_yky = 0.

   q_cos_xkx_minus_yky = charge
   q_sin_xkx_minus_yky = 0.

   do kk=1,nk

      kz = this % grid % kz(kk)
      ky = this % grid % ky(kk)

      q_cosxy_cosz = q_cos_xkx_plus_yky * cos_zkz(kz)
      q_sinxy_sinz = q_sin_xkx_plus_yky * sin_zkz(kz)
      q_sinxy_cosz = q_sin_xkx_plus_yky * cos_zkz(kz)
      q_cosxy_sinz = q_cos_xkx_plus_yky * sin_zkz(kz)

 !     write(*,*) 'ME:',q_cosxy_cosz,q_sinxy_sinz,q_sinxy_cosz,q_cosxy_sinz
!      write(*,*) 'ME: kx',this % grid % kx(kk),'ky',this % grid % ky(kk),'kz',kz

      ! ppp -> xkx + yky + zkz
      cos_ppp(kk) =  q_cosxy_cosz - q_sinxy_sinz  ! cos(X+Y+Z) = cos(X+Y)cos(Z) - sin(X+Y)sin(Z)
                                                                                      ! where X = xkx, Y=yky, Z=zkz
      sin_ppp(kk) =  q_sinxy_cosz + q_cosxy_sinz

      ! ppm -> xkx + yky - zkz
      if ( kz > 0 ) then
         cos_ppm(kk) =  q_cosxy_cosz + q_sinxy_sinz
         sin_ppm(kk) =  q_sinxy_cosz - q_cosxy_sinz
      else
         cos_ppm(kk) = 0.
         sin_ppm(kk) = 0.
      end if

!!! Now - the same for xkx - yky     
      q_cosxmy_cosz = q_cos_xkx_minus_yky * cos_zkz(kz)
      q_sinxmy_sinz = q_sin_xkx_minus_yky * sin_zkz(kz)
      q_sinxmy_cosz = q_sin_xkx_minus_yky * cos_zkz(kz)
      q_cosxmy_sinz = q_cos_xkx_minus_yky * sin_zkz(kz)

      ! pmp -> xkx - yky + zkz
      if ( ky > 0 ) then
         cos_pmp(kk) =  q_cosxmy_cosz - q_sinxmy_sinz  ! cos(X-Y+Z) = cos(X-Y)cos(Z) - sin(X-Y)sin(Z)
                                                                      ! where X = xkx, Y=yky, Z=zkz
         sin_pmp(kk) =  q_sinxmy_cosz + q_cosxmy_sinz
      else
         cos_pmp(kk) = 0.
         sin_pmp(kk) = 0.
      end if

      ! pmm -> xkx - yky - zkz
      if ( (kz > 0) .and. (ky > 0) ) then 
         cos_pmm(kk) =  q_cosxmy_cosz + q_sinxmy_sinz
         sin_pmm(kk) =  q_sinxmy_cosz - q_cosxmy_sinz
      else
         cos_pmm(kk) = 0.
         sin_pmm(kk) = 0.
      end if

      if( this % grid % ip(kk) /= 0) then ! next raw

         if( kk .lt. nk ) then
   
            kx = this % grid % kx(kk+1)
            ky = this % grid % ky(kk+1)
   
            q_cosx_cosy =  charge * cos_xkx(kx) * cos_yky(ky)
            q_sinx_siny =  charge * sin_xkx(kx) * sin_yky(ky)
            q_sinx_cosy =  charge * sin_xkx(kx) * cos_yky(ky) 
            q_cosx_siny =  charge * cos_xkx(kx) * sin_yky(ky)
   
            q_cos_xkx_plus_yky = q_cosx_cosy - q_sinx_siny
            q_sin_xkx_plus_yky = q_sinx_cosy + q_cosx_siny
            q_cos_xkx_minus_yky = q_cosx_cosy + q_sinx_siny
            q_sin_xkx_minus_yky = q_sinx_cosy - q_cosx_siny
   
     !       write(*,*) 'ME: ip',kk,kx,ky,cos_xkx(kx),cos_yky(ky)
   
         end if

      end if 


   end do

end subroutine

subroutine SinCosKR_nulify(this)  ! DOC
!DOC set all arrays to zero
!DOC Parameters:
   Type(TSinCosKR) :: this !DOC SinCosKR
   
   this % cos_ppp(:) = 0
   this % sin_ppp(:) = 0 

   this % cos_ppm(:) = 0 
   this % sin_ppm(:) = 0 

   this % cos_pmp(:) = 0 
   this % sin_pmp(:) = 0 

   this % cos_pmm(:) = 0 
   this % sin_pmm(:) = 0 

end subroutine


subroutine SinCosKR_copy(dst, src)  !DOC
!DOC Make a copy dst = src
!DOC Parameters:
   Type(TSinCosKR) :: dst !DOC destination
   Type(TSinCosKR),intent(in) :: src !DOC source
   
   dst % cos_ppp(:) = src % cos_ppp(:) 
   dst % sin_ppp(:) = src % sin_ppp(:) 

   dst % cos_ppm(:) = src % cos_ppm(:) 
   dst % sin_ppm(:) = src % sin_ppm(:) 

   dst % cos_pmp(:) = src % cos_pmp(:) 
   dst % sin_pmp(:) = src % sin_pmp(:) 

   dst % cos_pmm(:) = src % cos_pmm(:) 
   dst % sin_pmm(:) = src % sin_pmm(:) 

end subroutine


subroutine SinCosKR_mulScalar(this, multiplier )  !DOC
!DOC multiply by scalar:  this = this * multiplier
!DOC Parameters:
   Type(TSinCosKR),intent(inout) :: this !DOC SinCosKR
   real(8),intent(in) :: multiplier !DOC scalar multiplier
   
   this % cos_ppp(:) = this % cos_ppp(:) * multiplier 
   this % sin_ppp(:) = this % sin_ppp(:) * multiplier 

   this % cos_ppm(:) = this % cos_ppm(:) * multiplier 
   this % sin_ppm(:) = this % sin_ppm(:) * multiplier 

   this % cos_pmp(:) = this % cos_pmp(:) * multiplier 
   this % sin_pmp(:) = this % sin_pmp(:) * multiplier 

   this % cos_pmm(:) = this % cos_pmm(:) * multiplier 
   this % sin_pmm(:) = this % sin_pmm(:) * multiplier 

end subroutine


!!!!!!!!!!!!!!! SumSinCosKR !!!!!!!!!!!!!!!!!!! 

subroutine SumSinCosKR_alloc(this,grid) !DOC
!DOC allocate SumSinCosKR structure
!DOC Parameters:
   Type(TSumSinCosKR) :: this !DOC SumSinCosKR
   Type(TFourierGrid),intent(in),target :: grid !DOS KSpace grid

   integer :: nk

   this % grid => grid

   nk = grid % nk
   allocate( this % sumcos_ppp( nk ) )
   allocate( this % sumsin_ppp( nk ) )

   allocate( this % sumcos_ppm( nk ) )
   allocate( this % sumsin_ppm( nk ) )
 
   allocate( this % sumcos_pmp( nk ) )
   allocate( this % sumsin_pmp( nk ) )
   
   allocate( this % sumcos_pmm( nk ) )
   allocate( this % sumsin_pmm( nk ) )

   call SumSinCosKR_nulify(this)

end subroutine

subroutine SumSinCosKR_dealloc(this) !DOC
!DOC Deallocate SumSinCosKR
!DOC Parameters: 
   Type(TSumSinCosKR) :: this !DOC SumSinCosKR

   deallocate( this % sumcos_ppp )
   deallocate( this % sumsin_ppp )

   deallocate( this % sumcos_ppm )
   deallocate( this % sumsin_ppm )
 
   deallocate( this % sumcos_pmp )
   deallocate( this % sumsin_pmp )
   
   deallocate( this % sumcos_pmm )
   deallocate( this % sumsin_pmm )

   nullify(this % grid )

end subroutine

!!!!!!!!!!!!!!! OPERATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

subroutine SumSinCosKR_nulify(this) !DOC
!DOC Set SumSinCosKR to zero
!DOC Parameters:
  Type(TSumSinCosKR) :: this !DOC SumSinCosKR

   this % sumcos_ppp(:) = 0
   this % sumsin_ppp(:) = 0  

   this % sumcos_ppm(:) = 0
   this % sumsin_ppm(:) = 0  

   this % sumcos_pmp(:) = 0
   this % sumsin_pmp(:) = 0  

   this % sumcos_pmm(:) = 0
   this % sumsin_pmm(:) = 0  

end subroutine

subroutine SumSinCosKR_addAtom(this,sincos_kr) !DOC
!DOC add atom
!DOC Parameters:
  Type(TSumSinCosKR) :: this !DOC SumSinCosKR
  Type(TSinCosKR),intent(in) :: sincos_kr !DOC sincoskr for the new atom 

   this % sumcos_ppp(:) = this % sumcos_ppp(:) + sincos_kr % cos_ppp(:)
   this % sumsin_ppp(:) = this % sumsin_ppp(:) + sincos_kr % sin_ppp(:)

   this % sumcos_ppm(:) = this % sumcos_ppm(:) + sincos_kr % cos_ppm(:)
   this % sumsin_ppm(:) = this % sumsin_ppm(:) + sincos_kr % sin_ppm(:)

   this % sumcos_pmp(:) = this % sumcos_pmp(:) + sincos_kr % cos_pmp(:)
   this % sumsin_pmp(:) = this % sumsin_pmp(:) + sincos_kr % sin_pmp(:)

   this % sumcos_pmm(:) = this % sumcos_pmm(:) + sincos_kr % cos_pmm(:)
   this % sumsin_pmm(:) = this % sumsin_pmm(:) + sincos_kr % sin_pmm(:)

end subroutine

subroutine SumSinCosKR_addSum(this,sumsincos) !DOC
!DOC add sum (many atoms) this = this + sumsincos
!DOC Parameyers:
  Type(TSumSinCosKR) :: this !DOC SumSinCosKR
  Type(TSumSinCosKR),intent(in) :: sumsincos !DOC another SumSinCosKR for the new atoms

   this % sumcos_ppp(:) = this % sumcos_ppp(:) + sumsincos % sumcos_ppp(:)
   this % sumsin_ppp(:) = this % sumsin_ppp(:) + sumsincos % sumsin_ppp(:)

   this % sumcos_ppm(:) = this % sumcos_ppm(:) + sumsincos % sumcos_ppm(:)
   this % sumsin_ppm(:) = this % sumsin_ppm(:) + sumsincos % sumsin_ppm(:)

   this % sumcos_pmp(:) = this % sumcos_pmp(:) + sumsincos % sumcos_pmp(:)
   this % sumsin_pmp(:) = this % sumsin_pmp(:) + sumsincos % sumsin_pmp(:)

   this % sumcos_pmm(:) = this % sumcos_pmm(:) + sumsincos % sumcos_pmm(:)
   this % sumsin_pmm(:) = this % sumsin_pmm(:) + sumsincos % sumsin_pmm(:)

end subroutine

subroutine SumSinCosKR_subSum(this,sumsincos) !DOC
!DOC substract atoms: this = this - sumsincos
!DOC Parameters:
  Type(TSumSinCosKR) :: this !DOC this  SumSinCosKR
  Type(TSumSinCosKR),intent(in) :: sumsincos !DOC atoms to substract

   this % sumcos_ppp(:) = this % sumcos_ppp(:) - sumsincos % sumcos_ppp(:)
   this % sumsin_ppp(:) = this % sumsin_ppp(:) - sumsincos % sumsin_ppp(:)

   this % sumcos_ppm(:) = this % sumcos_ppm(:) - sumsincos % sumcos_ppm(:)
   this % sumsin_ppm(:) = this % sumsin_ppm(:) - sumsincos % sumsin_ppm(:)

   this % sumcos_pmp(:) = this % sumcos_pmp(:) - sumsincos % sumcos_pmp(:)
   this % sumsin_pmp(:) = this % sumsin_pmp(:) - sumsincos % sumsin_pmp(:)

   this % sumcos_pmm(:) = this % sumcos_pmm(:) - sumsincos % sumcos_pmm(:)
   this % sumsin_pmm(:) = this % sumsin_pmm(:) - sumsincos % sumsin_pmm(:)

end subroutine

subroutine SumSinCosKR_copy(dst, src)  !DOC
!DOC Copy: dst = src
!DOC Parameters:
   Type(TSumSinCosKR) :: dst !DOC destination
   Type(TSumSinCosKR),intent(in) :: src !DOC source
   
   dst % sumcos_ppp(:) = src % sumcos_ppp(:) 
   dst % sumsin_ppp(:) = src % sumsin_ppp(:) 

   dst % sumcos_ppm(:) = src % sumcos_ppm(:) 
   dst % sumsin_ppm(:) = src % sumsin_ppm(:) 

   dst % sumcos_pmp(:) = src % sumcos_pmp(:) 
   dst % sumsin_pmp(:) = src % sumsin_pmp(:) 

   dst % sumcos_pmm(:) = src % sumcos_pmm(:) 
   dst % sumsin_pmm(:) = src % sumsin_pmm(:) 

end subroutine


subroutine SumSinCosKR_mulScalar(this, multiplier )   !DOC
! multiply by scalar
!DOC Parameters:
   Type(TSumSinCosKR),intent(inout) :: this !DOC SumSinCosKR
   real(8),intent(in) :: multiplier !DOC scalar multiplier
   
   this % sumcos_ppp(:) = this % sumcos_ppp(:) * multiplier 
   this % sumsin_ppp(:) = this % sumsin_ppp(:) * multiplier 

   this % sumcos_ppm(:) = this % sumcos_ppm(:) * multiplier 
   this % sumsin_ppm(:) = this % sumsin_ppm(:) * multiplier 

   this % sumcos_pmp(:) = this % sumcos_pmp(:) * multiplier 
   this % sumsin_pmp(:) = this % sumsin_pmp(:) * multiplier 

   this % sumcos_pmm(:) = this % sumcos_pmm(:) * multiplier 
   this % sumsin_pmm(:) = this % sumsin_pmm(:) * multiplier 

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

End Module SumSinCosKR




