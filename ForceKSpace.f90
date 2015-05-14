Module ForceKSpace !DOC 
!DOC !FILE Module to calculate forces originating from KSpace Ewald sum 
use SumSinCosKR
use FourierGrid
use SystemSettings, only : TRealArrayPointer
!
! Forces in KSpace
! 
! F_p^C = 2q_p SUM_m>=0 (1 + sgn(kx)) beta(k_m)
!                  k+++.* ( C+++ sin(k+++ r_p) - S+++ cos(k+++ r_p) )
!      +   sgn(ky) k+-+.* ( C+-+ sin(k+-+ r_p) - S+-+ cos(k+-+ r_p) )
!      +   sgn(kz) k++-.* ( C++- sin(k++- r_p) - S++- cos(k++- r_p) )
!     + sgn(ky kz) k+--.* ( C+-- sin(k+-- r_p) - S+-- cos(k+-- r_p) )
!
! where k(sxsysz) = (sx*kx,sy*ky,sz*kz),   sx,sy,sz in {+1,-1} = {+,-}
!
! it can be rewiriten as: 
!
!  F_p^C = 2 SUM_m>=0 (1+sgn(kx)) beta(k_m) * k_m * (
!                 (+++).* ( C+++ q_p sin(k+++ r_p) - S+++ q_p cos(k+++ r_p) )
!      +   sgn(ky) (+-+).* ( C+-+ q_p sin(k+-+ r_p) - S+-+ q_p cos(k+-+ r_p) )
!      +   sgn(kz) (++-).* ( C++- q_p sin(k++- r_p) - S++- q_p cos(k++- r_p) )
!     + sgn(ky kz) (+--).* ( C+-- q_p sin(k+-- r_p) - S+-- q_p cos(k+-- r_p) )
!                   )
!               
!   where .* is a point-wise multiplication
!   note : q_p sin(k(sxsysz) r_p) is calculated in SinCosKR_fill
!
!
! F_p(tau)^LJ = 2 SUM m>=0 SUM_t beta_LJ(tau,t) 
!                  k+++ ( C^t+++ sin(k+++ r_p^tau) - S^t+++ cos(k+++ r_p^tau) )
!      +   sgn(ky) k+-+ ( C^t+-+ sin(k+-+ r_p^tau) - S^t+-+ cos(k+-+ r_p^tau) )
!      +   sgn(kz) k++- ( C^t++- sin(k++- r_p^tau) - S^t++- cos(k++- r_p^tau) )
!     + sgn(ky kz) k+-- ( C^t+-- sin(k+-- r_p^tau) - S^t+-- cos(k+-- r_p^tau) )
! 
!
!  we rewrite it as follows:
!
! F_p(tau)^LJ = 2 SUM_t SUM_m>=0 k * beta_LJ(tau,t) 
!                  (+++).* ( C^t+++ sin(k+++ r_p^tau) - S^t+++ cos(k+++ r_p^tau) )
!      +   sgn(ky) (+-+).* ( C^t+-+ sin(k+-+ r_p^tau) - S^t+-+ cos(k+-+ r_p^tau) )
!      +   sgn(kz) (++-).* ( C^t++- sin(k++- r_p^tau) - S^t++- cos(k++- r_p^tau) )
!     + sgn(ky kz) (+--).* ( C^t+-- sin(k+-- r_p^tau) - S^t+-- cos(k+-- r_p^tau) )
! 
!
!
!

implicit none

contains

! F_p^C = 2q_p SUM_m>=0 (1 + sgn(kx)) beta(k_m)
!                  k+++.* ( C+++ sin(k+++ r_p) - S+++ cos(k+++ r_p) )
!      +   sgn(ky) k+-+.* ( C+-+ sin(k+-+ r_p) - S+-+ cos(k+-+ r_p) )
!      +   sgn(kz) k++-.* ( C++- sin(k++- r_p) - S++- cos(k++- r_p) )
!     + sgn(ky kz) k+--.* ( C+-- sin(k+-- r_p) - S+-- cos(k+-- r_p) )
!
subroutine ForceKSpace_coulomb(beta,sumsincos_coulomb,sincos_kr_coulomb, fx_tot, fy_tot, fz_tot ) !DOC
!DOC Coulomb component of forces
   use constants, only : four_pi
!  fx,fy,fz = F_p^C(x,y,z)
!Parameters:
   real(8),dimension(:),intent(in) :: beta  !DOC in reality it is already (1+sgn(kx)) beta(k_m) 
   Type(TSumSinCosKR),intent(in),target :: sumsincos_coulomb !DOC SumSinCosKR structures for coulomb part: SUM_i sin(kR_i), SUM_i cos(kR_i)
   Type(TSinCosKR),intent(in),target :: sincos_kr_coulomb !DOC in this case it is q_p sin,cos(sx*x*kx + sy*y*ky +sz*zk*z) 
                                                   !DOC to be calculated with SinCosKR_fill... 
   real(8)  :: fx_tot,fy_tot,fz_tot !DOC output: components of the force
   real(8)  :: fx,fy,fz
   real(8)  :: A
 
   integer :: nk 
   integer :: kk 

   Type(TFourierGrid),pointer :: grid
   integer,dimension(:),pointer :: kx,ky,kz

   real(8),dimension(:),pointer :: C_ppp,S_ppp,C_pmp,S_pmp,C_ppm,S_ppm,C_pmm,S_pmm  ! C,S(sx,sy,sz) sx=+, sy,sz in {+,-}
   real(8),dimension(:),pointer :: sin_ppp,sin_ppm,sin_pmp,sin_pmm  ! sin(sx*x*kx + sy*y*ky + sz*z*kz )
   real(8),dimension(:),pointer :: cos_ppp,cos_ppm,cos_pmp,cos_pmm  ! cos(sx*x*kx + sy*y*ky + sz*z*kz )
   

   C_ppp => sumsincos_coulomb % sumcos_ppp
   C_pmp => sumsincos_coulomb % sumcos_pmp
   C_ppm => sumsincos_coulomb % sumcos_ppm
   C_pmm => sumsincos_coulomb % sumcos_pmm

   S_ppp => sumsincos_coulomb % sumsin_ppp
   S_pmp => sumsincos_coulomb % sumsin_pmp
   S_ppm => sumsincos_coulomb % sumsin_ppm
   S_pmm => sumsincos_coulomb % sumsin_pmm

   cos_ppp => sincos_kr_coulomb % cos_ppp
   cos_pmp => sincos_kr_coulomb % cos_pmp
   cos_ppm => sincos_kr_coulomb % cos_ppm
   cos_pmm => sincos_kr_coulomb % cos_pmm

   sin_ppp => sincos_kr_coulomb % sin_ppp
   sin_pmp => sincos_kr_coulomb % sin_pmp
   sin_ppm => sincos_kr_coulomb % sin_ppm
   sin_pmm => sincos_kr_coulomb % sin_pmm

 
   grid => sincos_kr_coulomb % grid

   nk = grid % nk

   kx => grid % kx
   ky => grid % ky  
   kz => grid % kz

   fx_tot = 0
   fy_tot = 0
   fz_tot = 0

   do kk=1,nk

!                 (+++).* ( C+++ q_p sin(k+++ r_p) - S+++ q_p cos(k+++ r_p) )
      A = C_ppp(kk) * sin_ppp(kk) - S_ppp(kk) * cos_ppp(kk)
      fx =  A !C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)
      fy =  A !C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)
      fz =  A !C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)

      if ( ky(kk) > 0 ) then
!         sgn(ky) (+-+).* ( C+-+ q_p sin(k+-+ r_p) - S+-+ q_p cos(k+-+ r_p) )
       A = C_pmp(kk) * sin_pmp(kk) - S_pmp(kk) * cos_pmp(kk)
         fx = fx +  A !C_pmp(kk) * sin_pmp(kk)  - S_pmp(kk) * cos_pmp(kk)
         fy = fy -  A !C_pmp(kk) * sin_pmp(kk)  + S_pmp(kk) * cos_pmp(kk)  ! consider the sign - for ky! (it is pmp )
         fz = fz +  A !C_pmp(kk) * sin_pmp(kk)  - S_pmp(kk) * cos_pmp(kk)
      end if 

      if ( kz(kk) > 0 ) then
!      +   sgn(kz) (++-).* ( C++- q_p sin(k++- r_p) - S++- q_p cos(k++- r_p) )
       A = C_ppm(kk) * sin_ppm(kk) - S_ppm(kk) * cos_ppm(kk)

         fx = fx +  A ! C_ppm(kk) * sin_ppm(kk)  - S_ppm(kk) * cos_ppm(kk)
         fy = fy +  A ! C_ppm(kk) * sin_ppm(kk)  - S_ppm(kk) * cos_ppm(kk)  
         fz = fz -  A ! C_ppm(kk) * sin_ppm(kk)  + S_ppm(kk) * cos_ppm(kk)  ! consider the sign - for kz! (it is ppm )
      end if 

      if ( (ky(kk) > 0) .and. (kz(kk) > 0) ) then
!     + sgn(ky kz) (+--).* ( C+-- q_p sin(k+-- r_p) - S+-- q_p cos(k+-- r_p) )

       A = C_pmm(kk) * sin_pmm(kk) - S_pmm(kk) * cos_pmm(kk)

         fx = fx +  A !C_pmm(kk) * sin_pmm(kk)  - S_pmm(kk) * cos_pmm(kk)
         fy = fy -  A !C_pmm(kk) * sin_pmm(kk)  + S_pmm(kk) * cos_pmm(kk)  ! consider the sign - for ky ! (it is pmm )
         fz = fz -  A !C_pmm(kk) * sin_pmm(kk)  + S_pmm(kk) * cos_pmm(kk)  ! consider the sign - for kz ! (it is pmm )
      end if 

  ! multiply by k_m * (1+sgn(kx)) beta(k_m):
      fx = fx * kx(kk) * beta(kk)  ! again: beta(kk) in our case is  indeed (1+sgn(kx)) beta(k_m) 
      fy = fy * ky(kk) * beta(kk) 
      fz = fz * kz(kk) * beta(kk)

      fx_tot = fx_tot + fx
      fy_tot = fy_tot + fy
      fz_tot = fz_tot + fz  

   end do  

! multiply the total sum by 2
   fx_tot = four_pi * fx_tot  ! why 4pi:  2 is the coefficient, before the sum, 2pi is the conversion from k to m 
   fy_tot = four_pi * fy_tot  ! ( k = 2pi m / L )
   fz_tot = four_pi * fz_tot

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!! LJ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ForceKSpace_LJ(beta_LJ,sumsincos_LJ,ntype,ityp,sincos_kr_ityp, fx_tot, fy_tot, fz_tot ) !DOC
!DOC Calculate the LJ component of forces
   use constants, only : four_pi
!  fx,fy,fz = F_p^C(x,y,z)
   Type(TRealArrayPointer),dimension(:),intent(in),target :: beta_LJ  !DOC 2D array (t1,t2) <--> beta_LJ( (t1-1)*Ntype + t2 )
   Type(TSumSinCosKR),dimension(:),intent(in),target :: sumsincos_LJ  !DOC sumsincos for all types (see SumSinCosKR)
   integer,intent(in) :: ntype !DOC number of types 
   integer,intent(in) :: ityp  !DOC type of the atom for which we are calculationg the forces
   Type(TSinCosKR),intent(in),target :: sincos_kr_ityp   !DOC sin,cos(sx*x*kx + sy*y*ky + sz*z*kz)  for that atom
   real(8) :: fx_tot,fy_tot,fz_tot !DOC output: the force components

   real(8) :: fx,fy,fz
   
   
   integer :: nk 
   integer :: kk 
   integer :: typ
   integer :: base_offset ! (ityp-1)*ntype


   Type(TFourierGrid),pointer :: grid
   integer,dimension(:),pointer :: kx,ky,kz


   real(8),dimension(:),pointer :: beta_LJ_current  ! pointer to betaLJ(typ,ityp)
   real(8),dimension(:),pointer :: C_ppp,S_ppp,C_pmp,S_pmp,C_ppm,S_ppm,C_pmm,S_pmm  ! C,S(sx,sy,sz) sx=+, sy,sz in {+,-}
   real(8),dimension(:),pointer :: sin_ppp,sin_ppm,sin_pmp,sin_pmm  ! sin(sx*x*kx + sy*y*ky + sz*z*kz )
   real(8),dimension(:),pointer :: cos_ppp,cos_ppm,cos_pmp,cos_pmm  ! cos(sx*x*kx + sy*y*ky + sz*z*kz )
   real(8) :: A
   

   cos_ppp => sincos_kr_ityp % cos_ppp
   cos_pmp => sincos_kr_ityp % cos_pmp
   cos_ppm => sincos_kr_ityp % cos_ppm
   cos_pmm => sincos_kr_ityp % cos_pmm

   sin_ppp => sincos_kr_ityp % sin_ppp
   sin_pmp => sincos_kr_ityp % sin_pmp
   sin_ppm => sincos_kr_ityp % sin_ppm
   sin_pmm => sincos_kr_ityp % sin_pmm

 
   grid => sincos_kr_ityp % grid

   nk = grid % nk

   kx => grid % kx
   ky => grid % ky  
   kz => grid % kz

   fx_tot = 0
   fy_tot = 0
   fz_tot = 0

   base_offset = (ityp - 1) * ntype  ! offset for beta_LJ for types (ityp,*)

   do typ=1,ntype

      C_ppp => sumsincos_LJ(typ) % sumcos_ppp
      C_pmp => sumsincos_LJ(typ) % sumcos_pmp
      C_ppm => sumsincos_LJ(typ) % sumcos_ppm
      C_pmm => sumsincos_LJ(typ) % sumcos_pmm
   
      S_ppp => sumsincos_LJ(typ) % sumsin_ppp
      S_pmp => sumsincos_LJ(typ) % sumsin_pmp
      S_ppm => sumsincos_LJ(typ) % sumsin_ppm
      S_pmm => sumsincos_LJ(typ) % sumsin_pmm

      beta_LJ_current => beta_LJ( base_offset + typ ) % ptr  ! beta_LJ(ityp,typ)

      do kk=1,nk   
         !     (+++).* ( C+++ q_p sin(k+++ r_p) - S+++ q_p cos(k+++ r_p) )
         A = C_ppp(kk) * sin_ppp(kk) - S_ppp(kk) * cos_ppp(kk) 
!         write(*,*) 'C_ppp',C_ppp(kk),'S_ppp',S_ppp(kk)
!         write(*,*) 'sin_ppp',sin_ppp(kk),'cos_ppp',cos_ppp(kk)

         fx =  A !C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)
         fy =  A !C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)
         fz =  A ! C_ppp(kk) * sin_ppp(kk)  - S_ppp(kk) * cos_ppp(kk)
  

   !     write(*,*) 'A_ppp',A 
         if ( ky(kk) > 0 ) then
            !   sgn(ky) (+-+).* ( C+-+ q_p sin(k+-+ r_p) - S+-+ q_p cos(k+-+ r_p) )
         A = C_pmp(kk) * sin_pmp(kk) - S_pmp(kk) * cos_pmp(kk) 

  !       write(*,*) 'A_pmp',A
            fx = fx +  A! C_pmp(kk) * sin_pmp(kk)  - S_pmp(kk) * cos_pmp(kk)
            fy = fy -  A! C_pmp(kk) * sin_pmp(kk)  + S_pmp(kk) * cos_pmp(kk)  ! consider the sign - for ky! (it is pmp )
            fz = fz +  A!C_pmp(kk) * sin_pmp(kk)  - S_pmp(kk) * cos_pmp(kk)
         end if 
   
         if ( kz(kk) > 0 ) then
!           write(*,*) 'MY ppm',kk
            !  +   sgn(kz) (++-).* ( C++- q_p sin(k++- r_p) - S++- q_p cos(k++- r_p) )
         A = C_ppm(kk) * sin_ppm(kk) - S_ppm(kk) * cos_ppm(kk) 
!         write(*,*) 'A_ppm',A

            fx = fx +  A !C_ppm(kk) * sin_ppm(kk)  - S_ppm(kk) * cos_ppm(kk)
            fy = fy +  A !C_ppm(kk) * sin_ppm(kk)  - S_ppm(kk) * cos_ppm(kk)  
            fz = fz -  A !C_ppm(kk) * sin_ppm(kk)  + S_ppm(kk) * cos_ppm(kk)  ! consider the sign - for kz! (it is ppm )
         end if 
   
         if ( (ky(kk) > 0) .and. (kz(kk) > 0) ) then
            !   + sgn(ky kz) (+--).* ( C+-- q_p sin(k+-- r_p) - S+-- q_p cos(k+-- r_p) )
         A = C_pmm(kk) * sin_pmm(kk) - S_pmm(kk) * cos_pmm(kk) 
 !        write(*,*) 'A_pmm',A

            fx = fx +  A !C_pmm(kk) * sin_pmm(kk)  - S_pmm(kk) * cos_pmm(kk)
            fy = fy -  A !C_pmm(kk) * sin_pmm(kk)  + S_pmm(kk) * cos_pmm(kk)  ! consider the sign - for ky ! (it is pmm )
            fz = fz -  A !C_pmm(kk) * sin_pmm(kk)  + S_pmm(kk) * cos_pmm(kk)  ! consider the sign - for kz ! (it is pmm )
         end if 
   
         ! multiply by k_m * beta_LJ(tau,t)(k_m):  (ok tau == ityp, t == typ )
         fx = fx * kx(kk) * beta_LJ_current(kk)  ! again: beta(kk) in our case is  indeed (1+sgn(kx)) beta(k_m) 
         fy = fy * ky(kk) * beta_LJ_current(kk) 
         fz = fz * kz(kk) * beta_LJ_current(kk)

         fx_tot = fx_tot + fx
         fy_tot = fy_tot + fy
         fz_tot = fz_tot + fz

     end do !  kk

   end do   ! typ

! multiply the total sum by 2
   fx_tot = four_pi * fx_tot ! 4pi = 2*2pi, 2 - sum prefactor, 2pi - koefficient in k to m conversion: k = 2pi m/L
   fy_tot = four_pi * fy_tot
   fz_tot = four_pi * fz_tot

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End Module ForceKSpace
