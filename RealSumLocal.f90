Module RealSumLocal !DOC
!DOC   auxilarly functions for particle-particle interactions
!DOC  to be used inside the Real Space Ewald Sums

!
!  
!    SUM q_iq_j erfc(alpha r) / r 
!  + SUM 4 eps_ij sigma_ij^6   exp(-a^2)/r^6 ( 1 + a^2 + a^4/2) 
!  + SUM 4 eps_ij sigma_ij^12  exp(-a^2)/r^12 ( 1 + a^2 + a^4/2 + a^6/6 + a^8/24 + a^10/120 )
!
! ( comment: (1+a^2+a^4/2) = C_6 
!            ( 1 + a^2 + a^4/2 + a^6/6 + a^8/24 + a^10/120 ) = C_12
!   definition of C see below    
! )
! 
!  a = alpha r
!  Forces:
!
!  F_p = SUM_ij F_R(i,j) (r_ij) (r_i - r_j) 
!  
!  FR^C = q_i q_j /r^2 ( erfc(a) + 2/sqrt(pi) exp(-a^2) ) 
!
! for w=2k+2
!  FR^w = w exp(-a^2) ( C_w(r) / r^w + alpha^w / (w/2)! ) / r^2
!
!  where C_w(r) = SUM_{l=0}^{w/2-1} a^2l / l!
!
!  FR^6 = 6 exp(-a^2) ( C_6(r)/r^6 + alpha^6 / 6 ) /r^2
!  FR^12 = 12 exp(-a^2) ( C_12(r)/r^12 + alpha^12/720 ) / r^2
!
!  FR^LJ(i,j) =4 epsilon( sigma^12 FR^12 - sigma^6 FR^6) =
!  = 4 epsilon exp(-a^2) * { 6S [ 2S (C12 + A^2/720)  - C6 - A/6 ] } / r^2
!
!  where S = (sigma/r)^6, A = (alpha r)^6 = a^6 

implicit none

 Type :: TRealSumLocal !DOC

   !DOC auxilarly variables
   real(8) :: r,r2 !DOC distance 
   real(8) :: a,a2,a6  !DOC  a = alpha * r, a2 = a^2, a6= a^6
   real(8) :: exp_minus_a2 !DOC exp(-a^2)   
   real(8) :: erfc_a   !DOC erfc(a) = erfc(alpha r )   
   real(8) :: C6,C12   !DOC C_w = sum_{l=0}^{w/2-1}  l^{2l} / l!  where w=6,12
                       !DOC used in LJ energy and forces
   real(8) :: four_epsilon !DOC 4eps
   real(8) :: sigma_over_r_6   !DOC (sigma/r)^6

   ! results of the calculations 
   real(8) :: potew  !DOC electrostatic interaction potential   q_iq_j * erfc(alpha r ) / r
   real(8) :: potlj12,potlj6 !DOC Lennard Jones potential compontnts:  U_w = A_ij exp(-alpha r)/r^w C_w, w=6,12
   real(8) :: pot   !DOC total potential potew + potlj12 - potlj6
   real(8) :: FR     !DOC Force radial component
                     !DOC  F_p = SUM_ij F_R(i,j) ( r_i - r_j) 
                     !DOC  where  FR = FR_C + FR_LJ  ( definitions see above, explanations in ewald.pdf ) 
       
End Type TRealSumLocal

contains


subroutine RealSumLocal_init_electrostatics( this, r2, sameMolecule ) !DOC
  use Functions, only : erfc_Luc,erf_Luc
  use parameters, only : alpha
!DOC init variables used for electrostatic energy and forces calculation
!DOC Parameters:
  Type(TRealSumLocal) :: this !DOC RelSumLocal
  real(8), intent(in) :: r2 !DOC r^2
  logical,intent(in) :: sameMolecule !DOC indicates that the particles are in the same molecule
                                     !DOC for compatibility with Luc Belloni code

  this % r = sqrt(r2)
!  write(*,*) 'ME: r=',this % r
  this % r2 = r2
  this % a = alpha * this % r
  this % a2 = this % a**2
  this % exp_minus_a2 = EXP(-this % a2)

  if ( sameMolecule ) then
     this % erfc_a = -erf_Luc(this % a, this % a2, this % exp_minus_a2 )
  else
     this % erfc_a = erfc_Luc(this % a,this % a2,this % exp_minus_a2) 
  end if

  this % pot = 0

end subroutine 

subroutine RealSumLocal_calc_electrostatics( this, qi, qj )  !DOC
  use parameters, only : xclb,two_alpha_over_sqrt_pi
!DOC calculate the electrostatic interactions
!DOC run init_electrostatics first
!DOC Parameters:
   Type(TRealSumLocal) :: this  !DOC RealSumLocal
   real(8),intent(in) :: qi, qj !DOC charges 

 
   this % potew = xclb*qi*qj* this % erfc_a / this % r
   this % pot =this % potew
     
!   write(*,*) 'xclb',xclb,'qi',qi,'qj',qj,'xclb*qiqj',xclb*qi*qj,'erfc',this % erfc_a,'r',this % r
 
   this % FR = xclb*qi*qj*( this % erfc_a /this % r + two_alpha_over_sqrt_pi * this % exp_minus_a2)/ this % r2     

!   write(*,*) 'xclb*qi*qj',xclb*qi*qj,'erfc(a):',this % erfc_a, 'r:',this % r, '2alp/sqrt(pi):',two_alpha_over_sqrt_pi
!   write(*,*) 'exp(-a^2):',this % exp_minus_a2,'r2:',this % r2,'FR:',this % FR
!    write(*,*) 'xclb=',xclb,'qi=',qi,'qj=',qj,'erfc_a=',this % erfc_a,'r=',this % r,'pot=',this % pot

end subroutine

subroutine RealSumLocal_init_LJ(this,sigma6) !DOC
!DOC note: LJ calculations should always be performed after the coulomb
!DOC at least: you should run RealSumLocal_init_coulomb before
!DOC Parameters:
   Type(TRealSumLocal) :: this !DOC RealSumLocal
   real(8),intent(in) :: sigma6 !DOC sigma^6

   real(8) :: a2  ! a^2 = (alpha r)^2
   real(8) :: r6  ! r^6

   a2 = this % a2

   r6 = (this % r2)**3
   this % sigma_over_r_6 = sigma6 / r6
!   write(*,*) 'sigma6=',sigma6,'r2=',this % r2
   this % C6 = 1.d0 + a2 * (1.d0 + a2/2.d0) !

   this % a6 = a2**3
   this % C12 = this % C6 + this % a6 * (1.d0 + a2*(1.d0 + a2/5.d0)/4.d0)/6.d0 ! C12 = 1 + a^2 + 0.5 a^4 + a^6/6 + a^8/24 + a^10/120

end subroutine

subroutine RealSumLocal_calc_LJ(this,four_epsilon) !DOC 
!DOC  Calculate LJ potentials and update FR (add LJ force radial component)
!DOC Parameters:
   Type(TRealSumLocal) :: this !DOC RealSumLocal
   real(8) :: four_epsilon !DOC 4eps

   real(8) :: sigma_over_r_6, exp_minus_a2,C6,C12,a6,r2

   sigma_over_r_6 = this % sigma_over_r_6 
   exp_minus_a2 = this % exp_minus_a2
   C6 = this % C6
   C12 = this % C12
   a6 = this % a6
   r2 = this % r2

   this % potlj6 = four_epsilon *  sigma_over_r_6 * exp_minus_a2 * C6    ! LJ6 potential, see the definition above
   this % potlj12 = four_epsilon *  sigma_over_r_6**2 * exp_minus_a2 * C12 ! LJ12 potential, see the definition above

   this % pot = this % pot + this % potlj12 - this % potlj6

!   write(*,*) 'pot6=',this % potlj6,'pot12=',this % potlj12,'potlj=', this % potlj12 - this % potlj6

!     ffr=ffr+xlj*6.d0*vv*ee2*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2
  ! LJ Force radial component, see the definition in the begining of the module, explanation ewald.pdf
   this % FR = this % FR + four_epsilon*6.d0 * sigma_over_r_6 * exp_minus_a2 *  &
                     (  2.d0 * sigma_over_r_6 * ( C12 + a6**2/720.d0 ) &
                        - (C6 + a6/6.d0) &
                     ) / r2
    
!   write(*,*)  'four_epsilon',four_epsilon,'sigma_over_r_6',sigma_over_r_6,'exp(-a^2)',exp_minus_a2,'a6',a6,'C6',C6,'C12',C12
                !
                ! 1/r^2(   4 eps * 12 sigma^12/r^12 e^-a^2 ( C12 + (alpha r)^12 / 6! )
                !        - 4 eps * 6 * sigma^6/r^6  e^-a^2 ( C6 + (alpha r)^6 / 3! )  )
                !
end subroutine 

!         if(r2 .ge. rmax2 ) cycle
!
!          r=SQRT(r2)
!         a=alpha * r
!         a2=a**2
!         exp_minus_a2=EXP(-a2)
!         erfc_a=erfc_Luc(a,a2,exp_minus_a2)          ! calcule erfc_Luc(x) en s'aidant de x**2 et exp(-x**2)
!         
! !        if(io==jo) er=-erf_Luc(alpr,alpr2,ee2)    ! retrancher le self eventuellement donc erfc_Luc-1=-erf  luc84p134
!                                                   ! io == jo <=> for the atoms of the same molecule 
!         ! erfc - 1 = - erf
!         ! Sergiievskyi 31 May 2014: 
!         ! > why not to remove the intermolecular interactions at all?
!         ! > after all: they are always constant for the rigid molecules, thus do not matter
!
!
!         potew = xclb*qi*qj*erfc_a/r
!         this % u_ew=this % u_ew + potew
!         pot=potew
!
!        
!        FR=xclb*qi*qj*(er/r+alp2pi*ee2)/r2     
!
!         if((ityp .gt. 0) .and. (jtyp .gt. 0) ) then            ! OO: LJ total  luc84p147
!
!            four_epsilon = LJTypes_four_epsilon(lj_types, ityp, jtyp )
!            sigma6 = LJTypes_sigma_power6(lj_types, ityp, jtyp )    
!            r6 = r2**3       
!
!            sigma_over_r_6 = sigma6 / r6
!            C6 = 1.d0 + a2 * (1.d0 + a2/2.d0) !
! 
!!           potlj6=-xlj*vv*ee2*c6            !  4*epsilon * (1+a^2 + 0.5a^4) exp(-a^2) (sigma/r)^6
!!           minus_potlj6 = -four_epsilon * sigma_over_r_6 * exp_minus_a2 * C6
!
!!           alpr6=alpr2**3                   !  a^6 
!            a6 = a2**3          
!
!!            c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0  
!            C12 = C6 + a6*(1.d0 + a2*(1.d0 + a2/5.d0)/4.d0)/6.d0 ! C12 = 1 + a^2 + 0.5 a^4 + a^6/6 + a^8/24 + a^10/120
!
!!            potlj12=xlj*vv**2*ee2*c12                              ! 4 epsilon * C12  * exp(-a^2) * (sigma/r)^12
!            potlj12 = four_epsilon * sigma_over_r_6**2 * exp_minus_a2 * C12
!
!!            uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
!            this % u_lj12 = this % u_lj12 + potlj12
!            this % minus_u_lj6 = this % minus_u_lj6 + minus_potlj6
!
!!            pot=pot+potlj12+potlj6
!            this % pot = this % pot + potlj12 + minus_potlj6
!
!!            ffr=ffr+xlj*6.d0*vv*ee2*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2
!            FR = FR + four_epsilon*6.d0 * sigma_over_r_6 * exp_minus_a2 *  &
!                              (  2.d0 * sigma_over_r_6 * ( C12 + a6**2/720.d0 ) &
!                                 - (C6 + a6/6.d0) &
!                              ) / r2
!                !
!                ! 1/r^2(   4 eps * 12 sigma^12/r^12 e^-a^2 ( C12 + (alpha r)^12 / 6! )
!                !        - 4 eps * 6 * sigma^6/r^6  e^-a^2 ( C6 + (alpha r)^6 / 3! )  )
!                !
!                
!             
!         end if
!
!              ! not same molecules
!!              if(io/=jo) then 
!              this % uu(i) = this % uu(i)+pot                  ! ne contient plus le self
!              this % uu(j) = this % uu(j)+pot                  ! mais contient l'intra!
!              fxij=xij*ffr                     ! NON, plus maintenant!
!              fx(i)=fx(i)+fxij
!              fx(j)=fx(j)-fxij                  
!              fyij=yij*ffr
!              fy(i)=fy(i)+fyij
!              fy(j)=fy(j)-fyij
!              fzij=zij*ffr
!              fz(i)=fz(i)+fzij
!              fz(j)=fz(j)-fzij
! !               end if ! io /= jo
! 



End Module RealSumLocal




