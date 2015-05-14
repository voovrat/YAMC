Module EwaldSumRealSpace !DOC 
!DOC !FILE Data structures and functions which deal with the real component of the Ewald sum
use AtomicData
use LJTypes

!
!    SUM q_iq_j erfc(alpha r) / r 
!  + SUM 4 eps_ij sigma_ij^6   exp(-a^2)/r^6 ( 1 + a^2 + a^4/2)
!  + SUM 4 eps_ij sigma_ij^12  exp(-a^2)/r^12 ( 1 + a^2 + a^4/2 + a^6/6 + a^8/24 + a^10/120 )
!
!  a = alpha r
!  Forces:
!
!  F_p = SUM_ij A_ij F_R(r_ij) (r_i - r_j) 
!  
!  FR^C = 1/r^2 ( erfc(a) + 2/sqrt(pi) exp(-a^2) ) 
!
! for w=2k+2
!  FR^w = w exp(-a^2) ( C_w(r) / r^w + alpha^w / (w/2)! ) / r^2
!
!  where C_w(r) = SUM_{l=0}^{w/2-1} a^2l / l!
!
!  FR^6 = 6 exp(-a^2) ( C_6(r)/r^6 + alpha^6 / 6 ) /r^2
!  FR^12 = 12 exp(-a^2) ( C_12(r)/r^12 + alpha^12/720 ) / r^2
!
!  FR^LJ =4 epsilon( sigma^12 FR^12 - sigma^6 FR^6) =
!  = 4 epsilon exp(-a^2) * { 6S [ 2S (C12 + A^2/720)  - C6 - A/6 ] } / r^2
!
!  where S = (sigma/r)^6, A = (alpha r)^6 = a^6 
implicit none

Type :: TEwaldSumRealSpace !DOC
!DOC Fields:
!    real(8) :: rmax,rmax2 now in parameters

   Type(TAtomicData),pointer :: atomic_data !DOC AtomicData (coordinates, charges)
   Type(TLJTypes),pointer :: lj_types !DOC LJ parametrs for each pair of atom types

   integer :: first_atom,last_atom !DOC the first and the last atoms (used for the partial sums which include only the atoms of the certain molecule)

   integer :: natom,nalloc !DOC number of atoms and the size of the allocated arrays

   real(8) ::  uu_ew, uu_lj6, uu_lj12  !DOC electrostatic, LJ6 and LJ12 components of the energy
   real(8) ::  uu_ew_intra,uu_lj6_intra,uu_lj12_intra !DOC components of the intra-molecular interaction (not used to my knowledge)

   real(8),dimension(:),pointer :: uu  !DOC energy per atom. U = sum_i uu(i)
   real(8),dimension(:),pointer :: uu_back  !DOC for the full sum uu_back => uu. Otherwise - it is the energy change per molecule
   real(8),dimension(:),pointer :: Fx,Fy,Fz  !DOC 1..Nnew: forces on the atoms of the  molecule (at new coordinates)
   real(8),dimension(:),pointer :: Fx_back,Fy_back,Fz_back  !DOC 1..Ntot: forces of the new coords on the all molecules in the system
                                    !DOC  for the full sum Fx,y,z_back => Fx,y,z
   
   real(8),dimension(:),pointer :: xmol,ymol,zmol  !DOC  1..Nmol: positions of atoms of the molecule (pointer to xx,yy,zz)
                                        !DOC for the full sum xmol,ymol,zmol => atomic_data % xx,yy,zz


   logical :: full_sum  !DOC logcal which indicates the full (or molecular) sum

End Type

contains 

!
! calculates SUM_{i =1}^{N_new} SUM_{j=1}^N 
!             (  q_i q_j erfc(alpha |r_new(i) - r(j)|) / |r_new(i) - r(j)| +
!               4epsilon sigma_ij^12 ( C12 / |r_new(i) - r(j) |^12 )
!               - 4epsilon sigma_ij^6( C6 / |r_new(i) - r(j)|^6 ) 
!             )
!
! 
subroutine EwaldSumRealSpace_alloc_full(this, atomic_data, lj_types ) !DOC
!DOC allocates the arrays for the full sum
!DOC Parameters:
   Type(TEwaldSumRealSpace) :: this !DOC EwaldSumRealSpace structure
   Type(TAtomicData),intent(in),target :: atomic_data !DOC coordinates and charges
   Type(TLJTypes),intent(in),target :: lj_types !DOC LJ parameters for each pair of atom types
   integer :: ntotal

   this % atomic_data => atomic_data 
   this % lj_types => lj_types 

   this % full_sum = .TRUE.

   ntotal = atomic_data % Natom
   this % natom = ntotal


   this % xmol => atomic_data % xx
   this % ymol => atomic_data % yy
   this % zmol => atomic_data % zz

   allocate(this % uu(ntotal))
   this % uu_back => this % uu 
 

   allocate(this % Fx(ntotal))
   allocate(this % Fy(ntotal))
   allocate(this % Fz(ntotal))   

   this % Fx_back => this % Fx
   this % Fy_back => this % Fy
   this % Fz_back => this % Fz
 
end subroutine

subroutine EwaldSumRealSpace_alloc_partial(this,atomic_data,lj_types) !DOC
!DOC allocate arrays for the partial (molecular) sum
   Type(TEwaldSumRealSpace) :: this !DOC EwaldSumRealSpace
   Type(TAtomicData),intent(in),target :: atomic_data !DOC AtomicData (coordinates and charges)
   Type(TLJTypes),intent(in),target :: lj_types !DOC LJ parameters for each pair of atom types

   integer :: ntotal

   this % atomic_data => atomic_data 
   this % lj_types => lj_types 

   this % full_sum = .FALSE.

   this % natom = 0
   ntotal = atomic_data % Natom;

   nullify( this % xmol )
   nullify( this % ymol )
   nullify( this % zmol )

   allocate(this % uu_back(ntotal))

   allocate(this % Fx_back(ntotal))
   allocate(this % Fy_back(ntotal))
   allocate(this % Fz_back(ntotal))
 
end subroutine

subroutine EwaldSumRealSpace_set_molecule(this, xmol, ymol, zmol, first_atom, last_atom ) !DOC
!DOC set the molecule (to the partial sum only) 
   use error
!DOC Parameters:
   implicit none
   Type(TEwaldSumRealSpace) :: this !DOC EwaldSumRealSpace
   real(8),dimension(:),intent(in),target  :: xmol,ymol,zmol !DOC coordinates of the molecule atoms (either pointer to AtomicData or xnew,ynew,znew arrays)
   integer,intent(in) :: first_atom, last_atom !DOC first_atom, last_atom are used to know charges and lj parameters

   integer :: nmolatom

   if( this % full_sum ) then
      write(error_message,*) 'EwaldSimRealSpace_set_molecule: function not valid for the full  sum!' 
      call error_throw(ERROR_WRONG_FUNCTION)
      return
   end if 


   this % xmol => xmol
   this % ymol => ymol
   this % zmol => zmol

   this % first_atom = first_atom
   this % last_atom = last_atom

   nmolatom = last_atom - first_atom + 1 

   if ( this % nalloc < nmolatom ) then

      if( this % nalloc > 0) then

         deallocate( this % Fx )
         deallocate( this % Fy )
         deallocate( this % Fz )

         deallocate( this % uu )

      end if 

      allocate(this % uu(nmolatom))

      allocate(this % Fx(nmolatom))
      allocate(this % Fy(nmolatom))
      allocate(this % Fz(nmolatom))   

      this % nalloc = nmolatom

   end if

   this % natom = nmolatom

end subroutine


subroutine EwaldSumRealSpace_dealloc_full(this) !DOC
!DOC Deallocate the full sum 
!DOC Parameters:
   use error
   Type(TEwaldSumRealSpace) :: this !DOC EwaldSumRealSpace structure

   if ( .not. this % full_sum ) then
     write(error_message,*) 'EwaldSumRealSpace_dealloc_full: cannot deallocate partial sum! use dealloc_partial instead! '
     call error_throw(ERROR_WRONG_FUNCTION)
     return
   end if 

   deallocate(this % uu)
   deallocate(this % Fx)
   deallocate(this % Fy)
   deallocate(this % Fz)

end subroutine

subroutine EwaldSumRealSpace_dealloc_partial(this) !DOC 
!DOC Deallocate the partial sum
   use error
   Type(TEwaldSumRealSpace) :: this !DOC EwaldSumRealSpace structure

   if ( this % full_sum ) then
     write(error_message,*) 'EwaldSumRealSpace_dealloc_partial: cannot deallocate full_sum! use dealloc_full instead! '
     call error_throw(ERROR_WRONG_FUNCTION)
     return
   end if 

   deallocate(this % Fx_back)
   deallocate(this % Fy_back)
   deallocate(this % Fz_back)
   deallocate(this % uu_back)

   if( this % nalloc > 0 ) then 
       deallocate(this % Fx)
      deallocate(this % Fy)
      deallocate(this % Fz)
      deallocate(this % uu)

  end if

end subroutine


subroutine EwaldSumRealSpace_calc(this,overlap,first,last,dont_stop_on_overlap) !DOC
!DOC Calculate the sum. The energy per atom uu and forces arrays will also be initialized
   use RealSumLocal  ! subroutines for calculation of energy and forces
   use parameters, only : rmax2,BoxLength
   use io, only : write_real_array
!
!DOC note: function works for both cases : full sum and sum of one molecule
!DOC       in case of the full sum xnew = xx, ynew = yy, znew = zz
!DOC       otherwise, other arrays xmol,ymol,zmol should be provided
!
!DOC for the full sum run EwaldSumRealSpace_calc(this,this % atomic_data % xx,this % atomic_data % yy, this % atomic_data % zz)
!DOC Parameters:
   Type(TEwaldSumRealSpace),target :: this !DOC EwaldSimRealSpace
   logical,intent(out) :: overlap !DOC logical output: indicates if the overlap occured during the sum calculation
   integer,intent(in),optional :: first,last  !DOC for the calculation of the molecule-molecule interactions, like u12. For the usual (even partial) sum calculations first=1, last=natom (or can be omited and thus set by default).
   logical,intent(in),optional :: dont_stop_on_overlap
   integer :: i,j,j0,j1
   integer :: ii  ! ii = i + first_atom - 1, where i=1,natom
                  ! for the full sum ii=i
   integer :: imol,jmol,ityp,jtyp
   real(8) :: q_i,q_j
   real(8) :: xij,yij,zij
   real(8) :: fxij,fyij,fzij
   real(8) :: r2  ! r^2
   real(8) :: qi,qj  ! charges 

   Type(TRealSumLocal) :: real_sum_local

   integer :: first_atom,last_atom
   integer,dimension(:),pointer :: molnum
   logical :: full_sum

   real(8),dimension(:),pointer :: charge
   real(8),dimension(:),pointer :: xx,yy,zz
   real(8),dimension(:),pointer :: xnew,ynew,znew

   Type(TLJTypes),pointer :: lj_types
   integer,dimension(:),pointer :: type_by_index
   real(8) :: four_epsilon, sigma6
   real(8) :: sigma_squared 
   real(8),dimension(:),pointer :: hard_core_angstr 

   integer :: base_offset,offset,ntype
   integer :: natom  ! number of atoms in partial sum (if full natom = ntotal )
   integer :: ntotal ! number of atoms in the system

   real(8) :: d_HS_i,d_HS_j,d_HS_12


   full_sum = this % full_sum

   molnum => this % atomic_data % molnum_by_atomnum 
   charge => this % atomic_data % charge
   lj_types => this % lj_types
   type_by_index => lj_types % type_by_index

   ntype = lj_types % NType

   xx => this % atomic_data % xx
   yy => this % atomic_data % yy
   zz => this % atomic_data % zz

   hard_core_angstr => this % atomic_data % hard_core_angstr 

   xnew => this % xmol
   ynew => this % ymol
   znew => this % zmol

   first_atom = this % first_atom
   last_atom = this % last_atom

   natom = this % natom
   ntotal = this % atomic_data % Natom

!   write(*,*) 'EwaldSumRealSpace_calc: '
!   write(*,*) 'x='
!   call write_real_array(0,xnew,natom)
!   write(*,*) 'y='
!   call write_real_array(0,ynew,natom)
!   write(*,*) 'z='
!   call write_real_array(0,znew,natom)


   overlap = .FALSE.

   this % uu(:) = 0
   this % uu_back(:) =0
 
   this % Fx(:) = 0
   this % Fx_back(:) = 0

   this % Fy(:) = 0
   this % Fy_back(:) = 0

   this % Fz(:) = 0
   this % Fz_back(:) = 0

   this % uu_ew = 0
   this % uu_lj6 = 0
   this % uu_lj12 = 0

   this % uu_ew_intra = 0
   this % uu_lj6_intra = 0
   this % uu_lj12_intra = 0

   do i=1,natom

      if(  full_sum ) then 
         j0 = i+1 
         ii = i ! ii - index of the molecule atom in the full list
                ! for the full sum it is the same of the inter-molecular index
      else
         j0 = 1
         ii = i + this % first_atom - 1  ! for the partial sum the full-sum index has an offset
      end if 

      j1 = ntotal

      if( present(first) ) j0 = first
      if( present(last) ) j1 = last
     

      imol = molnum( ii )
      ityp = type_by_index( ii )

      qi = charge(ii)
      d_HS_i = hard_core_angstr(ii) / BoxLength

      if ( ityp > 0 ) then
         base_offset = (ityp - 1) * ntype
      else
         base_offset = -1000 ! to get the error message if something is done with it (it should not be the case)
      end if

      do j=j0,j1
  
         jmol = molnum(j)
         jtyp = type_by_index( j  ) 


          if( (.not. full_sum) .and. ( imol == jmol ) ) cycle
   !       if( imol == jmol) cycle ! do not account inter-molecular interactions only in case of not-full sum !WHY?
!  because in full sum it is inconsistent with Ewald? 
!  maybe yes: ewald cutoff depends on BoxLength, so when we change it, we spliting between Rspace-kspace changes, so we need to recalculate everything 

     !    write(*,*) 'imol:',imol,' jmol:',jmol,'ityp:',ityp,' jtyp:',jtyp
 
         qj = charge(j)
         d_HS_j = hard_core_angstr(j) / BoxLength 
 
         d_HS_12 = (d_HS_i + d_HS_j) / 2

         if( (ityp>0) .and. (jtyp > 0) ) then
            offset = base_offset + jtyp 
            sigma_squared = lj_types % sigma2_tab( offset )
         end if 

         xij=xnew(i)-xx(j)
         yij=ynew(i)-yy(j)                ! donc compris entre -1 et 1
         zij=znew(i)-zz(j)
 
         if(xij.gt. 0.5) xij=xij-1.   ! nearliest neighbour. Remember: distance units are BoxLength
         if(xij.lt.-0.5) xij=xij+1.   ! i.e. Box is [-0.5:0.5]^3
 
         if(yij.gt. 0.5) yij=yij-1.
         if(yij.lt.-0.5) yij=yij+1.         ! prend moins de temps que anint!
 
         if(zij.gt. 0.5) zij=zij-1.
         if(zij.lt.-0.5) zij=zij+1.
 
         r2=xij**2+yij**2+zij**2
 
    !    write(*,*) 'xij=',xij,'yij=',yij,'zij=',zij,'r2=',r2

         
          
         if( (ityp>0) .and. (jtyp>0) .and. (imol /= jmol) .and. (r2<d_HS_12**2) )  then        ! overlap OO LJ, r<0.63*sigma:  on sort

!            write(*,*) 'sigma_squared=',sigma_squared,'0.4*sigma_squared',0.4*sigma_squared,'d_HS_12**2',d_HS_12**2
!            write(*,*) 'd_HS_i',d_HS_i * BoxLength,'d_HS_j',d_HS_j*BoxLength
!            write(*,*) 'xij=',xij,'yij=',yij,'zij=',zij
!
!            write(*,*) 'sigma_squared=',sigma_squared,' r2=',r2
!            write(*,*) 'OVERLAP: ii=',ii,'j=',j,'ityp',ityp
!            write(*,*) 'ii',xx(ii),yy(ii),zz(ii)
!            write(*,*) 'i(new): ',xnew(i),ynew(i),znew(i)
!            write(*,*) 'j:',xx(j),yy(j),zz(j) ,'jtyp',jtyp
 
            write(*,'(A,$)') 'O'            
             !stop
            overlap = .TRUE.

            if (present(dont_stop_on_overlap)) then 
              if ( .not. dont_stop_on_overlap) then
                return
              end if
            else
               return  
            end if

         end if 



         if(r2 .ge. rmax2 ) cycle

         call RealSumLocal_init_electrostatics( real_sum_local, r2, (imol == jmol) )
         call RealSumLocal_calc_electrostatics( real_sum_local, qi, qj )  ! potew, FR , pot = potew

         this % uu_ew = this % uu_ew + real_sum_local % potew
     
         if ( imol == jmol ) this % uu_ew_intra = this % uu_ew_intra + real_sum_local % potew
   

!          write(*,*) 'FR_C:',real_sum_local % FR 

         if( (ityp /=0) .and. (jtyp /= 0) ) then! if both atoms has LJ potential

             four_epsilon = lj_types % four_epsilon_tab( offset ) 
             sigma6 = lj_types % sigma6_tab( offset )   

             call RealSumLocal_init_LJ( real_sum_local, sigma6 ) 
             call RealSumLocal_calc_LJ( real_sum_local, four_epsilon ) ! potlj6, potlj12, pot += potlj12 - potlj6, FR += FR_LJ

             this % uu_lj6 = this % uu_lj6 + real_sum_local % potlj6
             this % uu_lj12 = this % uu_lj12 + real_sum_local % potlj12


             if ( imol == jmol ) then
                this % uu_lj6_intra = this % uu_lj6_intra + real_sum_local % potlj6
                this % uu_lj12_intra = this % uu_lj12_intra + real_sum_local % potlj12
             end if

!              write(*,*) 'FR_tot:',real_sum_local % FR
 
         end if
              ! not same molecules
               
         !write(*,*) 'ME : pot:',real_sum_local % pot,'ew:',real_sum_local % potew,'lj6:',&
          !                      real_sum_local % potlj6,'lj12:',real_sum_local % potlj12

         if( imol /= jmol ) then
            this % uu(i) = this % uu(i) + real_sum_local % pot                  ! ne contient plus le self
            this % uu_back(j) = this % uu_back(j)+ real_sum_local % pot                  ! mais contient l'intra!

            ! F_i = SUM_j FR(rij) ( r_i - r_j)  <-- fx,fy,fz
            ! [F_j]_M = SUM_(i in M) FR(rij) (r)j - r_i)  <-- fx_re

            fxij=xij * real_sum_local % FR                     ! NON, plus maintenant!

!            write(*,*) 'FR:',real_sum_local % FR,'fxij:',fxij,'x(i):',xx(i),'x(j):',xx(j),'xij',xij
            this % Fx(i) = this % Fx(i) + fxij
            this % Fx_back(j) = this % Fx_back(j)-fxij                  

            fyij=yij * real_sum_local % FR 
            this % Fy(i) = this % Fy(i) + fyij
            this % Fy_back(j)= this % Fy_back(j)-fyij   

            fzij=zij * real_sum_local % FR
            this % Fz(i) = this % Fz(i) + fzij
            this % Fz_back(j) = this % Fz_back(j)-fzij
         end if ! imol /= jmol

     end do  ! j

   end do ! i

end subroutine


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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module
