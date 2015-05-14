Module MonteCarloMove !DOC 
!
!DOC functions for choosing of new position and calculation of the 
!DOC  Molecule-independent functions! 
!
!DOC With minimal changes taken from MC code for H2O by Prof. Luc Belloni
!
implicit none

Type :: TMonteCarloMove  !DOC
!DOC Fields:
   ! old
   real(8) :: fxtotold,fytotold,fztotold,ftotold    !DOC old forces 
   real(8) :: couplxold,couplyold,couplzold,couplold !DOC old torques
   real(8) :: umaxold  !DOC lambda_f * F_tot * dr
   real(8) :: vmaxold !DOC lambda_c * tau_tot * d_angle
   real(8) :: cold,qold !DOC 

   real(8) :: vx,vy,vz !DOC 
   real(8) :: deltx,delty,deltz !DOC displacement
   real(8) :: angle !DOC 

   real(8) :: fxtotnew,fytotnew,fztotnew,ftotnew !DOC new forces
   real(8) :: couplxnew,couplynew,couplznew,couplnew !DOC new torq
   real(8) :: umaxnew,vmaxnew !DOC 
   real(8) :: cnew,qnew !DOC 

End Type TMonteCarloMove



contains

subroutine MonteCarloMove_init_old(this, fx_old,fy_old,fz_old,torq_x_old,torq_y_old,torq_z_old) !DOC 
!DOC initialize the move structure with old forces and torques
!DOC Parameters:
   Type(TMonteCarloMove) :: this !DOC MonteCarloMove structure

   real(8),intent(in) :: fx_old,fy_old,fz_old !DOC old force (of the molecue)
   real(8),intent(in) :: torq_x_old,torq_y_old,torq_z_old !DOC old torq (molecular)

   this % fxtotold = fx_old
   this % fytotold = fy_old
   this % fztotold = fz_old

   this % ftotold = SQRT( fx_old**2 + fy_old**2 + fz_old**2 ) 

   this % couplxold = torq_x_old
   this % couplyold = torq_y_old
   this % couplzold = torq_z_old

   this % couplold = SQRT( torq_x_old**2 + torq_y_old**2 + torq_z_old**2 )

end subroutine

subroutine MonteCarloMove_init_new(this, fx_new,fy_new,fz_new,torq_x_new,torq_y_new,torq_z_new) !DOC
!DOC Set new forces and torques
   use parameters , only : dr,d_angle,xlambda_f,xlambda_c
   use Functions, only : shi_sans_exp_f
!DOC Parameters:
   Type(TMonteCarloMove) :: this !DOC MonteCarloMove

   real(8),intent(in) :: fx_new, fy_new, fz_new !DOC new force
   real(8),intent(in) :: torq_x_new, torq_y_new, torq_z_new !DOC new torq

   real(8) :: umaxnew,e,t,cnew
   real(8) :: vmaxnew,qnew

   this % fxtotnew = fx_new
   this % fytotnew = fy_new
   this % fztotnew = fz_new

   this % ftotnew = SQRT( fx_new**2 + fy_new**2 + fz_new**2 ) 

   this % couplxnew = torq_x_new
   this % couplynew = torq_y_new
   this % couplznew = torq_z_new

   this % couplnew = SQRT( torq_x_new**2 + torq_y_new**2 + torq_z_new**2 )

   umaxnew=xlambda_f*this % ftotnew*dr
   IF(umaxnew.gt.1.d-2) then   ! a ne faire que si vraiment biaise
      e=EXP(-2.d0*umaxnew)
      t=(umaxnew*(1.d0+e)-(1.d0-e))/2.d0                 ! on gere les overflows luc70p197
      cnew=3.d0*t/umaxnew**3                         ! a exp(Umax) pres
   else
      umaxnew=0.
      cnew=1.d0
   endif

   this % umaxnew = umaxnew
   this % cnew = cnew

   vmaxnew=xlambda_c*this % couplnew*d_angle
   
   IF(vmaxnew.gt.1.d-2) then   ! a ne faire que si vraiment biaise
      qnew=shi_sans_exp_f(vmaxnew)/vmaxnew  ! je ne mets pas le facteur exp(vmax) pour eviter les overflows
   else
      vmaxnew=0.
      qnew=1.d0
   endif

   this % vmaxnew = vmaxnew
   this % qnew = qnew

end subroutine


subroutine MonteCarloMove_chooseDisplacement(this,deltx,delty,deltz)
   use parameters, only : dr,xlambda_f
   use MRandom
   use geometry, only : prod_vect

!DOC determine the displacement
!DOC Taken from MC code for H2O by Luc Belloni, with minimal changes
!DOC is universal for any molecules
!
!
!
!DOC     deplacement dans une sphere de rayon dr  avec proba exp(lamda.F.deltar)/Cold
!DOC        luc70p194 et luc84p160
!
!
!DOC Parameters:
  implicit none
 Type(TMonteCarloMove) :: this !DOC MonteCarloMove
 real(8),intent(out) :: deltx,delty,deltz !DOC result 

  real(8) :: fxtotold,fytotold,fztotold  ! F^M_old
  real(8) :: ftotold  ! abs(F_old^M)
  real(8) :: umax,umaxold  ! xlambda_f*ftotold*dr, used in acceptance
                                  ! (as i understand - to avoid the overflows) 
  real(8) :: e  ! exp(-2umax)
  real(8) :: sh,ch
  real(8) :: t,z,s,u,g,h,h1
  real(8) :: cold
  real(8) :: xksi ! random number
  real(8) :: xksishmax
  real(8) :: deltar
  real(8) :: xcostet,xsintet
  integer :: k
  real(8) :: vx,vy,vz  ! random vector
  real(8) :: v2        ! vx^2 + vy^2 + vz^2 
  real(8) :: wx,wy,wz  ! w = prod_vect( v,  F )
  real(8) :: w         ! |w|  

  fxtotold = this % fxtotold 
  fytotold = this % fytotold
  fztotold = this % fztotold

  ftotold = this % ftotold

 ! dr=dr_a/xl_a            ! car dr_a saisi dans mcluc2 est en unite A
! dr = dr_a / BoxLength
umax=xlambda_f*ftotold*dr


if(umax.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   !        d'abord choisir la norme deltar suivant deltar*sh(lambda*Fi*deltar)
   e=EXP(-2.d0*umax)
   t=(umax*(1.d0+e)-(1.d0-e))/2.d0                 ! on gere les overflows luc70p197
   cold=3.d0*t/umax**3                         ! a exp(Umax) pres
!   xksi=alea(1)
   xksi = 0.d0
   do while ( xksi .lt. 1d-10 )
      xksi = rand()
   end do

   !      IF(xksi*t<=0) PRINT*, 'log(xksi*t):',xksi,t

   z=umax+LOG(xksi*t)
   !       il faut trouver U=lambda*Fi*deltar tel que g(u)=UchU-shU=s donc u=g-1(s)
   !       ou plutot U tq h(u)=ln(g(U))=U+ln(g(u)/exp(u))=z=lns
   if(z<=LOG(4.d0)) then                  ! pt de depart luc84p169
      s=EXP(z)
      u=(3.d0*s)**(1.d0/3.d0)-s/10.d0
   else
      u=z-LOG(z/2.d0)
   endif

   do k=1,10                         ! qq NR  (10 NR iterations)
      e=EXP(-2.d0*u)
      ch=(1.d0+e)/2.d0
      sh=(1.d0-e)/2.d0
 
      if(u<1.d-5) sh=u-u**2

      g=u*ch-sh
      if(u<1.d-3) g=u**3/3.d0

      h=u+LOG(g)
      h1=u*sh/g
      if(k<10) u=u+(z-h)/h1            ! pas dans la derniere iteration, pour que u et e soient coherents

   end do

   !      IF(ABS(z-h)>1.d-6) PRINT*, 'attention, u pas converge a fond!',u,h,z
   deltar=u/(xlambda_f*ftotold)

   !       puis choisir x=costetha suivant exp(U*x)
   xksi = 1.d0   
   do while( (xksi .gt. (1.d0-1.d-10)) .or. (xksi .lt. 1.d-10))
      xksi= rand()
   end do

   xcostet=1.d0+LOG(xksi+(1.d0-xksi)*e)/u
   
   if(u<1.d-7) xcostet=1.d0-2.d0*(1.d0-xksi)
   !      IF(ABS(xcostet)>1) PRINT*, 'xcostet:',u,e,xksi,xcostet,1.d0-xcostet**2

   if ( xcostet > 1.d0 - 1d-14 ) xcostet = 1.d0 - 1d-14
   if ( xcostet < -1.d0 +1d-14 ) xcostet = -1.d0 + 1d-14

   xsintet=SQRT(1.d0-xcostet**2)
   !       vecteur au hasard en direction
   w=-1.d0
   
   do while(w <= 1.d-10 )

      v2=2.
      do WHILE(v2>1.)                ! vecteur v qq
         vx=2.*rand()-1.
         vy=2.*rand()-1.
         vz=2.*rand()-1.
         v2=vx**2+vy**2+vz**2
      end do

   !       vecteur w=v*F, perpendiculaire à F
      call prod_vect(vx,vy,vz,fxtotold,fytotold,fztotold,wx,wy,wz)

      w=SQRT(wx**2+wy**2+wz**2)
   end do  ! while w<=0

!       enfin deplacement suivant Fi et w
   deltx=deltar*(xcostet*fxtotold/ftotold+xsintet*wx/w)
   delty=deltar*(xcostet*fytotold/ftotold+xsintet*wy/w)
   deltz=deltar*(xcostet*fztotold/ftotold+xsintet*wz/w)

else ! (umax.le.1.d-2)  where umax=xlambda_f*ftotold*dr  
     !  si lambda=0 , on tire au hasard dans la sphere
         
   v2=2.
   do while(v2>1.)                ! vecteur v qq
      vx=2.*rand()-1.
      vy=2.*rand()-1.
      vz=2.*rand()-1.
      v2=vx**2+vy**2+vz**2
   end do

   deltx=dr*vx; delty=dr*vy; deltz=dr*vz
   umax=0.
   cold=1.d0

endif ! (umax.gt.1.d-2) 

this % deltx = deltx
this % delty = delty
this % deltz = deltz

this % cold = cold

umaxold=umax

this % umaxold = umaxold

end subroutine ! chooseDisplacement


subroutine MonteCarloMove_chooseRotation(this,rot) !DOC
  use MRandom, only : rand
  use parameters, only : d_angle, xlambda_c
  use geometry, only : prod_vect , TRotation
  use Functions, only : shi_sans_exp_f

!DOC determine the rotation
!DOC Taken from MC code for H2O by Luc Belloni, with minimal changes
!DOC is universal for any molecules
!DOC Parameters:
 implicit none
  Type(TMonteCarloMove) :: this !DOC MonteCarloMove
  Type(TRotation) :: rot !DOC rotation (in Luc's format, see geometry.f90) 
 
  real(8) :: couplxold,couplyold,couplzold  ! total torque 
  real(8) :: couplold ! abs(torque)
  real(8) :: vmax,vmaxold  ! xlambda_c*couplold*d_angle used in acceptance 
  real(8) :: vx,vy,vz  

  real(8) :: emax,shmax
  real(8) :: xksi ! random number 
  real(8) :: xksishmax
  real(8) :: proba
  real(8) :: e,v
  real(8) :: v2 ! = vx^2 + vy^2 + vz^2 
  real(8) :: wx,wy,wz  ! w = prod_vect( v, C)
  real(8) :: w  ! abs(w)
  real(8) :: sinx,cosx
  real(8) :: angle

  real(8) :: ca,ca1,sa
  real(8) :: st,ct
  real(8) :: cphi,sphi 

  real(8) :: qold 


  couplxold = this % couplxold
  couplyold = this % couplyold
  couplzold = this % couplzold

  couplold = this % couplold
!
!
!       rotation d'un angle alpha choisi dans -d_angle,+d_angle
!       par rapport a un vecteur V
!       avec proba exp(lamda*C.V*alpha)  ou C=couple
!       luc84p157,162,185
!
!
!call prod_vect(xx(i0+2)-xx(i0+1),yy(i0+2)-yy(i0+1),zz(i0+2)-zz(i0+1),fxold(2),fyold(2),fzold(2),couplx2,couply2,couplz2)
!call prod_vect(xx(i0+3)-xx(i0+1),yy(i0+3)-yy(i0+1),zz(i0+3)-zz(i0+1),fxold(3),fyold(3),fzold(3),couplx3,couply3,couplz3)

!couplxold=couplx2+couplx3                   ! couple total en ajoutant les 2 contributions des 2 H
!couplyold=couply2+couply3
!couplzold=couplz2+couplz3

!couplold=SQRT(couplxold**2+couplyold**2+couplzold**2)
!          d'abord sortir x=cos(angle entre V et C) entre 0 et 1

vmax=xlambda_c*couplold*d_angle

if(vmax.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   !       il faut tirer v entre 0 et vmax suivant shv/v
   !       pas possible, on tire donc d'abord suivant chv

   emax=EXP(-2.d0*vmax)
   shmax=(1.d0-emax)/2.d0  ! sans exp(vmax)
   xksi=1.; proba=0.

   do while(xksi>=proba)

      ! there was an error: NaN Produced in the next line   proba=(1.d0-e)/v/(1.d0+e)      
      ! it can be if v=0, then e=exp(-2v) = 1, (1-e)/v/(1+e) = NaN
      ! to avoid this, just excude too small v
      v=0.d0
      do while(dabs(v) .lt. 1.d-10) 
         xksi=rand()
         xksishmax=xksi*shmax
         v=vmax+LOG(xksishmax+SQRT(xksishmax**2+emax))   ! v suivant chv
         if(xksishmax**2<1.d-14*emax) v=xksishmax/SQRT(emax)
      end do 

      e=EXP(-2.d0*v)

      proba=(1.d0-e)/v/(1.d0+e)       ! accepte avec proba shv/v / chv
      if(v<1.d-5) proba=(2.d0*v-2.d0*v**2)/v/(1.d0+e)

      xksi=rand()
   end do

   cosx=v/vmax                     ! cos de l'angle entre C et V
   qold=shi_sans_exp_f(vmax)/vmax  ! je ne mets pas le facteur exp(vmax) pour eviter les overflows
   !        puis sortir l'angle de rotation alpha entre -d_angle et d_angle suivant exp(lamda*C.V*alpha)

   xksi = 1.d0
   do while( xksi .gt. (1d0 - 1d-10 ) ) 
      xksi=rand()
   end do

   angle=d_angle*(1.d0+LOG(xksi+(1.d0-xksi)*e)/v)

   if(v<1.d-7) angle=d_angle*(1.d0-2.d0*(1.d0-xksi))
   !    determiner l'orientation absolue de V qui fait l'angle x avec C
   !    vecteur au hasard en direction
   w=-1.d0
   do while(w <= 1.d-10 ) 
  
      v2=2.
      do WHILE( (v2>1.d0) .or. (v2<1d-10) )                ! vecteur v qq
         vx=2.d0*rand()-1.d0
         vy=2.d0*rand()-1.d0
         vz=2.d0*rand()-1.d0
         v2=vx**2+vy**2+vz**2
      end do

      !  vecteur w=v*C, perpendiculaire à C
      call prod_vect(vx,vy,vz,couplxold,couplyold,couplzold,wx,wy,wz)

      w=SQRT(wx**2+wy**2+wz**2)

   end do

   !IF(1.d0-cosx**2<0.) PRINT*, 'sinx:',vmax,v,cosx,1.d0-cosx**2
   sinx=SQRT(1.d0-cosx**2)
   vx=couplxold/couplold*cosx+wx/w*sinx              ! coord de l'axe de rotation V
   vy=couplyold/couplold*cosx+wy/w*sinx
   vz=couplzold/couplold*cosx+wz/w*sinx
   !PRINT*, cosx,sinx,couplxold*wx+couplyold*wy+couplzold*wz

else ! (vmax.le.1.d-2) where vmax = xlambda_c*couplold*d_angle
     ! si lambda=0 , on tire au hasard V

   v2=2.d0
   do WHILE( (v2>1.d0).or.(v2<1d-10) )                ! vecteur V qq
      vx=2.d0*rand()-1.d0
      vy=2.d0*rand()-1.d0
      vz=2.d0*rand()-1.d0
      v2=vx**2+vy**2+vz**2
   END do

   v=SQRT(v2)
   vx=vx/v; vy=vy/v; vz=vz/v
   angle=d_angle*(2.d0*rand()-1.d0)
   vmax=0.
   qold=1.d0

endif  ! (vmax.gt.1.d-2)

vmaxold=vmax
!                               je determine les angles de V dans le repere fixe
ca=COS(angle); sa=SIN(angle); ca1=1.d0-ca
ct=vz
!IF(1.d0-ct**2<0.) PRINT*, 'st:',cosx,sinx,couplzold,couplold,wz,w,vz,ct,1.d0-ct**2

st=SQRT(1.d0-ct**2)
IF(st>0.) THEN; cphi=vx/st; sphi=vy/st; endif

rot % cos_phi = cphi
rot % sin_phi = sphi
rot % cos_theta = ct 
rot % sin_theta = st 
rot % cos_angle = ca
rot % sin_angle = sa 

this % vx = vx
this % vy = vy
this % vz = vz

this % vmaxold = vmaxold

this % angle = angle

this % qold = qold

end subroutine ! chooseRotation


subroutine MonteCarloMove_accept_or_decline( this, dutot, accepted ) !DOC
!DOC accept or decline the move
   use parameters, only : xlambda_f , xlambda_c
   use MRandom, only : rand
!DOC Parameters:
   implicit none
   Type(TMonteCarloMove) :: this !DOC MonteCarloMove
   real(8),intent(in) :: dutot !DOC total change of the energy
   logical,intent(out) :: accepted !DOC output: accepted or declined

   real(8) :: vmaxold,vmaxnew
   real(8) :: umaxold,umaxnew
   real(8) :: cold,cnew,qold,qnew

   real(8) :: fxtotold,fytotold,fztotold,ftotold
   real(8) :: couplxold,couplyold,couplzold,couplold
   real(8) :: fxtotnew,fytotnew,fztotnew,ftotnew
   real(8) :: couplxnew,couplynew,couplznew,couplnew
   real(8) :: deltx,delty,deltz,vx,vy,vz
   real(8) :: angle

   real(8) :: xksi,factor

   real(8) :: exposu,exposf,exposfcorr,exposc,exposccorr

   vmaxold = this % vmaxold
   vmaxnew = this % vmaxnew
   
   umaxold = this % umaxold
   umaxnew = this % umaxnew

   cold = this % cold
   cnew = this % cnew

   qnew = this % qnew
   qold = this % qold

   fxtotold = this % fxtotold
   fytotold = this % fytotold
   fztotold = this % fztotold
   ftotold = this % ftotold

   fxtotnew = this % fxtotnew
   fytotnew = this % fytotnew
   fztotnew = this % fztotnew
   ftotnew = this % ftotnew

   couplxold = this % couplxold
   couplyold = this % couplyold
   couplzold = this % couplzold
   couplold = this % couplold

   couplxnew = this % couplxnew
   couplynew = this % couplynew
   couplznew = this % couplznew
   couplnew = this % couplnew
  
   deltx = this % deltx
   delty = this % delty
   deltz = this % deltz

   vx = this % vx
   vy = this % vy
   vz = this % vz

   angle = this % angle

   accepted = .FALSE.

!   IF(vmaxnew.gt.1.d-2) then   ! a ne faire que si vraiment biaise
!      qnew=shi_sans_exp_f(vmaxnew)/vmaxnew  ! je ne mets pas le facteur exp(vmax) pour eviter les overflows
!   else
!      vmaxnew=0.
!      qnew=1.d0
 !  endif
   
   !PRINT*, 'old'
   !PRINT*, fxold,fyold,fzold
   !PRINT*, fxtotold,fytotold,fztotold
   !PRINT*, couplxold,couplyold,couplzold
   !PRINT*, umaxold,vmaxold,cold,qold
   !PRINT*, 'new'
   !PRINT *, uinew
   !PRINT*, fxnew,fynew,fznew
   !PRINT*, fxtotnew,fytotnew,fztotnew
   !PRINT*, couplxnew,couplynew,couplznew
   !PRINT*, umaxnew,vmaxnew,cnew,qnew
   !PRINT*, dutotk,dutotk6,dutotk12
   !
   !        on accepte avec proba min(1,(CQ)old/(CQ)new*exp(-deltaU-lamda*(Fold+Fnew).deltar-lamda(Cold+Cnew).V*angle))
   !        ne pas oublier que C et Q ne contiennent pas les facteurs exp(umax) et exp(vmax)
   !        donc les rajouter dans l'exposant total  luc84p164
   !
   exposu=-dutot
   
   exposf=-xlambda_f*           &
    ((fxtotold+fxtotnew)*deltx+(fytotold+fytotnew)*delty+(fztotold+fztotnew)*deltz)    ! allen tildesley p224
   exposfcorr=umaxold-umaxnew               ! luc70p196
   
   exposc=-xlambda_c*angle*          &
    ((couplxold+couplxnew)*vx+(couplyold+couplynew)*vy+(couplzold+couplznew)*vz)    ! allen tildesley p224
   exposccorr=vmaxold-vmaxnew               ! luc70p196
   
   factor=cold*qold/(cnew*qnew)*EXP(exposu+exposf+exposfcorr+exposc+exposccorr)
   !PRINT*, 'facteur: ',factor
   !        on accepte avec la proba=factot
   xksi=rand()

!   write(*,*) 'ME: xksi',xksi,'factor',factor
!   write(*,*) 'cold',cold,'qold',qold,'cnew',cnew,'qnew',qnew
!   write(*,*) 'exposu',exposu,'exposf',exposf
!   write(*,*) 'exposfcorr',exposfcorr,'exposc',exposc,'exposccorr',exposccorr
   
   IF(xksi>factor) return                  ! NON!!! on n'accepte pas, et on sort!
!   IF(i_partie/=0) return          ! on n'accepte pas si test de vitesse partiel
   !
   !
   !        et donc OUI, on accepte le deplacement!
   !
   !
!   iaccep=1
   accepted = .TRUE.

end subroutine

subroutine MonteCarloMove_accept_or_decline_simple( dutot, accepted ) !DOC
!DOC without the force bias (used in exchange move)
   use MRandom, only : rand
!DOC Parameters:
   implicit none
   real(8),intent(in) :: dutot !DOC total change of energy
   logical,intent(out) :: accepted !DOC output: accepted or declined

   real(8) :: xksi,factor

   accepted = .FALSE.
   factor = exp(-dutot)
 
   xksi=rand()

!   write(*,*) 'accept_or_decline: du=',dutot,'factor=',factor,'xksi=',xksi

   IF(xksi>factor) return                  ! NON!!! on n'accepte pas, et on sort!

   ! dutot<0 ==> -dutot>0 ==> factor exp(-dutot) > 1 > xksi   (always accepted )
   ! dutot=+inf ==> exp(-dutot) = factor = 0 < xksi ( always declined)
   ! otherwise - acceptance proportional to exp(-beta Delta U)

   accepted = .TRUE.

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module MonteCarloMove
