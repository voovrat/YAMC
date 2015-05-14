Module Functions !DOC
!DOC !FILE Functions used in calculation of Ewald sums
!DOC !FILE Written by Luc Belloni in MC for water
! 

contains

!
!
!
! 
!  
! w=6:   
!     S_fourier = pi^4.5/3V SUM_m F_2(k_m) h_m^3 [ sqrt(pi) erfc(b) + (1/2b^3 - 1/b) exp(-b^2 ] 
!                = SUM_m beta6(k_m) F_2(k_m)
!
! w=12: 
!     S_fourier =  1/(945*120) pi^10.5/V *
!                  *  SUM_m<>0 F_2(k_m) h_m^9 (-16 sqrt(pi) erfc(b) + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) )
!                  = SUM_m beta12(k_m) F_2(k_m) 
! 
! (see above in the real-space cycle)
!
!  Here we calculate beta6 and beta12
! 
subroutine erfc_Luc_bet6_bet12(x,x2,e2,erfk,bet6,bet12) !DOC
!DOC calculate the beta6 and beta12 values
!DOC Parameters:
     IMPLICIT REAL(8) (a-h,o-z) 
!DOC  x, x2 --> b, b^2
!DOC  e2    --> exp(-b^2)
!DOC  erfk  --> erfc(b)
!
!DOC  b^2 = pi^2 h_m^2 / alpha  --> h_m = sqrt(alpha) b/pi
! 
!DOC          avec erfc_Luc a x<0.5, DL a x>8, integration numerique de x a 8 entre les 2
!DOC          luc84p150
!

pi=4.d0*ATAN(1.d0)
pi12=SQRT(pi)
pi15=pi**1.5
xmin=0.5d0
xmax=8.d0
if(x<=xmin) then
   erfk=erfc_Luc(x,x2,e2)

   ! h_m = alpha b/pi,  h_m^3 = alpha^3 b^3 / pi^3
   !
   bet6=pi15/3.d0*x**3*(pi12*erfk+(0.5d0/x2-1.d0)/x*e2)  ! pi^4.5/3 * h_m^3 * ( sqrt(pi)*erfc(b) + (1/2b^3 - 1/b)exp(-b^2) ) = 
                                                         ! alpha^3 pi^1.5/3 * b^3 (sqrt(pi) erfc(b) + (1/2b^3 -1b)exp(-b^2) )
                                                         ! as i understand, we do multiply by alpha afterwards, after returning form the subroutine... 


! beta12 =1/(945*120) *pi^10.5 {  h_m^9 (-16 sqrt(pi) erfc(b) + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) ) }
! again: h_m^9 = alpha^9 b^9 / pi^9
! and we multiply the result by alpha^9 in the main code
! 
!   so, we compute 1/(945*120) pi^1.5 { b^9 ( -16 sqrt(pi) erfc(b) + [ .... ]*exp(-b^2) ) }
!
! ok. 945*120 = 1080 * 105 
! i.e.  1/1080* -16/105  sqrt(pi) erfc(b) = 1/(945*120) (-16 sqrt(pi) erfc(b) )
!  1/1080 * 1/ b^9  = 1/(945*120) * 105 / b^9
!  1/1080 * 2/7 1/b^7 =  1/(945*120) * 105 * 2/7 / b^7 = 30 / b^7
!  1/1080 * 4/35 1/b^5 = 1/(945*120) * 105 * 4/35 / b^5 = 12/ b^5
!  1/1080 * 8/105 1/b^3 = 1/(945*120) * 105 * 8/105 / b^3 = 8/b^3
!  1/1080 * 16/105 1/b  = 1/(945*120) * 105 * 16/105 / b = 16/b
   bet12=pi15/1080.d0*x**9*(-16.d0*pi12/105.d0*erfk+           &
         (((((1.d0/x2-2.d0/7.d0)/x2+4.d0/35.d0)/x2-8.d0/105.d0)/x2+16.d0/105.d0))/x*e2)

elseif(x>=xmax) then
   erfk=e2/pi12/x*(1.d0-0.5d0/x2*(1.d0-1.5d0/x2*(1.d0-2.5d0/x2)))    !+105.d0/16.d0/b2**4)
   bet6=pi15*e2/4.d0/x**2*(1.d0-2.5d0/x2*(1.d0-3.5d0/x2*(1.d0-4.5d0/x2)))
   bet12=pi15*e2/240.d0/x**2*(1.d0-5.5d0/x2*(1.d0-6.5d0/x2*(1.d0-7.5d0/x2)))
else

! Sergiievskyi 31 May 2014:
!  > Why use the integration? 
!  > Why use different formulae for different intervals?
!  > Why not to tabulate the function (as they are position-independent one-argument functions )
!

   dx=1.d-3                ! au lieu de 10-4, ca devrait suffir
   n=(xmax-x)/dx+2
   dx=(xmax-x)/(n-1)
   som=dx*e2/2.d0-dx**2/12.d0*2.d0*x*e2
   som6=dx*e2/x**4/2.d0-dx**2/12.d0*(2.d0/x**3+4.d0/x**5)*e2
   som12=dx*e2/x**10/2.d0-dx**2/12.d0*(2.d0/x**9+10.d0/x**11)*e2

   do i=2,n
       y=x+(i-1)*dx
       e=EXP(-y**2)
       som=som+dx*e
       som6=som6+dx*e/y**4
       som12=som12+dx*e/y**10
   end do

   som=som-dx*e/2.d0+dx**2/12.d0*2.d0*y*e
   y2=y**2
   som=som+e/2.d0/y*(1.d0-0.5d0/y2*(1.d0-1.5d0/y2*(1.d0-2.5d0/y2)))   !+105.d0/16.d0/b12**4)
   erfk=2.d0/pi12*som
   som6=som6-dx*e/y**4/2.d0+dx**2/12.d0*(2.d0/y**3+4.d0/y**5)*e
   som6=som6+e/2.d0/y**5*(1.d0-2.5d0/y2*(1.d0-3.5d0/y2*(1.d0-4.5d0/y2)))
   bet6=pi15*x**3/2.d0*som6
   som12=som12-dx*e/y**10/2.d0+dx**2/12.d0*(2.d0/y**9+10.d0/y**11)*e
   som12=som12+e/2.d0/y**11*(1.d0-5.5d0/y2*(1.d0-6.5d0/y2*(1.d0-7.5d0/y2)))
   bet12=pi15*x**9/120.d0*som12

end if
!
end subroutine
!
!


function erfc_Luc(X,x2,ee2) !DOC 
!DOC calcualte erfc(x)
!DOC Parameters:
!DOC x --> b
!DOC x2 --> b^2
!DOC ee2 --> exp(-b^2)
!DOC Return value:
!DOC erfc(b)
      IMPLICIT REAL(8) (a-h,o-z)
      parameter (P=0.3275911d0,A1=0.254829592d0,        &
       A2=-0.284496736d0,A3=1.421413741d0,A4=-1.453152027d0,   &
       A5=1.061405429d0,pi05=1.7724538509055d0)
!DOC        calcule erfc_Luc(x)=2/racine(pi) integrale de x à infini de exp(-t**2)dt
!DOC        x2=x**2 et ee2=exp(-x**2)
!DOC        luc85p108
      if(x<0.455d0) then                ! DL de erf si x petit
      erfc_Luc=1.d0-2.d0/pi05*x*(1.d0-x2*(1.d0/3.d0-x2*(0.1d0-x2*(1./42.d0-x2/216.d0))))
      ELSEif(x>3.d0)  then              ! Devlpt asymptotique si x grand
      ax=(1.d0-0.5d0*(1.d0-1.5d0/x2)/x2)/(pi05*x)
      erfc_Luc=ax*ee2
          else                          ! Abramowitz Stegun sinon
      T=1.D0/(1.D0+P*x)
      AX=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
      erfc_Luc=ax*ee2
          end if
      RETURN
      END 
!
!
function erf_Luc(X,x2,ee2) !DOC
!DOC Calculate the erf(X)
!DOC Parameters:
!DOC  X --> b
!DOC  x2 --> b2
!DOC ee2 --> exp(-b^2)
!DOC Return value:
!DOC  erf(b)
      IMPLICIT REAL(8) (a-h,o-z)
      parameter (P=0.3275911d0,A1=0.254829592d0,        &
       A2=-0.284496736d0,A3=1.421413741d0,A4=-1.453152027d0,   &
       A5=1.061405429d0,pi05=1.7724538509055d0)
!DOC        calcule erf_Luc(x)=2/racine(pi) integrale de 0 à x de exp(-t**2)dt
!DOC        x2=x**2 et ee2=exp(-x**2)
!DOC        luc85p108
      if(x<0.455d0) then                ! DL de erf si x petit
      erf_Luc=2.d0/pi05*x*(1.d0-x2*(1.d0/3.d0-x2*(0.1d0-x2*(1./42.d0-x2/216.d0))))
      ELSEif(x>3.d0)  then              ! Devlpt asymptotique de erfc_Luc=1-erf si x grand
      ax=(1.d0-0.5d0*(1.d0-1.5d0/x2)/x2)/(pi05*x)
      erf_Luc=1.d0-ax*ee2
          else                          ! Abramowitz Stegun pour erfc_Luc sinon
      T=1.D0/(1.D0+P*x)
      AX=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
      erf_Luc=1.d0-ax*ee2
          end if
      RETURN
      END 
!
!      FUNCTION erfc_Luc_pas_utilisee(X,x2,ee2)
!      IMPLICIT REAL(8) (a-h,o-z)
!      parameter (pi05=1.7724538509055d0)
!!
!!         fonction qui donne erfc_Luc
!!         Numerical Recipes
!!         fonction gamma Q(1/2,x**2)
!!         developpement limite si x petit
!!         developpement fraction continue si x grand
!!         luc85p106
!!
!epsi=1.d-14
!epsiprime=1.d-30
!!x2=x**2
!!ee2=EXP(-x2)
!!
!if (x<=1.5d0) then                  ! DL
!som=1.d0
!a=1.d0
!do i=1,100
!a=2.d0/(2.d0*i+1.d0)*x2*a
!som=som+a
!if(a<som*epsi) exit
!end do
!!PRINT*, 'DL',i
!erf=ee2*2.d0/pi05*x*som
!erfc_Luc=1.d0-erf
!!
!else                                ! fraction continue  methode de Lorentz section 5.2 NR
!!
!b=x2+1.d0-0.5d0
!c=1.d0/epsiprime
!d=1.d0/b
!h=d
!do i=1,100
!a=-i*(i-0.5d0)
!b=b+2.d0
!d=a*d+b; if(ABS(d)<epsiprime) d=epsiprime
!c=b+a/c; if(ABS(c)<epsiprime) c=epsiprime
!d=1.d0/d
!dc=d*c
!h=h*dc
!if(ABS(dc-1.d0)<epsi) exit
!end do
!!print*, 'FC',i
!erfc_Luc=ee2/pi05*x*h
!end if
!!
!return
!end


function shi_sans_exp_f(x) !DOC
IMPLICIT REAL(8) (a-h,o-z)
!
!DOC        calcul auto de Shi(x)
!DOC        luc84p163
!DOC        en fait, on ne veut pas de facteur exp(x) qui peut occasionner un overflow
!DOC        donc donner plutot Shi(x)/exp(x)
!
dx0=0.2d0
IF(x<=1.d0) then          !  DL x<1
   x2=x**2
   shi_f=x
   fac=1.d0
   xn=x

   do i=1,4
      deuxi=2.d0*i
      xn=xn*x2
      fac=fac*(deuxi)*(deuxi+1.d0)
      shi_f=shi_f+xn/(deuxi+1.d0)/fac
   end do

   shi_sans_exp_f=EXP(-x)*shi_f
   return
   !
elseif(x>20.d0) then            ! D asympt x>20
   som=1.d0
   term=1.d0

   do i=1,6
      term=term*i/x
      som=som+term
   end do

   !shi_f=EXP(x)/2.d0/x*som
   shi_sans_exp_f=1.d0/2.d0/x*som  !!! c'est l'interet!
   return
   !
ELSEIF(x<5.d0) then             ! x entre 1 et 5
   x1=1.d0
   shi_f=1.057250875d0
ELSEIF(x<10.d0) then            ! x entre 1 et 10
   x1=5.d0
   shi_f=20.093211826d0
ELSEIF(x<15.d0) then            ! x entre 10 et 15
   x1=10.d0
   shi_f=1246.11449d0
else                            ! x entre 15 et 20
   x1=15.d0
   shi_f=117477.926d0
endif
!                           integration numerique a partir du point repertorie precedent
n=(x-x1)/dx0+2
dx=(x-x1)/(n-1)
shi_f=shi_f+dx*SINH(x1)/x1/2.+dx**2/12.*(COSH(x1)/x1-SINH(x1)/x1**2)

do i=2,n
   t=x1+(i-1)*dx
   s=SINH(t)
   shi_f=shi_f+dx*s/t
end do

shi_f=shi_f-dx*s/t/2.d0-dx**2/12.d0*(COSH(t)/t-s/t**2)
shi_sans_exp_f=EXP(-x)*shi_f
end
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module
