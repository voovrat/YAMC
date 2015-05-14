Module MCLuc !DOC 
!DOC !FILE Module which includes the procedures written by Luc Belloni for the MC for water. Is used for the test purposes only (to check that my results are the same as Luc's).
 


!Type :: TMCLuc
! 
!   integer :: n
!
!   real(8) :: BoxLength
!   real(8) :: sigma2,sigma_a,xlj
!  
!   real(8),dimension(:),allocatable :: bet,bet6,bet12
! 
!   Type(TFourierGrid),pointer :: grid 
!
!   real(8),dimension(:),allocatable :: xx,yy,zz
!   real(8),dimension(:),allocatable :: uu,fx,fy,fz
!
!End Type TMCLuc

 Type ::  MCLucInput  !DOC
!DOC The input parameters for the Luc procedures
!DOC Fields:

   real(8) :: BoxLength !DOC Length of the box
 
   real(8) :: sigma_a   !DOC sigma in angtroems
   real(8) :: xlj  !DOC 4epsilon (in kT)
   real(8) :: charg_h  !DOC charge of the hydrogen
   real(8) :: roh_a,theta_d !DOC distance OH and the angle HOH in water
 
   integer :: kmax !DOC maximal value for |k|
   real(8) :: rmax2 !DOC maximal value for r^2
   real(8) :: alpha !DOC alpha*L
   real(8) :: ssr,ssk !DOC  sr, sk
   
   real(8) :: temp  !DOC temperature
   real(8) :: cdiel_ext !DOC dielectrical permutivity
   real(8) :: dbjr_a !DOC Bjerum length in Angstroems

   real(8) :: dr_a !DOC maximal shift in angstroems
   real(8) :: d_angle !DOC maximal rotation in rad

   real(8) :: xlambda_f,xlambda_c !DOC constants for force and torq bias

End Type MCLucInput 


contains

subroutine MCLuc_init_input(input) !DOC 
!DOC initialize the MCLucInput structure with the standard values
   use parameters
!DOC Parameters:
   Type(MCLucInput) :: input !DOC output: MCLucInput structure

!   input % rhh_a = 1.6282461761048d0
   input % roh_a =   1.d0  ! 0.999992804974116d0 
   input % theta_d = 109.d0 ! 109.001493760613d0 
   input % sigma_a = 3.16d0 
   input % xlj = 4d0 * 0.15d0 / kT_kcal_mol
   input % charg_h = 0.4238d0

   input % BoxLength = BoxLength
   input % temp = temp
   input % cdiel_ext = external_permutivity
   input % dbjr_a = dbjr_a
   
   input % kmax = kmax
   input % rmax2 = rmax2 
   input % ssr = sr
   input % ssk = sk
   input % alpha = alpha 

   input % dr_a = dr_a
   input % d_angle = d_angle 
   
   input % xlambda_f = xlambda_f
   input % xlambda_c = xlambda_c


end subroutine

!     Last change:  LB    4 Feb 2014    3:46 pm
subroutine initmc_Luc(input,xx,yy,zz,n_h2o,systematic_correction) !DOC
!DOC run Luc's initmc function
USE lecture_fichier
!use parameters
use constants
!use io, only : rite_real_array, write_xyz_array,io_open,io_close
!use SystemSettings, only : SYSTEM_STRING_LENGTH
!use io
      implicit REAL(8) (a-h,o-z)
!DOC Parameters:     
      Type(MCLucInput) :: input !DOC MCLucInput structure
      real(8), dimension(:),intent(in) :: xx,yy,zz  !DOC atom coordinates 
      integer, intent(in) :: n_h2o !DOC number of h2o molecules
!      real(8),intent(in) :: sigma2
!      real(8),intent(in) :: xlj ! 4 epsilon
      logical,intent(in),optional :: systematic_correction !DOC use or not the Luc's systematic correction (I don't have it so for the comparison it should be .FALSE.)
!      real(8),intent(in),optional :: ssr,ssk   

      integer * 4 nittot,naccep,ncumul
      REAL * 8 iacc,iacc1,iacc2,iacc3,iacc_o

! DP maintenant      real * 4 alea
      character * 50 nomfic,nomaaa,nomlec
      common/nombre/n
!      common/positi/xx(10000),yy(10000),zz(10000)        ! positions

      common/xenergy/utotk,utotk6,utotk12,uu_lj6_r,uu_lj12_r
   common/dipole/summux,summuy,summuz

      common/deplac/dr_a,d_angle
      COMMON/forcebias/xlambda_f,xlambda_c

      common/iteri/nittot,naccep,ncumul
      common/numeri/p,num
      common/energi/uu(10000),vrnew(10000),uutot
      COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r
      COMMON/forces/fx(10000),fy(10000),fz(10000)
      common/moyenn/uutotm,duutot,evtest,devtes,mmu
      COMMON/nb_proj/mnmax,ialpmax
      common/cumul_1/iacc1(2000),iacc2(2000),iacc3(2000),    &  !gOO, OH, HH
      gr1(2000),gr2(2000),gr3(2000)
      PARAMETER(ialpmaxx=1226)                    ! 1226 pour nmax=6!
      PARAMETER(n_omega=8)
      COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)
      COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)
      common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190
      common/conc/conc_h2o,xl_a
      common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d

      COMMON/lennardjones/epsi_lj,xlj,ucut

      COMMON/ewald/temp,dbjr,charg_h,xxclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
!      COMMON/ewald/charg_h,xxclb,alp_a1,alp,kkmax,cdiel_ext

      COMMON/ewald_s/ssr,ssk                            ! pour bien sauvegarder sk!
      common/pression/vir,vir1,vir2,vir02,vir12,vir22,hypervir,hypervir2
      COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82
      COMMON/sauveaaa/nomaaa
      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
      common/gracan/activ,pinsdes,x_n,x_n2           ! pour GC
      PARAMETER(nkmax=4661,kmaxmax=20)
      COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
       sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
       sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
       cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
       cckz(0:kmaxmax),sskz(0:kmaxmax)
      COMMON/sommeenk_lj/bet6(nkmax),bet12(nkmax),                     &         ! pour Ewald LJ en k
       sumcos_o(nkmax),sumsin_o(nkmax),sumcos1_o(nkmax),sumsin1_o(nkmax),        &
       sumcos2_o(nkmax),sumsin2_o(nkmax),sumcos3_o(nkmax),sumsin3_o(nkmax)
      COMMON/numeriq/pasq,nptq,mq
      COMMON/spectre/sq1(500),sq2(500),sq3(500)
      COMMON/spectre1/sk1(500),sk2(500),sk3(500)
COMMON/inpout/inp
!
!        on initalise le programme Monte-Carlo pour H2O
!        luc belloni 84p121
!        on suit les forces en parallele luc84p133
!        finalement, pas necessaire luc84p137
!
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        pour la version FB, uu et fx,y,z contiennent la partie r SANS intra  luc84p187
!
!       on lit au clavier          luc84p92    et 85p130

!
!if(inp==0) inp=5
!!
!    1 format(1x,a)
!    2 format(1x,a,$)
!    3 format(1x,78(1h*))
!      pi=4.d0*ATAN(1.d0)
!      if(ichoixaaa.eq.0) ichoixaaa=1
!!       pour ecrire les valeurs de aaa en fichier
!      nomaaa='qq'
!      write(*,3)
!      write(*,3)
!      write(*,1) '          Monte-Carlo "luc" pour H2O"'
!      WRITE(*,1) '                   Modele SPCE'
!      WRITE(*,1) '                  LJ NON tronque!'
!      WRITE(*,1) '     utilise par des personnes jeunes, belles et dynamiques (et curieuses)'
!      write(*,3)
! 9999 write(*,3)
!      write(*,1) '                    MENU'
!      WRITE(*,1) '100--NNNEEEWWW!!! Clavier remplace par fichier d''entree'
!      write(*,1) '1--Parametres physiques'
!      write(*,1) '2--Positions de depart sur un reseau cfc'
!      write(*,1) '3--Lire en fichier les positions de depart'
!      write(*,1) '4--Parametres numeriques d''analyse'
!      write(*,1) '5--Lire en fichier les resulats deja accumules'
!      WRITE(*,1) '9--Variable pour estimer l''efficacite statistique'
!      WRITE(*,1) '10--Sauver les valeurs instantanees de la variable'
!      write(*,1) '11--"Lisser" des resultats accumules'
!      WRITE(*,1) '12--Cumuler 2 simulations'
!      write(*,2) '8--Demarrer les cycles Monte-Carlo  --->'
!      read(inp,*,err=9999) ic
!      write(*,3)
!!
!!
!      if(ic==100) then
!   !
!   !      on peut lire tous les paramètres en fichier           luc84p92
!   !
!         WRITE(*,2) 'Nom du fichier (pas de lecture si "qq") --->'
!         read(inp,'(A)') nomlec
!         if(nomlec=='stop') stop
!   
!         if(nomlec.ne.'qq'.and.nomlec.ne.'QQ') then
!              OPEN(8,FILE=nomlec,STATUS='old',ERR=9999)
!              inp=8
!         ELSE
!              CLOSE(8)
!              inp=5
!         END if ! not QQ
!
!      END if ! ic == 100
!!
!!
!!
!      goto(1000,2000,3000,4000,5000,9999,9999,8000,9600,9700,1100,1200),ic
!!        ic
!! 1000 : 1           1--Parametres physiques'
!! 2000 : 2           2--Positions de depart sur un reseau cfc'
!! 3000 : 3           3--Lire en fichier les positions de depart'
!! 4000 : 4           4--Parametres numeriques d''analyse'
!! 5000 : 5           5--Lire en fichier les resulats deja accumules'
!! 9999 : 6            -none-
!! 9999 : 7            -none-
!! 8000 : 8           8--Demarrer les cycles Monte-Carlo  --->'
!! 9600 : 9           9--Variable pour estimer l''efficacite statistique'
!! 9700 : 10          10--Sauver les valeurs instantanees de la variable'
!! 1100 : 11          11--"Lisser" des resultats accumules
!! 1200 : 12          12--Cumuler 2 simulations
!GOTO 9999
!!
!!        physique
!!
! 1000 write(*,*) 'Potentiel LJ TOTAL pour OO et charges partielles en O et H'
!write(*,3)
!write(*,*) 'La simulation Monte-Carlo va mettre n_h2o H2O dans un cube de taille L (qui sera prise comme unite de longueur)'
!write(*,3)
!write(*,2) 'Nombre de molecules d''eau n_h2o  --->'
!call lis_f(inp,n_h2o)
!write(*,2) 'Temperature (en K) --->'
!call lis_f(inp,temp)
!write(*,2) 'Concentration d''eau  (en M) --->'
!call lis_f(inp,conc_h2o)
!WRITE(*,3)
!WRITE(*,2) 'Diametre LJ de OO (en A) --->'
!call lis_f(inp,sigma_a)
!write(*,2) 'epsilon/kT LJ de OO --->'
!call lis_f(inp,epsi_lj)
!WRITE(*,2) 'Distance rOH (en A) et angle HOH (en degres) --->'
!call lis_f(inp,roh_a,theta_d)
!WRITE(*,2) 'Charge partielle de H --->'
!call lis_f(inp,charg_h)
!WRITE(*,3)
!WRITE(*,2) 'Coefficient K (en A-1) dans la decomposition Ewald (KL>>1) --->'
!call lis_f(inp,alp_a1)
!WRITE(*,2) 'On coupe les sommes en r et k de Ewald aux precisions sr,sk --->'
!call lis_f(inp,ssr,ssk)
!WRITE(*,2) 'Constante dielectrique exterieure --->'
!call lis_f(inp,cdiel_ext)
!WRITE(*,3)
!WRITE(*,2) 'NOUVEAU: Valeur maxi pour m,n dans le calcul des gmnl(r)  --->'
!call lis_f(inp,mnmax)
!if(mnmax==0) ialpmax=1
!if(mnmax==1) ialpmax=4
!if(mnmax==2) ialpmax=27
!if(mnmax==3) ialpmax=79
!if(mnmax==4) ialpmax=250
!if(mnmax==5) ialpmax=549
!if(mnmax==6) ialpmax=1226
!33 write(*,3)
!!        On calcule les parametres de la simulation avec L pris comme unite
!n=3*n_h2o
!avog=6.02214d-4
!!  C = N/ (NA*V) ==> V = N/ (C*NA)
!!  L = V^(1/3) 
!xl_a=(n_h2o/(conc_h2o*avog))**(1.d0/3.d0)   ! xl_a  - box length
!!  also: concentration is in mol/L = mol/dm^3
!!  so, answer is in dm, e.g.  L = ( N/ (C*NA))^(1/3) [dm]
!!  to convert from dm [10^-1 m ] to angstroms [10^-10 m] we multiply by  10^9 (or divide by 10^-9)
!!  this equivalently can be done inside the brackets (before taking ^(1/3)), which will give 10^-27 in the denominator
!!  combining 6.023e23 * 10^-27 we get 6.023e-4 (see the constant above)
!
!sigma=sigma_a/xl_a
!sigma2=sigma**2
!xlj=4.d0*epsi_lj                        ! donc v/kT=xlj*((sigma/r)**12-(sigma/r)**6)

   !rhh_a = input % rhh_a 
   roh_a = input % roh_a 
   theta_d = input % theta_d  
   sigma_a = input % sigma_a 
   xlj = input % xlj 
   charg_h = input % charg_h 

   xl_a = input % BoxLength 
   temp = input % temp 
   cdiel_ext = input % cdiel_ext 

   dbjr_a = input % dbjr_a 
   dbjr = dbjr_a / xl_a 
   
   kmax = input % kmax 
   rmax2 = input % rmax2 
   ssr = input % ssr 
   ssk = input % ssk
   alp = input % alpha
!   write(*,*) 'alp',alp,'input % alpha',input % alpha
   alp_a1 = alp / xl_a 

   sigma = sigma_a / xl_a
   sigma2 = sigma**2


   dr_a = input % dr_a 
   dr = dr_a / xl_a
   d_angle = input % d_angle 
   
   xlambda_f = input % xlambda_f
   xlambda_c = input % xlambda_c


rcutoff_a=1.d20                         ! sert a rien mais sera dans les fichiers acc
!roh_a = 1.d0
!theta_d = 109.d0
roh=roh_a/xl_a

!pi=4.d0*ATAN(1.d0)
!boltz=1.38065d-23
!elec=1.602176d-19
!epsi0=8.854187d-12   ! dialectrical permutivity of vacuum [ in SI ]
!dbjr_a=elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10    ! Bjerum length  en A
!dbjr=dbjr_a/xl_a
xxclb=charg_h**2*dbjr
theta=pi/180.d0*theta_d
theta2=theta/2.d0
cthet2=COS(theta2); sthet2=SIN(theta2)
rhh=2.d0*roh*sthet2

   write(*,*) 'xxclb',xxclb,'xl_a',xl_a,'alp',alp,'sigma2',sigma2,'charg_h',charg_h,'xlj',xlj,'cdiel_ext',cdiel_ext 

epsi_lj = xlj/4.0
ucut = 0.


!!u_self=xxclb*(-4.d0/roh+1.d0/rhh)
!!foh_self=-2.d0*xxclb/roh**2
!!fhh_self=+xxclb/rhh**2
!p=p_a/xl_a    ! p_a is read in section 4000 ( Parametres numeriques d''analyse )
!              ! it is spacing for g(r) 
!
!pasq=pasq_a1*xl_a ! also in 4000. It is spacing in q-space for diffusion spectra
!
!
!alp=alp_a1*xl_a   ! coefficient in the ewald decomposition (see above) 
pi2=2.d0*pi 
!pialp2=pi**2/alp**2
!alp2pi=2.d0*alp/dSQRT(pi)
!rmax=ssr/alp
!rmax2=rmax**2
!kmax=NINT(alp*ssk/pi)
!ssk_vrai=kmax*pi/alp
!xkmax2=kmax**2
!ucut=xlj*((sigma/rmax)**12-(sigma/rmax)**6)     ! si on s'amuse a revenir a LJ tronque shifte a rmax
!PRINT*, 'Nombre de sites en tout : ',n
!PRINT*, 'Taille de la boite L (en A) : ',xl_a
!PRINT*, 'Diametre en unite L : ',sigma
!PRINT*, 'Distance OH en unite L : ',roh
!      WRITE(*,*) 'KL (doit etre >>1) : ',alp
!      WRITE(*,*) 'rmax/L, kmax, sk reel pour Ewald= ',rmax,kmax,ssk_vrai
!if(kmax>kmaxmax) PRINT*, 'ATTENTION: kmax dépasse la taille des tableaux!!!!'
!PRINT*, 'Nombre de projections : ',ialpmax
!!
!!             attention: si on a change la concentration..., L a change
!!             et donc les positions absolues venant des memes positions relatives
!!             ont ete implicitement renormalisées
!!             bien pour les centres, faux pour les distances intra!!!
!!             donc corriger  luc81p31
!!
!if(ilec23==17) then
!do i_h2o=1,n_h2o
!i0=(i_h2o-1)*3+1
!do i=i0+1,i0+2
!xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
!rcs=SQRT(xcs**2+ycs**2+zcs**2)
!xx(i)=xx(i0)+roh*xcs/rcs
!yy(i)=yy(i0)+roh*ycs/rcs
!zz(i)=zz(i0)+roh*zcs/rcs
!end do
!end do
!end if
!!
!      ilec1=17
!      goto 9999
!!
!!        position initiale des particules sur un reseau cfc pour les coeurs
!!        et orientations au hasard pour les sites
!!        d'abord H2O pointe en haut puis 3 rotations aleatoires suivant les 3 axes
!!
! 2000 if(n==0) goto 9999
!      xn=n_h2o
!      m=(xn/4.)**(1./3.)+0.999
!      xl1=1./m
!      write(*,3)
!      write(*,*) 'La cellule cubique est decoupee en ',m**3,' mailles'
!      write(*,*) 'de cote ',xl1
!      write(*,*) 'Les Oxygènes vont etre distribuees sur les ',4*m**3,        &
!               ' sites du reseau cfc ainsi forme'
!      i=1
!        do 100 ix=1,m
!        do 100 iy=1,m
!        do 100 iz=1,m
!        xx(i)=-0.5+(ix-1)*xl1
!        yy(i)=-0.5+(iy-1)*xl1
!        zz(i)=-0.5+(iz-1)*xl1
!        i=i+3
!  100   if(i>n) goto 101
!        do 104 ix=1,m
!        do 104 iy=1,m
!        do 104 iz=1,m
!        xx(i)=-0.5+xl1/2.+(ix-1)*xl1
!        yy(i)=-0.5+(iy-1)*xl1
!        zz(i)=-0.5+xl1/2.+(iz-1)*xl1
!        i=i+3
!  104   if(i>n) goto 101
!        do 105 ix=1,m
!        do 105 iy=1,m
!        do 105 iz=1,m
!        xx(i)=-0.5+xl1/2.+(ix-1)*xl1
!        yy(i)=-0.5+xl1/2.+(iy-1)*xl1
!        zz(i)=-0.5+(iz-1)*xl1
!        i=i+3
!  105   if(i>n) goto 101
!        do 106 ix=1,m
!        do 106 iy=1,m
!        do 106 iz=1,m
!        xx(i)=-0.5+(ix-1)*xl1
!        yy(i)=-0.5+xl1/2.+(iy-1)*xl1
!        zz(i)=-0.5+xl1/2.+(iz-1)*xl1
!        i=i+3
!  106   if(i>n) goto 101
!  101 continue
!!       orientations au hasard
!PRINT*, 'Et orientations au hasard de la molécule H2O'
!!        initialisation du generateur de nombres aleatoires
!!      x=alea(0)
!do i=1,n_h2o
!vx1=sthet2                           ! d'abord pointe en haut
!vy1=0.
!vz1=cthet2
!vx2=-sthet2
!vy2=0.
!vz2=cthet2
!!angle=pi2*alea(1)               ! rotation aleatoire suivant x
!angle = pi2 
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,1)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,1)
!!angle=pi2*alea(1)               ! rotation aleatoire suivant y
!!angle = pi2
!angle = pi2
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,2)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,2)
!!angle=pi2*alea(1)               ! rotation aleatoire suivant z
!angle = pi2
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,3)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,3)
!i0=3*(i-1)+1
!xx(i0+1)=xx(i0)+roh*vx1
!yy(i0+1)=yy(i0)+roh*vy1
!zz(i0+1)=zz(i0)+roh*vz1
!xx(i0+2)=xx(i0)+roh*vx2
!yy(i0+2)=yy(i0)+roh*vy2
!zz(i0+2)=zz(i0)+roh*vz2
!end do
!!    attention, certains sites peuvent etre a l'exterieur de la boite
!!    maintenant, je m'en fiche! luc80p179
!      ilec23=17
!      goto 9999
!!
!!        lecture en fichier des positions
!!
! 3000 continue
!      write(*,2) 'nom du fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=3000,file=nomfic,STATUS='old')
!      read(7,*,end=3000) n
!      write(*,*) 'n= ',n
!      read(7,*) (xx(i),yy(i),zz(i),i=1,n)
!      close(7)
!!              il peut arriver que cette config (ancienne) n'ait pas autorise les sites H
!!              a sortir legerement de la boite, donc les ait deplaces
!!              de l'autre cote
!!              maintenant, on force les sites a rester autour de leur coeur O
!!              quitte a sortir de la boite
!do i0=0,n-1,3
!do j=2,3
!xr=xx(i0+j)-xx(i0+1)
!yyr=yy(i0+j)-yy(i0+1)
!zr=zz(i0+j)-zz(i0+1)
!xr=xr-NINT(xr)
!yyr=yyr-NINT(yyr)
!zr=zr-NINT(zr)
!xx(i0+j)=xx(i0+1)+xr
!yy(i0+j)=yy(i0+1)+yyr
!zz(i0+j)=zz(i0+1)+zr
!end do
!end do
!!
!n_h2o=n/3
!conc_h2o=n_h2o/xl_a**3/avog
!PRINT*, 'Nombre total de sites : ',n
!PRINT*, 'Nombre de H2O : ',n_h2o
!      write(*,*) 'Concentration de H2O (en M) : ',conc_h2o
!      end if
!!
!!             attention: si on a change la concentration..., L a change
!!             et donc les positions absolues venant des memes positions relatives
!!             ont ete implicitement renormalisées
!!             bien pour les centres, faux pour les distances intra!!!
!!             donc corriger  luc81p31
!!
!do i_h2o=1,n_h2o
!i0=(i_h2o-1)*3+1
!do i=i0+1,i0+2
!xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
!rcs=SQRT(xcs**2+ycs**2+zcs**2)
!xx(i)=xx(i0)+roh*xcs/rcs
!yy(i)=yy(i0)+roh*ycs/rcs
!zz(i)=zz(i0)+roh*zcs/rcs
!end do
!end do
!!
!      ilec23=17
!      goto 9999
!!
!!        parametres d'analyse pour la probabilite g(r)
!!
! 4000 if(ilec1.ne.17) goto 9999
!      write(*,2) 'pour g(r), nombre d''intervalles --->'
!      call lis_f(inp,num)
!      WRITE(*,2) 'et pas (en A) --->'
!      call lis_f(inp,p_a)
!      p=p_a/xl_a
!      write(*,*) 'Pas en unite L : ',p
!      WRITE(*,2) 'pour g(r,omega), pas en cos(beta) et pas en phi (en rd) --->'      ! g(r,omega)  luc85p176
!      call lis_f(inp,dcosbet,dphi)
!      WRITE(*,2) 'et coefficients correcteurs pour A=C=E=G et B=D=F=H  --->'
!      call lis_f(inp,coeff_a,coeff_b)
!      write(*,1) 'Pour le potentiel chimique, on fait des insertions de particules tests M fois au hasard'
!      write(*,2) 'Valeur de M --->'
!      call lis_f(inp,mmu)
!      WRITE(*,*) 'Pour le calcul des spectres de diffusion :'
!      WRITE(*,2) 'Nombre de points en q --->'
!      CALL lis_f(inp,nptq)
!      WRITE(*,2) 'Pas en q (A-1) --->'
!      call lis_f(inp,pasq_a1)
!      pasq=pasq_a1*xl_a
!      WRITE(*,2) 'Nombre d''orientations de q a chaque config --->'
!      call lis_f(inp,mq)
!  44  WRITE(*,2) 'On efface tout? (oui=1) --->'
!      READ(inp,*,ERR=44) ieff
!      if(ieff==1) then
!      iacc(1:ialpmax,1:num)=0.
!      iacc1(1:num)=0.; iacc2(1:num)=0.; iacc3(1:num)=0.
!      nittot=0
!      naccep=0
!      ncumul=0
!      uutotm=0.
!      duutot=0.
!      evtest=0.
!      devtes=0.
!      vir=0.
!      vir1=0.
!      vir2=0.
!      vir02=0.
!      vir12=0.
!      vir22=0.
!      hypervir=0.
!      hypervir2=0.
!      xmtotx=0.; xmtoty=0.; xmtotz=0.; xmtot2=0.; xmtot22=0.
!      do kq=1,nptq
!      sq1(kq)=0.
!      sq2(kq)=0.
!      sq3(kq)=0.
!      end do
!      iacc_o(:,1:num)=0.
!      end if
!      ilec4=17
!      goto 9999
!!
!!        lecture en fichier des resultats
!!
! 5000 continue
!      write(*,2) 'nom du fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=5000,status='old',file=nomfic)
!      read(7,*,err=5000) n_h2o
!      write(*,*) 'Nombre de H2O : ',n_h2o
!      read(7,*) temp
!      write(*,*) 'Temperature (en K) : ',temp
!      read(7,*) conc_h2o
!      write(*,*) 'Concentration (M) : ',conc_h2o
!      READ(7,*) sigma_a
!      WRITE(*,*) 'Diametre LJ OO (A) : ',sigma_a
!      READ(7,*) epsi_lj
!      WRITE(*,*) 'Epsilon/kT LJ OO : ',epsi_lj
!      READ(7,*) rcutoff_a
!      READ(7,*) roh_a,theta_d
!      WRITE(*,*) 'Distance OH (en A) et angle HOH (en degres) : ',roh_a,theta_d
!      READ(7,*) charg_h
!      WRITE(*,*) 'Charge partielle de H : ',charg_h
!      READ(7,*) alp_a1
!      PRINT*, 'Coefficient K Ewald (en A-1) : ',alp_a1
!      READ(7,*) ssr,ssk
!      PRINT*, 'Precisions r et k Ewald : ',ssr,ssk
!      READ(7,*) cdiel_ext
!      PRINT*, 'Constante dielectrique externe : ',cdiel_ext
!      READ(7,*) mnmax
!      PRINT*, 'Valeur max de m et n dans les projections : ',mnmax
!if(mnmax==0) ialpmax=1
!if(mnmax==1) ialpmax=4
!if(mnmax==2) ialpmax=27
!if(mnmax==3) ialpmax=79
!if(mnmax==4) ialpmax=250
!if(mnmax==5) ialpmax=549
!if(mnmax==6) ialpmax=1226
!PRINT*, 'Nombre de projections : ',ialpmax
!      read(7,*) num,p_a
!      write(*,*) 'Nombre d''intervalles et pas (en A) = ',num,p_a
!      read(7,*) nittot,naccep,ncumul
!      write(*,*) 'nombre d''iterations deja effectuees : ',nittot
!      write(*,*) 'nombre d''accumulations deja effectuees : ',ncumul
!      read(7,*) uutotm,duutot
!      read(7,*) vir,vir1,vir2,vir02,vir12,vir22
!      READ(7,*) hypervir,hypervir2
!      READ(7,*) xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
!      read(7,*) evtest,devtes,mmu
!      WRITE(*,*) 'Nombre d''insertions par configuration : ',mmu
!      do i=1,num
!      READ(7,*,END=5010) r,iacc1(i),iacc2(i),iacc3(i)    ! r ne sert pas
!      end do
!      do i=1,num
!      READ(7,*,END=5010) r,iacc(1:ialpmax,i)
!      end do
!      READ(7,*,END=5010) nptq,pasq_a1,mq
!      WRITE(*,*) 'En q: ',nptq,pasq_a1,mq
!      do kq=1,nptq
!      READ(7,*,END=5010) q,sq1(kq),sq2(kq),sq3(kq)
!      end do
!!         normalisation des Sij(q)
!      den=ncumul
!      n_o=n/3
!      n_h=n-n_o
!      do kq=1,nptq
!      sk1(kq)=sq1(kq)/(den*n_o)
!      sk2(kq)=sq2(kq)/(den*SQRT(n_o*n_h*1.))
!      sk3(kq)=sq3(kq)/(den*n_h)
!      end do
!      READ(7,*,END=5010) dcosbet,dphi                     ! g(r,omega) luc85p176
!      WRITE(*,*) 'pas en cosbeta et phi pour g(r,omega) : ',dcosbet,dphi
!      READ(7,*,END=5010) coeff_a,coeff_b
!      WRITE(*,*) 'coefficients correcteurs en A et en B : ',coeff_a,coeff_b
!      do i=1,num
!      READ(7,*,END=5010) r,iacc_o(1:n_omega,i)    ! r ne sert pas
!      end do
!!PRINT*, 'je divise par 2 les cas G et H'           ! juste pour rattraper les fichiers 48 et 49
!!iacc_o(7:8,1:num)=0.5d0*iacc_o(7:8,1:num)
!      GOTO 5011
! 5010 PRINT*, 'ATTENTION!!!!!! fin de fichier detectee!'
! 5011 continue
!      close(7)
!      ilec4=17
!      goto 33
!        end if
!      GOTO 9999
!!
!!      pour estimer l'efficacite statistique s, on a le choix entre differentes variables
!!
! 9600 WRITE(*,*) 'On a le choix entre les variables suivantes : '
!      WRITE(*,*) '1--Energie  2--viriel  '
!      WRITE(*,*) '6--exp(-vtest)  7--hyperviriel 8--M**2  9,10,11--Mx,y,z'
!      WRITE(*,*) '10000*alpha+i pour g_alpha(ri)'                 ! et non 1000*alpha... luc82p166
!      WRITE(*,*) 'Le choix actuel est : ',ichoixaaa
!      WRITE(*,2) 'Nouveau choix souhaite --->'
!      read(inp,*,ERR=9999) nent
!      ichoixaaa=nent
!      GOTO 9999
!!
!!      on peut sauvagarder en fichier les valeurs instantanees de la variable
!!
! 9700 WRITE(*,2) 'Nom du fichier de sauvegarde ("qq" pour ignorer) --->'
!      READ(inp,'(A)') nomaaa
!      GOTO 9999
!!
!!       on regroupe des intervalles pour lisser
!!
! 1100 if(ilec4.ne.17) goto 9999
!      write(*,*) 'Actuellement, le nombre d''intervalles vaut ',num
!      WRITE(*,*) 'et le pas (en A) ',p_a
!      write(*,2) 'On va les regrouper par paquets de combien ? --->'
!      call lis_f(inp,ipaq)
!      num1=num/ipaq
!      do k=1,num1
!      iacc(:,k)=iacc(:,(k-1)*ipaq+1)
!      iacc1(k)=iacc1((k-1)*ipaq+1)
!      iacc2(k)=iacc2((k-1)*ipaq+1)
!      iacc3(k)=iacc3((k-1)*ipaq+1)
!      iacc_o(:,k)=iacc_o(:,(k-1)*ipaq+1)
!       do l=2,ipaq
!       iacc(:,k)=iacc(:,k)+iacc(:,(k-1)*ipaq+l)
!       iacc1(k)=iacc1(k)+iacc1((k-1)*ipaq+l)
!       iacc2(k)=iacc2(k)+iacc2((k-1)*ipaq+l)
!       iacc3(k)=iacc3(k)+iacc3((k-1)*ipaq+l)
!       iacc_o(:,k)=iacc_o(:,k)+iacc_o(:,(k-1)*ipaq+l)
!       END do
!      END do
!      num=num1
!      p_a=p_a*float(ipaq)
!      p=p_a/xl_a
!      WRITE(*,*) 'Nouveaux nombre de points et pas: ',num,p_a
!      goto 9999
!!
!!        on cumule les donnees actuelles avec celles d'un autre fichier luc80p127
!!
! 1200 if(ilec4.ne.17) goto 9999
!      WRITE(*,1) 'On cumule les donnees actuelles avec celles lues dans un autre fichier'
!      WRITE(*,1) 'donnees = alpha * actuelles + beta * lues'
!      write(*,2) 'nom du nouveau fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=9999,status='old',file=nomfic)
!      WRITE(*,2) 'Parametres de melange alpha et beta --->'
!      call lis_f(inp,alpha,beta)
!      PRINT*, 'ATTENTION: VERifIER que les 2 fichiers correspondent au meme systeme!'
!      read(7,*,err=5000) n_h2o
!      write(*,*) 'Nombre de H2O : ',n_h2o
!      read(7,*) temp
!      write(*,*) 'Temperature (en K) : ',temp
!      read(7,*) conc_h2o
!      write(*,*) 'Concentration de H2O (en M) : ',conc_h2o
!      READ(7,*) sigma_a
!      WRITE(*,*) 'Diametre LJ OO (A) : ',sigma_a
!      READ(7,*) epsi_lj
!      WRITE(*,*) 'Epsilon/kT LJ OO : ',epsi_lj
!      READ(7,*) rcutoff_a
!      READ(7,*) roh_a,theta_d
!      WRITE(*,*) 'Distance OH (en A) et angle HOH (en degres) : ',roh_a,theta_d
!      READ(7,*) charg_h
!      WRITE(*,*) 'Charge partielle de H : ',charg_h
!      READ(7,*) alp_a1
!      PRINT*, 'Coefficient K Ewald (en A-1) : ',alp_a1
!      READ(7,*) ssr,ssk
!      PRINT*, 'Precisions r et k Ewald : ',ssr,ssk
!      READ(7,*) cdiel_ext
!      PRINT*, 'Constante dielectrique externe : ',cdiel_ext
!      READ(7,*) mnmax
!      PRINT*, 'Valeur max de m et n dans les projections : ',mnmax
!      read(7,*) num,p_a
!      write(*,*) 'Nombre d''intervalles et pas (en A) = ',num,p_a
!      read(7,*) nittot_1,naccep_1,ncumul_1
!      write(*,*) 'nombre d''iterations deja effectuees : ',nittot_1
!      write(*,*) 'nombre d''accumulations deja effectuees : ',ncumul_1
!      nittot=NINT(alpha*nittot+beta*nittot_1)
!      naccep=NINT(alpha*naccep+beta*naccep_1)
!      x=alpha*ncumul; x_1=beta*ncumul_1
!      ncumul=nint(alpha*ncumul+beta*ncumul_1)
!      x=x/ncumul; x_1=x_1/ncumul
!      read(7,*) uutotm_1,duutot_1
!      read(7,*) vir_1,vir1_1,vir2_1,vir02_1,vir12_1,vir22_1
!      READ(7,*) hypervir_1,hypervir2_1
!      READ(7,*,end=1210) xmtotx_1,xmtoty_1,xmtotz_1,xmtot2_1,xmtot22_1
!      uutotm=alpha*uutotm+beta*uutotm_1
!      duutot=alpha*duutot+beta*duutot_1
!      vir=alpha*vir+beta*vir_1
!      vir1=alpha*vir1+beta*vir1_1
!      vir2=alpha*vir2+beta*vir2_1
!      vir02=alpha*vir02+beta*vir02_1
!      vir12=alpha*vir12+beta*vir12_1
!      vir22=alpha*vir22+beta*vir22_1
!      hypervir=alpha*hypervir+beta*hypervir_1
!      hypervir2=alpha*hypervir2+beta*hypervir2_1
!      xmtotx=alpha*xmtotx+beta*xmtotx_1
!      xmtoty=alpha*xmtoty+beta*xmtoty_1
!      xmtotz=alpha*xmtotz+beta*xmtotz_1
!      xmtot2=alpha*xmtot2+beta*xmtot2_1
!      xmtot22=alpha*xmtot22+beta*xmtot22_1
!      read(7,*) evtest_1,devtes_1,mmu_1
!      WRITE(*,*) 'Nombre d''insertions par configuration : ',mmu_1
!!         attention: on peut envisager que mmu et mmu_1 different luc82p163
!!         rappel: evtest est deja normalise par mmu
!      if(mmu>0.or.mmu_1>0) then
!      xmu_moyen=x*mmu+x_1*mmu_1
!      evtest=(alpha*mmu*evtest+beta*mmu_1*evtest_1)/xmu_moyen
!      devtes=(alpha*mmu*devtes+beta*mmu_1*devtes_1)/xmu_moyen
!      end if
!!         puis on prend pour mmu la nouvelle valeur qui vient d'etre lue
!      mmu=mmu_1
!!       on se sert de gr comme tableau auxi
!      do i=1,num
!      read(7,*) r,gr(1:3,i)
!      end do
!      iacc1(1:num)=NINT(alpha*iacc1(1:num)+beta*gr(1,1:num))
!      iacc2(1:num)=NINT(alpha*iacc2(1:num)+beta*gr(2,1:num))
!      iacc3(1:num)=NINT(alpha*iacc3(1:num)+beta*gr(3,1:num))
!      do i=1,num
!      read(7,*) r,gr(1:ialpmax,i)
!      END do
!      iacc(1:ialpmax,1:num)=alpha*iacc(1:ialpmax,1:num)+beta*gr(1:ialpmax,1:num)
!      READ(7,*,END=1210) nptq,pasq_a1,mq
!      WRITE(*,*) 'En q: ',nptq,pasq_a1,mq
!      do kq=1,nptq
!      READ(7,*,END=1210) q,sk1(kq),sk2(kq),sk3(kq)
!      end do
!      sq1(1:nptq)=alpha*sq1(1:nptq)+beta*sk1(1:nptq)
!      sq2(1:nptq)=alpha*sq2(1:nptq)+beta*sk2(1:nptq)
!      sq3(1:nptq)=alpha*sq3(1:nptq)+beta*sk3(1:nptq)
!!         normalisation des Sij(q)
!      den=ncumul
!!      n_o=n/3
!!      n_h=n-n_o
!!      do kq=1,nptq
!!      sk1(kq)=sq1(kq)/(den*n_o)
!!      sk2(kq)=sq2(kq)/(den*SQRT(n_o*n_h*1.))
!!      sk3(kq)=sq3(kq)/(den*n_h)
!      end do
!      READ(7,*,END=1210) dcosbet,dphi
!      WRITE(*,*) 'pas en cosbeta et phi pour g(r,omega) : ',dcosbet,dphi
!!      READ(7,*,END=5010) coeff_a,coeff_b
!      READ(7,*) coeff_a,coeff_b
!      WRITE(*,*) 'coefficients correcteurs en A et en B : ',coeff_a,coeff_b
!      do i=1,num
!      read(7,*,END=1210) r,gr(1:n_omega,i)
!      END do
!      iacc_o(1:n_omega,1:num)=alpha*iacc_o(1:n_omega,1:num)+beta*gr(1:n_omega,1:num)
!      GOTO 1211
! 1210 PRINT*, 'ATTENTION: fin de fichier detectee!'
! 1211 continue
!      close(7)
!      GOTO 33
!        end if
!      GOTO 9999
!!
!!      on demarre les cycles si on a bien tout defini
!!      juste avant, on initialise les energies et les nbrs aleatoires
!!
! 8000 if(ilec1.ne.17.or.ilec23/=17.or.ilec4/=17) goto 9999
!
!        on initialise les energies et les forces en r SANS intra luc84p187
!        on calcule maintenant LJ total par Ewald aussi!
!
!     charge_h = 0.42380d0
!     xxclb = xclb * charge_h**2

      alp2pi=2.d0*alp/dSQRT(pi)
      n = n_h2o * 3

!      write(*,*) 'n=',n

      uu(1:n)=0.                         ! repere LJr+EWr
      fx(1:n)=0.; fy(1:n)=0.; fz(1:n)=0.
      uu_lj12=0.; uu_lj6=0.
      uu_ew=0.
      summux=0.;summuy=0.; summuz=0.
      
      !        pour test, je reviens a LJ tronque shifte a rcut=rmax  luc85p123   NON, plus maintenant
      !xlj0=xlj; xlj=0.
      
!       EWALD SUM
!
!       U = U_real + U_fourier - U_self
!
!       where 
!
!       U_real = SUM_i SUM_j<>i  q_i q_j erfc( alpha r_ij ) / r_ij 
!              (the first term. Other images are ommited, and as r_ij is taken the nearliest image)
!
!
!       U_fourier = 1 / 2V SUM_m<>0   4pi/k_m^2 SUM_i SUM_j q_i q_j exp(i k_m r_ij) exp( - k_m^2/4 alpha^2)
        
!       U_self = sqrt(alpha^2/pi) SUM_i q_i^2    


      ! calculation of potential (Real part of Ewald sum )

      do i=1,n
         iio=3*((i-1)/3)+1   ! equivalent to mod(i,3)==1 : 
                            ! i == 3k+1  --> (i-1)/3 = k,   iio = 3*k + 1 == i
                            ! i == 3k+2 --> (i-1)/3 = k     iio = 3*k + 1 <> i 
         chargi=+1.d0       ! i == 3k   --> (i-1)/3 = k-1   iio = 3(k-1) + 1 = 3k-2 <> i
          if(i==iio) chargi=-2.d0 ! 

         summux=summux+chargi*xx(i)                      ! rq: j'ignore ici que H peut sortir de la boite
         summuy=summuy+chargi*yy(i)
         summuz=summuz+chargi*zz(i)

         do j=i+1,n

            jo=3*((j-1)/3)+1
            chargj=+1.d0

            if(j==jo) chargj=-2.d0

            xij=xx(i)-xx(j)
            yij=yy(i)-yy(j)
            zij=zz(i)-zz(j)
!            xij=xij-anint(xij)    ! anint =  round  
!            yij=yij-anint(yij)    ! Coordinates are in "box" units (1=length of the box)
!            zij=zij-anint(zij)    ! i.e. x - round(x) = nearliest image
         if(xij.gt. 0.5) xij=xij-1.   ! nearliest neighbour. Remember: distance units are BoxLength
         if(xij.lt.-0.5) xij=xij+1.   ! i.e. Box is [-0.5:0.5]^3
 
         if(yij.gt. 0.5) yij=yij-1.
         if(yij.lt.-0.5) yij=yij+1.         ! prend moins de temps que anint!
 
         if(zij.gt. 0.5) zij=zij-1.
         if(zij.lt.-0.5) zij=zij+1.
                                  ! only one image is used ( efrc sum has only 1 term ) 

            r2=xij**2+yij**2+zij**2   ! distanse to the nearliest image


!            write(*,*) 'xij=',xij,'yij=',yij,'zij=',zij

            if(r2<rmax2) then                         ! rmax pour Ewald r   rmax = ssr / alp 
               r=SQRT(r2)
               alpr=alp*r
               alpr2=alpr**2
               ee2=EXP(-alpr2)
               er=erfc_Luc(alpr,alpr2,ee2)            ! calcule erfc_Luc(x) en s'aidant de x**2 et exp(-x**2)
 !               write(*,*) 'alpr=',alpr,'akpr2=',alpr2,'ee2=',ee2,'er=',er,'potew=',potew

            
               if(iio==jo) er=-erf_Luc(alpr,alpr2,ee2)    ! retrancher le self eventuellement donc erfc_Luc-1=-erf  luc84p134
                                                         ! iio == jo <=> for the atoms of the same molecule 
               ! erfc - 1 = - erf
               ! Sergiievskyi 31 May 2014: 
               ! > why not to remove the intermolecular interactiions at all?
               ! > after all: they are always constant for the rigid molecules, thus do not matter

!               write(*,*) 'er_after=',er



               potew=xxclb*chargi*chargj*er/r

  !             write(*,*) 'xclb',xclb,'qi',charge_h,'xxclb*chargi*chargj',xxclb*chargi*chargj,'er',er,'r',r
               uu_ew=uu_ew+potew
               pot=potew


               ! ffr = -1/r dU/dr
               ! for coulomb:  U_ij = q_i q_j erfc(alpha r)/ r
               ! dU/dr = q_i q_j ( - erfc(alpha r) / r^2 + derfc(alpha r) / dr * 1/r )
               ! erfc(alpha r ) = 1 - erf = 1 - 2 / sqrt(pi) int_0^{alpha r} e^-t^2 dt                
               ! -> derfc(alpha r) /dr =  alpha derfc(alpha r)/d(alpha r) = - 2 alpha/sqrt(pi) e^{-alpha^2 r^2} 
               !  
               ! ffr = q_i q_j ( erfc(alpha r) /r^3 +  {2 alpha /r^2 }exp(-alpha^2 r^2) ) 
               ffr=xxclb*chargi*chargj*(er/r+alp2pi*ee2)/r2     


!               write(*,*) 'xclbqiqj:',xxclb*chargi*chargj,'er:',er,'r:',r,'alp2pi:',alp2pi,'ee2:',ee2,'r2:',r2    
!               write(*,*) 'ffr_c:',ffr

               ! O-O ineractions : LJ ewald
               ! c.f. Karasawa, Goddard, JPC 1989, 93, 7320-7327
               !      Williams,   Acta Cryst A,  1971, 27, 452
               !    
               ! General case: potential 1/r^w
               ! 
               !     S =  S_real + S_fourier + S_{k=0} - S_self
               !
               ! where
               ! 
               !  S_real =  K  SUM_{sjn,r_sjn>0}  A_sj Gamma(w/2,a^2) r^-w 
               !  S_fourier = K  pi^{w-1.5}/V  SUM_{m<>0}  F_2(k_m) h_m^3 Gamma(3/2 - w/2,b^2) 
               !  S_{k=0}  =  K pi^1.5 alpha^{w - 3} / V * 2/(w-3) SUM_{sj} A_sj
               !  S_self = K 2 alpha^{w} / w * SUM_j A_jj 
               !
               !  K = 1 / (2 Gamma(w/2) )
               !  a^2  =  alpha^2 r_{sj;n}^2 
               !  b^2  =  pi^2 h_m^2 / alpha
               !  h_m  = k_m / 2pi = m / L
               !  r_{sj;n} = | r_s - r_j + n*L |
               !  Gamma(w/2; a^2) = int_{a^2}^inf t^{w/2-1} exp(-t) dt
               !  Gamma(w/2) = Gamma(w/2;0)
               !  F_2(k_m) = SUM_sj A_sj exp(i k_m (r_s - r_j ) )  
               !     for A_sj = q_s*q_j  F_2(k_m) = F(k_m) F(-k_m) = |F(k_m)|^2 
               !     where F(k_m) = SUM_j q_j exp(i k_m r_j)
               !
               ! Particular cases: w=1 --> usual Ewald sum
               !
               ! w=2k+2 (k integer): for that case
               !
               !  Gamma(w/2; a^2) = k! exp(-a^2)  SUM_{p=0}^k a^{2p}/p!
               !  Gamma(3/2 - w/2; b^2)  = (-2)^k / (2k-1)!! [  sqrt(pi) erfc(b)  - b^-1 exp(-b^2) ]
               !                            + exp(-b^2) / (2k-1)!! SUM_{p=1}^{k-1}  (-2)^(k-p) (2p-1)!! b^{-2p-1}
               !
               ! which gives:
               ! w=6:   
               !     S_real = 1/2 SUM_{sjn} A_sj ( 1 + a^2 + 0.5 a^4) exp(-a^2) r^-6
               !     S_fourier = pi^4.5/3V SUM_m F_2(k_m) h_m^3 [ sqrt(pi) erfc(b) + (1/2b^3 - 1/b) exp(-b^2 ] 
               !     S_{k=0} = 1/6V * pi^1.5 alpha^3 SUM_sj A_sj 
               !     S_self = alpha^6/12 SUM_j A_jj
               !
               ! w=12: 
               !     S_real = 1/2 SUM_{r_nsj>0} A_sj (a^10/120 + a^8/24 + a^6/6 + a^4/2 + a^2 + 1 ) exp(-a^2) r^-12
               !     S_fourier =  1/(945*120) pi^10.5/V *
               !                *  SUM_m<>0 F_2(k_m) h_m^9 (-16 sqrt(pi) erfc(b) + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) )
               !     S_{k=0} = pi^1.5 alpha^9 / 1080V * SUM_sj A_sj
               !     S_self =  alpha^12 / 1440 SUM_j A_jj 
  
 !              write(*,*) 'sigma2=',sigma2
!               write(*,*) 'i',i,'iio',iio,'j',j,'jo',jo 

             if(i==iio.and.j==jo) then            ! OO: LJ total  luc84p147

                  vv=(sigma2/r2)**3                ! vv = (sigma/r)^6
 
!                  write(*,*) 'sigma2=',sigma2,'r2=',r2
                   c6=1.d0+alpr2*(1.d0+alpr2/2.d0)  !   1 + a^2 + 0.5 a^4
                  potlj6=-xlj*vv*ee2*c6            !  4*epsilon * (1+a^2 + 0.5a^4) exp(-a^2) (sigma/r)^6
!                  write(*,*) 'xlj:',xlj,'vv:',vv,'ee2',ee2,'c6:',c6
                  alpr6=alpr2**3                   !  a^6 
                  c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0  ! C12 = 1 + a^2 + 0.5 a^4 + a^6/6 + a^8/24 + a^10/120
                  potlj12=xlj*vv**2*ee2*c12                              ! 4 epsilon * C12  * exp(-a^2) * (sigma/r)^12
                  uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
                  pot=pot+potlj12+potlj6
                  ffr=ffr+xlj*6.d0*vv*ee2*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2

  !                write(*,*) 'xlj',xlj,'vv',vv,'ee2',ee2,'c12',c12,'alpr6',alpr6,'c6',c6
 !                 write(*,*) 'ffr_tot:',ffr
!
 !                 write(*,*) 'vv=',vv,'c6=',c6,'potlj6=',potlj6,'potlj12=',potlj12,'c12=',c12
!                  write(*,*) 'potlj6=',potlj6,'potlj12=',potlj12,'potlj=',potlj12 + potlj6

                      !
                      ! 1/r^2(   4 eps * 12 sigma^12/r^12 e^-a^2 ( C12 + (alpha r)^12 / 6! )
                      !        - 4 eps * 6 * sigma^6/r^6  e^-a^2 ( C6 + (alpha r)^6 / 3! )  )
                      !
                      
                      
      !     pour test, LJ tronque shifte (xlj a ete mis a 0)         NON
      !            potlj6=-xlj0*vv
      !            potlj12=xlj0*vv**2
      !            uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
      !            pot=pot+potlj12+potlj6
      !            ffr=ffr+xlj0*6.d0*vv*(2.d0*vv-1.d0)/r2
                end if
 
!         write(*,*) 'LUC: pot:',pot,'ew:',potew,'lj6:',potlj6,'lj12:',potlj12

                ! not same molecules
                if(iio/=jo) then 
                   uu(i)=uu(i)+pot                  ! ne contient plus le self
                   uu(j)=uu(j)+pot                  ! mais contient l'intra!
                   fxij=xij*ffr                     ! NON, plus maintenant!

                   !write(*,*) 'ffr:',ffr,'fxij:',fxij,'x(i):',xx(i),'x(j):',xx(j),'xij',xij
                   fx(i)=fx(i)+fxij
                   fx(j)=fx(j)-fxij                  
                   fyij=yij*ffr
                   fy(i)=fy(i)+fyij
                   fy(j)=fy(j)-fyij
                   fzij=zij*ffr
                   fz(i)=fz(i)+fzij
                   fz(j)=fz(j)-fzij
                end if ! iio /= jo

             end if ! r2 < rmax2
         end do ! j
      end do ! i


      uu_ew_r=uu_ew
      uu_lj_r=uu_lj12+uu_lj6
      uu_lj6_r = uu_lj6
      uu_lj12_r = uu_lj12
      
!       PRINT*, 'LUC_REAL_SUM: uu_lj6:',uu_lj6,'uu_lj12:',uu_lj12,'uu_lj:',uu_lj_r,'ew:',uu_ew_r
  

!      write(*,*) 'uu_Luc:'
!      call write_real_array(0,uu,n)

!      write(*,*) 'forces Luc:'
!      call write_xyz_array(0,fx,fy,fz,n)

!      return
      !
      !        je prepare la somme en k pour Ewald  luc70p21
      !
      !
      !  U_fourier = 1/2V SUM_m<>0 4pi/k_m SUM_js q_jq_s exp(ik_m r_js) exp(-k_m^2/4alpha^2) 
      !
      !  Using that rho(k) = SUM q_j exp(i k r_j)
      !
      !  U_fourier = 1/2V SUM_k<>0 4pi/k^2 |rho(k)|^2 exp(-k^2/(4 alpha^2 )) 
      !
      pialp=pi/alp
      !pialp2 = pialp**2
      pialp2=pi**2/alp**2

      kk=0
      pi2=two_pi
        
 
      xkmax2 = kmax**2
 
!      write(*,*) 'KMAX = ',kmax
 
      do kx=0,kmax
     
 
        xkx2=kx**2
        kymax=SQRT(xkmax2-xkx2)
 
 !       write(*,*) 'KX:',kx,'Kymax=',kymax
      
        do ky=0,kymax
 
            xky2=ky**2
            kzmax=SQRT(xkmax2-xkx2-xky2)
      
  !          write(*,*) 'KY:',ky,'kzmax:',kzmax
      
            do kz=0,kzmax
                 xk2=xkx2+xky2+kz**2

               

                 if(xk2<0.1) cycle            ! k different de 0
                 kk=kk+1

   !              write(*,*) 'KZ:',kz,'xk2',xk2,'kk',kk
         
                 if(kk.gt.nkmax) then
                      WRITE(*,*) 'attention: on depasse la dimension en k',nkmax
                      return
!                      goto 9999
                 end if
     
                 ! pialp2 = pi^2 / alpha^2
                 pialpk2=pialp2*xk2    ! pi^2 / alpha^2 * k^2 
!                k_m = 2*pi*m/L 
!                exp( - k_m^2 / 4 alpha^2 ) = exp( - 4 pi^2 m^2 / L^2 / 4 alpha^2) = 
!                  = exp( - pi^2 * m^2 / L^2 /alpha^2 ) 
!                now, in our notation m <-> k, L=1
!                  exp ( -k_m^2 / 4 alpha^2 ) = exp( - pi^2 *k^2 / alpha^2 ) 
!                                     
                 eek=EXP(-pialpk2)     ! exp(-pi^2/alpha^2 * k^2) 
                 bet(kk)=eek/xk2/pi2   ! 1/2piL pour charge-charge
                                       ! bet(kk) = 1/2pi/k^2 * exp(-pi^2 / alpha^2 * k^2 ) 

                                     ! why 1/2pi: 1/2 is the sum prefactor
                                     ! 4pi/k_m^2 = 4pi/ (4pi^2 L^2/m^2) = 1/ (pi * m^2) = 1 / (pi* k^2) (in our notation)
                 xk=SQRT(xk2)          
                 pialpk=pialp*xk       ! pi/alpha * k


           ! calculation of beta6 and beta12 for LJ sums in Fourier space
           ! 
           !      erfk=erfc_Luc(pialpk,pialpk2,eek)
           !      bet6(kk)=pi45/3.d0*xk2*xk*(pi12*erfk+(0.5d0/pialpk2-1.d0)/pialpk*eek)
           !      bet12(kk)=pi105/1080.d0*xk2**4*xk*(-16.d0*pi12/105.d0*erfk+           &
           !         (((((1.d0/pialpk2-2.d0/7.d0)/pialpk2+4.d0/35.d0)/pialpk2-8.d0/105.d0)/pialpk2+16.d0/105.d0))/pialpk*eek)
           !      nouveau: on peut utiliser un DL ou/et une integration numerique partielle  luc84p150
           !
           ! ( for details:  see definition of erfc_Luc_bet6_bet12)

                 call erfc_Luc_bet6_bet12(pialpk,pialpk2,eek,erfk,x6,x12)
                 bet6(kk)=alp**3*x6   ! alpha^3 *   pi^1.5/3V * b^3 (sqrt(pi) erfc(b) + (1/2b^3 -1b)exp(-b^2) ) 
                 bet12(kk)=alp**9*x12 ! alpha^9 *  1/(945*120) *pi^1.5 *  b^9 * 
                                      !            ( -16 sqrt(pi) erfc(b) 
                                      !              + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) 
                                      !             ) 

      
                 if(kx>0) then                     ! pour compter -kx
                     bet(kk)=2.d0*bet(kk)
                     bet6(kk)=2.d0*bet6(kk)
                     bet12(kk)=2.d0*bet12(kk)
                 end if
      
                 if(ky.eq.0) then                           ! iq sert a savoir si l'un des ki est nul
                     if(kz.eq.0) iq(kk)=1
                     if(kz.gt.0) iq(kk)=3
                 else
                     if(kz.eq.0) iq(kk)=2
                     if(kz.gt.0) iq(kk)=4
                 end if
                 ! iq = 1  -->  ky=0, kz=0     
                 ! iq = 2  -->  ky>0, kz=0
                 ! iq = 3  -->  ky=0, kz>0
                 ! iq = 4  -->  ky>0, kz>0

 
                 ip(kk)=0
            end do  ! kz
      
            ip(kk)=1              ! repere le changement de ky au cas suivant
      
        end do ! ky
      
        ip(kk)=2     ! repere le changement de kx au cas suivant
      
      end do ! kx
      
      ip(kk)=3
      
      kkmax=kk
      
!      WRITE(*,*) 'nombre de k : ',kkmax
      

!    e.g. bet(kk) =  1/2pi/k^2 * exp(-pi^2 / alpha^2 * k^2 )


!      return 
      
      if(kkmax>nkmax) PRINT*, 'ATTENTION: dépasse la taille des tableaux!!!!'
      !
!        puis l'energie Ewald k  astuce luc70p35
!        il faut somme de Zi*exp(2ipi/L*ri*k)  luc70p35
!        et somme sur les i=O uniquement de exp()  luc84p148
      do kk=1,kkmax
         sumcos(kk)=0.
         sumsin(kk)=0.
         sumcos1(kk)=0.
         sumsin1(kk)=0.
         sumcos2(kk)=0.
         sumsin2(kk)=0.
         sumcos3(kk)=0.
         sumsin3(kk)=0.
         sumcos_o(kk)=0.
         sumsin_o(kk)=0.
         sumcos1_o(kk)=0.
         sumsin1_o(kk)=0.
         sumcos2_o(kk)=0.
         sumsin2_o(kk)=0.
         sumcos3_o(kk)=0.
         sumsin3_o(kk)=0.
      end do
!
      pi2 = two_pi

      do i=1,n                                   ! debut i
! do i=3,3
!
         iio=3*((i-1)/3)+1
         chi=+1.d0

         if(i==iio) chi=-2.d0

         cx=COS(pi2*xx(i))  ! cos ( 2 pi *1 /L * x)  where L=1
         sx=SIN(pi2*xx(i))  
         cy=COS(pi2*yy(i))
         sy=SIN(pi2*yy(i))
         cz=COS(pi2*zz(i))
         sz=SIN(pi2*zz(i))

         !write(*,*) 'x',xx(i),'cx:',cx,'y',yy(i),'cy',cy,'z',zz(i),'cz:',cz

         cckx(0)=1.d0
         sskx(0)=0.d0
         ccky(0)=1.d0
         ssky(0)=0.d0
         cckz(0)=1.d0
         sskz(0)=0.d0

          ! Table of the functions
          ! cos(2pi m /L x) = cos( x*k_m )
          ! sin(k_m x), cos(k_m y) ... 
         do k=1,kmax
            cckx(k)=cckx(k-1)*cx-sskx(k-1)*sx   
            sskx(k)=sskx(k-1)*cx+cckx(k-1)*sx
            ccky(k)=ccky(k-1)*cy-ssky(k-1)*sy                  ! tableaux suivant les 3D
            ssky(k)=ssky(k-1)*cy+ccky(k-1)*sy
            cckz(k)=cckz(k-1)*cz-sskz(k-1)*sz
            sskz(k)=sskz(k-1)*cz+cckz(k-1)*sz
         end do
 
!         write(tmpstr,'(AI3.3A)') 'cck',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!         call write_xyz_array(hfile,cckx,ccky,cckz,kmax)
!         call io_close(hfile)
!
!         write(tmpstr,'(AI3.3A)') 'ssk',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!         call write_xyz_array(hfile,sskx,ssky,sskz,kmax)
!         call io_close(hfile)

     

!        write(*,*) 'cckx:'
!        call write_xyz_array(0,cckx,ccky,cckz,kmax)
!        write(*,*) 'sskx:'
!        call write_xyz_array(0,sskx,ssky,sskz,kmax)

         kx=0
         ky=0
         kz=0
         ! chi - charge: for H chi=1, for O chi=-2
         cckxky=chi                               ! kx,ky   chi * cos (x*kx+ y*ky) = chi cos 0  
         sskxky=0.                                !         chi * sin(x*kx+y*ky) = chi shi 0 
         cckxkya=chi                              ! kx,-ky  chi * cos( 0 - 0)
         sskxkya=0.                               !         chi * sin( 0 - 0)
  !
         do kk=1,kkmax                           ! debut k

            kz=kz+1

            q1=cckxky*cckz(kz)
            q2=sskxky*sskz(kz)
            q3=sskxky*cckz(kz)
            q4=cckxky*sskz(kz)

!            write(*,*) 'LUCkx:',kx,'ky:',ky,'kz:',kz
        !    write(*,*) 'LUC:',q1,q2,q3,q4

            ! x*kx + y*ky + z*kz
            ! cos(x*kx + y*ky)cos(z*kz) - sin(x*kx + y*ky) sin(z*kz) = 
            ! = cos(x*kx + y*ky + z*kz )        
             ! sin(x*kx+y*ky)cos(z*kz) + cos(x*kx+y*ky)sin(z*kz) = sin(x*kx+y*ky+z*kz)                                               
            sumcos(kk)=sumcos(kk)+q1-q2       ! ky,kz
            sumsin(kk)=sumsin(kk)+q3+q4       

            ! kx+ky - kz
            sumcos1(kk)=sumcos1(kk)+q1+q2     ! ky,-kz  cos
            sumsin1(kk)=sumsin1(kk)+q3-q4     !         sin

            if(i==iio) then                 ! pour O, ajouter LJ
               sumcos_o(kk)=sumcos_o(kk)+q1-q2       ! ky,kz                     contient zi=-2!
               sumsin_o(kk)=sumsin_o(kk)+q3+q4
               sumcos1_o(kk)=sumcos1_o(kk)+q1+q2     ! ky,-kz
               sumsin1_o(kk)=sumsin1_o(kk)+q3-q4
            end if

            ! cckxkya = chi cos(kx-ky) 
            q1=cckxkya*cckz(kz)
            q2=sskxkya*sskz(kz)
            q3=sskxkya*cckz(kz)
            q4=cckxkya*sskz(kz)

            !  kx - ky + kz
            sumcos2(kk)=sumcos2(kk)+q1-q2         ! -ky,kz
            sumsin2(kk)=sumsin2(kk)+q3+q4

            ! kx - ky - kz
            sumcos3(kk)=sumcos3(kk)+q1+q2         ! -ky,-kz
            sumsin3(kk)=sumsin3(kk)+q3-q4

            if(i==iio) then                 ! pour O, ajouter LJ
               sumcos2_o(kk)=sumcos2_o(kk)+q1-q2       ! ky,kz
               sumsin2_o(kk)=sumsin2_o(kk)+q3+q4
               sumcos3_o(kk)=sumcos3_o(kk)+q1+q2     ! ky,-kz
               sumsin3_o(kk)=sumsin3_o(kk)+q3-q4
            end if

             
            ! ip - end of row indicator (filled in the cycle before)
            ! ip(kk) == 1  - end of row in kz
            ! ip(kk) == 2  - end of row in ky
            if(ip(kk).eq.1.or.ip(kk).eq.2) then       ! fin de boucle en kz et peut-etre en ky

               if(ip(kk).eq.1) then
                  ky=ky+1
               else
                  ky=0
                  kx=kx+1
               end if
      
               q1=cckx(kx)*ccky(ky)*chi
               q2=sskx(kx)*ssky(ky)*chi
               q3=sskx(kx)*ccky(ky)*chi
               q4=cckx(kx)*ssky(ky)*chi
      
               cckxky=q1-q2 ! chi*( cos(x*kx)cos(y*ky) - sin(x*kx)sin(y*ky) ) = chi * cos(x*kx+y*ky)
               sskxky=q3+q4 ! chi*( sin(x*kx)cos(y*ky) + cos(x*kx)sin(y*ky) ) = chi * sin(x*kx+y*ky)
               cckxkya=q1+q2! chi*( cos(x*kx)cos(y*ky) + sin(x*kx)sin(y*ky) ) = chi * cos(x*kx-y*ky)
               sskxkya=q3-q4! chi*( sin(x*kx)cos(y*ky) - cos(x*kx)sin(y*ky) ) = chi * sin(x*kx-y*ky)

!               write(*,*) 'LUC ip',kk,kx,ky,cckx(kx),ccky(ky)

               kz=-1
            end if
          end do                                  ! fin kk
  !
       end do                                   ! fin i
!                                                                attention, les _o contiennent zi=zO=-2
      sumcos_o(1:kkmax)=sumcos_o(1:kkmax)/(-2.d0)
      sumsin_o(1:kkmax)=sumsin_o(1:kkmax)/(-2.d0)
      sumcos1_o(1:kkmax)=sumcos1_o(1:kkmax)/(-2.d0)
      sumsin1_o(1:kkmax)=sumsin1_o(1:kkmax)/(-2.d0)
      sumcos2_o(1:kkmax)=sumcos2_o(1:kkmax)/(-2.d0)
      sumsin2_o(1:kkmax)=sumsin2_o(1:kkmax)/(-2.d0)
      sumcos3_o(1:kkmax)=sumcos3_o(1:kkmax)/(-2.d0)
      sumsin3_o(1:kkmax)=sumsin3_o(1:kkmax)/(-2.d0)
!


!      return

      utotk=0.
      utotk6=0.; utotk12=0.
      do kk=1,kkmax
         som=sumcos(kk)**2+sumsin(kk)**2

                 ! iq = 1  -->  ky=0, kz=0     
                 ! iq = 2  -->  ky>0, kz=0
                 ! iq = 3  -->  ky=0, kz>0
                 ! iq = 4  -->  ky>0, kz>0
         if(iq(kk).eq.3.or.iq(kk).eq.4) som=som+sumcos1(kk)**2+sumsin1(kk)**2        ! kz>0                           ! kz>0
         if(iq(kk).eq.2.or.iq(kk).eq.4) som=som+sumcos2(kk)**2+sumsin2(kk)**2        ! ky>0                           ! ky>0
         if(iq(kk).eq.4) som=som+sumcos3(kk)**2+sumsin3(kk)**2                       ! ky,kz>0                           ! ky,kz>0

         utotk=utotk+bet(kk)*som
         som_o=sumcos_o(kk)**2+sumsin_o(kk)**2

 
         if(iq(kk).eq.3.or.iq(kk).eq.4) som_o=som_o+sumcos1_o(kk)**2+sumsin1_o(kk)**2        ! kz>0                           ! kz>0
         if(iq(kk).eq.2.or.iq(kk).eq.4) som_o=som_o+sumcos2_o(kk)**2+sumsin2_o(kk)**2        ! ky>0                           ! ky>0
         if(iq(kk).eq.4) som_o=som_o+sumcos3_o(kk)**2+sumsin3_o(kk)**2                       ! ky,kz>0                           ! ky,kz>0
 
         utotk6=utotk6+bet6(kk)*som_o
         utotk12=utotk12+bet12(kk)*som_o
 
      end do

!      cdiel_ext = external_permutivity
!      write(*,*) 'cdiel_ext=',cdiel_ext
      uu_ew=uu_ew+xxclb*utotk
      ureste=pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)-alp2pi*n*xxclb
!      write(*,*) 'ureste=',ureste,'pi2',pi2,'xxclb',xxclb,'cdiel_ext',cdiel_ext,'alp2pi',alp2pi

      
      uu_ew=uu_ew+ureste
      uu_lj6=uu_lj6-xlj*sigma2**3*(utotk6+pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
      uu_lj12=uu_lj12+xlj*sigma2**6*(utotk12+pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)

!      write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew
!
!      write(*,*) 'xclb=',xclb,'xxclb=',xxclb
!      write(*,*) 'utotk=',xxclb*utotk,'utotk6=',utotk6,'utotk12=',utotk12
!      write(*,*) 'xlj:',xlj,'sigma6:',sigma2**3,'xlj*sigma6*utotk6',xlj*sigma2**3*utotk6
!      write(*,*) 'xlj:',xlj,'sigma12:',sigma2**6,'xlj*sigma12*utotk12',xlj*sigma2**6*utotk12


uu_lj6_k = xlj*sigma2**3 * utotk6
uu_lj12_k = xlj*sigma2**6 * utotk12
uu_lj_k = uu_lj12_k - uu_lj6_k

uu_ext = pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)

!write(*,*) 'alp2pi=',alp2pi,'alpha',alpha,'alpha/sqrt(pi)',alpha/sqrt(pi)
uu_tail_ew =  -alp2pi*n*xxclb  ! alp2pi = 2 alp / sqrt(pi)
! alp2pi * n *xxclb = xclb * q^2 6 alpha n_h2o / sqrt(pi) 
! SUM q_s = 2 n_h2o + n_h2o * 4 = 6 n_h2o
uu_tail_lj6 = xlj*sigma2**3*(pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
uu_tail_lj12 = xlj*sigma2**6*(pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)
uu_tail_lj = uu_tail_lj12 - uu_tail_lj6

call write_energy_LUC(n_h2o)

!write(*,*) '******************** LUC *******************************'
!write(*,*) 'U_REAL ew:',uu_ew_r,'lj:',uu_lj_r,'lj6:',uu_lj6_r,'lj12:',uu_lj12_r
!write(*,*) 'U_FOUR ew:',xxclb*utotk,'lj:',uu_lj_k,'lj6:',uu_lj6_k,'lj12:',uu_lj12_k
!write(*,*) 'U_TAIL ew:',uu_tail_ew,'lj:',uu_tail_lj,'lj6:',uu_tail_lj6,'lj12:',uu_tail_lj12
!write(*,*) 'U_EXT  ew:',uu_ext
!write(*,*) 'U_TOTL ew:',uu_ew,'lj:',uu_lj12 + uu_lj6,'lj6:',-uu_lj6,'lj12:',uu_lj12
!write(*,*) 'U_TOT tot:',uu_ew + uu_lj6 + uu_lj12



if( present(systematic_correction) .and. systematic_correction ) then 
   
   
   write(*,*) '********** LUC AFTER CORRECTION ******************'
   
   
   ssk_vrai=kmax*pi/alp
   !      return   
    
         !      je peux essayer des corrections systematiques  luc84p151
    
         qmax=2.d0*pi*kmax
   !      write(*,*) 'qmax=',qmax
   
!         roh = roh_a / xl_a
!         rhh = rhh_a / xl_a    
   
        sigma = SQRT(sigma2)
   
         uu_ew=uu_ew+n_h2o/pi*xxclb*alp*EXP(-ssk_vrai**2)/ssk_vrai*(6.d0-8.d0*SIN(qmax*roh)/(qmax*roh)&
               +2.d0*SIN(qmax*rhh)/(qmax*rhh))
         write(*,*) 'uu_ew=',uu_ew,'n_h2o=',n_h2o,'xxclb=',xxclb,'ssk_vrai=',ssk_vrai,'roh=',roh,'rhh=',rhh
         write(*,*) 'uu_ew=',uu_ew 
   
        uu_lj6=uu_lj6-epsi_lj*sigma**6*         &
            (2.d0*pi*n_h2o**2*alp**3*EXP(-ssr**2)/ssr+2.d0/SQRT(pi)*n_h2o*alp**6*EXP(-ssk_vrai**2)/ssk_vrai)
        write(*,*) 'uu_lj6=',uu_lj6,'epsi_lj=',epsi_lj,'sigma=',sigma
   
         !  
         ! SUM_sj A_sj * 2 pi alpha^3  
         uu_lj12=uu_lj12+epsi_lj*sigma**12*      &
            (pi/30.d0*n_h2o**2*alp**9*EXP(-ssr**2)/ssr+1.d0/30.d0/SQRT(pi)*n_h2o*alp**12*EXP(-ssk_vrai**2)/ssk_vrai)
         write(*,*) 'uu_lj12=',uu_lj12  

end if 

!      PRINT*, 'et apres corrections des queues en r et en k:'
!            write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew

      !      uself deja retranchee dans la boucle Ewald r
uu_lj=uu_lj12+uu_lj6
uutot=uu_lj+uu_ew

write(*,*) 'Energie totale actuelle : ',uutot


!this % bet => bet
!this % bet6 => bet6
!this % bet12 => this % bet12

!uu => this % uu
!fx => this % fx
!fy => this % fy
!fz => this % fz
!xx => this % xx
!yy => this % yy 
!zz => this % zz

!xl_a = this % BoxLength
!sigma2 = this % sigma2 
!sigma_a = sqrt(sigma2) * BoxLength
!xlj = this % xlj


      return 

      !xlj=xlj0       ! sert pour le test  NON, plus maintenant
      !
      !             pour les forces en k, on doit refaire une boucle en i et k
      !             luc84p133
      !             NON, plus maintenant, puisque fx,y,z ne contient que la partie en r luc84p187
      !
      !
      !          je definis les mnlmunu
      !

      ialp=0
      do mnmaxi=0,mnmax             ! boucle en max(m,n) a l'exterieur
      do m0=0,mnmax                      ! ainsi 1,1 sera avant 0,2 et 2,0

            mudeb=0
            mumax=m0
            ndeb=m0
    
            do n0=ndeb,mnmaxi
                
                if(m0/=mnmaxi.and.n0/=mnmaxi) cycle             ! m ou n egal a max(m,n)

                nudeb=0
                numax=n0
                
                do l0=ABS(m0-n0),m0+n0
                do mu0=mudeb,mumax

                      if((-1)**mu0==-1) cycle
                      if(m0==n0) nudeb=mu0

                      do nu0=nudeb,numax
  
                          if((-1)**nu0==-1) cycle
                          if(mu0==0.and.nu0==0.AND.(-1)**(m0+n0+l0)==-1) cycle
  
                          ialp=ialp+1
                          mm(ialp)=m0
                          nn(ialp)=n0
                          ll(ialp)=l0
                          mumu(ialp)=mu0
                          nunu(ialp)=nu0
            !              WRITE(*,*) ialp,m0,n0,l0,mu0,nu0
                          if(mu0>0.and.nu0>0.AND.(m0/=n0.OR.(mu0.ne.nu0).OR.(-1)**l0==1)) then
                             ialp=ialp+1                  ! - nu
                             mm(ialp)=m0
                             nn(ialp)=n0
                             ll(ialp)=l0
                             mumu(ialp)=mu0
                             nunu(ialp)=-nu0
               !              WRITE(*,*) ialp,m0,n0,l0,mu0,-nu0
                          end if
                      end do ! nu0

               end do ! mu0
               end do ! l0

            end do ! n0

      end do ! m0
      end do ! mnmaxi

      PRINT*, 'nombre de projections mnlmunu: ',ialp

      !
      !          je definis les mnmunukhi
      !
      ialp=0

      do mnmini=0,mnmax              ! boucle en min(m,n) a l'exterieur
      do m0=mnmini,mnmax                             ! ainsi 1,1 sera apres 0,2 et 2,0

         mudeb=0
         mumax=m0
         ndeb=m0

         do n0=ndeb,mnmax
           
              if(m0/=mnmini.and.n0/=mnmini) cycle             ! m ou n egal a min(m,n)
           
              nudeb=0
              numax=n0

              do khi0=-MIN(m0,n0),MIN(m0,n0)             ! khi et non l !
              do mu0=mudeb,mumax

                  if((-1)**mu0==-1) cycle

                  if(m0==n0) nudeb=mu0

                  do nu0=nudeb,numax

                     if((-1)**nu0==-1) cycle
                     if(mu0==0.and.nu0==0.AND.khi0<0) cycle

                     ialp=ialp+1
                     mm1(ialp)=m0
                     nn1(ialp)=n0
                     khi1(ialp)=khi0
                     mumu1(ialp)=mu0
                     nunu1(ialp)=nu0
      !              if(itest==0) WRITE(*,*) ialp,m0,n0,mu0,nu0,khi0
                     if(mu0>0.and.nu0>0.AND.(m0/=n0.OR.(mu0.ne.nu0).OR.khi0>=0)) then
                        ialp=ialp+1              ! -nu
                        mm1(ialp)=m0
                        nn1(ialp)=n0
                        khi1(ialp)=khi0
                        mumu1(ialp)=mu0
                        nunu1(ialp)=-nu0
          !              if(itest==0) WRITE(*,*) ialp,m0,n0,mu0,-nu0,khi0
                     end if

                  end do ! nu0

              end do ! mu0
              end do ! khi0

         end do ! n0

      end do ! m0
      end do ! mnmini

      PRINT*, 'nombre de projections mnmunukhi: ',ialp

!
!        initialisation du generateur de nombres aleatoires
!
!      x=alea(0)
!        on y va !
      return
!
end
!
subroutine rot_vect(x,y,z,x1,y1,z1,cc,ss,i) !DOC
IMPLICIT REAL(8) (a-h,o-z)
!
!DOC        fait la rotation de x,y,z en x1,y1,z1 d'angle cc=cos,ss=sin autour d'un (i) des 3 axes principaux
!DOC        attention: x1,y1,z1 peut etre a la meme place memoire que x,y,z
!DOC        luc80p176
!
if(i==1) then
x1=x
y2=y
y1=cc*y2-ss*z
z1=ss*y2+cc*z
end if
if(i==2) then
y1=y
x2=x
x1=cc*x2+ss*z
z1=-ss*x2+cc*z
end if
if(i==3) then
z1=z
x2=x
x1=cc*x2-ss*y
y1=ss*x2+cc*y
end if
end subroutine
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


!     Last change:  LB    4 Feb 2014    3:59 pm
subroutine movemc_Luc(xx,yy,zz,n_h2o, i_h2o,iaccep) !DOC
!DOC Luc's movemc function
!DOC Parameters:
!use parameters
use constants
use io, only : write_real_array, write_xyz_array,io_open,io_close
implicit REAL(8) (a-h,o-z)
! DP maintenant      real * 4 alea
common/nombre/n

real(8),dimension(:) :: xx,yy,zz !DOC coordinates of atoms
integer,intent(in) :: n_h2o !DOC number of h2o molecules
integer,intent(in) :: i_h2o !DOC number of molecule to move
integer,intent(out) :: iaccep !DOC whether the movement is accepted or not

integer :: hfile
character(256) :: tmpstr

common/dipole/summux,summuy,summuz
common/xenergy/utotk,utotk6,utotk12,uu_lj6_r,uu_lj12_r

integer :: n
!common/positi/xx(10000),yy(10000),zz(10000)        ! positions cellule centrale
common/energi/uu(10000),vrnew(10000),uutot
COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r
COMMON/forces/fx(10000),fy(10000),fz(10000)
common/conc/conc_h2o,xl_a
common/deplac/dr_a,d_angle
common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d
!common/diamet/sigma_a,rcutoff_a,rcut2,roh_a,theta_d

COMMON/lennardjones/epsi_lj,xlj,ucut                         ! xlj0 pour test
!COMMON/lennardjones/epsi_lj,ucut
COMMON/ewald/temp,dbjr,charg_h,xxclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
!COMMON/ewald/charg_h,xxclb,alp_a1,alp,kkmax,cdiel_ext


COMMON/forcebias/xlambda_f,xlambda_c
PARAMETER(nkmax=4661,kmaxmax=20)
COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
 sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
 sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
 cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
 cckz(0:kmaxmax),sskz(0:kmaxmax)
COMMON/sommeenk_lj/bet6(nkmax),bet12(nkmax),                     &         ! pour Ewald LJ en k
 sumcos_o(nkmax),sumsin_o(nkmax),sumcos1_o(nkmax),sumsin1_o(nkmax),        &
 sumcos2_o(nkmax),sumsin2_o(nkmax),sumcos3_o(nkmax),sumsin3_o(nkmax)
COMMON/sommeenkold/sumcosold(nkmax),sumsinold(nkmax),sumcos1old(nkmax),sumsin1old(nkmax),        &
 sumcos2old(nkmax),sumsin2old(nkmax),sumcos3old(nkmax),sumsin3old(nkmax)
COMMON/sommeenkold_lj/sumcosold_o(nkmax),sumsinold_o(nkmax),sumcos1old_o(nkmax),sumsin1old_o(nkmax),        &
 sumcos2old_o(nkmax),sumsin2old_o(nkmax),sumcos3old_o(nkmax),sumsin3old_o(nkmax)
COMMON/chrono/i_partie
real(8),DIMENSION(0:kmaxmax)::cckx_old_1,sskx_old_1,ccky_old_1,ssky_old_1,cckz_old_1,sskz_old_1,          &
                             cckx_old_2,sskx_old_2,ccky_old_2,ssky_old_2,cckz_old_2,sskz_old_2,          &
                             cckx_old_3,sskx_old_3,ccky_old_3,ssky_old_3,cckz_old_3,sskz_old_3
real(8),DIMENSION(0:kmaxmax)::cckx_new_1,sskx_new_1,ccky_new_1,ssky_new_1,cckz_new_1,sskz_new_1,          &
                             cckx_new_2,sskx_new_2,ccky_new_2,ssky_new_2,cckz_new_2,sskz_new_2,          &
                             cckx_new_3,sskx_new_3,ccky_new_3,ssky_new_3,cckz_new_3,sskz_new_3
DIMENSION xnew(3),ynew(3),znew(3),uinew(3)
DIMENSION fxold(3),fyold(3),fzold(3)
DIMENSION fxnew(3),fynew(3),fznew(3)
DIMENSION fxnew_r(3),fynew_r(3),fznew_r(3)
DIMENSION frxnew(10000),frynew(10000),frznew(10000)

!COMMON/ewald/temp,dbjr,charg_h,xclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
! temp,dbjr,kmax,rmax2 - in parameters
!charg_h = 0.42380d0
!xxclb = xclb * charg_h**2
!alp_a1 = alpha_angstr_m1  ! in parameters
!alp = alpha ! parameters
!kkmax = this % grid % nk  common
!cdiel_ext = external_permutivity ! parameters

!epsi_lj = xlj / 4

! should  be already initialized in inmc

!  They are in commons
!bet => this % bet
!bet6 => this % bet6
!bet12 => this % bet12

!uu => this % uu
!fx => this % fx
!fy => this % fy
!fz => this % fz
!xx => this % xx
!yy => this % yy 
!zz => this % zz

!write(*,*) 'i_h2o = ',i_h2o


!xl_a = BoxLength
!sigma2 = this % sigma2 
!sigma_a = sqrt(sigma2) * BoxLength
!xlj = this % xlj

!
!        sous-programme Monte-Carlo
!        on essaie le deplacement de la particule i
!        on l'accepte si... (iaccep=1)
!
!        Luc Belloni 22-11-90 puis 3-4-98 pour potentiel rectangle
!
!        version H2O luc84p122
!        donc i_h2o repere la molecule H2O
!
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        version Force-Bias  luc84p160
!
!        version "new" qui ruse en k pour reduire le nombre de tableaux et accelerer en k
!        luc84p183
!
!        j'actualise les forces en r pour eviter la boucle old r si deplt refuse  luc84p187
!
!           pour test, je reviens a LJ tronque qui utilise xlj0       NON, plus maintenant
!xlj=0.
!PRINT*, 'Quel i_h2o ? '
!READ(*,*) i_h2o


! ***********************************
!   write(*,*) 'xxclb',xxclb,'xl_a',xl_a,'alp',alp,'sigma2',sigma2,'charg_h',charg_h,'xlj',xlj
!   write(*,*) 'cdiel_ext',cdiel_ext,'dr_a',dr_a,'d_angle',d_angle 
!   write(*,*) 'xlambda_f',xlambda_f,'xlambda_c',xlambda_c
!   write(*,*) 'kmax',kmax,'kkmax',kkmax,'rmax2',rmax2
   


iaccep=0
i0=3*(i_h2o-1)                 ! les 3 atomes de i_h2o seront i0+1,+2,+3
iio=i0+1                        ! no du O
!pi=4.d0*ATAN(1.d0)
pi2=2.d0*pi
alp2pi=2.d0*alp/dSQRT(pi)
!
!
!        il faut commencer par calculer les forces sur les 3 atomes de la molecule i_h2o en position OLD
!
!        pour la partie en r, c'est facile, c'est contenu dans uu et fx,y,z !!!
!
!
fxold(1:3)=fx(i0+1:i0+3); fyold(1:3)=fy(i0+1:i0+3); fzold(1:3)=fz(i0+1:i0+3)        ! et leur force en r

!write(*,*) 'uu:'
!call write_real_array(0,uu(i0+1:i0+3),3)

!write(*,*) 'fxold real:'
!call write_xyz_array(0,fxold,fyold,fzold,3)

!
!
!        puis l'energie Ewald k  astuce luc70p35
!        il faut Zi*exp(2ipi/L*ri*k)  luc79p160
!        pour les 3 i O,H,H qui composent la molecule deplacee i_h2o
!
!        je fais les 3 i old a l'interieur de la meme boucle en kk !!!
!        luc84p181

! i_partie = 0 in mcluc_h2o_fb.f90. 
! And seems that it never changes, so it is always i_partie==0 (it looks like the 'kind of move' variable, and 0 is the most general case) 

!
! U_fourier = SUM_m  rho^2(k_m)  beta(k_m)   
!   * the same for beta6(k_m) and beta12(k_m), but then rho(k) is taken only over oxygens
!
! where rho(k_m) = SUM_j exp( i k_m r_j ) = SUM_j cos( k_m r_j) + i SUM_j sin( k_m r_j )
!       rho(k_m)^2 = [ SUM_j cos(k_m r_j) ]^2 + [ SUM_j sin(k_m r_j) ]^2
!
! rho(k_m)^2 = rho(-k_m)^2
!
! we can rewrite the sum for only k_m>0:
!
! SUM_m = 2 * SUM_{mx,my_mz>0}  ( rho(+++)^2 + rho(++-)^2 + rho(+-+)^2 + rho(+--)^2 )
!
!  where rho(s1,s2,s3) = C(s1,s2,s3)^2 + S(s1,s2,s3)^2
!  
!  C(s1,s2,s3) = SUM_j cos( s1 kx x_j + s2 ky y_j + s3 kz z_j )
!  S(s1,s2,s3) = SUM_j sin( s1 kx x_j + s2 ky y_j + s3 kz z_j )
!  
!  s1,s2,s3 \in {+,-}
!
!  now, let we move one molecule (3 atoms)
!
!  we should subtract rho(+++)^2,rho(++-)^2, rho(+-+)^2, rho(+--)^2 for these atoms at old coordinates 
!  and add the same values for the new coordinated
!
!  rho(+++)^2, rho(++-)^2 ...  should be computed at all k_m 
!  to do it C(+++) S(+++), C(++-) S(++-) ... should be computed
!    
!  Now  cos( x (kx + dkx) + y ky + zkz ) = 
!       cos( x kx + y ky + z kz) cos(x dkx )
!     - sin( x kx + y ky + z kz) sin(x dkx )

i_partie = 0

IF(i_partie==2.or.i_partie==0) then


!   write(*,*) 'xx.yy.zz(i0+1.2.3)'
!   call write_xyz_array(0,xx(i0+1:i0+3),yy(i0+1:i0+3),zz(i0+1:i0+3),3)

   sigma6=sigma2**3; sigma12=sigma6**2
   pi4=2.d0*pi2  ! pi2= 2pi, pi4=4pi
!
   cx_old_1=COS(pi2*xx(i0+1))
   sx_old_1=SIN(pi2*xx(i0+1))
   cy_old_1=COS(pi2*yy(i0+1))
   sy_old_1=SIN(pi2*yy(i0+1))
   cz_old_1=COS(pi2*zz(i0+1))
   sz_old_1=SIN(pi2*zz(i0+1))

   cx_old_2=COS(pi2*xx(i0+2))
   sx_old_2=SIN(pi2*xx(i0+2))
   cy_old_2=COS(pi2*yy(i0+2))
   sy_old_2=SIN(pi2*yy(i0+2))
   cz_old_2=COS(pi2*zz(i0+2))
   sz_old_2=SIN(pi2*zz(i0+2))

   cx_old_3=COS(pi2*xx(i0+3))
   sx_old_3=SIN(pi2*xx(i0+3))
   cy_old_3=COS(pi2*yy(i0+3))
   sy_old_3=SIN(pi2*yy(i0+3))
   cz_old_3=COS(pi2*zz(i0+3))
   sz_old_3=SIN(pi2*zz(i0+3))

   cckx_old_1(0)=1.
   sskx_old_1(0)=0.
   ccky_old_1(0)=1.
   ssky_old_1(0)=0.
   cckz_old_1(0)=1.
   sskz_old_1(0)=0.

   cckx_old_2(0)=1.
   sskx_old_2(0)=0.
   ccky_old_2(0)=1.
   ssky_old_2(0)=0.
   cckz_old_2(0)=1.
   sskz_old_2(0)=0.
   cckx_old_3(0)=1.

   sskx_old_3(0)=0.
   ccky_old_3(0)=1.
   ssky_old_3(0)=0.
   cckz_old_3(0)=1.
   sskz_old_3(0)=0.


   do k=1,kmax

      cckx_old_1(k)=cckx_old_1(k-1)*cx_old_1-sskx_old_1(k-1)*sx_old_1    ! cos(x (kx+1) ) = cos(x kx) cos(x) - sin(kx)sin(x)
      sskx_old_1(k)=sskx_old_1(k-1)*cx_old_1+cckx_old_1(k-1)*sx_old_1
      ccky_old_1(k)=ccky_old_1(k-1)*cy_old_1-ssky_old_1(k-1)*sy_old_1                  ! tableaux suivant les 3D
      ssky_old_1(k)=ssky_old_1(k-1)*cy_old_1+ccky_old_1(k-1)*sy_old_1
      cckz_old_1(k)=cckz_old_1(k-1)*cz_old_1-sskz_old_1(k-1)*sz_old_1
      sskz_old_1(k)=sskz_old_1(k-1)*cz_old_1+cckz_old_1(k-1)*sz_old_1

      cckx_old_2(k)=cckx_old_2(k-1)*cx_old_2-sskx_old_2(k-1)*sx_old_2
      sskx_old_2(k)=sskx_old_2(k-1)*cx_old_2+cckx_old_2(k-1)*sx_old_2
      ccky_old_2(k)=ccky_old_2(k-1)*cy_old_2-ssky_old_2(k-1)*sy_old_2                  ! tableaux suivant les 3D
      ssky_old_2(k)=ssky_old_2(k-1)*cy_old_2+ccky_old_2(k-1)*sy_old_2
      cckz_old_2(k)=cckz_old_2(k-1)*cz_old_2-sskz_old_2(k-1)*sz_old_2
      sskz_old_2(k)=sskz_old_2(k-1)*cz_old_2+cckz_old_2(k-1)*sz_old_2

      cckx_old_3(k)=cckx_old_3(k-1)*cx_old_3-sskx_old_3(k-1)*sx_old_3
      sskx_old_3(k)=sskx_old_3(k-1)*cx_old_3+cckx_old_3(k-1)*sx_old_3
      ccky_old_3(k)=ccky_old_3(k-1)*cy_old_3-ssky_old_3(k-1)*sy_old_3                  ! tableaux suivant les 3D
      ssky_old_3(k)=ssky_old_3(k-1)*cy_old_3+ccky_old_3(k-1)*sy_old_3
      cckz_old_3(k)=cckz_old_3(k-1)*cz_old_3-sskz_old_3(k-1)*sz_old_3
      sskz_old_3(k)=sskz_old_3(k-1)*cz_old_3+cckz_old_3(k-1)*sz_old_3
   end do
 

   cckxky_old_1=1.d0                               ! kx,ky
   sskxky_old_1=0.
   cckxkya_old_1=1.d0                              ! kx,-ky
   sskxkya_old_1=0.

   cckxky_old_2=1.d0                               ! kx,ky
   sskxky_old_2=0.
   cckxkya_old_2=1.d0                              ! kx,-ky
   sskxkya_old_2=0.

   cckxky_old_3=1.d0                               ! kx,ky
   sskxky_old_3=0.
   cckxkya_old_3=1.d0                              ! kx,-ky
   sskxkya_old_3=0.

   kx=0
   ky=0
   kz=0
!
!      debut de la somme en kk
!
   do kk=1,kkmax                           ! debut k

      kz=kz+1

      q1_old_1=cckxky_old_1*cckz_old_1(kz)  ! cos(x kx + y ky) cos( z kz )
      q2_old_1=sskxky_old_1*sskz_old_1(kz)  ! sin(x kx + y ky) sin( z kz )
      q3_old_1=sskxky_old_1*cckz_old_1(kz)  ! sin(x kx + y ky) cos( z kz )
      q4_old_1=cckxky_old_1*sskz_old_1(kz)  ! cos(x kx + y ky) sin( z kz )

      cosi_old_1=+q1_old_1-q2_old_1       ! ky,kz     ! c(x+y) c(z)  - s(x+y)s(z) = c(x+y+z)
      sini_old_1=+q3_old_1+q4_old_1                   ! s(x+y) c(z)  + c(x+y)s(z) = s(x+y+z)
      cosi1_old_1=+q1_old_1+q2_old_1     ! ky,-kz     ! c(x+y-z)
      sini1_old_1=+q3_old_1-q4_old_1                  ! s(x+y-z)

      q1_old_1=cckxkya_old_1*cckz_old_1(kz)         ! cos(x-y) cos(z) 
      q2_old_1=sskxkya_old_1*sskz_old_1(kz)         ! sin(x-y) sin(z)
      q3_old_1=sskxkya_old_1*cckz_old_1(kz)
      q4_old_1=cckxkya_old_1*sskz_old_1(kz)

      cosi2_old_1=+q1_old_1-q2_old_1         ! -ky,kz  ! c(x-y+z)
      sini2_old_1=+q3_old_1+q4_old_1                   ! s(x-y+z)
      cosi3_old_1=+q1_old_1+q2_old_1         ! -ky,-kz ! c(x-y-z)
      sini3_old_1=+q3_old_1-q4_old_1                   ! s(x-y-z)

    
      ! the same for the 2nd atom
      q1_old_2=cckxky_old_2*cckz_old_2(kz)
      q2_old_2=sskxky_old_2*sskz_old_2(kz)
      q3_old_2=sskxky_old_2*cckz_old_2(kz)
      q4_old_2=cckxky_old_2*sskz_old_2(kz)

      cosi_old_2=+q1_old_2-q2_old_2       ! ky,kz
      sini_old_2=+q3_old_2+q4_old_2
      cosi1_old_2=+q1_old_2+q2_old_2     ! ky,-kz
      sini1_old_2=+q3_old_2-q4_old_2

      q1_old_2=cckxkya_old_2*cckz_old_2(kz)
      q2_old_2=sskxkya_old_2*sskz_old_2(kz)
      q3_old_2=sskxkya_old_2*cckz_old_2(kz)
      q4_old_2=cckxkya_old_2*sskz_old_2(kz)

      cosi2_old_2=+q1_old_2-q2_old_2         ! -ky,kz
      sini2_old_2=+q3_old_2+q4_old_2
      cosi3_old_2=+q1_old_2+q2_old_2         ! -ky,-kz
      sini3_old_2=+q3_old_2-q4_old_2

      ! the same for the 3rd atom
      q1_old_3=cckxky_old_3*cckz_old_3(kz)
      q2_old_3=sskxky_old_3*sskz_old_3(kz)
      q3_old_3=sskxky_old_3*cckz_old_3(kz)
      q4_old_3=cckxky_old_3*sskz_old_3(kz)

      cosi_old_3=+q1_old_3-q2_old_3       ! ky,kz
      sini_old_3=+q3_old_3+q4_old_3
      cosi1_old_3=+q1_old_3+q2_old_3     ! ky,-kz
      sini1_old_3=+q3_old_3-q4_old_3

      q1_old_3=cckxkya_old_3*cckz_old_3(kz)
      q2_old_3=sskxkya_old_3*sskz_old_3(kz)
      q3_old_3=sskxkya_old_3*cckz_old_3(kz)
      q4_old_3=cckxkya_old_3*sskz_old_3(kz)

      cosi2_old_3=+q1_old_3-q2_old_3         ! -ky,kz
      sini2_old_3=+q3_old_3+q4_old_3
      cosi3_old_3=+q1_old_3+q2_old_3         ! -ky,-kz
      sini3_old_3=+q3_old_3-q4_old_3
!
      sumcosold(kk)=-2.d0*cosi_old_1+cosi_old_2+cosi_old_3 ! -2 is the oxygen charge, sumcosold = SUM_j q_j cos( k_m r_j )
      sumsinold(kk)=-2.d0*sini_old_1+sini_old_2+sini_old_3
  

      sumcosold_o(kk)=cosi_old_1   ! sum for oxygen only, for Lennard Jones ewald 
      sumsinold_o(kk)=sini_old_1   !

      ! sumcosold(sxsysz) = SUM_{j\in mol} q_i cos(sx x_j kx + sy y_j ky + sz z_j kz ) 
      ! sumsinold(sxsysz) = SUM_{j\in mol} q_i sin(sx x_j kx + sy y_j ky + sz z_j kz ) 

      cos_actuel=sumcos(kk)
      sin_actuel=sumsin(kk)

      wclb=xxclb*pi4*bet(kk)
      wlj=xlj*pi4*(sigma12*bet12(kk)-sigma6*bet6(kk))
   
      coeff=-2.d0*wclb*(-cosi_old_1*sin_actuel+sini_old_1*cos_actuel)      &
        +wlj*(-cosi_old_1*sumsin_o(kk)+sini_old_1*sumcos_o(kk))
  
!      write(*,*) 'cosi_old_1',cosi_old_1,'sin_actuel',sin_actuel 
     

      fxold(1)=fxold(1)+coeff*kx
      fyold(1)=fyold(1)+coeff*ky
      fzold(1)=fzold(1)+coeff*kz
   
      coeff=wclb*(-cosi_old_2*sin_actuel+sini_old_2*cos_actuel)
  
      fxold(2)=fxold(2)+coeff*kx
      fyold(2)=fyold(2)+coeff*ky
      fzold(2)=fzold(2)+coeff*kz
   
      coeff=wclb*(-cosi_old_3*sin_actuel+sini_old_3*cos_actuel)

      fxold(3)=fxold(3)+coeff*kx
      fyold(3)=fyold(3)+coeff*ky
      fzold(3)=fzold(3)+coeff*kz
!
      if (iq(kk).eq.3.or.iq(kk).eq.4) then        ! faire -kz  ppm

!          write(*,*) 'LUC ppm:',kk

         sumcos1old(kk)=-2.d0*cosi1_old_1+cosi1_old_2+cosi1_old_3
         sumsin1old(kk)=-2.d0*sini1_old_1+sini1_old_2+sini1_old_3

         sumcos1old_o(kk)=cosi1_old_1
         sumsin1old_o(kk)=sini1_old_1

         cos1_actuel=sumcos1(kk)
         sin1_actuel=sumsin1(kk)

         coeff=-2.d0*wclb*(-cosi1_old_1*sin1_actuel+sini1_old_1*cos1_actuel)      &
           +wlj*(-cosi1_old_1*sumsin1_o(kk)+sini1_old_1*sumcos1_o(kk))

         fxold(1)=fxold(1)+coeff*kx
         fyold(1)=fyold(1)+coeff*ky
         fzold(1)=fzold(1)-coeff*kz

         coeff=wclb*(-cosi1_old_2*sin1_actuel+sini1_old_2*cos1_actuel)

         fxold(2)=fxold(2)+coeff*kx
         fyold(2)=fyold(2)+coeff*ky
         fzold(2)=fzold(2)-coeff*kz

         coeff=wclb*(-cosi1_old_3*sin1_actuel+sini1_old_3*cos1_actuel)

         fxold(3)=fxold(3)+coeff*kx
         fyold(3)=fyold(3)+coeff*ky
         fzold(3)=fzold(3)-coeff*kz

      endif
!
      IF(iq(kk).eq.2.or.iq(kk).eq.4) then        ! faire -ky

         sumcos2old(kk)=-2.d0*cosi2_old_1+cosi2_old_2+cosi2_old_3
         sumsin2old(kk)=-2.d0*sini2_old_1+sini2_old_2+sini2_old_3

         sumcos2old_o(kk)=cosi2_old_1
         sumsin2old_o(kk)=sini2_old_1

         cos2_actuel=sumcos2(kk)
         sin2_actuel=sumsin2(kk)

         coeff=-2.d0*wclb*(-cosi2_old_1*sin2_actuel+sini2_old_1*cos2_actuel)      &
           +wlj*(-cosi2_old_1*sumsin2_o(kk)+sini2_old_1*sumcos2_o(kk))

         fxold(1)=fxold(1)+coeff*kx
         fyold(1)=fyold(1)-coeff*ky
         fzold(1)=fzold(1)+coeff*kz

         coeff=wclb*(-cosi2_old_2*sin2_actuel+sini2_old_2*cos2_actuel)

         fxold(2)=fxold(2)+coeff*kx
         fyold(2)=fyold(2)-coeff*ky
         fzold(2)=fzold(2)+coeff*kz

         coeff=wclb*(-cosi2_old_3*sin2_actuel+sini2_old_3*cos2_actuel)

         fxold(3)=fxold(3)+coeff*kx
         fyold(3)=fyold(3)-coeff*ky
         fzold(3)=fzold(3)+coeff*kz

      endif
!
      IF(iq(kk).eq.4) then        ! faire -ky,-kz

         sumcos3old(kk)=-2.d0*cosi3_old_1+cosi3_old_2+cosi3_old_3
         sumsin3old(kk)=-2.d0*sini3_old_1+sini3_old_2+sini3_old_3

         sumcos3old_o(kk)=cosi3_old_1
         sumsin3old_o(kk)=sini3_old_1

         cos3_actuel=sumcos3(kk)
         sin3_actuel=sumsin3(kk)

         coeff=-2.d0*wclb*(-cosi3_old_1*sin3_actuel+sini3_old_1*cos3_actuel)      &
           +wlj*(-cosi3_old_1*sumsin3_o(kk)+sini3_old_1*sumcos3_o(kk))

         fxold(1)=fxold(1)+coeff*kx
         fyold(1)=fyold(1)-coeff*ky
         fzold(1)=fzold(1)-coeff*kz

         coeff=wclb*(-cosi3_old_2*sin3_actuel+sini3_old_2*cos3_actuel)

         fxold(2)=fxold(2)+coeff*kx
         fyold(2)=fyold(2)-coeff*ky
         fzold(2)=fzold(2)-coeff*kz

         coeff=wclb*(-cosi3_old_3*sin3_actuel+sini3_old_3*cos3_actuel)

         fxold(3)=fxold(3)+coeff*kx
         fyold(3)=fyold(3)-coeff*ky
         fzold(3)=fzold(3)-coeff*kz

      endif
!
      if(ip(kk).eq.1.or.ip(kk).eq.2) then       ! fin de boucle en kz et peut-etre en ky

         if(ip(kk).eq.1) then
            ky=ky+1
         else
            ky=0
            kx=kx+1
         endif

         q1_old_1=cckx_old_1(kx)*ccky_old_1(ky)
         q2_old_1=sskx_old_1(kx)*ssky_old_1(ky)
         q3_old_1=sskx_old_1(kx)*ccky_old_1(ky)
         q4_old_1=cckx_old_1(kx)*ssky_old_1(ky)
 
         cckxky_old_1=q1_old_1-q2_old_1
         sskxky_old_1=q3_old_1+q4_old_1
         cckxkya_old_1=q1_old_1+q2_old_1
         sskxkya_old_1=q3_old_1-q4_old_1
 
         q1_old_2=cckx_old_2(kx)*ccky_old_2(ky)
         q2_old_2=sskx_old_2(kx)*ssky_old_2(ky)
         q3_old_2=sskx_old_2(kx)*ccky_old_2(ky)
         q4_old_2=cckx_old_2(kx)*ssky_old_2(ky)
 
         cckxky_old_2=q1_old_2-q2_old_2
         sskxky_old_2=q3_old_2+q4_old_2
         cckxkya_old_2=q1_old_2+q2_old_2
         sskxkya_old_2=q3_old_2-q4_old_2
 
         q1_old_3=cckx_old_3(kx)*ccky_old_3(ky)
         q2_old_3=sskx_old_3(kx)*ssky_old_3(ky)
         q3_old_3=sskx_old_3(kx)*ccky_old_3(ky)
         q4_old_3=cckx_old_3(kx)*ssky_old_3(ky)

         cckxky_old_3=q1_old_3-q2_old_3
         sskxky_old_3=q3_old_3+q4_old_3
         cckxkya_old_3=q1_old_3+q2_old_3
         sskxkya_old_3=q3_old_3-q4_old_3

         kz=-1
      endif

   end do  ! fin kk

endif !  (i_partie==2.or.i_partie==0)

!tmpstr = 'sumcosold.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile,sumcosold,kkmax)
!call io_close(hfile)

!tmpstr = 'sumsinold.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile,sumsinold,kkmax)
!call io_close(hfile)


!write(*,*) 'fxold real+fourier'
!call write_xyz_array(0,fxold,fyold,fzold,3)


!PRINT*, 'old en r et k '
!PRINT*, fxold,fyold,fzold
!
!        et enfin les forces en epsilon'
!
!write(*,*) 'n=',n
n = 3 * n_h2o

summux=0.;summuy=0.; summuz=0.
do j=1,n

   jo=3*((j-1)/3)+1

   chargj=+1.d0
   if(j==jo) chargj=-2.d0

   summux=summux+chargj*xx(j)                      ! rq: j'ignore ici que H peut sortir de la boite
   summuy=summuy+chargj*yy(j)
   summuz=summuz+chargj*zz(j)

end do


coeffx=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*summux
coeffy=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*summuy
coeffz=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*summuz

do i=1,3

   chi=+1.d0
   if(i==1) chi=-2.d0

   fxold(i)=fxold(i)+coeffx*chi
   fyold(i)=fyold(i)+coeffy*chi
   fzold(i)=fzold(i)+coeffz*chi
end do




!
!          fin calcul des 3 fi
!
fxtotold=SUM(fxold(1:3))                    ! force totale sur la molecule H2O
fytotold=SUM(fyold(1:3))
fztotold=SUM(fzold(1:3))
ftotold=SQRT(fxtotold**2+fytotold**2+fztotold**2)



!
!
!
!        deplacement dans une sphere de rayon dr  avec proba exp(lamda.F.deltar)/Cold
!        luc70p194 et luc84p160
!
!
!
dr=dr_a/xl_a            ! car dr_a saisi dans mcluc2 est en unite A

umax=xlambda_f*ftotold*dr

if(umax.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   !        d'abord choisir la norme deltar suivant deltar*sh(lambda*Fi*deltar)
   e=EXP(-2.d0*umax)
   t=(umax*(1.d0+e)-(1.d0-e))/2.d0                 ! on gere les overflows luc70p197
   cold=3.d0*t/umax**3                         ! a exp(Umax) pres
   xksi=alea(1)

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

   do k=1,10                         ! qq NR
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
   xksi=alea(1)
   xcostet=1.d0+LOG(xksi+(1.d0-xksi)*e)/u
   
   if(u<1.d-7) xcostet=1.d0-2.d0*(1.d0-xksi)
   !      IF(ABS(xcostet)>1) PRINT*, 'xcostet:',u,e,xksi,xcostet,1.d0-xcostet**2

   xsintet=SQRT(1.d0-xcostet**2)
   !       vecteur au hasard en direction
   w=-1.d0
   
   do while(w<=0.)

      v2=2.
      do WHILE(v2>1.)                ! vecteur v qq
         vx=2.*alea(1)-1.
         vy=2.*alea(1)-1.
         vz=2.*alea(1)-1.
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
        
!   write(*,*) 'umax.le.1.d-2'
 
   v2=2.
   do while(v2>1.)                ! vecteur v qq
      vx=2.*alea(1)-1.
      vy=2.*alea(1)-1.
      vz=2.*alea(1)-1.
      v2=vx**2+vy**2+vz**2
   end do

   deltx=dr*vx; delty=dr*vy; deltz=dr*vz
   umax=0.
   cold=1.d0

!   write(*,*) 'delta xyz',deltx,delty,deltz

endif ! (umax.gt.1.d-2) 

umaxold=umax
!
xnew(1)=xx(i0+1)+deltx                     ! deplacement de O
ynew(1)=yy(i0+1)+delty
znew(1)=zz(i0+1)+deltz

xnew(1)=xnew(1)-ANINT(xnew(1))         ! O qui doit rester dans la boite
ynew(1)=ynew(1)-ANINT(ynew(1))         ! again: coordinates in box units (1 is size of the box)
znew(1)=znew(1)-ANINT(znew(1))         ! => this is the nearliest neighbour law

!
!
!       rotation d'un angle alpha choisi dans -d_angle,+d_angle
!       par rapport a un vecteur V
!       avec proba exp(lamda*C.V*alpha)  ou C=couple
!       luc84p157,162,185
!
!
call prod_vect(xx(i0+2)-xx(i0+1),yy(i0+2)-yy(i0+1),zz(i0+2)-zz(i0+1),fxold(2),fyold(2),fzold(2),couplx2,couply2,couplz2)
call prod_vect(xx(i0+3)-xx(i0+1),yy(i0+3)-yy(i0+1),zz(i0+3)-zz(i0+1),fxold(3),fyold(3),fzold(3),couplx3,couply3,couplz3)

couplxold=couplx2+couplx3                   ! couple total en ajoutant les 2 contributions des 2 H
couplyold=couply2+couply3
couplzold=couplz2+couplz3


couplold=SQRT(couplxold**2+couplyold**2+couplzold**2)
!          d'abord sortir x=cos(angle entre V et C) entre 0 et 1


vmax=xlambda_c*couplold*d_angle

if(vmax.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   !       il faut tirer v entre 0 et vmax suivant shv/v
   !       pas possible, on tire donc d'abord suivant chv

   emax=EXP(-2.d0*vmax)
   shmax=(1.d0-emax)/2.d0  ! sans exp(vmax)
   xksi=1.; proba=0.

   do while(xksi>=proba)

      xksi=alea(1)
      xksishmax=xksi*shmax
      v=vmax+LOG(xksishmax+SQRT(xksishmax**2+emax))   ! v suivant chv
      if(xksishmax**2<1.d-14*emax) v=xksishmax/SQRT(emax)

      e=EXP(-2.d0*v)
      proba=(1.d0-e)/v/(1.d0+e)       ! accepte avec proba shv/v / chv
      if(v<1.d-5) proba=(2.d0*v-2.d0*v**2)/v/(1.d0+e)

      xksi=alea(1)
   end do

   cosx=v/vmax                     ! cos de l'angle entre C et V
   qold=shi_sans_exp_f(vmax)/vmax  ! je ne mets pas le facteur exp(vmax) pour eviter les overflows
   !        puis sortir l'angle de rotation alpha entre -d_angle et d_angle suivant exp(lamda*C.V*alpha)
   xksi=alea(1)
   angle=d_angle*(1.d0+LOG(xksi+(1.d0-xksi)*e)/v)

   if(v<1.d-7) angle=d_angle*(1.d0-2.d0*(1.d0-xksi))
   !    determiner l'orientation absolue de V qui fait l'angle x avec C
   !    vecteur au hasard en direction
   w=-1.d0
   do while(w<=0.)
  
      v2=2.
      do WHILE(v2>1.d0)                ! vecteur v qq
         vx=2.d0*alea(1)-1.d0
         vy=2.d0*alea(1)-1.d0
         vz=2.d0*alea(1)-1.d0
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
   do WHILE(v2>1.d0)                ! vecteur V qq
      vx=2.d0*alea(1)-1.d0
      vy=2.d0*alea(1)-1.d0
      vz=2.d0*alea(1)-1.d0
      v2=vx**2+vy**2+vz**2
   END do

   v=SQRT(v2)
   vx=vx/v; vy=vy/v; vz=vz/v
   angle=d_angle*(2.d0*alea(1)-1.d0)
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

!
do i=2,3                          ! rotation des 2 H      luc84p157
   xr=xx(i0+i)-xx(i0+1)
   yr=yy(i0+i)-yy(i0+1)
   zr=zz(i0+i)-zz(i0+1)

   xrnew=(cphi**2*st**2*ca1+ca)*xr+(sphi*cphi*st**2*ca1-ct*sa)*yr+st*(ct*ca1*cphi+sa*sphi)*zr
   yrnew=(cphi*sphi*st**2*ca1+ct*sa)*xr+(sphi**2*st**2*ca1+ca)*yr+st*(ct*ca1*sphi-sa*cphi)*zr
   zrnew=st*(ct*ca1*cphi-sa*sphi)*xr+st*(ct*ca1*sphi+sa*cphi)*yr+(st**2*ca+ct**2)*zr

   xnew(i)=xnew(1)+xrnew
   ynew(i)=ynew(1)+yrnew
   znew(i)=znew(1)+zrnew
end do


!PRINT*, i_h2o
!PRINT*, 'old'
!PRINT*, xx(i0+1:i0+3),yy(i0+1:i0+3),zz(i0+1:i0+3)
!PRINT*, fxtotold,fytotold,fztotold
!PRINT*, couplxold,couplyold,couplzold
!PRINT*, 'deplct'
!PRINT*, xcostet
!PRINT*, deltx,delty,deltz
!PRINT*, vx,vy,vz,angle
!PRINT*, 'new'
!PRINT*, xnew,ynew,znew
!PRINT*, i_h2o
!PRINT*, deltx,delty,deltz
!PRINT*, vx,vy,vz,angle,st,cphi,sphi
!
!
!            ca y est, la position new des 3 atomes de H2O est connue!
!            il faut alors calculer les nouvelles energies et forces!
!
!
vrnew(1:n)=0.                               ! pour tous les j sauf ceux de i_h2o
frxnew(1:n)=0.; frynew(1:n)=0.; frznew(1:n)=0.        ! et leur force r
uinew(1:3)=0.                               ! pour les 3 i de i_h2o
fxnew(1:3)=0.; fynew(1:3)=0.; fznew(1:3)=0.        ! et leur force totale (r+k+epsilon')
!
!             d'abord LJ+coulomb Ewald en r
!        je ne compte pas l'intra comme pour uu car c'est inchange lors du deplacement  luc84p164
!
IF(i_partie==3.or.i_partie==0) then
   potlj12new=0.; potlj6new=0.           ! servira pour la pression
   potewnew=0.
   do i=1,3                        ! les 3 atomes de la molecule deplacee

      if(i==1) chargi=-2.d0
      if(i>1) chargi=+1.d0

      do j=1,n                  ! tous les sites

         jo=3*((j-1)/3)+1
         if(jo==iio) cycle        ! j'ignore l'intra!    (comme pour uu)
 
         chargj=+1.d0
         if(j==jo) chargj=-2.d0
 
         xij=xnew(i)-xx(j)
         yij=ynew(i)-yy(j)                ! donc compris entre -1 et 1
         zij=znew(i)-zz(j)
 
         if(xij.gt.0.5) xij=xij-1.
         if(xij.lt.-0.5) xij=xij+1.
 
         if(yij.gt.0.5) yij=yij-1.
         if(yij.lt.-0.5) yij=yij+1.         ! prend moins de temps que anint!
 
         if(zij.gt.0.5) zij=zij-1.
         if(zij.lt.-0.5) zij=zij+1.
 
         r2=xij**2+yij**2+zij**2
 
         IF(i==1.and.j==jo.and.r2<0.4*sigma2) return        ! overlap OO LJ, r<0.63*sigma:  on sort

         IF(r2<rmax2) then

            r=SQRT(r2)
            alpr=alp*r
            alpr2=alpr**2
            ee2=EXP(-alpr2)
            er=erfc_Luc(alpr,alpr2,ee2)            ! calcule erfc_Luc(x) en s'aidant de x**2 et exp(-x**2)
   !         IF(jo==io) er=er-1.d0    ! retrancher le self eventuellement  luc84p134
            xclbzizj=xxclb*chargi*chargj
            potew=xclbzizj*er/r
            pot=potew
            potewnew=potewnew+potew
            ffr=xclbzizj*(er/r+alp2pi*ee2)/r2

            IF(i==1.and.j==jo) then  ! OO: LJ
               vv=(sigma2/r2)**3
               c6=1.d0+alpr2*(1.d0+alpr2/2.d0)
               coeff=xlj*vv*ee2
               potlj6=-coeff*c6
               alpr6=alpr2**3
               c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0
               potlj12=coeff*vv*c12
               potlj12new=potlj12new+potlj12; potlj6new=potlj6new+potlj6
               pot=pot+potlj12+potlj6
               ffr=ffr+6.d0*coeff*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2

               !     pour test, LJ tronque shifte (xlj a ete mis a 0)
               !            potlj6=-xlj0*vv
               !            potlj12=xlj0*vv**2
               !            potlj12new=potlj12new+potlj12; potlj6new=potlj6new+potlj6
               !            pot=pot+potlj12+potlj6
               !            ffr=ffr+xlj0*6.d0*vv*(2.d0*vv-1.d0)/r2
            endif

            vrnew(j)=vrnew(j)+pot               ! a conserver en cas de deplt accepte
            uinew(i)=uinew(i)+pot
   
            fxij=xij*ffr
            fxnew(i)=fxnew(i)+fxij
            frxnew(j)=frxnew(j)-fxij            ! idem
            fyij=yij*ffr
            fynew(i)=fynew(i)+fyij
            frynew(j)=frynew(j)-fyij
            fzij=zij*ffr
            fznew(i)=fznew(i)+fzij
            frznew(j)=frznew(j)-fzij
         endif ! r2<rmax2

      end do  ! j = 1,n
   end do ! i = 1,3

endif ! (i_partie==3.or.i_partie==0)

!PRINT*, 'new en r'
!PRINT*, fxnew,fynew,fznew
!        noter la contribution r de fnew pour remplacer eventuellement fx,y,z si deplt accepte

fxnew_r(1:3)=fxnew(1:3); fynew_r(1:3)=fynew(1:3); fznew_r(1:3)=fznew(1:3)

!
!
dutot=sum(uinew(1:3)-uu(i0+1:i0+3))        ! delta u en r pour la molecule i_h2o
     charge_h = 0.42380d0



!
!        puis l'energie Ewald k  astuce luc70p35
!        il faut Zi*exp(2ipi/L*ri*k)  luc79p160
!        pour les 3 i O,H,H qui composent la molecule deplacee i_h2o
!
!        je fais les 3 i new a l'interieur de la meme boucle en kk !!!
!        et delta Uk en meme temps
!        luc84p181,183
!

dutotk=0.
dutotk6=0.; dutotk12=0.
!
IF(i_partie==4.or.i_partie==0) then

   cx_new_1=COS(pi2*xnew(1))
   sx_new_1=SIN(pi2*xnew(1))
   cy_new_1=COS(pi2*ynew(1))
   sy_new_1=SIN(pi2*ynew(1))
   cz_new_1=COS(pi2*znew(1))
   sz_new_1=SIN(pi2*znew(1))

   cx_new_2=COS(pi2*xnew(2))
   sx_new_2=SIN(pi2*xnew(2))
   cy_new_2=COS(pi2*ynew(2))
   sy_new_2=SIN(pi2*ynew(2))
   cz_new_2=COS(pi2*znew(2))
   sz_new_2=SIN(pi2*znew(2))

   cx_new_3=COS(pi2*xnew(3))
   sx_new_3=SIN(pi2*xnew(3))
   cy_new_3=COS(pi2*ynew(3))
   sy_new_3=SIN(pi2*ynew(3))
   cz_new_3=COS(pi2*znew(3))
   sz_new_3=SIN(pi2*znew(3))

   cckx_new_1(0)=1.
   sskx_new_1(0)=0.
   ccky_new_1(0)=1.
   ssky_new_1(0)=0.
   cckz_new_1(0)=1.
   sskz_new_1(0)=0.

   cckx_new_2(0)=1.
   sskx_new_2(0)=0.
   ccky_new_2(0)=1.
   ssky_new_2(0)=0.
   cckz_new_2(0)=1.
   sskz_new_2(0)=0.

   cckx_new_3(0)=1.
   sskx_new_3(0)=0.
   ccky_new_3(0)=1.
   ssky_new_3(0)=0.
   cckz_new_3(0)=1.
   sskz_new_3(0)=0.

   do k=1,kmax
      cckx_new_1(k)=cckx_new_1(k-1)*cx_new_1-sskx_new_1(k-1)*sx_new_1
      sskx_new_1(k)=sskx_new_1(k-1)*cx_new_1+cckx_new_1(k-1)*sx_new_1
      ccky_new_1(k)=ccky_new_1(k-1)*cy_new_1-ssky_new_1(k-1)*sy_new_1                  ! tableaux suivant les 3D
      ssky_new_1(k)=ssky_new_1(k-1)*cy_new_1+ccky_new_1(k-1)*sy_new_1
      cckz_new_1(k)=cckz_new_1(k-1)*cz_new_1-sskz_new_1(k-1)*sz_new_1
      sskz_new_1(k)=sskz_new_1(k-1)*cz_new_1+cckz_new_1(k-1)*sz_new_1

      cckx_new_2(k)=cckx_new_2(k-1)*cx_new_2-sskx_new_2(k-1)*sx_new_2
      sskx_new_2(k)=sskx_new_2(k-1)*cx_new_2+cckx_new_2(k-1)*sx_new_2
      ccky_new_2(k)=ccky_new_2(k-1)*cy_new_2-ssky_new_2(k-1)*sy_new_2                  ! tableaux suivant les 3D
      ssky_new_2(k)=ssky_new_2(k-1)*cy_new_2+ccky_new_2(k-1)*sy_new_2
      cckz_new_2(k)=cckz_new_2(k-1)*cz_new_2-sskz_new_2(k-1)*sz_new_2
      sskz_new_2(k)=sskz_new_2(k-1)*cz_new_2+cckz_new_2(k-1)*sz_new_2

      cckx_new_3(k)=cckx_new_3(k-1)*cx_new_3-sskx_new_3(k-1)*sx_new_3
      sskx_new_3(k)=sskx_new_3(k-1)*cx_new_3+cckx_new_3(k-1)*sx_new_3
      ccky_new_3(k)=ccky_new_3(k-1)*cy_new_3-ssky_new_3(k-1)*sy_new_3                  ! tableaux suivant les 3D
      ssky_new_3(k)=ssky_new_3(k-1)*cy_new_3+ccky_new_3(k-1)*sy_new_3
      cckz_new_3(k)=cckz_new_3(k-1)*cz_new_3-sskz_new_3(k-1)*sz_new_3
      sskz_new_3(k)=sskz_new_3(k-1)*cz_new_3+cckz_new_3(k-1)*sz_new_3
   end do

   cckxky_new_1=1.d0                               ! kx,ky
   sskxky_new_1=0.
   cckxkya_new_1=1.d0                              ! kx,-ky
   sskxkya_new_1=0.

   cckxky_new_2=1.d0                               ! kx,ky
   sskxky_new_2=0.
   cckxkya_new_2=1.d0                              ! kx,-ky
   sskxkya_new_2=0.

   cckxky_new_3=1.d0                               ! kx,ky
   sskxky_new_3=0.
   cckxkya_new_3=1.d0                              ! kx,-ky
   sskxkya_new_3=0.

   kx=0
   ky=0
   kz=0
!
!      debut de la somme en kk
!
   do kk=1,kkmax                           ! debut k
      kz=kz+1
   
      q1_new_1=cckxky_new_1*cckz_new_1(kz)
      q2_new_1=sskxky_new_1*sskz_new_1(kz)
      q3_new_1=sskxky_new_1*cckz_new_1(kz)
      q4_new_1=cckxky_new_1*sskz_new_1(kz)

      cosi_new_1=+q1_new_1-q2_new_1       ! ky,kz
      sini_new_1=+q3_new_1+q4_new_1
      cosi1_new_1=+q1_new_1+q2_new_1     ! ky,-kz
      sini1_new_1=+q3_new_1-q4_new_1

      q1_new_1=cckxkya_new_1*cckz_new_1(kz)
      q2_new_1=sskxkya_new_1*sskz_new_1(kz)
      q3_new_1=sskxkya_new_1*cckz_new_1(kz)
      q4_new_1=cckxkya_new_1*sskz_new_1(kz)

      cosi2_new_1=+q1_new_1-q2_new_1         ! -ky,kz
      sini2_new_1=+q3_new_1+q4_new_1
      cosi3_new_1=+q1_new_1+q2_new_1         ! -ky,-kz
      sini3_new_1=+q3_new_1-q4_new_1

      q1_new_2=cckxky_new_2*cckz_new_2(kz)
      q2_new_2=sskxky_new_2*sskz_new_2(kz)
      q3_new_2=sskxky_new_2*cckz_new_2(kz)
      q4_new_2=cckxky_new_2*sskz_new_2(kz)

      cosi_new_2=+q1_new_2-q2_new_2       ! ky,kz
      sini_new_2=+q3_new_2+q4_new_2
      cosi1_new_2=+q1_new_2+q2_new_2     ! ky,-kz
      sini1_new_2=+q3_new_2-q4_new_2

      q1_new_2=cckxkya_new_2*cckz_new_2(kz)
      q2_new_2=sskxkya_new_2*sskz_new_2(kz)
      q3_new_2=sskxkya_new_2*cckz_new_2(kz)
      q4_new_2=cckxkya_new_2*sskz_new_2(kz)

      cosi2_new_2=+q1_new_2-q2_new_2         ! -ky,kz
      sini2_new_2=+q3_new_2+q4_new_2
      cosi3_new_2=+q1_new_2+q2_new_2         ! -ky,-kz
      sini3_new_2=+q3_new_2-q4_new_2

      q1_new_3=cckxky_new_3*cckz_new_3(kz)
      q2_new_3=sskxky_new_3*sskz_new_3(kz)
      q3_new_3=sskxky_new_3*cckz_new_3(kz)
      q4_new_3=cckxky_new_3*sskz_new_3(kz)

      cosi_new_3=+q1_new_3-q2_new_3       ! ky,kz
      sini_new_3=+q3_new_3+q4_new_3
      cosi1_new_3=+q1_new_3+q2_new_3     ! ky,-kz
      sini1_new_3=+q3_new_3-q4_new_3

      q1_new_3=cckxkya_new_3*cckz_new_3(kz)
      q2_new_3=sskxkya_new_3*sskz_new_3(kz)
      q3_new_3=sskxkya_new_3*cckz_new_3(kz)
      q4_new_3=cckxkya_new_3*sskz_new_3(kz)

      cosi2_new_3=+q1_new_3-q2_new_3         ! -ky,kz
      sini2_new_3=+q3_new_3+q4_new_3
      cosi3_new_3=+q1_new_3+q2_new_3         ! -ky,-kz
      sini3_new_3=+q3_new_3-q4_new_3
!
      dcos=-2.d0*cosi_new_1+cosi_new_2+cosi_new_3-sumcosold(kk)             ! delta M(h2o)
      dsin=-2.d0*sini_new_1+sini_new_2+sini_new_3-sumsinold(kk)

      sumcosold(kk)=dcos                             ! dorenavant, contient dif=new-old!
      sumsinold(kk)=dsin

      cos_actuel=sumcos(kk)
      sin_actuel=sumsin(kk)

      d_som=(2.d0*cos_actuel+dcos)*dcos+(2.d0*sin_actuel+dsin)*dsin         ! pour delta Uk
      dcos_o=cosi_new_1-sumcosold_o(kk)
      dsin_o=sini_new_1-sumsinold_o(kk)

      sumcosold_o(kk)=dcos_o                         ! idem, c'est dif!
      sumsinold_o(kk)=dsin_o

      cos_actuel_o=sumcos_o(kk)
      sin_actuel_o=sumsin_o(kk)

      d_som_o=(2.d0*cos_actuel_o+dcos_o)*dcos_o+(2.d0*sin_actuel_o+dsin_o)*dsin_o

      wclb=xxclb*pi4*bet(kk)
      wlj=xlj*pi4*(sigma12*bet12(kk)-sigma6*bet6(kk))

      coeff=-2.d0*wclb*(-cosi_new_1*(sin_actuel+dsin)+sini_new_1*(cos_actuel+dcos))      &
        +wlj*(-cosi_new_1*(sin_actuel_o+dsin_o)+sini_new_1*(cos_actuel_o+dcos_o))



      fxnew(1)=fxnew(1)+coeff*kx
      fynew(1)=fynew(1)+coeff*ky
      fznew(1)=fznew(1)+coeff*kz

      coeff=wclb*(-cosi_new_2*(sin_actuel+dsin)+sini_new_2*(cos_actuel+dcos))


      fxnew(2)=fxnew(2)+coeff*kx
      fynew(2)=fynew(2)+coeff*ky
      fznew(2)=fznew(2)+coeff*kz

      coeff=wclb*(-cosi_new_3*(sin_actuel+dsin)+sini_new_3*(cos_actuel+dcos))


      fxnew(3)=fxnew(3)+coeff*kx
      fynew(3)=fynew(3)+coeff*ky
      fznew(3)=fznew(3)+coeff*kz
!
      if(iq(kk).eq.3.or.iq(kk).eq.4) then        ! faire -kz

         dcos1=-2.d0*cosi1_new_1+cosi1_new_2+cosi1_new_3-sumcos1old(kk)
         dsin1=-2.d0*sini1_new_1+sini1_new_2+sini1_new_3-sumsin1old(kk)



         sumcos1old(kk)=dcos1                             ! dorenavant, contient dif=new-old!
         sumsin1old(kk)=dsin1

         cos1_actuel=sumcos1(kk)
         sin1_actuel=sumsin1(kk)

         d_som=d_som+(2.d0*cos1_actuel+dcos1)*dcos1+(2.d0*sin1_actuel+dsin1)*dsin1
         dcos1_o=cosi1_new_1-sumcos1old_o(kk)
         dsin1_o=sini1_new_1-sumsin1old_o(kk)

         sumcos1old_o(kk)=dcos1_o                         ! idem, c'est dif!
         sumsin1old_o(kk)=dsin1_o

         cos1_actuel_o=sumcos1_o(kk)
         sin1_actuel_o=sumsin1_o(kk)

         d_som_o=d_som_o+(2.d0*cos1_actuel_o+dcos1_o)*dcos1_o+(2.d0*sin1_actuel_o+dsin1_o)*dsin1_o

         coeff=-2.d0*wclb*(-cosi1_new_1*(sin1_actuel+dsin1)+sini1_new_1*(cos1_actuel+dcos1))      &
           +wlj*(-cosi1_new_1*(sin1_actuel_o+dsin1_o)+sini1_new_1*(cos1_actuel_o+dcos1_o))


         fxnew(1)=fxnew(1)+coeff*kx
         fynew(1)=fynew(1)+coeff*ky
         fznew(1)=fznew(1)-coeff*kz

         coeff=wclb*(-cosi1_new_2*(sin1_actuel+dsin1)+sini1_new_2*(cos1_actuel+dcos1))

         fxnew(2)=fxnew(2)+coeff*kx
         fynew(2)=fynew(2)+coeff*ky
         fznew(2)=fznew(2)-coeff*kz

         coeff=wclb*(-cosi1_new_3*(sin1_actuel+dsin1)+sini1_new_3*(cos1_actuel+dcos1))

         fxnew(3)=fxnew(3)+coeff*kx
         fynew(3)=fynew(3)+coeff*ky
         fznew(3)=fznew(3)-coeff*kz
      endif
!
      if(iq(kk).eq.2.or.iq(kk).eq.4) then        ! faire -ky
         dcos2=-2.d0*cosi2_new_1+cosi2_new_2+cosi2_new_3-sumcos2old(kk)
         dsin2=-2.d0*sini2_new_1+sini2_new_2+sini2_new_3-sumsin2old(kk)
      
         sumcos2old(kk)=dcos2                             ! dorenavant, contient dif=new-old!
         sumsin2old(kk)=dsin2

         cos2_actuel=sumcos2(kk)
         sin2_actuel=sumsin2(kk)

         d_som=d_som+(2.d0*cos2_actuel+dcos2)*dcos2+(2.d0*sin2_actuel+dsin2)*dsin2
         dcos2_o=cosi2_new_1-sumcos2old_o(kk)
         dsin2_o=sini2_new_1-sumsin2old_o(kk)

         sumcos2old_o(kk)=dcos2_o                         ! idem, c'est dif!
         sumsin2old_o(kk)=dsin2_o

         cos2_actuel_o=sumcos2_o(kk)
         sin2_actuel_o=sumsin2_o(kk)

         d_som_o=d_som_o+(2.d0*cos2_actuel_o+dcos2_o)*dcos2_o+(2.d0*sin2_actuel_o+dsin2_o)*dsin2_o

         coeff=-2.d0*wclb*(-cosi2_new_1*(sin2_actuel+dsin2)+sini2_new_1*(cos2_actuel+dcos2))      &
           +wlj*(-cosi2_new_1*(sin2_actuel_o+dsin2_o)+sini2_new_1*(cos2_actuel_o+dcos2_o))



         fxnew(1)=fxnew(1)+coeff*kx
         fynew(1)=fynew(1)-coeff*ky
         fznew(1)=fznew(1)+coeff*kz

         coeff=wclb*(-cosi2_new_2*(sin2_actuel+dsin2)+sini2_new_2*(cos2_actuel+dcos2))

         fxnew(2)=fxnew(2)+coeff*kx
         fynew(2)=fynew(2)-coeff*ky
         fznew(2)=fznew(2)+coeff*kz

         coeff=wclb*(-cosi2_new_3*(sin2_actuel+dsin2)+sini2_new_3*(cos2_actuel+dcos2))

         fxnew(3)=fxnew(3)+coeff*kx
         fynew(3)=fynew(3)-coeff*ky
         fznew(3)=fznew(3)+coeff*kz
      endif
!
      if(iq(kk).eq.4) then        ! faire -ky,-kz
         dcos3=-2.d0*cosi3_new_1+cosi3_new_2+cosi3_new_3-sumcos3old(kk)
         dsin3=-2.d0*sini3_new_1+sini3_new_2+sini3_new_3-sumsin3old(kk)

         sumcos3old(kk)=dcos3                             ! dorenavant, contient dif=new-old!
         sumsin3old(kk)=dsin3

         cos3_actuel=sumcos3(kk)
         sin3_actuel=sumsin3(kk)

         d_som=d_som+(2.d0*cos3_actuel+dcos3)*dcos3+(2.d0*sin3_actuel+dsin3)*dsin3
         dcos3_o=cosi3_new_1-sumcos3old_o(kk)
         dsin3_o=sini3_new_1-sumsin3old_o(kk)

         sumcos3old_o(kk)=dcos3_o                         ! idem, c'est dif!
         sumsin3old_o(kk)=dsin3_o

         cos3_actuel_o=sumcos3_o(kk)
         sin3_actuel_o=sumsin3_o(kk)

         d_som_o=d_som_o+(2.d0*cos3_actuel_o+dcos3_o)*dcos3_o+(2.d0*sin3_actuel_o+dsin3_o)*dsin3_o

         coeff=-2.d0*wclb*(-cosi3_new_1*(sin3_actuel+dsin3)+sini3_new_1*(cos3_actuel+dcos3))      &
           +wlj*(-cosi3_new_1*(sin3_actuel_o+dsin3_o)+sini3_new_1*(cos3_actuel_o+dcos3_o))



         fxnew(1)=fxnew(1)+coeff*kx
         fynew(1)=fynew(1)-coeff*ky
         fznew(1)=fznew(1)-coeff*kz

         coeff=wclb*(-cosi3_new_2*(sin3_actuel+dsin3)+sini3_new_2*(cos3_actuel+dcos3))

         fxnew(2)=fxnew(2)+coeff*kx
         fynew(2)=fynew(2)-coeff*ky
         fznew(2)=fznew(2)-coeff*kz

         coeff=wclb*(-cosi3_new_3*(sin3_actuel+dsin3)+sini3_new_3*(cos3_actuel+dcos3))

         fxnew(3)=fxnew(3)+coeff*kx
         fynew(3)=fynew(3)-coeff*ky
         fznew(3)=fznew(3)-coeff*kz
      endif
!
      dutotk=dutotk+bet(kk)*d_som
      dutotk6=dutotk6+bet6(kk)*d_som_o
      dutotk12=dutotk12+bet12(kk)*d_som_o
!
      if(ip(kk).eq.1.or.ip(kk).eq.2) then       ! fin de boucle en kz et peut-etre en ky

         IF(ip(kk).eq.1) then
            ky=ky+1
         else
            ky=0
            kx=kx+1
         endif

         q1_new_1=cckx_new_1(kx)*ccky_new_1(ky)
         q2_new_1=sskx_new_1(kx)*ssky_new_1(ky)
         q3_new_1=sskx_new_1(kx)*ccky_new_1(ky)
         q4_new_1=cckx_new_1(kx)*ssky_new_1(ky)

         cckxky_new_1=q1_new_1-q2_new_1
         sskxky_new_1=q3_new_1+q4_new_1
         cckxkya_new_1=q1_new_1+q2_new_1
         sskxkya_new_1=q3_new_1-q4_new_1

         q1_new_2=cckx_new_2(kx)*ccky_new_2(ky)
         q2_new_2=sskx_new_2(kx)*ssky_new_2(ky)
         q3_new_2=sskx_new_2(kx)*ccky_new_2(ky)
         q4_new_2=cckx_new_2(kx)*ssky_new_2(ky)

         cckxky_new_2=q1_new_2-q2_new_2
         sskxky_new_2=q3_new_2+q4_new_2
         cckxkya_new_2=q1_new_2+q2_new_2
         sskxkya_new_2=q3_new_2-q4_new_2

         q1_new_3=cckx_new_3(kx)*ccky_new_3(ky)
         q2_new_3=sskx_new_3(kx)*ssky_new_3(ky)
         q3_new_3=sskx_new_3(kx)*ccky_new_3(ky)
         q4_new_3=cckx_new_3(kx)*ssky_new_3(ky)

         cckxky_new_3=q1_new_3-q2_new_3
         sskxky_new_3=q3_new_3+q4_new_3
         cckxkya_new_3=q1_new_3+q2_new_3
         sskxkya_new_3=q3_new_3-q4_new_3

         kz=-1
      endif

   end do ! fin kk

endif ! (i_partie==4.or.i_partie==0)



dpotew=potewnew+xxclb*dutotk                     ! on n'a pas encore retranche old
dpotlj6=potlj6new-xlj*sigma6*dutotk6
dpotlj12=potlj12new+xlj*sigma12*dutotk12
dutot=dutot+xxclb*dutotk-xlj*sigma6*dutotk6+xlj*sigma12*dutotk12       ! la, par contre, old a ete retranche

!PRINT*, dutot
!PRINT*, 'new en r et en k'
!PRINT*, fxnew,fynew,fznew
!
!        et enfin les energie et forces en epsilon'
!
dsummux=0.;dsummuy=0.; dsummuz=0.      ! somme de zi*deltari POUR la molecule i_h2o
do i=1,3
   IF(i==1) chargi=-2.d0
   IF(i>1) chargi=+1.d0

   dsummux=dsummux+chargi*(xnew(i)-xx(i0+i))
   dsummuy=dsummuy+chargi*(ynew(i)-yy(i0+i))
   dsummuz=dsummuz+chargi*(znew(i)-zz(i0+i))
end do

dureste=pi2*xxclb/(2.d0*cdiel_ext+1.d0)*      &                   ! delta M**2   luc84p165
  ((2.d0*summux+dsummux)*dsummux+(2.d0*summuy+dsummuy)*dsummuy+(2.d0*summuz+dsummuz)*dsummuz)
dpotew=dpotew+dureste
dutot=dutot+dureste


coeffx=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux+dsummux)        ! c'est Mnew qui compte pour la force new
coeffy=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summuy+dsummuy)
coeffz=-2.d0*pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summuz+dsummuz)

do i=1,3
   chi=+1.d0
   IF(i==1) chi=-2.d0

   fxnew(i)=fxnew(i)+coeffx*chi
   fynew(i)=fynew(i)+coeffy*chi
   fznew(i)=fznew(i)+coeffz*chi
end do
!
!                  les 3 fi new sont maintenant totalement connues!
!
!
!                         on doit calculer Cnew et Qnew
!
fxtotnew=SUM(fxnew(1:3))                    ! force totale sur la molecule H2O
fytotnew=SUM(fynew(1:3))
fztotnew=SUM(fznew(1:3))
ftotnew=SQRT(fxtotnew**2+fytotnew**2+fztotnew**2)


umaxnew=xlambda_f*ftotnew*dr
IF(umaxnew.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   e=EXP(-2.d0*umaxnew)
   t=(umaxnew*(1.d0+e)-(1.d0-e))/2.d0                 ! on gere les overflows luc70p197
   cnew=3.d0*t/umaxnew**3                         ! a exp(Umax) pres
else
   umaxnew=0.
   cnew=1.d0
endif
!
call prod_vect(xnew(2)-xnew(1),ynew(2)-ynew(1),znew(2)-znew(1),fxnew(2),fynew(2),fznew(2),couplx2,couply2,couplz2)
call prod_vect(xnew(3)-xnew(1),ynew(3)-ynew(1),znew(3)-znew(1),fxnew(3),fynew(3),fznew(3),couplx3,couply3,couplz3)

couplxnew=couplx2+couplx3                   ! couple total
couplynew=couply2+couply3
couplznew=couplz2+couplz3
couplnew=SQRT(couplxnew**2+couplynew**2+couplznew**2)


vmaxnew=xlambda_c*couplnew*d_angle

IF(vmaxnew.gt.1.d-2) then   ! a ne faire que si vraiment biaise
   qnew=shi_sans_exp_f(vmaxnew)/vmaxnew  ! je ne mets pas le facteur exp(vmax) pour eviter les overflows
else
   vmaxnew=0.
   qnew=1.d0
endif

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
xksi=alea(1)

!write(*,*) 'LUC: xksi',xksi,'factor:',factor
!   write(*,*) 'cold',cold,'qold',qold,'cnew',cnew,'qnew',qnew
!   write(*,*) 'exposu',exposu,'exposf',exposf
!   write(*,*) 'exposfcorr',exposfcorr,'exposc',exposc,'exposccorr',exposccorr


! write(*,*) 'LUC: du:',dutot,'F:', fxtotnew,fytotnew,fztotnew,&
!                                'T:',couplxnew,couplynew,couplznew
 


IF(xksi>factor) return                  ! NON!!! on n'accepte pas, et on sort!
IF(i_partie/=0) return          ! on n'accepte pas si test de vitesse partiel
!
!
!        et donc OUI, on accepte le deplacement!
!
!
iaccep=1
!PRINT*, 'accepte'
!
uutot=uutot+dutot
uu(1:n)=uu(1:n)+vrnew(1:n)

fx(1:n)=fx(1:n)+frxnew(1:n); fy(1:n)=fy(1:n)+frynew(1:n); fz(1:n)=fz(1:n)+frznew(1:n)

uu(i0+1:i0+3)=uinew(1:3)
fx(i0+1:i0+3)=fxnew_r(1:3); fy(i0+1:i0+3)=fynew_r(1:3); fz(i0+1:i0+3)=fznew_r(1:3)
!
!             il faut retrancher l'energie et la force old r pour chaque j
!
potlj12old=0.; potlj6old=0.           ! servira pour la pression
potewold=0.
!        je ne compte pas l'intra car c'est inchange lors du deplacement  luc84p164
IF(i_partie==1.or.i_partie==0) then

   do i=i0+1,i0+3

      if(i==i0+1) chargi=-2.d0
      if(i>i0+1) chargi=+1.d0

      do j=1,n

         jo=3*((j-1)/3)+1
         IF(jo==iio) cycle        ! j'ignore l'intra!    (comme pour vrnew)

         chargj=+1.d0
         IF(j==jo) chargj=-2.d0

         xij=xx(i)-xx(j)
         yij=yy(i)-yy(j)                ! donc compris entre -1 et 1
         zij=zz(i)-zz(j)

         if(xij.gt.0.5) xij=xij-1.
         if(xij.lt.-0.5) xij=xij+1.
         if(yij.gt.0.5) yij=yij-1.
         if(yij.lt.-0.5) yij=yij+1.         ! prend moins de temps que anint!
         if(zij.gt.0.5) zij=zij-1.
         if(zij.lt.-0.5) zij=zij+1.

         r2=xij**2+yij**2+zij**2

         IF(r2<rmax2) then
            r=SQRT(r2)
            alpr=alp*r
            alpr2=alpr**2
            ee2=EXP(-alpr2)

            er=erfc_Luc(alpr,alpr2,ee2)            ! calcule erfc_Luc(x) en s'aidant de x**2 et exp(-x**2)
              !         IF(jo==io) er=er-1.d0    ! retrancher le self eventuellement  luc84p134

            xclbzizj=xxclb*chargi*chargj
            potew=xclbzizj*er/r
            pot=potew
            potewold=potewold+potew
            ffr=xclbzizj*(er/r+alp2pi*ee2)/r2

            IF(i==iio.and.j==jo) then   ! OO: LJ
               vv=(sigma2/r2)**3
               c6=1.d0+alpr2*(1.d0+alpr2/2.d0)
               coeff=xlj*vv*ee2
               potlj6=-coeff*c6
               alpr6=alpr2**3
               c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0
               potlj12=coeff*vv*c12
               potlj12old=potlj12old+potlj12; potlj6old=potlj6old+potlj6
               pot=pot+potlj12+potlj6
               ffr=ffr+6.d0*coeff*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2
               !     pour test, LJ tronque shifte (xlj a ete mis a 0)
               !            potlj6=-xlj0*vv
               !            potlj12=xlj0*vv**2
               !            potlj12old=potlj12old+potlj12; potlj6old=potlj6old+potlj6
               !            pot=pot+potlj12+potlj6
               !            ffr=ffr+xlj0*6.d0*vv*(2.d0*vv-1.d0)/r2
            endif

            uu(j)=uu(j)-pot
            fx(j)=fx(j)+xij*ffr            ! on retranche -xij*ffr
            fy(j)=fy(j)+yij*ffr
            fz(j)=fz(j)+zij*ffr

         endif !(r2<rmax2) 

      end do  ! j=1,n
   end do ! i=i0+1,i0+3
endif ! (i_partie==1.or.i_partie==0)

!PRINT*, 'old en r: '
!PRINT*, fxold,fyold,fzold
dpotew=dpotew-potewold
dpotlj6=dpotlj6-potlj6old
dpotlj12=dpotlj12-potlj12old
uu_lj12=uu_lj12+dpotlj12; uu_lj6=uu_lj6+dpotlj6
uu_ew=uu_ew+dpotew
!
IF(i_partie==6.or.i_partie==0) then
   !        et en k
   !        attention: on a mis la difference new-old dans sumcosold!  luc84p183
   do kk=1,kkmax
  
      sumcos(kk)=sumcos(kk)+sumcosold(kk)
      sumsin(kk)=sumsin(kk)+sumsinold(kk)
  
      sumcos_o(kk)=sumcos_o(kk)+sumcosold_o(kk)
      sumsin_o(kk)=sumsin_o(kk)+sumsinold_o(kk)
  
      IF(iq(kk).eq.3.or.iq(kk).eq.4) then            ! kz>0
         sumcos1(kk)=sumcos1(kk)+sumcos1old(kk)
         sumsin1(kk)=sumsin1(kk)+sumsin1old(kk)
         sumcos1_o(kk)=sumcos1_o(kk)+sumcos1old_o(kk)
         sumsin1_o(kk)=sumsin1_o(kk)+sumsin1old_o(kk)
      endif
 
      IF(iq(kk).eq.2.or.iq(kk).eq.4) then              ! ky>0
         sumcos2(kk)=sumcos2(kk)+sumcos2old(kk)
         sumsin2(kk)=sumsin2(kk)+sumsin2old(kk)
         sumcos2_o(kk)=sumcos2_o(kk)+sumcos2old_o(kk)
         sumsin2_o(kk)=sumsin2_o(kk)+sumsin2old_o(kk)
      endif                                            ! ky,kz>0

      IF(iq(kk).eq.4) then
         sumcos3(kk)=sumcos3(kk)+sumcos3old(kk)
         sumsin3(kk)=sumsin3(kk)+sumsin3old(kk)
         sumcos3_o(kk)=sumcos3_o(kk)+sumcos3old_o(kk)
         sumsin3_o(kk)=sumsin3_o(kk)+sumsin3old_o(kk)
      endif

   end do ! kk
endif ! i_partie==6.or.i_partie==0
!
xx(i0+1:i0+3)=xnew(1:3)
yy(i0+1:i0+3)=ynew(1:3)
zz(i0+1:i0+3)=znew(1:3)


!
return
end
!
!
subroutine prod_vect(x1,y1,z1,x2,y2,z2,x,y,z) !DOC
IMPLICIT REAL(8) (a-h,o-z)
!
!DOC        fait le produit vectoriel r = r1 * r2
!DOC        luc84 page 160
!DOC        le resultat peut-etre mis a la place du 1er vecteur!
!
x3=x1; y3=y1; z3=z1
x=y3*z2-z3*y2
y=z3*x2-x3*z2
z=x3*y2-y3*x2

end subroutine
!
!
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

function alea(n) !DOC
!DOC random number
!DOC Parmaters:
   use MRandom, only : rand,randomize
   real(8) :: alea
   integer,intent(in) :: n !DOC code. 0=initialize, 1-give the next random

   if ( n == 0 ) then
!      call randomize(0)
      alea = 0
   elseif ( n == 1 ) then
      alea = rand()  
   else
      alea = 0
   end if

end function

!     Last change:  LB   22 May 2014   11:52 pm
      subroutine piston(iaccep,press_a3,dlnvmax,xx,yy,zz,n_h2o,                        systematic_correction) !DOC
      use constants, only : pi
!DOC Luc's volume step 
!DOC Parameters:
      implicit REAL(8) (a-h,o-z)


  
      integer,intent(out) :: iaccep !DOC whether accepted or not
      real(8),intent(in) :: press_a3  !DOC pressure in kT/A^3
      real(8),intent(in) :: dlnvmax !DOC maximal value for ln V 
      real(8),dimension(:),intent(inout) :: xx,yy,zz !DOC coordinates. re-scale if needed
      integer,intent(in) :: n_h2o !DOC number of h2o molecules
      logical,intent(in),optional :: systematic_correction !DOC do the systematic correction or not

! DP maintenant      real * 4 alea
   common/xenergy/utotk_old,utotk6_old,utotk12_old,uu_lj6_r_old,uu_lj12_r_old
   common/dipole/summux,summuy,summuz

      common/nombre/n
!      common/positi/xx(10000),yy(10000),zz(10000)        ! positions cellule centrale
      common/numeri/p,num
      common/energi/uu_old(10000),uu(10000),uutot_old
      COMMON/energies/uu_lj_old,uu_lj12_old,uu_lj6_old,uu_ew_old,uu_ew_r_old
      COMMON/forces/fx_old(10000),fy_old(10000),fz_old(10000)
      common/conc/conc_h2o,xl_a_old
!      common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d
      common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d

      COMMON/lennardjones/epsi_lj,xlj,ucut                         ! xlj0 pour test
      COMMON/ewald/temp,dbjr,charg_h,xxclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
!COMMON/ewald/temp,dbjr,charg_h,xclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
!COMMON/ewald/charg_h,xxclb,alp_a1,alp,kkmax,cdiel_ext

      COMMON/ewald_s/ssr,ssk                            ! pour bien sauvegarder sk!
!      COMMON/isobare/press_pa,press_a3,rho1,rho2,vol1,vol2,vol_sph_cub(2000)
!      COMMON/chgt_vol/proba_vol,dlnvmax
      PARAMETER(nkmax=4661,kmaxmax=20)
      COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
       sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
       sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
       cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
       cckz(0:kmaxmax),sskz(0:kmaxmax)
      COMMON/sommeenk_lj/bet6(nkmax),bet12(nkmax),                     &         ! pour Ewald LJ en k
       sumcos_o(nkmax),sumsin_o(nkmax),sumcos1_o(nkmax),sumsin1_o(nkmax),        &
       sumcos2_o(nkmax),sumsin2_o(nkmax),sumcos3_o(nkmax),sumsin3_o(nkmax)
      COMMON/sommeenkold/sumcosold(nkmax),sumsinold(nkmax),sumcos1old(nkmax),sumsin1old(nkmax),        &
       sumcos2old(nkmax),sumsin2old(nkmax),sumcos3old(nkmax),sumsin3old(nkmax)
      COMMON/sommeenkold_lj/sumcosold_o(nkmax),sumsinold_o(nkmax),sumcos1old_o(nkmax),sumsin1old_o(nkmax),        &
       sumcos2old_o(nkmax),sumsin2old_o(nkmax),sumcos3old_o(nkmax),sumsin3old_o(nkmax)
DIMENSION fx(10000),fy(10000),fz(10000)
!
!        sous-programme Monte-Carlo NPT pour H2O
!
!        on envisage des changements de volume delta(lnV)   luc87p148
!        il faut changer V et les positions et recalculer U et tous les tableaux intermediaires comme dans inmc
!        mais auparavant, sauver les anciens tableaux pour les restaurer si changement de volume refuse
!
!!
!        version H2O luc84p122
!
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        version Force-Bias  luc84p160
!
!        version "new" qui ruse en k pour reduire le nombre de tableaux et accelerer en k
!        luc84p183
!
!        j'actualise les forces en r pour eviter la boucle old r si deplt refuse  luc84p187
!
!n_h2o=n/3
n = n_h2o * 3
!
!          on choisit dlnV uniformement entre -dlnVmax et +dlnVmax
!
dlnv=(2.d0*alea(1)-1.d0)*dlnvmax
xl_a=xl_a_old*EXP(dlnv/3.d0)
!                                     nouveaux parametres normalises
sigma=sigma_a/xl_a
sigma2=sigma**2
roh=roh_a/xl_a

!pi=4.d0*ATAN(1.d0)


dbjr=dbjr*xl_a_old/xl_a
xxclb=charg_h**2*dbjr
theta=pi/180.d0*theta_d
theta2=theta/2.d0
cthet2=COS(theta2); sthet2=SIN(theta2)
rhh=2.d0*roh*sthet2
      pi2=2.d0*pi
      alp2pi=2.d0*alp/dSQRT(pi)
      ssk_vrai=kmax*pi/alp
!          les si=ri/L sont inchanges pour les Oxygenes, par contre sOH doit etre renormalise par le nouveau L
do i_h2o=1,n_h2o
i0=(i_h2o-1)*3+1
do i=i0+1,i0+2
xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
rcs=dSQRT(xcs**2+ycs**2+zcs**2)
xx(i)=xx(i0)+roh*xcs/rcs
yy(i)=yy(i0)+roh*ycs/rcs
zz(i)=zz(i0)+roh*zcs/rcs
end do
end do
!
!        je sauvegarde les tableaux actuels
!
sumcosold(1:kkmax)=sumcos(1:kkmax)
sumsinold(1:kkmax)=sumsin(1:kkmax)
sumcos1old(1:kkmax)=sumcos1(1:kkmax)
sumsin1old(1:kkmax)=sumsin1(1:kkmax)
sumcos2old(1:kkmax)=sumcos2(1:kkmax)
sumsin2old(1:kkmax)=sumsin2(1:kkmax)
sumcos3old(1:kkmax)=sumcos3(1:kkmax)
sumsin3old(1:kkmax)=sumsin3(1:kkmax)
sumcosold_o(1:kkmax)=sumcos_o(1:kkmax)
sumsinold_o(1:kkmax)=sumsin_o(1:kkmax)
sumcos1old_o(1:kkmax)=sumcos1_o(1:kkmax)
sumsin1old_o(1:kkmax)=sumsin1_o(1:kkmax)
sumcos2old_o(1:kkmax)=sumcos2_o(1:kkmax)
sumsin2old_o(1:kkmax)=sumsin2_o(1:kkmax)
sumcos3old_o(1:kkmax)=sumcos3_o(1:kkmax)
sumsin3old_o(1:kkmax)=sumsin3_o(1:kkmax)
!
!        on calcule alors comme a la fin de inmc
!
!
!        on initialise les energies et les forces en r SANS intra luc84p187
!        on calcule maintenant LJ total par Ewald aussi!
!
      uu(1:n)=0.                         ! repere LJr+EWr
      fx(1:n)=0.; fy(1:n)=0.; fz(1:n)=0.
      uu_lj12=0.; uu_lj6=0.
      uu_ew=0.
summux=0.;summuy=0.; summuz=0.
!        pour test, je reviens a LJ tronque shifte a rcut=rmax  luc85p123   NON, plus maintenant
!xlj0=xlj; xlj=0.
        do i=1,n
        io=3*((i-1)/3)+1
        chargi=+1.d0
        IF(i==io) chargi=-2.d0
summux=summux+chargi*xx(i)                      ! rq: j'ignore ici que H peut sortir de la boite
summuy=summuy+chargi*yy(i)
summuz=summuz+chargi*zz(i)
        do j=i+1,n
        jo=3*((j-1)/3)+1
        chargj=+1.d0
        IF(j==jo) chargj=-2.d0
        xij=xx(i)-xx(j)
        yij=yy(i)-yy(j)
        zij=zz(i)-zz(j)
        xij=xij-anint(xij)
        yij=yij-anint(yij)
        zij=zij-anint(zij)
        r2=xij**2+yij**2+zij**2
  IF(i==io.and.j==jo.and.r2<0.4*sigma2) GOTO 1000        ! overlap OO LJ, r<0.63*sigma:  on sort
         IF(r2<rmax2) then                         ! rmax pour Ewald r
         r=SQRT(r2)
         alpr=alp*r
         alpr2=alpr**2
         ee2=EXP(-alpr2)
         er=erfc_Luc(alpr,alpr2,ee2)            ! calcule erfc(x) en s'aidant de x**2 et exp(-x**2)
         IF(io==jo) er=-erf_Luc(alpr,alpr2,ee2)    ! retrancher le self eventuellement donc erfc-1=-erf  luc84p134
         potew=xxclb*chargi*chargj*er/r
         uu_ew=uu_ew+potew
         pot=potew
         ffr=xxclb*chargi*chargj*(er/r+alp2pi*ee2)/r2
            IF(i==io.and.j==jo) then            ! OO: LJ total  luc84p147
            vv=(sigma2/r2)**3
            c6=1.d0+alpr2*(1.d0+alpr2/2.d0)
            potlj6=-xlj*vv*ee2*c6
            alpr6=alpr2**3
            c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0
            potlj12=xlj*vv**2*ee2*c12
            uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
            pot=pot+potlj12+potlj6
            ffr=ffr+xlj*6.d0*vv*ee2*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2
!     pour test, LJ tronque shifte (xlj a ete mis a 0)         NON
!            potlj6=-xlj0*vv
!            potlj12=xlj0*vv**2
!            uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
!            pot=pot+potlj12+potlj6
!            ffr=ffr+xlj0*6.d0*vv*(2.d0*vv-1.d0)/r2
            endif
           IF(io/=jo) then
         uu(i)=uu(i)+pot                  ! ne contient plus le self
         uu(j)=uu(j)+pot                  ! mais contient l'intra!
         fxij=xij*ffr                     ! NON, plus maintenant!
         fx(i)=fx(i)+fxij
         fx(j)=fx(j)-fxij
         fyij=yij*ffr
         fy(i)=fy(i)+fyij
         fy(j)=fy(j)-fyij
         fzij=zij*ffr
         fz(i)=fz(i)+fzij
         fz(j)=fz(j)-fzij
           endif
         endif
        END do
        END do
uu_ew_r=uu_ew
!PRINT*, 'uu_lj,r,uu_ew_r',uu_lj_r,uu_ew_r
uu_lj6_r = uu_lj6
uu_lj12_r = uu_lj12
!
!        puis l'energie Ewald k  astuce luc70p35
!        il faut somme de Zi*exp(2ipi/L*ri*k)  luc70p35
!        et somme sur les i=O uniquement de exp()  luc84p148
      do kk=1,kkmax
      sumcos(kk)=0.
      sumsin(kk)=0.
      sumcos1(kk)=0.
      sumsin1(kk)=0.
      sumcos2(kk)=0.
      sumsin2(kk)=0.
      sumcos3(kk)=0.
      sumsin3(kk)=0.
      sumcos_o(kk)=0.
      sumsin_o(kk)=0.
      sumcos1_o(kk)=0.
      sumsin1_o(kk)=0.
      sumcos2_o(kk)=0.
      sumsin2_o(kk)=0.
      sumcos3_o(kk)=0.
      sumsin3_o(kk)=0.
      end do
!
       do i=1,n                                   ! debut i
!
       io=3*((i-1)/3)+1
       chi=+1.d0
       IF(i==io) chi=-2.d0
       cx=COS(pi2*xx(i))
       sx=SIN(pi2*xx(i))
       cy=COS(pi2*yy(i))
       sy=SIN(pi2*yy(i))
       cz=COS(pi2*zz(i))
       sz=SIN(pi2*zz(i))
       cckx(0)=1.
       sskx(0)=0.
       ccky(0)=1.
       ssky(0)=0.
       cckz(0)=1.
       sskz(0)=0.
        do k=1,kmax
        cckx(k)=cckx(k-1)*cx-sskx(k-1)*sx
        sskx(k)=sskx(k-1)*cx+cckx(k-1)*sx
        ccky(k)=ccky(k-1)*cy-ssky(k-1)*sy                  ! tableaux suivant les 3D
        ssky(k)=ssky(k-1)*cy+ccky(k-1)*sy
        cckz(k)=cckz(k-1)*cz-sskz(k-1)*sz
        sskz(k)=sskz(k-1)*cz+cckz(k-1)*sz
        end do
      kx=0
      ky=0
      kz=0
      cckxky=chi                               ! kx,ky
      sskxky=0.
      cckxkya=chi                              ! kx,-ky
      sskxkya=0.
!
        do kk=1,kkmax                           ! debut k
        kz=kz+1
        q1=cckxky*cckz(kz)
        q2=sskxky*sskz(kz)
        q3=sskxky*cckz(kz)
        q4=cckxky*sskz(kz)
        sumcos(kk)=sumcos(kk)+q1-q2       ! ky,kz
        sumsin(kk)=sumsin(kk)+q3+q4
        sumcos1(kk)=sumcos1(kk)+q1+q2     ! ky,-kz
        sumsin1(kk)=sumsin1(kk)+q3-q4
           IF(i==io) then                 ! pour O, ajouter LJ
        sumcos_o(kk)=sumcos_o(kk)+q1-q2       ! ky,kz                     contient zi=-2!
        sumsin_o(kk)=sumsin_o(kk)+q3+q4
        sumcos1_o(kk)=sumcos1_o(kk)+q1+q2     ! ky,-kz
        sumsin1_o(kk)=sumsin1_o(kk)+q3-q4
           endif
        q1=cckxkya*cckz(kz)
        q2=sskxkya*sskz(kz)
        q3=sskxkya*cckz(kz)
        q4=cckxkya*sskz(kz)
        sumcos2(kk)=sumcos2(kk)+q1-q2         ! -ky,kz
        sumsin2(kk)=sumsin2(kk)+q3+q4
        sumcos3(kk)=sumcos3(kk)+q1+q2         ! -ky,-kz
        sumsin3(kk)=sumsin3(kk)+q3-q4
           IF(i==io) then                 ! pour O, ajouter LJ
        sumcos2_o(kk)=sumcos2_o(kk)+q1-q2       ! ky,kz
        sumsin2_o(kk)=sumsin2_o(kk)+q3+q4
        sumcos3_o(kk)=sumcos3_o(kk)+q1+q2     ! ky,-kz
        sumsin3_o(kk)=sumsin3_o(kk)+q3-q4
           endif
          IF(ip(kk).eq.1.or.ip(kk).eq.2) then       ! fin de boucle en kz et peut-etre en ky
        IF(ip(kk).eq.1) then
        ky=ky+1
        else
        ky=0
        kx=kx+1
        endif
        q1=cckx(kx)*ccky(ky)*chi
        q2=sskx(kx)*ssky(ky)*chi
        q3=sskx(kx)*ccky(ky)*chi
        q4=cckx(kx)*ssky(ky)*chi
        cckxky=q1-q2
        sskxky=q3+q4
        cckxkya=q1+q2
        sskxkya=q3-q4
        kz=-1
          endif
        end do                                  ! fin kk
!
       end do                                   ! fin i
!                                                                attention, les _o contiennent zi=zO=-2
      sumcos_o(1:kkmax)=sumcos_o(1:kkmax)/(-2.d0)
      sumsin_o(1:kkmax)=sumsin_o(1:kkmax)/(-2.d0)
      sumcos1_o(1:kkmax)=sumcos1_o(1:kkmax)/(-2.d0)
      sumsin1_o(1:kkmax)=sumsin1_o(1:kkmax)/(-2.d0)
      sumcos2_o(1:kkmax)=sumcos2_o(1:kkmax)/(-2.d0)
      sumsin2_o(1:kkmax)=sumsin2_o(1:kkmax)/(-2.d0)
      sumcos3_o(1:kkmax)=sumcos3_o(1:kkmax)/(-2.d0)
      sumsin3_o(1:kkmax)=sumsin3_o(1:kkmax)/(-2.d0)
!
      utotk=0.
      utotk6=0.; utotk12=0.
       do kk=1,kkmax
       som=sumcos(kk)**2+sumsin(kk)**2
       IF(iq(kk).eq.3.or.iq(kk).eq.4) som=som+sumcos1(kk)**2+sumsin1(kk)**2        ! kz>0                           ! kz>0
       IF(iq(kk).eq.2.or.iq(kk).eq.4) som=som+sumcos2(kk)**2+sumsin2(kk)**2        ! ky>0                           ! ky>0
       IF(iq(kk).eq.4) som=som+sumcos3(kk)**2+sumsin3(kk)**2                       ! ky,kz>0                           ! ky,kz>0
       utotk=utotk+bet(kk)*som
       som_o=sumcos_o(kk)**2+sumsin_o(kk)**2
       IF(iq(kk).eq.3.or.iq(kk).eq.4) som_o=som_o+sumcos1_o(kk)**2+sumsin1_o(kk)**2        ! kz>0                           ! kz>0
       IF(iq(kk).eq.2.or.iq(kk).eq.4) som_o=som_o+sumcos2_o(kk)**2+sumsin2_o(kk)**2        ! ky>0                           ! ky>0
       IF(iq(kk).eq.4) som_o=som_o+sumcos3_o(kk)**2+sumsin3_o(kk)**2                       ! ky,kz>0                           ! ky,kz>0
       utotk6=utotk6+bet6(kk)*som_o
       utotk12=utotk12+bet12(kk)*som_o
       end do
uu_ew=uu_ew+xxclb*utotk
ureste=pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)-alp2pi*n*xxclb
uu_ew=uu_ew+ureste
uu_lj6=uu_lj6-xlj*sigma2**3*(utotk6+pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
uu_lj12=uu_lj12+xlj*sigma2**6*(utotk12+pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)

!      write(*,*) 'utotk6',utotk6,'utotk12',utotk12
!      write(*,*) 'xlj',xlj,'sigma2',sigma2,'alp',alp,'n_h2o',n_h2o
!      write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew
!      je peux essayer des corrections systematiques  luc84p151

   if ( systematic_correction ) then
     qmax=2.d0*pi*kmax
      uu_ew=uu_ew+n_h2o/pi*xxclb*alp*EXP(-ssk_vrai**2)/ssk_vrai*(6.d0-8.d0*SIN(qmax*roh)/(qmax*roh)+2.d0*SIN(qmax*rhh)/(qmax*rhh))
      uu_lj6=uu_lj6-epsi_lj*sigma**6*         &
         (2.d0*pi*n_h2o**2*alp**3*EXP(-ssr**2)/ssr+2.d0/SQRT(pi)*n_h2o*alp**6*EXP(-ssk_vrai**2)/ssk_vrai)
      uu_lj12=uu_lj12+epsi_lj*sigma**12*      &
         (pi/30.d0*n_h2o**2*alp**9*EXP(-ssr**2)/ssr+1.d0/30.d0/SQRT(pi)*n_h2o*alp**12*EXP(-ssk_vrai**2)/ssk_vrai)
     
 !     PRINT*, 'et apres corrections des queues en r et en k:'
 !           write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew
     !!!       uself deja retranchee dans la boucle Ewald r
   end if

uu_lj=uu_lj12+uu_lj6
uutot=uu_lj+uu_ew
!      write(*,*) 'Energie totale actuelle : ',uutot
!
!
!       probabilite d'acceptation = (Vnew/Vold)**(N+1) * exp(-beta*P*(Vnew-Vold)-beta(Unew-Uold))
vold=xl_a_old**3
vnew=xl_a**3
hold=uutot_old+press_a3*vold
hnew=uutot+press_a3*vnew
proba=(vnew/vold)**(n_h2o+1)*EXP(-(hnew-hold))
IF(proba<1.d0) then
IF(alea(1)>proba) GOTO 1000
endif
!
!        on accepte le changement de volume!
!
!write(*,*) '------------> on accepte le changement de volume!'

p=p*xl_a_old/xl_a       ! contient p/Lactuel
xl_a_old=xl_a
avog=6.02214e-4
conc_h2o=n_h2o/xl_a**3/avog
uu_old(1:n)=uu(1:n)
uutot_old=uutot
uu_lj_old=uu_lj
uu_lj12_old=uu_lj12
uu_lj6_old=uu_lj6
uu_ew_old=uu_ew
uu_ew_r_old=uu_ew_r
fx_old(1:n)=fx(1:n)
fy_old(1:n)=fy(1:n)
fz_old(1:n)=fz(1:n)


utotk_old = utotk
utotk6_old = utotk6
utotk12_old = utotk12
uu_lj6_r_old = uu_lj6_r
uu_lj12_r_old = uu_lj12_r

iaccep=1
return
!
!        on refuse le changement de volume!
!
1000 continue
!
!write(*,*) '------------> on n''accepte pas le changement de volume!'
!write(*,*) 'uutot',uutot,'uutot_old',uutot_old

xl_a_new=xl_a
xl_a=xl_a_old
!                                     retour aux anciens parametres normalises
sigma=sigma_a/xl_a
sigma2=sigma**2
roh=roh_a/xl_a
dbjr=dbjr*xl_a_new/xl_a
xxclb=charg_h**2*dbjr
rhh=2.d0*roh*sthet2
!          les si=ri/L sont inchanges pour les Oxygenes, par contre sOH doit etre renormalise par le nouveau L
do i_h2o=1,n_h2o
i0=(i_h2o-1)*3+1
do i=i0+1,i0+2
xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
rcs=SQRT(xcs**2+ycs**2+zcs**2)
xx(i)=xx(i0)+roh*xcs/rcs
yy(i)=yy(i0)+roh*ycs/rcs
zz(i)=zz(i0)+roh*zcs/rcs
end do
end do
!
!        je recupere les tableaux actuels
!
sumcos(1:kkmax)=sumcosold(1:kkmax)
sumsin(1:kkmax)=sumsinold(1:kkmax)
sumcos1(1:kkmax)=sumcos1old(1:kkmax)
sumsin1(1:kkmax)=sumsin1old(1:kkmax)
sumcos2(1:kkmax)=sumcos2old(1:kkmax)
sumsin2(1:kkmax)=sumsin2old(1:kkmax)
sumcos3(1:kkmax)=sumcos3old(1:kkmax)
sumsin3(1:kkmax)=sumsin3old(1:kkmax)
sumcos_o(1:kkmax)=sumcosold_o(1:kkmax)
sumsin_o(1:kkmax)=sumsinold_o(1:kkmax)
sumcos1_o(1:kkmax)=sumcos1old_o(1:kkmax)
sumsin1_o(1:kkmax)=sumsin1old_o(1:kkmax)
sumcos2_o(1:kkmax)=sumcos2old_o(1:kkmax)
sumsin2_o(1:kkmax)=sumsin2old_o(1:kkmax)
sumcos3_o(1:kkmax)=sumcos3old_o(1:kkmax)
sumsin3_o(1:kkmax)=sumsin3old_o(1:kkmax)
iaccep=0
return
!
!
!
      end
!
!
subroutine write_energy_LUC(n_h2o) !DOC 
!DOC Write the energy calculated by Luc
   use constants, only : pi
!DOC Parameters:
!DOC :: n_h2o  - number of h2o molecules
   implicit REAL(8) (a-h,o-z)

   common/nombre/n
   COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r
   COMMON/lennardjones/epsi_lj,xlj,ucut
   common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d
   common/xenergy/utotk,utotk6,utotk12,uu_lj6_r,uu_lj12_r
   COMMON/ewald/temp,dbjr,charg_h,xxclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
   common/dipole/summux,summuy,summuz

   n = n_h2o * 3

   pi2 = 2.d0*pi
   alp2pi=2.d0*alp/dSQRT(pi)

   uu_lj_r = uu_lj12_r + uu_lj6_r

   uu_lj6_k = xlj*sigma2**3 * utotk6
   uu_lj12_k = xlj*sigma2**6 * utotk12
   uu_lj_k = uu_lj12_k - uu_lj6_k
   
   uu_ext = pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)
   
   !write(*,*) 'alp2pi=',alp2pi,'alpha',alpha,'alpha/sqrt(pi)',alpha/sqrt(pi)
   uu_tail_ew =  -alp2pi*n*xxclb  ! alp2pi = 2 alp / sqrt(pi)
   ! alp2pi * n *xxclb = xclb * q^2 6 alpha n_h2o / sqrt(pi) 
   ! SUM q_s = 2 n_h2o + n_h2o * 4 = 6 n_h2o
   uu_tail_lj6 = xlj*sigma2**3*(pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
   uu_tail_lj12 = xlj*sigma2**6*(pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)
   uu_tail_lj = uu_tail_lj12 - uu_tail_lj6
   
   write(*,*) '******************** LUC *******************************'
   write(*,*) 'U_REAL ew:',uu_ew_r,'lj:',uu_lj_r,'lj6:',-uu_lj6_r,'lj12:',uu_lj12_r
   write(*,*) 'U_FOUR ew:',xxclb*utotk,'lj:',uu_lj_k,'lj6:',uu_lj6_k,'lj12:',uu_lj12_k
   write(*,*) 'U_TAIL ew:',uu_tail_ew,'lj:',uu_tail_lj,'lj6:',uu_tail_lj6,'lj12:',uu_tail_lj12
   write(*,*) 'U_EXT  ew:',uu_ext
   write(*,*) 'U_TOTL ew:',uu_ew,'lj:',uu_lj12 + uu_lj6,'lj6:',-uu_lj6,'lj12:',uu_lj12
   write(*,*) 'U_TOT tot:',uu_ew + uu_lj6 + uu_lj12

end subroutine 

      subroutine initmc_Luc_OLD(xx,yy,zz,n_h2o,sigma2_in,xlj)
USE lecture_fichier
use parameters
use constants
use io, only : write_real_array, write_xyz_array,io_open,io_close
use SystemSettings, only : SYSTEM_STRING_LENGTH
!use io
      implicit REAL(8) (a-h,o-z)
      integer * 4 nittot,naccep,ncumul
      REAL * 8 iacc,iacc1,iacc2,iacc3,iacc_o
      real(8), dimension(:),intent(in) :: xx,yy,zz   
      integer, intent(in) :: n_h2o
      real(8),intent(in) :: sigma2_in
      real(8),intent(in) :: xlj ! 4 epsilon
      character(SYSTEM_STRING_LENGTH) :: tmpstr
      integer :: hfile 

! DP maintenant      real * 4 alea
      character * 50 nomfic,nomaaa,nomlec
!      common/nombre/n
!      common/positi/xx(10000),yy(10000),zz(10000)        ! positions
      common/iteri/nittot,naccep,ncumul
      common/numeri/p,num
      common/energi/uu(10000),vrnew(10000),uutot
      COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r
      COMMON/forces/fx(10000),fy(10000),fz(10000)
      common/moyenn/uutotm,duutot,evtest,devtes,mmu
      COMMON/nb_proj/mnmax,ialpmax
      common/cumul_1/iacc1(2000),iacc2(2000),iacc3(2000),    &  !gOO, OH, HH
      gr1(2000),gr2(2000),gr3(2000)
      PARAMETER(ialpmaxx=1226)                    ! 1226 pour nmax=6!
      PARAMETER(n_omega=8)
      COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)
      COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)
      common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190
      common/conc/conc_h2o,xl_a
!      common/diamet/sigma_a,rcutoff_a,rcut2,roh_a,theta_d
      common/diamet/sigma_a,sigma2,rcutoff_a,rcut2,roh_a,theta_d

!      COMMON/lennardjones/epsi_lj,xlj,ucut
!      COMMON/ewald/temp,dbjr,charg_h,xclb,alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
      COMMON/ewald_s/ssr,ssk                            ! pour bien sauvegarder sk!
      common/pression/vir,vir1,vir2,vir02,vir12,vir22,hypervir,hypervir2
      COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82
      COMMON/sauveaaa/nomaaa
      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
      common/gracan/activ,pinsdes,x_n,x_n2           ! pour GC
      PARAMETER(nkmax=4661,kmaxmax=20)
      COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
       sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
       sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
       cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
       cckz(0:kmaxmax),sskz(0:kmaxmax)
      COMMON/sommeenk_lj/bet6(nkmax),bet12(nkmax),                     &         ! pour Ewald LJ en k
       sumcos_o(nkmax),sumsin_o(nkmax),sumcos1_o(nkmax),sumsin1_o(nkmax),        &
       sumcos2_o(nkmax),sumsin2_o(nkmax),sumcos3_o(nkmax),sumsin3_o(nkmax)
      COMMON/numeriq/pasq,nptq,mq
      COMMON/spectre/sq1(500),sq2(500),sq3(500)
      COMMON/spectre1/sk1(500),sk2(500),sk3(500)
COMMON/inpout/inp
!
!        on initalise le programme Monte-Carlo pour H2O
!        luc belloni 84p121
!        on suit les forces en parallele luc84p133
!        finalement, pas necessaire luc84p137
!
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        pour la version FB, uu et fx,y,z contiennent la partie r SANS intra  luc84p187
!
!       on lit au clavier          luc84p92    et 85p130

!
!if(inp==0) inp=5
!!
!    1 format(1x,a)
!    2 format(1x,a,$)
!    3 format(1x,78(1h*))
!      pi=4.d0*ATAN(1.d0)
!      if(ichoixaaa.eq.0) ichoixaaa=1
!!       pour ecrire les valeurs de aaa en fichier
!      nomaaa='qq'
!      write(*,3)
!      write(*,3)
!      write(*,1) '          Monte-Carlo "luc" pour H2O"'
!      WRITE(*,1) '                   Modele SPCE'
!      WRITE(*,1) '                  LJ NON tronque!'
!      WRITE(*,1) '     utilise par des personnes jeunes, belles et dynamiques (et curieuses)'
!      write(*,3)
! 9999 write(*,3)
!      write(*,1) '                    MENU'
!      WRITE(*,1) '100--NNNEEEWWW!!! Clavier remplace par fichier d''entree'
!      write(*,1) '1--Parametres physiques'
!      write(*,1) '2--Positions de depart sur un reseau cfc'
!      write(*,1) '3--Lire en fichier les positions de depart'
!      write(*,1) '4--Parametres numeriques d''analyse'
!      write(*,1) '5--Lire en fichier les resulats deja accumules'
!      WRITE(*,1) '9--Variable pour estimer l''efficacite statistique'
!      WRITE(*,1) '10--Sauver les valeurs instantanees de la variable'
!      write(*,1) '11--"Lisser" des resultats accumules'
!      WRITE(*,1) '12--Cumuler 2 simulations'
!      write(*,2) '8--Demarrer les cycles Monte-Carlo  --->'
!      read(inp,*,err=9999) ic
!      write(*,3)
!!
!!
!      if(ic==100) then
!   !
!   !      on peut lire tous les paramètres en fichier           luc84p92
!   !
!         WRITE(*,2) 'Nom du fichier (pas de lecture si "qq") --->'
!         read(inp,'(A)') nomlec
!         if(nomlec=='stop') stop
!   
!         if(nomlec.ne.'qq'.and.nomlec.ne.'QQ') then
!              OPEN(8,FILE=nomlec,STATUS='old',ERR=9999)
!              inp=8
!         ELSE
!              CLOSE(8)
!              inp=5
!         END if ! not QQ
!
!      END if ! ic == 100
!!
!!
!!
!      goto(1000,2000,3000,4000,5000,9999,9999,8000,9600,9700,1100,1200),ic
!!        ic
!! 1000 : 1           1--Parametres physiques'
!! 2000 : 2           2--Positions de depart sur un reseau cfc'
!! 3000 : 3           3--Lire en fichier les positions de depart'
!! 4000 : 4           4--Parametres numeriques d''analyse'
!! 5000 : 5           5--Lire en fichier les resulats deja accumules'
!! 9999 : 6            -none-
!! 9999 : 7            -none-
!! 8000 : 8           8--Demarrer les cycles Monte-Carlo  --->'
!! 9600 : 9           9--Variable pour estimer l''efficacite statistique'
!! 9700 : 10          10--Sauver les valeurs instantanees de la variable'
!! 1100 : 11          11--"Lisser" des resultats accumules
!! 1200 : 12          12--Cumuler 2 simulations
!GOTO 9999
!!
!!        physique
!!
! 1000 write(*,*) 'Potentiel LJ TOTAL pour OO et charges partielles en O et H'
!write(*,3)
!write(*,*) 'La simulation Monte-Carlo va mettre n_h2o H2O dans un cube de taille L (qui sera prise comme unite de longueur)'
!write(*,3)
!write(*,2) 'Nombre de molecules d''eau n_h2o  --->'
!call lis_f(inp,n_h2o)
!write(*,2) 'Temperature (en K) --->'
!call lis_f(inp,temp)
!write(*,2) 'Concentration d''eau  (en M) --->'
!call lis_f(inp,conc_h2o)
!WRITE(*,3)
!WRITE(*,2) 'Diametre LJ de OO (en A) --->'
!call lis_f(inp,sigma_a)
!write(*,2) 'epsilon/kT LJ de OO --->'
!call lis_f(inp,epsi_lj)
!WRITE(*,2) 'Distance rOH (en A) et angle HOH (en degres) --->'
!call lis_f(inp,roh_a,theta_d)
!WRITE(*,2) 'Charge partielle de H --->'
!call lis_f(inp,charg_h)
!WRITE(*,3)
!WRITE(*,2) 'Coefficient K (en A-1) dans la decomposition Ewald (KL>>1) --->'
!call lis_f(inp,alp_a1)
!WRITE(*,2) 'On coupe les sommes en r et k de Ewald aux precisions sr,sk --->'
!call lis_f(inp,ssr,ssk)
!WRITE(*,2) 'Constante dielectrique exterieure --->'
!call lis_f(inp,cdiel_ext)
!WRITE(*,3)
!WRITE(*,2) 'NOUVEAU: Valeur maxi pour m,n dans le calcul des gmnl(r)  --->'
!call lis_f(inp,mnmax)
!if(mnmax==0) ialpmax=1
!if(mnmax==1) ialpmax=4
!if(mnmax==2) ialpmax=27
!if(mnmax==3) ialpmax=79
!if(mnmax==4) ialpmax=250
!if(mnmax==5) ialpmax=549
!if(mnmax==6) ialpmax=1226
!33 write(*,3)
!!        On calcule les parametres de la simulation avec L pris comme unite
!n=3*n_h2o
!avog=6.02214d-4
!!  C = N/ (NA*V) ==> V = N/ (C*NA)
!!  L = V^(1/3) 
!xl_a=(n_h2o/(conc_h2o*avog))**(1.d0/3.d0)   ! xl_a  - box length
!!  also: concentration is in mol/L = mol/dm^3
!!  so, answer is in dm, e.g.  L = ( N/ (C*NA))^(1/3) [dm]
!!  to convert from dm [10^-1 m ] to angstroms [10^-10 m] we multiply by  10^9 (or divide by 10^-9)
!!  this equivalently can be done inside the brackets (before taking ^(1/3)), which will give 10^-27 in the denominator
!!  combining 6.023e23 * 10^-27 we get 6.023e-4 (see the constant above)
!
!sigma=sigma_a/xl_a
!sigma2=sigma**2
!xlj=4.d0*epsi_lj                        ! donc v/kT=xlj*((sigma/r)**12-(sigma/r)**6)
!rcutoff_a=1.d20                         ! sert a rien mais sera dans les fichiers acc
!roh=roh_a/xl_a
!pi=4.d0*ATAN(1.d0)
!boltz=1.38065d-23
!elec=1.602176d-19
!epsi0=8.854187d-12   ! dialectrical permutivity of vacuum [ in SI ]
!dbjr_a=elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10    ! Bjerum length  en A
!dbjr=dbjr_a/xl_a
!xclb=charg_h**2*dbjr
!theta=pi/180.d0*theta_d
!theta2=theta/2.d0
!cthet2=COS(theta2); sthet2=SIN(theta2)
!rhh=2.d0*roh*sthet2
!!u_self=xclb*(-4.d0/roh+1.d0/rhh)
!!foh_self=-2.d0*xclb/roh**2
!!fhh_self=+xclb/rhh**2
!p=p_a/xl_a    ! p_a is read in section 4000 ( Parametres numeriques d''analyse )
!              ! it is spacing for g(r) 
!
!pasq=pasq_a1*xl_a ! also in 4000. It is spacing in q-space for diffusion spectra
!
!
!alp=alp_a1*xl_a   ! coefficient in the ewald decomposition (see above) 
!pi2=2.d0*pi 
!pialp2=pi**2/alp**2
!alp2pi=2.d0*alp/dSQRT(pi)
!rmax=ssr/alp
!rmax2=rmax**2
!kmax=NINT(alp*ssk/pi)
!ssk_vrai=kmax*pi/alp
!xkmax2=kmax**2
!ucut=xlj*((sigma/rmax)**12-(sigma/rmax)**6)     ! si on s'amuse a revenir a LJ tronque shifte a rmax
!PRINT*, 'Nombre de sites en tout : ',n
!PRINT*, 'Taille de la boite L (en A) : ',xl_a
!PRINT*, 'Diametre en unite L : ',sigma
!PRINT*, 'Distance OH en unite L : ',roh
!      WRITE(*,*) 'KL (doit etre >>1) : ',alp
!      WRITE(*,*) 'rmax/L, kmax, sk reel pour Ewald= ',rmax,kmax,ssk_vrai
!if(kmax>kmaxmax) PRINT*, 'ATTENTION: kmax dépasse la taille des tableaux!!!!'
!PRINT*, 'Nombre de projections : ',ialpmax
!!
!!             attention: si on a change la concentration..., L a change
!!             et donc les positions absolues venant des memes positions relatives
!!             ont ete implicitement renormalisées
!!             bien pour les centres, faux pour les distances intra!!!
!!             donc corriger  luc81p31
!!
!if(ilec23==17) then
!do i_h2o=1,n_h2o
!i0=(i_h2o-1)*3+1
!do i=i0+1,i0+2
!xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
!rcs=SQRT(xcs**2+ycs**2+zcs**2)
!xx(i)=xx(i0)+roh*xcs/rcs
!yy(i)=yy(i0)+roh*ycs/rcs
!zz(i)=zz(i0)+roh*zcs/rcs
!end do
!end do
!end if
!!
!      ilec1=17
!      goto 9999
!!
!!        position initiale des particules sur un reseau cfc pour les coeurs
!!        et orientations au hasard pour les sites
!!        d'abord H2O pointe en haut puis 3 rotations aleatoires suivant les 3 axes
!!
! 2000 if(n==0) goto 9999
!      xn=n_h2o
!      m=(xn/4.)**(1./3.)+0.999
!      xl1=1./m
!      write(*,3)
!      write(*,*) 'La cellule cubique est decoupee en ',m**3,' mailles'
!      write(*,*) 'de cote ',xl1
!      write(*,*) 'Les Oxygènes vont etre distribuees sur les ',4*m**3,        &
!               ' sites du reseau cfc ainsi forme'
!      i=1
!        do 100 ix=1,m
!        do 100 iy=1,m
!        do 100 iz=1,m
!        xx(i)=-0.5+(ix-1)*xl1
!        yy(i)=-0.5+(iy-1)*xl1
!        zz(i)=-0.5+(iz-1)*xl1
!        i=i+3
!  100   if(i>n) goto 101
!        do 104 ix=1,m
!        do 104 iy=1,m
!        do 104 iz=1,m
!        xx(i)=-0.5+xl1/2.+(ix-1)*xl1
!        yy(i)=-0.5+(iy-1)*xl1
!        zz(i)=-0.5+xl1/2.+(iz-1)*xl1
!        i=i+3
!  104   if(i>n) goto 101
!        do 105 ix=1,m
!        do 105 iy=1,m
!        do 105 iz=1,m
!        xx(i)=-0.5+xl1/2.+(ix-1)*xl1
!        yy(i)=-0.5+xl1/2.+(iy-1)*xl1
!        zz(i)=-0.5+(iz-1)*xl1
!        i=i+3
!  105   if(i>n) goto 101
!        do 106 ix=1,m
!        do 106 iy=1,m
!        do 106 iz=1,m
!        xx(i)=-0.5+(ix-1)*xl1
!        yy(i)=-0.5+xl1/2.+(iy-1)*xl1
!        zz(i)=-0.5+xl1/2.+(iz-1)*xl1
!        i=i+3
!  106   if(i>n) goto 101
!  101 continue
!!       orientations au hasard
!PRINT*, 'Et orientations au hasard de la molécule H2O'
!!        initialisation du generateur de nombres aleatoires
!!      x=alea(0)
!do i=1,n_h2o
!vx1=sthet2                           ! d'abord pointe en haut
!vy1=0.
!vz1=cthet2
!vx2=-sthet2
!vy2=0.
!vz2=cthet2
!!angle=pi2*alea(1)               ! rotation aleatoire suivant x
!angle = pi2 
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,1)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,1)
!!angle=pi2*alea(1)               ! rotation aleatoire suivant y
!!angle = pi2
!angle = pi2
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,2)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,2)
!!angle=pi2*alea(1)               ! rotation aleatoire suivant z
!angle = pi2
!ca=COS(angle)
!sa=SIN(angle)
!call rot_vect(vx1,vy1,vz1,vx1,vy1,vz1,ca,sa,3)
!call rot_vect(vx2,vy2,vz2,vx2,vy2,vz2,ca,sa,3)
!i0=3*(i-1)+1
!xx(i0+1)=xx(i0)+roh*vx1
!yy(i0+1)=yy(i0)+roh*vy1
!zz(i0+1)=zz(i0)+roh*vz1
!xx(i0+2)=xx(i0)+roh*vx2
!yy(i0+2)=yy(i0)+roh*vy2
!zz(i0+2)=zz(i0)+roh*vz2
!end do
!!    attention, certains sites peuvent etre a l'exterieur de la boite
!!    maintenant, je m'en fiche! luc80p179
!      ilec23=17
!      goto 9999
!!
!!        lecture en fichier des positions
!!
! 3000 continue
!      write(*,2) 'nom du fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=3000,file=nomfic,STATUS='old')
!      read(7,*,end=3000) n
!      write(*,*) 'n= ',n
!      read(7,*) (xx(i),yy(i),zz(i),i=1,n)
!      close(7)
!!              il peut arriver que cette config (ancienne) n'ait pas autorise les sites H
!!              a sortir legerement de la boite, donc les ait deplaces
!!              de l'autre cote
!!              maintenant, on force les sites a rester autour de leur coeur O
!!              quitte a sortir de la boite
!do i0=0,n-1,3
!do j=2,3
!xr=xx(i0+j)-xx(i0+1)
!yyr=yy(i0+j)-yy(i0+1)
!zr=zz(i0+j)-zz(i0+1)
!xr=xr-NINT(xr)
!yyr=yyr-NINT(yyr)
!zr=zr-NINT(zr)
!xx(i0+j)=xx(i0+1)+xr
!yy(i0+j)=yy(i0+1)+yyr
!zz(i0+j)=zz(i0+1)+zr
!end do
!end do
!!
!n_h2o=n/3
!conc_h2o=n_h2o/xl_a**3/avog
!PRINT*, 'Nombre total de sites : ',n
!PRINT*, 'Nombre de H2O : ',n_h2o
!      write(*,*) 'Concentration de H2O (en M) : ',conc_h2o
!      end if
!!
!!             attention: si on a change la concentration..., L a change
!!             et donc les positions absolues venant des memes positions relatives
!!             ont ete implicitement renormalisées
!!             bien pour les centres, faux pour les distances intra!!!
!!             donc corriger  luc81p31
!!
!do i_h2o=1,n_h2o
!i0=(i_h2o-1)*3+1
!do i=i0+1,i0+2
!xcs=xx(i)-xx(i0); ycs=yy(i)-yy(i0); zcs=zz(i)-zz(i0)
!rcs=SQRT(xcs**2+ycs**2+zcs**2)
!xx(i)=xx(i0)+roh*xcs/rcs
!yy(i)=yy(i0)+roh*ycs/rcs
!zz(i)=zz(i0)+roh*zcs/rcs
!end do
!end do
!!
!      ilec23=17
!      goto 9999
!!
!!        parametres d'analyse pour la probabilite g(r)
!!
! 4000 if(ilec1.ne.17) goto 9999
!      write(*,2) 'pour g(r), nombre d''intervalles --->'
!      call lis_f(inp,num)
!      WRITE(*,2) 'et pas (en A) --->'
!      call lis_f(inp,p_a)
!      p=p_a/xl_a
!      write(*,*) 'Pas en unite L : ',p
!      WRITE(*,2) 'pour g(r,omega), pas en cos(beta) et pas en phi (en rd) --->'      ! g(r,omega)  luc85p176
!      call lis_f(inp,dcosbet,dphi)
!      WRITE(*,2) 'et coefficients correcteurs pour A=C=E=G et B=D=F=H  --->'
!      call lis_f(inp,coeff_a,coeff_b)
!      write(*,1) 'Pour le potentiel chimique, on fait des insertions de particules tests M fois au hasard'
!      write(*,2) 'Valeur de M --->'
!      call lis_f(inp,mmu)
!      WRITE(*,*) 'Pour le calcul des spectres de diffusion :'
!      WRITE(*,2) 'Nombre de points en q --->'
!      CALL lis_f(inp,nptq)
!      WRITE(*,2) 'Pas en q (A-1) --->'
!      call lis_f(inp,pasq_a1)
!      pasq=pasq_a1*xl_a
!      WRITE(*,2) 'Nombre d''orientations de q a chaque config --->'
!      call lis_f(inp,mq)
!  44  WRITE(*,2) 'On efface tout? (oui=1) --->'
!      READ(inp,*,ERR=44) ieff
!      if(ieff==1) then
!      iacc(1:ialpmax,1:num)=0.
!      iacc1(1:num)=0.; iacc2(1:num)=0.; iacc3(1:num)=0.
!      nittot=0
!      naccep=0
!      ncumul=0
!      uutotm=0.
!      duutot=0.
!      evtest=0.
!      devtes=0.
!      vir=0.
!      vir1=0.
!      vir2=0.
!      vir02=0.
!      vir12=0.
!      vir22=0.
!      hypervir=0.
!      hypervir2=0.
!      xmtotx=0.; xmtoty=0.; xmtotz=0.; xmtot2=0.; xmtot22=0.
!      do kq=1,nptq
!      sq1(kq)=0.
!      sq2(kq)=0.
!      sq3(kq)=0.
!      end do
!      iacc_o(:,1:num)=0.
!      end if
!      ilec4=17
!      goto 9999
!!
!!        lecture en fichier des resultats
!!
! 5000 continue
!      write(*,2) 'nom du fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=5000,status='old',file=nomfic)
!      read(7,*,err=5000) n_h2o
!      write(*,*) 'Nombre de H2O : ',n_h2o
!      read(7,*) temp
!      write(*,*) 'Temperature (en K) : ',temp
!      read(7,*) conc_h2o
!      write(*,*) 'Concentration (M) : ',conc_h2o
!      READ(7,*) sigma_a
!      WRITE(*,*) 'Diametre LJ OO (A) : ',sigma_a
!      READ(7,*) epsi_lj
!      WRITE(*,*) 'Epsilon/kT LJ OO : ',epsi_lj
!      READ(7,*) rcutoff_a
!      READ(7,*) roh_a,theta_d
!      WRITE(*,*) 'Distance OH (en A) et angle HOH (en degres) : ',roh_a,theta_d
!      READ(7,*) charg_h
!      WRITE(*,*) 'Charge partielle de H : ',charg_h
!      READ(7,*) alp_a1
!      PRINT*, 'Coefficient K Ewald (en A-1) : ',alp_a1
!      READ(7,*) ssr,ssk
!      PRINT*, 'Precisions r et k Ewald : ',ssr,ssk
!      READ(7,*) cdiel_ext
!      PRINT*, 'Constante dielectrique externe : ',cdiel_ext
!      READ(7,*) mnmax
!      PRINT*, 'Valeur max de m et n dans les projections : ',mnmax
!if(mnmax==0) ialpmax=1
!if(mnmax==1) ialpmax=4
!if(mnmax==2) ialpmax=27
!if(mnmax==3) ialpmax=79
!if(mnmax==4) ialpmax=250
!if(mnmax==5) ialpmax=549
!if(mnmax==6) ialpmax=1226
!PRINT*, 'Nombre de projections : ',ialpmax
!      read(7,*) num,p_a
!      write(*,*) 'Nombre d''intervalles et pas (en A) = ',num,p_a
!      read(7,*) nittot,naccep,ncumul
!      write(*,*) 'nombre d''iterations deja effectuees : ',nittot
!      write(*,*) 'nombre d''accumulations deja effectuees : ',ncumul
!      read(7,*) uutotm,duutot
!      read(7,*) vir,vir1,vir2,vir02,vir12,vir22
!      READ(7,*) hypervir,hypervir2
!      READ(7,*) xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
!      read(7,*) evtest,devtes,mmu
!      WRITE(*,*) 'Nombre d''insertions par configuration : ',mmu
!      do i=1,num
!      READ(7,*,END=5010) r,iacc1(i),iacc2(i),iacc3(i)    ! r ne sert pas
!      end do
!      do i=1,num
!      READ(7,*,END=5010) r,iacc(1:ialpmax,i)
!      end do
!      READ(7,*,END=5010) nptq,pasq_a1,mq
!      WRITE(*,*) 'En q: ',nptq,pasq_a1,mq
!      do kq=1,nptq
!      READ(7,*,END=5010) q,sq1(kq),sq2(kq),sq3(kq)
!      end do
!!         normalisation des Sij(q)
!      den=ncumul
!      n_o=n/3
!      n_h=n-n_o
!      do kq=1,nptq
!      sk1(kq)=sq1(kq)/(den*n_o)
!      sk2(kq)=sq2(kq)/(den*SQRT(n_o*n_h*1.))
!      sk3(kq)=sq3(kq)/(den*n_h)
!      end do
!      READ(7,*,END=5010) dcosbet,dphi                     ! g(r,omega) luc85p176
!      WRITE(*,*) 'pas en cosbeta et phi pour g(r,omega) : ',dcosbet,dphi
!      READ(7,*,END=5010) coeff_a,coeff_b
!      WRITE(*,*) 'coefficients correcteurs en A et en B : ',coeff_a,coeff_b
!      do i=1,num
!      READ(7,*,END=5010) r,iacc_o(1:n_omega,i)    ! r ne sert pas
!      end do
!!PRINT*, 'je divise par 2 les cas G et H'           ! juste pour rattraper les fichiers 48 et 49
!!iacc_o(7:8,1:num)=0.5d0*iacc_o(7:8,1:num)
!      GOTO 5011
! 5010 PRINT*, 'ATTENTION!!!!!! fin de fichier detectee!'
! 5011 continue
!      close(7)
!      ilec4=17
!      goto 33
!        end if
!      GOTO 9999
!!
!!      pour estimer l'efficacite statistique s, on a le choix entre differentes variables
!!
! 9600 WRITE(*,*) 'On a le choix entre les variables suivantes : '
!      WRITE(*,*) '1--Energie  2--viriel  '
!      WRITE(*,*) '6--exp(-vtest)  7--hyperviriel 8--M**2  9,10,11--Mx,y,z'
!      WRITE(*,*) '10000*alpha+i pour g_alpha(ri)'                 ! et non 1000*alpha... luc82p166
!      WRITE(*,*) 'Le choix actuel est : ',ichoixaaa
!      WRITE(*,2) 'Nouveau choix souhaite --->'
!      read(inp,*,ERR=9999) nent
!      ichoixaaa=nent
!      GOTO 9999
!!
!!      on peut sauvagarder en fichier les valeurs instantanees de la variable
!!
! 9700 WRITE(*,2) 'Nom du fichier de sauvegarde ("qq" pour ignorer) --->'
!      READ(inp,'(A)') nomaaa
!      GOTO 9999
!!
!!       on regroupe des intervalles pour lisser
!!
! 1100 if(ilec4.ne.17) goto 9999
!      write(*,*) 'Actuellement, le nombre d''intervalles vaut ',num
!      WRITE(*,*) 'et le pas (en A) ',p_a
!      write(*,2) 'On va les regrouper par paquets de combien ? --->'
!      call lis_f(inp,ipaq)
!      num1=num/ipaq
!      do k=1,num1
!      iacc(:,k)=iacc(:,(k-1)*ipaq+1)
!      iacc1(k)=iacc1((k-1)*ipaq+1)
!      iacc2(k)=iacc2((k-1)*ipaq+1)
!      iacc3(k)=iacc3((k-1)*ipaq+1)
!      iacc_o(:,k)=iacc_o(:,(k-1)*ipaq+1)
!       do l=2,ipaq
!       iacc(:,k)=iacc(:,k)+iacc(:,(k-1)*ipaq+l)
!       iacc1(k)=iacc1(k)+iacc1((k-1)*ipaq+l)
!       iacc2(k)=iacc2(k)+iacc2((k-1)*ipaq+l)
!       iacc3(k)=iacc3(k)+iacc3((k-1)*ipaq+l)
!       iacc_o(:,k)=iacc_o(:,k)+iacc_o(:,(k-1)*ipaq+l)
!       END do
!      END do
!      num=num1
!      p_a=p_a*float(ipaq)
!      p=p_a/xl_a
!      WRITE(*,*) 'Nouveaux nombre de points et pas: ',num,p_a
!      goto 9999
!!
!!        on cumule les donnees actuelles avec celles d'un autre fichier luc80p127
!!
! 1200 if(ilec4.ne.17) goto 9999
!      WRITE(*,1) 'On cumule les donnees actuelles avec celles lues dans un autre fichier'
!      WRITE(*,1) 'donnees = alpha * actuelles + beta * lues'
!      write(*,2) 'nom du nouveau fichier (pas de lecture si "qq")  --->'
!      read(inp,'(A)') nomfic
!        if(nomfic.ne.'qq'.and.nomfic.ne.'QQ') then
!      open(7,err=9999,status='old',file=nomfic)
!      WRITE(*,2) 'Parametres de melange alpha et beta --->'
!      call lis_f(inp,alpha,beta)
!      PRINT*, 'ATTENTION: VERifIER que les 2 fichiers correspondent au meme systeme!'
!      read(7,*,err=5000) n_h2o
!      write(*,*) 'Nombre de H2O : ',n_h2o
!      read(7,*) temp
!      write(*,*) 'Temperature (en K) : ',temp
!      read(7,*) conc_h2o
!      write(*,*) 'Concentration de H2O (en M) : ',conc_h2o
!      READ(7,*) sigma_a
!      WRITE(*,*) 'Diametre LJ OO (A) : ',sigma_a
!      READ(7,*) epsi_lj
!      WRITE(*,*) 'Epsilon/kT LJ OO : ',epsi_lj
!      READ(7,*) rcutoff_a
!      READ(7,*) roh_a,theta_d
!      WRITE(*,*) 'Distance OH (en A) et angle HOH (en degres) : ',roh_a,theta_d
!      READ(7,*) charg_h
!      WRITE(*,*) 'Charge partielle de H : ',charg_h
!      READ(7,*) alp_a1
!      PRINT*, 'Coefficient K Ewald (en A-1) : ',alp_a1
!      READ(7,*) ssr,ssk
!      PRINT*, 'Precisions r et k Ewald : ',ssr,ssk
!      READ(7,*) cdiel_ext
!      PRINT*, 'Constante dielectrique externe : ',cdiel_ext
!      READ(7,*) mnmax
!      PRINT*, 'Valeur max de m et n dans les projections : ',mnmax
!      read(7,*) num,p_a
!      write(*,*) 'Nombre d''intervalles et pas (en A) = ',num,p_a
!      read(7,*) nittot_1,naccep_1,ncumul_1
!      write(*,*) 'nombre d''iterations deja effectuees : ',nittot_1
!      write(*,*) 'nombre d''accumulations deja effectuees : ',ncumul_1
!      nittot=NINT(alpha*nittot+beta*nittot_1)
!      naccep=NINT(alpha*naccep+beta*naccep_1)
!      x=alpha*ncumul; x_1=beta*ncumul_1
!      ncumul=nint(alpha*ncumul+beta*ncumul_1)
!      x=x/ncumul; x_1=x_1/ncumul
!      read(7,*) uutotm_1,duutot_1
!      read(7,*) vir_1,vir1_1,vir2_1,vir02_1,vir12_1,vir22_1
!      READ(7,*) hypervir_1,hypervir2_1
!      READ(7,*,end=1210) xmtotx_1,xmtoty_1,xmtotz_1,xmtot2_1,xmtot22_1
!      uutotm=alpha*uutotm+beta*uutotm_1
!      duutot=alpha*duutot+beta*duutot_1
!      vir=alpha*vir+beta*vir_1
!      vir1=alpha*vir1+beta*vir1_1
!      vir2=alpha*vir2+beta*vir2_1
!      vir02=alpha*vir02+beta*vir02_1
!      vir12=alpha*vir12+beta*vir12_1
!      vir22=alpha*vir22+beta*vir22_1
!      hypervir=alpha*hypervir+beta*hypervir_1
!      hypervir2=alpha*hypervir2+beta*hypervir2_1
!      xmtotx=alpha*xmtotx+beta*xmtotx_1
!      xmtoty=alpha*xmtoty+beta*xmtoty_1
!      xmtotz=alpha*xmtotz+beta*xmtotz_1
!      xmtot2=alpha*xmtot2+beta*xmtot2_1
!      xmtot22=alpha*xmtot22+beta*xmtot22_1
!      read(7,*) evtest_1,devtes_1,mmu_1
!      WRITE(*,*) 'Nombre d''insertions par configuration : ',mmu_1
!!         attention: on peut envisager que mmu et mmu_1 different luc82p163
!!         rappel: evtest est deja normalise par mmu
!      if(mmu>0.or.mmu_1>0) then
!      xmu_moyen=x*mmu+x_1*mmu_1
!      evtest=(alpha*mmu*evtest+beta*mmu_1*evtest_1)/xmu_moyen
!      devtes=(alpha*mmu*devtes+beta*mmu_1*devtes_1)/xmu_moyen
!      end if
!!         puis on prend pour mmu la nouvelle valeur qui vient d'etre lue
!      mmu=mmu_1
!!       on se sert de gr comme tableau auxi
!      do i=1,num
!      read(7,*) r,gr(1:3,i)
!      end do
!      iacc1(1:num)=NINT(alpha*iacc1(1:num)+beta*gr(1,1:num))
!      iacc2(1:num)=NINT(alpha*iacc2(1:num)+beta*gr(2,1:num))
!      iacc3(1:num)=NINT(alpha*iacc3(1:num)+beta*gr(3,1:num))
!      do i=1,num
!      read(7,*) r,gr(1:ialpmax,i)
!      END do
!      iacc(1:ialpmax,1:num)=alpha*iacc(1:ialpmax,1:num)+beta*gr(1:ialpmax,1:num)
!      READ(7,*,END=1210) nptq,pasq_a1,mq
!      WRITE(*,*) 'En q: ',nptq,pasq_a1,mq
!      do kq=1,nptq
!      READ(7,*,END=1210) q,sk1(kq),sk2(kq),sk3(kq)
!      end do
!      sq1(1:nptq)=alpha*sq1(1:nptq)+beta*sk1(1:nptq)
!      sq2(1:nptq)=alpha*sq2(1:nptq)+beta*sk2(1:nptq)
!      sq3(1:nptq)=alpha*sq3(1:nptq)+beta*sk3(1:nptq)
!!         normalisation des Sij(q)
!      den=ncumul
!!      n_o=n/3
!!      n_h=n-n_o
!!      do kq=1,nptq
!!      sk1(kq)=sq1(kq)/(den*n_o)
!!      sk2(kq)=sq2(kq)/(den*SQRT(n_o*n_h*1.))
!!      sk3(kq)=sq3(kq)/(den*n_h)
!      end do
!      READ(7,*,END=1210) dcosbet,dphi
!      WRITE(*,*) 'pas en cosbeta et phi pour g(r,omega) : ',dcosbet,dphi
!!      READ(7,*,END=5010) coeff_a,coeff_b
!      READ(7,*) coeff_a,coeff_b
!      WRITE(*,*) 'coefficients correcteurs en A et en B : ',coeff_a,coeff_b
!      do i=1,num
!      read(7,*,END=1210) r,gr(1:n_omega,i)
!      END do
!      iacc_o(1:n_omega,1:num)=alpha*iacc_o(1:n_omega,1:num)+beta*gr(1:n_omega,1:num)
!      GOTO 1211
! 1210 PRINT*, 'ATTENTION: fin de fichier detectee!'
! 1211 continue
!      close(7)
!      GOTO 33
!        end if
!      GOTO 9999
!!
!!      on demarre les cycles si on a bien tout defini
!!      juste avant, on initialise les energies et les nbrs aleatoires
!!
! 8000 if(ilec1.ne.17.or.ilec23/=17.or.ilec4/=17) goto 9999
!
!        on initialise les energies et les forces en r SANS intra luc84p187
!        on calcule maintenant LJ total par Ewald aussi!
!i
    sigma2 = sigma2_in

     charge_h = 0.42380d0
     xxclb = xclb * charge_h**2

      alp = alpha
      alp2pi=2.d0*alp/dSQRT(pi)
      n = n_h2o * 3

      uu(1:n)=0.                         ! repere LJr+EWr
      fx(1:n)=0.; fy(1:n)=0.; fz(1:n)=0.
      uu_lj12=0.; uu_lj6=0.
      uu_ew=0.
      summux=0.;summuy=0.; summuz=0.
      
      !        pour test, je reviens a LJ tronque shifte a rcut=rmax  luc85p123   NON, plus maintenant
      !xlj0=xlj; xlj=0.
      
!       EWALD SUM
!
!       U = U_real + U_fourier - U_self
!
!       where 
!
!       U_real = SUM_i SUM_j<>i  q_i q_j erfc( alpha r_ij ) / r_ij 
!              (the first term. Other images are ommited, and as r_ij is taken the nearliest image)
!
!
!       U_fourier = 1 / 2V SUM_m<>0   4pi/k_m^2 SUM_i SUM_j q_i q_j exp(i k_m r_ij) exp( - k_m^2/4 alpha^2)
        
!       U_self = sqrt(alpha^2/pi) SUM_i q_i^2    


      ! calculation of potential (Real part of Ewald sum )

      do i=1,n
         iio=3*((i-1)/3)+1   ! equivalent to mod(i,3)==1 : 
                            ! i == 3k+1  --> (i-1)/3 = k,   iio = 3*k + 1 == i
                            ! i == 3k+2 --> (i-1)/3 = k     iio = 3*k + 1 <> i 
         chargi=+1.d0       ! i == 3k   --> (i-1)/3 = k-1   iio = 3(k-1) + 1 = 3k-2 <> i
          if(i==iio) chargi=-2.d0 ! 

         summux=summux+chargi*xx(i)                      ! rq: j'ignore ici que H peut sortir de la boite
         summuy=summuy+chargi*yy(i)
         summuz=summuz+chargi*zz(i)

         do j=i+1,n

            jo=3*((j-1)/3)+1
            chargj=+1.d0

            if(j==jo) chargj=-2.d0

            xij=xx(i)-xx(j)
            yij=yy(i)-yy(j)
            zij=zz(i)-zz(j)
!            xij=xij-anint(xij)    ! anint =  round  
!            yij=yij-anint(yij)    ! Coordinates are in "box" units (1=length of the box)
!            zij=zij-anint(zij)    ! i.e. x - round(x) = nearliest image
         if(xij.gt. 0.5) xij=xij-1.   ! nearliest neighbour. Remember: distance units are BoxLength
         if(xij.lt.-0.5) xij=xij+1.   ! i.e. Box is [-0.5:0.5]^3
 
         if(yij.gt. 0.5) yij=yij-1.
         if(yij.lt.-0.5) yij=yij+1.         ! prend moins de temps que anint!
 
         if(zij.gt. 0.5) zij=zij-1.
         if(zij.lt.-0.5) zij=zij+1.
                                  ! only one image is used ( efrc sum has only 1 term ) 

            r2=xij**2+yij**2+zij**2   ! distanse to the nearliest image

!            write(*,*) 'xij=',xij,'yij=',yij,'zij=',zij

            if(r2<rmax2) then                         ! rmax pour Ewald r   rmax = ssr / alp 
               r=SQRT(r2)
               alpr=alp*r
               alpr2=alpr**2
               ee2=EXP(-alpr2)
               er=erfc_Luc_OLD(alpr,alpr2,ee2)            ! calcule erfc_Luc(x) en s'aidant de x**2 et exp(-x**2)
 !               write(*,*) 'alpr=',alpr,'akpr2=',alpr2,'ee2=',ee2,'er=',er,'potew=',potew

            
               if(iio==jo) er=-erf_Luc_OLD(alpr,alpr2,ee2)    ! retrancher le self eventuellement donc erfc_Luc-1=-erf  luc84p134
                                                         ! iio == jo <=> for the atoms of the same molecule 
               ! erfc - 1 = - erf
               ! Sergiievskyi 31 May 2014: 
               ! > why not to remove the intermolecular interactiions at all?
               ! > after all: they are always constant for the rigid molecules, thus do not matter

!               write(*,*) 'er_after=',er



               potew=xxclb*chargi*chargj*er/r

  !             write(*,*) 'xclb',xclb,'qi',charge_h,'xxclb*chargi*chargj',xxclb*chargi*chargj,'er',er,'r',r
               uu_ew=uu_ew+potew
               pot=potew


               ! ffr = -1/r dU/dr
               ! for coulomb:  U_ij = q_i q_j erfc(alpha r)/ r
               ! dU/dr = q_i q_j ( - erfc(alpha r) / r^2 + derfc(alpha r) / dr * 1/r )
               ! erfc(alpha r ) = 1 - erf = 1 - 2 / sqrt(pi) int_0^{alpha r} e^-t^2 dt                
               ! -> derfc(alpha r) /dr =  alpha derfc(alpha r)/d(alpha r) = - 2 alpha/sqrt(pi) e^{-alpha^2 r^2} 
               !  
               ! ffr = q_i q_j ( erfc(alpha r) /r^3 +  {2 alpha /r^2 }exp(-alpha^2 r^2) ) 
               ffr=xxclb*chargi*chargj*(er/r+alp2pi*ee2)/r2     


!               write(*,*) 'xclbqiqj:',xxclb*chargi*chargj,'er:',er,'r:',r,'alp2pi:',alp2pi,'ee2:',ee2,'r2:',r2    
!               write(*,*) 'ffr_c:',ffr

               ! O-O ineractions : LJ ewald
               ! c.f. Karasawa, Goddard, JPC 1989, 93, 7320-7327
               !      Williams,   Acta Cryst A,  1971, 27, 452
               !    
               ! General case: potential 1/r^w
               ! 
               !     S =  S_real + S_fourier + S_{k=0} - S_self
               !
               ! where
               ! 
               !  S_real =  K  SUM_{sjn,r_sjn>0}  A_sj Gamma(w/2,a^2) r^-w 
               !  S_fourier = K  pi^{w-1.5}/V  SUM_{m<>0}  F_2(k_m) h_m^3 Gamma(3/2 - w/2,b^2) 
               !  S_{k=0}  =  K pi^1.5 alpha^{w - 3} / V * 2/(w-3) SUM_{sj} A_sj
               !  S_self = K 2 alpha^{w} / w * SUM_j A_jj 
               !
               !  K = 1 / (2 Gamma(w/2) )
               !  a^2  =  alpha^2 r_{sj;n}^2 
               !  b^2  =  pi^2 h_m^2 / alpha
               !  h_m  = k_m / 2pi = m / L
               !  r_{sj;n} = | r_s - r_j + n*L |
               !  Gamma(w/2; a^2) = int_{a^2}^inf t^{w/2-1} exp(-t) dt
               !  Gamma(w/2) = Gamma(w/2;0)
               !  F_2(k_m) = SUM_sj A_sj exp(i k_m (r_s - r_j ) )  
               !     for A_sj = q_s*q_j  F_2(k_m) = F(k_m) F(-k_m) = |F(k_m)|^2 
               !     where F(k_m) = SUM_j q_j exp(i k_m r_j)
               !
               ! Particular cases: w=1 --> usual Ewald sum
               !
               ! w=2k+2 (k integer): for that case
               !
               !  Gamma(w/2; a^2) = k! exp(-a^2)  SUM_{p=0}^k a^{2p}/p!
               !  Gamma(3/2 - w/2; b^2)  = (-2)^k / (2k-1)!! [  sqrt(pi) erfc(b)  - b^-1 exp(-b^2) ]
               !                            + exp(-b^2) / (2k-1)!! SUM_{p=1}^{k-1}  (-2)^(k-p) (2p-1)!! b^{-2p-1}
               !
               ! which gives:
               ! w=6:   
               !     S_real = 1/2 SUM_{sjn} A_sj ( 1 + a^2 + 0.5 a^4) exp(-a^2) r^-6
               !     S_fourier = pi^4.5/3V SUM_m F_2(k_m) h_m^3 [ sqrt(pi) erfc(b) + (1/2b^3 - 1/b) exp(-b^2 ] 
               !     S_{k=0} = 1/6V * pi^1.5 alpha^3 SUM_sj A_sj 
               !     S_self = alpha^6/12 SUM_j A_jj
               !
               ! w=12: 
               !     S_real = 1/2 SUM_{r_nsj>0} A_sj (a^10/120 + a^8/24 + a^6/6 + a^4/2 + a^2 + 1 ) exp(-a^2) r^-12
               !     S_fourier =  1/(945*120) pi^10.5/V *
               !                *  SUM_m<>0 F_2(k_m) h_m^9 (-16 sqrt(pi) erfc(b) + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) )
               !     S_{k=0} = pi^1.5 alpha^9 / 1080V * SUM_sj A_sj
               !     S_self =  alpha^12 / 1440 SUM_j A_jj 
  
 !              write(*,*) 'sigma2=',sigma2
!               write(*,*) 'i',i,'iio',iio,'j',j,'jo',jo 

             if(i==iio.and.j==jo) then            ! OO: LJ total  luc84p147

                  vv=(sigma2/r2)**3                ! vv = (sigma/r)^6
 
!                  write(*,*) 'sigma2=',sigma2,'r2=',r2
                   c6=1.d0+alpr2*(1.d0+alpr2/2.d0)  !   1 + a^2 + 0.5 a^4
                  potlj6=-xlj*vv*ee2*c6            !  4*epsilon * (1+a^2 + 0.5a^4) exp(-a^2) (sigma/r)^6
!                  write(*,*) 'xlj:',xlj,'vv:',vv,'ee2',ee2,'c6:',c6
                  alpr6=alpr2**3                   !  a^6 
                  c12=c6+alpr6*(1.d0+alpr2*(1.d0+alpr2/5.d0)/4.d0)/6.d0  ! C12 = 1 + a^2 + 0.5 a^4 + a^6/6 + a^8/24 + a^10/120
                  potlj12=xlj*vv**2*ee2*c12                              ! 4 epsilon * C12  * exp(-a^2) * (sigma/r)^12
                  uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
                  pot=pot+potlj12+potlj6
                  ffr=ffr+xlj*6.d0*vv*ee2*(2.d0*vv*(c12+alpr6**2/720.d0)-(c6+alpr6/6.d0))/r2

  !                write(*,*) 'xlj',xlj,'vv',vv,'ee2',ee2,'c12',c12,'alpr6',alpr6,'c6',c6
 !                 write(*,*) 'ffr_tot:',ffr
!
 !                 write(*,*) 'vv=',vv,'c6=',c6,'potlj6=',potlj6,'potlj12=',potlj12,'c12=',c12
!                  write(*,*) 'potlj6=',potlj6,'potlj12=',potlj12,'potlj=',potlj12 + potlj6

                      !
                      ! 1/r^2(   4 eps * 12 sigma^12/r^12 e^-a^2 ( C12 + (alpha r)^12 / 6! )
                      !        - 4 eps * 6 * sigma^6/r^6  e^-a^2 ( C6 + (alpha r)^6 / 3! )  )
                      !
                      
                      
      !     pour test, LJ tronque shifte (xlj a ete mis a 0)         NON
      !            potlj6=-xlj0*vv
      !            potlj12=xlj0*vv**2
      !            uu_lj12=uu_lj12+potlj12; uu_lj6=uu_lj6+potlj6
      !            pot=pot+potlj12+potlj6
      !            ffr=ffr+xlj0*6.d0*vv*(2.d0*vv-1.d0)/r2
                end if
 
!         write(*,*) 'LUC: pot:',pot,'ew:',potew,'lj6:',potlj6,'lj12:',potlj12

                ! not same molecules
                if(iio/=jo) then 
                   uu(i)=uu(i)+pot                  ! ne contient plus le self
                   uu(j)=uu(j)+pot                  ! mais contient l'intra!
                   fxij=xij*ffr                     ! NON, plus maintenant!

                   !write(*,*) 'ffr:',ffr,'fxij:',fxij,'x(i):',xx(i),'x(j):',xx(j),'xij',xij
                   fx(i)=fx(i)+fxij
                   fx(j)=fx(j)-fxij                  
                   fyij=yij*ffr
                   fy(i)=fy(i)+fyij
                   fy(j)=fy(j)-fyij
                   fzij=zij*ffr
                   fz(i)=fz(i)+fzij
                   fz(j)=fz(j)-fzij
                end if ! iio /= jo

             end if ! r2 < rmax2
         end do ! j
      end do ! i


      uu_ew_r=uu_ew
      uu_lj_r=uu_lj12+uu_lj6
      uu_lj6_r = uu_lj6
      uu_lj12_r = uu_lj12
      
!       PRINT*, 'LUC_REAL_SUM: uu_lj6:',uu_lj6,'uu_lj12:',uu_lj12,'uu_lj:',uu_lj_r,'ew:',uu_ew_r
  

!      write(*,*) 'uu_Luc:'
!      call write_real_array(0,uu,n)

!      write(*,*) 'forces Luc:'
!      call write_xyz_array(0,fx,fy,fz,n)

!      return
      !
      !        je prepare la somme en k pour Ewald  luc70p21
      !
      !
      !  U_fourier = 1/2V SUM_m<>0 4pi/k_m SUM_js q_jq_s exp(ik_m r_js) exp(-k_m^2/4alpha^2) 
      !
      !  Using that rho(k) = SUM q_j exp(i k r_j)
      !
      !  U_fourier = 1/2V SUM_k<>0 4pi/k^2 |rho(k)|^2 exp(-k^2/(4 alpha^2 )) 
      !
      pialp=pi/alp
      !pialp2 = pialp**2
      pialp2=pi**2/alp**2

      kk=0
      pi2=two_pi
        
 
      xkmax2 = kmax**2
 
!      write(*,*) 'KMAX = ',kmax
 
      do kx=0,kmax
     
 
        xkx2=kx**2
        kymax=SQRT(xkmax2-xkx2)
 
 !       write(*,*) 'KX:',kx,'Kymax=',kymax
      
        do ky=0,kymax
 
            xky2=ky**2
            kzmax=SQRT(xkmax2-xkx2-xky2)
      
  !          write(*,*) 'KY:',ky,'kzmax:',kzmax
      
            do kz=0,kzmax
                 xk2=xkx2+xky2+kz**2

               

                 if(xk2<0.1) cycle            ! k different de 0
                 kk=kk+1

   !              write(*,*) 'KZ:',kz,'xk2',xk2,'kk',kk
         
                 if(kk.gt.nkmax) then
                      WRITE(*,*) 'attention: on depasse la dimension en k',nkmax
                      return
!                      goto 9999
                 end if
     
                 ! pialp2 = pi^2 / alpha^2
                 pialpk2=pialp2*xk2    ! pi^2 / alpha^2 * k^2 
!                k_m = 2*pi*m/L 
!                exp( - k_m^2 / 4 alpha^2 ) = exp( - 4 pi^2 m^2 / L^2 / 4 alpha^2) = 
!                  = exp( - pi^2 * m^2 / L^2 /alpha^2 ) 
!                now, in our notation m <-> k, L=1
!                  exp ( -k_m^2 / 4 alpha^2 ) = exp( - pi^2 *k^2 / alpha^2 ) 
!                                     
                 eek=EXP(-pialpk2)     ! exp(-pi^2/alpha^2 * k^2) 
                 bet(kk)=eek/xk2/pi2   ! 1/2piL pour charge-charge
                                       ! bet(kk) = 1/2pi/k^2 * exp(-pi^2 / alpha^2 * k^2 ) 

                                     ! why 1/2pi: 1/2 is the sum prefactor
                                     ! 4pi/k_m^2 = 4pi/ (4pi^2 L^2/m^2) = 1/ (pi * m^2) = 1 / (pi* k^2) (in our notation)
                 xk=SQRT(xk2)          
                 pialpk=pialp*xk       ! pi/alpha * k


           ! calculation of beta6 and beta12 for LJ sums in Fourier space
           ! 
           !      erfk=erfc_Luc(pialpk,pialpk2,eek)
           !      bet6(kk)=pi45/3.d0*xk2*xk*(pi12*erfk+(0.5d0/pialpk2-1.d0)/pialpk*eek)
           !      bet12(kk)=pi105/1080.d0*xk2**4*xk*(-16.d0*pi12/105.d0*erfk+           &
           !         (((((1.d0/pialpk2-2.d0/7.d0)/pialpk2+4.d0/35.d0)/pialpk2-8.d0/105.d0)/pialpk2+16.d0/105.d0))/pialpk*eek)
           !      nouveau: on peut utiliser un DL ou/et une integration numerique partielle  luc84p150
           !
           ! ( for details:  see definition of erfc_Luc_bet6_bet12)

                 call erfc_Luc_bet6_bet12_OLD(pialpk,pialpk2,eek,erfk,x6,x12)
                 bet6(kk)=alp**3*x6   ! alpha^3 *   pi^1.5/3V * b^3 (sqrt(pi) erfc(b) + (1/2b^3 -1b)exp(-b^2) ) 
                 bet12(kk)=alp**9*x12 ! alpha^9 *  1/(945*120) *pi^1.5 *  b^9 * 
                                      !            ( -16 sqrt(pi) erfc(b) 
                                      !              + [16/b - 8/b^3 + 12/b^5 - 30/b^7 + 105/b^9 ]*exp(-b^2) 
                                      !             ) 

      
                 if(kx>0) then                     ! pour compter -kx
                     bet(kk)=2.d0*bet(kk)
                     bet6(kk)=2.d0*bet6(kk)
                     bet12(kk)=2.d0*bet12(kk)
                 end if
      
                 if(ky.eq.0) then                           ! iq sert a savoir si l'un des ki est nul
                     if(kz.eq.0) iq(kk)=1
                     if(kz.gt.0) iq(kk)=3
                 else
                     if(kz.eq.0) iq(kk)=2
                     if(kz.gt.0) iq(kk)=4
                 end if
                 ! iq = 1  -->  ky=0, kz=0     
                 ! iq = 2  -->  ky>0, kz=0
                 ! iq = 3  -->  ky=0, kz>0
                 ! iq = 4  -->  ky>0, kz>0

 
                 ip(kk)=0
            end do  ! kz
      
            ip(kk)=1              ! repere le changement de ky au cas suivant
      
        end do ! ky
      
        ip(kk)=2     ! repere le changement de kx au cas suivant
      
      end do ! kx
      
      ip(kk)=3
      
      kkmax=kk
      
!      WRITE(*,*) 'nombre de k : ',kkmax
      

!    e.g. bet(kk) =  1/2pi/k^2 * exp(-pi^2 / alpha^2 * k^2 )


!      return 
      
      if(kkmax>nkmax) PRINT*, 'ATTENTION: dépasse la taille des tableaux!!!!'
      !
!        puis l'energie Ewald k  astuce luc70p35
!        il faut somme de Zi*exp(2ipi/L*ri*k)  luc70p35
!        et somme sur les i=O uniquement de exp()  luc84p148
      do kk=1,kkmax
         sumcos(kk)=0.
         sumsin(kk)=0.
         sumcos1(kk)=0.
         sumsin1(kk)=0.
         sumcos2(kk)=0.
         sumsin2(kk)=0.
         sumcos3(kk)=0.
         sumsin3(kk)=0.
         sumcos_o(kk)=0.
         sumsin_o(kk)=0.
         sumcos1_o(kk)=0.
         sumsin1_o(kk)=0.
         sumcos2_o(kk)=0.
         sumsin2_o(kk)=0.
         sumcos3_o(kk)=0.
         sumsin3_o(kk)=0.
      end do
!
      pi2 = two_pi

      do i=1,n                                   ! debut i
! do i=3,3
!
         iio=3*((i-1)/3)+1
         chi=+1.d0

         if(i==iio) chi=-2.d0

         cx=COS(pi2*xx(i))  ! cos ( 2 pi *1 /L * x)  where L=1
         sx=SIN(pi2*xx(i))  
         cy=COS(pi2*yy(i))
         sy=SIN(pi2*yy(i))
         cz=COS(pi2*zz(i))
         sz=SIN(pi2*zz(i))

         !write(*,*) 'x',xx(i),'cx:',cx,'y',yy(i),'cy',cy,'z',zz(i),'cz:',cz

         cckx(0)=1.
         sskx(0)=0.
         ccky(0)=1.
         ssky(0)=0.
         cckz(0)=1.
         sskz(0)=0.

          ! Table of the functions
          ! cos(2pi m /L x) = cos( x*k_m )
          ! sin(k_m x), cos(k_m y) ... 
         do k=1,kmax
            cckx(k)=cckx(k-1)*cx-sskx(k-1)*sx   
            sskx(k)=sskx(k-1)*cx+cckx(k-1)*sx
            ccky(k)=ccky(k-1)*cy-ssky(k-1)*sy                  ! tableaux suivant les 3D
            ssky(k)=ssky(k-1)*cy+ccky(k-1)*sy
            cckz(k)=cckz(k-1)*cz-sskz(k-1)*sz
            sskz(k)=sskz(k-1)*cz+cckz(k-1)*sz
         end do
 
!         write(tmpstr,'(AI3.3A)') 'cck',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!         call write_xyz_array(hfile,cckx,ccky,cckz,kmax)
!         call io_close(hfile)
!
!         write(tmpstr,'(AI3.3A)') 'ssk',i,'.txt'
!         hfile = io_open(tmpstr,'w')
!         call write_xyz_array(hfile,sskx,ssky,sskz,kmax)
!         call io_close(hfile)

     

!        write(*,*) 'cckx:'
!        call write_xyz_array(0,cckx,ccky,cckz,kmax)
!        write(*,*) 'sskx:'
!        call write_xyz_array(0,sskx,ssky,sskz,kmax)

         kx=0
         ky=0
         kz=0
         ! chi - charge: for H chi=1, for O chi=-2
         cckxky=chi                               ! kx,ky   chi * cos (x*kx+ y*ky) = chi cos 0  
         sskxky=0.                                !         chi * sin(x*kx+y*ky) = chi shi 0 
         cckxkya=chi                              ! kx,-ky  chi * cos( 0 - 0)
         sskxkya=0.                               !         chi * sin( 0 - 0)
  !
         do kk=1,kkmax                           ! debut k

            kz=kz+1

            q1=cckxky*cckz(kz)
            q2=sskxky*sskz(kz)
            q3=sskxky*cckz(kz)
            q4=cckxky*sskz(kz)

!            write(*,*) 'LUCkx:',kx,'ky:',ky,'kz:',kz
        !    write(*,*) 'LUC:',q1,q2,q3,q4

            ! x*kx + y*ky + z*kz
            ! cos(x*kx + y*ky)cos(z*kz) - sin(x*kx + y*ky) sin(z*kz) = 
            ! = cos(x*kx + y*ky + z*kz )        
             ! sin(x*kx+y*ky)cos(z*kz) + cos(x*kx+y*ky)sin(z*kz) = sin(x*kx+y*ky+z*kz)                                               
            sumcos(kk)=sumcos(kk)+q1-q2       ! ky,kz
            sumsin(kk)=sumsin(kk)+q3+q4       

            ! kx+ky - kz
            sumcos1(kk)=sumcos1(kk)+q1+q2     ! ky,-kz  cos
            sumsin1(kk)=sumsin1(kk)+q3-q4     !         sin

            if(i==iio) then                 ! pour O, ajouter LJ
               sumcos_o(kk)=sumcos_o(kk)+q1-q2       ! ky,kz                     contient zi=-2!
               sumsin_o(kk)=sumsin_o(kk)+q3+q4
               sumcos1_o(kk)=sumcos1_o(kk)+q1+q2     ! ky,-kz
               sumsin1_o(kk)=sumsin1_o(kk)+q3-q4
            end if

            ! cckxkya = chi cos(kx-ky) 
            q1=cckxkya*cckz(kz)
            q2=sskxkya*sskz(kz)
            q3=sskxkya*cckz(kz)
            q4=cckxkya*sskz(kz)

            !  kx - ky + kz
            sumcos2(kk)=sumcos2(kk)+q1-q2         ! -ky,kz
            sumsin2(kk)=sumsin2(kk)+q3+q4

            ! kx - ky - kz
            sumcos3(kk)=sumcos3(kk)+q1+q2         ! -ky,-kz
            sumsin3(kk)=sumsin3(kk)+q3-q4

            if(i==iio) then                 ! pour O, ajouter LJ
               sumcos2_o(kk)=sumcos2_o(kk)+q1-q2       ! ky,kz
               sumsin2_o(kk)=sumsin2_o(kk)+q3+q4
               sumcos3_o(kk)=sumcos3_o(kk)+q1+q2     ! ky,-kz
               sumsin3_o(kk)=sumsin3_o(kk)+q3-q4
            end if

             
            ! ip - end of row indicator (filled in the cycle before)
            ! ip(kk) == 1  - end of row in kz
            ! ip(kk) == 2  - end of row in ky
            if(ip(kk).eq.1.or.ip(kk).eq.2) then       ! fin de boucle en kz et peut-etre en ky

               if(ip(kk).eq.1) then
                  ky=ky+1
               else
                  ky=0
                  kx=kx+1
               end if
      
               q1=cckx(kx)*ccky(ky)*chi
               q2=sskx(kx)*ssky(ky)*chi
               q3=sskx(kx)*ccky(ky)*chi
               q4=cckx(kx)*ssky(ky)*chi
      
               cckxky=q1-q2 ! chi*( cos(x*kx)cos(y*ky) - sin(x*kx)sin(y*ky) ) = chi * cos(x*kx+y*ky)
               sskxky=q3+q4 ! chi*( sin(x*kx)cos(y*ky) + cos(x*kx)sin(y*ky) ) = chi * sin(x*kx+y*ky)
               cckxkya=q1+q2! chi*( cos(x*kx)cos(y*ky) + sin(x*kx)sin(y*ky) ) = chi * cos(x*kx-y*ky)
               sskxkya=q3-q4! chi*( sin(x*kx)cos(y*ky) - cos(x*kx)sin(y*ky) ) = chi * sin(x*kx-y*ky)

!               write(*,*) 'LUC ip',kk,kx,ky,cckx(kx),ccky(ky)

               kz=-1
            end if
          end do                                  ! fin kk
  !
       end do                                   ! fin i
!                                                                attention, les _o contiennent zi=zO=-2
      sumcos_o(1:kkmax)=sumcos_o(1:kkmax)/(-2.d0)
      sumsin_o(1:kkmax)=sumsin_o(1:kkmax)/(-2.d0)
      sumcos1_o(1:kkmax)=sumcos1_o(1:kkmax)/(-2.d0)
      sumsin1_o(1:kkmax)=sumsin1_o(1:kkmax)/(-2.d0)
      sumcos2_o(1:kkmax)=sumcos2_o(1:kkmax)/(-2.d0)
      sumsin2_o(1:kkmax)=sumsin2_o(1:kkmax)/(-2.d0)
      sumcos3_o(1:kkmax)=sumcos3_o(1:kkmax)/(-2.d0)
      sumsin3_o(1:kkmax)=sumsin3_o(1:kkmax)/(-2.d0)
!

      tmpstr = 'sumcos_o.txt'
      hfile = io_open(tmpstr,'w')
      call write_real_array(hfile,sumcos_o,kkmax)
      call io_close(hfile)

      tmpstr='bet_Luc.txt'
      hfile=io_open(tmpstr,'w')
      call write_real_array(hfile,bet,kkmax)
      call io_close(hfile)


      tmpstr='sumcos_ppp_Luc.txt'
      hfile=io_open(tmpstr,'w')
      call write_real_array(hfile,sumcos,kkmax)
      call io_close(hfile)

      tmpstr='sumsin_ppp_Luc.txt'
      hfile=io_open(tmpstr,'w')
      call write_real_array(hfile,sumsin,kkmax)
      call io_close(hfile)

      tmpstr='bet6.txt'
      hfile=io_open(tmpstr,'w')
      call write_real_array(hfile,bet6,kkmax)
      call io_close(hfile)

      tmpstr='bet12.txt'
      hfile=io_open(tmpstr,'w')
      call write_real_array(hfile,bet12,kkmax)
      call io_close(hfile)      
 
!      return

      utotk=0.
      utotk6=0.; utotk12=0.
      do kk=1,kkmax
         som=sumcos(kk)**2+sumsin(kk)**2

                 ! iq = 1  -->  ky=0, kz=0     
                 ! iq = 2  -->  ky>0, kz=0
                 ! iq = 3  -->  ky=0, kz>0
                 ! iq = 4  -->  ky>0, kz>0
         if(iq(kk).eq.3.or.iq(kk).eq.4) som=som+sumcos1(kk)**2+sumsin1(kk)**2        ! kz>0                           ! kz>0
         if(iq(kk).eq.2.or.iq(kk).eq.4) som=som+sumcos2(kk)**2+sumsin2(kk)**2        ! ky>0                           ! ky>0
         if(iq(kk).eq.4) som=som+sumcos3(kk)**2+sumsin3(kk)**2                       ! ky,kz>0                           ! ky,kz>0

         utotk=utotk+bet(kk)*som
         som_o=sumcos_o(kk)**2+sumsin_o(kk)**2

!         write(*,*) 'som=',som,'bet=',bet(kk)   
 
         if(iq(kk).eq.3.or.iq(kk).eq.4) som_o=som_o+sumcos1_o(kk)**2+sumsin1_o(kk)**2        ! kz>0                           ! kz>0
         if(iq(kk).eq.2.or.iq(kk).eq.4) som_o=som_o+sumcos2_o(kk)**2+sumsin2_o(kk)**2        ! ky>0                           ! ky>0
         if(iq(kk).eq.4) som_o=som_o+sumcos3_o(kk)**2+sumsin3_o(kk)**2                       ! ky,kz>0                           ! ky,kz>0
 
!         write(*,*) 'som_o:',som_o,'bet6:',bet6(kk)
         utotk6=utotk6+bet6(kk)*som_o
         utotk12=utotk12+bet12(kk)*som_o
 
      end do

      cdiel_ext = external_permutivity
!      write(*,*) 'cdiel_ext=',cdiel_ext
      uu_ew=uu_ew+xxclb*utotk
      ureste=pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)-alp2pi*n*xxclb
!      write(*,*) 'ureste=',ureste
      uu_ew=uu_ew+ureste
      uu_lj6=uu_lj6-xlj*sigma2**3*(utotk6+pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
      uu_lj12=uu_lj12+xlj*sigma2**6*(utotk12+pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)

!      write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew
!
!      write(*,*) 'xclb=',xclb,'xxclb=',xxclb
!      write(*,*) 'utotk=',xxclb*utotk,'utotk6=',utotk6,'utotk12=',utotk12
!      write(*,*) 'xlj:',xlj,'sigma6:',sigma2**3,'xlj*sigma6*utotk6',xlj*sigma2**3*utotk6
!      write(*,*) 'xlj:',xlj,'sigma12:',sigma2**6,'xlj*sigma12*utotk12',xlj*sigma2**6*utotk12


uu_lj6_k = xlj*sigma2**3 * utotk6
uu_lj12_k = xlj*sigma2**6 * utotk12
uu_lj_k = uu_lj12_k - uu_lj6_k

uu_ext = pi2*xxclb/(2.d0*cdiel_ext+1.d0)*(summux**2+summuy**2+summuz**2)

!write(*,*) 'alp2pi=',alp2pi,'alpha',alpha,'alpha/sqrt(pi)',alpha/sqrt(pi)
uu_tail_ew =  -alp2pi*n*xxclb  ! alp2pi = 2 alp / sqrt(pi)
! alp2pi * n *xxclb = xclb * q^2 6 alpha n_h2o / sqrt(pi) 
! SUM q_s = 2 n_h2o + n_h2o * 4 = 6 n_h2o
uu_tail_lj6 = xlj*sigma2**3*(pi**1.5*alp**3/6.d0*n_h2o**2-alp**6/12.d0*n_h2o)
uu_tail_lj12 = xlj*sigma2**6*(pi**1.5*alp**9/1080.d0*n_h2o**2-alp**12/1440.d0*n_h2o)
uu_tail_lj = uu_tail_lj12 - uu_tail_lj6

write(*,*) '******************** LUC *******************************'
write(*,*) 'U_REAL ew:',uu_ew_r,'lj:',uu_lj_r,'lj6:',uu_lj6_r,'lj12:',uu_lj12_r
write(*,*) 'U_FOUR ew:',xxclb*utotk,'lj:',uu_lj_k,'lj6:',uu_lj6_k,'lj12:',uu_lj12_k
write(*,*) 'U_TAIL ew:',uu_tail_ew,'lj:',uu_tail_lj,'lj6:',uu_tail_lj6,'lj12:',uu_tail_lj12
write(*,*) 'U_EXT  ew:',uu_ext
write(*,*) 'U_TOTL ew:',uu_ew,'lj:',uu_lj12 + uu_lj6,'lj6:',-uu_lj6,'lj12:',uu_lj12
write(*,*) 'U_TOT tot:',uu_ew + uu_lj6 + uu_lj12

return


ssk_vrai=kmax*pi/alp
!      return   
 
      !      je peux essayer des corrections systematiques  luc84p151
 
      qmax=2.d0*pi*kmax
      uu_ew=uu_ew+n_h2o/pi*xclb*alp*EXP(-ssk_vrai**2)/ssk_vrai*(6.d0-8.d0*SIN(qmax*roh)/(qmax*roh)+2.d0*SIN(qmax*rhh)/(qmax*rhh))
  
      uu_lj6=uu_lj6-epsi_lj*sigma**6*         &
         (2.d0*pi*n_h2o**2*alp**3*EXP(-ssr**2)/ssr+2.d0/SQRT(pi)*n_h2o*alp**6*EXP(-ssk_vrai**2)/ssk_vrai)
      !  
      ! SUM_sj A_sj * 2 pi alpha^3  
 
      uu_lj12=uu_lj12+epsi_lj*sigma**12*      &
         (pi/30.d0*n_h2o**2*alp**9*EXP(-ssr**2)/ssr+1.d0/30.d0/SQRT(pi)*n_h2o*alp**12*EXP(-ssk_vrai**2)/ssk_vrai)
  
      PRINT*, 'et apres corrections des queues en r et en k:'
            write(*,*) 'Energies LJ12,LJ6,coulomb actuelles : ',uu_lj12,uu_lj6,uu_ew

      !      uself deja retranchee dans la boucle Ewald r
      uu_lj=uu_lj12+uu_lj6
      uutot=uu_lj+uu_ew

      write(*,*) 'Energie totale actuelle : ',uutot

      return 

      !xlj=xlj0       ! sert pour le test  NON, plus maintenant
      !
      !             pour les forces en k, on doit refaire une boucle en i et k
      !             luc84p133
      !             NON, plus maintenant, puisque fx,y,z ne contient que la partie en r luc84p187
      !
      !
      !          je definis les mnlmunu
      !

      ialp=0
      do mnmaxi=0,mnmax             ! boucle en max(m,n) a l'exterieur
      do m0=0,mnmax                      ! ainsi 1,1 sera avant 0,2 et 2,0

            mudeb=0
            mumax=m0
            ndeb=m0
    
            do n0=ndeb,mnmaxi
                
                if(m0/=mnmaxi.and.n0/=mnmaxi) cycle             ! m ou n egal a max(m,n)

                nudeb=0
                numax=n0
                
                do l0=ABS(m0-n0),m0+n0
                do mu0=mudeb,mumax

                      if((-1)**mu0==-1) cycle
                      if(m0==n0) nudeb=mu0

                      do nu0=nudeb,numax
  
                          if((-1)**nu0==-1) cycle
                          if(mu0==0.and.nu0==0.AND.(-1)**(m0+n0+l0)==-1) cycle
  
                          ialp=ialp+1
                          mm(ialp)=m0
                          nn(ialp)=n0
                          ll(ialp)=l0
                          mumu(ialp)=mu0
                          nunu(ialp)=nu0
            !              WRITE(*,*) ialp,m0,n0,l0,mu0,nu0
                          if(mu0>0.and.nu0>0.AND.(m0/=n0.OR.(mu0.ne.nu0).OR.(-1)**l0==1)) then
                             ialp=ialp+1                  ! - nu
                             mm(ialp)=m0
                             nn(ialp)=n0
                             ll(ialp)=l0
                             mumu(ialp)=mu0
                             nunu(ialp)=-nu0
               !              WRITE(*,*) ialp,m0,n0,l0,mu0,-nu0
                          end if
                      end do ! nu0

               end do ! mu0
               end do ! l0

            end do ! n0

      end do ! m0
      end do ! mnmaxi

      PRINT*, 'nombre de projections mnlmunu: ',ialp

      !
      !          je definis les mnmunukhi
      !
      ialp=0

      do mnmini=0,mnmax              ! boucle en min(m,n) a l'exterieur
      do m0=mnmini,mnmax                             ! ainsi 1,1 sera apres 0,2 et 2,0

         mudeb=0
         mumax=m0
         ndeb=m0

         do n0=ndeb,mnmax
           
              if(m0/=mnmini.and.n0/=mnmini) cycle             ! m ou n egal a min(m,n)
           
              nudeb=0
              numax=n0

              do khi0=-MIN(m0,n0),MIN(m0,n0)             ! khi et non l !
              do mu0=mudeb,mumax

                  if((-1)**mu0==-1) cycle

                  if(m0==n0) nudeb=mu0

                  do nu0=nudeb,numax

                     if((-1)**nu0==-1) cycle
                     if(mu0==0.and.nu0==0.AND.khi0<0) cycle

                     ialp=ialp+1
                     mm1(ialp)=m0
                     nn1(ialp)=n0
                     khi1(ialp)=khi0
                     mumu1(ialp)=mu0
                     nunu1(ialp)=nu0
      !              if(itest==0) WRITE(*,*) ialp,m0,n0,mu0,nu0,khi0
                     if(mu0>0.and.nu0>0.AND.(m0/=n0.OR.(mu0.ne.nu0).OR.khi0>=0)) then
                        ialp=ialp+1              ! -nu
                        mm1(ialp)=m0
                        nn1(ialp)=n0
                        khi1(ialp)=khi0
                        mumu1(ialp)=mu0
                        nunu1(ialp)=-nu0
          !              if(itest==0) WRITE(*,*) ialp,m0,n0,mu0,-nu0,khi0
                     end if

                  end do ! nu0

              end do ! mu0
              end do ! khi0

         end do ! n0

      end do ! m0
      end do ! mnmini

      PRINT*, 'nombre de projections mnmunukhi: ',ialp

!
!        initialisation du generateur de nombres aleatoires
!
!      x=alea(0)
!        on y va !
      return
!
end
!
subroutine rot_vect_OLD(x,y,z,x1,y1,z1,cc,ss,i) !DOC
IMPLICIT REAL(8) (a-h,o-z)
!
!DOC        fait la rotation de x,y,z en x1,y1,z1 d'angle cc=cos,ss=sin autour d'un (i) des 3 axes principaux
!DOC        attention: x1,y1,z1 peut etre a la meme place memoire que x,y,z
!DOC        luc80p176
!
if(i==1) then
x1=x
y2=y
y1=cc*y2-ss*z
z1=ss*y2+cc*z
end if
if(i==2) then
y1=y
x2=x
x1=cc*x2+ss*z
z1=-ss*x2+cc*z
end if
if(i==3) then
z1=z
x2=x
x1=cc*x2-ss*y
y1=ss*x2+cc*y
end if
end subroutine
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
subroutine erfc_Luc_bet6_bet12_OLD(x,x2,e2,erfk,bet6,bet12) !DOC
      IMPLICIT REAL(8) (a-h,o-z)
!DOC  x, x2 --> b, b^2
!DOC  e2    --> exp(-b^2)
!DOC  erfk  --> erfc(b)
!
!DOC  b^2 = pi^2 h_m^2 / alpha  --> h_m = sqrt(alpha) b/pi
! 
!          avec erfc_Luc a x<0.5, DL a x>8, integration numerique de x a 8 entre les 2
!          luc84p150
!

pi=4.d0*ATAN(1.d0)
pi12=SQRT(pi)
pi15=pi**1.5
xmin=0.5d0
xmax=8.d0
if(x<=xmin) then
   erfk=erfc_Luc_OLD(x,x2,e2)

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
      FUNCTION erfc_Luc_OLD(X,x2,ee2)
      IMPLICIT REAL(8) (a-h,o-z)
      parameter (P=0.3275911d0,A1=0.254829592d0,        &
       A2=-0.284496736d0,A3=1.421413741d0,A4=-1.453152027d0,   &
       A5=1.061405429d0,pi05=1.7724538509055d0)
!        calcule erfc_Luc(x)=2/racine(pi) integrale de x à infini de exp(-t**2)dt
!        x2=x**2 et ee2=exp(-x**2)
!        luc85p108
      if(x<0.455d0) then                ! DL de erf si x petit
      erfc_Luc_OLD=1.d0-2.d0/pi05*x*(1.d0-x2*(1.d0/3.d0-x2*(0.1d0-x2*(1./42.d0-x2/216.d0))))
      ELSEif(x>3.d0)  then              ! Devlpt asymptotique si x grand
      ax=(1.d0-0.5d0*(1.d0-1.5d0/x2)/x2)/(pi05*x)
      erfc_Luc_OLD=ax*ee2
          else                          ! Abramowitz Stegun sinon
      T=1.D0/(1.D0+P*x)
      AX=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
      erfc_Luc_OLD=ax*ee2
          end if
      RETURN
      END 
!
!
      FUNCTION erf_Luc_OLD(X,x2,ee2)
      IMPLICIT REAL(8) (a-h,o-z)
      parameter (P=0.3275911d0,A1=0.254829592d0,        &
       A2=-0.284496736d0,A3=1.421413741d0,A4=-1.453152027d0,   &
       A5=1.061405429d0,pi05=1.7724538509055d0)
!        calcule erf_Luc(x)=2/racine(pi) integrale de 0 à x de exp(-t**2)dt
!        x2=x**2 et ee2=exp(-x**2)
!        luc85p108
      if(x<0.455d0) then                ! DL de erf si x petit
      erf_Luc_OLD=2.d0/pi05*x*(1.d0-x2*(1.d0/3.d0-x2*(0.1d0-x2*(1./42.d0-x2/216.d0))))
      ELSEif(x>3.d0)  then              ! Devlpt asymptotique de erfc_Luc=1-erf si x grand
      ax=(1.d0-0.5d0*(1.d0-1.5d0/x2)/x2)/(pi05*x)
      erf_Luc_OLD=1.d0-ax*ee2
          else                          ! Abramowitz Stegun pour erfc_Luc sinon
      T=1.D0/(1.D0+P*x)
      AX=T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))
      erf_Luc_OLD=1.d0-ax*ee2
          end if
      RETURN
      END 

End Module MCLuc
