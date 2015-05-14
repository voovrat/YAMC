Module MCAccumLuc !DOC 
use AtomicData
use LJTypes
use composition
use Molecule
use MoleculeHandler
!DOC !FILE Module which contains the functions necessary for the projection calculations
!DOC  Is the interface to the functions written by Luc Belloni for the projection calculation
!DOC 
!DOC  NOTE: quantities like pressure, compressibility etc are calculated incorrectly
!DOC NOTE: density is not updates and thus is calculated incorrectly (always 0.0332891 particles/nm^3).
!DOC to get a correct assymptote, the projections should be re-normalized to the real inverse density mean(V^-1) which can be get from frames.dat file 


real(8) :: p_a

contains
!
!subroutine accumu_set_zeros
!        implicit REAL(8) (a-h,o-z)
!  
!  integer * 4 nittot,naccep,ncumul
!    common/moyenn/uutotm,duutot,evtest,devtes,mmu
!      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
!
!
!ncumul = 0
!!
!duutotm = 0.
!uutotm = 0.
!ncumul = 0 
!
! xmtotx = 0.
! xmtoty = 0.
! xmtotz = 0.
!
! xmtot2 = 0.
! xmtot22 = 0.
!
!
!
!end subroutine


subroutine accumu_init(comp,box_length,MAX_MN) !DOC
!DOC initialize the accumulation 
     use parameters, only : kT_kcal_mol, pressure, alpha, sr, sk,external_permutivity
     use constants, only : boltz, epsi0, elec
     use Molecule
!      implicit REAL(8) (a-h,o-z)
     implicit none
!DOC Parameters:
      integer :: MAX_MN !DOC maximum value for m,n 

!      Type(TAtomicData) :: atomic_data
!      Type(TLJTypes) :: lj_types  ! are bad because values are in length units
      Type(TComposition) :: comp !DOC system composition
      real(8) :: box_length !DOC box length
 

     Type(TMolecule),pointer :: water_ptr,Na_ptr,Cl_ptr
      integer :: hwater,hNa,hCl

!!!!!!!!!!! Commons
      REAL * 8 iacc,iacc_o

      integer(4) :: ncstmax,ncrsmax, nomax
      PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O

       real(8) :: nittot,naccep
       integer :: ncumul
       common/iteri/nittot,naccep,ncumul   ! note : ncumul is updated in accumu

      integer(4) :: ncst_ion,ncst_o,ncst_h,ncrs,ij_tab,ialpdeb,ialptot
      COMMON/melange/ncst_ion,ncst_o,ncst_h,ncrs,ij_tab(nomax,nomax),ialpdeb(ncrsmax),ialptot(ncrsmax)

      integer(4) :: n_atomes, nmol, nion_tot, n_h2o, nmol_tot
      common/nombre/n_atomes,nmol(ncstmax),nion_tot,n_h2o,nmol_tot

      real(8) :: xx,yy,zz
      integer :: i_espece
      common/positi/xx(10000),yy(10000),zz(10000),i_espece(10000)        ! positions

      real(8) :: p
      integer :: num
      common/numeri/p,num

      integer :: mnmax,ialpmax
      COMMON/nb_proj/mnmax,ialpmax

      integer :: ialpmaxx,n_omega
      PARAMETER(ialpmaxx=(nomax-1)*nomax/2+(nomax-1)*16+1226)                    ! max pour nmax=6!
      PARAMETER(n_omega=8)

      integer :: mm,nn,ll,mumu,nunu
      COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)

      integer :: mm1, nn1, khi1, mumu1, nunu1
      COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)

      real(8) :: gr, grint
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)

      real(8) :: gr_o,dcosbet, dphi, coeff_a, coeff_b
      common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190

      real(8) :: diam_a,sigma_a, sigma, sigma2, rcutoff_a, rcut2, roh_a, theta_d 
      common/diamet/diam_a(ncstmax),sigma_a(ncrsmax),sigma(ncrsmax),sigma2(ncrsmax),rcutoff_a,rcut2,roh_a,theta_d


      real(8) :: conc_h2o, xl_a
      common/conc/conc_h2o,xl_a

      real(8) :: aaa,aaa_bis
      integer :: ichoixaaa
      COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82

      real(8) :: temp, dbjr, charge, alp_a1, alp, rmax2 , cdiel_ext      
      integer :: kmax,kkmax
      COMMON/ewald/temp,dbjr,charge(ncstmax),alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext

      real(8) :: epsi_lj,xlj,ucut
      COMMON/lennardjones/epsi_lj(ncstmax),xlj(ncrsmax),ucut

      real(8) :: pasq
      integer :: nptq, mq
      COMMON/numeriq/pasq,nptq,mq
      real(8) :: pasq_a1

      real(8) :: press_pa, press_a3, rho1, rho2, vol1,vol2, vol_sph_cub
      COMMON/isobare/press_pa,press_a3,rho1,rho2,vol1,vol2,vol_sph_cub(2000)


  integer :: mmu
  real(8) :: uutotm, duutot, evtest, devtes
    common/moyenn/uutotm,duutot,evtest,devtes,mmu

      real(8) ::      xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22



   integer :: nkmax,kmaxmax
   PARAMETER(nkmax=4661,kmaxmax=20)


      real(8) :: roh,rhh,rohh

     integer :: m0,n0,l0,mu0,nu0,khi0
     integer :: mnmaxdei,mnmaxdej,mnmaxi,mnmini,mumax,numax
     integer :: mudeb,ndeb,nudeb
     integer :: n 
  !!!!!!!!!!!!!

   real(8) :: avog    
   real(8) :: alp2pi
   real(8) :: cthet2,sthet2
   real(8) :: dbjr_a
   integer :: i,j,ialp,icst,ij

!   real(8) :: p_a

   real(8) :: pi,pi2,pialp2 ! pi2 = 2*pi

   real(8) :: rmax

   real(8) :: ssr,ssk,ssk_vrai
   real(8) :: theta,theta2,xkmax2
! LUC
avog=6.02214d-4  ! conversion between the molar concentration and particle/A^3
cdiel_ext = external_permutivity
temp = 300.d0
!boltz=1.38065d-23

! SVP
ncumul = 0

duutot = 0.
uutotm = 0.
ncumul = 0 

 xmtotx = 0.
 xmtoty = 0.
 xmtotz = 0.

 xmtot2 = 0.
 xmtot22 = 0.


charge(1) = 1 ! Na
charge(2) = -1 ! Cl
charge(3) = -0.8476 ! O
charge(4) = 0.4238  ! H

ncst_ion = 2  ! number of ionic species
theta_d = 109.47
roh_a = 1
press_pa = pressure
p_a = 0.05  ! dr for the projection calculation
num = 1000  ! (maximum) number of the g(r) and projection bins
pasq_a1 = 0.05 ! dk for the spectrum in A^{-1}

alp = alpha
ssr = sr
ssk = sk

dcosbet = 0.1  ! don't know if it is enough or not
dphi = 0.1
coeff_a = 1.d0
coeff_b = 1.d0

mnmax = MAX_MN ! maximum number in projection summation  should be 1..6


ichoixaaa = 0 ! the "choice" of the user, which quantity to display. I think we don't need it

rho1=0.; rho2=0.; vol1=0.; vol2=0.; vol_sph_cub(1:num)=0. 



! in my case i have:  water(1),Na(2),Cl(3) in composition
! in Luc case : Na(1),Cl(2),water(3)
! i.e.

   nmol(1) = comp % mol_numbers(2)  ! Na
   nmol(2) = comp % mol_numbers(3)  ! Cl
   nmol(3) = comp % mol_numbers(1)  ! water

  

   conc_h2o = 0.0332891 / avog  ! = 55.278

   hWater = comp % mol_types(1)
   hNa = comp % mol_types(2)
   hCl = comp % mol_types(3)

   call MoleculeHandler_getMolecule(hWater,water_ptr)
   call MoleculeHandler_getMolecule(hNa,Na_ptr)
   call MoleculeHandler_getMolecule(hCl,Cl_ptr)

! LUC


ncst_o=ncst_ion+1               ! nbre d'especes moleculaires   ! nsct_o = 3
ncst_h=ncst_ion+2               ! nbre d'especes atomiques      ! ncst_h = 4
nion_tot=SUM(nmol(1:ncst_ion))                   ! nbre d'ion

n_h2o=nmol(ncst_o)          !nmol(3) = n_h2o              ! nbre de H2O
nmol(ncst_h)=2*n_h2o        !nmol(4) = 2*n_h2o = number of hydrogens 

nmol_tot=nion_tot+n_h2o                       ! nbre de molecules
n_atomes=nion_tot+3*n_h2o                               ! nbre total d'atomes
n=n_atomes
!xl_a=(n_h2o/(conc_h2o*avog))**(1.d0/3.d0)  ! initial length of the box. 
                                           ! we don't need it (i think) but let it be
xl_a = box_length

diam_a(ncst_h)=0.                          ! sigma of hyfrogen is zero

! SVP
   ! Na Luc(1)  SVP(2) 
   diam_a(1) = Na_ptr % sigma(1)
   epsi_lj(1) = Na_ptr % epsilon(1) / kT_kcal_mol 
   ! Cl Luc(2) SVP(3)
   diam_a(2) = Cl_ptr % sigma(1)
   epsi_lj(2) = Cl_ptr % epsilon(1) / kT_kcal_mol 
   ! O  Luc(3)  SVP(1)
   diam_a(3) = water_ptr % sigma(1)
   epsi_lj(3) = water_ptr % epsilon(1) / kT_kcal_mol
!  

rcutoff_a=1.d20                         ! sert a rien mais sera dans les fichiers acc
epsi_lj(ncst_h)=0.
roh=roh_a/xl_a
pi=4.d0*ATAN(1.d0)
!boltz=1.38065d-23  ! already in module constant
!elec=1.602176d-19
!epsi0=8.854187d-12

dbjr_a=elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10    ! en A
dbjr=dbjr_a/xl_a

theta=pi/180.d0*theta_d
theta2=theta/2.d0
cthet2=COS(theta2); sthet2=SIN(theta2)
rhh=2.d0*roh*sthet2
charge(ncst_o)=-2.d0*charge(ncst_h)
press_a3=press_pa*1.d-30/(boltz*temp)            ! avec ca, beta*P*V=press_a3*L3 ou L en A
p=p_a/xl_a
pasq=pasq_a1*xl_a
!u_self=xclb*(-4.d0/roh+1.d0/rhh)
!foh_self=-2.d0*xclb/roh**2
!fhh_self=+xclb/rhh**2
!alp=alp_a1*xl_a
      pi2=2.d0*pi
      pialp2=pi**2/alp**2
      alp2pi=2.d0*alp/dSQRT(pi)
      rmax=ssr/alp
      rmax2=rmax**2
      kmax=NINT(alp*ssk/pi)
      ssk_vrai=kmax*pi/alp
      xkmax2=kmax**2
!ucut=xlj*((sigma/rmax)**12-(sigma/rmax)**6)     ! si on s'amuse a revenir a LJ tronque shifte a rmax
PRINT*, 'Nombre de sites en tout : ',n
PRINT*, 'Taille ACTUELLE de la boite L (en A) : ',xl_a
PRINT*, 'soit une concentration (en M) : ',conc_h2o
PRINT*, 'Distance OH en unite L : ',roh
      WRITE(*,*) 'K (en A-1) : ',alp/xl_a
      WRITE(*,*) 'rmax/L, kmax, sk reel pour Ewald= ',rmax,kmax,ssk_vrai
IF(kmax>kmaxmax) PRINT*, 'ATTENTION: kmax dépasse la taille des tableaux!!!!'


ij=0
ialp=1
do i=1,ncst_o  ! 1..3
 do j=i,ncst_o ! 1..3
 ij=ij+1
 ij_tab(i,j)=ij
 ij_tab(j,i)=ij
 sigma_a(ij)=(diam_a(i)+diam_a(j))/2.d0
 sigma(ij)=sigma_a(ij)/xl_a
 sigma2(ij)=sigma(ij)**2
 xlj(ij)=4.d0*SQRT(epsi_lj(i)*epsi_lj(j))
 ialpdeb(ij)=ialp                                ! start index for the projections for the pair ij
 IF(j<ncst_o) ialpmax=1                          ! ion-ion
 IF(i<ncst_o.and.j==ncst_o) then                 ! ion-H2O
IF(mnmax==0) ialpmax=1
IF(mnmax==1) ialpmax=2
IF(mnmax==2) ialpmax=4
IF(mnmax==3) ialpmax=6
IF(mnmax==4) ialpmax=9
IF(mnmax==5) ialpmax=12
IF(mnmax==6) ialpmax=16
 endif
 IF(i==ncst_o) then                              ! H2O-H2O
IF(mnmax==0) ialpmax=1
IF(mnmax==1) ialpmax=4
IF(mnmax==2) ialpmax=27
IF(mnmax==3) ialpmax=79
IF(mnmax==4) ialpmax=250
IF(mnmax==5) ialpmax=549
IF(mnmax==6) ialpmax=1226
 endif
 ialptot(ij)=ialpmax  ! ij = ij_tab(i,j)   ialptot(ij) = number of projections for the pair ij
 ialp=ialp+ialpmax
 end do
end do
ncrs=ij                         ! nombre de paires ij mol-mol
ialpmax=ialp-1                  ! nombre total de projections
!
PRINT*, 'Diametres ij en unite L : ',sigma(1:ncrs)
PRINT*, 'Nombre de projections pour chaque paire : ',ialptot(1:ncrs)
PRINT*, 'qui commencent dans la liste totale a : ',ialpdeb(1:ncrs)
PRINT*, 'Nombre total : ',ialpmax
!
!
i=0                                           ! i_espece(i) donne le type de i
do icst=1,ncst_ion
i_espece(i+1:i+nmol(icst))=icst
i=i+nmol(icst)
end do
i_espece(i+1:n:3)=ncst_o
i_espece(i+2:n:3)=ncst_h
i_espece(i+3:n:3)=ncst_h
!


ialp=0
ij=0
  do i=1,ncst_o
  mnmaxdei=0; IF(i==ncst_o) mnmaxdei=mnmax
    do j=i,ncst_o
    PRINT*, 'i=',i,', j=',j
    ij=ij+1
    mnmaxdej=0; IF(j==ncst_o) mnmaxdej=mnmax
        do mnmaxi=0,MAX(mnmaxdei,mnmaxdej)              ! boucle en max(m,n) a l'exterieur
      do m0=0,MIN(mnmaxdei,mnmaxi)                      ! ainsi 1,1 sera avant 0,2 et 2,0
      mudeb=0
      mumax=0; IF(i==ncst_o) mumax=m0
      ndeb=0; IF(i==j) ndeb=m0
        do n0=ndeb,MIN(mnmaxdej,mnmaxi)
        IF(m0/=mnmaxi.and.n0/=mnmaxi) cycle             ! m ou n egal a max(m,n)
        nudeb=0
        numax=0; IF(j==ncst_o) numax=n0
          do l0=ABS(m0-n0),m0+n0
            do mu0=mudeb,mumax
            IF((-1)**mu0==-1) cycle
            IF(i==j.and.m0==n0) nudeb=mu0
              do nu0=nudeb,numax
              IF((-1)**nu0==-1) cycle
              IF(mu0==0.and.nu0==0.AND.(-1)**(m0+n0+l0)==-1) cycle
              ialp=ialp+1
              mm(ialp)=m0
              nn(ialp)=n0
              ll(ialp)=l0
              mumu(ialp)=mu0
              nunu(ialp)=nu0
!              WRITE(*,*) ialp,m0,n0,l0,mu0,nu0
                IF(mu0>0.and.nu0>0.and.i==ncst_o.and.j==ncst_o.AND.     &         ! -nu
                          (i/=j.or.m0/=n0.OR.(mu0.ne.nu0).OR.(-1)**l0==1)) THEN
                ialp=ialp+1
              mm(ialp)=m0
              nn(ialp)=n0
              ll(ialp)=l0
              mumu(ialp)=mu0
              nunu(ialp)=-nu0
!              WRITE(*,*) ialp,m0,n0,l0,mu0,-nu0
                endif
              end do
            end do
          end do
        end do
      end do
        end do
    end do
  end do
PRINT*, 'nombre de projections mnlmunu: ',ialp
!
!          je definis les mnmunukhi
!
ialp=0
ij=0
  do i=1,ncst_o
  mnmaxdei=0; IF(i==ncst_o) mnmaxdei=mnmax
    do j=i,ncst_o
    PRINT*, 'i=',i,', j=',j
    ij=ij+1
    mnmaxdej=0; IF(j==ncst_o) mnmaxdej=mnmax
        do mnmini=0,min(mnmaxdei,mnmaxdej)              ! boucle en min(m,n) a l'exterieur
      do m0=mnmini,mnmaxdei                             ! ainsi 1,1 sera apres 0,2 et 2,0
      mudeb=0
      mumax=0; IF(i==ncst_o) mumax=m0
      ndeb=mnmini; IF(i==j) ndeb=MAX(mnmini,m0)
        do n0=ndeb,mnmaxdej
        IF(m0/=mnmini.and.n0/=mnmini) cycle             ! m ou n egal a min(m,n)
        nudeb=0
        numax=0; IF(j==ncst_o) numax=n0
          do khi0=-MIN(m0,n0),MIN(m0,n0)             ! khi et non l !
            do mu0=mudeb,mumax
            IF((-1)**mu0==-1) cycle
            IF(i==j.and.m0==n0) nudeb=mu0
              do nu0=nudeb,numax
              IF((-1)**nu0==-1) cycle
              IF(mu0==0.and.nu0==0.AND.khi0<0) cycle
              ialp=ialp+1
              mm1(ialp)=m0
              nn1(ialp)=n0
              khi1(ialp)=khi0
              mumu1(ialp)=mu0
              nunu1(ialp)=nu0
!              WRITE(*,*) ialp,m0,n0,mu0,nu0,khi0
                IF(mu0>0.and.nu0>0.and.i==ncst_o.and.j==ncst_o.AND.     &
                          (i/=j.or.m0/=n0.OR.(mu0.ne.nu0).OR.khi0>=0)) THEN
                ialp=ialp+1
              mm1(ialp)=m0
              nn1(ialp)=n0
              khi1(ialp)=khi0
              mumu1(ialp)=mu0
              nunu1(ialp)=-nu0
!              WRITE(*,*) ialp,m0,n0,mu0,-nu0,khi0
                endif
              end do
            end do
          end do
        end do
      end do
        end do
    end do
  end do
PRINT*, 'nombre de projections mnmunukhi: ',ialp


end subroutine 


subroutine accumu_set_frame(comp,atomic_data,box_length,kspace_sum,rspace_sum,ext_sum,inp_nittot,inp_naccep) !DOC
!DOC initializes the arrays used in projection calculation with the data read from file
!DOC Parameters:
   use AtomicData
   use composition
   use EwaldSumKSpace
   use EwaldSumRealSpace
   use EwaldSumExternal
   use EwaldSumTails
 
   use parameters, only : external_permutivity

   implicit none

      real(8) :: p 
      integer :: num  
      common/numeri/p,num


     integer :: ncstmax,ncrsmax,nomax
     PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O

      integer :: mnmax,ialpmax
      COMMON/nb_proj/mnmax,ialpmax


   Type(TComposition) :: comp !DOC system composition
   Type(TAtomicData) :: atomic_data !DOC coordinates of atoms
   real(8),intent(in) :: box_length !DOC box length
   Type(TEwaldSumKSpace) :: kspace_sum !DOC EwaldSumKSpace structure
   Type(TEwaldSumRealSpace) :: rspace_sum !DOC EwaldSumRealSpace structure
   Type(TEwaldSumExternal) :: ext_sum !DOC EwaldSumExternal
   real(8) :: inp_nittot, inp_naccep !DOC total number of iterations and accepted number of iterations


   integer :: nion,nwater

!     implicit REAL(8) (a-h,o-z)
   real(8) :: xx,yy,zz
   integer :: i_espece
   common/positi/xx(10000),yy(10000),zz(10000),i_espece(10000)        ! positions

   real(8) :: temp, dbjr, charge, alp_a1, alp, rmax2, cdiel_ext
   integer :: kmax, kkmax
   COMMON/ewald/temp,dbjr,charge(ncstmax),alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext

   real(8) :: uu_lj, uu_lj12, uu_lj6, uu_ew, uu_ew_r
   COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r

   real(8) :: nittot,naccep
   integer :: ncumul
   common/iteri/nittot,naccep,ncumul   ! note : ncumul is updated in accumu


   real(8) :: uu,vrnew,uutot ! uu,vrnew - not used,  uutot - total energy
   common/energi/uu(10000),vrnew(10000),uutot


   integer :: nkmax,kmaxmax
   PARAMETER(nkmax=4661,kmaxmax=20)

   real(8) :: bet, sumcos, sumsin, sumcos1, sumsin1, sumcos2, sumsin2, sumcos3,sumsin3,cckx,ccky,cckz,sskx,ssky,sskz
   integer :: ip,iq
   COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                         &
       sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),          &
       sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),        &
       cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
       cckz(0:kmaxmax),sskz(0:kmaxmax)  ! corresponds to the electrostatics only

   real(8) :: u_ext
   real(8) :: u_tail_ew,u_tail_lj6,u_tail_lj12

      real(8) :: conc_h2o, xl_a
      common/conc/conc_h2o,xl_a


 cdiel_ext = external_permutivity 

!!!!!!!!!!!!!!!!!

  xl_a = box_length

  p = p_a / xl_a

   nittot = inp_nittot
   naccep = inp_naccep 

   nion = comp % mol_numbers(2)
   nwater = comp % mol_numbers(1)
   
! in my representation- water is the first one, in Luc's - water is the last one
xx(1:2*nion) = atomic_data % xx( 3*nwater+1:3*nwater+2*nion) 
yy(1:2*nion) = atomic_data % yy( 3*nwater+1:3*nwater+2*nion) 
zz(1:2*nion) = atomic_data % zz( 3*nwater+1:3*nwater+2*nion) 

xx(2*nion+1:2*nion + 3*nwater) = atomic_data % xx(1:3*nwater)
yy(2*nion+1:2*nion + 3*nwater) = atomic_data % yy(1:3*nwater)
zz(2*nion+1:2*nion + 3*nwater) = atomic_data % zz(1:3*nwater)



  kmax = kspace_sum % grid % kmax
  kkmax = kspace_sum % grid % nk
  bet(1:kkmax) = kspace_sum % beta(1:kkmax)/dbjr

  ip(1:kkmax) = kspace_sum % grid % ip(1:kkmax)
  iq(1:kkmax) = kspace_sum % grid % iq(1:kkmax)



  !sumcos  ! kx ky  --> ppp
  !sumcos1 ! ky, -kz -> ppm
  !sumcos2 ! -ky, kz -> pmp
  !sumcos3 ! -ky -kz -> pmm
   sumcos(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumcos_ppp(1:kkmax)
   sumsin(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumsin_ppp(1:kkmax)
  
   sumcos1(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumcos_ppm(1:kkmax)
   sumsin1(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumsin_ppm(1:kkmax)

   sumcos2(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumcos_pmp(1:kkmax)
   sumsin2(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumsin_pmp(1:kkmax)

   sumcos3(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumcos_pmm(1:kkmax)
   sumsin3(1:kkmax) = kspace_sum % rho_squared_total % sumsincos_coulomb % sumsin_pmm(1:kkmax)

  !      cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &  seems we don't need them) 
  !     cckz(0:kmaxmax),sskz(0:kmaxmax)  ! corresponds to the electrostatics only

   u_ext = EwaldSumExternal_calc_energy(  ext_sum )
   u_tail_ew = ewald_sum_coulomb_tail( atomic_data % charge , atomic_data % natom )
   call ewald_sum_lj_tails( comp, rspace_sum % lj_types, u_tail_lj6, u_tail_lj12 )

   uu_ew = kspace_sum % energy_coulomb + rspace_sum % uu_ew +  u_ext + u_tail_ew
   uu_lj6 = -( kspace_sum % energy_lj + rspace_sum % uu_lj6 + u_tail_lj6 )
   uu_lj12 = rspace_sum % uu_lj12 + u_tail_lj12

   uutot = uu_ew + uu_lj6 + uu_lj12
!!!!!!!!!!!!!!!!!!


end subroutine


!     Last change:  LB   16 Jul 2014    5:11 pm
!     Last change:  LB   13 Jan 2012    3:13 pm
      subroutine accumu !DOC
!DOC do the accumulation 
      use Functions, only : erf_Luc, erfc_Luc
   !   implicit REAL(8) (a-h,o-z)
     implicit none
      REAL * 8 iacc_ss
     integer :: ncstmax,ncrsmax,nomax
     PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O


      integer :: ncst_ion, ncst_o, ncst_h, ncrs, ij_tab, ialpdeb, ialptot ! accum_init
      COMMON/melange/ncst_ion,ncst_o,ncst_h,ncrs,ij_tab(nomax,nomax),ialpdeb(ncrsmax),ialptot(ncrsmax)

      integer :: n,nmol,nion_tot,n_h2o,nmol_tot ! accum_init  
      common/nombre/n,nmol(ncstmax),nion_tot,n_h2o,nmol_tot

      real(8) :: xx,yy,zz
      integer :: i_espece
      common/positi/xx(10000),yy(10000),zz(10000),i_espece(10000)        ! positions

      real(8) :: nittot,naccep
      integer :: ncumul ! set to zero in accum_set_zero,  accumulated here
      common/iteri/nittot,naccep,ncumul

      real(8) :: p  ! dr, accum_init
      integer :: num  ! number of samples in g(r)
      common/numeri/p,num

      real(8) :: uu,vrnew,uutot ! uu,vrnew - not used,  uutot - set in set_frame
      common/energi/uu(10000),vrnew(10000),uutot

      real(8) :: uu_lj, uu_lj12, uu_lj6, uu_ew, uu_ew_r
      COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r

      real(8) :: uutotm, duutot, evtest, devtes
      integer :: mmu
      common/moyenn/uutotm,duutot,evtest,devtes,mmu


      integer :: mnmax, ialpmax
      COMMON/nb_proj/mnmax,ialpmax

      real(8) :: gr_ss
      common/cumul_1/iacc_ss(ncstmax,2000),gr_ss(ncstmax,2000)     ! gsite-H : ion-H, OH et HH

      real(8) :: diam_a, sigma_a, sigma, sigma2, rcutoff_a, rcut2, roh_a, theta_d
      common/diamet/diam_a(ncstmax),sigma_a(ncrsmax),sigma(ncrsmax),sigma2(ncrsmax),rcutoff_a,rcut2,roh_a,theta_d

      real(8) :: epsi_lj, xlj, ucut
      COMMON/lennardjones/epsi_lj(ncstmax),xlj(ncrsmax),ucut

      real(8) :: temp, dbjr, charge, alp_a1, alp, rmax2,  cdiel_ext 
      integer :: kmax, kkmax
      COMMON/ewald/temp,dbjr,charge(ncstmax),alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext

      real(8) :: conc_h2o, xl_a
      common/conc/conc_h2o,xl_a

      real(8) :: press_pa, press_a3, rho1, rho2, vol1,vol2, vol_sph_cub
      COMMON/isobare/press_pa,press_a3,rho1,rho2,vol1,vol2,vol_sph_cub(2000)

      real(8) :: vir, vir1,vir2,vir02,vir12,vir22,hypervir, hypervir2
      common/pression/vir,vir1,vir2,vir02,vir12,vir22,hypervir,hypervir2

      real(8) :: aaa, aaa_bis
       integer :: ichoixaaa
      COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82

      real(8) :: xmtotx, xmtoty, xmtotz, xmtot2, xmtot22
      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22


      integer :: nkmax, kmaxmax
      PARAMETER(nkmax=4661,kmaxmax=20)

      real(8) :: bet, sumcos, sumsin, sumcos1, sumsin1, sumcos2, sumsin2, sumcos3, sumsin3,sskx,ssky,sskz,cckx,ccky,cckz
      integer :: ip,iq
      COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
       sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
       sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
       cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
       cckz(0:kmaxmax),sskz(0:kmaxmax)

      real(8) :: sumcost,sumsint, sumcos1t, sumsin1t, sumcos2t, sumsin2t, sumcos3t, sumsin3t
      COMMON/sommeenktest/                  &         ! pour Ewald en k
       sumcost(nkmax),sumsint(nkmax),sumcos1t(nkmax),sumsin1t(nkmax),       &
       sumcos2t(nkmax),sumsin2t(nkmax),sumcos3t(nkmax),sumsin3t(nkmax)

      real(8) :: pasq
      integer :: nptq, mq
      COMMON/numeriq/pasq,nptq,mq

      real(8) :: sq1, sq2, sq3, sk1,sk2,sk3
      COMMON/spectre/sq1(500),sq2(500),sq3(500)
      COMMON/spectre1/sk1(500),sk2(500),sk3(500)

      real(8) :: alp2, alp2pi, alpr, alpr2
      real(8) :: ee2,er,ffr,ffr1
      real(8) :: cckxky,cckxkya, cosi, cosi1,cosi2, cosi3, cx, cy, cz
      real(8) :: sskxky,sskxkya, sini, sini1,sini2, sini3, sx, sy, sz
      real(8) :: chargi,chargj, chi 
      real(8) :: dal,dfl  ! da/dl df/dl
      real(8) :: xkdi,xkdi1,xkdi2,xkdi3
      real(8) :: dvol

      real(8) :: xi,xj,xij
      real(8) :: yi,yj,yij
      real(8) :: zi,zj,zij    
      real(8) :: dxi,dyi,dzi,dxj,dyj,dzj 

!!!!!!!!!!!!!!!
      integer :: i,j,io,jo,ies,jes,k,kk,kx,ky,kz
      real(8) :: pi,pi2 
      real(8) :: q1,q2,q3,q4
      real(8) :: r,r2,rmax2_p,rijdij

      real(8) :: som,sumvir,sumhypervir
      real(8) :: summux,summuy,summuz

      real(8) :: vol,vol_old,vol_new
      real(8) :: xl_a_12 ! sqrt(xl_a)
      real(8) :: xmcarre


!
!        on accumule la densit‚ de probabilit‚ g(r)
!        … partir de la position des particules
!
!        luc belloni 4-4-98
!
!        nouveau: j'essaie de calculer y(r) en meme temps que muexc par la methode d'insertion luc78p139
!
!        projections gmnl(r)  luc79p163
!
!        je suis plusieurs variables statistiques en meme temps  luc82p82
!
!        version H2O luc84p124
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        version isobare  luc87p144
!
!        version electrolyte luc87p188
!
      ncumul=ncumul+1
      pi=4.d0*ATAN(1.d0)
      pi2=2.*pi
      alp2pi=2.*alp/SQRT(pi)
      alp2=alp**2
!
rmax2_p=(num*p)**2
dfl=0.          ! pour P :  df/dL=d(LV)/dL=Somme fi.di    luc84p140
dal=0.          ! pour khi : dA/dL=d(L**2*df/dL)/dL
summux=0.;summuy=0.; summuz=0.
!
!               en r, je construis les g(r) et j'accumule P et Khi Ewald r
!
        do 1 i=1,n
!
        ies=i_espece(i)
        io=i; IF(i>nion_tot) io=3*((i-nion_tot-1)/3)+1+nion_tot          ! io repere le site central de la molecule portant i
        chargi=charge(ies)
summux=summux+chargi*xx(i)                  ! dipole total M
summuy=summuy+chargi*yy(i)
summuz=summuz+chargi*zz(i)
        xi=xx(i)
        yi=yy(i)
        zi=zz(i)
        dxi=xi-xx(io)         ! di=ri-rO
        dyi=yi-yy(io)
        dzi=zi-zz(io)
!
         do 1 j=i+1,n
!
        jes=i_espece(j)
        jo=j; IF(j>nion_tot) jo=3*((j-nion_tot-1)/3)+1+nion_tot
        chargj=charge(jes)
        xj=xx(j)
        yj=yy(j)
        zj=zz(j)
        dxj=xj-xx(jo)         ! dj=rj-rO
        dyj=yj-yy(jo)
        dzj=zj-zz(jo)
         xij=xi-xj
         yij=yi-yj
         zij=zi-zj
        if(xij.gt.0.5) xij=xij-1.
        if(xij.lt.-0.5) xij=xij+1.
        if(yij.gt.0.5) yij=yij-1.
        if(yij.lt.-0.5) yij=yij+1.
        if(zij.gt.0.5) zij=zij-1.
        if(zij.lt.-0.5) zij=zij+1.
         r2=xij**2+yij**2+zij**2
         IF(r2.gt.rmax2_p.and.r2>rmax2) cycle
         r=sqrt(r2)
         k=int(r/p)+1
         IF(k<=num.and.jo/=io) then                ! g(r) site-site
         IF(ies<ncst_h.and.jes==ncst_h) iacc_ss(ies,k)=iacc_ss(ies,k)+1.
         IF(jes<ncst_h.and.ies==ncst_h) iacc_ss(jes,k)=iacc_ss(jes,k)+1.
         IF(ies==ncst_h.and.jes==ncst_h) iacc_ss(ies,k)=iacc_ss(ies,k)+1.
         endif
!      pour le viriel et l'hyperviriel Ewald r   luc84p140
         IF(r2<rmax2) then
         alpr=alp*r
         alpr2=alpr**2
         ee2=EXP(-alpr2)
         er=erfc_Luc(alpr,alpr2,ee2)            ! calcule erfc(x) en s'aidant de x**2 et exp(-x**2)
         IF(io==jo) er=-erf_Luc(alpr,alpr2,ee2)    ! retrancher le self eventuellement  luc84p134
         ffr=dbjr*chargi*chargj*(er/r+alp2pi*ee2)/r2
         rijdij=xij*(dxi-dxj)+yij*(dyi-dyj)+zij*(dzi-dzj)
         dfl=dfl+ffr*rijdij
         ffr1=dbjr*chargi*chargj*(3.d0*er/r/r2+alp2pi*ee2*(3.d0/r2+2.d0*alp2))/r2
         dal=dal+ffr1*rijdij**2-ffr*((dxi-dxj)**2+(dyi-dyj)**2+(dzi-dzj)**2)
         endif
!
    1 continue
!
!             boucle en i et k pour P et Khi
!             luc84p140
!
      do kk=1,kkmax
      sumcost(kk)=0.                   ! va contenir T=somme zi.ei.(k.di)
      sumsint(kk)=0.
      sumcos1t(kk)=0.
      sumsin1t(kk)=0.
      sumcos2t(kk)=0.
      sumsin2t(kk)=0.
      sumcos3t(kk)=0.
      sumsin3t(kk)=0.
      end do
!
       do i=1,n                                   ! debut i
!
       ies=i_espece(i)
       io=i; IF(i>nion_tot) io=3*((i-nion_tot-1)/3)+1+nion_tot          ! io repere le site central de la molecule portant i
       IF(i==io) cycle              ! ne concerne que les H!
       chi=charge(ies)
       dxi=xx(i)-xx(io)
       dyi=yy(i)-yy(io)
       dzi=zz(i)-zz(io)
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
        cosi=+q1-q2       ! ky,kz
        sini=+q3+q4
        cosi1=+q1+q2     ! ky,-kz
        sini1=+q3-q4
        xkdi=kx*dxi+ky*dyi+kz*dzi
        sumcost(kk)=sumcost(kk)+cosi*xkdi       ! ky,kz
        sumsint(kk)=sumsint(kk)+sini*xkdi
        xkdi1=kx*dxi+ky*dyi-kz*dzi
        sumcos1t(kk)=sumcos1t(kk)+cosi1*xkdi1     ! ky,-kz
        sumsin1t(kk)=sumsin1t(kk)+sini1*xkdi1
        q1=cckxkya*cckz(kz)
        q2=sskxkya*sskz(kz)
        q3=sskxkya*cckz(kz)
        q4=cckxkya*sskz(kz)
        cosi2=+q1-q2         ! -ky,kz
        sini2=+q3+q4
        cosi3=+q1+q2         ! -ky,-kz
        sini3=+q3-q4
        xkdi2=kx*dxi-ky*dyi+kz*dzi
        sumcos2t(kk)=sumcos2t(kk)+cosi2*xkdi2       ! -ky,kz
        sumsin2t(kk)=sumsin2t(kk)+sini2*xkdi2
        xkdi3=kx*dxi-ky*dyi-kz*dzi
        sumcos3t(kk)=sumcos3t(kk)+cosi3*xkdi3       ! -ky,-kz
        sumsin3t(kk)=sumsin3t(kk)+sini3*xkdi3
        dfl=dfl+dbjr*2.d0*pi2*bet(kk)*(-cosi*sumsin(kk)+sini*sumcos(kk))*xkdi
        dal=dal-dbjr*2.d0*pi2**2*bet(kk)*(cosi*sumcos(kk)+sini*sumsin(kk))*xkdi**2
        IF(iq(kk).eq.3.or.iq(kk).eq.4) then        ! faire -kz
        dfl=dfl+dbjr*2.d0*pi2*bet(kk)*(-cosi1*sumsin1(kk)+sini1*sumcos1(kk))*xkdi1
        dal=dal-dbjr*2.d0*pi2**2*bet(kk)*(cosi1*sumcos1(kk)+sini1*sumsin1(kk))*xkdi1**2
        endif
        IF(iq(kk).eq.2.or.iq(kk).eq.4) then        ! faire -ky
        dfl=dfl+dbjr*2.d0*pi2*bet(kk)*(-cosi2*sumsin2(kk)+sini2*sumcos2(kk))*xkdi2
        dal=dal-dbjr*2.d0*pi2**2*bet(kk)*(cosi2*sumcos2(kk)+sini2*sumsin2(kk))*xkdi2**2
        endif
        IF(iq(kk).eq.4) then        ! faire -ky,-kz
        dfl=dfl+dbjr*2.d0*pi2*bet(kk)*(-cosi3*sumsin3(kk)+sini3*sumcos3(kk))*xkdi3
        dal=dal-dbjr*2.d0*pi2**2*bet(kk)*(cosi3*sumcos3(kk)+sini3*sumsin3(kk))*xkdi3**2
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
!
       do kk=1,kkmax
       som=sumcost(kk)**2+sumsint(kk)**2
       IF(iq(kk).eq.3.or.iq(kk).eq.4) som=som+sumcos1t(kk)**2+sumsin1t(kk)**2        ! kz>0                           ! kz>0
       IF(iq(kk).eq.2.or.iq(kk).eq.4) som=som+sumcos2t(kk)**2+sumsin2t(kk)**2        ! ky>0                           ! ky>0
       IF(iq(kk).eq.4) som=som+sumcos3t(kk)**2+sumsin3t(kk)**2                       ! ky,kz>0                           ! ky,kz>0
       dal=dal+dbjr*2.d0*pi2**2*bet(kk)*som
       end do
!
      xmcarre=summux**2+summuy**2+summuz**2
      dfl=dfl-2.d0*pi2*dbjr/(2.d0*cdiel_ext+1.d0)*xmcarre
      dal=dal+2.d0*pi2*dbjr/(2.d0*cdiel_ext+1.d0)*xmcarre
!
      sumvir=12.d0*uu_lj12+6.d0*uu_lj6+uu_ew-dfl
      sumhypervir=144.d0*uu_lj12+36.d0*uu_lj6+uu_ew-3.d0*dfl+dal
!       en isobare, V varie et c'est W/V qu'il faut moyenner! luc87p145
!       en fait, il faut moyenner P/kT tot = N/V+W/V
!       pour hypervir, je moyenne KHI/V aussi   luc87p146
      vol=xl_a**3
      sumvir=nmol_tot/vol+1.d0/3.d0*sumvir/vol
      sumhypervir=sumhypervir/vol
!
!          attention: il faut moyenner M**2/V donc M'**2/L ou M'=somme zi*si (si=ri/L)   luc87p147
      xl_a_12=SQRT(xl_a)
      xmtotx=xmtotx+summux  /xl_a_12
      xmtoty=xmtoty+summuy  /xl_a_12
      xmtotz=xmtotz+summuz  /xl_a_12
      xmtot2=xmtot2+xmcarre  /xl_a
      xmtot22=xmtot22+xmcarre**2  /xl_a**2
!        energie
      uutotm=uutotm+uutot
      duutot=duutot+(uutot+press_a3*vol)**2             ! accumule en fait (U+PV)**2 pour Cp  luc87p145
!        pression
      vir=vir+sumvir
      vir02=vir02+sumvir**2
!        hyperviriel
      hypervir=hypervir+sumhypervir
      hypervir2=hypervir2+sumhypervir**2
!        volume
      rho1=rho1+conc_h2o
      rho2=rho2+conc_h2o**2
      vol1=vol1+vol
      vol2=vol2+vol**2
      vol_old=0.                                   ! il faut moyenner le volume, intersection avec le cube
       do i=1,num
       vol_new=volume_sph_cub(i*p)                 ! c'est en unite V
       dvol=vol_new-vol_old
       vol_sph_cub(i)=vol_sph_cub(i)+vol*dvol
       vol_old=vol_new
       end do
      IF(ichoixaaa.EQ.1) aaa=uutot/nmol_tot
      IF(ichoixaaa.EQ.2) aaa=sumvir
      IF(ichoixaaa.EQ.3) aaa=conc_h2o
      IF(ichoixaaa.EQ.4) aaa=vol
      IF(ichoixaaa.EQ.7) aaa=sumhypervir/9.d0
      IF(ichoixaaa.EQ.8) aaa=xmcarre  *xl_a**2
      IF(ichoixaaa.EQ.9) aaa=summux  *xl_a
      IF(ichoixaaa.EQ.10) aaa=summuy  *xl_a
      IF(ichoixaaa.EQ.11) aaa=summuz  *xl_a
      aaa_bis(1)=uutot/nmol_tot
      aaa_bis(2)=sumvir
      aaa_bis(3)=conc_h2o
      aaa_bis(4)=vol
      aaa_bis(7)=sumhypervir/9.d0
      aaa_bis(8)=xmcarre  *xl_a**2
      aaa_bis(9)=summux  *xl_a
      aaa_bis(10)=summuy  *xl_a
!
!        j'accumule les gmnmunukhi(r) dans un autre sous-programme
!

      call accumu_mnmunukhi
!
      aaa_bis(80)=aaa           ! je mets la variable aaa d'avant a la fin du tableau aaa_bis
!
      return
end subroutine


!     Last change:  LB   16 Jul 2014    5:11 pm
!     Last change:  LB   13 Jan 2012    3:13 pm
      subroutine accumu_mnmunukhi !DOC
!DOC do the accumulation for projections
!      implicit REAL(8) (a-h,o-z)
      implicit none

      REAL * 8 iacc,iacc_o

      integer :: ncstmax, ncrsmax, nomax
      PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O

      
      integer :: ncst_ion, ncst_o, ncst_h, ncrs, ij_tab, ialpdeb, ialptot
      COMMON/melange/ncst_ion,ncst_o,ncst_h,ncrs,ij_tab(nomax,nomax),ialpdeb(ncrsmax),ialptot(ncrsmax)

      integer :: n_atomes, nmol, nion_tot, n_h2o, nmol_tot
      common/nombre/n_atomes,nmol(ncstmax),nion_tot,n_h2o,nmol_tot

      real(8) :: xx,yy,zz
      integer :: i_espece
      common/positi/xx(10000),yy(10000),zz(10000),i_espece(10000)        ! positions

      real(8) :: p
      integer :: num
      common/numeri/p,num

      integer :: mnmax,ialpmax
      COMMON/nb_proj/mnmax,ialpmax

      integer :: ialpmaxx,n_omega
      PARAMETER(ialpmaxx=(nomax-1)*nomax/2+(nomax-1)*16+1226)                    ! max pour nmax=6!
      PARAMETER(n_omega=8)

      integer :: mm,nn,ll,mumu,nunu
      COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)

      integer :: mm1,nn1,khi1,mumu1, nunu1
      COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)

      real(8) :: gr, grint
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)

      real(8) :: gr_o, dcosbet, dphi, coeff_a, coeff_b 
      common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190

      real(8) :: diam_a, sigma_a, sigma, sigma2, rcutoff_a, rcut2, roh_a, theta_d
      common/diamet/diam_a(ncstmax),sigma_a(ncrsmax),sigma(ncrsmax),sigma2(ncrsmax),rcutoff_a,rcut2,roh_a,theta_d

      real(8) :: conc_h2o, xl_a, aaa, aaa_bis
      integer :: ichoixaaa
      common/conc/conc_h2o,xl_a
      COMMON/effica/aaa,ichoixaaa,aaa_bis(80)                         ! differentes aaa dans aaa_bis luc82p82

      integer,parameter :: max_mnmax = 6
      
      integer,DIMENSION(100) ::  k_bis
      real(8),dimension(-1:2*max_mnmax+1) ::  rac
      real(8) :: rac2 

     REAL(8), DIMENSION(0:max_mnmax,-max_mnmax:max_mnmax,-max_mnmax:max_mnmax):: a,b,c,d,fi,gi,fj,gj
      
      REAL(8),DIMENSION(3,3):: ri,rj,rij,rmati,rmatj


      integer :: io , ih1, ih2, k
      integer ::  i_mol, j_mol, ialp, ialp_aaa, ies
      integer :: ij


      real(8) ::  pi,theta,theta2
      real(8) :: cthet2, sthet2, cosdphi2, sindphi2
      real(8) :: den

      real(8) :: roh,rhh,rohh

      integer :: jes,jo,jh1,jh2
      integer :: k_aaa, k_bis_max

      integer :: kk,kr,l,l1,m,m1,m1min,n,mu,nu,khi
     
      real(8) :: phi,r,r2,rij1_norm,rmax2_p
      real(8) :: sisj,titj
      real(8) :: xh1i,xh1j, xh2i,xh2j, xij, xi, xj
      real(8) :: yh1i,yh1j, yh2i,yh2j, yij, yi, yj
      real(8) :: zh1i,zh1j, zh2i,zh2j, zij, zi, zj

!
!interface
!FUNCTION vector(r1,r2)
!implicit REAL(8) (a-h,o-z)
!REAL(8), DIMENSION(1:3):: r1,r2,vector
!end function
!END interface
!
!
!
!        projections gmnmunukhi(r) accumulees dans iacc  pour H2O
!        luc85p23
!        et g(r,omega) pour qqs omegas caracteristiques luc85p176
!




cosdphi2=COS(dphi/2.d0)
sindphi2=sin(dphi/2.d0)
!
roh=roh_a/xl_a
pi=4.d0*ATAN(1.d0)
theta=pi/180.d0*theta_d
theta2=theta/2.d0
cthet2=COS(theta2); sthet2=SIN(theta2)
rhh=2.d0*roh*sthet2  ! H -H 
!                         sth2
rohh=2.d0*roh*cthet2  !  H-- * -- H
!                         \  |  /
!                          \ |<--- cth2
                      !     \|/
                      !      O 
!  rohh : double distance from the oxygen to the H-H line



rac2=sqrt(2.d0)
rac(-1)=0.                    ! pour eviter des pbs

do k=0,2*mnmax+1
  rac(k)=SQRT(dble(k))                           ! racine(k) en tableau
end do
do l=1,mnmax
  do m=-l,l
    do m1=-l+1,l-1
      a(l,m,m1)=rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l-m1))              ! coeff a,b,c,d lmm'
      b(l,m,m1)=rac(l+m)*rac(l+m-1)/(rac(2)*rac(l+m1)*rac(l-m1))
    end do
    m1=l
    c(l,m,m1)=rac(2)*rac(l+m)*rac(l-m)/(rac(l+m1)*rac(l+m1-1))
    d(l,m,m1)=rac(l+m)*rac(l+m-1)/(rac(l+m1)*rac(l+m1-1))
  end do
end do



k_aaa=0
IF(ichoixaaa>100) then                  ! on veut suivre g_alpha(ri)
aaa=0.
ialp_aaa=ichoixaaa/10000                  ! et non /1000  luc82p166
k_aaa=ichoixaaa-10000*ialp_aaa
endif
                                        ! on suit plein de variables avec aaa_bis
do k=1,9
kk=nint(k*sigma(ncrs)/p)                ! pour h2o-h2o pour l'instant
IF(kk>num) exit
k_bis(k)=kk
end do
k_bis_max=k                ! et le dernier point g(num)
k_bis(k_bis_max)=num
aaa_bis(11:80)=0.

rmax2_p=(num*p)**2
!
! pour test, ne regarder une projection particuliere que d'une seule paire
!PRINT*, 'Test de gmnlmunu'
!PRINT*, 'pour quel m,n,l,mu,nu ?'
!READ(*,*) m0,n0,l0,mu0,nu0
!PRINT*, 'pour quel i,j ?'
!READ(*,*) i0,j0
!               en r, je construis les gmnmunukhi(r)
!




do i_mol=1,nmol_tot             ! molecule i_mol
   !
   !IF(i_mol<=ncst_ion) then                 ! c'est un ion
   if(i_mol<=nion_tot) then
   ies=i_espece(i_mol)
   io=i_mol
   else                                     ! c'est une H2O
   ies=ncst_o
   io=nion_tot+3*(i_mol-nion_tot-1)+1
   ih1=io+1
   ih2=io+2
   endif
   !


  !          construire les coord des vecteurs lies a i_h2o dans repere fixe
  IF(ies==ncst_o) then
        ri(3,1)=(xx(ih1)+xx(ih2)-2.d0*xx(io))/rohh        ! axe principal z suivant O - H1+H2
        ri(3,2)=(yy(ih1)+yy(ih2)-2.d0*yy(io))/rohh
        ri(3,3)=(zz(ih1)+zz(ih2)-2.d0*zz(io))/rohh
        ri(1,1)=(xx(ih2)-xx(ih1))/rhh                      ! axe x suivant H1H2
        ri(1,2)=(yy(ih2)-yy(ih1))/rhh
        ri(1,3)=(zz(ih2)-zz(ih1))/rhh
        ri(2,1:3)=vector(ri(3,1:3),ri(1,1:3))                  ! 3eme axe, perpendiculaire aux premiers
  endif
!


  do j_mol=i_mol+1,nmol_tot
        !
        !IF(i_h2o/=i0.or.j_h2o/=j0) cycle   ! pour test
        !PRINT*, 'i fixe',((ri(k,l),l=1,3),k=1,3)
        !IF(j_mol<=ncst_ion) then                 ! c'est un ion
        if(j_mol<=nion_tot) then
           jes=i_espece(j_mol)
           jo=j_mol
        else                                     ! c'est une H2O
           jes=ncst_o
           jo=nion_tot+3*(j_mol-nion_tot-1)+1
           jh1=jo+1
           jh2=jo+2
        endif
      !
        xij=xx(jo)-xx(io)                               ! vecteur Oi - Oj
        yij=yy(jo)-yy(io)
        zij=zz(jo)-zz(io)

        if ( mnmax < 0) then
            write(*,*) 'jmol',j_mol,'A: mnmax',mnmax
            stop
        end if

        if(xij.gt.0.5) xij=xij-1.
        if(xij.lt.-0.5) xij=xij+1.
        if(yij.gt.0.5) yij=yij-1.
        if(yij.lt.-0.5) yij=yij+1.
        if(zij.gt.0.5) zij=zij-1.
        if(zij.lt.-0.5) zij=zij+1.

        r2=xij**2+yij**2+zij**2
        IF(r2.gt.rmax2_p) cycle
        r=sqrt(r2)

        if ( r < 1e-9 ) then
           write(*,*) 'r=',r,'io=',io,'jo=',jo
           cycle
        end if


        kr=int(r/p)+1
!         construire les coord des vecteurs du repere intermoleculaire ij dans repere fixe
        rij(3,1)=xij/r; rij(3,2)=yij/r; rij(3,3)=zij/r                ! vecteur unitaire rij
        IF(ABS(rij(3,1))<0.99) then                                ! puis x*rij (sauf si rij suivant x)
        rij(1,1:3)=vector((/1.d0,0.d0,0.d0/),rij(3,1:3))
        else                                              ! dans ce cas y*rij
        rij(1,1:3)=vector((/0.d0,1.d0,0.d0/),rij(3,1:3))
        endif
        rij1_norm=SQRT(SUM(rij(1,:)**2))
        rij(1,:)=rij(1,:)/rij1_norm
        rij(2,1:3)=vector(rij(3,1:3),rij(1,1:3))
!PRINT*, 'rij fixe',((rij(k,l),l=1,3),k=1,3)
!          construire les coord des vecteurs lies a j_h2o dans repere fixe
        IF(jes==ncst_o) then
           rj(3,1)=(xx(jh1)+xx(jh2)-2.d0*xx(jo))/rohh        ! axe principal z suivant O - H1+H2
           rj(3,2)=(yy(jh1)+yy(jh2)-2.d0*yy(jo))/rohh
           rj(3,3)=(zz(jh1)+zz(jh2)-2.d0*zz(jo))/rohh
           rj(1,1)=(xx(jh2)-xx(jh1))/rhh                      ! axe x suivant H1H2
           rj(1,2)=(yy(jh2)-yy(jh1))/rhh
           rj(1,3)=(zz(jh2)-zz(jh1))/rhh
           rj(2,1:3)=vector(rj(3,1:3),rj(1,1:3))               ! 3eme axe, perpendiculaire aux premiers
        endif
!PRINT*, 'j fixe',((rj(k,l),l=1,3),k=1,3)
!
!            calcul de g(r,omega) pour qq omega  luc85p176,186,189
!        compter 1/2 quand l'angle entre HHi et HHj peut etre 0 ou pi (A,C,E), pi/2 ou -pi/2 (B,D,F)
!        compter 1/2 quand l'angle entre les normales a i et j peu etre 0 ou pi (G), pi/2 ou -pi/2 (H)
!        compter 1/2 quand on renverse i et j en meme temps (A,B)
!        compter 1/4 quand un des 4 H pointe vers l'autre O (G,H)
!


        IF(ies==ncst_o.and.jes==ncst_o) then        ! h2o-h2o seulement pour l'instant
           xi=DOT_PRODUCT(ri(3,:),rij(3,:))
           xj=DOT_PRODUCT(rj(3,:),rij(3,:))
           sisj=ABS(DOT_PRODUCT(ri(1,:),rj(1,:)))
           titj=ABS(DOT_PRODUCT(ri(2,:),rj(2,:)))

           IF(abs(xi)>1.d0-dcosbet.and.ABS(xj)>1.d0-dcosbet) then            ! A-F
              IF(xi*xj>0.) then
                  IF(sisj>cosdphi2) iacc_o(1,kr)=iacc_o(1,kr)+0.25d0                 ! A
                  IF(sisj<sindphi2) iacc_o(2,kr)=iacc_o(2,kr)+0.25d0                 ! B
              ELSEIF(xi>0.) then
                  IF(sisj>cosdphi2) iacc_o(3,kr)=iacc_o(3,kr)+0.5d0                 ! C
                  IF(sisj<sindphi2) iacc_o(4,kr)=iacc_o(4,kr)+0.5d0                 ! D
              else
                  IF(sisj>cosdphi2) iacc_o(5,kr)=iacc_o(5,kr)+0.5d0                 ! E
                  IF(sisj<sindphi2) iacc_o(6,kr)=iacc_o(6,kr)+0.5d0                 ! F
              endif
!ELSEIF((ABS(xi-costetra)<dcosbet/2.d0.and.xj>1.d0-dcosbet).or.(xi<-1.d0+dcosbet.and.ABS(xi+costetra)<dcosbet/2.d0)) then
!IF(titj>cosdphi2) iacc_o(7,kr)=iacc_o(7,kr)+0.25d0                 ! G
!IF(titj<sindphi2) iacc_o(8,kr)=iacc_o(8,kr)+0.25d0                 ! H
           endif

           xh1i=DOT_PRODUCT((/xx(ih1)-xx(io),yy(ih1)-yy(io),zz(ih1)-zz(io)/),rij(3,:))/roh
           xh2i=DOT_PRODUCT((/xx(ih2)-xx(io),yy(ih2)-yy(io),zz(ih2)-zz(io)/),rij(3,:))/roh
           xh1j=DOT_PRODUCT((/xx(jh1)-xx(jo),yy(jh1)-yy(jo),zz(jh1)-zz(jo)/),rij(3,:))/roh
           xh2j=DOT_PRODUCT((/xx(jh2)-xx(jo),yy(jh2)-yy(jo),zz(jh2)-zz(jo)/),rij(3,:))/roh

           IF((MAX(xh1i,xh2i)>1.d0-dcosbet.and.xj>1.d0-dcosbet).OR.(xi<-1.d0+dcosbet.and.MIN(xh1j,xh2j)<-1.d0+dcosbet)) then
               IF(titj>cosdphi2) iacc_o(7,kr)=iacc_o(7,kr)+0.125d0                 ! G
               IF(titj<sindphi2) iacc_o(8,kr)=iacc_o(8,kr)+0.125d0                 ! H
           endif
       endif  ! ies == ncst_o .AND. jes==ncst_o


       fi(0,0,0)=1.d0; gi(0,0,0)=0.          ! facile!
       IF(ies==ncst_o) then          ! pas la peine d'aller plus loin si ion spherique!
!         coord des vecteurs de i dans le repere intermol
          do k=1,3
             do l=1,3
                rmati(k,l)=DOT_PRODUCT(rij(k,:),ri(l,:))
             end do
          end do
!PRINT*, 'i interm',((rmati(k,l),l=1,3),k=1,3)
!         matrice Messiah Rlmm' pour l=1  = Flmm' + i Glmm'
          fi(1,-1,-1)=(rmati(2,2)+rmati(1,1))/2.d0
          fi(1,-1,0)=rmati(1,3)/rac2
          fi(1,-1,+1)=(rmati(2,2)-rmati(1,1))/2.d0
          fi(1,0,-1)=rmati(3,1)/rac2
          fi(1,0,0)=rmati(3,3)
          fi(1,0,+1)=-rmati(3,1)/rac2
          fi(1,+1,-1)=(rmati(2,2)-rmati(1,1))/2.d0
          fi(1,+1,0)=-rmati(1,3)/rac2
          fi(1,+1,+1)=(rmati(2,2)+rmati(1,1))/2.d0
          gi(1,-1,-1)=(rmati(2,1)-rmati(1,2))/2.d0
          gi(1,-1,0)=rmati(2,3)/rac2
          gi(1,-1,+1)=(-rmati(2,1)-rmati(1,2))/2.d0
          gi(1,0,-1)=-rmati(3,2)/rac2
          gi(1,0,0)=0.
          gi(1,0,+1)=-rmati(3,2)/rac2
          gi(1,+1,-1)=(rmati(2,1)+rmati(1,2))/2.d0
          gi(1,+1,0)=rmati(2,3)/rac2
          gi(1,+1,+1)=(rmati(1,2)-rmati(2,1))/2.d0
!        matrice Messiah Rlmm' pour l>1 par recurrence Choi... luc85p17
          do l=2,mnmax
             l1=l-1
             do m=-l,l

                m1min=0; 
                IF(m>0) m1min=1

                do m1=m1min,l-1
                   fi(l,m,m1)=a(l,m,m1)*(fi(1,0,0)*fi(l1,m,m1)-gi(1,0,0)*gi(l1,m,m1))+         &
                             b(l,m,m1)*(fi(1,+1,0)*fi(l1,m-1,m1)-gi(1,+1,0)*gi(l1,m-1,m1))+   &
                             b(l,-m,m1)*(fi(1,-1,0)*fi(l1,m+1,m1)-gi(1,-1,0)*gi(l1,m+1,m1))
                   gi(l,m,m1)=a(l,m,m1)*(fi(1,0,0)*gi(l1,m,m1)+gi(1,0,0)*fi(l1,m,m1))+         &
                             b(l,m,m1)*(fi(1,+1,0)*gi(l1,m-1,m1)+gi(1,+1,0)*fi(l1,m-1,m1))+   &
                             b(l,-m,m1)*(fi(1,-1,0)*gi(l1,m+1,m1)+gi(1,-1,0)*fi(l1,m+1,m1))
                   fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                   gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
                end do
                m1=l
                fi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*fi(l1,m,m1-1)-gi(1,0,+1)*gi(l1,m,m1-1))+         &
                          d(l,m,m1)*(fi(1,+1,+1)*fi(l1,m-1,m1-1)-gi(1,+1,+1)*gi(l1,m-1,m1-1))+   &
                          d(l,-m,m1)*(fi(1,-1,+1)*fi(l1,m+1,m1-1)-gi(1,-1,+1)*gi(l1,m+1,m1-1))
                gi(l,m,m1)=c(l,m,m1)*(fi(1,0,+1)*gi(l1,m,m1-1)+gi(1,0,+1)*fi(l1,m,m1-1))+         &
                          d(l,m,m1)*(fi(1,+1,+1)*gi(l1,m-1,m1-1)+gi(1,+1,+1)*fi(l1,m-1,m1-1))+   &
                          d(l,-m,m1)*(fi(1,-1,+1)*gi(l1,m+1,m1-1)+gi(1,-1,+1)*fi(l1,m+1,m1-1))
                fi(l,-m,-m1)=(-1)**(m+m1)*fi(l,m,m1)
                gi(l,-m,-m1)=-(-1)**(m+m1)*gi(l,m,m1)
             end do  ! m
          end do ! l
       endif ! ies == ncst_o
       !
       !           idem pour j
       !
       Fj(0,0,0)=1.d0; gj(0,0,0)=0.          ! facile!
       IF(jes==ncst_o) then
       !         coord des vecteurs de j dans le repere intermol
          do k=1,3
             do l=1,3
                rmatj(k,l)=DOT_PRODUCT(rij(k,:),rj(l,:))
             end do
          end do
          !PRINT*, 'j interm',((rmatj(k,l),l=1,3),k=1,3)
          !         matrice Messiah Rlmm' pour l=1  = Flmm' + i Glmm'
          fj(1,-1,-1)=(rmatj(2,2)+rmatj(1,1))/2.d0
          fj(1,-1,0)=rmatj(1,3)/rac2
          fj(1,-1,+1)=(rmatj(2,2)-rmatj(1,1))/2.d0
          fj(1,0,-1)=rmatj(3,1)/rac2
          fj(1,0,0)=rmatj(3,3)
          fj(1,0,+1)=-rmatj(3,1)/rac2
          fj(1,+1,-1)=(rmatj(2,2)-rmatj(1,1))/2.d0
          fj(1,+1,0)=-rmatj(1,3)/rac2
          fj(1,+1,+1)=(rmatj(2,2)+rmatj(1,1))/2.d0
          gj(1,-1,-1)=(rmatj(2,1)-rmatj(1,2))/2.d0
          gj(1,-1,0)=rmatj(2,3)/rac2
          gj(1,-1,+1)=(-rmatj(2,1)-rmatj(1,2))/2.d0
          gj(1,0,-1)=-rmatj(3,2)/rac2
          gj(1,0,0)=0.
          gj(1,0,+1)=-rmatj(3,2)/rac2
          gj(1,+1,-1)=(rmatj(2,1)+rmatj(1,2))/2.d0
          gj(1,+1,0)=rmatj(2,3)/rac2
          gj(1,+1,+1)=(rmatj(1,2)-rmatj(2,1))/2.d0
          !        matrice Messiah Rlmm' pour l>1 par recurrence Choi... luc85p17

          do l=2,mnmax
             l1=l-1
             do m=-l,l
                m1min=0; 
                IF(m>0) m1min=1
                do m1=m1min,l-1
                   fj(l,m,m1)=a(l,m,m1)*(fj(1,0,0)*fj(l1,m,m1)-gj(1,0,0)*gj(l1,m,m1))+         &
                             b(l,m,m1)*(fj(1,+1,0)*fj(l1,m-1,m1)-gj(1,+1,0)*gj(l1,m-1,m1))+   &
                             b(l,-m,m1)*(fj(1,-1,0)*fj(l1,m+1,m1)-gj(1,-1,0)*gj(l1,m+1,m1))
                   gj(l,m,m1)=a(l,m,m1)*(fj(1,0,0)*gj(l1,m,m1)+gj(1,0,0)*fj(l1,m,m1))+         &
                             b(l,m,m1)*(fj(1,+1,0)*gj(l1,m-1,m1)+gj(1,+1,0)*fj(l1,m-1,m1))+   &
                             b(l,-m,m1)*(fj(1,-1,0)*gj(l1,m+1,m1)+gj(1,-1,0)*fj(l1,m+1,m1))
                   fj(l,-m,-m1)=(-1)**(m+m1)*fj(l,m,m1)
                   gj(l,-m,-m1)=-(-1)**(m+m1)*gj(l,m,m1)
                end do ! m1
                m1=l
                fj(l,m,m1)=c(l,m,m1)*(fj(1,0,+1)*fj(l1,m,m1-1)-gj(1,0,+1)*gj(l1,m,m1-1))+         &
                           d(l,m,m1)*(fj(1,+1,+1)*fj(l1,m-1,m1-1)-gj(1,+1,+1)*gj(l1,m-1,m1-1))+   &
                           d(l,-m,m1)*(fj(1,-1,+1)*fj(l1,m+1,m1-1)-gj(1,-1,+1)*gj(l1,m+1,m1-1))
            
                if  ( (l > mnmax) .OR. ( m > mnmax).OR.( m1 > mnmax) .OR. (l1 > mnmax) ) then
                   write(*,*) 'l',l,'l1',l1,'m',m,'m1',m1
                   write(*,*) 'c(l,m,m1)',c(l,m,m1)
                   write(*,*) 'gj(l1,m,m1-1)',gj(l1,m,m1-1)
                   write(*,*) 'fj(l1,m,m1-1)',fj(l1,m,m1-1)
                   write(*,*) 'd(l,m,m1)',d(l,m,m1)
                   write(*,*) 'd(l,-m,m1)',d(l,-m,m1)
                   write(*,*) 'fj(l1,m+1,m1-1)',fj(l1,m+1,m1-1)
                end if
                gj(l,m,m1)=c(l,m,m1)*(fj(1,0,+1)*gj(l1,m,m1-1)+gj(1,0,+1)*fj(l1,m,m1-1))+         &
                          d(l,m,m1)*(fj(1,+1,+1)*gj(l1,m-1,m1-1)+gj(1,+1,+1)*fj(l1,m-1,m1-1))+   &
                          d(l,-m,m1)*(fj(1,-1,+1)*gj(l1,m+1,m1-1)+gj(1,-1,+1)*fj(l1,m+1,m1-1))
               !write(*,*) 'Ugh'
               
                fj(l,-m,-m1)=(-1)**(m+m1)*fj(l,m,m1)
               
               !write(*,*) 'Ogh'
                gj(l,-m,-m1)=-(-1)**(m+m1)*gj(l,m,m1)
             end do ! m
          end do ! l
       endif ! jes == ncst_o
       !
       !        construire alors Reelle(Rmkhimu(i)*Rn-khinu(j)) pour tous les ialp de la liste
       !
       !phitest1=0.

       ij=ij_tab(ies,jes)

       do ialp=ialpdeb(ij),ialpdeb(ij)+ialptot(ij)-1
          m=mm1(ialp); n=nn1(ialp); mu=mumu1(ialp); nu=nunu1(ialp); khi=khi1(ialp)
          phi=fi(m,khi,mu)*fj(n,-khi,nu)-gi(m,khi,mu)*gj(n,-khi,nu)

          IF(ies==jes.AND.(m/=n.or.mu/=nu)) phi=(phi+(-1)**(m+n)*(fi(n,khi,nu)*fj(m,-khi,mu)-gi(n,khi,nu)*gj(m,-khi,mu)))/2.d0
          ! pour test   rappel: maintenant, c'est i0,j0 et j0,i0  luc86p66
          !IF(m==m0.and.n==n0.and.mu==mu0.and.nu==nu0) then
          !PRINT*, khi,phi
          !     s3j=symbol3j(m,n,l0,khi,-khi,0)
          !     x=s3j*phi
          !IF(m==n.and.nu==-mu) x=s3j*(phi+fi(n,khi,nu)*fj(m,-khi,mu)-gi(n,khi,nu)*gj(m,-khi,mu))/2.d0   NON, plus maintenant
          !     IF(((mu==0.and.nu==0).OR.(m==n.and.mu==-nu)).and.khi>0) x=2.d0*x     ! khi et -khi
          !phitest1=phitest1+SQRT((2.d0*m0+1.d0)*(2.d0*n0+1.d0))*x
          !endif
          !  fin test

          phi=rac(2*m+1)*rac(2*n+1)*phi

          if ( isnan(phi) ) then
              write(*,*) 'phi is NAN','ialp=',ialp,'imol=',i_mol,'jmol=',j_mol
              stop
 
          end if

          iacc(ialp,kr)=iacc(ialp,kr)+phi

          IF(ialp==ialp_aaa.and.kr==k_aaa) aaa=aaa+phi
          IF(ialp<=7) then
             do kk=1,k_bis_max
                IF(kr==k_bis(kk)) aaa_bis(10*ialp+kk)=aaa_bis(10*ialp+kk)+phi
             end do
          endif
       end do ! ialp


!
    END do                    ! fin j_mol
!
END do               ! fin i_mol
!
!PRINT*, 'phimnlmunu(ij) fute = ',phitest1
!                                              on normalise aaa pour que ce soit g
      IF(ichoixaaa>100) then
      den=n_h2o**2/2.*4.*pi/3.*((k_aaa*p)**3-((k_aaa-1)*p)**3)
      aaa=aaa/den
      endif
      do kk=1,k_bis_max
      den=n_h2o**2/2.*4.*pi/3.*((k_bis(kk)*p)**3-((k_bis(kk)-1)*p)**3)
      aaa_bis(10*(/1,2,3,4,5,6,7/)+kk)=aaa_bis(10*(/1,2,3,4,5,6,7/)+kk)/den
      end do
!
!
      aaa_bis(80)=aaa           ! je mets la variable aaa d'avant a la fin du tableau aaa_bis
!
      return
      end
!
!
function vector(r1,r2) !DOC
IMPLICIT REAL(8) (a-h,o-z)
DIMENSION r1(3),r2(3),vector(3)
!
!DOC        fait le produit vectoriel r = r1 * r2
!DOC        luc85 page 23
!
x3=r1(1); y3=r1(2); z3=r1(3)              ! r peut eventuellement remplacer en memoire r1
vector(1)=y3*r2(3)-z3*r2(2)
vector(2)=z3*r2(1)-x3*r2(3)
vector(3)=x3*r2(2)-y3*r2(1)
end
!
function harm_sph(m,mu,mup,beta) !DOC
implicit real(8) (a-h,o-z)
COMMON/facto/fac(0:101)
!
!DOC          pour harmonique spherique Rm,mu,mup(Omega=omega,beta,phi)
!DOC                                    =exp(-i*omega)*r(m)mu,mup*exp(-i*phi)
!DOC          calcul l'element mu,mup de la matrice r(m) en fonction de l'angle beta
!DOC          formule de Wigner dans Messiah eq.72 p922 betement
!DOC          luc72p143
!DOC          si mu ou mup nul, formule de recurrence stable
!
IF(ABS(mu)>m.or.ABS(mup)>m) THEN; harm_sph=0.; RETURN; endif
!
IF(mu==0.or.mup==0) THEN                      ! mu=0 ou mup=0
!
mu0=mu; beta0=beta
IF(mu==0) THEN; mu0=mup; beta0=-beta; ENDIF       ! je mets le 0 en second
!          si mu negatif, ca vaut (-1)**mu * valeur pour -mu
x=1.d0;
IF(mu0<0) THEN; x=(-1)**mu0; mu0=-mu0; endif
cc=COS(beta0)                            ! plutôt a partir d'une formule de recurrence stable
pm1=0.                                  ! des polynomes de Legendre associes Pl,m
pm=1.d0                                 ! luc73p96
do l=mu0+1,m
pm2=pm1
pm1=pm
pm=(cc*dble(2*l-1)*pm1-dble(l+mu0-1)*pm2)/DBLE(l-mu0)
end do
harm_sph=x*(-1)**mu0*SQRT(fac(m-mu0)/fac(m+mu0))*fac(2*mu0)/(2.d0**mu0*fac(mu0))*SIN(beta0)**mu0*pm
!
        else                  !   donc mu et mup non nuls, utiliser betement la formule de Wigner
!
harm_sph=0.
cc=COS(0.5d0*beta)
ss=SIN(0.5d0*beta)
do it=MAX(0,mu-mup),MIN(m+mu,m-mup)
harm_sph=harm_sph+(-1)**it/(fac(m+mu-it)*fac(m-mup-it)*fac(it)*fac(it-mu+mup))*  &
                           cc**(2*m+mu-mup-2*it)*ss**(2*it-mu+mup)
end do
harm_sph=SQRT(fac(m+mu)*fac(m-mu)*fac(m+mup)*fac(m-mup))*harm_sph
        endif
end function


!     Last change:  LB   16 Jul 2014    5:11 pm
!     Last change:  LB    2 Mar 2012    4:05 pm
      subroutine resul !DOC
!DOC get the results of accumulation (compute the projections)

      implicit REAL(8) (a-h,o-z)
      real * 8 nittot,naccep
      integer * 4 ncumul
      REAL * 8 iacc,iacc_ss,iacc_o
      PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O
      COMMON/melange/ncst_ion,ncst_o,ncst_h,ncrs,ij_tab(nomax,nomax),ialpdeb(ncrsmax),ialptot(ncrsmax)
      common/nombre/n,nmol(ncstmax),nion_tot,n_h2o,nmol_tot
      common/positi/xx(10000),yy(10000),zz(10000),i_espece(10000)        ! positions
      common/iteri/nittot,naccep,ncumul
      common/numeri/p,num
      common/energi/uu(10000),vrnew(10000),uutot
      COMMON/energies/uu_lj,uu_lj12,uu_lj6,uu_ew,uu_ew_r
      common/moyenn/uutotm,duutot,evtest,devtes,mmu
      COMMON/nb_proj/mnmax,ialpmax
      common/cumul_1/iacc_ss(ncstmax,2000),gr_ss(ncstmax,2000)     ! gsite-H : ion-H, OH et HH
      PARAMETER(ialpmaxx=(nomax-1)*nomax/2+(nomax-1)*16+1226)                    ! max pour nmax=6!
      PARAMETER(n_omega=8)
      COMMON/proj_l/mm(ialpmaxx),nn(ialpmaxx),ll(ialpmaxx),mumu(ialpmaxx),nunu(ialpmaxx)
      COMMON/proj_khi/mm1(ialpmaxx),nn1(ialpmaxx),khi1(ialpmaxx),mumu1(ialpmaxx),nunu1(ialpmaxx)
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)
      common/cumul_omega/iacc_o(n_omega,2000),gr_o(n_omega,2000),dcosbet,dphi,coeff_a,coeff_b           ! pour g(r,omega) luc85p176,190
      common/conc/conc_h2o,xl_a
      common/diamet/diam_a(ncstmax),sigma_a(ncrsmax),sigma(ncrsmax),sigma2(ncrsmax),rcutoff_a,rcut2,roh_a,theta_d
      COMMON/lennardjones/epsi_lj(ncstmax),xlj(ncrsmax),ucut
      COMMON/ewald/temp,dbjr,charge(ncstmax),alp_a1,alp,rmax2,kmax,kkmax,cdiel_ext
      COMMON/isobare/press_pa,press_a3,rho1,rho2,vol1,vol2,vol_sph_cub(2000)
      common/pression/vir,vir1,vir2,vir02,vir12,vir22,hypervir,hypervir2
      COMMON/cte_diel/xmtotx,xmtoty,xmtotz,xmtot2,xmtot22
!
!        resultats cumul‚s
!
!        luc belloni 4-4-98
!
!        version H2O
!        LJ total, non tronque! calcule avec Ewald  luc84p144
!
!        version isobare luc87p144
!
!        version electrolyte  luc87p188
!
    3 format(1x,78(1h*))
      write(*,3)
      write(*,*)'Nombre total de cycles par particule realises: ',nittot
      xn=n
      xn_h2o=n_h2o
!        pourcentage de deplacements acceptes
      x=naccep
      y=nittot
      alp1=x/(y*nmol_tot)
      write(*,*) 'Fraction de deplacements acceptes: ',alp1
      write(*,*) 'Nombre d''accumulations: ',ncumul
      write(*,3)
!        energie
PRINT*, 'Configuration actuelle:'
PRINT*, 'Concentration h2o (en M): ',conc_h2o
PRINT*, 'Taille de cellule L (en A): ',xl_a
PRINT*, 'Energies LJ12 et LJ6: ',uu_lj12,uu_lj6
      write(*,*) 'Energie EW : ',uu_ew
      write(*,*) 'Energie totale (kT) : ',uutot
!        moyenne des <Ri**2> et <ri**2> actuelles
      sum0=0.
      do i=1,n
      sum0=sum0+xx(i)**2+yy(i)**2+zz(i)**2
      end do

      write(*,*) 'sum0',sum0,'xn',xn
      sum0=SQRT(sum0/xn)
      WRITE(*,*) 'Actuellement, <ri> = ',sum0
      write(*,3)
!        fonction de distribution
!        on calcule g(r) a partir de iacc
!        il faut diviser par le nombre ideal
!        rho=N/L**3 = N puisque les longueurs sont en unite L
      pi=4.d0*ATAN(1.d0)
      x=ncumul
      xno=n_h2o
      xnh=nmol(ncst_h)
!
                      if(ncumul.ne.0) then
!
!           conc et volume
       conc_moy=rho1/ncumul                                                            ! <N/V>
       dconc=SQRT(MAX(0.,rho2/ncumul-conc_moy**2))
       write(*,18) 'Concentration H2O moyenne (M) : ',conc_moy,' +/- ',dconc
       vol_moy=vol1/ncumul
       WRITE(*,18) 'ou N/<V> (M) : ',xn_h2o/vol_moy/6.02214e-4                         ! N/<V>
       dvol2=vol2/ncumul-vol_moy**2
       WRITE(*,17) 'concentratios ioniques moyennes (M) : ',nmol(1:ncst_ion)*conc_moy/n_h2o
       compress_isobar=xn_h2o*dvol2/vol_moy**2
      WRITE(*,18) 'Compressibilite "isobare"= ',compress_isobar
      rhokt=conc_moy*1000.d0*6.02214d23*1.38065d-23*temp/1.d6       ! rho*kT en Mpa
      WRITE(*,18) 'soit en Mpa-1 = ',compress_isobar/rhokt
!           on normalise gOO,gOH,gHH
!           attention: en NPT, il faut diviser par <rho> et non rho_actuel
!           donc diviser les expressions NVT par V_actuel*<1/V>  luc87p145
!           et comme l'intersection de la calotte spherique avec le cube peut changer (a r>L/2),
!           on a accumule ce volume dans vol_sph_cub(:) luc87p152
!           en fait, il faudrait moyenner iacc/vol_sph_cub! luc87p153
avog=6.02214d-4
vol_inv_moy=conc_moy/xn_h2o*avog       ! <1/V> en A-3
fac_inv=xl_a**3*vol_inv_moy
        do k=1,num
!       xk=(k-1)*p
!       xk1=xk+p
!       vk=xk1**3-xk**3
!        IF(xk1>0.5d0) vk=(volume_sph_cub(xk1)-volume_sph_cub(xk))/(4.d0*pi/3.d0)
        gr_ss(1:ncst_h,k)=iacc_ss(1:ncst_h,k)/vol_sph_cub(k)/(x*xnh*nmol(1:ncst_h)*vol_inv_moy)
        gr_ss(ncst_h,k)=gr_ss(ncst_h,k)*2.d0
        END do
!           et idem pour les g(r,omega)  luc85p176,190
        do k=1,num
!        xk=(k-1)*p
!        xk1=xk+p
!        vk=xk1**3-xk**3
!        IF(xk1>0.5d0) vk=(volume_sph_cub(xk1)-volume_sph_cub(xk))/(4.d0*pi/3.d0)

        gr_o(:,k)=iacc_o(:,k)/vol_sph_cub(k)/(x*xno/2.d0*conc_moy*avog)/(dcosbet**2*dphi)*(8.d0*pi) /   &
            (/coeff_a,coeff_b,coeff_a,coeff_b,coeff_a,coeff_b,coeff_a,coeff_b/)   /fac_inv      ! luc85p190
        END do
!           on normalise gmnmunukhi et on transforme dans la foulee en gmnlmunu
call factoriel
gr=0.
do ies=1,ncst_o
 do jes=ies,ncst_o
 ij=ij_tab(ies,jes)

!   write(*,*) 'ies',ies,'jes',jes,'nmol(ies)',nmol(ies),'nmol(jes)',nmol(jes),'x'
    do ialp1=ialpdeb(ij),ialpdeb(ij)+ialptot(ij)-1                               ! toutes les proj en khipour ies,jes

       do k=1,num
!        xk=(k-1)*p
!        xk1=xk+p
!        vk=xk1**3-xk**3
!        IF(xk1>0.5d0) vk=(volume_sph_cub(xk1)-volume_sph_cub(xk))/(4.d0*pi/3.d0)
!        grint(k)=iacc(ialp1,k)/den_oo/vk  /fac_inv           ! sert de tableau auxi
          r = (k-1)*p_a
          dv = 4d0/3d0*pi*(3d0*r**2*p_a + 3d0*r*p_a**2 + p_a**3)
          grint(k)=iacc(ialp1,k)/dv/(ncumul*nmol(ies)*nmol(jes)*vol_inv_moy)           ! sert de tableau auxi
          IF(ies==jes) grint(k)=grint(k)*2.d0
       END do ! k

       do ialp=ialpdeb(ij),ialpdeb(ij)+ialptot(ij)-1                      ! je cherche les proj en l concernees

         IF(mm(ialp)==mm1(ialp1).and.nn(ialp)==nn1(ialp1).and.mumu(ialp)==mumu1(ialp1).and.nunu(ialp)==nunu1(ialp1)) then
            s3j=symbol3j(mm(ialp),nn(ialp),ll(ialp),khi1(ialp1),-khi1(ialp1),0)
            x=(2.d0*ll(ialp)+1.d0)*s3j

            IF(((mumu(ialp)==0.and.nunu(ialp)==0).OR.(ies==jes.and.mm(ialp)==nn(ialp).and.mumu(ialp)==-nunu(ialp)))        &
                              .and.khi1(ialp1)>0) x=2.d0*x     ! khi et -khi

            gr(ialp,1:num)=gr(ialp,1:num)+x*grint(1:num)

          endif
        end do ! ialp

    END do! ialp1
 end do
end do
!        grint(i) est le nombre de voisins dans les i premiers intervalles
!        vaut donc Nv(i*dr)
      den=x*xn_h2o/2.d0
      grint(1)=iacc(ialpdeb(ncrs),1)/den
        do 4 k=2,num
    4   grint(k)=grint(k-1)+iacc(ialpdeb(ncrs),k)/den
!
!        on calcule l'energie interne=energie micros moyenne
!        et la capacite calorifique=fluctuations en energie micro
!        attention ici, en NPT, Cp vient des fluctuations de U+PV
!        uutotm a accumule U alors que duutot a accumule (U+PV)**2
      u_moy=uutotm/ncumul
      cp=duutot/ncumul-(u_moy+press_a3*vol_moy)**2
   17 FORMAT(1x,a,5f12.6)
   18 FORMAT(1x,a,f14.6,a,f12.6)
      write(*,18) 'Energie interne U/kT : ',u_moy
      write(*,18) 'ou, U/NkT :',u_moy/nmol_tot
      xktkj=1.38065d-23*temp*6.02214d23/1000.d0         ! pour passer de kT/part a kJ/mol
      write(*,18) 'soit en kJ/mol :',u_moy/nmol_tot*xktkj
      write(*,18) 'Enthalpie H/kT : ',u_moy+press_a3*vol_moy,' +/- ',SQRT(cp)
      write(*,18) 'ou, H/NkT :',(u_moy+press_a3*vol_moy)/nmol_tot,' +/- ',SQRT(cp)/nmol_tot
      write(*,18) 'Capacite calorifique Cp/Nk :',cp/nmol_tot
      write(*,18) 'soit en kJ/mol/degre :',cp/nmol_tot*xktkj/temp
!       PAS de correction ici car LJ total!
!        pression du viriel
!        attention: vir contient maintenant N/V+W/V (et non 3W) en A-3
!        donc p/<rho>kT=vir/<rho>
      pi=4.d0*ATAN(1.d0)
      conc_tot_moy=conc_moy/xn_h2o*nmol_tot
      pvir=vir/ncumul  / (conc_tot_moy*avog)
      dpvir=SQRT(MAX(0.,vir02/ncumul-(vir/ncumul)**2)) / (conc_tot_moy*avog)
      write(*,18)'Pression viriel total P/(<rho_tot>kT) = ',pvir,' +/- ',dpvir
      rhokt=conc_tot_moy*1000.d0*6.02214d23*1.38065d-23*temp/1.d6       ! rho*kT en Mpa
      write(*,18)'soit en Mpa = ',pvir*rhokt,' +/- ',dpvir*rhokt
!        hyperviriel et compressibilite luc78p189
!        hypervir contient deja 9KHI/V
!        a priori, la formule NVT va donner khi=infini avec fluctuations NPT de P!  luc87p146
      hvir=1./9.*hypervir/ncumul
      dhvir=SQRT(MAX(0.,hypervir2/ncumul-(hypervir/ncumul)**2))/9.
      WRITE(*,18) 'Hyperviriel <X/V> = ',hvir,' +/- ',dhvir
      compress=1./(pvir+hvir/(conc_tot_moy*avog)-nmol_tot*dpvir**2)
      dcompress=compress**2*SQRT(dpvir**2+dhvir**2/(conc_tot_moy*avog)**2+nmol_tot**2*dpvir**4)         ! luc82p88 et corrige luc85p31!
      WRITE(*,18) 'Compressibilite ("infinie" en NPT!) = ',compress,' +/- ',dcompress
      WRITE(*,18) 'soit en Mpa-1 = ',compress/rhokt,' +/- ',dcompress/rhokt
!        on calcule le potentiel chimique ainsi que g(0) d'apres les fluctuations
!      xmu=-log(evtest/ncumul)
!      dxmu=SQRT(devtes*ncumul/evtest**2-1.d0)
!      write(*,18)'Potentiel chimique d''exces = ln(gamma) :',xmu,' +/- ',dxmu
!      WRITE(*,18) 'Apres correction, ln(gamma) = ',xmu+2.*utail+xmushift
      write(*,3)
!
!         constante dielectrique
!
!         soit par fluctuations du moment dipolaire total
!         xmtotx a moyenne M'x/racine(L)
!         xmtot2 a moyenne M'**2/L=M**2/V
      WRITE(*,17) 'Moment dipolaire total / racine(V) moyen : ',             &
        xmtotx/ncumul,xmtoty/ncumul,xmtotz/ncumul
      WRITE(*,17) 'Moment dipolaire quadratique /racine(V) moyen : ',SQRT(xmtot2/ncumul)
      xmu_h2o=2.d0*charge(ncst_h)*roh_a*COS(theta_d/180.d0*pi/2.d0)
      g1=(xmtot2/ncumul)/(xn_h2o*vol_inv_moy*xmu_h2o**2)
      dg1=SQRT(MAX(0.,xmtot22/ncumul-(xmtot2/ncumul)**2))/(xn_h2o*vol_inv_moy*xmu_h2o**2)              ! avec fluctuations luc82p88
      WRITE(*,18) 'Parametre g=<M**2/V>/<N/V>*Mu**2 : ',g1,' +/- ',dg1
      y=4.*pi/9.*xn_h2o*vol_inv_moy*xmu_h2o**2*(dbjr*xl_a)
      epsilon1=1.+1./(1./(3.*g1*y)-1./(2.*cdiel_ext+1.))
      depsilon1=(epsilon1-1.)**2*dg1/(3.*y*g1**2)
!        soit par integrale de g110 , 3eme proj ici
      g2=0.
      do k=1,num
!      dv=((k*p)**3-((k-1)*p)**3)/3.
      dv=(volume_sph_cub(k*p)-volume_sph_cub((k-1)*p))/(4.d0*pi)
      dv=vol_sph_cub(k)/ncumul/(4.d0*pi)                           ! j'utilise <V(k)>
      g2=g2+gr(ialpdeb(ncrs)+2,k)*dv
      end do
!      g2=1.d0-4.d0*pi/SQRT(3.d0)*xn_h2o*vol_inv_moy*xl_a**3*g2
      g2=1.d0-4.d0*pi/SQRT(3.d0)*xn_h2o*vol_inv_moy*g2
      WRITE(*,18) 'Parametre g par integrale de g110 : ',g2
!      epsilon2=1.+1./(1./(3.*g2*y)-1./(2.*cdiel_ext+1.))   Kirkwood  luc83p19
      epsilon2=1.+1./(1./(3.*g2*y)-1./(2.*cdiel_ext+1.))
!        ou j'utilise la formule de Caillol  luc80p201
      y2=9.*y*g2
      e=epsilon1         ! essai
      ep=cdiel_ext
!      v=4.*pi/3.*(k*p)**3
      v=volume_sph_cub(num*p)
        do l=1,20        ! 20 cycles NR
        fe=(e-1.)*(2.*e+1.)/e+2.*(e-1.)**2/e*(ep-e)/(2.*ep+e)*v
        e2=e**2
        fpe=(1./e2+2.)+2.*((1.-1./e2)*(ep-e)/(2.*ep+e)-               &
                           (e-1.)**2/e*3.*ep/(2.*ep+e)**2)*v
        e=e+(y2-fe)/fpe
        end do
      epsilon3=e
      WRITE(*,18) 'Constante dielectrique : ',epsilon1,' +/- ',depsilon1
      WRITE(*,18) 'ou : ',epsilon2,' ou ',epsilon3
           endif
      return
      end
!
    function volume_sph_cub(r) !DOC
      implicit REAL(8) (a-h,o-z)
!
!DOC       calcule le volume egal a l'intersection d'une sphere de rayon r et d'un cube de cote L
!DOC       tout en unite L
!DOC       formule de Caillol
!DOC       ou ma formule, luc82p197
!
      pi=4.d0*ATAN(1.d0)
      rl=r/0.5d0                ! tout en unite L/2=1/2
      rl2=rl**2; rl3=rl*rl2
      IF(rl<1.d0) then                          ! r<L/2
      volume_sph_cub=4.d0*pi/3.d0*rl3
      ELSEIF(rl<SQRT(2.d0)) then            ! L/2<r<L/2*racine(2)
      volume_sph_cub=-8.d0*pi/3.d0*rl3+6.d0*pi*rl2-2.d0*pi
      else                                     ! L/2*racine(2)<r<L/2*racine(3)
      xl1=SQRT(rl2-1.d0)
      xl2=SQRT(rl2-2.d0)
      v=1.d0-xl2+pi/4.d0*(4.d0/3.d0*rl3+3.d0*rl2-1.d0)+     &
      (1.d0-3.d0*rl2)*asin(1.d0/xl1)-rl3*ATAN(rl*xl2)
      volume_sph_cub=8.d0-8.d0*v
      endif
      volume_sph_cub=volume_sph_cub/8.d0
      end
!
      subroutine factoriel !DOC
      implicit REAL(8) (a-h,o-z)
      COMMON/facto/fac(0:101)
      fac(0)=1.d0
      do k=1,101
      fac(k)=dble(k)*fac(k-1)
      end do
      return
      end subroutine
!
      function itriangle(m,n,l) !DOC
!DOC        nul sauf si |m-n|<l<m+n
!DOC        rq: ne depend pas de l'ordre des 3 entiers
      itriangle=0
      IF(l.ge.ABS(m-n).and.l.le.m+n) itriangle=1
      end
!
      function delta(m,n,l) !DOC
      implicit REAL(8) (a-h,o-z)
      COMMON/facto/fac(0:101)
      delta=fac(m+n-l)*fac(n+l-m)*fac(l+m-n)/fac(m+n+l+1)
      end
!
      function symbol3j(m,n,l,mu,nu,lu) !DOC
      implicit REAL(8) (a-h,o-z)
      COMMON/facto/fac(0:101)
!
!DOC        symbole 3j
!DOC        Messiah page 910 eq.21
!
      IF(itriangle(m,n,l).eq.0.or.mu+nu+lu.NE.0.or.                                     &
         ABS(mu).gt.m.or.abs(nu).gt.n.or.abs(lu).gt.l) then
      symbol3j=0.
      else
      som=0.
      do it=MAX(0,n-l-mu,m-l+nu),MIN(m+n-l,m-mu,n+nu)
      som=som+(-1)**it/(fac(it)*fac(l-n+it+mu)*fac(l-m+it-nu)*                          &
                fac(m+n-l-it)*fac(m-it-mu)*fac(n-it+nu))
      end do
      symbol3j=(-1)**(m-n-lu)*SQRT(delta(m,n,l))*                                       &
     SQRT(fac(m+mu)*fac(m-mu)*fac(n+nu)*fac(n-nu)*fac(l+lu)*fac(l-lu))                 &
                 *som
      IF(mu==0.and.nu==0.and.lu==0.and.2*((m+n+l)/2)/=m+n+l) symbol3j=0.
      endif

end function
!

subroutine write_projections(hfile) !DOC
!DOC write projections to the file
      use io, only : io_open,io_close
      use SystemSettings, only: SYSTEM_STRING_LENGTH
      implicit none
 !DOC Parameters:
      integer :: hfile !DOC file handler

      real(8) :: p 
      integer :: num  
      common/numeri/p,num

      integer :: mnmax,ialpmax
      COMMON/nb_proj/mnmax,ialpmax

      integer(4) :: ncstmax,ncrsmax, nomax
      PARAMETER(ncstmax=4,ncrsmax=6,nomax=ncstmax-1)               ! composants +,-,O,H    croises ++,+-,+H2O,--,-H2O,H2Oh2O

      integer :: ialpmaxx   
      PARAMETER(ialpmaxx=(nomax-1)*nomax/2+(nomax-1)*16+1226)                    ! max pour nmax=6!
  
      real(8) :: iacc,gr,grint 
      common/cumul/iacc(ialpmaxx,2000),gr(ialpmaxx,2000),grint(2000)
  
      integer :: k,ialp

      integer :: hfile_acc, hfile_int
      character(SYSTEM_STRING_LENGTH) :: tmpstr

!      tmpstr = 'grint.dat'
!      hfile_int = io_open(tmpstr,'w')
!      tmpstr = 'iacc.dat'
!      hfile_acc = io_open(tmpstr,'w')

      do k=1,num

         ! k = int(r/p) + 1 --> k-1 = int(r/p), r = p(k-1)  
         write(hfile,'(F10.6,$)')  p_a*(k-1)

         do ialp = 1,ialpmax

            write(hfile,'(F20.10,$)')  gr(ialp,k)
   !         write(hfile_int,'(F20.10,$)')  grint(k)
   !         write(hfile_acc,'(F20.10,$)')  iacc(ialp,k)
         end do 

         do ialp = 1,ialpmax
            write(hfile,'(F20.10,$)') iacc(ialp,k)
         end do

         write(hfile,*)
  !       write(hfile_int,*)
  !       write(hfile_acc,*)
      end do

   !   call io_close(hfile_int)
   !   call io_close(hfile_acc)
     

end subroutine




End Module MCAccumLuc
