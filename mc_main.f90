program runmc_main !DOC !PROG
!DOC !PROG !FILE Main Program: runmc
use RhoSquared
use SumSinCosKR
use parameters
use MoleculeTable
use AtomicData
use LJTypes
use MC , only : initmc, movemc,clearmc,RealSpace_TotalSum,KSpace_TotalSum,mc_calc_total_sums
use MC , only : U_tail_lj6,U_tail_lj12,U_tail_coulomb, U_total_ew,U_total_lj,U_total,U_ext
use io, only : write_integer_array,write_real_array,write_xyz_array,io_open,io_close,write_real_matrix
!use MCLuc, only : initmc_Luc,movemc_Luc,MCLuc_init_input,MCLucInput,piston ,write_energy_LUC
use runmc, only : run_mc
use MRandom


!DOC !PROG The program uses the input files which describe the system to produce the trajectories of the particles.

!DOC !PROG 

!DOC !PROG Usage:      runmc  parameters.prm system.composition  input.moltab [ extra_parameters ]
!DOC !PROG Arguments:
!DOC !PROG ::    parameters.prm - file, which include main parameters of the simulation
!DOC !PROG ::    system.composition - file, which describes the numbers of molecules of each kind in the system
!DOC !PROG ::    input.moltab   - initial configuration of the system
!DOC !PROG ::    extra_parameters  - coma separated list of   pair   parameter1=value1,parameter2=value2

implicit none

character(SYSTEM_STRING_LENGTH) :: parameters_file,composition_file,moltab_file
character(SYSTEM_STRING_LENGTH) :: tmpstr
character(SYSTEM_STRING_LENGTH) :: extra_parameters_string
logical :: extra_parameters = .FALSE.

Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TAtomicData) :: atomic_data
Type(TLJTypes) :: lj_types

integer :: i,j
integer :: hfile,nk
integer :: MaxMolecule,MaxAtom
integer :: imol

logical :: accepted

real(8) :: rnd
real(8) :: sigma2,xlj
real(8),dimension(:),pointer :: pa,pb,pc
Type(TSumSinCosKR),pointer :: Coul,LJ
real(8) :: uu_lj

real(8),dimension(:),allocatable :: xx,yy,zz ! for test movemc

call parse_command_line
call read_input_files


call randomize(rnd_seed)

if( param_out_file(1:5) /= 'none ' ) then
   hfile = io_open( param_out_file ,'w')
   call parameters_write(hfile)
   call io_close(hfile)
end if

 
!call test_initmc
call initmc(comp, mol_tab, atomic_data, lj_types )


!write(*,*) '************** BEROFE ALL ******************'
!call write_energy_ME

!call movemc(

!accepted = .FALSE.
!imol = 1

!do while ( .not. accepted )
       
!   call movemc( imol, mol_tab,  accepted )
!   write(*,*) 'A:',accepted
!end do

!write(*,*) '***************************  AFTER THE MOVE U+dU *******'

!call write_energy_ME


!call mc_calc_total_sums(mol_tab,atomic_data,lj_types)

!call initmc(comp,mol_tab,atomic_data, lj_types)

!write(*,*) '**************************** AFTER RECALC **************'

!call write_energy_ME

!tmpstr = 'sumsincos_total_RECALC.rho2'
!call RhoSquared_save_to_file( KSpace_TotalSum % rho_squared_total, tmpstr )

!tmpstr = 'beta_LJ.dat'
!hfile = io_open(tmpstr,'w')
!call write_real_matrix(hfile, KSpace_TotalSum % beta_LJ, KSpace_TotalSum % grid % nk, lj_types % NType ** 2)
!call io_close(hfile)

call run_mc(mol_tab,max_cycle)

!call test_mcvol
  

!   write(*,*) accepted
!end do

!call clearmc

call deallocate_all


contains 


subroutine write_energy_ME !DOC 
!DOC write the energy components (DEBUG ONLY)

write(*,*) '***************  ME **************** '
uu_lj = RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6 
write(*,*) 'U_REAL ew:',RealSpace_TotalSum % uu_ew,'lj:',uu_lj,'lj6:',&
                        RealSpace_TotalSum % uu_lj6,'lj12:',RealSpace_TotalSum % uu_lj12
write(*,*) 'U_FOUR ew:',KSpace_TotalSum % energy_coulomb,'lj:',KSpace_TotalSum % energy_LJ
write(*,*) 'U_TAIL ew:',U_tail_coulomb,'lj:',U_tail_lj12 - U_tail_lj6,'lj6:',U_tail_lj6,'lj12:',U_tail_lj12
write(*,*) 'U_EXT  ew:',U_ext
write(*,*) 'U_TOTL ew:',U_total_ew,'lj:',U_total_lj
write(*,*) 'U_TOT tot:',U_total

end subroutine


!subroutine test_initmc
!!Type(MCLucInput) :: input
!!write(*,*) 'MaxAtom=',MaxAtom
!
!!call write_real_array(0,atomic_data % epsilon,MaxAtom)
!
!!call write_real_array(0,lj_types % four_epsilon_tab,1)
!
!call initmc(comp, mol_tab, atomic_data, lj_types )
!
!!write(*,*) 'U_REAL=',SUM(RealSpace_TotalSum % uu)
!!write(*,*) 'F_REAL=',SUM(RealSpace_TotalSum % Fx),SUM(RealSpace_TotalSum % Fy),SUM(RealSpace_TotalSum % Fz)
!call write_energy_ME
!
!sigma2 = lj_types % sigma2_tab(1)
!xlj = lj_types % four_epsilon_tab(1)
!
!!write(*,*) 'sigma2=',sigma2,'xlj=',xlj
!
!call MCLuc_init_input(input)
!call initmc_Luc(input,atomic_data % xx,atomic_data % yy,atomic_data % zz,MaxMolecule,.FALSE.)
!
!
!!write(*,*) 'uu_lj6=',RealSpace_TotalSum % uu_lj6,'uu_lj12 =',RealSpace_TotalSum % uu_lj12,&
!!            'uu_lj=',RealSpace_TotalSum % uu_lj12 - RealSpace_TotalSum % uu_lj6, &
!!             'uu_ew = ',RealSpace_TotalSum % uu_ew
!
!nk = KSpace_TotalSum % grid % nk
!Coul => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
!
!tmpstr = 'sumcos_ppp.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, Coul % sumcos_ppp, nk)
!call io_close(hfile)
!
!tmpstr = 'sumsin_ppp.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, Coul % sumsin_ppp, nk)
!call io_close(hfile)
!
!tmpstr = 'beta.txt'
!hfile=io_open(tmpstr,'w')
!call write_real_array(hfile,KSpace_TotalSum % beta,nk)
!call io_close(hfile)
!
!LJ => KSpace_TotalSum % rho_squared_total % sumsincos_LJ(1)
!
!tmpstr = 'sumcos_LJ_ppp.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, LJ % sumcos_ppp, nk )
!call io_close(hfile)
!
!tmpstr = 'beta6.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, KSpace_TotalSum % beta6, nk)
!call io_close(hfile)
!
!tmpstr = 'beta12.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, KSpace_TotalSum % beta12, nk)
!call io_close(hfile)
!
!tmpstr = 'betaLJ.txt'
!hfile = io_open(tmpstr,'w')
!call write_real_array(hfile, KSpace_TotalSum % beta_LJ(1) % ptr, nk)
!call io_close(hfile)
!
!
!
!
!
!!do i=1,MaxAtom
!
! !        write(tmpstr,'(AI3.3A)') 'sumcosC',i,'.txt'
!  !       hfile = io_open(tmpstr,'w')
!
!!         pa => KSpace_TotalSum % rho_squared_total % 
!         
! !        call write_xyz_array(hfile,cckx,ccky,cckz,kmax)
!  !       call io_close(hfile)
!
!   !      write(tmpstr,'(AI3.3A)') 'sumsinC',i,'.txt'
!   !      hfile = io_open(tmpstr,'w')
!   !      call write_xyz_array(hfile,sskx,ssky,sskz,kmax)
!   !      call io_close(hfile)
!!end do
!
!!write(*,*) 'uu='
!!call write_real_array(0,RealSpace_TotalSum % uu, MaxAtom )
!
!!write(*,*) 'uu_back='
!!call write_real_array(0,RealSpace_TotalSum % uu_back, MaxAtom )
!
!
!!write(*,*) 'Fx='
!!call write_xyz_array(0,RealSpace_TotalSum % Fx, RealSpace_TotalSum % Fy, RealSpace_TotalSum % Fz, MaxAtom)
!
!!tmpstr = 'xk2.txt'
!!hfile = io_open(tmpstr,'w')
!!call write_real_array(hfile,KSpace_TotalSum % grid % xk2, KSpace_TotalSum % grid % nk )
!!call io_close(hfile)
!
!!write(*,*) 'Fx_back='
!!call write_xyz_array(0,RealSpace_TotalSum % Fx_back, RealSpace_TotalSum % Fy_back, &
!!                       RealSpace_TotalSum % Fz_back, MaxAtom)
!
!
!!do i=1,100
!
!
!
!
!end subroutine
!
!
!subroutine test_mcvol
!   use MC, only : RealSpace_TotalSum, KSpace_TotalSum,mc_total_energy,mcvol 
!   use RhoSquared 
!   use SumSinCosKR
!   use parameters, only : log_max_volume_scaling,pressure_angstr_kT,BoxLength
!   
!
!   Type(TSumSinCosKR),pointer :: sscew,ssclj
!
!   integer :: i_h2o
!   integer :: iaccep
!   real(8) :: dx,dy,dz,S
! 
!   integer :: nk
!
!   real(8) uu,vrnew,uutot
!
!   real(8) :: conc_h2o,xl_a
!
!common/conc/conc_h2o,xl_a
!common/energi/uu(10000),vrnew(10000),uutot
!
!  real(8) fx,fy,fz
!COMMON/forces/fx(10000),fy(10000),fz(10000)
!
!  real(8),parameter :: charg_h = 0.4238d0
!!  PARAMETER(nkmax=4661,kmaxmax=20)
! integer, parameter :: nkmax = 4661, kmaxmax = 20
!
!
! real(8) bet,                     &         ! pour Ewald en k
! sumcos,sumsin,sumcos1,sumsin1,        &
! sumcos2,sumsin2,sumcos3,sumsin3,      &
! cckx,sskx,ccky,ssky,    &
! cckz,sskz
! integer ip,iq
!
!COMMON/sommeenk/bet(nkmax),ip(nkmax),iq(nkmax),                     &         ! pour Ewald en k
! sumcos(nkmax),sumsin(nkmax),sumcos1(nkmax),sumsin1(nkmax),        &
! sumcos2(nkmax),sumsin2(nkmax),sumcos3(nkmax),sumsin3(nkmax),      &
! cckx(0:kmaxmax),sskx(0:kmaxmax),ccky(0:kmaxmax),ssky(0:kmaxmax),    &
! cckz(0:kmaxmax),sskz(0:kmaxmax)
!
!
!real(8) bet6,bet12,                     &         ! pour Ewald LJ en k
! sumcos_o,sumsin_o,sumcos1_o,sumsin1_o,        &
! sumcos2_o,sumsin2_o,sumcos3_o,sumsin3_o
!
!
!COMMON/sommeenk_lj/bet6(nkmax),bet12(nkmax),                     &         ! pour Ewald LJ en k
! sumcos_o(nkmax),sumsin_o(nkmax),sumcos1_o(nkmax),sumsin1_o(nkmax),        &
! sumcos2_o(nkmax),sumsin2_o(nkmax),sumcos3_o(nkmax),sumsin3_o(nkmax)
!
!   allocate(xx(MaxAtom))
!   allocate(yy(MaxAtom))
!   allocate(zz(MaxAtom))
!
!   xx(:) = atomic_data % xx(:)
!   yy(:) = atomic_data % yy(:)
!   zz(:) = atomic_data % zz(:)
! 
!   tmpstr = 'my_xyz_before.txt'
!   hfile = io_open(tmpstr,'w')
!   call write_xyz_array(hfile,atomic_data % xx, atomic_data % yy, atomic_data % zz, MaxAtom )
!   call io_close(hfile)
!
!
!   tmpstr = 'luc_xyz_before.txt'
!   hfile = io_open(tmpstr,'w')
!   call write_xyz_array(hfile, xx, yy,  zz, MaxAtom )
!   call io_close(hfile)
!
! 
!
!   do i=1,115
!
!      ! **************** VOLUME *******************
!
!
! !     write(*,*) '************ MCVOL : ME ****************'
! !     write(*,*) 'BoxLength:',BoxLength,'Energy:',mc_total_energy()
!      call randomize(i)
!      call mcvol( mol_tab, accepted )
!  !    write(*,*) 'Accepted:',accepted,'BoxLength:',BoxLength,'Energy:',mc_total_energy()
!
!!      write(*,*) '******* BEFORE CALC TOTAL SUMS *********'
!!      call write_energy_ME 
!!      call mc_calc_total_sums_dbg(mol_tab,atomic_data,lj_types)
!!     write(*,*) '******* AFTER CALC_TOTAL_SUMS *********'
!!      call write_energy_ME
!
!
! !     write(*,*) '********** MCVOL : LUC *****************'
! !     write(*,*) 'xl_a:',xl_a,'uutot:',uutot
!      call randomize(i)
!      call piston(iaccep,pressure_angstr_kT,log_max_volume_scaling,xx,yy,zz,MaxMolecule,.FALSE.)
!   !   write(*,*) 'iaccep:',iaccep,'xl_a:',xl_a,'uutot:',uutot
!
!    !  write(*,*)
!    !  write(*,*) '-----------------------------------------------------------------------'
!    !  write(*,*) 
!
!!      call write_energy_ME
!!      call write_energy_LUC(MaxMolecule)
!
!
! !     exit 
!
!
!      call randomize(i)
!      i_h2o = random(MaxMolecule) + 1
!
!      call randomize(i)
!      call movemc(i_h2o, mol_tab, accepted )
!   
!      
!      call randomize(i) 
!     call movemc_Luc(xx,yy,zz,MaxMolecule, i_h2o,iaccep)
!   
!!      call movemc_Luc(xx,yy,zz,MaxMolecule,sigma2,xlj , i_h2o,iaccep)
!
!!      write(*,*) 'i_h20:',i_h2o,'accepted:',accepted,'iaccep',iaccep
!
!      S = 0
!      do j=1,MaxAtom
!
!        dx = atomic_data % xx(j) - xx(j)
!        dy = atomic_data % yy(j) - yy(j)
!        dz = atomic_data % zz(j) - zz(j)
!
!        if( dx > 0.5 )   dx = dx - 1
!        if( dx < -0.5 )  dx = dx + 1
!
!        if( dy > 0.5 )   dy = dy - 1
!        if( dy < -0.5 )  dy = dy + 1
!
!        if( dz > 0.5 )   dz = dz - 1
!        if( dz < -0.5 )  dz = dz + 1
!
!
!         S = S + abs(dx) + abs(dy) + abs(dz) 
!      end do
!
!      write(*,*) 'Step',i,'imol',i_h2o,'A:',accepted,':',iaccep,'S=',S,'L=',BoxLength
! 
!      if ( S < 1.0d-8 ) then
!
!          nk = KSpace_TotalSum % grid % nk
!
!          sscew => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
!          ssclj => KSpace_TotalSum % rho_squared_total % sumsincos_LJ(1)
!
!          xx(:) = atomic_data % xx(:)
!          yy(:) = atomic_data % yy(:)
!          zz(:) = atomic_data % zz(:) 
!
!          fx(1:MaxAtom) = RealSpace_TotalSum % Fx(1:MaxAtom)
!          fy(1:MaxAtom) = RealSpace_TotalSum % Fy(1:MaxAtom)
!          fz(1:MaxAtom) = RealSpace_TotalSum % Fz(1:MaxAtom)
!
!          uu(1:MaxAtom) = RealSpace_TotalSum % uu(1:MaxAtom)
!
!
!!          write(*,*) 'sumcos1'
!!          call write_real_array(0,sumcos1*charg_h,5)
!!
!!          write(*,*) 'sumcos2'
!!          call write_real_array(0,sumcos2,5)
!!
!!          write(*,*) 'C_ppm'
!!          call write_real_array(0,sscew % sumcos_ppm,5) 
!!
!!          write(*,*) 'C_pmp'
!!          call write_real_array(0,sscew % sumcos_pmp,5) 
!
!          sumcos(1:nk) = sscew % sumcos_ppp(1:nk) / charg_h
!          sumsin(1:nk) = sscew % sumsin_ppp(1:nk) / charg_h
!      
!          sumcos1(1:nk) = sscew % sumcos_ppm(1:nk) / charg_h
!          sumsin1(1:nk) = sscew % sumsin_ppm(1:nk) / charg_h
!
!
!          sumcos2(1:nk) = sscew % sumcos_pmp(1:nk) / charg_h
!          sumsin2(1:nk) = sscew % sumsin_pmp(1:nk) / charg_h
!          
!          sumcos3(1:nk) = sscew % sumcos_pmm(1:nk) / charg_h
!          sumsin3(1:nk) = sscew % sumsin_pmm(1:nk) / charg_h
! 
!          sumcos_o(1:nk) =  ssclj % sumcos_ppp(1:nk)
!          sumsin_o(1:nk) =  ssclj % sumsin_ppp(1:nk)
!      
!          sumcos1_o(1:nk) = ssclj % sumcos_ppm(1:nk)
!          sumsin1_o(1:nk) = ssclj % sumsin_ppm(1:nk)
!
!          sumcos2_o(1:nk) = ssclj % sumcos_pmp(1:nk)
!          sumsin2_o(1:nk) = ssclj % sumsin_pmp(1:nk) 
!          
!          sumcos3_o(1:nk) = ssclj % sumcos_pmm(1:nk)
!          sumsin3_o(1:nk) = ssclj % sumsin_pmm(1:nk)
!        
!!          sumcos_o,sumsin_o,sumcos1_o,sumsin1_o,        &
!!          sumcos2_o,sumsin2_o,sumcos3_o,sumsin3_o
!
!
!
!      end if
!
!!      write(*,*) 'After the move:=='
!!
!!      call write_energy_ME 
!!      call write_energy_LUC(MaxMolecule)
!
!!      write(*,*) '==> </After>'
! 
!   end do 
!
!
!   tmpstr = 'my_xyz.txt'
!   hfile = io_open(tmpstr,'w')
!   call write_xyz_array(hfile,atomic_data % xx, atomic_data % yy, atomic_data % zz, MaxAtom )
!   call io_close(hfile)
!
!
!   tmpstr = 'luc_xyz.txt'
!   hfile = io_open(tmpstr,'w')
!   call write_xyz_array(hfile, xx, yy,  zz, MaxAtom )
!   call io_close(hfile)
!
!
!
!end subroutine


subroutine parse_command_line !DOC
!DOC read command line arguments

if ( iargc() < 3 ) then

   write(*,*)  'Usage:  runmc  parameters.prm system.composition  input.moltab [ extra_parameters ]'
   write(*,*)  '   system.composition - number and types of molecules in the system'
   write(*,*)  '   input.moltab - in binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   extra_parameters - coma separated  prm1=val1,prm2=val2,...   NO SPACES ALLOWED! '

   stop
end if

call getarg(1,parameters_file)
call getarg(2,composition_file)
call getarg(3,moltab_file)


if ( iargc() .ge. 4 ) then
    call getarg(4,extra_parameters_string)
    extra_parameters = .TRUE.
else
    extra_parameters = .FALSE.
end if



end subroutine


subroutine read_input_files !DOC
!DOC read input files

!write(*,*)  'input=',input_file(1:20), 'output=',output_file(1:20), 'nbytes_xyz=',nbytes_xyz, 'nbytes_ang=',nbytes_ang

! read the system composition
call Composition_nulify(comp)
call Composition_read_from_file(comp, composition_file )

MaxAtom = Composition_count_atoms(comp)
MaxMolecule = SUM( comp % mol_numbers(1:comp % n_types) )

! read and initialize the  parameters

if ( .not.extra_parameters ) then
     call Parameters_init(parameters_file,MaxMolecule)  ! MAxMolecule is used for the BoxLength calculation (from density, if required)
else
     call Parameters_init(parameters_file,MaxMolecule,extra_parameters_string)
end if


! read the molecule table
call MoleculeTable_nulify(mol_tab)

hfile = io_open(moltab_file,'b')
call MoleculeTable_load_binary( mol_tab, hfile, comp, input_nbytes_xyz, input_nbytes_ang )
call io_close(hfile)

! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )

end subroutine


subroutine deallocate_all !DOC 
!DOC deallocate everything

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)

call LJTypes_dealloc(lj_types)


end subroutine

end program runmc_main
