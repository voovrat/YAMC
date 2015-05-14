Module Parameters !DOC
use constants, only : boltz, elec, epsi0,kcal_mol
use SystemSettings, only : SYSTEM_STRING_LENGTH

!DOC !FILE parameters of the simulation
!DOC  everything which is independent of anything is in constants
!DOC everything else is here (e.g. which depends on kT, box Length etc.)

  implicit none

  logical :: Parameters_initialized = .FALSE. !DOC indicates that the parameters was initialized. Consider also Parameters_areConsistent(AtomicData)
                                              ! because some constants are BoxLengthDependent


! Parameters which have to be read from file/ keyboard
  logical :: Parameters_areRead = .FALSE. !DOC indicates that the parameters are read from the parameters file
   
  real(8) :: sr = 4.d0 !DOC rmax = sr / alpha
                !DOC err_r =  Q (sr/ alpha L^3)^(1/2) exp(-sr^2) / sr^2   
                !DOC  typical value sr = 3..4

  real(8) :: sk = 4.d0 !DOC kmax = s(L*alpha) / pi
                !DOC err_k = Q (sk/2 alpha L^3)^(1/2) exp(-sk^2) / sk^2
                !DOC typical value sk = 3..4
 
!  real(8) :: rmax = 0.5  ! to be read from input file   
!   rmax = sr / alpha
!  real(8) :: kmax
!  integer :: kmax = 4
  real(8) :: alpha = 0.d0  !DOC alpha * L  , by default  2*sr 
                        !DOC rmax =  sr / alpha < L/2 ==>  alpha*L > 2*sr
                        !DOC  alpha*L = alpha [L^-1] = 2*sr corresponds to the cutoff L/2 (default, if alpha=0)

!  real(8) :: alpha_angstr_m1 = 1  ! alpha [Angstr^-1] 
  real(8) :: temp =298.15d0 !DOC temperature, K
  real(8) :: external_permutivity = 50.d0  !DOC epsilon_ext
  real(8) :: BoxLength = 0.d0  !DOC if box_length = 0 --> to be recalculated from density
  real(8) :: density = 33.33d0 !DOC /particles/nm^3 
  real(8) :: dr_a = 0.4d0  !DOC max. dispacement in Angstroems  
  real(8) :: d_angle_degree = 40.d0 !DOC max rotation  
  real(8) :: xlambda_c =0.5d0  !DOC coupling constant for torque ( note: normally xlambda_f = xlambda_c AND xlambda_f = 0.5 )
  real(8) :: xlambda_f = 0.5d0  !DOC coupling constant for  force

  integer :: rnd_seed = 0 !DOC random seed 
  

  ! simulation
  real(8) :: pressure   !DOC in Pascals, i.e. J/m^3 

  real(8) :: pressure_step_multiplier = 0.d0  !DOC make NPT step at average each nmol * pressure_step_multiplier steps
                                           !DOC set 0 for NVT simulation
  real(8) :: max_volume_scaling = 1.05d0  !DOC  vnew = vold * max_volume_scaling ** lambda, -1 < lambda < 1


  integer :: max_cycle = 1000000   !DOC number of cycles to do 
  
  real(8) :: n_store_traj_interval = 100.d0   !DOC save trajectory frames each nmol * n_store_traj_interval steps 
                                           !DOC ( each n_store_traj_interval_cycles)
 
  real(8) :: n_store_energy_interval = 100.d0     !DOC save energy interval
 
 ! io
  integer :: input_nbytes_xyz=2, input_nbytes_ang = 2 !DOC  number of bytes to store values in input moltab file
  integer :: output_nbytes_xyz=2, output_nbytes_ang = 2 !DOC number of bytes to store values in the trajectory files

  character(SYSTEM_STRING_LENGTH) :: traj_file = 'traj.moltab' !DOC traj output file
  character(SYSTEM_STRING_LENGTH) :: energy_file = 'energy.dat'!DOC energy output file
  character(SYSTEM_STRING_LENGTH) :: boxlength_file = 'boxlength.dat' !DOC boxlength output file
  character(SYSTEM_STRING_LENGTH) :: frames_file = 'frames.dat' !DOC frames output file
  character(SYSTEM_STRING_LENGTH) :: param_out_file = 'parameters.out' !DOC parameters output file
 
  character(SYSTEM_STRING_LENGTH) :: freq_file = 'none' !DOC frequences input file
 
! parameters which can be re-calculated
  real(8) :: rmax         !DOC  rmax = sr / alpha
  integer :: kmax         !DOC  kmax = sk*(alpha L ) / pi
!   real(8) :: alpha        ! in BoxLength^-1
  real(8) :: kT_kcal_mol  !DOC boltz * temp

  real(8) :: rmax2 !DOC rmax^2

  real(8) :: dbjr_a   !DOC bjerum length in Angtroem:  elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10  
  real(8) :: dbjr   !DOC dbjr/BoxLength
                    ! note: we use kT as energy units, and positron charge as charges, and BoxLengths as distance units
                    ! e.g. U_C [kT] = dbjr * q_i q_j / r 
  real(8) :: xclb   !DOC =dbjr. coulomb potential prefactor. Just more clear name.
  real(8) :: alpha_over_sqrt_pi      !DOC alpha/sqrt(pi), use in EwaldSumTails
  real(8) :: two_alpha_over_sqrt_pi  !DOC 2alpha/sqrt(pi), use in real coulomb forces

  real(8) :: kext  !DOC  for external sum:  2 pi / (2 external_permutivity + 1)
  real(8) :: minus_two_kext !DOC - 2 * kext = -4pi/(2 external_permutivity + 1), used in external sum forces  
 
  real(8) :: dr  !DOC maximal displacement = dr_a / BoxLength 
  real(8) :: d_angle !DOC maximal rotation in radians
  
  real(8) :: log_max_volume_scaling !DOC log(max_volume_scaling), used in mcvol  
  real(8) :: pressure_angstr_kT  !DOC pressure in kT/A^3


contains 

subroutine read_parameter(nam,val) !DOC
!DOC initialize the parameter using the text pair nam and val
    use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC Parameters:
    character(SYSTEM_STRING_LENGTH),intent(in) :: nam,val !DOC nam and val strings read from the input file as nam=val

   ! if( nam == 'kmax ' ) read(val,*) kmax
   ! if( nam == 'rmax ' ) read(val,*) rmax
    if( nam == 'sr' ) read(val,*) sr
    if( nam == 'sk' ) read(val,*) sk
    if( nam == 'alphaL ' ) read(val,*) alpha
    if( nam == 'temp ' ) read(val,*) temp
    if( nam == 'external_permutivity ' ) read(val,*) external_permutivity  ! epsilon_ext
    if( nam == 'BoxLength ') read(val,*) BoxLength   ! if box_length = 0 --> to be recalculated from density
    if( nam == 'density ')  read(val,*) density ! particles/nm^3 
    if( nam == 'dr_a  ' )  read(val,*) dr_a  ! max. dispacement in Angstroems  
    if( nam == 'd_angle_degree ') read(val,*) d_angle_degree ! max rotation  
    if( nam == 'xlambda_c ') read(val,*) xlambda_c  ! coupling constant for torque ( note: normally xlambda_f = xlambda_c AND xlambda_f = 0.5 )
    if( nam == 'xlambda_f ') read(val,*)  xlambda_f  ! coupling constant for  force

    if( nam == 'input_nbytes_xyz ') read(val,*) input_nbytes_xyz
    if( nam == 'input_nbytes_ang ') read(val,*) input_nbytes_ang

    if( nam == 'output_nbytes_xyz ') read(val,*) output_nbytes_xyz
    if( nam == 'output_nbytes_ang ') read(val,*) output_nbytes_ang

    if( nam == 'pressure_step_multiplier ' ) read(val,*) pressure_step_multiplier
    if( nam == 'max_volume_scaling ' ) read(val,*) max_volume_scaling
    if( nam == 'n_store_traj_interval ') read(val,*) n_store_traj_interval
    if( nam == 'n_store_energy_interval ') read(val,*) n_store_energy_interval
    if( nam == 'max_cycle ') read(val,*) max_cycle
  
    if( nam == 'pressure ') read(val,*) pressure

    if( nam == 'traj_file ') traj_file = val
    if( nam == 'energy_file ') energy_file = val    
    if( nam == 'boxlength_file ') boxlength_file = val
    if( nam == 'frames_file ') frames_file = val

    if( nam == 'freq_file ') freq_file = val
    if( nam == 'param_out_file ') param_out_file = val


    if( nam == 'random_seed ') read(val,*) rnd_seed
 

end subroutine


subroutine parameters_write(h) !DOC 
!DOC write the parameters in the text format to the file 
  use string, only : str_get_next_pos
!DOC Parameters:
  integer,intent(in) :: h !DOC file handler
   ! NOTE: YOU CAN USE STDIN STDOUT AND STDERR from SystemSettings
  write(h,*)  'sr =', sr 
  write(h,*)  'sk =', sk
  write(h,*)  'alphaL = ',alpha   
  write(h,*)  'temp = ',temp 
  write(h,*)  'external_permutivity = ',external_permutivity 
  write(h,*)  'BoxLength = ',BoxLength 
  write(h,*)  'density = ',density  
  write(h,*)  'dr_a = ',dr_a 
  write(h,*)  'd_angle_degree = ',d_angle_degree
  write(h,*)  'xlambda_c = ',xlambda_c 
  write(h,*)  'xlambda_f = ',xlambda_f
  write(h,*)  'pressure = ',pressure   
  write(h,*)  'pressure_step_multiplier = ',pressure_step_multiplier
  write(h,*)  'max_volume_scaling = ',max_volume_scaling 
  write(h,*)  'max_cycle = ',max_cycle  
  write(h,*)  'n_store_traj_interval = ',n_store_traj_interval  
  write(h,*)  'n_store_energy_interval = ',n_store_energy_interval
  write(h,*)  'input_nbytes_xyz = ', input_nbytes_xyz
  write(h,*)  'input_nbytes_ang = ', input_nbytes_ang 
  write(h,*)  'output_nbytes_xyz = ',output_nbytes_xyz
  write(h,*)  'output_nbytes_ang = ',output_nbytes_ang
  write(h,*)  'rmax = ',rmax         !  rmax = sr / alpha
  write(h,*)  'kmax = ',kmax         !  kmax = sk*(alpha L ) / pi
  write(h,*)  'kT_kcal_mol = ',kT_kcal_mol  ! boltz * temp
  write(h,*)  'rmax2 = ',rmax2 ! rmax^2
  write(h,*)  'dbjr_a = ',dbjr_a   ! bjerum length in Angtroem:  elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10  
  write(h,*)  'dbjr = ',dbjr   ! dbjr/BoxLength
  write(h,*)  'xclb = ',xclb   ! =dbjr. coulomb potential prefactor. Just more clear name.
  write(h,*)  'alpha_over_sqrt_pi = ',alpha_over_sqrt_pi      ! alpha/sqrt(pi), use in EwaldSumTails
  write(h,*)  'two_alpha_over_sqrt_pi = ',two_alpha_over_sqrt_pi  ! 2alpha/sqrt(pi), use in real coulomb forces
  write(h,*)  'kext = ',kext  !  for external sum:  2 pi / (2 external_permutivity + 1)
  write(h,*)  'minus_two_kext = ',minus_two_kext ! - 2 * kext = -4pi/(2 external_permutivity + 1), used in external sum forces  
  write(h,*)  'dr = ',dr  ! dr_a / BoxLength 
  write(h,*)  'd_angle = ',d_angle
  write(h,*)  'log_max_volume_scaling = ',log_max_volume_scaling ! log(max_volume_scaling), used in mcvol  
  write(h,*)  'pressure_angstr_kT = ',pressure_angstr_kT  ! pressure in kT/A^3
  
  write(h,*)  'traj_file = ',traj_file(1:str_get_next_pos(traj_file,1,' ')) 
  write(h,*)  'energy_file = ',energy_file(1:str_get_next_pos(energy_file,1,' '))
  write(h,*)  'boxlength_file = ',boxlength_file(1:str_get_next_pos(boxlength_file,1,' '))
  write(h,*)  'frames_file = ',frames_file(1:str_get_next_pos(frames_file,1,' '))
  write(h,*)  'freq_file = ',freq_file(1:str_get_next_pos(freq_file,1,' '))
  write(h,*)  'param_out_file = ',param_out_file(1:str_get_next_pos(param_out_file,1,' ')) 

  write(h,*)  'random_seed = ',rnd_seed

   

end subroutine



subroutine Parameters_read_string(str) !DOC 
!DOC read the parameter nam=val pair from the string  
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use string, only : str_get_next_pos
   use error
!DOC Parameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: str !DOC input string

  ! character(SYSTEM_STRING_LENGTH) :: tmp
   character(SYSTEM_STRING_LENGTH) :: nam,val
   
   integer :: ibeg,iend,ilen
   logical :: finish

   !tmp(:) = str(:)

   !call str_subs(tmp,',',' ')
    
   ibeg = 1
   iend = 1
   finish = .FALSE.
   do while( .not. finish )

     iend = str_get_next_pos(str,ibeg,'=') - 1
     ilen = iend - ibeg + 1

     if( iend < 0 ) then
         write(error_message,*) 'Parameters_read_string: expected = after the offset',ibeg
         call error_throw(ERROR_PARAMETER)
         exit
     end if

     nam(:) = ' '
     nam(1:ilen) = str(ibeg:iend)

     ibeg = iend + 2
     iend = str_get_next_pos(str,ibeg,',') - 1

     if(iend < 0 ) then
        iend = str_get_next_pos(str,ibeg,' ') -1
        finish = .TRUE.
     end if

     ilen = iend - ibeg + 1

     val(:) = ' '
     val(1:ilen) = str(ibeg:iend)

     call read_parameter(nam,val)

     ibeg = iend + 2

   end do

end subroutine


subroutine Parameters_read(fname) !DOC
!DOC read parameters fule
   use io, only : io_open,io_close
   use error, only : error_code
   use string, only : str_subs,str_isempty
   use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC Parameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: fname !DOC file name

   character(SYSTEM_STRING_LENGTH) :: line
   character(SYSTEM_STRING_LENGTH) :: nam,val
   integer :: stat
   integer :: hfile

   hfile = io_open(fname,'r')
   if( error_code /= 0 ) return

   do while(.TRUE.)

      read(hfile,'(A)',iostat=stat) line

      if ( stat < 0 )  exit  ! end of file
      if ( str_isempty(line)) cycle  ! LINE IS EMPTY

      call str_subs(line,'=',' ')

      read(line,*) nam,val

      if ( nam(1:1) == '#' ) cycle  ! comment  

      call read_parameter(nam,val)

   end do


   call io_close(hfile)

   Parameters_areRead = .TRUE.
   
end subroutine

pure function BoxLength_from_density( density, nmol ) !DOC
!DOC calculate the box length from the density
    real(8) :: BoxLength_from_density
!DOC Parameters:
    real(8),intent(in) :: density !DOC in particles/nm^3
    integer,intent(in) :: nmol  !DOC number of molecules
!DOC Return value: 
  !DOC BoxLength

    real(8) :: volume

    volume = nmol / density   ! in nm^3
    BoxLength_from_density = volume**(1.d0/3.d0) * 10.d0 ! in angstroems

end function

subroutine Parameters_init(filename,nmol,extra_parameters_string) !DOC 
!DOC initialize the parameters of simulation
   use error
   use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC Parameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: filename !DOC name of the parameters file
   integer,intent(in) :: nmol !DOC number of molecules
   character(SYSTEM_STRING_LENGTH),intent(in),optional :: extra_parameters_string !DOC extra parameters in format nam1=val1,nam2=val2,...

   call Parameters_read(filename)

   if ( present(extra_parameters_string) ) then
       call Parameters_read_string(extra_parameters_string)
   end if

   if ( abs(BoxLength) < 1d-9 ) then

      if ( abs(density) < 1d-9 ) then
          write(error_message,*) 'Parameters_init: both: BoxLength and density are zero!'
          call error_throw(ERROR_PARAMETER) 
          return
      end if

      BoxLength = BoxLength_from_density( density, nmol )

    end if
   
    call Parameters_recalc( BoxLength )

end subroutine 


subroutine Parameters_recalc(box_length)   !DOC
!DOC recalculate the parameters for the new boxlength,temp etc 
   use error
   use constants, only : pi,elec,epsi0,boltz
   ! actually, only atomic_data % BoxLength is used
   ! but this manifestates, that atomic_data should be initialized before the parameters 
!DOC Parameters:
   real(8),intent(in) :: box_length !DOC if BoxLength = 0 --> recalculate from density and Nmol
                                        !  otherwise Nmol is unused

   real(8) :: Volume

   if(.not. Parameters_areRead ) then

      write(error_message,*) 'Parameters_init: Parameters are not read! run Parameters_read first'
      call error_throw(ERROR_INITIALIZATION)
      return

   end if
 
   BoxLength = box_length

 

!  alpha_angstr_m1 * A^-1 = alpha * BoxLength^-1 -->
!  alpha = alpha_angstr_m1 * A^-1 / BoxLength^-1 = alpha_angstr_m1 * BoxLength [A] * [A^-1]
 !  alpha = alpha_angstr_m1 * BoxLength 

   if ( alpha < 1d-9 ) then
       alpha = 2.d0 * sr  ! corresponds to rmax = L/2
   end if

   rmax = sr / alpha
   kmax = nint( sk * alpha / pi )

   kT_kcal_mol = boltz * temp / kcal_mol

!   BoxLength = atomic_data % BoxLength
   dbjr_a = elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10
   dbjr = dbjr_a / BoxLength
   xclb = dbjr    

   dr = dr_a / BoxLength
   d_angle = d_angle_degree * pi / 180d0

   alpha_over_sqrt_pi = alpha / dsqrt(pi)
   two_alpha_over_sqrt_pi = 2.d0 * alpha_over_sqrt_pi  

   kext = 2.d0 * xclb * pi / ( 2.d0 * external_permutivity + 1.d0)  
   minus_two_kext = -2*kext

   rmax2 = rmax**2

   log_max_volume_scaling = dlog( max_volume_scaling )

   ! pressure [ kT/A^3] = pressure [ J/m^3] / kT[J] * A[m]^3  = pressure [J/m^3] / kT * 10^-30
   pressure_angstr_kT = pressure / ( boltz * temp ) * 1.d-30

end subroutine 



End Module Parameters
