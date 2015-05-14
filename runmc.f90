Module runmc !DOC
!DOC !FILE Module which performs the MC cycles. Is called from the mc_main.f90
use MoleculeTable, only : TMoleculeTable
use MC, only : movemc, mcvol, mc_xchange
use BiasedRandom

implicit none

integer :: ifirst = 1 
real(8) :: n_mv_steps =0.d0, n_mv_accepted=0.d0, n_xchg_accepted=0d0, n_xchg_steps=0d0,n_vol_steps=0.d0, n_vol_accepted=0.d0 


Type(TBiasedRandom) :: move_freq, xchange_freq  ! arrays for choosing the molecules according to the biased probability
real(8) :: xchange_probability  ! probability of the xchange step. 
                                ! is equal to SUM(xchange_relative_prob) / (SUM(xchange_relative_prob) + SUM(move_relative_prob) )


!!!  "Frequences" represent the probabilities ot move one or other molecule
!!! they are stored in the freq_file given in parameters)
!!! if the file is "none" --> equal frequences are used

contains



subroutine init_freq(nmol) !DOC
!DOC initialize the frequences (read them from file or just make flat if no file given)
!DOC Parameters:
    use parameters, only : freq_file

    integer,intent(in) :: nmol !DOC number of molecules

    integer :: hfile
    real(8),dimension(:),allocatable :: mv_num, xchange_num

    real(8) :: Smv,Sxch

    allocate(mv_num(nmol))
    allocate(xchange_num(nmol))

    if( freq_file(1:5) == 'none ' ) then

        mv_num(:) = 1.d0
        xchange_num(:) = 0.d0

    else

       call read_freq_file(freq_file,mv_num,xchange_num, nmol )

    end if

    call BiasedRandom_alloc(move_freq,nmol)
    call BiasedRandom_alloc(xchange_freq,nmol)

    call BiasedRandom_init(move_freq, mv_num )
    call BiasedRandom_init(xchange_freq, xchange_num )

    Smv  = SUM(mv_num)
    Sxch = SUM(xchange_num)

    xchange_probability = Sxch / ( Smv + Sxch )

    deallocate(mv_num)
    deallocate(xchange_num)

end subroutine 


subroutine read_freq_file( fname, mv_num, xchange_num, nmol  ) !DOC
!DOC read frequences from file
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use io, only : io_open, io_close
   use string
   use error

!DOC # Frequences file format: 
!DOC # first non-comment line - number of intervals
!DOC # next lines:
!DOC # column 1:  interval in format  first-last
!DOC # column 2:  relative probability to choose the molecule in this interval
!DOC # column 3:  realtive probability to xchange the molecules from this interval
!DOC Parameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: fname !DOC name of the file
   real(8),dimension(:) :: mv_num, xchange_num !DOC output: relative probabilities of movement and exchange
   integer,intent(in) :: nmol !DOC number of molecules

   integer :: hfile
   character(SYSTEM_STRING_LENGTH) :: tmpstr

   integer :: nlin ! number of line in the file
   integer :: stat
   integer :: nrec ! number of records in file
   integer :: i,j
   integer :: first,last
   real(8) :: mv_prob,xchange_prob

   hfile = io_open(fname, 'r')
  
  nlin = 0
  ! skip the comments and empty lines in the file
   do while(.TRUE.)    

     read(hfile,'(A)',iostat=stat) tmpstr   
     nlin = nlin + 1

     if( stat /= 0 ) then
         write(error_message,*) 'read_freq_file: Unexpected end of file ',fname(1:20),'at line ',nlin
         call error_throw(ERROR_IO)
         return
     end if

     if ( .not.str_isempty( tmpstr ) .AND. &
        ( str_get_next_pos( tmpstr, 1, '#') .le. 0 ) ) exit

   end do
  
   ! read the number of records
   read(tmpstr,*) nrec

   ! initialize mc_num, xchange_num arrays with -1
   mv_num(:) = -1.d0
   xchange_num(:) = -1.d0

   ! read the records
   do i=1,nrec

      read(hfile,'(A)',iostat=stat) tmpstr
      nlin = nlin + 1     
 
      if( stat /= 0 ) then
         write(error_message,*) 'read_freq_file: Unexpected end of file ',fname(1:20),' at line ',nlin
         call error_throw(ERROR_IO)
         return
     end if

     call str_subs(tmpstr,'-',' ')

     read(tmpstr,*) first,last, mv_prob, xchange_prob 

     if ( last > nmol ) then
        write(*,*) 'read_freq_file: File ',fname(1:20),' is inconsistent with the composition'
        call error_throw(ERROR_PARAMETER)
       
     end if

     do j=first,last

         if( mv_num(j) > 0 ) then
              write(error_message,*) 'read_freq_file: Crossing intervals in ',fname(1:20),' at line ',nlin
              call error_throw(ERROR_PARAMETER)
              return
         end if

         mv_num(j) = mv_prob
         xchange_num(j) = xchange_prob

     end do

   end do

   call io_close(hfile)

   do i=1,nmol
      
      if( mv_num(i) < 0 ) then
          write(error_message,*) 'read_freq_file: not all frequences defined (',fname(1:20),' is inconsistent with composition?)'
          call error_throw(ERROR_PARAMETER)
          return
      end if

   end do

end subroutine


!!!!!! RUN MC 

subroutine run_mc( mol_table, ncycle ) !DOC 
!DOC run the MC cycles
   use MRandom
   use parameters, only : pressure_step_multiplier, n_store_traj_interval,n_store_energy_interval
   use io, only : io_open,io_close,write_xyz_array,write_real_array
   use MC, only : RealSpace_TotalSum, KSpace_TotalSum
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use SumSinCosKR, only : TSumSinCosKR
!DOC Parameters:
   Type(TMoleculeTable),intent(inout) :: mol_table !DOC MoleculeTable
   integer,intent(in) :: ncycle !DOC max number of cycles

   integer :: i,j,nmol,imol
   integer :: imol1,imol2  ! for xchange 
   real(8) :: rnd
   logical :: accepted
   integer :: hfile
   integer :: natom,nk
   real(8),dimension(:),pointer :: xx,yy,zz
   character(SYSTEM_STRING_LENGTH) :: filename

   Type(TSumSinCosKR),pointer :: coul,lj

   nmol = mol_table % nmol
   call init_freq(nmol)

   do i=ifirst,ncycle
 
      rnd = rand() * ( pressure_step_multiplier + 1)
      if ( (pressure_step_multiplier > 0) .and. (rnd > pressure_step_multiplier) ) then 
                ! pressure
        call mcvol(mol_table,accepted)

        if( accepted ) n_vol_accepted = n_vol_accepted + 1
        n_vol_steps = n_vol_steps + 1
!        write(*,*) 'VOL'
!        write(*,'(A,$)'),'v' 

        if(accepted) then
           write(*,'(A,$)') 'v'
           call store_boxlength(i)
        else
           write(*,'(A,$)') '_'
        end if       


      else


        ! move or xchange step
        write(*,'(A,$)') '.'

        do j=1,nmol

           rnd = rand()
         !  write(*,*) 'rnd',rnd,'xchange_prob',xchange_probability
         !  xchange_probability=1.d0
           if ( rnd < xchange_probability) then
         
              n_xchg_steps = n_xchg_steps + 1
              ! xchange step
              imol1 = BiasedRandom_choose(xchange_freq)
              imol2 = imol1

        !      write(*,*) 'xchange_freq % N',xchange_freq % N
        !      call write_xyz_array(0,xchange_freq % freq_prob, xchange_freq % freq_beg,&
        !                             xchange_freq % freq_end, xchange_freq % N) 

         !     return
              do while( imol2 == imol1 ) 
                 imol2 = BiasedRandom_choose(xchange_freq)
!                 write(*,*) 'imol1',imol1,'imol2',imol2
              end do
              
              call mc_xchange(imol1, imol2, mol_table, accepted ) 
              
              if(accepted) n_xchg_accepted = n_xchg_accepted + 1
           else
            ! move step
!              imol = aint( rand() * nmol + 1 )
              imol = BiasedRandom_choose(move_freq)
 !             write(*,*) 'imol',imol
              call movemc(imol, mol_table,accepted)

           end if

           if( accepted ) n_mv_accepted = n_mv_accepted + 1
           n_mv_steps = n_mv_steps + 1


        end do

       

      end if 



      if( mod( n_mv_steps / nmol + n_vol_steps, n_store_energy_interval ) < 1d-10 ) then
         call store_energy(i)
         write(*,'(A,$)') 'e'
      end if

      if ( mod( ( n_mv_steps) / nmol + n_vol_steps , n_store_traj_interval ) < 1d-10 ) then
         call store_traj_frame(mol_table,i)
         write(*,*) i
      end if
 
   end do
 
   call BiasedRandom_dealloc(move_freq)
   call BiasedRandom_dealloc(xchange_freq)
 

!   natom = RealSpace_TotalSum % atomic_data % natom
!  
!   xx => RealSpace_TotalSum % atomic_data % xx
!   yy => RealSpace_TotalSum % atomic_data % yy
!   zz => RealSpace_TotalSum % atomic_data % zz
!
!
!   filename = 'xyz.dat'
!   hfile = io_open(filename,'w');
!   call write_xyz_array(hfile,xx,yy,zz,natom)
!   call io_close(hfile)
!
!   filename = 'forces.dat'
!   hfile = io_open(filename,'w');
!   call write_xyz_array(hfile,RealSpace_TotalSum % Fx,RealSpace_TotalSum % Fy,RealSpace_TotalSum % Fz,natom)
!   call io_close(hfile)
!
!   coul => KSpace_TotalSum % rho_squared_total % sumsincos_coulomb
!   lj => KSpace_TotalSum % rho_squared_total % sumsincos_lj(1)
!
!   nk = KSpace_TotalSum % grid % nk
!
!   filename = 'sumcos_ppp.txt'
!   hfile = io_open(filename,'w')
!   call write_real_array(hfile,coul % sumcos_ppp, nk)
!   call io_close(hfile)
!
!   filename = 'sumsin_ppp.txt'
!   hfile = io_open(filename,'w')
!   call write_real_array(hfile,coul % sumsin_ppp, nk)
!   call io_close(hfile)
!    
!   filename = 'sumcos_ppp_lj.txt'
!   hfile = io_open(filename,'w')
!   call write_real_array(hfile,lj % sumcos_ppp, nk)
!   call io_close(hfile)
!
!   filename = 'sumsin_ppp_lj.txt'
!   hfile = io_open(filename,'w')
!   call write_real_array(hfile,lj % sumsin_ppp, nk)
!   call io_close(hfile)

end subroutine

subroutine store_energy(iframe) !DOC
!DOC save the energies to the energy file
  use io,only : io_open,io_close
  use SystemSettings, only : SEEK_END
  use MC, only : mc_total_energy,RealSpace_TotalSum,KSpace_TotalSum,U_total_ew,U_total_lj 
  use parameters, only : energy_file
!DOC Parameters:
  integer,intent(in) :: iframe !DOC current frame
  integer :: hfile   

! energy
   hfile = io_open(energy_file,'w')  
   call fseek(hfile,0,SEEK_END)
   write(hfile,*) iframe,mc_total_energy(),U_total_ew,U_total_lj,&
                  RealSpace_TotalSum % uu_ew,RealSpace_TotalSum % uu_lj6,RealSpace_TotalSum % uu_lj12, &
                  KSpace_TotalSum % energy_coulomb,KSpace_TotalSum % energy_lj
   call io_close(hfile)

end subroutine

subroutine store_boxlength(iframe) !DOC 
!DOC save current boxlength to the box length file
  use io,only : io_open,io_close
  use SystemSettings, only : SEEK_END
  use parameters, only : BoxLength,boxlength_file
!DOC Parameters:
  integer,intent(in) :: iframe !DOC current frame number
  integer :: hfile   

  ! boxlength
   hfile = io_open(boxlength_file,'w')  
   call fseek(hfile,0,SEEK_END)
   write(hfile,*) iframe,BoxLength
   call io_close(hfile)

end subroutine

subroutine store_traj_frame(mol_tab,framecount) !DOC
!DOC save current molecule positions to the trajectory file
   use MoleculeTable, only : TMoleculeTable,MoleculeTable_calcPositionsOrientations, MoleculeTable_save_binary
   use MC, only : RealSpace_TotalSum,mc_total_energy
   use EwaldSumRealSpace, only : TEwaldSumRealSpace
   use parameters, only : BoxLength,output_nbytes_xyz,output_nbytes_ang, traj_file,frames_file
   use SystemSettings, only : SEEK_END
   use io, only : io_open,io_close
!DOC Parameters:
   Type(TMoleculeTable),intent(inout) :: mol_tab !DOC MoleculeTable
   integer,intent(in) :: framecount !DOC current frame

   real(8),dimension(:),pointer :: xx,yy,zz
   integer :: hfile
  


   ! trajectory

   xx => RealSpace_TotalSum % atomic_data % xx
   yy => RealSpace_TotalSum % atomic_data % yy
   zz => RealSpace_TotalSum % atomic_data % zz

   hfile = io_open(traj_file,'b')
   call fseek(hfile,0,SEEK_END)

   call MoleculeTable_calcPositionsOrientations( mol_tab, xx, yy, zz, BoxLength )
   call MoleculeTable_save_binary(mol_tab,hfile,output_nbytes_xyz,output_nbytes_ang)

   call io_close(hfile)

   ! boxlength
!   hfile = io_open(boxlength_file,'w')  
!   call fseek(hfile,0,SEEK_END)
!   write(hfile,*) BoxLength
!   call io_close(hfile)

   ! framecount
   hfile = io_open(frames_file,'w')  
   call fseek(hfile,0,SEEK_END)
   write(hfile,'(I12,F15.10,F12.0,F12.0,F12.0,F12.0,F12.0,F12.0,I12)') framecount,BoxLength,n_mv_steps,n_mv_accepted,&
                     n_vol_steps,n_vol_accepted, n_xchg_steps, n_xchg_accepted, time()
   call io_close(hfile)

end subroutine

subroutine read_traj_frame(hfile,stat)
  use parameters, only : BoxLength

  integer,intent(in) :: hfile
  integer,intent(out) :: stat

  read(hfile,*,iostat=stat) ifirst,BoxLength,n_mv_steps,n_mv_accepted,n_vol_steps,n_vol_accepted

end subroutine

subroutine read_last_frame(frames_file,nframes) !DOC 
!DOC read the last frame from the frames.dat file
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use io, only : io_open,io_close
   use error
   use parameters, only : BoxLength
!DOC PArameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: frames_file !DOC frames file
   integer  :: nframes !DOC output: number of frames

   integer :: hfile,stat
  

   hfile = io_open(frames_file,'r')

   nframes = 0
   do while(.TRUE.)

      call read_traj_frame(hfile,stat) 
      if ( stat /= 0) exit

      nframes = nframes + 1
   end do

   call io_close(hfile)

   ifirst = ifirst + 1
 
!   write(*,*) 'ifirst',ifirst,'BoxLength',BoxLength,'nframes=',nframes

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module runmc
