program mcrdf !DOC !PROG
!DOC !PROG !FILE calculate the Radial distribution functions between the atoms
use parameters, only : Parameters_init, &
                      traj_file, boxlength_file, frames_file, &
                      output_nbytes_xyz, output_nbytes_ang 
use MoleculeTable
use io, only : io_open,io_close
use composition 
use SystemSettings, only : SYSTEM_STRING_LENGTH, TRealArrayPointer

!DOC !PROG Usage:  mcrdf2 parameters.prm composition output_prefix [dr Rmax  [mol_labels [nskip[-maxfram]  ]  ]   ] '
!DOC !PROG
!DOC !PROG               *****   SIMPLIFIED VERSION  ****
!DOC !PROG
!DOC !PROG Arguments::
!DOC !PROG  ::  parameters.prm - parameters of the simulation. in format prm = val at each line '
!DOC !PROG                     they should have AT LEAST such fields: '
!DOC !PROG                        frames_file = ...'
!DOC !PROG                        traj_file = ... '
!DOC !PROG                        output_nbytes_xyz = ... '
!DOC !PROG                        output_nbytes_arg = ... '
!DOC !PROG  ::  composition - number and types of molecules in the system'
!DOC !PROG  ::  traj  - trajectory binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
!DOC !PROG  ::  frames.dat - information about the boxlength at each frame in traj '
!DOC !PROG  ::  output_prefix - prefix for output files'
!DOC !PROG                  output_files are:  output_prefixN1_N2.dat, where N1_N2 are numbers or labels of species'
!DOC !PROG                  ouput files are text  three-coulomn files. '
!DOC !PROG                  The columns are :  r  g(r) ' 
!DOC !PROG 
!DOC !PROG  ::  dr,Rmax - bin size and size for Rdf in angstroems. Default  dr=0.1  Rmax=12. '
!DOC !PROG    ::        Note: Nomrally Rmax should be less than min(BoxLength)/2 '   
!DOC !PROG   :: mol_labels  - optional coma separated labels used to produce the output files. If no labels given, numbers are used'
!DOC !PROG   :: nskip - number of frames to skip before the start of counting g(r) '
!DOC !PROG  :: maxfram - the number of frame to stop the accumulation


implicit none

character(SYSTEM_STRING_LENGTH) :: parameters_file,composition_file
character(SYSTEM_STRING_LENGTH) :: output_file
character(SYSTEM_STRING_LENGTH) :: tmpstr

logical :: extra_parameters = .FALSE.
character(SYSTEM_STRING_LENGTH) :: extra_parameters_string


character(SYSTEM_STRING_LENGTH) :: labels_string
logical :: has_labels = .FALSE.
integer,dimension(:),allocatable :: label_first_char, label_last_char


Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
 
real(8) :: dr = 0.1  ! angstr
real(8) :: Rmax = 12 ! angstr 
integer :: N   ! = Rmax / dr + 1

integer :: nskip ! number of frames to skip
integer :: maxframe ! collect up to the maxframe

integer :: i,j
integer :: htraj,hfram,hout
integer :: MaxMolecule

integer :: ntypes
!real(8),dimension(:) :: r_arr 
Type(TRealArrayPointer),dimension(:),allocatable :: cnt ! ntypes x ntypes rdfs
                                              ! i*ntypes + j is the rdf(ij)
Type(TRealArrayPointer),dimension(:),allocatable :: cnt2  !  sum_ij 1/4pi r_ij^2
Type(TRealArrayPointer),dimension(:),allocatable :: cnt3  !  1/ V_frame sum_ij 1/4pi_r_ij^2  (see the difference in RDF.pdf)
Type(TRealArrayPointer),dimension(:),allocatable :: cnt4 

real(8),dimension(:),allocatable :: N_total


real(8) :: avg_vol,avg_inv_vol

integer :: N_frame,N_vol


call parse_command_line
call read_input_files

call count_distances ! averages are counted here

!call count_volumes  ! is not good for "broken" files with wrong step numeration !!!
call write_averages

call save_rdfs
!call test_mcvol
write(*,*) 'N_total:'
do i=1,ntypes
   do j=i,ntypes
      write(*,*) N_total( (i-1)*ntypes + j )
   end do
end do
  

!   write(*,*) accepted
!end do

!call clearmc

call deallocate_all


contains 

subroutine count_distances !DOC 
!DOC accumulate the total counts at each distance
   use error
   use io, only : write_real_array
   use constants, only : four_pi

   integer :: ityp,jtyp
   integer :: offset, base_offset, r_offset
   real(8) :: x1,y1,z1,x2,y2,z2
   real(8) :: dx,dy,dz,r2,r,Rmax2 
   integer :: istep
   real(8) :: box_length,n_mv_steps,n_mv_accepted,n_vol_steps,n_vol_accepted
   integer :: ios
   integer :: nx,ny,nz,Nbox
   

   real(8) :: dV,r0,w0,w1
   real(8) :: V_frame
  
   ntypes = comp % n_types
   
   allocate(cnt(ntypes**2))
   allocate(cnt2(ntypes**2))
   allocate(cnt3(ntypes**2))
   allocate(cnt4(ntypes**2))
   allocate(N_total(ntypes**2))  
 
   N = aint( (Rmax - 1d-10) / dr )+1 +1 ! yet another +1 for linear interpolation
   
   base_offset = 0 
   
   do i=1,ntypes
      do j=i,ntypes
   
         offset = base_offset + j
         allocate( cnt( offset ) % ptr(N) )
         cnt( offset ) % ptr(:) = 0.d0

         allocate( cnt2( offset ) % ptr(N) )
         cnt2( offset ) % ptr(:) = 0d0 
 
         allocate( cnt3( offset ) % ptr(N) )
         cnt3( offset ) % ptr(:) = 0d0

         allocate( cnt4( offset) % ptr(N) )
         cnt4( offset ) % ptr(:) = 0d0         

         N_total(offset) = 0.d0
 
      end do
      base_offset = base_offset + ntypes
   end do
  
!   do i=0,N-1
!      r_arr(i+1) = i * dr
!   end do
 
   htraj = io_open(traj_file,'b')
   hfram = io_open(frames_file,'r')

   call error_set_catch(ERROR_IO) 
   
   Rmax2 = Rmax**2

   avg_vol = 0d0
   avg_inv_vol = 0d0

   N_frame = 0
   do while (.TRUE.)
  

      call MoleculeTable_load_binary(mol_tab, htraj, comp, output_nbytes_xyz, output_nbytes_ang )
!     MoleculeTable_load_binary(this,fid,comp, nbytes_xyz,nbytes_ang)

      if( error_code /= 0 ) exit  ! the frames are over

      N_frame = N_frame + 1

      if( maxframe > 0 .and. N_frame > maxframe ) exit
  

      read(hfram,*,iostat=ios) istep, box_length, n_mv_steps,n_mv_accepted, n_vol_steps, n_vol_accepted
    
      if( ios < 0) exit ! end of file      
    

      if( N_frame .ge. nskip ) then

      write(*,*) 'Nframe:',N_frame

      else
          write(*,*) ' Frame ',N_frame, ' ---> skiped'
          cycle
      end if

 
      V_frame = box_length**3
  

      avg_vol = avg_vol + V_frame
      avg_inv_vol = avg_inv_vol  + 1d0 / V_frame
      ! if Rmax > box_length/2 ==>  consider periodic images 
      if(Rmax < box_length/2) then
         Nbox = 0
      else
       !        L    L/2  L/2     L
       ! |---------|----*----|---------|
       !                *-------)
       !                    R
       !                     |--)
       !                     R-L/2  
       !      minimal integer lager than R-L/2 
        Nbox = AINT ( (Rmax-box_length/2) / box_length ) + 1
      end if

!      Nbox = 0

      do i=1,MaxMolecule

         ityp = mol_tab % mol_type(i)

         x1 = mol_tab % x(i)
         y1 = mol_tab % y(i)
         z1 = mol_tab % z(i)

         do j = 1,MaxMolecule

            jtyp = mol_tab % mol_type(j)

            x2 = mol_tab % x(j)
            y2 = mol_tab % y(j)
            z2 = mol_tab % z(j)

            if( x1 - x2 >  0.5d0 ) x2 = x2 + 1.d0
            if( x1 - x2 < -0.5d0 ) x2 = x2 - 1.d0  
            
            if( y1 - y2 >  0.5d0 ) y2 = y2 + 1.d0
            if( y1 - y2 < -0.5d0 ) y2 = y2 - 1.d0

            if( z1 - z2 >  0.5d0 ) z2 = z2 + 1.d0
            if( z1 - z2 < -0.5d0 ) z2 = z2 - 1.d0 

            do nx = -Nbox,Nbox
            do ny = -Nbox,Nbox
            do nz = -Nbox,Nbox

               if ( i == j ) then
                   if ( (nx == 0) .and. (ny == 0) .and. (nz == 0) ) cycle
               end if

               dx = (x1 - x2 + nx) * box_length
               dy = (y1 - y2 + ny) * box_length
               dz = (z1 - z2 + nz) * box_length

               r2 = dx**2 + dy**2 + dz**2
               if( r2 > Rmax2 - 1d-10 ) cycle
   
!               N_total = N_total + 1 ! total number of distances
 


               r = sqrt(r2) 
   
               r_offset = aint( r / dr ) + 1
     
               if ( r_offset > N ) then
                   write(*,*) 'r_offset',r_offset,'N',N,'r',r,'Rmax',Rmax
               end if
               ! in ideal gas any distance < Rmax is equivally probable
               ! probability to find the offset  [r0;r0+dr] is proportional to the volume of segment
               ! in comparison with the volume of the sphere 
               ! r0 = (r_offset - 1)*dr
               ! dV = 4*pi*r0**2 * dr + 4*pi*r0*dr**2 + 4/3*pi*dr**3
               !
               ! P_id([r;r+dr]) =  dV / (4/3 pi Rmax^3)
               ! in our case the same probability is calculated as 
               ! P([r;r+dr] = cnt([r;r+dr]) / N_total
               ! which gives:
               ! g(r) = P/P_idd = cnt(r) / N_total *  Vmax / dV(r) 
               ! where Vmax = 4/3 pi Rmax^3
                
               if ( ityp < jtyp ) then
                   offset = (ityp - 1) * ntypes + jtyp  
               else
                   offset = (jtyp - 1) * ntypes + ityp
               end if

               N_total( offset ) = N_total( offset ) + 1.d0   
   !            write(*,*) 'r_offset:',r_offset,'N:',N 
           !    cnt( offset ) % ptr(r_offset ) = cnt( offset ) % ptr(r_offset) + 1
           ! to improve accuracy, we use the linear interpolation between the r_offset and r_offset + 1
                r0 = (r_offset - 1) * dr;
                w0 = (r0+dr - r) /dr;
                w1 = (r-r0) / dr     
   
                cnt( offset ) % ptr(r_offset) = cnt(offset) % ptr(r_offset) + w0
                cnt( offset ) % ptr(r_offset+1) = cnt(offset) % ptr(r_offset + 1) + w1

                cnt2( offset ) % ptr(r_offset ) = cnt2( offset ) % ptr( r_offset ) + 1d0 /  r2 

                cnt3( offset ) % ptr(r_offset ) = cnt3( offset ) % ptr( r_offset ) + 1d0 / r2 / V_frame

                cnt4( offset ) % ptr(r_offset) = cnt4(offset) % ptr(r_offset) + 1

            end do  ! nz
            end do  ! ny
            end do  ! nx
         end do ! j
      end do  ! i
   
   end do ! write .TRUE. (reading the frames )


   call error_clear_catch(ERROR_IO)

   

   call io_close(htraj)  
   call io_close(hfram)

   avg_vol = avg_vol / N_frame
   avg_inv_vol = avg_inv_vol / N_frame

!   call write_real_array(0,cnt(1) % ptr, N)
  

end subroutine


subroutine count_volumes !DOC
!DOC calculate the mean volume and mean inverse volume

   real(8) :: sum_vol,sum_inv_vol
   real(8) :: box_length,vol
   integer :: hfile
   integer :: istep, istep_old
   integer :: N_step, N_vol
   integer :: ios

   hfile = io_open(boxlength_file,'r')

   istep_old = 0
   do while (.TRUE.)

      read(hfile,*,iostat=ios) istep,box_length
      
      if( ios /= 0) exit

      vol = box_length**3
      N_step = istep - istep_old
      N_vol =  N_vol + N_step
      
      sum_vol = sum_vol + vol * N_step 
      sum_inv_vol = sum_inv_vol + (1.d0 / vol) * N_step 
      
      istep_old = istep

   end do

   call io_close(hfile)
 
   avg_vol = sum_vol / N_vol
   avg_inv_vol = sum_inv_vol / N_vol

end subroutine

subroutine write_averages !DOC
!DOC write average values
 
   write(*,*) 'Average_Volume:',avg_vol,' [Angstr^3]'
   write(*,*) 'Average_Inverse_Volume:',avg_inv_vol,'[Angstr^-3]'
   write(*,*) 'Average_density  [particles/Angstr^3]:  '
   write(*,*) ' N/<V>:', MaxMolecule/avg_vol
   write(*,*) 'N*<1/V>',MaxMolecule * avg_inv_vol  

   write(*,*) 'N_frame:',N_frame
   write(*,*) 'N_total:',N_total

end subroutine


subroutine save_rdfs !DOC 
!DOC save rdfs to file
  use constants, only : pi,four_pi
  use string, only : str_get_next_pos

  integer :: hfile
  integer :: k
  real(8) :: SmallBall
  real(8) :: four_pi_dr
  real(8) :: rdf1,rdf2,rdf3
  real(8) :: dV
  real(8) :: r
  integer :: strlen,strlen2
  integer :: offset
  real(8) :: Vmax
  real(8) :: rho_i,rho_j
  real(8) :: factor

  ! rdf(r) = count(r) / N_frame / dV  *  ( 1 / rho )
  !  dV \approx 4 pi r^2 dr 
  !  dV = 4/3 pi ( r+dr)^3  - 4/3 pi r^3 =  4 pi r^2 dr + 4 pi r dr^2 + 4/3 pi dr^3 =
  !                                         4 pi * dr ( r^2 + r * dr)  + 4/3 pi *dr^3  
  ! 
  ! two possibilities:    1/rho =  <V> / N 
  !                       1/rho =  1 / ( N * <1/V> )

   four_pi_dr = four_pi * dr
   SmallBall = 4.d0 / 3.d0 * pi * dr**3
   Vmax = 4.d0  /3.d0 * pi * Rmax**3

   strlen = str_get_next_pos(output_file,1,' ')-1

   do i=1,ntypes
      do j=i,ntypes
   

         
         offset = (i-1)*ntypes + j

         rho_i = comp % mol_numbers(i) * avg_inv_vol
         rho_j = comp % mol_numbers(j) * avg_inv_vol

         write(*,*) 'rho_i',rho_i,'rho_j',rho_j 


         if ( has_labels ) then
            write(tmpstr,'(A,A,A,A,A)') output_file(1:strlen),&
                            labels_string(label_first_char(i):label_last_char(i)),'_', &
                            labels_string(label_first_char(j):label_last_char(j)),'.dat'
           
         else
            write(tmpstr,'(A,I2.2,A,I2.2,A)') output_file(1:strlen),i,'_',j,'.dat'
         end if

         strlen2 = str_get_next_pos(tmpstr,1,' ')
         write(*,*) 'Saving ',tmpstr(1:strlen2)

         hfile = io_open(tmpstr,'w')
         r = 0.d0
         do k=1,N
              
            dV =  four_pi_dr * ( r**2 + r * dr ) + SmallBall

            ! see count_ditances: g(r) = cnt(r)/N_total * Vmax / dV

            if ( i == j ) then
               factor = 1.d0
            else
               factor = 2.d0
            end if

            if ( N_total(offset) > 0 ) then
              ! rdf = cnt / Nfram * V / dV  / Nmol1 / Nmol2
               rdf1 = cnt4( offset ) % ptr(k) / N_frame / avg_inv_vol / comp % mol_numbers(i) / comp % mol_numbers(j) / dV / factor
            else
               rdf1 = 0d0
            end if

!            cnt2( (i-1) * ntypes + j) % ptr(k) 

!            rdf1 = cnt(i*ntypes + j) % ptr(k) / N_frame *  avg_vol / comp % mol_numbers(j)
!            rdf2 = cnt(i*ntypes + j) % ptr(k) / N_frame *  1.d0 /  ( comp % mol_numbers(j)  * avg_inv_vol ) 
!             
            rdf2 = 1d0 / (rho_i * rho_j * dr ) * cnt2( offset ) % ptr(k) / N_frame / four_pi / avg_vol
            rdf3 = 1d0 / (rho_i * rho_j * dr ) * cnt3( offset ) % ptr(k) / N_frame / four_pi

            !write(hfile,*)  r+dr/2, cnt( offset ) % ptr(k) , rdf1, cnt2(offset) % ptr(k), rdf2, cnt3(offset) % ptr(k) ,rdf3, cnt4( offset) % ptr(k)
         write(hfile,*)  r+dr/2, rdf1 !, cnt2(offset) % ptr(k), rdf2, cnt3(offset) % ptr(k) ,rdf3, cnt4( offset) % ptr(k)


            r = r + dr
         end do
         call io_close(hfile)

      end do
   end do
   

end subroutine


subroutine parse_command_line !DOC 
!DOC read command line arguments
  use string
  integer :: pos 
if ( iargc() < 3 ) then

   write(*,*)  'Usage:  mcrdf parameters.prm composition output_prefix [dr Rmax  [mol_labels [nskip[-maxfram]  ]  ]   ] '
   write(*,*)
   write(*,*)  '                ****  SIMPLIFIED VERSION **** '
   write(*,*)
   write(*,*)  '   parameters.prm - parameters of the simulation. in format prm = val at each line '
   write(*,*)  '                    they should have AT LEAST such fields: '
   write(*,*)  '                       frames_file = ...'
   write(*,*)  '                       traj_file = ... '
   write(*,*)  '                       output_nbytes_xyz = ... '
   write(*,*)  '                       output_nbytes_arg = ... '
   write(*,*)  '   composition - number and types of molecules in the system'
   write(*,*)  '   traj  - trajectory binary format, (fixed point):   Nmolecule records: (x,y,z,theta,phi,psi)  '
   write(*,*)  '   frames.dat - information about the boxlength at each frame in traj '
   write(*,*)  '   output_prefix - prefix for output files'
   write(*,*)  '                 output_files are:  output_prefixN1_N2.dat, where N1_N2 are numbers or labels of species'
   write(*,*)  '                 ouput files are text  three-coulomn files. '
   write(*,*)  '                 The columns are :  r g(r)' 
  write(*,*)  ' '
   write(*,*)  '   dr,Rmax - bin size and size for Rdf in angstroems. Default  dr=0.1  Rmax=12. '
   write(*,*)  '           Note: Nomrally Rmax should be less than min(BoxLength)/2 '   
   write(*,*)  '   mol_labels  - optional coma separated labels used to produce the output files. If no labels given, numbers are used'
   write(*,*)  '   nskip - number of frames to skip before the start of counting g(r) '

   stop
end if

call getarg(1,parameters_file)
call getarg(2,composition_file)
call getarg(3,output_file)



if ( iargc() .ge. 4 ) then
    call getarg(4,tmpstr)
    read(tmpstr,*) dr
end if

if (iargc() .ge. 5 ) then
   call getarg(5,tmpstr)
   read(tmpstr,*) Rmax
end if

if ( iargc() .ge. 6 ) then
   call getarg(6,labels_string)
   has_labels = .TRUE.
else
   has_labels = .FALSE.
end if

if ( iargc() .ge. 7 ) then

   call getarg(7,tmpstr)

   pos = str_get_next_pos(tmpstr,1,'-')
   if ( pos > 0 ) then
     call str_subs(tmpstr,'-',' ')
     read(tmpstr,*)  nskip,maxframe 
   else
      read(tmpstr,*)  nskip
      maxframe = -1
   end if
else
   nskip = 0
end if


end subroutine


subroutine read_input_files !DOC
!DOC read the input files
   use string, only : str_subs,str_get_next_pos
   use error

   integer :: i,pos,old_pos
!write(*,*)  'input=',input_file(1:20), 'output=',output_file(1:20), 'nbytes_xyz=',nbytes_xyz, 'nbytes_ang=',nbytes_ang

! read the system composition
call Composition_nulify(comp)
call Composition_read_from_file(comp, composition_file )

if( has_labels ) then

  call str_subs(labels_string,',',' ') 

  old_pos = 1

  allocate(label_first_char( comp % n_types ) )
  allocate(label_last_char( comp % n_types ) )

  do i = 1,comp % n_types

      if( labels_string(old_pos:old_pos) == ' ' ) then
          write(error_message,*) 'read_input_files: not enough labels (',&
                                   comp % n_types,'expected ):',labels_string(1:old_pos)
          call error_throw( ERROR_PARAMETER)
          return
      end if
 
      pos = str_get_next_pos(labels_string,old_pos,' ')
      
      label_first_char(i) = old_pos
      label_last_char(i) = pos - 1

      old_pos = pos + 1      

  end do


end if


!MaxAtom = Composition_count_atoms(comp)
MaxMolecule = SUM( comp % mol_numbers(1:comp % n_types) )

! read and initialize the  parameters

if ( .not.extra_parameters ) then
     call Parameters_init(parameters_file,MaxMolecule)  ! MAxMolecule is used for the BoxLength calculation (from density, if required)
else
     call Parameters_init(parameters_file,MaxMolecule,extra_parameters_string)
end if


! read the molecule table
call MoleculeTable_nulify(mol_tab)

!hfile = io_open(moltab_file,'b')
!call MoleculeTable_load_binary( mol_tab, hfile, comp, input_nbytes_xyz, input_nbytes_ang )
!call io_close(hfile)

! initialize the atomic_data
!call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
!call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
!call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )

end subroutine


subroutine deallocate_all
! deallocate everything

!call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)


if ( has_labels ) then
   deallocate( label_first_char )
   deallocate( label_last_char )
end if
!call LJTypes_dealloc(lj_types)


end subroutine

end program mcrdf
