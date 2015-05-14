program mc_rmsd !DOC !PROG
!DOC !PROG !FILE calculate the displacement from the original position for a given set of molecules

use SystemSettings, only : SYSTEM_STRING_LENGTH
use composition
use MoleculeTable
use io, only : io_open, io_close
use string

!DOC !PROG Usage:  mc_rmsd  system.composition  traj.moltraj frames.dat interval  >  output.dat'
!DOC !PROG Arguments:
!DOC !PROG ::  interval - in format num1-num2, where num1 and num2 are the numbers of the first and the last molecules for rmsd ' 
!DOC !PROG ::   output.dat - text file with columns (one per molecule), displacement from the initial position' 


implicit none

character(SYSTEM_STRING_LENGTH) :: composition_file,frames_file,traj_file,interval_str
character(SYSTEM_STRING_LENGTH) :: tmpstr

integer :: nbytes_xyz,nbytes_ang

Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
!Type(TAtomicData) :: atomic_data
integer :: i,ii
integer :: htraj,hfram
real(8) :: box_length, half_box

integer :: ifram

integer :: first_mol,last_mol,nmol
integer :: stat
logical :: first_frame


real(8),dimension(:),allocatable :: x0,y0,z0, xprev,yprev,zprev
integer,dimension(:), allocatable :: nx,ny,nz ! number of the turning the box boundaries
real(8) :: dx,dy,dz,dist, xcurr,ycurr,zcurr



if ( iargc() < 4 ) then

   write(*,*)  '  calculates the displacement(frame) for a given molecules'
   write(*,*)  'Usage:  mc_rmsd  system.composition  traj.moltraj frames.dat interval  >  output.dat'

   write(*,*)  '   interval - in format num1-num2, where num1 and num2 are the numbers of the first and the last molecules for rmsd ' 
   write(*,*)  '   output.dat - text file with columns (one per molecule), displacement from the initial position' 

   stop
end if

call getarg(1,composition_file)
call getarg(2,traj_file)
call getarg(3,frames_file)
call getarg(4,interval_str)
!call getarg(5,output_file)

call str_subs(interval_str,'-',' ')

read(interval_str,*)  first_mol,last_mol


!write(*,*)  'input=',input_file(1:20), 'output=',output_file(1:20), 'nbytes_xyz=',nbytes_xyz, 'nbytes_ang=',nbytes_ang

! read the system composition
call Composition_nulify(comp)
call Composition_read_from_file(comp, composition_file )

! read the molecule table
call MoleculeTable_nulify(mol_tab)

htraj = io_open(traj_file,'b')
hfram = io_open(frames_file,'r')

first_frame = .TRUE.

nmol = last_mol-first_mol + 1

allocate(x0(nmol))
allocate(y0(nmol))
allocate(z0(nmol))

allocate(xprev(nmol))
allocate(yprev(nmol))
allocate(zprev(nmol))

allocate(nx(nmol))
allocate(ny(nmol))
allocate(nz(nmol))

nx(:) = 0
ny(:) = 0
nz(:) = 0

do while(.TRUE.)

   read(hfram,'(A)',iostat=stat) tmpstr
   
   if( stat /= 0 ) exit

   read(tmpstr,*)  ifram,box_length

   call MoleculeTable_load_binary(mol_tab, htraj, comp, 2,2) 

   write(*,'(I10)',advance='no') ifram

   do i=first_mol,last_mol

       ii = i - first_mol + 1

       if ( first_frame ) then

           x0(ii) = mol_tab % x(i) * box_length
           y0(ii) = mol_tab % y(i) * box_length
           z0(ii) = mol_tab % z(i) * box_length
          
           xprev(ii) = x0(ii)
           yprev(ii) = y0(ii)
           zprev(ii) = z0(ii)
 
       end if
 
       half_box = box_length * 0.5
      

       xcurr = mol_tab % x(i)  * box_length
       ycurr = mol_tab % y(i)  * box_length
       zcurr = mol_tab % z(i)  * box_length

       if ( xcurr - xprev(ii) > half_box ) nx(ii) = nx(ii) - 1  ! if become very big --> turn the lowest boundary
       if ( xcurr - xprev(ii) < -half_box) nx(ii) = nx(ii) + 1  ! very small - turn the upper boundary
 
       if ( ycurr - yprev(ii) > half_box ) ny(ii) = ny(ii) - 1
       if ( ycurr - yprev(ii) < -half_box) ny(ii) = ny(ii) + 1
        
       if ( zcurr - zprev(ii) > half_box ) nz(ii) = nz(ii) - 1
       if ( zcurr - zprev(ii) < -half_box) nz(ii) = nz(ii) + 1
 

       dx = dabs( xcurr + nx(ii) * box_length - x0(ii))
       dy = dabs( ycurr + ny(ii) * box_length - y0(ii))
       dz = dabs( zcurr + nz(ii) * box_length - z0(ii))

       xprev(ii) = xcurr
       yprev(ii) = ycurr
       zprev(ii) = zcurr

!       if ( dx > half_box ) dx = dx - box_length
!       if ( dy > half_box ) dy = dy - box_length
!       if ( dz > half_box ) dz = dz - box_length

       dist = sqrt(dx**2 + dy**2 + dz**2)
       
       write(*,'(F10.6)',advance='no') dist


   end do   
   
   write(*,*)
   first_frame = .FALSE.

end do


call io_close(htraj)
call io_close(hfram)


deallocate(x0)
deallocate(y0)
deallocate(z0)

deallocate(xprev)
deallocate(yprev)
deallocate(zprev)

deallocate(nx)
deallocate(ny)
deallocate(nz)

call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)





end program mc_rmsd
