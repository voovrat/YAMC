program mc_mean_force !DOC !PROG
!DOC !PROG !FILE calculate the mean forces for specific molecules
use parameters , only : BoxLength, traj_file, frames_file, Parameters_init, Parameters_recalc
use parameters , only : kT_kcal_mol, output_nbytes_xyz, output_nbytes_ang
use MoleculeTable
use AtomicData
use geometry
use string
!use LJTypes

!DOC !PROG mc_mean_force: calculate the mean force projection between the molecules'
!DOC !PROG
!DOC !PROG Usage:  mc_mean_force  parameters.prm system.composition forces.ftraj  rangeA rangeB maxR dr output.dat
!DOC !PROG Arguments:
!DOC !PROG ::  system.composition - number and types of molecules in the system'
!DOC !PROG :: extra_parameters - coma separated  prm1=val1,prm2=val2,...   NO SPACES ALLOWED! '
!DOC !PROG :: forces.ftraj - forces trajectory file. In float (24+8) format. See FloatingPoint.f90 for details'
!DOC !PROG  :: rangeA,rangeB  - ranges for the 1st and 2nd molecules'
!DOC !PROG                     in format  num1[-num2]'
!DOC !PROG                     the MeanForce will be calculated for all pairs A-B where A is in range A, B is in rangeB'
!DOC !PROG  :: maxR, dr   -  samples for MeanForce will be [0:dr:maxR]'
!DOC !PROG  :: output.dat - four  columns:  r  sum(f12)  N(r) sum(f12)/N(r) where N(r) is number of AB pairs at distance r'
!DOC !PROG

!use EwaldSumKSpace
!use EwaldSumRealSpace
!use EwaldSumExternal
!use RhoSquared
!use SumSinCosKR

!use MC , only : initmc, movemc,clearmc,RealSpace_TotalSum,KSpace_TotalSum,mc_calc_total_sums_dbg,mc_calc_total_sums
!use MC , only : U_tail_lj6,U_tail_lj12,U_tail_coulomb, U_total_ew,U_total_lj,U_total,U_ext
!use io, only : write_integer_array,write_real_array,write_xyz_array,io_open,io_close,write_real_matrix
!use MCLuc, only : initmc_Luc,movemc_Luc,MCLuc_init_input,MCLucInput,piston ,write_energy_LUC
!use runmc, only : run_mc,ifirst
!use MRandom

implicit none

character(SYSTEM_STRING_LENGTH) :: parameters_file,composition_file,forces_file,output_file
character(SYSTEM_STRING_LENGTH) :: tmpstr
character(SYSTEM_STRING_LENGTH) :: rangeA,rangeB   ! ranges for the 1st and 2nd molecules
                                                   ! in format  num1[-num2]
                                                   ! this means that the MeanForce will be calculated for all pairs A-B where A is in range A, B is in rangeB

integer :: A_first,A_last, B_first,B_last ! first and last molecule indeces for the ranges

real(8) :: maxR, dr
integer :: nsampl ! AINT(maxR/dr)

integer,dimension(:),allocatable :: NN
real(8),dimension(:),allocatable :: FF

character(SYSTEM_STRING_LENGTH) :: extra_parameters_string
logical :: extra_parameters = .FALSE.

Type(TComposition) :: comp
Type(TMoleculeTable) :: mol_tab
Type(TAtomicData),target :: atomic_data
!Type(TLJTypes) :: lj_types

integer :: i,j
integer :: hfile,nk
integer :: MaxMolecule,MaxAtom
integer :: hfram,stat,htraj,hforce
integer :: nfram

!Type(TEwaldSumKSpace),target :: kspace_sum
!Type(TEwaldSumRealSpace) :: real_sum
!Type(TEwaldSumExternal) :: ext_sum

real(8),dimension(:),allocatable :: fx,fy,fz

real(8),dimension(:),allocatable :: fxmol,fymol,fzmol


integer,dimension(:),allocatable :: molnum_by_atomnum,first_atom,last_atom
real(8),dimension(:),allocatable,target :: atom_mass ! masses of the atoms


!real(8) :: fx_tot,fy_tot,fz_tot  ! total force (the sum of all forces) 
                                 ! is subtracted before applying the virial equation

!real(8) :: x_center,y_center,z_center ! center of mass of the box



call parse_command_line
call read_input_files
!call allocate_ewald

allocate( fx(MaxAtom) )
allocate( fy(MaxAtom) )
allocate( fz(MaxAtom) )


allocate(fxmol(MaxMolecule))
allocate(fymol(MaxMolecule))
allocate(fzmol(MaxMolecule))

allocate(molnum_by_atomnum(MaxAtom))
allocate(first_atom(MaxMolecule))
allocate(last_atom(MaxMolecule))
allocate(atom_mass(MaxAtom))


nsampl = aint( maxR / dr) + 1

allocate( NN(nsampl) )
allocate( FF(nsampl) )

FF(1:nsampl) = 0.d0
NN(1:nsampl) = 0


call fill_atom_indeces(comp,molnum_by_atomnum,first_atom,last_atom,atom_mass)

hforce = io_open(forces_file,'b')

! first frame is already read
nfram = 1
do while ( .TRUE. )
 
  write(*,*) 'Fram:',nfram
 
   call read_forces(hforce,fx,fy,fz,MaxAtom)
   call calc_mol_forces(comp,fx,fy,fz,fxmol,fymol,fzmol)
   call calc_mean_force

   call load_frame
   
   if( stat /= 0) exit


  nfram = nfram + 1

end do

call io_close(hforce)

call save_results



call deallocate_all


contains 


subroutine calc_mean_force !DOC
!DOC calculate the mean force at the given frame
  integer :: i

  integer :: first1, last1, natom1, first2, last2, natom2
  real(8),dimension(:),pointer :: xx1, yy1, zz1, mass1, xx2, yy2, zz2, mass2
  real(8) :: x_center1, y_center1, z_center1, x_center2, y_center2, z_center2
  real(8) :: dx,dy,dz, halfBox
  real(8) :: fx1,fy1,fz1,fx2,fy2,fz2,dfx,dfy,dfz, F12
  real(8) :: R12 
 
  integer :: offset
  

  halfBox = 0.5*BoxLength

  do i=A_first,A_last

     first1 = first_atom(i)
     last1 = last_atom(i)

     xx1 => atomic_data % xx(first1:last1)
     yy1 => atomic_data % yy(first1:last1)
     zz1 => atomic_data % zz(first1:last1)    
     mass1 => atom_mass(first1:last1)
     natom1 = last1 - first1 + 1

     call center_of_mass(xx1,yy1,zz1,mass1,natom1, x_center1, y_center1, z_center1)
         
     x_center1 = x_center1 * BoxLength
     y_center1 = y_center1 * BoxLength
     z_center1 = z_center1 * BoxLength

     fx1 = fxmol(i)
     fy1 = fymol(i)
     fz1 = fzmol(i)


     do j=B_first,B_last

        if( i==j ) cycle
 
        first2 = first_atom(j)
        last2 = last_atom(j)
   
        xx2 => atomic_data % xx(first2:last2)
        yy2 => atomic_data % yy(first2:last2)
        zz2 => atomic_data % zz(first2:last2)    
        mass2 => atom_mass(first2:last2)
        natom2 = last2 - first2 + 1
   
        call center_of_mass(xx2,yy2,zz2,mass2,natom2, x_center2, y_center2, z_center2)
            
        x_center2 = x_center2 * BoxLength
        y_center2 = y_center2 * BoxLength
        z_center2 = z_center2 * BoxLength
   
        
        dx = x_center2 - x_center1
        dy = y_center2 - y_center1
        dz = z_center2 - z_center1

        if( dx < -halfBox ) dx = dx + BoxLength
        if( dx > halfBox )  dx = dx - BoxLength

        if( dy < -halfBox ) dy = dy + BoxLength
        if( dy > halfBox )  dy = dy - BoxLength
       
        if( dz < -halfBox ) dz = dz + BoxLength
        if( dz > halfBox )  dz = dz - BoxLength

        ! now:  (F2 - F1)  * (r2-r1) / | r2 - r1| 

        R12 = sqrt(dx**2 + dy**2 + dz**2 );

        if( R12 >= maxR ) cycle

        fx2 = fxmol(j)
        fy2 = fymol(j)
        fz2 = fzmol(j)        

        dfx = fx2 - fx1
        dfy = fy2 - fy1
        dfz = fz2 - fz1
         
         ! projection to the r12
        F12 = (dfx * dx + dfy * dy + dfz * dz ) / R12

        offset = nint(R12 / dr ) + 1

        if ( offset > nsampl) then
           write(*,*) 'offset',offset,'nsampl',nsampl,'r',R12
           stop
        end if 

        NN(offset) = NN(offset) + 1
        FF(offset) = FF(offset) + F12

!        write(*,*) 'offset',offset,'NN',NN(offset),'F12',F12,'FF',FF(offset)

     end do
  end do


end subroutine

subroutine save_results !DOC 
!DOC write the results to file

   integer :: hout
   integer :: i
   integer :: offset

   real(8) :: r

   hout = io_open(output_file,'w')
   r = 0
   do while( r < maxR ) 
        
      offset = nint(r/dr) + 1

      if( NN(offset) > 0) then
         write(hout,*)  r,FF(offset),NN(offset), FF(offset)/NN(offset) 
      else
         write(hout,*)  r,FF(offset), NN(offset), 0.d0
      end if

      r = r + dr
   end do


   call io_close(hout)

end subroutine

subroutine parse_command_line !DOC 
!DOC read the command line arguments
   use string

   integer :: separator_pos


if ( iargc() < 8 ) then

   write(*,*)  '  mc_mean_force: calculate the mean force projection between the molecules'
   write(*,*)
   write(*,*)  'Usage:  mc_mean_force  parameters.prm system.composition forces.ftraj  rangeA rangeB maxR dr output.dat'
   write(*,*)  '   system.composition - number and types of molecules in the system'
   write(*,*)  '   extra_parameters - coma separated  prm1=val1,prm2=val2,...   NO SPACES ALLOWED! '
   write(*,*)  '   forces.ftraj - forces trajectory file. In float (24+8) format. See FloatingPoint.f90 for details'
   write(*,*)  '   rangeA, rangeB  - ranges for the 1st and 2nd molecules'
   write(*,*)  '                      in format  num1[-num2]'
   write(*,*)  '                      the MeanForce will be calculated for all pairs A-B where A is in range A, B is in rangeB'
   write(*,*)  '   maxR, dr   -  samples for MeanForce will be [0:dr:maxR]'
   write(*,*)  '   output.dat - four  columns:  r  sum(f12)  N(r) sum(f12)/N(r) where N(r) is number of AB pairs at distance r'

   stop
end if

call getarg(1,parameters_file)
call getarg(2,composition_file)
call getarg(3,forces_file)

call getarg(4,rangeA)
call read_range(rangeA,A_first,A_last)

call getarg(5,rangeB)
call read_range(rangeB,B_first,B_last)

call getarg(6,tmpstr)
read(tmpstr,*) maxR

call getarg(7,tmpstr)
read(tmpstr,*) dr
call getarg(8,output_file)

    extra_parameters = .FALSE.

end subroutine

subroutine read_range(str,first,last) !DOC
!DOC extract the first and last ranges from the string of format first-last
   use string, only : str_get_next_pos
   character(SYSTEM_STRING_LENGTH),intent(in) :: str
   integer :: first,last
 
   integer :: sep_pos 
  
   sep_pos = str_get_next_pos(str,1,'-')
   if (sep_pos .gt. 0) then
     ! from-to
!     write(*,*) str(1:20),sep_pos,'[',str(1:sep_pos-1),']',str(sep_pos+1:SYSTEM_STRING_LENGTH),']'

      read(str(1:sep_pos-1),*) first
     read(str(sep_pos+1:SYSTEM_STRING_LENGTH),*) last  

  !  write(*,*) 'first',first,'last',last

   else ! one number, from=to
  
     read(str,*) first
     last = first
  
   end if
  
end subroutine


subroutine read_input_files !DOC 
!DOC read input files
   use runmc, only : read_traj_frame
   use SystemSettings, only : SEEK_SET
   use io, only : io_open,io_close


   integer :: nframes
   integer :: block_size


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

hfram = io_open(frames_file,'r')

! read the frames file 
call read_traj_frame(hfram,stat)
! recalc parameters for the new boxlength
call Parameters_recalc(BoxLength)

! read the molecule table
call MoleculeTable_nulify(mol_tab)
htraj = io_open(traj_file,'b') 
call MoleculeTable_load_binary( mol_tab, htraj, comp, output_nbytes_xyz, output_nbytes_ang )
!call io_close(hfile)

! initialize the atomic_data
call  AtomicData_alloc(atomic_data, MaxAtom, BoxLength, .TRUE. )

! generate the atomic data
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
!call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )

end subroutine



subroutine load_frame !DOC
!DOC load next frame from file
   use runmc, only : read_traj_frame

real(8) :: old_box_length
!hfram = io_open(frames_file,'r')

old_box_length = BoxLength

! read the frames file 
call read_traj_frame(hfram,stat)
! recalc parameters for the new boxlength

if (stat /=0 ) return


call Parameters_recalc(BoxLength)

! read the molecule table
!call MoleculeTable_nulify(mol_tab)
!htraj = io_open(traj_file,'b') 
call MoleculeTable_load_binary( mol_tab, htraj, comp, output_nbytes_xyz, output_nbytes_ang )
call MoleculeTable_fillAtomicData( mol_tab, atomic_data, BoxLength, kT_kcal_mol )

! prepare lj_types arrays
!call LJTypes_fill(lj_types , atomic_data % sigma , atomic_data % epsilon, MaxAtom )
!call LJTypes_scale_sigma(lj_types, old_box_length / BoxLength )

end subroutine

subroutine read_forces(hforce,fx,fy,fz,MaxAtom) !DOC
!DOC read the forces from the forces trajectory file
!DOC Parameters:
   use FloatingPoint, only : read_float
 
   integer,intent(in) :: hforce !DOC file handler
   real(8),dimension(:)  :: fx,fy,fz !DOC output: forces for each atom
    integer,intent(in) :: MaxAtom !DOC number of atoms

   integer :: i
  

   do i=1,MaxAtom

      fx(i) = read_float(hforce)  ! forces are in kT/Angstr 
      fy(i) = read_float(hforce)
      fz(i) = read_float(hforce)
  

   end do

end subroutine

subroutine fill_atom_indeces(comp,molnum_by_atomnum,first_atom,last_atom,atom_mass) !DOC
!DOC fill molnum_by_atomnum, first_atom, last_atom, atom_mass arrays, using the composition
!DOC Parameters: 
    use MoleculeHandler 
    use Molecule

    Type(TComposition),intent(in) :: comp !DOC composition (input)
    integer,dimension(:) :: molnum_by_atomnum !DOC molecule index by atom index (output)
    integer,dimension(:) ::  first_atom,last_atom  !DOC first and last atom indeces for each molecule (output) 
    real(8),dimension(:) :: atom_mass !DOC mass of each atom (output)
    


    integer :: i,j,mol_typ,nmolatom , nmol, imol
    integer :: first,last
   
    Type(TMolecule),pointer :: mol_ptr  


    imol=1
    first=1
    do i=1,comp % n_types 

        mol_typ = comp % mol_types(i)
        call MoleculeHandler_getMolecule(mol_typ,mol_ptr)
          
        nmolatom = mol_ptr % NAtoms

        nmol = comp % mol_numbers(i)

        do j=1,nmol
            last = first + nmolatom-1

            molnum_by_atomnum(first:last) = imol
            first_atom(imol) = first
            last_atom(imol) = last
            atom_mass(first:last)  = mol_ptr % mass(1:last-first+1)

  !           write(*,*) 'first',first,'last',last,'MaxAtom',MaxAtom,'MaxMolecule',MaxMolecule,'imol',imol 
 
 !           fxmol(imol) = SUM(fx(first:last))
 !           fymol(imol) = SUM(fy(first:last))
 !           fzmol(imol) = SUM(fz(first:last))

            imol = imol + 1
            first = first + nmolatom
        end do 

   end do    


end subroutine


subroutine calc_mol_forces(comp,fx,fy,fz,fxmol,fymol,fzmol) !DOC
!DOC calculate the forces acting on each molecule
    use MoleculeHandler 
    use Molecule

!DOC Parameters:
    Type(TComposition),intent(in) :: comp !DOC composition
    real(8),dimension(:),intent(in) :: fx,fy,fz !DOC atomic forces
    real(8),dimension(:) :: fxmol,fymol,fzmol !DOC molecular forces (output)

    integer :: i,j,mol_typ,nmolatom , nmol, imol
    integer :: first,last
   
    Type(TMolecule),pointer :: mol_ptr  


    imol=1
    first=1
    do i=1,comp % n_types 

        mol_typ = comp % mol_types(i)
        call MoleculeHandler_getMolecule(mol_typ,mol_ptr)
          
        nmolatom = mol_ptr % NAtoms

        nmol = comp % mol_numbers(i)

        do j=1,nmol
            last = first + nmolatom-1

            fxmol(imol) = SUM(fx(first:last))
            fymol(imol) = SUM(fy(first:last))
            fzmol(imol) = SUM(fz(first:last))

            imol = imol + 1
            first = first + nmolatom
        end do 

   end do    


end subroutine

subroutine deallocate_all !DOC
!DOC deallocate everything


!call deallocate_ewald

deallocate( fx )
deallocate( fy )
deallocate( fz )

deallocate( fxmol )
deallocate( fymol )
deallocate( fzmol )

deallocate( NN )
deallocate( FF )

deallocate( molnum_by_atomnum )
deallocate( first_atom )
deallocate( last_atom )
deallocate( atom_mass )

call AtomicData_dealloc(atomic_data)
call MoleculeTable_dealloc(mol_tab)

call Composition_dealloc_molecules(comp)
call Composition_dealloc(comp)

!call LJTypes_dealloc(lj_types)


end subroutine


end program mc_mean_force
