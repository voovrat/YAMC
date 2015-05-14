Module MoleculeTable !DOC 
!DOC !FILE Molecule table contains the coordinates and orientations of the molecules and the indeces of molecule type by atom types and first and last atom indeces in each molecule
use Molecule
use AtomicData, only : TAtomicData
use composition

   implicit none
   
   Type :: TMoleculeTable  !DOC
!DOC Fields
      integer, dimension(:),allocatable :: mol_type !DOC molecule types
      integer, dimension(:),allocatable :: first_atom,last_atom   !DOC first atom indeces of the molecules
      real(8), dimension(:),allocatable :: x,y,z,theta,phi,psi  !DOC theta phi psi: 3 rotatiions, psi over Oz, theta over Ox, phi over Oz   
      integer :: nalloc !DOC allocated size
      integer :: nmol !DOC number of molecules

      integer :: MaxMolAtom    !DOC  number of atoms in the largest molecule

   END TYPE TMoleculeTable

   
   Type :: TMolTypeSpool !DOC
      integer, dimension(:), allocatable :: limits  !DOC array of n_type elements 
                             !DOC limits for the mol_types
                             !DOC at given position, the molecule type is taken randomly
                             !DOC but the probability to take each molecule should be proportional to
                             !DOC the number of molecule of that type left in the spool
                             !DOC so, we take an interval of length n_left, and divide it to the sub-intervals 
                             !DOC which are equal to the numbers of molecules left
                             !DOC then we take a random number between 1 and n_left, and check, in which interval it is
                             !DOC
                             !DOC limits is an array of n_types numbers, containing the upper boundaries of the intervals
      integer :: n_left !DOC  number of molecules left in the spool
      integer :: n_types !DOC number of molecule types
   END TYPE TMolTypeSpool


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!  MolTypeSpool 

subroutine MolTypeSpool_alloc(this, mol_numbers , n_types ) !DOC
!DOC allocate the MolTypeSpool
!DOC Parameters:
   implicit none
   Type(TMolTypeSpool)   :: this !MolTypeSpool
   integer,dimension(:) :: mol_numbers !Numbers of molecules of each type
   integer :: n_types !DOC number of types
   integer :: i , n_curr

   allocate( this % limits(n_types) )

   this % n_types = n_types 
 
   n_curr = 0
   do i=1,n_types
      n_curr = n_curr + mol_numbers(i)  
      this % limits(i) = n_curr 
   end do

   this % n_left = n_curr
end subroutine

subroutine MolTypeSpool_dealloc(this) !DOC
!DOC deallocate the MolTypeSpool
    implicit none
    Type(TMolTypeSpool) :: this

    deallocate( this % limits )

end subroutine

integer function MolTypeSpool_peekMolecule( this ) !DOC
!DOC Peak the molecule (by chance)
    use MRandom, only : random

    implicit none
!DOC Parameters:
    Type(TMolTypeSpool) :: this !DOC MolTypeSpool
    integer :: N_rand
    integer :: i,j     

    N_rand = random(this % n_left )

    ! find the interval in which N_rand is
    do i=1,this % n_types 
       if ( N_rand .lt. this % limits(i) ) exit
    end do

    MolTypeSpool_peekMolecule = i

    ! update the interval limits after this interval 
    ! we take one molecule from i-th interval, so interval bound (i..n_types) move one molecule left
    do j=i, this % n_types
       this % limits(j) = this % limits(j) - 1
    end do  

    this % n_left = this % n_left - 1

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MoleculeTable 

subroutine MoleculeTable_nulify(this) !DOC
!DOC Set molecule table to the zero state
!DOC Parameters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable
 
   this % nalloc = 0
   this % nmol = 0

end subroutine 

subroutine MoleculeTable_alloc( this, nalloc ) !DOC
!DOC Allocate the molecule table arrays
!DOC Parameters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable
   integer :: nalloc !DOC size of arrays

   allocate( this % x(nalloc)  )
   allocate( this % y(nalloc) )
   allocate( this % z(nalloc) )
   allocate( this % theta(nalloc) )
   allocate( this % phi(nalloc) )
   allocate( this % psi(nalloc) )

   allocate( this % mol_type(nalloc) )
   allocate( this % first_atom(nalloc) )
   allocate( this % last_atom(nalloc) )
  
   this % nalloc = nalloc
   this % nmol = 0

end subroutine

subroutine MoleculeTable_copy( dst ,src ) !DOC
!DOC copy the molecule table
!DOC Prameters:
   Type(TMoleculeTable) :: dst,src !DOC destination and source tables
   integer :: nmol

   nmol = src % nmol
   if (dst % nalloc <  nmol ) then

      if ( dst % nalloc > 0 ) call MoleculeTable_dealloc(dst)
     
      call MoleculeTable_alloc( dst, nmol )

   end if

   dst % nmol = nmol

   dst % x(1:nmol) = src % x(1:nmol)
   dst % y(1:nmol) = src % y(1:nmol)
   dst % z(1:nmol) = src % z(1:nmol)
   
   dst % theta(1:nmol) = src % theta(1:nmol)
   dst % phi(1:nmol) = src % phi(1:nmol)
   dst % psi(1:nmol) = src % psi(1:nmol)

   dst % mol_type(1:nmol) = src % mol_type(1:nmol)
   dst % first_atom(1:nmol) = src % first_atom(1:nmol)
   dst % last_atom(1:nmol) = src % last_atom(1:nmol)

   dst % MaxMolAtom = src % MaxMolAtom

end subroutine


subroutine MoleculeTable_save_binary(this,fid,nbytes_xyz,nbytes_ang) !DOC
!DOC save the MoleculeTable in binary format (see trajectory file)
!DOC Parameters:
   use io, only : write_little_endian
   use constants, only : two_pi
   Type(TMoleculeTable),intent(in) :: this !DOC MoleculeTable
   integer,intent(in) :: fid !DOC file handler
   integer,intent(in) :: nbytes_xyz  !DOC number of bytes per coordinate sampe (fixed point). allowed values: 1,2,3
   integer,intent(in) :: nbytes_ang  !DOC number of bytes per angular sample. allowed values: 1,2,3

   integer :: i
   integer :: xi,yi,zi
   integer :: itheta,iphi,ipsi

   integer :: Pow_xyz,Pow_ang 

   Pow_xyz = 1
   Pow_ang = 1

   do i=1,nbytes_xyz
       Pow_xyz = Pow_xyz * 256
   end do

   do i=1,nbytes_ang
      Pow_ang = Pow_ang * 256 
   end do

   do i=1,this % nmol
      
      xi = nint( (this % x(i) + 0.5 ) * Pow_xyz )
      yi = nint( (this % y(i) + 0.5 ) * Pow_xyz )
      zi = nint( (this % z(i) + 0.5 ) * Pow_xyz )

      itheta = nint( (this % theta(i)) / two_pi * Pow_ang )
      iphi = nint( (this % phi(i)) / two_pi * Pow_ang )
      ipsi = nint( (this % psi(i)) / two_pi * Pow_ang )  

      call write_little_endian(fid, xi, nbytes_xyz )
      call write_little_endian(fid, yi, nbytes_xyz )
      call write_little_endian(fid, zi, nbytes_xyz )

      call write_little_endian( fid, itheta, nbytes_ang )
      call write_little_endian( fid, iphi, nbytes_ang ) 
      call write_little_endian( fid, ipsi, nbytes_ang )   

   end do
   
end subroutine

subroutine MoleculeTable_save_text(this,fid,BoxLength) !DOC 
!DOC Save the molecule table in the text format
   use io, only : write_little_endian
   use constants, only : two_pi
!DOC Parameters
   Type(TMoleculeTable),intent(in) :: this !DOC MoleculeTable
   integer,intent(in) :: fid !DOC file handler
   real(8),intent(in) :: BoxLength !DOC BoxLength

   integer :: i

   do i=1,this % nmol
      
      write(fid,*)  this % x(i) * BoxLength, this % y(i) * BoxLength, this % z(i) * BoxLength, &
                    this % theta(i), this % phi(i), this % psi(i) 
   end do
   
end subroutine



subroutine MoleculeTable_load_binary(this,fid,comp, nbytes_xyz,nbytes_ang) !DOC 
!DOC load molecule table from the binary format
   use io, only : read_little_endian 
   use constants, only : two_pi

!DOC Parameters:
  Type(TMoleculeTable) :: this !DOC MoleculeTable
   integer,intent(in) :: fid !DOC file handler
   Type(TComposition),intent(in) :: comp !DOC composition

   integer,intent(in) :: nbytes_xyz  !DOC number of bytes per coordinate sampe (fixed point). allowed values: 1,2,3
   integer,intent(in) :: nbytes_ang  !DOC number of bytes per angular sample. allowed values: 1,2,3

   integer :: i,ityp,ii
   integer :: xi,yi,zi
   integer :: itheta,iphi,ipsi

   real(8) :: Pow_xyz,Pow_ang 

   integer :: ntotal

   ntotal = SUM( comp % mol_numbers(1:comp % n_types) )

   if( this % nalloc < ntotal ) then 

      if ( this % nalloc > 0) then
          call MoleculeTable_dealloc(this )
      end if 

      call MoleculeTable_alloc( this, ntotal )

      this % nalloc = ntotal
  
   end if 

   this % nmol = ntotal

   ii=1
   do ityp=1,comp % n_types 

      do i=1,comp % mol_numbers(ityp)

         this % mol_type(ii) = comp % mol_types( ityp )
         ii = ii + 1        

      end do
   end do

   Pow_xyz = 1
   Pow_ang = 1

   do i=1,nbytes_xyz
       Pow_xyz = Pow_xyz * 256
   end do

   do i=1,nbytes_ang
      Pow_ang = Pow_ang * 256 
   end do

   do i=1,this % nmol
 
      call read_little_endian(fid, nbytes_xyz, xi )
      call read_little_endian(fid, nbytes_xyz, yi )
      call read_little_endian(fid, nbytes_xyz, zi )

      call read_little_endian( fid,  nbytes_ang, itheta )
      call read_little_endian( fid,  nbytes_ang, iphi ) 
      call read_little_endian( fid,  nbytes_ang, ipsi )   
     
      this % x(i) = xi / Pow_xyz - 0.5
      this % y(i) = yi / Pow_xyz - 0.5
      this % z(i) = zi / Pow_xyz - 0.5

      this % theta(i) = itheta / Pow_ang * two_pi
      this % phi(i) = iphi / Pow_ang * two_pi
      this % psi(i) = ipsi / Pow_ang * two_pi

   end do
 
end subroutine

subroutine MoleculeTable_load_text(this,fid,comp,BoxLength) !DOC
!DOC Load the MoleculeTable from the text file
!DOC Parameters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable
   integer,intent(in) :: fid !DOC file handler
   Type(TComposition),intent(in) :: comp !DOC composition

   real(8),intent(in) :: BoxLength !DOC BoxLength
   integer :: ntotal
   integer :: ityp
   integer :: i,ii 

   ntotal = SUM( comp % mol_numbers(1:comp % n_types) )

   if( this % nalloc < ntotal ) then 

      if ( this % nalloc > 0) then
          call MoleculeTable_dealloc(this )
      end if 

      call MoleculeTable_alloc( this, ntotal )

      this % nalloc = ntotal
  
   end if 

   this % nmol = ntotal

   write(*,*) 'Load Text: nmol',this % nmol

   ii=1
   do ityp=1,comp % n_types 

      do i=1,comp % mol_numbers(ityp)

         this % mol_type(ii) = comp % mol_types( ityp )
         ii = ii + 1        

      end do
   end do

   
   do i=1,this % nmol
 

      read(fid,*)  this % x(i), this % y(i), this % z(i), this % theta(i), this % phi(i), this % psi(i)

      this % x(i) = this % x(i) / BoxLength 
      this % y(i) = this % y(i) / BoxLength
      this % z(i) = this % z(i) / BoxLength 

   end do
 
end subroutine

subroutine MoleculeTable_dealloc(this) !DOC
!DOC deallocate the molecule table
!DOC Parmeters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable

   deallocate( this % x)
   deallocate( this % y)
   deallocate( this % z )
   deallocate( this % theta )
   deallocate( this % phi )
   deallocate( this % psi )
 
   deallocate( this % mol_type )
   deallocate( this % first_atom )
   deallocate( this % last_atom )

   this % nalloc = 0
   this % nmol = 0
end subroutine

subroutine  MoleculeTable_placeCube(this,x_shft,y_shft,z_shft, i, m) !DOC
      !DOC places m particles in a 1x1x1 cube  
      implicit none
!DOC Parameters
      Type(TMoleculeTable) :: this !DOC MoleculeTable
      real(8) :: x_shft,y_shft,z_shft   !DOC dispacement of the cube
      integer,intent(inout) :: i        !DOC  Total particle counter. Updated in the subroutine .
      integer :: m  !DOC m^3 particles have to be placed

      real(8) :: dx  ! dx = 1/m
      integer :: ix,iy,iz

      if ( i>this % nmol) return

      dx = 1.0/m

      do ix=1,m
      do iy=1,m
      do iz=1,m

          this % x(i) = -0.5 + x_shft + (ix-1) * dx
          this % y(i) = -0.5 + y_shft + (iy-1) * dx
          this % z(i) = -0.5 + z_shft + (iz-1) * dx

          i=i+1
          if (i>this % nmol) return  ! exceed the total number to be placed in all cubes

      end do
      end do
      end do

end subroutine


subroutine MoleculeTable_placeMoleculesToGrid(this, mol_types,  mol_numbers, n_types ) !DOC 
 !DOC place molecules into the MoleculeTable
   use MRandom, only :  rand
   use constants, only : pi, two_pi
  
   implicit none
!DOC Parameters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable
   integer,dimension(:) :: mol_types    !DOC types of the molecules. Array of n_types elements
                           !DOC (actually, these should be their handlers, loaded to memory
                           !DOC  this function does not use this, but to convrt them to atom coors it is necessary)
   integer,dimension(:) :: mol_numbers  !DOC array of n_types elements. How many molecules of each type will be placed into the box  
   integer :: n_types  !DOC number of mol_types and mol_numbers

!VARIABLES
   real(8) :: dx ! grid step
   real(8) :: O = 0.0

   real(8) :: xn  ! real(n) 
   integer :: m ! number of molecules in one sub-box in one direction
   integer :: i,ii
!   real(8) :: vx1,vy1,vz1,vx2,vy2,vz2 ! hydrogen atom positions in the intramolecular coordinates
   Type(TMolTypeSpool) :: spool  ! to check types of molecules randomly
   integer :: n_total ! total number of molecules
   integer :: interval ! current peeked interval

   real(8),dimension(:),allocatable :: xx,yy,zz

 ! Now: we want to put the molecules of all times sequentially
 ! to do this we create the list of counters - current position of the last molecule of each type
 ! counter(i) is the current position of the next molecule of type i to be places
   integer,dimension(:),allocatable :: counters

! calculate total number of molecules

  allocate(counters(n_types) )
 

  counters(1) = 1

  do i=1,n_types-1

     counters(i+1) = counters(i) + mol_numbers(i)

  end do

  call MolTypeSpool_alloc(spool, mol_numbers, n_types) 
 
   n_total = spool % n_left 

   write(*,*) 'n_total=',n_total

   if( this % nalloc .lt. n_total ) then
  
      if( this % nalloc > 0) then 
         write(*,*) 'moltable dealloc'
         call MoleculeTable_dealloc(this)
      end if
      write(*,*) 'moltable alloc ',n_total
      call MoleculeTable_alloc(this,n_total)
   end if

  allocate(xx(n_total))
  allocate(yy(n_total))
  allocate(zz(n_total))

   xn=n_total
    m=(xn/4.)**(1./3.)+0.999
    dx = 1./m

! m * dx = 1
!
!      dx
!   |------|
!     
!     
!   1--2--1--2--1
!   |..   |     |
!   4--3--4--3--4
!   |     |     |
!   1--2- 1--2--1 
!   |     |     |
!   4--3--4--3--4
!   |     |     |
!   1--2--1--2--1
!
!    z=0  z=0.5/m z=1/m  z=1.5/m .... 
!   12121         12121
!          43434          43434
!   12121         12121
!          43434          43434
!   12121         12121 
!
!   1,2,3,4 - cycles below
!
!   1,3 are in the same plane
!   2,4 are at dz =  dx/2  above
!
   O = 0.0  ! but real(8)! 
   i=1   
   this % nmol = n_total  ! used in the subroutine placeCube
                                      !x_shft y_shft z_shft   counter(inout)   cube_length
   call MoleculeTable_placeCube( this,    O,       O,      O,        i,          m)
   call MoleculeTable_placeCube( this,   dx/2,     O,    dx/2,        i,          m)
   call MoleculeTable_placeCube( this,   dx/2, dx/2,       O,        i,          m)
   call MoleculeTable_placeCube( this,     O, dx/2,    dx/2,        i,          m)


 ! random orientations. Just stored in the molecule array + randomly chosen types 

  do i=1,n_total
     interval = MolTypeSpool_peekMolecule( spool )
 
     ii = counters(interval)

   !  write(*,*) interval,ii
 
     this % mol_type( ii ) = mol_types(interval)

     this % theta( ii) = rand() * pi
     this % phi(ii ) = rand() * two_pi
     this % psi(ii ) = rand() * two_pi

     xx(ii) = this % x(i)
     yy(ii) = this % y(i)
     zz(ii) = this % z(i)

     counters(interval) = counters(interval) + 1
  end do

  this % x(1:n_total) = xx(1:n_total)
  this % y(1:n_total) = yy(1:n_total)
  this % z(1:n_total) = zz(1:n_total)
  

  call MolTypeSpool_dealloc(spool) 

  deallocate(counters)

  deallocate(xx)
  deallocate(yy)
  deallocate(zz)

end subroutine



!   Type :: TAtomicData
!       real(8), dimension(:), allocatable :: xx,yy,zz
!       real(8), dimension(:), allocatable :: sigma,epsilon,charge
!       character(4),dimension(:),allocatable :: atomnames
!       integer,dimension(:),allocatable :: molnum_by_atomnum
!       integer :: natom, nalloc
!   End Type TAtomicData
subroutine MoleculeTable_fillAtomicData( this, atomic_data, BoxLength, kT_kcal_mol ) !DOC
!DOC convert the molecule coordinates in molecule table to the atom coordinates in the AtomicData
!DOC atomicData should be pre-allocated
!DOC atomnames are optional and filled only if fillAtomNames=.TRUE. AND corresponding molecules have atomnames 
   use error
   use geometry, only : rot_vect
   use MoleculeHandler
   use Molecule 
!   use parameters, only : BoxLength, kT_kcal_mol 
 ! BoxLength is used to convert molecular units (anstroems) to internal units (BoxLengths)
 ! k_B*T, for conversion from molecular units (usually kcal/mol) to internal units 

   implicit none
! PARAMETERS
   Type(TMoleculeTable) :: this !DOC MoleculeTable
   Type(TAtomicData) :: atomic_data  !DOC AtomicData (output)
 ! ***********************************************
                                     !  DO NOT USE intent(out) WITH STRUCTURES!!!!
                                     ! IT CORRUPTS THE DATA ~!!!!!!
                                     ! ***********************************************

!   real(8),intent(in)    :: kT           ! k_B*T, for conversion from molecular units (usually kcal/mol) to internal units 

   real(8),intent(in) :: BoxLength    !DOC can be taken from parameters. But if not - it is parameters independent
   real(8),intent(in) :: kT_kcal_mol !DOC kt in kcal/mol. sometimes can be set to 1, if for example, only coordinates matter and epsilon is irrelevant

! VARIABLES
   integer :: imol, imolatom
   Type(TMolecule),pointer ::  mol_ptr     ! pointer to the molecule  of given type returned by MoleculeHandler
   integer :: n_mol_atoms
   real(8) :: cos_theta,sin_theta,cos_phi,sin_phi,cos_psi,sin_psi    ! sin and cos of theta,phi,psi (sequentially), use by rot_vect
   real(8) :: xmolatom,ymolatom,zmolatom ! atom coordinates in molecular coordinates
   !real(8)    :: BoxLength    ! BoxLength is used to convert molecular units (anstroems) to internal units (BoxLengths)
   integer :: natom ! number of atoms 
 
   !BoxLength = atomic_data % BoxLength
   atomic_data % BoxLength = BoxLength

   this % MaxMolAtom = 0

   natom = 0
   do imol=1,this % nmol

      call MoleculeHandler_getMolecule( this % mol_type(imol), mol_ptr )
   
      if ( natom + mol_ptr % Natoms > atomic_data % nalloc ) then

          !write(*,*) 'natoms=',mol_ptr % Natoms, ' nalloc=', atomic_data % nalloc

          write(error_message,*) 'MoleculeTable_fillAtomicData: cannot place molecule #',imol,': atom array limit exceeded'
          call error_throw(ERROR_LIMITS)
          return
      end if    

      if ( this % MaxMolAtom < mol_ptr % Natoms ) then

         this % MaxMolAtom = mol_ptr % Natoms
 
      end if 

      this % first_atom(imol) = natom + 1
      this % last_atom(imol) = natom +mol_ptr % Natoms

      cos_theta = cos(this % theta(imol) )
      sin_theta = sin(this % theta(imol) )

      cos_phi = cos(this % phi(imol) )
      sin_phi = sin(this % phi(imol) )

      cos_psi = cos(this % psi(imol) )
      sin_psi = sin(this % psi(imol) )
   
      do imolatom = 1,mol_ptr % Natoms

          ! Rotate in a following direction: by psi over Oz, by theta over Oy, by phi over Oz again

          !  by psi over Oz
          call rot_vect(mol_ptr % x(imolatom), mol_ptr % y(imolatom), mol_ptr % z(imolatom), &
                       xmolatom, ymolatom, zmolatom, &
                       cos_psi,sin_psi,3) 
          
          ! by theta over Oy
          call rot_vect(xmolatom, ymolatom, zmolatom, &
                       xmolatom, ymolatom, zmolatom, &
                       cos_theta,sin_theta,2)  ! note: subroutine works well even when xin and xout are the same variables
          
          ! by phi over Oz again
          call rot_vect(xmolatom, ymolatom, zmolatom, &
                        xmolatom, ymolatom, zmolatom, &
                        cos_phi,sin_phi,3)  ! note: subroutine works well even when xin and xout are the same variables
 
          ! convert from angstroems to BoxLengths and place to the arrays
          natom = natom + 1

          atomic_data % xx(natom) = xmolatom / BoxLength + this % x(imol)
          atomic_data % yy(natom) = ymolatom / BoxLength + this % y(imol)
          atomic_data % zz(natom) = zmolatom / BoxLength + this % z(imol)

          ! get charges
          atomic_data % charge(natom) = mol_ptr % charge(imolatom)
          ! convert sigma to internal units
          atomic_data % sigma(natom) = mol_ptr % sigma(imolatom) / BoxLength

          atomic_data % hard_core_angstr(natom) = mol_ptr % hard_core(imolatom)

!          write(*,*) 'natom:',natom,'imolatom:',imolatom,'mol_ptr.eps(imolatom)=',mol_ptr % epsilon(imolatom),'kT',kT_kcal_mol
          atomic_data % epsilon(natom) = mol_ptr % epsilon(imolatom) / kT_kcal_mol

          ! save atomnames if any
          if ( mol_ptr % hasNames .and. atomic_data % hasAtomNames ) then
             atomic_data % atomnames(natom) = mol_ptr % atomnames(imolatom)
          end if 

          atomic_data % molnum_by_atomnum(natom) = imol
 
      end do ! imolatom
      
   end do ! imol

   atomic_data % natom = natom

end subroutine

subroutine MoleculeTable_calcPositionsOrientations( this, xx, yy, zz, BoxLength ) !DOC
!DOC Extract the positions and the orientations of the molecules from the coordinates of the atoms
   use MoleculeHandler, only : MoleculeHandler_getMolecule
   use geometry, only : xyz_to_angles,center_of_mass
   ! re-calculate the molecule positions and orientations using the atomic coordinates
   ! the sizes of xx,yy,zz should correspond to the current MoleculeTable
   implicit none
!DOC Parameters:
   Type(TMoleculeTable) :: this !DOC MoleculeTable
   real(8),dimension(:),intent(in),target :: xx,yy,zz !DOC coordinates of the atoms
   real(8),intent(in) :: BoxLength !DOC box length

   integer :: i 
   Type(TMolecule),pointer :: mol_ptr
   real(8) :: theta, phi, psi   
   integer :: first_atom,last_atom,natom
   real(8),dimension(:),allocatable :: xx_new,yy_new,zz_new
   real(8) :: x_center,y_center,z_center

   allocate(xx_new(this % MaxMolAtom ) )
   allocate(yy_new(this % MaxMolAtom ) )
   allocate(zz_new(this % MaxMolAtom ) )

  ! write(*,*) 'nmol=',this % nmol

   do i=1,this % nmol
      
      call MoleculeHandler_getMolecule( this % mol_type(i), mol_ptr )
      
      first_atom = this % first_atom(i)
      last_atom = this % last_atom(i)
      natom = mol_ptr % Natoms

   !   write(*,*) 'first_atom:',first_atom,'last_atom',last_atom

      xx_new(1:natom) = xx(first_atom:last_atom)
      yy_new(1:natom) = yy(first_atom:last_atom)
      zz_new(1:natom) = zz(first_atom:last_atom)

      call center_of_mass(xx_new,yy_new,zz_new, mol_ptr % mass, natom,x_center,y_center,z_center )

      xx_new(1:natom) = ( xx_new(1:natom) - x_center ) * BoxLength 
      yy_new(1:natom) = ( yy_new(1:natom) - y_center ) * BoxLength
      zz_new(1:natom) = ( zz_new(1:natom) - z_center ) * BoxLength

      call xyz_to_angles( mol_ptr % x, mol_ptr % y, mol_ptr % z, &
                          xx_new,yy_new,zz_new, natom, theta, phi, psi )

      this % x(i) = x_center
      this % y(i) = y_center 
      this % z(i) = z_center

      this % theta(i) = theta
      this % phi(i) = phi
      this % psi(i) = psi
 
   end do

   deallocate( xx_new )
   deallocate( yy_new )
   deallocate( zz_new ) 

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module MoleculeTable
