Module MRandom !DOC
!DOC !FILE functions to work with random numbers random number

   implicit none

   integer, dimension(:), allocatable :: MRandom_seed_array
   integer :: MRandom_seed_size = 0
   real :: MRandom_current_random_number

   integer,parameter :: RND_STANDARD = 1
   integer,parameter :: RND_CONSTANT = 2
   real(8) :: rnd_constant_const = 0.5
   real(8) :: MRandom_rand_epsilon = 1d-14

   integer,parameter :: RND_FILE = 3
   integer :: rnd_file_hfile = 0  ! to be set 

   integer :: random_method = RND_STANDARD
 
     

 
CONTAINS

   subroutine randomize(seed) !DOC
!DOC initialize the random number generator with seed
      integer :: seed
      integer :: i      

  !    write(*,*) 'randomize(',seed,'), seed_size=',MRandom_seed_size

      if ( MRandom_seed_size.gt.0 ) then
!         write(*,*) 'Deallocating the last array'
         deallocate(MRandom_seed_array)
      end if

      call random_seed(size=MRandom_seed_size)

   !   write(*,*) 'New seed_size=',MRandom_seed_size 

      allocate(MRandom_seed_array(MRandom_seed_size))

    !  write(*,*) 'allocating new array'

       do i=1,MRandom_seed_size
          MRandom_seed_array(i) = seed;
     !     write(*,*) 'array(',i,')=',seed
       end do

       call random_seed(put=MRandom_seed_array)

      ! write(*,*) 'Pushing the array to random_seed'

   end subroutine

   subroutine set_rand_epsilon( eps  ) !DOC 
!DOC set the minimal random value. The randoms will be generated in [eps;1-eps]  (to avoid zero results which sometimes cause errors)
!DOC Parameters:
       real(8) :: eps !DOC epsilon
 
       MRandom_rand_epsilon = eps

   end subroutine 

    function rand() !DOC
!DOC get the random value from (0;1)
       real(8) :: rand 
       integer :: ios

        if ( random_method == RND_STANDARD ) then      
                       
           MRandom_current_random_number = 0.d0
           do while( ( MRandom_current_random_number .lt. MRandom_rand_epsilon) .or. &
                     ( MRandom_current_random_number .gt. 1.d0 - MRandom_rand_epsilon ) & 
                   )                  
               call random_number(MRandom_current_random_number)
           end do
   !        write(*,*) 'MRandom_current_random_number',MRandom_current_random_number
           rand = MRandom_current_random_number 
        elseif ( random_method == RND_CONSTANT ) then

           rand = rnd_constant_const

        elseif ( random_method == RND_FILE ) then

           read(rnd_file_hfile,*,iostat=ios)
           if( ios < 0) then 
               call fseek(rnd_file_hfile, 0, 0, ios)
           end if 

        end if 

   end function      

   integer function random(N) !DOC
!DOC get the integer random value from 0 to N-1
        integer :: N
        real(8) :: rnd
       
        do while ( .TRUE. )  
   
           rnd = rand()
           if ( rnd < 1 ) exit

        end do
       
        random = aint(rnd * N )
   end function random


END MODULE MRandom
