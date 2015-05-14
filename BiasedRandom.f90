Module BiasedRandom !DOC 
!DOC !FILE BiasedRandom - generate the random numbers 1..N , where the probability of i is proportional to prob(i)

Type :: TBiasedRandom !DOC
!DOC Fields:
   integer :: N  !DOC number of intervals 
   real(8),dimension(:),allocatable :: freq_prob  !DOC probabilities to move the particle:  freq_prob = freq_num / SUM(freq_num)
   real(8),dimension(:),allocatable :: freq_beg, freq_end !DOC begins and ends of the intervals for each molecule
 
End Type TBiasedRandom

contains

subroutine BiasedRandom_alloc(this,N) !DOC
!DOC  Allocate the BiasedRandom structure 
!DOC Parameters:
   Type(TBiasedRandom) :: this !DOC BiasedRandom structure
   integer,intent(in) :: N  !DOC number of the elements to be allocated

!   allocate(this % freq_num(N) )
   allocate(this % freq_prob(N) )
   allocate(this % freq_beg(N) )
   allocate(this % freq_end(N) )

   this % N = N

end subroutine

subroutine BiasedRandom_dealloc(this) !DOC
!DOC Deallocate the BiasedRandom structure
!DOC Parameters:
   Type(TBiasedRandom) :: this !DOC BiasedRandom structure

!   deallocate( this % freq_num )
   deallocate( this % freq_prob )
   deallocate( this % freq_beg )
   deallocate( this % freq_end ) 

   this % N = 0

end subroutine


subroutine BiasedRandom_init(this, freq_num ) !DOC
!DOC initializes the frequency array in the BiasedRandom structure
!DOC Parameters:
  Type(TBiasedRandom) :: this !DOC BiasedRandom structure 
  real(8),dimension(:),intent(in) :: freq_num !DOC frequences distribution (maybe not normalized to 1 )
 
  integer :: i
  real(8) :: S

  S = SUM(freq_num(:))
  
  if ( S == 0 ) then
     this % freq_prob(:) = 1.d0 / this % N 
  else
     this % freq_prob(:) = freq_num(:) / S
  end if

  this % freq_beg(1) = 0d0
  this % freq_end(1) = this % freq_prob(1)
  
  do i=2,this % N

     this % freq_beg(i) = this % freq_end(i-1)
     this % freq_end(i) = this % freq_beg(i) + this % freq_prob(i)

  end do

  if ( S > 0) then 
     do i=1,this % N

       if ( freq_num(i) == 0 ) this % freq_end(i) = -1d0 ! to avoid some problems with zero size intervals. 0 means 0

     end do
  end if 

end subroutine

function  BiasedRandom_choose(this) !DOC
!DOC chooses the random number using the given frequences defined in BiasedRandom structure
!DOC Parameters:
   use MRandom, only : rand
   integer :: BiasedRandom_choose 

   Type(TBiasedRandom) :: this  !DOC BiasedRandom structure
!DOC Return value:
!DOC       random integer number from 1..N with probability proportional to the given frequences array  
   real(8) :: rnd
   integer :: i,interval

   rnd = rand()

   interval = 1

   do i=1,this % N
   
      if ( ( rnd > this % freq_beg(i) ) .AND. ( rnd < this % freq_end(i) ) ) then

          interval = i
          exit
      end if

   end do

   BiasedRandom_choose = interval

end function



End Module BiasedRandom
