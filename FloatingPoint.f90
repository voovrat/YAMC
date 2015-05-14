Module FloatingPoint !DOC
!DOC !FILE  Module to work with the floating point numbers (construct from mantisse and exponent, store to and restore from file)
implicit none
!DOC Why in general to use this module, not the standard one functions?
!DOC Because, the standard functions can be system- or hardware- dependent.
!DOC This module is completely system independent
!
!DOC subroutines to extract mantisse and exponent
!DOC in representation   a = sign*m*2^n, where  1 <= m < 2 
!DOC
!DOC  Format to save:
!DOC   mm = (|m|-1)*(2^(mw-1)-1)  +  (m<0)*2^(mw-1)
!DOC   nn = |n| + (n<0)*2^(nw-1)
!DOC
!DOC  where mw stands for "mantisse width", number of bits for mantisse (normally 24)
!DOC        nw - "exponent width", number of bits for exponent (normally 8)
!DOC
!DOC  The values mm,nn are stored in the file in the little-endian format.
!DOC  I did not optimize for non-whole byte widths, so if e.g. mv=18,nw=6 it will use 4 bytes anyway

contains

function myint(x) !DOC
!DOC myint is the maximal integer number smaller than x
!DOC Parameters:
   integer :: myint
   real(8),intent(in) :: x !DOC real number x
!DOC Return value:
!DOC   maximal integer smaller than x

   if( x>0) then 
      myint = aint(x)
   else
      myint = aint(x - 1.d0 + 1.d-14 )
   endif

end function


subroutine mantisse_exponent_real(a,m,n) !DOC
!DOC Get mantisse and exponent for the real number
   use constants, only : LN2
!DOC Parameters:
   real(8),intent(in) :: a  !DOC real number
   real(8) :: m !DOC output: mantisse(real)
   integer :: n !DOC output: exponent(integer)

   if( dabs(a) < 1e-14) then
      m = sign(1.d0,a)
      n = -1000
      return
   end if

   if( dabs(a) > 1.e14 ) then
      m = sign(1.d0,a)
      n = 1000
      return
   end if       
   
   ! 2^n < |a| < 2^(n+1) ==>  n < |a| < n+1
   n = myint( dlog( dabs(a)) / LN2 )
   m = a / 2.d0**n

end subroutine

subroutine mantisse_exponent_bits(m,n,mantisse_width,exp_width,mm,nn) !DOC
!DOC convert mantisse and exponent to the integer numbers
!  
!DOC   a = sign * m * 2^n
!DOC   m,n taken from mantisse_exponent_real  
!
!DOC Parameters: 
  real(8),intent(in) :: m  !DOC mantisse (real)
  integer,intent(in) :: n  !DOC exponent (integer)
  integer,intent(in) :: mantisse_width,exp_width !DOC number of bits in mantisse and exponent
  integer :: mm,nn !output: integer numbers which represent the mantisse and exponent

 
  integer :: m_scale ! 2^(mantisse_width-1)
  integer :: n_scale ! 2^(exp_width-1)

  integer :: sgn_m,sgn_n

  m_scale = 2**(mantisse_width-1) 
  
  n_scale = 2**(exp_width-1) 

  if( m > 0) then
    sgn_m = 0
  else
    sgn_m = 1
  end if 

  if( n > 0) then
     sgn_n = 0
  else
     sgn_n = 1
  end if

  if( abs(n) .ge. n_scale ) then ! under/overflow
    mm = sgn_m * m_scale  
    nn = sgn_n * n_scale  + n_scale - 1
    return
  end if 

  mm =   sgn_m * m_scale + aint( ( dabs(m) - 1.d0 ) * (m_scale - 1) )
  nn =   sgn_n * n_scale + abs(n)

end subroutine

function construct_float(mm,nn,mantisse_width,exp_width) !DOC
!DOC construct the real number from integer values representing mantisse and exponent
   real(8) :: construct_float
 !DOC Parameters:
   integer,intent(in) :: mm,nn  !DOC integer values representing mantisse and exponent
   integer,intent(in) :: mantisse_width,exp_width !DOC number of bits in mantisse and exponent

   integer :: mm1,nn1

   integer :: m_scale,n_scale
   integer :: sign_m,sign_n
   
   real(8) :: x
   real(8) :: m 
 
   m_scale = 2**(mantisse_width-1) 
   n_scale = 2**(exp_width-1) 

   if( mm .ge. m_scale ) then
      mm1 = mm - m_scale
      sign_m = -1
   else
      mm1 = mm
      sign_m = 1
   end if

   if( nn .ge. n_scale ) then
       nn1 = nn - n_scale
       sign_n = -1
   else
       nn1 = nn
       sign_n = 1
   end if

!   write(*,*) 'nn',nn1,'mm',mm1,'n_scale',n_scale,'m_scale',m_scale,'sign_m',sign_m,'sign_n',sign_n
 
   if( nn1 == n_scale-1) then

       if( sign_n == -1) then
           construct_float = 0.d0
       else   
           construct_float = 1.d14
       end if

       return
   end if


  m = mm1 / ( m_scale - 1.d0) + 1.d0
!  write(*,*) 'mm1',mm1,'m_scale-1',m_scale-1.d0,'m',m
  x = sign_m * m  * 2.d0**(sign_n * nn1 )
  construct_float = x 

end function

subroutine write_float(hfile,x,mantisse_width_in,exp_width_in) !DOC
!DOC write real number to file
!DOC Parameters:
   use io, only : write_little_endian
   integer,intent(in) :: hfile !DOC file handler
   real(8),intent(in) :: x !DOC real number

   integer,intent(in),optional :: mantisse_width_in,exp_width_in !DOC number of bits in mantisse and exponent

   real(8) :: m
   integer :: n
   integer :: mm,nn
   
   integer :: mantisse_width,exp_width, mantisse_bytes, exp_bytes

   if(present(mantisse_width_in)) then
      mantisse_width = mantisse_width_in
   else
      mantisse_width = 24
   end if

   if(present(exp_width_in)) then
      exp_width = exp_width_in
   else
      exp_width = 8
   end if

   mantisse_bytes = aint( (mantisse_width+0.999d0) / 8 )
   exp_bytes = aint( (exp_width + 0.999d0) / 8 )
 
!   write(*,*) 'mb',mantisse_bytes,'nb',exp_bytes

 
   call mantisse_exponent_real(x,m,n)
   call mantisse_exponent_bits(m,n,mantisse_width,exp_width,mm,nn)

!   write(*,*) 'x',x,'m',m,'n',n,'mm',mm,'nn',nn

   call write_little_endian(hfile,mm,mantisse_bytes)
   call write_little_endian(hfile,nn,exp_bytes)

end subroutine

function read_float(hfile,mantisse_width_in,exp_width_in) !DOC
!DOC read real number from file
   use io, only : read_little_endian
!DOC Parameters:
   real(8) :: read_float 

   integer,intent(in) :: hfile !DOC file handler
   integer,intent(in),optional :: mantisse_width_in,exp_width_in !DOC number of bits in mantisse and exponent
!DOC Return value:
!DOC  real number read from file
   real(8) :: m
   integer :: n
   integer :: mm,nn
   
   integer :: mantisse_width,exp_width, mantisse_bytes, exp_bytes

   if(present(mantisse_width_in)) then
      mantisse_width = mantisse_width_in
   else
      mantisse_width = 24
   end if

   if(present(exp_width_in)) then
      exp_width = exp_width_in
   else
      exp_width = 8
   end if

   mantisse_bytes = aint( (mantisse_width+0.999d0) / 8 )
   exp_bytes = aint( (exp_width + 0.999d0) / 8 )

   call read_little_endian(hfile,mantisse_bytes,mm)
   call read_little_endian(hfile,exp_bytes,nn)

!   write(*,*) 'mm',mm,'nn',nn

   read_float = construct_float(mm,nn,mantisse_width,exp_width)

end function



End Module FloatingPoint
