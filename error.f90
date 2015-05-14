Module Error !DOC
use SystemSettings, only : SYSTEM_STRING_LENGTH
!DOC !FILE The module contains the functions and data structures to deal with run-time errors
!DOC how do errors work:
!DOC if there is an error in the program, the function sets error_message to some value, and run throw with some error code
!DOC 
!DOC throw checks, if the error is in catch list.
!DOC if it is, that means that the programmer thought about the possibility of that error, thus functions just returns
!DOC otherwise, the executions stops

! example:
!
! function one_over_x(x)
!   use Error  
!  real :: x
!
!   if(x==0) then
!      error_message = '1/x: argument is zero!'
!      call error_throw(ERROR_PARAMETER)
!      return
!   end if 
!
!   one_over_x = 1/x
!
! end function
!
! program test
!   
!   real :: x=0,y
!
!   call error_set_catch(ERROR_PARAMETER)
!    
!   y = 1/x
!
!   if ( error_code == ERROR_PARAMETER ) then
!        write(*,*) 'y is infinity. this is not a problem for us'
!   end if
!
!   call error_clear_catch(ERROR_PARAMETER)
!
! end program
!

    implicit none
! public:
!DOC error codes:
    integer, parameter :: ERROR_IO = 1   !DOC input/output error
    integer, parameter :: ERROR_LIMITS = 2   !DOC error with limits (sizes of arrays)
    integer, parameter :: ERROR_PARAMETER = 3 !DOC incorrect parameter given
    integer, parameter :: ERROR_INITIALIZATION = 4 !DOC the function was called before the initialization 
    integer, parameter :: ERROR_WRONG_FUNCTION = 5 !DOC the function was called for incorrect data 

! private:

! catch list 
    integer, parameter :: MAX_CATCH_ERROR = 100
    integer, dimension(MAX_CATCH_ERROR) :: error_catch_list
    integer :: error_last_catch = 0    

!public: 
    character(SYSTEM_STRING_LENGTH) :: error_message
    integer :: error_code = 0
  
CONTAINS
 
subroutine error_set_catch(err_code)  !DOC
!DOC the function which is used to set the error catch for the specific error code, which means that this code will not cause the stop of the program
   implicit none
!DOC Parameters:
   integer,intent(in) :: err_code !DOC error code

   if(error_last_catch.ge.MAX_CATCH_ERROR) then
       error_message='catch_error: error stack is full!'
       call error_throw(ERROR_LIMITS)
       return
   end if

   error_last_catch = error_last_catch + 1

   error_catch_list(error_last_catch) = err_code

end subroutine

subroutine error_clear_catch(err_code) !DOC
   implicit none
!DOC the function clears the catch for the error, that means that the errors with err_code will cause the stop of the program
!DOC Parameters:
   integer,intent(in) :: err_code   !DOC error code

   if(error_catch_list(error_last_catch) /= err_code) then
      error_message = 'error_clear_catch: last err_code does not match'
      call error_throw(ERROR_PARAMETER)
      return 
   end if 

   if(error_last_catch.le.0) then
      error_message = 'error_clear_catch: error_catch_list is empty'
      call error_throw(ERROR_PARAMETER)
      return 
   end if 

   if(error_code == err_code) then  ! ok, error_code is global variable err_code - local
                                    ! if the error was caught, and we reached this statement -- this means the error is already elaborated, 
                                    ! the program is not in error status any more
      error_code = 0
   end if

   error_last_catch = error_last_catch - 1

end subroutine 

subroutine error_throw(err_code) !DOC
   implicit none
!DOC the subroutine is called then some errorneus situation occurs.
!DOC Parameters:
   integer,intent(in) :: err_code  !DOC   err_code describes the situation (see error codes above)
   integer :: i

   error_code = err_code

   do i=1,error_last_catch

      if(error_catch_list(i) == err_code ) return

   end do

   write(*,*) 'ERROR: ',error_message
   stop

end subroutine


END MODULE Error 
