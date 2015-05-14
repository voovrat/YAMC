Module String !DOC 
!DOC !FILE Functions to work with the strings
use SystemSettings, only : SYSTEM_STRING_LENGTH

contains

pure function str_get_next_pos(str,offset,ch) !DOC 
!DOC get the position of the first occurance of the symbol ch in the string str starting from the position offset
!DOC if not found then -1 is returned
   integer :: str_get_next_pos
!DOC Parameters:
   character(SYSTEM_STRING_LENGTH),intent(in) :: str !DOC string to search
   integer,intent(in) :: offset !DOC start position
   character,intent(in) :: ch !DOC character

   integer :: i

   i = offset

   do while( (i < SYSTEM_STRING_LENGTH ) .and. ( str(i:i) /= ch ) )

      i = i + 1    

   end do

   if ( i == SYSTEM_STRING_LENGTH ) i = -1

   str_get_next_pos = i


end function


pure function str_isempty(str) !DOC
!DOC check if the string is empty (contains only spaces)
   logical :: str_isempty

   character(SYSTEM_STRING_LENGTH),intent(in) :: str

   integer :: i
 
   do i=1,SYSTEM_STRING_LENGTH

      if( str(i:i) /= ' ') then
         str_isempty = .FALSE.
         return
      end if

   end do

   str_isempty = .TRUE.
   

end function

subroutine str_subs(str, ch_old, ch_new) !DOC
!DOC substitute in the string str the symbol ch_old with the symbol ch_new
   character(SYSTEM_STRING_LENGTH) :: str
   character :: ch_old, ch_new
  
   integer :: i

   do i=1,SYSTEM_STRING_LENGTH

      if ( str(i:i) == ch_old ) then
         str(i:i) = ch_new
      end if

   end do
   
end subroutine


End Module String
