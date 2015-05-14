Module io !DOC
!DOC !FILE input/output interface 

   implicit none

   integer, parameter :: IO_MAX_HANDLER = 100
   integer, parameter :: IO_RESERVED_HANDLER = 10

! private
    logical, dimension( IO_MAX_HANDLER ) ::  io_free_handlers
    logical :: io_initialized = .FALSE.

CONTAINS


subroutine io_init !DOC
!DOC initialize the io module
   implicit none
   integer :: i

   do i=1,IO_RESERVED_HANDLER
      io_free_handlers(i) = .FALSE.
   end do 

   do i=IO_RESERVED_HANDLER + 1,IO_MAX_HANDLER
      io_free_handlers(i) = .TRUE.
   end do
    
   io_initialized = .TRUE.

end subroutine

subroutine close_all !DOC
!DOC close all the opened files
   use error, only : error_code
   implicit none
   integer :: h
   
   do h=IO_RESERVED_HANDLER+1,IO_MAX_HANDLER

     if(.not.io_free_handlers(h)) then
          call io_close(h)
          if(error_code/=0) return
      end if

   end do

end subroutine


integer function io_open(filename,mode) !DOC
!DOC opens the file and  return file handler
!DOC Parameters:
    use SystemSettings, only : SYSTEM_STRING_LENGTH
    use error

    implicit none
    character(SYSTEM_STRING_LENGTH) :: filename !DOC name of the fule
    character :: mode !DOC mode: r - read, w - write
    integer :: h
    integer :: stat

    if ( .not.io_initialized ) then
       call io_init
    end if

    do h=IO_RESERVED_HANDLER+1,IO_MAX_HANDLER
       if ( io_free_handlers(h) ) then
        
          if ( mode == 'r' ) then

             open(h,file=filename,iostat=stat)
             read(h,*,iostat=stat)
             close(h)

             if(stat==0) open(h,file=filename,iostat=stat)

          elseif (mode == 'w' ) then

             open(h,file=filename,form='formatted',iostat=stat)

          elseif (mode == 'b' ) then

             open(h,file=filename,form='unformatted',access='stream',iostat=stat)

          else 
             write(error_message,*)  ' io_open : unknown mode :  ',mode
             call error_throw(ERROR_PARAMETER) 
             io_open = -ERROR_PARAMETER
             return
          end if

          if( stat /= 0) then
             
             write(error_message,*)   'io_open: File cannot be opened:',filename(1:80)
             call error_throw(ERROR_IO)
             io_open = -ERROR_IO
             return
          end if
  
          io_free_handlers(h) = .FALSE.
          io_open = h;

          return
       end if ! free handler

    end do

    
    write(error_message,*) 'io_open: Free handler not found!!!'
    call error_throw(ERROR_LIMITS)
    io_open = -ERROR_LIMITS;
    return  

end function io_open
 
subroutine io_close(hfile) !DOC
!DOC close the file opened with io_open
!DOC Parameters: 
   use error
   implicit none
   integer :: hfile !DOC file handler

    if ( .not.io_initialized ) then
       call io_init
    end if

    if(io_free_handlers(hfile)) then
       write(error_message,*) 'io_close: File is not open',hfile
       call error_throw(ERROR_PARAMETER)
       return
    end if

    if( (hfile.le.IO_RESERVED_HANDLER) .or. (hfile.gt.IO_MAX_HANDLER) ) then
       write(error_message,*) 'io_close: invalid file handler', hfile
       call error_throw(ERROR_PARAMETER)
       return
    end if
 
    close(hfile)
    io_free_handlers(hfile) = .TRUE.

end subroutine

function io_count_file_lines(filename) !DOC
!DOC Count number of lines in the file
   use SystemSettings, only : SYSTEM_STRING_LENGTH
   use error
   integer :: io_count_file_lines

!DOC Parameters:
   character(SYSTEM_STRING_LENGTH) :: filename  !name fo the file
!DOC Return value:
!DOC number of lines in the file
   integer :: hfile  
   integer :: ios
   integer :: n

   hfile = io_open(filename,'r');   

   if(error_code /=0) return

   n = 0
   do while(.true.)

      read(hfile,*,iostat=ios)

      if(ios>0) then
         write(error_message,*) 'io_count_file_lines: read error'
         call error_throw(ERROR_IO)
         io_count_file_lines = -ERROR_IO
         call io_close(hfile)
         return
      elseif(ios<0) then ! end of file reached
          exit
      end if 

      n=n+1
  
   end do

   io_count_file_lines = n 
   call io_close(hfile)
 
end function

subroutine write_little_endian( fid, val, n) !DOC
!DOC Write the integer value to the file
!DOC Parameter:
  implicit none
  integer,intent(in) :: fid !DOC file handler
  integer,intent(in) :: val !DOC the integer value
  integer,intent(in) :: n !DOC numer of bytes to be written
  integer :: vval,nn

  integer :: last_byte
  integer(1) :: ch

  vval = val 
  nn = n

  do while(nn>0)

    last_byte = mod( vval, 256 )
    if ( last_byte .ge. 128 ) then

       ch = last_byte - 256

    else
       ch = last_byte
    end if

    write(fid) ch  
    nn = nn - 1
    vval = vval / 256
  end do


end subroutine

subroutine read_little_endian(fid,n,val) !DOC
!DOC read integer value form file
!DOC Parameters:
   use error 

   integer,intent(in) :: fid !DOC file descriptor
   integer,intent(in) :: n !DOC number of bytes to read
   integer,intent(out) :: val !DOC output: value
   integer(1) :: ch

   integer :: nn,last_byte
   integer :: power  
   integer :: ios

   nn = n
   val = 0
   power = 1
   do while (nn>0)

      read(fid,iostat=ios) ch

      if ( ios /= 0) then
          error_message = 'read_little_endian : IO Error'
          call error_throw(ERROR_IO)
          return
      end if

      if( ch<0) then
         last_byte = ch + 256
      else
         last_byte = ch
      end if

      val = val +  last_byte * power
      power = power * 256
      nn = nn - 1
   end do

end subroutine

subroutine write_real_array(hFile,arr,n) !DOC
!DOC write the array to file or on the screen
!DOC Parameters:
   integer,intent(in) :: hFile !DOC file handler. 0 means screen
   real(8),dimension(:),intent(in) :: arr !DOC array
   integer,intent(in) :: n !DOC number of elements in the array
   
   integer :: i

   do i=1,n

      if ( hFile == 0 ) then
         write(*,*) arr(i)
      else
         write(hFile,*) arr(i)
      end if

   end do

end subroutine

subroutine write_xyz_array(hFile,arr_x,arr_y,arr_z,n) !DOC
!DOC write 3 arrays representing x,y,z coordinates
!DOC Parameters:
    integer,intent(in) :: hFile !DOC file handler (0 means screen)
   real(8),dimension(:),intent(in) :: arr_x,arr_y,arr_z !DOC arrays
   integer,intent(in) :: n !DOC number of elements
   
   integer :: i

   do i=1,n

      if ( hFile == 0 ) then
         write(*,*) arr_x(i),arr_y(i),arr_z(i)
      else
         write(hFile,*) arr_x(i),arr_y(i),arr_z(i)
      end if

   end do

end subroutine


subroutine write_real_matrix(hfile, matrix, m, n )!DOC
!DOC write the matrix m*n to the file
    use SystemSettings, only : TRealArrayPointer
!DOC Parameters:
    integer,intent(in) :: hfile
    Type(TRealArrayPointer),dimension(:),intent(in) :: matrix  !DOC the matrix: array of arrays
    integer,intent(in) :: m,n   !DOC size of the matrix
 
    integer :: i,j

    do i=1,m

        do j=1,n

            write(hfile,'(F20.10,$)')  matrix(j) % ptr(i)

        end do

        write(hfile,*)   

    end do

end subroutine

subroutine write_integer_array(hFile,arr,n) !DOC
!DOC write real array 
!DOC Parameters:
   integer,intent(in) :: hFile !DOC file handler
   integer,dimension(:),intent(in) :: arr !DOC array
   integer,intent(in) :: n !DOC number of elements
   
   integer :: i

   do i=1,n

      if ( hFile == 0 ) then
         write(*,*) arr(i)
      else
         write(hFile,*) arr(i)
      end if

   end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE io

