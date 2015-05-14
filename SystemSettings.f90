Module SystemSettings !DOC
!DOC !FILE System-dependent constants

   integer, parameter :: SYSTEM_STRING_LENGTH = 256 !DOC Standart string length

   integer, parameter :: STDERR =0 !DOC standard error output file
   integer, parameter ::  STDIN =5 !DOC standard input file
   integer, parameter ::  STDOUT =6 !DOC standard output file
   integer, parameter :: SEEK_SET = 0 !DOC relative to the beginging of the file
   integer, parameter ::  SEEK_CURR = 1 !DOC relative to the current position
   integer, parameter ::  SEEK_END = 2    !DOC relative to the end of the file

   Type :: TRealArrayPointer !DOC
   !DOC the structure is used to have pointers to the array pointers
       real(8),dimension(:),pointer :: ptr !DOC pointer
   End Type

END MODULE SystemSettings
