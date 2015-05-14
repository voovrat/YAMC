Module MoleculeHandler !DOC 
!DOC !FILE the global storage for the molecule structures. Each molecule has it's identifier in MoleculeHandler
use Molecule, only : TMolecule

!
!DOC each molecule can be accessed using its internal "name" (handler)
!DOC these names are available for all parts of the program

! constants

     implicit none
     integer,parameter :: MAX_MOLECULE_HANDLER = 100   

!   
     logical :: MoleculeHandler_inited = .FALSE.
     logical, dimension(MAX_MOLECULE_HANDLER) :: MoleculeHandler_isFree ! private
     Type(TMolecule), dimension(MAX_MOLECULE_HANDLER),target :: MoleculeHandler_mol_list
     
CONTAINS

subroutine MoleculeHandler_init !DOC
!DOC Initialize the handler array
  
   MoleculeHandler_isFree(:) = .TRUE.
   MoleculeHandler_inited = .TRUE.

end subroutine

subroutine MoleculeHandler_getMolecule(h,mol_ptr) !DOC 
!DOC get the pointer to the molecule by its handler
   use error
!DOC Parameters: 
   integer,intent(in) :: h !DOC molecule handler
   Type(TMolecule),pointer :: mol_ptr !DOC output: pointer to the molecule

    if(.not.MoleculeHandler_inited) then
        call MoleculeHandler_init        
    end if
 

   nullify(mol_ptr)

   if(MoleculeHandler_isFree(h)) then
      write(error_message,*) 'MoleculeHandler_getMolecule: Handler ',h,'is not occupied'
      call error_throw(ERROR_PARAMETER)
      return
   end if

   mol_ptr => MoleculeHandler_mol_list(h)

end subroutine


subroutine MoleculeHandler_occupy(h,mol_ptr)  !DOC
!DOC try to occupy the handler, associates mol_ptr with the  handler h
!DOC Parameters:
    use error
    integer,intent(in) :: h !DOC molecule handler
    Type(TMolecule),pointer :: mol_ptr    !DOC output: pointer to the handler h (you should allocate it after occupation)

    if(.not.MoleculeHandler_inited) then
        call MoleculeHandler_init
    end if
 
    nullify(mol_ptr)

    if( ( h .le. 0 ).or.( h .gt. MAX_MOLECULE_HANDLER ) ) then

       write(error_message,*) 'MoleculeHandler_occupy: invalid handler:',h
       call error_throw(ERROR_PARAMETER)
       return
    end if

    if(.not.MoleculeHandler_isFree(h)) then
 
      write(error_message,*) 'MoleculeHandler_occupy: Handler',h,' is not free. Run MoleculeHandler_release first!'
      call error_throw(ERROR_PARAMETER)
      return 
    end if
    
    MoleculeHandler_isFree(h) = .FALSE.
    mol_ptr => MoleculeHandler_mol_list(h)

end subroutine   

subroutine MoleculeHandler_release(h) !DOC 
!DOC release the handler h
   use error
!DOC Parameters:
   integer,intent(in) :: h !DOC molecule handler

   if (.not.MoleculeHandler_inited) then 
      call MoleculeHandler_init
   end if
   

   if( ( h .le. 0 ).or.( h .gt. MAX_MOLECULE_HANDLER ) ) then

       write(error_message,*) 'MoleculeHandler_release: invalid handler:',h
       call error_throw(ERROR_PARAMETER)
       return
    end if

    if(MoleculeHandler_isFree(h)) then
 
      write(error_message,*) 'MoleculeHandler_released: Handler',h,' is free. Invalid (double?) release!'
      call error_throw(ERROR_PARAMETER)
      return 
    end if
    
    MoleculeHandler_isFree(h) = .TRUE.

end subroutine

subroutine MoleculeHandler_clear_all   !DOC
!DOC release and deallocates all molecules
   use Molecule, only : Molecule_deallocate
   implicit none
   integer :: i

   if(.not.MoleculeHandler_inited) then
      call MoleculeHandler_init
      return
   end if

   do i=1,MAX_MOLECULE_HANDLER
      
      if(.not.MoleculeHandler_isFree(i) ) then

          call Molecule_deallocate(MoleculeHandler_mol_list(i))
          MoleculeHandler_isFree(i) = .TRUE.
      end if

   end do


end subroutine

integer function MoleculeHandler_getFreeHandler() !DOC
!DOC get the free handler
    use error
    integer :: i

    if(.not.MoleculeHandler_inited) then
        call MoleculeHandler_init
    end if
 
    do i=1,MAX_MOLECULE_HANDLER
       
       if(MoleculeHandler_isFree(i)) then
           MoleculeHandler_getFreeHandler = i
           return
       end if
    end do

    error_message = 'MoleculeHandler_getFreeHandler: No free handlers available'
    call error_throw(ERROR_LIMITS)
    MoleculeHandler_getFreeHandler = -ERROR_LIMITS

end function



END MODULE MoleculeHandler
