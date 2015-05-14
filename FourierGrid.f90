Module FourierGrid !DOC
!
!DOC !FILE  Grid in fourier space
!DOC  includes all k such that |kk|<kmax
!

implicit none


Type :: TFourierGrid !DOC
!DOC Fields:
   integer :: kmax  !DOC limit for |k| < kmax
   integer :: nk     !DOC number of grid points  
   integer :: nalloc !DOC size of allocated arrays

   integer,dimension(:), allocatable :: ip  !DOC next raw indicator
         !DOC ip(kk) = 0  -- no charnges
         !DOC ip(kk) = 1  -- next ky
         !DOC ip(kk) = 2  -- next kz
         !DOC ip(kk) = 3  -- end of array

   integer,dimension(:), allocatable :: iq  !DOC zero indicator (if kx,ky,or kz == 0 )
    !DOC iq(kk) = 1  -->  ky=0, kz=0     
    !DOC iq(kk) = 2  -->  ky>0, kz=0
    !DOC iq(kk) = 3  -->  ky=0, kz>0
    !DOC iq(kk) = 4  -->  ky>0, kz>0

   integer,dimension(:),allocatable :: kx,ky,kz   !DOC indeed these are mx,my,mz in my notation, where k = 2pi m/L
   real(8),dimension(:),allocatable :: xk,xk2     !DOC |k|, k^2

End Type TFourierGrid


contains


subroutine FourierGrid_calcNalloc(kmax,nalloc) !DOC 
!DOC calculate number of grid points for a given kmax
!DOC Parameters:
   integer,intent(in) :: kmax ! maximal value for |k|
   integer,intent(out) :: nalloc  ! output: number of k values

   integer :: kx,ky,kz
   real(8) :: xkmax2,xkx2,xky2,xkz2  ! kmax^2, kx^2, ky^2, kz^2
   real(8) :: xk2                    ! |k|^2
   integer :: kymax,kzmax

   integer :: kk

   xkmax2 = kmax**2

   kk=0 
   do kx=0,kmax
      
      xkx2=kx**2
      kymax=SQRT(xkmax2-xkx2)
      
      do ky=0,kymax
      
         xky2=ky**2
         kzmax=SQRT(xkmax2-xkx2-xky2)
       
         do kz=0,kzmax

             xk2=xkx2+xky2+kz**2
             if(xk2<0.1) cycle            ! k different de 0

             kk=kk+1

         end do  ! kz
      end do ! ky
   end do ! kx

   nalloc = kk


end subroutine


subroutine FourierGrid_alloc(this,nalloc) !DOC
!DOC allocate the KSpace grid
!DOC Parameters:
   Type(TFourierGrid) :: this !DOC FourierGrid structure
   integer,intent(in) :: nalloc !DOC number of k values

   allocate( this % kx(nalloc) )
   allocate( this % ky(nalloc) )
   allocate( this % kz(nalloc) )
   
   allocate( this % ip(nalloc) )
   allocate( this % iq(nalloc) )

   allocate( this % xk(nalloc) )
   allocate( this % xk2(nalloc) )

   this % nalloc = nalloc
   this % nk = 0

end subroutine

subroutine FourierGrid_dealloc(this) !DOC
!DOC deallocate the KSpace grid
!DOC Parameters:
   Type(TFourierGrid) :: this !DOC FourierGrid structure

   deallocate( this % kx )
   deallocate( this % ky )
   deallocate( this % kz )
   
   deallocate( this % ip )
   deallocate( this % iq )
  
   this % nalloc = 0
   this % nk = 0   

end subroutine

subroutine FourierGrid_init(this,kmax) !DOC
!DOC initialize the KSpace grid (fill the kx,ky,kz,k^2,ip,iq arrays)
   use error
!DOC Parameters:
   Type(TFourierGrid) :: this !DOC FourierGrid structure
   integer,intent(in) :: kmax !DOC maximum value for |k|
   
   integer :: kx,ky,kz
   real(8) :: xkmax2,xkx2,xky2,xkz2  ! kmax^2, kx^2, ky^2, kz^2  
   real(8) :: xk2                    ! |k|^2 
                                     ! note: prefix 'x' indicates that the values are real
   integer :: kymax,kzmax

   integer :: kk


   this % kmax = kmax
   xkmax2 = kmax**2

   kk=0 
   do kx=0,kmax
      
      xkx2=kx**2
      kymax=SQRT(xkmax2-xkx2)
      
      do ky=0,kymax
      
         xky2=ky**2
         kzmax=SQRT(xkmax2-xkx2-xky2)
       
         do kz=0,kzmax

             xk2=xkx2+xky2+kz**2
             if(xk2<0.1) cycle            ! k different de 0

             kk=kk+1
             if( kk > this % nalloc ) then
                 write(error_message,*) 'FourierGrid_init: k exceeds allocated limit ',this % nalloc
                 call error_throw(ERROR_LIMITS)
                 return
             end if

             this % kx(kk) = kx
             this % ky(kk) = ky
             this % kz(kk) = kz

             this % xk2(kk) = xk2
             this % xk(kk) = sqrt(xk2)

             if(ky.eq.0) then                           ! iq sert a savoir si l'un des ki est nul
                if(kz.eq.0) this % iq(kk)=1
                if(kz.gt.0) this % iq(kk)=3
             else
                if(kz.eq.0) this % iq(kk)=2
                if(kz.gt.0) this % iq(kk)=4
             end if

             this % ip(kk)=0

         end do  ! kz

         this % ip(kk) = 1
      end do ! ky

      this % ip(kk) = 2 
   end do ! kx
  
   this % ip(kk) = 3

   this % nk = kk

end subroutine


End Module FourierGrid
