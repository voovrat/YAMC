Module matrix3x3 !DOC
!
!DOC !FILE Operations with matrices 3x3
!DOC ( they are particullary intersting, because can be used for the coordinate-transformations)
!DOC  Also, for these matrices one knows the explicit relations for determinant and inverse matrices
! 

implicit none

contains


subroutine matrix3x3_mul(A,B,C) !DOC
!DOC matrix multiplication C=A*B
!DOC Parameters:
   real(8),dimension(3,3),intent(in) :: A,B !DOC multiplicands
   real(8),dimension(3,3) :: C !DOC result

   integer :: i,j,k
   real(8) :: S

   do i=1,3
      do j=1,3

      S = 0
      do k=1,3
         S = S + A(i,k) * B(k,j)
      end do
      C(i,j) = S

      end do
   end do

end subroutine 

!  
! 
pure function matrix3x3_det(A) !DOC
!DOC calculate matrix determinant
!DOC Parameters:
    real(8) :: matrix3x3_det
    real(8),dimension(3,3),intent(in) :: A !DOC the matrix
!DOC Return value:
!DOC determinant det(A)

    matrix3x3_det =  A(1,1)*A(2,2)*A(3,3)  + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) & 
                   - A(1,1)*A(2,3)*A(3,2)  - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)  

end function


 
subroutine matrix3x3_adj(A,B) !DOC
!DOC auxilarly function to calculate the inverse matrix
!DOC  adj(A) = det(A) * A^-1
!  adj(A)_ij = (-1)^{i+j}*det(Minor_ji)
!  can be obtained using the explicit formula and symbolic calculations
!  see alse http://de.wikipedia.org/wiki/Inverse_Matrix for explicit forumla
!DOC Parameters:
  real(8),dimension(3,3),intent(in) :: A !DOC input matrix
  real(8),dimension(3,3) :: B !DOC output matrix

  B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
  B(2,1) =  - A(2,1)*A(3,3) + A(2,3)*A(3,1)
  B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
  B(1,2) =  - A(1,2)*A(3,3) + A(1,3)*A(3,2)
  B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
  B(3,2) =  - A(1,1)*A(3,2) + A(1,2)*A(3,1)
  B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
  B(2,3) =  - A(1,1)*A(2,3) + A(1,3)*A(2,1)
  B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

end subroutine


subroutine matrix3x3_inv(A,B) !DOC
!DOC calculate matrix inverse
! A^-1 = inv(A) = adj(A)/det(A)   
!DOC Parameters: 
  real(8),dimension(3,3),intent(in) :: A !DOC input matrix
   real(8),dimension(3,3) :: B !DOC B=A^-1

   real(8) :: D
   integer :: i,j

   D = matrix3x3_det(A)
   call matrix3x3_adj(A,B)

   do i=1,3
      do j=1,3
         B(i,j) = B(i,j) / D
      end do
   end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module matrix3x3
