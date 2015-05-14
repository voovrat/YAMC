Module EwaldSumExternal !DOC
!DOC !FILE The module contains the data structures and the functions to deal with the external permutivity component of the ewald sums, i.e.
!DOC !FILE   2 pi / (2 external_permutivity + 1) (mu_x^2 + mu_y^2 + mu_z^2)
!DOC !FILE where mu = SUM_j q_j r_j 

implicit none

Type :: TEwaldSumExternal !DOC
!DOC Fields:
   real(8) :: mu_x,mu_y,mu_z  !DOC mu = SUM_j q_j r_j
   real(8),dimension(:),pointer :: xx,yy,zz,charge !DOC pointers to the coordinates and charges arrays (usually stored in AtomicData)
   integer :: natom !DOC number of atoms

End Type TEwaldSumExternal

contains

subroutine EwaldSumExternal_init(this,xx,yy,zz,charge,natom) !DOC
!DOC initialize the EwaldSumExternal structure
!DOC Parameters:
   Type(TEwaldSumExternal) :: this !DOC EwaldSumExternal structure
   real(8),dimension(:),intent(in),target :: xx,yy,zz,charge !DOC coordinates and charges
   integer,intent(in) :: natom !DOC number of atoms

   this % xx => xx
   this % yy => yy
   this % zz => zz
   this % charge => charge

   this % natom = natom

end subroutine

subroutine EwaldSumExternal_calc_mu(this) !DOC
!DOC calculate the moment mu=(mu_x,mu_y,mu_z)
   Type(TEwaldSumExternal) :: this !DOC EwaldSumExternal structure

   integer :: i

   this % mu_x = 0
   this % mu_y = 0
   this % mu_z = 0

   do i=1,this % natom

      this % mu_x = this % mu_x + this % charge(i) * this % xx(i)
      this % mu_y = this % mu_y + this % charge(i) * this % yy(i)
      this % mu_z = this % mu_z + this % charge(i) * this % zz(i)

   end do

end subroutine

pure function EwaldSumExternal_calc_energy(this) !DOC
!DOC calculate the energy E = kext * (mu_x^2 + mu_y^2 + mu_z^2)
!DOC where kext =  2 pi / (2 external_permutivity + 1)
   use parameters, only : kext
   real(8) :: EwaldSumExternal_calc_energy
!DOC Parameters:
   Type(TEwaldSumExternal),intent(in) :: this !DOC EwaldSumExternal structure
!DOC Return value: 
!DOC  external permutivity component of the ewald sum
   
    EwaldSumExternal_calc_energy = kext * (this % mu_x**2 + this % mu_y**2 + this %  mu_z**2 ) 

end function

pure function EwaldSumExternal_calc_dU(this, dmu_x,dmu_y,dmu_z) !DOC
  use parameters, only : kext
  real(8) :: EwaldSumExternal_calc_dU
!DOC Calculate the energy difference dU = kext*[(mu+dmu)^2 - mu^2]
!DOC Parameters:
  Type(TEwaldSumExternal),intent(in) :: this !DOC EwaldSumExternal structure
  real(8),intent(in) :: dmu_x, dmu_y, dmu_z !DOC x,y,z components of the dmu vector

  EwaldSumExternal_calc_dU = kext *(  &
                             ( 2 * this % mu_x  + dmu_x) * dmu_x  &
                           + ( 2 * this % mu_y  + dmu_y) * dmu_y  &
                           + ( 2 * this % mu_z  + dmu_z) * dmu_z  )  

end function


subroutine external_sum_calc_forces( mu_x, mu_y, mu_z, charge, natom, fx,fy,fz ) !DOC
   use parameters, only : minus_two_kext
!DOC Calculate the forces F = dU/dr
!DOC Parameters:
   real(8),intent(in) :: mu_x,mu_y,mu_z !DOC components of mu vector
   real(8),dimension(:),intent(in) ::  charge !DOC charges (array)
   integer,intent(in) :: natom !DOC number of atoms (length of charge array)

   real(8),dimension(:) :: fx,fy,fz !DOC Output arguments: forces (arrays)

   ! minus_two_kext = -2 * (2 pi / ( 2 epsilon_ext + 1 ))


   fx(1:natom) = minus_two_kext * charge(1:natom) * mu_x
   fy(1:natom) = minus_two_kext * charge(1:natom) * mu_y
   fz(1:natom) = minus_two_kext * charge(1:natom) * mu_z

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module EwaldSumExternal
