Module EwaldSumTails !DOC
!DOC !FILE This module contains the functions to compute the coordinate-independent parts of the Ewald sum
implicit none

!DOC The position-independent "tails" of the ewald sum:
!DOC
!DOC  Coulomb Term: 
!DOC  alpha / sqrt(pi) SUM q_i^2
!DOC
!DOC  LJ6 Term: 
!DOC  1/6 pi^1.5 alpha^3 / V * SUM_ij A_ij   -  alpha^6/12 SUM_j A_jj
!DOC
!DOC  LJ12 Term:
!DOC  1/1080 V * pi^1.5 alpha^9 SUM_ij A_ij  - alpha^12/1440 SUM_j A_jj

contains

function ewald_sum_coulomb_tail( charge , natom ) !DOC
   use parameters, only : alpha_over_sqrt_pi,xclb
!DOC computes the position-independent coulomb term
!DOC Parameters:
   real(8) :: ewald_sum_coulomb_tail 
   real(8), dimension(:), intent(in) :: charge !DOC charges array
   integer,intent(in) :: natom      !DOC number of atoms
!DOC Return value:
!DOC  alpha / sqrt(pi) SUM q_i^2

   real(8) :: S
   integer :: i
    
   S = 0 
   do i=1,natom
      S = S + charge(i)**2
   end do

   ewald_sum_coulomb_tail = -xclb * alpha_over_sqrt_pi * S 

end function

!  LJ6 Term: 
!  1/6 pi^1.5 alpha^3 / V * SUM_ij A_ij   -  alpha^6/12 SUM_j A_jj
!
!  LJ12 Term:
!  1/1080V * pi^1.5 alpha^9 SUM_ij A_ij  - alpha^12/1440 SUM_j A_jj
subroutine ewald_sum_lj_tails( comp, lj_types, lj6_tail, lj12_tail ) !DOC
    use composition, only : TComposition
    use LJTypes, only : TLJTypes
    use parameters, only :  alpha
    use constants, only : pi
!DOC Compute the position-independet terms of the Ewald LJ sums
!DOC Parameters:
    Type(TComposition),intent(in) :: comp !DOC composition of the system
    Type(TLJTypes),intent(in) :: lj_types !DOC LJ parameters for each pair of atom types
    real(8),intent(out) :: lj6_tail !DOC ouput: LJ6 term
    real(8),intent(out) :: lj12_tail !DOC output: LJ12 term

    real(8) :: alpha3,alpha6,alpha9,alpha12   ! alpha^3,alpha^6,alpha^9,alpha^12
    real(8) :: pi_sqrt_pi   ! pi^(1.5) = pi sqrt(pi)

    real(8) :: sum_ij_lj6, sum_j_lj6  ! SUM_ij A_ij, SUM_j A_jj for LJ6
    real(8) :: sum_ij_lj12, sum_j_lj12 ! SUM_ij A_ij, SUM_j_A_jj for LJ12

    integer :: ntype
    integer :: t,t1,t2
    integer :: N1N2 ! N_t1 * N_t2
    real(8) :: local_sum_lj6,local_sum_lj12

    integer :: base_offset, offset  ! base_offset = (t1-1) * ntype, offset = base_offset + t2

    alpha3 = alpha**3
    alpha6 = alpha3**2
    alpha9 = alpha3*alpha6
    alpha12 = alpha6**2
 
    pi_sqrt_pi = pi**1.5

    ! SUM_ij A_ij = SUM_t1t2 N_t1 N_t2 A_t1t2 
    ! = SUM_t1 SUM_t2>=t1 (2-delta_t1t2) N_t1 N_t2 A_t1t2
    ntype = lj_types % ntype
    base_offset = 0
    sum_ij_lj6 = 0
    sum_ij_lj12 = 0
    do t1=1,ntype
       do t2=t1,ntype

          offset = base_offset + t2

          N1N2 = comp % mol_numbers(t1) * comp % mol_numbers(t2)
          local_sum_lj6 = lj_types % LJ6Tab( offset ) * N1N2
          local_sum_lj12 = lj_types % LJ12Tab( offset ) * N1N2

          if( t1 /= t2 ) then
              local_sum_lj6 = local_sum_lj6 + local_sum_lj6  ! (2-delta_t1t2) * A_t1t2 *N1N2
              local_sum_lj12 = local_sum_lj12 + local_sum_lj12  ! (2-delta_t1t2) * A_t1t2 *N1N2
          end if
         
          sum_ij_lj6 = sum_ij_lj6 + local_sum_lj6
          sum_ij_lj12 = sum_ij_lj12 + local_sum_lj12

       end do
 
       base_offset = base_offset + ntype 

    end do
 
    ! SUM_j A_jj = SUM_t N_t A_tt
    base_offset = 0
    sum_j_lj6 = 0
    sum_j_lj12 = 0 
    do t = 1,ntype
       
       offset = base_offset + t
       local_sum_lj6 = lj_types % LJ6Tab( offset ) * comp % mol_numbers( t )
       local_sum_lj12 = lj_types % LJ12Tab( offset ) * comp % mol_numbers( t )

       sum_j_lj6 = sum_j_lj6 + local_sum_lj6
       sum_j_lj12 = sum_j_lj12 + local_sum_lj12
 
       base_offset = base_offset + ntype
    end do
    
    ! LJ6:  1/6V pi^1.5 alpha^3 * SUM_ij A_ij   -  alpha^6 / 12 SUM_j A_jj
    ! LJ12: 1/1080V pi^1.5 alpha^9 SUM_ij A_ij  - alpha^12 /1440 SUM_j A_jj
    lj6_tail = pi_sqrt_pi * alpha3 / 6.0 *  sum_ij_lj6     - alpha6 / 12.0 * sum_j_lj6
    lj12_tail = pi_sqrt_pi * alpha9 / 1080.0 * sum_ij_lj12 - alpha12 /1440.0 * sum_j_lj12

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module EwaldSumTails
