Module geometry !DOC
!DOC !FILE Data structures and functions for geometrical calculations

! in Luc's convention
! i even do not know, how it works
Type :: TRotation !DOC
!DOC Rotation matrix in Luc's convention (i don't know how it works)
!DOC used in choosing the rotation of the molecule
!DOC Fields:
   real(8) :: sin_phi,cos_phi !DOC  sin(phi),cos(phi)
   real(8) :: sin_theta,cos_theta !DOC sin(theta), cos(theta)
   real(8) :: sin_angle,cos_angle  !DOC sin(angle), cos(angle)

End Type

Type :: TRotMatrix !DOC
!DOC Rotational matix 3x3
!DOC Fields:
   real(8) :: xx,xy,xz  !DOC 
   real(8) :: yx,yy,yz  !DOC 
   real(8) :: zx,zy,zz  !DOC

End Type


CONTAINS

subroutine rot_vect(x,y,z,x1,y1,z1,cc,ss,i ) !DOC
IMPLICIT REAL(8) (a-h,o-z)
!DOC rotate the coordinates (x,y,z) by angle alpha along one of 3 axes: x,y,z
!DOC        fait la rotation de x,y,z en x1,y1,z1 d'angle cc=cos,ss=sin autour d'un (i) des 3 axes principaux
!DOC        attention: x1,y1,z1 peut etre a la meme place memoire que x,y,z
!DOC        luc80p176
!
!DOC Taken from the MC code for H2O by Luc Belloni, no changes
!
!
!DOC Parameters:
!DOC :: x,y,z - coordinates of the first vector
!DOC :: x1,y1,z1 - output coordinates
!DOC :: cc,ss - cos(alpha),sin(alpha)
!DOC :: i - number of the axis: 1 is x, 2 is y, 3 is z
!

if(i==1) then
x1=x
y2=y
y1=cc*y2-ss*z
z1=ss*y2+cc*z
end if
if(i==2) then
y1=y
x2=x
x1=cc*x2+ss*z
z1=-ss*x2+cc*z
end if
if(i==3) then
z1=z
x2=x
x1=cc*x2-ss*y
y1=ss*x2+cc*y
end if
end subroutine

subroutine prod_vect(x1,y1,z1,x2,y2,z2,xprod,yprod,zprod) !DOC
!IMPLICIT REAL(8) (a-h,o-z)
!DOC vector product
!DOC Parameters:
implicit none
real(8),intent(in) :: x1,y1,z1  !DOC first vector
real(8),intent(in) ::  x2,y2,z2 !DOC second vector
real(8)    :: xprod,yprod,zprod !DOC output: result

real(8) :: x1_,y1_,z1_,x2_,y2_,z2_ ! tmp vars for the case x1==xprod etc... 

x1_ = x1
y1_ = y1
z1_ = z1

x2_ = x2
y2_ = y2
z2_ = z2

!
!        fait le produit vectoriel r = r1 * r2
!        luc84 page 160
!        le resultat peut-etre mis a la place du 1er vecteur!
!
!
! Taken from the MC code for H2O by Luc Belloni, no changes
!
!
!x3=x1; y3=y1; z3=z1
!x = y3*z2-z3*y2
!y=z3*x2-x3*z2
!z=x3*y2-y3*x2
xprod = y1_ * z2_ - z1_ * y2_
yprod = z1_ * x2_ - x1_ * z2_
zprod = x1_ * y2_ - y1_ * x2_

end subroutine

subroutine fill_rot_matrix(rot,rot_matrix) !DOC
!DOC convert rotation in Luc's notation to the rotational matrix
!DOC Parameters:
  implicit none
  Type(TRotation),intent(in) :: rot !DOC rotational structure in Luc's notation
  Type(TRotMatrix) :: rot_matrix !DOC 3x3 rotational matrix 

  real(8) :: cphi, sphi, st, ct, ca, sa, ca1

   cphi = rot % cos_phi
   sphi = rot % sin_phi

   st = rot % sin_theta
   ct = rot % cos_theta

   ca = rot % cos_angle
   sa = rot % sin_angle

   ca1 = 1 - ca

  rot_matrix % xx = (cphi**2*st**2*ca1+ca)
  rot_matrix % xy = (sphi*cphi*st**2*ca1-ct*sa)
  rot_matrix % xz = st*(ct*ca1*cphi+sa*sphi)
 
  rot_matrix % yx = (cphi*sphi*st**2*ca1+ct*sa)
  rot_matrix % yy = (sphi**2*st**2*ca1+ca )
  rot_matrix % yz = st*(ct*ca1*sphi-sa*cphi)

  rot_matrix % zx = st*(ct*ca1*cphi-sa*sphi)
  rot_matrix % zy = st*(ct*ca1*sphi+sa*cphi)
  rot_matrix % zz = (st**2*ca+ct**2)
!   xrnew=(cphi**2*st**2*ca1+ca)*xr+(sphi*cphi*st**2*ca1-ct*sa)*yr+st*(ct*ca1*cphi+sa*sphi)*zr
!   yrnew=(cphi*sphi*st**2*ca1+ct*sa)*xr+(sphi**2*st**2*ca1+ca)*yr+st*(ct*ca1*sphi-sa*cphi)*zr
!   zrnew=st*(ct*ca1*cphi-sa*sphi)*xr+st*(ct*ca1*sphi+sa*cphi)*yr+(st**2*ca+ct**2)*zr

end subroutine


subroutine rotate_vect(R,xr,yr,zr, xrnew, yrnew, zrnew ) !DOC
  implicit none
!DOC Rotate vector using the rotational matirx R
!DOC Parameters:

  Type(TRotMatrix),intent(in) :: R !DOC rotational matrix
  real(8),intent(in) :: xr,yr,zr  !DOC input coordinates
  real(8),intent(out) :: xrnew, yrnew, zrnew !DOC coordinates after rotation

  xrnew = R % xx * xr + R % xy * yr + R % xz * zr
  yrnew = R % yx * xr + R % yy * yr + R % yz * zr
  zrnew = R % zx * xr + R % zy * yr + R % zz * zr
!   cphi = rot % cos_phi
!   sphi = rot % sin_phi
!
!   st = rot % sin_theta
!   ct = rot % cos_theta
!
!   ca = rot % cos_alpha
!   sa = rot % sin_alpha
!
!   ca1 = 1 - ca
!
!   xrnew=(cphi**2*st**2*ca1+ca)*xr+(sphi*cphi*st**2*ca1-ct*sa)*yr+st*(ct*ca1*cphi+sa*sphi)*zr
!   yrnew=(cphi*sphi*st**2*ca1+ct*sa)*xr+(sphi**2*st**2*ca1+ca)*yr+st*(ct*ca1*sphi-sa*cphi)*zr
!   zrnew=st*(ct*ca1*cphi-sa*sphi)*xr+st*(ct*ca1*sphi+sa*cphi)*yr+(st**2*ca+ct**2)*zr

end subroutine 

subroutine center_of_mass(xx,yy,zz,mass,natom, x_center, y_center, z_center) !DOC
!DOC compute center of mass of the molecule
!DOC Parameters:
   implicit none  
   real(8),dimension(:),intent(in) :: xx,yy,zz,mass !DOC coordinates and masses of atoms
   integer,intent(in) :: natom !DOC number of atoms
   real(8),intent(out) :: x_center,y_center,z_center !DOC output: coordinates of the center of mass

   integer :: i 
   real(8) :: sum_mx,sum_my,sum_mz  ! SUM m(i)x,y,z(i)
   real(8) :: totalMass

   sum_mx = 0
   sum_my = 0
   sum_mz = 0
   totalMass = 0

   do i=1,natom

      sum_mx = sum_mx + mass(i) * xx(i)
      sum_my = sum_my + mass(i) * yy(i)
      sum_mz = sum_mz + mass(i) * zz(i)

      totalMass = totalMass + mass(i)

   end do

   x_center = sum_mx / totalMass
   y_center = sum_my / totalMass
   z_center = sum_mz / totalMass

end subroutine


function atan2_two_pi( y , x ) !DOC
!DOC compute arctangent from y and x and convert it to the angle from 0 to 2pi
    use constants, only : two_pi
    implicit none
    real(8) :: atan2_two_pi
    real(8) :: y,x,R
     
    R = atan2(y,x)
    if( R < 0 ) R = R + two_pi
    atan2_two_pi = R

end function
  




subroutine ZYZRotation_matrix_to_angles( R, theta, phi, psi ) !DOC
!DOC Convert the rotational matrix to angles, theta phi,psi
     use constants, only : two_pi
!DOC Parameters:  
   implicit none
     real(8),dimension(3,3),intent(in) :: R !DOC Rotation matrix
     real(8) :: theta,phi,psi !DOC output: angles

     real(8) :: cos_theta

!DOC Conventional rotation in my case (not Luc's) is:
!DOC 1) rotate over Oz by 0<psi<2pi
!DOC 2) rotate over Oy by 0<theta<pi
!DOC 3) rotate over Oz by 0<phi<2pi
!
!DOC This rotation is used in MoleculeTable for example
!DOC I refer it as "ZYZ rotation"
!DOC do not mix this with the "Luc's rotation" which is different
! 
!DOC ZYZRotation can be defined either by angles (theta,phi,psi) or by the rotation matrix
! 
!
!DOC (note, that direction (clockwise or counter-clockwise is also important)
! As I use the Luc's rot_vect function (see above), we have:
!
! Ox rotation: 
!x1=x
!y2=y
!y1=cc*y2-ss*z
!z1=ss*y2+cc*z
!
!  which gives the matrix   ( 1   0   0   )
!                           ( 0  cos -sin )  ( counter-clockwise)
!                           ( 0  sin cos  ) 
!
!  Oy rotation:
!  
!y1=y
!x2=x
!x1=cc*x2+ss*z
!z1=-ss*x2+cc*z
!
! matrix:   ( cos  0  sin )
!           ( 0    1   0  )   (clockwise)
!           (-sin  0  cos )
!
! Oz rotation:
! 
!z1=z
!x2=x
!x1=cc*x2-ss*y
!y1=ss*x2+cc*y
!
!  matrix:  (  cos -sin 0 )
!           (  sin cos  0 )  ( counter-clockwise)
!           (  0    0   1 ) 
!
!  I will stick to these orientations 
!
! The rotation matrix for this rotation is:
!
! Rot= A * B * C 
!            ( cos_phi -sin_phi 0 )       ( cos_theta  0  sin_theta )      ( cos_psi -sin_psi  0 )
!  where A = ( sin_phi cos_phi  0 )   B = ( 0          1     0      )  C = ( sin_psi  cos_psi  0 )
!            ( 0        0       1 )       ( -sin_theta 0  cos_theta )      (   0       0       1 )
!
! this results in a following matrix (checked with symmatrix_mul):
!
!  R(1,1) = cos_phi*cos_psi*cos_theta - sin_phi*sin_psi
!  R(2,1) = cos_psi*cos_theta*sin_phi + cos_phi*sin_psi
!  R(3,1) =  - cos_psi*sin_theta
!  R(1,2) =  - cos_phi*cos_theta*sin_psi - cos_psi*sin_phi
!  R(2,2) =  - cos_theta*sin_phi*sin_psi + cos_phi*cos_psi
!  R(3,2) = sin_psi*sin_theta
!  R(1,3) = cos_phi*sin_theta
!  R(2,3) = sin_phi*sin_theta
!  R(3,3) = cos_theta
!
!  which allows to calculate the angles:
!    theta = acos( R(3,3) )
!    sin_theta = sqrt( 1 - R(3,3)**2 )
 
     cos_theta = R(3,3)

     if ( abs( cos_theta - 1d0) < 1d-8 ) then !rotation over Oy is missing, B == I
!     IF (theta == 0 (i.e. cos_theta == 1) ) THEN rotation over Oy is missing, B == I
!       both rotations which left are over Oz --> i.e. one can be simply omitted, e.g. psi = 0
!       then from the first two lines
!       R(1,1) = cos_phi * 1 * 1 - sin_phi  * 0  --> cos_phi = R(1,1)
!       R(1,2) = - cos_phi * 1 * 0  - 1 * sin_phi --> sin_phi = -R(1,2) 
!       phi = atan2( sin_phi, cos_phi)  
        theta = 0d0
        psi = 0d0

        !cos_phi = R(1,1)
        !sin_phi = -R(1,2) 
        phi = atan2_two_pi( -R(1,2), R(1,1) ) ! = atan2_two_pi( sin_phi, cos_phi )
     else ! non-zero theta. use elements of R to extract phi and psi (see above): 
        theta = acos( cos_theta )    ! 0<theta<pi
        !sin_theta = sqrt( 1 - cos_theta**2 )

        ! sin_psi = R(3,2) / sin_theta     !  R(3,2) = sin_psi*sin_theta
        ! cos_psi = - R(3,1) / sin_theta   !  R(3,1) =  - cos_psi*sin_theta
        psi = atan2_two_pi( R(3,2), -R(3,1) ) ! = atan2_two_pi( sin_psi, cos_psi )

        ! sin_phi = R(2,3) / sin_theta   !  R(2,3) = sin_phi*sin_theta
        ! cos_phi = R(1,3) / sin_theta   !  R(1,3) = cos_phi*sin_theta
        phi = atan2_two_pi( R(2,3), R(1,3) ) ! =atan2_two_pi( sin_phi, cos_phi )
     end if

end subroutine 


subroutine rotation_matrix( X_before, X_after, R ) !DOC
!DOC Calculate the rotation matrix which converts X_before to X_after
!DOC i.e. X_after = R*X_before, R = X_after*X_before^-1
!DOC Parameters:
   use matrix3x3, only : matrix3x3_mul, matrix3x3_inv
   implicit none
   real(8),dimension(3,3),intent(in) :: X_before,X_after !DOC coordinates of 3 points before and after rotation
   real(8),dimension(3,3) :: R  !DOC output:  rotation_matrix

   real(8),dimension(3,3) :: invX_before
   ! R*X_before = X_after --> R = X_after*X_before^-1

   call matrix3x3_inv( X_before, invX_before )
   call matrix3x3_mul( X_after,  invX_before, R )

end subroutine 


subroutine xyz_to_angles( xx,yy,zz, xx_new,yy_new,zz_new,natom, theta, phi, psi ) !DOC
!DOC compute the rotation matrix for the given coordinates of atoms
!    use io, only : write_xyz_array
!DOC Parameters:
    implicit none
    real(8), dimension(:),intent(in) :: xx,yy,zz !DOC inital coordinates
    real(8), dimension(:),intent(in) :: xx_new, yy_new, zz_new !DOC new coordinates of atoms  
                                                                !DOC should be centered by center of mass (both: before & after )
    integer,intent(in) :: natom !DOC number of atoms
    real(8) :: theta,phi,psi !DOC output: angles which describe the rotation

    real(8),dimension(3,3) :: X,X_new,R
    real(8) :: xprod,yprod,zprod
    integer :: i,i1,i2

    ! find the first non-zero atom
    i1 = 0
    do i=1,natom
       if( abs(xx(i)) + abs(yy(i)) + abs(zz(i)) > 1d-8 ) then
           i1 = i
           exit
       end if
    end do     

    if ( i1 == 0 ) then ! not found. Normally this can be only in one case: natom == 1
       theta = 0 
       phi = 0
       psi = 0
       return
    end if

    ! find the second, non-zero and non-collinear with first atom
    i2 = 0
    do i=i1+1,natom

       if( abs(xx(i)) + abs(yy(i)) + abs(zz(i)) < 1d-8 ) cycle ! zero atom

       ! if non-zero check parralelity using the vector product 
       ! (it is zero for parralell vectors ) 
       call prod_vect( xx(i1), yy(i1), zz(i1), xx(i), yy(i), zz(i), xprod, yprod, zprod )
  
       if( abs(xprod) + abs(yprod) + abs(zprod) > 1d-8 ) then  ! the second atom is found
          i2 = i
          exit
       end if

    end do

    if ( i2 == 0 ) then ! if the second atom not found (for example 2atom or linear molecule )
                        ! use the procedure for two_atom molecules
        call xyz_to_angles_two_atoms(xx(i1),yy(i1),zz(i1),xx_new(i1),yy_new(i1),zz_new(i1), theta, phi, psi )
        return
    end if

    ! else: the full molecule with 2 non-zero atoms.
    !       however, we need 3 atoms for our procedure.
    !       the third dummy "atom" can be found as vector product of two first
    !       for old positions we already have it: xprod,yprod,zprod.
    !       so, we have the first matrix:

    ! old xyz matrix:
    ! first atom
    X(1,1) = xx(i1)
    X(2,1) = yy(i1)
    X(3,1) = zz(i1)

    ! second atom
    X(1,2) = xx(i2)
    X(2,2) = yy(i2)
    X(3,2) = zz(i2)
 
    ! third pseudo-atom
    X(1,3) = xprod
    X(2,3) = yprod
    X(3,3) = zprod

    ! new xyz matrix
    ! first atom
    X_new(1,1) = xx_new(i1)
    X_new(2,1) = yy_new(i1) 
    X_new(3,1) = zz_new(i1)

    ! second atom
    X_new(1,2) = xx_new(i2)
    X_new(2,2) = yy_new(i2)
    X_new(3,2) = zz_new(i2)

    ! third pseudo-atom
    call prod_vect(xx_new(i1),yy_new(i1),zz_new(i1),xx_new(i2),yy_new(i2),zz_new(i2),xprod,yprod,zprod)
    X_new(1,3) = xprod
    X_new(2,3) = yprod
    X_new(3,3) = zprod

       ! calculate rotation matrix
    call rotation_matrix(X,X_new,R)

    ! extract the angles from the rotation matrix
    call  ZYZRotation_matrix_to_angles( R, theta, phi, psi )

end subroutine 

pure function frac(x)
   real(8) :: frac
   real(8),intent(in) :: x

   frac = x - floor(x)

end function

subroutine xyz_to_angles_two_atoms(x,y,z,x_new,y_new,z_new, theta, phi, psi ) !DOC 
!DOC special case of the previous subroutine for the case of 2atom molecule
!DOC one atom is expected to be in the origin
!DOC only coordinates of non_zero atom 
    use constants, only : pi,two_pi
    implicit none
    real(8),intent(in) :: x,y,z,x_new,y_new,z_new !DOC coordinates of atoms  
                                                    !DOC should be centered by center of mass
    real(8) :: theta,phi,psi !DOC output angles

    real(8) :: Dxz, alpha
    real(8) :: a, beta, beta_plus_phi
        psi = 0
        ! only two angles, R = A * B where A and B are defined above (see ZYZRotation_matrix_to_angles )
        ! in that case we have: 
        !    (x)   ( cos_phi*cos_theta  - sin_phi cos_phi*sin_theta  ) (x)    (x_new)
        ! R *(y) = ( cos_theta*sin_phi    cos_phi  sin_phi*sin_theta ) (y)  = (y_new) 
        !    (z)   (  - sin_theta          0   cos_theta             ) (z)    (z_new)
        !
        !  which gives
        !   z_new = -x*sin_theta + z*cos_theta
        !
        !  dividing by  Dxz = sqrt(x^2 + z^2 ) :
        !  z_new / Dxz =  cos_alpha * sin_theta + sin_alpha * cos_theta = sin( alpha + theta ) 
        !   where cos_alpha = -x / Dxz,  sin_alpha = z / Dxz
        !     alpha = atan2( z / Dxz, -x / Dxz ) 
        !   theta = asin( z_new / Dxz ) - alpha        
        !    considering that 0 <= theta <= pi : 
        !   theta = frac ( (asin(z_new /Dxz) - alpha )/pi ) * pi

            Dxz = sqrt( x**2 + z**2 )
            alpha = atan2( z, -x )  ! = atat2( z/Dxz, -x/Dxz) 
            theta = frac( (asin(z_new/Dxz) - alpha) / pi ) * pi
        !
        !  for the first two lines we have 2 eqs: 
        !  
        !   x * cos_phi * cos_theta - y * sin_phi + z * cos_phi*sin_theta = x_new   (1)
        !   x * sin_phi * cos_theta + y * cos_phi + z * sin_phi*sin_theta = y_new   (2)
        ! 
        ! defining  a = x*cos_theta + z*sin_theta   we have:
        !   a * cos_phi - y * sin_phi = x_new
        !   a * sin_phi + y * cos_phi = y_new
        !  
           a = x * cos(theta) + z*sin(theta) ! this is not the most efficient implementation, but for the sake of clarity... 

        ! dividing both by Day = sqrt( a^2 + y^2 ) we have
        !  cos_beta * cos_phi - sin_beta * sin_phi = x_new / Day = cos( beta + phi )
        !  cos_beta * sin_phi + sin_beta * sin_phi = y_new / Day = sin( beta + phi )
        ! 
        ! where cos_beta = a / Day,  sin_beta = y / Day
        ! 
        ! which gives 
        !         beta =  atan2( y/Day , a/Day )
        !      beta_plus_phi = atan2( y_new/Day, x_new/Day)
        !      phi = frac( (beta_plus_phi - beta) / 2pi) * 2pi

!           Day = sqrt( a**2 + y**2 )
           beta = atan2(y,a) ! =  atan2( y/Day, a/Day )
           beta_plus_phi = atan2( y_new, x_new ) ! =atan2( y_new/Day, x_new/Day )
           phi = frac( (beta_plus_phi - beta) / two_pi ) * two_pi

end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

END MODULE geometry

