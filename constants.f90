Module Constants !DOC
!DOC !FILE The module contains the constants used in the program


    real(8), parameter :: pi=3.14159265358979323846d0   !DOC pi
    real(8), parameter :: two_pi = 2.d0*pi      !DOC 2*pi
    real(8), parameter :: four_pi = 4.d0*pi     !DOC 4*pi

    real(8), parameter :: boltz=1.38065d-23     !DOC  Boltzmann constant in Jouls
    real(8), parameter :: elec=1.602176d-19    !DOC charge of the electron in Coulombs
    real(8), parameter :: epsi0=8.854187d-12   !DOC dialectrical permutivity of vacuum [ in SI ]

    real(8), parameter :: kcal_mol=6.9477d-21   !DOC kcal/mol in Joul

    real(8) ,parameter :: LN2 = 0.6931471805599453d0 !DOC natural logarithm of 2

!    real(8), parameter :: dbjr_a = elec**2/(4.d0*pi*epsi0*boltz*temp)/1.d-10  ! Bjerum Length in Angstr
                               ! also : U = 1/(4 pi epsi0


END MODULE Constants
