Module periodic_table !DOC
!DOC !FILE this module contains the masses for the elements in periodic table
!DOC !FILE actually, it is not used, since the mass now is given in mol file
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Periodic Table related data definitions
!> \par History
!>      none
!> \author JGH
!> modified by Maximilien Levesque, Dec. 15, 2011
! It now comes as a simple toolbox for MDFT code
! 
!   VP Sergiievskyi
!    - KIND=dp removed, element_by_name
! *****************************************************************************
!  USE precision_kinds,only : dp
  IMPLICIT NONE
  PUBLIC :: init_periodic_table, ptable, element, nelem, atom 
! *****************************************************************************
  Type :: element !DOC
!DOC Fields:
     CHARACTER ( LEN = 2 ) :: symbol !DOC
     CHARACTER ( LEN = 14 ) :: name !DOC
     INTEGER :: number !DOC
     REAL (8) :: amass !DOC average mass in formula units
     REAL (8) :: mass !DOC mass of most abundant isotope in formula units
     REAL (8) :: covalent_radius !DOC in Angstroms
     REAL (8) :: vdw_radius !DOC in Angstroms
     INTEGER :: e_conv ( 0:3 ) !DOC ?
     REAL (8) :: heat_of_formation !DOC in kcal/mol
     REAL (4) :: eht_param ( 0:3 ) !DOC in eV
     REAL (8) :: gyrom_ratio !DOC in Mhz/Tesla
     INTEGER :: gyrom_ratio_isotope !DOC isotope number corresponding with gyrom_ratio
  END TYPE element
  type atom
     character(len=4) :: symbol
     character(len=2) :: element
     character(len=14) :: name
     INTEGER :: number
     REAL (8) :: amass ! average mass in formula units
     REAL (8) :: mass ! mass of most abundant isotope in formula units
     REAL (8) :: covalent_radius ! in Angstroms
     REAL (8) :: vdw_radius ! in Angstroms
     INTEGER :: e_conv ( 0:3 )
  end type atom
  INTEGER, PARAMETER :: nelem = 106
  TYPE ( element ) :: ptable ( 0:nelem ) !DOC periodic table array
  REAL(4), PARAMETER :: z=0.0, z1=2.015940, z2=2.015650
  INTEGER :: table_is_initialized=0
CONTAINS
! *****************************************************************************
!> \brief Initialization of Periodic Table related data
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
 function element_by_name(el) !DOC
!DOC return the number of element by its name
   integer :: element_by_name 
!DOC Parameters:
   character(2),intent(in) :: el !DOC name of the element

   integer :: i

   do i=1,nelem
      write(*,*) '"',el,'"=="',ptable(i) % symbol,'"?'  
      if ( el == ptable(i) % symbol  ) then
         element_by_name = i
         return
      end if

   end do

   element_by_name = -1

end function


SUBROUTINE init_periodic_table() !DOC
!DOC initialize the table
  if (table_is_initialized.eq.1) return
! Dummy
  ptable(0) % symbol = 'X '
  ptable(0) % name = 'Dummy'
  ptable(0) % number = 0
  ptable(0) % amass = z1
  ptable(0) % mass = z2
  ptable(0) % covalent_radius = z
  ptable(0) % vdw_radius = z
  ptable(0) % e_conv(0:3) = (/ 0, 0, 0, 0 /)
  ptable(0) % eht_param(0:3) = (/ z, z, z, z/)
! Hydrogen
  ptable(1) % symbol = 'H '
  ptable(1) % name = 'Hydrogen'
  ptable(1) % number = 1
  ptable(1) % amass = 1.00797
  ptable(1) % mass = 1.007825
  ptable(1) % covalent_radius = 0.32
  ptable(1) % vdw_radius = 1.2
  ptable(1) % e_conv(0:3) = (/ 1, 0, 0, 0 /)
  ptable(1) % eht_param(0:3) = (/ -13.60, z, z, z/)
! Deuterium
!   ptable(1) % symbol = 'D '
!   ptable(1) % name = 'Deuterium'
!   ptable(1) % number = 107
!   ptable(1) % amass = 2.0 * 1.00797
!   ptable(1) % mass = 2.0 * 1.007825
!   ptable(1) % covalent_radius = 0.32
!   ptable(1) % vdw_radius = 1.2
!   ptable(1) % e_conv(0:3) = (/ 1, 0, 0, 0 /)
!   ptable(1) % eht_param(0:3) = (/ -13.60, z, z, z/)
! Helium
  ptable(2) % symbol = 'He'
  ptable(2) % name = 'Helium'
  ptable(2) % number = 2
  ptable(2) % amass = 4.00260
  ptable(2) % mass = 4.00260
  ptable(2) % covalent_radius = 0.9300
  ptable(2) % vdw_radius = z
  ptable(2) % e_conv(0:3) = (/ 2, 0, 0, 0 /)
  ptable(2) % eht_param(0:3) = (/ -23.40, z, z, z/)
! Lithium
  ptable(3) % symbol = 'Li'
  ptable(3) % name = 'Lithium'
  ptable(3) % number = 3
  ptable(3) % amass = 6.93900
  ptable(3) % mass = 7.01600
  ptable(3) % covalent_radius = 1.2300
  ptable(3) % vdw_radius = z
  ptable(3) % e_conv(0:3) = (/ 3, 0, 0, 0 /)
  ptable(3) % eht_param(0:3) = (/ -5.40, -3.50, z, z/)
! Beryllium
  ptable(4) % symbol = 'Be'
  ptable(4) % name = 'Beryllium'
  ptable(4) % number = 4
  ptable(4) % amass = 9.01220
  ptable(4) % mass = 9.01218
  ptable(4) % covalent_radius = 0.9000
  ptable(4) % vdw_radius = z
  ptable(4) % e_conv(0:3) = (/ 4, 0, 0, 0 /)
  ptable(4) % eht_param(0:3) = (/ -10.00, -6.00, z, z/)
! Boron
  ptable(5) % symbol = 'B '
  ptable(5) % name = 'Boron'
  ptable(5) % number = 5
  ptable(5) % amass = 10.81100
  ptable(5) % mass = 11.00931
  ptable(5) % covalent_radius = 0.8200
  ptable(5) % vdw_radius = z
  ptable(5) % e_conv(0:3) = (/ 4, 1, 0, 0 /)
  ptable(5) % eht_param(0:3) = (/ -15.20, -8.50, z, z/)
! Carbon
  ptable(6) % symbol = 'C '
  ptable(6) % name = 'Carbon'
  ptable(6) % number = 6
  ptable(6) % amass = 12.01115
  ptable(6) % mass = 12.0000
  ptable(6) % covalent_radius = 0.7700
  ptable(6) % vdw_radius = 1.7           ! obtained from www.webelement.com
  ptable(6) % e_conv(0:3) = (/ 4, 2, 0, 0 /)
  ptable(6) % eht_param(0:3) = (/ -21.40, -11.40, z, z/)
! Nitrogen
  ptable(7) % symbol = 'N '
  ptable(7) % name = 'Nitrogen'
  ptable(7) % number = 7
  ptable(7) % amass = 14.00670
  ptable(7) % mass = 14.00307
  ptable(7) % covalent_radius = 0.7500
  ptable(7) % vdw_radius = 1.5
  ptable(7) % e_conv(0:3) = (/ 4, 3, 0, 0 /)
  ptable(7) % eht_param(0:3) = (/ -26.00, -13.40, z, z/)
! Oxygen
  ptable(8) % symbol = 'O '
  ptable(8) % name = 'Oxygen'
  ptable(8) % number = 8
  ptable(8) % amass = 15.99940
  ptable(8) % mass = 15.99491
  ptable(8) % covalent_radius = 0.7300
  ptable(8) % vdw_radius = 1.40
  ptable(8) % e_conv(0:3) = (/ 4, 4, 0, 0 /)
  ptable(8) % eht_param(0:3) = (/ -32.30, -14.80, z, z/)
! Fluorine
  ptable(9) % symbol = 'F '
  ptable(9) % name = 'Fluorine'
  ptable(9) % number = 9
  ptable(9) % amass = 18.99840
  ptable(9) % mass = 18.99840
  ptable(9) % covalent_radius = 0.7200
  ptable(9) % vdw_radius = 1.35
  ptable(9) % e_conv(0:3) = (/ 4, 5, 0, 0 /)
  ptable(9) % eht_param(0:3) = (/ -40.00, -18.10, z, z/)
! Neon
  ptable(10) % symbol = 'Ne'
  ptable(10) % name = 'Neon'
  ptable(10) % number = 10
  ptable(10) % amass = 20.18300
  ptable(10) % mass = 19.99244
  ptable(10) % covalent_radius = 0.7100
  ptable(10) % vdw_radius = z
  ptable(10) % e_conv(0:3) = (/ 4, 6, 0, 0 /)
  ptable(10) % eht_param(0:3) = (/ -43.20, -20.00, z, z/)
! Sodium
  ptable(11) % symbol = 'Na'
  ptable(11) % name = 'Sodium'
  ptable(11) % number = 11
  ptable(11) % amass = 22.98980
  ptable(11) % mass = 22.9898
  ptable(11) % covalent_radius = 1.5400
  ptable(11) % vdw_radius = z
  ptable(11) % e_conv(0:3) = (/ 5, 6, 0, 0 /)
  ptable(11) % eht_param(0:3) = (/ -5.10, -3.00, z, z/)
! Magnesium
  ptable(12) % symbol = 'Mg'
  ptable(12) % name = 'Magnesium'
  ptable(12) % number = 12
  ptable(12) % amass = 24.31200
  ptable(12) % mass = 23.98504
  ptable(12) % covalent_radius = 1.3600
  ptable(12) % vdw_radius = z
  ptable(12) % e_conv(0:3) = (/ 6, 6, 0, 0 /)
  ptable(12) % eht_param(0:3) = (/ -9.00, -4.50, z, z/)
! Aluminium
  ptable(13) % symbol = 'Al'
  ptable(13) % name = 'Aluminium'
  ptable(13) % number = 13
  ptable(13) % amass = 26.98153
  ptable(13) % mass = 26.98153
  ptable(13) % covalent_radius = 1.1800
  ptable(13) % vdw_radius = z
  ptable(13) % e_conv(0:3) = (/ 6, 7, 0, 0 /)
  ptable(13) % eht_param(0:3) = (/ -12.30, -6.50, z, z/)
! Silicon
  ptable(14) % symbol = 'Si'
  ptable(14) % name = 'Silicon'
  ptable(14) % number = 14
  ptable(14) % amass = 28.08600
  ptable(14) % mass = 27.97693
  ptable(14) % covalent_radius = 1.1100
  ptable(14) % vdw_radius = z
  ptable(14) % e_conv(0:3) = (/ 6, 8, 0, 0 /)
  ptable(14) % eht_param(0:3) = (/ -17.30, -9.20, z, z/)
! Phosphorus
  ptable(15) % symbol = 'P '
  ptable(15) % name = 'Phosphorus'
  ptable(15) % number = 15
  ptable(15) % amass = 30.97380
  ptable(15) % mass = 30.97376
  ptable(15) % covalent_radius = 1.0600
  ptable(15) % vdw_radius = 1.9
  ptable(15) % e_conv(0:3) = (/ 6, 9, 0, 0 /)
  ptable(15) % eht_param(0:3) = (/ -18.60, -14.00, z, z/)
! Sulfur
  ptable(16) % symbol = 'S '
  ptable(16) % name = 'Sulfur'
  ptable(16) % number = 16
  ptable(16) % amass = 32.06400
  ptable(16) % mass = 31.97207
  ptable(16) % covalent_radius = 1.0200
  ptable(16) % vdw_radius = 1.85
  ptable(16) % e_conv(0:3) = (/ 6, 10, 0, 0 /)
  ptable(16) % eht_param(0:3) = (/ -20.00, -11.00, z, z/)
! Chlorine
  ptable(17) % symbol = 'Cl'
  ptable(17) % name = 'Chlorine'
  ptable(17) % number = 17
  ptable(17) % amass = 35.45300
  ptable(17) % mass = 34.96885
  ptable(17) % covalent_radius = 0.9900
  ptable(17) % vdw_radius = 1.80
  ptable(17) % e_conv(0:3) = (/ 6, 11, 0, 0 /)
  ptable(17) % eht_param(0:3) = (/ -26.30, -14.20, z, z/)
! Argon
  ptable(18) % symbol = 'Ar'
  ptable(18) % name = 'Argon'
  ptable(18) % number = 18
  ptable(18) % amass = 39.94800
  ptable(18) % mass = 39.94800
  ptable(18) % covalent_radius = 0.9800
  ptable(18) % vdw_radius = 3.83
  ptable(18) % e_conv(0:3) = (/ 6, 12, 0, 0 /)
  ptable(18) % eht_param(0:3) = (/ z, z, z, z/)
! Potassium
  ptable(19) % symbol = 'K '
  ptable(19) % name = 'Potassium'
  ptable(19) % number = 19
  ptable(19) % amass = 39.10200
  ptable(19) % mass = 38.96371
  ptable(19) % covalent_radius = 2.0300
  ptable(19) % vdw_radius = z
  ptable(19) % e_conv(0:3) = (/ 7, 12, 0, 0 /)
  ptable(19) % eht_param(0:3) = (/ -4.34, -2.73, z, z/)
! Calcium
  ptable(20) % symbol = 'Ca'
  ptable(20) % name = 'Calcium'
  ptable(20) % number = 20
  ptable(20) % amass = 40.08000
  ptable(20) % mass = 39.96259
  ptable(20) % covalent_radius = 1.7400
  ptable(20) % vdw_radius = z
  ptable(20) % e_conv(0:3) = (/ 8, 12, 0, 0 /)
  ptable(20) % eht_param(0:3) = (/ -7.00, -4.00, z, z/)
! Scandium
  ptable(21) % symbol = 'Sc'
  ptable(21) % name = 'Scandium'
  ptable(21) % number = 21
  ptable(21) % amass = 44.95600
  ptable(21) % mass = 44.95592
  ptable(21) % covalent_radius = 1.4400
  ptable(21) % vdw_radius = z
  ptable(21) % e_conv(0:3) = (/ 8, 12, 1, 0 /)
  ptable(21) % eht_param(0:3) = (/ -8.87, -2.75, -8.51, z/)
! Titanium
  ptable(22) % symbol = 'Ti'
  ptable(22) % name = 'Titanium'
  ptable(22) % number = 22
  ptable(22) % amass = 47.90000
  ptable(22) % mass = 48.00000
  ptable(22) % covalent_radius = 1.3200
  ptable(22) % vdw_radius = z
  ptable(22) % e_conv(0:3) = (/ 8, 12, 2, 0 /)
  ptable(22) % eht_param(0:3) = (/ -8.97, -5.44, -10.81, z/)
! Vanadium
  ptable(23) % symbol = 'V '
  ptable(23) % name = 'Vanadium'
  ptable(23) % number = 23
  ptable(23) % amass = 50.94200
  ptable(23) % mass = 50.94400
  ptable(23) % covalent_radius = 1.2200
  ptable(23) % vdw_radius = z
  ptable(23) % e_conv(0:3) = (/ 8, 12, 3, 0 /)
  ptable(23) % eht_param(0:3) = (/ -8.81, -5.52, -11.00, z/)
! Chromium
  ptable(24) % symbol = 'Cr'
  ptable(24) % name = 'Chromium'
  ptable(24) % number = 24
  ptable(24) % amass = 51.99600
  ptable(24) % mass = 51.94050
  ptable(24) % covalent_radius = 1.1800
  ptable(24) % vdw_radius = z
  ptable(24) % e_conv(0:3) = (/ 7, 12, 5, 0 /)
  ptable(24) % eht_param(0:3) = (/ -8.66, -5.24, -11.22, z/)
! Manganese
  ptable(25) % symbol = 'Mn'
  ptable(25) % name = 'Manganese'
  ptable(25) % number = 25
  ptable(25) % amass = 54.93800
  ptable(25) % mass = 54.93810
  ptable(25) % covalent_radius = 1.1700
  ptable(25) % vdw_radius = z
  ptable(25) % e_conv(0:3) = (/ 8, 12, 5, 0 /)
  ptable(25) % eht_param(0:3) = (/ -9.75, -5.89, -11.67, z/)
! Iron
  ptable(26) % symbol = 'Fe'
  ptable(26) % name = 'Iron'
  ptable(26) % number = 26
  ptable(26) % amass = 55.84700
  ptable(26) % mass = 55.93490
  ptable(26) % covalent_radius = 1.1700
  ptable(26) % vdw_radius = z
  ptable(26) % e_conv(0:3) = (/ 8, 12, 6, 0 /)
  ptable(26) % eht_param(0:3) = (/ -9.10, -5.32, -12.60, z/)
! Cobalt
  ptable(27) % symbol = 'Co'
  ptable(27) % name = 'Cobalt'
  ptable(27) % number = 27
  ptable(27) % amass = 58.93300
  ptable(27) % mass = 58.93320
  ptable(27) % covalent_radius = 1.1600
  ptable(27) % vdw_radius = z
  ptable(27) % e_conv(0:3) = (/ 8, 12, 7, 0 /)
  ptable(27) % eht_param(0:3) = (/ -9.21, -5.29, -13.18, z/)
! Nickel
  ptable(28) % symbol = 'Ni'
  ptable(28) % name = 'Nickel'
  ptable(28) % number = 28
  ptable(28) % amass = 58.71000
  ptable(28) % mass = 57.93530
  ptable(28) % covalent_radius = 1.1500
  ptable(28) % vdw_radius = z
  ptable(28) % e_conv(0:3) = (/ 8, 12, 8, 0 /)
  ptable(28) % eht_param(0:3) = (/ -9.17, -5.15, -13.49, z/)
! Copper
  ptable(29) % symbol = 'Cu'
  ptable(29) % name = 'Copper'
  ptable(29) % number = 29
  ptable(29) % amass = 63.54000
  ptable(29) % mass = 63.20000
  ptable(29) % covalent_radius = 1.1700
  ptable(29) % vdw_radius = z
  ptable(29) % e_conv(0:3) = (/ 7, 12, 10, 0 /)
  ptable(29) % eht_param(0:3) = (/ -11.40, -6.06, -14.00, z/)
! Zinc
  ptable(30) % symbol = 'Zn'
  ptable(30) % name = 'Zinc'
  ptable(30) % number = 30
  ptable(30) % amass = 65.37000
  ptable(30) % mass = 63.92910
  ptable(30) % covalent_radius = 1.2500
  ptable(30) % vdw_radius = z
  ptable(30) % e_conv(0:3) = (/ 8, 12, 10, 0 /)
  ptable(30) % eht_param(0:3) = (/ -12.41, -6.53, z, z/)
! Gallium
  ptable(31) % symbol = 'Ga'
  ptable(31) % name = 'Gallium'
  ptable(31) % number = 31
  ptable(31) % amass = 69.72000
  ptable(31) % mass = 68.92570
  ptable(31) % covalent_radius = 1.2600
  ptable(31) % vdw_radius = z
  ptable(31) % e_conv(0:3) = (/ 8, 13, 10, 0 /)
  ptable(31) % eht_param(0:3) = (/ -14.58, -6.75, z, z/)
! Germanium
  ptable(32) % symbol = 'Ge'
  ptable(32) % name = 'Germanium'
  ptable(32) % number = 32
  ptable(32) % amass = 72.59000
  ptable(32) % mass = 73.92140
  ptable(32) % covalent_radius = 1.2200
  ptable(32) % vdw_radius = z
  ptable(32) % e_conv(0:3) = (/ 8, 14, 10, 0 /)
  ptable(32) % eht_param(0:3) = (/ -16.00, -9.00, z, z/)
! Arsenic
  ptable(33) % symbol = 'As'
  ptable(33) % name = 'Arsenic'
  ptable(33) % number = 33
  ptable(33) % amass = 74.92200
  ptable(33) % mass = 74.92160
  ptable(33) % covalent_radius = 1.2000
  ptable(33) % vdw_radius = 2.0
  ptable(33) % e_conv(0:3) = (/ 8, 15, 10, 0 /)
  ptable(33) % eht_param(0:3) = (/ -16.22, -12.16, z, z/)
! Selenium
  ptable(34) % symbol = 'Se'
  ptable(34) % name = 'Selenium'
  ptable(34) % number = 34
  ptable(34) % amass = 78.96000
  ptable(34) % mass = 79.91650
  ptable(34) % covalent_radius = 1.1600
  ptable(34) % vdw_radius = 2.00
  ptable(34) % e_conv(0:3) = (/ 8, 16, 10, 0 /)
  ptable(34) % eht_param(0:3) = (/ -20.50, -14.40, z, z/)
! Bromine
  ptable(35) % symbol = 'Br'
  ptable(35) % name = 'Bromine'
  ptable(35) % number = 35
  ptable(35) % amass = 79.90900
  ptable(35) % mass = 78.91830
  ptable(35) % covalent_radius = 1.1400
  ptable(35) % vdw_radius = 1.95
  ptable(35) % e_conv(0:3) = (/ 8, 17, 10, 0 /)
  ptable(35) % eht_param(0:3) = (/ -22.07, -13.10, z, z/)
! Krypton
  ptable(36) % symbol = 'Kr'
  ptable(36) % name = 'Krypton'
  ptable(36) % number = 36
  ptable(36) % amass = 83.80000
  ptable(36) % mass = 84.00000
  ptable(36) % covalent_radius = 1.1200
  ptable(36) % vdw_radius = z
  ptable(36) % e_conv(0:3) = (/ 8, 18, 10, 0 /)
  ptable(36) % eht_param(0:3) = (/ z, z, z, z/)
! Rubidium
  ptable(37) % symbol = 'Rb'
  ptable(37) % name = 'Rubidium'
  ptable(37) % number = 37
  ptable(37) % amass = 85.47000
  ptable(37) % mass = 84.91170
  ptable(37) % covalent_radius = 2.1600
  ptable(37) % vdw_radius = z
  ptable(37) % e_conv(0:3) = (/ 9, 18, 10, 0 /)
  ptable(37) % eht_param(0:3) = (/ -4.18, -2.60, z, z/)
! Strontium
  ptable(38) % symbol = 'Sr'
  ptable(38) % name = 'Strontium'
  ptable(38) % number = 38
  ptable(38) % amass = 87.62000
  ptable(38) % mass = 87.90560
  ptable(38) % covalent_radius = 1.9100
  ptable(38) % vdw_radius = z
  ptable(38) % e_conv(0:3) = (/ 10, 18, 10, 0 /)
  ptable(38) % eht_param(0:3) = (/ -6.62, -3.92, z, z/)
! Yttrium
  ptable(39) % symbol = 'Y '
  ptable(39) % name = 'Yttrium'
  ptable(39) % number = 39
  ptable(39) % amass = 88.90500
  ptable(39) % mass = 88.90590
  ptable(39) % covalent_radius = 1.6200
  ptable(39) % vdw_radius = z
  ptable(39) % e_conv(0:3) = (/ 10, 18, 11, 0 /)
  ptable(39) % eht_param(0:3) = (/ z, z, z, z/)
! Zirconium
  ptable(40) % symbol = 'Zr'
  ptable(40) % name = 'Zirconium'
  ptable(40) % number = 40
  ptable(40) % amass = 91.22000
  ptable(40) % mass = 89.90430
  ptable(40) % covalent_radius = 1.4500
  ptable(40) % vdw_radius = z
  ptable(40) % e_conv(0:3) = (/ 10, 18, 12, 0 /)
  ptable(40) % eht_param(0:3) = (/ -8.00, -5.40, -10.20, z/)
! Niobium
  ptable(41) % symbol = 'Nb'
  ptable(41) % name = 'Niobium'
  ptable(41) % number = 41
  ptable(41) % amass = 92.90600
  ptable(41) % mass = 92.90600
  ptable(41) % covalent_radius = 1.3400
  ptable(41) % vdw_radius = z
  ptable(41) % e_conv(0:3) = (/ 9, 18, 14, 0 /)
  ptable(41) % eht_param(0:3) = (/ -10.10, -6.86, -12.10, z/)
! Molybdenum
  ptable(42) % symbol = 'Mo'
  ptable(42) % name = 'Molybdenum'
  ptable(42) % number = 42
  ptable(42) % amass = 95.94000
  ptable(42) % mass = 97.90550
  ptable(42) % covalent_radius = 1.3000
  ptable(42) % vdw_radius = z
  ptable(42) % e_conv(0:3) = (/ 9, 18, 15, 0 /)
  ptable(42) % eht_param(0:3) = (/ -8.34, -5.25, -10.50, z/)
! Technetium
  ptable(43) % symbol = 'Tc'
  ptable(43) % name = 'Technetium'
  ptable(43) % number = 43
  ptable(43) % amass = 98.90600
  ptable(43) % mass = 98.90600
  ptable(43) % covalent_radius = 1.2700
  ptable(43) % vdw_radius = z
  ptable(43) % e_conv(0:3) = (/ 9, 18, 16, 0 /)
  ptable(43) % eht_param(0:3) = (/ -10.07, -5.40, -12.82, z/)
! Ruthenium
  ptable(44) % symbol = 'Ru'
  ptable(44) % name = 'Ruthenium'
  ptable(44) % number = 44
  ptable(44) % amass = 101.07000
  ptable(44) % mass = 101.90370
  ptable(44) % covalent_radius = 1.2500
  ptable(44) % vdw_radius = z
  ptable(44) % e_conv(0:3) = (/ 9, 18, 17, 0 /)
  ptable(44) % eht_param(0:3) = (/ -10.40, -6.87, -14.90, z/)
! Rhodium
  ptable(45) % symbol = 'Rh'
  ptable(45) % name = 'Rhodium'
  ptable(45) % number = 45
  ptable(45) % amass = 102.90500
  ptable(45) % mass = 102.90480
  ptable(45) % covalent_radius = 1.2500
  ptable(45) % vdw_radius = z
  ptable(45) % e_conv(0:3) = (/ 9, 18, 18, 0 /)
  ptable(45) % eht_param(0:3) = (/ -8.09, -4.57, -12.50, z/)
! Palladium
  ptable(46) % symbol = 'Pd'
  ptable(46) % name = 'Palladium'
  ptable(46) % number = 46
  ptable(46) % amass = 106.40000
  ptable(46) % mass = 105.90320
  ptable(46) % covalent_radius = 1.2800
  ptable(46) % vdw_radius = z
  ptable(46) % e_conv(0:3) = (/ 8, 18, 20, 0 /)
  ptable(46) % eht_param(0:3) = (/ -7.32, -3.75, -12.02, z/)
! Silver
  ptable(47) % symbol = 'Ag'
  ptable(47) % name = 'Silver'
  ptable(47) % number = 47
  ptable(47) % amass = 107.87000
  ptable(47) % mass = 106.90509
  ptable(47) % covalent_radius = 1.3400
  ptable(47) % vdw_radius = z
  ptable(47) % e_conv(0:3) = (/ 9, 18, 20, 0 /)
  ptable(47) % eht_param(0:3) = (/ z, z, z, z/)
! Cadmium
  ptable(48) % symbol = 'Cd'
  ptable(48) % name = 'Cadmium'
  ptable(48) % number = 48
  ptable(48) % amass = 112.40000
  ptable(48) % mass = 113.90360
  ptable(48) % covalent_radius = 1.4800
  ptable(48) % vdw_radius = z
  ptable(48) % e_conv(0:3) = (/ 10, 18, 20, 0 /)
  ptable(48) % eht_param(0:3) = (/ z, z, z, z/)
! Indium
  ptable(49) % symbol = 'In'
  ptable(49) % name = 'Indium'
  ptable(49) % number = 49
  ptable(49) % amass = 114.82000
  ptable(49) % mass = 114.90410
  ptable(49) % covalent_radius = 1.4400
  ptable(49) % vdw_radius = z
  ptable(49) % e_conv(0:3) = (/ 10, 19, 20, 0 /)
  ptable(49) % eht_param(0:3) = (/ -12.60, -6.19, z, z/)
! Tin
  ptable(50) % symbol = 'Sn'
  ptable(50) % name = 'Tin'
  ptable(50) % number = 50
  ptable(50) % amass = 118.69000
  ptable(50) % mass = 120.00000
  ptable(50) % covalent_radius = 1.4100
  ptable(50) % vdw_radius = z
  ptable(50) % e_conv(0:3) = (/ 10, 20, 20, 0 /)
  ptable(50) % eht_param(0:3) = (/ -16.16, -8.32, z, z/)
! Antimony
  ptable(51) % symbol = 'Sb'
  ptable(51) % name = 'Antimony'
  ptable(51) % number = 51
  ptable(51) % amass = 121.75000
  ptable(51) % mass = 120.90380
  ptable(51) % covalent_radius = 1.4000
  ptable(51) % vdw_radius = 2.2
  ptable(51) % e_conv(0:3) = (/ 10, 21, 20, 0 /)
  ptable(51) % eht_param(0:3) = (/ -18.80, -11.70, z, z/)
! Tellurium
  ptable(52) % symbol = 'Te'
  ptable(52) % name = 'Tellurium'
  ptable(52) % number = 52
  ptable(52) % amass = 127.60000
  ptable(52) % mass = 129.90670
  ptable(52) % covalent_radius = 1.3600
  ptable(52) % vdw_radius = 2.20
  ptable(52) % e_conv(0:3) = (/ 10, 22, 20, 0 /)
  ptable(52) % eht_param(0:3) = (/ -20.80, -13.20, z, z/)
! Iodine
  ptable(53) % symbol = 'I '
  ptable(53) % name = 'Iodine'
  ptable(53) % number = 53
  ptable(53) % amass = 126.90440
  ptable(53) % mass = 126.90440
  ptable(53) % covalent_radius = 1.3300
  ptable(53) % vdw_radius = 2.15
  ptable(53) % e_conv(0:3) = (/ 10, 23, 20, 0 /)
  ptable(53) % eht_param(0:3) = (/ -18.00, -12.70, z, z/)
! Xenon
  ptable(54) % symbol = 'Xe'
  ptable(54) % name = 'Xenon'
  ptable(54) % number = 54
  ptable(54) % amass = 131.30000
  ptable(54) % mass = 131.90420
  ptable(54) % covalent_radius = 1.3100
  ptable(54) % vdw_radius = z
  ptable(54) % e_conv(0:3) = (/ 10, 24, 20, 0 /)
  ptable(54) % eht_param(0:3) = (/ z, z, z, z/)
! Cesium
  ptable(55) % symbol = 'Cs'
  ptable(55) % name = 'Cesium'
  ptable(55) % number = 55
  ptable(55) % amass = 132.90500
  ptable(55) % mass = 132.90510
  ptable(55) % covalent_radius = 2.3500
  ptable(55) % vdw_radius = z
  ptable(55) % e_conv(0:3) = (/ 11, 24, 20, 0 /)
  ptable(55) % eht_param(0:3) = (/ -3.88, -2.49, z, z/)
! Barium
  ptable(56) % symbol = 'Ba'
  ptable(56) % name = 'Barium'
  ptable(56) % number = 56
  ptable(56) % amass = 137.34000
  ptable(56) % mass = 137.90500
  ptable(56) % covalent_radius = 1.9800
  ptable(56) % vdw_radius = z
  ptable(56) % e_conv(0:3) = (/ 12, 24, 20, 0 /)
  ptable(56) % eht_param(0:3) = (/ z, z, z, z/)
! Lantanum
  ptable(57) % symbol = 'La'
  ptable(57) % name = 'Lantanum'
  ptable(57) % number = 57
  ptable(57) % amass = 138.91000
  ptable(57) % mass = 138.90610
  ptable(57) % covalent_radius = 1.6900
  ptable(57) % vdw_radius = z
  ptable(57) % e_conv(0:3) = (/ 12, 24, 21, 0 /)
  ptable(57) % eht_param(0:3) = (/ -7.67, -5.01, -8.21, z/)
! Cerium
  ptable(58) % symbol = 'Ce'
  ptable(58) % name = 'Cerium'
  ptable(58) % number = 58
  ptable(58) % amass = 140.12000
  ptable(58) % mass = 139.90530
  ptable(58) % covalent_radius = 1.6500
  ptable(58) % vdw_radius = z
  ptable(58) % e_conv(0:3) = (/ 12, 24, 20, 2 /)
  ptable(58) % eht_param(0:3) = (/ z, z, z, z/)
! Praseodymium
  ptable(59) % symbol = 'Pr'
  ptable(59) % name = 'Praseodymium'
  ptable(59) % number = 59
  ptable(59) % amass = 140.90700
  ptable(59) % mass = 140.90740
  ptable(59) % covalent_radius = 1.6500
  ptable(59) % vdw_radius = z
  ptable(59) % e_conv(0:3) = (/ 12, 24, 20, 3 /)
  ptable(59) % eht_param(0:3) = (/ z, z, z, z/)
! Neodymium
  ptable(60) % symbol = 'Nd'
  ptable(60) % name = 'Neodymium'
  ptable(60) % number = 60
  ptable(60) % amass = 144.24000
  ptable(60) % mass = 141.90750
  ptable(60) % covalent_radius = 1.6400
  ptable(60) % vdw_radius = z
  ptable(60) % e_conv(0:3) = (/ 12, 24, 20, 4 /)
  ptable(60) % eht_param(0:3) = (/ z, z, z, z/)
! Promethium
  ptable(61) % symbol = 'Pm'
  ptable(61) % name = 'Promethium'
  ptable(61) % number = 61
  ptable(61) % amass = 144.91300
  ptable(61) % mass = 144.91300
  ptable(61) % covalent_radius = 1.6300
  ptable(61) % vdw_radius = z
  ptable(61) % e_conv(0:3) = (/ 12, 24, 20, 5 /)
  ptable(61) % eht_param(0:3) = (/ z, z, z, z/)
! Samarium
  ptable(62) % symbol = 'Sm'
  ptable(62) % name = 'Samarium'
  ptable(62) % number = 62
  ptable(62) % amass = 150.35000
  ptable(62) % mass = 151.91950
  ptable(62) % covalent_radius = 1.6200
  ptable(62) % vdw_radius = z
  ptable(62) % e_conv(0:3) = (/ 12, 24, 20, 6 /)
  ptable(62) % eht_param(0:3) = (/ -4.86, -4.86, -6.06, &
          -11.28/)
! Europium
  ptable(63) % symbol = 'Eu'
  ptable(63) % name = 'Europium'
  ptable(63) % number = 63
  ptable(63) % amass = 151.96000
  ptable(63) % mass = 152.92090
  ptable(63) % covalent_radius = 1.8500
  ptable(63) % vdw_radius = z
  ptable(63) % e_conv(0:3) = (/ 12, 24, 20, 7 /)
  ptable(63) % eht_param(0:3) = (/ z, z, z, z/)
! Gadolinium
  ptable(64) % symbol = 'Gd'
  ptable(64) % name = 'Gadolinium'
  ptable(64) % number = 64
  ptable(64) % amass = 157.25000
  ptable(64) % mass = 157.92410
  ptable(64) % covalent_radius = 1.6100
  ptable(64) % vdw_radius = z
  ptable(64) % e_conv(0:3) = (/ 12, 24, 21, 7 /)
  ptable(64) % eht_param(0:3) = (/ z, z, z, z/)
! Terbium
  ptable(65) % symbol = 'Tb'
  ptable(65) % name = 'Terbium'
  ptable(65) % number = 65
  ptable(65) % amass = 158.92400
  ptable(65) % mass = 158.92500
  ptable(65) % covalent_radius = 1.5900
  ptable(65) % vdw_radius = z
  ptable(65) % e_conv(0:3) = (/ 12, 24, 20, 9 /)
  ptable(65) % eht_param(0:3) = (/ z, z, z, z/)
! Dysprosium
  ptable(66) % symbol = 'Dy'
  ptable(66) % name = 'Dysprosium'
  ptable(66) % number = 66
  ptable(66) % amass = 162.50000
  ptable(66) % mass = 163.92880
  ptable(66) % covalent_radius = 1.5900
  ptable(66) % vdw_radius = z
  ptable(66) % e_conv(0:3) = (/ 12, 24, 20, 10 /)
  ptable(66) % eht_param(0:3) = (/ z, z, z, z/)
! Holmium
  ptable(67) % symbol = 'Ho'
  ptable(67) % name = 'Holmium'
  ptable(67) % number = 67
  ptable(67) % amass = 164.93000
  ptable(67) % mass = 164.93000
  ptable(67) % covalent_radius = 1.5800
  ptable(67) % vdw_radius = z
  ptable(67) % e_conv(0:3) = (/ 12, 24, 20, 11 /)
  ptable(67) % eht_param(0:3) = (/ z, z, z, z/)
! Erbium
  ptable(68) % symbol = 'Er'
  ptable(68) % name = 'Erbium'
  ptable(68) % number = 68
  ptable(68) % amass = 167.26000
  ptable(68) % mass = 165.93040
  ptable(68) % covalent_radius = 1.5700
  ptable(68) % vdw_radius = z
  ptable(68) % e_conv(0:3) = (/ 12, 24, 20, 12 /)
  ptable(68) % eht_param(0:3) = (/ z, z, z, z/)
! Thulium
  ptable(69) % symbol = 'Tm'
  ptable(69) % name = 'Thulium'
  ptable(69) % number = 69
  ptable(69) % amass = 168.93400
  ptable(69) % mass = 168.93440
  ptable(69) % covalent_radius = 1.5600
  ptable(69) % vdw_radius = z
  ptable(69) % e_conv(0:3) = (/ 12, 24, 20, 13 /)
  ptable(69) % eht_param(0:3) = (/ z, z, z, z/)
! Ytterbium
  ptable(70) % symbol = 'Yb'
  ptable(70) % name = 'Ytterbium'
  ptable(70) % number = 70
  ptable(70) % amass = 173.04000
  ptable(70) % mass = 173.93900
  ptable(70) % covalent_radius = 1.5600
  ptable(70) % vdw_radius = z
  ptable(70) % e_conv(0:3) = (/ 12, 24, 20, 14 /)
  ptable(70) % eht_param(0:3) = (/ -5.35, -5.35, -5.21, &
          -13.86/)
! Lutetium
  ptable(71) % symbol = 'Lu'
  ptable(71) % name = 'Lutetium'
  ptable(71) % number = 71
  ptable(71) % amass = 174.97000
  ptable(71) % mass = 174.94090
  ptable(71) % covalent_radius = 1.5600
  ptable(71) % vdw_radius = z
  ptable(71) % e_conv(0:3) = (/ 12, 24, 21, 14 /)
  ptable(71) % eht_param(0:3) = (/ -6.05, -6.05, -5.12, &
          -22.40/)
! Hafnium
  ptable(72) % symbol = 'Hf'
  ptable(72) % name = 'Hafnium'
  ptable(72) % number = 72
  ptable(72) % amass = 178.49000
  ptable(72) % mass = 179.94680
  ptable(72) % covalent_radius = 1.4400
  ptable(72) % vdw_radius = z
  ptable(72) % e_conv(0:3) = (/ 12, 24, 22, 14 /)
  ptable(72) % eht_param(0:3) = (/ z, z, z, z/)
! Tantalum
  ptable(73) % symbol = 'Ta'
  ptable(73) % name = 'Tantalum'
  ptable(73) % number = 73
  ptable(73) % amass = 180.94800
  ptable(73) % mass = 180.94800
  ptable(73) % covalent_radius = 1.3400
  ptable(73) % vdw_radius = z
  ptable(73) % e_conv(0:3) = (/ 12, 24, 23, 14 /)
  ptable(73) % eht_param(0:3) = (/ -10.10, -6.86, -12.10, z/)
! Tungsten
  ptable(74) % symbol = 'W '
  ptable(74) % name = 'Tungsten'
  ptable(74) % number = 74
  ptable(74) % amass = 183.85000
  ptable(74) % mass = 183.95100
  ptable(74) % covalent_radius = 1.3000
  ptable(74) % vdw_radius = z
  ptable(74) % e_conv(0:3) = (/ 12, 24, 24, 14 /)
  ptable(74) % eht_param(0:3) = (/ -8.26, -5.17, -10.37, z/)
! Rhenium
  ptable(75) % symbol = 'Re'
  ptable(75) % name = 'Rhenium'
  ptable(75) % number = 75
  ptable(75) % amass = 186.20000
  ptable(75) % mass = 186.95600
  ptable(75) % covalent_radius = 1.2800
  ptable(75) % vdw_radius = z
  ptable(75) % e_conv(0:3) = (/ 12, 24, 25, 14 /)
  ptable(75) % eht_param(0:3) = (/ -9.36, -5.96, -12.66, z/)
! Osmium
  ptable(76) % symbol = 'Os'
  ptable(76) % name = 'Osmium'
  ptable(76) % number = 76
  ptable(76) % amass = 190.20000
  ptable(76) % mass = 192.00000
  ptable(76) % covalent_radius = 1.2600
  ptable(76) % vdw_radius = z
  ptable(76) % e_conv(0:3) = (/ 12, 24, 26, 14 /)
  ptable(76) % eht_param(0:3) = (/ -8.17, -4.81, -11.84, z/)
! Iridium
  ptable(77) % symbol = 'Ir'
  ptable(77) % name = 'Iridium'
  ptable(77) % number = 77
  ptable(77) % amass = 192.20000
  ptable(77) % mass = 192.96330
  ptable(77) % covalent_radius = 1.2700
  ptable(77) % vdw_radius = z
  ptable(77) % e_conv(0:3) = (/ 12, 24, 27, 14 /)
  ptable(77) % eht_param(0:3) = (/ -11.36, -4.50, -12.17, z/)
! Platinum
  ptable(78) % symbol = 'Pt'
  ptable(78) % name = 'Platinum'
  ptable(78) % number = 78
  ptable(78) % amass = 195.09000
  ptable(78) % mass = 194.96480
  ptable(78) % covalent_radius = 1.3000
  ptable(78) % vdw_radius = z
  ptable(78) % e_conv(0:3) = (/ 11, 24, 29, 14 /)
  ptable(78) % eht_param(0:3) = (/ -9.077, -5.475, -12.59, &
          z/)
! Gold
  ptable(79) % symbol = 'Au'
  ptable(79) % name = 'Gold'
  ptable(79) % number = 79
  ptable(79) % amass = 196.96700
  ptable(79) % mass = 196.96660
  ptable(79) % covalent_radius = 1.3400
  ptable(79) % vdw_radius = z
  ptable(79) % e_conv(0:3) = (/ 11, 24, 30, 14 /)
  ptable(79) % eht_param(0:3) = (/ -10.92, -5.55, -15.076, &
          z/)
! Mercury
  ptable(80) % symbol = 'Hg'
  ptable(80) % name = 'Mercury'
  ptable(80) % number = 80
  ptable(80) % amass = 200.59000
  ptable(80) % mass = 201.97060
  ptable(80) % covalent_radius = 1.4900
  ptable(80) % vdw_radius = z
  ptable(80) % e_conv(0:3) = (/ 12, 24, 30, 14 /)
  ptable(80) % eht_param(0:3) = (/ -13.68, -8.47, -17.50, z/)
! Thallium
  ptable(81) % symbol = 'Tl'
  ptable(81) % name = 'Thallium'
  ptable(81) % number = 81
  ptable(81) % amass = 204.37000
  ptable(81) % mass = 204.97450
  ptable(81) % covalent_radius = 1.4800
  ptable(81) % vdw_radius = z
  ptable(81) % e_conv(0:3) = (/ 12, 25, 30, 14 /)
  ptable(81) % eht_param(0:3) = (/ -11.60, -5.80, z, z/)
! Lead
  ptable(82) % symbol = 'Pb'
  ptable(82) % name = 'Lead'
  ptable(82) % number = 82
  ptable(82) % amass = 207.19000
  ptable(82) % mass = 207.97660
  ptable(82) % covalent_radius = 1.4700
  ptable(82) % vdw_radius = z
  ptable(82) % e_conv(0:3) = (/ 12, 26, 30, 14 /)
  ptable(82) % eht_param(0:3) = (/ -15.70, -8.00, z, z/)
! Bismuth
  ptable(83) % symbol = 'Bi'
  ptable(83) % name = 'Bismuth'
  ptable(83) % number = 83
  ptable(83) % amass = 208.98000
  ptable(83) % mass = 208.98040
  ptable(83) % covalent_radius = 1.4600
  ptable(83) % vdw_radius = z
  ptable(83) % e_conv(0:3) = (/ 12, 27, 30, 14 /)
  ptable(83) % eht_param(0:3) = (/ -15.19, -7.79, z, z/)
! Polonium
  ptable(84) % symbol = 'Po'
  ptable(84) % name = 'Polonium'
  ptable(84) % number = 84
  ptable(84) % amass = 209.98290
  ptable(84) % mass = 209.98290
  ptable(84) % covalent_radius = 1.4600
  ptable(84) % vdw_radius = z
  ptable(84) % e_conv(0:3) = (/ 12, 28, 30, 14 /)
  ptable(84) % eht_param(0:3) = (/ z, z, z, z/)
! Astatine
  ptable(85) % symbol = 'At'
  ptable(85) % name = 'Astatine'
  ptable(85) % number = 85
  ptable(85) % amass = 209.98700
  ptable(85) % mass = 209.98700
  ptable(85) % covalent_radius = 1.4500
  ptable(85) % vdw_radius = z
  ptable(85) % e_conv(0:3) = (/ 12, 29, 30, 14 /)
  ptable(85) % eht_param(0:3) = (/ z, z, z, z/)
! Radon
  ptable(86) % symbol = 'Rn'
  ptable(86) % name = 'Radon'
  ptable(86) % number = 86
  ptable(86) % amass = 222.01750
  ptable(86) % mass = 222.01750
  ptable(86) % covalent_radius = z
  ptable(86) % vdw_radius = z
  ptable(86) % e_conv(0:3) = (/ 12, 30, 30, 14 /)
  ptable(86) % eht_param(0:3) = (/ z, z, z, z/)
! Francium
  ptable(87) % symbol = 'Fr'
  ptable(87) % name = 'Francium'
  ptable(87) % number = 87
  ptable(87) % amass = 223.01980
  ptable(87) % mass = 223.01980
  ptable(87) % covalent_radius = z
  ptable(87) % vdw_radius = z
  ptable(87) % e_conv(0:3) = (/ 13, 30, 30, 14 /)
  ptable(87) % eht_param(0:3) = (/ z, z, z, z/)
! Radium
  ptable(88) % symbol = 'Ra'
  ptable(88) % name = 'Radium'
  ptable(88) % number = 88
  ptable(88) % amass = 226.02540
  ptable(88) % mass = 226.02540
  ptable(88) % covalent_radius = z
  ptable(88) % vdw_radius = z
  ptable(88) % e_conv(0:3) = (/ 14, 30, 30, 14 /)
  ptable(88) % eht_param(0:3) = (/ z, z, z, z/)
! Actinium
  ptable(89) % symbol = 'Ac'
  ptable(89) % name = 'Actinium'
  ptable(89) % number = 89
  ptable(89) % amass = 227.02780
  ptable(89) % mass = 227.02780
  ptable(89) % covalent_radius = z
  ptable(89) % vdw_radius = z
  ptable(89) % e_conv(0:3) = (/ 14, 30, 31, 14 /)
  ptable(89) % eht_param(0:3) = (/ z, z, z, z/)
! Thorium
  ptable(90) % symbol = 'Th'
  ptable(90) % name = 'Thorium'
  ptable(90) % number = 90
  ptable(90) % amass = 232.03810
  ptable(90) % mass = 232.03810
  ptable(90) % covalent_radius = 1.6500
  ptable(90) % vdw_radius = z
  ptable(90) % e_conv(0:3) = (/ 14, 30, 32, 14 /)
  ptable(90) % eht_param(0:3) = (/ -5.39, -5.39, -10.11, &
          -9.64/)
! Proctactinium
  ptable(91) % symbol = 'Pa'
  ptable(91) % name = 'Proctactinium'
  ptable(91) % number = 91
  ptable(91) % amass = 231.03590
  ptable(91) % mass = 231.03590
  ptable(91) % covalent_radius = z
  ptable(91) % vdw_radius = z
  ptable(91) % e_conv(0:3) = (/ 14, 30, 31, 16 /)
  ptable(91) % eht_param(0:3) = (/ z, z, z, z/)
! Uranium
  ptable(92) % symbol = 'U '
  ptable(92) % name = 'Uranium'
  ptable(92) % number = 92
  ptable(92) % amass = 238.05080
  ptable(92) % mass = 238.05080
  ptable(92) % covalent_radius = 1.4200
  ptable(92) % vdw_radius = z
  ptable(92) % e_conv(0:3) = (/ 14, 30, 31, 17 /)
  ptable(92) % eht_param(0:3) = (/ -5.50, -5.50, -9.19, &
          -10.62/)
! Neptunium
  ptable(93) % symbol = 'Np'
  ptable(93) % name = 'Neptunium'
  ptable(93) % number = 93
  ptable(93) % amass = 237.04820
  ptable(93) % mass = 237.04820
  ptable(93) % covalent_radius = z
  ptable(93) % vdw_radius = z
  ptable(93) % e_conv(0:3) = (/ 14, 30, 31, 18 /)
  ptable(93) % eht_param(0:3) = (/ z, z, z, z/)
! Plutonium
  ptable(94) % symbol = 'Pu'
  ptable(94) % name = 'Plutonium'
  ptable(94) % number = 94
  ptable(94) % amass = 244.0640
  ptable(94) % mass = 244.0640
  ptable(94) % covalent_radius = z
  ptable(94) % vdw_radius = z
  ptable(94) % e_conv(0:3) = (/ 14, 30, 30, 20 /)
  ptable(94) % eht_param(0:3) = (/ z, z, z, z/)
! Americum
  ptable(95) % symbol = 'Am'
  ptable(95) % name = 'Americum'
  ptable(95) % number = 95
  ptable(95) % amass = 243.0614
  ptable(95) % mass = 243.0614
  ptable(95) % covalent_radius = z
  ptable(95) % vdw_radius = z
  ptable(95) % e_conv(0:3) = (/ 14, 30, 30, 21 /)
  ptable(95) % eht_param(0:3) = (/ z, z, z, z/)
! Curium
  ptable(96) % symbol = 'Cm'
  ptable(96) % name = 'Curium'
  ptable(96) % number = 96
  ptable(96) % amass = 247.0700
  ptable(96) % mass = 247.0700
  ptable(96) % covalent_radius = z
  ptable(96) % vdw_radius = z
  ptable(96) % e_conv(0:3) = (/ 14, 30, 31, 21 /)
  ptable(96) % eht_param(0:3) = (/ z, z, z, z/)
! Berkelium
  ptable(97) % symbol = 'Bk'
  ptable(97) % name = 'Berkelium'
  ptable(97) % number = 97
  ptable(97) % amass = 251.0800
  ptable(97) % mass = 251.0800
  ptable(97) % covalent_radius = z
  ptable(97) % vdw_radius = z
  ptable(97) % e_conv(0:3) = (/ 14, 30, 31, 22 /)
  ptable(97) % eht_param(0:3) = (/ z, z, z, z/)
! Californium
  ptable(98) % symbol = 'Cf'
  ptable(98) % name = 'Californium'
  ptable(98) % number = 98
  ptable(98) % amass = 252.0820
  ptable(98) % mass = 252.0820
  ptable(98) % covalent_radius = z
  ptable(98) % vdw_radius = z
  ptable(98) % e_conv(0:3) = (/ 14, 30, 30, 24 /)
  ptable(98) % eht_param(0:3) = (/ z, z, z, z/)
! Einsteinium
  ptable(99) % symbol = 'Es'
  ptable(99) % name = 'Einsteinium'
  ptable(99) % number = 99
  ptable(99) % amass = 252.0829
  ptable(99) % mass = 252.0829
  ptable(99) % covalent_radius = z
  ptable(99) % vdw_radius = z
  ptable(99) % e_conv(0:3) = (/ 14, 30, 30, 25 /)
  ptable(99) % eht_param(0:3) = (/ z, z, z, z/)
! Fermium
  ptable(100) % symbol = 'Fm'
  ptable(100) % name = 'Fermium'
  ptable(100) % number = 100
  ptable(100) % amass = 257.0950
  ptable(100) % mass = 257.0950
  ptable(100) % covalent_radius = z
  ptable(100) % vdw_radius = z
  ptable(100) % e_conv(0:3) = (/ 14, 30, 30, 26 /)
  ptable(100) % eht_param(0:3) = (/ z, z, z, z/)
! Mendelevium
  ptable(101) % symbol = 'Md'
  ptable(101) % name = 'Mendelevium'
  ptable(101) % number = 101
  ptable(101) % amass = 256.0000
  ptable(101) % mass = 256.0000
  ptable(101) % covalent_radius = z
  ptable(101) % vdw_radius = z
  ptable(101) % e_conv(0:3) = (/ 14, 30, 30, 27 /)
  ptable(101) % eht_param(0:3) = (/ z, z, z, z/)
! Nobelium
  ptable(102) % symbol = 'No'
  ptable(102) % name = 'Nobelium'
  ptable(102) % number = 102
  ptable(102) % amass = 254.0000
  ptable(102) % mass = 254.0000
  ptable(102) % covalent_radius = z
  ptable(102) % vdw_radius = z
  ptable(102) % e_conv(0:3) = (/ 14, 30, 30, 28 /)
  ptable(102) % eht_param(0:3) = (/ z, z, z, z/)
! Lawrencium
  ptable(103) % symbol = 'Lr'
  ptable(103) % name = 'Lawrencium'
  ptable(103) % number = 103
  ptable(103) % amass = 257.0000
  ptable(103) % mass = 257.0000
  ptable(103) % covalent_radius = z
  ptable(103) % vdw_radius = z
  ptable(103) % e_conv(0:3) = (/ 14, 30, 31, 28 /)
  ptable(103) % eht_param(0:3) = (/ z, z, z, z/)
! Unnilquadium
  ptable(104) % symbol = 'Uq'
  ptable(104) % name = 'Unnilquadium'
  ptable(104) % number = 104
  ptable(104) % amass = 261.0000
  ptable(104) % mass = 261.0000
  ptable(104) % covalent_radius = z
  ptable(104) % vdw_radius = z
  ptable(104) % e_conv(0:3) = (/ 14, 30, 32, 28 /)
  ptable(104) % eht_param(0:3) = (/ z, z, z, z/)
! Unnilpentium
  ptable(105) % symbol = 'Up'
  ptable(105) % name = 'Unnilpentium'
  ptable(105) % number = 105
  ptable(105) % amass = 262.0000
  ptable(105) % mass = 262.0000
  ptable(105) % covalent_radius = z
  ptable(105) % vdw_radius = z
  ptable(105) % e_conv(0:3) = (/ 14, 30, 33, 28 /)
  ptable(105) % eht_param(0:3) = (/ z, z, z, z/)
! Unnilhexium
  ptable(106) % symbol = 'Uh'
  ptable(106) % name = 'Unnilhexium'
  ptable(106) % number = 106
  ptable(106) % amass = 263.0000
  ptable(106) % mass = 263.0000
  ptable(106) % covalent_radius = z
  ptable(106) % vdw_radius = z
  ptable(106) % e_conv(0:3) = (/ 14, 30, 34, 28 /)
  ptable(106) % eht_param(0:3) = (/ z, z, z, z/)
! Initialize heat of formation
  CALL init_eheat
! Initialize gyromagnetic ratio 
  CALL init_gratio
  table_is_initialized=1
END SUBROUTINE init_periodic_table
! *****************************************************************************
SUBROUTINE init_eheat 
! All values in kcal/mol
! Dummy
  ptable(0) % heat_of_formation = z
! Hydrogen
  ptable(1) % heat_of_formation = 52.102
! Helium
  ptable(2) % heat_of_formation = z
! Lithium
  ptable(3) % heat_of_formation = 38.410
! Beryllium
  ptable(4) % heat_of_formation = 76.960
! Boron
  ptable(5) % heat_of_formation = 135.700
! Carbon
  ptable(6) % heat_of_formation = 170.890
! Nitrogen
  ptable(7) % heat_of_formation = 113.000
! Oxygen
  ptable(8) % heat_of_formation = 59.559
! Fluorine
  ptable(9) % heat_of_formation = 18.890
! Neon
  ptable(10) % heat_of_formation = z
! Sodium
  ptable(11) % heat_of_formation = 25.850
! Magnesium
  ptable(12) % heat_of_formation = 35.000
! Aluminium
  ptable(13) % heat_of_formation = 79.490
! Silicon
  ptable(14) % heat_of_formation = 108.390
! Phosphorus
  ptable(15) % heat_of_formation = 75.570
! Sulfur
  ptable(16) % heat_of_formation = 66.400
! Chlorine
  ptable(17) % heat_of_formation = 28.990
! Argon
  ptable(18) % heat_of_formation = z
! Potassium
  ptable(19) % heat_of_formation = 21.420
! Calcium
  ptable(20) % heat_of_formation = 42.600
! Scandium
  ptable(21) % heat_of_formation = 90.300
! Titanium
  ptable(22) % heat_of_formation = 112.300
! Vanadium
  ptable(23) % heat_of_formation = 122.900
! Chromium
  ptable(24) % heat_of_formation = 95.000
! Manganese
  ptable(25) % heat_of_formation = 67.700
! Iron
  ptable(26) % heat_of_formation = 99.300
! Cobalt
  ptable(27) % heat_of_formation = 102.400
! Nickel
  ptable(28) % heat_of_formation = 102.800
! Copper
  ptable(29) % heat_of_formation = 80.700
! Zinc
  ptable(30) % heat_of_formation = 31.170
! Gallium
  ptable(31) % heat_of_formation = 65.400
! Germanium
  ptable(32) % heat_of_formation = 89.500
! Arsenic
  ptable(33) % heat_of_formation = 72.300
! Selenium
  ptable(34) % heat_of_formation = 54.300
! Bromine
  ptable(35) % heat_of_formation = 26.740
! Krypton
  ptable(36) % heat_of_formation = z
! Rubidium
  ptable(37) % heat_of_formation = 19.600
! Strontium
  ptable(38) % heat_of_formation = 39.100
! Yttrium
  ptable(39) % heat_of_formation = 101.500
! Zirconium
  ptable(40) % heat_of_formation = 145.500
! Niobium
  ptable(41) % heat_of_formation = 172.400
! Molybdenum
  ptable(42) % heat_of_formation = 157.300
! Technetium
  ptable(43) % heat_of_formation = z
! Ruthenium
  ptable(44) % heat_of_formation = 155.500
! Rhodium
  ptable(45) % heat_of_formation = 133.000
! Palladium
  ptable(46) % heat_of_formation = 90.000
! Silver
  ptable(47) % heat_of_formation = 68.100
! Cadmium
  ptable(48) % heat_of_formation = 26.720
! Indium
  ptable(49) % heat_of_formation = 58.000
! Tin
  ptable(50) % heat_of_formation = 72.200
! Antimony
  ptable(51) % heat_of_formation = 63.200
! Tellurium
  ptable(52) % heat_of_formation = 47.000
! Iodine
  ptable(53) % heat_of_formation = 25.517
! Xenon
  ptable(54) % heat_of_formation = z
! Cesium
  ptable(55) % heat_of_formation = 18.700
! Barium
  ptable(56) % heat_of_formation = 42.500
! Lantanum
  ptable(57) % heat_of_formation = z
! Cerium
  ptable(58) % heat_of_formation = 101.300
! Praseodymium
  ptable(59) % heat_of_formation = z
! Neodymium
  ptable(60) % heat_of_formation = z
! Promethium
  ptable(61) % heat_of_formation = z
! Samarium
  ptable(62) % heat_of_formation = 49.400
! Europium
  ptable(63) % heat_of_formation = z
! Gadolinium
  ptable(64) % heat_of_formation = z
! Terbium
  ptable(65) % heat_of_formation = z
! Dysprosium
  ptable(66) % heat_of_formation = z
! Holmium
  ptable(67) % heat_of_formation = z
! Erbium
  ptable(68) % heat_of_formation = 75.800
! Thulium
  ptable(69) % heat_of_formation = z
! Ytterbium
  ptable(70) % heat_of_formation = 36.350
! Lutetium
  ptable(71) % heat_of_formation = z
! Hafnium
  ptable(72) % heat_of_formation = 148.000
! Tantalum
  ptable(73) % heat_of_formation = 186.900
! Tungsten
  ptable(74) % heat_of_formation = 203.100
! Rhenium
  ptable(75) % heat_of_formation = 185.000
! Osmium
  ptable(76) % heat_of_formation = 188.000
! Iridium
  ptable(77) % heat_of_formation = 160.000
! Platinum
  ptable(78) % heat_of_formation = 135.200
! Gold
  ptable(79) % heat_of_formation = 88.000
! Mercury
  ptable(80) % heat_of_formation = 14.690
! Thallium
  ptable(81) % heat_of_formation = 43.550
! Lead
  ptable(82) % heat_of_formation = 46.620
! Bismuth
  ptable(83) % heat_of_formation = 50.100
! Polonium
  ptable(84) % heat_of_formation = z
! Astatine
  ptable(85) % heat_of_formation = z
! Radon
! from Radon no parametrisation in dynamo
  ptable(86) % heat_of_formation = z
! Francium
  ptable(87) % heat_of_formation = z
! Radium
  ptable(88) % heat_of_formation = z
! Actinium
  ptable(89) % heat_of_formation = z
! Thorium
  ptable(90) % heat_of_formation = z
! Proctactinium
  ptable(91) % heat_of_formation = z
! Uranium
  ptable(92) % heat_of_formation = z
! Neptunium
  ptable(93) % heat_of_formation = z
! Plutonium
  ptable(94) % heat_of_formation = z
! Americum
  ptable(95) % heat_of_formation = z
! Curium
  ptable(96) % heat_of_formation = z
! Berkelium
  ptable(97) % heat_of_formation = z
! Californium
  ptable(98) % heat_of_formation = z
! Einsteinium
  ptable(99) % heat_of_formation = z
! Fermium
  ptable(100) % heat_of_formation = z
! Mendelevium
  ptable(101) % heat_of_formation = z
! Nobelium
  ptable(102) % heat_of_formation = z
! Lawrencium
  ptable(103) % heat_of_formation = z
! Unnilquadium
  ptable(104) % heat_of_formation = z
! Unnilpentium
  ptable(105) % heat_of_formation = z
! Unnilhexium
  ptable(106) % heat_of_formation = z
END SUBROUTINE init_eheat 
! *****************************************************************************
SUBROUTINE init_gratio
! D. M. Granty and R. K. Harris, Encyclopedia of Nuclear 
! Magnetic Resonance vol. 5 (Wiley, Chichester UK, 1996)
!
! gyrom_ratio values in Mhz/Tesla
!
! Dummy
  ptable(0) % gyrom_ratio = 0.0 
  ptable(0) % gyrom_ratio_isotope = 0
! Hydrogen
  ptable(1) % gyrom_ratio = 42.5774690577
  ptable(1) % gyrom_ratio_isotope = 1
! Helium
  ptable(2) % gyrom_ratio = -32.4360299810
  ptable(2) % gyrom_ratio_isotope = 3
! Lithium
  ptable(3) % gyrom_ratio = 16.5484555869
  ptable(3) % gyrom_ratio_isotope = 7
! Beryllium
  ptable(4) % gyrom_ratio = -5.9836942827
  ptable(4) % gyrom_ratio_isotope = 9
! Boron
  ptable(5) % gyrom_ratio = 13.6629814024
  ptable(5) % gyrom_ratio_isotope = 11
! Carbon
  ptable(6) % gyrom_ratio = 10.7083965713
  ptable(6) % gyrom_ratio_isotope = 13
! Nitrogen
  ptable(7) % gyrom_ratio = 3.0777051853
  ptable(7) % gyrom_ratio_isotope = 14
! Oxygen
  ptable(8) % gyrom_ratio = -5.7742686593
  ptable(8) % gyrom_ratio_isotope = 17
! Fluorine
  ptable(9) % gyrom_ratio = 40.0775701637
  ptable(9) % gyrom_ratio_isotope = 19
! Neon
  ptable(10) % gyrom_ratio = -3.3630712715
  ptable(10) % gyrom_ratio_isotope = 21
! Sodium
  ptable(11) % gyrom_ratio = 11.2695216738
  ptable(11) % gyrom_ratio_isotope = 23
! Magnesium
  ptable(12) % gyrom_ratio = -2.6083426159
  ptable(12) % gyrom_ratio_isotope = 25
! Aluminium
  ptable(13) % gyrom_ratio = 11.1030809358
  ptable(13) % gyrom_ratio_isotope = 27
! Silicon
  ptable(14) % gyrom_ratio = -8.4654514231
  ptable(14) % gyrom_ratio_isotope = 29
! Phosphorus
  ptable(15) % gyrom_ratio = 17.2514409015
  ptable(15) % gyrom_ratio_isotope = 31
! Sulfur
  ptable(16) % gyrom_ratio = 3.2717242919
  ptable(16) % gyrom_ratio_isotope = 33
! Chlorine
  ptable(17) % gyrom_ratio = 4.1765408335
  ptable(17) % gyrom_ratio_isotope = 35
! Argon
  ptable(18) % gyrom_ratio = 0.0
  ptable(18) % gyrom_ratio_isotope = 0
! Potassium
  ptable(19) % gyrom_ratio = 1.9895335549
  ptable(19) % gyrom_ratio_isotope = 39
! Calcium
  ptable(20) % gyrom_ratio = -2.8696734409
  ptable(20) % gyrom_ratio_isotope = 43
! Scandium
  ptable(21) % gyrom_ratio = 10.3590726388
  ptable(21) % gyrom_ratio_isotope = 45
! Titanium
  ptable(22) % gyrom_ratio = -2.4040354154
  ptable(22) % gyrom_ratio_isotope = 47
! Vanadium
  ptable(23) % gyrom_ratio = 11.2132801367
  ptable(23) % gyrom_ratio_isotope = 51
! Chromium
  ptable(24) % gyrom_ratio = -2.4115156977
  ptable(24) % gyrom_ratio_isotope = 53
! Manganese
  ptable(25) % gyrom_ratio = 10.5762511769
  ptable(25) % gyrom_ratio_isotope = 55
! Iron
  ptable(26) % gyrom_ratio = 1.3815642187
  ptable(26) % gyrom_ratio_isotope = 57
! Cobalt
  ptable(27) % gyrom_ratio = 10.0776909966
  ptable(27) % gyrom_ratio_isotope = 59
! Nickel
  ptable(28) % gyrom_ratio = -3.8114425772
  ptable(28) % gyrom_ratio_isotope = 61
! Copper
  ptable(29) % gyrom_ratio = 11.3187637358
  ptable(29) % gyrom_ratio_isotope = 63
! Zinc
  ptable(30) % gyrom_ratio = 2.6685318322
  ptable(30) % gyrom_ratio_isotope = 67
! Gallium
  ptable(31) % gyrom_ratio = 10.2477560110
  ptable(31) % gyrom_ratio_isotope = 69
! Germanium
  ptable(32) % gyrom_ratio = -1.4897384913
  ptable(32) % gyrom_ratio_isotope = 73
! Arsenic
  ptable(33) % gyrom_ratio = 7.3150206071
  ptable(33) % gyrom_ratio_isotope = 75
! Selenium
  ptable(34) % gyrom_ratio = 8.1573046941
  ptable(34) % gyrom_ratio_isotope = 77
! Bromine
  ptable(35) % gyrom_ratio = 10.7041503174
  ptable(35) % gyrom_ratio_isotope = 79
! Krypton
  ptable(36) % gyrom_ratio = -1.6442297171
  ptable(36) % gyrom_ratio_isotope = 83
! Rubidium
  ptable(37) % gyrom_ratio = 4.1264181673
  ptable(37) % gyrom_ratio_isotope = 85
! Strontium
  ptable(38) % gyrom_ratio = -1.8524642249
  ptable(38) % gyrom_ratio_isotope = 87
! Yttrium
  ptable(39) % gyrom_ratio = -2.0949232525
  ptable(39) % gyrom_ratio_isotope = 89
! Zirconium
  ptable(40) % gyrom_ratio = -3.9747832953
  ptable(40) % gyrom_ratio_isotope = 91
! Niobium
  ptable(41) % gyrom_ratio = 10.4523417326
  ptable(41) % gyrom_ratio_isotope = 93
! Molybdenum
  ptable(42) % gyrom_ratio = -2.7868030535
  ptable(42) % gyrom_ratio_isotope = 95
! Technetium
  ptable(43) % gyrom_ratio = 9.6225078593
  ptable(43) % gyrom_ratio_isotope = 99
! Ruthenium
  ptable(44) % gyrom_ratio = -2.1915635664
  ptable(44) % gyrom_ratio_isotope = 101
! Rhodium
  ptable(45) % gyrom_ratio = -1.3477240581
  ptable(45) % gyrom_ratio_isotope = 103
! Palladium
  ptable(46) % gyrom_ratio = -1.9576058000
  ptable(46) % gyrom_ratio_isotope = 105
! Silver
  ptable(47) % gyrom_ratio = -1.7330669824
  ptable(47) % gyrom_ratio_isotope = 107
! Cadmium
  ptable(48) % gyrom_ratio = -9.0691469715
  ptable(48) % gyrom_ratio_isotope = 111
! Indium
  ptable(49) % gyrom_ratio = 9.3856853040
  ptable(49) % gyrom_ratio_isotope = 115
! Tin
  ptable(50) % gyrom_ratio = -15.9659464261
  ptable(50) % gyrom_ratio_isotope = 119
! Antimony
  ptable(51) % gyrom_ratio = 10.2551487581
  ptable(51) % gyrom_ratio_isotope = 121
! Tellurium
  ptable(52) % gyrom_ratio = -13.5454231953
  ptable(52) % gyrom_ratio_isotope = 125
! Iodine
  ptable(53) % gyrom_ratio = 8.5777718410
  ptable(53) % gyrom_ratio_isotope = 127
! Xenon
  ptable(54) % gyrom_ratio = -11.8603902888
  ptable(54) % gyrom_ratio_isotope = 129
! Cesium
  ptable(55) % gyrom_ratio = 5.6233482338
  ptable(55) % gyrom_ratio_isotope = 133
! Barium
  ptable(56) % gyrom_ratio = 4.7634278693
  ptable(56) % gyrom_ratio_isotope = 137
! Lantanum
  ptable(57) % gyrom_ratio = 6.0611483090
  ptable(57) % gyrom_ratio_isotope = 139
! Cerium
  ptable(58) % gyrom_ratio = 0.0
  ptable(58) % gyrom_ratio_isotope = 0
! Praseodymium
  ptable(59) % gyrom_ratio = 13.0359039238
  ptable(59) % gyrom_ratio_isotope = 141
! Neodymium
  ptable(60) % gyrom_ratio = -2.3188875208
  ptable(60) % gyrom_ratio_isotope = 143
! Promethium
  ptable(61) % gyrom_ratio = 5.7502680939
  ptable(61) % gyrom_ratio_isotope = 147
! Samarium
  ptable(62) % gyrom_ratio = -1.7745776155
  ptable(62) % gyrom_ratio_isotope = 147
! Europium
  ptable(63) % gyrom_ratio = 4.6742215237
  ptable(63) % gyrom_ratio_isotope = 153
! Gadolinium
  ptable(64) % gyrom_ratio = -1.7139395822
  ptable(64) % gyrom_ratio_isotope = 157
! Terbium
  ptable(65) % gyrom_ratio = 10.2352543902
  ptable(65) % gyrom_ratio_isotope = 159
! Dysprosium
  ptable(66) % gyrom_ratio = 2.0515072165
  ptable(66) % gyrom_ratio_isotope = 163
! Holmium
  ptable(67) % gyrom_ratio = 9.0877472505
  ptable(67) % gyrom_ratio_isotope = 165
! Erbium
  ptable(68) % gyrom_ratio = -1.2279917944
  ptable(68) % gyrom_ratio_isotope = 167
! Thulium
  ptable(69) % gyrom_ratio = -3.5300566378
  ptable(69) % gyrom_ratio_isotope = 169
! Ytterbium
  ptable(70) % gyrom_ratio = -2.0729931338
  ptable(70) % gyrom_ratio_isotope = 173
! Lutetium
  ptable(71) % gyrom_ratio = 4.8625018213
  ptable(71) % gyrom_ratio_isotope = 175
! Hafnium
  ptable(72) % gyrom_ratio = 1.7284226820
  ptable(72) % gyrom_ratio_isotope = 177
! Tantalum
  ptable(73) % gyrom_ratio = 5.1626680440
  ptable(73) % gyrom_ratio_isotope = 181
! Tungsten
  ptable(74) % gyrom_ratio = 1.7956502074
  ptable(74) % gyrom_ratio_isotope = 183
! Rhenium
  ptable(75) % gyrom_ratio = 9.8169951998
  ptable(75) % gyrom_ratio_isotope = 187
! Osmium
  ptable(76) % gyrom_ratio = 3.3536015524
  ptable(76) % gyrom_ratio_isotope = 189
! Iridium
  ptable(77) % gyrom_ratio = 0.8319028875
  ptable(77) % gyrom_ratio_isotope = 193
! Platinum
  ptable(78) % gyrom_ratio = 9.2922613524
  ptable(78) % gyrom_ratio_isotope = 195
! Gold
  ptable(79) % gyrom_ratio = 0.7528983738
  ptable(79) % gyrom_ratio_isotope = 197
! Mercury
  ptable(80) % gyrom_ratio = 7.7123168633
  ptable(80) % gyrom_ratio_isotope = 199
! Thallium
  ptable(81) % gyrom_ratio = 24.9748814221
  ptable(81) % gyrom_ratio_isotope = 205
! Lead
  ptable(82) % gyrom_ratio = 8.8815779373
  ptable(82) % gyrom_ratio_isotope = 207
! Bismuth
  ptable(83) % gyrom_ratio = 6.9630287603
  ptable(83) % gyrom_ratio_isotope = 209
! Polonium
  ptable(84) % gyrom_ratio = 11.7774657888
  ptable(84) % gyrom_ratio_isotope = 209
! Astatine
  ptable(85) % gyrom_ratio = 0.0
  ptable(85) % gyrom_ratio_isotope = 0
! Radon
  ptable(86) % gyrom_ratio = 0.0
  ptable(86) % gyrom_ratio_isotope = 0
! Francium
  ptable(87) % gyrom_ratio = 0.0
  ptable(87) % gyrom_ratio_isotope = 0
! Radium
  ptable(88) % gyrom_ratio = 0.0
  ptable(88) % gyrom_ratio_isotope = 0
! Actinium
  ptable(89) % gyrom_ratio = 5.5704230082
  ptable(89) % gyrom_ratio_isotope = 227
! Thorium
  ptable(90) % gyrom_ratio = 0.6366197724
  ptable(90) % gyrom_ratio_isotope = 229
! Proctactinium
  ptable(91) % gyrom_ratio = 5.1088736732
  ptable(91) % gyrom_ratio_isotope = 231
! Uranium
  ptable(92) % gyrom_ratio = -0.8276057041
  ptable(92) % gyrom_ratio_isotope = 235
! Neptunium
  ptable(93) % gyrom_ratio = 4.9338032358
  ptable(93) % gyrom_ratio_isotope = 237
! Plutonium
  ptable(94) % gyrom_ratio = 1.5469860469
  ptable(94) % gyrom_ratio_isotope = 239
! Americum
  ptable(95) % gyrom_ratio = 2.4509861236
  ptable(95) % gyrom_ratio_isotope = 243
! Curium
  ptable(96) % gyrom_ratio = 0.3183098862
  ptable(96) % gyrom_ratio_isotope = 247
! Berkelium
  ptable(97) % gyrom_ratio = 0.0
  ptable(97) % gyrom_ratio_isotope = 0
! Californium
  ptable(98) % gyrom_ratio = 0.0
  ptable(98) % gyrom_ratio_isotope = 0
! Einsteinium
  ptable(99) % gyrom_ratio = 0.0
  ptable(99) % gyrom_ratio_isotope = 0
! Fermium
  ptable(100) % gyrom_ratio = 0.0
  ptable(100) % gyrom_ratio_isotope = 0
! Mendelevium
  ptable(101) % gyrom_ratio = 0.0
  ptable(101) % gyrom_ratio_isotope = 0
! Nobelium
  ptable(102) % gyrom_ratio = 0.0
  ptable(102) % gyrom_ratio_isotope = 0
! Lawrencium
  ptable(103) % gyrom_ratio = 0.0
  ptable(103) % gyrom_ratio_isotope = 0
! Unnilquadium
  ptable(104) % gyrom_ratio = 0.0
  ptable(104) % gyrom_ratio_isotope = 0
! Unnilpentium
  ptable(105) % gyrom_ratio = 0.0
  ptable(105) % gyrom_ratio_isotope = 0
! Unnilhexium
  ptable(106) % gyrom_ratio = 0.0
  ptable(106) % gyrom_ratio_isotope = 0
END SUBROUTINE init_gratio
END MODULE periodic_table
