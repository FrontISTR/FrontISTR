!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides data for gauss quadrature
module gauss_integration
  use hecmw
  implicit none
  real(kind=kreal) :: XG(3, 3) !< abscissa of gauss points
  real(kind=kreal) :: WGT(3, 3) !< wieght of gauss points
  !****************************
  !* Gauss Integration Table **
  !****************************
  !** 1st ***
  data XG (1,1)/0.0/
  data WGT(1,1)/2.0/
  !** 2nd ***
  data XG(2,1),XG(2,2)/-0.577350269189626,0.577350269189626/
  data WGT(2,1),WGT(2,2)/1.0,1.0/
  !** 3rd ***
  data XG(3,1),XG(3,2),XG(3,3)/  &
    -0.7745966692,               &
    0.0,                         &
    0.7745966692/
  data WGT(3,1),WGT(3,2),WGT(3,3)/ &
    0.5555555555,                  &
    0.8888888888,                  &
    0.5555555555/
  ! end of this module
end module gauss_integration
