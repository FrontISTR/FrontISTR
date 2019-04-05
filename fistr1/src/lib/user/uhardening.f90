!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This function provides user interface to define hardening
!>  tangent
!> This subroutine calculate isotropic hardening tangent
function uhardening( matl, pstrain )
  use hecmw
  implicit none
  real( kind=kreal ), intent(in) :: matl(:) !< material proerties
  real( kind=kreal ), intent(in) :: pstrain !< plastic strain
  real(kind=kreal) :: uhardening

  uhardening = 0.d0
end function

!> This subroutine calculate kinematic hardening tangent
function ukhardening( matl, pstrain )
  use hecmw
  implicit none
  real( kind=kreal ), intent(in) :: matl(:) !< material proerties
  real( kind=kreal ), intent(in) :: pstrain !< plastic strain
  real(kind=kreal) :: ukhardening

  ukhardening = 0.d0
end function

!> This subroutine calculate current yield value
function uCurrYield( matl, pstrain )
  use hecmw
  implicit none
  real( kind=kreal ), intent(in) :: matl(:) !< material proerties
  real( kind=kreal ), intent(in) :: pstrain !< plastic strain
  real(kind=kreal) :: uCurrYield

  uCurrYield = 0.d0
end function
