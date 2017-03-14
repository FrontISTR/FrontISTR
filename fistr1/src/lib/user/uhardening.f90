!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This function provides user interface to define hardening
!>  tangent
!> This subroutine calculate isotropic hardening tangent
function uhardening( matl, pstrain )
    use hecmw
    implicit none
    real( KIND=kreal ), INTENT(IN) :: matl(:) !< material proerties
    real( KIND=kreal ), INTENT(IN) :: pstrain !< plastic strain
    real(kind=kreal) :: uhardening

    uhardening = 0.d0
end function

!> This subroutine calculate kinematic hardening tangent
function ukhardening( matl, pstrain )
    use hecmw
    implicit none
    real( KIND=kreal ), INTENT(IN) :: matl(:) !< material proerties
    real( KIND=kreal ), INTENT(IN) :: pstrain !< plastic strain
    real(kind=kreal) :: ukhardening

    ukhardening = 0.d0
end function

!> This subroutine calculate current yield value
function uCurrYield( matl, pstrain )
    use hecmw
    implicit none
    real( KIND=kreal ), INTENT(IN) :: matl(:) !< material proerties
    real( KIND=kreal ), INTENT(IN) :: pstrain !< plastic strain
    real(kind=kreal) :: uCurrYield

    uCurrYield = 0.d0
end function
