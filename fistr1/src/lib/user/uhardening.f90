!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!>  \brief   This function provides user interface to define hardening
!>  tangent 
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/01/12
!>  \version    0.00
!======================================================================!

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
