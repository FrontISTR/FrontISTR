!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by X. YUAN                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!                                                                      !
!> \brief  This module contains functions for interpolation in 4 node 
!!    tetrahedron element (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/15
!>  \version    0.00
!======================================================================!

module shape_tet4n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_tet4n(volcoord,func)
      real(kind=kreal), intent(in) :: volcoord(3)
      real(kind=kreal) :: func(4)
      func(2:4) = volcoord(1:3)
      func(1)   = 1.d0-volcoord(1)-volcoord(2)-volcoord(3)
    end subroutine

    subroutine ShapeDeriv_tet4n(func)
      real(kind=kreal), intent(out) :: func(4,3)
      func(1,1) = -1.d0
      func(2,1) = 1.d0
      func(3,1) = 0.d0
      func(4,1) = 0.d0

      func(1,2) = -1.d0
      func(2,2) = 0.d0
      func(3,2) = 1.d0
      func(4,2) = 0.d0

      func(1,3) = -1.d0
      func(2,3) = 0.d0
      func(3,3) = 0.d0
      func(4,3) = 1.d0
    end subroutine

end module
