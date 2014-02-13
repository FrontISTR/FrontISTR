!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
!> \brief  This module contains functions for interpolation in 3 nodes 
!!     line element   (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/27
!>  \version    0.00
!======================================================================!

module shape_line3n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_line3n(lcoord,func)
      real(kind=kreal), intent(in) :: lcoord(1)
      real(kind=kreal) :: func(3)
      func(1) =-0.5d0*(1.d0-lcoord(1))*lcoord(1)
      func(2) = 0.5d0*(1.d0+lcoord(1))*lcoord(1)
      func(3) = 1.0d0-lcoord(1)*lcoord(1)
    end subroutine

    subroutine ShapeDeriv_line3n(lcoord,func)
      real(kind=kreal), intent(in) :: lcoord(1)
      real(kind=kreal) :: func(3,1)
      func(1,1) = lcoord(1)-0.5d0
      func(2,1) = lcoord(1)+0.5d0
      func(3,1) =-2.d0*lcoord(1)
    end subroutine

    subroutine Shape2ndDeriv_line3n(func)
      real(kind=kreal) :: func(3,1,1)
      func(1,1,1) = 1.d0
      func(2,1,1) = 1.d0
      func(3,1,1) = -2.d0
    end subroutine

end module
