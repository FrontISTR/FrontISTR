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
!> \brief  This module contains functions for interpolation in 8 node
!!      hexahedral element  (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/15
!>  \version    0.00
!======================================================================!

module shape_hex8n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_hex8n(localcoord,func)
      real(kind=kreal) :: localcoord(3)
      real(kind=kreal) :: func(8)
      func(1) = 0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(2))*(1.d0-localcoord(3))
      func(2) = 0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(2))*(1.d0-localcoord(3))
      func(3) = 0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(2))*(1.d0-localcoord(3))
      func(4) = 0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(2))*(1.d0-localcoord(3))
      func(5) = 0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(2))*(1.d0+localcoord(3))
      func(6) = 0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(2))*(1.d0+localcoord(3))
      func(7) = 0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(2))*(1.d0+localcoord(3))
      func(8) = 0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(2))*(1.d0+localcoord(3))
    end subroutine

    subroutine ShapeDeriv_hex8n(localcoord, func)
      real(kind=kreal) :: localcoord(3)
      real(kind=kreal) :: func(8,3)
      func(1,1) = -0.125d0*(1.d0-localcoord(2))*(1.d0-localcoord(3))
      func(2,1) =  0.125d0*(1.d0-localcoord(2))*(1.d0-localcoord(3))
      func(3,1) =  0.125d0*(1.d0+localcoord(2))*(1.d0-localcoord(3))
      func(4,1) = -0.125d0*(1.d0+localcoord(2))*(1.d0-localcoord(3))
      func(5,1) = -0.125d0*(1.d0-localcoord(2))*(1.d0+localcoord(3))
      func(6,1) =  0.125d0*(1.d0-localcoord(2))*(1.d0+localcoord(3))
      func(7,1) =  0.125d0*(1.d0+localcoord(2))*(1.d0+localcoord(3))
      func(8,1) = -0.125d0*(1.d0+localcoord(2))*(1.d0+localcoord(3))

      func(1,2) = -0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(3))
      func(2,2) = -0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(3))
      func(3,2) =  0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(3))
      func(4,2) =  0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(3))
      func(5,2) = -0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(3))
      func(6,2) = -0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(3))
      func(7,2) =  0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(3))
      func(8,2) =  0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(3))

      func(1,3) = -0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(2))
      func(2,3) = -0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(2))
      func(3,3) = -0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(2))
      func(4,3) = -0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(2))
      func(5,3) =  0.125d0*(1.d0-localcoord(1))*(1.d0-localcoord(2))
      func(6,3) =  0.125d0*(1.d0+localcoord(1))*(1.d0-localcoord(2))
      func(7,3) =  0.125d0*(1.d0+localcoord(1))*(1.d0+localcoord(2))
      func(8,3) =  0.125d0*(1.d0-localcoord(1))*(1.d0+localcoord(2))
    end subroutine

end module
