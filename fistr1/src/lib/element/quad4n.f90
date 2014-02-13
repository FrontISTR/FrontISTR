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
!!   qudrilateral element (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/05/08
!>  \version    0.00
!======================================================================!

MODULE shape_quad4n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_quad4n(lcoord,func)
      real(kind=kreal), intent(in) :: lcoord(2)
      real(kind=kreal) :: func(4)
      func(1) = 0.25d0*(1.d0-lcoord(1))*(1.d0-lcoord(2))
      func(2) = 0.25d0*(1.d0+lcoord(1))*(1.d0-lcoord(2))
      func(3) = 0.25d0*(1.d0+lcoord(1))*(1.d0+lcoord(2))
      func(4) = 0.25d0*(1.d0-lcoord(1))*(1.d0+lcoord(2))
    end subroutine

    subroutine ShapeDeriv_quad4n(lcoord,func)
      real(kind=kreal), intent(in) :: lcoord(2)
      real(kind=kreal) :: func(4,2)
      func(1,1) = -0.25d0*(1.d0-lcoord(2))
      func(2,1) =  0.25d0*(1.d0-lcoord(2))
      func(3,1) =  0.25d0*(1.d0+lcoord(2))
      func(4,1) = -0.25d0*(1.d0+lcoord(2))

      func(1,2) = -0.25d0*(1.d0-lcoord(1))
      func(2,2) = -0.25d0*(1.d0+lcoord(1))
      func(3,2) =  0.25d0*(1.d0+lcoord(1))
      func(4,2) =  0.25d0*(1.d0-lcoord(1))
    end subroutine
    
    subroutine Shape2ndDeriv_quad4n(func)
      real(kind=kreal) :: func(4,2,2)
      func(:,1,1) = 0.d0
      func(1,1,2) = 0.25d0
      func(2,1,2) = -0.25d0
      func(3,1,2) = 0.25d0
      func(4,1,2) = -0.25d0

      func(1,2,1) = 0.25d0
      func(2,2,1) = -0.25d0
      func(3,2,1) = 0.25d0
      func(4,2,1) = -0.25d0
      func(:,2,2) = 0.d0
    end subroutine


      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE NodalNaturalCoord_quad4n(nncoord)
!####################################################################
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal), INTENT(OUT) :: nncoord(4, 2)
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      nncoord(1, 1) = -1.0D0
      nncoord(2, 1) =  1.0D0
      nncoord(3, 1) =  1.0D0
      nncoord(4, 1) = -1.0D0
      ! eta-coordinate at a node in a local element
      nncoord(1, 2) = -1.0D0
      nncoord(2, 2) = -1.0D0
      nncoord(3, 2) =  1.0D0
      nncoord(4, 2) =  1.0D0
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE NodalNaturalCoord_quad4n
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


END MODULE
