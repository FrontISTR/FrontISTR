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
!> \brief  This module contains functions for interpolation in 3 node 
!!  trianglar element (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/05/12
!>  \version    0.00
!======================================================================!

module shape_tri3n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_tri3n(areacoord,func)
      real(kind=kreal), intent(in) :: areacoord(2)
      real(kind=kreal) :: func(3)
      func(1:2) = areacoord(1:2)
      func(3)   = 1.d0-areacoord(1)-areacoord(2)
    end subroutine

    subroutine ShapeDeriv_tri3n(func)
      real(kind=kreal) :: func(3,2)
      func(1,1) = 1.d0
      func(2,1) = 0.d0
      func(3,1) = -1.d0

      func(1,2) = 0.d0
      func(2,2) = 1.d0
      func(3,2) = -1.d0
    end subroutine
    
    subroutine Shape2ndDeriv_tri3n(func)
      real(kind=kreal) :: func(3,2,2)
      func(:,:,:) = 0.d0
    end subroutine


      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE NodalNaturalCoord_tri3n(nncoord)
!####################################################################
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal), INTENT(OUT) :: nncoord(3, 2)
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      nncoord(1, 1) =  1.0D0
      nncoord(2, 1) =  0.0D0
      nncoord(3, 1) =  0.0D0
      ! eta-coordinate at a node in a local element
      nncoord(1, 2) =  0.0D0
      nncoord(2, 2) =  1.0D0
      nncoord(3, 2) =  0.0D0
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE NodalNaturalCoord_tri3n
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


end module
