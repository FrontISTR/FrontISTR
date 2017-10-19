!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains functions for interpolation in 4 node
!!   qudrilateral element (Langrange  interpolation)
module shape_quad4n
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
  subroutine NodalNaturalCoord_quad4n(nncoord)
    !####################################################################

    implicit none

    !--------------------------------------------------------------------

    real(kind = kreal), intent(out) :: nncoord(4, 2)

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

    return

    !####################################################################
  end subroutine NodalNaturalCoord_quad4n
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


end module
