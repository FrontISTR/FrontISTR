!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains functions for interpolation in 3 node
!!  trianglar element (Langrange  interpolation)
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
  subroutine NodalNaturalCoord_tri3n(nncoord)
    !####################################################################

    implicit none

    !--------------------------------------------------------------------

    real(kind = kreal), intent(out) :: nncoord(3, 2)

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

    return

    !####################################################################
  end subroutine NodalNaturalCoord_tri3n
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


end module
