!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains functions for interpolation in 9 node
!!   quadrilateral element
module shape_quad9n
  !####################################################################

  implicit none

  !--------------------------------------------------------------------

  integer, parameter, private :: kreal = kind( 0.0D0 )

  !--------------------------------------------------------------------

contains


  !####################################################################
  subroutine ShapeFunc_quad9n(lcoord, func)
    !####################################################################

    real(kind = kreal), intent(in)  :: lcoord(2)
    real(kind = kreal), intent(out) :: func(9)

    !--------------------------------------------------------------------

    real(kind = kreal) :: xi(9), eta(9)
    real(kind = kreal) :: n_xi(9), n_eta(9)

    integer :: na

    !--------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    xi(1)  = -1.0D0
    xi(2)  =  1.0D0
    xi(3)  =  1.0D0
    xi(4)  = -1.0D0
    xi(5)  =  0.0D0
    xi(6)  =  1.0D0
    xi(7)  =  0.0D0
    xi(8)  = -1.0D0
    xi(9)  =  0.0D0
    ! eta-coordinate at a node in a local element
    eta(1) = -1.0D0
    eta(2) = -1.0D0
    eta(3) =  1.0D0
    eta(4) =  1.0D0
    eta(5) = -1.0D0
    eta(6) =  0.0D0
    eta(7) =  1.0D0
    eta(8) =  0.0D0
    eta(9) =  0.0D0

    !--------------------------------------------------------------------

    do na = 1, 9

      n_xi(na)  = ( 0.5D0*xi(na) *lcoord(1) )    &
        *( 1.0D0+xi(na) *lcoord(1) )   &
        +( 1.0D0-xi(na) *xi(na)  )      &
        *( 1.0D0-lcoord(1)*lcoord(1) )
      n_eta(na) = ( 0.5D0*eta(na)*lcoord(2) )    &
        *( 1.0D0+eta(na)*lcoord(2) )   &
        +( 1.0D0-eta(na)*eta(na) )      &
        *( 1.0D0-lcoord(2)*lcoord(2) )

    end do

    !--------------------------------------------------------------------

    do na = 1, 9

      func(na) = n_xi(na)*n_eta(na)

    end do

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine ShapeFunc_quad9n
  !####################################################################

  !####################################################################
  subroutine ShapeDeriv_quad9n(lcoord, func)
    !####################################################################

    real(kind = kreal), intent(in)  :: lcoord(2)
    real(kind = kreal), intent(out) :: func(9, 2)

    !--------------------------------------------------------------------

    real(kind = kreal) :: xi(9), eta(9)
    real(kind = kreal) :: n_xi(9), n_eta(9)
    real(kind = kreal) :: dn_xi(9), dn_eta(9)

    integer :: na

    !--------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    xi(1)  = -1.0D0
    xi(2)  =  1.0D0
    xi(3)  =  1.0D0
    xi(4)  = -1.0D0
    xi(5)  =  0.0D0
    xi(6)  =  1.0D0
    xi(7)  =  0.0D0
    xi(8)  = -1.0D0
    xi(9)  =  0.0D0
    ! eta-coordinate at a node in a local element
    eta(1) = -1.0D0
    eta(2) = -1.0D0
    eta(3) =  1.0D0
    eta(4) =  1.0D0
    eta(5) = -1.0D0
    eta(6) =  0.0D0
    eta(7) =  1.0D0
    eta(8) =  0.0D0
    eta(9) =  0.0D0

    !--------------------------------------------------------------------

    do na = 1, 9

      n_xi(na)  = ( 0.5D0*xi(na) *lcoord(1) )    &
        *( 1.0D0+xi(na) *lcoord(1) )   &
        +( 1.0D0-xi(na) *xi(na)  )      &
        *( 1.0D0-lcoord(1)*lcoord(1) )
      n_eta(na) = ( 0.5D0*eta(na)*lcoord(2) )    &
        *( 1.0D0+eta(na)*lcoord(2) )   &
        +( 1.0D0-eta(na)*eta(na) )      &
        *( 1.0D0-lcoord(2)*lcoord(2) )

      dn_xi(na)  = ( 0.5D0*xi(na)  )            &
        *( 1.0D0+xi(na) *lcoord(1) ) &
        +( 0.5D0*xi(na) *lcoord(1) )  &
        *xi(na)                      &
        +( 1.0D0-xi(na) *xi(na)  )    &
        *( -2.0D0*lcoord(1) )
      dn_eta(na) = ( 0.5D0*eta(na) )            &
        *( 1.0D0+eta(na)*lcoord(2) ) &
        +( 0.5D0*eta(na)*lcoord(2) )  &
        *eta(na)                     &
        +( 1.0D0-eta(na)*eta(na) )    &
        *( -2.0D0*lcoord(2) )

    end do

    !--------------------------------------------------------------------

    do na = 1, 9

      func(na, 1) = dn_xi(na)*n_eta(na)
      func(na, 2) = n_xi(na) *dn_eta(na)

    end do

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine ShapeDeriv_quad9n
  !####################################################################

  !####################################################################
  subroutine NodalNaturalCoord_quad9n(nncoord)
    !####################################################################

    real(kind = kreal), intent(out) :: nncoord(9, 2)

    !--------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    nncoord(1, 1)  = -1.0D0
    nncoord(2, 1)  =  1.0D0
    nncoord(3, 1)  =  1.0D0
    nncoord(4, 1)  = -1.0D0
    nncoord(5, 1)  =  0.0D0
    nncoord(6, 1)  =  1.0D0
    nncoord(7, 1)  =  0.0D0
    nncoord(8, 1)  = -1.0D0
    nncoord(9, 1)  =  0.0D0
    ! eta-coordinate at a node in a local element
    nncoord(1, 2) = -1.0D0
    nncoord(2, 2) = -1.0D0
    nncoord(3, 2) =  1.0D0
    nncoord(4, 2) =  1.0D0
    nncoord(5, 2) = -1.0D0
    nncoord(6, 2) =  0.0D0
    nncoord(7, 2) =  1.0D0
    nncoord(8, 2) =  0.0D0
    nncoord(9, 2) =  0.0D0

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine NodalNaturalCoord_quad9n
  !####################################################################

end module shape_quad9n
