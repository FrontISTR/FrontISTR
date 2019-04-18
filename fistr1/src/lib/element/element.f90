!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module encapsulate the basic functions of all elements
!!      provide by this software.
!!
!>      After providing a specified elemental type, this module provide
!!      functions to fetch
!!-#      Space dimension
!!-#      Number of elements' nodes
!!-#      Suggested number of quadrature points
!!-#      Quadrature position and weight
!!-#      Shape function, first and second derivatives of shape
!!        functions in natural coordinate
!!-#      Elemental surfaces (Nodes' number)
!!
!>      Those would be supposed to adopt this module include
!!-       Those need to calculate shape function, shape derivatives
!!        in global or user-defined coordinate system.
!!-       Those need to do quadrature calculation.
!!-       Those need to fetch the information of sub-surfaces ( when
!!        calculate external surface load, e.g.)
!!-       And something else ...
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/16
!>  \version    0.00
!!
!!
!       If you wish introduce new elements with new geometry or/and
!    new shape functions, you need do the followings
!!
!!
!!-    Introduce new element ID corresponding to your element in this module.
!!-    Provide corresponding shape function and shape derivative in a new
!!    module and include this new module here.
!!-    Add the information of your new element into functions listed above.
!!-    If new quadrature method needed, do some modification in MODULE
!!    Quadrature.
!!-    If you introduce new surface geometry and do contact calculation also,
!!    you may need do some modification in MODULE mSurfElement.
module elementInfo

  use shape_line2n
  use shape_line3n
  use shape_tri3n
  use shape_tri6n
  use shape_quad4n
  use shape_quad8n
  use shape_quad9n
  use shape_hex8n
  use shape_hex20n
  use shape_tet4n
  use shape_tet10n
  use shape_prism6n
  use shape_prism15n
  implicit none

  integer, parameter, private :: kreal = kind(0.0d0)

  !-------------------------------------------
  !     Fllowing ID of element types
  !-------------------------------------------
  integer, parameter :: fe_unknown  = -1

  integer, parameter :: fe_line2n   = 111
  integer, parameter :: fe_line3n   = 112
  integer, parameter :: fe_tri3n    = 231
  integer, parameter :: fe_tri6n    = 232
  integer, parameter :: fe_tri6nc   = 2322
  integer, parameter :: fe_quad4n   = 241
  integer, parameter :: fe_quad8n   = 242
  integer, parameter :: fe_truss    = 301
  integer, parameter :: fe_tet4n    = 341
  integer, parameter :: fe_tet4n_pipi = 3414
  integer, parameter :: fe_tet10n   = 342
  integer, parameter :: fe_tet10nc  = 3422
  integer, parameter :: fe_prism6n  = 351
  integer, parameter :: fe_prism15n = 352
  integer, parameter :: fe_hex8n    = 361
  integer, parameter :: fe_hex20n   = 362
  integer, parameter :: fe_hex27n   = 363

  integer, parameter :: fe_beam2n   = 611
  integer, parameter :: fe_beam3n   = 612
  integer, parameter :: fe_beam341  = 641

  integer, parameter :: fe_tri6n_shell  = 732
  integer, parameter :: fe_dsg3_shell   = 733
  integer, parameter :: fe_mitc3_shell  = 731
  integer, parameter :: fe_mitc4_shell  = 741
  integer, parameter :: fe_mitc8_shell  = 742
  integer, parameter :: fe_mitc9_shell  = 743

  integer, parameter :: fe_mitc3_shell361  = 761
  integer, parameter :: fe_mitc4_shell361  = 781

  integer, parameter :: fe_tri3n_patch    = 1031
  integer, parameter :: fe_tri6n_patch    = 1032
  integer, parameter :: fe_quad4n_patch   = 1041
  integer, parameter :: fe_quad8n_patch   = 1042
  ! ---------------------------------------------

contains

  !************************************
  !    Following geometric information
  !************************************
  !> Obtain the space dimension of the element
  integer(kind=kind(2)) function getSpaceDimension( etype )
    integer, intent(in) :: etype    !< element type

    select case( etype)
      case (fe_line2n, fe_line3n)
        getSpaceDimension = 1
      case (fe_tri3n, fe_tri6n, fe_tri6nc, fe_quad4n, fe_quad8n)
        getSpaceDimension = 2
      case default
        getSpaceDimension = 3
    end select
  end function

  !> Obtain number of nodes of the element
  integer(kind=kind(2)) function getNumberOfNodes( etype )
    integer, intent(in) :: etype     !< element type

    select case (etype)
      case (fe_line2n, fe_beam2n)
        getNumberOfNodes = 2
      case (fe_line3n, fe_beam3n)
        getNumberOfNodes = 3
      case (fe_tri3n, fe_mitc3_shell, fe_mitc3_shell361, fe_tri3n_patch )
        getNumberOfNodes = 3
      case ( fe_tri6n, fe_tri6nc, fe_tri6n_shell, fe_tri6n_patch )
        getNumberOfNodes = 6
      case ( fe_quad4n, fe_mitc4_shell, fe_mitc4_shell361, fe_quad4n_patch )
        getNumberOfnodes = 4
      case ( fe_quad8n, fe_mitc8_shell, fe_quad8n_patch)
        getNumberOfNodes = 8
      case ( fe_mitc9_shell )
        getNumberOfNodes = 9
      case ( fe_tet4n, fe_tet4n_pipi, fe_beam341 )
        getNumberOfNodes = 4
      case ( fe_tet10n, fe_tet10nc )
        getNumberOfNodes = 10
      case ( fe_prism6n )
        getNumberOfNodes = 6
      case ( fe_prism15n )
        getNumberOfNodes = 15
      case ( fe_hex8n )
        getNumberOfNodes = 8
      case ( fe_hex20n )
        getNumberOfNodes = 20
      case default
        getNumberOfNodes = -1
        ! error message
    end select
  end function

  !> Obtain number of sub-surface
  integer(kind=kind(2)) function getNumberOfSubface( etype )
    integer, intent(in) :: etype    !< element type

    select case (etype)
      case (fe_line2n, fe_line3n)
        getNumberOfSubface = 2
      case (fe_tri3n, fe_mitc3_shell, fe_tri6n, fe_tri6nc, fe_tri6n_shell, fe_mitc3_shell361  )
        getNumberOfSubface = 3
      case ( fe_quad4n, fe_mitc4_shell, fe_quad8n, fe_mitc8_shell, fe_mitc9_shell, fe_mitc4_shell361  )
        getNumberOfSubface = 4
      case ( fe_tet4n, fe_tet4n_pipi, fe_tet10n, fe_tet10nc, fe_beam341 )
        getNumberOfSubface = 4
      case ( fe_prism6n, fe_prism15n )
        getNumberOfSubface = 5
      case ( fe_hex8n, fe_hex20n)
        getNumberOfSubface = 8
      case ( fe_tri3n_patch, fe_tri6n_patch, fe_quad4n_patch, fe_quad8n_patch )
        getNumberOfSubface = 1
      case default
        getNumberOfSubface = -1
        ! error message
    end select
  end function

  !> Find the definition of surface of the element
  subroutine getSubFace( intype, innumber, outtype, nodes )
    integer, intent(in)    :: intype   !< input element type, with space dimension n
    integer, intent(in)    :: innumber !< number of sub-surface
    integer, intent(out)   :: outtype  !< output element type, with space dimension n-1
    integer, intent(out)   :: nodes(:) !< output face nodes' ID

    if( innumber>getNumberOfSubface( intype ) ) stop "Error in getting subface"
    select case ( intype )
      case (fe_tet4n, fe_tet4n_pipi, fe_beam341)
        outtype = fe_tri3n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3
          case (2)
            nodes(1)=4; nodes(2)=2; nodes(3)=1
          case (3)
            nodes(1)=4; nodes(2)=3; nodes(3)=2
          case (4)
            nodes(1)=4; nodes(2)=1; nodes(3)=3
        end select
      case (fe_tet10n)
        outtype = fe_tri6n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3
            nodes(4)=5; nodes(5)=6; nodes(6)=7
          case (2)
            nodes(1)=4; nodes(2)=2; nodes(3)=1
            nodes(4)=9; nodes(5)=5; nodes(6)=8
          case (3)
            nodes(1)=4;  nodes(2)=3; nodes(3)=2
            nodes(4)=10; nodes(5)=6; nodes(6)=9
          case (4)
            nodes(1)=4; nodes(2)=1; nodes(3)=3
            nodes(4)=8; nodes(5)=7; nodes(6)=10
        end select
      case (fe_tet10nc)
        outtype = fe_tri6nc
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3
            nodes(4)=5; nodes(5)=6; nodes(6)=7
          case (2)
            nodes(1)=4; nodes(2)=2; nodes(3)=1
            nodes(4)=9; nodes(5)=5; nodes(6)=8
          case (3)
            nodes(1)=4;  nodes(2)=3; nodes(3)=2
            nodes(4)=10; nodes(5)=6; nodes(6)=9
          case (4)
            nodes(1)=4; nodes(2)=1; nodes(3)=3
            nodes(4)=8; nodes(5)=7; nodes(6)=10
        end select
      case ( fe_hex8n )
        outtype = fe_quad4n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3; nodes(4)=4
          case (2)
            nodes(1)=8; nodes(2)=7; nodes(3)=6; nodes(4)=5
          case (3)
            nodes(1)=5; nodes(2)=6; nodes(3)=2; nodes(4)=1
          case (4)
            nodes(1)=6; nodes(2)=7; nodes(3)=3; nodes(4)=2
          case (5)
            nodes(1)=7; nodes(2)=8; nodes(3)=4; nodes(4)=3
          case (6)
            nodes(1)=8; nodes(2)=5; nodes(3)=1; nodes(4)=4
          case default
            ! error
        end select
      case (fe_hex20n)
        outtype = fe_quad8n
        select case ( innumber )
          case (1)
            nodes(1)=1;  nodes(2)=2;  nodes(3)=3;  nodes(4)=4
            nodes(5)=9;  nodes(6)=10; nodes(7)=11; nodes(8)=12
          case (2)
            nodes(1)=8;  nodes(2)=7;  nodes(3)=6;  nodes(4)=5
            nodes(5)=15; nodes(6)=14; nodes(7)=13; nodes(8)=16
          case (3)
            nodes(1)=5;  nodes(2)=6;  nodes(3)=2;  nodes(4)=1
            nodes(5)=13; nodes(6)=18; nodes(7)=9;  nodes(8)=17
          case (4)
            nodes(1)=6;  nodes(2)=7;  nodes(3)=3;  nodes(4)=2
            nodes(5)=14; nodes(6)=19; nodes(7)=10; nodes(8)=18
          case (5)
            nodes(1)=7;  nodes(2)=8;  nodes(3)=4;  nodes(4)=3
            nodes(5)=15; nodes(6)=20; nodes(7)=11; nodes(8)=19
          case (6)
            nodes(1)=8;  nodes(2)=5;  nodes(3)=1;  nodes(4)=4
            nodes(5)=16; nodes(6)=17; nodes(7)=12; nodes(8)=20
          case default
            ! error
        end select
      case (fe_prism6n)
        select case ( innumber )
          case (1)
            outtype = fe_tri3n
            nodes(1)=1; nodes(2)=2; nodes(3)=3
          case (2)
            outtype = fe_tri3n
            nodes(1)=6; nodes(2)=5; nodes(3)=4
          case (3)
            outtype = fe_quad4n
            nodes(1)=4; nodes(2)=5; nodes(3)=2; nodes(4)=1
          case (4)
            outtype = fe_quad4n
            nodes(1)=5; nodes(2)=6; nodes(3)=3; nodes(4)=2
          case (5)
            outtype = fe_quad4n
            nodes(1)=6; nodes(2)=4; nodes(3)=1; nodes(4)=3
        end select
      case (fe_prism15n)
        select case ( innumber )
          case (1)
            outtype = fe_tri6n
            nodes(1)=1; nodes(2)=2; nodes(3)=3
            nodes(4)=7; nodes(5)=8; nodes(6)=9
          case (2)
            outtype = fe_tri6n
            nodes(1)=6;  nodes(2)=5;  nodes(3)=4
            nodes(4)=11; nodes(5)=10; nodes(6)=12
          case (3)
            outtype = fe_quad8n
            nodes(1)=4;  nodes(2)=5;  nodes(3)=2;  nodes(4)=1
            nodes(5)=10; nodes(6)=14; nodes(7)=7;  nodes(8)=13
          case (4)
            outtype = fe_quad8n
            nodes(1)=5;  nodes(2)=6;  nodes(3)=3;  nodes(4)=2
            nodes(5)=11; nodes(6)=15; nodes(7)=8;  nodes(8)=14
          case (5)
            outtype = fe_quad8n
            nodes(1)=6;  nodes(2)=4;  nodes(3)=1;  nodes(4)=3
            nodes(5)=12; nodes(6)=13; nodes(7)=9;  nodes(8)=15
        end select
      case ( fe_tri3n, fe_mitc3_shell )
        outtype = fe_line2n
        select case (innumber )
          case (1)
            nodes(1) = 1; nodes(2)=2
          case (2)
            nodes(1) = 2; nodes(2)=3
          case (3)
            nodes(1) = 3; nodes(2)=1
        end select
      case ( fe_tri6n, fe_tri6nc, fe_tri6n_shell )
        outtype = fe_line3n
        select case (innumber )
          case (1)
            nodes(1) = 1; nodes(2)=2;  nodes(3)=4
          case (2)
            nodes(1) = 2; nodes(2)=3;  nodes(3)=5
          case (3)
            nodes(1) = 3; nodes(2)=1;  nodes(3)=6
        end select
      case ( fe_quad4n, fe_mitc4_shell )
        outtype = fe_line2n
        select case (innumber )
          case (1)
            nodes(1) = 1; nodes(2)=2
          case (2)
            nodes(1) = 2; nodes(2)=3
          case (3)
            nodes(1) = 3; nodes(2)=4
          case (4)
            nodes(1) = 4; nodes(2)=1
        end select
      case ( fe_quad8n, fe_mitc8_shell, fe_mitc9_shell )
        outtype = fe_line3n
        select case (innumber )
          case (1)
            nodes(1) = 1; nodes(2)=2;  nodes(3)=5
          case (2)
            nodes(1) = 2; nodes(2)=3;  nodes(3)=6
          case (3)
            nodes(1) = 3; nodes(2)=4;  nodes(3)=7
          case (4)
            nodes(1) = 4; nodes(2)=1;  nodes(3)=8
        end select
      case (fe_mitc3_shell361)
        select case ( innumber )
          case (1)
            outtype = fe_tri3n
            nodes(1)=1; nodes(2)=2; nodes(3)=3
          case (2)
            outtype = fe_tri3n
            nodes(1)=6; nodes(2)=5; nodes(3)=4
          case (3)
            outtype = fe_quad4n
            nodes(1)=4; nodes(2)=5; nodes(3)=2; nodes(4)=1
          case (4)
            outtype = fe_quad4n
            nodes(1)=5; nodes(2)=6; nodes(3)=3; nodes(4)=2
          case (5)
            outtype = fe_quad4n
            nodes(1)=6; nodes(2)=4; nodes(3)=1; nodes(4)=3
        end select
      case ( fe_mitc4_shell361 )
        outtype = fe_quad4n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3; nodes(4)=4
          case (2)
            nodes(1)=8; nodes(2)=7; nodes(3)=6; nodes(4)=5
          case (3)
            nodes(1)=5; nodes(2)=6; nodes(3)=2; nodes(4)=1
          case (4)
            nodes(1)=6; nodes(2)=7; nodes(3)=3; nodes(4)=2
          case (5)
            nodes(1)=7; nodes(2)=8; nodes(3)=4; nodes(4)=3
          case (6)
            nodes(1)=8; nodes(2)=5; nodes(3)=1; nodes(4)=4
          case default
            ! error
        end select
      case ( fe_tri3n_patch )
        outtype = fe_tri3n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3
          case default
            !error
        end select
      case ( fe_tri6n_patch )
        outtype = fe_tri6n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3
            nodes(4)=4; nodes(5)=5; nodes(6)=6
          case default
            !error
        end select
      case ( fe_quad4n_patch )
        outtype = fe_quad4n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3; nodes(4)=4
          case default
            !error
        end select
      case ( fe_quad8n_patch )
        outtype = fe_quad8n
        select case ( innumber )
          case (1)
            nodes(1)=1; nodes(2)=2; nodes(3)=3; nodes(4)=4
            nodes(5)=5; nodes(6)=6; nodes(7)=7; nodes(8)=8
          case default
            !error
        end select
      case default
        outtype = fe_unknown
        stop "element type not defined-sbs"
        ! error message
    end select
  end subroutine

  !> Obtains the number of quadrature points of the element
  integer function NumOfQuadPoints( fetype )
    integer, intent(in) :: fetype         !< element type
    select case (fetype)
      case (fe_line2n, fe_tri3n, fe_tet4n, fe_beam2n , fe_beam341, fe_truss )
        NumOfQuadPoints = 1
      case ( fe_tri6n )
        NumOfQuadPoints = 3
      case ( fe_tri6nc )
        NumOfQuadPoints = 4
      case (fe_line3n )
        NumOfQuadPoints = 2
      case ( fe_quad4n, fe_mitc4_shell, fe_mitc4_shell361 )
        NumOfQuadPoints = 4
      case ( fe_quad8n, fe_mitc9_shell )
        NumOfQuadPoints = 9
      case ( fe_hex8n )
        NumOfQuadPoints = 8
      case ( fe_hex20n, fe_mitc8_shell )
        NumOfQuadPoints = 27
      case ( fe_prism6n )
        NumOfQuadPoints = 2
      case ( fe_mitc3_shell, fe_mitc3_shell361 )
        NumOfQuadPoints = 3
      case ( fe_prism15n, fe_tri6n_shell )
        NumOfQuadPoints = 9
      case ( fe_tet10n, fe_tet4n_pipi )
        NumOfQuadPoints = 4
      case ( fe_tet10nc )
        NumOfQuadPoints = 12
      case default
        NumOfQuadPoints = -1
        ! error message
        stop "element type not defined-np"
    end select
  end function

  !> Fetch the coordinate of gauss point
  subroutine getQuadPoint( fetype, np, pos )
    use Quadrature
    integer, intent(in)           :: fetype    !< element type
    integer, intent(in)           :: np        !< number of curr quadrature point
    real(kind=kreal), intent(out) :: pos(:)    !< natural coord of curr quadrature point

    if( np<1 .or. np>NumOfQuadPoints(fetype) ) then
      ! error
    endif

    select case (fetype)
      case (fe_tri3n)
        pos(1:2)=gauss2d4(:,np)
      case ( fe_tri6n, fe_mitc3_shell )
        pos(1:2)=gauss2d5(:,np)
      case (fe_tri6nc )
        pos(1:2)=gauss2d6(:,np)
      case ( fe_quad4n, fe_mitc4_shell )
        pos(1:2)=gauss2d2(:,np)
      case ( fe_quad8n, fe_mitc9_shell )
        pos(1:2)=gauss2d3(:,np)
      case ( fe_hex8n, fe_mitc4_shell361 )
        pos(1:3)=gauss3d2(:,np)
      case ( fe_hex20n, fe_mitc8_shell )
        pos(1:3)=gauss3d3(:,np)
      case ( fe_prism6n, fe_mitc3_shell361 )
        pos(1:3)=gauss3d7(:,np)
      case ( fe_prism15n, fe_tri6n_shell )
        pos(1:3)=gauss3d8(:,np)
      case ( fe_tet4n, fe_beam341 )
        pos(1:3)=gauss3d4(:,np)
      case ( fe_tet10n, fe_tet4n_pipi )
        pos(1:3)=gauss3d5(:,np)
      case ( fe_tet10nc )
        pos(1:3)=np
      case ( fe_line2n )
        pos(1:1)=gauss1d1(:,np)
      case ( fe_line3n )
        pos(1:1)=gauss1d2(:,np)
      case default
        ! error message
        stop "element type not defined-qp"
    end select
  end subroutine

  !> Fetch the weight value in given gauss point
  real(kind=kreal) function getWeight( fetype, np )
    use Quadrature
    integer, intent(in)           :: fetype       !< element type
    integer, intent(in)           :: np           !< number of curr quadrature point
    if( np<1 .or. np>NumOfQuadPoints(fetype) ) then
      ! error
    endif

    select case (fetype)
      case (fe_tri3n)
        getWeight = weight2d4(1)
      case ( fe_tri6n, fe_mitc3_shell )
        getWeight = weight2d5(np)
      case ( fe_quad4n, fe_mitc4_shell )
        getWeight = weight2d2(np)
      case ( fe_quad8n, fe_mitc9_shell )
        getWeight = weight2d3(np)
      case ( fe_hex8n, fe_mitc4_shell361 )
        getWeight = weight3d2(np)
      case ( fe_hex20n)
        getWeight = weight3d3(np)
      case ( fe_prism6n, fe_mitc3_shell361 )
        getWeight = weight3d7(np)
      case ( fe_prism15n )
        getWeight = weight3d8(np)
      case ( fe_tet4n, fe_beam341 )
        getWeight = weight3d4(1)
      case ( fe_tet10n, fe_tet4n_pipi )
        getWeight = weight3d5(np)
      case ( fe_line2n )
        getWeight = weight1d1(1)
      case ( fe_line3n )
        getWeight = weight1d2(np)
      case default
        getWeight = 0.d0
        ! error message
    end select
  end function

  !************************************
  !    Following shape function information
  !************************************
  !> Calculate deivatives of shape fucntion in natural coordiante system
  subroutine getShapeDeriv( fetype, localcoord, shapederiv )
    integer, intent(in)           :: fetype             !< input element type
    real(kind=kreal), intent(in)  :: localcoord(:)      !< natural points
    real(kind=kreal), intent(out) :: shapederiv(:,:)    !< deivative of shape function

    select case (fetype)
      case ( fe_tri3n, fe_mitc3_shell )
        !error check
        call ShapeDeriv_tri3n(shapederiv(1:3,1:2))
      case (fe_tri6n)
        !error check
        call ShapeDeriv_tri6n(localcoord,shapederiv(1:6,1:2) )
      case ( fe_quad4n, fe_mitc4_shell )
        !error check
        call ShapeDeriv_quad4n(localcoord,shapederiv(1:4,1:2))
      case (fe_quad8n)
        !error check
        call ShapeDeriv_quad8n(localcoord,shapederiv(1:8,1:2))
      case ( fe_mitc9_shell )
        !error check
        call ShapeDeriv_quad9n(localcoord,shapederiv(1:9,1:2))
      case (fe_hex8n, fe_mitc4_shell361)
        ! error check
        call ShapeDeriv_hex8n(localcoord,shapederiv(1:8,1:3))
      case (fe_hex20n)
        ! error check
        call ShapeDeriv_hex20n(localcoord, shapederiv(1:20,1:3))
      case (fe_prism6n, fe_mitc3_shell361)
        call ShapeDeriv_prism6n(localcoord,shapederiv(1:6,1:3))
      case (fe_prism15n)
        call ShapeDeriv_prism15n(localcoord,shapederiv(1:15,1:3))
      case (fe_tet4n, fe_tet4n_pipi, fe_beam341)
        ! error check
        call ShapeDeriv_tet4n(shapederiv(1:4,1:3))
      case (fe_tet10n)
        ! error check
        call ShapeDeriv_tet10n(localcoord,shapederiv(1:10,1:3))
      case default
        ! error message
        stop "Element type not defined-sde"
    end select
  end subroutine

  !> Calculate the 2nd derivative of shape function in natural coodinate system
  subroutine getShape2ndDeriv( fetype, localcoord, shapederiv )
    integer, intent(in)           :: fetype             !< elemental type
    real(kind=kreal), intent(in)  :: localcoord(:)      !< natural points
    real(kind=kreal), intent(out) :: shapederiv(:,:,:)  !< 2nd order shape derivatives

    select case (fetype)
      case ( fe_tri3n, fe_mitc3_shell )
        !error check
        call Shape2ndDeriv_tri3n(shapederiv(1:3,1:2,1:2))
      case (fe_tri6n)
        !error check
        call Shape2ndDeriv_tri6n(shapederiv(1:6,1:2,1:2))
      case ( fe_quad4n, fe_mitc4_shell )
        !error check
        call Shape2ndDeriv_quad4n(shapederiv(1:4,1:2,1:2))
      case (fe_quad8n)
        !error check
        call Shape2ndDeriv_quad8n(localcoord,shapederiv(1:8,1:2,1:2))
      case default
        ! error message
        stop "Cannot calculate second derivatives of shape function"
    end select
  end subroutine

  !> Calculate the shape function in natural coodinate system
  subroutine getShapeFunc( fetype, localcoord, func )
    integer, intent(in)           :: fetype            !< input element type
    real(kind=kreal), intent(in)  :: localcoord(:)     !< natural points
    real(kind=kreal), intent(out) :: func(:)           !< shape function

    select case (fetype)
      case ( fe_tri3n, fe_mitc3_shell )
        !error check
        call ShapeFunc_tri3n(localcoord,func(1:3))
      case (fe_tri6n)
        !error check
        call ShapeFunc_tri6n(localcoord,func(1:6))
      case ( fe_quad4n, fe_mitc4_shell )
        !error check
        call ShapeFunc_quad4n(localcoord,func(1:4))
      case (fe_quad8n)
        !error check
        call ShapeFunc_quad8n(localcoord,func(1:8))
      case (fe_hex8n, fe_mitc4_shell361)
        ! error check
        call ShapeFunc_hex8n(localcoord,func(1:8))
      case ( fe_mitc9_shell )
        !error check
        call ShapeFunc_quad9n(localcoord,func(1:9))
      case (fe_hex20n)
        ! error check
        call ShapeFunc_hex20n(localcoord,func(1:20))
      case (fe_prism6n, fe_mitc3_shell361)
        call ShapeFunc_prism6n(localcoord,func(1:6))
      case (fe_prism15n)
        call ShapeFunc_prism15n(localcoord,func(1:15))
      case (fe_tet4n, fe_tet4n_pipi, fe_beam341)
        ! error check
        call ShapeFunc_tet4n(localcoord,func(1:4))
      case (fe_tet10n)
        ! error check
        call ShapeFunc_tet10n(localcoord,func(1:10))
      case (fe_line2n)
        !error check
        call ShapeFunc_line2n(localcoord,func(1:2))
      case (fe_line3n)
        !error check
        call ShapeFunc_line3n(localcoord,func(1:3))
      case default
        stop "Element type not defined-sf"
        ! error message
    end select
  end subroutine


  ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
  !####################################################################
  subroutine getNodalNaturalCoord(fetype, nncoord)
    !####################################################################

    integer, intent(in)             :: fetype
    real(kind = kreal), intent(out) :: nncoord(:, :)

    !--------------------------------------------------------------------

    select case( fetype )
      case( fe_tri3n, fe_mitc3_shell, fe_mitc3_shell361 )

        !error check
        call NodalNaturalCoord_tri3n( nncoord(1:3, 1:2) )

      case( fe_quad4n, fe_mitc4_shell, fe_mitc4_shell361 )

        !error check
        call NodalNaturalCoord_quad4n( nncoord(1:4, 1:2) )

      case( fe_mitc9_shell )

        !error check
        call NodalNaturalCoord_quad9n( nncoord(1:9, 1:2) )

      case default

        ! error message
        stop "Element type not defined-sde"

    end select

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine getNodalNaturalCoord
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


  !> Calculate shape derivative in global coordinate system
  subroutine getGlobalDeriv( fetype, nn, localcoord, elecoord, det, gderiv )
    integer, intent(in)           :: fetype          !< element type
    integer, intent(in)           :: nn              !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(:)   !< curr position with natural coord
    real(kind=kreal), intent(in)  :: elecoord(:,:)   !< nodal coord of curr element
    real(kind=kreal), intent(out) :: det             !< nodal coord of curr element
    real(kind=kreal), intent(out) :: gderiv(:,:)     !< shape deivative in global coordinate system

    real(kind=kreal) :: DUM, XJ(3,3), XJI(3,3), deriv(nn,3)
    integer          :: nspace, i

    nspace = getSpaceDimension( fetype )
    call getShapeDeriv( fetype, localCoord(:), deriv(1:nn,:) )

    if( nspace==2 ) then
      XJ(1:2,1:2)=matmul( elecoord(1:2,1:nn), deriv(1:nn,1:2) )
      DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)
      if( det==0.d0 ) stop "Math error in GetGlobalDeriv! Determinant==0.0"
      DUM=1.d0/DET
      XJI(1,1)= XJ(2,2)*DUM
      XJI(1,2)=-XJ(1,2)*DUM
      XJI(2,1)=-XJ(2,1)*DUM
      XJI(2,2)= XJ(1,1)*DUM
    else
      !  JACOBI MATRIX
      XJ(1:3,1:3)= matmul( elecoord(1:3,1:nn), deriv(1:nn,1:3) )
      !DETERMINANT OF JACOBIAN
      DET=XJ(1,1)*XJ(2,2)*XJ(3,3)                                             &
        +XJ(2,1)*XJ(3,2)*XJ(1,3)                                             &
        +XJ(3,1)*XJ(1,2)*XJ(2,3)                                             &
        -XJ(3,1)*XJ(2,2)*XJ(1,3)                                             &
        -XJ(2,1)*XJ(1,2)*XJ(3,3)                                             &
        -XJ(1,1)*XJ(3,2)*XJ(2,3)
      if( det==0.d0 ) stop "Math error in GetGlobalDeriv! Determinant==0.0"
      ! INVERSION OF JACOBIAN
      DUM=1.d0/DET
      XJI(1,1)=DUM*( XJ(2,2)*XJ(3,3)-XJ(3,2)*XJ(2,3) )
      XJI(1,2)=DUM*(-XJ(1,2)*XJ(3,3)+XJ(3,2)*XJ(1,3) )
      XJI(1,3)=DUM*( XJ(1,2)*XJ(2,3)-XJ(2,2)*XJ(1,3) )
      XJI(2,1)=DUM*(-XJ(2,1)*XJ(3,3)+XJ(3,1)*XJ(2,3) )
      XJI(2,2)=DUM*( XJ(1,1)*XJ(3,3)-XJ(3,1)*XJ(1,3) )
      XJI(2,3)=DUM*(-XJ(1,1)*XJ(2,3)+XJ(2,1)*XJ(1,3) )
      XJI(3,1)=DUM*( XJ(2,1)*XJ(3,2)-XJ(3,1)*XJ(2,2) )
      XJI(3,2)=DUM*(-XJ(1,1)*XJ(3,2)+XJ(3,1)*XJ(1,2) )
      XJI(3,3)=DUM*( XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2) )
    endif

    gderiv(1:nn,1:nspace)=matmul( deriv(1:nn,1:nspace), XJI(1:nspace,1:nspace) )
  end subroutine

  !> Calculate shape derivative in global coordinate system
  real(kind=kreal) function getDeterminant( fetype, nn, localcoord, elecoord )
    integer, intent(in)           :: fetype          !< element type
    integer, intent(in)           :: nn              !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(:)   !< curr position with natural coord
    real(kind=kreal), intent(in)  :: elecoord(:,:)   !< nodal coord of curr element

    real(kind=kreal) :: XJ(3,3), deriv(nn,3)
    integer          :: nspace

    nspace = getSpaceDimension( fetype )
    call getShapeDeriv( fetype, localCoord(:), deriv(1:nn,:) )

    if( nspace==2 ) then
      XJ(1:2,1:2)=matmul( elecoord(1:2,1:nn), deriv(1:nn,1:2) )
      getDeterminant=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)
    else
      XJ(1:3,1:3)= matmul( elecoord(1:3,1:nn), deriv(1:nn,1:3) )
      getDeterminant=XJ(1,1)*XJ(2,2)*XJ(3,3)                                             &
        +XJ(2,1)*XJ(3,2)*XJ(1,3)                                             &
        +XJ(3,1)*XJ(1,2)*XJ(2,3)                                             &
        -XJ(3,1)*XJ(2,2)*XJ(1,3)                                             &
        -XJ(2,1)*XJ(1,2)*XJ(3,3)                                             &
        -XJ(1,1)*XJ(3,2)*XJ(2,3)
    endif

  end function

  !> calculate Jacobian matrix, its determinant and inverse
  subroutine getJacobian( fetype, nn, localcoord, elecoord, det, jacobian, inverse )
    integer, intent(in)           :: fetype          !< element type
    integer, intent(in)           :: nn              !< number of element nodes
    real(kind=kreal), intent(in)  :: localcoord(:)   !< curr position with natural coord
    real(kind=kreal), intent(in)  :: elecoord(:,:)   !< nodal coord of curr element
    real(kind=kreal), intent(out) :: det             !< nodal coord of curr element
    real(kind=kreal), intent(out) :: jacobian(:,:)   !< jacobian
    real(kind=kreal), intent(out) :: inverse(:,:)    !< inverse of jacobian

    real(kind=kreal) :: dum, deriv(NN,3)
    integer          :: nspace

    nspace = getSpaceDimension( fetype )
    call getShapeDeriv( fetype, localCoord(:), deriv(1:NN,:) )

    if( nspace==2 ) then
      jacobian(1:2,1:2)=matmul( elecoord(1:2,1:NN), deriv(1:NN,1:2) )
      det=jacobian(1,1)*jacobian(2,2)-jacobian(2,1)*jacobian(1,2)
      if( det==0.d0 ) stop "Math error in getJacobain! Determinant==0.0"
      dum=1.0/det
      inverse(1,1)= jacobian(2,2)*dum
      inverse(1,2)=-jacobian(1,2)*dum
      inverse(2,1)=-jacobian(2,1)*dum
      inverse(2,2)= jacobian(1,1)*dum
    else
      !  JACOBI MATRIX
      jacobian(1:3,1:3)= matmul( elecoord(1:3,1:NN), deriv(1:NN,1:3) )
      !DETERMINANT OF JACOBIAN
      det=jacobian(1,1)*jacobian(2,2)*jacobian(3,3)                                         &
        +jacobian(2,1)*jacobian(3,2)*jacobian(1,3)                                         &
        +jacobian(3,1)*jacobian(1,2)*jacobian(2,3)                                         &
        -jacobian(3,1)*jacobian(2,2)*jacobian(1,3)                                         &
        -jacobian(2,1)*jacobian(1,2)*jacobian(3,3)                                         &
        -jacobian(1,1)*jacobian(3,2)*jacobian(2,3)
      if( det==0.d0 ) stop "Math error in getJacobain! Determinant==0.0"
      ! INVERSION OF JACOBIAN
      dum=1.d0/det
      inverse(1,1)=DUM*( jacobian(2,2)*jacobian(3,3)-jacobian(3,2)*jacobian(2,3) )
      inverse(1,2)=DUM*(-jacobian(1,2)*jacobian(3,3)+jacobian(3,2)*jacobian(1,3) )
      inverse(1,3)=DUM*( jacobian(1,2)*jacobian(2,3)-jacobian(2,2)*jacobian(1,3) )
      inverse(2,1)=DUM*(-jacobian(2,1)*jacobian(3,3)+jacobian(3,1)*jacobian(2,3) )
      inverse(2,2)=DUM*( jacobian(1,1)*jacobian(3,3)-jacobian(3,1)*jacobian(1,3) )
      inverse(2,3)=DUM*(-jacobian(1,1)*jacobian(2,3)+jacobian(2,1)*jacobian(1,3) )
      inverse(3,1)=DUM*( jacobian(2,1)*jacobian(3,2)-jacobian(3,1)*jacobian(2,2) )
      inverse(3,2)=DUM*(-jacobian(1,1)*jacobian(3,2)+jacobian(3,1)*jacobian(1,2) )
      inverse(3,3)=DUM*( jacobian(1,1)*jacobian(2,2)-jacobian(2,1)*jacobian(1,2) )
    endif
  end subroutine

  !> Calculate  normal of 3d-surface
  function SurfaceNormal( fetype, nn, localcoord, elecoord ) result( normal )
    integer, intent(in)           :: fetype            !< type of surface element
    integer, intent(in)           :: nn                !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(2)     !< position
    real(kind=kreal), intent(in)  :: elecoord(3,nn)    !< nodes coordinates of element
    real(kind=kreal) :: normal(3)
    real(kind=kreal) :: deriv(nn,2), gderiv(3,2)

    select case (fetype)
      case (fe_tri3n)
        !error check
        call ShapeDeriv_tri3n(deriv(1:3,1:2))
      case (fe_tri6n)
        !error check
        call ShapeDeriv_tri6n(localcoord,deriv(1:6,1:2))
      case (fe_quad4n)
        !error check
        call ShapeDeriv_quad4n(localcoord,deriv(1:4,1:2))
      case (fe_quad8n)
        !error check
        call ShapeDeriv_quad8n(localcoord,deriv(1:8,1:2))
      case default
        ! error message
        normal =0.d0
        return
    end select

    gderiv = matmul( elecoord, deriv )
    normal(1) = gderiv(2,1)*gderiv(3,2) - gderiv(3,1)*gderiv(2,2)
    normal(2) = gderiv(3,1)*gderiv(1,2) - gderiv(1,1)*gderiv(3,2)
    normal(3) = gderiv(1,1)*gderiv(2,2) - gderiv(2,1)*gderiv(1,2)
    !  normal = normal/dsqrt(dot_product(normal, normal))
  end function

  !> Calculate normal of 2d-edge
  function EdgeNormal( fetype, nn, localcoord, elecoord ) result( normal )
    integer, intent(in)           :: fetype            !< type of surface element
    integer, intent(in)           :: nn                !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(1)     !< position
    real(kind=kreal), intent(in)  :: elecoord(2,nn)    !< nodes coordinates of element
    real(kind=kreal) :: normal(2)
    real(kind=kreal) :: deriv(nn,1), gderiv(2,1)

    select case (fetype)
      case (fe_line2n)
        !error check
        call ShapeDeriv_line2n(deriv(1:nn,:))
      case (fe_line3n)
        !error check
        call ShapeDeriv_line3n(localcoord,deriv(1:nn,:))
      case default
        ! error message
        normal =0.d0
        return
    end select

    gderiv = matmul( elecoord, deriv )
    normal(1) = -gderiv(2,1)
    normal(2) = gderiv(1,1)
    !  normal = normal/dsqrt(dot_product(normal, normal))
  end function

  !> Calculate base vector of tangent space of 3d surface
  subroutine TangentBase( fetype, nn, localcoord, elecoord, tangent )
    integer, intent(in)           :: fetype            !< type of surface element
    integer, intent(in)           :: nn                !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(2)     !< position
    real(kind=kreal), intent(in)  :: elecoord(3,nn)    !< nodes coordinates of element
    real(kind=kreal), intent(out) :: tangent(3,2)      !< two tangent vectors
    real(kind=kreal) :: deriv(nn,2)

    select case (fetype)
      case (fe_tri3n)
        !error check
        call ShapeDeriv_tri3n(deriv(1:3,1:2))
      case (fe_tri6n)
        !error check
        call ShapeDeriv_tri6n(localcoord,deriv(1:6,1:2))
      case (fe_tri6nc)
        !error check
        call ShapeDeriv_tri6n(localcoord,deriv(1:6,1:2))
      case (fe_quad4n)
        !error check
        call ShapeDeriv_quad4n(localcoord,deriv(1:4,1:2))
      case (fe_quad8n)
        !error check
        call ShapeDeriv_quad8n(localcoord,deriv(1:8,1:2))
      case default
        ! error message
        tangent =0.d0
        return
    end select

    tangent = matmul( elecoord, deriv )
  end subroutine TangentBase

  !> Calculate curvature tensor at a point along 3d surface
  subroutine Curvature( fetype, nn, localcoord, elecoord, l2ndderiv, normal, curv )
    integer, intent(in)           :: fetype            !< type of surface element
    integer, intent(in)           :: nn                !< number of elemental nodes
    real(kind=kreal), intent(in)  :: localcoord(2)     !< position
    real(kind=kreal), intent(in)  :: elecoord(3,nn)    !< nodes coordinates of element
    real(kind=kreal), intent(out) :: l2ndderiv(3,2,2)  !< 2nd derivative of shape function
    real(kind=kreal), intent(in), optional  :: normal(3)     !< noraml direction of surface
    real(kind=kreal), intent(out), optional :: curv(2,2)     !< curvature tensor
    real(kind=kreal) :: deriv2(nn,2,2)

    select case (fetype)
      case (fe_tri3n)
        !error check
        call Shape2ndDeriv_tri3n(deriv2(1:3,1:2,1:2))
      case (fe_tri6n)
        !error check
        call Shape2ndDeriv_tri6n(deriv2(1:6,1:2,1:2))
      case (fe_tri6nc)
        !error check
        call Shape2ndDeriv_tri6n(deriv2(1:6,1:2,1:2))
        ! deriv2=0.d0
      case (fe_quad4n)
        !error check
        call Shape2ndDeriv_quad4n(deriv2(1:4,1:2,1:2))
      case (fe_quad8n)
        !error check
        call Shape2ndDeriv_quad8n(localcoord,deriv2(1:8,1:2,1:2))
      case default
        ! error message
        stop "Cannot calculate second derivatives of shape function"
    end select

    l2ndderiv(1:3,1,1) = matmul( elecoord(1:3,1:nn), deriv2(1:nn,1,1) )
    l2ndderiv(1:3,1,2) = matmul( elecoord(1:3,1:nn), deriv2(1:nn,1,2) )
    l2ndderiv(1:3,2,1) = matmul( elecoord(1:3,1:nn), deriv2(1:nn,2,1) )
    l2ndderiv(1:3,2,2) = matmul( elecoord(1:3,1:nn), deriv2(1:nn,2,2) )
    if( present(curv) ) then
      curv(1,1) = dot_product( l2ndderiv(:,1,1), normal(:) )
      curv(1,2) = dot_product( l2ndderiv(:,1,2), normal(:) )
      curv(2,1) = dot_product( l2ndderiv(:,2,1), normal(:) )
      curv(2,2) = dot_product( l2ndderiv(:,2,2), normal(:) )
    endif
  end subroutine Curvature

  !> Return natural coordinate of the center of surface element
  subroutine getElementCenter( fetype, localcoord )
    integer, intent(in)           :: fetype            !< type of surface element
    real(kind=kreal), intent(out) :: localcoord(2)     !< center coordinate

    select case (fetype)
      case (fe_tri3n, fe_tri6n, fe_tri6nc)
        localcoord(:) = 1.d0/3.d0
      case (fe_quad4n, fe_quad8n)
        localcoord(:) = 0.d0
      case default
        ! error message
        localcoord(:) = 0.d0
    end select
  end subroutine getElementCenter

  !> if a point is inside a surface element
  !> -1: No; 0: Yes; >0: Node's (vertex) number
  integer function isInsideElement( fetype, localcoord, clearance )
    integer, intent(in)              :: fetype          !< type of surface element
    real(kind=kreal), intent(inout)  :: localcoord(2)   !< natural coord
    real(kind=kreal), optional       :: clearance       !< clearance used for judgement
    real(kind=kreal) :: clr, coord3

    clr = 1.d-6
    if( present(clearance) ) clr = clearance
    if( dabs(localcoord(1))<clr ) localcoord(1)=0.d0
    if( dabs(localcoord(2))<clr ) localcoord(2)=0.d0
    if( dabs(dabs(localcoord(1))-1.d0)<clr )    &
      localcoord(1)=sign(1.d0,localcoord(1))
    if( dabs(dabs(localcoord(2))-1.d0)<clr )    &
      localcoord(2)=sign(1.d0,localcoord(2))
    isInsideElement = -1
    select case (fetype)
      case (fe_tri3n, fe_tri6n, fe_tri6nc)
        !error check
        coord3 = 1.d0-(localcoord(1)+localcoord(2))
        if( dabs(coord3)<clr ) coord3=0.d0
        if( localcoord(1)>=0.d0 .and. localcoord(1)<=1.d0 .and.   &
            localcoord(2)>=0.d0 .and. localcoord(2)<=1.d0 .and.   &
            coord3>=0.d0 .and. coord3<=1.d0 )  then
          isInsideElement = 0
          if( localcoord(1)==1.d0 ) then
            isInsideElement = 1
          elseif( localcoord(2)==1.d0 ) then
            isInsideElement = 2
          elseif( coord3==1.d0 ) then
            isInsideElement = 3
          elseif( coord3==0.d0 ) then
            isInsideElement = 12
          elseif( localcoord(1)==0.d0 ) then
            isInsideElement = 23
          elseif( localcoord(2)==0.d0 ) then
            isInsideElement = 31
          endif
        endif
      case (fe_quad4n, fe_quad8n)
        !error check
        if( all(dabs(localcoord)<=1.d0) ) then
          isInsideElement = 0
          if( localcoord(1)==-1.d0 .and. localcoord(2)==-1.d0 ) then
            isInsideElement = 1
          elseif( localcoord(1)==1.d0 .and. localcoord(2)==-1.d0 ) then
            isInsideElement = 2
          elseif( localcoord(1)==1.d0 .and. localcoord(2)==1.d0 ) then
            isInsideElement = 3
          elseif( localcoord(1)==-1.d0 .and. localcoord(2)==1.d0 ) then
            isInsideElement = 4
          elseif( localcoord(2)==-1.d0 ) then
            isInsideElement = 12
          elseif( localcoord(1)==1.d0 ) then
            isInsideElement = 23
          elseif( localcoord(2)==1.d0 ) then
            isInsideElement = 34
          elseif( localcoord(1)==-1.d0 ) then
            isInsideElement = 41
          endif
        endif
    end select
  end function isInsideElement

  !> Get the natural coord of a vertex node
  subroutine getVertexCoord( fetype, cnode, localcoord )
    integer, intent(in)            :: fetype          !< type of surface element
    integer, intent(in)            :: cnode           !< current node
    real(kind=kreal), intent(out)  :: localcoord(2)   !< natural coord

    select case (fetype)
      case (fe_tri3n, fe_tri6n, fe_tri6nc)
        if( cnode==1 ) then
          localcoord(1) =1.d0
          localcoord(2) =0.d0
        elseif( cnode==2 ) then
          localcoord(1) =0.d0
          localcoord(2) =1.d0
        else
          localcoord(1) =0.d0
          localcoord(2) =0.d0
        endif
      case (fe_quad4n, fe_quad8n)
        if( cnode==1 ) then
          localcoord(1) =-1.d0
          localcoord(2) =-1.d0
        elseif( cnode==2 ) then
          localcoord(1) =1.d0
          localcoord(2) =-1.d0
        elseif( cnode==3 ) then
          localcoord(1) =1.d0
          localcoord(2) =1.d0
        else
          localcoord(1) =-1.d0
          localcoord(2) =1.d0
        endif
    end select
  end subroutine

  !> This subroutine extrapolate a point value into elemental nodes
  subroutine extrapolateValue( lpos, fetype, nnode, pvalue, ndvalue )
    real(kind=kreal), intent(in)  :: lpos(:)        !< poisition of value given
    integer, intent(in)           :: fetype         !< element type
    integer, intent(in)           :: nnode          !< number of element node
    real(kind=kreal), intent(in)  :: pvalue(:)      !< value to be extropolated
    real(kind=kreal), intent(out) :: ndvalue(:,:)   !< equivalent nodal value

    integer         :: i
    real(kind=kreal) :: shapefunc(nnode)
    call getShapeFunc( fetype, lpos, shapefunc )
    do i=1,nnode
      ndvalue(i,:) = shapefunc(i)*pvalue(:)
    enddo
  end subroutine

  !> This subroutine interapolate element nodes value into a point value
  subroutine interapolateValue( lpos, fetype, nnode, pvalue, ndvalue )
    real(kind=kreal), intent(in)  :: lpos(:)        !< poisition of value given
    integer, intent(in)           :: fetype         !< element type
    integer, intent(in)           :: nnode          !< number of element node
    real(kind=kreal), intent(out) :: pvalue(:)      !< value to be extropolated
    real(kind=kreal), intent(in)  :: ndvalue(:,:)   !< equivalent nodal value

    integer         :: i
    real(kind=kreal) :: shapefunc(nnode)
    call getShapeFunc( fetype, lpos, shapefunc )
    pvalue(:) = 0
    do i=1,nnode
      pvalue(:) = pvalue(:)+ shapefunc(i)*ndvalue(i,:)
    enddo
  end subroutine

  !> This subroutine extroplate value in quadrature point to element nodes
  subroutine Gauss2Node( fetype, gaussv, nodev )
    integer, intent(in)           :: fetype       !< element type
    real(kind=kreal), intent(in)  :: gaussv(:,:)  !< values in quadrature points
    real(kind=kreal), intent(out) :: nodev(:,:)   !< values in nodes

    integer :: i, ngauss, nnode
    real(kind=kreal) :: localcoord(3), func(100)
    ngauss = NumOfQuadPoints( fetype )
    nnode = getNumberOfNodes( fetype )
    ! error checking
    select case (fetype)
      case (fe_tri3n)
        !error check
        forall(i=1:nnode)
          nodev(i,:) = gaussv(1,:)
        end forall
      case (fe_tri6n)
        !error check
        !   func(1:6) = ShapeFunc_tri6n(localcoord)
      case (fe_quad4n)
        !error check
        !  nodev(:,:) = gaussv(1,:)
      case (fe_quad8n)
        !error check
        call ShapeFunc_quad8n(localcoord,func(1:8))
      case (fe_hex8n, fe_mitc4_shell361)
        ! error check
        call ShapeFunc_hex8n(localcoord,func(1:8))
      case (fe_hex20n)
        ! error check
        call ShapeFunc_hex20n(localcoord,func(1:20))
      case (fe_prism6n, fe_mitc3_shell361)
        forall(i=1:3)
          nodev(i,:) = gaussv(1,:)
        end forall
        forall(i=1:3)
          nodev(i+3,:) = gaussv(2,:)
        end forall
      case (fe_prism15n)
        call ShapeFunc_prism15n(localcoord,func(1:15))
      case (fe_tet4n, fe_tet4n_pipi, fe_beam341)
        ! error check
        forall(i=1:nnode)
          nodev(i,:) = gaussv(1,:)
        end forall
      case (fe_tet10n)
        ! error check
        call ShapeFunc_tet10n(localcoord,func(1:10))
      case default
        stop "Element type not defined"
        ! error message
    end select
  end subroutine

  !> This function calculates reference length at a point in surface
  real(kind=kreal) function getReferenceLength( fetype, nn, localcoord, elecoord )
    integer, intent(in)         :: fetype           !< surface element type
    integer, intent(in)         :: nn               !< number of elemental nodes
    real(kind=kreal),intent(in) :: localcoord(2)    !< natural coordinates
    real(kind=kreal),intent(in) :: elecoord(3,nn)   !< nodes coordinates of surface element
    real(kind=kreal) :: detJxy, detJyz, detJxz, detJ
    detJxy = getDeterminant( fetype, nn, localcoord, elecoord(1:2,1:nn) )
    detJyz = getDeterminant( fetype, nn, localcoord, elecoord(2:3,1:nn) )
    detJxz = getDeterminant( fetype, nn, localcoord, elecoord(1:3:2,1:nn) )
    detJ = dsqrt( detJxy **2 + detJyz **2 + detJxz **2 )
    getReferenceLength = dsqrt( detJ )
  end function getReferenceLength


end module
