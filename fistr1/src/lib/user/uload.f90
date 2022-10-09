!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined loading
!>  tangent
module mULoad
  use hecmw
  implicit none

  !> Structure for user defines load. User may need to fill in it
  !> according to specified loads
  type tULoad
    integer, pointer :: nodeID(:)=>null()    !< nodes' ID
    integer, pointer :: dof(:)=>null()       !< dof to be loaded
    ! == add further definitions here ==
  end type

  type(tULoad), pointer, save :: uloads(:)=>null()

contains

  !> This subroutine read in variables needs to define user-defined external loads
  integer function ureadload( fname )
    character(len=*), intent(in)    :: fname   !< input file name
    ureadload = 0
  end function

  !> This subroutine take consider of user-defined external loading
  subroutine uloading( cstep, factor, exForce )
    integer, intent(in)             :: cstep      !< current step number
    real(kind=kreal), intent(in)    :: factor     !< loading factor of current step
    real(kind=kreal), intent(inout) :: exForce(:) !< external force

  end subroutine

  !> This subroutine take consider of user-defined external loading
  subroutine uResidual( cstep, factor, residual )
    integer, intent(in)             :: cstep        !< current step number
    real(kind=kreal), intent(in)    :: factor       !< loading factor of current step
    real(kind=kreal), intent(inout) :: residual(:)  !< residual

  end subroutine

end module

