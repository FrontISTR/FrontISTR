!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
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
      ! == add futher defintiions here ==
   end type

   type(tULoad), pointer, save :: uloads(:)=>null()

contains

!> This suborutine read in variables needs to define user-defined external loads
   integer function ureadload( fname )
   character(len=*), intent(in)    :: fname   !< input file name
     ureadload = 0
   end function

!> This subroutine take consider of user-defined external loading
   subroutine uloading( cstep, factor, exForce )
     integer, INTENT(IN)             :: cstep      !< current step number
     REAL(KIND=kreal), INTENT(IN)    :: factor     !< loading factor of current step
     REAL(KIND=kreal), INTENT(INOUT) :: exForce(:) !< external force

   end subroutine

!> This subroutine take consider of user-defined external loading
   subroutine uResidual( cstep, factor, residual )
     integer, INTENT(IN)             :: cstep        !< current step number
     REAL(KIND=kreal), INTENT(IN)    :: factor       !< loading factor of current step
     REAL(KIND=kreal), INTENT(INOUT) :: residual(:)  !< residual

   end subroutine

end module

