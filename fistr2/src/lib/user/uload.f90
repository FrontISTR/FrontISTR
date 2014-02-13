!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!>  \brief   This subroutine read in used-defined loading
!>  tangent 
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/01/30
!>  \version    0.00
!======================================================================!
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
   subroutine uloading( cstep, factor, force )
     integer, INTENT(IN)             :: cstep     !< current step number
     REAL(KIND=kreal), INTENT(IN)    :: factor    !< loading factor of current step
     REAL(KIND=kreal), INTENT(INOUT) :: force(:)  !< external force
   end subroutine

!> This subroutine take consider of user-defined external loading
   subroutine uResidual( cstep, factor )
     integer, INTENT(IN)             :: cstep        !< current step number
     REAL(KIND=kreal), INTENT(IN)    :: factor       !< loading factor of current step
	
   end subroutine

end module

