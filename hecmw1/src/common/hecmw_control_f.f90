!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuaki Sakane (RIST)                          !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_control
    use hecmw_util
    implicit none

    contains

    subroutine hecmw_ctrl_get_control_file(name_ID, filename)
        character(len=HECMW_NAME_LEN) :: name_ID
        character(len=HECMW_FILENAME_LEN) :: filename
        integer(kind=kint) :: ierr

        call hecmw_ctrl_get_control_file_if(name_ID, filename, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_ctrl_get_control_file

end module hecmw_control

