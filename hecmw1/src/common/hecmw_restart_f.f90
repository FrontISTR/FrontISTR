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

module hecmw_restart
    use hecmw_util
    implicit none

    contains

!C=============================================================================
!C Write restart data to file
!C=============================================================================
    subroutine hecmw_restart_add_int(src, n_data)
        integer(kind=kint),dimension(:) :: src
        integer(kind=kint) :: n_data,ierr

        call hecmw_restart_add_int_if(src, n_data, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_add_int


    subroutine hecmw_restart_add_real(src, n_data)
        real(kind=kreal),dimension(:) :: src
        integer(kind=kint) :: n_data,ierr

        call hecmw_restart_add_real_if(src, n_data, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_add_real


    subroutine hecmw_restart_write()
        integer(kind=kint) :: ierr

        call hecmw_restart_write_if(ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_write


    subroutine hecmw_restart_write_by_name(name_ID)
        integer(kind=kint) :: ierr
        character(len=HECMW_NAME_LEN) :: name_ID 

        call hecmw_restart_write_by_name_if(name_ID, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_write_by_name


!C=============================================================================
!C Read restart data from file
!C=============================================================================
    subroutine hecmw_restart_open()
        integer(kind=kint) :: ierr

        call hecmw_restart_open_if(ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_open


    subroutine hecmw_restart_open_by_name(name_ID)
        integer(kind=kint) :: ierr
        character(len=HECMW_NAME_LEN) :: name_ID

        call hecmw_restart_open_by_name_if(name_ID, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_open_by_name


    subroutine hecmw_restart_close()
        integer(kind=kint) :: ierr

        call hecmw_restart_close_if(ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_close


    subroutine hecmw_restart_read_int(dst)
        integer(kind=kint) :: ierr  
        integer(kind=kint),dimension(:) :: dst

        call hecmw_restart_read_int_if(dst, ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_read_int


    subroutine hecmw_restart_read_real(dst)
        integer(kind=kint) :: ierr  
        real(kind=kreal),dimension(:) :: dst

        call hecmw_restart_read_real_if(dst, ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_restart_read_real

end module hecmw_restart

