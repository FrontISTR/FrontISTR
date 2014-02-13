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

module hecmw_logging
    use hecmw_util
    implicit none

    integer(kind=kint),parameter :: HECMW_LOG_NONE  = 0
    integer(kind=kint),parameter :: HECMW_LOG_ERROR = 1
    integer(kind=kint),parameter :: HECMW_LOG_WARN  = 2
    integer(kind=kint),parameter :: HECMW_LOG_INFO  = 4
    integer(kind=kint),parameter :: HECMW_LOG_DEBUG = 8
    integer(kind=kint),parameter :: HECMW_LOG_ALL   = 15

    contains

    subroutine hecmw_log(loglv, msg)
        integer(kind=kint) :: loglv
        character(len=HECMW_MSG_LEN) :: msg

        call hecmw_log_if(loglv, msg)
    end subroutine hecmw_log


    subroutine hecmw_setloglv(loglv)
        integer(kind=kint) :: loglv

        call hecmw_setloglv_if(loglv)
    end subroutine hecmw_setloglv
    

    subroutine hecmw_log_set_enable(from, to, true_or_false )
        integer(kind=kint) :: from, to, true_or_false

        call hecmw_log_set_enable_if(from, to, true_or_false)
    end subroutine hecmw_log_set_enable


    subroutine hecmw_openlog(logfile, loglv, options, id, ierror)
        character(len=HECMW_FILENAME_LEN) :: logfile
        integer(kind=kint) :: loglv,options,id,ierror

        call hecmw_openlog_if(logfile, loglv, options, id, ierror)
    end subroutine hecmw_openlog


    subroutine hecmw_closelog(id, ierror)
        integer(kind=kint) :: id,ierror

        call hecmw_closelog_if(id, ierror)
    end subroutine hecmw_closelog
end module hecmw_logging
