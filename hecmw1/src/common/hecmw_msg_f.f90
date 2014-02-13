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

module hecmw_msg
!      use hecmw_util
      use hecmw_msgno

      contains

      subroutine hecmw_strmsg(msgno, msg)
          integer(kind=kint) :: msgno
          character(len=HECMW_MSG_LEN) :: msg

          call hecmw_strmsg_if(msgno, msg)
      end subroutine hecmw_strmsg
end module hecmw_msg
