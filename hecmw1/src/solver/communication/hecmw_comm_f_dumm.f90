 !======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.1                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module m_hecmw_comm_f

contains
!C DUMMY SUBROUTINE
!C
!C***
!C*** hecmw_barrier
!C***
!C
      subroutine hecmw_barrier (hecMESH)
      use hecmw_util
      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_barrier

      subroutine hecmw_scatterv_DP(sbuf, sc, disp, rbuf, rc, root, comm)
      use hecmw_util
      double precision   :: sbuf(:) !send buffer
      integer(kind=kint) :: sc      !send counts
      integer(kind=kint) :: disp    !displacement
      double precision   :: rbuf(:) !receive buffer
      integer(kind=kint) :: rc      !receive counts
      integer(kind=kint) :: root
      integer(kind=kint) :: comm

      rbuf(1) = sbuf(1)

      end subroutine hecmw_scatterv_DP

!C
!C***
!C*** hecmw_allREDUCE
!C***
!C
      subroutine hecmw_allreduce_DP(VAL,VALM,n,hec_op,comm )
      use hecmw_util
      integer(kind=kint) :: i,n, hec_op, comm
      real(kind=kreal) :: VAL(*)
      real(kind=kreal) :: VALM(*)
      do i=1,n
          VALM(i) = VAL(i)
      end do
      end subroutine hecmw_allreduce_DP

      subroutine hecmw_allreduce_DP1(VAL,VALM,hec_op,comm )
      use hecmw_util
      integer(kind=kint) ::  hec_op, comm
      real(kind=kreal) :: VAL
      real(kind=kreal) :: VALM
      VALM = VAL
      end subroutine hecmw_allreduce_DP1

!C
!C***
!C*** hecmw_allREDUCE_R
!C***
!C
      subroutine hecmw_allreduce_R (hecMESH, VAL, n, ntag)
      use hecmw_util
      integer(kind=kint):: n, ntag
      real(kind=kreal), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_allreduce_R

      subroutine hecmw_allreduce_R1 (hecMESH, s, ntag)
      use hecmw_util
      integer(kind=kint)::  ntag
      real(kind=kreal) ::s
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_allreduce_R1
!C
!C***
!C*** hecmw_allREDUCE_I
!C***
!C
      subroutine hecmw_allreduce_I (hecMESH, VAL, n, ntag)
      use hecmw_util
      integer(kind=kint):: n, ntag
      integer(kind=kint), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_allreduce_I

      subroutine hecmw_allreduce_I1 (hecMESH, s, ntag)
      use hecmw_util
      integer(kind=kint)::  ntag, s
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_allreduce_I1
!C
!C***
!C*** hecmw_bcast_R
!C***
!C
      subroutine hecmw_bcast_R (hecMESH, VAL, n, nbase)
      use hecmw_util
      integer(kind=kint):: n, nbase
      real(kind=kreal), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_R

      subroutine hecmw_bcast_R1 (hecMESH, VAL, nbase)
      use hecmw_util
      integer(kind=kint):: nbase
      real(kind=kreal) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_R1
!C
!C***
!C*** hecmw_bcast_I
!C***
!C

      subroutine hecmw_bcast_I (hecMESH, VAL, n, nbase)
      use hecmw_util
      integer(kind=kint):: n, nbase
      integer(kind=kint), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_I

      subroutine hecmw_bcast_I1 (hecMESH, VAL, nbase)
      use hecmw_util
      integer(kind=kint)::  nbase
      integer(kind=kint):: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_I1
!C
!C***
!C*** hecmw_bcast_C
!C***
!C
      subroutine hecmw_bcast_C (hecMESH, VAL, n, nn, nbase)
      use hecmw_util
      integer(kind=kint):: n, nn, nbase
      character(len=n) :: VAL(nn)
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_bcast_C

!C
!C***
!C*** hecmw_update_1_R
!C***
!C
!C    1-DOF, REAL
!C
      subroutine hecmw_update_1_R (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      real(kind=kreal), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_1_R

!C
!C***
!C*** hecmw_update_2_R
!C***
!C
!C    2-DOF, REAL
!C
      subroutine hecmw_update_2_R (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      real(kind=kreal), dimension(2*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_2_R

!C
!C***
!C*** hecmw_update_3_R
!C***
!C
!C    3-DOF, REAL
!C
      subroutine hecmw_update_3_R (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      real(kind=kreal), dimension(3*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_3_R

!C
!C***
!C*** hecmw_update_m_R
!C***
!C
!C    m-DOF, REAL
!C
      subroutine hecmw_update_m_R (hecMESH, VAL, n, m)
      use hecmw_util
      integer(kind=kint):: n, m
      real(kind=kreal), dimension(m*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_m_R

!C
!C***
!C*** hecmw_update_1_I
!C***
!C
!C    1-DOF, INTEGER
!C
      subroutine hecmw_update_1_I (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      integer(kind=kint), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_1_I

!C
!C***
!C*** hecmw_update_2_I
!C***
!C
!C    2-DOF, INTEGER
!C
      subroutine hecmw_update_2_I (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      integer(kind=kint), dimension(2*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_2_I

!C
!C***
!C*** hecmw_update_3_I
!C***
!C
!C    3-DOF, INTEGER
!C
      subroutine hecmw_update_3_I (hecMESH, VAL, n)
      use hecmw_util
      integer(kind=kint):: n
      integer(kind=kint), dimension(3*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_3_I

!C
!C***
!C*** hecmw_update_m_I
!C***
!C
!C    m-DOF, REAL
!C
      subroutine hecmw_update_m_I (hecMESH, VAL, n, m)
      use hecmw_util
      integer(kind=kint):: n, m
      integer(kind=kint), dimension(m*n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      end subroutine hecmw_update_m_I

end module m_hecmw_comm_f



