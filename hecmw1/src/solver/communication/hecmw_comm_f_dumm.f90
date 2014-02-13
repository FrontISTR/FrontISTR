 !======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
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
      implicit none
      double precision   :: sbuf(:) !send buffer
      integer(kind=kint) :: sc      !send counts
      integer(kind=kint) :: disp    !displacement
      double precision   :: rbuf(:) !receive buffer
      integer(kind=kint) :: rc      !receive counts
      integer(kind=kint) :: root
      integer(kind=kint) :: comm

      rbuf(1) = sbuf(1)

      end subroutine hecmw_scatterv_DP

      subroutine hecmw_scatter_int_1(sbuf, rval, root, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sbuf(*) !send buffer
      integer(kind=kint) :: rval    !receive buffer
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      rval=sbuf(1)
      end subroutine hecmw_scatter_int_1

      subroutine hecmw_gather_int_1(sval, rbuf, root, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sval    !send buffer
      integer(kind=kint) :: rbuf(*) !receive buffer
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      rbuf(1)=sval
      end subroutine hecmw_gather_int_1

      subroutine hecmw_allgather_int_1(sval, rbuf, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sval    !send buffer
      integer(kind=kint) :: rbuf(*) !receive buffer
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      rbuf(1)=sval
      end subroutine hecmw_allgather_int_1

      subroutine hecmw_scatterv_real(sbuf, scs, disp, &
     &     rbuf, rc, root, comm)
      use hecmw_util
      implicit none
      real(kind=kreal)   :: sbuf(*) !send buffer
      integer(kind=kint) :: scs(*)  !send counts
      integer(kind=kint) :: disp(*) !displacement
      real(kind=kreal)   :: rbuf(*) !receive buffer
      integer(kind=kint) :: rc      !receive counts
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      rbuf(1:rc)=sbuf(1:rc)
      end subroutine hecmw_scatterv_real

      subroutine hecmw_gatherv_real(sbuf, sc, &
     &     rbuf, rcs, disp, root, comm)
      use hecmw_util
      implicit none
      real(kind=kreal)   :: sbuf(*) !send buffer
      integer(kind=kint) :: sc      !send counts
      real(kind=kreal)   :: rbuf(*) !receive buffer
      integer(kind=kint) :: rcs(*)  !receive counts
      integer(kind=kint) :: disp(*) !displacement
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      rbuf(1:sc)=sbuf(1:sc)
      end subroutine hecmw_gatherv_real

      subroutine hecmw_allreduce_int_1(sval, rval, op, comm)
      use hecmw_util
      implicit none
      integer(kind=kint)   :: sval !send buffer
      integer(kind=kint)   :: rval !receive buffer
      integer(kind=kint):: op, comm, ierr
      rval=sval
      end subroutine hecmw_allreduce_int_1

      subroutine hecmw_isend_int(sbuf, sc, dest, tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sbuf(*)
      integer(kind=kint) :: sc
      integer(kind=kint) :: dest
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      return
      end subroutine hecmw_isend_int
      
      subroutine hecmw_isend_r(sbuf, sc, dest, tag, comm, req)
      use hecmw_util
      implicit none
      real(kind=kreal) :: sbuf(*)
      integer(kind=kint) :: sc
      integer(kind=kint) :: dest
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      return
      end subroutine hecmw_isend_r
	  
      subroutine hecmw_irecv_int(rbuf, rc, source, &
     &     tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: rbuf(*)
      integer(kind=kint) :: rc
      integer(kind=kint) :: source
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      integer(kind=kint) :: ierr
      end subroutine hecmw_irecv_int
      
      subroutine hecmw_irecv_r(rbuf, rc, source, &
     &     tag, comm, req)
      use hecmw_util
      implicit none
      real(kind=kreal) :: rbuf(*)
      integer(kind=kint) :: rc
      integer(kind=kint) :: source
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      integer(kind=kint) :: ierr
      end subroutine hecmw_irecv_r

      subroutine hecmw_irecv(rbuf, rc, source, tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: rbuf(*)
      integer(kind=kint) :: rc
      integer(kind=kint) :: source
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      return
      end subroutine hecmw_irecv

      subroutine hecmw_waitall(cnt, reqs, stats)
      use hecmw_util
      implicit none
      integer(kind=kint) :: cnt
      integer(kind=kint) :: reqs(*)
      integer(kind=kint) :: stats(HECMW_STATUS_SIZE,*)
      return
      end subroutine hecmw_waitall
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
      
      subroutine hecmw_bcast_R_comm (VAL, n, nbase, comm)
      use hecmw_util
      integer(kind=kint):: n, nbase
      real(kind=kreal), dimension(n) :: VAL
      integer(kind=kint):: comm
      end subroutine hecmw_bcast_R_comm

      subroutine hecmw_bcast_R1 (hecMESH, VAL, nbase)
      use hecmw_util
      integer(kind=kint):: nbase
      real(kind=kreal) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_R1
      
      subroutine hecmw_bcast_R1_comm (VAL, nbase, comm)
      use hecmw_util
      integer(kind=kint):: nbase
      real(kind=kreal):: VAL
      integer(kind=kint):: comm
      end subroutine hecmw_bcast_R1_comm
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
      
      subroutine hecmw_bcast_I_comm (VAL, n, nbase, comm)
      use hecmw_util
      integer(kind=kint):: n, nbase
      integer(kind=kint), dimension(n) :: VAL
      integer(kind=kint):: comm
      end subroutine hecmw_bcast_I_comm

      subroutine hecmw_bcast_I1 (hecMESH, VAL, nbase)
      use hecmw_util
      integer(kind=kint)::  nbase
      integer(kind=kint):: VAL
      type (hecmwST_local_mesh) :: hecMESH
      end subroutine hecmw_bcast_I1

      subroutine hecmw_bcast_I1_comm (VAL, nbase, comm)
      use hecmw_util
      integer(kind=kint):: nbase
      integer(kind=kint):: VAL
      integer(kind=kint):: comm
      end subroutine hecmw_bcast_I1_comm
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
      
      subroutine hecmw_bcast_C_comm (VAL, n, nn, nbase, comm)
      use hecmw_util
      integer(kind=kint):: n, nn, nbase
      character(len=n) :: VAL(nn)
      integer(kind=kint):: comm

      end subroutine hecmw_bcast_C_comm

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



