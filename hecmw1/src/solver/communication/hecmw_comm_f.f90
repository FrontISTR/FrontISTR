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

!C
!C***
!C*** hecmw_barrier
!C***
!C
      subroutine hecmw_barrier (hecMESH)
      use hecmw_util
      implicit none
      integer(kind=kint):: ierr
      type (hecmwST_local_mesh) :: hecMESH

      call MPI_BARRIER (hecMESH%MPI_COMM, ierr)
      end subroutine hecmw_barrier

      subroutine hecmw_scatterv_DP(sbuf, sc, disp, rbuf, rc, root, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sc      !send counts
      double precision   :: sbuf(sc) !send buffer
      integer(kind=kint) :: disp    !displacement
      integer(kind=kint) :: rc      !receive counts
      double precision   :: rbuf(rc) !receive buffer
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr

      CALL MPI_scatterv( sbuf, sc, disp, MPI_DOUBLE_PRECISION, &
                         rbuf, rc, MPI_DOUBLE_PRECISION, &
                         root, comm, ierr )

      end subroutine hecmw_scatterv_DP

      function hecmw_operation_hec2mpi(operation)
      use hecmw_util
      implicit none
      integer(kind=kint) :: hecmw_operation_hec2mpi
      integer(kind=kint) :: operation
      if (operation == HECMW_SUM) then
        hecmw_operation_hec2mpi = MPI_SUM
      elseif (operation == HECMW_PROD) then
        hecmw_operation_hec2mpi = MPI_PROD
      elseif (operation == HECMW_MAX) then
        hecmw_operation_hec2mpi = MPI_MAX
      elseif (operation == HECMW_MIN) then
        hecmw_operation_hec2mpi = MPI_MIN
      endif
      end function hecmw_operation_hec2mpi

      subroutine hecmw_scatter_int_1(sbuf, rval, root, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sbuf(*) !send buffer
      integer(kind=kint) :: rval !receive value
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      CALL MPI_scatter( sbuf, 1, MPI_INTEGER, &
     &     rval, 1, MPI_INTEGER, root, comm, ierr )
      end subroutine hecmw_scatter_int_1

      subroutine hecmw_gather_int_1(sval, rbuf, root, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sval !send buffer
      integer(kind=kint) :: rbuf(*) !receive buffer
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      CALL MPI_gather( sval, 1, MPI_INTEGER, &
     &     rbuf, 1, MPI_INTEGER, root, comm, ierr )
      end subroutine hecmw_gather_int_1

      subroutine hecmw_allgather_int_1(sval, rbuf, comm)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sval !send buffer
      integer(kind=kint) :: rbuf(*) !receive buffer
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      CALL MPI_allgather( sval, 1, MPI_INTEGER, &
     &     rbuf, 1, MPI_INTEGER, comm, ierr )
      end subroutine hecmw_allgather_int_1

      subroutine hecmw_scatterv_real(sbuf, scs, disp, &
     &     rbuf, rc, root, comm)
      use hecmw_util
      implicit none
      real(kind=kreal)   :: sbuf(*) !send buffer
      integer(kind=kint) :: scs(*) !send counts
      integer(kind=kint) :: disp(*) !displacement
      real(kind=kreal)   :: rbuf(*) !receive buffer
      integer(kind=kint) :: rc  !receive counts
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      CALL MPI_scatterv( sbuf, scs, disp, MPI_REAL8, &
     &     rbuf, rc, MPI_REAL8, root, comm, ierr )
      end subroutine hecmw_scatterv_real

      subroutine hecmw_gatherv_real(sbuf, sc, &
     &     rbuf, rcs, disp, root, comm)
      use hecmw_util
      implicit none
      real(kind=kreal)   :: sbuf(*) !send buffer
      integer(kind=kint) :: sc  !send counts
      real(kind=kreal)   :: rbuf(*) !receive buffer
      integer(kind=kint) :: rcs(*) !receive counts
      integer(kind=kint) :: disp(*) !displacement
      integer(kind=kint) :: root
      integer(kind=kint) :: comm
      integer(kind=kint) :: ierr
      CALL MPI_gatherv( sbuf, sc, MPI_REAL8, &
     &     rbuf, rcs, disp, MPI_REAL8, root, comm, ierr )
      end subroutine hecmw_gatherv_real

      subroutine hecmw_allreduce_int_1(sval, rval, op, comm)
      use hecmw_util
      implicit none
      integer(kind=kint)   :: sval !send val
      integer(kind=kint)   :: rval !receive val
      integer(kind=kint):: op, comm, ierr
      call MPI_ALLREDUCE(sval, rval, 1, MPI_INTEGER, &
     &     hecmw_operation_hec2mpi(op), comm, ierr)
      end subroutine hecmw_allreduce_int_1

      subroutine hecmw_isend_int(sbuf, sc, dest, &
     &     tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sbuf(*)
      integer(kind=kint) :: sc
      integer(kind=kint) :: dest
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      integer(kind=kint) :: ierr
      call MPI_ISEND(sbuf, sc, MPI_INTEGER, &
     &     dest, tag, comm, req, ierr)
      end subroutine hecmw_isend_int
      
      subroutine hecmw_isend_r(sbuf, sc, dest, &
     &     tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: sc
      double precision, dimension(sc) :: sbuf
      integer(kind=kint) :: dest
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      integer(kind=kint) :: ierr
      call MPI_ISEND(sbuf, sc, MPI_DOUBLE_PRECISION, &
     &     dest, tag, comm, req, ierr)
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
      call MPI_IRECV(rbuf, rc, MPI_INTEGER, &
     &     source, tag, comm, req, ierr)
      end subroutine hecmw_irecv_int
      
      subroutine hecmw_irecv_r(rbuf, rc, source, &
     &     tag, comm, req)
      use hecmw_util
      implicit none
      integer(kind=kint) :: rc
      double precision, dimension(rc) :: rbuf
      integer(kind=kint) :: source
      integer(kind=kint) :: tag
      integer(kind=kint) :: comm
      integer(kind=kint) :: req
      integer(kind=kint) :: ierr
      call MPI_IRECV(rbuf, rc, MPI_DOUBLE_PRECISION, &
     &     source, tag, comm, req, ierr)
      end subroutine hecmw_irecv_r

      subroutine hecmw_waitall(cnt, reqs, stats)
      use hecmw_util
      implicit none
      integer(kind=kint) :: cnt
      integer(kind=kint) :: reqs(*)
      integer(kind=kint) :: stats(MPI_STATUS_SIZE,*)
      integer(kind=kint) :: ierr
      call MPI_WAITALL(cnt, reqs, stats, ierr)
      end subroutine hecmw_waitall
!C
!C***
!C*** hecmw_allREDUCE
!C***
!C
      subroutine hecmw_allreduce_DP(VAL,VALM,n,hec_op,comm )
      use hecmw_util
      implicit none
      integer(kind=kint) :: n, hec_op,op, comm, ierr
      double precision, dimension(n) :: VAL
      double precision, dimension(n) :: VALM

      select case( hec_op )
      case ( hecmw_sum )
          op = MPI_SUM
      case ( hecmw_prod )
          op = MPI_PROD
      case ( hecmw_max )
          op = MPI_MAX
      case ( hecmw_min )
          op = MPI_MIN
      end select
      call MPI_allREDUCE(VAL,VALM,n,MPI_DOUBLE_PRECISION,op,comm,ierr)

      end subroutine hecmw_allREDUCE_DP

      subroutine hecmw_allreduce_DP1(s1,s2,hec_op,comm )
      use hecmw_util
      implicit none
      integer(kind=kint) ::  hec_op, comm
      double precision :: s1, s2
      double precision, dimension(1) :: VAL
      double precision, dimension(1) :: VALM
      VAL(1) = s1
      VALM(1) = s2
      call hecmw_allreduce_DP(VAL,VALM,1,hec_op,comm )
      s1 = VAL(1)
      s2 = VALM(1)
      end subroutine hecmw_allreduce_DP1
!C
!C***
!C*** hecmw_allREDUCE_R
!C***
!C
      subroutine hecmw_allreduce_R (hecMESH, VAL, n, ntag)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, ntag, ierr
      real(kind=kreal), dimension(n) :: VAL
      real(kind=kreal), dimension(:), allocatable :: VALM
      type (hecmwST_local_mesh) :: hecMESH

      allocate (VALM(n))
      VALM= 0.d0
      if (ntag .eq. hecmw_sum) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_DOUBLE_PRECISION, MPI_SUM,              &
     &        hecMESH%MPI_COMM, ierr)
      endif

      if (ntag .eq. hecmw_max) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_DOUBLE_PRECISION, MPI_MAX,              &
     &        hecMESH%MPI_COMM, ierr)
      endif

      if (ntag .eq. hecmw_min) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_DOUBLE_PRECISION, MPI_MIN,              &
     &        hecMESH%MPI_COMM, ierr)
      endif

      VAL= VALM
      deallocate (VALM)

      end subroutine hecmw_allreduce_R

      subroutine hecmw_allreduce_R1 (hecMESH, s, ntag)
      use hecmw_util
      implicit none
      integer(kind=kint):: ntag
      real(kind=kreal) :: s
      real(kind=kreal), dimension(1) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      VAL(1) = s
      call hecmw_allreduce_R(hecMESH, VAL, 1, ntag )
      s = VAL(1)
      end subroutine hecmw_allreduce_R1

!C
!C***
!C*** hecmw_allREDUCE_I
!C***
!C
      subroutine hecmw_allreduce_I(hecMESH, VAL, n, ntag)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, ntag, ierr
      integer(kind=kint), dimension(n) :: VAL
      integer(kind=kint), dimension(:), allocatable :: VALM
      type (hecmwST_local_mesh) :: hecMESH

      allocate (VALM(n))
      VALM= 0
      if (ntag .eq. hecmw_sum) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_INTEGER, MPI_SUM,                       &
     &        hecMESH%MPI_COMM, ierr)
      endif

      if (ntag .eq. hecmw_max) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_INTEGER, MPI_MAX,                       &
     &        hecMESH%MPI_COMM, ierr)
      endif

      if (ntag .eq. hecmw_min) then
        call MPI_allREDUCE                                              &
     &       (VAL, VALM, n, MPI_INTEGER, MPI_MIN,                       &
     &        hecMESH%MPI_COMM, ierr)
      endif


      VAL= VALM
      deallocate (VALM)
      end subroutine hecmw_allreduce_I

      subroutine hecmw_allreduce_I1 (hecMESH, s, ntag)
      use hecmw_util
      implicit none
      integer(kind=kint)::  ntag, s
      integer(kind=kint), dimension(1) :: VAL
      type (hecmwST_local_mesh) :: hecMESH

      VAL(1) = s
      call hecmw_allreduce_I(hecMESH, VAL, 1, ntag )
      s = VAL(1)
      end subroutine hecmw_allreduce_I1

!C
!C***
!C*** hecmw_bcast_R
!C***
!C
      subroutine hecmw_bcast_R (hecMESH, VAL, n, nbase)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nbase, ierr
      real(kind=kreal), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      call MPI_BCAST (VAL, n, MPI_DOUBLE_PRECISION, nbase, hecMESH%MPI_COMM, ierr)
      end subroutine hecmw_bcast_R
      
      subroutine hecmw_bcast_R_comm (VAL, n, nbase, comm)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nbase, ierr
      real(kind=kreal), dimension(n) :: VAL
      integer(kind=kint):: comm
      call MPI_BCAST (VAL, n, MPI_DOUBLE_PRECISION, nbase, comm, ierr)
      end subroutine hecmw_bcast_R_comm

      subroutine hecmw_bcast_R1 (hecMESH, s, nbase)
      use hecmw_util
      implicit none
      integer(kind=kint):: nbase, ierr
      real(kind=kreal) :: s
      real(kind=kreal), dimension(1) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      VAL(1)=s
      call MPI_BCAST (VAL, 1, MPI_DOUBLE_PRECISION, nbase, hecMESH%MPI_COMM, ierr)
      s = VAL(1)
      end subroutine hecmw_bcast_R1
      
      subroutine hecmw_bcast_R1_comm (s, nbase, comm)
      use hecmw_util
      implicit none
      integer(kind=kint):: nbase, ierr
      real(kind=kreal) :: s
      real(kind=kreal), dimension(1) :: VAL
      integer(kind=kint):: comm
      VAL(1)=s
      call MPI_BCAST (VAL, 1, MPI_DOUBLE_PRECISION, nbase, comm, ierr)
      s = VAL(1)
      end subroutine hecmw_bcast_R1_comm
!C
!C***
!C*** hecmw_bcast_I
!C***
!C
      subroutine hecmw_bcast_I (hecMESH, VAL, n, nbase)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nbase, ierr
      integer(kind=kint), dimension(n) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      call MPI_BCAST (VAL, n, MPI_INTEGER, nbase, hecMESH%MPI_COMM, ierr)
      end subroutine hecmw_bcast_I
      
      subroutine hecmw_bcast_I_comm (VAL, n, nbase, comm)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nbase, ierr
      integer(kind=kint), dimension(n) :: VAL
      integer(kind=kint):: comm
      call MPI_BCAST (VAL, n, MPI_INTEGER, nbase, comm, ierr)
      end subroutine hecmw_bcast_I_comm

      subroutine hecmw_bcast_I1 (hecMESH, s, nbase)
      use hecmw_util
      implicit none
      integer(kind=kint):: nbase, ierr,s
      integer(kind=kint), dimension(1) :: VAL
      type (hecmwST_local_mesh) :: hecMESH
      VAL(1) = s
      call MPI_BCAST (VAL, 1, MPI_INTEGER, nbase, hecMESH%MPI_COMM, ierr)
      s = VAL(1)
      end subroutine hecmw_bcast_I1
      
      subroutine hecmw_bcast_I1_comm (s, nbase, comm)
      use hecmw_util
      implicit none
      integer(kind=kint):: nbase, ierr,s
      integer(kind=kint), dimension(1) :: VAL
      integer(kind=kint):: comm
      VAL(1) = s
      call MPI_BCAST (VAL, 1, MPI_INTEGER, nbase, comm, ierr)
      s = VAL(1)
      end subroutine hecmw_bcast_I1_comm
!C
!C***
!C*** hecmw_bcast_C
!C***
!C
      subroutine hecmw_bcast_C (hecMESH, VAL, n, nn, nbase)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nn, nbase, ierr
      character(len=n) :: VAL(nn)
      type (hecmwST_local_mesh) :: hecMESH

      call MPI_BCAST (VAL, n*nn, MPI_CHARACTER, nbase, hecMESH%MPI_COMM,&
     &                                                 ierr)
      end subroutine hecmw_bcast_C
      
      subroutine hecmw_bcast_C_comm (VAL, n, nn, nbase, comm)
      use hecmw_util
      implicit none
      integer(kind=kint):: n, nn, nbase, ierr
      character(len=n) :: VAL(nn)
      integer(kind=kint):: comm

      call MPI_BCAST (VAL, n*nn, MPI_CHARACTER, nbase, comm,&
     &                                                 ierr)
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
      use  hecmw_solver_SR_11

      implicit none
      integer(kind=kint):: n, ns, nr
      real(kind=kreal), dimension(n) :: VAL
      real(kind=kreal), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(ns), WR(nr))
      call hecmw_solve_SEND_RECV_11                                     &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_22

      implicit none
      integer(kind=kint):: n, ns, nr
      real(kind=kreal), dimension(2*n) :: VAL
      real(kind=kreal), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(2*ns), WR(2*nr))
      call hecmw_solve_SEND_RECV_22                                     &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_33

      implicit none
      integer(kind=kint):: n, ns, nr
      real(kind=kreal), dimension(3*n) :: VAL
      real(kind=kreal), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(3*ns), WR(3*nr))
      call hecmw_solve_SEND_RECV_33                                     &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_mm

      implicit none
      integer(kind=kint):: n, m, ns, nr
      real(kind=kreal), dimension(m*n) :: VAL
      real(kind=kreal), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(m*ns), WR(m*nr))
      call hecmw_solve_SEND_RECV_mm                                     &
     &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_11i

      implicit none
      integer(kind=kint):: n, ns, nr
      integer(kind=kint), dimension(n) :: VAL
      integer(kind=kint), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(ns), WR(nr))
      call hecmw_solve_SEND_RECV_11i                                    &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_22i

      implicit none
      integer(kind=kint):: n, ns, nr
      integer(kind=kint), dimension(2*n) :: VAL
      integer(kind=kint), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(2*ns), WR(2*nr))
      call hecmw_solve_SEND_RECV_22i                                    &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_33i

      implicit none
      integer(kind=kint):: n, ns, nr
      integer(kind=kint), dimension(3*n) :: VAL
      integer(kind=kint), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(3*ns), WR(3*nr))
      call hecmw_solve_SEND_RECV_33i                                    &
     &   ( n, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

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
      use  hecmw_solver_SR_mmi

      implicit none
      integer(kind=kint):: n, m, ns, nr
      integer(kind=kint), dimension(m*n) :: VAL
      integer(kind=kint), dimension(:), allocatable :: WS, WR
      type (hecmwST_local_mesh) :: hecMESH

      if( hecMESH%n_neighbor_pe == 0 ) return

      ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
      nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

      allocate (WS(m*ns), WR(m*nr))
      call hecmw_solve_SEND_RECV_mmi                                    &
     &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, VAL , hecMESH%MPI_COMM, hecMESH%my_rank)
      deallocate (WS, WR)

      end subroutine hecmw_update_m_I

end module m_hecmw_comm_f


