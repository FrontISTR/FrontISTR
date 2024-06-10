!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

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
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BARRIER (hecMESH%MPI_COMM, ierr)
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_scatterv( sbuf, sc, disp, MPI_DOUBLE_PRECISION, &
      rbuf, rc, MPI_DOUBLE_PRECISION, &
      root, comm, ierr )
#else
    rbuf(1) = sbuf(1)
#endif
  end subroutine hecmw_scatterv_DP

#ifndef HECMW_SERIAL
  function hecmw_operation_hec2mpi(operation)
    use hecmw_util
    implicit none
    integer(kind=kint) :: hecmw_operation_hec2mpi
    integer(kind=kint) :: operation
    hecmw_operation_hec2mpi = -1
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
#endif

  subroutine hecmw_scatter_int_1(sbuf, rval, root, comm)
    use hecmw_util
    implicit none
    integer(kind=kint) :: sbuf(*) !send buffer
    integer(kind=kint) :: rval !receive value
    integer(kind=kint) :: root
    integer(kind=kint) :: comm
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_scatter( sbuf, 1, MPI_INTEGER, &
      &     rval, 1, MPI_INTEGER, root, comm, ierr )
#else
    rval=sbuf(1)
#endif
  end subroutine hecmw_scatter_int_1

  subroutine hecmw_gather_int_1(sval, rbuf, root, comm)
    use hecmw_util
    implicit none
    integer(kind=kint) :: sval !send buffer
    integer(kind=kint) :: rbuf(*) !receive buffer
    integer(kind=kint) :: root
    integer(kind=kint) :: comm
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_gather( sval, 1, MPI_INTEGER, &
      &     rbuf, 1, MPI_INTEGER, root, comm, ierr )
#else
    rbuf(1)=sval
#endif
  end subroutine hecmw_gather_int_1

  subroutine hecmw_allgather_int_1(sval, rbuf, comm)
    use hecmw_util
    implicit none
    integer(kind=kint) :: sval !send buffer
    integer(kind=kint) :: rbuf(*) !receive buffer
    integer(kind=kint) :: comm
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_allgather( sval, 1, MPI_INTEGER, &
      &     rbuf, 1, MPI_INTEGER, comm, ierr )
#else
    rbuf(1)=sval
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_scatterv( sbuf, scs, disp, MPI_REAL8, &
      &     rbuf, rc, MPI_REAL8, root, comm, ierr )
#else
    rbuf(1:rc)=sbuf(1:rc)
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_gatherv( sbuf, sc, MPI_REAL8, &
      &     rbuf, rcs, disp, MPI_REAL8, root, comm, ierr )
#else
    rbuf(1:sc)=sbuf(1:sc)
#endif
  end subroutine hecmw_gatherv_real

  subroutine hecmw_gatherv_int(sbuf, sc, &
      &     rbuf, rcs, disp, root, comm)
    use hecmw_util
    implicit none
    integer(kind=kint) :: sbuf(*) !send buffer
    integer(kind=kint) :: sc  !send counts
    integer(kind=kint) :: rbuf(*) !receive buffer
    integer(kind=kint) :: rcs(*) !receive counts
    integer(kind=kint) :: disp(*) !displacement
    integer(kind=kint) :: root
    integer(kind=kint) :: comm
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_gatherv( sbuf, sc, MPI_INTEGER, &
      &     rbuf, rcs, disp, MPI_INTEGER, root, comm, ierr )
#else
    rbuf(1:sc)=sbuf(1:sc)
#endif
  end subroutine hecmw_gatherv_int

  subroutine hecmw_allreduce_int_1(sval, rval, op, comm)
    use hecmw_util
    implicit none
    integer(kind=kint)   :: sval !send val
    integer(kind=kint)   :: rval !receive val
    integer(kind=kint):: op, comm, ierr
#ifndef HECMW_SERIAL
    call MPI_ALLREDUCE(sval, rval, 1, MPI_INTEGER, &
      &     hecmw_operation_hec2mpi(op), comm, ierr)
#else
    rval=sval
#endif
  end subroutine hecmw_allreduce_int_1

  subroutine hecmw_alltoall_int(sbuf, sc, rbuf, rc, comm)
    use hecmw_util
    implicit none
    integer(kind=kint) :: sbuf(*)
    integer(kind=kint) :: sc
    integer(kind=kint) :: rbuf(*)
    integer(kind=kint) :: rc
    integer(kind=kint) :: comm
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_alltoall( sbuf, sc, MPI_INTEGER, &
      &     rbuf, rc, MPI_INTEGER, comm, ierr )
#else
    rbuf(1:sc)=sbuf(1:sc)
#endif
  end subroutine hecmw_alltoall_int

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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_ISEND(sbuf, sc, MPI_INTEGER, &
      &     dest, tag, comm, req, ierr)
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_ISEND(sbuf, sc, MPI_DOUBLE_PRECISION, &
      &     dest, tag, comm, req, ierr)
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_IRECV(rbuf, rc, MPI_INTEGER, &
      &     source, tag, comm, req, ierr)
#endif
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
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_IRECV(rbuf, rc, MPI_DOUBLE_PRECISION, &
      &     source, tag, comm, req, ierr)
#endif
  end subroutine hecmw_irecv_r

  subroutine hecmw_waitall(cnt, reqs, stats)
    use hecmw_util
    implicit none
    integer(kind=kint) :: cnt
    integer(kind=kint) :: reqs(*)
    integer(kind=kint) :: stats(HECMW_STATUS_SIZE,*)
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_WAITALL(cnt, reqs, stats, ierr)
#endif
  end subroutine hecmw_waitall

  subroutine hecmw_recv_int(rbuf, rc, source, &
      &     tag, comm, stat)
    use hecmw_util
    implicit none
    integer(kind=kint) :: rbuf(*)
    integer(kind=kint) :: rc
    integer(kind=kint) :: source
    integer(kind=kint) :: tag
    integer(kind=kint) :: comm
    integer(kind=kint) :: stat(HECMW_STATUS_SIZE)
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_RECV(rbuf, rc, MPI_INTEGER, &
      &     source, tag, comm, stat, ierr)
#endif
  end subroutine hecmw_recv_int

  subroutine hecmw_recv_r(rbuf, rc, source, &
      &     tag, comm, stat)
    use hecmw_util
    implicit none
    integer(kind=kint) :: rc
    double precision, dimension(rc) :: rbuf
    integer(kind=kint) :: source
    integer(kind=kint) :: tag
    integer(kind=kint) :: comm
    integer(kind=kint) :: stat(HECMW_STATUS_SIZE)
#ifndef HECMW_SERIAL
    integer(kind=kint) :: ierr
    call MPI_RECV(rbuf, rc, MPI_DOUBLE_PRECISION, &
      &     source, tag, comm, stat, ierr)
#endif
  end subroutine hecmw_recv_r
  !C
  !C***
  !C*** hecmw_allREDUCE
  !C***
  !C
  subroutine hecmw_allreduce_DP(val,VALM,n,hec_op,comm )
    use hecmw_util
    implicit none
    integer(kind=kint) :: n, hec_op,op, comm, ierr
    double precision, dimension(n) :: val
    double precision, dimension(n) :: VALM
#ifndef HECMW_SERIAL
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
    call MPI_allREDUCE(val,VALM,n,MPI_DOUBLE_PRECISION,op,comm,ierr)
#else
    integer(kind=kint) :: i
    do i=1,n
      VALM(i) = val(i)
    end do
#endif
  end subroutine hecmw_allREDUCE_DP

  subroutine hecmw_allreduce_DP1(s1,s2,hec_op,comm )
    use hecmw_util
    implicit none
    integer(kind=kint) ::  hec_op, comm
    double precision :: s1, s2
    double precision, dimension(1) :: val
    double precision, dimension(1) :: VALM
#ifndef HECMW_SERIAL
    val(1) = s1
    VALM(1) = s2
    call hecmw_allreduce_DP(val,VALM,1,hec_op,comm )
    s1 = val(1)
    s2 = VALM(1)
#else
    s2 = s1
#endif
  end subroutine hecmw_allreduce_DP1
  !C
  !C***
  !C*** hecmw_allREDUCE_R
  !C***
  !C
  subroutine hecmw_allreduce_R (hecMESH, val, n, ntag)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, ntag
    real(kind=kreal), dimension(n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    real(kind=kreal), dimension(:), allocatable :: VALM

    allocate (VALM(n))
    VALM= 0.d0
    if (ntag .eq. hecmw_sum) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_DOUBLE_PRECISION, MPI_SUM,              &
        &        hecMESH%MPI_COMM, ierr)
    endif

    if (ntag .eq. hecmw_max) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_DOUBLE_PRECISION, MPI_MAX,              &
        &        hecMESH%MPI_COMM, ierr)
    endif

    if (ntag .eq. hecmw_min) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_DOUBLE_PRECISION, MPI_MIN,              &
        &        hecMESH%MPI_COMM, ierr)
    endif

    val= VALM
    deallocate (VALM)
#endif
  end subroutine hecmw_allreduce_R

  subroutine hecmw_allreduce_R1 (hecMESH, s, ntag)
    use hecmw_util
    implicit none
    integer(kind=kint):: ntag
    real(kind=kreal) :: s
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    real(kind=kreal), dimension(1) :: val
    val(1) = s
    call hecmw_allreduce_R(hecMESH, val, 1, ntag )
    s = val(1)
#endif
  end subroutine hecmw_allreduce_R1

  !C
  !C***
  !C*** hecmw_allREDUCE_I
  !C***
  !C
  subroutine hecmw_allreduce_I(hecMESH, val, n, ntag)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, ntag
    integer(kind=kint), dimension(n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    integer(kind=kint), dimension(:), allocatable :: VALM

    allocate (VALM(n))
    VALM= 0
    if (ntag .eq. hecmw_sum) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_INTEGER, MPI_SUM,                       &
        &        hecMESH%MPI_COMM, ierr)
    endif

    if (ntag .eq. hecmw_max) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_INTEGER, MPI_MAX,                       &
        &        hecMESH%MPI_COMM, ierr)
    endif

    if (ntag .eq. hecmw_min) then
      call MPI_allREDUCE                                              &
        &       (val, VALM, n, MPI_INTEGER, MPI_MIN,                       &
        &        hecMESH%MPI_COMM, ierr)
    endif


    val= VALM
    deallocate (VALM)
#endif
  end subroutine hecmw_allreduce_I

  subroutine hecmw_allreduce_I1 (hecMESH, s, ntag)
    use hecmw_util
    implicit none
    integer(kind=kint)::  ntag, s
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint), dimension(1) :: val

    val(1) = s
    call hecmw_allreduce_I(hecMESH, val, 1, ntag )
    s = val(1)
#endif
  end subroutine hecmw_allreduce_I1

  !C
  !C***
  !C*** hecmw_bcast_R
  !C***
  !C
  subroutine hecmw_bcast_R (hecMESH, val, n, nbase)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nbase
    real(kind=kreal), dimension(n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n, MPI_DOUBLE_PRECISION, nbase, hecMESH%MPI_COMM, ierr)
#endif
  end subroutine hecmw_bcast_R

  subroutine hecmw_bcast_R_comm (val, n, nbase, comm)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nbase
    real(kind=kreal), dimension(n) :: val
    integer(kind=kint):: comm
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n, MPI_DOUBLE_PRECISION, nbase, comm, ierr)
#endif
  end subroutine hecmw_bcast_R_comm

  subroutine hecmw_bcast_R1 (hecMESH, s, nbase)
    use hecmw_util
    implicit none
    integer(kind=kint):: nbase, ierr
    real(kind=kreal) :: s
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    real(kind=kreal), dimension(1) :: val
    val(1)=s
    call MPI_BCAST (val, 1, MPI_DOUBLE_PRECISION, nbase, hecMESH%MPI_COMM, ierr)
    s = val(1)
#endif
  end subroutine hecmw_bcast_R1

  subroutine hecmw_bcast_R1_comm (s, nbase, comm)
    use hecmw_util
    implicit none
    integer(kind=kint):: nbase
    real(kind=kreal) :: s
    integer(kind=kint):: comm
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    real(kind=kreal), dimension(1) :: val
    val(1)=s
    call MPI_BCAST (val, 1, MPI_DOUBLE_PRECISION, nbase, comm, ierr)
    s = val(1)
#endif
  end subroutine hecmw_bcast_R1_comm
  !C
  !C***
  !C*** hecmw_bcast_I
  !C***
  !C
  subroutine hecmw_bcast_I (hecMESH, val, n, nbase)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nbase
    integer(kind=kint), dimension(n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n, MPI_INTEGER, nbase, hecMESH%MPI_COMM, ierr)
#endif
  end subroutine hecmw_bcast_I

  subroutine hecmw_bcast_I_comm (val, n, nbase, comm)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nbase
    integer(kind=kint), dimension(n) :: val
    integer(kind=kint):: comm
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n, MPI_INTEGER, nbase, comm, ierr)
#endif
  end subroutine hecmw_bcast_I_comm

  subroutine hecmw_bcast_I1 (hecMESH, s, nbase)
    use hecmw_util
    implicit none
    integer(kind=kint):: nbase, s
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    integer(kind=kint), dimension(1) :: val
    val(1) = s
    call MPI_BCAST (val, 1, MPI_INTEGER, nbase, hecMESH%MPI_COMM, ierr)
    s = val(1)
#endif
  end subroutine hecmw_bcast_I1

  subroutine hecmw_bcast_I1_comm (s, nbase, comm)
    use hecmw_util
    implicit none
    integer(kind=kint):: nbase, s
    integer(kind=kint):: comm
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    integer(kind=kint), dimension(1) :: val
    val(1) = s
    call MPI_BCAST (val, 1, MPI_INTEGER, nbase, comm, ierr)
    s = val(1)
#endif
  end subroutine hecmw_bcast_I1_comm
  !C
  !C***
  !C*** hecmw_bcast_C
  !C***
  !C
  subroutine hecmw_bcast_C (hecMESH, val, n, nn, nbase)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nn, nbase
    character(len=n) :: val(nn)
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n*nn, MPI_CHARACTER, nbase, hecMESH%MPI_COMM,&
      &                                                 ierr)
#endif
  end subroutine hecmw_bcast_C

  subroutine hecmw_bcast_C_comm (val, n, nn, nbase, comm)
    use hecmw_util
    implicit none
    integer(kind=kint):: n, nn, nbase
    character(len=n) :: val(nn)
    integer(kind=kint):: comm
#ifndef HECMW_SERIAL
    integer(kind=kint):: ierr
    call MPI_BCAST (val, n*nn, MPI_CHARACTER, nbase, comm,&
      &                                                 ierr)
#endif
  end subroutine hecmw_bcast_C_comm

  !C
  !C***
  !C*** hecmw_assemble_R
  !C***
  subroutine hecmw_assemble_R (hecMESH, val, n, m)
    use hecmw_util
    use  hecmw_solver_SR

    implicit none
    integer(kind=kint):: n, m
    real(kind=kreal), dimension(m*n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ns, nr
    real(kind=kreal), dimension(:), allocatable :: WS, WR

    if( hecMESH%n_neighbor_pe == 0 ) return

    ns = hecMESH%import_index(hecMESH%n_neighbor_pe)
    nr = hecMESH%export_index(hecMESH%n_neighbor_pe)

    allocate (WS(m*ns), WR(m*nr))
    call hecmw_solve_REV_SEND_RECV                                    &
      &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
      &     hecMESH%import_index, hecMESH%import_item,                   &
      &     hecMESH%export_index, hecMESH%export_item,                   &
      &     WS, WR, val , hecMESH%MPI_COMM, hecMESH%my_rank)
    deallocate (WS, WR)
#endif
  end subroutine hecmw_assemble_R

  !C
  !C***
  !C*** hecmw_assemble_I
  !C***
  subroutine hecmw_assemble_I (hecMESH, val, n, m)
    use hecmw_util
    use  hecmw_solver_SR_i

    implicit none
    integer(kind=kint):: n, m
    integer(kind=kint), dimension(m*n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ns, nr
    integer(kind=kint), dimension(:), allocatable :: WS, WR

    if( hecMESH%n_neighbor_pe == 0 ) return

    ns = hecMESH%import_index(hecMESH%n_neighbor_pe)
    nr = hecMESH%export_index(hecMESH%n_neighbor_pe)

    allocate (WS(m*ns), WR(m*nr))
    call hecmw_solve_REV_SEND_RECV_i                                    &
      &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
      &     hecMESH%import_index, hecMESH%import_item,                   &
      &     hecMESH%export_index, hecMESH%export_item,                   &
      &     WS, WR, val , hecMESH%MPI_COMM, hecMESH%my_rank)
    deallocate (WS, WR)
#endif
  end subroutine hecmw_assemble_I

  !C
  !C***
  !C*** hecmw_update_R
  !C***
  !C
  subroutine hecmw_update_R (hecMESH, val, n, m)
    use hecmw_util
    use  hecmw_solver_SR

    implicit none
    integer(kind=kint):: n, m
    real(kind=kreal), dimension(m*n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ns, nr
    real(kind=kreal), dimension(:), allocatable :: WS, WR

    if( hecMESH%n_neighbor_pe == 0 ) return

    ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
    nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

    allocate (WS(m*ns), WR(m*nr))
    call hecmw_solve_SEND_RECV                                     &
      &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
      &     hecMESH%import_index, hecMESH%import_item,                   &
      &     hecMESH%export_index, hecMESH%export_item,                   &
      &     WS, WR, val , hecMESH%MPI_COMM, hecMESH%my_rank)
    deallocate (WS, WR)
#endif
  end subroutine hecmw_update_R

  !C
  !C***
  !C*** hecmw_update_R_async
  !C***
  !C
  !C    REAL
  !C
  subroutine hecmw_update_R_async (hecMESH, val, n, m, ireq)
    use hecmw_util
    use  hecmw_solver_SR
    implicit none
    integer(kind=kint):: n, m, ireq
    real(kind=kreal), dimension(m*n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    if( hecMESH%n_neighbor_pe == 0 ) return

    call hecmw_solve_ISEND_IRECV                                   &
      &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
      &     hecMESH%import_index, hecMESH%import_item,                   &
      &     hecMESH%export_index, hecMESH%export_item,                   &
      &     val , hecMESH%MPI_COMM, hecMESH%my_rank, ireq)
#endif
  end subroutine hecmw_update_R_async

  !C
  !C***
  !C*** hecmw_update_R_wait
  !C***
  !C
  !C    REAL
  !C
  subroutine hecmw_update_R_wait (hecMESH, ireq)
    use hecmw_util
    use  hecmw_solver_SR
    implicit none
    integer(kind=kint):: ireq
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    if( hecMESH%n_neighbor_pe == 0 ) return

    call hecmw_solve_ISEND_IRECV_WAIT(hecMESH%n_dof, ireq )
#endif
  end subroutine hecmw_update_R_wait


  !C
  !C***
  !C*** hecmw_update_I
  !C***
  !C
  subroutine hecmw_update_I (hecMESH, val, n, m)
    use hecmw_util
    use  hecmw_solver_SR_i

    implicit none
    integer(kind=kint):: n, m
    integer(kind=kint), dimension(m*n) :: val
    type (hecmwST_local_mesh) :: hecMESH
#ifndef HECMW_SERIAL
    integer(kind=kint):: ns, nr
    integer(kind=kint), dimension(:), allocatable :: WS, WR

    if( hecMESH%n_neighbor_pe == 0 ) return

    ns = hecMESH%export_index(hecMESH%n_neighbor_pe)
    nr = hecMESH%import_index(hecMESH%n_neighbor_pe)

    allocate (WS(m*ns), WR(m*nr))
    call hecmw_solve_SEND_RECV_i                                    &
      &   ( n, m, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,            &
      &     hecMESH%import_index, hecMESH%import_item,                   &
      &     hecMESH%export_index, hecMESH%export_item,                   &
      &     WS, WR, val , hecMESH%MPI_COMM, hecMESH%my_rank)
    deallocate (WS, WR)
#endif
  end subroutine hecmw_update_I

end module m_hecmw_comm_f
