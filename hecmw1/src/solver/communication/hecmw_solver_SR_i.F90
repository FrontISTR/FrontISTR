!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_SR_i
!C***
!C
module hecmw_solver_SR_i
contains
  !C
  !C*** SOLVER_SEND_RECV
  !C
  subroutine  HECMW_SOLVE_SEND_RECV_i                             &
      &                ( N, m, NEIBPETOT, NEIBPE,                        &
      &                  STACK_IMPORT, NOD_IMPORT,                       &
      &                  STACK_EXPORT, NOD_EXPORT, WS, WR, X,            &
      &                  SOLVER_COMM,my_rank)

    use hecmw_util
    implicit none
    !      include  'mpif.h'
    !      include  'hecmw_config_f.h'

    integer(kind=kint )                , intent(in)   ::  N, m
    integer(kind=kint )                , intent(in)   ::  NEIBPETOT
    integer(kind=kint ), pointer :: NEIBPE      (:)
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    integer(kind=kint ), pointer :: STACK_EXPORT(:)
    integer(kind=kint ), pointer :: NOD_EXPORT  (:)
    integer(kind=kint ), dimension(:  ), intent(inout):: WS
    integer(kind=kint ), dimension(:  ), intent(inout):: WR
    integer(kind=kint ), dimension(:  ), intent(inout):: X
    integer(kind=kint )                , intent(in)   ::SOLVER_COMM
    integer(kind=kint )                , intent(in)   :: my_rank

#ifndef HECMW_SERIAL
    integer(kind=kint ), dimension(:,:), allocatable :: sta1
    integer(kind=kint ), dimension(:,:), allocatable :: sta2
    integer(kind=kint ), dimension(:  ), allocatable :: req1
    integer(kind=kint ), dimension(:  ), allocatable :: req2

    integer(kind=kint ), save :: NFLAG
    data NFLAG/0/
    ! local valiables
    integer(kind=kint ) :: neib,istart,inum,k,ii,kk,ierr,nreq1,nreq2
    !C
    !C-- INIT.
    allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
    allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
    allocate (req1(NEIBPETOT))
    allocate (req2(NEIBPETOT))

    !C
    !C-- SEND
    nreq1=0
    do neib= 1, NEIBPETOT
      istart= STACK_EXPORT(neib-1)
      inum  = STACK_EXPORT(neib  ) - istart
      if (inum==0) cycle
      nreq1=nreq1+1
      do k= istart+1, istart+inum
        ii   = m*NOD_EXPORT(k)
        do kk= 1, m
          WS(m*k-kk+1)= X(ii-kk+1)
        enddo
      enddo

      call MPI_ISEND (WS(m*istart+1), m*inum, MPI_INTEGER,            &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req1(nreq1), ierr)
    enddo

    !C
    !C-- RECEIVE
    nreq2=0
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      if (inum==0) cycle
      nreq2=nreq2+1
      call MPI_IRECV (WR(m*istart+1), m*inum, MPI_INTEGER,            &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req2(nreq2), ierr)
    enddo

    call MPI_WAITALL (nreq2, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = m*NOD_IMPORT(k)
        do kk= 1, m
          X(ii-kk+1)= WR(m*k-kk+1)
        enddo
      enddo
    enddo

    call MPI_WAITALL (nreq1, req1, sta1, ierr)

    deallocate (sta1, sta2, req1, req2)
#endif
  end subroutine hecmw_solve_send_recv_i

  !C
  !C*** SOLVER_REVERSE_SEND_RECV_i
  !C
  subroutine  HECMW_SOLVE_REV_SEND_RECV_i                               &
      &                ( N, M, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  WS, WR, X, SOLVER_COMM,my_rank)

    use hecmw_util
    implicit none

    integer(kind=kint )                , intent(in)   ::  N, M
    integer(kind=kint )                , intent(in)   ::  NEIBPETOT
    integer(kind=kint ), pointer :: NEIBPE      (:)
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    integer(kind=kint ), pointer :: STACK_EXPORT(:)
    integer(kind=kint ), pointer :: NOD_EXPORT  (:)
    integer(kind=kint ), dimension(:  ), intent(inout):: WS
    integer(kind=kint ), dimension(:  ), intent(inout):: WR
    integer(kind=kint ), dimension(:  ), intent(inout):: X
    integer(kind=kint )                , intent(in)   ::SOLVER_COMM
    integer(kind=kint )                , intent(in)   :: my_rank

#ifndef HECMW_SERIAL
    integer(kind=kint ), dimension(:,:), allocatable :: sta1
    integer(kind=kint ), dimension(:,:), allocatable :: sta2
    integer(kind=kint ), dimension(:  ), allocatable :: req1
    integer(kind=kint ), dimension(:  ), allocatable :: req2

    ! local valiables
    integer(kind=kint ) :: neib,istart,inum,k,kk,ii,ierr,nreq1,nreq2
    !C
    !C-- INIT.
    allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
    allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
    allocate (req1(NEIBPETOT))
    allocate (req2(NEIBPETOT))

    !C
    !C-- SEND
    nreq1=0
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      if (inum==0) cycle
      nreq1=nreq1+1
      do k= istart+1, istart+inum
        ii   = NOD_IMPORT(k)
        do kk = 1, M
          WS(M*(k-1)+kk)= X(M*(ii-1)+kk)
        enddo
      enddo

      call MPI_ISEND (WS(M*istart+1), M*inum,MPI_INTEGER,    &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req1(nreq1), ierr)
    enddo

    !C
    !C-- RECEIVE
    nreq2=0
    do neib= 1, NEIBPETOT
      istart= STACK_EXPORT(neib-1)
      inum  = STACK_EXPORT(neib  ) - istart
      if (inum==0) cycle
      nreq2=nreq2+1
      call MPI_IRECV (WR(M*istart+1), M*inum, MPI_INTEGER,   &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req2(nreq2), ierr)
    enddo

    call MPI_WAITALL (nreq2, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= STACK_EXPORT(neib-1)
      inum  = STACK_EXPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = NOD_EXPORT(k)
        do kk = 1, M
          X(M*(ii-1)+kk)= X(M*(ii-1)+kk)+WR(M*(k-1)+kk)
        enddo
      enddo
    enddo

    call MPI_WAITALL (nreq1, req1, sta1, ierr)
    deallocate (sta1, sta2, req1, req2)
#endif
  end subroutine HECMW_SOLVE_REV_SEND_RECV_i

end module     hecmw_solver_SR_i
