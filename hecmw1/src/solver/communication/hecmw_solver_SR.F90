!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_SR_mm
!C***
!C
module hecmw_solver_SR
  use hecmw_util

#ifndef HECMW_SERIAL
  type async_buf
    integer(kind=kint )   ::  NEIBPETOT = 0
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    real   (kind=kreal), pointer :: WS(:)
    real   (kind=kreal), pointer :: WR(:)
    real   (kind=kreal), pointer :: X(:)
    type(MPI_REQUEST), pointer :: req1(:)
    type(MPI_REQUEST), pointer :: req2(:)
    integer(kind=kint )   ::  nreq1
    integer(kind=kint )   ::  nreq2
  end type async_buf

  integer(kind=kint), parameter :: MAX_NREQ = 8
  type(async_buf), save :: abuf(MAX_NREQ)
#endif

contains
  !C
  !C*** SOLVER_SEND_RECV
  !C
  subroutine  HECMW_SOLVE_SEND_RECV                              &
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
    real   (kind=kreal), dimension(:  ), intent(inout):: WS
    real   (kind=kreal), dimension(:  ), intent(inout):: WR
    real   (kind=kreal), dimension(:  ), intent(inout):: X
    integer(kind=kint )                , intent(in)   ::SOLVER_COMM
    integer(kind=kint )                , intent(in)   :: my_rank

#ifndef HECMW_SERIAL
    type(MPI_STATUS), dimension(:), allocatable :: sta1
    type(MPI_STATUS), dimension(:), allocatable :: sta2
    type(MPI_REQUEST), dimension(:), allocatable :: req1
    type(MPI_REQUEST), dimension(:), allocatable :: req2

    integer(kind=kint ), save :: NFLAG
    data NFLAG/0/
    ! local valiables
    integer(kind=kint ) :: neib,istart,inum,k,kk,ii,ierr,nreq1,nreq2
    type(MPI_COMM) :: mpicomm
    mpicomm%mpi_val = SOLVER_COMM
    !C
    !C-- INIT.
    allocate (sta1(NEIBPETOT))
    allocate (sta2(NEIBPETOT))
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
        ii   = NOD_EXPORT(k)
        do kk = 1, m
          WS(m*(k-1)+kk)= X(m*(ii-1)+kk)
        enddo
      enddo

      call MPI_ISEND (WS(m*istart+1), m*inum,MPI_DOUBLE_PRECISION,    &
        &                  NEIBPE(neib), 0, mpicomm, req1(nreq1), ierr)
    enddo

    !C
    !C-- RECEIVE
    nreq2=0
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      if (inum==0) cycle
      nreq2=nreq2+1
      call MPI_IRECV (WR(m*istart+1), m*inum, MPI_DOUBLE_PRECISION,   &
        &                  NEIBPE(neib), 0, mpicomm, req2(nreq2), ierr)
    enddo

    call MPI_WAITALL (nreq2, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = NOD_IMPORT(k)
        do kk= 1, m
          X(m*(ii-1)+kk)= WR(m*(k-1)+kk)
        enddo
      enddo
    enddo

    call MPI_WAITALL (nreq1, req1, sta1, ierr)

    deallocate (sta1, sta2, req1, req2)
#endif
  end subroutine hecmw_solve_send_recv

  !C
  !C*** SOLVER_ISEND_IRECV
  !C
  subroutine  HECMW_SOLVE_ISEND_IRECV                              &
      &                ( N, M, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  X, SOLVER_COMM,my_rank,ireq)
    implicit none
    integer(kind=kint )                , intent(in)   ::  N, M
    integer(kind=kint )                , intent(in)   ::  NEIBPETOT
    integer(kind=kint ), pointer :: NEIBPE      (:)
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    integer(kind=kint ), pointer :: STACK_EXPORT(:)
    integer(kind=kint ), pointer :: NOD_EXPORT  (:)
    real   (kind=kreal), target, intent(inout):: X (:)
    integer(kind=kint )        , intent(in)   ::SOLVER_COMM
    integer(kind=kint )        , intent(in)   :: my_rank
    integer(kind=kint )        , intent(out)  :: ireq

#ifndef HECMW_SERIAL
    ! local valiables
    real   (kind=kreal), pointer :: WS(:)
    real   (kind=kreal), pointer :: WR(:)
    type(MPI_REQUEST), pointer :: req1(:)
    type(MPI_REQUEST), pointer :: req2(:)
    integer(kind=kint ) :: neib,istart,inum,k,kk,ii,ierr,i,nreq1,nreq2
    type(MPI_COMM) :: mpicomm
    mpicomm%mpi_val = SOLVER_COMM
    !C
    !C-- INIT.
    allocate (WS(M*STACK_EXPORT(NEIBPETOT)))
    allocate (WR(M*STACK_IMPORT(NEIBPETOT)))
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
        ii   = NOD_EXPORT(k)
        do kk = 1, m
          WS(m*(k-1)+kk)= X(m*(ii-1)+kk)
        enddo
      enddo
      call MPI_ISEND (WS(m*istart+1), m*inum,MPI_DOUBLE_PRECISION,    &
        &                  NEIBPE(neib), 0, mpicomm, req1(nreq1), ierr)
    enddo
    !C
    !C-- RECEIVE
    nreq2=0
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      if (inum==0) cycle
      nreq2=nreq2+1
      call MPI_IRECV (WR(m*istart+1), m*inum, MPI_DOUBLE_PRECISION,   &
        &                  NEIBPE(neib), 0, mpicomm, req2(nreq2), ierr)
    enddo
    !C
    !C-- Find empty abuf
    ireq = 0
    do i = 1, MAX_NREQ
      if (abuf(i)%NEIBPETOT == 0) then
        ireq = i
        exit
      endif
    enddo
    if (ireq == 0) then
      stop 'Error: HECMW_SOLVE_ISEND_IRECV: exceeded maximum num of requests'
    endif
    !C
    !C-- Save in abuf
    abuf(ireq)%NEIBPETOT   =  NEIBPETOT
    abuf(ireq)%STACK_IMPORT=> STACK_IMPORT
    abuf(ireq)%NOD_IMPORT  => NOD_IMPORT
    abuf(ireq)%WS          => WS
    abuf(ireq)%WR          => WR
    abuf(ireq)%X           => X
    abuf(ireq)%req1        => req1
    abuf(ireq)%req2        => req2
    abuf(ireq)%nreq1       =  nreq1
    abuf(ireq)%nreq2       =  nreq2
#endif
  end subroutine HECMW_SOLVE_ISEND_IRECV

  !C
  !C*** SOLVER_ISEND_IRECV_WAIT
  !C
  subroutine  HECMW_SOLVE_ISEND_IRECV_WAIT( m, ireq )
    use hecmw_util
    implicit none
    integer(kind=kint ), intent(in)   :: m, ireq

#ifndef HECMW_SERIAL
    ! local valiables
    integer(kind=kint ) ::  NEIBPETOT
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    real   (kind=kreal), pointer :: WS(:)
    real   (kind=kreal), pointer :: WR(:)
    real   (kind=kreal), pointer :: X (:)
    type(MPI_REQUEST), pointer :: req1(:)
    type(MPI_REQUEST), pointer :: req2(:)
    type(MPI_STATUS), dimension(:), allocatable :: sta1
    type(MPI_STATUS), dimension(:), allocatable :: sta2
    integer(kind=kint ) :: neib,istart,inum,k,j,ii,ierr,nreq1,nreq2
    !C-- Check ireq
    if (ireq < 0 .or. ireq > MAX_NREQ) then
      stop 'ERROR: HECMW_SOLVE_ISEND_IRECV_WAIT: invalid ireq'
    endif
    !C-- Restore from abuf
    NEIBPETOT   =  abuf(ireq)%NEIBPETOT
    STACK_IMPORT=> abuf(ireq)%STACK_IMPORT
    NOD_IMPORT  => abuf(ireq)%NOD_IMPORT
    WS          => abuf(ireq)%WS
    WR          => abuf(ireq)%WR
    X           => abuf(ireq)%X
    req1        => abuf(ireq)%req1
    req2        => abuf(ireq)%req2
    nreq1       =  abuf(ireq)%nreq1
    nreq2       =  abuf(ireq)%nreq2
    !C-- Free abuf
    abuf(ireq)%NEIBPETOT = 0
    !C-- Allocate
    allocate (sta1(NEIBPETOT))
    allocate (sta2(NEIBPETOT))
    !C-- Wait irecv
    call MPI_WAITALL (nreq2, req2, sta2, ierr)
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = NOD_IMPORT(k)
        do j = 1, m
          X(m*(ii-1)+j)= WR(m*(k-1)+j)
        enddo
      enddo
    enddo
    !C-- Wait isend
    call MPI_WAITALL (nreq1, req1, sta1, ierr)
    !C-- Deallocate
    deallocate (sta1, sta2)
    deallocate (req1, req2)
    deallocate (WS, WR)
#endif
  end subroutine HECMW_SOLVE_ISEND_IRECV_WAIT

  !C
  !C*** SOLVER_REVERSE_SEND_RECV
  !C
  subroutine  HECMW_SOLVE_REV_SEND_RECV                               &
      &                ( N, M, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  WS, WR, X, SOLVER_COMM,my_rank)

    implicit none

    integer(kind=kint )                , intent(in)   ::  N, M
    integer(kind=kint )                , intent(in)   ::  NEIBPETOT
    integer(kind=kint ), pointer :: NEIBPE      (:)
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    integer(kind=kint ), pointer :: STACK_EXPORT(:)
    integer(kind=kint ), pointer :: NOD_EXPORT  (:)
    real   (kind=kreal), dimension(:  ), intent(inout):: WS
    real   (kind=kreal), dimension(:  ), intent(inout):: WR
    real   (kind=kreal), dimension(:  ), intent(inout):: X
    integer(kind=kint )                , intent(in)   ::SOLVER_COMM
    integer(kind=kint )                , intent(in)   :: my_rank

#ifndef HECMW_SERIAL
    type(MPI_STATUS), dimension(:), allocatable :: sta1
    type(MPI_STATUS), dimension(:), allocatable :: sta2
    type(MPI_REQUEST), dimension(:), allocatable :: req1
    type(MPI_REQUEST), dimension(:), allocatable :: req2

    ! local valiables
    integer(kind=kint ) :: neib,istart,inum,k,kk,ii,ierr,nreq1,nreq2
    type(MPI_COMM) :: mpicomm
    mpicomm%mpi_val = SOLVER_COMM
    !C
    !C-- INIT.
    allocate (sta1(NEIBPETOT))
    allocate (sta2(NEIBPETOT))
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

      call MPI_ISEND (WS(M*istart+1), M*inum,MPI_DOUBLE_PRECISION,    &
        &                  NEIBPE(neib), 0, mpicomm, req1(nreq1), ierr)
    enddo

    !C
    !C-- RECEIVE
    nreq2=0
    do neib= 1, NEIBPETOT
      istart= STACK_EXPORT(neib-1)
      inum  = STACK_EXPORT(neib  ) - istart
      if (inum==0) cycle
      nreq2=nreq2+1
      call MPI_IRECV (WR(M*istart+1), M*inum, MPI_DOUBLE_PRECISION,   &
        &                  NEIBPE(neib), 0, mpicomm, req2(nreq2), ierr)
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
  end subroutine HECMW_SOLVE_REV_SEND_RECV

end module     hecmw_solver_SR
