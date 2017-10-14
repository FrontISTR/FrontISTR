!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_SR_33
!C***
!C
module hecmw_solver_SR_33
  use hecmw_util

  private

  public :: HECMW_SOLVE_SEND_RECV_33
  public :: HECMW_SOLVE_ISEND_IRECV_33
  public :: HECMW_SOLVE_ISEND_IRECV_33_WAIT

#ifndef HECMW_SERIAL
  type async_buf
    integer(kind=kint )   ::  NEIBPETOT = 0
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    real   (kind=kreal), pointer :: WS(:)
    real   (kind=kreal), pointer :: WR(:)
    real   (kind=kreal), pointer :: X(:)
    integer(kind=kint ), pointer :: req1(:)
    integer(kind=kint ), pointer :: req2(:)
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
  subroutine  HECMW_SOLVE_SEND_RECV_33                              &
      &                ( N, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  WS, WR, X, SOLVER_COMM,my_rank)

    implicit none

    integer(kind=kint )                , intent(in)   ::  N
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
    integer(kind=kint ), dimension(:,:), allocatable :: sta1
    integer(kind=kint ), dimension(:,:), allocatable :: sta2
    integer(kind=kint ), dimension(:  ), allocatable :: req1
    integer(kind=kint ), dimension(:  ), allocatable :: req2

    ! local valiables
    integer(kind=kint ) :: neib,istart,inum,k,ii,ierr,nreq1,nreq2
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
        ii   = 3*NOD_EXPORT(k)
        WS(3*k-2)= X(ii-2)
        WS(3*k-1)= X(ii-1)
        WS(3*k  )= X(ii  )
      enddo

      call MPI_ISEND (WS(3*istart+1), 3*inum,MPI_DOUBLE_PRECISION,    &
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
      call MPI_IRECV (WR(3*istart+1), 3*inum, MPI_DOUBLE_PRECISION,   &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req2(nreq2), ierr)
    enddo

    call MPI_WAITALL (nreq2, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = 3*NOD_IMPORT(k)
        X(ii-2)= WR(3*k-2)
        X(ii-1)= WR(3*k-1)
        X(ii  )= WR(3*k  )
      enddo
    enddo

    call MPI_WAITALL (nreq1, req1, sta1, ierr)
    deallocate (sta1, sta2, req1, req2)
#endif
  end subroutine hecmw_solve_send_recv_33

  !C
  !C*** SOLVER_ISEND_IRECV
  !C
  subroutine  HECMW_SOLVE_ISEND_IRECV_33                              &
      &                ( N, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  X, SOLVER_COMM,my_rank,ireq)
    implicit none
    integer(kind=kint )                , intent(in)   ::  N
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
    integer(kind=kint ), pointer :: req1(:)
    integer(kind=kint ), pointer :: req2(:)
    integer(kind=kint ) :: neib,istart,inum,k,ii,ierr,i,nreq1,nreq2
    !C
    !C-- INIT.
    allocate (WS(3*STACK_EXPORT(NEIBPETOT)))
    allocate (WR(3*STACK_IMPORT(NEIBPETOT)))
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
        ii   = 3*NOD_EXPORT(k)
        WS(3*k-2)= X(ii-2)
        WS(3*k-1)= X(ii-1)
        WS(3*k  )= X(ii  )
      enddo
      call MPI_ISEND (WS(3*istart+1), 3*inum,MPI_DOUBLE_PRECISION,    &
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
      call MPI_IRECV (WR(3*istart+1), 3*inum, MPI_DOUBLE_PRECISION,   &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req2(nreq2), ierr)
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
      stop 'Error: hecmw_solve_isend_irecv_33: exceeded maximum num of requests'
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
  end subroutine hecmw_solve_isend_irecv_33

  !C
  !C*** SOLVER_ISEND_IRECV_WAIT
  !C
  subroutine  HECMW_SOLVE_ISEND_IRECV_33_WAIT( ireq )
    implicit none
    integer(kind=kint ), intent(in)   :: ireq

#ifndef HECMW_SERIAL
    ! local valiables
    integer(kind=kint ) ::  NEIBPETOT
    integer(kind=kint ), pointer :: STACK_IMPORT(:)
    integer(kind=kint ), pointer :: NOD_IMPORT  (:)
    real   (kind=kreal), pointer :: WS(:)
    real   (kind=kreal), pointer :: WR(:)
    real   (kind=kreal), pointer :: X (:)
    integer(kind=kint ), pointer :: req1(:)
    integer(kind=kint ), pointer :: req2(:)
    integer(kind=kint ), dimension(:,:), allocatable :: sta1
    integer(kind=kint ), dimension(:,:), allocatable :: sta2
    integer(kind=kint ) :: neib,istart,inum,k,ii,ierr,nreq1,nreq2
    !C-- Check ireq
    if (ireq < 0 .or. ireq > MAX_NREQ) then
      stop 'ERROR: hecmw_solve_isend_irecv_33_wait: invalid ireq'
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
    allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
    allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
    !C-- Wait irecv
    call MPI_WAITALL (nreq2, req2, sta2, ierr)
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii   = 3*NOD_IMPORT(k)
        X(ii-2)= WR(3*k-2)
        X(ii-1)= WR(3*k-1)
        X(ii  )= WR(3*k  )
      enddo
    enddo
    !C-- Wait isend
    call MPI_WAITALL (nreq1, req1, sta1, ierr)
    !C-- Deallocate
    deallocate (sta1, sta2)
    deallocate (req1, req2)
    deallocate (WS, WR)
#endif
  end subroutine hecmw_solve_isend_irecv_33_wait

end module     hecmw_solver_SR_33
