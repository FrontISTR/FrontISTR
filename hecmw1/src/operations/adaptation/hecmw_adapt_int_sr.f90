!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

module hecmw_adapt_INT_SR
contains
  !C
  !C***
  !C*** hecmw_adapt_INT_SEND_RECV
  !C***
  !C
  subroutine  hecmw_adapt_INT_SEND_RECV                             &
      &                ( N,   NEIBPETOT,NEIBPE,STACK_IMPORT, NOD_IMPORT, &
      &                                        STACK_EXPORT, NOD_EXPORT, &
      &                  WS, WR, X, SOLVER_COMM,my_rank, NB, m)

    use hecmw_util
    implicit real*8 (A-H,O-Z)

    integer(kind=kint)                , intent(in)   ::  N, m
    integer(kind=kint)                , intent(in)   ::  NEIBPETOT
    integer(kind=kint), pointer :: NEIBPE      (:)
    integer(kind=kint), pointer :: STACK_IMPORT(:)
    integer(kind=kint), pointer :: NOD_IMPORT  (:)
    integer(kind=kint), pointer :: STACK_EXPORT(:)
    integer(kind=kint), pointer :: NOD_EXPORT  (:)
    integer(kind=kint), dimension(NB*m), intent(inout):: WS
    integer(kind=kint), dimension(NB*m), intent(inout):: WR
    integer(kind=kint), dimension(NB*N), intent(inout):: X
    integer(kind=kint)                 , intent(in)   ::SOLVER_COMM
    integer(kind=kint)                 , intent(in)   :: my_rank

    integer(kind=kint ), dimension(:,:), save, allocatable :: sta1
    integer(kind=kint ), dimension(:,:), save, allocatable :: sta2
    integer(kind=kint ), dimension(:  ), save, allocatable :: req1
    integer(kind=kint ), dimension(:  ), save, allocatable :: req2
    integer(kind=kint ), save :: NFLAG
    data NFLAG/0/

    !C
    !C-- INIT.
    if (NFLAG.eq.0) then
      allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (req1(NEIBPETOT))
      allocate (req2(NEIBPETOT))
      NFLAG= 1
    endif

    !C
    !C-- SEND
    do neib= 1, NEIBPETOT
      istart= STACK_EXPORT(neib-1)
      inum  = STACK_EXPORT(neib  ) - istart

      do k= istart+1, istart+inum
        ii= NB*NOD_EXPORT(k) - NB
        ik= NB*k             - NB
        do j= 1, NB
          WS(ik+j)= X(ii+j)
        enddo
      enddo
      call MPI_ISEND (WS(NB*istart+1), NB*inum, MPI_INTEGER,          &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req1(neib), ierr)
    enddo

    !C
    !C-- RECEIVE
    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      call MPI_IRECV (WR(NB*istart+1), NB*inum, MPI_INTEGER,          &
        &                  NEIBPE(neib), 0, SOLVER_COMM, req2(neib), ierr)
    enddo

    call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)

    do neib= 1, NEIBPETOT
      istart= STACK_IMPORT(neib-1)
      inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
        ii= NB*NOD_IMPORT(k) - NB
        ik= NB*k             - NB
        do j= 1, NB
          X(ii+j)= WR(ik+j)
        enddo
      enddo
    enddo

    call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

  end subroutine hecmw_adapt_INT_SEND_RECV
end module     hecmw_adapt_INT_SR



