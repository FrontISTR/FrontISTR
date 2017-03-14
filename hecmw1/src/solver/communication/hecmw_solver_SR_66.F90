!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_SR_66
!C***
!C
      module hecmw_solver_SR_66
      contains
!C
!C*** SOLVER_SEND_RECV
!C
      subroutine  HECMW_SOLVE_SEND_RECV_66                              &
     &                ( N, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
     &                                        STACK_EXPORT, NOD_EXPORT, &
     &                  WS, WR, X, SOLVER_COMM,my_rank)

      use hecmw_util
      implicit none
!      include  'mpif.h'
!      include  'hecmw_config_f.h'

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

      integer(kind=kint ), save :: NFLAG
      data NFLAG/0/

      ! local valiables
      integer(kind=kint ) :: neib,istart,inum,k,ii,ierr
!C
!C-- INIT.
      allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
      allocate (req1(NEIBPETOT))
      allocate (req2(NEIBPETOT))

!C
!C-- SEND
      do neib= 1, NEIBPETOT
        istart= STACK_EXPORT(neib-1)
        inum  = STACK_EXPORT(neib  ) - istart
        do k= istart+1, istart+inum
               ii   = 6*NOD_EXPORT(k)
           WS(6*k-5)= X(ii-5)
           WS(6*k-4)= X(ii-4)
           WS(6*k-3)= X(ii-3)
           WS(6*k-2)= X(ii-2)
           WS(6*k-1)= X(ii-1)
           WS(6*k  )= X(ii  )
        enddo

        call MPI_ISEND (WS(6*istart+1), 6*inum,MPI_DOUBLE_PRECISION,    &
     &                  NEIBPE(neib), 0, SOLVER_COMM, req1(neib), ierr)
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        istart= STACK_IMPORT(neib-1)
        inum  = STACK_IMPORT(neib  ) - istart
        call MPI_IRECV (WR(6*istart+1), 6*inum, MPI_DOUBLE_PRECISION,   &
     &                  NEIBPE(neib), 0, SOLVER_COMM, req2(neib), ierr)
      enddo

      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)

      do neib= 1, NEIBPETOT
        istart= STACK_IMPORT(neib-1)
        inum  = STACK_IMPORT(neib  ) - istart
      do k= istart+1, istart+inum
          ii   = 6*NOD_IMPORT(k)
        X(ii-5)= WR(6*k-5)
        X(ii-4)= WR(6*k-4)
        X(ii-3)= WR(6*k-3)
        X(ii-2)= WR(6*k-2)
        X(ii-1)= WR(6*k-1)
        X(ii  )= WR(6*k  )
      enddo
      enddo

      call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)
      deallocate (sta1, sta2, req1, req2)
#endif
      end subroutine hecmw_solve_send_recv_66
      end module     hecmw_solver_SR_66
