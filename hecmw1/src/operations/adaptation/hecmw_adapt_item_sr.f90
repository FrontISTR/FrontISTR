!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

      module hecmw_adapt_ITEM_SR
      contains
!C
!C***
!C*** hecmw_adapt_ITEM_SEND_RECV
!C***
!C
!C    exchange IMPORT/EXPORT item information
!C    form communication table
!C
      subroutine  hecmw_adapt_ITEM_SEND_RECV                            &
     &                ( N, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT, &
     &                                        STACK_EXPORT, NOD_EXPORT, &
     &                  WS, WR, SOLVER_COMM, my_rank, NTAB)

      use hecmw_util

      integer(kind=kint ), intent(in)   :: N, NEIBPETOT
      integer(kind=kint ), pointer      :: NEIBPE      (:)
      integer(kind=kint ), pointer      :: STACK_IMPORT(:)
      integer(kind=kint ), pointer      :: STACK_EXPORT(:)
      integer(kind=kint ), pointer      :: NOD_IMPORT(:), NOD_EXPORT(:)

      integer(kind=kint ), dimension(N*NTAB) :: WS, WR
      integer(kind=kint ),                   :: SOLVER_COMM, my_rank

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
        istart= STACK_IMPORT(neib-1)
        inum  = STACK_IMPORT(neib  ) - istart
        do k= istart+1, istart+inum
          is= (k-1)*NTAB
          do jj= 1, NTAB
            WS(is+jj)= NOD_IMPORT(is+jj)
          enddo
        enddo
        is= istart*NTAB
        call MPI_ISEND (WS(is+1), NTAB*inum, MPI_INTEGER,               &
     &                  NEIBPE(neib), 0, SOLVER_COMM,                   &
     &                  req1(neib), ierr)
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        istart= STACK_EXPORT(neib-1)
        inum  = STACK_EXPORT(neib  ) - istart
        is= istart*NTAB
        call MPI_IRECV (WR(is+1), NTAB*inum, MPI_INTEGER,               &
     &                  NEIBPE(neib), 0, SOLVER_COMM,                   &
     &                  req2(neib), ierr)
      enddo

      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
   
      do neib= 1, NEIBPETOT
        istart= STACK_EXPORT(neib-1)
        inum  = STACK_EXPORT(neib  ) - istart
      do k= istart+1, istart+inum
        is= (k-1)*NTAB
        do jj= 1, NTAB
          NOD_EXPORT(is+jj)= WR(is+jj)
        enddo
      enddo
      enddo

      call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

      end subroutine hecmw_adapt_ITEM_SEND_RECV
      end module     hecmw_adapt_ITEM_SR
      


