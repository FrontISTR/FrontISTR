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

      module hecmw_adapt_STACK_SR
      contains
!C
!C***
!C*** hecmw_adapt_STACK_SEND_RECV
!C***
!C
!C    exchange IMPORT/EXPORT item NUMBER information
!C
      subroutine  hecmw_adapt_STACK_SEND_RECV                           &
     &              ( NEIBPETOT, NEIBPE, STACK_IMPORT, STACK_EXPORT,    &
     &                SOLVER_COMM,my_rank)

      use hecmw_util
      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint ), intent(in) :: NEIBPETOT
      integer(kind=kint ), pointer    :: NEIBPE      (:)
      integer(kind=kint ), pointer    :: STACK_IMPORT(:)
      integer(kind=kint ), pointer    :: STACK_EXPORT(:)
      integer(kind=kint )             :: SOLVER_COMM, my_rank

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
        num= STACK_IMPORT(neib)
        call MPI_ISEND (num, 1, MPI_INTEGER, NEIBPE(neib), 0,           &
     &                  SOLVER_COMM, req1(neib), ierr)
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        call MPI_IRECV (STACK_EXPORT(neib), 1, MPI_INTEGER,             &
     &                  NEIBPE(neib), 0, SOLVER_COMM, req2(neib), ierr)
      enddo

      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
      call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

      end subroutine hecmw_adapt_STACK_SEND_RECV
      end module     hecmw_adapt_STACK_SR



