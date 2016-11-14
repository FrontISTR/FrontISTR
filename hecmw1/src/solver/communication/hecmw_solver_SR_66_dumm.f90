!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!C DUMMY SUBROUTINE
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
      implicit REAL*8 (A-H,O-Z)
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


      end subroutine hecmw_solve_send_recv_66
      end module     hecmw_solver_SR_66



