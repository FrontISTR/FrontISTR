!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_misc

      contains

!C
!C***
!C*** hecmw_innerProduct_I
!C***
!C
      subroutine hecmw_innerProduct_I (hecMESH, ndof, X, Y, sum, COMMtime )
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: ndof
      integer(kind=kint) :: X(:), Y(:)
      integer(kind=kint)        :: sum
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i
      real(kind=kreal) :: START_TIME, END_TIME

      sum = 0
      do i = 1, hecMESH%nn_internal * ndof
         sum = sum + X(i)*Y(i)
      end do

      START_TIME= HECMW_WTIME()
      call hecmw_allreduce_I1 (hecMESH, sum, hecmw_sum)
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      end subroutine hecmw_innerProduct_I

!C
!C***
!C*** hecmw_innerProduct_R
!C***
!C
      subroutine hecmw_innerProduct_R (hecMESH, ndof, X, Y, sum, COMMtime )
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: ndof
      real(kind=kreal) :: X(:), Y(:)
      real(kind=kreal)          :: sum
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i
      real(kind=kreal) :: START_TIME, END_TIME

      sum = 0.0d0
      do i = 1, hecMESH%nn_internal * ndof
         sum = sum + X(i)*Y(i)
      end do

      START_TIME= HECMW_WTIME()
      call hecmw_allreduce_R1 (hecMESH, sum, hecmw_sum)
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      end subroutine hecmw_innerProduct_R

!C
!C***
!C*** hecmw_innerProduct_R_nocomm
!C***
!C
      subroutine hecmw_innerProduct_R_nocomm (hecMESH, ndof, X, Y, sum)
      use hecmw_util

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: ndof
      real(kind=kreal) :: X(:), Y(:)
      real(kind=kreal)          :: sum

      integer(kind=kint) :: i

      sum = 0.0d0
      do i = 1, hecMESH%nn_internal * ndof
         sum = sum + X(i)*Y(i)
      end do

      end subroutine hecmw_innerProduct_R_nocomm

!C
!C***
!C*** hecmw_time_statistics
!C***
!C
      subroutine hecmw_time_statistics (hecMESH, time, &
           t_max, t_min, t_avg, t_sd)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      type (hecmwST_local_mesh), intent(in) :: hecMESH
      real(kind=kreal), intent(in) :: time
      real(kind=kreal), intent(out) :: t_max
      real(kind=kreal), intent(out), optional :: t_min, t_avg, t_sd
      real(kind=kreal) :: t2_avg
      integer(kind=kint) :: nprocs

      nprocs = hecmw_comm_get_size()

      t_max = time
      call hecmw_allreduce_R1(hecMESH, t_max, hecmw_max)

      if (.not. present(t_min)) return
      t_min = time
      call hecmw_allreduce_R1(hecMESH, t_min, hecmw_min)

      if (.not. present(t_avg)) return
      t_avg = time
      call hecmw_allreduce_R1(hecMESH, t_avg, hecmw_sum)
      t_avg = t_avg / nprocs

      if (.not. present(t_sd)) return
      t2_avg = time*time
      call hecmw_allreduce_R1(hecMESH, t2_avg, hecmw_sum)
      t2_avg = t2_avg / nprocs

      t_sd = dsqrt(t2_avg - t_avg*t_avg)
      end subroutine hecmw_time_statistics

end module hecmw_solver_misc
