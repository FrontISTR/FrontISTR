module hecmw_solver_misc

      contains

!C
!C***
!C*** hecmw_innerProduct_33I
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
      real(kind=kreal) :: COMMtime

      integer(kind=kint) :: i
      real(kind=kreal) :: START_TIME, END_TIME

      sum = 0
      do i = 1, hecMESH%nn_internal * ndof
         sum = sum + X(i)*Y(i)
      end do

      START_TIME= HECMW_WTIME()
      call hecmw_allreduce_I1 (hecMESH, sum, hecmw_sum)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

      end subroutine hecmw_innerProduct_I

!C
!C***
!C*** hecmw_innerProduct_33R
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
      real(kind=kreal) :: COMMtime

      integer(kind=kint) :: i
      real(kind=kreal) :: START_TIME, END_TIME

      sum = 0.0d0
      do i = 1, hecMESH%nn_internal * ndof
         sum = sum + X(i)*Y(i)
      end do

      START_TIME= HECMW_WTIME()
      call hecmw_allreduce_R1 (hecMESH, sum, hecmw_sum)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

      end subroutine hecmw_innerProduct_R

end module hecmw_solver_misc