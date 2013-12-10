    module hecmw_solver_misc_11

      contains

!C
!C***
!C*** hecmw_matvec_11
!C***
!C
      subroutine hecmw_matvec_11 (hecMESH, hecMAT, X, Y, n, COMMtime)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      integer(kind=kint):: n
      real(kind=kreal), dimension(n) :: X, Y
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: START_TIME, END_TIME
      integer(kind=kint) :: i, j, jS, jE, in
      real(kind=kreal) :: YV

      if (present(COMMtime)) then
      START_TIME= HECMW_WTIME()
      call hecmw_update_1_R (hecMESH, X, hecMAT%NP)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      else
      call hecmw_update_1_R (hecMESH, X, hecMAT%NP)
      endif

      do i= 1, hecMAT%N
        YV= hecMAT%D(i) * X(i)
        jS= hecMAT%indexL(i-1) + 1
        jE= hecMAT%indexL(i  )
        do j= jS, jE
          in= hecMAT%itemL(j)
          YV= YV + hecMAT%AL(j) * X(in)
        enddo
        jS= hecMAT%indexU(i-1) + 1
        jE= hecMAT%indexU(i  )
        do j= jS, jE
          in= hecMAT%itemU(j)
          YV= YV + hecMAT%AU(j) * X(in)
        enddo
        Y(i)= YV
      enddo

      end subroutine hecmw_matvec_11

!C
!C***
!C*** hecmw_matresid_11
!C***
!C
      subroutine hecmw_matresid_11 (hecMESH, hecMAT, X, B, R, COMMtime)
      use hecmw_util

      implicit none
      real(kind=kreal) :: X(:), B(:), R(:)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i

      if (present(COMMtime)) then
      call hecmw_matvec_11 (hecMESH, hecMAT, X, R, hecMAT%NP, COMMtime)
      else
      call hecmw_matvec_11 (hecMESH, hecMAT, X, R, hecMAT%NP)
      endif
      do i = 1, hecMAT%N
        R(i) = B(i) - R(i)
      enddo

      end subroutine hecmw_matresid_11

!C
!C***
!C*** hecmw_rel_resid_L2_11
!C***
!C
      function hecmw_rel_resid_L2_11 (hecMESH, hecMAT, COMMtime)
      use hecmw_util
      use hecmw_solver_misc

      implicit none
      real(kind=kreal) :: hecmw_rel_resid_L2_11
      type ( hecmwST_local_mesh ), intent(in) :: hecMESH
      type ( hecmwST_matrix     ), intent(in) :: hecMAT
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: r(hecMAT%NDOF*hecMAT%NP)
      real(kind=kreal) :: bnorm2, rnorm2

      if (present(COMMtime)) then
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, COMMtime)
      call hecmw_matresid_11(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, COMMtime)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, COMMtime)
      else
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2)
      call hecmw_matresid_11(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2)
      endif
      hecmw_rel_resid_L2_11 = sqrt(rnorm2 / bnorm2)

      end function hecmw_rel_resid_L2_11

end module hecmw_solver_misc_11
