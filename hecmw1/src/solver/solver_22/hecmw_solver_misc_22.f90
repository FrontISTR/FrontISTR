module hecmw_solver_misc_22

      contains

!C
!C***
!C*** hecmw_matvec_22
!C***
!C
      subroutine hecmw_matvec_22 (hecMESH, hecMAT, X, Y, n, COMMtime)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      integer(kind=kint):: n
      real(kind=kreal), dimension(2*n) :: X, Y
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: START_TIME, END_TIME
      integer(kind=kint) :: i, j, jS, jE, in
      real(kind=kreal) :: YV1, YV2, X1, X2

      if (present(COMMtime)) then
      START_TIME= HECMW_WTIME()
      call hecmw_update_2_R (hecMESH, X, hecMAT%NP)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      else
      call hecmw_update_2_R (hecMESH, X, hecMAT%NP)
      endif

      do i= 1, hecMAT%N
        X1= X(2*i-1)
        X2= X(2*i  )
        YV1= hecMAT%D(4*i-3)*X1 + hecMAT%D(4*i-2)*X2
        YV2= hecMAT%D(4*i-1)*X1 + hecMAT%D(4*i  )*X2

        jS= hecMAT%indexL(i-1) + 1
        jE= hecMAT%indexL(i  )
        do j= jS, jE
          in  = hecMAT%itemL(j)
          X1 = X(2*in-1)
          X2 = X(2*in  )
          YV1= YV1 + hecMAT%AL(4*j-3)*X1 + hecMAT%AL(4*j-2)*X2
          YV2= YV2 + hecMAT%AL(4*j-1)*X1 + hecMAT%AL(4*j  )*X2
        enddo
        jS= hecMAT%indexU(i-1) + 1
        jE= hecMAT%indexU(i  )
        do j= jS, jE
          in  = hecMAT%itemU(j)
          X1 = X(2*in-1)
          X2 = X(2*in  )
          YV1= YV1 + hecMAT%AU(4*j-3)*X1 + hecMAT%AU(4*j-2)*X2
          YV2= YV2 + hecMAT%AU(4*j-1)*X1 + hecMAT%AU(4*j  )*X2
        enddo
        Y(2*i-1)= YV1
        Y(2*i  )= YV2
      enddo

      end subroutine hecmw_matvec_22

!C
!C***
!C*** hecmw_matresid_22
!C***
!C
      subroutine hecmw_matresid_22 (hecMESH, hecMAT, X, B, R, COMMtime)
      use hecmw_util

      implicit none
      real(kind=kreal) :: X(:), B(:), R(:)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i

      if (present(COMMtime)) then
      call hecmw_matvec_22 (hecMESH, hecMAT, X, R, hecMAT%NP, COMMtime)
      else
      call hecmw_matvec_22 (hecMESH, hecMAT, X, R, hecMAT%NP)
      endif
      do i = 1, hecMAT%N * 2
        R(i) = B(i) - R(i)
      enddo

      end subroutine hecmw_matresid_22

!C
!C***
!C*** hecmw_rel_resid_L2_22
!C***
!C
      function hecmw_rel_resid_L2_22 (hecMESH, hecMAT, COMMtime)
      use hecmw_util
      use hecmw_solver_misc

      implicit none
      real(kind=kreal) :: hecmw_rel_resid_L2_22
      type ( hecmwST_local_mesh ), intent(in) :: hecMESH
      type ( hecmwST_matrix     ), intent(in) :: hecMAT
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: r(hecMAT%NDOF*hecMAT%NP)
      real(kind=kreal) :: bnorm2, rnorm2

      if (present(COMMtime)) then
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, COMMtime)
      call hecmw_matresid_22(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, COMMtime)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, COMMtime)
      else
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2)
      call hecmw_matresid_22(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2)
      endif
      hecmw_rel_resid_L2_22 = sqrt(rnorm2 / bnorm2)

      end function hecmw_rel_resid_L2_22

end module hecmw_solver_misc_22
