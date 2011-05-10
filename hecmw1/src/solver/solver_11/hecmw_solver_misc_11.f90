    module hecmw_solver_misc_11

      contains

!C
!C***
!C*** hecmw_matvec_11
!C***
!C
      subroutine hecmw_matvec_11 (hecMESH, hecMAT, X, Y, n)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      integer(kind=kint):: n
      real(kind=kreal), dimension(n) :: X, Y
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: i, j, jS, jE, in
      real(kind=kreal) :: YV

      call hecmw_update_1_R (hecMESH, X, hecMAT%NP)

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
      subroutine hecmw_matresid_11 (hecMESH, hecMAT, X, B, R)
      use hecmw_util

      implicit none
      real(kind=kreal) :: X(:), B(:), R(:)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: i

      call hecmw_matvec_11 (hecMESH, hecMAT, X, R, hecMAT%NP)
      do i = 1, hecMAT%N
        R(i) = B(i) - R(i)
      enddo

      end subroutine hecmw_matresid_11

end module hecmw_solver_misc_11

