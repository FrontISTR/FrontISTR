!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Dense tile kernels of the multifrontal solver, symmetric (LDLt) mode.
!>
!> A tile is a column major block a(h,w) with leading dimension h. Diagonal tiles hold the unit
!> lower triangle of L below the diagonal and D on the diagonal; their strict upper part is not
!> referenced. The BLAS3 path is used when HECMW_WITH_LAPACK is defined, otherwise plain loops.
module hecmw_mf_kernel
  use hecmw_util
  implicit none

  private
  public :: hecmw_mf_kernel_ldlt
  public :: hecmw_mf_kernel_trsm
  public :: hecmw_mf_kernel_scale
  public :: hecmw_mf_kernel_gemm
  public :: hecmw_mf_kernel_trsv
  public :: hecmw_mf_kernel_trsv_t
  public :: hecmw_mf_kernel_gemv
  public :: hecmw_mf_kernel_gemv_t

contains

  !> In place LDLt of the lower triangle of a diagonal tile without pivoting; info is the index
  !> of the first pivot that is not positive, 0 on success.
  subroutine hecmw_mf_kernel_ldlt(n, a, info)
    implicit none
    integer(kind=kint), intent(in) :: n
    real(kind=kreal), intent(inout) :: a(n,n)
    integer(kind=kint), intent(out) :: info
    integer(kind=kint) :: i, j, k
    real(kind=kreal) :: d, l

    info = 0
    do k = 1, n
      d = a(k,k)
      if (d <= 0.0d0) then
        info = k
        return
      endif
      do j = k+1, n
        l = a(j,k) / d
        do i = j, n
          a(i,j) = a(i,j) - a(i,k)*l
        enddo
      enddo
      do i = k+1, n
        a(i,k) = a(i,k) / d
      enddo
    enddo
  end subroutine hecmw_mf_kernel_ldlt

  !> a(m,n) <- a * L^-T * D^-1 with L, D taken from the factored diagonal tile l(n,n).
  subroutine hecmw_mf_kernel_trsm(m, n, l, a)
    implicit none
    integer(kind=kint), intent(in) :: m, n
    real(kind=kreal), intent(in) :: l(n,n)
    real(kind=kreal), intent(inout) :: a(m,n)
    integer(kind=kint) :: i, j
#ifdef HECMW_WITH_LAPACK
    external :: dtrsm

    call dtrsm('R', 'L', 'T', 'U', m, n, 1.0d0, l, n, a, m)
    do j = 1, n
      do i = 1, m
        a(i,j) = a(i,j) / l(j,j)
      enddo
    enddo
#else
    integer(kind=kint) :: k

    do j = 1, n
      do k = 1, j-1
        do i = 1, m
          a(i,j) = a(i,j) - a(i,k)*l(j,k)
        enddo
      enddo
      do i = 1, m
        a(i,j) = a(i,j) / l(j,j)
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_trsm

  !> w(m,n) <- a * D with D the diagonal of the factored tile l(n,n).
  subroutine hecmw_mf_kernel_scale(m, n, l, a, w)
    implicit none
    integer(kind=kint), intent(in) :: m, n
    real(kind=kreal), intent(in) :: l(n,n)
    real(kind=kreal), intent(in) :: a(m,n)
    real(kind=kreal), intent(out) :: w(m,n)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = 1, m
        w(i,j) = a(i,j) * l(j,j)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_scale

  !> c(m,n) <- c - w(m,k) * b(n,k)^T
  subroutine hecmw_mf_kernel_gemm(m, n, k, w, b, c)
    implicit none
    integer(kind=kint), intent(in) :: m, n, k
    real(kind=kreal), intent(in) :: w(m,k)
    real(kind=kreal), intent(in) :: b(n,k)
    real(kind=kreal), intent(inout) :: c(m,n)
#ifdef HECMW_WITH_LAPACK
    external :: dgemm

    call dgemm('N', 'T', m, n, k, -1.0d0, w, m, b, n, 1.0d0, c, m)
#else
    integer(kind=kint) :: i, j, p

    do j = 1, n
      do p = 1, k
        do i = 1, m
          c(i,j) = c(i,j) - w(i,p)*b(j,p)
        enddo
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_gemm

  !> x <- L^-1 x with the unit lower triangle of the diagonal tile l(n,n).
  subroutine hecmw_mf_kernel_trsv(n, l, x)
    implicit none
    integer(kind=kint), intent(in) :: n
    real(kind=kreal), intent(in) :: l(n,n)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = j+1, n
        x(i) = x(i) - l(i,j)*x(j)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_trsv

  !> x <- L^-T x with the unit lower triangle of the diagonal tile l(n,n).
  subroutine hecmw_mf_kernel_trsv_t(n, l, x)
    implicit none
    integer(kind=kint), intent(in) :: n
    real(kind=kreal), intent(in) :: l(n,n)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: i, j

    do j = n, 1, -1
      do i = j+1, n
        x(j) = x(j) - l(i,j)*x(i)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_trsv_t

  !> y(m) <- y - a(m,n) * x(n)
  subroutine hecmw_mf_kernel_gemv(m, n, a, x, y)
    implicit none
    integer(kind=kint), intent(in) :: m, n
    real(kind=kreal), intent(in) :: a(m,n)
    real(kind=kreal), intent(in) :: x(n)
    real(kind=kreal), intent(inout) :: y(m)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = 1, m
        y(i) = y(i) - a(i,j)*x(j)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_gemv

  !> y(n) <- y - a(m,n)^T * x(m)
  subroutine hecmw_mf_kernel_gemv_t(m, n, a, x, y)
    implicit none
    integer(kind=kint), intent(in) :: m, n
    real(kind=kreal), intent(in) :: a(m,n)
    real(kind=kreal), intent(in) :: x(m)
    real(kind=kreal), intent(inout) :: y(n)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = 1, m
        y(j) = y(j) - a(i,j)*x(i)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_gemv_t

end module hecmw_mf_kernel
