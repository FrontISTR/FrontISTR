!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : coarsest-level solver
!!
!! Abstract coarsest-solver interface (setup_structure / setup_values / solve);
!! the first backend is a redundant dense symmetric LDL^T factorization (LAPACK
!! dsytrf/dsytrs).  dsytrf (not dpotrf) is used so the same backend handles the
!! symmetric-indefinite coarse matrices of the future u-p extension.  In the MPI
!! phase the coarse matrix will be Allgather'd and factored redundantly on every
!! rank; here (sequential) the level matrix is simply densified in place.
module hecmw_precond_SAAMG_coarse
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_comm, only: hecmw_saamg_abort, hecmw_saamg_check_alloc, &
       hecmw_saamg_require_lapack
  implicit none

  private
  public :: hecmwST_saamg_coarse
  public :: hecmw_saamg_coarse_setup
  public :: hecmw_saamg_coarse_solve
  public :: hecmw_saamg_coarse_free

  !> Dense coarsest solver state: symmetric LDL^T (dsytrf) or, for a non-symmetric
  !! coarse operator (friction contact), general LU (dgetrf).
  type hecmwST_saamg_coarse
    integer(kind=kint) :: n = 0
    logical            :: symmetric = .true.     !< .true.=dsytrf LDL^T, .false.=dgetrf LU
    real(kind=kreal),   allocatable :: ldl(:,:)  !< n x n, dsytrf/dgetrf factors
    integer(kind=kint), allocatable :: ipiv(:)
    integer(kind=kint) :: n_neg = 0              !< # negative eigenvalues (inertia; LDL^T only)
  end type hecmwST_saamg_coarse

contains

  !> Densify Ac and factor it: symmetric -> LDL^T (dsytrf), else -> LU (dgetrf, for a
  !! non-symmetric coarse operator).  bcsr_to_dense yields the full matrix, so LU is
  !! correct for non-symmetric A (dsytrf would use only the lower triangle).
  subroutine hecmw_saamg_coarse_setup(Ac, symmetric, cs)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: Ac
    logical,                  intent(in)  :: symmetric
    type(hecmwST_saamg_coarse), intent(out) :: cs
    real(kind=kreal), allocatable :: work(:)
    real(kind=kreal) :: wq(1)
    integer(kind=kint) :: n, lwork, info, astat
#ifdef HECMW_WITH_LAPACK
    external :: dsytrf, dgetrf

    n = Ac%n
    call hecmw_saamg_coarse_free(cs)
    cs%n = n; cs%symmetric = symmetric
    allocate(cs%ldl(n,n), cs%ipiv(n), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'coarse_setup dense factor (ldl: n*n)')
    call hecmw_saamg_to_dense_blk(Ac, cs%ldl)   ! densify directly from block storage

    if (symmetric) then
      call dsytrf('L', n, cs%ldl, n, cs%ipiv, wq, -1, info); lwork = max(int(wq(1)), 1)
      allocate(work(lwork))
      call dsytrf('L', n, cs%ldl, n, cs%ipiv, work, lwork, info)
      deallocate(work)
      if (info /= 0) then
        write(*,'(a,i0)') 'hecmw_saamg_coarse_setup: dsytrf info=', info
        call hecmw_saamg_abort('coarse_setup: dense LDL^T factorization failed (dsytrf)')
      end if
      cs%n_neg = ldlt_n_negative(cs)
    else
      call dgetrf(n, n, cs%ldl, n, cs%ipiv, info)
      if (info /= 0) then
        write(*,'(a,i0)') 'hecmw_saamg_coarse_setup: dgetrf info=', info
        call hecmw_saamg_abort('coarse_setup: dense LU factorization failed (dgetrf)')
      end if
      cs%n_neg = 0
    end if
#else
    call hecmw_saamg_require_lapack('dense coarsest factorization (dsytrf/dgetrf)')
#endif
  end subroutine hecmw_saamg_coarse_setup

  !> Number of negative eigenvalues (inertia) from the Bunch-Kaufman LDL^T factors
  !! produced by dsytrf('L'): walk the block-diagonal D (1x1 and 2x2 pivots).  For
  !! an SPD coarsest operator this is 0; a positive count signals an indefinite
  !! (non-SPD) operator.
  function ldlt_n_negative(cs) result(nneg)
    type(hecmwST_saamg_coarse), intent(in) :: cs
    integer(kind=kint) :: nneg, k
    real(kind=kreal) :: a, b, c, det, tr
    nneg = 0
    k = 1
    do while (k <= cs%n)
      if (cs%ipiv(k) > 0) then            ! 1x1 pivot
        if (cs%ldl(k,k) < 0.0d0) nneg = nneg + 1
        k = k + 1
      else                                ! 2x2 pivot on (k, k+1)
        a = cs%ldl(k,k); b = cs%ldl(k+1,k); c = cs%ldl(k+1,k+1)
        det = a*c - b*b; tr = a + c
        if (det < 0.0d0) then
          nneg = nneg + 1                  ! one positive, one negative
        else if (tr < 0.0d0) then
          nneg = nneg + 2                  ! both negative
        end if
        k = k + 2
      end if
    end do
  end function ldlt_n_negative

  !> Solve  Ac x = b  using the stored LDL^T factors.
  subroutine hecmw_saamg_coarse_solve(cs, b, x)
    implicit none
    type(hecmwST_saamg_coarse), intent(in)  :: cs
    real(kind=kreal),         intent(in)  :: b(:)
    real(kind=kreal),         intent(out) :: x(:)
    integer(kind=kint) :: info
#ifdef HECMW_WITH_LAPACK
    external :: dsytrs, dgetrs
    x(1:cs%n) = b(1:cs%n)
    if (cs%symmetric) then
      call dsytrs('L', cs%n, 1, cs%ldl, cs%n, cs%ipiv, x, cs%n, info)
    else
      call dgetrs('N', cs%n, 1, cs%ldl, cs%n, cs%ipiv, x, cs%n, info)
    end if
    if (info /= 0) then
      write(*,'(a,i0)') 'hecmw_saamg_coarse_solve: trs info=', info
      call hecmw_saamg_abort('coarse_solve: dense triangular solve failed')
    end if
#else
    x(1:cs%n) = b(1:cs%n)
    call hecmw_saamg_require_lapack('dense coarsest solve (dsytrs/dgetrs)')
#endif
  end subroutine hecmw_saamg_coarse_solve

  subroutine hecmw_saamg_coarse_free(cs)
    implicit none
    type(hecmwST_saamg_coarse), intent(inout) :: cs
    if (allocated(cs%ldl))  deallocate(cs%ldl)
    if (allocated(cs%ipiv)) deallocate(cs%ipiv)
    cs%n = 0
  end subroutine hecmw_saamg_coarse_free

end module hecmw_precond_SAAMG_coarse
