!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 10 (block-CSR kernels)
!
! Verifies the block-CSR storage and kernels against the scalar kernels, on
! blocked operators:
!   * sync_block: block storage densifies to the same matrix as the scalar CSR;
!   * matvec_blk == scalar matvec (square nb=3 and rectangular nb=3,mb=2);
!   * spgemm_blk == scalar SpGEMM  (A*P, rectangular blocks);
!   * transpose_blk == scalar transpose (P, swaps block dims);
!   * galerkin_blk == scalar Galerkin (Ac = P^T A P, 2x2 coarse blocks).
!-------------------------------------------------------------------------------
program test_step10
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  implicit none

  integer(kind=kint) :: npass, nfail
  integer(kind=kint), parameter :: nnode = 4, nb = 3, m = 2
  integer(kind=kint), parameter :: n = nnode*nb, ncc = nnode*m
  type(hecmwST_saamg_bcsr) :: A, P
  npass = 0; nfail = 0

  call build_A(A)
  call build_P(P)

  call test_sync(A, P, npass, nfail)
  call test_matvec(A, P, npass, nfail)
  call test_spgemm(A, P, npass, nfail)
  call test_transpose(P, npass, nfail)
  call test_galerkin(A, P, npass, nfail)
  call test_triple(A, P, npass, nfail)

  call hecmw_saamg_bcsr_free(A); call hecmw_saamg_bcsr_free(P)
  write(*,'(/,a,i0,a,i0,a)') '=== Step 10: ', npass, ' passed, ', nfail, ' failed ==='
  if (nfail > 0) stop 1

contains

  subroutine check(name, cond, npass, nfail)
    character(len=*),   intent(in)    :: name
    logical,            intent(in)    :: cond
    integer(kind=kint), intent(inout) :: npass, nfail
    if (cond) then
      npass = npass + 1; write(*,'(a,a)') '  [PASS] ', name
    else
      nfail = nfail + 1; write(*,'(a,a)') '  [FAIL] ', name
    end if
  end subroutine check

  function lcg(seed) result(r)
    integer(kind=8), intent(inout) :: seed
    real(kind=kreal) :: r
    seed = seed*6364136223846793005_8 + 1442695040888963407_8
    r = real(seed,kreal)/9.223372036854776d18    ! ~ uniform in (-1,1)
  end function lcg

  !> SPD block-tridiagonal operator (nnode nodes, nb dofs each).
  subroutine build_A(A)
    type(hecmwST_saamg_bcsr), intent(out) :: A
    real(kind=kreal) :: d(n,n)
    integer(kind=kint) :: i, j, ib
    d = 0.0d0
    do ib = 1, nnode
      ! diagonal block: SPD (I*10 + small symmetric perturbation)
      do i = 1, nb; do j = 1, nb
        if (i == j) then; d((ib-1)*nb+i,(ib-1)*nb+j) = 10.0d0
        else; d((ib-1)*nb+i,(ib-1)*nb+j) = 0.5d0; end if
      end do; end do
      ! off-diagonal coupling to next node: -I
      if (ib < nnode) then
        do i = 1, nb
          d((ib-1)*nb+i, ib*nb+i) = -1.0d0
          d(ib*nb+i, (ib-1)*nb+i) = -1.0d0
        end do
      end if
    end do
    call hecmw_saamg_bcsr_from_dense(d, n, nb, A)        ! mb defaults to nb
  end subroutine build_A

  !> Random tentative-like prolongator P (nnode fine nodes -> nnode coarse nodes
  !! is overkill; use 2 coarse nodes), nb x m blocks.  One nonzero block per fine
  !! node, mapping fine node ib to coarse node (ib-1)/2+1.
  subroutine build_P(P)
    type(hecmwST_saamg_bcsr), intent(out) :: P
    integer(kind=kint) :: ti(n*m), tj(n*m), nt, ib, cn, i, c
    real(kind=kreal) :: tv(n*m)
    integer(kind=8) :: seed
    seed = 12345_8
    nt = 0
    do ib = 1, nnode
      cn = (ib-1)/2 + 1                ! coarse node (2 fine -> 1 coarse)
      do i = 1, nb
        do c = 1, m
          nt = nt + 1
          ti(nt) = (ib-1)*nb + i
          tj(nt) = (cn-1)*m + c
          tv(nt) = lcg(seed)
        end do
      end do
    end do
    call hecmw_saamg_bcsr_from_triplets(n, (nnode/2)*m, nb, ti, tj, tv, nt, P, mb=m)
  end subroutine build_P

  !-----------------------------------------------------------------------------
  subroutine test_sync(A, P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: A, P
    integer(kind=kint),    intent(inout) :: npass, nfail
    real(kind=kreal), allocatable :: ds(:,:), db(:,:)
    allocate(ds(A%n,A%ncol), db(A%n,A%ncol))
    call hecmw_saamg_bcsr_to_dense(A, ds)
    call hecmw_saamg_to_dense_blk(A, db)
    call check('sync_block A: block==scalar', maxval(abs(ds-db)) < 1.0d-15, npass, nfail)
    call check('A block dims', A%nb==nb .and. A%mb==nb .and. A%nbrow==nnode, npass, nfail)
    deallocate(ds, db)
    allocate(ds(P%n,P%ncol), db(P%n,P%ncol))
    call hecmw_saamg_bcsr_to_dense(P, ds)
    call hecmw_saamg_to_dense_blk(P, db)
    call check('sync_block P: block==scalar', maxval(abs(ds-db)) < 1.0d-15, npass, nfail)
    call check('P block dims', P%nb==nb .and. P%mb==m .and. P%nbcol==nnode/2, npass, nfail)
    deallocate(ds, db)
  end subroutine test_sync

  subroutine test_matvec(A, P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: A, P
    integer(kind=kint),    intent(inout) :: npass, nfail
    real(kind=kreal) :: x(n), y1(n), y2(n), xc(ncc/2*2), yc1(n), yc2(n)
    integer(kind=8) :: seed
    integer(kind=kint) :: i
    seed = 7_8
    do i = 1, n; x(i) = lcg(seed); end do
    call hecmw_saamg_matvec(A, x, y1)
    call hecmw_saamg_matvec_blk(A, x, y2)
    call check('matvec_blk == matvec (square)', maxval(abs(y1-y2)) < 1.0d-13, npass, nfail)
    do i = 1, (nnode/2)*m; xc(i) = lcg(seed); end do
    call hecmw_saamg_matvec(P, xc, yc1)
    call hecmw_saamg_matvec_blk(P, xc, yc2)
    call check('matvec_blk == matvec (rectangular)', maxval(abs(yc1-yc2)) < 1.0d-13, npass, nfail)
  end subroutine test_matvec

  subroutine test_spgemm(A, P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: A, P
    integer(kind=kint),    intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: C1, C2
    real(kind=kreal), allocatable :: d1(:,:), d2(:,:)
    call hecmw_saamg_spgemm(A, P, C1)
    call hecmw_saamg_spgemm_blk(A, P, C2)
    call check('spgemm_blk dims', C2%n==n .and. C2%ncol==P%ncol .and. &
         C2%nb==nb .and. C2%mb==m, npass, nfail)
    allocate(d1(n,P%ncol), d2(n,P%ncol))
    call hecmw_saamg_bcsr_to_dense(C1, d1)
    call hecmw_saamg_to_dense_blk(C2, d2)
    call check('spgemm_blk == scalar SpGEMM', maxval(abs(d1-d2)) < 1.0d-12, npass, nfail)
    deallocate(d1, d2)
    call hecmw_saamg_bcsr_free(C1); call hecmw_saamg_bcsr_free(C2)
  end subroutine test_spgemm

  subroutine test_transpose(P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: P
    integer(kind=kint),    intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: T1, T2
    real(kind=kreal), allocatable :: d1(:,:), d2(:,:)
    call hecmw_saamg_bcsr_transpose(P, T1)
    call hecmw_saamg_transpose_blk(P, T2)
    call check('transpose_blk dims', T2%n==P%ncol .and. T2%ncol==P%n .and. &
         T2%nb==m .and. T2%mb==nb, npass, nfail)
    allocate(d1(P%ncol,P%n), d2(P%ncol,P%n))
    call hecmw_saamg_bcsr_to_dense(T1, d1)
    call hecmw_saamg_to_dense_blk(T2, d2)
    call check('transpose_blk == scalar transpose', maxval(abs(d1-d2)) < 1.0d-15, npass, nfail)
    deallocate(d1, d2)
    call hecmw_saamg_bcsr_free(T1); call hecmw_saamg_bcsr_free(T2)
  end subroutine test_transpose

  subroutine test_galerkin(A, P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: A, P
    integer(kind=kint),    intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: Ac1, Ac2
    real(kind=kreal), allocatable :: d1(:,:), d2(:,:)
    integer(kind=kint) :: ncoarse
    ncoarse = P%ncol
    call hecmw_saamg_galerkin_local(A, P, m, Ac1)
    call hecmw_saamg_galerkin_blk(A, P, m, Ac2)
    call check('galerkin_blk dims', Ac2%n==ncoarse .and. Ac2%nb==m .and. Ac2%mb==m, &
         npass, nfail)
    allocate(d1(ncoarse,ncoarse), d2(ncoarse,ncoarse))
    call hecmw_saamg_bcsr_to_dense(Ac1, d1)
    call hecmw_saamg_to_dense_blk(Ac2, d2)
    call check('galerkin_blk == scalar Galerkin', maxval(abs(d1-d2)) < 1.0d-11, npass, nfail)
    deallocate(d1, d2)
    call hecmw_saamg_bcsr_free(Ac1); call hecmw_saamg_bcsr_free(Ac2)
  end subroutine test_galerkin

  !> Fused triple product Ac = P^T A P (via triple_blk on Pt=P^T, A, P) must equal
  !! the two-SpGEMM Galerkin, and transpose_blk_rows must equal a full transpose.
  subroutine test_triple(A, P, npass, nfail)
    type(hecmwST_saamg_bcsr), intent(in)    :: A, P
    integer(kind=kint),    intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: Pt, Pt2, Ac1, Ac2
    real(kind=kreal), allocatable :: d1(:,:), d2(:,:)
    integer(kind=kint) :: ncoarse
    ncoarse = P%ncol
    ! transpose_blk_rows on all rows == full transpose
    call hecmw_saamg_transpose_blk(P, Pt)
    call hecmw_saamg_transpose_blk_rows(P, P%nbrow, Pt2)
    allocate(d1(P%ncol,P%n), d2(P%ncol,P%n))
    call hecmw_saamg_to_dense_blk(Pt, d1); call hecmw_saamg_to_dense_blk(Pt2, d2)
    call check('transpose_blk_rows(all) == transpose', maxval(abs(d1-d2)) < 1.0d-15, npass, nfail)
    deallocate(d1, d2)
    ! fused triple product Pt*A*P == galerkin
    call hecmw_saamg_galerkin_local(A, P, m, Ac1)
    call hecmw_saamg_triple_blk(Pt, A, P, Ac2)
    call check('triple_blk dims', Ac2%n==ncoarse .and. Ac2%nb==m .and. Ac2%mb==m, npass, nfail)
    allocate(d1(ncoarse,ncoarse), d2(ncoarse,ncoarse))
    call hecmw_saamg_bcsr_to_dense(Ac1, d1)
    call hecmw_saamg_to_dense_blk(Ac2, d2)
    call check('triple_blk == galerkin (fused P^T A P)', maxval(abs(d1-d2)) < 1.0d-11, npass, nfail)
    deallocate(d1, d2)
    call hecmw_saamg_bcsr_free(Pt); call hecmw_saamg_bcsr_free(Pt2)
    call hecmw_saamg_bcsr_free(Ac1); call hecmw_saamg_bcsr_free(Ac2)
  end subroutine test_triple

end program test_step10
