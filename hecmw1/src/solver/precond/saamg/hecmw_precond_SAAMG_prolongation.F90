!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : prolongation
!!
!! Two operations:
!!   * tentative prolongator P-hat and coarse near-kernel B via per-aggregate
!!     thin Householder QR (LAPACK dgeqrf/dorgqr): B_k = Q_k R_k, Q_k -> P-hat
!!     block column k, R_k -> B_coarse block row k.  By construction
!!     P-hat * B_coarse = B_fine and P-hat^T P-hat = I.
!!   * prolongator smoothing  P = (I - omega D^{-1} A) P-hat, computed column by
!!     column reusing the tested matvec and block-diagonal solve (a SpGEMM-based
!!     path is a later optimization).
!!
!! Coarse layout: aggregate k becomes coarse node k with m dofs; coarse dof
!! (k-1)*m + c.  Excluded nodes (aggr=0) get an all-zero P-hat row.
module hecmw_precond_SAAMG_prolongation
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_smoother
  use hecmw_precond_SAAMG_comm, only: hecmw_saamg_abort, hecmw_saamg_require_lapack
  implicit none

  private
  public :: hecmw_saamg_tentative
  public :: hecmw_saamg_smooth_prolongator

contains

  !> Build P-hat (n x naggr*m block-CSR) and B_coarse (naggr*m x m) from B_fine.
  subroutine hecmw_saamg_tentative(bfine, n, m, nb, aggr, naggr, phat, bcoarse)
    implicit none
    integer(kind=kint),    intent(in)  :: n, m, nb, naggr
    real(kind=kreal),      intent(in)  :: bfine(n, m)
    integer(kind=kint),    intent(in)  :: aggr(:)
    type(hecmwST_saamg_bcsr), intent(out) :: phat
    real(kind=kreal), allocatable, intent(out) :: bcoarse(:,:)
    integer(kind=kint), allocatable :: memcnt(:), aggptr(:), aggnodes(:), pos(:)
    integer(kind=kint), allocatable :: bi(:), bj(:)
    real(kind=kreal),   allocatable :: bv(:,:), bk(:,:), tau(:), work(:)
    integer(kind=kint) :: nnode, ncoarse, inode, k, l, ii, c, r, node
    integer(kind=kint) :: nk, rows, nt, maxt, lwork, info, rwork
    real(kind=kreal) :: wq(1)
#ifdef HECMW_WITH_LAPACK
    external :: dgeqrf, dorgqr

    nnode = n / nb
    ncoarse = naggr * m
    allocate(bcoarse(ncoarse, m)); bcoarse = 0.0d0

    ! member node list per aggregate (CSR-like)
    allocate(memcnt(naggr)); memcnt = 0
    do inode = 1, nnode
      if (aggr(inode) > 0) memcnt(aggr(inode)) = memcnt(aggr(inode)) + 1
    end do
    allocate(aggptr(naggr+1)); aggptr(1) = 0
    do k = 1, naggr
      aggptr(k+1) = aggptr(k) + memcnt(k)
    end do
    allocate(aggnodes(aggptr(naggr+1)), pos(naggr))
    do k = 1, naggr
      pos(k) = aggptr(k)
    end do
    do inode = 1, nnode
      k = aggr(inode)
      if (k > 0) then
        pos(k) = pos(k) + 1
        aggnodes(pos(k)) = inode
      end if
    end do

    ! block-triplet buffer (one nb x m block per assigned node)
    maxt = aggptr(naggr+1)
    allocate(bi(maxt), bj(maxt), bv(nb*m, maxt)); nt = 0

    do k = 1, naggr
      nk = memcnt(k); rows = nk * nb
      if (rows < m) then
        write(*,'(a,i0,a,i0,a,i0)') 'hecmw_saamg_tentative: aggregate ', k, &
             ' has rows=', rows, ' < m=', m
        call hecmw_saamg_abort('tentative: aggregate has fewer rows than near-kernel size m')
      end if
      allocate(bk(rows, m), tau(m))
      do l = 0, nk-1
        node = aggnodes(aggptr(k)+l+1)
        do ii = 1, nb
          do c = 1, m
            bk(l*nb+ii, c) = bfine((node-1)*nb+ii, c)
          end do
        end do
      end do

      ! workspace query (max of dgeqrf / dorgqr)
      call dgeqrf(rows, m, bk, rows, tau, wq, -1, info); lwork = int(wq(1))
      call dorgqr(rows, m, m, bk, rows, tau, wq, -1, info); rwork = int(wq(1))
      lwork = max(lwork, rwork, 1)
      allocate(work(lwork))

      call dgeqrf(rows, m, bk, rows, tau, work, lwork, info)
      if (info /= 0) then; write(*,*) 'dgeqrf info=', info
        call hecmw_saamg_abort('tentative: per-aggregate QR failed (dgeqrf)'); end if
      ! Extract R into B_coarse block row k.  R may be (near-)singular when the
      ! aggregate is geometrically degenerate (e.g. collinear nodes -> a rotation
      ! axis is unrepresented): that aggregate then carries fewer than m
      ! independent near-kernel modes.  This is admissible -- dorgqr still returns
      ! an orthonormal Q (all m columns), so P-hat^T P-hat = I and the exact
      ! identity P-hat * B_coarse = B_fine both hold; only B_coarse is rank
      ! deficient, which the next-coarser level tolerates.
      do r = 1, m
        do c = r, m
          bcoarse((k-1)*m+r, c) = bk(r, c)
        end do
      end do

      call dorgqr(rows, m, m, bk, rows, tau, work, lwork, info)
      if (info /= 0) then; write(*,*) 'dorgqr info=', info
        call hecmw_saamg_abort('tentative: per-aggregate Q formation failed (dorgqr)'); end if

      ! emit one P-hat block per member node (rows = node dofs, cols = coarse node k)
      do l = 0, nk-1
        node = aggnodes(aggptr(k)+l+1)
        nt = nt + 1; bi(nt) = node; bj(nt) = k
        do c = 1, m
          do ii = 1, nb
            bv((c-1)*nb+ii, nt) = bk(l*nb+ii, c)
          end do
        end do
      end do

      deallocate(bk, tau, work)
    end do

    call hecmw_saamg_bcsr_from_block_triplets(nnode, naggr, nb, m, bi, bj, bv, nt, phat)

    deallocate(memcnt, aggptr, aggnodes, pos, bi, bj, bv)
#else
    call hecmw_saamg_require_lapack('per-aggregate tentative QR (dgeqrf/dorgqr)')
#endif
  end subroutine hecmw_saamg_tentative

  !> Smoothed prolongator  P = P-hat - omega * D^{-1} (A P-hat), via SpGEMM.
  !! Computes A*P-hat once as a sparse product (instead of one full A-matvec per
  !! coarse column, which is O(ncoarse * nnz_A)), applies the block-diagonal
  !! D^{-1} row-wise, then forms the sparse difference by summing triplets.
  subroutine hecmw_saamg_smooth_prolongator(A, D, omega, phat, p)
    implicit none
    type(hecmwST_saamg_bcsr),  intent(in)  :: A
    type(hecmwST_saamg_blockdiag), intent(in) :: D
    real(kind=kreal),       intent(in)  :: omega
    type(hecmwST_saamg_bcsr),  intent(in)  :: phat
    type(hecmwST_saamg_bcsr),  intent(out) :: p
    type(hecmwST_saamg_bcsr) :: ap
    integer(kind=kint), allocatable :: bi(:), bj(:)
    real(kind=kreal),   allocatable :: bv(:,:)
    integer(kind=kint) :: nb, m, bs, nrow, i, t, nt

    nb = phat%nb; m = phat%mb; bs = nb*m
    ! P has only the INTERNAL fine rows (A%nbrow): for the distributed operator A
    ! is rectangular (nint x NP) and phat is the halo-extended P-hat (NP rows);
    ! only the first A%nbrow block rows of phat are kept.  Sequential: A%nbrow ==
    ! phat%nbrow, so this is a no-op restriction.
    nrow = A%nbrow
    call hecmw_saamg_spgemm_blk(A, phat, ap)              ! ap = A P-hat   (block, nint rows)
    call hecmw_saamg_blockdiag_apply_bcsr_inplace(D, ap)  ! ap := D^{-1} A P-hat (in place)

    ! P = P-hat - omega * ap : sum the two block sets (from_block_triplets sums dups)
    allocate(bi(phat%nnzb + ap%nnzb), bj(phat%nnzb + ap%nnzb), bv(bs, phat%nnzb + ap%nnzb))
    nt = 0
    do i = 1, nrow                                     ! + P-hat (internal rows only)
      do t = phat%browptr(i), phat%browptr(i+1)-1
        nt = nt + 1; bi(nt) = i; bj(nt) = phat%bcol(t)
        bv(1:bs, nt) = phat%bval((t-1)*bs+1:t*bs)
      end do
    end do
    do i = 1, ap%nbrow                                 ! - omega D^{-1} A P-hat
      do t = ap%browptr(i), ap%browptr(i+1)-1
        nt = nt + 1; bi(nt) = i; bj(nt) = ap%bcol(t)
        bv(1:bs, nt) = -omega * ap%bval((t-1)*bs+1:t*bs)
      end do
    end do

    call hecmw_saamg_bcsr_from_block_triplets(nrow, phat%nbcol, nb, m, bi, bj, bv, nt, p)

    deallocate(bi, bj, bv)
    call hecmw_saamg_bcsr_free(ap)
  end subroutine hecmw_saamg_smooth_prolongator

end module hecmw_precond_SAAMG_prolongation
