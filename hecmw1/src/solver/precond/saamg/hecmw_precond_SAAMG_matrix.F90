!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : internal block-CSR matrix
!!
!! Mesh-independent block-CSR matrix used at every level of the SA-AMG hierarchy.
!! Each stored nonzero is an nb x mb block (column major), with nb = row-block size
!! (dofs per row node) and mb = column-block size (dofs per column node).  For a
!! square operator nb == mb; for a prolongator P the row block is the fine node
!! size and the column block is the coarse near-kernel size m (nb /= mb).  matvec,
!! SpGEMM, transpose and the Galerkin product all run as block algorithms, so the
!! node-block structure is exploited directly (dense per-block GEMM/GEMV).  Indices
!! are 1-based to match Fortran/FrontISTR conventions.
module hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_comm, only: hecmwST_saamg_comm, hecmw_saamg_comm_update, &
       hecmw_saamg_abort, hecmw_saamg_check_alloc
  implicit none

  private
  public :: hecmwST_saamg_bcsr
  public :: hecmw_saamg_bcsr_free
  public :: hecmw_saamg_bcsr_copy
  public :: hecmw_saamg_bcsr_move
  public :: hecmw_saamg_bcsr_from_dense
  public :: hecmw_saamg_bcsr_from_triplets
  public :: hecmw_saamg_bcsr_from_block_triplets
  public :: hecmw_saamg_bcsr_transpose
  public :: hecmw_saamg_bcsr_to_dense
  public :: hecmw_saamg_matvec
  public :: hecmw_saamg_matvec_d
  public :: hecmw_saamg_spgemm
  public :: hecmw_saamg_galerkin_local
  public :: hecmw_saamg_is_symmetric
  ! block-CSR kernels (also reachable directly; the names above delegate to these)
  public :: hecmw_saamg_to_dense_blk
  public :: hecmw_saamg_matvec_blk
  public :: hecmw_saamg_transpose_blk
  public :: hecmw_saamg_spgemm_blk
  public :: hecmw_saamg_galerkin_blk
  public :: hecmw_saamg_triple_blk
  public :: hecmw_saamg_transpose_blk_rows

  !> Block-CSR sparse matrix (n rows x ncol columns; square when ncol == n).
  type hecmwST_saamg_bcsr
    integer(kind=kint) :: n     = 0  !< number of rows ( = nbrow * nb )
    integer(kind=kint) :: ncol  = 0  !< number of columns ( = nbcol * mb )
    integer(kind=kint) :: nb    = 0  !< row-block size (dofs per row node)
    integer(kind=kint) :: mb    = 0  !< column-block size (dofs per col node; = nb if square)
    integer(kind=kint) :: nnz   = 0  !< stored entries ( = nnzb * nb * mb )
    integer(kind=kint) :: nbrow = 0  !< number of block rows ( = n / nb )
    integer(kind=kint) :: nbcol = 0  !< number of block columns ( = ncol / mb )
    integer(kind=kint) :: nnzb  = 0  !< number of stored nonzero blocks
    integer(kind=kint), allocatable :: browptr(:)  !< block-row pointer, size nbrow+1
    integer(kind=kint), allocatable :: bcol(:)     !< block-column index, size nnzb
    real(kind=kreal),   allocatable :: bval(:)     !< nb*mb*nnzb, each block column major
  end type hecmwST_saamg_bcsr

contains

  !> Release the storage held by a hecmwST_saamg_bcsr
  subroutine hecmw_saamg_bcsr_free(A)
    implicit none
    type(hecmwST_saamg_bcsr), intent(inout) :: A
    if (allocated(A%browptr)) deallocate(A%browptr)
    if (allocated(A%bcol))    deallocate(A%bcol)
    if (allocated(A%bval))    deallocate(A%bval)
    A%n = 0; A%ncol = 0; A%nb = 0; A%mb = 0; A%nnz = 0
    A%nbrow = 0; A%nbcol = 0; A%nnzb = 0
  end subroutine hecmw_saamg_bcsr_free

  !> Copy:  B = A.
  subroutine hecmw_saamg_bcsr_copy(A, B)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_bcsr), intent(out) :: B
    integer(kind=kint) :: astat
    call hecmw_saamg_bcsr_free(B)
    B%n = A%n; B%ncol = A%ncol; B%nb = A%nb; B%mb = A%mb; B%nnz = A%nnz
    B%nbrow = A%nbrow; B%nbcol = A%nbcol; B%nnzb = A%nnzb
    if (allocated(A%browptr)) then
      allocate(B%browptr(A%nbrow+1), B%bcol(A%nnzb), B%bval(A%nb*A%mb*A%nnzb), stat=astat)
      call hecmw_saamg_check_alloc(astat, 'bcsr_copy (browptr/bcol/bval)')
      B%browptr = A%browptr; B%bcol = A%bcol; B%bval = A%bval
    end if
  end subroutine hecmw_saamg_bcsr_copy

  !> Move:  B = A, transferring A's storage (move_alloc, no copy); A is emptied.
  !! Use when the source is no longer needed -- avoids holding the operator twice.
  subroutine hecmw_saamg_bcsr_move(A, B)
    implicit none
    type(hecmwST_saamg_bcsr), intent(inout) :: A
    type(hecmwST_saamg_bcsr), intent(out)   :: B
    call hecmw_saamg_bcsr_free(B)
    B%n = A%n; B%ncol = A%ncol; B%nb = A%nb; B%mb = A%mb; B%nnz = A%nnz
    B%nbrow = A%nbrow; B%nbcol = A%nbcol; B%nnzb = A%nnzb
    if (allocated(A%browptr)) call move_alloc(A%browptr, B%browptr)
    if (allocated(A%bcol))    call move_alloc(A%bcol,    B%bcol)
    if (allocated(A%bval))    call move_alloc(A%bval,    B%bval)
    call hecmw_saamg_bcsr_free(A)
  end subroutine hecmw_saamg_bcsr_move

  !> Build a block-CSR matrix from a dense n x n array.  A block is stored when any
  !! of its entries is nonzero (the full nb x mb block is kept, in-block zeros and
  !! all).  Intended for unit tests and small reference matrices, not production.
  subroutine hecmw_saamg_bcsr_from_dense(dense, n, nb, A, mb)
    implicit none
    integer(kind=kint), intent(in)  :: n, nb
    real(kind=kreal),   intent(in)  :: dense(n,n)
    type(hecmwST_saamg_bcsr),     intent(out) :: A
    integer(kind=kint), optional, intent(in) :: mb   !< column-block size (default nb)
    real(kind=kreal), allocatable :: blk(:)
    integer(kind=kint) :: mb_, nbrow, nbcol, bs, i, j, ii, jj, kout, rbase, cbase
    logical :: any_nz

    mb_ = nb; if (present(mb)) mb_ = mb
    call hecmw_saamg_bcsr_free(A)
    A%n = n; A%ncol = n; A%nb = nb; A%mb = mb_
    nbrow = n / nb; nbcol = n / mb_; bs = nb*mb_
    A%nbrow = nbrow; A%nbcol = nbcol
    allocate(A%browptr(nbrow+1), blk(bs))

    ! pass 1: count nonzero blocks per block row
    A%browptr(1) = 1; A%nnzb = 0
    do i = 1, nbrow
      do j = 1, nbcol
        if (block_nonzero(i, j)) A%nnzb = A%nnzb + 1
      end do
      A%browptr(i+1) = A%nnzb + 1
    end do
    allocate(A%bcol(A%nnzb), A%bval(bs*A%nnzb)); A%nnz = bs*A%nnzb

    ! pass 2: copy each nonzero block
    kout = 0
    do i = 1, nbrow
      rbase = (i-1)*nb
      do j = 1, nbcol
        cbase = (j-1)*mb_; any_nz = .false.
        do jj = 1, mb_
          do ii = 1, nb
            blk((jj-1)*nb+ii) = dense(rbase+ii, cbase+jj)
            if (blk((jj-1)*nb+ii) /= 0.0d0) any_nz = .true.
          end do
        end do
        if (any_nz) then
          kout = kout + 1; A%bcol(kout) = j
          A%bval((kout-1)*bs+1:kout*bs) = blk(1:bs)
        end if
      end do
    end do
    deallocate(blk)
  contains
    logical function block_nonzero(ib, jb) result(nz)
      integer(kind=kint), intent(in) :: ib, jb
      integer(kind=kint) :: a, b
      nz = .false.
      do b = 1, mb_
        do a = 1, nb
          if (dense((ib-1)*nb+a, (jb-1)*mb_+b) /= 0.0d0) then; nz = .true.; return; end if
        end do
      end do
    end function block_nonzero
  end subroutine hecmw_saamg_bcsr_from_dense

  !> Sparse matrix-vector product : y = A x.
  subroutine hecmw_saamg_matvec(A, x, y)
    implicit none
    type(hecmwST_saamg_bcsr),   intent(in)  :: A
    real(kind=kreal), intent(in)  :: x(:)
    real(kind=kreal), intent(out) :: y(:)
    call hecmw_saamg_matvec_blk(A, x, y)
  end subroutine hecmw_saamg_matvec

  !> Distributed matvec  y = A x : refresh x's halo region from the owning ranks
  !! (via the level comm table), then run the local product.  x must have length
  !! nb*nnode (internal + halo room); y has length A%n = nb*nint (internal rows).
  !! Reduces to hecmw_saamg_matvec when the comm table has no neighbors (serial).
  subroutine hecmw_saamg_matvec_d(cmt, A, x, y)
    implicit none
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    type(hecmwST_saamg_bcsr),  intent(in)    :: A
    real(kind=kreal),       intent(inout) :: x(:)
    real(kind=kreal),       intent(out)   :: y(:)
    call hecmw_saamg_comm_update(cmt, A%nb, x)
    call hecmw_saamg_matvec(A, x, y)
  end subroutine hecmw_saamg_matvec_d

  !> Relative asymmetry  ||A - A^T||_F / ||A||_F  (debug / verification helper).
  !! Returns 0 for a structurally and numerically symmetric matrix.  Square only.
  function hecmw_saamg_is_symmetric(A) result(rel)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in) :: A
    real(kind=kreal) :: rel
    real(kind=kreal) :: nrm, dif, aij, aji
    integer(kind=kint) :: i, j, t, tt, ii, jj, m, boff, foff
    logical :: found

    m = A%nb                                   ! square operator: nb == mb
    nrm = 0.0d0; dif = 0.0d0
    do i = 1, A%nbrow
      do t = A%browptr(i), A%browptr(i+1)-1
        j = A%bcol(t); boff = (t-1)*m*m
        ! locate the transposed block A(j,i)
        found = .false.; foff = 0
        do tt = A%browptr(j), A%browptr(j+1)-1
          if (A%bcol(tt) == i) then; found = .true.; foff = (tt-1)*m*m; exit; end if
        end do
        do jj = 1, m
          do ii = 1, m
            aij = A%bval(boff + (jj-1)*m + ii)
            nrm = nrm + aij*aij
            aji = 0.0d0
            if (found) aji = A%bval(foff + (ii-1)*m + jj)   ! A(j,i) entry (jj,ii)
            dif = dif + (aij - aji)**2
          end do
        end do
      end do
    end do
    if (nrm > 0.0d0) then
      rel = sqrt(dif) / sqrt(nrm)
    else
      rel = 0.0d0
    end if
  end function hecmw_saamg_is_symmetric

  !> Assemble a block-CSR (nrow x ncol) from a scalar triplet list (i,j,v):
  !! each entry is scattered into its block (block row (i-1)/nb+1, block col
  !! (j-1)/mb+1, in-block position), duplicate entries summed.  Used by the
  !! distributed operator assembly (routed scalar triplets -> owned block rows).
  subroutine hecmw_saamg_bcsr_from_triplets(nrow, ncol, nb, ti, tj, tv, nt, A, mb)
    implicit none
    integer(kind=kint),    intent(in)  :: nrow, ncol, nb, nt
    integer(kind=kint),    intent(in)  :: ti(:), tj(:)
    real(kind=kreal),      intent(in)  :: tv(:)
    type(hecmwST_saamg_bcsr), intent(out) :: A
    integer(kind=kint), optional, intent(in) :: mb   !< column-block size (default nb)
    integer(kind=kint), allocatable :: rcnt(:), rpos(:), jtmp(:), otmp(:), mark(:), touched(:)
    real(kind=kreal),   allocatable :: vtmp(:)
    integer(kind=kint) :: mb_, nbrow, nbcol, bs, t, i, br, rr, bc, cc, j, p, q, ntouch, kout, blk, astat

    mb_ = nb; if (present(mb)) mb_ = mb
    call hecmw_saamg_bcsr_free(A)
    A%n = nrow; A%ncol = ncol; A%nb = nb; A%mb = mb_
    nbrow = nrow / nb; nbcol = ncol / mb_; bs = nb*mb_
    A%nbrow = nbrow; A%nbcol = nbcol

    ! bucket triplets by block row, storing (block col, in-block offset, value)
    allocate(rcnt(nbrow+1)); rcnt = 0
    do t = 1, nt
      br = (ti(t)-1)/nb + 1; rcnt(br+1) = rcnt(br+1) + 1
    end do
    do i = 1, nbrow
      rcnt(i+1) = rcnt(i+1) + rcnt(i)
    end do
    allocate(jtmp(nt), otmp(nt), vtmp(nt), rpos(nbrow), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'bcsr_from_triplets (jtmp/otmp/vtmp)')
    do i = 1, nbrow
      rpos(i) = rcnt(i)
    end do
    do t = 1, nt
      br = (ti(t)-1)/nb + 1; rr = ti(t) - (br-1)*nb
      bc = (tj(t)-1)/mb_ + 1; cc = tj(t) - (bc-1)*mb_
      rpos(br) = rpos(br) + 1
      jtmp(rpos(br)) = bc; otmp(rpos(br)) = (cc-1)*nb + rr; vtmp(rpos(br)) = tv(t)
    end do

    ! count distinct block columns per block row
    allocate(mark(nbcol), touched(nbcol)); mark = 0
    A%nnzb = 0
    do i = 1, nbrow
      ntouch = 0
      do p = rcnt(i)+1, rcnt(i+1)
        if (mark(jtmp(p)) == 0) then
          ntouch = ntouch + 1; touched(ntouch) = jtmp(p); mark(jtmp(p)) = 1
        end if
      end do
      A%nnzb = A%nnzb + ntouch
      do q = 1, ntouch
        mark(touched(q)) = 0
      end do
    end do

    allocate(A%browptr(nbrow+1), A%bcol(A%nnzb), A%bval(bs*A%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'from_triplets (bcol/bval)')
    A%bval = 0.0d0; A%browptr(1) = 1; kout = 0
    do i = 1, nbrow
      ntouch = 0
      do p = rcnt(i)+1, rcnt(i+1)
        j = jtmp(p)
        if (mark(j) == 0) then
          ntouch = ntouch + 1; touched(ntouch) = j; mark(j) = ntouch
          kout = kout + 1; A%bcol(kout) = j
        end if
        blk = kout - ntouch + mark(j)
        A%bval((blk-1)*bs + otmp(p)) = A%bval((blk-1)*bs + otmp(p)) + vtmp(p)
      end do
      do q = 1, ntouch
        mark(touched(q)) = 0
      end do
      A%browptr(i+1) = kout + 1
    end do
    A%nnz = bs*A%nnzb

    deallocate(rcnt, rpos, jtmp, otmp, vtmp, mark, touched)
  end subroutine hecmw_saamg_bcsr_from_triplets

  !> Assemble a block-CSR (nbrow x nbcol block grid, nb x mb blocks) from a list of
  !! block triplets (bi, bj, block); duplicate (bi,bj) blocks are summed.  Each block
  !! bv(:,t) is the nb*mb column-major entries of block t.  The block-native assembler
  !! used by the operator producers (adapter, tentative, smoothed prolongator).
  subroutine hecmw_saamg_bcsr_from_block_triplets(nbrow, nbcol, nb, mb, bi, bj, bv, nt, A)
    implicit none
    integer(kind=kint),    intent(in)  :: nbrow, nbcol, nb, mb, nt
    integer(kind=kint),    intent(in)  :: bi(:), bj(:)
    real(kind=kreal),      intent(in)  :: bv(:,:)   ! (nb*mb, nt)
    type(hecmwST_saamg_bcsr), intent(out) :: A
    integer(kind=kint), allocatable :: rcnt(:), rpos(:), btmp(:), mark(:), touched(:)
    real(kind=kreal),   allocatable :: vtmp(:,:)
    integer(kind=kint) :: t, i, j, p, q, r, ntouch, kout, bs, astat

    bs = nb*mb
    call hecmw_saamg_bcsr_free(A)
    A%n = nbrow*nb; A%ncol = nbcol*mb; A%nb = nb; A%mb = mb
    A%nbrow = nbrow; A%nbcol = nbcol

    ! bucket triplets by block row (unsorted, with duplicates)
    allocate(rcnt(nbrow+1)); rcnt = 0
    do t = 1, nt
      rcnt(bi(t)+1) = rcnt(bi(t)+1) + 1
    end do
    do i = 1, nbrow
      rcnt(i+1) = rcnt(i+1) + rcnt(i)
    end do
    allocate(btmp(nt), vtmp(bs, nt), rpos(nbrow), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'bcsr_from_block_triplets (btmp/vtmp)')
    do i = 1, nbrow
      rpos(i) = rcnt(i)
    end do
    do t = 1, nt
      i = bi(t)
      rpos(i) = rpos(i) + 1
      btmp(rpos(i)) = bj(t)
      vtmp(1:bs, rpos(i)) = bv(1:bs, t)
    end do

    ! count distinct block columns per row
    allocate(mark(nbcol), touched(nbcol)); mark = 0
    A%nnzb = 0
    do i = 1, nbrow
      ntouch = 0
      do p = rcnt(i)+1, rcnt(i+1)
        if (mark(btmp(p)) == 0) then
          ntouch = ntouch + 1; touched(ntouch) = btmp(p); mark(btmp(p)) = 1
        end if
      end do
      A%nnzb = A%nnzb + ntouch
      do q = 1, ntouch
        mark(touched(q)) = 0
      end do
    end do

    allocate(A%browptr(nbrow+1), A%bcol(A%nnzb), A%bval(bs*A%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'from_triplets (bcol/bval)')
    A%bval = 0.0d0; A%browptr(1) = 1; kout = 0
    do i = 1, nbrow
      ntouch = 0
      do p = rcnt(i)+1, rcnt(i+1)
        j = btmp(p)
        if (mark(j) == 0) then
          ntouch = ntouch + 1; touched(ntouch) = j; mark(j) = ntouch
          kout = kout + 1
          A%bcol(kout) = j
          A%bval((kout-1)*bs+1:kout*bs) = vtmp(1:bs, p)
        else
          r = kout - ntouch + mark(j)
          A%bval((r-1)*bs+1:r*bs) = A%bval((r-1)*bs+1:r*bs) + vtmp(1:bs, p)
        end if
      end do
      do q = 1, ntouch
        mark(touched(q)) = 0
      end do
      A%browptr(i+1) = kout + 1
    end do
    A%nnz = bs*A%nnzb

    deallocate(rcnt, rpos, btmp, vtmp, mark, touched)
  end subroutine hecmw_saamg_bcsr_from_block_triplets

  !> Transpose: At = A^T (At is ncol x n).
  subroutine hecmw_saamg_bcsr_transpose(A, At)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_bcsr), intent(out) :: At
    call hecmw_saamg_transpose_blk(A, At)
  end subroutine hecmw_saamg_bcsr_transpose

  !> Sparse matrix-matrix product  C = A * B.  Requires A%mb == B%nb.
  subroutine hecmw_saamg_spgemm(A, B, C)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A, B
    type(hecmwST_saamg_bcsr), intent(out) :: C
    call hecmw_saamg_spgemm_blk(A, B, C)
  end subroutine hecmw_saamg_spgemm

  !> Galerkin coarse operator  Ac = P^T A P.  mblk = coarse block size.
  subroutine hecmw_saamg_galerkin_local(A, P, mblk, Ac)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    integer(kind=kint),    intent(in)  :: mblk
    type(hecmwST_saamg_bcsr), intent(out) :: Ac
    call hecmw_saamg_galerkin_blk(A, P, mblk, Ac)
  end subroutine hecmw_saamg_galerkin_local

  !> Densify a matrix into a (n x ncol) array (verification helper).
  subroutine hecmw_saamg_bcsr_to_dense(A, dense)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    real(kind=kreal),      intent(out) :: dense(:,:)
    call hecmw_saamg_to_dense_blk(A, dense)
  end subroutine hecmw_saamg_bcsr_to_dense

  !-----------------------------------------------------------------------------
  ! block-CSR kernels
  !-----------------------------------------------------------------------------

  !> Densify from the block-CSR storage.
  subroutine hecmw_saamg_to_dense_blk(A, dense)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    real(kind=kreal),      intent(out) :: dense(:,:)
    integer(kind=kint) :: i, t, j, ii, jj, rbase, cbase, boff
    dense = 0.0d0
    do i = 1, A%nbrow
      rbase = (i-1)*A%nb
      do t = A%browptr(i), A%browptr(i+1)-1
        j = A%bcol(t); cbase = (j-1)*A%mb; boff = (t-1)*A%nb*A%mb
        do jj = 1, A%mb
          do ii = 1, A%nb
            dense(rbase+ii, cbase+jj) = dense(rbase+ii, cbase+jj) + A%bval(boff+(jj-1)*A%nb+ii)
          end do
        end do
      end do
    end do
  end subroutine hecmw_saamg_to_dense_blk

  !> Block matvec  y = A x  (x indexed by column blocks, y by row blocks).
  subroutine hecmw_saamg_matvec_blk(A, x, y)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    real(kind=kreal),      intent(in)  :: x(:)
    real(kind=kreal),      intent(out) :: y(:)
    integer(kind=kint) :: i, t, j, ii, jj, rbase, cbase, boff, nb, mb
    real(kind=kreal) :: xj
    nb = A%nb; mb = A%mb
    ! GPU (OpenACC) vs CPU (OpenMP) selected by #ifdef _OPENACC (mutually exclusive,
    ! so a build with both -acc and -mp is safe -- matches upstream las_33.F90).
    ! Arrays are unified (-gpu=mem:managed) so no data clauses; inner block-GEMV loops
    ! stay seq (accumulation, bitwise-stable).
#ifdef _OPENACC
    !$acc parallel loop gang vector private(t,j,ii,jj,rbase,cbase,boff,xj)
#else
    !$omp parallel do default(none) shared(A, x, y, nb, mb) private(i,t,j,ii,jj,rbase,cbase,boff,xj)
#endif
    do i = 1, A%nbrow
      rbase = (i-1)*nb
      y(rbase+1:rbase+nb) = 0.0d0
      !$acc loop seq
      do t = A%browptr(i), A%browptr(i+1)-1
        j = A%bcol(t); cbase = (j-1)*mb; boff = (t-1)*nb*mb
        !$acc loop seq
        do jj = 1, mb
          xj = x(cbase+jj)
          !$acc loop seq
          do ii = 1, nb
            y(rbase+ii) = y(rbase+ii) + A%bval(boff+(jj-1)*nb+ii) * xj
          end do
        end do
      end do
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif
  end subroutine hecmw_saamg_matvec_blk

  !> Block transpose  At = A^T : swap block dims and transpose each block.
  subroutine hecmw_saamg_transpose_blk(A, At)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_bcsr), intent(out) :: At
    integer(kind=kint), allocatable :: pos(:)
    integer(kind=kint) :: i, t, j, ii, jj, nb, mb, p, boff, aoff, astat
    nb = A%nb; mb = A%mb
    call hecmw_saamg_bcsr_free(At)
    At%n = A%ncol; At%ncol = A%n; At%nb = mb; At%mb = nb
    At%nbrow = A%nbcol; At%nbcol = A%nbrow; At%nnzb = A%nnzb; At%nnz = A%nnz
    allocate(At%browptr(At%nbrow+1), At%bcol(At%nnzb), At%bval(nb*mb*At%nnzb), pos(At%nbrow), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'transpose_blk (At)')
    At%browptr = 0
    do t = 1, A%nnzb
      At%browptr(A%bcol(t)+1) = At%browptr(A%bcol(t)+1) + 1
    end do
    At%browptr(1) = 1
    do j = 1, At%nbrow
      At%browptr(j+1) = At%browptr(j+1) + At%browptr(j)
    end do
    do j = 1, At%nbrow
      pos(j) = At%browptr(j)
    end do
    do i = 1, A%nbrow
      do t = A%browptr(i), A%browptr(i+1)-1
        j = A%bcol(t); p = pos(j)
        At%bcol(p) = i
        boff = (p-1)*mb*nb; aoff = (t-1)*nb*mb
        do jj = 1, mb            ! A block (ii,jj) -> At block (jj,ii)
          do ii = 1, nb
            At%bval(boff + (ii-1)*mb + jj) = A%bval(aoff + (jj-1)*nb + ii)
          end do
        end do
        pos(j) = pos(j) + 1
      end do
    end do
    deallocate(pos)
  end subroutine hecmw_saamg_transpose_blk

  !> Block SpGEMM  C = A * B : block Gustavson with dense block GEMM accumulation.
  !! Requires A%mb == B%nb (inner block dimension).
  subroutine hecmw_saamg_spgemm_blk(A, B, C)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A, B
    type(hecmwST_saamg_bcsr), intent(out) :: C
    integer(kind=kint), allocatable :: mark(:), clist(:)
    real(kind=kreal),   allocatable :: spa(:,:)
    integer(kind=kint) :: rb, kb, cb, ncb, i, ta, kk, tb, l, ii, jj, pp, cnt, nci, t, pos
    integer(kind=kint) :: aoff, boff, coff, astat
    real(kind=kreal) :: bkl

    if (A%mb /= B%nb .or. A%nbcol /= B%nbrow) then
      write(*,*) 'hecmw_saamg_spgemm_blk: inner dim mismatch (mb,nbcol)=', A%mb, A%nbcol, &
           ' vs (nb,nbrow)=', B%nb, B%nbrow
      call hecmw_saamg_abort('spgemm_blk: inner dimension mismatch')
    end if
    rb = A%nb; kb = A%mb; cb = B%mb; ncb = B%nbcol
    call hecmw_saamg_bcsr_free(C)
    C%n = A%n; C%ncol = B%ncol; C%nb = rb; C%mb = cb
    C%nbrow = A%nbrow; C%nbcol = ncb
    allocate(C%browptr(A%nbrow+1), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'spgemm_blk (browptr)')

    ! pass 1: count nonzero blocks per row (rows independent; per-thread mark).  The
    ! mark(l)/=i trick avoids re-zeroing: mark(l) holds the last row that touched l
    ! (distinct per row), so it works unchanged with a thread-private mark.
    !$omp parallel default(none) shared(A, B, C, ncb) private(mark, i, cnt, ta, kk, tb, l) if(A%nbrow > 256)
    allocate(mark(ncb)); mark = 0
    !$omp do
    do i = 1, A%nbrow
      cnt = 0
      do ta = A%browptr(i), A%browptr(i+1)-1
        kk = A%bcol(ta)
        do tb = B%browptr(kk), B%browptr(kk+1)-1
          l = B%bcol(tb)
          if (mark(l) /= i) then; mark(l) = i; cnt = cnt + 1; end if
        end do
      end do
      C%browptr(i+1) = cnt
    end do
    !$omp end do
    deallocate(mark)
    !$omp end parallel

    C%browptr(1) = 1                                  ! prefix-sum counts -> row pointers
    do i = 1, A%nbrow
      C%browptr(i+1) = C%browptr(i) + C%browptr(i+1)
    end do
    C%nnzb = C%browptr(A%nbrow+1) - 1
    allocate(C%bcol(C%nnzb), C%bval(rb*cb*C%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'spgemm_blk (C bcol/bval)')
    C%nnz = rb*cb*C%nnzb

    ! pass 2: accumulate block products (rows independent; each row writes its own
    ! browptr(i)..browptr(i+1)-1 slice; per-thread mark/clist/spa accumulator)
    !$omp parallel default(none) shared(A, B, C, rb, kb, cb, ncb) if(A%nbrow > 256) &
    !$omp&   private(mark, clist, spa, i, nci, ta, kk, aoff, tb, l, boff, jj, pp, bkl, ii, t, pos, coff)
    allocate(mark(ncb), clist(ncb), spa(rb*cb, ncb)); mark = 0
    !$omp do
    do i = 1, A%nbrow
      nci = 0
      do ta = A%browptr(i), A%browptr(i+1)-1
        kk = A%bcol(ta); aoff = (ta-1)*rb*kb
        do tb = B%browptr(kk), B%browptr(kk+1)-1
          l = B%bcol(tb); boff = (tb-1)*kb*cb
          if (mark(l) /= i) then
            mark(l) = i; nci = nci + 1; clist(nci) = l; spa(1:rb*cb, l) = 0.0d0
          end if
          do jj = 1, cb
            do pp = 1, kb
              bkl = B%bval(boff + (jj-1)*kb + pp)
              if (bkl == 0.0d0) cycle
              do ii = 1, rb
                spa((jj-1)*rb+ii, l) = spa((jj-1)*rb+ii, l) + A%bval(aoff+(pp-1)*rb+ii) * bkl
              end do
            end do
          end do
        end do
      end do
      pos = C%browptr(i) - 1
      do t = 1, nci
        pos = pos + 1
        l = clist(t)
        C%bcol(pos) = l
        coff = (pos-1)*rb*cb
        C%bval(coff+1:coff+rb*cb) = spa(1:rb*cb, l)
      end do
    end do
    !$omp end do
    deallocate(mark, clist, spa)
    !$omp end parallel
  end subroutine hecmw_saamg_spgemm_blk

  !> Block Galerkin  Ac = P^T A P  (block transpose + two block SpGEMM).
  !! mblk = coarse block size, tagged on Ac%nb / Ac%mb.
  subroutine hecmw_saamg_galerkin_blk(A, P, mblk, Ac)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    integer(kind=kint),    intent(in)  :: mblk
    type(hecmwST_saamg_bcsr), intent(out) :: Ac
    type(hecmwST_saamg_bcsr) :: Pt, C
    call hecmw_saamg_transpose_blk(P, Pt)   ! Pt = P^T
    call hecmw_saamg_spgemm_blk(A, P, C)    ! C  = A P
    call hecmw_saamg_spgemm_blk(Pt, C, Ac)  ! Ac = P^T C
    Ac%nb = mblk; Ac%mb = mblk
    call hecmw_saamg_bcsr_free(Pt)
    call hecmw_saamg_bcsr_free(C)
  end subroutine hecmw_saamg_galerkin_blk

  !> Transpose only the first nbrow_keep block rows of A : At = (A[1:nbrow_keep])^T.
  !! Used to form the restriction P^T from the internal rows of the halo-extended
  !! prolongator without first copying out those rows.
  subroutine hecmw_saamg_transpose_blk_rows(A, nbrow_keep, At)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    integer(kind=kint),    intent(in)  :: nbrow_keep
    type(hecmwST_saamg_bcsr), intent(out) :: At
    integer(kind=kint), allocatable :: pos(:)
    integer(kind=kint) :: i, t, j, ii, jj, nb, mb, p, boff, aoff, nnzb, astat
    nb = A%nb; mb = A%mb
    nnzb = A%browptr(nbrow_keep+1) - 1
    call hecmw_saamg_bcsr_free(At)
    At%n = A%ncol; At%ncol = nbrow_keep*nb; At%nb = mb; At%mb = nb
    At%nbrow = A%nbcol; At%nbcol = nbrow_keep; At%nnzb = nnzb; At%nnz = nnzb*nb*mb
    allocate(At%browptr(At%nbrow+1), At%bcol(nnzb), At%bval(nb*mb*nnzb), pos(At%nbrow), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'transpose_blk_rows (At)')
    At%browptr = 0
    do t = 1, nnzb
      At%browptr(A%bcol(t)+1) = At%browptr(A%bcol(t)+1) + 1
    end do
    At%browptr(1) = 1
    do j = 1, At%nbrow
      At%browptr(j+1) = At%browptr(j+1) + At%browptr(j)
    end do
    do j = 1, At%nbrow
      pos(j) = At%browptr(j)
    end do
    do i = 1, nbrow_keep
      do t = A%browptr(i), A%browptr(i+1)-1
        j = A%bcol(t); p = pos(j)
        At%bcol(p) = i
        boff = (p-1)*mb*nb; aoff = (t-1)*nb*mb
        do jj = 1, mb
          do ii = 1, nb
            At%bval(boff + (ii-1)*mb + jj) = A%bval(aoff + (jj-1)*nb + ii)
          end do
        end do
        pos(j) = pos(j) + 1
      end do
    end do
    deallocate(pos)
  end subroutine hecmw_saamg_transpose_blk_rows

  !> FUSED block triple product  Ac = L * A * R, computed WITHOUT materializing the
  !! intermediate A*R (or L*A): the i -> j(L) -> k(A) -> l(R) loop accumulates the
  !! block products L(i,j) A(j,k) R(k,l) directly into Ac(i,l) (modelled on the
  !! FrontISTR T^t K T product).  Only a small per-(i,j,k) temporary L*A is held, so
  !! the large A*R intermediate is never formed -- critical for the level-1 Galerkin
  !! P^T A P on huge meshes.  Requires L%mb==A%nb, A%mb==R%nb (and matching counts).
  subroutine hecmw_saamg_triple_blk(L, A, R, Ac)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)  :: L, A, R
    type(hecmwST_saamg_bcsr), intent(out) :: Ac
    integer(kind=kint), allocatable :: mark(:), clist(:)
    real(kind=kreal),   allocatable :: spa(:,:), la(:)
    integer(kind=kint) :: rb, k1, k2, cb, ncr, i, jt, j, kt, k, lt, ll, p, q, rr, c
    integer(kind=kint) :: cnt, nci, t, pos, loff, aoff, roff, coff, astat
    real(kind=kreal) :: av, rv

    if (L%mb /= A%nb .or. A%mb /= R%nb .or. L%nbcol /= A%nbrow .or. A%nbcol /= R%nbrow) then
      call hecmw_saamg_abort('triple_blk: dimension mismatch')
    end if
    rb = L%nb; k1 = L%mb; k2 = A%mb; cb = R%mb; ncr = R%nbcol
    call hecmw_saamg_bcsr_free(Ac)
    Ac%n = L%n; Ac%ncol = R%ncol; Ac%nb = rb; Ac%mb = cb
    Ac%nbrow = L%nbrow; Ac%nbcol = ncr
    allocate(Ac%browptr(L%nbrow+1), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'triple_blk (browptr)')

    ! pass 1: count nonzero blocks per row of Ac (rows independent; per-thread mark)
    !$omp parallel default(none) shared(L, A, R, Ac, ncr) private(mark, i, cnt, jt, j, kt, k, lt, ll) if(L%nbrow > 256)
    allocate(mark(ncr)); mark = 0
    !$omp do
    do i = 1, L%nbrow
      cnt = 0
      do jt = L%browptr(i), L%browptr(i+1)-1
        j = L%bcol(jt)
        do kt = A%browptr(j), A%browptr(j+1)-1
          k = A%bcol(kt)
          do lt = R%browptr(k), R%browptr(k+1)-1
            ll = R%bcol(lt)
            if (mark(ll) /= i) then; mark(ll) = i; cnt = cnt + 1; end if
          end do
        end do
      end do
      Ac%browptr(i+1) = cnt
    end do
    !$omp end do
    deallocate(mark)
    !$omp end parallel

    Ac%browptr(1) = 1                                 ! prefix-sum counts -> row pointers
    do i = 1, L%nbrow
      Ac%browptr(i+1) = Ac%browptr(i) + Ac%browptr(i+1)
    end do
    Ac%nnzb = Ac%browptr(L%nbrow+1) - 1
    allocate(Ac%bcol(Ac%nnzb), Ac%bval(rb*cb*Ac%nnzb), stat=astat); Ac%nnz = rb*cb*Ac%nnzb
    call hecmw_saamg_check_alloc(astat, 'triple_blk (Ac bcol/bval)')

    ! pass 2: accumulate L(i,j) A(j,k) R(k,ll) into Ac(i,ll) (rows independent;
    ! each row writes its own slice; per-thread mark/clist/spa/la scratch)
    !$omp parallel default(none) shared(L, A, R, Ac, rb, k1, k2, cb, ncr) if(L%nbrow > 256) &
    !$omp&   private(mark, clist, spa, la, i, nci, jt, j, loff, kt, k, aoff, q, p, av, rr, &
    !$omp&           lt, ll, roff, c, rv, t, pos, coff)
    allocate(mark(ncr), clist(ncr), spa(rb*cb, ncr), la(rb*k2)); mark = 0
    !$omp do
    do i = 1, L%nbrow
      nci = 0
      do jt = L%browptr(i), L%browptr(i+1)-1
        j = L%bcol(jt); loff = (jt-1)*rb*k1
        do kt = A%browptr(j), A%browptr(j+1)-1
          k = A%bcol(kt); aoff = (kt-1)*k1*k2
          ! la = L(i,j) * A(j,k)   (rb x k1) * (k1 x k2) = rb x k2
          la(1:rb*k2) = 0.0d0
          do q = 1, k2
            do p = 1, k1
              av = A%bval(aoff + (q-1)*k1 + p)
              if (av == 0.0d0) cycle
              do rr = 1, rb
                la((q-1)*rb+rr) = la((q-1)*rb+rr) + L%bval(loff+(p-1)*rb+rr) * av
              end do
            end do
          end do
          ! spa(:,ll) += la * R(k,ll)   (rb x k2) * (k2 x cb) = rb x cb
          do lt = R%browptr(k), R%browptr(k+1)-1
            ll = R%bcol(lt); roff = (lt-1)*k2*cb
            if (mark(ll) /= i) then
              mark(ll) = i; nci = nci + 1; clist(nci) = ll; spa(1:rb*cb, ll) = 0.0d0
            end if
            do c = 1, cb
              do q = 1, k2
                rv = R%bval(roff + (c-1)*k2 + q)
                if (rv == 0.0d0) cycle
                do rr = 1, rb
                  spa((c-1)*rb+rr, ll) = spa((c-1)*rb+rr, ll) + la((q-1)*rb+rr) * rv
                end do
              end do
            end do
          end do
        end do
      end do
      pos = Ac%browptr(i) - 1
      do t = 1, nci
        pos = pos + 1
        ll = clist(t); Ac%bcol(pos) = ll
        coff = (pos-1)*rb*cb
        Ac%bval(coff+1:coff+rb*cb) = spa(1:rb*cb, ll)
      end do
    end do
    !$omp end do
    deallocate(mark, clist, spa, la)
    !$omp end parallel
  end subroutine hecmw_saamg_triple_blk

end module hecmw_precond_SAAMG_matrix
