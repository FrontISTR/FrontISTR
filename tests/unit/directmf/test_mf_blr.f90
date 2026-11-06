!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of the BLR compression of hecmw_mf_numeric on synthetic block matrices
!> over grid graphs. Links against the hecmw library and LAPACK; see CMakeLists.txt.
!>
!> Matrix kinds ('smooth' blocks are outer products of per-node vectors, so the fronts are
!> low rank and the compression triggers; 'random' blocks stay near full rank and mostly
!> exercise the full rank fallback):
!>   smoothspd  symmetric positive definite, rank-1 couplings           (LDLt, no pivoting)
!>   smoothsad  smoothspd with zero diagonal multiplier nodes           (LDLt, delays)
!>   smoothlu   structurally symmetric, unsymmetric rank-1 couplings    (LU)
!>   smoothnsd  smoothlu with zero diagonal multiplier nodes            (LU, delays)
!>   randspd    random symmetric positive definite                      (LDLt)
!>   randoff    random unsymmetric, small diagonal                      (LU, row swaps)
!> For each case:
!>   (a) BLR with eps = 1e-30 keeps the full rank solution: bitwise when only zero tiles compress
!>       (rank 0 drops them exactly), else with the residual within 10 times the full rank one
!>       (whether the QR of a rank deficient tile leaves trailing columns below eps times the
!>       leading one depends on the BLAS kernel, and a tile compressed at a positive rank is exact
!>       only to rounding)
!>   (b) the residual of the eps sweep stays below C*eps (C prints with the case)
!>   (c) iterative refinement with the BLR factor reaches 1e-10
!>   (d) a repeated factorization and solve at the largest eps is bitwise reproducible
!> and separately:
!>   (e) the compression of a single tile follows eps, which (a) cannot tell on tiles of exact low
!>       rank: a zero tile compresses at rank 0, and a rank 2 tile with a 1e-12 perturbation
!>       (which keeps every trailing QR column away from zero under any BLAS) stays full rank at
!>       eps = 1e-30 and compresses at rank 2 at eps = 1e-10
program test_mf_blr
  use hecmw_util
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_numeric
  use hecmw_mf_kernel, only: hecmw_mf_kernel_compress
  implicit none

  integer(kind=kint) :: nfail, ntest
  integer(kind=8) :: seed
  real(kind=kreal), parameter :: irtol = 1.0d-10

  nfail = 0
  ntest = 0
  seed = 27182

  call run_case('smoothspd_grid3d_d3_metis_t16', 'smoothspd', 3, 7, 6, 5, 3, 'metis', 0, 16, 0, 1.0d-2, 50.0d0, .true., .true.)
  call run_case('smoothspd_grid2d_d1_id_t8', 'smoothspd', 2, 20, 17, 1, 1, 'identity', 0, 8, 0, 1.0d-2, 50.0d0, .true., .false.)
  call run_case('smoothsad_grid2d_d3_relax_t8', 'smoothsad', 2, 12, 10, 1, 3, 'metis', 5, 8, 2, 1.0d-2, 50.0d0, .true., .false.)
  call run_case('smoothlu_grid3d_d3_metis_t16', 'smoothlu', 3, 7, 6, 5, 3, 'metis', 0, 16, 0, 1.0d-2, 50.0d0, .true., .true.)
  call run_case('smoothlu_grid2d_d3_relax_t8', 'smoothlu', 2, 12, 10, 1, 3, 'metis', 5, 8, 0, 1.0d-2, 50.0d0, .true., .false.)
  call run_case('smoothnsd_grid2d_d1_metis_t8', 'smoothnsd', 2, 12, 10, 1, 1, 'metis', 0, 8, 3, 1.0d-2, 50.0d0, .true., .false.)
  call run_case('randspd_grid3d_d3_metis_t16', 'randspd', 3, 7, 6, 5, 3, 'metis', 0, 16, 0, 1.0d-2, 100.0d0, .false., .false.)
  call run_case('randoff_grid2d_d3_metis_t16', 'randoff', 2, 12, 10, 1, 3, 'metis', 0, 16, 0, 1.0d-6, 5.0d3, .false., .false.)
  call run_compress_eps('compress_eps')

  write(*,'(a,i0,a,i0)') 'test_mf_blr: failed ', nfail, ' of ', ntest
  if (nfail > 0) stop 1

contains

  subroutine run_case(name, kind, dim, nx, ny, nz, ndof, ordering, relax, tile, mstride, eps1, cres, want_lr, clus)
    character(len=*), intent(in) :: name, kind, ordering
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof, relax, tile, mstride
    real(kind=kreal), intent(in) :: eps1, cres
    logical, intent(in) :: want_lr, clus
    real(kind=kreal) :: epss(3)
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    type(hecmwST_mf_factor) :: fct
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint), allocatable :: xadj(:), adjncy(:), perm(:), invp(:), nd(:)
    logical, allocatable :: ismult(:)
    real(kind=kreal), allocatable :: r(:), x0(:), x1(:), dd(:)
    integer(kind=kint) :: n, i, ie, it, ierr
    integer(kind=8) :: nlr
    real(kind=kreal) :: resid0, resid, rir, bnrm
    logical :: ok

    ntest = ntest + 1
    ok = .true.
    epss = (/ eps1, min(eps1, 1.0d-6), 1.0d-10 /)
    call make_graph(dim, nx, ny, nz, mstride, n, xadj, adjncy, ismult)
    allocate(nd(n), perm(n), invp(n))
    nd(1:n) = ndof
    call hecmw_mf_graph_set(graph, n, nd, xadj, adjncy)
    call make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)
    if (ordering == 'metis') then
      call hecmw_ordering_gen(n, n + (xadj(n+1)-1)/2, xadj, adjncy, perm, invp, 0, 0)
    else
      do i = 1, n
        perm(i) = i
      enddo
    endif
    call hecmw_mf_symbolic_build(graph, perm, relax, sym)
    if (clus) then
      call hecmw_mf_symbolic_cluster(graph, sym, perm)
      call hecmw_mf_symbolic_finalize(sym)
      call hecmw_mf_symbolic_build(graph, perm, relax, sym)
    endif
    call hecmw_mf_numeric_init(graph, sym, tile, .true., MF_BLR_ETA, fct)
    allocate(x0(n*ndof), x1(n*ndof), dd(n*ndof), r(n*ndof))
    bnrm = sqrt(dot_product(hecMAT%B, hecMAT%B))

    ! full rank reference
    resid0 = 0.0d0
    fct%blr = .false.
    call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    if (ierr /= 0) then
      write(*,'(a,a,a,i0)') '  ', name, ': full rank factorization failed, ierr = ', ierr
      ok = .false.
    else
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x0)
      call residual(hecMAT, x0, r)
      resid0 = sqrt(dot_product(r, r)) / bnrm
    endif

    ! (a) eps = 1e-30: bitwise the full rank solution unless a tile compressed at a positive
    ! rank, then the residual at the full rank level
    fct%blr = .true.
    fct%eps = 1.0d-30
    if (ok) then
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      if (ierr /= 0) ok = .false.
    endif
    if (ok) then
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x1)
      if (fct%blr_rank_max == 0) then
        if (.not. all(x1 == x0)) then
          write(*,'(a,a,a)') '  ', name, ': eps -> 0 solution is not bitwise the full rank one'
          ok = .false.
        endif
      else
        call residual(hecMAT, x1, r)
        resid = sqrt(dot_product(r, r)) / bnrm
        if (.not. (resid <= 10.0d0*resid0)) then
          write(*,'(a,a,a,1pe12.4,a,i0)') '  ', name, ': eps -> 0 residual = ', resid, &
            ' with tiles of rank up to ', fct%blr_rank_max
          ok = .false.
        endif
      endif
    endif

    ! (b) eps sweep: residual below cres*eps (never above the full rank level by much)
    nlr = 0
    resid = 0.0d0
    do ie = 1, size(epss)
      if (.not. ok) exit
      fct%eps = epss(ie)
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      if (ierr /= 0) then
        write(*,'(a,a,a,1pe9.2,a,i0)') '  ', name, ': eps = ', epss(ie), ' factorization failed, ierr = ', ierr
        ok = .false.
        exit
      endif
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x1)
      call residual(hecMAT, x1, r)
      resid = sqrt(dot_product(r, r)) / bnrm
      if (.not. (resid <= max(cres*epss(ie), 10.0d0*resid0))) then
        write(*,'(a,a,a,1pe9.2,a,1pe12.4)') '  ', name, ': eps = ', epss(ie), ' residual = ', resid
        ok = .false.
      endif
      if (ie == 1) nlr = fct%blr_tiles_lr
    enddo
    if (ok .and. want_lr .and. nlr == 0) then
      write(*,'(a,a,a)') '  ', name, ': no tile was compressed at the largest eps'
      ok = .false.
    endif
    if (ok) then
      if (fct%lu .neqv. (kind == 'smoothlu' .or. kind == 'smoothnsd' .or. kind == 'randoff')) then
        write(*,'(a,a,a,l1)') '  ', name, ': unexpected mode, lu = ', fct%lu
        ok = .false.
      endif
    endif

    ! (c) iterative refinement at the largest eps reaches irtol
    if (ok) then
      fct%eps = epss(1)
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x1)
      rir = 0.0d0
      do it = 1, 40
        call residual(hecMAT, x1, r)
        rir = sqrt(dot_product(r, r)) / bnrm
        if (rir <= irtol) exit
        call hecmw_mf_numeric_solve(sym, fct, r, dd)
        x1 = x1 + dd
      enddo
      if (.not. (rir <= irtol)) then
        write(*,'(a,a,a,1pe12.4)') '  ', name, ': refinement stalled at residual = ', rir
        ok = .false.
      endif
    endif

    ! (d) bitwise reproducibility of the BLR factorization and solve, then of a repeated
    ! factorization with the rank reuse enabled (a tile skipped from the previous full rank
    ! outcome takes the same path the repeated compression would have taken)
    if (ok) then
      fct%eps = epss(1)
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x1)
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, dd)
      if (.not. all(dd == x1)) then
        write(*,'(a,a,a)') '  ', name, ': repeated BLR run is not bitwise reproducible'
        ok = .false.
      endif
      fct%blr_reuse = .true.
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, dd)
      fct%blr_reuse = .false.
      if (.not. all(dd == x1)) then
        write(*,'(a,a,a)') '  ', name, ': BLR rerun with rank reuse changed the solution'
        ok = .false.
      endif
    endif

    if (ok) then
      write(*,'(a,a,a,i0,a,l1,a,i0,a,i0,a,i0,a,1pe10.2,a,1pe10.2)') 'PASS ', name, ': ndof_tot=', n*ndof, &
        ' lu=', fct%lu, ' delayed=', fct%n_delay, ' lr_tiles=', nlr, ' of ', fct%blr_tiles, &
        ' resid(1e-2)=', resid, ' fr_resid=', resid0
    else
      write(*,'(a,a)') 'FAIL ', name
      nfail = nfail + 1
    endif
    call hecmw_mf_numeric_finalize(fct)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(hecMAT%D, hecMAT%AL, hecMAT%AU, hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
    deallocate(hecMAT%B, hecMAT%X)
    deallocate(xadj, adjncy, perm, invp, nd, ismult)
    deallocate(x0, x1, dd, r)
  end subroutine run_case

  !> (e) the compression of a single tile follows eps
  subroutine run_compress_eps(name)
    character(len=*), intent(in) :: name
    real(kind=kreal) :: a0(12,8), a(12,8), v(8,8)
    integer(kind=kint) :: i, j, rk(3)

    ntest = ntest + 1
    do j = 1, 8
      do i = 1, 12
        a0(i,j) = 1.0d0 + real(i, kind=kreal)/12.0d0*cos(real(j, kind=kreal)) + 1.0d-12*(2.0d0*rand01() - 1.0d0)
      enddo
    enddo
    a(:,:) = 0.0d0
    call hecmw_mf_kernel_compress(12, 8, 12, a, 1.0d-30, 8, v, rk(1))
    a(:,:) = a0(:,:)
    call hecmw_mf_kernel_compress(12, 8, 12, a, 1.0d-30, 8, v, rk(2))
    a(:,:) = a0(:,:)
    call hecmw_mf_kernel_compress(12, 8, 12, a, 1.0d-10, 8, v, rk(3))
    if (rk(1) == 0 .and. rk(2) == -1 .and. rk(3) == 2) then
      write(*,'(a,a)') 'PASS ', name
    else
      write(*,'(a,a,a,i0,a,i0,a,i0)') 'FAIL ', name, ': rank of the zero tile ', rk(1), ', at eps = 1e-30 ', rk(2), &
        ', at eps = 1e-10 ', rk(3)
      nfail = nfail + 1
    endif
  end subroutine run_compress_eps

  !> Grid graph (chain, 2D 9-point, 3D 27-point); with mstride > 0 every mstride-th grid node
  !> gets a multiplier node attached by a single edge (ismult).
  subroutine make_graph(dim, nx, ny, nz, mstride, n, xadj, adjncy, ismult)
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, mstride
    integer(kind=kint), intent(out) :: n
    integer(kind=kint), allocatable, intent(out) :: xadj(:), adjncy(:)
    logical, allocatable, intent(out) :: ismult(:)
    integer(kind=kint), allocatable :: mult(:)
    integer(kind=kint) :: ix, iy, iz, jx, jy, jz, i, j, cnt, my, mz, pass, ng, nm

    my = ny
    mz = nz
    if (dim < 2) my = 1
    if (dim < 3) mz = 1
    ng = nx*my*mz
    allocate(mult(ng))
    mult(1:ng) = 0
    nm = 0
    if (mstride > 0) then
      do i = 1, ng, mstride
        nm = nm + 1
        mult(i) = ng + nm
      enddo
    endif
    n = ng + nm
    allocate(xadj(n+1), ismult(n))
    ismult(1:ng) = .false.
    ismult(ng+1:n) = .true.
    do pass = 1, 2
      cnt = 0
      xadj(1) = 1
      do iz = 1, mz
        do iy = 1, my
          do ix = 1, nx
            i = ix + nx*((iy-1) + my*(iz-1))
            do jz = max(1,iz-1), min(mz,iz+1)
              do jy = max(1,iy-1), min(my,iy+1)
                do jx = max(1,ix-1), min(nx,ix+1)
                  j = jx + nx*((jy-1) + my*(jz-1))
                  if (j == i) cycle
                  cnt = cnt + 1
                  if (pass == 2) adjncy(cnt) = j
                enddo
              enddo
            enddo
            if (mult(i) > 0) then
              cnt = cnt + 1
              if (pass == 2) adjncy(cnt) = mult(i)
            endif
            xadj(i+1) = cnt + 1
          enddo
        enddo
      enddo
      do i = 1, ng
        if (mult(i) == 0) cycle
        cnt = cnt + 1
        if (pass == 2) adjncy(cnt) = i
        xadj(mult(i)+1) = cnt + 1
      enddo
      if (pass == 1) allocate(adjncy(cnt))
    enddo
    deallocate(mult)
  end subroutine make_graph

  !> Block matrix in the row major block layout of hecmwST_matrix. The 'smooth' kinds build
  !> every coupling block as an outer product of per-node vectors, so that the frontal
  !> matrices carry low rank off-diagonal tiles; the 'rand' kinds use random blocks. Diagonal
  !> blocks are diagonally dominant (SPD by Gershgorin) except at multiplier nodes (zero) and
  !> for randoff (small random diagonal).
  subroutine make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: xadj(:), adjncy(:)
    logical, intent(in) :: ismult(:)
    character(len=*), intent(in) :: kind
    type(hecmwST_matrix), intent(out) :: hecMAT
    integer(kind=kint) :: i, j, p, nl, nu, a, b, ku, kl, nd2
    integer(kind=kint), allocatable :: upos(:)
    real(kind=kreal), allocatable :: rowsum(:), uvec(:,:), vvec(:,:)
    real(kind=kreal) :: v
    logical :: smooth, lu

    smooth = (kind(1:6) == 'smooth')
    lu = (kind == 'smoothlu' .or. kind == 'smoothnsd' .or. kind == 'randoff')
    nd2 = ndof*ndof
    hecMAT%N = n
    hecMAT%NP = n
    hecMAT%NDOF = ndof
    hecMAT%symmetric = .not. lu
    allocate(hecMAT%indexL(0:n), hecMAT%indexU(0:n))
    allocate(hecMAT%itemL(xadj(n+1)-1), hecMAT%itemU(xadj(n+1)-1))
    hecMAT%indexL(0) = 0
    hecMAT%indexU(0) = 0
    nl = 0
    nu = 0
    do i = 1, n
      do p = xadj(i), xadj(i+1)-1
        j = adjncy(p)
        if (j < i) then
          nl = nl + 1
          hecMAT%itemL(nl) = j
        else
          nu = nu + 1
          hecMAT%itemU(nu) = j
        endif
      enddo
      hecMAT%indexL(i) = nl
      hecMAT%indexU(i) = nu
    enddo
    hecMAT%NPL = nl
    hecMAT%NPU = nu
    allocate(hecMAT%D(n*nd2), hecMAT%AL(nl*nd2), hecMAT%AU(nu*nd2), hecMAT%B(n*ndof), hecMAT%X(n*ndof))
    allocate(rowsum(n*ndof), upos(n), uvec(ndof, n), vvec(ndof, n))
    rowsum(:) = 0.0d0
    do i = 1, n
      do a = 1, ndof
        uvec(a, i) = 2.0d0*rand01() - 1.0d0
        vvec(a, i) = 2.0d0*rand01() - 1.0d0
      enddo
    enddo

    do i = 1, n
      do ku = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(ku)
        do b = 1, ndof
          do a = 1, ndof
            if (smooth .and. .not. (ismult(i) .or. ismult(j))) then
              ! rank-1 coupling; multiplier couplings stay full rank (else singular)
              v = uvec(a, i) * vvec(b, j)
            else
              v = 2.0d0*rand01() - 1.0d0
            endif
            hecMAT%AU((ku-1)*nd2 + (a-1)*ndof + b) = v
          enddo
        enddo
      enddo
    enddo
    upos(1:n) = hecMAT%indexU(0:n-1) + 1
    do j = 1, n
      do kl = hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i = hecMAT%itemL(kl)
        ku = upos(i)
        upos(i) = upos(i) + 1
        do a = 1, ndof
          do b = 1, ndof
            hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b) = hecMAT%AU((ku-1)*nd2 + (b-1)*ndof + a)
          enddo
        enddo
      enddo
    enddo
    if (lu) then
      ! unsymmetric values: an independent outer product (or random blocks) in the lower part
      do j = 1, n
        do kl = hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
          i = hecMAT%itemL(kl)
          do b = 1, ndof
            do a = 1, ndof
              if (smooth .and. .not. (ismult(i) .or. ismult(j))) then
                v = uvec(a, j) * vvec(b, i)
              else
                v = 2.0d0*rand01() - 1.0d0
              endif
              hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b) = v
            enddo
          enddo
        enddo
      enddo
    endif

    do i = 1, n
      do kl = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        do a = 1, ndof
          do b = 1, ndof
            rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b))
          enddo
        enddo
      enddo
      do ku = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        do a = 1, ndof
          do b = 1, ndof
            rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(hecMAT%AU((ku-1)*nd2 + (a-1)*ndof + b))
          enddo
        enddo
      enddo
      do a = 1, ndof
        do b = a+1, ndof
          v = 2.0d0*rand01() - 1.0d0
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = v
          hecMAT%D((i-1)*nd2 + (b-1)*ndof + a) = v
          rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(v)
          rowsum((i-1)*ndof + b) = rowsum((i-1)*ndof + b) + abs(v)
        enddo
      enddo
      do a = 1, ndof
        hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = rowsum((i-1)*ndof + a) + 1.0d0 + rand01()
      enddo
      if (kind == 'randoff') then
        do a = 1, ndof
          do b = 1, ndof
            if (a /= b) hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = 2.0d0*rand01() - 1.0d0
          enddo
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = 0.1d0*(2.0d0*rand01() - 1.0d0)
        enddo
      endif
      if ((kind == 'smoothsad' .or. kind == 'smoothnsd') .and. ismult(i)) then
        hecMAT%D((i-1)*nd2 + 1:i*nd2) = 0.0d0
      endif
    enddo

    do i = 1, n*ndof
      hecMAT%B(i) = 2.0d0*rand01() - 1.0d0
      hecMAT%X(i) = 0.0d0
    enddo
    deallocate(rowsum, upos, uvec, vvec)
  end subroutine make_matrix

  !> r = B - K X with the full (D, AL, AU) block storage.
  subroutine residual(hecMAT, x, r)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), intent(in) :: x(:)
    real(kind=kreal), intent(out) :: r(:)
    integer(kind=kint) :: i, j, k, a, b, ndof, nd2

    ndof = hecMAT%NDOF
    nd2 = ndof*ndof
    r(:) = hecMAT%B(:)
    do i = 1, hecMAT%NP
      do a = 1, ndof
        do b = 1, ndof
          r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%D((i-1)*nd2 + (a-1)*ndof + b)*x((i-1)*ndof + b)
        enddo
      enddo
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        j = hecMAT%itemL(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AL((k-1)*nd2 + (a-1)*ndof + b)*x((j-1)*ndof + b)
          enddo
        enddo
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AU((k-1)*nd2 + (a-1)*ndof + b)*x((j-1)*ndof + b)
          enddo
        enddo
      enddo
    enddo
  end subroutine residual

  function rand01() result(v)
    real(kind=kreal) :: v

    seed = mod(1103515245_8*seed + 12345_8, 2147483648_8)
    v = real(seed, kind=kreal) / 2147483648.0d0
  end function rand01

end program test_mf_blr
