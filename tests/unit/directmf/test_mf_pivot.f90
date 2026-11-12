!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of the pivoting of hecmw_mf_numeric on synthetic symmetric indefinite block
!> matrices assembled as hecmwST_matrix over grid graphs. Links against the hecmw library
!> and LAPACK; see CMakeLists.txt.
!>
!> Each case builds an SPD matrix as test_mf_numeric does and then makes it indefinite:
!>   signflip  random signs on the diagonal entries (diagonally dominant, pivots exchanged)
!>   saddle    multiplier nodes with zero diagonal block attached to grid nodes (2x2 and delays)
!>   shift     K - s I with s inside the spectrum
!>   zerodiag  every diagonal entry zero (delays propagate to the root)
!>   spd       unchanged (no exchange, no delay expected)
!> The factorization must give a relative residual below 1e-10 (1e-9 for zerodiag, whose
!> factors reach the growth bound 1/u of the threshold) and the inertia of dsytrf.
program test_mf_pivot
  use hecmw_util
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_numeric
  implicit none

  integer(kind=kint) :: nfail, ntest
  integer(kind=8) :: seed
  real(kind=kreal), parameter :: tol = 1.0d-10

  nfail = 0
  ntest = 0
  seed = 2718

  call run_case('spd_grid3d_d3_metis_t16',     'spd',      3,  7,  6, 5, 3, 'metis',    0, 16, 0)
  call run_case('spd_grid2d_d3_relax_t4096',   'spd',      2, 20, 17, 1, 3, 'metis',    5, 4096, 0)
  call run_case('signflip_grid2d_d1_id_t8',    'signflip', 2, 12, 10, 1, 1, 'identity', 0, 8, 0)
  call run_case('signflip_grid3d_d3_metis_t16','signflip', 3,  7,  6, 5, 3, 'metis',    0, 16, 0)
  call run_case('signflip_grid2d_d3_relax_t4096','signflip',2, 20, 17, 1, 3, 'metis',   5, 4096, 0)
  call run_case('saddle_grid2d_d1_metis_t8',   'saddle',   2, 12, 10, 1, 1, 'metis',    0, 8, 3)
  call run_case('saddle_grid2d_d3_relax_t16',  'saddle',   2, 12, 10, 1, 3, 'metis',    5, 16, 2)
  call run_case('saddle_grid3d_d3_metis_t4096','saddle',   3,  6,  5, 5, 3, 'metis',    0, 4096, 4)
  call run_case('saddle_chain_d2_relax_t8',    'saddle',   1, 40,  1, 1, 2, 'identity', 8, 8, 2)
  call run_case('shift_grid2d_d3_metis_t8',    'shift',    2, 12, 10, 1, 3, 'metis',    0, 8, 0)
  call run_case('shift_grid3d_d1_rand_t16',    'shift',    3,  7,  6, 5, 1, 'random',   0, 16, 0)
  call run_case('shift_grid3d_d3_relax_t4096', 'shift',    3,  6,  5, 5, 3, 'metis',    5, 4096, 0)
  call run_case('zerodiag_chain_d1_id_t8',     'zerodiag', 1,  6,  1, 1, 1, 'identity', 0, 8, 0)
  call run_case('zerodiag_grid2d_d1_id_t4096', 'zerodiag', 2,  4,  3, 1, 1, 'identity', 0, 4096, 0)
  call run_case('zerodiag_grid2d_d1_id_t8',    'zerodiag', 2,  8,  7, 1, 1, 'identity', 0, 8, 0)
  call run_case('zerodiag_grid2d_d3_metis_t8', 'zerodiag', 2, 12, 10, 1, 3, 'metis',    0, 8, 0)
  call run_case('zerodiag_grid3d_d3_relax_t4096','zerodiag',3, 6,  5, 5, 3, 'metis',    5, 4096, 0)

  write(*,'(a,i0,a,i0)') 'test_mf_pivot: failed ', nfail, ' of ', ntest
  if (nfail > 0) stop 1

contains

  subroutine run_case(name, kind, dim, nx, ny, nz, ndof, ordering, relax, tile, mstride)
    character(len=*), intent(in) :: name, kind, ordering
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof, relax, tile, mstride
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    type(hecmwST_mf_factor) :: fct
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint), allocatable :: xadj(:), adjncy(:), perm(:), invp(:), nd(:)
    logical, allocatable :: ismult(:)
    real(kind=kreal), allocatable :: r(:)
    integer(kind=kint) :: n, i, ierr, npos_ref, nneg_ref
    real(kind=kreal) :: resid, resid_ref
    logical :: ok

    ntest = ntest + 1
    call make_graph(dim, nx, ny, nz, mstride, n, xadj, adjncy, ismult)
    allocate(nd(n), perm(n), invp(n))
    nd(1:n) = ndof
    call hecmw_mf_graph_set(graph, n, nd, xadj, adjncy)
    call make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)

    select case (ordering)
    case ('identity')
      do i = 1, n
        perm(i) = i
      enddo
    case ('random')
      call shuffle(n, perm)
    case ('metis')
      call hecmw_ordering_gen(n, n + (xadj(n+1)-1)/2, xadj, adjncy, perm, invp, 0, 0)
    end select

    call hecmw_mf_symbolic_build(graph, perm, relax, sym)
    call hecmw_mf_numeric_init(graph, sym, tile, .false., 0, fct)
    call inertia_ref(hecMAT, npos_ref, nneg_ref, resid_ref)
    ok = .true.
    resid = 0.0d0

    call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    if (ierr /= 0) then
      write(*,'(a,a,a,i0)') '  ', name, ': factorization failed, ierr = ', ierr
      ok = .false.
    else
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, hecMAT%X)
      allocate(r(n*ndof))
      call residual(hecMAT, r)
      resid = sqrt(dot_product(r, r)) / sqrt(dot_product(hecMAT%B, hecMAT%B))
      deallocate(r)
      if (.not. (resid <= tol) .and. .not. (kind == 'zerodiag' .and. resid <= 10.0d0*tol)) then
        write(*,'(a,a,a,1pe12.4)') '  ', name, ': relative residual = ', resid
        ok = .false.
      endif
      if (fct%n_pos /= npos_ref .or. fct%n_neg /= nneg_ref) then
        write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a)') '  ', name, ': inertia (', fct%n_pos, ',', fct%n_neg, &
          ') vs dsytrf (', npos_ref, ',', nneg_ref, ')'
        ok = .false.
      endif
      if (kind == 'spd' .and. (fct%n_delay /= 0 .or. fct%n_swap /= 0 .or. fct%n_2x2 /= 0)) then
        write(*,'(a,a,a)') '  ', name, ': pivoting happened on an SPD matrix'
        ok = .false.
      endif
      if (kind == 'spd' .and. fct%factor_words_act /= fct%factor_words) then
        write(*,'(a,a,a)') '  ', name, ': factor size differs from the estimate on an SPD matrix'
        ok = .false.
      endif
    endif

    if (ok) then
      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe10.2,a,1pe10.2,a)') 'PASS ', name, ': ndof_tot=', n*ndof, &
        ' nsuper=', sym%nsuper, ' neg=', fct%n_neg, ' 2x2=', fct%n_2x2, ' swaps=', fct%n_swap, &
        ' delayed=', fct%n_delay, ' growth=', fct%max_growth, ' resid=', resid, ' (dsytrs ', resid_ref, ')'
    else
      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,1pe10.2,a)') 'FAIL ', name, ': nsuper=', sym%nsuper, ' 2x2=', fct%n_2x2, &
        ' swaps=', fct%n_swap, ' delayed=', fct%n_delay, ' growth=', fct%max_growth, ' (dsytrs resid ', resid_ref, ')'
      nfail = nfail + 1
    endif

    call hecmw_mf_numeric_finalize(fct)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(hecMAT%D, hecMAT%AL, hecMAT%AU, hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
    deallocate(hecMAT%B, hecMAT%X)
    deallocate(xadj, adjncy, perm, invp, nd, ismult)
  end subroutine run_case

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

  !> Symmetric block matrix with random off-diagonal blocks and diagonally dominant diagonal
  !> blocks (SPD by Gershgorin) in the row major block layout of hecmwST_matrix, then modified
  !> according to kind (see the program header).
  subroutine make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: xadj(:), adjncy(:)
    logical, intent(in) :: ismult(:)
    character(len=*), intent(in) :: kind
    type(hecmwST_matrix), intent(out) :: hecMAT
    integer(kind=kint) :: i, j, p, nl, nu, a, b, ku, kl, nd2
    integer(kind=kint), allocatable :: upos(:)
    real(kind=kreal), allocatable :: rowsum(:)
    real(kind=kreal) :: v, s

    nd2 = ndof*ndof
    hecMAT%N = n
    hecMAT%NP = n
    hecMAT%NDOF = ndof
    hecMAT%symmetric = .true.
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
    allocate(rowsum(n*ndof), upos(n))
    rowsum(:) = 0.0d0

    do ku = 1, nu
      do a = 1, nd2
        hecMAT%AU((ku-1)*nd2 + a) = 2.0d0*rand01() - 1.0d0
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
    enddo

    select case (kind)
    case ('signflip')
      do i = 1, n
        do a = 1, ndof
          if (rand01() < 0.5d0) hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = -hecMAT%D((i-1)*nd2 + (a-1)*ndof + a)
        enddo
      enddo
    case ('saddle')
      do i = 1, n
        if (ismult(i)) hecMAT%D((i-1)*nd2 + 1:i*nd2) = 0.0d0
      enddo
    case ('shift')
      s = 0.0d0
      do i = 1, n
        do a = 1, ndof
          s = s + hecMAT%D((i-1)*nd2 + (a-1)*ndof + a)
        enddo
      enddo
      s = s / (n*ndof)
      do i = 1, n
        do a = 1, ndof
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) - s
        enddo
      enddo
    case ('zerodiag')
      do i = 1, n
        do a = 1, ndof
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = 0.0d0
        enddo
      enddo
    end select

    do i = 1, n*ndof
      hecMAT%B(i) = 2.0d0*rand01() - 1.0d0
      hecMAT%X(i) = 0.0d0
    enddo
    deallocate(rowsum, upos)
  end subroutine make_matrix

  !> Inertia of the dense matrix by the Bunch-Kaufman factorization of LAPACK (dsytrf, 'L') and
  !> the relative residual of the solution by dsytrs for the same right hand side.
  subroutine inertia_ref(hecMAT, npos, nneg, resid)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: npos, nneg
    real(kind=kreal), intent(out) :: resid
    real(kind=kreal), allocatable :: a(:,:), work(:), r(:)
    integer(kind=kint), allocatable :: ipiv(:)
    integer(kind=kint) :: n, ndof, nd2, i, j, k, ia, ib, info, lwork
    real(kind=kreal) :: wq(1), d11, d21, d22, det
    external :: dsytrf, dsytrs

    ndof = hecMAT%NDOF
    nd2 = ndof*ndof
    n = hecMAT%NP*ndof
    allocate(a(n,n), ipiv(n))
    a(:,:) = 0.0d0
    do i = 1, hecMAT%NP
      do ia = 1, ndof
        do ib = 1, ndof
          a((i-1)*ndof + ia, (i-1)*ndof + ib) = hecMAT%D((i-1)*nd2 + (ia-1)*ndof + ib)
        enddo
      enddo
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        j = hecMAT%itemL(k)
        do ia = 1, ndof
          do ib = 1, ndof
            a((i-1)*ndof + ia, (j-1)*ndof + ib) = hecMAT%AL((k-1)*nd2 + (ia-1)*ndof + ib)
          enddo
        enddo
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        do ia = 1, ndof
          do ib = 1, ndof
            a((i-1)*ndof + ia, (j-1)*ndof + ib) = hecMAT%AU((k-1)*nd2 + (ia-1)*ndof + ib)
          enddo
        enddo
      enddo
    enddo
    call dsytrf('L', n, a, n, ipiv, wq, -1, info)
    lwork = max(int(wq(1)), 1)
    allocate(work(lwork))
    call dsytrf('L', n, a, n, ipiv, work, lwork, info)
    deallocate(work)
    if (info /= 0) then
      write(*,'(a,i0)') '  dsytrf info = ', info
    endif
    hecMAT%X(1:n) = hecMAT%B(1:n)
    call dsytrs('L', n, 1, a, n, ipiv, hecMAT%X, n, info)
    allocate(r(n))
    call residual(hecMAT, r)
    resid = sqrt(dot_product(r, r)) / sqrt(dot_product(hecMAT%B, hecMAT%B))
    deallocate(r)
    hecMAT%X(1:n) = 0.0d0
    npos = 0
    nneg = 0
    k = 1
    do while (k <= n)
      if (ipiv(k) > 0) then
        if (a(k,k) > 0.0d0) then
          npos = npos + 1
        else
          nneg = nneg + 1
        endif
        k = k + 1
      else
        d11 = a(k,k)
        d21 = a(k+1,k)
        d22 = a(k+1,k+1)
        det = d11*d22 - d21*d21
        if (det < 0.0d0) then
          npos = npos + 1
          nneg = nneg + 1
        else if (d11 + d22 > 0.0d0) then
          npos = npos + 2
        else
          nneg = nneg + 2
        endif
        k = k + 2
      endif
    enddo
    deallocate(a, ipiv)
  end subroutine inertia_ref

  !> r = B - K X with the full (D, AL, AU) block storage.
  subroutine residual(hecMAT, r)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), intent(out) :: r(:)
    integer(kind=kint) :: i, j, k, a, b, ndof, nd2

    ndof = hecMAT%NDOF
    nd2 = ndof*ndof
    r(:) = hecMAT%B(:)
    do i = 1, hecMAT%NP
      do a = 1, ndof
        do b = 1, ndof
          r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%D((i-1)*nd2 + (a-1)*ndof + b)*hecMAT%X((i-1)*ndof + b)
        enddo
      enddo
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        j = hecMAT%itemL(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AL((k-1)*nd2 + (a-1)*ndof + b)*hecMAT%X((j-1)*ndof + b)
          enddo
        enddo
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AU((k-1)*nd2 + (a-1)*ndof + b)*hecMAT%X((j-1)*ndof + b)
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

  subroutine shuffle(n, perm)
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(out) :: perm(:)
    integer(kind=kint) :: i, j, t

    do i = 1, n
      perm(i) = i
    enddo
    do i = n, 2, -1
      j = 1 + int(rand01()*i, kind=kint)
      if (j > i) j = i
      t = perm(i)
      perm(i) = perm(j)
      perm(j) = t
    enddo
  end subroutine shuffle

end program test_mf_pivot
