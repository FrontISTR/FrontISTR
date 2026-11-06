!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of the LU mode of hecmw_mf_numeric on synthetic block matrices with a
!> symmetric structure and unsymmetric values, assembled as hecmwST_matrix over grid graphs.
!> Links against the hecmw library and LAPACK; see CMakeLists.txt.
!>
!> Each case builds an SPD matrix as test_mf_pivot does and then modifies it:
!>   perturb   small random unsymmetric perturbation of the upper part (few or no exchanges)
!>   offdom    independent random lower and upper parts, small diagonal (row exchanges)
!>   nsaddle   unsymmetric couplings and multiplier nodes with zero diagonal block (delays)
!>   zerodiag  unsymmetric couplings and every diagonal entry zero (delays reach the root)
!>   symlu     unchanged (symmetric), solved in LDLt mode and, with the symmetric flag off,
!>             in LU mode; both solutions must agree
!> The mode chosen by the factorization must match the kind, the relative residual must be
!> below 1e-10 (1e-9 for zerodiag, as in test_mf_pivot) and is compared with dgetrs.
program test_mf_lu
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
  seed = 3141

  call run_case('perturb_grid2d_d1_id_t8',     'perturb',  2, 12, 10, 1, 1, 'identity', 0, 8, 0)
  call run_case('perturb_grid3d_d3_metis_t16', 'perturb',  3,  7,  6, 5, 3, 'metis',    0, 16, 0)
  call run_case('perturb_grid2d_d3_relax_t4096','perturb', 2, 20, 17, 1, 3, 'metis',    5, 4096, 0)
  call run_case('offdom_grid2d_d1_id_t8',      'offdom',   2, 12, 10, 1, 1, 'identity', 0, 8, 0)
  call run_case('offdom_grid3d_d3_metis_t16',  'offdom',   3,  7,  6, 5, 3, 'metis',    0, 16, 0)
  call run_case('offdom_grid3d_d1_rand_t4096', 'offdom',   3,  7,  6, 5, 1, 'random',   0, 4096, 0)
  call run_case('offdom_grid2d_d3_relax_t8',   'offdom',   2, 20, 17, 1, 3, 'metis',    5, 8, 0)
  call run_case('nsaddle_grid2d_d1_metis_t8',  'nsaddle',  2, 12, 10, 1, 1, 'metis',    0, 8, 3)
  call run_case('nsaddle_grid2d_d3_relax_t16', 'nsaddle',  2, 12, 10, 1, 3, 'metis',    5, 16, 2)
  call run_case('nsaddle_grid3d_d3_metis_t4096','nsaddle', 3,  6,  5, 5, 3, 'metis',    0, 4096, 4)
  call run_case('nsaddle_chain_d2_relax_t8',   'nsaddle',  1, 40,  1, 1, 2, 'identity', 8, 8, 2)
  call run_case('zerodiag_chain_d1_id_t8',     'zerodiag', 1,  6,  1, 1, 1, 'identity', 0, 8, 0)
  call run_case('zerodiag_grid2d_d1_id_t4096', 'zerodiag', 2,  4,  3, 1, 1, 'identity', 0, 4096, 0)
  call run_case('zerodiag_grid2d_d3_metis_t8', 'zerodiag', 2, 12, 10, 1, 3, 'metis',    0, 8, 0)
  call run_case('zerodiag_grid3d_d3_relax_t4096','zerodiag',3, 6,  5, 5, 3, 'metis',    5, 4096, 0)
  call run_case('symlu_grid3d_d3_metis_t16',   'symlu',    3,  7,  6, 5, 3, 'metis',    0, 16, 0)
  call run_case('symlu_grid2d_d3_relax_t4096', 'symlu',    2, 20, 17, 1, 3, 'metis',    5, 4096, 0)

  write(*,'(a,i0,a,i0)') 'test_mf_lu: failed ', nfail, ' of ', ntest
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
    real(kind=kreal), allocatable :: r(:), xsym(:)
    integer(kind=kint) :: n, i, ierr
    real(kind=kreal) :: resid, resid_ref, diff
    logical :: ok, want_lu

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
    ! the asymmetry scan is log output only; enable it so the check below can verify that
    ! the generated unsymmetric matrices are numerically unsymmetric
    fct%scan = .true.
    call residual_ref(hecMAT, resid_ref)
    ok = .true.
    resid = 0.0d0
    diff = 0.0d0
    want_lu = (kind /= 'symlu')

    if (kind == 'symlu') then
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
      if (ierr /= 0 .or. fct%lu) then
        write(*,'(a,a,a,i0,a,l1)') '  ', name, ': LDLt run failed, ierr = ', ierr, ', lu = ', fct%lu
        ok = .false.
      else
        allocate(xsym(n*ndof))
        call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, xsym)
      endif
      hecMAT%symmetric = .false.
    endif

    if (ok) then
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
        if (.not. fct%lu) then
          write(*,'(a,a,a,1pe12.4)') '  ', name, ': LDLt mode was chosen, asymmetry = ', fct%asym
          ok = .false.
        endif
        if (want_lu .and. fct%asym <= 1.0d-12) then
          write(*,'(a,a,a,1pe12.4)') '  ', name, ': asymmetry not detected, ', fct%asym
          ok = .false.
        endif
        if (kind == 'symlu') then
          diff = maxval(abs(hecMAT%X - xsym)) / maxval(abs(xsym))
          if (.not. (diff <= tol)) then
            write(*,'(a,a,a,1pe12.4)') '  ', name, ': LU and LDLt solutions differ by ', diff
            ok = .false.
          endif
          ! the estimate counts one face; the LU mode stores the L and U faces
          if (fct%n_swap /= 0 .or. fct%n_delay /= 0 .or. fct%factor_words_act /= 2*fct%factor_words) then
            write(*,'(a,a,a)') '  ', name, ': pivoting happened or the factor size differs from the estimate'
            ok = .false.
          endif
        endif
      endif
    endif

    if (ok) then
      write(*,'(a,a,a,i0,a,i0,a,1pe9.2,a,i0,a,i0,a,i0,a,1pe10.2,a,1pe10.2,a)') 'PASS ', name, ': ndof_tot=', n*ndof, &
        ' nsuper=', sym%nsuper, ' asym=', fct%asym, ' swaps=', fct%n_swap, ' delayed=', fct%n_delay, &
        ' growth=', fct%max_growth, ' resid=', resid, ' (dgetrs ', resid_ref, ')'
    else
      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,1pe10.2,a)') 'FAIL ', name, ': nsuper=', sym%nsuper, ' swaps=', fct%n_swap, &
        ' delayed=', fct%n_delay, ' growth=', fct%max_growth, ' (dgetrs resid ', resid_ref, ')'
      nfail = nfail + 1
    endif

    call hecmw_mf_numeric_finalize(fct)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    if (allocated(xsym)) deallocate(xsym)
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
    real(kind=kreal) :: v

    nd2 = ndof*ndof
    hecMAT%N = n
    hecMAT%NP = n
    hecMAT%NDOF = ndof
    hecMAT%symmetric = (kind == 'symlu')
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
    case ('perturb')
      do ku = 1, nu
        do a = 1, nd2
          hecMAT%AU((ku-1)*nd2 + a) = hecMAT%AU((ku-1)*nd2 + a) + 1.0d-3*(2.0d0*rand01() - 1.0d0)
        enddo
      enddo
      do i = 1, n
        do a = 1, ndof
          do b = a+1, ndof
            hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) + 1.0d-3*(2.0d0*rand01() - 1.0d0)
          enddo
        enddo
      enddo
    case ('offdom', 'nsaddle', 'zerodiag')
      do kl = 1, nl
        do a = 1, nd2
          hecMAT%AL((kl-1)*nd2 + a) = 2.0d0*rand01() - 1.0d0
        enddo
      enddo
      do i = 1, n
        do a = 1, ndof
          do b = 1, ndof
            if (a == b) cycle
            hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = 2.0d0*rand01() - 1.0d0
          enddo
        enddo
      enddo
      do i = 1, n
        do a = 1, ndof
          if (kind == 'offdom') then
            hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = 0.1d0*(2.0d0*rand01() - 1.0d0)
          else if (kind == 'zerodiag') then
            hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = 0.0d0
          endif
        enddo
        if (kind == 'nsaddle' .and. ismult(i)) hecMAT%D((i-1)*nd2 + 1:i*nd2) = 0.0d0
      enddo
    end select

    do i = 1, n*ndof
      hecMAT%B(i) = 2.0d0*rand01() - 1.0d0
      hecMAT%X(i) = 0.0d0
    enddo
    deallocate(rowsum, upos)
  end subroutine make_matrix

  !> Relative residual of the solution of the dense matrix by dgetrf/dgetrs for the same right
  !> hand side.
  subroutine residual_ref(hecMAT, resid)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(out) :: resid
    real(kind=kreal), allocatable :: a(:,:), r(:)
    integer(kind=kint), allocatable :: ipiv(:)
    integer(kind=kint) :: n, ndof, nd2, i, j, k, ia, ib, info
    external :: dgetrf, dgetrs

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
    call dgetrf(n, n, a, n, ipiv, info)
    if (info /= 0) then
      write(*,'(a,i0)') '  dgetrf info = ', info
    endif
    hecMAT%X(1:n) = hecMAT%B(1:n)
    call dgetrs('N', n, 1, a, n, ipiv, hecMAT%X, n, info)
    allocate(r(n))
    call residual(hecMAT, r)
    resid = sqrt(dot_product(r, r)) / sqrt(dot_product(hecMAT%B, hecMAT%B))
    deallocate(r)
    hecMAT%X(1:n) = 0.0d0
    deallocate(a, ipiv)
  end subroutine residual_ref

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

end program test_mf_lu
