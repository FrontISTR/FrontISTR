!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of hecmw_mf_numeric on synthetic SPD block matrices assembled as
!> hecmwST_matrix over grid graphs. Links against the hecmw library; see CMakeLists.txt.
!>
!> Each case builds the symbolic structure (relax 0 or >0), the numeric factor with a given
!> tile size, solves for a random right hand side, and checks the relative residual and that the
!> stored factor entries match the symbolic estimate.
program test_mf_numeric
  use hecmw_util
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_numeric
  implicit none

  integer(kind=kint) :: nfail, ntest
  integer(kind=8) :: seed
  real(kind=kreal), parameter :: tol = 1.0d-12

  nfail = 0
  ntest = 0
  seed = 4711

  call run_case('grid2d_d1_id_t8',       2, 12, 10, 1, 1, 'identity', 0, 8)
  call run_case('grid2d_d1_metis_t4096', 2, 12, 10, 1, 1, 'metis',    0, 4096)
  call run_case('grid2d_d3_metis_t8',    2, 20, 17, 1, 3, 'metis',    0, 8)
  call run_case('grid2d_d3_relax_t8',    2, 20, 17, 1, 3, 'metis',    5, 8)
  call run_case('grid2d_d3_relax_t4096', 2, 20, 17, 1, 3, 'metis',    5, 4096)
  call run_case('grid3d_d1_rand_t16',    3,  7,  6, 5, 1, 'random',   0, 16)
  call run_case('grid3d_d3_metis_t16',   3,  9,  8, 7, 3, 'metis',    0, 16)
  call run_case('grid3d_d3_relax_t16',   3,  9,  8, 7, 3, 'metis',    5, 16)
  call run_case('grid3d_d3_relax_t4096', 3,  9,  8, 7, 3, 'metis',    5, 4096)
  call run_case('chain_d2_relax_t8',     1, 50,  1, 1, 2, 'identity', 8, 8)
  call run_case('single_node_d3',        1,  1,  1, 1, 3, 'identity', 0, 8)

  write(*,'(a,i0,a,i0)') 'test_mf_numeric: failed ', nfail, ' of ', ntest
  if (nfail > 0) stop 1

contains

  subroutine run_case(name, dim, nx, ny, nz, ndof, ordering, relax, tile)
    character(len=*), intent(in) :: name, ordering
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof, relax, tile
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    type(hecmwST_mf_factor) :: fct
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint), allocatable :: xadj(:), adjncy(:), perm(:), invp(:), nd(:)
    real(kind=kreal), allocatable :: r(:)
    integer(kind=kint) :: n, i, ierr
    real(kind=kreal) :: resid
    logical :: ok

    ntest = ntest + 1
    call make_grid(dim, nx, ny, nz, n, xadj, adjncy)
    allocate(nd(n), perm(n), invp(n))
    nd(1:n) = ndof
    call hecmw_mf_graph_set(graph, n, nd, xadj, adjncy)
    call make_matrix(n, ndof, xadj, adjncy, hecMAT)

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
    ok = .true.
    if (fct%factor_nnz /= sym%factor_nnz) then
      write(*,'(a,a,a,i0,a,i0)') '  ', name, ': stored factor entries ', fct%factor_nnz, ' vs symbolic ', sym%factor_nnz
      ok = .false.
    endif

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
      if (resid > tol) then
        write(*,'(a,a,a,1pe12.4)') '  ', name, ': relative residual = ', resid
        ok = .false.
      endif
    endif

    if (ok) then
      write(*,'(a,a,a,i0,a,i0,a,i0,a,1pe10.2)') 'PASS ', name, ': ndof_tot=', n*ndof, ' nsuper=', sym%nsuper, &
        ' tiles/front(max)=', fct%max_tiles, ' resid=', resid
    else
      write(*,'(a,a)') 'FAIL ', name
      nfail = nfail + 1
    endif

    call hecmw_mf_numeric_finalize(fct)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(hecMAT%D, hecMAT%AL, hecMAT%AU, hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
    deallocate(hecMAT%B, hecMAT%X)
    deallocate(xadj, adjncy, perm, invp, nd)
  end subroutine run_case

  subroutine make_grid(dim, nx, ny, nz, n, xadj, adjncy)
    integer(kind=kint), intent(in) :: dim, nx, ny, nz
    integer(kind=kint), intent(out) :: n
    integer(kind=kint), allocatable, intent(out) :: xadj(:), adjncy(:)
    integer(kind=kint) :: ix, iy, iz, jx, jy, jz, i, j, cnt, my, mz, pass

    my = ny
    mz = nz
    if (dim < 2) my = 1
    if (dim < 3) mz = 1
    n = nx*my*mz
    allocate(xadj(n+1))
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
            xadj(i+1) = cnt + 1
          enddo
        enddo
      enddo
      if (pass == 1) allocate(adjncy(cnt))
    enddo
  end subroutine make_grid

  !> Symmetric block matrix with random off-diagonal blocks and diagonally dominant diagonal
  !> blocks (SPD by Gershgorin), stored in the row major block layout of hecmwST_matrix.
  subroutine make_matrix(n, ndof, xadj, adjncy, hecMAT)
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: xadj(:), adjncy(:)
    type(hecmwST_matrix), intent(out) :: hecMAT
    integer(kind=kint) :: i, j, p, nl, nu, a, b, ku, kl, nd2
    integer(kind=kint), allocatable :: upos(:)
    real(kind=kreal), allocatable :: rowsum(:)
    real(kind=kreal) :: v

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

    ! upper blocks random in (-1,1); the lower block of (j,i) is the transpose of the upper block of (i,j)
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
    do i = 1, n*ndof
      hecMAT%B(i) = 2.0d0*rand01() - 1.0d0
      hecMAT%X(i) = 0.0d0
    enddo
    deallocate(rowsum, upos)
  end subroutine make_matrix

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

end program test_mf_numeric
