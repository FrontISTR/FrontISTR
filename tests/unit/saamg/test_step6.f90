!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 6 (coarsest LDL^T, V-cycle, PCG)
!
! Builds the full hierarchy and uses one V-cycle as a CG preconditioner.  The
! key check: on a series-refined 2D grid the SA-AMG-PCG
! iteration count stays ~constant (mesh independence), while plain CG grows.
!-------------------------------------------------------------------------------
program test_step6
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_param
  use hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_core
  implicit none

  integer(kind=kint) :: npass, nfail
  integer(kind=kint) :: sizes(3), it_pre(3), it_un(3), s
  logical :: all_conv
  npass = 0; nfail = 0
  sizes = [16, 32, 64]

  all_conv = .true.
  write(*,'(a)') '   grid     n   levels  PCG-iters  plainCG-iters'
  do s = 1, 3
    call run_one(sizes(s), it_pre(s), it_un(s), all_conv)
  end do

  call check('all SA-AMG-PCG solves converged', all_conv, npass, nfail)
  call check('PCG iters bounded (<= 40)', maxval(it_pre) <= 40, npass, nfail)
  call check('PCG mesh-independent (64 vs 16 within +8)', &
       it_pre(3) - it_pre(1) <= 8, npass, nfail)
  call check('preconditioner clearly helps (64: plain >> PCG)', &
       it_un(3) > 3*it_pre(3), npass, nfail)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 6: ', npass, ' passed, ', nfail, ' failed ==='
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

  subroutine build_grid(nx, ny, A)
    integer(kind=kint), intent(in)  :: nx, ny
    type(hecmwST_saamg_bcsr), intent(out) :: A
    real(kind=kreal), allocatable :: dense(:,:)
    integer(kind=kint) :: ix, iy, id, n
    n = nx*ny
    allocate(dense(n,n)); dense = 0.0d0
    do iy = 1, ny
      do ix = 1, nx
        id = (iy-1)*nx + ix
        dense(id,id) = 4.0d0
        if (ix > 1)  dense(id, id-1)  = -1.0d0
        if (ix < nx) dense(id, id+1)  = -1.0d0
        if (iy > 1)  dense(id, id-nx) = -1.0d0
        if (iy < ny) dense(id, id+nx) = -1.0d0
      end do
    end do
    call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)
    deallocate(dense)
  end subroutine build_grid

  subroutine run_one(ng, iter_pre, iter_un, all_conv)
    integer(kind=kint), intent(in)    :: ng
    integer(kind=kint), intent(out)   :: iter_pre, iter_un
    logical,            intent(inout) :: all_conv
    type(hecmwST_saamg_bcsr) :: A, Awork
    type(hecmwST_saamg_params) :: prm
    type(hecmwST_saamg_hier) :: dh
    type(hecmwST_saamg_comm) :: cmt
    real(kind=kreal), allocatable :: b(:), x(:), bfine(:,:)
    integer(kind=kint) :: n, i
    integer(kind=8) :: seed
    logical :: conv

    call build_grid(ng, ng, A)
    n = A%n
    allocate(bfine(n,1)); bfine = 1.0d0       ! constant near-kernel
    call hecmw_saamg_comm_init_serial(cmt, A%nbrow, A%nb)
    call hecmw_saamg_bcsr_copy(A, Awork)       ! setup takes ownership (move_in)
    call hecmw_saamg_setup(Awork, cmt, bfine, 1, prm, dh, move_in=.true.)

    allocate(b(n), x(n))
    seed = 2024_8
    do i = 1, n
      seed = seed*6364136223846793005_8 + 1442695040888963407_8
      b(i) = real(seed,kreal)/9.223372036854776d18
    end do

    call pcg(A, dh, b, x, 1.0d-8, 300, .true.,  iter_pre, conv)
    if (.not. conv) all_conv = .false.
    call pcg(A, dh, b, x, 1.0d-8, 2000, .false., iter_un, conv)

    write(*,'(i7,i7,i7,i10,i12)') ng, n, dh%nlevel, iter_pre, iter_un

    deallocate(b, x, bfine)
    call hecmw_saamg_free(dh)
    call hecmw_saamg_comm_free(cmt)
    call hecmw_saamg_bcsr_free(A)
  end subroutine run_one

  !> (Preconditioned) CG.  use_pre=.false. -> plain CG.
  subroutine pcg(A, dh, b, x, tol, maxit, use_pre, iter, conv)
    type(hecmwST_saamg_bcsr),       intent(in)    :: A
    type(hecmwST_saamg_hier), intent(inout) :: dh
    real(kind=kreal),            intent(in)    :: b(:), tol
    real(kind=kreal),            intent(out)   :: x(:)
    integer(kind=kint),          intent(in)    :: maxit
    logical,                     intent(in)    :: use_pre
    integer(kind=kint),          intent(out)   :: iter
    logical,                     intent(out)   :: conv
    real(kind=kreal), allocatable :: r(:), z(:), p(:), ap(:)
    real(kind=kreal) :: rz, rz1, alpha, beta, bnrm, pap
    integer(kind=kint) :: n
    n = A%n
    allocate(r(n), z(n), p(n), ap(n))
    x = 0.0d0; r = b
    bnrm = sqrt(dot_product(b,b))
    conv = .false.
    rz1 = 1.0d0; p = 0.0d0
    do iter = 1, maxit
      if (sqrt(dot_product(r,r))/bnrm < tol) then; conv = .true.; exit; end if
      if (use_pre) then
        z = 0.0d0
        call hecmw_saamg_apply(dh, r, z)
      else
        z = r
      end if
      rz = dot_product(r, z)
      if (iter == 1) then
        p = z
      else
        beta = rz / rz1
        p = z + beta * p
      end if
      call hecmw_saamg_matvec(A, p, ap)
      pap = dot_product(p, ap)
      alpha = rz / pap
      x = x + alpha * p
      r = r - alpha * ap
      rz1 = rz
    end do
    if (.not. conv .and. sqrt(dot_product(r,r))/bnrm < tol) conv = .true.
    deallocate(r, z, p, ap)
  end subroutine pcg

end program test_step6
