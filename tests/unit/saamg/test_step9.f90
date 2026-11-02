!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : F5res (numeric-only refresh / setup reuse)
!
! Builds a hierarchy on A, then on a matrix A2 with the SAME sparsity but
! different values: refresh(hier, A2) must reproduce a full rebuild on A2 (same
! aggregation + P-hat are reused, only numeric values recomputed) and still
! precondition A2.
!-------------------------------------------------------------------------------
program test_step9
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_param
  use hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_core
  implicit none

  integer(kind=kint), parameter :: ng = 32
  type(hecmwST_saamg_bcsr) :: A, A2, Awork
  type(hecmwST_saamg_params) :: prm
  type(hecmwST_saamg_hier) :: dh, dhfull
  type(hecmwST_saamg_comm) :: cmt, cmtfull
  real(kind=kreal), allocatable :: bfine(:,:), r(:), z_ref(:), z_rf(:), x(:)
  integer(kind=kint) :: n, npass, nfail, iter
  real(kind=kreal) :: relres
  logical :: conv

  npass = 0; nfail = 0
  call build_grid(ng, 4.0d0, A)        ! diagonal 4
  call build_grid(ng, 4.5d0, A2)       ! same pattern, different values (still SPD)
  n = A%n
  allocate(bfine(n,1)); bfine = 1.0d0

  ! build on A, then refresh to A2 (drive the distributed hierarchy directly)
  call hecmw_saamg_comm_init_serial(cmt, A%nbrow, A%nb)
  call hecmw_saamg_bcsr_copy(A, Awork)
  call hecmw_saamg_setup(Awork, cmt, bfine, 1, prm, dh, move_in=.true.)
  call hecmw_saamg_refresh(dh, A2)

  ! full rebuild on A2
  call hecmw_saamg_comm_init_serial(cmtfull, A2%nbrow, A2%nb)
  call hecmw_saamg_bcsr_copy(A2, Awork)
  call hecmw_saamg_setup(Awork, cmtfull, bfine, 1, prm, dhfull, move_in=.true.)

  ! the two preconditioners must be identical
  allocate(r(n), z_ref(n), z_rf(n), x(n))
  call seeded(r, n)
  z_rf = 0.0d0;  call hecmw_saamg_apply(dh,     r, z_rf)
  z_ref = 0.0d0; call hecmw_saamg_apply(dhfull, r, z_ref)
  call check('refresh == full rebuild (same A2)', &
       maxval(abs(z_rf - z_ref)) < 1.0d-12, npass, nfail)
  write(*,'(a,es10.3)') '   max|z_refresh - z_full| = ', maxval(abs(z_rf - z_ref))

  ! refreshed preconditioner still solves A2
  call pcg(A2, dh, r, x, 1.0d-8, 200, iter, conv, relres)
  call check('refreshed PCG converges on A2', conv, npass, nfail)
  write(*,'(a,i0)') '   refreshed PCG iters on A2 = ', iter

  call hecmw_saamg_free(dh);  call hecmw_saamg_comm_free(cmt)
  call hecmw_saamg_free(dhfull); call hecmw_saamg_comm_free(cmtfull)
  call hecmw_saamg_bcsr_free(A); call hecmw_saamg_bcsr_free(A2)
  deallocate(bfine, r, z_ref, z_rf, x)

  write(*,'(/,a,i0,a,i0,a)') '=== F5res: ', npass, ' passed, ', nfail, ' failed ==='
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

  subroutine build_grid(nxy, diag, A)
    integer(kind=kint), intent(in)  :: nxy
    real(kind=kreal),   intent(in)  :: diag
    type(hecmwST_saamg_bcsr), intent(out) :: A
    real(kind=kreal), allocatable :: d(:,:)
    integer(kind=kint) :: ix, iy, id, n
    n = nxy*nxy
    allocate(d(n,n)); d = 0.0d0
    do iy = 1, nxy
      do ix = 1, nxy
        id = (iy-1)*nxy + ix
        d(id,id) = diag
        if (ix > 1)   d(id, id-1)   = -1.0d0
        if (ix < nxy) d(id, id+1)   = -1.0d0
        if (iy > 1)   d(id, id-nxy) = -1.0d0
        if (iy < nxy) d(id, id+nxy) = -1.0d0
      end do
    end do
    call hecmw_saamg_bcsr_from_dense(d, n, 1, A)
    deallocate(d)
  end subroutine build_grid

  subroutine seeded(r, n)
    real(kind=kreal),   intent(out) :: r(:)
    integer(kind=kint), intent(in)  :: n
    integer(kind=8) :: s
    integer(kind=kint) :: i
    s = 13579_8
    do i = 1, n
      s = s*6364136223846793005_8 + 1442695040888963407_8
      r(i) = real(s,kreal)/9.223372036854776d18
    end do
  end subroutine seeded

  subroutine pcg(A, dh, b, x, tol, maxit, iter, conv, relres)
    type(hecmwST_saamg_bcsr),       intent(in)    :: A
    type(hecmwST_saamg_hier), intent(inout) :: dh
    real(kind=kreal),            intent(in)    :: b(:), tol
    real(kind=kreal),            intent(out)   :: x(:), relres
    integer(kind=kint),          intent(in)    :: maxit
    integer(kind=kint),          intent(out)   :: iter
    logical,                     intent(out)   :: conv
    real(kind=kreal), allocatable :: rr(:), z(:), p(:), ap(:)
    real(kind=kreal) :: rz, rz1, al, be, bn, pap
    integer(kind=kint) :: nn
    nn = A%n; allocate(rr(nn), z(nn), p(nn), ap(nn))
    x = 0.0d0; rr = b; bn = sqrt(dot_product(b,b)); conv = .false.; rz1 = 1.0d0; p = 0.0d0
    do iter = 1, maxit
      relres = sqrt(dot_product(rr,rr))/bn
      if (relres < tol) then; conv = .true.; exit; end if
      z = 0.0d0; call hecmw_saamg_apply(dh, rr, z)
      rz = dot_product(rr, z)
      if (iter == 1) then; p = z; else; be = rz/rz1; p = z + be*p; end if
      call hecmw_saamg_matvec(A, p, ap); pap = dot_product(p, ap); al = rz/pap
      x = x + al*p; rr = rr - al*ap; rz1 = rz
    end do
    deallocate(rr, z, p, ap)
  end subroutine pcg

end program test_step9
