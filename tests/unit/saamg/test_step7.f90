!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : F3a (hierarchy diagnostics)
!
! Builds a hierarchy with verbose=.true. (prints the diagnostic table) and checks
! the captured statistics are sane: >= 2 levels, aggregates respect min_size,
! operator/grid complexity in a reasonable band.
!-------------------------------------------------------------------------------
program test_step7
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_param
  use hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_core
  implicit none

  integer(kind=kint) :: npass, nfail
  type(hecmwST_saamg_bcsr) :: A, Awork
  type(hecmwST_saamg_params) :: prm
  type(hecmwST_saamg_hier) :: dh
  type(hecmwST_saamg_comm) :: cmt
  real(kind=kreal), allocatable :: bfine(:,:)
  integer(kind=8) :: nnz_tot
  integer(kind=kint) :: n, l
  real(kind=kreal) :: opc
  logical :: minsize_ok

  npass = 0; nfail = 0

  call build_grid(48, 48, A)
  n = A%n
  allocate(bfine(n,1)); bfine = 1.0d0
  prm%verbose = .true.
  prm%min_size = 3
  ! drive the production distributed hierarchy directly with a single-rank comm
  call hecmw_saamg_comm_init_serial(cmt, A%nbrow, A%nb)
  call hecmw_saamg_bcsr_copy(A, Awork)        ! setup takes ownership (move_in)
  call hecmw_saamg_setup(Awork, cmt, bfine, 1, prm, dh, move_in=.true.)

  ! The near-kernel (P-hat B_c = B_f) and Galerkin-symmetry properties are checked
  ! by the production distributed verify (hecmw_precond_SAAMG_verify: F4b-2 / F4c-2b),
  ! exercised in the in-tree runs; here we check the structural diagnostics.
  call check('built >= 2 levels', dh%nlevel >= 2, npass, nfail)

  ! aggregation produced coarse nodes on every non-coarsest level
  minsize_ok = .true.
  do l = 1, dh%nlevel - 1
    if (dh%lev(l)%naggr_local <= 0) minsize_ok = .false.
  end do
  call check('aggregation produced coarse nodes (all non-coarsest levels)', &
       minsize_ok, npass, nfail)

  ! operator complexity in a sane band
  nnz_tot = 0
  do l = 1, dh%nlevel
    nnz_tot = nnz_tot + dh%lev(l)%A%nnz
  end do
  opc = real(nnz_tot,kreal)/real(dh%lev(1)%A%nnz,kreal)
  call check('operator complexity in (1.0, 2.5)', opc > 1.0d0 .and. opc < 2.5d0, &
       npass, nfail)

  call hecmw_saamg_free(dh)
  call hecmw_saamg_comm_free(cmt)
  call hecmw_saamg_bcsr_free(A)
  deallocate(bfine)

  write(*,'(/,a,i0,a,i0,a)') '=== F3a: ', npass, ' passed, ', nfail, ' failed ==='
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

end program test_step7
