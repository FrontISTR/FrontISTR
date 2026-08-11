!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : distributed self-checks
!!
!! Diagnostic routines run only under the `verify` option at setup.  Each rebuilds
!! a stage of the distributed pipeline (halo matvec / aggregation+coarse comm table
!! / tentative & smoothed prolongator / numeric refresh) and prints a global,
!! partition-invariant quantity, so a 1-rank and a K-rank run can be compared to
!! confirm the MPI construction.  Separated from the production backend so the
!! preconditioner lifecycle stays free of diagnostic code.
module hecmw_precond_SAAMG_verify
  use hecmw_util
  use hecmw_solver_misc,             only: hecmw_innerProduct_R
  use hecmw_precond_SAAMG_adapt
  use hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_core
  use hecmw_precond_SAAMG_param,     only: hecmwST_saamg_params
  implicit none

  private
  public :: hecmw_saamg_verify_matvec
  public :: hecmw_saamg_verify_coarsen
  public :: hecmw_saamg_verify_prolong
  public :: hecmw_saamg_verify_smoothp
  public :: hecmw_saamg_verify_refresh

contains

  !> Verify the distributed numeric refresh: refreshing with the SAME operator A
  !! must reproduce the freshly-built preconditioner (refresh reuses aggregation /
  !! coarse comm table / tentative P-hat, recomputes the value-dependent parts).
  !! Applies M^{-1} to a partition-invariant seed before and after refresh and
  !! reports the global relative difference, which must be ~machine precision.
  subroutine hecmw_saamg_verify_refresh(dh, A, cmt, B, m, prm, n)
    implicit none
    type(hecmwST_saamg_hier), intent(inout) :: dh
    type(hecmwST_saamg_bcsr),      intent(in)    :: A
    type(hecmwST_saamg_comm),      intent(in)    :: cmt
    real(kind=kreal),            intent(in)    :: B(:,:)
    integer(kind=kint),          intent(in)    :: m
    type(hecmwST_saamg_params),    intent(in)    :: prm
    integer(kind=kint),          intent(in)    :: n
    type(hecmwST_saamg_bcsr) :: Ap
    type(hecmwST_saamg_hier) :: dh2
    real(kind=kreal), allocatable :: r(:), z1(:), z2(:)
    real(kind=kreal) :: d2, n2, rel_same, rel_pert
    integer(kind=kint) :: i, node, dof, gid, ir, k, rn, cn, nb1, e
    real(kind=kreal) :: fr, fc

    ! The finest-level matvec hook (if any) points at the REAL hecMAT operator,
    ! which the perturbation Ap cannot reach.  Disable it for the duration of this
    ! self-check so every apply uses dh's own (refreshable) level operators -- the
    ! comparison must be refresh-vs-build of the SAME operator, not hook-vs-Ap.
    dh%hook_disabled = .true.

    nb1 = dh%lev(1)%nb_l
    allocate(r(n), z1(n), z2(n))
    do i = 1, n
      node = (i-1)/nb1 + 1; dof = mod(i-1, nb1) + 1
      gid = node; if (allocated(cmt%gnode)) gid = cmt%gnode(node)
      r(i) = sin(0.17d0*real(gid, kreal) + 0.31d0*real(dof, kreal))
    end do

    ! (a) same A: refresh must reproduce the full build, bit-for-bit
    call hecmw_saamg_apply(dh, r, z1)
    call hecmw_saamg_refresh(dh, A)
    call hecmw_saamg_apply(dh, r, z2)
    d2 = 0.0d0; n2 = 0.0d0
    do i = 1, n; d2 = d2 + (z1(i)-z2(i))**2; n2 = n2 + z1(i)**2; end do
    call hecmw_saamg_comm_allreduce_sum_r(cmt, d2)
    call hecmw_saamg_comm_allreduce_sum_r(cmt, n2)
    rel_same = sqrt(d2/max(n2,1.0d-300))

    ! (b) value perturbation A' = Df A Df (symmetric congruence: same pattern, still
    ! SPD).  A refresh from A' must match a FRESH full build from A' (S4 routing
    ! re-sends only values).  Df(node) = 1 + 0.3*sin(gid) > 0.
    call hecmw_saamg_bcsr_copy(A, Ap)
    do ir = 1, Ap%nbrow                         ! block row = node
      gid = ir; if (allocated(cmt%gnode)) gid = cmt%gnode(ir)
      fr = 1.0d0 + 0.3d0*sin(0.5d0*real(gid,kreal))
      do k = Ap%browptr(ir), Ap%browptr(ir+1)-1
        cn = Ap%bcol(k)
        gid = cn; if (allocated(cmt%gnode)) gid = cmt%gnode(cn)
        fc = 1.0d0 + 0.3d0*sin(0.5d0*real(gid,kreal))
        do e = 1, Ap%nb*Ap%mb                   ! scale the whole nb x mb block
          Ap%bval((k-1)*Ap%nb*Ap%mb+e) = Ap%bval((k-1)*Ap%nb*Ap%mb+e) * fr * fc
        end do
      end do
    end do
    call hecmw_saamg_setup(Ap, cmt, B, m, prm, dh2)   ! fresh reference build
    call hecmw_saamg_refresh(dh, Ap)                  ! refresh existing hierarchy
    call hecmw_saamg_apply(dh2, r, z1)
    call hecmw_saamg_apply(dh,  r, z2)
    d2 = 0.0d0; n2 = 0.0d0
    do i = 1, n; d2 = d2 + (z1(i)-z2(i))**2; n2 = n2 + z1(i)**2; end do
    call hecmw_saamg_comm_allreduce_sum_r(cmt, d2)
    call hecmw_saamg_comm_allreduce_sum_r(cmt, n2)
    rel_pert = sqrt(d2/max(n2,1.0d-300))
    call hecmw_saamg_free(dh2)
    call hecmw_saamg_refresh(dh, A)   ! restore the hierarchy to the real A
    dh%hook_disabled = .false.             ! re-enable the finest-level matvec hook

    if (cmt%my_rank == 0) then
      write(*,'(a,es10.3)') &
        '#### SA-AMG refresh verify (same A):    ||z_refresh - z_build||/||z|| = ', rel_same
      write(*,'(a,es10.3)') &
        '#### SA-AMG refresh verify (perturb A): ||z_refresh - z_rebuild||/||z|| = ', rel_pert
    end if
    call hecmw_saamg_bcsr_free(Ap)
    deallocate(r, z1, z2)
  end subroutine hecmw_saamg_verify_refresh

  !> F4a self-check: the distributed (halo-aware) matvec applied to a globally
  !! consistent vector x (defined from global node ids) must give the same global
  !! ||A x||^2 regardless of how the mesh is partitioned.  Compare the printed
  !! value between a 1-rank run and a K-rank run: equality proves the comm table
  !! and halo exchange are correct (dropping halo columns would change it).
  subroutine hecmw_saamg_verify_matvec(hecMAT, hecMESH, cmt, ndof)
    implicit none
    type(hecmwST_matrix),     intent(in) :: hecMAT
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_saamg_comm),   intent(in) :: cmt
    integer(kind=kint),       intent(in) :: ndof
    type(hecmwST_saamg_bcsr) :: Afull
    real(kind=kreal), allocatable :: x(:), y(:)
    real(kind=kreal) :: nax2, nx2
    integer(kind=kint) :: ni, np, i, d, gid

    call hecmw_saamg_from_hecmat(hecMAT, Afull, include_halo=.true.)
    ni = hecMESH%nn_internal; np = hecMESH%n_node
    allocate(x(np*ndof), y(ni*ndof)); x = 0.0d0
    do i = 1, ni
      gid = hecMESH%global_node_ID(i)
      do d = 1, ndof
        x(ndof*(i-1)+d) = sin(0.1d0*real(gid,kreal) + 0.37d0*real(d,kreal))
      end do
    end do
    call hecmw_saamg_matvec_d(cmt, Afull, x, y)            ! halo exchange + local product
    call hecmw_innerProduct_R(hecMESH, ndof, y, y, nax2)   ! internal-only sum + allreduce
    call hecmw_innerProduct_R(hecMESH, ndof, x, x, nx2)
    if (hecMESH%my_rank == 0) write(*,'(a,es22.15,a,es22.15)') &
         '#### SA-AMG F4a verify: global ||A x||^2= ', nax2, '  ||x||^2= ', nx2
    deallocate(x, y)
    call hecmw_saamg_bcsr_free(Afull)
  end subroutine hecmw_saamg_verify_matvec

  !> F4b self-check: build the uncoupled aggregation and the coarse comm table on
  !! the distributed finest operator, verify the coarse table is symmetric (each
  !! neighbor's coarse export to me equals my coarse import from it), and report
  !! the global coarse-node count.  The global coarse-node count is
  !! partition-independent-ish (uncoupled aggregation differs only at boundaries);
  !! the symmetry flag must be true on every rank.
  subroutine hecmw_saamg_verify_coarsen(hecMAT, cmt, m, prm)
    implicit none
    type(hecmwST_matrix),     intent(in) :: hecMAT
    type(hecmwST_saamg_comm),   intent(in) :: cmt
    integer(kind=kint),       intent(in) :: m
    type(hecmwST_saamg_params), intent(in) :: prm
    type(hecmwST_saamg_bcsr) :: Afull
    type(hecmwST_saamg_comm) :: cmt_c
    integer(kind=kint), allocatable :: aggr(:), agg_loc(:), buf(:)
    integer(kind=kint) :: naggr_local, ncnode, nprocs, tot_aggr, tot_halo, okint
    logical :: ok

    call hecmw_saamg_from_hecmat(hecMAT, Afull, include_halo=.true.)
    call hecmw_saamg_coarsen_struct(Afull, cmt, prm%min_size, prm%max_size, &
         prm%theta, m, aggr, naggr_local, agg_loc, ncnode, cmt_c, agg_order=prm%agg_order)
    call hecmw_saamg_verify_commtable(cmt_c, ok)

    nprocs = hecmw_saamg_comm_size(cmt)
    allocate(buf(nprocs))
    call hecmw_saamg_comm_allgather_int(cmt, naggr_local,        buf); tot_aggr = sum(buf(1:nprocs))
    call hecmw_saamg_comm_allgather_int(cmt, ncnode-naggr_local, buf); tot_halo = sum(buf(1:nprocs))
    okint = 0; if (ok) okint = 1
    call hecmw_saamg_comm_allgather_int(cmt, okint, buf)
    ok = (minval(buf(1:nprocs)) == 1)
    if (cmt%my_rank == 0) write(*,'(a,i0,a,i0,a,l1)') &
         '#### SA-AMG F4b verify: total coarse nodes= ', tot_aggr, &
         '  coarse-halo refs= ', tot_halo, '  commtable symmetric= ', ok

    deallocate(buf, aggr, agg_loc)
    call hecmw_saamg_comm_free(cmt_c)
    call hecmw_saamg_bcsr_free(Afull)
  end subroutine hecmw_saamg_verify_coarsen

  !> F4b-2 self-check: build the halo-extended tentative prolongator P-hat_ext on
  !! the distributed operator and confirm the exact QR identity P-hat_ext * B_c =
  !! B_fine holds on ALL rows (internal AND halo).  The halo rows exercise the
  !! P-hat block exchange and the coarse near-kernel exchange (coarse comm table).
  !! The near-kernel is built with zero centroid so it is globally consistent
  !! across ranks (owner and neighbor agree on a shared node's B_fine).
  subroutine hecmw_saamg_verify_prolong(hecMAT, hecMESH, cmt, ndof, m, prm)
    implicit none
    type(hecmwST_matrix),     intent(in) :: hecMAT
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_saamg_comm),   intent(in) :: cmt
    integer(kind=kint),       intent(in) :: ndof, m
    type(hecmwST_saamg_params), intent(in) :: prm
    type(hecmwST_saamg_bcsr) :: Afull, phat_ext
    type(hecmwST_saamg_comm) :: cmt_c
    integer(kind=kint), allocatable :: aggr(:), agg_loc(:)
    real(kind=kreal),   allocatable :: Bext(:,:), bcoarse(:,:), bc(:), y(:)
    integer(kind=kint) :: naggr_local, ncnode, np, ni, nrow, i, base, c, j
    real(kind=kreal) :: xx, yy, zz, resid

    call hecmw_saamg_from_hecmat(hecMAT, Afull, include_halo=.true.)
    call hecmw_saamg_coarsen_struct(Afull, cmt, prm%min_size, prm%max_size, &
         prm%theta, m, aggr, naggr_local, agg_loc, ncnode, cmt_c, agg_order=prm%agg_order)

    np = hecMESH%n_node; ni = hecMESH%nn_internal; nrow = np*ndof
    ! rigid-body near-kernel over ALL nodes, ZERO centroid (rank-consistent)
    allocate(Bext(nrow, m)); Bext = 0.0d0
    do i = 1, np
      xx = hecMESH%node(3*i-2); yy = hecMESH%node(3*i-1); zz = hecMESH%node(3*i)
      base = ndof*(i-1)
      Bext(base+1,1)=1.0d0; Bext(base+1,5)= zz; Bext(base+1,6)=-yy
      Bext(base+2,2)=1.0d0; Bext(base+2,4)=-zz; Bext(base+2,6)= xx
      Bext(base+3,3)=1.0d0; Bext(base+3,4)= yy; Bext(base+3,5)=-xx
      if (ndof == 6) then
        Bext(base+4,4)=1.0d0; Bext(base+5,5)=1.0d0; Bext(base+6,6)=1.0d0
      end if
    end do

    call hecmw_saamg_tentative_ext(cmt, ndof, m, Bext(1:ni*ndof, :), aggr, &
         naggr_local, agg_loc, ncnode, phat_ext, bcoarse)

    allocate(bc(ncnode*m), y(nrow))
    resid = 0.0d0
    do c = 1, m
      bc = 0.0d0
      do j = 1, naggr_local*m
        bc(j) = bcoarse(j, c)
      end do
      call hecmw_saamg_comm_update(cmt_c, m, bc)     ! fill halo coarse rows from owners
      call hecmw_saamg_matvec(phat_ext, bc, y)       ! y = P-hat_ext * bc (all rows)
      do i = 1, nrow
        resid = max(resid, abs(y(i) - Bext(i, c)))
      end do
    end do
    call hecmw_saamg_comm_allreduce_max_r(cmt, resid)
    if (cmt%my_rank == 0) write(*,'(a,es10.3)') &
         '#### SA-AMG F4b-2 verify: max |P-hat_ext B_c - B_fine| (all rows)= ', resid

    deallocate(aggr, agg_loc, Bext, bcoarse, bc, y)
    call hecmw_saamg_comm_free(cmt_c)
    call hecmw_saamg_bcsr_free(phat_ext)
    call hecmw_saamg_bcsr_free(Afull)
  end subroutine hecmw_saamg_verify_prolong

  !> F4c-1 self-check: build the distributed SMOOTHED prolongator P and confirm
  !! the global near-kernel residual ||P B_c - B_f||^2 (over internal rows) is
  !! partition-independent.  Smoothing does not preserve B exactly (P B_c = B_f -
  !! omega D^{-1} A B_f), so the residual is nonzero, but it is the SAME global
  !! value on 1 / K ranks iff the smoothed P, the shared lambda_max (omega) and
  !! the coarse near-kernel exchange are all correct.  Compare the printed value
  !! between a 1-rank and a K-rank run.
  subroutine hecmw_saamg_verify_smoothp(hecMAT, hecMESH, cmt, ndof, m, prm)
    implicit none
    type(hecmwST_matrix),     intent(in) :: hecMAT
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_saamg_comm),   intent(in) :: cmt
    integer(kind=kint),       intent(in) :: ndof, m
    type(hecmwST_saamg_params), intent(in) :: prm
    type(hecmwST_saamg_bcsr) :: Afull, P, Ac_global
    type(hecmwST_saamg_comm) :: cmt_c
    integer(kind=kint), allocatable :: aggr(:), agg_loc(:), chalo_gid(:)
    real(kind=kreal),   allocatable :: Bext(:,:), bcoarse(:,:), bc(:), y(:)
    integer(kind=kint) :: naggr_local, ncnode, np, ni, nir, i, base, c, j, my_off
    real(kind=kreal) :: xx, yy, zz, omega, res2, prows_mism, asym

    call hecmw_saamg_from_hecmat(hecMAT, Afull, include_halo=.true.)
    call hecmw_saamg_coarsen_struct(Afull, cmt, prm%min_size, prm%max_size, &
         prm%theta, m, aggr, naggr_local, agg_loc, ncnode, cmt_c, my_off, chalo_gid, &
         agg_order=prm%agg_order)

    np = hecMESH%n_node; ni = hecMESH%nn_internal; nir = ni*ndof
    allocate(Bext(np*ndof, m)); Bext = 0.0d0
    do i = 1, np
      xx = hecMESH%node(3*i-2); yy = hecMESH%node(3*i-1); zz = hecMESH%node(3*i)
      base = ndof*(i-1)
      Bext(base+1,1)=1.0d0; Bext(base+1,5)= zz; Bext(base+1,6)=-yy
      Bext(base+2,2)=1.0d0; Bext(base+2,4)=-zz; Bext(base+2,6)= xx
      Bext(base+3,3)=1.0d0; Bext(base+3,4)= yy; Bext(base+3,5)=-xx
      if (ndof == 6) then
        Bext(base+4,4)=1.0d0; Bext(base+5,5)=1.0d0; Bext(base+6,6)=1.0d0
      end if
    end do

    call hecmw_saamg_prolongator(Afull, cmt, m, prm%safety, prm%lanczos_iter, &
         Bext(1:nir, :), aggr, naggr_local, agg_loc, ncnode, P, bcoarse, omega)

    allocate(bc(ncnode*m), y(nir))
    res2 = 0.0d0
    do c = 1, m
      bc = 0.0d0
      do j = 1, naggr_local*m
        bc(j) = bcoarse(j, c)
      end do
      call hecmw_saamg_comm_update(cmt_c, m, bc)
      call hecmw_saamg_matvec(P, bc, y)              ! y = P bc (internal rows)
      do i = 1, nir
        res2 = res2 + (y(i) - Bext(i, c))**2
      end do
    end do
    call hecmw_saamg_comm_allreduce_sum_r(cmt, res2)
    if (cmt%my_rank == 0) write(*,'(a,es22.15,a,es12.5)') &
         '#### SA-AMG F4c-1 verify: global ||P B_c - B_f||^2= ', res2, '  omega= ', omega

    ! F4c-2a: verify the smoothed-P halo row exchange (round-trip checksum)
    call hecmw_saamg_verify_prows(cmt, m, P, naggr_local, my_off, chalo_gid, prows_mism)
    if (cmt%my_rank == 0) write(*,'(a,es10.3)') &
         '#### SA-AMG F4c-2a verify: max P halo-row exchange mismatch= ', prows_mism

    ! F4c-2b: build the global coarse operator (redundant assembly) and check it
    ! is symmetric -- the partition-invariant correctness test of the Galerkin.
    call hecmw_saamg_galerkin_global(Afull, cmt, P, naggr_local, ncnode, &
         my_off, chalo_gid, m, Ac_global)
    asym = hecmw_saamg_is_symmetric(Ac_global)
    if (cmt%my_rank == 0) write(*,'(a,i0,a,es10.3)') &
         '#### SA-AMG F4c-2b verify: global coarse op size= ', Ac_global%n/m, &
         '  relative asymmetry= ', asym

    ! S1: distributed Galerkin (row-routed) gathered == redundant Galerkin; cmt_op symmetric
    block
      real(kind=kreal) :: g_reldiff; logical :: g_sym
      call hecmw_saamg_verify_galerkin(Afull, cmt, P, naggr_local, ncnode, my_off, &
           chalo_gid, m, g_reldiff, g_sym)
      if (cmt%my_rank == 0) write(*,'(a,es10.3,a,l1)') &
           '#### SA-AMG S1 verify: |Ac_dist - Ac_global|/|Ac| = ', g_reldiff, &
           '   cmt_op symmetric= ', g_sym
    end block

    deallocate(aggr, agg_loc, chalo_gid, Bext, bcoarse, bc, y)
    call hecmw_saamg_comm_free(cmt_c)
    call hecmw_saamg_bcsr_free(P)
    call hecmw_saamg_bcsr_free(Ac_global)
    call hecmw_saamg_bcsr_free(Afull)
  end subroutine hecmw_saamg_verify_smoothp

end module hecmw_precond_SAAMG_verify
