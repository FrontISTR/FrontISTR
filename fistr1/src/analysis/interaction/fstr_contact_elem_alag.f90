!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Alag method implementations for contact element calculations
module m_fstr_contact_elem_alag
  use hecmw
  use elementInfo
  use mContactDef
  use mSurfElement
  use m_fstr_contact_geom
  use m_fstr_contact_elem_common
  use m_fstr_contact_smoothing
  implicit none

  public :: getContactStiffness_Alag
  public :: getContactNodalForce_Alag
  public :: getTiedStiffness_Alag
  public :: getTiedNodalForce_Alag
  public :: updateContactMultiplier_Alag
  public :: get_unique_map
  public :: getIntGap
  public :: build_group_tangent_basis
  public :: getTangentSlip
  public :: group_return_mapping
  public :: resolve_lambda_cur
  public :: getContactStiffness_Alag_SurfSurf
  public :: getContactNodalForce_Alag_SurfSurf

contains

  subroutine getContactStiffness_Alag(cstate, tSurf, ele, mu, mut, fcoeff, symm, stiff, force, smoothing_type, edisp, iter, &
      slvpos)

    type(tContactState), intent(inout) :: cstate       !< contact state (inout for projection info)
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: ele(:,:)        !< coord of surface element
    real(kind=kreal), intent(in)    :: mu, mut         !< penalty parameters
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient
    logical, intent(in)             :: symm            !< freeze the friction cone radius at the multiplier
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force direction
    integer(kind=kint), optional, intent(in) :: smoothing_type  !< kcsNONE or kcsNAGATA
    real(kind=kreal), optional, intent(in) :: edisp(:)  !< displacement increment for friction evaluation
    integer(kind=kint), intent(in) :: iter    !< NR iteration number (for tangent switching)
    real(kind=kreal), intent(in)    :: slvpos(3)       !< slave node position (coord+disp)

    integer          :: i, j, nnode
    real(kind=kreal) :: Bn(size(tSurf%nodes)*3+3), Ht(2,size(tSurf%nodes)*3+3), Gt(2,size(tSurf%nodes)*3+3)
    real(kind=kreal) :: metric(2,2)
    real(kind=kreal) :: A(2,2)       !< 2D local tangent operator
    real(kind=kreal) :: alpha_proj, that_dir(2)
    real(kind=kreal) :: K_fric(size(tSurf%nodes)*3+3,size(tSurf%nodes)*3+3)  !< friction stiffness
    real(kind=kreal) :: tmp_vec(2)
    real(kind=kreal) :: dummy_force(size(tSurf%nodes)*3+3)  !< dummy for computeFrictionForce_ALag
    real(kind=kreal) :: eval_disp(size(tSurf%nodes)*3+3)   !< displacement for trial friction evaluation
    real(kind=kreal) :: curpos(size(tSurf%nodes)*3+3)      !< current positions (coord+disp+ddisp)
    real(kind=kreal) :: lam_cone     !< normal force the friction cone radius is built on
    real(kind=kreal) :: Htt(size(tSurf%nodes)*3+3)  !< Ht^T * tdir: direction of the slip force
    real(kind=kreal) :: tdir(2)      !< slip direction the coupling block is linearised about
    real(kind=kreal) :: invmetric(2,2), det, norm_lamt

    nnode = size(tSurf%nodes)

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(cstate, tSurf, ele, Bn, metric, Ht, Gt, smoothing_type)

    ! Normal stiffness: stiff = mu * Bn * Bn^T
    do j = 1, nnode*3+3
      do i = 1, nnode*3+3
        stiff(i,j) = mu * Bn(i) * Bn(j)
      enddo
    enddo
    force(1:nnode*3+3) = Bn(:)

    ! frictional component
    if( fcoeff /= 0.d0 ) then
      ! Evaluate trial friction at current displacement for consistent tangent
      if( present(edisp) ) then
        eval_disp(1:nnode*3+3) = edisp(1:nnode*3+3)
      else
        eval_disp = 0.0d0
      endif
      ! Radius of the friction cone.  By default it follows the normal force this element
      ! actually applies, lambda_n + mu*g_n clipped at 0, which is the value
      ! getContactNodalForce_Alag distributes; the multiplier alone lags that force by the
      ! penalty term within a substep.  With !CONTACT_ALGO, FRICTION_CONE=FROZEN (symm) the
      ! radius stays at the multiplier of the last augmentation, which keeps the friction
      ! terms symmetric and leaves the Coulomb condition to the augmentation loop.
      if( symm ) then
        lam_cone = cstate%multiplier(1)
      else
        curpos(1:3) = slvpos(1:3) + eval_disp(1:3)
        do j = 1, nnode
          curpos(j*3+1:j*3+3) = ele(1:3,j) + eval_disp(j*3+1:j*3+3)
        enddo
        lam_cone = max( 0.d0, cstate%multiplier(1) + mu*dot_product( Bn(1:nnode*3+3), curpos(1:nnode*3+3) ) )
      endif

      call computeFrictionForce_ALag(cstate, fcoeff, lam_cone, metric, &
                                      Ht, Gt, eval_disp, nnode*3+3, dummy_force, &
                                      mut, alpha=alpha_proj, that=that_dir)

      ! Friction tangent operator A in 2D metric space:  K_fric = Ht^T * A * Ht
      if( lam_cone <= 0.0d0 .or. alpha_proj <= 1.0d-20 ) then
        ! No normal contact force: no friction contribution
        A = 0.0d0
      else if( alpha_proj >= 0.999d0 ) then
        ! Stick: A = mu_t * M (exact)
        A(1,1) = mut * metric(1,1)
        A(1,2) = mut * metric(1,2)
        A(2,1) = mut * metric(2,1)
        A(2,2) = mut * metric(2,2)
      else
        ! Slip: switch tangent by NR iteration count for stability.
        !   iter <= 2: A = alpha*mu_t*M  (stable, no directional correction)
        !   iter >= 3: A = alpha*mu_t*(M - t_hat x t_hat)  (consistent tangent)
        ! The first two NR steps use M to let t_hat stabilize; after that the consistent tangent is used to regain quadratic convergence.
        if( iter <= 2 ) then
          A(1,1) = alpha_proj * mut * metric(1,1)
          A(1,2) = alpha_proj * mut * metric(1,2)
          A(2,1) = alpha_proj * mut * metric(2,1)
          A(2,2) = alpha_proj * mut * metric(2,2)
        else
          A(1,1) = alpha_proj * mut * (metric(1,1) - that_dir(1)*that_dir(1))
          A(1,2) = alpha_proj * mut * (metric(1,2) - that_dir(1)*that_dir(2))
          A(2,1) = alpha_proj * mut * (metric(2,1) - that_dir(2)*that_dir(1))
          A(2,2) = alpha_proj * mut * (metric(2,2) - that_dir(2)*that_dir(2))
        endif
      endif

      ! Compute friction stiffness: K_fric = Ht^T * A * Ht (consistent tangent)
      do j = 1, nnode*3+3
        tmp_vec = matmul(A, Ht(1:2,j))
        do i = 1, nnode*3+3
          K_fric(i,j) = dot_product(Ht(1:2,i), tmp_vec)
        enddo
      enddo

      ! Coupling block from the radius following the normal force.  On the slip branch the
      ! friction force is f_t = R*tdir with R = fcoeff*lam_cone, and R varies with u through
      ! g_n, so d(Ht^T f_t)/du gains  Ht^T tdir * dR/du = fcoeff*mu * (Ht^T tdir) (x) Bn.
      ! A stuck node does not use the radius (f_t is the full trial), so the block belongs to
      ! the slip branch only, and it vanishes where the clip at 0 is active (lam_cone = 0).
      ! Rows are a slip direction and columns a normal map, so the block is unsymmetric: the
      ! caller has to set the linear solver up for a general matrix (see fstr_Newton_contactALag).
      if( .not.symm .and. lam_cone > 0.0d0 .and. alpha_proj > 1.0d-20 .and. alpha_proj < 0.999d0 ) then
        ! The direction is the tangential multiplier, not the trial direction that_dir, so
        ! that it stays fixed inside the augmentation step; at the fixed point of the
        ! augmentation the two coincide and the tangent is still the consistent one.
        tdir(1:2) = that_dir(1:2)
        det = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
        if( abs(det) > 1.0d-20 ) then
          invmetric(1,1) =  metric(2,2)/det
          invmetric(2,2) =  metric(1,1)/det
          invmetric(1,2) = -metric(1,2)/det
          invmetric(2,1) = -metric(2,1)/det
          tmp_vec(1:2) = matmul( invmetric(1:2,1:2), cstate%multiplier(2:3) )
          norm_lamt = dsqrt( dot_product( cstate%multiplier(2:3), tmp_vec(1:2) ) )
          if( norm_lamt > 1.0d-20 ) tdir(1:2) = cstate%multiplier(2:3) / norm_lamt
        endif
        Htt(1:nnode*3+3) = matmul( transpose(Ht(1:2,1:nnode*3+3)), tdir(1:2) )
        do j = 1, nnode*3+3
          do i = 1, nnode*3+3
            K_fric(i,j) = K_fric(i,j) + fcoeff * mu * Htt(i) * Bn(j)
          enddo
        enddo
      endif

      stiff(1:nnode*3+3,1:nnode*3+3) = stiff(1:nnode*3+3,1:nnode*3+3) + K_fric(1:nnode*3+3,1:nnode*3+3)
    endif

  end subroutine getContactStiffness_Alag

  subroutine get_unique_map(sSurf, maplist, master_idxs, unique_count)
    type(tContactSurf)  :: sSurf !< surface element structure
    integer(kind=kint), intent(out)  :: maplist(:), master_idxs(:)  !< IP->group / group->masterID, sized by the caller
    integer(kind=kint), intent(out)  :: unique_count
    integer(kind=kint)  :: tmp(MAX_N_INTP)
    integer(kind=kint)  :: i, j, n_intp, ctsurf
    logical :: found

    n_intp = sSurf%n_intp
    maplist = 0

    unique_count = 0

    do i = 1, n_intp
      if( sSurf%states(i)%state == CONTACTFREE ) cycle
      ctsurf = sSurf%states(i)%surface
      found = .false.
      ! Search existing groups by master surface index
      do j = 1, unique_count
        if (tmp(j) == ctsurf) then
          maplist(i) = j
          found = .true.
          exit
        endif
      enddo
      if (.not. found) then
        unique_count = unique_count + 1
        tmp(unique_count) = ctsurf
        maplist(i) = unique_count
      endif
    enddo

    master_idxs(1:unique_count) = tmp(1:unique_count)

  end subroutine get_unique_map

  !> \brief Compute the per-node mortar constraint quantities of one slave segment.
  !!
  !! The caller passes the unique_count/maplist/master_idxs of get_unique_map and sizes the
  !! outputs with them. For each group g and each slave-surf node a, it returns the
  !! decomposition that drives the residual/stiffness/augmentation:
  !!   Snode(g,a)      = sum_{IP in g} N_s(a) * weight                       (node tributary area)
  !!   ANnode(g,a,:)   = sum_{IP in g} N_s(a) * [N_s(b)|-N_m(k)] * weight*dir (per-node constraint accumulator)
  !!   Nsnode(g,a,:)   = ANnode(g,a,:) / Snode(g,a)                          (per-node averaged constraint grad)
  !!   gapwnode(g,a)   = ANnode(g,a,:) . curr_pos                            (per-node weighted gap)
  !! Summing over a re-collapses the group totals on a flat, uniform contact; on curved geometry
  !! they differ (one constraint per slave node).
  subroutine getIntGap(slave_surf, master, coord, disp, ddisp, &
     unique_count, maplist, master_idxs, &
     Snode, Nsnode, gapwnode)
    type(tContactSurf)  :: slave_surf
    type(tSurfElement) :: master(:)
    real(kind=kreal), intent(in)                :: coord(:), disp(:), ddisp(:)
    integer(kind=kint), intent(in)              :: unique_count
    integer(kind=kint), intent(in)              :: maplist(:), master_idxs(:)
    real(kind=kreal), intent(out)               :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:)  !< (g,a) / (g,a,24) per node in group

    integer(kind=kint) :: i, j, g, a, nnode_s, nnode_m, etype, slave, n_intp, ctsurf, nd
    integer(kind=kint) :: ndLocal(l_max_surface_node+1)
    real(kind=kreal)   :: snode_pos(3,4), weight(MAX_N_INTP)
    real(kind=kreal)   :: ncoord(2), shapefunc_s(4), shapefunc_m(4), direction(3)
    real(kind=kreal)   :: curr_pos(24)

    nnode_s = size(slave_surf%nodes)

    Snode = 0.d0
    Nsnode = 0.d0
    gapwnode = 0.d0

    ! Slave node positions at start of substep (coord + disp), used for IP weights.
    snode_pos = 0.d0
    do i = 1, nnode_s
      slave = slave_surf%nodes(i)
      snode_pos(:,i) = coord(3*slave-2:3*slave) + disp(3*slave-2:3*slave)
    enddo
    n_intp = slave_surf%n_intp
    weight = 0.d0
    call get_intp_weights(slave_surf%etype, nnode_s, n_intp, snode_pos, weight(1:n_intp))

    ! Accumulate the per-node-within-group constraint accumulator (Nsnode source) and area.
    do i = 1, n_intp
      if( slave_surf%states(i)%state == CONTACTFREE ) cycle
      ctsurf = slave_surf%states(i)%surface
      etype = master(ctsurf)%etype
      nnode_m = size(master(ctsurf)%nodes)
      direction = slave_surf%states(i)%direction(1:3)
      call getIntPoint4ss(slave_surf%etype, i, ncoord, n_intp, shapefunc_s)
      call getShapeFunc(etype, slave_surf%states(i)%lpos(1:2), shapefunc_m)
      g = maplist(i)
      ! Per-node-within-group: weight the whole IP constraint by the slave node shape function N_s(a).
      do a = 1, nnode_s
        Snode(g,a) = Snode(g,a) + shapefunc_s(a)*weight(i)
        do j = 1, nnode_s
          Nsnode(g,a,3*j-2:3*j) = Nsnode(g,a,3*j-2:3*j) &
            + shapefunc_s(a)*shapefunc_s(j)*weight(i)*direction(1:3)
        enddo
        do j = nnode_s+1, nnode_s+nnode_m
          Nsnode(g,a,3*j-2:3*j) = Nsnode(g,a,3*j-2:3*j) &
            - shapefunc_s(a)*shapefunc_m(j-nnode_s)*weight(i)*direction(1:3)
        enddo
      enddo
    enddo

    ! Per-node Nsnode (=ANnode/Snode) and weighted gap, at end of substep.
    do g = 1, unique_count
      ctsurf = master_idxs(g)
      nnode_m = size(master(ctsurf)%nodes)
      ndLocal(1:nnode_s) = slave_surf%nodes(1:nnode_s)
      ndLocal(nnode_s+1:nnode_s+nnode_m) = master(ctsurf)%nodes(1:nnode_m)
      do j = 1, nnode_s + nnode_m
        nd = ndLocal(j)
        curr_pos(3*j-2:3*j) = coord(3*nd-2:3*nd) + disp(3*nd-2:3*nd) + ddisp(3*nd-2:3*nd)
      enddo
      do a = 1, nnode_s
        ! gapwnode uses the un-normalized accumulator (= ANnode . curr_pos); the residual/aug add mu*gapwnode.
        gapwnode(g,a) = dot_product(Nsnode(g,a,1:(nnode_s+nnode_m)*3), curr_pos(1:(nnode_s+nnode_m)*3))
        if( Snode(g,a) > 0.d0 ) then
          Nsnode(g,a,1:(nnode_s+nnode_m)*3) = Nsnode(g,a,1:(nnode_s+nnode_m)*3) / Snode(g,a)
        else
          Nsnode(g,a,1:(nnode_s+nnode_m)*3) = 0.d0
        endif
      enddo
    enddo

  end subroutine getIntGap

  !> \brief Build an orthonormal tangent basis (t1,t2) as the orthogonal complement
  !!        of a group-representative slave inward normal n_hat (Householder / I - n(x)n).
  !!        The master geometry is NOT used.
  !!
  !! Picks the smallest-magnitude Cartesian axis e_k, removes its normal part to form
  !! a seed, normalizes to t1, then t2 = n_hat x t1. The result is orthonormal so the
  !! tangent metric is the identity (no metric inverse machinery needed).
  subroutine build_group_tangent_basis(n_hat, t1, t2)
    real(kind=kreal), intent(in)  :: n_hat(3)   !< unit slave inward normal (already normalized)
    real(kind=kreal), intent(out) :: t1(3), t2(3)

    integer(kind=kint) :: k
    real(kind=kreal)   :: v(3), vn, dotk

    ! choose the axis least aligned with n_hat for numerical stability
    k = 1
    if( abs(n_hat(2)) < abs(n_hat(k)) ) k = 2
    if( abs(n_hat(3)) < abs(n_hat(k)) ) k = 3

    ! v = e_k - (e_k . n_hat) n_hat  (projection of e_k onto the tangent plane)
    dotk = n_hat(k)
    v(1:3) = -dotk * n_hat(1:3)
    v(k)   = v(k) + 1.0d0

    vn = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
    t1(1:3) = v(1:3) / vn

    ! t2 = n_hat x t1  (completes the right-handed orthonormal frame)
    t2(1) = n_hat(2)*t1(3) - n_hat(3)*t1(2)
    t2(2) = n_hat(3)*t1(1) - n_hat(1)*t1(3)
    t2(3) = n_hat(1)*t1(2) - n_hat(2)*t1(1)
  end subroutine build_group_tangent_basis

  !> \brief Aggregate, per master group g and slave-surf node a of one mortar slave segment,
  !!        the mortar-weighted 3D relative displacement and the slave-normal accumulator.
  !!
  !! Mirrors getIntGap's per-IP loop (same weights, shape functions and direction field) but
  !! accumulates the mortar relative displacement du_rel(IP) = N_s . (slave ddisp) - N_m . (master ddisp).
  !! The caller passes the maplist/master_idxs/unique_count from getIntGap so the group->IP
  !! mapping matches. Each IP slip/normal is distributed to slave node a via shapefunc_s(a),
  !! so the caller can build a per-node frame, Dxi and cone radius. Outputs per group g and node a:
  !!   Sigma_node(g,a,:) = sum_IP N_s(a) * weight * du_rel(IP)   (per-node mortar slip)
  !!   nacc_node(g,a,:)  = sum_IP N_s(a) * weight * direction    (per-node slave normal)
  subroutine getTangentSlip(slave_surf, master, coord, disp, ddisp, &
     unique_count, maplist, master_idxs, Sigma_node, nacc_node)
    type(tContactSurf)  :: slave_surf
    type(tSurfElement) :: master(:)
    real(kind=kreal), intent(in)                :: coord(:), disp(:), ddisp(:)
    integer(kind=kint), intent(in)              :: unique_count
    integer(kind=kint), intent(in)              :: maplist(:), master_idxs(:)
    real(kind=kreal), intent(out)               :: Sigma_node(:,:,:)  !< (unique_count,nnode_s,3) per-node weighted 3D relative disp
    real(kind=kreal), intent(out)               :: nacc_node(:,:,:)   !< (unique_count,nnode_s,3) per-node weighted slave normal

    integer(kind=kint) :: i, j, g, a, nnode_s, nnode_m, etype, slave, n_intp, ctsurf, nd
    real(kind=kreal)   :: snode_pos(3,4), weight(MAX_N_INTP)
    real(kind=kreal)   :: ncoord(2), shapefunc_s(4), shapefunc_m(4), direction(3)
    real(kind=kreal)   :: du_rel(3), du_master(3)

    nnode_s = size(slave_surf%nodes)
    Sigma_node = 0.d0
    nacc_node  = 0.d0

    ! Slave node positions at start of substep (coord + disp), used for IP weights
    ! (identical to getIntGap so the weights match the normal aggregation exactly).
    snode_pos = 0.d0
    do i = 1, nnode_s
      slave = slave_surf%nodes(i)
      snode_pos(:,i) = coord(3*slave-2:3*slave) + disp(3*slave-2:3*slave)
    enddo
    n_intp = slave_surf%n_intp
    weight = 0.d0
    call get_intp_weights(slave_surf%etype, nnode_s, n_intp, snode_pos, weight(1:n_intp))

    do i = 1, n_intp
      if( slave_surf%states(i)%state == CONTACTFREE ) cycle
      ctsurf = slave_surf%states(i)%surface
      etype = master(ctsurf)%etype
      nnode_m = size(master(ctsurf)%nodes)
      direction = slave_surf%states(i)%direction(1:3)
      call getIntPoint4ss(slave_surf%etype, i, ncoord, n_intp, shapefunc_s)
      call getShapeFunc(etype, slave_surf%states(i)%lpos(1:2), shapefunc_m)
      g = maplist(i)

      ! Mortar relative displacement = Tm_mortar . edisp (edisp = ddisp).
      ! slave block: + N_s(j) * ddisp(slave_j) ; master block: - N_m(j) * ddisp(master_j).
      du_rel(1:3) = 0.d0
      do j = 1, nnode_s
        nd = slave_surf%nodes(j)
        du_rel(1:3) = du_rel(1:3) + shapefunc_s(j) * ddisp(3*nd-2:3*nd)
      enddo
      du_master(1:3) = 0.d0
      do j = 1, nnode_m
        nd = master(ctsurf)%nodes(j)
        du_master(1:3) = du_master(1:3) + shapefunc_m(j) * ddisp(3*nd-2:3*nd)
      enddo
      du_rel(1:3) = du_rel(1:3) - du_master(1:3)

      ! Per-node: weight the IP slip/normal by the slave node shape function N_s(a).
      do a = 1, nnode_s
        Sigma_node(g, a, 1:3) = Sigma_node(g, a, 1:3) + shapefunc_s(a) * weight(i) * du_rel(1:3)
        nacc_node(g, a, 1:3)  = nacc_node(g, a, 1:3)  + shapefunc_s(a) * weight(i) * direction(1:3)
      enddo
    enddo
  end subroutine getTangentSlip

  !> \brief Group-level Coulomb return mapping for the tangent multiplier (metric = I).
  !!
  !! Same 2D return mapping as computeFrictionForce_ALag, on the orthonormal basis (metric = I):
  !!   trial = lam_t_in + rho_t*Dxi, radius = fcoeff*lam_n
  !!   ||trial|| <= radius -> STICK (lam_t_out = trial), else SLIP (projected onto the cone).
  !!   lam_n <= 0 gives lam_t_out = 0 with the state kept.
  !! that(2) (= trial/||trial||) is returned for the slip-tangent stiffness.
  subroutine group_return_mapping(lam_t_in, rho_t, Dxi, fcoeff, lam_n, eps_fric_band, &
                                  lam_t_out, fric_state, alpha, that, update_state)
    real(kind=kreal),   intent(in)    :: lam_t_in(2)  !< incoming (warm-start) tangent multiplier, group frame
    real(kind=kreal),   intent(in)    :: rho_t        !< tangential penalty rho_t = tPenalty*refStiff
    real(kind=kreal),   intent(in)    :: Dxi(2)       !< group-frame slip increment (area-integrated)
    real(kind=kreal),   intent(in)    :: fcoeff       !< friction coefficient
    real(kind=kreal),   intent(in)    :: lam_n        !< group normal multiplier (area-integrated)
    real(kind=kreal),   intent(in)    :: eps_fric_band !< hysteresis half-band (0 = no band = legacy behavior)
    real(kind=kreal),   intent(out)   :: lam_t_out(2) !< updated tangent multiplier, group frame
    integer(kind=kint), intent(inout) :: fric_state   !< CONTACTSTICK/CONTACTSLIP (updated only when update_state)
    real(kind=kreal),   intent(out)   :: alpha        !< projection ratio (for the tangent stiffness)
    real(kind=kreal),   intent(out)   :: that(2)      !< trial/||trial|| (for the tangent stiffness)
    logical, optional,  intent(in)    :: update_state !< .true. at the augmentation commits the state; .false. (default) reads it

    real(kind=kreal) :: trial(2), norm_trial, radius
    logical          :: is_stick, do_update

    do_update = .false.
    if( present(update_state) ) do_update = update_state

    trial(1:2) = lam_t_in(1:2) + rho_t * Dxi(1:2)
    norm_trial = sqrt( trial(1)*trial(1) + trial(2)*trial(2) )   ! metric = I
    radius     = fcoeff * lam_n

    if( norm_trial > 1.0d-20 ) then
      that(1:2) = trial(1:2) / norm_trial
    else
      that(1:2) = 0.0d0
    endif

    if( lam_n <= 0.0d0 ) then
      ! No normal force: cone collapses, no friction. Leave fric_state as-is.
      lam_t_out(1:2) = 0.0d0
      alpha = 0.0d0
    else if( do_update ) then
      ! Augmentation: this is the only place where the stick/slip state is (re)decided.
      ! Hysteresis band: read the warm-start (previous) state and switch only when
      ! ||trial|| crosses the asymmetric thresholds. Inside the band the previous
      ! branch is kept, so marginal flips at the cone boundary cannot drive a limit
      ! cycle across augmentations.
      is_stick = ( fric_state == CONTACTSTICK )
      if( is_stick ) then
        if( norm_trial > (1.0d0 + eps_fric_band) * radius ) is_stick = .false.   ! genuine slip onset
      else
        if( norm_trial <= (1.0d0 - eps_fric_band) * radius ) is_stick = .true.    ! genuine stick recovery
      endif

      if( is_stick ) then
        ! Stick branch: keep full trial multiplier.
        lam_t_out(1:2) = trial(1:2)
        fric_state = CONTACTSTICK
        alpha = 1.0d0
      else
        ! Slip branch: project onto cone surface.
        alpha = radius / norm_trial
        lam_t_out(1:2) = alpha * trial(1:2)
        fric_state = CONTACTSLIP
      endif
    else
      ! Inner Newton-Raphson (read-only): pin the state frozen at the previous
      ! augmentation, skipping the band re-judgement. fric_state is not written.
      ! Same isolation as NTS-ALAG (update_multiplier=.false.): the stick/slip
      ! branch cannot flip on the live trial inside the inner NR, so the residual
      ! stays smooth and the per-iteration branch-flip limit cycle is broken.
      if( fric_state == CONTACTSTICK ) then
        ! Frozen STICK: keep full trial multiplier (no cone projection).
        lam_t_out(1:2) = trial(1:2)
        alpha = 1.0d0
      else
        ! Frozen SLIP: project onto the live cone surface with min(1, radius/||trial||), the same
        ! expression computeFrictionForce_ALag uses on the node-to-surface side.  The projection
        ! never scales a trial up: a node whose trial has come back inside the cone keeps the full
        ! trial force, which is what the alpha >= 0.999 branch of the tangent is linearised about.
        ! A node whose state was frozen SLIP at an augmentation where it carried no normal force
        ! (the lam_n <= 0 branch leaves fric_state untouched) can come back with a zero trial, and
        ! radius/0 would turn the whole residual into NaN.  No trial force means no friction force,
        ! which is what alpha = 0 gives.
        if( norm_trial > 1.0d-20 ) then
          alpha = min( 1.0d0, radius / norm_trial )
        else
          alpha = 0.0d0
        endif
        lam_t_out(1:2) = alpha * trial(1:2)
      endif
    endif
  end subroutine group_return_mapping

  !> Mortar: resolve the current per-node lambda of each active group of one slave surf.
  !> master_idxs is sorted ascending and merged against the ascending begin/working
  !> buffers with the rule: working hit -> working / else begin hit -> begin / else 0.
  subroutine resolve_lambda_cur( surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                                 lambda_node, lam_t_cur, fric_state_cur )
    type(tContactSurf), intent(in)  :: surf
    integer(kind=kint), intent(in)  :: master_idxs(:)   !< group->masterID (get_unique_map output, unsorted)
    integer(kind=kint), intent(in)  :: unique_count
    integer(kind=kint), intent(in)  :: nnode_s          !< number of slave-surf nodes
    integer(kind=kint), intent(out) :: sorted_idx(:)    !< ascending rank r -> original group g
    real(kind=kreal),   intent(out) :: lambda_node(:,:) !< (nnode_s, unique_count) per-node current lambda_n
    ! Optional friction warm-start: same reference rule as lambda_n (working -> begin ->
    ! default), riding the same merge.
    real(kind=kreal),   intent(out), optional :: lam_t_cur(:,:,:)    !< (2, nnode_s, unique_count) per-node tangent multiplier
    integer(kind=kint), intent(out), optional :: fric_state_cur(:,:) !< (nnode_s, unique_count) per-node friction state
    integer(kind=kint) :: r, j, tmp, ib, iw, g, mid
    logical            :: do_fric

    ! argsort master_idxs ascending (unique_count <= 27, insertion sort)
    do r = 1, unique_count
      sorted_idx(r) = r
    enddo
    do r = 2, unique_count
      tmp = sorted_idx(r)
      j = r - 1
      do while( j >= 1 )
        if( master_idxs(sorted_idx(j)) <= master_idxs(tmp) ) exit
        sorted_idx(j+1) = sorted_idx(j)
        j = j - 1
      enddo
      sorted_idx(j+1) = tmp
    enddo

    do_fric = present(lam_t_cur) .and. present(fric_state_cur)

    ! 2-pointer merge over ascending masters / ascending begin / ascending working
    ib = 1; iw = 1
    do r = 1, unique_count
      mid = master_idxs(sorted_idx(r))
      do while( ib <= surf%lam_begin_n .and. surf%lam_begin_id(ib) < mid ); ib = ib + 1; enddo
      do while( iw <= surf%lam_work_n  .and. surf%lam_work_id(iw)  < mid ); iw = iw + 1; enddo
      g = sorted_idx(r)
      if( iw <= surf%lam_work_n .and. surf%lam_work_id(iw) == mid ) then
        lambda_node(1:nnode_s,g) = surf%lam_work_val(1:nnode_s,iw)
        if( do_fric ) then
          lam_t_cur(1:2,1:nnode_s,g)  = surf%lam_work_t(1:2,1:nnode_s,iw)
          fric_state_cur(1:nnode_s,g) = surf%lam_work_fstate(1:nnode_s,iw)
        endif
      else if( ib <= surf%lam_begin_n .and. surf%lam_begin_id(ib) == mid ) then
        lambda_node(1:nnode_s,g) = surf%lam_begin_val(1:nnode_s,ib)
        if( do_fric ) then
          lam_t_cur(1:2,1:nnode_s,g)  = surf%lam_begin_t(1:2,1:nnode_s,ib)
          fric_state_cur(1:nnode_s,g) = surf%lam_begin_fstate(1:nnode_s,ib)
        endif
      else
        lambda_node(1:nnode_s,g) = 0.d0
        if( do_fric ) then
          lam_t_cur(1:2,1:nnode_s,g)  = 0.d0
          fric_state_cur(1:nnode_s,g) = CONTACTSTICK
        endif
      endif
    enddo
  end subroutine resolve_lambda_cur

  !> \brief Mortar (SURF-SURF) ALag: contact stiffness of one slave segment.
  !!
  !! Builds the element stiffness of every (slave-surf node a, master group g) constraint of
  !! one slave segment and hands the blocks to the caller:
  !!   stiff_n(:,:,a,g) : normal,   mu * Snode(g,a) * Nsnode(g,a) (x) Nsnode(g,a)
  !!   stiff_t(:,:,a,g) : friction, the consistent tangent of the per-node return mapping
  !!                      (plus the coupling block of a cone radius following the normal force)
  !! active_n / active_t mark the blocks that take part; an inactive block is not computed and
  !! stays zero. A block is ordered like the element vector [slave-surf nodes | master nodes of
  !! group g], so the caller builds ndLocal from master_idxs(g) and assembles the block as it
  !! stands. Normal and friction are kept apart because the caller assembles them in two passes.
  !! The caller sizes the blocks with the unique_count of get_unique_map and passes the friction
  !! pair only when fcoeff /= 0.
  subroutine getContactStiffness_Alag_SurfSurf( slave_surf, master, coord, disp, ddisp, &
      mu, mut, fcoeff, symm, eps_fric_band, unique_count, maplist, master_idxs, &
      stiff_n, active_n, stiff_t, active_t )
    type(tContactSurf), intent(in)  :: slave_surf       !< slave segment
    type(tSurfElement), intent(in)  :: master(:)        !< master surface elements
    real(kind=kreal), intent(in)    :: coord(:)         !< mesh coordinate
    real(kind=kreal), intent(in)    :: disp(:)          !< disp till current step
    real(kind=kreal), intent(in)    :: ddisp(:)         !< disp till current substep
    real(kind=kreal), intent(in)    :: mu, mut          !< penalty parameters
    real(kind=kreal), intent(in)    :: fcoeff           !< friction coefficient
    logical, intent(in)             :: symm             !< symmetricalize (cone radius frozen at the multiplier)
    real(kind=kreal), intent(in)    :: eps_fric_band    !< hysteresis half-band of the return mapping
    integer(kind=kint), intent(in)  :: unique_count     !< number of master groups of this segment
    integer(kind=kint), intent(in)  :: maplist(:)       !< integration point -> group
    integer(kind=kint), intent(in)  :: master_idxs(:)   !< group -> master surface index
    real(kind=kreal),   intent(out) :: stiff_n(:,:,:,:) !< (24,24,node,group) normal stiffness
    logical,            intent(out) :: active_n(:,:)    !< (node,group) block to assemble
    real(kind=kreal),   intent(out), optional :: stiff_t(:,:,:,:) !< (24,24,node,group) friction stiffness
    logical,            intent(out), optional :: active_t(:,:)    !< (node,group) block to assemble

    integer(kind=kint) :: g, a, j, k, nnode_m, nnode_s
    integer(kind=kint), allocatable :: sorted_idx(:)
    real(kind=kreal),   allocatable :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:), lambda_node(:,:)
    real(kind=kreal) :: Ns(24)
    ! --- friction consistent tangent ---
    real(kind=kreal),   allocatable :: lam_t_cur(:,:,:), nacc_node(:,:,:), Sigma_node(:,:,:)
    integer(kind=kint), allocatable :: fric_state_cur(:,:)
    real(kind=kreal) :: nhat(3), t1(3), t2(3), nrm, Dxi(2)
    real(kind=kreal) :: alpha, that(2), lam_t_new(2), Amat(2,2), T3d(3,2), M3(3,3)
    real(kind=kreal) :: Wb(l_max_surface_node+1)
    real(kind=kreal) :: lam_cone, that3d(3)
    integer(kind=kint) :: na, nb, fstate

    nnode_s = size(slave_surf%nodes)
    allocate(Snode(unique_count,nnode_s), Nsnode(unique_count,nnode_s,24), gapwnode(unique_count,nnode_s))
    call getIntGap(slave_surf, master, coord, disp, ddisp, &
                   unique_count, maplist, master_idxs, &
                   Snode, Nsnode, gapwnode)

    allocate(sorted_idx(unique_count), lambda_node(nnode_s,unique_count))
    stiff_n = 0.d0
    active_n = .false.
    if( fcoeff /= 0.d0 ) then
      allocate(lam_t_cur(2,nnode_s,unique_count), fric_state_cur(nnode_s,unique_count))
      call resolve_lambda_cur(slave_surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                              lambda_node, lam_t_cur, fric_state_cur)
    else
      call resolve_lambda_cur(slave_surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                              lambda_node)
    endif

    ! ===== Normal stiffness =====
    do g = 1, unique_count
      nnode_m = size(master(master_idxs(g))%nodes)
      do a = 1, nnode_s
        ! ALag contact condition per node: augmented force must be positive
        if( lambda_node(a,g)+mu*gapwnode(g,a) < 0.d0 ) cycle
        active_n(a,g) = .true.
        Ns = 0.d0
        Ns(1:(nnode_s+nnode_m)*3) = Nsnode(g, a, 1:(nnode_s+nnode_m)*3)
        do j = 1, (nnode_s+nnode_m)*3
          do k = 1, (nnode_s+nnode_m)*3
            stiff_n(j,k,a,g) = mu*Snode(g,a)*Ns(j)*Ns(k)
          enddo
        enddo
      enddo
    enddo

    if( fcoeff /= 0.d0 ) then
      ! ===== Friction consistent tangent (per slave node) =====
      ! Linearization of the per-node friction residual at the same live slip state:
      !   K_a(b,c) = Snode(g,a) * Wbar(a,b) * Wbar(a,c) * M3_a,  M3_a = T3d_a * A_a * T3d_a^T
      ! with Wbar(a,b) = Nsnode(g,a,b).nhat_a, the same map as the residual back-distribution.
      stiff_t = 0.d0
      active_t = .false.
      allocate(Sigma_node(unique_count,nnode_s,3), nacc_node(unique_count,nnode_s,3))
      call getTangentSlip(slave_surf, master, coord, disp, ddisp, &
                          unique_count, maplist, master_idxs, Sigma_node, nacc_node)
      do g = 1, unique_count
        nnode_m = size(master(master_idxs(g))%nodes)
        do a = 1, nnode_s
          if( lambda_node(a,g) <= 0.d0 ) cycle   ! no per-node normal force -> no friction
          ! Radius of the friction cone.  With FRICTION_CONE=FROZEN it stays at the multiplier
          ! of the last augmentation, which keeps the friction terms symmetric and leaves the
          ! Coulomb condition to the augmentation loop; with !CONTACT_ALGO, FRICTION_CONE=FOLLOW it
          ! follows the normal force this node actually applies, lambda_node + rho_n*gapwnode,
          ! the same expression the residual distributes as nrlforce.
          if( symm ) then
            lam_cone = lambda_node(a,g)
          else
            lam_cone = lambda_node(a,g) + mu*gapwnode(g,a)
          endif
          nrm = sqrt( nacc_node(g,a,1)**2 + nacc_node(g,a,2)**2 + nacc_node(g,a,3)**2 )
          if( nrm < 1.d-30 ) cycle
          nhat(1:3) = nacc_node(g,a,1:3) / nrm
          call build_group_tangent_basis(nhat, t1, t2)
          ! Live per-node slip projection and read-only return mapping (writes only OUT args).
          Dxi(1) = dot_product(t1(1:3), Sigma_node(g,a,1:3))
          Dxi(2) = dot_product(t2(1:3), Sigma_node(g,a,1:3))
          fstate = fric_state_cur(a,g)
          call group_return_mapping(lam_t_cur(1:2,a,g), mut, Dxi, fcoeff, lam_cone, &
                                    eps_fric_band, lam_t_new, fstate, alpha, that)
          ! 2D tangent operator A (same construction as getContactStiffness_Alag).
          if( alpha <= 1.0d-20 ) then
            Amat = 0.d0
          else if( alpha >= 0.999d0 ) then
            Amat = 0.d0
            Amat(1,1) = mut
            Amat(2,2) = mut
          else
            Amat(1,1) = alpha * mut * (1.0d0 - that(1)*that(1))
            Amat(1,2) = alpha * mut * (-that(1)*that(2))
            Amat(2,1) = alpha * mut * (-that(2)*that(1))
            Amat(2,2) = alpha * mut * (1.0d0 - that(2)*that(2))
          endif
          T3d(1:3,1) = t1(1:3)
          T3d(1:3,2) = t2(1:3)
          M3 = matmul( matmul(T3d, Amat), transpose(T3d) )

          ! Per-node averaged mortar weight of each node (= ANnode/Snode, recovered via nhat_a).
          do na = 1, nnode_s + nnode_m
            Wb(na) = dot_product(Nsnode(g,a,3*na-2:3*na), nhat(1:3))
          enddo
          active_t(a,g) = .true.
          do nb = 1, nnode_s + nnode_m
            do na = 1, nnode_s + nnode_m
              do k = 1, 3
                do j = 1, 3
                  stiff_t(3*na-3+j, 3*nb-3+k, a, g) = Snode(g,a) * Wb(na) * Wb(nb) * M3(j,k)
                enddo
              enddo
            enddo
          enddo
          ! Coupling block of a cone radius that follows the normal force.  On the slip branch
          ! the friction force is f_t = R*that3d with R = fcoeff*lam_cone, and R varies with u
          ! through gapwnode:  d(gapwnode(g,a))/du = Snode(g,a)*Nsnode(g,a,:), the map the normal
          ! stiffness uses, so the residual -Wbar(a,b)*f_t gains
          !   K_a(b,c) += fcoeff*rho_n*Snode(g,a) * Wbar(a,b)*that3d (x) Nsnode(g,a,c).
          ! Rows are a slip direction and columns a normal map, so the block is unsymmetric and
          ! the solver is set up for a general matrix (fstr_is_contactALag_symmetric).  A stuck
          ! node does not use the radius (f_t is the full trial), hence the slip-branch window,
          ! the same one the consistent tangent above uses.
          if( .not.symm .and. alpha > 1.0d-20 .and. alpha < 0.999d0 ) then
            that3d(1:3) = that(1)*t1(1:3) + that(2)*t2(1:3)
            do nb = 1, nnode_s + nnode_m
              do na = 1, nnode_s + nnode_m
                do k = 1, 3
                  do j = 1, 3
                    stiff_t(3*na-3+j, 3*nb-3+k, a, g) = stiff_t(3*na-3+j, 3*nb-3+k, a, g) &
                      + fcoeff * mu * Snode(g,a) * Wb(na) * that3d(j) * Nsnode(g,a,3*nb-3+k)
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      deallocate(Sigma_node, nacc_node, lam_t_cur, fric_state_cur)
    endif

    deallocate(sorted_idx)
    deallocate(Snode, Nsnode, gapwnode, lambda_node)
  end subroutine getContactStiffness_Alag_SurfSurf

  !> \brief Mortar (SURF-SURF) ALag: contact nodal force of one slave segment.
  !!
  !! Builds the element force of every (slave-surf node a, master group g) constraint of one
  !! slave segment, in the same (a,g) layout as getContactStiffness_Alag_SurfSurf:
  !!   ctNForce(:,a,g) : normal,   -(lambda_node + mu*gapwnode) * Nsnode(g,a)   for the residual,
  !!                               -lambda_node * Nsnode(g,a)                   for the output
  !!   ctTForce(:,a,g) : friction, the traction of the per-node return mapping distributed
  !!                     through the per-node mortar weight Wbar(a,j) = Nsnode(g,a,j).nhat_a
  !!                     (kctForOutput keeps the frozen multiplier instead of the live trial)
  !! The force is signed as the residual contribution, so the caller only adds it up.
  !! active_n / active_t, the block ordering and the sizing of the outputs are as in
  !! getContactStiffness_Alag_SurfSurf.
  subroutine getContactNodalForce_Alag_SurfSurf( purpose, slave_surf, master, coord, disp, ddisp, &
      mu, mut, fcoeff, symm, eps_fric_band, unique_count, maplist, master_idxs, &
      ctNForce, active_n, ctTForce, active_t )
    integer(kind=kint), intent(in)  :: purpose          !< kctForResidual or kctForOutput
    type(tContactSurf), intent(in)  :: slave_surf       !< slave segment
    type(tSurfElement), intent(in)  :: master(:)        !< master surface elements
    real(kind=kreal), intent(in)    :: coord(:)         !< mesh coordinate
    real(kind=kreal), intent(in)    :: disp(:)          !< disp till current step
    real(kind=kreal), intent(in)    :: ddisp(:)         !< disp till current substep
    real(kind=kreal), intent(in)    :: mu, mut          !< penalty parameters
    real(kind=kreal), intent(in)    :: fcoeff           !< friction coefficient
    logical, intent(in)             :: symm             !< symmetricalize (cone radius frozen at the multiplier)
    real(kind=kreal), intent(in)    :: eps_fric_band    !< hysteresis half-band of the return mapping
    integer(kind=kint), intent(in)  :: unique_count     !< number of master groups of this segment
    integer(kind=kint), intent(in)  :: maplist(:)       !< integration point -> group
    integer(kind=kint), intent(in)  :: master_idxs(:)   !< group -> master surface index
    real(kind=kreal),   intent(out) :: ctNForce(:,:,:)  !< (24,node,group) normal force vector
    logical,            intent(out) :: active_n(:,:)    !< (node,group) vector to assemble
    real(kind=kreal),   intent(out), optional :: ctTForce(:,:,:) !< (24,node,group) friction force vector
    logical,            intent(out), optional :: active_t(:,:)   !< (node,group) vector to assemble

    integer(kind=kint) :: g, a, j, nnode_m, nnode_s
    integer(kind=kint), allocatable :: sorted_idx(:)
    real(kind=kreal),   allocatable :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:), lambda_node(:,:)
    real(kind=kreal) :: nrlforce
    real(kind=kreal) :: Ns(24)
    ! --- friction force back-distribution (live return mapping) ---
    real(kind=kreal),   allocatable :: lam_t_cur(:,:,:), nacc_node(:,:,:), Sigma_node(:,:,:)
    integer(kind=kint), allocatable :: fric_state_cur(:,:)
    real(kind=kreal) :: nhat(3), t1(3), t2(3), nrm, fvec(3), Wbar
    real(kind=kreal) :: Dxi(2), alpha, that(2), lam_t_new(2), lam_cone
    integer(kind=kint) :: fstate

    nnode_s = size(slave_surf%nodes)
    allocate(Snode(unique_count,nnode_s), Nsnode(unique_count,nnode_s,24), gapwnode(unique_count,nnode_s))
    call getIntGap(slave_surf, master, coord, disp, ddisp, &
                   unique_count, maplist, master_idxs, &
                   Snode, Nsnode, gapwnode)

    allocate(sorted_idx(unique_count), lambda_node(nnode_s,unique_count))
    ctNForce = 0.d0
    active_n = .false.
    if( fcoeff /= 0.d0 ) then
      allocate(lam_t_cur(2,nnode_s,unique_count), fric_state_cur(nnode_s,unique_count))
      call resolve_lambda_cur(slave_surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                              lambda_node, lam_t_cur, fric_state_cur)
    else
      call resolve_lambda_cur(slave_surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                              lambda_node)
    endif

    ! ===== Normal force: per-node back-distribution =====
    ! nrlforce_a = lambda_node(a,g) + mu*gapwnode(g,a) for the residual, lambda_node(a,g) for output.
    do g = 1, unique_count
      nnode_m = size(master(master_idxs(g))%nodes)
      do a = 1, nnode_s
        nrlforce = lambda_node(a,g) + mu*gapwnode(g,a)
        ! ALag contact condition per node: augmented force must be positive
        if( nrlforce < 0.d0 ) cycle
        active_n(a,g) = .true.
        Ns = 0.d0
        Ns(1:(nnode_s+nnode_m)*3) = Nsnode(g, a, 1:(nnode_s+nnode_m)*3)
        do j = 1, nnode_s + nnode_m
          if( purpose == kctForResidual ) then
            ctNForce(3*j-2:3*j,a,g) = -nrlforce*Ns(3*j-2:3*j)
          else if ( purpose == kctForOutput ) then
            ! Output: multiplier only (converges to true contact force)
            ctNForce(3*j-2:3*j,a,g) = -lambda_node(a,g)*Ns(3*j-2:3*j)
          end if
        enddo
      enddo
    enddo

    if( fcoeff /= 0.d0 ) then
      ! ===== Friction force (per slave node): live return mapping, back-distributed =====
      ! trial = lam_t_warm(a) + rho_t*Dxi_live(a), projected onto the cone of radius lambda_node(a,g).
      ! The return mapping is read-only here (the augmentation update is the sole writer of the
      ! lambda_t / fric_state buffers). The resulting traction is distributed through the per-node
      ! mortar weight Wbar(a,j) = Nsnode(g,a,j).nhat_a, mirroring the normal back-distribution.
      ! Output (kctForOutput) keeps the frozen multiplier.
      ctTForce = 0.d0
      active_t = .false.
      allocate(Sigma_node(unique_count,nnode_s,3), nacc_node(unique_count,nnode_s,3))
      call getTangentSlip(slave_surf, master, coord, disp, ddisp, &
                          unique_count, maplist, master_idxs, Sigma_node, nacc_node)
      do g = 1, unique_count
        nnode_m = size(master(master_idxs(g))%nodes)
        do a = 1, nnode_s
          if( lambda_node(a,g) <= 0.d0 ) cycle   ! no per-node normal force -> no friction
          ! Cone radius: the frozen multiplier with FRICTION_CONE=FROZEN, the normal force this
          ! node just applied above (nrlforce = lambda_node + rho_n*gapwnode) with FRICTION_CONE=FOLLOW,
          ! the same radius the tangent uses (see getContactStiffness_Alag_SurfSurf).  A negative
          ! lam_cone reaches group_return_mapping as lam_n <= 0 and gives zero friction, which is
          ! what the normal back-distribution above does with a negative nrlforce too.
          if( symm ) then
            lam_cone = lambda_node(a,g)
          else
            lam_cone = lambda_node(a,g) + mu*gapwnode(g,a)
          endif
          nrm = sqrt( nacc_node(g,a,1)**2 + nacc_node(g,a,2)**2 + nacc_node(g,a,3)**2 )
          if( nrm < 1.d-30 ) cycle
          nhat(1:3) = nacc_node(g,a,1:3) / nrm
          call build_group_tangent_basis(nhat, t1, t2)
          if( purpose == kctForResidual ) then
            ! Live trial: project the live per-node mortar slip onto the per-node frame and return-map.
            Dxi(1) = dot_product(t1(1:3), Sigma_node(g,a,1:3))
            Dxi(2) = dot_product(t2(1:3), Sigma_node(g,a,1:3))
            fstate = fric_state_cur(a,g)
            call group_return_mapping(lam_t_cur(1:2,a,g), mut, Dxi, fcoeff, lam_cone, &
                                      eps_fric_band, lam_t_new, fstate, alpha, that)
            fvec(1:3) = lam_t_new(1)*t1(1:3) + lam_t_new(2)*t2(1:3)
          else
            ! Output: frozen multiplier only (converges to the true friction force).
            fvec(1:3) = lam_t_cur(1,a,g)*t1(1:3) + lam_t_cur(2,a,g)*t2(1:3)
          end if

          active_t(a,g) = .true.
          do j = 1, nnode_s + nnode_m
            ! per-node averaged mortar weight of node j (= ANnode/Snode, recovered via nhat_a)
            Wbar = dot_product(Nsnode(g,a,3*j-2:3*j), nhat(1:3))
            ctTForce(3*j-2:3*j,a,g) = -fvec(1:3) * Wbar
          enddo
        enddo
      enddo

      deallocate(Sigma_node, nacc_node, lam_t_cur, fric_state_cur)
    endif

    deallocate(sorted_idx)
    deallocate(Snode, Nsnode, gapwnode, lambda_node)
  end subroutine getContactNodalForce_Alag_SurfSurf

  subroutine getContactNodalForce_Alag(ctState,tSurf,ndCoord,ndDu,mu,mut,fcoeff,symm,lagrange,ctNForce,ctTForce,cflag, &
      smoothing_type)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint) :: nnode !< number of nodes of master segment
    integer(kind=kint) :: j
    real(kind=kreal), intent(in) :: mu, mut !< penalty parameters
    real(kind=kreal)   :: fcoeff !< friction coefficient
    logical, intent(in) :: symm  !< freeze the friction cone radius at the multiplier
    real(kind=kreal)   :: lagrange !< not used for ALagrange (kept for interface compatibility)
    real(kind=kreal)   :: ndCoord(:), ndDu(:) !< nodal coordinates (coord+disp+ddisp); nodal displacement increment (ddisp)
    real(kind=kreal)   :: ctNForce(:) !< contact normal force vector
    real(kind=kreal)   :: ctTForce(:) !< contact tangential force vector
    logical            :: cflag  !< not used for ALagrange (kept for interface compatibility)
    integer(kind=kint), optional, intent(in) :: smoothing_type  !< kcsNONE or kcsNAGATA

    real(kind=kreal)   :: normal(3) !< normal vector at target point
    real(kind=kreal)   :: Bn(3*l_max_elem_node+3) !< normal distribution vector
    real(kind=kreal)   :: Ht(2,3*l_max_elem_node+3), Gt(2,3*l_max_elem_node+3) !< tangent and covariant maps
    real(kind=kreal)   :: elemcrd(3, l_max_elem_node) !< master node coords (coord+disp, for computeContactMaps_ALag)
    real(kind=kreal)   :: edisp(3*l_max_elem_node+3) !< displacement increment
    real(kind=kreal)   :: dgn, nrlforce !< normal gap; normal force
    real(kind=kreal)   :: lam_cone !< normal force the friction cone radius is built on
    real(kind=kreal)   :: metric(2,2)
    integer(kind=kint) :: edof  !< element vector size (nnode*3+3)

    nnode = size(tSurf%nodes)
    edof = nnode*3+3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    normal(1:3) = ctState%direction(1:3)

    ! Prepare elemcrd = ndCoord - ndDu (i.e., coord + disp) for computeContactMaps_ALag
    do j = 1, nnode
      elemcrd(1:3, j) = ndCoord(j*3+1:j*3+3) - ndDu(j*3+1:j*3+3)
    enddo

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(ctState, tSurf, elemcrd(:,1:nnode), &
                                  Bn, metric, Ht, Gt, smoothing_type)

    ! Normal gap: dgn = Bn^T * ndCoord (using normal distribution vector)
    dgn = dot_product( Bn(1:edof), ndCoord(1:edof) )

    ! Normal force: multiplier + penalty * gap
    nrlforce = ctState%multiplier(1) + mu*dgn

    ! Distribute normal force using Bn: ctNForce = -nrlforce * Bn
    ctNForce(1:edof) = -nrlforce * Bn(1:edof)

    ! Lagrange row (not used in ALagrange, set to 0)
    ctNForce((nnode+1)*3+1) = 0.d0

    if( fcoeff == 0.d0 ) return

    ! --- Tangent component ---

    ! Prepare edisp from ndDu
    edisp(1:3) = ndDu(1:3)  ! slave
    do j = 1, nnode
      edisp(j*3+1:j*3+3) = ndDu(j*3+1:j*3+3)  ! master nodes
    enddo

    ! Compute friction force using common routine.  With FRICTION_CONE=FOLLOW the cone radius
    ! is bounded by the normal force just distributed above rather than by the multiplier
    ! alone, the same radius the tangent uses (see getContactStiffness_Alag).
    if( symm ) then
      lam_cone = ctState%multiplier(1)
    else
      lam_cone = max( 0.d0, nrlforce )
    endif
    call computeFrictionForce_ALag(ctState, fcoeff, lam_cone, metric, &
                                    Ht, Gt, edisp, edof, ctTForce, &
                                    mut)

    ! Lagrange row (not used in ALagrange, set to 0)
    ctTForce((nnode+1)*3+1) = 0.d0

  end subroutine getContactNodalForce_Alag

  subroutine updateContactMultiplier_Alag(ctState,ndLocal,coord,disp,ddisp,&
     &  mu,mut,fcoeff,tSurf,lgnt,ctchanged,ctNForce,ctTForce,jump_ratio,smoothing_type)

    type(tContactState), intent(inout)   :: ctState             !< contact state
    integer(kind=kint), intent(in)       :: ndLocal(:)          !< global node numbers (slave + master)
    real(kind=kreal), intent(in)         :: coord(:)            !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)             !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)            !< disp till current substep
    real(kind=kreal), intent(in)         :: mu, mut             !< penalty parameters
    real(kind=kreal), intent(in)         :: fcoeff              !< friction coefficient
    type(tSurfElement), intent(in)       :: tSurf               !< surface element structure (with vertex_normals for Nagata)
    real(kind=kreal), intent(inout)      :: lgnt(2)             !< convergence metrics
    logical, intent(inout)               :: ctchanged           !< contact state changed flag
    real(kind=kreal), intent(out)        :: ctNForce(:)         !< contact normal force vector
    real(kind=kreal), intent(out)        :: ctTForce(:)         !< contact tangential force vector
    real(kind=kreal), intent(out)        :: jump_ratio          !< stick trial / slip limit force ratio
    integer(kind=kint), optional, intent(in) :: smoothing_type  !< kcsNONE or kcsNAGATA

    integer(kind=kint)  :: nnode !< number of master nodes
    integer(kind=kint)  :: slave, j
    real(kind=kreal)    :: Bn(3*l_max_elem_node+3) !< normal distribution vector
    real(kind=kreal)    :: Ht(2,3*l_max_elem_node+3), Gt(2,3*l_max_elem_node+3) !< tangent and covariant maps
    real(kind=kreal)    :: elemcrd(3,l_max_elem_node) !< master node coords (coord+disp, for computeContactMaps_ALag)
    real(kind=kreal)    :: curpos(3*l_max_elem_node+3) !< current positions (coord+disp+ddisp)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3) !< displacement increment
    real(kind=kreal)    :: dgn, nrlforce !< normal gap; normal force
    real(kind=kreal)    :: metric(2,2)
    real(kind=kreal)    :: dxy(2)        !< for convergence check
    integer(kind=kint)  :: edof  !< element vector size (nnode*3+3)

    nnode = size(ndLocal) - 1
    slave = ndLocal(1)
    edof = nnode*3+3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Prepare elemcrd (coord+disp) and current positions (coord+disp+ddisp)
    curpos(1:3) = coord(3*slave-2:3*slave) + disp(3*slave-2:3*slave) + ddisp(3*slave-2:3*slave)
    edisp(1:3) = ddisp(3*slave-2:3*slave)
    do j = 1, nnode
      elemcrd(1:3,j) = coord(3*ndLocal(j+1)-2:3*ndLocal(j+1)) + disp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
      curpos(j*3+1:j*3+3) = elemcrd(1:3,j) + ddisp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
      edisp(j*3+1:j*3+3) = ddisp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
    enddo

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(ctState, tSurf, elemcrd(:,1:nnode), &
                                  Bn, metric, Ht, Gt, smoothing_type)

    ! Normal gap: dgn = Bn^T * curpos (using normal distribution vector)
    dgn = dot_product( Bn(1:edof), curpos(1:edof) )

    ! Update multiplier and working distance
    ctState%wkdist = -dgn
    ctState%multiplier(1) = ctState%multiplier(1) - mu*ctState%wkdist
    ctState%distance = ctState%wkdist
    lgnt(1) = lgnt(1) - ctState%wkdist

    ! Normal force: use updated multiplier
    nrlforce = ctState%multiplier(1)

    ! Distribute normal force using Bn: ctNForce = -nrlforce * Bn
    ctNForce(1:edof) = -nrlforce * Bn(1:edof)

    if( fcoeff == 0.d0 ) return

    ! --- Tangent component ---

    ! The multiplier has just absorbed mu*g_n above, so it already is the normal force this
    ! configuration applies and serves as the cone radius for either FRICTION_CONE setting.
    call computeFrictionForce_ALag(ctState, fcoeff, ctState%multiplier(1), metric, &
                                    Ht, Gt, edisp, edof, ctTForce, &
                                    mut, &
                                    update_multiplier=.true., slave_id=slave, ctchanged=ctchanged, &
                                    jump_ratio=jump_ratio)

    ! Tangent displacement for convergence check: use Gt to project curpos directly
    dxy = matmul( Gt(:,1:edof), curpos(1:edof) )
    lgnt(2) = lgnt(2) + dsqrt( dxy(1)*dxy(1) + dxy(2)*dxy(2) )

  end subroutine updateContactMultiplier_Alag

  subroutine getTiedStiffness_Alag(cstate, tSurf, mu, stiff, force)

    type(tContactState), intent(in) :: cstate          !< contact state
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: mu              !< penalty parameter
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force direction

    integer          :: i, j, nnode, edof
    real(kind=kreal) :: Tm(3, 3*(l_max_surface_node+1))
    real(kind=kreal) :: Tt(3, 3*(l_max_surface_node+1))  !< unused, required by computeTm_Tt interface

    nnode = size(tSurf%nodes)
    edof = nnode*3+3

    stiff = 0.d0

    ! Use common mapping routine to compute Tm
    call computeTm_Tt(cstate, tSurf, 0.0d0, Tm, Tt)

    ! Tied stiffness: stiff = mu * Tm^T * Tm
    do j = 1, edof
      do i = 1, edof
        stiff(i,j) = mu * dot_product(Tm(1:3,i), Tm(1:3,j))
      enddo
    enddo
    force(1:edof) = 0.d0  ! not used for tied (3-direction constraint)

  end subroutine getTiedStiffness_Alag

  subroutine getTiedNodalForce_Alag(ctState,tSurf,ndu,mu,ctNForce,ctTForce)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    real(kind=kreal)   :: ndu(:) !< nodal total displacement (disp+ddisp)
    real(kind=kreal), intent(in) :: mu !< penalty parameter
    real(kind=kreal)   :: ctNForce(:)  !< contact force vector
    real(kind=kreal)   :: ctTForce(:)  !< contact force vector (not used for tied)

    integer(kind=kint) :: edof
    real(kind=kreal)   :: Tm(3, 3*(l_max_surface_node+1))  !< relative displacement mapping matrix
    real(kind=kreal)   :: Tt(3, 3*(l_max_surface_node+1))  !< unused, required by computeTm_Tt interface
    real(kind=kreal)   :: dg(3) !< gap vector
    real(kind=kreal)   :: nrlforce(3) !< nodal force (3 components)

    nnode = size(tSurf%nodes)
    edof = nnode*3+3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Use common mapping routine to compute Tm
    call computeTm_Tt(ctState, tSurf, 0.0d0, Tm, Tt)

    ! Gap vector: dg = Tm * ndu (3-component relative displacement)
    dg(1:3) = matmul(Tm(1:3, 1:edof), ndu(1:edof))

    ! Force: multiplier + penalty * gap (3 components)
    nrlforce(1:3) = ctState%multiplier(1:3) + mu*dg(1:3)

    ! Distribute force: ctNForce = -Tm^T * nrlforce
    ctNForce(1:edof) = -matmul(transpose(Tm(1:3, 1:edof)), nrlforce(1:3))

    ! Lagrange row (not used in ALagrange, set to 0)
    ctNForce((nnode+1)*3+1) = 0.d0
    ctTForce((nnode+1)*3+1) = 0.d0

  end subroutine getTiedNodalForce_Alag

end module m_fstr_contact_elem_alag
