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
    integer(kind=kint), allocatable, intent(out)  :: maplist(:), master_idxs(:)
    integer(kind=kint), intent(out)  :: unique_count
    integer(kind=kint), allocatable :: tmp(:)
    integer(kind=kint)  :: i, j, n_intp, ctsurf
    logical :: found

    n_intp = sSurf%n_intp
    allocate(maplist(n_intp))
    allocate(tmp(n_intp))
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

    allocate(master_idxs(unique_count))
    master_idxs = tmp(1:unique_count)

  end subroutine get_unique_map

  !> \brief Compute group-level normal distribution (Ns), area (S), and current gap (g)
  !!        for one mortar slave segment.
  !!
  !! For each unique master surface in contact with this slave segment, accumulates
  !!   AN_g = sum_ip (shapefunc * weight * direction)
  !!   S_g  = sum_ip (weight)
  !! and returns Ns_g = AN_g / S_g, gap_g = Ns_g . curr_pos.
  !!
  !! Caller owns the returned allocatables. It also returns the per-node-within-group
  !! decomposition that drives the residual/stiffness/augmentation: for group g and node a,
  !!   Snode(g,a)      = sum_{IP in g} N_s(a) * weight                       (node tributary area)
  !!   ANnode(g,a,:)   = sum_{IP in g} N_s(a) * [N_s(b)|-N_m(k)] * weight*dir (per-node constraint accumulator)
  !!   Nsnode(g,a,:)   = ANnode(g,a,:) / Snode(g,a)                          (per-node averaged constraint grad)
  !!   gapwnode(g,a)   = ANnode(g,a,:) . curr_pos                            (per-node weighted gap)
  !! Summing over a re-collapses the group totals on a flat, uniform contact; on curved geometry
  !! they differ (one constraint per slave node). The group S / Ns_list / integrated_gaps are also returned.
  subroutine getIntGap(slave_surf, master, coord, disp, ddisp, &
     unique_count, maplist, master_idxs, S, Ns_list, integrated_gaps, &
     Snode, Nsnode, gapwnode)
    type(tContactSurf)  :: slave_surf
    type(tSurfElement) :: master(:)
    real(kind=kreal), intent(in)                :: coord(:), disp(:), ddisp(:)
    integer(kind=kint), intent(out)             :: unique_count
    integer(kind=kint), allocatable, intent(out):: maplist(:), master_idxs(:)
    real(kind=kreal), allocatable, intent(out)  :: S(:), Ns_list(:, :), integrated_gaps(:)
    real(kind=kreal), allocatable, intent(out)  :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:)  !< (g,a) / (g,a,24) per-node-within-group

    integer(kind=kint) :: i, j, g, a, nnode_s, nnode_m, etype, slave, n_intp, ctsurf, nd
    integer(kind=kint) :: ndLocal(l_max_surface_node+1)
    real(kind=kreal)   :: snode_pos(3,4), weight(MAX_N_INTP)
    real(kind=kreal)   :: ncoord(2), shapefunc_s(4), shapefunc_m(4), direction(3)
    real(kind=kreal)   :: curr_pos(24)
    real(kind=kreal), allocatable :: AN(:,:)

    call get_unique_map(slave_surf, maplist, master_idxs, unique_count)

    nnode_s = size(slave_surf%nodes)

    allocate(AN(unique_count, 24))
    allocate(S(unique_count))
    allocate(Ns_list(unique_count, 24))
    allocate(integrated_gaps(unique_count))
    allocate(Snode(unique_count, nnode_s))
    allocate(Nsnode(unique_count, nnode_s, 24))
    allocate(gapwnode(unique_count, nnode_s))
    AN = 0.d0
    S = 0.d0
    Ns_list = 0.d0
    integrated_gaps = 0.d0
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

    ! Accumulate AN and S per group, and the per-node-within-group constraint accumulator (Nsnode source).
    do i = 1, n_intp
      if( slave_surf%states(i)%state == CONTACTFREE ) cycle
      ctsurf = slave_surf%states(i)%surface
      etype = master(ctsurf)%etype
      nnode_m = size(master(ctsurf)%nodes)
      direction = slave_surf%states(i)%direction(1:3)
      call getIntPoint4ss(slave_surf%etype, i, ncoord, n_intp, shapefunc_s)
      call getShapeFunc(etype, slave_surf%states(i)%lpos(1:2), shapefunc_m)
      g = maplist(i)
      do j = 1, nnode_s
        AN(g, 3*j-2:3*j) = AN(g, 3*j-2:3*j) + shapefunc_s(j)*weight(i)*direction(1:3)
      enddo
      do j = nnode_s+1, nnode_s+nnode_m
        AN(g, 3*j-2:3*j) = AN(g, 3*j-2:3*j) - shapefunc_m(j-nnode_s)*weight(i)*direction(1:3)
      enddo
      S(g) = S(g) + weight(i)
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

    ! Per-group Ns and current gap; per-node Nsnode (=ANnode/Snode) and weighted gap, all at end of substep.
    do g = 1, unique_count
      ctsurf = master_idxs(g)
      nnode_m = size(master(ctsurf)%nodes)
      ndLocal(1:nnode_s) = slave_surf%nodes(1:nnode_s)
      ndLocal(nnode_s+1:nnode_s+nnode_m) = master(ctsurf)%nodes(1:nnode_m)
      Ns_list(g, 1:(nnode_s+nnode_m)*3) = AN(g, 1:(nnode_s+nnode_m)*3) / S(g)
      do j = 1, nnode_s + nnode_m
        nd = ndLocal(j)
        curr_pos(3*j-2:3*j) = coord(3*nd-2:3*nd) + disp(3*nd-2:3*nd) + ddisp(3*nd-2:3*nd)
      enddo
      integrated_gaps(g) = dot_product(Ns_list(g, 1:(nnode_s+nnode_m)*3),curr_pos(1:(nnode_s+nnode_m)*3))
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

    deallocate(AN)
  end subroutine getIntGap

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
