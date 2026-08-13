!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact processing at assembly level (all pairs in one tContact object)
!>
module m_fstr_contact_assembly
  use hecmw
  use m_fstr
  use mContactDef
  use m_fstr_contact_element
  use m_fstr_contact_elem_common, only: computeTm_Tt
  use m_fstr_contact_interference
  use m_fstr_contact_elem_alag
  use m_fstr_contact_damping
  implicit none

  public :: calc_contact_pair_refStiff

contains

  !> Calculate reference stiffness for one contact pair
  subroutine calc_contact_pair_refStiff(contact, diag, ndof, hecMESH)
    type(tContact), intent(inout)        :: contact   !< contact pair
    real(kind=kreal), intent(in)         :: diag(:)   !< diagonal vector (size = ndof * np)
    integer(kind=kint), intent(in)       :: ndof      !< degrees of freedom
    type(hecmwST_local_mesh), intent(in) :: hecMESH   !< mesh
    
    integer(kind=kint) :: i, j, k, slave_node, master_node, nnode, ctsurf
    integer(kind=kint) :: idx_start, idx_end, n_slave
    integer(kind=kint) :: cgrp, ic, iss, outtype, fnodes(100)
    real(kind=kreal)   :: maxv
    real(kind=kreal)   :: A_rep, slave_reflen_sum
    real(kind=kreal)   :: elem(3, l_max_surface_node), r0(2)

    maxv = 0.0d0

    ! Loop over slave nodes
    do j = 1, size(contact%slave)
      slave_node = contact%slave(j)
      ! The mortar refStiff must be partition-invariant: skip GHOST(external) rows, whose
      ! diagonal is not fully assembled here. Every slave node is internal on exactly one
      ! rank, so the allreduce-MAX below still sees every fully-assembled diagonal.
      if( contact%method == CONTACTS2S .and. slave_node > hecMESH%nn_internal ) cycle
      idx_start = ndof * (slave_node - 1) + 1
      idx_end = ndof * slave_node
      maxv = max(maxv, maxval(diag(idx_start:idx_end)))
    enddo
    
    ! Loop over master surfaces and nodes
    if( contact%method == CONTACTS2S ) then
      ! Enumerate the master faces from the surface group instead of contact%master: with
      ! !PARTITION, CONTACT_OWNER=SLAVE a rank that owns no slave node takes no master surface
      ! at all (fstr_contact_init), so the diagonals of the master nodes it owns would never
      ! enter the max and refStiff would depend on the partition. The surface group items are
      ! present wherever the element is, so taking internal rows only on every rank lets the
      ! allreduce-MAX reproduce the serial value.
      cgrp = contact%surf_id2
      if( cgrp > 0 ) then
        do i = hecMESH%surf_group%grp_index(cgrp-1)+1, hecMESH%surf_group%grp_index(cgrp)
          ic = hecMESH%surf_group%grp_item(2*i-1)
          call getSubFace( hecMESH%elem_type(ic), hecMESH%surf_group%grp_item(2*i), outtype, fnodes )
          nnode = getNumberOfNodes( outtype )
          iss = hecMESH%elem_node_index(ic-1)
          do j = 1, nnode
            master_node = hecMESH%elem_node_item( iss + fnodes(j) )
            if( master_node > hecMESH%nn_internal ) cycle
            idx_start = ndof * (master_node - 1) + 1
            idx_end = ndof * master_node
            maxv = max(maxv, maxval(diag(idx_start:idx_end)))
          enddo
        enddo
      endif
    else
      do ctsurf = 1, size(contact%master)
        nnode = size(contact%master(ctsurf)%nodes)
        do j = 1, nnode
          master_node = contact%master(ctsurf)%nodes(j)
          idx_start = ndof * (master_node - 1) + 1
          idx_end = ndof * master_node
          maxv = max(maxv, maxval(diag(idx_start:idx_end)))
        enddo
      enddo
    endif
    
    ! Parallel reduction
    call hecmw_allREDUCE_R1(hecMESH, maxv, hecmw_max)
    
    ! Set reference stiffness for this contact pair
    contact%refStiff = maxv

    ! Mortar dimensional correction: the mortar penalty multiplies refStiff by the contact
    ! area, so refStiff is divided by a representative tributary area A_rep = (slave reflen)^2
    ! to make mu = nPenalty*refStiff a pressure density. The same nPenalty then gives the same
    ! effective stiffness as the NODE-SURF per-node form. reflen is taken from the slave element
    ! (the constraint is integrated on the slave surface). NODE-SURF pairs are left untouched.
    if( contact%method == CONTACTS2S ) then
      ! Local slave reference-length sum and slave-face count (0 if this rank owns none).
      slave_reflen_sum = 0.0d0
      n_slave = 0
      if( associated(contact%slave_surf) ) n_slave = size(contact%slave_surf)
      do j = 1, n_slave
        nnode = size(contact%slave_surf(j)%nodes)
        do k = 1, nnode
          ctsurf = contact%slave_surf(j)%nodes(k)
          elem(1:3,k) = hecMESH%node(3*ctsurf-2:3*ctsurf)
        enddo
        call getElementCenter( contact%slave_surf(j)%etype, r0 )
        slave_reflen_sum = slave_reflen_sum + &
          getReferenceLength( contact%slave_surf(j)%etype, nnode, r0, elem )
      enddo
      ! Reduced outside the slave_surf>0 gate (ranks owning no slave face must still join the
      ! collective) but inside the method gate. Slave faces are owned by exactly one rank, so the
      ! SUM reproduces the serial value on every rank.
      call hecmw_allREDUCE_R1(hecMESH, slave_reflen_sum, hecmw_sum)
      call hecmw_allreduce_I1(hecMESH, n_slave, hecmw_sum)
      A_rep = 0.0d0
      if( n_slave > 0 ) then
        A_rep = ( slave_reflen_sum / dble(n_slave) ) ** 2
        if( A_rep > 0.0d0 ) contact%refStiff = contact%refStiff / A_rep
      endif
    endif

    ! Report penalty settings
    if (hecmw_comm_get_rank() == 0) then
      write(*,'(A,A,A,1pE12.3,A,1pE12.3,A,1pE12.3)') "  Contact [", &
        trim(contact%pair_name), "] set penalty: normal & tied ", &
        contact%nPenalty * contact%refStiff, ", tangential ", &
        contact%tPenalty * contact%refStiff, ", refStiff ", contact%refStiff
    endif

  end subroutine calc_contact_pair_refStiff

  !> \brief Assemble contact nodal force into residual vector (conMAT%B).
  subroutine assemble_contact_force_residual(nnode,ndLocal,id_lagrange,ctNForce,ctTForce,conMAT)
    integer(kind=kint), intent(in)       :: nnode              !< number of master nodes
    integer(kind=kint), intent(in)       :: ndLocal(nnode + 1) !< global node numbers
    integer(kind=kint), intent(in)       :: id_lagrange        !< Lagrange multiplier index (0 if none)
    real(kind=kreal), intent(in)         :: ctNForce((nnode+1)*3+1) !< normal contact force vector
    real(kind=kreal), intent(in)         :: ctTForce((nnode+1)*3+1) !< tangential contact force vector
    type(hecmwST_matrix), intent(inout)  :: conMAT             !< contact matrix

    integer(kind=kint) :: i, inod, idx

    do i = 1, nnode + 1
      inod = ndLocal(i)
      idx = (inod-1)*3+1
      conMAT%B(idx:idx+2) = conMAT%B(idx:idx+2) + ctNForce((i-1)*3+1:(i-1)*3+3) + ctTForce((i-1)*3+1:(i-1)*3+3)
    enddo

    ! Accumulate: several contributions can target the same Lagrange row (the caller
    ! zero-clears conMAT%B before the contact assembly).
    if( id_lagrange > 0 ) then
      conMAT%B(conMAT%NP*conMAT%NDOF+id_lagrange) = &
      conMAT%B(conMAT%NP*conMAT%NDOF+id_lagrange) + ctNForce((nnode+1)*3+1) + ctTForce((nnode+1)*3+1)
    endif

  end subroutine assemble_contact_force_residual

  !> \brief Accumulate contact nodal force into output arrays (CONT_NFORCE/CONT_FRIC).
  subroutine assemble_contact_force_output(nnode,ndLocal,ctNForce,ctTForce,cont_nforce,cont_fric)
    integer(kind=kint), intent(in)          :: nnode              !< number of master nodes
    integer(kind=kint), intent(in)          :: ndLocal(nnode + 1) !< global node numbers
    real(kind=kreal), intent(in)            :: ctNForce((nnode+1)*3+1) !< normal contact force vector
    real(kind=kreal), intent(in)            :: ctTForce((nnode+1)*3+1) !< tangential contact force vector
    real(kind=kreal), pointer, intent(inout) :: cont_nforce(:)    !< output normal force
    real(kind=kreal), pointer, optional, intent(inout) :: cont_fric(:) !< output friction force

    integer(kind=kint) :: i, inod, idx

    do i = 1, nnode + 1
      inod = ndLocal(i)
      idx = (inod-1)*3+1
      cont_nforce(idx:idx+2) = cont_nforce(idx:idx+2) + ctNForce((i-1)*3+1:(i-1)*3+3)
      if( present(cont_fric) ) cont_fric(idx:idx+2) = cont_fric(idx:idx+2) + ctTForce((i-1)*3+1:(i-1)*3+3)
    enddo

  end subroutine assemble_contact_force_output

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_contact_multiplier( ctAlgo, contact, coord, disp, ddisp, fcoeff, &
    hecMESH, hecLagMAT, gnt, ctchanged )
    integer(kind=kint), intent(in)       :: ctAlgo         !< contact algorithm
    type( tContact ), intent(inout)      :: contact        !< contact info
    real(kind=kreal), intent(in)         :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)         :: fcoeff         !< frictional coeff
    type(hecmwST_local_mesh), intent(in) :: hecMESH        !< mesh for allreduce
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< Lagrange matrix
    real(kind=kreal), intent(out)        :: gnt(2)         !< convergency information
    logical, intent(inout)               :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave, etype, master
    integer(kind=kint)  :: nn, i, g, cnt
    real(kind=kreal)    :: lgnt(2)
    integer(kind=kint)  :: ndLocal(l_max_elem_node+1)
    real(kind=kreal)    :: ctNForce(l_max_elem_node*3+3)
    real(kind=kreal)    :: ctTForce(l_max_elem_node*3+3)
    real(kind=kreal)    :: max_jump_ratio, jump_ratio_local
    real(kind=kreal)    :: mut_old, mut_new, threthold
    ! --- mortar (SURF-SURF) locals ---
    integer(kind=kint)  :: unique_count, r, a, nnode_s
    integer(kind=kint), allocatable :: maplist(:), master_idxs(:), sorted_idx(:)
    real(kind=kreal),   allocatable :: S(:), Ns_list(:,:), integrated_gaps(:), lambda_cur(:)
    ! per-node-within-group quantities; the per-node lambda_n drives the normal path
    real(kind=kreal),   allocatable :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:), lambda_node(:,:)
    real(kind=kreal)    :: mu, lambda_new
    ! --- friction: per-node slip/normal, per-node basis, return mapping ---
    real(kind=kreal),   allocatable :: Sigma_node(:,:,:), nacc_node(:,:,:)
    real(kind=kreal),   allocatable :: lam_t_cur(:,:,:)
    integer(kind=kint), allocatable :: fric_state_cur(:,:)
    real(kind=kreal)    :: nhat(3), t1(3), t2(3), nrm, Dxi(2)
    real(kind=kreal)    :: rho_t, alpha, that(2), lam_t_new(2)
    integer(kind=kint)  :: fstate

    cnt = 0
    lgnt(:) = 0.d0
    max_jump_ratio = 0.0d0
    if( contact%method == CONTACTS2S ) then
      ! ===== mortar multiplier update (per slave segment) =====
      mu = contact%nPenalty * contact%refStiff
      rho_t = contact%tPenalty * contact%refStiff   ! tangential penalty, used only if fcoeff/=0
      do i = 1, size(contact%slave_surf)
        if( contact%slave_surf(i)%state == CONTACTFREE ) cycle

        call getIntGap(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                       unique_count, maplist, master_idxs, S, Ns_list, integrated_gaps, &
                       Snode, Nsnode, gapwnode)

        nnode_s = size(contact%slave_surf(i)%nodes)
        allocate(sorted_idx(unique_count), lambda_cur(unique_count), lambda_node(nnode_s,unique_count))
        if( fcoeff /= 0.d0 ) then
          ! Resolve the tangent warm-start (working -> begin -> 0/STICK) before the working
          ! buffer is rebuilt below. lambda_n resolution is identical to the fcoeff=0 path.
          allocate(lam_t_cur(2,nnode_s,unique_count), fric_state_cur(nnode_s,unique_count))
          call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                  lambda_cur, lambda_node, lam_t_cur, fric_state_cur)
        else
          call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                  lambda_cur, lambda_node)
        endif

        ! Per-node augmented update: lambda_node(a,g) += mu*gapwnode(g,a), clamped at 0.
        ! lambda_cur(g) (the node sum) is refreshed for the friction cone radius below.
        do g = 1, unique_count
          do a = 1, nnode_s
            lambda_new = lambda_node(a,g) + (mu * gapwnode(g,a))
            if( lambda_new < 0.d0 ) lambda_new = 0.d0
            lambda_node(a,g) = lambda_new
          enddo
          lambda_cur(g) = sum( lambda_node(1:nnode_s,g) )

          ! Convergence tracking (group gap)
          lgnt(1) = lgnt(1) + integrated_gaps(g)
        enddo

        ! Rebuild the working buffer from the active master set (ascending); BEGIN left it empty.
        contact%slave_surf(i)%lam_work_n = unique_count
        do r = 1, unique_count
          contact%slave_surf(i)%lam_work_id(r)  = master_idxs(sorted_idx(r))
          contact%slave_surf(i)%lam_work_val(1:nnode_s,r) = lambda_node(1:nnode_s,sorted_idx(r))
        enddo
        cnt = cnt + 1

        ! --- friction tangent update: per slave node, project the mortar slip onto the
        !     per-node tangent frame, return-map with cone radius fcoeff*lambda_node(a,g),
        !     and write lambda_t / fric_state into the working buffer. ---
        if( fcoeff /= 0.d0 ) then
          allocate(Sigma_node(unique_count,nnode_s,3), nacc_node(unique_count,nnode_s,3))
          call getTangentSlip(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                              unique_count, maplist, master_idxs, Sigma_node, nacc_node)
          do g = 1, unique_count
            do a = 1, nnode_s
              nrm = sqrt( nacc_node(g,a,1)**2 + nacc_node(g,a,2)**2 + nacc_node(g,a,3)**2 )
              ! Project the per-node slip onto the per-node orthonormal frame, then Coulomb
              ! return-map. rho_t*Dxi matches the per-node mu*gapwnode area weighting (both
              ! node-tributary integrated), so the averaged back-distribution cancels the area.
              if( nrm < 1.d-30 ) then
                Dxi(1:2) = 0.d0
                nhat(1:3) = 0.d0
              else
                nhat(1:3) = nacc_node(g,a,1:3) / nrm
                call build_group_tangent_basis(nhat, t1, t2)
                Dxi(1) = dot_product(t1(1:3), Sigma_node(g,a,1:3))
                Dxi(2) = dot_product(t2(1:3), Sigma_node(g,a,1:3))
              endif
              fstate = fric_state_cur(a,g)
              call group_return_mapping(lam_t_cur(1:2,a,g), rho_t, Dxi, fcoeff, lambda_node(a,g), &
                                        fstrPR%eps_fric_band, lam_t_new, fstate, alpha, that, &
                                        update_state=.true.)
              lam_t_cur(1:2,a,g)  = lam_t_new(1:2)
              fric_state_cur(a,g) = fstate
            enddo
          enddo
          ! Write the tangent working buffer parallel to the rebuilt lambda_n (ascending master order).
          do r = 1, unique_count
            contact%slave_surf(i)%lam_work_t(1:2,1:nnode_s,r)  = lam_t_cur(1:2,1:nnode_s,sorted_idx(r))
            contact%slave_surf(i)%lam_work_fstate(1:nnode_s,r) = fric_state_cur(1:nnode_s,sorted_idx(r))
          enddo
          deallocate(Sigma_node, nacc_node, lam_t_cur, fric_state_cur)
        endif

        deallocate(maplist, master_idxs, S, Ns_list, integrated_gaps, sorted_idx, lambda_cur)
        deallocate(Snode, Nsnode, gapwnode, lambda_node)
      enddo
    else
      ! ===== NODE-SURF multiplier update (per slave node) =====
      do i = 1, size(contact%slave)
        if(.not. is_contact_active(contact%states(i)%state)) cycle   ! only STICK/SLIP

        slave = contact%slave(i)
        master = contact%states(i)%surface
        nn = size(contact%master(master)%nodes)
        etype = contact%master(master)%etype

        ndLocal(1) = slave
        ndLocal(2:nn+1) = contact%master(master)%nodes(1:nn)

        ! Update multiplier and calculate forces
        call updateContactMultiplier_Alag(contact%states(i), ndLocal(1:nn+1), coord, disp, ddisp, &
          contact%nPenalty * contact%refStiff, contact%tPenalty * contact%refStiff, &
          fcoeff, contact%master(master), lgnt, ctchanged, ctNForce, ctTForce, jump_ratio_local, contact%smoothing)

        ! Track maximum jump ratio
        max_jump_ratio = max(max_jump_ratio, jump_ratio_local)

        cnt = cnt + 1
      enddo
    endif

    if(cnt > 0) lgnt(:) = lgnt(:) / cnt
    gnt = gnt + lgnt
    
    call hecmw_allREDUCE_R1(hecMESH, max_jump_ratio, hecmw_max)
    
    ! Adjust tPenalty
    threthold = 100.d0
    if (max_jump_ratio > threthold) then
      mut_old = contact%tPenalty * contact%refStiff
      contact%tPenalty = contact%tPenalty * max(1.d0/dsqrt(threthold), 1.0d0/dsqrt(max_jump_ratio))
      mut_new = contact%tPenalty * contact%refStiff
      if (hecmw_comm_get_rank() == 0) then
        write(*,'(A,A,A,1pE12.3,A,1pE12.3,A)') "  Contact [", trim(contact%pair_name), &
          "] tangential penalty adjusted: ", mut_old, " -> ", mut_new, " (friction jump)"
      endif
    endif
      
  end subroutine update_contact_multiplier

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_tied_multiplier( contact, disp, ddisp, ctchanged )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    logical, intent(inout)            :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave, etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: dg(3), dgmax
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3)
    real(kind=kreal)    :: mu                          !< penalty from contact

    ! Calculate penalty from contact structure
    mu = contact%nPenalty * contact%refStiff

    do i= 1, size(contact%slave)
      if( .not. is_contact_active(contact%states(i)%state) ) cycle   ! only STICK/SLIP
      slave = contact%slave(i)
      edisp(1:3) = disp(3*slave-2:3*slave)+ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        edisp(3*j+1:3*j+3) = disp(3*iSS-2:3*iSS)+ddisp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      dg(1:3) = edisp(1:3)
      do j=1,nn
        dg(1:3) = dg(1:3)-shapefunc(j)*edisp(3*j+1:3*j+3)
      enddo

      contact%states(i)%multiplier(1:3) = contact%states(i)%multiplier(1:3) + mu*dg(1:3)

      ! check if tied constraint converged
      dgmax = 0.d0
      do j=1,(nn+1)*3
        dgmax = dgmax + dabs(edisp(j))
      enddo
      dgmax = dgmax/dble((nn+1)*3)
      do j=1,3
        if( dabs(dg(j))/dmax1(1.d0,dgmax) > 1.d-3 ) ctchanged = .true.
      enddo

    enddo
  end subroutine

  subroutine update_contact_TangentForce( contact )
    type( tContact ), intent(inout)   :: contact        !< contact info

    integer(kind=kint)  :: i

    do i= 1, size(contact%slave)
      if( .not. is_contact_active(contact%states(i)%state) ) then
        contact%states(i)%tangentForce(1:3) = 0.d0
        contact%states(i)%tangentForce_trial(1:3) = 0.d0
        contact%states(i)%tangentForce_final(1:3) = 0.d0
      else
        contact%states(i)%tangentForce(1:3) = contact%states(i)%tangentForce_final(1:3)
      end if
      contact%states(i)%tangentForce1(1:3) = contact%states(i)%tangentForce(1:3)
    enddo
  end subroutine update_contact_TangentForce

  !>\brief This subroutine calculates contact stiffness for each contact pair
  !! and assembles it into global stiffness matrix
  subroutine calcu_contact_stiffness_NodeSurf( ctAlgo, contact, coord, disp, ddisp, iter, lagrange_array, &
    conMAT, hecLagMAT)
    integer(kind=kint), intent(in)             :: ctAlgo          !< contact analysis algorithm
    type(tContact), intent(inout)              :: contact         !< contact info
    real(kind=kreal), intent(in)               :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)               :: disp(:)         !< displacement
    real(kind=kreal), intent(in)               :: ddisp(:)        !< displacement increment
    integer(kind=kint), intent(in)             :: iter            !< iteration number
    real(kind=kreal), intent(in)               :: lagrange_array(:) !< Lagrange multiplier array
    type(hecmwST_matrix), intent(inout)        :: conMAT          !< contact stiffness matrix
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT     !< Lagrange matrix

    integer(kind=kint) :: ctsurf, nnode, ndLocal(21), etype
    integer(kind=kint) :: j, k, algtype, id_lagrange
    real(kind=kreal)   :: lagrange
    real(kind=kreal)   :: stiffness((l_max_surface_node+1)*3+1, (l_max_surface_node+1)*3+1)
    real(kind=kreal)   :: elecoord(3, l_max_surface_node)  !< master node coordinates
    real(kind=kreal)   :: eledisp(l_max_surface_node*3+3)  !< element displacement increment for friction
    real(kind=kreal)   :: force(l_max_surface_node*3+3)    !< contact force direction
    logical            :: is_contact_active_flag, is_damping_active_flag

    algtype = contact%algtype

    do j = 1, size(contact%slave)

      ! stick or sliding contact is active
      is_contact_active_flag = is_contact_active(contact%states(j)%state)
      ! damping is active
      is_damping_active_flag = contact%states(j)%state == CONTACTNEAR .and. &
        &  is_damping_enabled(contact)

      if( .not. is_contact_active_flag .and. .not. is_damping_active_flag ) cycle

      ctsurf = contact%states(j)%surface
      etype = contact%master(ctsurf)%etype
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)

      ! Prepare master node coordinates for ALagrange (deformed configuration)
      do k = 1, nnode
        elecoord(1:3, k) = coord(3*ndLocal(k+1)-2:3*ndLocal(k+1)) + disp(3*ndLocal(k+1)-2:3*ndLocal(k+1))
      enddo

      if( is_contact_active_flag ) then

        if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then

          if( ctAlgo == kcaSLagrange ) then
            id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
            id_lagrange = id_lagrange + 1
            lagrange = lagrange_array(id_lagrange)
            call getContactStiffness_Slag(contact%states(j), contact%master(ctsurf), iter, &
              contact%tPenalty, contact%fcoeff, lagrange, stiffness, smoothing_type=contact%smoothing)

            ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
            call hecmw_mat_ass_contactlag(nnode, ndLocal, id_lagrange, contact%fcoeff, stiffness, conMAT, hecLagMAT)

          else if( ctAlgo == kcaALagrange ) then
            ! Build element displacement increment for consistent tangent evaluation
            eledisp(1:3) = ddisp(3*ndLocal(1)-2:3*ndLocal(1))
            do k = 1, nnode
              eledisp(k*3+1:k*3+3) = ddisp(3*ndLocal(k+1)-2:3*ndLocal(k+1))
            enddo
            call getContactStiffness_Alag(contact%states(j), contact%master(ctsurf), elecoord(:,1:nnode), &
              contact%nPenalty * contact%refStiff, contact%tPenalty * contact%refStiff, &
              contact%fcoeff, contact%symmetric, stiffness, force, &
              smoothing_type=contact%smoothing, edisp=eledisp(1:nnode*3+3), iter=iter, &
              slvpos=coord(3*ndLocal(1)-2:3*ndLocal(1)) + disp(3*ndLocal(1)-2:3*ndLocal(1)))

            ! Assemble contact stiffness matrix into global stiffness matrix
            call hecmw_mat_ass_elem(conMAT, nnode+1, ndLocal, stiffness)

          end if

        else if( algtype == CONTACTTIED ) then

          if( ctAlgo == kcaSLagrange ) then
            id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
            do k = 1, 3
              id_lagrange = id_lagrange + 1
              lagrange = lagrange_array(id_lagrange)

              call getTiedStiffness_Slag(contact%states(j), contact%master(ctsurf), k, stiffness, &
                contact%smoothing)
              ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
              call hecmw_mat_ass_contactlag(nnode, ndLocal, id_lagrange, 0.d0, stiffness, conMAT, hecLagMAT)
            enddo

          else if( ctAlgo == kcaALagrange ) then
            call getTiedStiffness_Alag(contact%states(j), contact%master(ctsurf), &
              contact%nPenalty * contact%refStiff, stiffness, force)

            ! Assemble contact stiffness matrix into global stiffness matrix
            call hecmw_mat_ass_elem(conMAT, nnode+1, ndLocal, stiffness)

          end if

        endif

      else if( is_damping_active_flag ) then
        call getDampingStiffness(contact%states(j), contact%master(ctsurf), &
          contact%damp_alpha * contact%refStiff, contact%damp_gact, &
          stiffness, smoothing_type=contact%smoothing)

        ! Assemble full damping stiffness for slave+master contact element
        call hecmw_mat_ass_elem(conMAT, nnode+1, ndLocal, stiffness)

      endif

    enddo

  end subroutine calcu_contact_stiffness_NodeSurf

  subroutine calcu_contact_stiffness_SurfSurf( ctAlgo, contact, coord, disp, ddisp, hecMAT )
    integer(kind=kint), intent(in)             :: ctAlgo          !< contact analysis algorithm
    type(tContact), intent(inout)              :: contact         !< contact info
    real(kind=kreal), intent(in)               :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)               :: disp(:)         !< displacement
    real(kind=kreal), intent(in)               :: ddisp(:)        !< displacement increment (NR)
    type(hecmwST_matrix), intent(inout)        :: hecMAT          !< global stiffness matrix

    integer(kind=kint) :: i, g, a, j, k, nnode_m, nnode_s, unique_count, ctsurf
    integer(kind=kint) :: ndLocal(l_max_surface_node+1)
    integer(kind=kint), allocatable :: maplist(:), master_idxs(:), sorted_idx(:)
    real(kind=kreal),   allocatable :: S(:), Ns_list(:,:), integrated_gaps(:), lambda_cur(:)
    real(kind=kreal),   allocatable :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:), lambda_node(:,:)
    real(kind=kreal) :: mu
    real(kind=kreal) :: stiff(24, 24), Ns(24)
    ! --- friction consistent tangent ---
    real(kind=kreal),   allocatable :: lam_t_cur(:,:,:), nacc_node(:,:,:), Sigma_node(:,:,:)
    integer(kind=kint), allocatable :: fric_state_cur(:,:)
    real(kind=kreal) :: rho_t, nhat(3), t1(3), t2(3), nrm, Dxi(2)
    real(kind=kreal) :: alpha, that(2), lam_t_new(2), Amat(2,2), T3d(3,2), M3(3,3)
    real(kind=kreal) :: Wb(l_max_surface_node+1), stiff_f(24,24)
    integer(kind=kint) :: na, nb, fstate

    mu = contact%nPenalty * contact%refStiff

    do i = 1, size(contact%slave_surf)
      if( contact%slave_surf(i)%state == CONTACTFREE ) cycle
      if( ctAlgo /= kcaALagrange ) cycle

      call getIntGap(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                     unique_count, maplist, master_idxs, S, Ns_list, integrated_gaps, &
                     Snode, Nsnode, gapwnode)

      nnode_s = size(contact%slave_surf(i)%nodes)
      allocate(sorted_idx(unique_count), lambda_cur(unique_count), lambda_node(nnode_s,unique_count))
      if( contact%fcoeff /= 0.d0 ) then
        allocate(lam_t_cur(2,nnode_s,unique_count), fric_state_cur(nnode_s,unique_count))
        call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                lambda_cur, lambda_node, lam_t_cur, fric_state_cur)
      else
        call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                lambda_cur, lambda_node)
      endif

      ! ===== Normal stiffness: per-node rank-1 sum, mu*Snode(g,a)*Nsnode(g,a)(x)Nsnode(g,a) =====
      do g = 1, unique_count
        ctsurf = master_idxs(g)
        nnode_m = size(contact%master(ctsurf)%nodes)
        ndLocal(1:nnode_s) = contact%slave_surf(i)%nodes(1:nnode_s)
        ndLocal(nnode_s+1:nnode_s+nnode_m) = contact%master(ctsurf)%nodes(1:nnode_m)
        do a = 1, nnode_s
          ! ALag contact condition per node: augmented force must be positive
          if( lambda_node(a,g)+mu*gapwnode(g,a) < 0.d0 ) cycle
          Ns = 0.d0
          Ns(1:(nnode_s+nnode_m)*3) = Nsnode(g, a, 1:(nnode_s+nnode_m)*3)
          do j = 1, (nnode_s+nnode_m)*3
            do k = 1, (nnode_s+nnode_m)*3
              stiff(j,k) = mu*Snode(g,a)*Ns(j)*Ns(k)
            enddo
          enddo
          ! The friction consistent tangent is assembled in the fcoeff-guarded block below.
          call hecmw_mat_ass_elem(hecMAT, nnode_s+nnode_m, ndLocal, stiff)
        enddo
      enddo

      ! ===== Friction consistent tangent (per slave node) =====
      ! Linearization of the per-node friction residual at the same live slip state:
      !   K_a(b,c) = Snode(g,a) * Wbar(a,b) * Wbar(a,c) * M3_a,  M3_a = T3d_a * A_a * T3d_a^T
      ! with Wbar(a,b) = Nsnode(g,a,b).nhat_a, the same map as the residual back-distribution.
      if( contact%fcoeff /= 0.d0 ) then
        rho_t = contact%tPenalty * contact%refStiff
        allocate(Sigma_node(unique_count,nnode_s,3), nacc_node(unique_count,nnode_s,3))
        call getTangentSlip(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                            unique_count, maplist, master_idxs, Sigma_node, nacc_node)
        do g = 1, unique_count
          ctsurf = master_idxs(g)
          nnode_m = size(contact%master(ctsurf)%nodes)
          ndLocal(1:nnode_s) = contact%slave_surf(i)%nodes(1:nnode_s)
          ndLocal(nnode_s+1:nnode_s+nnode_m) = contact%master(ctsurf)%nodes(1:nnode_m)
          do a = 1, nnode_s
            if( lambda_node(a,g) <= 0.d0 ) cycle   ! no per-node normal force -> no friction
            nrm = sqrt( nacc_node(g,a,1)**2 + nacc_node(g,a,2)**2 + nacc_node(g,a,3)**2 )
            if( nrm < 1.d-30 ) cycle
            nhat(1:3) = nacc_node(g,a,1:3) / nrm
            call build_group_tangent_basis(nhat, t1, t2)
            ! Live per-node slip projection and read-only return mapping (writes only OUT args).
            Dxi(1) = dot_product(t1(1:3), Sigma_node(g,a,1:3))
            Dxi(2) = dot_product(t2(1:3), Sigma_node(g,a,1:3))
            fstate = fric_state_cur(a,g)
            call group_return_mapping(lam_t_cur(1:2,a,g), rho_t, Dxi, contact%fcoeff, lambda_node(a,g), &
                                      fstrPR%eps_fric_band, lam_t_new, fstate, alpha, that)
            ! 2D tangent operator A (same construction as getContactStiffness_Alag).
            if( alpha <= 1.0d-20 ) then
              Amat = 0.d0
            else if( alpha >= 0.999d0 ) then
              Amat = 0.d0
              Amat(1,1) = rho_t
              Amat(2,2) = rho_t
            else
              Amat(1,1) = alpha * rho_t * (1.0d0 - that(1)*that(1))
              Amat(1,2) = alpha * rho_t * (-that(1)*that(2))
              Amat(2,1) = alpha * rho_t * (-that(2)*that(1))
              Amat(2,2) = alpha * rho_t * (1.0d0 - that(2)*that(2))
            endif
            ! M3 = T3d * A * T3d^T (3x3), T3d = [t1 t2]
            T3d(1:3,1) = t1(1:3)
            T3d(1:3,2) = t2(1:3)
            M3 = matmul( matmul(T3d, Amat), transpose(T3d) )

            ! Per-node averaged mortar weight of each node (= ANnode/Snode, recovered via nhat_a).
            do na = 1, nnode_s + nnode_m
              Wb(na) = dot_product(Nsnode(g,a,3*na-2:3*na), nhat(1:3))
            enddo
            ! K_a(b,c) = Snode(g,a) * Wbar(a,b) * Wbar(a,c) * M3_a
            stiff_f = 0.d0
            do nb = 1, nnode_s + nnode_m
              do na = 1, nnode_s + nnode_m
                do k = 1, 3
                  do j = 1, 3
                    stiff_f(3*na-3+j, 3*nb-3+k) = Snode(g,a) * Wb(na) * Wb(nb) * M3(j,k)
                  enddo
                enddo
              enddo
            enddo
            call hecmw_mat_ass_elem(hecMAT, nnode_s+nnode_m, ndLocal, stiff_f)
          enddo
        enddo
        deallocate(Sigma_node, nacc_node, lam_t_cur, fric_state_cur)
      endif

      deallocate(maplist, master_idxs, S, Ns_list, integrated_gaps, sorted_idx, lambda_cur)
      deallocate(Snode, Nsnode, gapwnode, lambda_node)
    enddo

  end subroutine calcu_contact_stiffness_SurfSurf

  !>\brief This subroutine calculates contact nodal force for each contact pair
  !! and assembles it into contact matrix and/or force arrays.
  !! When purpose == kctForResidual, forces are assembled into conMAT%B.
  !! When purpose == kctForOutput, forces are stored in CONT_NFORCE/CONT_FRIC using multiplier only (no penalty).
  subroutine calcu_contact_ndforce_NodeSurf( purpose, ctAlgo, contact, coord, disp, ddisp, lagrange_array, &
    conMAT, CONT_NFORCE, CONT_FRIC, hecLagMAT )
    integer(kind=kint), intent(in)       :: purpose         !< kctForResidual or kctForOutput
    integer(kind=kint), intent(in)       :: ctAlgo          !< contact analysis algorithm
    type( tContact ), intent(inout)      :: contact         !< contact info
    real(kind=kreal), intent(in)         :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)         !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)        !< disp till current substep
    real(kind=kreal), intent(in)         :: lagrange_array(:) !< Lagrange multiplier array
    type(hecmwST_matrix), intent(inout)  :: conMAT          !< contact matrix
    real(kind=kreal), pointer            :: CONT_NFORCE(:)  !< contact normal force
    real(kind=kreal), pointer            :: CONT_FRIC(:)    !< contact friction force
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT  !< Lagrange matrix

    integer(kind=kint) :: ctsurf, nnode, ndLocal(21)
    integer(kind=kint) :: j, k, algtype, id_lagrange
    real(kind=kreal)   :: ndCoord(21*3)
    real(kind=kreal)   :: ndu(21*3), ndDu(21*3)
    real(kind=kreal)   :: lagrange
    real(kind=kreal)   :: ctNForce(21*3+1)
    real(kind=kreal)   :: ctTForce(21*3+1)
    real(kind=kreal)   :: mu_n, mu_t
    logical            :: if_flag
    logical            :: is_contact_active_flag, is_damping_active_flag
    real(kind=kreal)   :: ctime, etime
    integer(kind=kint) :: if_type

    algtype = contact%algtype
    if_flag = (contact%if_type /= 0)
    if(if_flag)then
      ctime = contact%ctime
      etime = contact%if_etime
      if_type = contact%if_type
    end if

    do j = 1, size(contact%slave)

      ! stick or sliding contact is active
      is_contact_active_flag = is_contact_active(contact%states(j)%state)
      ! damping is active (residual only)
      is_damping_active_flag = (purpose == kctForResidual) .and. &
        contact%states(j)%state == CONTACTNEAR .and. is_damping_enabled(contact)

      if( .not. is_contact_active_flag .and. .not. is_damping_active_flag ) cycle

      ctsurf = contact%states(j)%surface
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)
      do k = 1, nnode+1
        ndDu((k-1)*3+1:(k-1)*3+3) = ddisp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
        ndu((k-1)*3+1:(k-1)*3+3) = disp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
        ndCoord((k-1)*3+1:(k-1)*3+3) = coord((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
      enddo

      if( is_contact_active_flag ) then

        if(if_flag) call set_shrink_factor(ctime, contact%states(j), etime, if_type)

        ! --- Determine penalty parameters: zero for output (multiplier-only)
        if( ctAlgo == kcaALagrange .and. purpose == kctForOutput ) then
          mu_n = 0.0d0
          mu_t = 0.0d0
        else
          mu_n = contact%nPenalty * contact%refStiff
          mu_t = contact%tPenalty * contact%refStiff
        endif

        if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
          ! Obtain contact nodal force vector of contact pair
          if(if_flag) call get_shrink_elemact_surf(contact%states(j),ndCoord, nnode)

          if( ctAlgo == kcaSLagrange ) then
            id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
            id_lagrange = id_lagrange + 1
            lagrange = lagrange_array(id_lagrange)
            call getContactNodalForce_Slag(contact%states(j),contact%master(ctsurf),ndCoord,ndDu,    &
              contact%tPenalty,contact%fcoeff,lagrange,ctNForce,ctTForce,.true.,contact%smoothing)

          else if( ctAlgo == kcaALagrange ) then
            id_lagrange = 0
            lagrange = 0.d0
            call getContactNodalForce_Alag(contact%states(j),contact%master(ctsurf),ndCoord,ndDu,    &
              mu_n, mu_t, contact%fcoeff,contact%symmetric,lagrange,ctNForce,ctTForce,.true.,contact%smoothing)

          end if

          ! Assemble contact force
          if( purpose == kctForResidual ) then
            call assemble_contact_force_residual(nnode,ndLocal,id_lagrange,ctNForce,ctTForce,conMAT)
          else
            call assemble_contact_force_output(nnode,ndLocal,ctNForce,ctTForce,CONT_NFORCE,CONT_FRIC)
          endif

        else if( algtype == CONTACTTIED ) then

          if( ctAlgo == kcaSLagrange ) then
            id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
            do k=1,3
              id_lagrange = id_lagrange + 1
              lagrange = lagrange_array(id_lagrange)
              contact%states(j)%multiplier(k) = lagrange

              call getTiedNodalForce_Slag(contact%states(j),contact%master(ctsurf),k,ndu, &
              &  lagrange,ctNForce,ctTForce,contact%smoothing)
              if( purpose == kctForResidual ) then
                call assemble_contact_force_residual(nnode,ndLocal,id_lagrange,ctNForce,ctTForce,conMAT)
              else
                call assemble_contact_force_output(nnode,ndLocal,ctNForce,ctTForce,CONT_NFORCE)
              endif
            end do

          else if( ctAlgo == kcaALagrange ) then
            id_lagrange = 0
            call getTiedNodalForce_Alag(contact%states(j),contact%master(ctsurf),ndu,    &
              mu_n, ctNForce,ctTForce)
            if( purpose == kctForResidual ) then
              call assemble_contact_force_residual(nnode,ndLocal,id_lagrange,ctNForce,ctTForce,conMAT)
            else
              call assemble_contact_force_output(nnode,ndLocal,ctNForce,ctTForce,CONT_NFORCE)
            endif

          end if

        endif

      else if( is_damping_active_flag ) then

        call getDampingNodalForce(contact%states(j), contact%master(ctsurf), ndDu, &
          contact%damp_alpha * contact%refStiff, contact%damp_gact, &
          ctNForce, ctTForce, smoothing_type=contact%smoothing)

        ! Assemble damping force for slave+master contact element
        call assemble_contact_force_residual(nnode,ndLocal,0,ctNForce,ctTForce,conMAT)

      endif

    enddo

  end subroutine calcu_contact_ndforce_NodeSurf

  !> \brief Compute contact nodal normal force for output from the stored contact
  !! multiplier, for the explicit dynamic method.
  !!
  !! Explicit counterpart of the kctForOutput path of calcu_contact_ndforce_NodeSurf.
  !! The forward-increment Lagrange corrector stores the converged contact normal
  !! force in states(:)%multiplier(1) but has no Lagrange matrix, so the multiplier
  !! is taken directly from the contact state instead of hecLagMAT%Lagrange. The same
  !! element routine (getContactNodalForce_Slag) and assembly (assemble_contact_force_output)
  !! as the implicit/static output path are reused, so CONT_NFORCE is produced identically.
  !! The explicit corrector stores the converged tangential force in
  !! states(:)%tangentForce_final; distribute it with the same relative-displacement
  !! mapping used by the implicit formulation.
  subroutine calcu_contact_ndforce_exp( contact, coord, disp, ddisp, CONT_NFORCE, CONT_FRIC )
    type( tContact ), intent(inout)      :: contact         !< contact info
    real(kind=kreal), intent(in)         :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)         !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)        !< disp increment of current substep
    real(kind=kreal), pointer            :: CONT_NFORCE(:)  !< contact normal force (output)
    real(kind=kreal), pointer            :: CONT_FRIC(:)    !< contact friction force (output)

    integer(kind=kint) :: ctsurf, nnode, ndLocal(21), j, k
    real(kind=kreal)   :: ndCoord(21*3), ndu(21*3), ndDu(21*3)
    real(kind=kreal)   :: ctNForce(21*3+1), ctTForce(21*3+1)
    real(kind=kreal)   :: Tm(3,3*(l_max_surface_node+1)), Tt(3,3*(l_max_surface_node+1))

    do j = 1, size(contact%slave)
      if( .not. is_contact_active(contact%states(j)%state) ) cycle

      ctsurf = contact%states(j)%surface
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)
      do k = 1, nnode+1
        ndDu((k-1)*3+1:(k-1)*3+3) = ddisp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
        ndu((k-1)*3+1:(k-1)*3+3) = disp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
        ndCoord((k-1)*3+1:(k-1)*3+3) = coord((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
      enddo

      call getContactNodalForce_Slag( contact%states(j), contact%master(ctsurf), ndCoord, ndDu, &
        0.d0, 0.d0, contact%states(j)%multiplier(1), ctNForce, ctTForce, .false., contact%smoothing )

      call computeTm_Tt( contact%states(j), contact%master(ctsurf), contact%fcoeff, &
        Tm, Tt, contact%smoothing )
      ctTForce(:) = 0.d0
      ctTForce(1:3*(nnode+1)) = -matmul(transpose(Tm(1:3,1:3*(nnode+1))), &
        contact%states(j)%tangentForce_final)

      call assemble_contact_force_output( nnode, ndLocal, ctNForce, ctTForce, CONT_NFORCE, CONT_FRIC )
    enddo

  end subroutine calcu_contact_ndforce_exp

  subroutine calcu_contact_ndforce_SurfSurf( purpose, ctAlgo, contact, coord, disp, ddisp, &
    conMAT, CONT_NFORCE, CONT_FRIC )
    integer(kind=kint), intent(in)       :: purpose         !< kctForResidual or kctForOutput
    integer(kind=kint), intent(in)       :: ctAlgo          !< contact analysis algorithm
    type( tContact ), intent(inout)      :: contact         !< contact info
    real(kind=kreal), intent(in)         :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)         !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)        !< disp till current substep
    type(hecmwST_matrix), intent(inout)  :: conMAT          !< contact matrix
    real(kind=kreal), pointer            :: CONT_NFORCE(:)  !< contact normal force
    real(kind=kreal), pointer            :: CONT_FRIC(:)    !< contact friction force

    integer(kind=kint) :: i, g, a, j, nd, nnode_m, nnode_s, unique_count, ctsurf
    integer(kind=kint) :: ndLocal(l_max_surface_node+1)
    integer(kind=kint), allocatable :: maplist(:), master_idxs(:), sorted_idx(:)
    real(kind=kreal),   allocatable :: S(:), Ns_list(:,:), integrated_gaps(:), lambda_cur(:)
    real(kind=kreal),   allocatable :: Snode(:,:), Nsnode(:,:,:), gapwnode(:,:), lambda_node(:,:)
    real(kind=kreal) :: mu, nrlforce
    real(kind=kreal) :: Ns(24)
    ! --- friction force back-distribution (live return mapping) ---
    real(kind=kreal),   allocatable :: lam_t_cur(:,:,:), nacc_node(:,:,:), Sigma_node(:,:,:)
    integer(kind=kint), allocatable :: fric_state_cur(:,:)
    real(kind=kreal) :: nhat(3), t1(3), t2(3), nrm, fvec(3), Wbar, fk(3)
    real(kind=kreal) :: rho_t, Dxi(2), alpha, that(2), lam_t_new(2)
    integer(kind=kint) :: fstate

    mu = contact%nPenalty * contact%refStiff

    do i = 1, size(contact%slave_surf)
      if( contact%slave_surf(i)%state == CONTACTFREE ) cycle
      if( ctAlgo /= kcaALagrange ) cycle

      call getIntGap(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                     unique_count, maplist, master_idxs, S, Ns_list, integrated_gaps, &
                     Snode, Nsnode, gapwnode)

      nnode_s = size(contact%slave_surf(i)%nodes)
      allocate(sorted_idx(unique_count), lambda_cur(unique_count), lambda_node(nnode_s,unique_count))
      if( contact%fcoeff /= 0.d0 ) then
        allocate(lam_t_cur(2,nnode_s,unique_count), fric_state_cur(nnode_s,unique_count))
        call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                lambda_cur, lambda_node, lam_t_cur, fric_state_cur)
      else
        call resolve_lambda_cur(contact%slave_surf(i), master_idxs, unique_count, nnode_s, sorted_idx, &
                                lambda_cur, lambda_node)
      endif

      ! ===== Normal force: per-node back-distribution =====
      ! nrlforce_a = lambda_node(a,g) + mu*gapwnode(g,a) for the residual, lambda_node(a,g) for output.
      do g = 1, unique_count
        ctsurf = master_idxs(g)
        nnode_m = size(contact%master(ctsurf)%nodes)
        ndLocal(1:nnode_s) = contact%slave_surf(i)%nodes(1:nnode_s)
        ndLocal(nnode_s+1:nnode_s+nnode_m) = contact%master(ctsurf)%nodes(1:nnode_m)
        do a = 1, nnode_s
          nrlforce = lambda_node(a,g) + mu*gapwnode(g,a)
          ! ALag contact condition per node: augmented force must be positive
          if( nrlforce < 0.d0 ) cycle
          Ns = 0.d0
          Ns(1:(nnode_s+nnode_m)*3) = Nsnode(g, a, 1:(nnode_s+nnode_m)*3)
          do j = 1, nnode_s + nnode_m
            nd = ndLocal(j)
            if( purpose == kctForResidual ) then
              conMAT%B(3*nd-2:3*nd) = conMAT%B(3*nd-2:3*nd) - nrlforce*Ns(3*j-2:3*j)
            else if ( purpose == kctForOutput ) then
              ! Output: multiplier only (converges to true contact force)
              CONT_NFORCE(3*nd-2:3*nd) = CONT_NFORCE(3*nd-2:3*nd) - lambda_node(a,g)*Ns(3*j-2:3*j)
            end if
          enddo
        enddo
      enddo

      ! ===== Friction force (per slave node): live return mapping, back-distributed =====
      ! trial = lam_t_warm(a) + rho_t*Dxi_live(a), projected onto the cone of radius lambda_node(a,g).
      ! The return mapping is read-only here (the augmentation update is the sole writer of the
      ! lambda_t / fric_state buffers). The resulting traction is distributed through the per-node
      ! mortar weight Wbar(a,j) = Nsnode(g,a,j).nhat_a, mirroring the normal back-distribution.
      ! Output (kctForOutput) keeps the frozen multiplier.
      if( contact%fcoeff /= 0.d0 ) then
        rho_t = contact%tPenalty * contact%refStiff
        allocate(Sigma_node(unique_count,nnode_s,3), nacc_node(unique_count,nnode_s,3))
        call getTangentSlip(contact%slave_surf(i), contact%master, coord, disp, ddisp, &
                            unique_count, maplist, master_idxs, Sigma_node, nacc_node)
        do g = 1, unique_count
          ctsurf = master_idxs(g)
          nnode_m = size(contact%master(ctsurf)%nodes)
          ndLocal(1:nnode_s) = contact%slave_surf(i)%nodes(1:nnode_s)
          ndLocal(nnode_s+1:nnode_s+nnode_m) = contact%master(ctsurf)%nodes(1:nnode_m)
          do a = 1, nnode_s
            if( lambda_node(a,g) <= 0.d0 ) cycle   ! no per-node normal force -> no friction
            nrm = sqrt( nacc_node(g,a,1)**2 + nacc_node(g,a,2)**2 + nacc_node(g,a,3)**2 )
            if( nrm < 1.d-30 ) cycle
            nhat(1:3) = nacc_node(g,a,1:3) / nrm
            call build_group_tangent_basis(nhat, t1, t2)
            if( purpose == kctForResidual ) then
              ! Live trial: project the live per-node mortar slip onto the per-node frame and return-map.
              Dxi(1) = dot_product(t1(1:3), Sigma_node(g,a,1:3))
              Dxi(2) = dot_product(t2(1:3), Sigma_node(g,a,1:3))
              fstate = fric_state_cur(a,g)
              call group_return_mapping(lam_t_cur(1:2,a,g), rho_t, Dxi, contact%fcoeff, lambda_node(a,g), &
                                        fstrPR%eps_fric_band, lam_t_new, fstate, alpha, that)
              fvec(1:3) = lam_t_new(1)*t1(1:3) + lam_t_new(2)*t2(1:3)
            else
              ! Output: frozen multiplier only (converges to the true friction force).
              fvec(1:3) = lam_t_cur(1,a,g)*t1(1:3) + lam_t_cur(2,a,g)*t2(1:3)
            end if

            do j = 1, nnode_s + nnode_m
              nd = ndLocal(j)
              ! per-node averaged mortar weight of node j (= ANnode/Snode, recovered via nhat_a)
              Wbar = dot_product(Nsnode(g,a,3*j-2:3*j), nhat(1:3))
              fk(1:3) = fvec(1:3) * Wbar
              if( purpose == kctForResidual ) then
                conMAT%B(3*nd-2:3*nd) = conMAT%B(3*nd-2:3*nd) - fk(1:3)
              else if( purpose == kctForOutput ) then
                CONT_FRIC(3*nd-2:3*nd) = CONT_FRIC(3*nd-2:3*nd) - fk(1:3)
              end if
            enddo
          enddo
        enddo
        deallocate(Sigma_node, nacc_node, lam_t_cur, fric_state_cur)
      endif

      deallocate(maplist, master_idxs, S, Ns_list, integrated_gaps, sorted_idx, lambda_cur)
      deallocate(Snode, Nsnode, gapwnode, lambda_node)
    enddo

  end subroutine calcu_contact_ndforce_SurfSurf

  !> Mortar: resolve the current per-node lambda of each active group of one slave surf.
  !> master_idxs is sorted ascending and merged against the ascending begin/working
  !> buffers with the rule: working hit -> working / else begin hit -> begin / else 0.
  !> lambda_cur(g) is the node sum of lambda_node(:,g).
  subroutine resolve_lambda_cur( surf, master_idxs, unique_count, nnode_s, sorted_idx, &
                                 lambda_cur, lambda_node, lam_t_cur, fric_state_cur )
    type(tContactSurf), intent(in)  :: surf
    integer(kind=kint), intent(in)  :: master_idxs(:)   !< group->masterID (get_unique_map output, unsorted)
    integer(kind=kint), intent(in)  :: unique_count
    integer(kind=kint), intent(in)  :: nnode_s          !< number of slave-surf nodes
    integer(kind=kint), intent(out) :: sorted_idx(:)    !< ascending rank r -> original group g
    real(kind=kreal),   intent(out) :: lambda_cur(:)    !< group-order g current lambda (= sum_a lambda_node)
    real(kind=kreal),   intent(out) :: lambda_node(:,:) !< (nnode_s, unique_count) per-node current lambda_n
    ! Optional friction warm-start: same reference rule as lambda_n (working -> begin ->
    ! default), riding the same merge. The lambda_n logic is unchanged.
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
      lambda_cur(g) = sum( lambda_node(1:nnode_s,g) )   ! group lambda_n
    enddo
  end subroutine resolve_lambda_cur

end module m_fstr_contact_assembly
