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
  use m_fstr_contact_interference
  use m_fstr_contact_elem_alag
  implicit none

  public :: calc_contact_pair_refStiff

contains

  !> Calculate reference stiffness for one contact pair
  subroutine calc_contact_pair_refStiff(contact, diag, ndof, hecMESH)
    type(tContact), intent(inout)        :: contact   !< contact pair
    real(kind=kreal), intent(in)         :: diag(:)   !< diagonal vector (size = ndof * np)
    integer(kind=kint), intent(in)       :: ndof      !< degrees of freedom
    type(hecmwST_local_mesh), intent(in) :: hecMESH   !< mesh
    
    integer(kind=kint) :: j, slave_node, master_node, nnode, ctsurf
    integer(kind=kint) :: idx_start, idx_end
    real(kind=kreal)   :: maxv
    
    maxv = 0.0d0
    
    ! Loop over slave nodes
    do j = 1, size(contact%slave)
      slave_node = contact%slave(j)
      idx_start = ndof * (slave_node - 1) + 1
      idx_end = ndof * slave_node
      maxv = max(maxv, maxval(diag(idx_start:idx_end)))
    enddo
    
    ! Loop over master surfaces and nodes
    do ctsurf = 1, size(contact%master)
      nnode = size(contact%master(ctsurf)%nodes)
      do j = 1, nnode
        master_node = contact%master(ctsurf)%nodes(j)
        idx_start = ndof * (master_node - 1) + 1
        idx_end = ndof * master_node
        maxv = max(maxv, maxval(diag(idx_start:idx_end)))
      enddo
    enddo
    
    ! Parallel reduction
    call hecmw_allREDUCE_R1(hecMESH, maxv, hecmw_max)
    
    ! Set reference stiffness for this contact pair
    contact%refStiff = maxv

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

    if( id_lagrange > 0 ) then
      conMAT%B(conMAT%NP*conMAT%NDOF+id_lagrange) = ctNForce((nnode+1)*3+1) + ctTForce((nnode+1)*3+1)
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
    integer(kind=kint)  :: nn, i, cnt
    real(kind=kreal)    :: lgnt(2)
    integer(kind=kint)  :: ndLocal(l_max_elem_node+1)
    real(kind=kreal)    :: ctNForce(l_max_elem_node*3+3)
    real(kind=kreal)    :: ctTForce(l_max_elem_node*3+3)
    real(kind=kreal)    :: max_jump_ratio, jump_ratio_local
    real(kind=kreal)    :: mut_old, mut_new, threthold

    cnt = 0
    lgnt(:) = 0.d0
    max_jump_ratio = 0.0d0
    
    do i = 1, size(contact%slave)
      if(contact%states(i)%state == CONTACTFREE) cycle   ! not in contact
      
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
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
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
      if( contact%states(i)%state==CONTACTFREE ) then
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

    algtype = contact%algtype

    do j = 1, size(contact%slave)

      if( contact%states(j)%state == CONTACTFREE ) cycle

      ctsurf = contact%states(j)%surface
      etype = contact%master(ctsurf)%etype
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)

      ! Prepare master node coordinates for ALagrange (deformed configuration)
      do k = 1, nnode
        elecoord(1:3, k) = coord(3*ndLocal(k+1)-2:3*ndLocal(k+1)) + disp(3*ndLocal(k+1)-2:3*ndLocal(k+1))
      enddo

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
              smoothing_type=contact%smoothing, edisp=eledisp(1:nnode*3+3))

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

    enddo

  end subroutine calcu_contact_stiffness_NodeSurf

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

      if( contact%states(j)%state == CONTACTFREE ) cycle
      if(if_flag) call set_shrink_factor(ctime, contact%states(j), etime, if_type)

      ctsurf = contact%states(j)%surface
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)
      do k = 1, nnode+1
        ndDu((k-1)*3+1:(k-1)*3+3) = ddisp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
        ndu((k-1)*3+1:(k-1)*3+3) = disp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
        ndCoord((k-1)*3+1:(k-1)*3+3) = coord((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
      enddo

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
            mu_n, mu_t, contact%fcoeff,lagrange,ctNForce,ctTForce,.true.,contact%smoothing)

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

    enddo

  end subroutine calcu_contact_ndforce_NodeSurf

end module m_fstr_contact_assembly
