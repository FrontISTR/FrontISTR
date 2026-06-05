!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides a function to deal with prescribed displacement.

module m_fstr_AddBC
  implicit none
contains

  !>  Add Essential Boundary Conditions
  !------------------------------------------------------------------------------------------*
  subroutine fstr_AddBC(cstep,hecMESH,hecMAT,fstrSOLID,fstrPARAM,hecLagMAT,iter,conMAT,RHSvector)
    !------------------------------------------------------------------------------------------*
    use m_fstr
    use mContact
    use m_static_LIB_1d
    use m_static_LIB_shell, only: ShellRelativeRotationVector
    use m_utilities
    use m_fstr_TimeInc
    use m_fstr_NodalKinematics, only: fstr_mark_finite_rotation_nodes
    integer, intent(in)                  :: cstep !< current step
    type(hecmwST_local_mesh)             :: hecMESH !< hecmw mesh
    type(hecmwST_matrix)                 :: hecMAT !< hecmw matrix
    type(fstr_solid)                     :: fstrSOLID !< fstr_solid
    type(fstr_param)                     :: fstrPARAM !< analysis control parameters
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint)                   :: iter !< NR iterations
    type(hecmwST_matrix), optional       :: conMAT !< hecmw matrix for contact only
    real(kind=kreal), optional           :: RHSvector(:) !< only Right Hand Side vector

    integer(kind=kint) :: ig0, ig, ig1, ig_node, ityp, idofS, idofE, idof, idofS1, idofE1
    integer(kind=kint) :: idof1, iS0, iE0, iS1, iE1, ik, ik1, in
    real(kind=kreal)   :: RHS0, RHS, factor, factor0
    integer(kind=kint) :: ndof, grpid, istot
    integer(kind=kint) :: rot_dof_offset, rot_idofS, rot_idofE
    integer(kind=kint), allocatable :: shell_node_mode(:)
    logical :: target_found(3), use_shell_rotation_bc, node_in_group
    real(kind=kreal) :: theta_old(3), theta_target(3), theta_increment(3)

    !for rotation
    integer(kind=kint) :: n_rot, rid, jj_n_amp
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), cdiff0(3)
    real(kind=kreal)   :: cdisp(3), cddisp(3)

    !
    ndof = hecMAT%NDOF
    if( ndof >= 6 ) then
      allocate( shell_node_mode(hecMESH%n_node) )
      call fstr_mark_finite_rotation_nodes( hecMESH, fstrSOLID, ndof, shell_node_mode )
    endif

    n_rot = fstrSOLID%BOUNDARY_ngrp_rot
    if( n_rot > 0 ) call fstr_RotInfo_init(n_rot, rinfo)

    !   ----- Prescibed displacement Boundary Conditions
    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( iter>1 ) then
        factor=0.d0
      else
        jj_n_amp = fstrSOLID%BOUNDARY_ngrp_amp(ig0)
        if( jj_n_amp <= 0 ) then  ! Amplitude not defined
          factor0 = fstrSOLID%FACTOR(1)
          factor = fstrSOLID%FACTOR(2)
        else
          call table_amp(hecMESH,fstrSOLID,cstep,jj_n_amp,fstr_get_time(),factor0)
          call table_amp(hecMESH,fstrSOLID,cstep,jj_n_amp,fstr_get_time()+fstr_get_timeinc(),factor)
        endif
        factor = factor - factor0
        if(fstrSOLID%step_ctrl(cstep)%solution==stepVisco)then
          factor = 0.d0
          if(factor0 < 1.d-10) factor = 1.d0
        endif
      endif

      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      RHS0 = fstrSOLID%BOUNDARY_ngrp_val(ig0)
      !
      ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      idofS = ityp/10
      idofE = ityp - idofS*10
      !
      istot = fstrSOLID%BOUNDARY_ngrp_istot(ig0)
      !
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )

      if( fstrSOLID%BOUNDARY_ngrp_rotID(ig0) > 0 ) then ! setup rotation information
        rid = fstrSOLID%BOUNDARY_ngrp_rotID(ig0)
        if( .not. rinfo%conds(rid)%active ) then
          rinfo%conds(rid)%active = .true.
          rinfo%conds(rid)%center_ngrp_id = fstrSOLID%BOUNDARY_ngrp_centerID(ig0)
          rinfo%conds(rid)%torque_ngrp_id = ig
        endif
        if( istot == 0 ) then
          RHS= RHS0*factor
        else
          ! not implemented yet
          write(*,*) 'Error: rotational boundary cannot be specified with total value'
          call hecmw_abort( hecmw_comm_get_comm() )
        endif
        do idof=idofS,idofE
          if( idof>ndof ) then
            rinfo%conds(rid)%vec(idof-ndof) = RHS
          else
            rinfo%conds(rid)%vec(idof) = RHS
          endif
        enddo
        cycle
      endif

      !
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        !
        do idof = idofS, idofE
          use_shell_rotation_bc = .false.
          if( istot /= 0 .and. allocated(shell_node_mode) ) then
            rot_dof_offset = -1
            if( shell_node_mode(in) == 1 .and. ndof >= 6 .and. idof >= 4 .and. idof <= 6 ) then
              rot_dof_offset = 3
            endif
            if( rot_dof_offset >= 0 ) then
              theta_old(1:3) = fstrSOLID%unode_bak(ndof*(in-1)+rot_dof_offset+1: &
                ndof*(in-1)+rot_dof_offset+3)
              theta_target(1:3) = theta_old(1:3)
              target_found(1:3) = .false.
              do ig1 = 1, fstrSOLID%BOUNDARY_ngrp_tot
                if( .not. fstr_isBoundaryActive( fstrSOLID, fstrSOLID%BOUNDARY_ngrp_GRPID(ig1), cstep ) ) cycle
                if( fstrSOLID%BOUNDARY_ngrp_rotID(ig1) > 0 ) cycle
                if( fstrSOLID%BOUNDARY_ngrp_istot(ig1) == 0 ) cycle
                ig_node = fstrSOLID%BOUNDARY_ngrp_ID(ig1)
                iS1 = hecMESH%node_group%grp_index(ig_node-1) + 1
                iE1 = hecMESH%node_group%grp_index(ig_node  )
                node_in_group = .false.
                do ik1 = iS1, iE1
                  if( hecMESH%node_group%grp_item(ik1) == in ) then
                    node_in_group = .true.
                    exit
                  endif
                enddo
                if( .not. node_in_group ) cycle
                ityp = fstrSOLID%BOUNDARY_ngrp_type(ig1)
                idofS1 = ityp/10
                idofE1 = ityp - idofS1*10
                rot_idofS = rot_dof_offset + 1
                rot_idofE = rot_dof_offset + 3
                do idof1 = max(idofS1, rot_idofS), min(idofE1, rot_idofE)
                  theta_target(idof1-rot_dof_offset) = fstrSOLID%BOUNDARY_ngrp_val(ig1)
                  target_found(idof1-rot_dof_offset) = .true.
                enddo
              enddo
              if( all(target_found(1:3)) ) then
                ! finite rotation boundary increment
                theta_target(1:3) = theta_old(1:3) + factor*(theta_target(1:3) - theta_old(1:3))
                call ShellRelativeRotationVector( theta_old, theta_target, theta_increment )
                RHS = theta_increment(idof-rot_dof_offset)
                use_shell_rotation_bc = .true.
              endif
            endif
          endif
          if( .not. use_shell_rotation_bc ) then
            if( istot == 0 ) then
              RHS = RHS0*factor
            else
              RHS = (RHS0 - fstrSOLID%unode_bak(ndof*(in-1)+idof))*factor
            endif
          endif
          if(present(RHSvector)) then
            RHSvector(ndof*(in-1)+idof) = RHS
            ! write(6,*) 'BC: ', ndof*(in-1)+idof, RHS
            cycle
          endif
          if(present(conMAT)) then
            call hecmw_mat_ass_bc(hecMAT, in, idof, RHS, conMAT)
          else
            call hecmw_mat_ass_bc(hecMAT, in, idof, RHS)
          endif
          if( fstr_is_contact_active() .and. fstrPARAM%solution_type == kstSTATIC   &
              .and. fstrPARAM%contact_algo == kcaSLagrange ) then
            if(present(conMAT)) then
              call hecmw_mat_ass_bc_contactlag(conMAT,hecLagMAT,in,idof,RHS)
            else
              call hecmw_mat_ass_bc_contactlag(hecMAT,hecLagMAT,in,idof,RHS)
            endif
          endif

        enddo
      enddo
    enddo

    !Apply rotational boundary condition
    !need to fix!!
    do rid = 1, n_rot
      if( .not. rinfo%conds(rid)%active ) cycle
      cdiff = 0.d0
      cdiff0 = 0.d0
      cddisp = 0.d0

      if( factor > 0.d0 ) then
        ig = rinfo%conds(rid)%center_ngrp_id
        do idof = 1, ndof
          ccoord(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, hecMESH%node)
          cdisp(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, fstrSOLID%unode)
          cddisp(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, hecMAT%B)
        enddo
        ccoord(1:ndof) = ccoord(1:ndof) + cdisp(1:ndof)
      endif

      ig = rinfo%conds(rid)%torque_ngrp_id
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        if( factor > 0.d0 ) then
          cdiff0(1:ndof) = hecMESH%node(ndof*(in-1)+1:ndof*in)+fstrSOLID%unode(ndof*(in-1)+1:ndof*in)-ccoord(1:ndof)
          cdiff(1:ndof) = cdiff0(1:ndof)
          call rotate_3dvector_by_Rodrigues_formula(rinfo%conds(rid)%vec(1:ndof),cdiff(1:ndof))
        endif
        do idof = 1, ndof
          RHS = cdiff(idof)-cdiff0(idof)+cddisp(idof)
          if(present(RHSvector)) then
            RHSvector(ndof*(in-1)+idof) = RHS
            ! write(6,*) 'BC(rot): ', ndof*(in-1)+idof, RHS
            cycle
          endif
          if(present(conMAT)) then
            call hecmw_mat_ass_bc(hecMAT, in, idof, RHS, conMAT)
          else
            call hecmw_mat_ass_bc(hecMAT, in, idof, RHS)
          endif
          if( fstr_is_contact_active() .and. fstrPARAM%solution_type == kstSTATIC   &
              .and. fstrPARAM%contact_algo == kcaSLagrange ) then
            if(present(conMAT)) then
              call hecmw_mat_ass_bc_contactlag(conMAT,hecLagMAT,in,idof,RHS)
            else
              call hecmw_mat_ass_bc_contactlag(hecMAT,hecLagMAT,in,idof,RHS)
            endif
          endif
        enddo
      enddo
    enddo
    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)
    if( allocated(shell_node_mode) ) deallocate( shell_node_mode )

    !
    !   ------ Truss element Diagonal Modification
    call truss_diag_modify(hecMAT,hecMESH)

    !
    !   ------ Equation boundary conditions
    do ig0=1,fstrSOLID%n_fix_mpc
      if( fstrSOLID%mpc_const(ig0) == 0.d0 ) cycle
      ! we need to confirm if it is active in curr step here
      RHS = fstrSOLID%mpc_const(ig0)*factor
      hecMESH%mpc%mpc_const(ig0) = RHS
    enddo

  end subroutine fstr_AddBC

end module m_fstr_AddBC
