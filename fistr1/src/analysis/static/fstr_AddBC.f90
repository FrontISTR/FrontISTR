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
  subroutine fstr_AddBC(cstep,hecMESH,hecMAT,fstrSOLID,fstrPARAM,hecLagMAT,iter,conMAT)
    !------------------------------------------------------------------------------------------*
    use m_fstr
    use mContact
    use m_static_LIB_1d
    use m_utilities
    integer, intent(in)                  :: cstep !< current step
    type(hecmwST_local_mesh)             :: hecMESH !< hecmw mesh
    type(hecmwST_matrix)                 :: hecMAT !< hecmw matrix
    type(fstr_solid)                     :: fstrSOLID !< fstr_solid
    type(fstr_param)                     :: fstrPARAM !< analysis control parameters
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint)                   :: iter !< NR iterations
    type(hecmwST_matrix), optional       :: conMAT !< hecmw matrix for contact only

    integer(kind=kint) :: ig0, ig, ityp, idofS, idofE, idof, iS0, iE0, ik, in
    real(kind=kreal)   :: RHS0, RHS, factor
    integer(kind=kint) :: ndof, grpid, istot

    !for rotation
    integer(kind=kint) :: n_rot, rid
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), cdiff0(3)
    real(kind=kreal)   :: cdisp(3), cddisp(3)

    !
    ndof = hecMAT%NDOF
    factor = fstrSOLID%FACTOR(2)-fstrSOLID%FACTOR(1)

    if( cstep<=fstrSOLID%nstep_tot .and. fstrSOLID%step_ctrl(cstep)%solution==stepVisco ) then
      factor = 0.d0
      if( fstrSOLID%FACTOR(1) < 1.d-10 ) factor = 1.d0
    endif
    if( iter>1 ) factor=0.d0

    n_rot = fstrSOLID%BOUNDARY_ngrp_rot
    if( n_rot > 0 ) call fstr_RotInfo_init(n_rot, rinfo)

    !   ----- Prescibed displacement Boundary Conditions
    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
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
          if( istot == 0 ) then
            RHS = RHS0*factor
          else
            RHS = (RHS0 - fstrSOLID%unode_bak(ndof*(in-1)+idof))*factor
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
