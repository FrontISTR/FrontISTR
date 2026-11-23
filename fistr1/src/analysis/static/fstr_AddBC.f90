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
  subroutine fstr_AddBC(cstep,hecMESH,hecMAT,fstrSOLID,fstrPARAM,iter,hecEBC,RHSvector)
    !------------------------------------------------------------------------------------------*
    use m_fstr
    use mContact
    use m_static_LIB_1d
    use m_utilities
    use m_fstr_TimeInc
    use hecmw_ebc_defer
    integer, intent(in)                  :: cstep !< current step
    type(hecmwST_local_mesh)             :: hecMESH !< hecmw mesh
    type(hecmwST_matrix)                 :: hecMAT !< hecmw matrix
    type(fstr_solid)                     :: fstrSOLID !< fstr_solid
    type(fstr_param)                     :: fstrPARAM !< analysis control parameters
    integer(kind=kint)                   :: iter !< NR iterations
    type(hecmwST_ebc)                    :: hecEBC !< prescribed displacements to be imposed after the MPC processing
    real(kind=kreal), optional           :: RHSvector(:) !< only Right Hand Side vector

    integer(kind=kint) :: ig0, ig, ityp, idofS, idofE, idof, iS0, iE0, ik, in
    real(kind=kreal)   :: RHS0, RHS, factor, factor0
    integer(kind=kint) :: ndof, grpid, istot

    !for rotation
    integer(kind=kint) :: n_rot, rid, jj_n_amp, ntdof
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), cdiff0(3)
    real(kind=kreal)   :: cdisp(3), cddisp(3)

    !
    ndof = hecMAT%NDOF
    ntdof = min(ndof,3)

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
          factor0 = 1.d0
          call fstr_get_amplitude(hecMESH,fstrSOLID,cstep,jj_n_amp,fstr_get_time(),factor0)
          factor = 1.d0
          call fstr_get_amplitude(hecMESH,fstrSOLID,cstep,jj_n_amp,fstr_get_time()+fstr_get_timeinc(),factor)
        endif
        factor = factor - factor0
        if(fstrSOLID%step_ctrl(cstep)%amp_default_type==stepAmpStep)then
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
          rinfo%conds(rid)%vec(mod(idof-1,3)+1) = RHS
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
          call hecmw_ebc_set(hecEBC, in, idof, RHS)
          if(present(RHSvector)) then
            RHSvector(ndof*(in-1)+idof) = RHS
            ! write(6,*) 'BC: ', ndof*(in-1)+idof, RHS
            cycle
          endif
          if( fstr_is_contact_active() .and. fstrPARAM%solution_type == kstSTATIC   &
              .and. fstrPARAM%contact_algo == kcaSLagrange ) then
            call hecmw_mat_ass_bc_contactlag(hecMAT,in,idof,RHS)
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
        do idof = 1, ntdof
          ccoord(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, 3, idof, hecMESH%node)
          cdisp(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, fstrSOLID%unode)
          cddisp(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, hecEBC%val)
        enddo
        ccoord(1:ntdof) = ccoord(1:ntdof) + cdisp(1:ntdof)
      endif

      ig = rinfo%conds(rid)%torque_ngrp_id
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        if( factor > 0.d0 ) then
          cdiff0(1:ntdof) = hecMESH%node(3*(in-1)+1:3*(in-1)+ntdof) &
            + fstrSOLID%unode(ndof*(in-1)+1:ndof*(in-1)+ntdof) - ccoord(1:ntdof)
          cdiff(1:3) = cdiff0(1:3)
          call rotate_3dvector_by_Rodrigues_formula(rinfo%conds(rid)%vec(1:3),cdiff(1:3))
        endif
        do idof = 1, ntdof
          RHS = cdiff(idof)-cdiff0(idof)+cddisp(idof)
          call hecmw_ebc_set(hecEBC, in, idof, RHS)
          if(present(RHSvector)) then
            RHSvector(ndof*(in-1)+idof) = RHS
            ! write(6,*) 'BC(rot): ', ndof*(in-1)+idof, RHS
            cycle
          endif
          if( fstr_is_contact_active() .and. fstrPARAM%solution_type == kstSTATIC   &
              .and. fstrPARAM%contact_algo == kcaSLagrange ) then
            call hecmw_mat_ass_bc_contactlag(hecMAT,in,idof,RHS)
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
