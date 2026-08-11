!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains functions to set velocity boundary condition in dynamic analysis
module m_dynamic_mat_ass_bc_vl
contains

  !C***
  !> Gather rotational velocity information (!VELOCITY, ROT_CENTER=) into rinfo.
  !> For each velocity BC with ROT_CENTER, the given DOF values are the angular
  !> velocity vector components (rad/s), multiplied by the amplitude at t_curr.
  !C***
  subroutine gather_rotvel_info(hecMESH, fstrSOLID, fstrDYNAMIC, t_curr, ndof, rinfo, cstep)
    use m_fstr
    use m_table_dyn
    implicit none
    type(hecmwST_local_mesh)             :: hecMESH
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_dynamic)                   :: fstrDYNAMIC
    real(kind=kreal), intent(in)         :: t_curr
    integer(kind=kint), intent(in)       :: ndof
    type(tRotInfo)                       :: rinfo
    integer(kind=kint), optional         :: cstep

    integer(kind=kint) :: ig0, ig, ityp, idofS, idofE, idof, rid, grpid, flag_u
    real(kind=kreal)   :: RHS, f_t

    flag_u = 2
    do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
      if( fstrSOLID%VELOCITY_ngrp_rotID(ig0) <= 0 ) cycle
      grpid = fstrSOLID%VELOCITY_ngrp_GRPID(ig0)
      if( present(cstep) ) then
        if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      endif

      RHS = fstrSOLID%VELOCITY_ngrp_val(ig0)
      call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, f_t, flag_u)
      RHS = RHS * f_t

      rid = fstrSOLID%VELOCITY_ngrp_rotID(ig0)
      ig  = fstrSOLID%VELOCITY_ngrp_ID(ig0)
      if( .not. rinfo%conds(rid)%active ) then
        rinfo%conds(rid)%active = .true.
        rinfo%conds(rid)%center_ngrp_id = fstrSOLID%VELOCITY_ngrp_centerID(ig0)
        rinfo%conds(rid)%torque_ngrp_id = ig
      endif

      ityp = fstrSOLID%VELOCITY_ngrp_type(ig0)
      idofS = ityp/10
      idofE = ityp - idofS*10
      do idof = idofS, idofE
        if( idof > ndof ) then
          rinfo%conds(rid)%vec(idof-ndof) = RHS
        else
          rinfo%conds(rid)%vec(idof) = RHS
        endif
      enddo
    enddo
  end subroutine gather_rotvel_info

  !C***
  !> Return the current-configuration coordinate of the rotation center node group.
  !C***
  subroutine get_rotcenter_coord(hecMESH, fstrSOLID, ig, ndof, ccoord)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh)       :: hecMESH
    type(fstr_solid)               :: fstrSOLID
    integer(kind=kint), intent(in) :: ig, ndof
    real(kind=kreal), intent(out)  :: ccoord(3)

    integer(kind=kint) :: idof
    real(kind=kreal)   :: cdisp(3)

    ccoord = 0.d0
    cdisp  = 0.d0
    do idof = 1, ndof
      ccoord(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, hecMESH%node)
      cdisp(idof)  = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, fstrSOLID%unode)
    enddo
    ccoord(1:ndof) = ccoord(1:ndof) + cdisp(1:ndof)
  end subroutine get_rotcenter_coord

  !C***
  !> This subrouitne set velocity boundary condition in dynamic analysis
  !C***

  subroutine DYNAMIC_MAT_ASS_BC_VL(cstep, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, hecLagMAT, t_curr, iter, conMAT)
    use m_fstr
    use m_table_dyn
    use mContact
    use m_utilities

    implicit none
    integer(kind=kint) :: cstep
    type(hecmwST_matrix)                 :: hecMAT
    type(hecmwST_local_mesh)             :: hecMESH
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_param)                     :: fstrPARAM !< analysis control parameters
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    real(kind=kreal)                     :: t_curr
    type(hecmwST_matrix), optional       :: conMAT

    integer, optional :: iter

    integer(kind=kint) :: ig0, ig, ityp, NDOF, iS0, iE0, ik, in, idofS, idofE, idof
    integer(kind=kint) :: flag_u, grpid
    real(kind=kreal)   :: b2, b3, b4, c1
    real(kind=kreal)   :: RHS, RHS0, f_t

    !for rotation
    integer(kind=kint) :: n_rot, rid
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), vnode(3), omega(3)

    if( fstrSOLID%VELOCITY_type == kbcInitial )return

    flag_u = 2

    if(dabs(fstrDYNAMIC%gamma) .lt. 1.0e-20) then
      if( hecMESH%my_rank == 0 ) then
        write(imsg,*) 'stop due to fstrDYNAMIC%gamma = 0'
      end if
      call hecmw_abort( hecmw_comm_get_comm())
    end if

    b2 = fstrDYNAMIC%t_delta   &
      *(fstrDYNAMIC%gamma-fstrDYNAMIC%beta)/fstrDYNAMIC%gamma
    b3 = fstrDYNAMIC%t_delta**2  &
      *(fstrDYNAMIC%gamma-2.0*fstrDYNAMIC%beta)    &
      /(2.0*fstrDYNAMIC%gamma)
    b4 = fstrDYNAMIC%t_delta*fstrDYNAMIC%beta/fstrDYNAMIC%gamma
    c1 = 2.0*fstrDYNAMIC%t_delta

    NDOF = hecMAT%NDOF

    n_rot = fstrSOLID%VELOCITY_ngrp_rot
    if( n_rot > 0 ) then
      call fstr_RotInfo_init(n_rot, rinfo)
      call gather_rotvel_info(hecMESH, fstrSOLID, fstrDYNAMIC, t_curr, NDOF, rinfo, cstep)
    endif

    !C=============================C
    !C-- implicit dynamic analysis
    !C=============================C
    if( fstrDYNAMIC%idx_eqa == 1 ) then

      do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
        if( fstrSOLID%VELOCITY_ngrp_rotID(ig0) > 0 ) cycle
        ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
        grpid = fstrSOLID%VELOCITY_ngrp_GRPID(ig0)
        if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
        RHS  = fstrSOLID%VELOCITY_ngrp_val(ig0)

        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, f_t, flag_u)
        RHS = RHS * f_t
        RHS0 = RHS

        ityp = fstrSOLID%VELOCITY_ngrp_type(ig0)

        idofS = ityp/10
        idofE = ityp - idofS*10

        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          do idof = idofS, idofE

            if( present(iter) ) then   ! increment
              if( iter>1 ) then
                RHS = 0.d0
              else
                RHS =              &
                  + b2*fstrDYNAMIC%VEL (NDOF*in-(NDOF-idof),1)     &
                  + b3*fstrDYNAMIC%ACC (NDOF*in-(NDOF-idof),1)     &
                  + b4*RHS0
              endif
            else
              RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),1)     &
                + b2*fstrDYNAMIC%VEL (NDOF*in-(NDOF-idof),1)     &
                + b3*fstrDYNAMIC%ACC (NDOF*in-(NDOF-idof),1)     &
                + b4*RHS0
            endif
            if(present(conMAT)) then
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS, conMAT)
            else
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS)
            endif
            if( fstr_is_contact_active() .and. fstrPARAM%contact_algo == kcaSLagrange  &
                .and. fstrPARAM%nlgeom .and. fstrDYNAMIC%idx_resp == 1 ) then
              if(present(conMAT)) then
                call hecmw_mat_ass_bc_contactlag(conMAT,hecLagMAT,in,idof,RHS)
              else
                call hecmw_mat_ass_bc_contactlag(hecMAT,hecLagMAT,in,idof,RHS)
              endif
            endif
          enddo
        enddo
      enddo

      !C-- rotational velocity boundary condition (implicit)
      do rid = 1, n_rot
        if( .not. rinfo%conds(rid)%active ) cycle
        omega(1:3) = rinfo%conds(rid)%vec(1:3)
        call get_rotcenter_coord(hecMESH, fstrSOLID, rinfo%conds(rid)%center_ngrp_id, NDOF, ccoord)
        ig = rinfo%conds(rid)%torque_ngrp_id
        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          cdiff = 0.d0
          cdiff(1:NDOF) = hecMESH%node(NDOF*(in-1)+1:NDOF*in) &
            + fstrSOLID%unode(NDOF*(in-1)+1:NDOF*in) - ccoord(1:NDOF)
          call cross_product(omega, cdiff, vnode)   ! v = omega x (x - x_center)
          do idof = 1, NDOF
            if( present(iter) ) then
              if( iter>1 ) then
                RHS = 0.d0
              else
                RHS = b2*fstrDYNAMIC%VEL(NDOF*in-(NDOF-idof),1)   &
                  + b3*fstrDYNAMIC%ACC(NDOF*in-(NDOF-idof),1)     &
                  + b4*vnode(idof)
              endif
            else
              RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),1)       &
                + b2*fstrDYNAMIC%VEL(NDOF*in-(NDOF-idof),1)       &
                + b3*fstrDYNAMIC%ACC(NDOF*in-(NDOF-idof),1)       &
                + b4*vnode(idof)
            endif
            if(present(conMAT)) then
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS, conMAT)
            else
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS)
            endif
            if( fstr_is_contact_active() .and. fstrPARAM%contact_algo == kcaSLagrange  &
                .and. fstrPARAM%nlgeom .and. fstrDYNAMIC%idx_resp == 1 ) then
              if(present(conMAT)) then
                call hecmw_mat_ass_bc_contactlag(conMAT,hecLagMAT,in,idof,RHS)
              else
                call hecmw_mat_ass_bc_contactlag(hecMAT,hecLagMAT,in,idof,RHS)
              endif
            endif
          enddo
        enddo
      enddo
      !C
      !C-- end of implicit dynamic analysis
      !C

      !C=============================C
      !C-- explicit dynamic analysis
      !C=============================C
    else if( fstrDYNAMIC%idx_eqa == 11 ) then
      !C
      do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
        if( fstrSOLID%VELOCITY_ngrp_rotID(ig0) > 0 ) cycle
        ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
        RHS  = fstrSOLID%VELOCITY_ngrp_val(ig0)

        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, f_t, flag_u)
        RHS = RHS * f_t
        RHS0 = RHS

        ityp = fstrSOLID%VELOCITY_ngrp_type(ig0)

        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        idofS = ityp/10
        idofE = ityp - idofS*10

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          do idof = idofS, idofE
            RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),3)    &
              + c1*RHS0
            hecMAT%B        (NDOF*in-(NDOF-idof)) = RHS
            fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof)) = 1.0d0
          end do
        enddo
      enddo

      !C-- rotational velocity boundary condition (explicit)
      do rid = 1, n_rot
        if( .not. rinfo%conds(rid)%active ) cycle
        omega(1:3) = rinfo%conds(rid)%vec(1:3)
        call get_rotcenter_coord(hecMESH, fstrSOLID, rinfo%conds(rid)%center_ngrp_id, NDOF, ccoord)
        ig = rinfo%conds(rid)%torque_ngrp_id
        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          cdiff = 0.d0
          cdiff(1:NDOF) = hecMESH%node(NDOF*(in-1)+1:NDOF*in) &
            + fstrSOLID%unode(NDOF*(in-1)+1:NDOF*in) - ccoord(1:NDOF)
          call cross_product(omega, cdiff, vnode)
          do idof = 1, NDOF
            RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),3) + c1*vnode(idof)
            hecMAT%B        (NDOF*in-(NDOF-idof)) = RHS
            fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof)) = 1.0d0
          enddo
        enddo
      enddo
      !C
      !C-- end of explicit dynamic analysis
      !C
    end if
    !
    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)

    return
  end subroutine DYNAMIC_MAT_ASS_BC_VL


  !C***
  !> This function sets initial condition of velocity
  !C***
  subroutine DYNAMIC_BC_INIT_VL(hecMESH, hecMAT, fstrSOLID ,fstrDYNAMIC, t_curr)
    use m_fstr
    use m_table_dyn
    use m_utilities

    implicit none
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC

    integer(kind=kint) :: NDOF, ig0, ig, ityp, iS0, iE0, ik, in, idofS, idofE, idof

    integer(kind=kint) :: flag_u, grpid
    real(kind=kreal)   :: RHS, f_t, t_curr

    !for rotation
    integer(kind=kint) :: n_rot, rid
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), vnode(3), omega(3)

    if( fstrSOLID%VELOCITY_type == kbcTransit )return

    flag_u = 2
    NDOF = hecMAT%NDOF

    n_rot = fstrSOLID%VELOCITY_ngrp_rot
    if( n_rot > 0 ) then
      call fstr_RotInfo_init(n_rot, rinfo)
      call gather_rotvel_info(hecMESH, fstrSOLID, fstrDYNAMIC, t_curr, NDOF, rinfo, 1)
    endif

    do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
      if( fstrSOLID%VELOCITY_ngrp_rotID(ig0) > 0 ) cycle
      ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
      RHS  = fstrSOLID%VELOCITY_ngrp_val(ig0)
      grpid = fstrSOLID%VELOCITY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, 1 ) ) cycle

      call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, f_t, flag_u)
      RHS = RHS * f_t

      ityp = fstrSOLID%VELOCITY_ngrp_type(ig0)

      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      idofS = ityp/10
      idofE = ityp - idofS*10

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)

        do idof = idofS, idofE
          fstrDYNAMIC%VEL (NDOF*in-(NDOF-idof),1) = RHS
        end do
      enddo
    enddo

    !C-- rotational velocity initial condition
    do rid = 1, n_rot
      if( .not. rinfo%conds(rid)%active ) cycle
      omega(1:3) = rinfo%conds(rid)%vec(1:3)
      call get_rotcenter_coord(hecMESH, fstrSOLID, rinfo%conds(rid)%center_ngrp_id, NDOF, ccoord)
      ig = rinfo%conds(rid)%torque_ngrp_id
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        cdiff = 0.d0
        cdiff(1:NDOF) = hecMESH%node(NDOF*(in-1)+1:NDOF*in) &
          + fstrSOLID%unode(NDOF*(in-1)+1:NDOF*in) - ccoord(1:NDOF)
        call cross_product(omega, cdiff, vnode)
        do idof = 1, NDOF
          fstrDYNAMIC%VEL (NDOF*in-(NDOF-idof),1) = vnode(idof)
        enddo
      enddo
    enddo

    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)

    return
  end subroutine DYNAMIC_BC_INIT_VL

  subroutine DYNAMIC_EXPLICIT_ASS_VL(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, t_curr, iter)
    use m_fstr
    use m_table_dyn
    use mContact
    use m_utilities

    implicit none
    type(hecmwST_matrix)                 :: hecMAT
    type(hecmwST_local_mesh)             :: hecMESH
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_dynamic)                   :: fstrDYNAMIC
    real(kind=kreal)                     :: t_curr
    integer, optional :: iter

    integer(kind=kint) :: ig0, ig, ityp, NDOF, iS0, iE0, ik, in, idofS, idofE, idof
    integer(kind=kint) :: flag_u
    real(kind=kreal)   :: b2, b3, b4, c1
    real(kind=kreal)   :: RHS, RHS0, f_t

    !for rotation
    integer(kind=kint) :: n_rot, rid
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: ccoord(3), cdiff(3), vnode(3), omega(3)

    if( fstrSOLID%VELOCITY_type == kbcInitial )return

    flag_u = 2

    if(dabs(fstrDYNAMIC%gamma) .lt. 1.0e-20) then
      if( hecMESH%my_rank == 0 ) then
        write(imsg,*) 'stop due to fstrDYNAMIC%gamma = 0'
      end if
      call hecmw_abort( hecmw_comm_get_comm())
    end if

    b2 = fstrDYNAMIC%t_delta   &
      *(fstrDYNAMIC%gamma-fstrDYNAMIC%beta)/fstrDYNAMIC%gamma
    b3 = fstrDYNAMIC%t_delta**2  &
      *(fstrDYNAMIC%gamma-2.0*fstrDYNAMIC%beta)    &
      /(2.0*fstrDYNAMIC%gamma)
    b4 = fstrDYNAMIC%t_delta*fstrDYNAMIC%beta/fstrDYNAMIC%gamma
    c1 = 2.0*fstrDYNAMIC%t_delta

    NDOF = hecMAT%NDOF

    n_rot = fstrSOLID%VELOCITY_ngrp_rot
    if( n_rot > 0 ) then
      call fstr_RotInfo_init(n_rot, rinfo)
      call gather_rotvel_info(hecMESH, fstrSOLID, fstrDYNAMIC, t_curr, NDOF, rinfo)
    endif

      !C
      do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
        if( fstrSOLID%VELOCITY_ngrp_rotID(ig0) > 0 ) cycle
        ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
        RHS  = fstrSOLID%VELOCITY_ngrp_val(ig0)

        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, f_t, flag_u)
        RHS = RHS * f_t
        RHS0 = RHS

        ityp = fstrSOLID%VELOCITY_ngrp_type(ig0)

        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        idofS = ityp/10
        idofE = ityp - idofS*10

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          do idof = idofS, idofE
            RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),3)    &
              + c1*RHS0
            hecMAT%B(NDOF*in-(NDOF-idof)) = RHS* fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof))
       !     fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof)) = 1.0d0
          end do
        enddo
      enddo

      !C-- rotational velocity boundary condition (explicit)
      do rid = 1, n_rot
        if( .not. rinfo%conds(rid)%active ) cycle
        omega(1:3) = rinfo%conds(rid)%vec(1:3)
        call get_rotcenter_coord(hecMESH, fstrSOLID, rinfo%conds(rid)%center_ngrp_id, NDOF, ccoord)
        ig = rinfo%conds(rid)%torque_ngrp_id
        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          cdiff = 0.d0
          cdiff(1:NDOF) = hecMESH%node(NDOF*(in-1)+1:NDOF*in) &
            + fstrSOLID%unode(NDOF*(in-1)+1:NDOF*in) - ccoord(1:NDOF)
          call cross_product(omega, cdiff, vnode)
          do idof = 1, NDOF
            RHS = fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),3) + c1*vnode(idof)
            hecMAT%B(NDOF*in-(NDOF-idof)) = RHS* fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof))
          enddo
        enddo
      enddo

    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)

  end subroutine DYNAMIC_EXPLICIT_ASS_VL

end module m_dynamic_mat_ass_bc_vl
