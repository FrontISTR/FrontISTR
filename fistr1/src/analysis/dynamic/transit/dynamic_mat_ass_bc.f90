!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains functions to set displacement boundary condition in dynamic analysis

module m_dynamic_mat_ass_bc

contains


  !>  This subroutine setup disp bundary condition
  subroutine DYNAMIC_MAT_ASS_BC(hecMESH, hecMAT, fstrSOLID ,fstrDYNAMIC, fstrPARAM, fstrMAT, iter, conMAT)
    use m_fstr
    use m_table_dyn
    use fstr_matrix_con_contact
    use m_addContactStiffness
    use mContact
    use m_utilities

    implicit none
    type(hecmwST_matrix)                 :: hecMAT
    type(hecmwST_local_mesh)             :: hecMESH
    type(fstr_solid)                     :: fstrSOLID
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_param)                     :: fstrPARAM !< analysis control parameters
    type(fstrST_matrix_contact_lagrange) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, optional                    :: iter
    type(hecmwST_matrix), optional       :: conMAT

    integer(kind=kint) :: ig0, ig, ityp, NDOF, iS0, iE0, ik, in, idofS, idofE, idof

    integer(kind=kint) :: flag_u
    real(kind=kreal)   :: RHS, f_t, f_t1

    !for rotation
    integer(kind=kint) :: n_rot, rid, n_nodes
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: theta, normal(3), direc(3), ccoord(3), cdiff(3), cdiff0(3)
    real(kind=kreal)   :: cdisp(3), cddisp(3)
    !
    ndof = hecMAT%NDOF
    n_rot = fstrSOLID%BOUNDARY_ngrp_rot
    if( n_rot > 0 ) call fstr_RotInfo_init(n_rot, rinfo)
    fstrSOLID%REACTION = 0.d0

    flag_u = 1
    !C=============================C
    !C-- implicit dynamic analysis
    !C=============================C
    if( fstrDYNAMIC%idx_eqa == 1 ) then

      do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        RHS  = fstrSOLID%BOUNDARY_ngrp_val(ig0)

        if( present(iter) ) then
          if( iter>1 ) then
            RHS=0.d0
          else
            fstrDYNAMIC%i_step = fstrDYNAMIC%i_step-1
            fstrDYNAMIC%t_curr = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
            call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t1, flag_u)
            fstrDYNAMIC%i_step = fstrDYNAMIC%i_step+1
            fstrDYNAMIC%t_curr = fstrDYNAMIC%t_curr + fstrDYNAMIC%t_delta
            call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
            RHS = RHS * (f_t-f_t1)
          endif
        else
          call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
          RHS = RHS * f_t
        endif

        ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)
        idofS = ityp/10
        idofE = ityp - idofS*10

        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )

        if( fstrSOLID%BOUNDARY_ngrp_rotID(ig0) > 0 ) then ! setup rotation information
          rid = fstrSOLID%BOUNDARY_ngrp_rotID(ig0)
          if( .not. rinfo%conds(rid)%active ) then
            rinfo%conds(rid)%active = .true.
            rinfo%conds(rid)%center_ngrp_id = fstrSOLID%BOUNDARY_ngrp_centerID(ig0)
            rinfo%conds(rid)%torque_ngrp_id = ig
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

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)

          do idof = idofS, idofE
            if(present(conMAT)) then
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS, conMAT)
            else
              call hecmw_mat_ass_bc(hecMAT, in, idof, RHS)
            endif
            if( fstr_is_contact_active() .and. fstrPARAM%contact_algo == kcaSLagrange  &
                .and. fstrPARAM%nlgeom .and. fstrDYNAMIC%idx_resp == 1 )  then
              if(present(conMAT)) then
                call fstr_mat_ass_bc_contact(conMAT,fstrMAT,in,idof,RHS)
              else
                call fstr_mat_ass_bc_contact(hecMAT,fstrMAT,in,idof,RHS)
              endif
            endif

            !for output reaction force
            fstrSOLID%REACTION(ndof*(in-1)+idof) = fstrSOLID%QFORCE(ndof*(in-1)+idof)
          enddo
        enddo

      enddo

      !Apply rotational boundary condition
      do rid = 1, n_rot
        if( .not. rinfo%conds(rid)%active ) cycle
        cdiff = 0.d0
        cdiff0 = 0.d0
        cddisp = 0.d0

        if( f_t > 0.d0 ) then
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
          if( f_t > 0.d0 ) then
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
                call fstr_mat_ass_bc_contact(conMAT,fstrMAT,in,idof,RHS)
              else
                call fstr_mat_ass_bc_contact(hecMAT,fstrMAT,in,idof,RHS)
              endif
            endif

            !for output reaction force
            fstrSOLID%REACTION(ndof*(in-1)+idof) = fstrSOLID%QFORCE(ndof*(in-1)+idof)
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
      NDOF = hecMAT%NDOF
      do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        RHS  = fstrSOLID%BOUNDARY_ngrp_val(ig0)

        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
        RHS = RHS * f_t

        ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)

        iS0 = hecMESH%node_group%grp_index(ig-1) + 1
        iE0 = hecMESH%node_group%grp_index(ig  )
        idofS = ityp/10
        idofE = ityp - idofS*10

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)

          do idof = idofS, idofE
            hecMAT%B        (NDOF*in-(NDOF-idof)) = RHS
            fstrDYNAMIC%VEC1(NDOF*in-(NDOF-idof)) = 1.0d0

            !for output reaction force
            fstrSOLID%REACTION(NDOF*(in-1)+idof) = fstrSOLID%QFORCE(NDOF*(in-1)+idof)
          end do
        enddo
      enddo
      !C
      !C-- end of explicit dynamic analysis
      !C
    end if

    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)

  end subroutine DYNAMIC_MAT_ASS_BC


  !C***
  !> This subroutine setup initial condition of displacement
  !C***
  subroutine DYNAMIC_BC_INIT(hecMESH, hecMAT, fstrSOLID ,fstrDYNAMIC)
    use m_fstr
    use m_table_dyn

    implicit none
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC

    integer(kind=kint) :: NDOF, ig0, ig, ityp, iS0, iE0, ik, in, idofS, idofE, idof
    integer(kind=kint) :: flag_u
    real(kind=kreal)   :: RHS, f_t

    flag_u = 1
    NDOF = hecMAT%NDOF

    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      RHS  = fstrSOLID%BOUNDARY_ngrp_val(ig0)


      call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
      RHS = RHS * f_t

      ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)

      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      idofS = ityp/10
      idofE = ityp - idofS*10

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)

        do idof = idofS, idofE
          fstrDYNAMIC%DISP(NDOF*in-(NDOF-idof),1) = RHS
        end do
      enddo
    enddo

    return
  end subroutine DYNAMIC_BC_INIT

end module m_dynamic_mat_ass_bc
