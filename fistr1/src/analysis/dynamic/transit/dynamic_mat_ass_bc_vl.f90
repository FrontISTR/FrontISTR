!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains functions to set velocity boundary condition in dynamic analysis
module m_dynamic_mat_ass_bc_vl
contains

  !C***
  !> This subrouitne set velocity boundary condition in dynamic analysis
  !C***

  subroutine DYNAMIC_MAT_ASS_BC_VL(cstep, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, hecLagMAT, t_curr, iter, conMAT)
    use m_fstr
    use m_table_dyn
    use mContact

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

    !C=============================C
    !C-- implicit dynamic analysis
    !C=============================C
    if( fstrDYNAMIC%idx_eqa == 1 ) then

      do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
        ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
        grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
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

            !for output reaction force
            fstrSOLID%REACTION(NDOF*(in-1)+idof) = fstrSOLID%QFORCE(NDOF*(in-1)+idof)
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

            !for output reaction force
            fstrSOLID%REACTION(NDOF*(in-1)+idof) = fstrSOLID%QFORCE(NDOF*(in-1)+idof)
          end do
        enddo
      enddo
      !C
      !C-- end of explicit dynamic analysis
      !C
    end if
    !
    return
  end subroutine DYNAMIC_MAT_ASS_BC_VL


  !C***
  !> This function sets initial condition of velocity
  !C***
  subroutine DYNAMIC_BC_INIT_VL(hecMESH, hecMAT, fstrSOLID ,fstrDYNAMIC, t_curr)
    use m_fstr
    use m_table_dyn

    implicit none
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC

    integer(kind=kint) :: NDOF, ig0, ig, ityp, iS0, iE0, ik, in, idofS, idofE, idof

    integer(kind=kint) :: flag_u
    real(kind=kreal)   :: RHS, f_t, t_curr

    if( fstrSOLID%VELOCITY_type == kbcTransit )return

    flag_u = 2
    NDOF = hecMAT%NDOF

    do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
      ig   = fstrSOLID%VELOCITY_ngrp_ID(ig0)
      RHS  = fstrSOLID%VELOCITY_ngrp_val(ig0)


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

    return
  end subroutine DYNAMIC_BC_INIT_VL

  subroutine DYNAMIC_EXPLICIT_ASS_VL(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, t_curr, iter)
    use m_fstr
    use m_table_dyn
    use mContact

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

      !C
      do ig0 = 1, fstrSOLID%VELOCITY_ngrp_tot
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

            !for output reaction force
            fstrSOLID%REACTION(NDOF*(in-1)+idof) = fstrSOLID%QFORCE(NDOF*(in-1)+idof)
          end do
        enddo
      enddo

  end subroutine DYNAMIC_EXPLICIT_ASS_VL

end module m_dynamic_mat_ass_bc_vl
