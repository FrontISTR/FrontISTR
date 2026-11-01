!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize, allocate and finalize
!! fstr_solid, and to read the cards of solid analysis into it. Boundary
!! conditions are appended, and the contact pairs already read are checked and
!! initialized, independently of the input file format.
!! header_name is the card name shown in error messages. A blank name of
!! amplitude or rotation center means that it is not given.
module m_fstr_solid_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , elem_grp_name_to_id_ex &
    , surf_grp_name_to_id_ex &
    , hecmw_streqr &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_ctrl_err_stop &
    , amp_name_to_id &
    , dload_grp_name_to_id_ex &
    , fstr_setup_util_err_stop
  use fstr_ctrl_common
  use fstr_ctrl_static
  use fstr_ctrl_dynamic
  use fstr_ctrl_material
  use fstr_ctrl_util_f
  implicit none

contains

  !> Initializer of structure fstr_solid
  subroutine fstr_solid_init( hecMESH, fstrSOLID )
    use m_fstr
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_solid)                :: fstrSOLID

    fstrSOLID%file_type  = kbcfFSTR

    fstrSOLID%BOUNDARY_ngrp_tot = 0
    fstrSOLID%BOUNDARY_ngrp_rot = 0
    fstrSOLID%CLOAD_ngrp_tot    = 0
    fstrSOLID%CLOAD_ngrp_rot    = 0
    fstrSOLID%DLOAD_ngrp_tot    = 0
    fstrSOLID%DLOAD_follow      = 1
    fstrSOLID%TEMP_ngrp_tot     = 0
    fstrSOLID%SPRING_ngrp_tot   = 0
    fstrSOLID%TEMP_irres        = 0
    fstrSOLID%TEMP_tstep        = 1
    fstrSOLID%TEMP_interval     = 1
    fstrSOLID%TEMP_rtype        = 1
    fstrSOLID%TEMP_factor       = 1.d0
    fstrSOLID%VELOCITY_ngrp_tot = 0
    fstrSOLID%VELOCITY_ngrp_rot = 0
    fstrSOLID%ACCELERATION_ngrp_tot = 0
    fstrSOLID%COUPLE_ngrp_tot   = 0

    fstrSOLID%restart_nout= 0
    fstrSOLID%is_smoothing_active = .false.
    fstrSOLID%has_finite_rotation_kinematics = .false.
    fstrSOLID%finite_rotation_state_ready = .false.

    fstrSOLID%elemopt361 = 0
    fstrSOLID%AutoINC_stat = 0
    fstrSOLID%CutBack_stat = 0
    fstrSOLID%NRstat_i(:) = 0
    fstrSOLID%NRstat_r(:) = 0.d0

  end subroutine fstr_solid_init

  !> Initializer of structure fstr_solid
  subroutine fstr_solid_alloc( hecMESH, fstrSOLID )
    use m_fstr
    use m_fstr_FiniteRotationKinematics, only: fstr_has_finite_rotation_kinematics, &
      fstr_mark_finite_rotation_nodes
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_solid)                :: fstrSOLID

    integer :: ndof, ntotal, ierror

    ndof=hecMESH%n_dof
    ntotal=ndof*hecMESH%n_node
    fstrSOLID%has_finite_rotation_kinematics = &
      fstr_has_finite_rotation_kinematics( hecMESH, fstrSOLID )
    fstrSOLID%finite_rotation_state_ready = .false.

    allocate ( fstrSOLID%GL( ntotal )          ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, GL>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%GL0( ntotal )          ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, GL0>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%EFORCE( ntotal )      ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, EFORCE>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    !        allocate ( fstrSOLID%TOTAL_DISP( ntotal )  ,STAT=ierror )
    !            if( ierror /= 0 ) then
    !              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TOTAL_DISP>'
    !              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
    !              call flush(idbg)
    !              call hecmw_abort( hecmw_comm_get_comm())
    !            end if
    allocate ( fstrSOLID%unode( ntotal )  ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, unode>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%unode_bak( ntotal )  ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, unode>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%dunode( ntotal )  ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, dunode>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%ddunode( ntotal )  ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, ddunode>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( fstrSOLID%has_finite_rotation_kinematics ) then
      allocate ( fstrSOLID%shell_node_mode( hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_node_mode>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      call fstr_mark_finite_rotation_nodes( hecMESH, fstrSOLID, ndof, fstrSOLID%shell_node_mode )
      allocate ( fstrSOLID%shell_rot_state( hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_rot_state>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_ref_triad( 9*hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_ref_triad>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_triad( 9*hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_triad>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_triad_bak( 9*hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_triad_bak>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_dtriad( 9*hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_dtriad>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_drill( hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_drill>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_drill_bak( hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_drill_bak>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate ( fstrSOLID%shell_ddrill( hecMESH%n_node )  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, shell_ddrill>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    end if
    allocate ( fstrSOLID%QFORCE( ntotal )      ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, QFORCE>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%DFORCE( ntotal )      ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, DFORCE>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate ( fstrSOLID%QFORCE_bak( ntotal )      ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <FSTR_SOLID, QFORCE_bak>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if

    fstrSOLID%GL(:)=0.d0
    fstrSOLID%GL0(:)=0.d0
    !        fstrSOLID%TOTAL_DISP(:)=0.d0
    fstrSOLID%unode(:)      = 0.d0
    fstrSOLID%unode_bak(:)  = 0.d0
    fstrSOLID%dunode(:)     = 0.d0
    fstrSOLID%ddunode(:)    = 0.d0
    if( associated(fstrSOLID%shell_rot_state) ) fstrSOLID%shell_rot_state(:) = 0
    if( associated(fstrSOLID%shell_ref_triad) ) fstrSOLID%shell_ref_triad(:) = 0.d0
    if( associated(fstrSOLID%shell_triad) )     fstrSOLID%shell_triad(:)     = 0.d0
    if( associated(fstrSOLID%shell_triad_bak) ) fstrSOLID%shell_triad_bak(:) = 0.d0
    if( associated(fstrSOLID%shell_dtriad) )    fstrSOLID%shell_dtriad(:)    = 0.d0
    if( associated(fstrSOLID%shell_drill) )     fstrSOLID%shell_drill(:)     = 0.d0
    if( associated(fstrSOLID%shell_drill_bak) ) fstrSOLID%shell_drill_bak(:) = 0.d0
    if( associated(fstrSOLID%shell_ddrill) )    fstrSOLID%shell_ddrill(:)    = 0.d0
    fstrSOLID%QFORCE(:)     = 0.d0
    fstrSOLID%QFORCE_bak(:) = 0.d0
    ! DFORCE is assembled only by implicit dynamic analysis, but the convergence
    ! check reads it in every analysis type.
    fstrSOLID%DFORCE(:)     = 0.d0
    fstrSOLID%FACTOR( 1:2 ) = 0.d0

    ! for MPC
    fstrSOLID%n_fix_mpc = hecMESH%mpc%n_mpc
    if( fstrSOLID%n_fix_mpc>0 ) then
      allocate( fstrSOLID%mpc_const( fstrSOLID%n_fix_mpc ) )
      fstrSOLID%mpc_const(:) = hecMESH%mpc%mpc_const(:)
    endif

    ! initialize for linear static problems
    fstrSOLID%FACTOR(2)=1.d0
    fstrSOLID%FACTOR(1)=0.d0
  end subroutine fstr_solid_alloc

  subroutine fstr_smoothed_element_init( hecMESH, fstrSOLID )
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_solid)                :: fstrSOLID

    logical, allocatable :: is_selem_list(:)
    integer :: i, isect

    do isect=1,hecMESH%section%n_sect
      if( fstrSOLID%sections(isect)%elemopt341 == kel341SESNS ) fstrSOLID%is_smoothing_active = .true.
    end do
    if( .not. fstrSOLID%is_smoothing_active ) return

    allocate(is_selem_list(hecMESH%n_elem), stat=i)
    if( i /= 0 ) then
      write(*,*) 'Allocation error: is_selem_list'
      return
    endif
    is_selem_list(:) = .false.

    do i=1,hecMESH%n_elem
      isect= hecMESH%section_ID(i)
      if( hecMESH%elem_type(i) /= fe_tet4n ) cycle
      if( fstrSOLID%sections(isect)%elemopt341 == kel341SESNS ) is_selem_list(i) = .true.
    enddo

    call hecmw_create_smoothing_element_connectivity(hecMESH,is_selem_list)

    deallocate(is_selem_list)

  end subroutine

  subroutine fstr_smoothed_element_calcmaxcon( hecMESH, fstrSOLID )
    use m_static_LIB_C3D4SESNS
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_solid)                :: fstrSOLID

    integer :: i, nodlocal(fstrSOLID%max_ncon), iiS, nn, con_stf

    if( fstrSOLID%max_ncon_stf > 20 ) fstrSOLID%max_ncon_stf = 20   

    do i=1,hecMESH%elem_type_index(hecMESH%n_elem_type)
      if( hecMESH%elem_type(i) /= fe_nodesmooth_tet4n .and. hecMESH%elem_type(i) /= fe_edgesmooth_tet4n ) cycle
      iiS = hecMESH%elem_node_index(i-1)
      nn = hecMESH%elem_node_index(i) - iiS
      nodlocal(1:nn) = hecMESH%elem_node_item(iiS+1:iiS+nn)
      con_stf = Return_nn_comp_C3D4_SESNS(nn, nodlocal)
      if( con_stf > fstrSOLID%max_ncon_stf ) fstrSOLID%max_ncon_stf = con_stf
    enddo
  end subroutine

  !> Initialize elements info in static calculation
  subroutine fstr_element_init( hecMESH, fstrSOLID, solution_type )
    use elementInfo
    use mMechGauss
    use m_fstr_FiniteRotationKinematics, only: fstr_uses_finite_rotation_kinematics
    use m_fstr
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_solid)                :: fstrSOLID
    integer(kind=kint), intent(in)  :: solution_type

    integer :: i, j, ng, isect, ndof, id, nn, n_elem, nthick

    if( hecMESH%n_elem <=0 ) then
      stop "no element defined!"
    endif

    fstrSOLID%maxn_gauss = 0
    fstrSOLID%max_ncon   = 0

    ! elemopt341 = kel341ES
    call fstr_smoothed_element_init( hecMESH, fstrSOLID )

    ! number of elements
    n_elem = hecMESH%elem_type_index(hecMESH%n_elem_type)
    allocate( fstrSOLID%elements(n_elem) )
    
    do i=1,n_elem
      fstrSOLID%elements(i)%elemact_flag = kELACT_UNDEFINED
      if( solution_type == kstHEAT) cycle !fstrSOLID is used only for elemact element in heat transfer analysis

      fstrSOLID%elements(i)%etype = hecMESH%elem_type(i)
      if( hecMESH%elem_type(i)==301 ) fstrSOLID%elements(i)%etype=111
      if (hecmw_is_etype_link(fstrSOLID%elements(i)%etype)) cycle
      if (hecmw_is_etype_patch(fstrSOLID%elements(i)%etype)) cycle
      ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
      if( ng > fstrSOLID%maxn_gauss ) fstrSOLID%maxn_gauss = ng
      if(ng>0) allocate( fstrSOLID%elements(i)%gausses( ng ) )

      isect= hecMESH%section_ID(i)
      ndof = getSpaceDimension( fstrSOLID%elements(i)%etype )
      if (ndof == 2) then              ! why do this???
        id=hecMESH%section%sect_opt(isect)
        if( id==0 ) then
          fstrSOLID%elements(i)%iset=1
        else if( id==1) then
          fstrSOLID%elements(i)%iset=0
        else if( id==2) then
          fstrSOLID%elements(i)%iset=2
        endif
      endif

      if( isect<0 .or. isect>hecMESH%section%n_sect )   &
        stop "Error in element's section definition"
      id = hecMESH%section%sect_mat_ID_item(isect)
      fstrSOLID%materials(id)%cdsys_ID = hecMESH%section%sect_orien_ID(isect)
      nn = hecMESH%elem_node_index(i)-hecMESH%elem_node_index(i-1)
      do j=1,ng
        fstrSOLID%elements(i)%gausses(j)%pMaterial => fstrSOLID%materials(id)
        call fstr_init_gauss( fstrSOLID%elements(i)%gausses( j )  )
      enddo
      nthick = 0
      if( fstr_uses_finite_rotation_kinematics( fstrSOLID%elements(i)%etype, nn, &
          fstrSOLID%materials(id) ) .or. &
          fstr_uses_plane_stress_j2_shell(fstrSOLID%elements(i)%etype, nn, &
          fstrSOLID%materials(id)) ) then
        nthick = fstr_shell_num_thickness_points( fstrSOLID%elements(i)%etype )
      endif
      if( nthick > 0 ) call fstr_init_shell_layer_gausses( fstrSOLID%elements(i), ng, &
        fstrSOLID%materials(id)%totallyr, nthick )

      allocate(fstrSOLID%elements(i)%equiForces(nn*ndof))
      fstrSOLID%elements(i)%equiForces = 0.0d0
      if( nn > fstrSOLID%max_ncon ) fstrSOLID%max_ncon = nn

      if( hecMESH%elem_type(i)==361 ) then
        if( fstrSOLID%sections(isect)%elemopt361==kel361IC ) then
          allocate( fstrSOLID%elements(i)%aux(3,3) )
          fstrSOLID%elements(i)%aux = 0.0d0
        endif
        fstrSOLID%elements(i)%p(1) = 0.0d0
      endif

    enddo

    fstrSOLID%max_ncon_stf = fstrSOLID%max_ncon
    if( fstrSOLID%is_smoothing_active ) call fstr_smoothed_element_calcmaxcon( hecMESH, fstrSOLID )

    call hecmw_allreduce_I1(hecMESH,fstrSOLID%maxn_gauss,HECMW_MAX)
  end subroutine

  recursive subroutine fstr_solid_physic_val_finalize( phys )
    use m_fstr
    implicit none
    type(fstr_solid_physic_val) :: phys
    integer :: i

    if(associated(phys%STRESS)) deallocate(phys%STRESS)
    if(associated(phys%STRAIN)) deallocate(phys%STRAIN)
    if(associated(phys%MISES)) deallocate(phys%MISES)

    if(associated(phys%PSTRESS)) deallocate(phys%PSTRESS)
    if(associated(phys%PSTRAIN)) deallocate(phys%PSTRAIN)
    if(associated(phys%PSTRESS_VECT)) deallocate(phys%PSTRESS_VECT)
    if(associated(phys%PSTRAIN_VECT)) deallocate(phys%PSTRAIN_VECT)

    if(associated(phys%ESTRESS)) deallocate(phys%ESTRESS)
    if(associated(phys%ESTRAIN)) deallocate(phys%ESTRAIN)
    if(associated(phys%EMISES)) deallocate(phys%EMISES)
    if(associated(phys%EPLSTRAIN)) deallocate(phys%EPLSTRAIN)

    if(associated(phys%EPSTRESS)) deallocate(phys%EPSTRESS)
    if(associated(phys%EPSTRAIN)) deallocate(phys%EPSTRAIN)
    if(associated(phys%EPSTRESS_VECT)) deallocate(phys%EPSTRESS_VECT)
    if(associated(phys%EPSTRAIN_VECT)) deallocate(phys%EPSTRAIN_VECT)
    if(associated(phys%ENQM)) deallocate(phys%ENQM)

    if(associated(phys%LAYER)) then
      do i=1, size(phys%LAYER)
        call fstr_solid_physic_val_finalize(phys%LAYER(i))
      end do
      deallocate(phys%LAYER)
    end if
    if(associated(phys%PLUS)) then
      call fstr_solid_physic_val_finalize(phys%PLUS)
      deallocate(phys%PLUS)
    end if
    if(associated(phys%MINUS)) then
      call fstr_solid_physic_val_finalize(phys%MINUS)
      deallocate(phys%MINUS)
    end if
  end subroutine

  !> Finalizer of fstr_solid
  subroutine fstr_solid_finalize( fstrSOLID )
    use mMechGauss, only: fstr_finalize_shell_layer_gausses
    implicit none
    type(fstr_solid) :: fstrSOLID
    integer :: i, j, ierror
    if( associated(fstrSOLID%materials) ) then
      do j=1,size(fstrSOLID%materials)
        call finalizeMaterial(fstrSOLID%materials(j))
      enddo
      deallocate( fstrSOLID%materials )
    endif
    if( .not. associated(fstrSOLID%elements ) ) return
    do i=1,size(fstrSOLID%elements)
      if( associated(fstrSOLID%elements(i)%gausses) ) then
        do j=1,size(fstrSOLID%elements(i)%gausses)
          call fstr_finalize_gauss(fstrSOLID%elements(i)%gausses(j))
        enddo
        deallocate( fstrSOLID%elements(i)%gausses )
      endif
      call fstr_finalize_shell_layer_gausses( fstrSOLID%elements(i) )
      if(associated(fstrSOLID%elements(i)%equiForces) ) then
        deallocate(fstrSOLID%elements(i)%equiForces)
      endif
      if( associated(fstrSOLID%elements(i)%aux) ) then
        deallocate(fstrSOLID%elements(i)%aux)
      endif
    enddo

    deallocate( fstrSOLID%elements )
    if( associated( fstrSOLID%mpc_const ) ) then
      deallocate( fstrSOLID%mpc_const )
    endif
    call free_stepInfo( fstrSOLID%step_ctrl_restart )
    if( associated(fstrSOLID%step_ctrl) ) then
      do i=1,size(fstrSOLID%step_ctrl)
        call free_stepInfo( fstrSOLID%step_ctrl(i) )
      enddo
      deallocate( fstrSOLID%step_ctrl )
    endif
    if(associated(fstrSOLID%output_ctrl) ) then
      do i=1,size(fstrSOLID%output_ctrl)
        if( fstrSOLID%output_ctrl(i)%filenum==IUTB )  &
          close(fstrSOLID%output_ctrl(i)%filenum)
      enddo
      deallocate(fstrSOLID%output_ctrl)
    endif
    if( associated( fstrSOLID%sections ) ) then
      deallocate( fstrSOLID%sections )
    endif

    if( associated(fstrSOLID%GL) ) then
      deallocate(fstrSOLID%GL               ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, GL>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%EFORCE) ) then
      deallocate(fstrSOLID%EFORCE           ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, EFORCE>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%unode) ) then
      deallocate(fstrSOLID%unode       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, unode>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%unode_bak) ) then
      deallocate(fstrSOLID%unode_bak   ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, unode_bak>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%dunode) ) then
      deallocate(fstrSOLID%dunode       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, dunode>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%ddunode) ) then
      deallocate(fstrSOLID%ddunode       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, ddunode>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_node_mode) ) then
      deallocate(fstrSOLID%shell_node_mode       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_node_mode>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_rot_state) ) then
      deallocate(fstrSOLID%shell_rot_state       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_rot_state>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_ref_triad) ) then
      deallocate(fstrSOLID%shell_ref_triad       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_ref_triad>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_triad) ) then
      deallocate(fstrSOLID%shell_triad       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_triad>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_triad_bak) ) then
      deallocate(fstrSOLID%shell_triad_bak       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_triad_bak>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_dtriad) ) then
      deallocate(fstrSOLID%shell_dtriad       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_dtriad>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_drill) ) then
      deallocate(fstrSOLID%shell_drill       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_drill>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_drill_bak) ) then
      deallocate(fstrSOLID%shell_drill_bak       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_drill_bak>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%shell_ddrill) ) then
      deallocate(fstrSOLID%shell_ddrill       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, shell_ddrill>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%QFORCE) ) then
      deallocate(fstrSOLID%QFORCE           ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, QFORCE>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%DFORCE) ) then
      deallocate(fstrSOLID%DFORCE           ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, DFORCE>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%temperature) ) then
      deallocate(fstrSOLID%temperature       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, temperature>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%last_temp) ) then
      deallocate(fstrSOLID%last_temp       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, reftemp>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrSOLID%temp_bak) ) then
      deallocate(fstrSOLID%temp_bak       ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, reftemp>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif

    ! Allocated in fstr_solid_append_boundary
    if( associated(fstrSOLID%BOUNDARY_ngrp_GRPID) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_GRPID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_GRPID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_ID) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_ID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_ID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_type) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_type, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_type>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_val) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_val, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_val>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_amp) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_amp, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_amp>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_istot) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_istot, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_istot>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_rotID) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_rotID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_rotID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%BOUNDARY_ngrp_centerID) ) then
       deallocate(fstrSOLID%BOUNDARY_ngrp_centerID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, BOUNDARY_ngrp_centerID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif

    if(associated(fstrSOLID%VELOCITY_ngrp_GRPID)) deallocate(fstrSOLID%VELOCITY_ngrp_GRPID)
    if(associated(fstrSOLID%VELOCITY_ngrp_ID)) deallocate(fstrSOLID%VELOCITY_ngrp_ID)
    if(associated(fstrSOLID%VELOCITY_ngrp_type)) deallocate(fstrSOLID%VELOCITY_ngrp_type)
    if(associated(fstrSOLID%VELOCITY_ngrp_amp)) deallocate(fstrSOLID%VELOCITY_ngrp_amp)
    if(associated(fstrSOLID%VELOCITY_ngrp_val)) deallocate(fstrSOLID%VELOCITY_ngrp_val)
    if(associated(fstrSOLID%VELOCITY_ngrp_rotID)) deallocate(fstrSOLID%VELOCITY_ngrp_rotID)
    if(associated(fstrSOLID%VELOCITY_ngrp_centerID)) deallocate(fstrSOLID%VELOCITY_ngrp_centerID)

    if(associated(fstrSOLID%ACCELERATION_ngrp_GRPID)) deallocate(fstrSOLID%ACCELERATION_ngrp_GRPID)
    if(associated(fstrSOLID%ACCELERATION_ngrp_ID)) deallocate(fstrSOLID%ACCELERATION_ngrp_ID)
    if(associated(fstrSOLID%ACCELERATION_ngrp_type)) deallocate(fstrSOLID%ACCELERATION_ngrp_type)
    if(associated(fstrSOLID%ACCELERATION_ngrp_amp)) deallocate(fstrSOLID%ACCELERATION_ngrp_amp)
    if(associated(fstrSOLID%ACCELERATION_ngrp_val)) deallocate(fstrSOLID%ACCELERATION_ngrp_val)

    ! Allocated in fstr_solid_append_cload
    if( associated(fstrSOLID%CLOAD_ngrp_GRPID) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_GRPID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_GRPID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_ID) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_ID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_ID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_DOF) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_DOF, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_DOF>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_val) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_val, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_val>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_amp) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_amp, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_amp>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_rotID) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_rotID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_rotID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif
    if( associated(fstrSOLID%CLOAD_ngrp_centerID) ) then
       deallocate(fstrSOLID%CLOAD_ngrp_centerID, stat=ierror)
       if( ierror /= 0 ) then
          write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, CLOAD_ngrp_centerID>'
          call flush(idbg)
          call hecmw_abort( hecmw_comm_get_comm())
       end if
    endif

    if(associated(fstrSOLID%DLOAD_ngrp_GRPID)) deallocate(fstrSOLID%DLOAD_ngrp_GRPID)
    if(associated(fstrSOLID%DLOAD_ngrp_ID)) deallocate(fstrSOLID%DLOAD_ngrp_ID)
    if(associated(fstrSOLID%DLOAD_ngrp_LID)) deallocate(fstrSOLID%DLOAD_ngrp_LID)
    if(associated(fstrSOLID%DLOAD_ngrp_amp)) deallocate(fstrSOLID%DLOAD_ngrp_amp)
    if(associated(fstrSOLID%DLOAD_ngrp_params)) deallocate(fstrSOLID%DLOAD_ngrp_params)

    if(associated(fstrSOLID%TEMP_ngrp_GRPID)) deallocate(fstrSOLID%TEMP_ngrp_GRPID)
    if(associated(fstrSOLID%TEMP_ngrp_ID)) deallocate(fstrSOLID%TEMP_ngrp_ID)
    if(associated(fstrSOLID%TEMP_ngrp_val)) deallocate(fstrSOLID%TEMP_ngrp_val)

    if(associated(fstrSOLID%SPRING_ngrp_GRPID)) deallocate(fstrSOLID%SPRING_ngrp_GRPID)
    if(associated(fstrSOLID%SPRING_ngrp_ID)) deallocate(fstrSOLID%SPRING_ngrp_ID)
    if(associated(fstrSOLID%SPRING_ngrp_DOF)) deallocate(fstrSOLID%SPRING_ngrp_DOF)
    if(associated(fstrSOLID%SPRING_ngrp_amp)) deallocate(fstrSOLID%SPRING_ngrp_amp)
    if(associated(fstrSOLID%SPRING_ngrp_val)) deallocate(fstrSOLID%SPRING_ngrp_val)
    if(associated(fstrSOLID%SPRING_incremental)) deallocate(fstrSOLID%SPRING_incremental)

    if(associated(fstrSOLID%COUPLE_ngrp_ID)) deallocate(fstrSOLID%COUPLE_ngrp_ID)

!    if(associated(fstrSOLID%STRESS)) deallocate(fstrSOLID%STRESS)
!    if(associated(fstrSOLID%STRAIN)) deallocate(fstrSOLID%STRAIN)
!    if(associated(fstrSOLID%MISES)) deallocate(fstrSOLID%MISES)
!    if(associated(fstrSOLID%PSTRESS)) deallocate(fstrSOLID%PSTRESS)
!    if(associated(fstrSOLID%PSTRAIN)) deallocate(fstrSOLID%PSTRAIN)
!    if(associated(fstrSOLID%PSTRESS_VECT)) deallocate(fstrSOLID%PSTRESS_VECT)
!    if(associated(fstrSOLID%PSTRAIN_VECT)) deallocate(fstrSOLID%PSTRAIN_VECT)

!    if(associated(fstrSOLID%ESTRESS)) deallocate(fstrSOLID%ESTRESS)
!    if(associated(fstrSOLID%ESTRAIN)) deallocate(fstrSOLID%ESTRAIN)
!    if(associated(fstrSOLID%EMISES)) deallocate(fstrSOLID%EMISES)
!    if(associated(fstrSOLID%EPLSTRAIN)) deallocate(fstrSOLID%EPLSTRAIN)

!    if(associated(fstrSOLID%EPSTRESS)) deallocate(fstrSOLID%EPLSTRAIN)
!    if(associated(fstrSOLID%EPSTRAIN)) deallocate(fstrSOLID%EPLSTRAIN)
!    if(associated(fstrSOLID%EPSTRESS_VECT)) deallocate(fstrSOLID%EPLSTRAIN)
!    if(associated(fstrSOLID%EPSTRAIN_VECT)) deallocate(fstrSOLID%EPLSTRAIN)

    if(associated(fstrSOLID%TNSTRAIN)) deallocate(fstrSOLID%TNSTRAIN)
    if(associated(fstrSOLID%TESTRAIN)) deallocate(fstrSOLID%TESTRAIN)

    if(associated(fstrSOLID%YIELD_RATIO)) deallocate(fstrSOLID%YIELD_RATIO)

!    if(associated(fstrSOLID%ENQM)) deallocate(fstrSOLID%ENQM)
    if(associated(fstrSOLID%REACTION)) deallocate(fstrSOLID%REACTION)

    if(associated(fstrSOLID%CONT_NFORCE)) deallocate(fstrSOLID%CONT_NFORCE)
    if(associated(fstrSOLID%CONT_FRIC)) deallocate(fstrSOLID%CONT_FRIC)
    if(associated(fstrSOLID%CONT_RELVEL)) deallocate(fstrSOLID%CONT_RELVEL)
    if(associated(fstrSOLID%CONT_STATE)) deallocate(fstrSOLID%CONT_STATE)
    if(associated(fstrSOLID%CONT_SGRP_ID)) deallocate(fstrSOLID%CONT_SGRP_ID)
    if(associated(fstrSOLID%CONT_AREA)) deallocate(fstrSOLID%CONT_AREA)
    if(associated(fstrSOLID%CONT_NTRAC)) deallocate(fstrSOLID%CONT_NTRAC)
    if(associated(fstrSOLID%CONT_FTRAC)) deallocate(fstrSOLID%CONT_FTRAC)
    if(associated(fstrSOLID%EMBED_NFORCE)) deallocate(fstrSOLID%EMBED_NFORCE)

    if(associated(fstrSOLID%SOLID)) call fstr_solid_physic_val_finalize(fstrSOLID%SOLID)
    if(associated(fstrSOLID%SHELL)) call fstr_solid_physic_val_finalize(fstrSOLID%SHELL)
    if(associated(fstrSOLID%BEAM)) call fstr_solid_physic_val_finalize(fstrSOLID%BEAM)

    if(associated(fstrSOLID%is_rot)) deallocate(fstrSOLID%is_rot)
    if(associated(fstrSOLID%GL0)) deallocate(fstrSOLID%GL0)
    if(associated(fstrSOLID%QFORCE_bak)) deallocate(fstrSOLID%QFORCE_bak)

    if(associated(fstrSOLID%contacts)) then
!      do i=1, size(fstrSOLID%contacts)
!        call fstr_contact_finalize(fstrSOLID%contacts(i))
!      end do
      deallocate(fstrSOLID%contacts)
    end if
    if(associated(fstrSOLID%embeds)) then
!      do i=1, size(fstrSOLID%embeds)
!        call fstr_contact_finalize(fstrSOLID%embeds(i))
!      end do
      deallocate(fstrSOLID%embeds)
    end if
    if(associated(fstrSOLID%unode_bkup)) deallocate(fstrSOLID%unode_bkup)
    if(associated(fstrSOLID%QFORCE_bkup)) deallocate(fstrSOLID%QFORCE_bkup)
    if(associated(fstrSOLID%DFORCE_bkup)) deallocate(fstrSOLID%DFORCE_bkup)
    if(associated(fstrSOLID%last_temp_bkup)) deallocate(fstrSOLID%last_temp_bkup)
    if(associated(fstrSOLID%shell_triad_bkup)) deallocate(fstrSOLID%shell_triad_bkup)
    if(associated(fstrSOLID%shell_drill_bkup)) deallocate(fstrSOLID%shell_drill_bkup)
    if(associated(fstrSOLID%elements_bkup)) then
      do i=1,size(fstrSOLID%elements_bkup)
        if( associated(fstrSOLID%elements_bkup(i)%gausses) ) then
          do j=1,size(fstrSOLID%elements_bkup(i)%gausses)
            call fstr_finalize_gauss(fstrSOLID%elements_bkup(i)%gausses(j))
          enddo
          deallocate( fstrSOLID%elements_bkup(i)%gausses )
        endif
        call fstr_finalize_shell_layer_gausses( fstrSOLID%elements_bkup(i) )
        if(associated(fstrSOLID%elements_bkup(i)%equiForces) ) then
          deallocate(fstrSOLID%elements_bkup(i)%equiForces)
        endif
        if( associated(fstrSOLID%elements_bkup(i)%aux) ) then
          deallocate(fstrSOLID%elements_bkup(i)%aux)
        endif
      enddo
      deallocate(fstrSOLID%elements_bkup)
    end if
    if(associated(fstrSOLID%contacts_bkup)) then
!      do i=1, size(fstrSOLID%contacts_bkup)
!        call fstr_contact_finalize(fstrSOLID%contacts_bkup(i))
!      end do
      deallocate(fstrSOLID%contacts_bkup)
    end if
    if(associated(fstrSOLID%embeds_bkup)) then
!      do i=1, size(fstrSOLID%embeds_bkup)
!        call fstr_contact_finalize(fstrSOLID%embeds_bkup(i))
!      end do
      deallocate(fstrSOLID%embeds_bkup)
    end if

  end subroutine

  !-----------------------------------------------------------------------------!
  !> Read in !STEP                                                             !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_STEP( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P
    character(HECMW_NAME_LEN) :: amp
    integer(kind=kint) :: amp_id

    integer(kind=kint) :: rcode, iproc

    amp = ' '
    rcode = fstr_ctrl_get_STEP( ctrl, amp, iproc )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    call amp_name_to_id( P%MESH, '!STEP', amp, amp_id )
    !    P%SOLID%NLSTATIC_ngrp_amp = amp_id;

  end subroutine fstr_setup_STEP

  !> Fix the convergence criteria of a step from !STEP and the !CONVERG_PARAM it refers to
  subroutine fstr_setup_step_convergence( stepinfo, cpname, cnvparam, istep )
    implicit none
    type( step_info ), intent(inout)  :: stepinfo       !< step info
    character(len=*), intent(in)      :: cpname         !< CONVERGPARAM of !STEP
    type( tParamConverg ), intent(in) :: cnvparam(0:)   !< !CONVERG_PARAM read so far; 0 gives none
    integer(kind=kint), intent(in)    :: istep          !< step number

    integer(kind=kint) :: i, id

    id = 0
    if( len_trim(cpname) > 0 ) then
      do i=1,ubound(cnvparam,1)
        if( hecmw_streqr( cnvparam(i)%name, cpname ) ) then
          id = i; exit
        endif
      enddo
      if( id == 0 ) then
        write(*,*) '### Error: CONVERG_PARAM referred to by STEP is not defined before it : ', trim(cpname)
        write(ILOG,*) '### Error: CONVERG_PARAM referred to by STEP is not defined before it : ', trim(cpname)
        stop HECMW_EXIT_INPUT
      endif
    endif

    call setup_stepInfo_converg( stepinfo, cnvparam(id) )

    if( .not. any( stepinfo%cnv_check(kcnvResidual:kcnvCorrection, kcnvTranslation:kcnvRotation, :) ) ) then
      write(*,*) '### Error: RESIDUAL and CORRECTION are all SKIP in step : ', istep
      write(ILOG,*) '### Error: RESIDUAL and CORRECTION are all SKIP in step : ', istep
      stop HECMW_EXIT_INPUT
    endif

  end subroutine fstr_setup_step_convergence

  !> Links step_ctrl(istep), whose step is already read, to !TIME_POINTS, !AUTOINC_PARAM and
  !> !CONVERG_PARAM given by name, and keeps the !DYNAMIC time increment for a fixed increment.
  subroutine fstr_setup_step_links( fstrPARAM, fstrSOLID, fstrDYNAMIC, istep, tp_name, ainc_name, cnv_name )
    implicit none
    type(fstr_param)               :: fstrPARAM
    type(fstr_solid)               :: fstrSOLID
    type(fstr_dynamic)             :: fstrDYNAMIC
    integer(kind=kint), intent(in) :: istep          !< step number
    character(len=*), intent(in)   :: tp_name        !< name of !TIME_POINTS
    character(len=*), intent(in)   :: ainc_name      !< name of !AUTOINC_PARAM
    character(len=*), intent(in)   :: cnv_name       !< name of !CONVERG_PARAM

    integer(kind=kint) :: i

    call fstr_setup_step_convergence( fstrSOLID%step_ctrl(istep), cnv_name, fstrPARAM%cnvparam, istep )
    ! For DYNAMIC fixed-increment: keep the !DYNAMIC time increment while preserving !STEP duration.
    ! fstr_ctrl_get_ISTEP unconditionally sets initdt=1/num_substep which is wrong for DYNAMIC.
    ! Only override initdt/mindt/maxdt; keep elapsetime and num_substep as-is.
    if( fstrPARAM%solution_type==kstDYNAMIC .and. &
      & fstrSOLID%step_ctrl(istep)%inc_type == stepFixedInc ) then
      fstrSOLID%step_ctrl(istep)%initdt      = fstrDYNAMIC%t_delta
      fstrSOLID%step_ctrl(istep)%mindt       = fstrDYNAMIC%t_delta
      fstrSOLID%step_ctrl(istep)%maxdt       = fstrDYNAMIC%t_delta
    endif
    if( associated(fstrPARAM%timepoints) ) then
      do i=1,size(fstrPARAM%timepoints)
        if( hecmw_streqr( fstrPARAM%timepoints(i)%name, tp_name ) ) then
          fstrSOLID%step_ctrl(istep)%timepoint_id = i; exit
        endif
      enddo
    endif
    if( associated(fstrPARAM%ainc) ) then
      do i=1,size(fstrPARAM%ainc)-1
        if( hecmw_streqr( fstrPARAM%ainc(i)%name, ainc_name ) ) then
          fstrSOLID%step_ctrl(istep)%AincParam_id = i; exit
        endif
      enddo
    endif
  end subroutine fstr_setup_step_links

  !> Read in !WRITE after fstr_setup_alloc and store the output controls.
  !> fg_result and fg_visual are set by the caller before fstr_setup_alloc.
  subroutine fstr_setup_WRITE_store( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    type(t_output_ctrl) :: outctrl
    integer(kind=kint) :: islog, result, visual, femap, c_output, ierror
    character(len=HECMW_FILENAME_LEN) :: logfileNAME

    fstrSOLID => P%SOLID
    write( logfileNAME, '(i5,''.log'')' ) myrank

    call fstr_ctrl_get_output( ctrl, outctrl, islog, result, visual, femap )
    if( islog == 1 ) then
      c_output=1
      outctrl%filename = trim(logfileNAME)
      outctrl%filenum = ILOG
      call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
    endif
    if( femap == 1 ) then
      c_output=2
      write( outctrl%filename, '(a,i0,a)') 'utable.',myrank,'.dat'
      outctrl%filenum = IUTB
      call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
      open( unit=outctrl%filenum, file=outctrl%filename, status='REPLACE', iostat=ierror )
      if( ierror /= 0 ) then
        write(*,*) 'Warning: cannot open output file: ', trim(outctrl%filename)
      endif
      fstrSOLID%output_ctrl(c_output)%outinfo%grp_id = 1
    endif
    if( result == 1 ) then
      c_output=3
      call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
    endif
    if( visual == 1 ) then
      c_output=4
      call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
    endif
  end subroutine fstr_setup_WRITE_store

  !-----------------------------------------------------------------------------!
  !> Read in !COUPLE                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_COUPLE( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P
    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint) :: n

    if( P%SOLID%file_type /= kbcfFSTR ) return

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    rcode = fstr_ctrl_get_COUPLE( ctrl,           &
      P%PARAM%fg_couple_type,       &
      P%PARAM%fg_couple_first,      &
      P%PARAM%fg_couple_window,     &
      grp_id_name, HECMW_NAME_LEN  )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_solid_append_couple( P%MESH, P%SOLID, P%PARAM, '!COUPLE', n, grp_id_name )

    deallocate( grp_id_name )

  end subroutine fstr_setup_COUPLE

  !> Read in !ELEMENT_ACTIVATION
  subroutine fstr_setup_ELEMENT_ACTIVATION( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid, mode, measure, state
    real(kind=kreal)   :: eps
    real(kind=kreal), pointer :: thlow(:), thup(:)

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n), thlow(n), thup(n) )
    amp = ' '
    eps = 1.d-3
    rcode = fstr_ctrl_get_ELEMENT_ACTIVATION( ctrl, amp, eps, grp_id_name, mode, measure, state, thlow, thup )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_solid_append_element_activation( P%MESH, P%SOLID, '!ELEMENT_ACTIVATION', gid, &
      amp, eps, measure, state, n, grp_id_name, thlow, thup )

    deallocate( grp_id_name )
  end subroutine fstr_setup_ELEMENT_ACTIVATION

  !-----------------------------------------------------------------------------!
  !> Read in !STATIC(old)                                                           !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_STATIC( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P
    integer(kind=kint) :: rcode

    integer :: nout, nout_monit,node_monit_1 ,elem_monit_1 ,intg_monit_1
    integer :: ipt, idx_elpl, iout_list(6)
    real(kind=kreal) :: sig_y0, h_dash
    real(kind=kreal) :: DT, ETIME, EPS ! dummy
    integer :: ITMAX ! dummy

    if( counter > 1 ) then
      write(*,*)
    endif

    ipt = 0
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ', 'INFINITESIMAL,NLGEOM,INFINITE ', 0, 'P', ipt  )/=0 )  &
      return
    if( ipt == 2 ) P%PARAM%nlgeom = .true.

    ! for backward compatibility
    if( ipt == 3 ) then
      write(*,*) "Warning : !STATIC : parameter 'TYPE=INFINITE' is deprecated." &
           & //  " Please use the replacement parameter 'TYPE=INFINITESIMAL'"
    endif

    rcode = fstr_ctrl_get_STATIC( ctrl, &
      DT, ETIME, ITMAX, EPS, P%SOLID%restart_nout, &
      idx_elpl, &
      iout_list, &
      sig_y0, h_dash, &
      nout, nout_monit, node_monit_1, &
      elem_monit_1, intg_monit_1 )

    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_STATIC

  !-----------------------------------------------------------------------------!
  !> Read in !BOUNDARY                                                    !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_BOUNDARY( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp, rotc_name(1)
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: dof_ids (:)
    integer(kind=kint),pointer :: dof_ide (:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint) :: n

    integer(kind=kint) :: gid, istot

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

    istot = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TOTAL ', '# ', 0, 'E', istot )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    ! get center of torque load
    rotc_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'ROT_CENTER ', '# ', 0, 'S', rotc_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    ! ENTIRE -----------------------------------------------
    P%SOLID%file_type = kbcfFSTR

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_id_name(n) )
    allocate( dof_ids (n) )
    allocate( dof_ide (n) )
    allocate( val_ptr(n) )

    amp = ' '
    val_ptr = 0.0d0
    if( n > 0 ) then
      rcode = fstr_ctrl_get_BOUNDARY( ctrl, amp, grp_id_name, HECMW_NAME_LEN, dof_ids, dof_ide, val_ptr)
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif
    call fstr_solid_append_boundary( P%MESH, P%SOLID, '!BOUNDARY', gid, istot, rotc_name(1), &
      amp, n, grp_id_name, dof_ids, dof_ide, val_ptr )

    deallocate( grp_id_name )
    deallocate( dof_ids )
    deallocate( dof_ide )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( dof_ids )
    nullify( dof_ide )
    nullify( val_ptr )

  end subroutine fstr_setup_BOUNDARY

  !-----------------------------------------------------------------------------!
  !> Read in !CLOAD                                                           !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_CLOAD( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp, rotc_name(1)
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint),pointer :: id_ptr(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid

    if( P%SOLID%file_type /= kbcfFSTR ) return
    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    ! get center of torque load
    rotc_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'ROT_CENTER ', '# ', 0, 'S', rotc_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_id_name(n))
    allocate( id_ptr(n) )
    allocate( val_ptr(n) )
    amp = ' '
    id_ptr = 0
    val_ptr = 0.0d0
    if( n > 0 ) then
      rcode = fstr_ctrl_get_CLOAD( ctrl, amp, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif
    call fstr_solid_append_cload( P%MESH, P%SOLID, '!CLOAD', gid, rotc_name(1), &
      amp, n, grp_id_name, id_ptr, val_ptr )

    deallocate( grp_id_name )
    deallocate( id_ptr )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( id_ptr )
    nullify( val_ptr )

  end subroutine fstr_setup_CLOAD

  !> Read in !DLOAD
  subroutine fstr_setup_DLOAD( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    integer(kind=kint) :: follow
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: new_params(:,:)
    integer(kind=kint),pointer :: lid_ptr(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid

    if( P%SOLID%file_type /= kbcfFSTR ) return

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( lid_ptr(n) )
    allocate( new_params(0:6,n))
    new_params = 0
    amp = ' '
    follow = P%SOLID%DLOAD_follow
    if( .not. P%PARAM%nlgeom ) follow = 0
    rcode = fstr_ctrl_get_DLOAD( ctrl, amp, follow, &
      grp_id_name, HECMW_NAME_LEN,    &
      lid_ptr, new_params )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    call fstr_solid_append_dload( P%MESH, P%SOLID, '!DLOAD', gid, follow, &
      amp, n, grp_id_name, lid_ptr, new_params )
    deallocate( grp_id_name )
    deallocate( lid_ptr )
    deallocate( new_params )
    nullify( grp_id_name )
    nullify( lid_ptr )
    nullify( new_params )
  end subroutine fstr_setup_DLOAD

  !-----------------------------------------------------------------------------!
  !> Read in !TEMPERATURE                                                  !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_TEMPERATURE( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, gid
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint) :: n

    if( P%SOLID%file_type /= kbcfFSTR ) return

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_id_name(n))
    allocate( val_ptr(n) )
    val_ptr = 0.0d0

    rcode = fstr_ctrl_get_TEMPERATURE( ctrl,      &
      P%SOLID%TEMP_irres,           &
      P%SOLID%TEMP_tstep,           &
      P%SOLID%TEMP_interval,        &
      P%SOLID%TEMP_rtype,           &
      grp_id_name, HECMW_NAME_LEN,  &
      val_ptr )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( P%SOLID%TEMP_irres == 0 ) then
      call fstr_solid_append_temperature( P%MESH, P%SOLID, '!TEMPERATURE', gid, n, grp_id_name, val_ptr )
    endif

    deallocate( grp_id_name )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( val_ptr )

  end subroutine fstr_setup_TEMPERATURE

  !-----------------------------------------------------------------------------!
  !> Read in !SPRING                                                            !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_SPRING( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, incremental
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint),pointer :: id_ptr(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid

    if( P%SOLID%file_type /= kbcfFSTR ) return
    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )
    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( id_ptr(n) )
    allocate( val_ptr(n) )
    amp = ' '
    incremental = 0
    id_ptr = 0
    val_ptr = 0.0d0
    rcode = fstr_ctrl_get_SPRING( ctrl, amp, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr, incremental )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    call fstr_solid_append_spring( P%MESH, P%SOLID, '!SPRING', gid, incremental, &
      amp, n, grp_id_name, id_ptr, val_ptr )

    deallocate( grp_id_name )
    deallocate( id_ptr )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( id_ptr )
    nullify( val_ptr )

  end subroutine fstr_setup_SPRING

  !-----------------------------------------------------------------------------!
  !> Read in !VELOCITY                                                         !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_VELOCITY( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    integer(kind=kint) :: vType
    character(HECMW_NAME_LEN) :: amp, rotc_name(1)
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: dof_ids (:)
    integer(kind=kint),pointer :: dof_ide (:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',  0, 'I', gid  )

    ! get center of rotation ( same syntax as !BOUNDARY,ROT_CENTER= )
    rotc_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'ROT_CENTER ', '# ', 0, 'S', rotc_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_id_name(n))
    allocate( dof_ids (n))
    allocate( dof_ide (n))
    allocate( val_ptr(n) )

    vType = kbcTransit
    amp = ''
    val_ptr = 0.0d0
    if( n > 0 ) then
      rcode = fstr_ctrl_get_VELOCITY( ctrl,  vType, amp,   &
        grp_id_name, HECMW_NAME_LEN,  &
        dof_ids, dof_ide, val_ptr )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif
    call fstr_solid_append_velocity( P%MESH, P%SOLID, P%DYN, '!VELOCITY', gid, vType, rotc_name(1), &
      amp, n, grp_id_name, dof_ids, dof_ide, val_ptr )

    deallocate( grp_id_name )
    deallocate( dof_ids )
    deallocate( dof_ide )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( dof_ids )
    nullify( dof_ide )
    nullify( val_ptr )

  end subroutine fstr_setup_VELOCITY

  !-----------------------------------------------------------------------------!
  !> Read in !ACCELERATION                                                     !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_ACCELERATION( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    integer(kind=kint) :: aType
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: dof_ids (:)
    integer(kind=kint),pointer :: dof_ide (:)
    real(kind=kreal),pointer :: val_ptr(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: gid

    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',  0, 'I', gid  )

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( dof_ids (n))
    allocate( dof_ide (n))
    allocate( val_ptr(n))

    amp = ' '
    val_ptr = 0.0d0
    rcode = fstr_ctrl_get_ACCELERATION( ctrl,  aType, amp,   &
      grp_id_name, HECMW_NAME_LEN,  &
      dof_ids, dof_ide,  val_ptr)
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    call fstr_solid_append_acceleration( P%MESH, P%SOLID, P%DYN, '!ACCELERATION', gid, aType, &
      amp, n, grp_id_name, dof_ids, dof_ide, val_ptr )

    deallocate( grp_id_name )
    deallocate( dof_ids )
    deallocate( dof_ide )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( dof_ids )
    nullify( dof_ide )
    nullify( val_ptr )
  end subroutine fstr_setup_ACCELERATION

  subroutine fstr_setup_solid_nastran( ctrl, hecMESH, fstrSOLID )
    implicit none
    integer(kind=kint) :: ctrl
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid        ) :: fstrSOLID
    write(ILOG,*) '### Error : In !BOUNDARY, TYPE=NASTRAN is not supported.'
    call fstr_abort( HECMW_EXIT_MODEL )
  end subroutine fstr_setup_solid_nastran

  !-----------------------------------------------------------------------------!
  !> Read in !CONTACT                                                           !
  !-----------------------------------------------------------------------------!
  !> c_contact is the number of the pairs already read, and the pairs of this card are added to it.
  subroutine fstr_setup_CONTACT( ctrl, c_contact, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_contact
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: n, k
    real(kind=kreal) :: ee, pp, rho, alpha
    character(len=HECMW_FILENAME_LEN) :: mName

    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( .not. fstr_ctrl_get_CONTACT( ctrl, n, fstrSOLID%contacts(c_contact+1:c_contact+n)   &
        ,ee, pp, rho, alpha, P%PARAM%contact_algo, mName, k ) ) then
      write(*,*) '### Error: Fail in read in contact condition : ', c_contact
      write(ILOG,*) '### Error: Fail in read in contact condition : ', c_contact
      stop HECMW_EXIT_INPUT
    endif
    ! initialize contact condition (cdotp and mut are now obsolete - penalty is in tContact)
    if( rho>0.d0 ) cgn = rho
    if( alpha>0.d0 ) cgt = alpha
    call fstr_solid_init_contacts( P%MESH, fstrSOLID, fstrPARAM, c_contact, n, mName, k )
    c_contact = c_contact+n
  end subroutine fstr_setup_CONTACT

  !-----------------------------------------------------------------------------!
  !> Read in !EMBED                                                             !
  !-----------------------------------------------------------------------------!
  !> c_embed is the number of the pairs already read, and the pairs of this card are added to it.
  subroutine fstr_setup_EMBED( ctrl, c_embed, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_embed
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: n, k
    character(len=HECMW_FILENAME_LEN) :: mName

    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( .not. fstr_ctrl_get_EMBED( ctrl, n, fstrSOLID%embeds(c_embed+1:c_embed+n), mName, k ) ) then
      write(*,*) '### Error: Fail in read in embed condition : ', c_embed
      write(ILOG,*) '### Error: Fail in read in embed condition : ', c_embed
      stop HECMW_EXIT_INPUT
    endif
    call fstr_solid_init_embeds( P%MESH, fstrSOLID, fstrPARAM, c_embed, n, mName, k )
    c_embed = c_embed+n
  end subroutine fstr_setup_EMBED

  !-----------------------------------------------------------------------------!
  !> Read in !ISTEP                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_ISTEP( ctrl, c_istep, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_istep
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: i
    character(len=HECMW_FILENAME_LEN) :: mName, mName2, mName3

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID

    if( .not. fstr_ctrl_get_ISTEP( ctrl, hecMESH, fstrSOLID%step_ctrl(c_istep), mName, mName2, mName3 ) ) then
      write(*,*) '### Error: Fail in read in step definition : ' , c_istep
      write(ILOG,*) '### Error: Fail in read in step definition : ', c_istep
      stop HECMW_EXIT_INPUT
    endif
    call fstr_setup_step_convergence( fstrSOLID%step_ctrl(c_istep), mName3, fstrPARAM%cnvparam, c_istep )
    if( associated(fstrPARAM%timepoints) ) then
      do i=1,size(fstrPARAM%timepoints)
        if( hecmw_streqr( fstrPARAM%timepoints(i)%name, mName ) ) then
          fstrSOLID%step_ctrl(c_istep)%timepoint_id = i; exit
        endif
      enddo
    endif
    if( associated(fstrPARAM%ainc) ) then
      do i=1,size(fstrPARAM%ainc)
        if( hecmw_streqr( fstrPARAM%ainc(i)%name, mName2 ) ) then
          fstrSOLID%step_ctrl(c_istep)%AincParam_id = i; exit
        endif
      enddo
    endif
  end subroutine fstr_setup_ISTEP

  !-----------------------------------------------------------------------------!
  !> Read in !STEP                                                              !
  !-----------------------------------------------------------------------------!
  !> Read in !STEP of !VERSION 1 or later after fstr_setup_alloc and store the step.
  subroutine fstr_setup_STEP_store( ctrl, c_istep, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_istep
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_dynamic), pointer :: fstrDYNAMIC
    character(len=HECMW_FILENAME_LEN) :: mName, mName2, mName3

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID
    fstrDYNAMIC => P%DYN

    if( .not. fstr_ctrl_get_ISTEP( ctrl, hecMESH, fstrSOLID%step_ctrl(c_istep), mName, mName2, mName3 ) ) then
      write(*,*) '### Error: Fail in read in step definition : ' , c_istep
      write(ILOG,*) '### Error: Fail in read in step definition : ', c_istep
      stop HECMW_EXIT_INPUT
    endif
    call fstr_setup_step_links( fstrPARAM, fstrSOLID, fstrDYNAMIC, c_istep, mName, mName2, mName3 )
  end subroutine fstr_setup_STEP_store

  !-----------------------------------------------------------------------------!
  !> Read in !SECTION                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_SECTION( ctrl, c_section, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_section
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_solid), pointer :: fstrSOLID

    hecMESH => P%MESH
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_SECTION( ctrl, hecMESH, fstrSOLID%sections )/=0 ) then
      write(*,*) '### Error: Fail in read in SECTION definition : ' , c_section
      write(ILOG,*) '### Error: Fail in read in SECTION definition : ', c_section
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_SECTION

  !-----------------------------------------------------------------------------!
  !> Read in !ELEMOPT                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_ELEMOPT( ctrl, c_elemopt, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_elemopt
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID

    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_ELEMOPT( ctrl, fstrSOLID%elemopt361 )/=0 ) then
      write(*,*) '### Error: Fail in read in ELEMOPT definition : ' , c_elemopt
      write(ILOG,*) '### Error: Fail in read in ELEMOPT definition : ', c_elemopt
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_ELEMOPT

  !-----------------------------------------------------------------------------!
  !> Read in !MATERIAL                                                          !
  !-----------------------------------------------------------------------------!
  !> P%STATE%cid is set to the material of this card, and P%STATE%cache to the material expected next.
  subroutine fstr_setup_MATERIAL( ctrl, c_material, P )
    use mMaterial
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: c_material
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: i, cid, cache
    character(len=HECMW_FILENAME_LEN) :: mName

    hecMESH => P%MESH
    fstrSOLID => P%SOLID
    cache = P%STATE%cache

    if( fstr_ctrl_get_MATERIAL( ctrl, mName )/=0 ) then
      write(*,*) '### Error: Fail in read in material definition : ' , c_material
      write(ILOG,*) '### Error: Fail in read in material definition : ', c_material
      stop HECMW_EXIT_INPUT
    endif
    cid = 0
    if(cache < hecMESH%material%n_mat) then
      if(hecmw_streqr( hecMESH%material%mat_name(cache), mName ))then
        cid = cache
        cache = cache + 1
      endif
    endif
    if(cid == 0)then
      do i=1,hecMESH%material%n_mat
        if( hecmw_streqr( hecMESH%material%mat_name(i), mName ) ) then
          cid = i
          cache = i + 1
          exit
        endif
      enddo
    endif
    if(cid == 0)then
      write(*,*) '### Error: Fail in read in material definition : ' , c_material
      write(ILOG,*) '### Error: Fail in read in material definition : ', c_material
      stop HECMW_EXIT_INPUT
    endif
    fstrSOLID%materials(cid)%name = hecMESH%material%mat_name(cid)
    if(c_material>hecMESH%material%n_mat) call initMaterial( fstrSOLID%materials(cid) )
    P%STATE%cid = cid
    P%STATE%cache = cache
  end subroutine fstr_setup_MATERIAL

  !> Material of the last !MATERIAL, for a material card. Stops at a material card given before the first !MATERIAL.
  integer(kind=kint) function fstr_setup_current_material( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: header_name

    if( P%STATE%cid == 0 ) then
      rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )
      call fstr_setup_util_err_stop( '### Error: '//trim(header_name)// &
        ' : given before the first !MATERIAL; put it after !MATERIAL, NAME=<material name>' )
    endif
    fstr_setup_current_material = P%STATE%cid
  end function fstr_setup_current_material

  !-----------------------------------------------------------------------------!
  !> Read in !ELASTIC                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_ELASTIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_ELASTICITY( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables,   &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in elasticity definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in elasticity definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_ELASTIC

  !-----------------------------------------------------------------------------!
  !> Read in !PLASTIC                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_PLASTIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_PLASTICITY( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables,   &
        fstrSOLID%materials(cid)%table,   &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_PLASTIC

  !-----------------------------------------------------------------------------!
  !> Read in !HYPERELASTIC                                                      !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_HYPERELASTIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_HYPERELASTIC( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables )/=0 ) then
      write(*,*) '### Error: Fail in read in elasticity definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in elasticity definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_HYPERELASTIC

  !-----------------------------------------------------------------------------!
  !> Read in !VISCOELASTIC                                                      !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_VISCOELASTIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_VISCOELASTICITY( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_VISCOELASTIC

  !-----------------------------------------------------------------------------!
  !> Read in !TRS                                                               !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_TRS( ctrl, P )
    use mMaterial
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstrSOLID%materials(cid)%mtype/=VISCOELASTIC ) then
      write(*,*) '### WARNING: TRS can only be defined for viscoelastic material! It is ignored! '
      write(ILOG,*) '### WARNING: TRS can only be defined for viscoelastic material! It is ignored! '
    else
      if( fstr_ctrl_get_TRS( ctrl, fstrSOLID%materials(cid)%mtype, fstrSOLID%materials(cid)%variables)/=0 ) then
        write(*,*) '### Error: Fail in read in TRS definition : ' , cid
        write(ILOG,*) '### Error: Fail in read in TRS definition : ', cid
        stop HECMW_EXIT_INPUT
      endif
    endif
  end subroutine fstr_setup_TRS

  !-----------------------------------------------------------------------------!
  !> Read in !CREEP                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_CREEP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_VISCOPLASTICITY( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_CREEP

  !-----------------------------------------------------------------------------!
  !> Read in !DENSITY                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_DENSITY( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_DENSITY( ctrl, fstrSOLID%materials(cid)%variables )/=0 ) then
      write(*,*) '### Error: Fail in read in density definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in density definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_DENSITY

  !-----------------------------------------------------------------------------!
  !> Read in !EXPANSION_COEF                                                    !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_EXPANSION_COEFF( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_EXPANSION_COEFF( ctrl, fstrSOLID%materials(cid)%variables, &
        fstrSOLID%materials(cid)%dict)/=0 )  then
      write(*,*) '### Error: Fail in read in expansion coefficient definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in expansion coefficient definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_EXPANSION_COEFF

  !-----------------------------------------------------------------------------!
  !> Read in !DAMPING                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_DAMPING( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_RAYLEIGH_DAMPING( ctrl, fstrSOLID%materials(cid)%variables, &
        fstrSOLID%materials(cid)%is_elem_Rayleigh_damping)/=0 )  then
      write(*,*) '### Error: Fail in read in damping definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in damping definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_DAMPING

  !-----------------------------------------------------------------------------!
  !> Read in !FLUID                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_FLUID( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_FLUID( ctrl,                                 &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables,   &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in fluid definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in fluid definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_FLUID

  !-----------------------------------------------------------------------------!
  !> Read in !SPRING_D                                                          !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_SPRING_D( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_SPRING_D( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables_i,  &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in spring_d definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in spring_d definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_SPRING_D

  !-----------------------------------------------------------------------------!
  !> Read in !SPRING_A                                                          !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_SPRING_A( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_SPRING_A( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables_i,  &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in spring_a definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in spring_a definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_SPRING_A

  !-----------------------------------------------------------------------------!
  !> Read in !DASHPOT_D                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_DASHPOT_D( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_DASHPOT_D( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables_i,  &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in spring_d definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in spring_d definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_DASHPOT_D

  !-----------------------------------------------------------------------------!
  !> Read in !DASHPOT_A                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_DASHPOT_A( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_DASHPOT_A( ctrl,                                        &
        fstrSOLID%materials(cid)%mtype,       &
        fstrSOLID%materials(cid)%nlgeom_flag, &
        fstrSOLID%materials(cid)%variables_i,  &
        fstrSOLID%materials(cid)%dict)/=0 ) then
      write(*,*) '### Error: Fail in read in spring_a definition : ' , cid
      write(ILOG,*) '### Error: Fail in read in spring_a definition : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_DASHPOT_A

  !-----------------------------------------------------------------------------!
  !> Read in !USER_MATERIAL                                                     !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_USER_MATERIAL( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: cid

    cid = fstr_setup_current_material( ctrl, P )
    fstrSOLID => P%SOLID

    if( fstr_ctrl_get_USERMATERIAL( ctrl, fstrSOLID%materials(cid)%mtype,   &
        fstrSOLID%materials(cid)%nlgeom_flag, fstrSOLID%materials(cid)%nfstatus, &
        fstrSOLID%materials(cid)%variables(101:) )/=0 ) then
      write(*,*) '### Error: Fail in read in user defined material : ' , cid
      write(ILOG,*) '### Error: Fail in read in user defined material : ', cid
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_USER_MATERIAL

  !-----------------------------------------------------------------------------!
  !> Read in !OUTPUT_RES                                                        !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_OUTPUT_RES( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: c_output, i

    hecMESH => P%MESH
    fstrSOLID => P%SOLID

    c_output=3
    if( .not. fstr_ctrl_get_outitem( ctrl, hecMESH, fstrSOLID%output_ctrl(c_output)%outinfo ) ) then
      write(*,*) '### Error: Fail in read in node output definition : ' , c_output
      write(ILOG,*) '### Error: Fail in read in node output definition : ', c_output
      stop HECMW_EXIT_INPUT
    endif
    if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name /= 'ALL' ) then
      c_output=2
      do i=1,hecMESH%node_group%n_grp
        if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name == hecMESH%node_group%grp_name(i) ) then
          fstrSOLID%output_ctrl(c_output)%outinfo%grp_id = i; exit
        endif
      enddo
    endif
  end subroutine fstr_setup_OUTPUT_RES

  !-----------------------------------------------------------------------------!
  !> Read in !OUTPUT_VIS                                                        !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_OUTPUT_VIS( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: c_output, i

    hecMESH => P%MESH
    fstrSOLID => P%SOLID

    c_output=4
    if( .not. fstr_ctrl_get_outitem( ctrl, hecMESH, fstrSOLID%output_ctrl(c_output)%outinfo ) ) then
      write(*,*) '### Error: Fail in read in element output definition : ' , c_output
      write(ILOG,*) '### Error: Fail in read in element output definition : ', c_output
      stop HECMW_EXIT_INPUT
    endif
    if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name /= 'ALL' ) then
      c_output=2
      do i=1,hecMESH%node_group%n_grp
        if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name == hecMESH%node_group%grp_name(i) ) then
          fstrSOLID%output_ctrl(c_output)%outinfo%grp_id = i; exit
        endif
      enddo
    endif
  end subroutine fstr_setup_OUTPUT_VIS

  !-----------------------------------------------------------------------------!
  !> Read in !ULOAD                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_ULOAD( ctrl )
    implicit none
    integer(kind=kint) :: ctrl

    if( fstr_ctrl_get_USERLOAD( ctrl )/=0 ) then
      write(*,*) '### Error: Fail in read in ULOAD definition : '
      write(ILOG,*) '### Error: Fail in read in ULOAD definition : '
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_ULOAD

  !-----------------------------------------------------------------------------!
  !> Append !BOUNDARY                                                           !
  !-----------------------------------------------------------------------------!
  !> With n = 0, only the group of rotation center is counted.
  subroutine fstr_solid_append_boundary( hecMESH, fstrSOLID, header_name, gid, istot, rotc_name, &
      amp_name, n, grp_names, dof_ids, dof_ide, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: istot
    character(len=HECMW_NAME_LEN), intent(in) :: rotc_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dof_ids(:)
    integer(kind=kint), intent(in)            :: dof_ide(:)
    real(kind=kreal), intent(in)              :: vals(:)

    character(len=HECMW_NAME_LEN) :: rotc(1)
    integer(kind=kint) :: amp_id, rotc_id(1), n_rotc
    integer(kind=kint) :: i, old_size, new_size

    rotc_id = -1
    n_rotc = -1
    if( rotc_name /= ' ' ) then
      if( istot /= 0 ) call fstr_setup_util_err_stop( 'fstr control file error : '//header_name// &
        ' : rotational boundary cannot be specified with total value' )
      fstrSOLID%BOUNDARY_ngrp_rot = fstrSOLID%BOUNDARY_ngrp_rot + 1
      n_rotc = fstrSOLID%BOUNDARY_ngrp_rot
      rotc(1) = rotc_name
      call node_grp_name_to_id_ex( hecMESH, header_name//',ROT_CENTER=', 1, rotc, rotc_id )
    endif

    if( n == 0 ) return
    old_size = fstrSOLID%BOUNDARY_ngrp_tot
    new_size = old_size + n
    fstrSOLID%BOUNDARY_ngrp_tot = new_size
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_ID,    old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_type,  old_size, new_size )
    call hecmw_expand_real_array    (fstrSOLID%BOUNDARY_ngrp_val,   old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_amp,   old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_istot, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_rotID, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%BOUNDARY_ngrp_centerID, old_size, new_size )

    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    fstrSOLID%BOUNDARY_ngrp_GRPID(old_size+1:new_size) = gid
    call node_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrSOLID%BOUNDARY_ngrp_ID(old_size+1:))
    fstrSOLID%BOUNDARY_ngrp_istot(old_size+1:new_size) = istot

    ! set up information about rotation ( default value is set if ROT_CENTER is not given.)
    fstrSOLID%BOUNDARY_ngrp_rotID(old_size+1:) = n_rotc
    fstrSOLID%BOUNDARY_ngrp_centerID(old_size+1:) = rotc_id(1)

    do i = 1, n
      if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
        call fstr_setup_util_err_stop( 'fstr control file error : '//header_name// &
          ' : range of dof_ids and dof_ide is from 1 to 6' )
      end if
      fstrSOLID%BOUNDARY_ngrp_val(old_size+i) = vals(i)
      fstrSOLID%BOUNDARY_ngrp_type(old_size+i) = 10 * dof_ids(i) + dof_ide(i)
      fstrSOLID%BOUNDARY_ngrp_amp(old_size+i) = amp_id
    end do

  end subroutine fstr_solid_append_boundary

  !-----------------------------------------------------------------------------!
  !> Append !CLOAD                                                              !
  !-----------------------------------------------------------------------------!
  !> With n = 0, only the group of rotation center is counted.
  subroutine fstr_solid_append_cload( hecMESH, fstrSOLID, header_name, gid, rotc_name, &
      amp_name, n, grp_names, dofs, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    character(len=HECMW_NAME_LEN), intent(in) :: rotc_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dofs(:)
    real(kind=kreal), intent(in)              :: vals(:)

    character(len=HECMW_NAME_LEN) :: rotc(1)
    integer(kind=kint) :: amp_id, rotc_id(1), n_rotc
    integer(kind=kint) :: i, old_size, new_size
    integer(kind=kint) :: ngrp_id, nnode_grp
    character(len=HECMW_MSG_LEN) :: msg(2)

    rotc_id = -1
    n_rotc = -1
    if( rotc_name /= ' ' ) then
      fstrSOLID%CLOAD_ngrp_rot = fstrSOLID%CLOAD_ngrp_rot + 1
      n_rotc = fstrSOLID%CLOAD_ngrp_rot
      rotc(1) = rotc_name
      call node_grp_name_to_id_ex( hecMESH, header_name//',ROT_CENTER=', 1, rotc, rotc_id )
    endif

    if( n == 0 ) return
    old_size = fstrSOLID%CLOAD_ngrp_tot
    new_size = old_size + n
    fstrSOLID%CLOAD_ngrp_tot = new_size
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_ID,  old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_DOF, old_size, new_size )
    call hecmw_expand_real_array    ( fstrSOLID%CLOAD_ngrp_val, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_amp, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_rotID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%CLOAD_ngrp_centerID, old_size, new_size )

    ! set up information about torque load ( default value is set if ROT_CENTER is not given.)
    fstrSOLID%CLOAD_ngrp_rotID(old_size+1:) = n_rotc
    fstrSOLID%CLOAD_ngrp_centerID(old_size+1:) = rotc_id(1)

    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    do i=1,n
      fstrSOLID%CLOAD_ngrp_amp(old_size+i) = amp_id
      fstrSOLID%CLOAD_ngrp_DOF(old_size+i) = dofs(i)
      fstrSOLID%CLOAD_ngrp_val(old_size+i) = vals(i)
    end do
    fstrSOLID%CLOAD_ngrp_GRPID(old_size+1:new_size) = gid
    call node_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrSOLID%CLOAD_ngrp_ID(old_size+1:))

    ! the refinement adds nodes to the group, so the load given to each node of the group
    ! sums up to more than intended before the refinement
    if( hecMESH%n_refine > 0 ) then
      do i=1,n
        ngrp_id = fstrSOLID%CLOAD_ngrp_ID(old_size+i)
        nnode_grp = hecmw_ngrp_get_number(hecMESH, ngrp_id)
        if( nnode_grp > 1 .and. hecMESH%my_rank == 0 ) then
          msg(1) = ' ### Warning: '//header_name//' on a refined mesh'
          write(msg(2),'(a,a,a,i0,a)') '   node group "', trim(hecMESH%node_group%grp_name(ngrp_id)), &
            '" now has ', nnode_grp, ' nodes; the same load is applied to each of them'
          write(*,'(/a/a/)') trim(msg(1)), trim(msg(2))
          write(ILOG,'(/a/a/)') trim(msg(1)), trim(msg(2))
        endif
      enddo
    endif

  end subroutine fstr_solid_append_cload

  !-----------------------------------------------------------------------------!
  !> Reset !DLOAD                                                        !
  !-----------------------------------------------------------------------------!

  subroutine fstr_expand_dload_array( array, old_size, new_size )
    implicit none
    real(kind=kreal), pointer :: array(:,:)
    integer(kind=kint) :: old_size, new_size, i, j
    real(kind=kreal), pointer :: temp(:,:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(0:6, old_size))
      temp = array
      deallocate(array)
      allocate(array(0:6, new_size))
      array = 0
      do i=1,old_size
        do j=0,6
          array(j,i) = temp(j,i)
        end do
      end do
      deallocate(temp)
    else
      allocate(array(0:6, new_size))
      array = 0
    end if
  end subroutine fstr_expand_dload_array

  !-----------------------------------------------------------------------------!
  !> Append !DLOAD                                                              !
  !-----------------------------------------------------------------------------!
  !> follow is one value shared by all the DLOAD conditions; each call with n > 0 overwrites it.
  subroutine fstr_solid_append_dload( hecMESH, fstrSOLID, header_name, gid, follow, &
      amp_name, n, grp_names, lids, params )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: follow
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: lids(:)
    real(kind=kreal), intent(in)              :: params(0:,:)

    logical, allocatable :: fg_surface(:)
    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%DLOAD_ngrp_tot
    new_size = old_size + n
    fstrSOLID%DLOAD_ngrp_tot = new_size
    call hecmw_expand_integer_array ( fstrSOLID%DLOAD_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%DLOAD_ngrp_ID,  old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%DLOAD_ngrp_LID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%DLOAD_ngrp_amp, old_size, new_size )
    call fstr_expand_dload_array ( fstrSOLID%DLOAD_ngrp_params, old_size, new_size )

    allocate( fg_surface(n) )
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    fstrSOLID%DLOAD_follow = follow
    do i=1,n
      fstrSOLID%DLOAD_ngrp_amp(old_size+i) = amp_id
      fstrSOLID%DLOAD_ngrp_LID(old_size+i) = lids(i)
      do j=0, 6
        fstrSOLID%DLOAD_ngrp_params(j,old_size+i) = params(j,i)
      end do
      fg_surface(i) =  ( lids(i) == 100 )
    end do
    fstrSOLID%DLOAD_ngrp_GRPID(old_size+1:new_size) = gid
    call dload_grp_name_to_id_ex( hecMESH, n, grp_names, fg_surface, fstrSOLID%DLOAD_ngrp_ID(old_size+1:))
    deallocate( fg_surface )

  end subroutine fstr_solid_append_dload

  !-----------------------------------------------------------------------------!
  !> Append !TEMPERATURE                                                        !
  !-----------------------------------------------------------------------------!
  !> Not used for the temperature read from a result (TEMP_irres > 0).
  subroutine fstr_solid_append_temperature( hecMESH, fstrSOLID, header_name, gid, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: i, old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%TEMP_ngrp_tot
    new_size = old_size + n
    call hecmw_expand_integer_array ( fstrSOLID%TEMP_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%TEMP_ngrp_ID, old_size, new_size )
    call hecmw_expand_real_array    ( fstrSOLID%TEMP_ngrp_val,old_size, new_size )

    do i = 1, n
      fstrSOLID%TEMP_ngrp_val(old_size+i) = vals(i)
    enddo
    fstrSOLID%TEMP_ngrp_GRPID(old_size+1:new_size) = gid
    fstrSOLID%TEMP_ngrp_tot = new_size
    call node_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrSOLID%TEMP_ngrp_ID(old_size+1:))

  end subroutine fstr_solid_append_temperature

  !-----------------------------------------------------------------------------!
  !> Append !SPRING                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_solid_append_spring( hecMESH, fstrSOLID, header_name, gid, incremental, &
      amp_name, n, grp_names, dofs, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: incremental
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dofs(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%SPRING_ngrp_tot
    new_size = old_size + n
    fstrSOLID%SPRING_ngrp_tot = new_size
    call hecmw_expand_integer_array ( fstrSOLID%SPRING_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%SPRING_ngrp_ID,  old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%SPRING_ngrp_DOF, old_size, new_size )
    call hecmw_expand_real_array    ( fstrSOLID%SPRING_ngrp_val, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%SPRING_ngrp_amp, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%SPRING_incremental, old_size, new_size )

    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    do i=1,n
      fstrSOLID%SPRING_ngrp_amp(old_size+i) = amp_id
      fstrSOLID%SPRING_ngrp_DOF(old_size+i) = dofs(i)
      fstrSOLID%SPRING_ngrp_val(old_size+i) = vals(i)
    end do
    fstrSOLID%SPRING_ngrp_GRPID(old_size+1:new_size) = gid
    call node_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrSOLID%SPRING_ngrp_ID(old_size+1:))
    fstrSOLID%SPRING_incremental(old_size+1:new_size) = incremental

  end subroutine fstr_solid_append_spring

  !-----------------------------------------------------------------------------!
  !> Append !VELOCITY                                                           !
  !-----------------------------------------------------------------------------!
  !> With n = 0, only the group of rotation center is counted.
  !> vtype is one value shared by all the VELOCITY conditions; each call with n > 0 overwrites it.
  subroutine fstr_solid_append_velocity( hecMESH, fstrSOLID, fstrDYNAMIC, header_name, gid, vtype, &
      rotc_name, amp_name, n, grp_names, dof_ids, dof_ide, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    type(fstr_dynamic)                        :: fstrDYNAMIC
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: vtype
    character(len=HECMW_NAME_LEN), intent(in) :: rotc_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dof_ids(:)
    integer(kind=kint), intent(in)            :: dof_ide(:)
    real(kind=kreal), intent(in)              :: vals(:)

    character(len=HECMW_NAME_LEN) :: rotc(1)
    integer(kind=kint) :: amp_id, rotc_id(1), n_rotc
    integer(kind=kint) :: i, j, old_size, new_size

    rotc_id = -1
    n_rotc = -1
    if( rotc_name /= ' ' ) then
      fstrSOLID%VELOCITY_ngrp_rot = fstrSOLID%VELOCITY_ngrp_rot + 1
      n_rotc = fstrSOLID%VELOCITY_ngrp_rot
      rotc(1) = rotc_name
      call node_grp_name_to_id_ex( hecMESH, header_name//',ROT_CENTER=', 1, rotc, rotc_id )
    endif

    if( n == 0 ) return
    old_size = fstrSOLID%VELOCITY_ngrp_tot
    new_size = old_size + n
    fstrSOLID%VELOCITY_ngrp_tot = new_size

    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_ID  , old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_type, old_size, new_size )
    call hecmw_expand_real_array    (fstrSOLID%VELOCITY_ngrp_val , old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_amp , old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_rotID, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%VELOCITY_ngrp_centerID, old_size, new_size )

    fstrSOLID%VELOCITY_type = vtype
    if( vtype == kbcInitial ) fstrDYNAMIC%VarInitialize = .true.
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    call node_grp_name_to_id_ex( hecMESH, header_name, &
      n, grp_names, fstrSOLID%VELOCITY_ngrp_ID(old_size+1:))
    fstrSOLID%VELOCITY_ngrp_GRPID(old_size+1:new_size) = gid

    ! set up information about rotation ( default value is set if ROT_CENTER is not given.)
    fstrSOLID%VELOCITY_ngrp_rotID(old_size+1:) = n_rotc
    fstrSOLID%VELOCITY_ngrp_centerID(old_size+1:) = rotc_id(1)

    j = old_size+1
    do i = 1, n
      if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
        call fstr_setup_util_err_stop( 'fstr control file error : '//header_name// &
          ' : range of dof_ids and dof_ide is from 1 to 6' )
      end if
      fstrSOLID%VELOCITY_ngrp_type(j) = 10 * dof_ids(i) + dof_ide(i)
      fstrSOLID%VELOCITY_ngrp_amp(j) = amp_id
      fstrSOLID%VELOCITY_ngrp_val(old_size+i) = vals(i)
      j = j+1
    end do

  end subroutine fstr_solid_append_velocity

  !-----------------------------------------------------------------------------!
  !> Append !ACCELERATION                                                       !
  !-----------------------------------------------------------------------------!
  !> atype is one value shared by all the ACCELERATION conditions; each call with n > 0 overwrites it.
  subroutine fstr_solid_append_acceleration( hecMESH, fstrSOLID, fstrDYNAMIC, header_name, gid, atype, &
      amp_name, n, grp_names, dof_ids, dof_ide, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    type(fstr_dynamic)                        :: fstrDYNAMIC
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: atype
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dof_ids(:)
    integer(kind=kint), intent(in)            :: dof_ide(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%ACCELERATION_ngrp_tot
    new_size = old_size + n
    fstrSOLID%ACCELERATION_ngrp_tot = new_size

    call hecmw_expand_integer_array (fstrSOLID%ACCELERATION_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%ACCELERATION_ngrp_ID  , old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%ACCELERATION_ngrp_type, old_size, new_size )
    call hecmw_expand_real_array    (fstrSOLID%ACCELERATION_ngrp_val , old_size, new_size )
    call hecmw_expand_integer_array (fstrSOLID%ACCELERATION_ngrp_amp , old_size, new_size )

    fstrSOLID%ACCELERATION_type = atype
    if( atype == kbcInitial )fstrDYNAMIC%VarInitialize = .true.
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    call node_grp_name_to_id_ex( hecMESH, header_name, &
      n, grp_names, fstrSOLID%ACCELERATION_ngrp_ID(old_size+1:))
    fstrSOLID%ACCELERATION_ngrp_GRPID(old_size+1:new_size) = gid

    j = old_size+1
    do i = 1, n
      if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
        call fstr_setup_util_err_stop( 'fstr control file error : '//header_name// &
          ' : range of dof_ids and dof_ide is from 1 to 6' )
      end if
      fstrSOLID%ACCELERATION_ngrp_type(j) = 10 * dof_ids(i) + dof_ide(i)
      fstrSOLID%ACCELERATION_ngrp_amp(j) = amp_id
      fstrSOLID%ACCELERATION_ngrp_val(old_size+i) = vals(i)
      j = j+1
    end do

  end subroutine fstr_solid_append_acceleration

  !-----------------------------------------------------------------------------!
  !> Append !COUPLE                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_solid_append_couple( hecMESH, fstrSOLID, fstrPARAM, header_name, n, grp_names )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    type(fstr_param)                          :: fstrPARAM
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)

    integer(kind=kint) :: old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%COUPLE_ngrp_tot
    new_size = old_size + n
    fstrSOLID%COUPLE_ngrp_tot = new_size

    call hecmw_expand_integer_array ( fstrSOLID%COUPLE_ngrp_ID,  old_size, new_size )

    call surf_grp_name_to_id_ex( hecMESH, header_name, &
      n, grp_names, fstrSOLID%COUPLE_ngrp_ID(old_size+1:))

    fstrPARAM%fg_couple = 1

  end subroutine fstr_solid_append_couple

  !-----------------------------------------------------------------------------!
  !> Append !ELEMENT_ACTIVATION                                                 !
  !-----------------------------------------------------------------------------!
  subroutine fstr_solid_append_element_activation( hecMESH, fstrSOLID, header_name, gid, &
      amp_name, eps, measure, state, n, grp_names, thlow, thup )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    real(kind=kreal), intent(in)              :: eps
    integer(kind=kint), intent(in)            :: measure
    integer(kind=kint), intent(in)            :: state
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: thlow(:)
    real(kind=kreal), intent(in)              :: thup(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, old_size, new_size

    if( n == 0 ) return
    old_size = fstrSOLID%elemact%ELEMACT_egrp_tot
    new_size = old_size + n
    fstrSOLID%elemact%ELEMACT_egrp_tot = new_size

    call hecmw_expand_integer_array ( fstrSOLID%elemact%ELEMACT_egrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%elemact%ELEMACT_egrp_ID,  old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%elemact%ELEMACT_egrp_amp, old_size, new_size )
    call hecmw_expand_real_array ( fstrSOLID%elemact%ELEMACT_egrp_eps, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%elemact%ELEMACT_egrp_depends, old_size, new_size )
    call hecmw_expand_real_array ( fstrSOLID%elemact%ELEMACT_egrp_ts_lower, old_size, new_size )
    call hecmw_expand_real_array ( fstrSOLID%elemact%ELEMACT_egrp_ts_upper, old_size, new_size )
    call hecmw_expand_integer_array ( fstrSOLID%elemact%ELEMACT_egrp_state, old_size, new_size )

    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )
    do i=1,n
      fstrSOLID%elemact%ELEMACT_egrp_amp(old_size+i) = amp_id
      fstrSOLID%elemact%ELEMACT_egrp_eps(old_size+i) = eps
    end do
    fstrSOLID%elemact%ELEMACT_egrp_GRPID(old_size+1:new_size) = gid
    fstrSOLID%elemact%ELEMACT_egrp_depends(old_size+1:new_size) = measure
    fstrSOLID%elemact%ELEMACT_egrp_ts_lower(old_size+1:new_size) = thlow(1:n)
    fstrSOLID%elemact%ELEMACT_egrp_ts_upper(old_size+1:new_size) = thup(1:n)
    fstrSOLID%elemact%ELEMACT_egrp_state(old_size+1:new_size) = state

    call elem_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrSOLID%elemact%ELEMACT_egrp_ID(old_size+1:))

  end subroutine fstr_solid_append_element_activation

  !-----------------------------------------------------------------------------!
  !> Initialize !CONTACT                                                        !
  !-----------------------------------------------------------------------------!
  !> Checks and initializes contacts(is+1:is+n), whose pairs are already read. cparam_name is the
  !> name of !CONTACT_PARAM; the default parameters are used if it is not found.
  subroutine fstr_solid_init_contacts( hecMESH, fstrSOLID, fstrPARAM, is, n, cparam_name, smoothing )
    implicit none
    type(hecmwST_local_mesh), pointer         :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    type(fstr_param)                          :: fstrPARAM
    integer(kind=kint), intent(in)            :: is
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: cparam_name
    integer(kind=kint), intent(in)            :: smoothing

    integer(kind=kint) :: i, cparam_id
    integer(kind=kint) :: ss_nslav, ss_nmas_vis, ss_nmas_glob, ss_bad   ! MORTAR=YES visibility guard scratch
    logical :: isOK

    cparam_id = 0
    do i=1,size(fstrPARAM%contactparam)-1
      if( hecmw_streqr( fstrPARAM%contactparam(i)%name, cparam_name ) ) then
        cparam_id = i; exit
      endif
    enddo
    do i=1,n
      fstrSOLID%contacts(is+i)%smoothing = smoothing
      if( .not. fstr_contact_check( fstrSOLID%contacts(is+i), hecMESH ) ) then
        write(*,*) '### Error: Inconsistence in contact and surface definition : ' , i+is
        write(ILOG,*) '### Error: Inconsistence in contact and surface definition : ', i+is
        stop HECMW_EXIT_MODEL
      else
        isOK = fstr_contact_init( fstrSOLID%contacts(is+i), hecMESH, fstrPARAM%contactparam(cparam_id))
        ! !CONTACT_ALGO, FRICTION_CONE=FOLLOW gives up the symmetry of the ALagrange
        ! friction terms; see getContactStiffness_Alag
        fstrSOLID%contacts(is+i)%symmetric = .not. fstrPARAM%fric_cone_follow
        fstrSOLID%contacts(is+i)%eps_fric_band = fstrPARAM%eps_fric_band
      endif
      if( fstrSOLID%contacts(is+i)%method == CONTACTS2S ) then
        ! Only the static path is wired for MORTAR=YES. The implicit dynamic solver calls
        ! the shared scan / matrix connectivity / stiffness routines, so the pair would be
        ! dispatched to the mortar code, but it never calls fstr_begin_lambda_txn /
        ! fstr_commit_lambda_txn and the multipliers would not advance over the substeps.
        ! The explicit solver applies the NODE-SURF scan to every pair regardless of the
        ! method. Both would be silent, so this must be an error.
        ! ELEMCHECK and NZPROF run no analysis at all (fstr_input_precheck returns before
        ! the solver dispatch and never touches the contact pairs), so they are let through.
        if( fstrPARAM%solution_type /= kstSTATIC .and. &
          & fstrPARAM%solution_type /= kstPRECHECK .and. &
          & fstrPARAM%solution_type /= kstNZPROF ) then
          write(*,*)    '### Error: MORTAR=YES is supported only in static analysis : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES is supported only in static analysis : ', i+is
          write(*,*)    '           Specify !SOLUTION, TYPE=STATIC or TYPE=NLSTATIC.'
          write(ILOG,*) '           Specify !SOLUTION, TYPE=STATIC or TYPE=NLSTATIC.'
          stop HECMW_EXIT_MODEL
        endif
        ! Restart is not supported. The restart file keeps the node-level contact state
        ! only (contacts(i)%states over contacts(i)%slave), so the segment states in
        ! slave_surf and the converged multipliers are not saved. A resumed run would
        ! start with every segment free and without the warm start, and the initial scan
        ! is skipped on restart, so the difference would be silent.
        if( fstrSOLID%restart_nout < 0 ) then
          write(*,*)    '### Error: MORTAR=YES does not support restart : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES does not support restart : ', i+is
          write(*,*)    '           The contact state of mortar pairs is not saved in the restart file.'
          write(ILOG,*) '           The contact state of mortar pairs is not saved in the restart file.'
          stop HECMW_EXIT_MODEL
        endif
        if( fstrPARAM%contact_algo == kcaSLagrange ) then
          write(*,*)    '### Error: MORTAR=YES is not supported with !CONTACT_ALGO TYPE=SLAGRANGE : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES is not supported with !CONTACT_ALGO TYPE=SLAGRANGE : ', i+is
          stop HECMW_EXIT_MODEL
        endif
        if( fstrSOLID%contacts(is+i)%surf_id1_sgrp <= 0 ) then
          write(*,*)    '### Error: MORTAR=YES requires !CONTACT_PAIR TYPE=SURF-SURF : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES requires !CONTACT_PAIR TYPE=SURF-SURF : ', i+is
          stop HECMW_EXIT_MODEL
        endif
        ! Only finite sliding is implemented. With any other INTERACTION the per-integration
        ! point tracking is skipped, and with it the activation of the slave segment, so the
        ! assembly sees every segment as free and no contact force is produced at all
        ! (a silent no-op). INTERACTION defaults to SSLID, so this must be an error.
        if( fstrSOLID%contacts(is+i)%algtype /= CONTACTFSLID ) then
          write(*,*)    '### Error: MORTAR=YES is supported only with INTERACTION=FSLID : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES is supported only with INTERACTION=FSLID : ', i+is
          write(*,*)    '           Specify INTERACTION=FSLID on !CONTACT (the default is SSLID).'
          write(ILOG,*) '           Specify INTERACTION=FSLID on !CONTACT (the default is SSLID).'
          stop HECMW_EXIT_MODEL
        endif
        if( fstrSOLID%contacts(is+i)%smoothing /= kcsNONE ) then
          write(*,*)    '### Error: MORTAR=YES is not supported with SMOOTHING= : ', i+is
          write(ILOG,*) '### Error: MORTAR=YES is not supported with SMOOTHING= : ', i+is
          stop HECMW_EXIT_MODEL
        endif
        ! Visibility guard: the mortar integral of a slave segment needs every master face
        ! it can project onto, so a rank owning slave faces must see the entire master
        ! surface (owned + ghost), which !PARTITION, CONTACT_OWNER=SLAVE provides. Verified
        ! by comparing the visible master count with the global owned-only count.
        ! Serial (np=1) skips the block, so the numerical path is unchanged.
        if( hecmw_comm_get_size() > 1 ) then
          ss_nslav     = size( fstrSOLID%contacts(is+i)%slave_surf )
          ss_nmas_vis  = size( fstrSOLID%contacts(is+i)%master )
          ss_nmas_glob = fstrSOLID%contacts(is+i)%n_master_owned
          call hecmw_allreduce_I1( hecMESH, ss_nmas_glob, hecmw_sum )
          ! a slave-owning rank must see every master face; flag incomplete ghost
          if( ss_nslav > 0 .and. ss_nmas_vis /= ss_nmas_glob ) then
            ss_bad = 1
          else
            ss_bad = 0
          endif
          call hecmw_allreduce_I1( hecMESH, ss_bad, hecmw_max )
          if( ss_bad /= 0 ) then
            write(*,*)    '### Error: MORTAR=YES contact pair has incomplete master ghost : ', i+is
            write(ILOG,*) '### Error: MORTAR=YES contact pair has incomplete master ghost : ', i+is
            write(*,*)    '           A rank owning slave faces cannot see the full master surface.'
            write(ILOG,*) '           A rank owning slave faces cannot see the full master surface.'
            write(*,*)    '           Re-partition the mesh with !PARTITION, CONTACT_OWNER=SLAVE,'
            write(ILOG,*) '           Re-partition the mesh with !PARTITION, CONTACT_OWNER=SLAVE,'
            write(*,*)    '           or run serial (np=1).'
            write(ILOG,*) '           or run serial (np=1).'
            stop HECMW_EXIT_MODEL
          endif
        endif
      endif
    enddo

  end subroutine fstr_solid_init_contacts

  !-----------------------------------------------------------------------------!
  !> Initialize !EMBED                                                          !
  !-----------------------------------------------------------------------------!
  !> Checks and initializes embeds(is+1:is+n), whose pairs are already read. cparam_name is the
  !> name of !CONTACT_PARAM; the default parameters are used if it is not found.
  subroutine fstr_solid_init_embeds( hecMESH, fstrSOLID, fstrPARAM, is, n, cparam_name, smoothing )
    implicit none
    type(hecmwST_local_mesh), pointer         :: hecMESH
    type(fstr_solid)                          :: fstrSOLID
    type(fstr_param)                          :: fstrPARAM
    integer(kind=kint), intent(in)            :: is
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: cparam_name
    integer(kind=kint), intent(in)            :: smoothing

    integer(kind=kint) :: i, cparam_id
    logical :: isOK

    cparam_id = 0
    do i=1,size(fstrPARAM%contactparam)-1
      if( hecmw_streqr( fstrPARAM%contactparam(i)%name, cparam_name ) ) then
        cparam_id = i; exit
      endif
    enddo
    do i=1,n
      fstrSOLID%embeds(is+i)%smoothing = smoothing
      if( .not. fstr_contact_check( fstrSOLID%embeds(is+i), hecMESH ) ) then
        write(*,*) '### Error: Inconsistence in contact and surface definition : ' , i+is
        write(ILOG,*) '### Error: Inconsistence in contact and surface definition : ', i+is
        stop HECMW_EXIT_MODEL
      else
        isOK = fstr_embed_init( fstrSOLID%embeds(is+i), hecMESH, fstrPARAM%contactparam(cparam_id))
      endif
    enddo

  end subroutine fstr_solid_init_embeds

end module m_fstr_solid_setup
