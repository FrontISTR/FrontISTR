!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize, allocate and finalize
!! fstr_dynamic, and to read !DYNAMIC into it.
module m_fstr_dynamic_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_ctrl_err_stop
  use fstr_ctrl_dynamic
  implicit none

contains

  !> Initial setting of dynamic calculation
  subroutine fstr_dynamic_init( fstrDYNAMIC )
    use m_fstr
    type(fstr_dynamic) :: fstrDYNAMIC
    fstrDYNAMIC%idx_eqa  = 1
    fstrDYNAMIC%idx_resp = 1
    fstrDYNAMIC%n_step   = 1
    fstrDYNAMIC%t_start  = 0.0
    fstrDYNAMIC%t_curr   = 0.0d0
    fstrDYNAMIC%t_end    = 1.0
    fstrDYNAMIC%t_delta  = 1.0
    fstrDYNAMIC%gamma    = 0.5
    fstrDYNAMIC%beta     = 0.25
    fstrDYNAMIC%alpha    = 0.0
    fstrDYNAMIC%idx_mas  = kMassLumped
    fstrDYNAMIC%idx_dmp  = 1
    fstrDYNAMIC%ray_m    = 0.0
    fstrDYNAMIC%ray_k    = 0.0
    fstrDYNAMIC%restart_nout = 0
    fstrDYNAMIC%nout         = 100
    fstrDYNAMIC%ngrp_monit   = 0
    fstrDYNAMIC%nout_monit   = 1
    fstrDYNAMIC%iout_list(1) = 0
    fstrDYNAMIC%iout_list(2) = 0
    fstrDYNAMIC%iout_list(3) = 0
    fstrDYNAMIC%iout_list(4) = 0
    fstrDYNAMIC%iout_list(5) = 0
    fstrDYNAMIC%iout_list(6) = 0

  end subroutine fstr_dynamic_init

  !> Initial setting of dynamic calculation
  subroutine fstr_dynamic_alloc( hecMESH, fstrDYNAMIC )
    use m_fstr
    type(hecmwST_local_mesh),target :: hecMESH
    type(fstr_dynamic) :: fstrDYNAMIC

    integer :: ierror, ndof,nnod

    ndof=hecMESH%n_dof
    nnod=hecMESH%n_node
    if(fstrDYNAMIC%idx_eqa == 11) then
      allocate( fstrDYNAMIC%DISP(ndof*nnod,3)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate( fstrDYNAMIC%VEL (ndof*nnod,1)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate( fstrDYNAMIC%ACC (ndof*nnod,1)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    else
      allocate( fstrDYNAMIC%DISP(ndof*nnod,2)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate( fstrDYNAMIC%VEL (ndof*nnod,2)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      allocate( fstrDYNAMIC%ACC (ndof*nnod,2)  ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif


    allocate( fstrDYNAMIC%VEC1(ndof*nnod)    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC1>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate( fstrDYNAMIC%VEC2(ndof*nnod)    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC2>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    allocate( fstrDYNAMIC%VEC3(ndof*nnod)    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC3>'
      write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( fstrDYNAMIC%idx_eqa == 2 ) then
      allocate( fstrDYNAMIC%HHT_FORCE(ndof*nnod)    ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, HHT_FORCE>'
        write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrDYNAMIC%HHT_FORCE(:) = 0.0d0
    endif

  end subroutine fstr_dynamic_alloc

  !> Finalizer of fstr_dynamic
  subroutine fstr_dynamic_finalize( fstrDYNAMIC )
    type(fstr_dynamic) :: fstrDYNAMIC

    integer :: ierror
    if( associated(fstrDYNAMIC%DISP) )        &
      deallocate( fstrDYNAMIC%DISP    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( associated(fstrDYNAMIC%VEL) )        &
      deallocate( fstrDYNAMIC%VEL     ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( associated(fstrDYNAMIC%ACC) )        &
      deallocate( fstrDYNAMIC%ACC     ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( associated(fstrDYNAMIC%HHT_FORCE) ) then
      deallocate( fstrDYNAMIC%HHT_FORCE    ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, HHT_FORCE>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(fstrDYNAMIC%VEC1) )        &
      deallocate( fstrDYNAMIC%VEC1    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC1>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( associated(fstrDYNAMIC%VEC2) )        &
      deallocate( fstrDYNAMIC%VEC2    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC2>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    if( associated(fstrDYNAMIC%VEC3) )        &
      deallocate( fstrDYNAMIC%VEC3    ,stat=ierror )
    if( ierror /= 0 ) then
      write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC3>'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm())
    end if

  end subroutine

  !-----------------------------------------------------------------------------!
  !> Read in !DYNAMIC                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_DYNAMIC( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P
    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: grp_id_name(1)
    integer(kind=kint) :: grp_id(1)
    logical :: is_alpha_given

    rcode = fstr_ctrl_get_DYNAMIC( ctrl, &
      P%PARAM%nlgeom,  &
      P%DYN%idx_eqa, &
      P%DYN%idx_resp,&
      P%DYN%n_step,  &
      P%DYN%t_start, &
      P%DYN%t_end,   &
      P%DYN%t_delta, &
      P%DYN%gamma,   &
      P%DYN%beta,    &
      P%DYN%alpha,   &
      P%DYN%idx_mas, &
      P%DYN%idx_dmp, &
      P%DYN%ray_m,   &
      P%DYN%ray_k,   &
      P%DYN%nout,    &
      grp_id_name(1), HECMW_NAME_LEN,  &
      P%DYN%nout_monit,  &
      P%DYN%iout_list )

    if( rcode /= 0) call fstr_ctrl_err_stop

    if( P%DYN%idx_eqa /= 1 .and. P%DYN%idx_eqa /= 2 .and. P%DYN%idx_eqa /= 11 ) then
      write(*,*) '### Error: !DYNAMIC time integration type must be 1 (Newmark-beta), 2 (HHT) or 11 (explicit)'
      write(ILOG,*) '### Error: !DYNAMIC time integration type must be 1 (Newmark-beta), 2 (HHT) or 11 (explicit)'
      call flush(6)
      call flush(ILOG)
      call fstr_ctrl_err_stop
    endif
    if( P%DYN%idx_eqa == 2 .and. P%DYN%idx_resp /= 1 ) then
      write(*,*) '### Error: HHT method is available only for time history analysis'
      write(ILOG,*) '### Error: HHT method is available only for time history analysis'
      call flush(6)
      call flush(ILOG)
      call fstr_ctrl_err_stop
    endif

    is_alpha_given = ( P%DYN%alpha < huge(0.0d0) )
    if( P%DYN%idx_eqa == 2 ) then
      if( .not. is_alpha_given ) then
        write(*,*) '### Warning: !DYNAMIC HHT parameter alpha is not given; Newmark-beta method is used'
        write(ILOG,*) '### Warning: !DYNAMIC HHT parameter alpha is not given; Newmark-beta method is used'
        P%DYN%idx_eqa = 1
        P%DYN%alpha = 0.0d0
      else
        if( P%DYN%alpha < -1.0d0/3.0d0 .or. P%DYN%alpha > 0.0d0 ) then
          write(*,*) '### Warning: !DYNAMIC HHT parameter alpha must be in [-1/3, 0]; alpha = 0 is used'
          write(ILOG,*) '### Warning: !DYNAMIC HHT parameter alpha must be in [-1/3, 0]; alpha = 0 is used'
          P%DYN%alpha = 0.0d0
        endif
        ! gamma and beta given with alpha are not used
        P%DYN%beta  = 0.25d0*(1.0d0-P%DYN%alpha)**2
        P%DYN%gamma = 0.5d0-P%DYN%alpha
        if( myrank == 0 ) write(IMSG,'(a,3(1x,1pe12.5))') &
          'INFO: HHT method with alpha, beta, gamma =', P%DYN%alpha, P%DYN%beta, P%DYN%gamma
      endif
    else
      if( is_alpha_given ) then
        write(*,*) '### Warning: !DYNAMIC HHT parameter alpha is ignored except for the HHT method (2)'
        write(ILOG,*) '### Warning: !DYNAMIC HHT parameter alpha is ignored except for the HHT method (2)'
      endif
      P%DYN%alpha = 0.0d0
    endif

    if( P%DYN%idx_mas /= kMassLumped .and. P%DYN%idx_mas /= kMassConsistent ) then
      write(*,*) '### Error: !DYNAMIC mass matrix type must be 1 (lumped) or 2 (consistent)'
      write(ILOG,*) '### Error: !DYNAMIC mass matrix type must be 1 (lumped) or 2 (consistent)'
      call flush(6)
      call flush(ILOG)
      call fstr_ctrl_err_stop
    endif
    if( P%DYN%idx_mas == kMassConsistent .and. &
        (P%DYN%idx_eqa == 11 .or. P%DYN%idx_resp /= 1) ) then
      write(*,*) '### Error: consistent mass matrix is available only for implicit transient analysis'
      write(ILOG,*) '### Error: consistent mass matrix is available only for implicit transient analysis'
      call flush(6)
      call flush(ILOG)
      call fstr_ctrl_err_stop
    endif

    if (P%DYN%idx_resp == 1) then
      call node_grp_name_to_id_ex( P%MESH, '!DYNAMIC', 1, grp_id_name, grp_id)
      P%DYN%ngrp_monit = grp_id(1)
    else
      read(grp_id_name,*) P%DYN%ngrp_monit
    endif

  end subroutine fstr_setup_DYNAMIC

end module m_fstr_dynamic_setup
