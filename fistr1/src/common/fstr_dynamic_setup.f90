!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize, allocate and finalize
!! fstr_dynamic.
module m_fstr_dynamic_setup
  use m_fstr
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

  !> Finalizer of fstr_solid
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
    if( associated(fstrDYNAMIC%HHT_FORCE) ) then
      deallocate( fstrDYNAMIC%HHT_FORCE    ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, HHT_FORCE>'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif

  end subroutine

end module m_fstr_dynamic_setup
