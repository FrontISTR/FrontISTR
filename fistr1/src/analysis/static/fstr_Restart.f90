!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to read in and write out
!>         restart files

module m_fstr_Restart
  use m_utilities
  use m_fstr
  use mMechGauss, only: tGaussStatus
  use hecmw, only: hecmw_abort, hecmw_comm_get_comm
  implicit none

contains

  !> Stop restart processing with a rank-zero diagnostic.
  subroutine fstr_restart_abort(message)
    character(len=*), intent(in) :: message

    if( myrank == 0 ) write(*,*) 'ERROR: ', trim(message)
    call hecmw_abort(hecmw_comm_get_comm())
  end subroutine fstr_restart_abort

  !> Validate the allocated nodal arrays required by the version-6 shell block.
  subroutine fstr_validate_shell_nodal_state(fstrSOLID, nnode)
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint), intent(in) :: nnode

    if( .not. associated(fstrSOLID%shell_node_mode) ) then
      call fstr_restart_abort('shell_node_mode is not allocated.')
    else if( size(fstrSOLID%shell_node_mode) /= nnode ) then
      call fstr_restart_abort('shell_node_mode size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_rot_state) ) then
      call fstr_restart_abort('shell_rot_state is not allocated.')
    else if( size(fstrSOLID%shell_rot_state) /= nnode ) then
      call fstr_restart_abort('shell_rot_state size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_ref_triad) ) then
      call fstr_restart_abort('shell_ref_triad is not allocated.')
    else if( size(fstrSOLID%shell_ref_triad) /= 9*nnode ) then
      call fstr_restart_abort('shell_ref_triad size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_triad) ) then
      call fstr_restart_abort('shell_triad is not allocated.')
    else if( size(fstrSOLID%shell_triad) /= 9*nnode ) then
      call fstr_restart_abort('shell_triad size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_triad_bak) ) then
      call fstr_restart_abort('shell_triad_bak is not allocated.')
    else if( size(fstrSOLID%shell_triad_bak) /= 9*nnode ) then
      call fstr_restart_abort('shell_triad_bak size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_dtriad) ) then
      call fstr_restart_abort('shell_dtriad is not allocated.')
    else if( size(fstrSOLID%shell_dtriad) /= 9*nnode ) then
      call fstr_restart_abort('shell_dtriad size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_drill) ) then
      call fstr_restart_abort('shell_drill is not allocated.')
    else if( size(fstrSOLID%shell_drill) /= nnode ) then
      call fstr_restart_abort('shell_drill size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_drill_bak) ) then
      call fstr_restart_abort('shell_drill_bak is not allocated.')
    else if( size(fstrSOLID%shell_drill_bak) /= nnode ) then
      call fstr_restart_abort('shell_drill_bak size does not match the mesh.')
    endif
    if( .not. associated(fstrSOLID%shell_ddrill) ) then
      call fstr_restart_abort('shell_ddrill is not allocated.')
    else if( size(fstrSOLID%shell_ddrill) /= nnode ) then
      call fstr_restart_abort('shell_ddrill size does not match the mesh.')
    endif
  end subroutine fstr_validate_shell_nodal_state

  !> Return true when the model contains shell state that cannot be represented by restart version 5.
  logical function fstr_requires_shell_restart_v6(fstrSOLID)
    type(fstr_solid), intent(in) :: fstrSOLID
    integer :: i

    fstr_requires_shell_restart_v6 = fstrSOLID%has_finite_rotation_kinematics
    if( fstr_requires_shell_restart_v6 ) return
    if( .not. associated(fstrSOLID%elements) ) return
    do i = 1, size(fstrSOLID%elements)
      if( associated(fstrSOLID%elements(i)%shell_layer_gausses) ) then
        fstr_requires_shell_restart_v6 = .true.
        return
      endif
    end do
  end function fstr_requires_shell_restart_v6

  !> Abort instead of silently discarding finite-rotation or layer shell state.
  subroutine fstr_require_shell_restart_version(fstrSOLID, restart_version)
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint), intent(in) :: restart_version

    if( restart_version >= 6 ) return
    if( .not. fstr_requires_shell_restart_v6(fstrSOLID) ) return
    call fstr_restart_abort( &
      'restart VERSION >= 6 is required for finite-rotation or layered shell state.')
  end subroutine fstr_require_shell_restart_version

  !> Store one complete shell layer/thickness Gauss state.
  subroutine fstr_write_restart_shell_gauss(gauss)
    type(tGaussStatus), intent(in) :: gauss
    integer(kind=kint) :: nif(2)
    real(kind=kreal) :: scalars(4)

    nif = 0
    if( associated(gauss%istatus) ) nif(1) = size(gauss%istatus)
    if( associated(gauss%fstatus) ) nif(2) = size(gauss%fstatus)
    scalars = (/ gauss%plstrain, gauss%strain_energy, gauss%strain_energy_bak, gauss%plpotential /)
    call hecmw_restart_add_int(nif, size(nif))
    call hecmw_restart_add_real(gauss%strain, size(gauss%strain))
    call hecmw_restart_add_real(gauss%strain_bak, size(gauss%strain_bak))
    call hecmw_restart_add_real(gauss%stress, size(gauss%stress))
    call hecmw_restart_add_real(gauss%stress_bak, size(gauss%stress_bak))
    call hecmw_restart_add_real(gauss%nqm, size(gauss%nqm))
    call hecmw_restart_add_real(gauss%strain_out, size(gauss%strain_out))
    call hecmw_restart_add_real(gauss%stress_out, size(gauss%stress_out))
    call hecmw_restart_add_real(scalars, size(scalars))
    if( nif(1) > 0 ) call hecmw_restart_add_int(gauss%istatus, nif(1))
    if( nif(2) > 0 ) call hecmw_restart_add_real(gauss%fstatus, nif(2))
  end subroutine fstr_write_restart_shell_gauss

  !> Restore one complete shell layer/thickness Gauss state.
  subroutine fstr_read_restart_shell_gauss(gauss)
    type(tGaussStatus), intent(inout) :: gauss
    integer(kind=kint) :: nif(2), current_nif(2)
    real(kind=kreal) :: scalars(4)

    call hecmw_restart_read_int(nif)
    current_nif = 0
    if( associated(gauss%istatus) ) current_nif(1) = size(gauss%istatus)
    if( associated(gauss%fstatus) ) current_nif(2) = size(gauss%fstatus)
    if( any(current_nif /= nif) ) then
      call fstr_restart_abort('shell material-state dimensions do not match the restart file.')
    endif
    call hecmw_restart_read_real(gauss%strain)
    call hecmw_restart_read_real(gauss%strain_bak)
    call hecmw_restart_read_real(gauss%stress)
    call hecmw_restart_read_real(gauss%stress_bak)
    call hecmw_restart_read_real(gauss%nqm)
    call hecmw_restart_read_real(gauss%strain_out)
    call hecmw_restart_read_real(gauss%stress_out)
    call hecmw_restart_read_real(scalars)
    gauss%plstrain = scalars(1)
    gauss%strain_energy = scalars(2)
    gauss%strain_energy_bak = scalars(3)
    gauss%plpotential = scalars(4)
    if( nif(1) > 0 ) call hecmw_restart_read_int(gauss%istatus)
    if( nif(2) > 0 ) call hecmw_restart_read_real(gauss%fstatus)
  end subroutine fstr_read_restart_shell_gauss

  !> Write committed nodal shell kinematics and all layer/thickness material histories.
  subroutine fstr_write_restart_shell_state(hecMESH, fstrSOLID)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint) :: nodal_meta(2), layer_meta(3)
    integer :: i, j

    nodal_meta = 0
    if( fstrSOLID%has_finite_rotation_kinematics ) then
      call fstr_validate_shell_nodal_state(fstrSOLID, int(hecMESH%n_node, kind=kint))
      nodal_meta = (/ 1_kint, int(hecMESH%n_node, kind=kint) /)
    endif
    call hecmw_restart_add_int(nodal_meta, size(nodal_meta))
    if( nodal_meta(1) == 1 ) then
      call hecmw_restart_add_int(fstrSOLID%shell_node_mode, size(fstrSOLID%shell_node_mode))
      call hecmw_restart_add_int(fstrSOLID%shell_rot_state, size(fstrSOLID%shell_rot_state))
      call hecmw_restart_add_real(fstrSOLID%shell_ref_triad, size(fstrSOLID%shell_ref_triad))
      call hecmw_restart_add_real(fstrSOLID%shell_triad, size(fstrSOLID%shell_triad))
      call hecmw_restart_add_real(fstrSOLID%shell_triad_bak, size(fstrSOLID%shell_triad_bak))
      call hecmw_restart_add_real(fstrSOLID%shell_drill, size(fstrSOLID%shell_drill))
      call hecmw_restart_add_real(fstrSOLID%shell_drill_bak, size(fstrSOLID%shell_drill_bak))
    endif

    do i = 1, hecMESH%n_elem
      layer_meta = 0
      if( associated(fstrSOLID%elements(i)%shell_layer_gausses) ) then
        layer_meta(1) = fstrSOLID%elements(i)%shell_nlayer
        layer_meta(2) = fstrSOLID%elements(i)%shell_nthick
        layer_meta(3) = size(fstrSOLID%elements(i)%shell_layer_gausses)
      endif
      call hecmw_restart_add_int(layer_meta, size(layer_meta))
      do j = 1, layer_meta(3)
        call fstr_write_restart_shell_gauss(fstrSOLID%elements(i)%shell_layer_gausses(j))
      enddo
    enddo
  end subroutine fstr_write_restart_shell_state

  !> Restore committed nodal shell kinematics and all layer/thickness material histories.
  subroutine fstr_read_restart_shell_state(hecMESH, fstrSOLID)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(inout) :: fstrSOLID
    integer(kind=kint) :: nodal_meta(2), layer_meta(3)
    integer :: i, j

    call hecmw_restart_read_int(nodal_meta)
    if( nodal_meta(1) == 1 ) then
      if( .not. fstrSOLID%has_finite_rotation_kinematics ) then
        call fstr_restart_abort('restart contains finite-rotation shell state for an incompatible model.')
      endif
      if( nodal_meta(2) /= hecMESH%n_node ) then
        call fstr_restart_abort('finite-rotation shell restart node count does not match the mesh.')
      endif
      call fstr_validate_shell_nodal_state(fstrSOLID, int(hecMESH%n_node, kind=kint))
      call hecmw_restart_read_int(fstrSOLID%shell_node_mode)
      call hecmw_restart_read_int(fstrSOLID%shell_rot_state)
      call hecmw_restart_read_real(fstrSOLID%shell_ref_triad)
      call hecmw_restart_read_real(fstrSOLID%shell_triad)
      call hecmw_restart_read_real(fstrSOLID%shell_triad_bak)
      call hecmw_restart_read_real(fstrSOLID%shell_drill)
      call hecmw_restart_read_real(fstrSOLID%shell_drill_bak)
      fstrSOLID%shell_dtriad = fstrSOLID%shell_triad
      fstrSOLID%shell_ddrill = fstrSOLID%shell_drill
      fstrSOLID%finite_rotation_state_ready = .true.
    else if( fstrSOLID%has_finite_rotation_kinematics ) then
      call fstr_restart_abort('restart file has no finite-rotation shell nodal state.')
    endif

    do i = 1, hecMESH%n_elem
      call hecmw_restart_read_int(layer_meta)
      if( layer_meta(3) > 0 ) then
        if( .not. associated(fstrSOLID%elements(i)%shell_layer_gausses) ) then
          call fstr_restart_abort( &
            'restart contains shell layer state for an unallocated element.')
        endif
        if( layer_meta(1) /= fstrSOLID%elements(i)%shell_nlayer .or. &
            layer_meta(2) /= fstrSOLID%elements(i)%shell_nthick .or. &
            layer_meta(3) /= size(fstrSOLID%elements(i)%shell_layer_gausses) ) then
          call fstr_restart_abort( &
            'shell layer restart dimensions do not match the current model.')
        endif
      else if( associated(fstrSOLID%elements(i)%shell_layer_gausses) ) then
        call fstr_restart_abort( &
          'restart file has no shell layer state required by the current model.')
      endif
      do j = 1, layer_meta(3)
        call fstr_read_restart_shell_gauss(fstrSOLID%elements(i)%shell_layer_gausses(j))
      enddo
    enddo
  end subroutine fstr_read_restart_shell_state

  !> Read in restart file
  !----------------------------------------------------------------------*
  subroutine fstr_read_restart(cstep,substep,step_count,ctime,dtime,hecMESH,fstrSOLID,fstrPARAM,contactNode)
    !----------------------------------------------------------------------*
    integer, intent(out)                  :: cstep       !< current step
    integer, intent(out)                  :: substep     !< current sub step
    integer, intent(out)                  :: step_count  !< current through step
    real(kind=kreal), intent(out)         :: ctime       !< current time
    real(kind=kreal), intent(out)         :: dtime       !< current time increment
    integer, intent(out)                  :: contactNode !< total number of contact nodes
    type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
    type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(3),nif(2),istat(1),nload_prev(1),naux(2)
    real(kind=kreal) :: times(3)

    call hecmw_restart_open()

    call hecmw_restart_read_int(restrt_step)
    if( fstrPARAM%restart_version >= 5 ) then
      if( myrank == 0 ) write(*,*) 'Reading restart file as new format(>=ver5.0)'
      call hecmw_restart_read_real(times)
      call hecmw_restart_read_int(fstrSOLID%NRstat_i)
      call hecmw_restart_read_real(fstrSOLID%NRstat_r)
      call hecmw_restart_read_int(istat)
    else
      if( myrank == 0 ) write(*,*) 'Reading restart file as old format(<ver5.0)'
    endif
    call hecmw_restart_read_int(nload_prev) !load info at previous step
    if( nload_prev(1)>0 ) then
      allocate(fstrSOLID%step_ctrl_restart%Load(nload_prev(1)))
      call hecmw_restart_read_int(fstrSOLID%step_ctrl_restart%Load)
    endif

    call hecmw_restart_read_real(fstrSOLID%unode)
    call hecmw_restart_read_real(fstrSOLID%unode_bak)
    call hecmw_restart_read_real(fstrSOLID%QFORCE)
    call fstr_require_shell_restart_version(fstrSOLID, fstrPARAM%restart_version)
    if( fstrPARAM%restart_version >= 6 ) call fstr_read_restart_shell_state(hecMESH, fstrSOLID)

    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      do j= 1, size(fstrSOLID%elements(i)%gausses)
        call hecmw_restart_read_int(nif)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain_bak)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress_bak)
        if( nif(1)>0 ) call hecmw_restart_read_int(fstrSOLID%elements(i)%gausses(j)%istatus)
        if( nif(2)>0 ) call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%fstatus)
      enddo
      call hecmw_restart_read_int(naux)
      do j= 1, naux(2)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%aux(:,j))
      enddo
    enddo

    if( associated( fstrSOLID%contacts ) ) then
      call hecmw_restart_read_int(nif)
      contactNode = nif(1)
      do i= 1, fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          call hecmw_restart_read_int(nif)
          fstrSOLID%contacts(i)%states(j)%surface = nif(1)
          fstrSOLID%contacts(i)%states(j)%state = nif(2)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%lpos)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%direction)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%multiplier)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%tangentForce_trial)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%tangentForce_final)
          fstrSOLID%contacts(i)%states(j)%tangentForce(1:3)  = fstrSOLID%contacts(i)%states(j)%tangentForce_final(1:3)
          fstrSOLID%contacts(i)%states(j)%tangentForce1(1:3) = fstrSOLID%contacts(i)%states(j)%tangentForce_final(1:3)
        enddo
      enddo
    endif

    call hecmw_restart_close()

    cstep = restrt_step(1)
    substep = restrt_step(2) + 1
    step_count = restrt_step(3)
    if( fstrPARAM%restart_version >= 5 ) then
      ctime = times(1)
      dtime = times(2)
      fstrSOLID%AutoINC_stat = istat(1)
      if( dabs(times(1)-times(3)) < 1.d-10 ) then
        cstep = cstep + 1
        substep = 1
      endif
      do i=1,size(fstrSOLID%step_ctrl)  !shift start time
        fstrSOLID%step_ctrl(i)%starttime = fstrSOLID%step_ctrl(i)%starttime + times(3)
      end do
    else
      ctime = fstrSOLID%step_ctrl(cstep)%starttime
      ctime = ctime + dble(substep-1)*fstrSOLID%step_ctrl(cstep)%initdt
      dtime = fstrSOLID%step_ctrl(cstep)%initdt
      if( dabs(ctime-fstrSOLID%step_ctrl(cstep)%starttime-fstrSOLID%step_ctrl(cstep)%elapsetime) < 1.d-10 ) then
        cstep = cstep + 1
        substep = 1
      endif
    endif

  end subroutine fstr_read_restart

  !> write out restart file
  !----------------------------------------------------------------------*
  subroutine fstr_write_restart(cstep,cstep_ext,substep,step_count,ctime,dtime,hecMESH,  &
      &  fstrSOLID,fstrPARAM,is_StepFinished,contactNode)
    !----------------------------------------------------------------------*
    integer, intent(in)                   :: cstep       !< current step (internal step id)
    integer, intent(in)                   :: cstep_ext   !< current step (external step id)
    integer, intent(in)                   :: substep    !< current sub step
    integer, intent(in)                   :: step_count !< current through step
    real(kind=kreal), intent(in)          :: ctime       !< current time
    real(kind=kreal), intent(in)          :: dtime       !< current time increment
    logical, intent(in)                   :: is_StepFinished
    integer, intent(in)                   :: contactNode!< number of contact nodes
    type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
    type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(3),nif(2),istat(1),nload_prev(1),naux(2)
    real(kind=kreal) :: times(3)

    restrt_step(1) = cstep_ext
    restrt_step(2) = substep
    restrt_step(3) = step_count
    times(1) = ctime
    times(2) = dtime
    if( is_StepFinished ) then
      times(3) = ctime
    else
      times(3) = fstrSOLID%step_ctrl(cstep)%starttime
    end if
    istat(1) = fstrSOLID%AutoINC_stat
    call hecmw_restart_add_int(restrt_step,size(restrt_step))
    if( fstrPARAM%restart_version >= 5 ) then
      call hecmw_restart_add_real(times,size(times))
      call hecmw_restart_add_int(fstrSOLID%NRstat_i,size(fstrSOLID%NRstat_i))
      call hecmw_restart_add_real(fstrSOLID%NRstat_r,size(fstrSOLID%NRstat_r))
      call hecmw_restart_add_int(istat,1)
    endif
    nload_prev(1) = 0
    if( is_StepFinished ) then
      if( associated(fstrSOLID%step_ctrl(cstep)%Load) ) nload_prev(1) = size(fstrSOLID%step_ctrl(cstep)%Load)
      call hecmw_restart_add_int(nload_prev,1)
      if( nload_prev(1)>0 ) call hecmw_restart_add_int(fstrSOLID%step_ctrl(cstep)%Load,nload_prev(1))
    else
      if( cstep>1 ) then
        if( associated(fstrSOLID%step_ctrl(cstep-1)%Load) ) nload_prev(1) = size(fstrSOLID%step_ctrl(cstep-1)%Load)
        call hecmw_restart_add_int(nload_prev,1)
        if( nload_prev(1)>0 ) call hecmw_restart_add_int(fstrSOLID%step_ctrl(cstep-1)%Load,nload_prev(1))
      else
        if( associated(fstrSOLID%step_ctrl_restart%Load) ) nload_prev(1) = size(fstrSOLID%step_ctrl_restart%Load)
        call hecmw_restart_add_int(nload_prev,1)
        if( nload_prev(1)>0 ) call hecmw_restart_add_int(fstrSOLID%step_ctrl_restart%Load,nload_prev(1))
      endif
    end if

    call hecmw_restart_add_real(fstrSOLID%unode,size(fstrSOLID%unode))
    call hecmw_restart_add_real(fstrSOLID%unode_bak,size(fstrSOLID%unode_bak))
    call hecmw_restart_add_real(fstrSOLID%QFORCE,size(fstrSOLID%QFORCE))
    call fstr_require_shell_restart_version(fstrSOLID, fstrPARAM%restart_version)
    if( fstrPARAM%restart_version >= 6 ) call fstr_write_restart_shell_state(hecMESH, fstrSOLID)

    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      do j= 1, size(fstrSOLID%elements(i)%gausses)
        nif = 0
        if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) ) nif(1)=size(fstrSOLID%elements(i)%gausses(j)%istatus)
        if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) ) nif(2)=size(fstrSOLID%elements(i)%gausses(j)%fstatus)
        call hecmw_restart_add_int(nif,size(nif))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain,size(fstrSOLID%elements(i)%gausses(j)%strain))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain_bak,size(fstrSOLID%elements(i)%gausses(j)%strain_bak))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress,size(fstrSOLID%elements(i)%gausses(j)%stress))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress_bak,size(fstrSOLID%elements(i)%gausses(j)%stress_bak))
        if( nif(1)>0 ) then
          call hecmw_restart_add_int(fstrSOLID%elements(i)%gausses(j)%istatus,size(fstrSOLID%elements(i)%gausses(j)%istatus))
        endif
        if( nif(2)>0 ) then
          call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%fstatus,size(fstrSOLID%elements(i)%gausses(j)%fstatus))
        endif
      enddo
      naux = 0
      if( associated(fstrSOLID%elements(i)%aux) ) naux=shape(fstrSOLID%elements(i)%aux)
      call hecmw_restart_add_int(naux,size(naux))
      do j= 1, naux(2)
        call hecmw_restart_add_real(fstrSOLID%elements(i)%aux(:,j),naux(1))
      enddo
    enddo

    if( associated( fstrSOLID%contacts ) ) then
      nif(1) = contactNode
      call hecmw_restart_add_int(nif,size(nif))
      do i= 1, fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          nif(1) = fstrSOLID%contacts(i)%states(j)%surface
          nif(2) = fstrSOLID%contacts(i)%states(j)%state
          call hecmw_restart_add_int(nif,size(nif))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%lpos,size(fstrSOLID%contacts(i)%states(j)%lpos))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%direction,size(fstrSOLID%contacts(i)%states(j)%direction))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%multiplier,size(fstrSOLID%contacts(i)%states(j)%multiplier))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%tangentForce_trial, &
            size(fstrSOLID%contacts(i)%states(j)%tangentForce_trial))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%tangentForce_final, &
            size(fstrSOLID%contacts(i)%states(j)%tangentForce_final))
        enddo
      enddo
    endif

    call hecmw_restart_write()

  end subroutine fstr_write_restart

  !> Read in restart file for nonlinear dynamic analysis
  !----------------------------------------------------------------------*
  subroutine fstr_read_restart_dyna_nl(cstep,substep,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,contactNode,step_count)
    !----------------------------------------------------------------------*
    integer, intent(out)                  :: cstep       !< current step
    integer, intent(out)                  :: substep     !< current substep
    integer, intent(out), optional        :: contactNode !< number of contact nodes
    integer, intent(out), optional        :: step_count  !< total step count for result output
    type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
    type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
    type ( fstr_dynamic), intent(inout)   :: fstrDYNAMIC
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(3),nif(2),istat(1),nload_prev(1),naux(2),dyna_int(1)
    real(kind=kreal) :: times(3),dyna_data(1)

    call hecmw_restart_open()

    !--- Common header (same format as static restart) ---
    call hecmw_restart_read_int(restrt_step)
    if( fstrPARAM%restart_version >= 5 ) then
      if( myrank == 0 ) write(*,*) 'Reading dynamic restart file as new format(>=ver5.0)'
      call hecmw_restart_read_real(times)
      call hecmw_restart_read_int(fstrSOLID%NRstat_i)
      call hecmw_restart_read_real(fstrSOLID%NRstat_r)
      call hecmw_restart_read_int(istat)
    else
      if( myrank == 0 ) write(*,*) 'Reading dynamic restart file as old format(<ver5.0)'
    endif
    call hecmw_restart_read_int(nload_prev)
    if( nload_prev(1)>0 ) then
      allocate(fstrSOLID%step_ctrl_restart%Load(nload_prev(1)))
      call hecmw_restart_read_int(fstrSOLID%step_ctrl_restart%Load)
    endif

    call hecmw_restart_read_real(fstrSOLID%unode)
    call hecmw_restart_read_real(fstrSOLID%unode_bak)
    call hecmw_restart_read_real(fstrSOLID%QFORCE)
    call fstr_require_shell_restart_version(fstrSOLID, fstrPARAM%restart_version)
    if( fstrPARAM%restart_version >= 6 ) call fstr_read_restart_shell_state(hecMESH, fstrSOLID)

    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      do j= 1, size(fstrSOLID%elements(i)%gausses)
        call hecmw_restart_read_int(nif)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain_bak)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress_bak)
        if( nif(1)>0 ) call hecmw_restart_read_int(fstrSOLID%elements(i)%gausses(j)%istatus)
        if( nif(2)>0 ) call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%fstatus)
      enddo
      call hecmw_restart_read_int(naux)
      do j= 1, naux(2)
        call hecmw_restart_read_real(fstrSOLID%elements(i)%aux(:,j))
      enddo
    enddo

    if(present(contactNode)) then
      call hecmw_restart_read_int(nif)
      contactNode = nif(1)
      do i= 1, fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          call hecmw_restart_read_int(nif)
          fstrSOLID%contacts(i)%states(j)%surface = nif(1)
          fstrSOLID%contacts(i)%states(j)%state = nif(2)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%lpos)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%direction)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%multiplier)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%tangentForce_trial)
          call hecmw_restart_read_real(fstrSOLID%contacts(i)%states(j)%tangentForce_final)
          fstrSOLID%contacts(i)%states(j)%tangentForce(1:3)  = fstrSOLID%contacts(i)%states(j)%tangentForce_final(1:3)
          fstrSOLID%contacts(i)%states(j)%tangentForce1(1:3) = fstrSOLID%contacts(i)%states(j)%tangentForce_final(1:3)
        enddo
      enddo
    endif

    !--- Dynamic-specific data ---
    call hecmw_restart_read_int(dyna_int)
    fstrDYNAMIC%idx_eqa = dyna_int(1)
    call hecmw_restart_read_real(dyna_data)
    fstrDYNAMIC%strainEnergy = dyna_data(1)
    if( fstrDYNAMIC%idx_eqa == 1 ) then
      call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
      call hecmw_restart_read_real(fstrDYNAMIC%VEL(:,1))
      call hecmw_restart_read_real(fstrDYNAMIC%ACC(:,1))
    else
      call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
      call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,3))
    endif
    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      call hecmw_restart_read_real(fstrSOLID%elements(i)%equiForces)
    enddo

    call hecmw_restart_close()

    !--- Restore step info (same logic as static restart) ---
    cstep = restrt_step(1)
    substep = restrt_step(2) + 1
    if( present(step_count) ) step_count = restrt_step(3)
    if( fstrPARAM%restart_version >= 5 ) then
      fstrDYNAMIC%t_curr = times(1)
      fstrDYNAMIC%t_delta = times(2)
      fstrSOLID%AutoINC_stat = istat(1)
      if( dabs(times(1)-times(3)) < 1.d-10 ) then
        cstep = cstep + 1
        substep = 1
      endif
      do i=1,size(fstrSOLID%step_ctrl)
        fstrSOLID%step_ctrl(i)%starttime = fstrSOLID%step_ctrl(i)%starttime + times(3)
      end do
    else
      fstrDYNAMIC%t_curr = fstrSOLID%step_ctrl(cstep)%starttime
      fstrDYNAMIC%t_curr = fstrDYNAMIC%t_curr + dble(substep-1)*fstrSOLID%step_ctrl(cstep)%initdt
      fstrDYNAMIC%t_delta = fstrSOLID%step_ctrl(cstep)%initdt
      if( dabs(fstrDYNAMIC%t_curr-fstrSOLID%step_ctrl(cstep)%starttime &
        -fstrSOLID%step_ctrl(cstep)%elapsetime) < 1.d-10 ) then
        cstep = cstep + 1
        substep = 1
      endif
    endif

  end subroutine fstr_read_restart_dyna_nl

  !> write out restart file for nonlinear dynamic analysis
  !----------------------------------------------------------------------*
  subroutine fstr_write_restart_dyna_nl(cstep,substep,hecMESH,fstrSOLID, &
    fstrDYNAMIC,fstrPARAM,is_StepFinished,contactNode,step_count)
    !----------------------------------------------------------------------*
    integer, intent(in)                   :: cstep      !< current step
    integer, intent(in)                   :: substep    !< current substep
    logical, intent(in)                   :: is_StepFinished !< whether the step has finished
    integer, intent(in), optional         :: contactNode!< number of contact nodes
    integer, intent(in), optional         :: step_count !< total step count for result output
    type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
    type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
    type ( fstr_dynamic), intent(in)      :: fstrDYNAMIC
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(3),nif(2),istat(1),nload_prev(1),naux(2),dyna_int(1)
    real(kind=kreal) :: times(3),dyna_data(1)

    !--- Common header (same format as static restart) ---
    restrt_step(1) = cstep
    restrt_step(2) = substep
    if( present(step_count) ) then
      restrt_step(3) = step_count
    else
      restrt_step(3) = 0
    endif
    call hecmw_restart_add_int(restrt_step,size(restrt_step))
    if( fstrPARAM%restart_version >= 5 ) then
      times(1) = fstrDYNAMIC%t_curr
      times(2) = fstrDYNAMIC%t_delta
      if( is_StepFinished ) then
        times(3) = fstrDYNAMIC%t_curr
      else
        times(3) = fstrSOLID%step_ctrl(cstep)%starttime
      end if
      call hecmw_restart_add_real(times,size(times))
      call hecmw_restart_add_int(fstrSOLID%NRstat_i,size(fstrSOLID%NRstat_i))
      call hecmw_restart_add_real(fstrSOLID%NRstat_r,size(fstrSOLID%NRstat_r))
      istat(1) = fstrSOLID%AutoINC_stat
      call hecmw_restart_add_int(istat,1)
    endif
    nload_prev(1) = 0
    call hecmw_restart_add_int(nload_prev,1)

    call hecmw_restart_add_real(fstrSOLID%unode,size(fstrSOLID%unode))
    call hecmw_restart_add_real(fstrSOLID%unode_bak,size(fstrSOLID%unode_bak))
    call hecmw_restart_add_real(fstrSOLID%QFORCE,size(fstrSOLID%QFORCE))
    call fstr_require_shell_restart_version(fstrSOLID, fstrPARAM%restart_version)
    if( fstrPARAM%restart_version >= 6 ) call fstr_write_restart_shell_state(hecMESH, fstrSOLID)

    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      do j= 1, size(fstrSOLID%elements(i)%gausses)
        nif = 0
        if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) ) nif(1)=size(fstrSOLID%elements(i)%gausses(j)%istatus)
        if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) ) nif(2)=size(fstrSOLID%elements(i)%gausses(j)%fstatus)
        call hecmw_restart_add_int(nif,size(nif))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain,size(fstrSOLID%elements(i)%gausses(j)%strain))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain_bak,size(fstrSOLID%elements(i)%gausses(j)%strain_bak))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress,size(fstrSOLID%elements(i)%gausses(j)%stress))
        call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress_bak,size(fstrSOLID%elements(i)%gausses(j)%stress_bak))
        if( nif(1)>0 ) then
          call hecmw_restart_add_int(fstrSOLID%elements(i)%gausses(j)%istatus,size(fstrSOLID%elements(i)%gausses(j)%istatus))
        endif
        if( nif(2)>0 ) then
          call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%fstatus,size(fstrSOLID%elements(i)%gausses(j)%fstatus))
        endif
      enddo
      naux = 0
      if( associated(fstrSOLID%elements(i)%aux) ) naux=shape(fstrSOLID%elements(i)%aux)
      call hecmw_restart_add_int(naux,size(naux))
      do j= 1, naux(2)
        call hecmw_restart_add_real(fstrSOLID%elements(i)%aux(:,j),naux(1))
      enddo
    enddo

    if(present(contactNode)) then
      nif(1) = contactNode
      call hecmw_restart_add_int(nif,size(nif))
      do i= 1, fstrSOLID%n_contacts
        do j= 1, size(fstrSOLID%contacts(i)%slave)
          nif(1) = fstrSOLID%contacts(i)%states(j)%surface
          nif(2) = fstrSOLID%contacts(i)%states(j)%state
          call hecmw_restart_add_int(nif,size(nif))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%lpos,size(fstrSOLID%contacts(i)%states(j)%lpos))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%direction,size(fstrSOLID%contacts(i)%states(j)%direction))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%multiplier,size(fstrSOLID%contacts(i)%states(j)%multiplier))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%tangentForce_trial, &
            size(fstrSOLID%contacts(i)%states(j)%tangentForce_trial))
          call hecmw_restart_add_real(fstrSOLID%contacts(i)%states(j)%tangentForce_final, &
            size(fstrSOLID%contacts(i)%states(j)%tangentForce_final))
        enddo
      enddo
    endif

    !--- Dynamic-specific data ---
    dyna_int(1) = fstrDYNAMIC%idx_eqa
    call hecmw_restart_add_int(dyna_int,1)
    dyna_data(1) = fstrDYNAMIC%strainEnergy
    call hecmw_restart_add_real(dyna_data,size(dyna_data))
    if( fstrDYNAMIC%idx_eqa == 1 ) then
      call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
      call hecmw_restart_add_real(fstrDYNAMIC%VEL(:,1),size(fstrDYNAMIC%VEL(:,1)))
      call hecmw_restart_add_real(fstrDYNAMIC%ACC(:,1),size(fstrDYNAMIC%ACC(:,1)))
    else
      call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
      call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,3),size(fstrDYNAMIC%DISP(:,3)))
    endif
    do i= 1, hecMESH%n_elem
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      call hecmw_restart_add_real(fstrSOLID%elements(i)%equiForces,size(fstrSOLID%elements(i)%equiForces))
    enddo

    call hecmw_restart_write()

  end subroutine fstr_write_restart_dyna_nl

end module m_fstr_restart
