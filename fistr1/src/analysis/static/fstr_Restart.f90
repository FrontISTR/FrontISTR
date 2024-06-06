!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to read in and write out
!>         restart files

module m_fstr_Restart
  use m_utilities
  use m_fstr
  implicit none

contains

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
  subroutine fstr_read_restart_dyna_nl(cstep,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,contactNode)
    !----------------------------------------------------------------------*
    integer, intent(out)                  :: cstep       !< current step
    integer, intent(out), optional        :: contactNode !< number of contact nodes
    type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
    type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
    type ( fstr_dynamic), intent(inout)   :: fstrDYNAMIC
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(1),nif(2),naux(2)
    real(kind=kreal) :: data(2)

    call hecmw_restart_open()

    call hecmw_restart_read_int(restrt_step)
    cstep = restrt_step(1)
    call hecmw_restart_read_real(fstrSOLID%unode)
    call hecmw_restart_read_real(fstrSOLID%QFORCE)

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
        enddo
      enddo
    endif

    call hecmw_restart_read_int(restrt_step)
    fstrDYNAMIC%idx_eqa = restrt_step(1)
    call hecmw_restart_read_real(data)
    fstrDYNAMIC%t_curr = data(1)
    fstrDYNAMIC%strainEnergy = data(2)
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

  end subroutine fstr_read_restart_dyna_nl

  !> write out restart file for nonlinear dynamic analysis
  !----------------------------------------------------------------------*
  subroutine fstr_write_restart_dyna_nl(cstep,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,contactNode)
    !----------------------------------------------------------------------*
    integer, intent(in)                   :: cstep      !< current step
    integer, intent(in), optional         :: contactNode!< number of contact nodes
    type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
    type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
    type ( fstr_dynamic), intent(in)      :: fstrDYNAMIC
    type(fstr_param), intent(in)          :: fstrPARAM

    integer :: i,j,restrt_step(1),nif(2),naux(2)
    real(kind=kreal) :: data(2)

    restrt_step(1) = cstep
    call hecmw_restart_add_int(restrt_step,size(restrt_step))
    call hecmw_restart_add_real(fstrSOLID%unode,size(fstrSOLID%unode))
    call hecmw_restart_add_real(fstrSOLID%QFORCE,size(fstrSOLID%QFORCE))

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

    restrt_step(1) = fstrDYNAMIC%idx_eqa
    call hecmw_restart_add_int(restrt_step,size(restrt_step))
    data(1) = fstrDYNAMIC%t_curr
    data(2) = fstrDYNAMIC%strainEnergy
    call hecmw_restart_add_real(data,size(data))
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
