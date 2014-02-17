!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato(Advancesoft), X. YUAN(AdavanceSoft)    !
!                       Z. Sun(ASTOM)                                  !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to read in and write out
!>         restart fiels
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/12/26
!>  \version    0.00
!>  \author     Z. Sun(ASTOM)
!>  \date       2011/11   
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Restart
   use m_utilities
   use m_fstr
   implicit none

   contains

!> Read in restart file
!----------------------------------------------------------------------*
      subroutine fstr_read_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM,contactNode)
!----------------------------------------------------------------------*
      integer, intent(out)                  :: cstep       !< current step
      integer, intent(out)                  :: substep     !< current sub step
      integer, intent(out)                  :: step_count  !< current through step
      integer, intent(out), optional        :: contactNode !< total number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
      type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
      type(fstr_param), intent(in)          :: fstrPARAM

      integer :: i,j,restrt_step(3),nif(2)

      call hecmw_restart_open()

      call hecmw_restart_read_int(restrt_step)
      call hecmw_restart_read_real(fstrSOLID%unode)
      call hecmw_restart_read_real(fstrSOLID%QFORCE)

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            call hecmw_restart_read_int(nif)
            call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain)
            call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress)
            if( nif(1)>0 ) call hecmw_restart_read_int(fstrSOLID%elements(i)%gausses(j)%istatus)
            if( nif(2)>0 ) call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%fstatus)
        enddo
      enddo

      if(present(contactNode)) then
        call hecmw_restart_read_int(nif)
        contactNode = nif(1)
        do i= 1, size(fstrSOLID%contacts)
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
      if(fstrPARAM%solution_type==kstNLSTATIC) then
        if(restrt_step(2)==fstrSOLID%step_ctrl(cstep)%num_substep) then
          cstep = cstep + 1
          substep = 1
        endif 
      endif

      end subroutine fstr_read_restart

!> write out restart file
!----------------------------------------------------------------------*
      subroutine fstr_write_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM,contactNode)
!----------------------------------------------------------------------*
      integer, intent(in)                   :: cstep      !< current step
      integer, intent(in)                   :: substep    !< current sub step
      integer, intent(in)                   :: step_count !< current through step
      integer, intent(in),optional          :: contactNode!< number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
      type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
      type(fstr_param), intent(in)          :: fstrPARAM

      integer :: i,j,restrt_step(3),nif(2)

      restrt_step(1) = cstep
      restrt_step(2) = substep
	  restrt_step(3) = step_count
      call hecmw_restart_add_int(restrt_step,size(restrt_step))
      call hecmw_restart_add_real(fstrSOLID%unode,size(fstrSOLID%unode))
      call hecmw_restart_add_real(fstrSOLID%QFORCE,size(fstrSOLID%QFORCE))

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            nif = 0
            if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) ) nif(1)=size(fstrSOLID%elements(i)%gausses(j)%istatus)
            if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) ) nif(2)=size(fstrSOLID%elements(i)%gausses(j)%fstatus)
            call hecmw_restart_add_int(nif,size(nif))
            call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain,size(fstrSOLID%elements(i)%gausses(j)%strain))
            call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress,size(fstrSOLID%elements(i)%gausses(j)%stress))
            if( nif(1)>0 ) then
              call hecmw_restart_add_int(fstrSOLID%elements(i)%gausses(j)%istatus,size(fstrSOLID%elements(i)%gausses(j)%istatus))
            endif
            if( nif(2)>0 ) then
              call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%fstatus,size(fstrSOLID%elements(i)%gausses(j)%fstatus))
            endif
        enddo
      enddo

      if(present(contactNode)) then
        nif(1) = contactNode
        call hecmw_restart_add_int(nif,size(nif))
        do i= 1, size(fstrSOLID%contacts)
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

      integer :: i,j,restrt_step(1),nif(2)
      real(kind=kreal) :: data(2)

      call hecmw_restart_open()

      call hecmw_restart_read_int(restrt_step)
      cstep = restrt_step(1)
      call hecmw_restart_read_real(fstrSOLID%unode)
      call hecmw_restart_read_real(fstrSOLID%QFORCE)

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            call hecmw_restart_read_int(nif)
            call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%strain)
            call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%stress)
            if( nif(1)>0 ) call hecmw_restart_read_int(fstrSOLID%elements(i)%gausses(j)%istatus)
            if( nif(2)>0 ) call hecmw_restart_read_real(fstrSOLID%elements(i)%gausses(j)%fstatus)
        enddo
      enddo

      if(present(contactNode)) then
        call hecmw_restart_read_int(nif)
        contactNode = nif(1)
        do i= 1, size(fstrSOLID%contacts)
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

      integer :: i,j,restrt_step(1),nif(2)
      real(kind=kreal) :: data(2)

      restrt_step(1) = cstep
      call hecmw_restart_add_int(restrt_step,size(restrt_step))
      call hecmw_restart_add_real(fstrSOLID%unode,size(fstrSOLID%unode))
      call hecmw_restart_add_real(fstrSOLID%QFORCE,size(fstrSOLID%QFORCE))

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            nif = 0
            if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) ) nif(1)=size(fstrSOLID%elements(i)%gausses(j)%istatus)
            if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) ) nif(2)=size(fstrSOLID%elements(i)%gausses(j)%fstatus)
            call hecmw_restart_add_int(nif,size(nif))
            call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%strain,size(fstrSOLID%elements(i)%gausses(j)%strain))
            call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%stress,size(fstrSOLID%elements(i)%gausses(j)%stress))
            if( nif(1)>0 ) then
              call hecmw_restart_add_int(fstrSOLID%elements(i)%gausses(j)%istatus,size(fstrSOLID%elements(i)%gausses(j)%istatus))
            endif
            if( nif(2)>0 ) then
              call hecmw_restart_add_real(fstrSOLID%elements(i)%gausses(j)%fstatus,size(fstrSOLID%elements(i)%gausses(j)%fstatus))
            endif
        enddo
      enddo

      if(present(contactNode)) then
        nif(1) = contactNode
        call hecmw_restart_add_int(nif,size(nif))
        do i= 1, size(fstrSOLID%contacts)
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
        call hecmw_restart_add_real(fstrSOLID%elements(i)%equiForces,size(fstrSOLID%elements(i)%equiForces))
      enddo

      call hecmw_restart_write()

      end subroutine fstr_write_restart_dyna_nl

!> read in restart file for linear dynamic analysis
!----------------------------------------------------------------------*
      subroutine fstr_read_restart_dyna_linear(cstep,fstrDYNAMIC)
!----------------------------------------------------------------------*
      implicit none
      integer, intent(out)                  :: cstep       !< current step
      type ( fstr_dynamic), intent(inout)   :: fstrDYNAMIC

      integer :: restrt_step(1)

      call hecmw_restart_open()
      call hecmw_restart_read_int(restrt_step)
      cstep = restrt_step(1)
      if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
        call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
        call hecmw_restart_read_real(fstrDYNAMIC%VEL (:,1))
        call hecmw_restart_read_real(fstrDYNAMIC%ACC (:,1))
      else
        call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
        call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,3))
      endif
      call hecmw_restart_close()
      end subroutine fstr_read_restart_dyna_linear

!> write out restart file for linear dynamic analysis
!----------------------------------------------------------------------*
      subroutine fstr_write_restart_dyna_linear(cstep,fstrDYNAMIC)
!----------------------------------------------------------------------*
      implicit none
      integer, intent(in)                   :: cstep      !< current step
      type ( fstr_dynamic), intent(in)      :: fstrDYNAMIC

      integer :: restrt_step(1)

      restrt_step(1) = cstep
      call hecmw_restart_add_int(restrt_step,size(restrt_step))
      if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
        call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
        call hecmw_restart_add_real(fstrDYNAMIC%VEL (:,1),size(fstrDYNAMIC%VEL (:,1)))
        call hecmw_restart_add_real(fstrDYNAMIC%ACC (:,1),size(fstrDYNAMIC%ACC (:,1)))
      else
        call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
        call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,3),size(fstrDYNAMIC%DISP(:,3)))
      endif
      call hecmw_restart_write()
      end subroutine fstr_write_restart_dyna_linear

end module m_fstr_restart
