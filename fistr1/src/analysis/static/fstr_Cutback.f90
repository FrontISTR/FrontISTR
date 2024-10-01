!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to deal with cutback

module m_fstr_Cutback
  use m_fstr
  use elementInfo
  use mMechGauss
  use mContactDef

  implicit none

  logical, private :: is_cutback_active = .false.

contains

  logical function fstr_cutback_active()
    fstr_cutback_active = is_cutback_active
  end function

  !> Initializer of cutback variables
  subroutine fstr_cutback_init( hecMESH, fstrSOLID, fstrPARAM )
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_param)         :: fstrPARAM
    type(fstr_solid)         :: fstrSOLID

    integer(kind=kint) :: istep, i, j
    integer(kind=kint) :: ng
    integer(kind=kint) :: ncont, nstate

    do istep=1,fstrSOLID%nstep_tot
      if( fstrSOLID%step_ctrl(istep)%inc_type == stepAutoInc ) is_cutback_active = .true.
    end do

    if( .not. is_cutback_active ) return

    !allocate variables to store analysis state
    !(1)nodal values
    allocate(fstrSOLID%unode_bkup(size(fstrSOLID%unode)))
    allocate(fstrSOLID%QFORCE_bkup(size(fstrSOLID%QFORCE)))
    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      allocate(fstrSOLID%last_temp_bkup(size(fstrSOLID%last_temp)))
    endif

    !(2)elemental values
    allocate(fstrSOLID%elements_bkup(size(fstrSOLID%elements)))
    do i=1,size(fstrSOLID%elements)
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
      allocate( fstrSOLID%elements_bkup(i)%gausses(ng) )
      do j=1,ng
        fstrSOLID%elements_bkup(i)%gausses(j)%pMaterial => fstrSOLID%elements(i)%gausses(j)%pMaterial
        call fstr_init_gauss( fstrSOLID%elements_bkup(i)%gausses(j) )
      end do
      if( associated( fstrSOLID%elements(i)%aux ) ) then
        allocate( fstrSOLID%elements_bkup(i)%aux(3,3) )
      endif
    end do

    !(3)contact values
    ncont = fstrSOLID%n_contacts
    if( associated(fstrSOLID%contacts) ) then
      allocate(fstrSOLID%contacts_bkup(ncont))
      do i=1,ncont
        nstate = size(fstrSOLID%contacts(i)%states)
        allocate(fstrSOLID%contacts_bkup(i)%states(nstate))
      end do
    end if

    !(4)embed values
    ncont = fstrSOLID%n_embeds
    if( associated(fstrSOLID%embeds) ) then
      allocate(fstrSOLID%embeds_bkup(ncont))
      do i=1,ncont
        nstate = size(fstrSOLID%embeds(i)%states)
        allocate(fstrSOLID%embeds_bkup(i)%states(nstate))
      end do
    end if
    
  end subroutine

  !> Finalizer of cutback variables
  subroutine fstr_cutback_finalize( fstrSOLID )
    type(fstr_solid) :: fstrSOLID

    integer(kind=kint) :: i, j
    integer(kind=kint) :: ng
    integer(kind=kint) :: ncont

    if( .not. is_cutback_active ) return

    !deallocate variables to store analysis state
    !(1)nodal values
    if( associated(fstrSOLID%unode_bkup) ) deallocate(fstrSOLID%unode_bkup)
    if( associated(fstrSOLID%QFORCE_bkup) ) deallocate(fstrSOLID%QFORCE_bkup)
    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      if( associated(fstrSOLID%last_temp_bkup) ) deallocate(fstrSOLID%last_temp_bkup)
    endif

    !(2)elemental values
    do i=1,size(fstrSOLID%elements)
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
      do j=1,ng
        call fstr_finalize_gauss( fstrSOLID%elements_bkup(i)%gausses(j) )
      end do
      deallocate( fstrSOLID%elements_bkup(i)%gausses )
      if( associated( fstrSOLID%elements_bkup(i)%aux ) ) then
        deallocate( fstrSOLID%elements_bkup(i)%aux )
      endif
    end do
    deallocate(fstrSOLID%elements_bkup)

    !(3)contact values
    ncont = fstrSOLID%n_contacts
    if( associated(fstrSOLID%contacts) ) then
      do i=1,ncont
        deallocate(fstrSOLID%contacts_bkup(i)%states)
      end do
      deallocate(fstrSOLID%contacts_bkup)
    end if

    !(4)embed values
    ncont = fstrSOLID%n_embeds
    if( associated(fstrSOLID%embeds) ) then
      do i=1,ncont
        deallocate(fstrSOLID%embeds_bkup(i)%states)
      end do
      deallocate(fstrSOLID%embeds_bkup)
    end if

        
  end subroutine

  !> Save analysis status
  subroutine fstr_cutback_save( fstrSOLID, infoCTChange, infoCTChange_bak )
    type(fstr_solid), intent(inout)              :: fstrSOLID
    type(fstr_info_contactChange), intent(inout) :: infoCTChange !< contact change info
    type(fstr_info_contactChange), intent(inout) :: infoCTChange_bak !< contact change info

    integer(kind=kint) :: i, j
    integer(kind=kint) :: ng
    integer(kind=kint) :: ncont, nstate

    if( .not. is_cutback_active ) return

    !(1)nodal values
    do i=1,size(fstrSOLID%unode)
      fstrSOLID%unode_bkup(i) = fstrSOLID%unode(i)
    end do
    do i=1,size(fstrSOLID%QFORCE)
      fstrSOLID%QFORCE_bkup(i) = fstrSOLID%QFORCE(i)
    end do
    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      do i=1,size(fstrSOLID%last_temp)
        fstrSOLID%last_temp_bkup(i) = fstrSOLID%last_temp(i)
      end do
    endif

    !(2)elemental values
    do i=1,size(fstrSOLID%elements)
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
      do j=1,ng
        call fstr_copy_gauss( fstrSOLID%elements(i)%gausses(j), fstrSOLID%elements_bkup(i)%gausses(j) )
      end do
      if( associated( fstrSOLID%elements(i)%aux ) ) then
        fstrSOLID%elements_bkup(i)%aux(:,:) = fstrSOLID%elements(i)%aux(:,:)
      endif
    end do

    !(3)contact values
    ncont = fstrSOLID%n_contacts
    if( associated(fstrSOLID%contacts) ) then
      do i=1,ncont
        nstate = size(fstrSOLID%contacts(i)%states)
        do j=1,nstate
          call contact_state_copy(fstrSOLID%contacts(i)%states(j), fstrSOLID%contacts_bkup(i)%states(j))
        enddo
      end do
      infoCTChange_bak = infoCTChange
    end if

    !(4)embed values
    ncont = fstrSOLID%n_embeds
    if( associated(fstrSOLID%embeds) ) then
      do i=1,ncont
        nstate = size(fstrSOLID%embeds(i)%states)
        do j=1,nstate
          call contact_state_copy(fstrSOLID%embeds(i)%states(j), fstrSOLID%embeds_bkup(i)%states(j))
        enddo
      end do
      infoCTChange_bak = infoCTChange
    end if
  end subroutine

  !> Load analysis status
  subroutine fstr_cutback_load( fstrSOLID, infoCTChange, infoCTChange_bak )
    type(fstr_solid), intent(inout)              :: fstrSOLID
    type(fstr_info_contactChange), intent(inout) :: infoCTChange !< contact change info
    type(fstr_info_contactChange), intent(inout) :: infoCTChange_bak !< contact change info

    integer(kind=kint) :: i, j
    integer(kind=kint) :: ng
    integer(kind=kint) :: ncont, nstate

    if( .not. is_cutback_active ) return

    !(1)nodal values
    do i=1,size(fstrSOLID%unode)
      fstrSOLID%unode(i) = fstrSOLID%unode_bkup(i)
    end do
    do i=1,size(fstrSOLID%QFORCE)
      fstrSOLID%QFORCE(i) = fstrSOLID%QFORCE_bkup(i)
    end do
    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      do i=1,size(fstrSOLID%last_temp)
        fstrSOLID%last_temp(i) = fstrSOLID%last_temp_bkup(i)
      end do
    endif

    !(2)elemental values
    do i=1,size(fstrSOLID%elements)
      if (hecmw_is_etype_link( fstrSOLID%elements(i)%etype )) cycle
      if (hecmw_is_etype_patch( fstrSOLID%elements(i)%etype )) cycle
      ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
      do j=1,ng
        call fstr_copy_gauss( fstrSOLID%elements_bkup(i)%gausses(j), fstrSOLID%elements(i)%gausses(j) )
      end do
      if( associated( fstrSOLID%elements(i)%aux ) ) then
        fstrSOLID%elements(i)%aux(:,:) = fstrSOLID%elements_bkup(i)%aux(:,:)
      endif
    end do

    !(3)contact values
    ncont = fstrSOLID%n_contacts
    if( associated(fstrSOLID%contacts) ) then
      do i=1,ncont
        nstate = size(fstrSOLID%contacts(i)%states)
        do j=1,nstate
          call contact_state_copy(fstrSOLID%contacts_bkup(i)%states(j), fstrSOLID%contacts(i)%states(j))
        enddo
      end do
      infoCTChange = infoCTChange_bak
    end if

    !(4)embed values
    ncont = fstrSOLID%n_embeds
    if( associated(fstrSOLID%embeds) ) then
      do i=1,ncont
        nstate = size(fstrSOLID%embeds(i)%states)
        do j=1,nstate
          call contact_state_copy(fstrSOLID%embeds_bkup(i)%states(j), fstrSOLID%embeds(i)%states(j))
        enddo
      end do
      infoCTChange = infoCTChange_bak
    end if

  end subroutine

end module m_fstr_Cutback
