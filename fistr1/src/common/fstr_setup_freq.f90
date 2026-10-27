!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to append loads of frequency
!! response analysis to fstr_freqanalysis, independently of the input file format.
!! header_name is the card name shown in error messages. loadcase is
!! kFLOADCASE_RE (real part) or kFLOADCASE_IM (imaginary part).
module fstr_setup_freq
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array
  use fstr_setup_util, only &
    : nodesurf_grp_name_to_id_ex &
    , fstr_setup_util_err_stop
  implicit none

contains

  !-----------------------------------------------------------------------------!
  !> Append !FLOAD                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_freq_append_fload( hecMESH, fstrFREQ, header_name, gid, loadcase, &
      n, grp_names, dofs, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_freqanalysis)                   :: fstrFREQ
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: gid
    integer(kind=kint), intent(in)            :: loadcase
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dofs(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: i, old_size, new_size

    if( n == 0 ) return
    old_size = fstrFREQ%FLOAD_ngrp_tot
    new_size = old_size + n

    !expand data array
    fstrFREQ%FLOAD_ngrp_tot = new_size
    call hecmw_expand_integer_array( fstrFREQ%FLOAD_ngrp_GRPID, old_size, new_size )
    call hecmw_expand_integer_array( fstrFREQ%FLOAD_ngrp_ID,    old_size, new_size )
    call hecmw_expand_integer_array( fstrFREQ%FLOAD_ngrp_TYPE,  old_size, new_size )
    call hecmw_expand_integer_array( fstrFREQ%FLOAD_ngrp_DOF,   old_size, new_size )
    call hecmw_expand_real_array   ( fstrFREQ%FLOAD_ngrp_valre, old_size, new_size )
    call hecmw_expand_real_array   ( fstrFREQ%FLOAD_ngrp_valim, old_size, new_size )

    !fill bc data
    if(loadcase == kFLOADCASE_RE) then
      do i = 1, n
        fstrFREQ%FLOAD_ngrp_DOF(old_size+i) = dofs(i)
        fstrFREQ%FLOAD_ngrp_valre(old_size+i) = vals(i)
      enddo
    else if(loadcase == kFLOADCASE_IM) then
      do i = 1, n
        fstrFREQ%FLOAD_ngrp_DOF(old_size+i) = dofs(i)
        fstrFREQ%FLOAD_ngrp_valim(old_size+i) = vals(i)
      enddo
    end if
    fstrFREQ%FLOAD_ngrp_GRPID(old_size+1:new_size) = gid
    call nodesurf_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, &
         fstrFREQ%FLOAD_ngrp_ID(old_size+1:), fstrFREQ%FLOAD_ngrp_TYPE(old_size+1:))

  end subroutine fstr_freq_append_fload

  !-----------------------------------------------------------------------------!
  !> Append !FACCELERATION                                                      !
  !-----------------------------------------------------------------------------!
  subroutine fstr_freq_append_facceleration( hecMESH, fstrFREQ, header_name, loadcase, &
      n, grp_names, dof_ids, dof_ide, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_freqanalysis)                   :: fstrFREQ
    character(len=*), intent(in)              :: header_name
    integer(kind=kint), intent(in)            :: loadcase
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: dof_ids(:)
    integer(kind=kint), intent(in)            :: dof_ide(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: i, old_size, new_size

    if( n == 0 ) return
    old_size = fstrFREQ%FACC_ngrp_tot
    new_size = old_size + n

    !the expanded part is zero, so the part not given by loadcase stays zero
    fstrFREQ%FACC_ngrp_tot = new_size
    call hecmw_expand_integer_array( fstrFREQ%FACC_ngrp_ID,    old_size, new_size )
    call hecmw_expand_integer_array( fstrFREQ%FACC_ngrp_type,  old_size, new_size )
    call hecmw_expand_real_array   ( fstrFREQ%FACC_ngrp_valre, old_size, new_size )
    call hecmw_expand_real_array   ( fstrFREQ%FACC_ngrp_valim, old_size, new_size )

    !node ID or node group name -> node group ID
    call node_grp_name_to_id_ex( hecMESH, header_name, n, grp_names, fstrFREQ%FACC_ngrp_ID(old_size+1:) )

    do i = 1, n
      !check the range of the DOFs as !ACCELERATION does
      if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < dof_ids(i)).or.(6 < dof_ide(i)) ) then
        call fstr_setup_util_err_stop( 'fstr control file error : '//header_name// &
          ' : range of dof_ids and dof_ide is from 1 to 6' )
      end if
      !encode the DOF range as !ACCELERATION does
      fstrFREQ%FACC_ngrp_type(old_size+i) = 10 * dof_ids(i) + dof_ide(i)
      !store the amplitude in the part selected by loadcase
      if( loadcase == kFLOADCASE_RE ) then
        fstrFREQ%FACC_ngrp_valre(old_size+i) = vals(i)
      else
        fstrFREQ%FACC_ngrp_valim(old_size+i) = vals(i)
      end if
    end do

  end subroutine fstr_freq_append_facceleration

end module fstr_setup_freq
