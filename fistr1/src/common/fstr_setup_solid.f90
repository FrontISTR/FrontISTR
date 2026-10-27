!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to append boundary conditions of
!! solid analysis to fstr_solid, independently of the input file format.
!! header_name is the card name shown in error messages. A blank name of
!! amplitude or rotation center means that it is not given.
module fstr_setup_solid
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array
  use fstr_setup_util, only &
    : amp_name_to_id &
    , fstr_setup_util_err_stop
  implicit none

contains

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

end module fstr_setup_solid
