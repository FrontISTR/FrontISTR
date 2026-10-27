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
    , dload_grp_name_to_id_ex &
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

end module fstr_setup_solid
