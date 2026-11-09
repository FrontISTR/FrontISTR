!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to read the cards of frequency
!! response analysis into fstr_freqanalysis. Loads are appended independently
!! of the input file format.
!! header_name is the card name shown in error messages. loadcase is
!! kFLOADCASE_RE (real part) or kFLOADCASE_IM (imaginary part).
module m_fstr_freq_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_ctrl_err_stop &
    , nodesurf_grp_name_to_id_ex &
    , fstr_setup_util_err_stop
  use fstr_ctrl_util_f
  implicit none

contains

  subroutine fstr_freqanalysis_finalize( f )
    implicit none
    type( fstr_freqanalysis ), intent(inout) :: f
    if(associated(f%FLOAD_ngrp_GRPID)) deallocate(f%FLOAD_ngrp_GRPID)
    if(associated(f%FLOAD_ngrp_ID)) deallocate(f%FLOAD_ngrp_ID)
    if(associated(f%FLOAD_ngrp_TYPE)) deallocate(f%FLOAD_ngrp_TYPE)
    if(associated(f%FLOAD_ngrp_DOF)) deallocate(f%FLOAD_ngrp_DOF)
    if(associated(f%FLOAD_ngrp_valre)) deallocate(f%FLOAD_ngrp_valre)
    if(associated(f%FLOAD_ngrp_valim)) deallocate(f%FLOAD_ngrp_valim)
    if(associated(f%FACC_ngrp_ID)) deallocate(f%FACC_ngrp_ID)
    if(associated(f%FACC_ngrp_type)) deallocate(f%FACC_ngrp_type)
    if(associated(f%FACC_ngrp_valre)) deallocate(f%FACC_ngrp_valre)
    if(associated(f%FACC_ngrp_valim)) deallocate(f%FACC_ngrp_valim)
  end subroutine fstr_freqanalysis_finalize

  subroutine fstr_freqanalysis_data_finalize( f )
    implicit none
    type( fstr_freqanalysis_data ), intent(inout) :: f
    if(associated(f%eigOmega)) deallocate(f%eigOmega)
    if(associated(f%eigVector)) deallocate(f%eigVector)
  end subroutine fstr_freqanalysis_data_finalize

  !-----------------------------------------------------------------------------!
  !> Read in !FLOAD                                                             !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_FLOAD( ctrl, counter, P )
    implicit none
  !---- args
    integer(kind=kint)   :: ctrl
    integer(kind=kint)   :: counter
    type(fstr_param_pack) :: P
  !---- vals
    integer(kind=kint)                  :: rcode
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal), pointer           :: val_ptr(:)
    integer(kind=kint), pointer        :: id_ptr(:)
    integer(kind=kint)                  :: n
    integer(kind=kint)                  :: gid, loadcase
  !---- body

    if( P%SOLID%file_type /= kbcfFSTR) return

    !read grpid
    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',  0, 'I', gid )
    !read loadcase (real=1:default, img=2)
    loadcase = kFLOADCASE_RE
    rcode = fstr_ctrl_get_param_ex( ctrl, 'LOAD CASE ', '# ', 0, 'I', loadcase)
    !write(*,*) "loadcase=", loadcase
    !pause

    !read the num of dataline
    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    !fill bc data
    allocate( grp_id_name(n) )
    allocate( id_ptr(n) )
    allocate( val_ptr(n) )
    id_ptr  = 0
    val_ptr = 0.0d0
    rcode = fstr_ctrl_get_FLOAD( ctrl, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr)
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if(loadcase /= kFLOADCASE_RE .and. loadcase /= kFLOADCASE_IM) then
      !error
      write(*,*)    "Error this load set is not defined!"
      write(ilog,*) "Error this load set is not defined!"
      stop HECMW_EXIT_MODEL
    end if
    call fstr_freq_append_fload( P%MESH, P%FREQ, '!FLOAD', gid, loadcase, &
      n, grp_id_name, id_ptr, val_ptr )

    deallocate( grp_id_name )
    deallocate( id_ptr )
    deallocate( val_ptr )
    nullify( grp_id_name )
    nullify( id_ptr )
    nullify( val_ptr )
    return

    contains

    function fstr_ctrl_get_FLOAD(ctrl, node_id, node_id_len, dof_id, value)
      integer(kind=kint)                    :: ctrl
      character(len=HECMW_NAME_LEN)        :: node_id(:)  !Node group name
      integer(kind=kint), pointer          :: dof_id(:)
      integer(kind=kint)                    :: node_id_len
      real(kind=kreal), pointer             :: value(:)
      integer(kind=kint)                    :: fstr_ctrl_get_FLOAD !return value
      character(len=HECMW_NAME_LEN)        :: data_fmt, ss

      write(ss,*) node_id_len
      write(data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IR '

      fstr_ctrl_get_FLOAD = fstr_ctrl_get_data_array_ex(ctrl, data_fmt, node_id, dof_id, value)
    end function

  end subroutine

  !-----------------------------------------------------------------------------!
  !> Read in !FACCELERATION                                                     !
  !-----------------------------------------------------------------------------!
  !> Data line: NODE_ID, DOF_idS, DOF_idE, Value (same as !ACCELERATION).
  !> Value is the acceleration amplitude of DOF_idS..DOF_idE, stored as the real
  !> part (LOAD CASE=1, default) or the imaginary part (LOAD CASE=2).
  !> counter is not used; it is kept for the same interface as the other readers.
  subroutine fstr_setup_FACCELERATION( ctrl, counter, P )
    implicit none
  !---- args
    integer(kind=kint)    :: ctrl
    integer(kind=kint)    :: counter
    type(fstr_param_pack) :: P
  !---- vals
    integer(kind=kint)                 :: rcode
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint), pointer        :: dof_ids(:)
    integer(kind=kint), pointer        :: dof_ide(:)
    real(kind=kreal), pointer          :: val_ptr(:)
    integer(kind=kint)                 :: n
    integer(kind=kint)                 :: loadcase
  !---- body

    !only the FrontISTR format is supported
    if( P%SOLID%file_type /= kbcfFSTR) return

    !read loadcase (real=1:default, img=2)
    loadcase = kFLOADCASE_RE
    rcode = fstr_ctrl_get_param_ex( ctrl, 'LOAD CASE ', '# ', 0, 'I', loadcase)
    if( loadcase /= kFLOADCASE_RE .and. loadcase /= kFLOADCASE_IM ) then
      write(*,*)    'fstr control file error : !FACCELERATION : LOAD CASE must be 1 or 2'
      write(ilog,*) 'fstr control file error : !FACCELERATION : LOAD CASE must be 1 or 2'
      stop HECMW_EXIT_INPUT
    end if

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    !read the data lines into temporary arrays
    allocate( grp_id_name(n) )
    allocate( dof_ids(n) )
    allocate( dof_ide(n) )
    allocate( val_ptr(n) )
    val_ptr = 0.0d0
    rcode = fstr_ctrl_get_FACCELERATION( ctrl, grp_id_name, HECMW_NAME_LEN, dof_ids, dof_ide, val_ptr )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    !append the data lines of this header to the entries read so far
    call fstr_freq_append_facceleration( P%MESH, P%FREQ, '!FACCELERATION', loadcase, &
      n, grp_id_name, dof_ids, dof_ide, val_ptr )

    deallocate( grp_id_name )
    deallocate( dof_ids )
    deallocate( dof_ide )
    deallocate( val_ptr )
    return

    contains

    !read all data lines in the format of NODE_ID, DOF_idS, DOF_idE, Value
    function fstr_ctrl_get_FACCELERATION(ctrl, node_id, node_id_len, dof_ids, dof_ide, value)
      integer(kind=kint)            :: ctrl
      character(len=HECMW_NAME_LEN) :: node_id(:)
      integer(kind=kint)            :: node_id_len
      integer(kind=kint), pointer   :: dof_ids(:)
      integer(kind=kint), pointer   :: dof_ide(:)
      real(kind=kreal), pointer     :: value(:)
      integer(kind=kint)            :: fstr_ctrl_get_FACCELERATION
      character(len=HECMW_NAME_LEN) :: data_fmt, ss

      !format: string of node_id_len characters, two integers and a real
      write(ss,*) node_id_len
      write(data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IIr '

      fstr_ctrl_get_FACCELERATION = fstr_ctrl_get_data_array_ex(ctrl, data_fmt, node_id, dof_ids, dof_ide, value)
    end function

  end subroutine

  !-----------------------------------------------------------------------------!
  !> Read in !EIGENREAD                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_eigenread( ctrl, counter, P )
  !---- args
    integer(kind=kint)    :: ctrl
    integer(kind=kint)    :: counter
    type(fstr_param_pack) :: P
  !---- vals
    integer(kind=kint)                :: filename_len
    character(len=HECMW_NAME_LEN) :: datafmt, ss
  !---- body

    filename_len = HECMW_FILENAME_LEN
    write(ss,*) filename_len
    write(datafmt, '(a,a,a)') 'F', trim(adjustl(ss)), ' '

    if( fstr_ctrl_get_data_ex( ctrl, 1, datafmt, P%FREQ%eigenlog_filename ) /= 0) return
    if( fstr_ctrl_get_data_ex( ctrl, 2, 'ii ', P%FREQ%start_mode, P%FREQ%end_mode ) /= 0) return

    return

  end subroutine

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

end module m_fstr_freq_setup
