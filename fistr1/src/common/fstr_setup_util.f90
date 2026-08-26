!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains auxiliary functions in calculation setup

module fstr_setup_util
  use m_fstr
  use hecmw
  use fstr_ctrl_util_f

contains
  !------------------------------------------------------------------------------

  subroutine fstr_strupr( s )
    implicit none
    character(*) :: s
    integer :: i, n, a

    n = len_trim(s)
    do i = 1, n
      a = iachar(s(i:i))
      if( a >= iachar('a') .and. a <= iachar('z')) then
        s(i:i) = achar(a - 32)
      end if
    end do
  end subroutine fstr_strupr

  !------------------------------------------------------------------------------

  subroutine fstr_ctrl_err_stop
    implicit none
    character(len=256) :: msg

    call fstr_ctrl_get_err_msg( msg, 256 )
    write(*,*) msg
    write(imsg,*) msg
    call hecmw_abort( hecmw_comm_get_comm())
  end subroutine fstr_ctrl_err_stop

  !------------------------------------------------------------------------------

  subroutine fstr_setup_util_err_stop( msg )
    implicit none
    character(*) :: msg

    write(*,*) msg
    write(imsg,*) msg
    call hecmw_abort( hecmw_comm_get_comm())
  end subroutine fstr_setup_util_err_stop

  subroutine append_node_grp_from_surf_grp( hecMESH, sgrp_id, ngrp_id )
    use hecmw_setup_util, only : append_new_group
    use hecmw_array_util, only : hecmw_qsort_int_array, hecmw_uniq_int_array
    implicit none
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    integer(kind=kint), intent(in) :: sgrp_id
    integer(kind=kint), intent(out) :: ngrp_id
    integer(kind=kint) :: is, ie, nnode, i, ic, isurf, ic_type, stype, nn, j0, j, ndup, new_nnode
    integer(kind=kint) :: snode(20)
    integer(kind=kint), allocatable :: node(:)
    character(len=HECMW_NAME_LEN) :: grp_name
    is= hecMESH%surf_group%grp_index(sgrp_id-1) + 1
    ie= hecMESH%surf_group%grp_index(sgrp_id  )
    ! count num of nodes on surface incl duplication
    nnode = 0
    do i=is,ie
      ic   = hecMESH%surf_group%grp_item(2*i-1)
      isurf = hecMESH%surf_group%grp_item(2*i)
      ic_type = hecMESH%elem_type(ic)
      call getSubFace( ic_type, isurf, stype, snode )
      nnode = nnode + getNumberOfNodes( stype )
    enddo
    ! extract nodes on surface incl duplication
    allocate( node(nnode) )
    nnode = 0
    do i=is,ie
      ic   = hecMESH%surf_group%grp_item(2*i-1)
      isurf = hecMESH%surf_group%grp_item(2*i)
      ic_type = hecMESH%elem_type(ic)
      call getSubFace( ic_type, isurf, stype, snode )
      nn = getNumberOfNodes( stype )
      j0 = hecMESH%elem_node_index(ic-1)
      do j=1,nn
        node(nnode+j) = hecMESH%elem_node_item(j0+snode(j))
      enddo
      nnode = nnode + nn
    enddo
    ! sort and uniq node list
    call hecmw_qsort_int_array(node, 1, nnode)
    call hecmw_uniq_int_array(node, 1, nnode, ndup)
    new_nnode = nnode - ndup
    ! append node group
    write( grp_name, '(a,a)') 'FSTR_S2N_',trim(hecMESH%surf_group%grp_name(sgrp_id))
    call append_new_group(hecMESH, 'node_grp', grp_name, new_nnode, node, ngrp_id)
    deallocate(node)
  end subroutine append_node_grp_from_surf_grp

  subroutine append_intersection_node_grp( hecMESH, ngrp_id1, ngrp_id2 )
    use hecmw_setup_util, only : append_new_group
    use hecmw_array_util, only : hecmw_qsort_int_array
    implicit none
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    integer(kind=kint), intent(in) :: ngrp_id1, ngrp_id2
    integer(kind=kint) :: nnode1, nnode2, nnode, is, i, nisect, ngrp_id
    integer(kind=kint), allocatable :: node(:), isect(:)
    character(len=HECMW_NAME_LEN) :: grp_name
    nnode1 = hecMESH%node_group%grp_index(ngrp_id1) - hecMESH%node_group%grp_index(ngrp_id1-1)
    nnode2 = hecMESH%node_group%grp_index(ngrp_id2) - hecMESH%node_group%grp_index(ngrp_id2-1)
    nnode = nnode1 + nnode2
    allocate( node(nnode) )
    is= hecMESH%node_group%grp_index(ngrp_id1-1)
    do i=1,nnode1
      node(i) = hecMESH%node_group%grp_item(is+i)
    enddo
    is= hecMESH%node_group%grp_index(ngrp_id2-1)
    do i=1,nnode2
      node(nnode1+i) = hecMESH%node_group%grp_item(is+i)
    enddo
    call hecmw_qsort_int_array(node, 1, nnode)
    allocate( isect(nnode) )
    nisect = 0
    do i=1,nnode-1
      if( node(i) == node(i+1) ) then
        nisect = nisect + 1
        isect(nisect) = node(i)
      endif
    enddo
    write( grp_name, '(a,a,a,a)') &
         'FSTR_ISCT_',trim(hecMESH%node_group%grp_name(ngrp_id1)),'_AND_',trim(hecMESH%node_group%grp_name(ngrp_id2))
    call append_new_group(hecMESH, 'node_grp', grp_name, nisect, isect, ngrp_id)
    deallocate(node)
    deallocate(isect)
  end subroutine append_intersection_node_grp

  !------------------------------------------------------------------------------
  ! JP-3
  ! JP-4
  ! type_name : 'node', 'element'
  ! name : group name
  ! local_id : local id (set only when return value > 0)
  ! return : -1 if name is not a number
  !          0 if name is a number and a node with ID=name is not in myrank
  !          >0 if name is a number and a node with ID=name is in myrank

  function get_local_member_index( hecMESH, type_name, name, local_id )
    use hecmw_setup_util, only: hecmw_str2index
    implicit none
    integer(kind=kint) :: get_local_member_index
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: type_name
    character(len=*) :: name
    integer(kind=kint) :: local_id
    integer(kind=kint) :: i, n, no, fg
    integer(kind=kint),pointer :: global_item(:)

    if( .not. hecmw_str2index(name, no) ) then
      get_local_member_index = -1
      return
    end if

    if( type_name == 'node' ) then
      fg = 1
      n  =  hecMESH%n_node
      global_item => hecMESH%global_node_ID
    else if( type_name == 'element' ) then
      fg = 2
      n  =  hecMESH%n_elem
      global_item => hecMESH%global_elem_ID
    else
      stop 'assert in get_local_member_index: unknown type_name'
    end if

    do i = 1, n
      if( no == global_item(i)) then
        local_id = i
        get_local_member_index = local_id
        return
      end if
    end do
    local_id = 0
    get_local_member_index = 0
    return
  end function get_local_member_index

  !-----------------------------------------------------------------------------!
  !

  function get_sorted_local_member_index( hecMESH, hecPARAM, type_name, name, local_id )
    use hecmw_setup_util, only: hecmw_str2index
    use hecmw_array_util, only: hecmw_bsearch_int_array
    implicit none
    integer(kind=kint) :: get_sorted_local_member_index
    type (hecmwST_local_mesh),target :: hecMESH
    type(fstr_param), target         :: hecPARAM
    character(len=*) :: type_name
    character(len=*) :: name
    integer(kind=kint) :: local_id, idx
    integer(kind=kint) :: n, no, fg

    if( .not. hecmw_str2index(name, no) ) then
      get_sorted_local_member_index = -1
      return
    end if

    if( type_name == 'node' ) then
      fg = 1
      n  =  hecMESH%nn_internal
      !   item => hecMESH%global_node_ID
      ! else if( type_name == 'element' ) then
      !   fg = 2
      !   n  =  hecMESH%n_elem
      !   item => hecMESH%global_elem_ID
    else
      stop 'assert in get_sorted_local_member_index: unknown type_name'
    end if

    call hecmw_bsearch_int_array(hecPARAM%global_local_ID(1,:), 1, n, no, idx)
    if(idx > 0)then
      get_sorted_local_member_index = hecPARAM%global_local_ID(2,idx)
      local_id = get_sorted_local_member_index
      return
    endif

    get_sorted_local_member_index = 0
    return
  end function get_sorted_local_member_index
  !-----------------------------------------------------------------------------!

  !Find node/surf group from name or nodeid

  subroutine nodesurf_grp_name_to_id_ex(hecMESH, header_name, n, grp_id_name, grp_ID, grp_TYPE)
    use hecmw_setup_util, only : append_single_group, hecmw_str2index, hecmw_streqr
    use m_fstr
    implicit none
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*)                 :: header_name
    integer(kind=kint) :: n
    character(len=HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint) :: grp_ID(:)
    integer(kind=kint) :: grp_TYPE(:)

    integer(kind=kint) :: i, id
    integer(kind=kint) :: no, no_count, exist_n
    integer(kind=kint),pointer :: no_list(:)
    character(HECMW_NAME_LEN) :: name
    character(len=256) :: msg

    allocate( no_list( n ))
    no_count = 0
    do i = 1, n
      if( hecmw_str2index( grp_id_name(i), no )) then
        no_count = no_count + 1
        no_list(no_count) = no
        grp_ID(i)  = hecMESH%node_group%n_grp + no_count
        grp_TYPE(i) = kFLOADTYPE_NODE
      else
        !Find node group
        grp_ID(i) = -1
        do id = 1, hecMESH%node_group%n_grp
          if (hecmw_streqr(hecMESH%node_group%grp_name(id), grp_id_name(i))) then
            grp_ID(i)   = id
            grp_TYPE(i) = kFLOADTYPE_NODE
            exit
          end if
        end do
        !Find surf group
        if (grp_ID(i) == -1) then
          do id = 1, hecMESH%surf_group%n_grp
            if (hecmw_streqr(hecMESH%surf_group%grp_name(id), grp_id_name(i))) then
              grp_ID(i)   = id
              grp_TYPE(i) = kFLOADTYPE_SURF
              exit
            end if
          end do
        end if

        !not found => exit
        if( grp_ID(i) == -1 ) then
          write(msg,*) '### Error: ', header_name,' : Node group "',grp_id_name(i),'" does not exist.'
          call fstr_setup_util_err_stop(msg)
        end if
      end if
    end do
    if( no_count > 0 ) then
      name = 'node_grp'
      exist_n = append_single_group( hecMESH, name, no_count, no_list )
    end if

    deallocate( no_list )

  end subroutine nodesurf_grp_name_to_id_ex

  !------------------------------------------------------------------------------

  subroutine dload_grp_name_to_id_ex( hecMESH, n, grp_id_name, fg_surface, grp_ID )
    use hecmw_setup_util, only : append_single_group, hecmw_str2index, hecmw_streqr
    implicit none
    type (hecmwST_local_mesh),target :: hecMESH
    integer(kind=kint) :: n
    integer(kind=kint),save :: casha = 1, cashb = 1
    character(HECMW_NAME_LEN) :: grp_id_name(:)
    logical :: fg_surface(:)
    integer(kind=kint) :: grp_ID(:)
    integer(kind=kint) :: i, id
    integer(kind=kint) :: no, no_count, exist_n
    integer(kind=kint),pointer :: no_list(:)
    character(HECMW_NAME_LEN) :: name
    character(len=256) :: msg

    allocate( no_list( n ))
    no_count = 0
    do i = 1, n
      if( fg_surface(i) ) then
        grp_ID(i) = -1
        if(casha < hecMESH%surf_group%n_grp)then
          if(hecmw_streqr(hecMESH%surf_group%grp_name(casha), grp_id_name(i))) then
            grp_ID(i) = casha
            casha = casha + 1
            cycle
          end if
        endif
        do id = 1, hecMESH%surf_group%n_grp
          if(hecmw_streqr(hecMESH%surf_group%grp_name(id), grp_id_name(i))) then
            grp_ID(i) = id
            casha = id + 1
            exit
          end if
        end do
        if( grp_ID(i) == -1 ) then
          write(msg,*) '### Error: !DLOAD : Surface group "',&
            grp_id_name(i),'" does not exist.'
          call fstr_setup_util_err_stop(msg)
        end if
      else
        if( hecmw_str2index( grp_id_name(i), no )) then
          no_count = no_count + 1
          no_list(no_count) = no
          grp_ID(i) = hecMESH%elem_group%n_grp + no_count
        else
          grp_ID(i) = -1
          if(cashb < hecMESH%surf_group%n_grp)then
            if(hecmw_streqr(hecMESH%surf_group%grp_name(cashb), grp_id_name(i))) then
              grp_ID(i) = cashb
              cashb = cashb + 1
              cycle
            end if
          endif
          do id = 1, hecMESH%elem_group%n_grp
            if(hecmw_streqr(hecMESH%elem_group%grp_name(id), grp_id_name(i))) then
              grp_ID(i) = id
              cashb = cashb + 1
              exit
            end if
          end do
          if( grp_ID(i) == -1 ) then
            write(msg,*) '### Error: !DLOAD : Element group "',&
              grp_id_name(i),'" does not exist.'
            call fstr_setup_util_err_stop(msg)
          end if
        end if
      end if
    end do

    if( no_count > 0 ) then
      name = 'elem_grp'
      exist_n = append_single_group( hecMESH, name, no_count, no_list )
      !   if( exist_n < no_count ) then
      !     write(*,*) '### Warning: !DLOAD : following elements are not exist'
      !     if( hecMESH%my_rank == 0 ) then
      !       write(imsg,*) '### Warning: !DLOAD : following elements are not exist'
      !     end if
      !     do i=1, no_count
      !       if( no_list(i)<0 ) then
      !         write(*,*) -no_list(i)
      !         if( hecMESH%my_rank == 0 ) then
      !           write(imsg,*) -no_list(i)
      !         endif
      !       end if
      !     end do
      !   end if
    end if

    deallocate( no_list )
  end subroutine dload_grp_name_to_id_ex

  !------------------------------------------------------------------------------
  ! JP-7

  !> Append new amplitude table at the end of existing amplitude tables
  subroutine append_new_amplitude( amp, name, type_def, type_time, type_val, np, val, table )
    use hecmw_setup_util, only &
      : hecmw_expand_index_array &
      , hecmw_expand_char_array &
      , hecmw_expand_integer_array &
      , hecmw_expand_real_array &
      , hecmw_streqr
    type( hecmwST_amplitude ),     intent(inout) :: amp       !< amplitude table structure
    character(len=HECMW_NAME_LEN), intent(in)    :: name      !< name of amplitude table
    integer(kind=kint),            intent(in)    :: type_def  !< HECMW_AMP_TYPEDEF_TABULAR
    integer(kind=kint),            intent(in)    :: type_time !< HECMW_AMP_TYPETIME_STEP
    integer(kind=kint),            intent(in)    :: type_val  !< HECMW_AMP_TYPEVAL_{RELATIVE|ABSOLUTE}
    integer(kind=kint),            intent(in)    :: np        !< number of table items
    real(kind=kreal),              intent(in)    :: val(:)    !< values of the table
    real(kind=kreal),              intent(in)    :: table(:)  !< time points of the table

    ! type(fstr_str_arr) :: amp_name
    integer(kind=kint) :: n_amp, new_size, old_size, i

    do i=1,amp%n_amp
      if( hecmw_streqr(amp%amp_name(i), name) ) then
        write(*,*) 'Error: AMPLITUDE with NAME=',trim(name),' already exists'
        call fstr_ctrl_err_stop
      endif
    enddo

    n_amp = amp%n_amp
    new_size = n_amp+1
    amp%n_amp = new_size
    call hecmw_expand_index_array( amp%amp_index, n_amp+1, new_size+1 )
    ! amp_name%s => amp%amp_name
    ! call hecmw_expand_name_array( amp_name, n_amp, new_size )
    ! amp%amp_name => amp_name%s
    call hecmw_expand_char_array( amp%amp_name, n_amp, new_size )
    call hecmw_expand_integer_array( amp%amp_type_definition, n_amp, new_size )
    call hecmw_expand_integer_array( amp%amp_type_time, n_amp, new_size )
    call hecmw_expand_integer_array( amp%amp_type_value, n_amp, new_size )
    old_size = amp%amp_index( n_amp )
    new_size = old_size+np
    call hecmw_expand_real_array( amp%amp_val, old_size, new_size )
    call hecmw_expand_real_array( amp%amp_table, old_size, new_size )

    amp%amp_index(amp%n_amp) = amp%amp_index(amp%n_amp-1)+np
    amp%amp_name(amp%n_amp) = name
    amp%amp_type_definition(amp%n_amp) = type_def
    amp%amp_type_time(amp%n_amp) = type_time
    amp%amp_type_value(amp%n_amp) = type_val
    do i=1,np
      amp%amp_val(old_size+i) = val(i)
      amp%amp_table(old_size+i) = table(i)
    enddo
  end subroutine append_new_amplitude


  subroutine amp_name_to_id( hecMESH, header_name, aname, id )
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    character(len=*)          :: header_name
    character(len=HECMW_NAME_LEN)::aname
    integer(kind=kint) :: id
    character(len=256) :: msg

    id = 0
    if(  aname .eq. ' ' )  return
    call get_amp_id( hecMESH, aname, id )
    if( id == 0 ) then
      write(msg,*) '### Error: ', header_name,' : Amplitude group "',&
        aname,'" does not exist.'
      call fstr_setup_util_err_stop(msg)
    end if
  end subroutine amp_name_to_id


  !GET AMPLITUDE INDEX

  subroutine get_amp_id( hecMESH, aname, id )
    use hecmw_setup_util, only: hecmw_streqr
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    character(len=HECMW_NAME_LEN)::aname
    integer(kind=kint) :: id

    integer(kind=kint) :: i

    id = 0
    if(  aname .eq. ' ' )  return

    do i = 1, hecMESH%amp%n_amp
      if( hecmw_streqr(hecMESH%amp%amp_name(i), aname)) then
        id = i
        return
      end if
    end do
  end subroutine get_amp_id

  !------------------------------------------------------------------------------

  subroutine reallocate_integer( array, n )
    implicit none
    integer(kind=kint),pointer :: array(:)
    integer(kind=kint) :: n;

    if( associated( array )) deallocate(array)
    allocate( array(n));
  end subroutine reallocate_integer

  subroutine reallocate_real( array, n )
    implicit none
    real(kind=kreal),pointer :: array(:)
    integer(kind=kint) :: n;

    if( associated( array )) deallocate(array)
    allocate( array(n));
  end subroutine reallocate_real

  !-----------------------------------------------------------------------------!
  ! FSTR_SETUP_VISUALIZE                                                        !
  ! 1) Seeking header to 'WRITE'                                                !
  ! 2) If parameter 'VISUAL' exists, then 'hecmw_vis.ini' is opened.            !
  ! 3) All following lines under the header are written to the opened file      !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_visualize( ctrl, hecMESH )
    implicit none
    integer(kind=kint) :: ctrl
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: rcode
    character(HECMW_FILENAME_LEN) :: vis_filename = 'hecmw_vis.ini'
    logical :: is_exit

    rcode = fstr_ctrl_seek_header( ctrl, '!VISUAL ' )
    if(rcode == 0) return

    if(hecMESH%my_rank == 0)then
      call fstr_setup_visualize_main( ctrl, vis_filename )
    endif

    call hecmw_barrier( hecMESH )

    inquire(file = vis_filename, EXIST = is_exit)

    if(.not. is_exit)then
      call fstr_setup_visualize_main( ctrl, vis_filename )
    endif
  end subroutine fstr_setup_visualize

  subroutine fstr_setup_visualize_main( ctrl, vis_filename )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: rcode
    integer(kind=kint) :: i, start_n, end_n
    character(HECMW_FILENAME_LEN) :: vis_filename
    integer(kind=kint), parameter :: buffsize = 127
    character( buffsize ) :: buff
    character( buffsize ) :: head
    character( buffsize ) :: msg

    start_n = fstr_ctrl_get_c_h_pos( ctrl )
    end_n = fstr_ctrl_get_rec_number( ctrl )

    open ( IFVS, file = trim(vis_filename), status = 'replace', err = 1000)
    do i=start_n, end_n
      rcode = fstr_ctrl_get_line( ctrl, i, buff, buffsize )
      if( rcode /= 0 ) exit
      read( buff, *) head
      if( head == '!END') exit
      write( IFVS, '(a)') buff
    end do
    close( IFVS );

    return

    1000    write(msg,*) 'Error: cannot create file:"', trim(vis_filename), '" for visualization'
    call fstr_setup_util_err_stop(msg)
  end subroutine fstr_setup_visualize_main

  !******************************************************************************

end module fstr_setup_util
