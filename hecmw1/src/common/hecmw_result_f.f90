!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module hecmw_result
  use hecmw_util
  use hecmw_etype
  implicit none

  public :: hecmwST_result_data
  public :: hecmw_nullify_result_data
  public :: hecmw_result_copy_c2f
  public :: hecmw_result_copy_f2c
  public :: hecmw_result_init
  public :: hecmw_result_add
  public :: hecmw_result_write_by_name
  public :: hecmw_result_write_st_by_name
  public :: hecmw_result_write_by_addfname
  public :: hecmw_result_finalize
  public :: hecmw_result_free
  public :: hecmw_result_read_by_name
  public :: hecmw_result_checkfile_by_name
  private :: put_node_component
  private :: put_elem_component
  private :: refine_result
  private :: get_node_component
  private :: get_elem_component

  type hecmwST_result_data
    integer(kind=kint) :: ng_component
    integer(kind=kint) :: nn_component
    integer(kind=kint) :: ne_component
    integer(kind=kint),pointer :: ng_dof(:)
    integer(kind=kint),pointer :: nn_dof(:)
    integer(kind=kint),pointer :: ne_dof(:)
    character(len=HECMW_NAME_LEN),pointer :: global_label(:)
    character(len=HECMW_NAME_LEN),pointer :: node_label(:)
    character(len=HECMW_NAME_LEN),pointer :: elem_label(:)
    real(kind=kreal),pointer :: global_val_item(:)
    real(kind=kreal),pointer :: node_val_item(:)
    real(kind=kreal),pointer :: elem_val_item(:)
  end type hecmwST_result_data

  private
  character(len=HECMW_NAME_LEN) :: sname,vname

  ! for PERFORMANCE TUNING (BLOCK REORDERED NODE ID)
  logical :: tuning_block_reorder_on
  integer(kind=kint), pointer :: nid_new2org(:) ! logical node ID reorder table (reordered => original)
  integer(kind=kint), pointer :: org_order_global_node_ID(:)
  integer(kind=kint) :: nnode, nelem
  logical :: MPC_exist
  integer(kind=kint) :: nelem_wo_MPC = 0
  integer(kind=kint), allocatable :: eid_wo_MPC(:)
  integer(kind=kint), allocatable :: elemID_wo_MPC(:)

contains

  !C=============================================================================
  !C nullify pointer
  !C=============================================================================

  subroutine hecmw_nullify_result_data( P )
    type( hecmwST_result_data ) :: P
    nullify( P%ng_dof )
    nullify( P%nn_dof )
    nullify( P%ne_dof )
    nullify( P%global_label )
    nullify( P%node_label )
    nullify( P%elem_label )
    nullify( P%global_val_item )
    nullify( P%node_val_item )
    nullify( P%elem_val_item )
  end subroutine hecmw_nullify_result_data

  !C=============================================================================
  !C Write result data to file
  !C=============================================================================

  subroutine hecmw_result_init(hecMESH, i_step, header, comment)
    type(hecmwST_local_mesh):: hecMESH
    integer(kind=kint) :: nnode, nelem, i_step, ierr
    integer :: i, org_id
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN) :: comment

    integer(kind=kint) :: itype, iS, iE, ic_type, icel

    MPC_exist = .false.
    do itype= 1, hecMESH%n_elem_type
      ic_type = hecMESH%elem_type_item(itype)
      if (hecmw_is_etype_patch(ic_type)) MPC_exist = .true.
      if (hecmw_is_etype_link(ic_type)) MPC_exist = .true.
    end do

    nnode = hecMESH%n_node
    nelem = hecMESH%n_elem

    ! TUNING BLOCK NODE REORDER
    tuning_block_reorder_on = hecMESH%tuning_block_reorder_on

    allocate(nid_new2org(nnode))
    nid_new2org(:) = hecMESH%tuning_block_reorder_new2old(:)

    allocate(org_order_global_node_ID(nnode))
    do i = 1, nnode
      org_id = nid_new2org(i)
      org_order_global_node_ID(org_id) = hecMESH%global_node_ID(i)
    end do

    if( MPC_exist ) then

      if( nelem_wo_MPC == 0 ) then
        allocate(eid_wo_MPC(nelem))
        allocate(elemID_wo_MPC(nelem))
        eid_wo_MPC(:) = 0
        elemID_wo_MPC(:) = 0

        nelem_wo_MPC = 0
        do itype= 1, hecMESH%n_elem_type
          iS= hecMESH%elem_type_index(itype-1) + 1
          iE= hecMESH%elem_type_index(itype  )
          ic_type= hecMESH%elem_type_item(itype)

          if (hecmw_is_etype_patch(ic_type)) cycle
          if (hecmw_is_etype_link(ic_type)) cycle

          do icel= iS, iE
            nelem_wo_MPC = nelem_wo_MPC + 1
            elemID_wo_MPC(nelem_wo_MPC) = hecMESH%global_elem_ID(icel)
            eid_wo_MPC(nelem_wo_MPC) = icel
          end do
        end do
      end if

      call hecmw_result_init_if(nnode, nelem_wo_MPC, org_order_global_node_ID, elemID_wo_MPC, i_step, header, comment, ierr)
    else
      call hecmw_result_init_if(nnode, nelem, org_order_global_node_ID, hecMESH%global_elem_ID, i_step, header, comment, ierr)
    end if
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_init


  subroutine hecmw_result_add(dtype, n_dof, label, data)
    integer(kind=kint) :: dtype, n_dof, ierr
    character(len=HECMW_NAME_LEN) :: label
    real(kind=kreal) :: data(:)
    real(kind=kreal), allocatable :: org_order_data(:)
    integer :: i, j, ofset_org, ofset_new, org_id
    integer(kind=kint) :: icel
    real(kind=kreal), pointer :: data_wo_MPC(:)

    if (dtype == 1) then ! node
      if (tuning_block_reorder_on) then ! reorderd node ID
        allocate(org_order_data(size(data)))

        do i = 1, nnode
          org_id = nid_new2org(i)
          ofset_org = (org_id - 1) * n_dof
          ofset_new = (i      - 1) * n_dof
          do j = 1, n_dof
            org_order_data(ofset_org + j) = data(ofset_new + j)
          end do
        end do

        call hecmw_result_add_if(dtype, n_dof, label, org_order_data, ierr)
        if (ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
        deallocate(org_order_data)

      else ! original node ID
        call hecmw_result_add_if(dtype, n_dof, label, data, ierr)
        if (ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
      end if
    else if( dtype == 2 .and. MPC_exist ) then !element output without patch element

      allocate(data_wo_MPC(n_dof*nelem_wo_MPC))
      data_wo_MPC(:) = 0.d0

      do i= 1, nelem_wo_MPC
        icel = eid_wo_MPC(i)
        data_wo_MPC(n_dof*(i-1)+1:n_dof*i) = data(n_dof*(icel-1)+1:n_dof*icel)
      end do

      call hecmw_result_add_if(dtype, n_dof, label, data_wo_MPC, ierr)

      deallocate(data_wo_MPC)
      if (ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    else
      return ! NEVER COME HERE (not node, not elem)
    end if

  end subroutine hecmw_result_add


  subroutine hecmw_result_write_by_name(name_ID)
    integer(kind=kint) :: ierr
    character(len=HECMW_NAME_LEN) :: name_ID

    call hecmw_result_write_by_name_if(name_ID, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_write_by_name


  subroutine hecmw_result_finalize()
    integer(kind=kint) :: ierr

    call hecmw_result_finalize_if(ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_finalize


  subroutine hecmw_result_write_st_by_name(name_ID, result_data)
    integer(kind=kint) :: ierr
    type(hecmwST_result_data):: result_data
    character(len=HECMW_NAME_LEN):: name_ID

    call hecmw_result_write_st_init_if(ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    call hecmw_result_copy_f2c(result_data, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    call hecmw_result_write_st_by_name_if(name_ID, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    call hecmw_result_write_st_finalize_if(ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_write_st_by_name


  subroutine hecmw_result_write_by_addfname(name_ID, addfname)
    integer(kind=kint) :: ierr
    character(len=HECMW_NAME_LEN) :: name_ID, addfname

    call hecmw_result_write_by_addfname_if(name_ID, addfname, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_write_by_addfname


  subroutine  hecmw_result_copy_f2c( result_data, ierr )
    type(hecmwST_result_data), intent(in)    :: result_data
    integer(kind=kint),        intent(inout) :: ierr

    call  put_global_component( result_data, ierr )
    if( ierr /= 0 )  return
    call  put_node_component( result_data, ierr )
    if( ierr /= 0 )  return
    call  put_elem_component( result_data, ierr )
    if( ierr /= 0 )  return
  end subroutine  hecmw_result_copy_f2c


  subroutine  put_global_component( result_data, ierr )
    type(hecmwST_result_data), intent(in)    :: result_data
    integer(kind=kint),        intent(inout) :: ierr

    sname = "hecmwST_result_data"

    vname = "ng_component"
    call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%ng_component, ierr )
    if( ierr /= 0 )  return

    if( result_data%ng_component /= 0 )  then
      vname = "ng_dof"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%ng_dof, ierr )
      if( ierr /= 0 )  return

      vname = "global_label"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%global_label, ierr )
      if( ierr /= 0 )  return

      vname = "global_val_item"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%global_val_item, ierr )
      if( ierr /= 0 )  return
    endif
  end subroutine  put_global_component

  subroutine  put_node_component( result_data, ierr )
    type(hecmwST_result_data), intent(in)    :: result_data
    integer(kind=kint),        intent(inout) :: ierr

    sname = "hecmwST_result_data"

    vname = "nn_component"
    call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%nn_component, ierr )
    if( ierr /= 0 )  return

    if( result_data%nn_component /= 0 )  then
      vname = "nn_dof"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%nn_dof, ierr )
      if( ierr /= 0 )  return

      vname = "node_label"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%node_label, ierr )
      if( ierr /= 0 )  return

      vname = "node_val_item"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%node_val_item, ierr )
      if( ierr /= 0 )  return
    endif
  end subroutine  put_node_component

  subroutine  put_elem_component( result_data, ierr )
    type(hecmwST_result_data), intent(in)    :: result_data
    integer(kind=kint),        intent(inout) :: ierr

    sname = "hecmwST_result_data"

    vname = "ne_component"
    call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%ne_component, ierr )
    if( ierr /= 0 )  return

    if( result_data%ne_component /= 0 )  then
      vname = "ne_dof"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%ne_dof, ierr )
      if( ierr /= 0 )  return

      vname = "elem_label"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%elem_label, ierr )
      if( ierr /= 0 )  return

      vname = "elem_val_item"
      call  hecmw_result_copy_f2c_set_if( sname, vname, result_data%elem_val_item, ierr )
      if( ierr /= 0 )  return
    endif
  end subroutine  put_elem_component

  !C=============================================================================
  !C Read result data from file
  !C=============================================================================

  subroutine hecmw_result_checkfile_by_name(name_ID, i_step, ierr)
    character(len=HECMW_NAME_LEN), intent(in) :: name_ID
    integer(kind=kint), intent(in) :: i_step
    integer(kind=kint), intent(out) :: ierr

    call hecmw_result_checkfile_by_name_if(name_ID, i_step, ierr)
  end subroutine hecmw_result_checkfile_by_name


  subroutine hecmw_result_read_by_name(hecMESH, name_ID, i_step, result)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    character(len=HECMW_NAME_LEN), intent(in) :: name_ID
    integer(kind=kint), intent(in) :: i_step
    type(hecmwST_result_data), intent(inout) :: result
    integer(kind=kint) :: n_node, n_elem, ierr

    call hecmw_result_read_by_name_if(name_ID, i_step, n_node, n_elem, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_result_copy_c2f(result, n_node, n_elem, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_result_read_finalize_if(ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call refine_result(hecMESH, n_node, result, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_read_by_name


  subroutine refine_result(hecMESH, n_node, result, ierr)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: n_node
    type(hecmwST_result_data), intent(inout) :: result
    integer(kind=kint), intent(out) :: ierr
    real(kind=kreal), pointer :: tmp_val(:)
    integer(kind=kint) :: iref, i, j, k, is, ie, js, je, i0
    integer(kind=kint) :: jj, j0, nn_comp_tot, nn, n_node_ref
    ierr = 0
    if(n_node == hecMESH%n_node) return
    if(n_node > hecMESH%n_node) then
      write(*,*) 'ERROR: result needs to be coarsened; not implemented yet'
      ierr = 1
      return
    else
      !write(0,*) 'DEBUG: result needs to be refined'
      nn_comp_tot = 0
      do i = 1, result%nn_component
        nn_comp_tot = nn_comp_tot + result%nn_dof(i)
      enddo
      do iref = 1, hecMESH%n_refine
        is = hecMESH%refine_origin%index(iref-1)
        ie = hecMESH%refine_origin%index(iref)
        n_node_ref = ie - is
        if(n_node >= n_node_ref) cycle
        !write(0,*) 'DEBUG: start refining result; step=',iref
        allocate(tmp_val(n_node_ref * nn_comp_tot))
        tmp_val = 0.d0
        do i = 1, n_node_ref
          js = hecMESH%refine_origin%item_index(is+i-1)
          je = hecMESH%refine_origin%item_index(is+i)
          nn = je - js
          i0 = (i-1)*nn_comp_tot
          do j = js+1, je
            jj = hecMESH%refine_origin%item_item(j)
            j0 = (jj-1)*nn_comp_tot
            do k = 1, nn_comp_tot
              tmp_val(i0+k) = tmp_val(i0+k) + result%node_val_item(j0+k) / nn
            enddo
          enddo
        enddo
        deallocate(result%node_val_item)
        result%node_val_item => tmp_val
        !write(0,*) 'DEBUG: end refining result; step=',iref
      enddo
      !write(0,*) 'DEBUG: refining result done'
    endif
  end subroutine refine_result


  subroutine hecmw_result_copy_c2f(result, n_node, n_elem, ierr)
    integer(kind=kint) :: n_node, n_elem, ierr
    type(hecmwST_result_data) :: result

    call get_global_component(result, n_node, ierr)
    if(ierr /= 0) return
    call get_node_component(result, n_node, ierr)
    if(ierr /= 0) return
    call get_elem_component(result, n_elem, ierr)
    if(ierr /= 0) return
  end subroutine hecmw_result_copy_c2f


  subroutine get_global_component(result, n_global, ierr)
    integer(kind=kint) :: n_global, ierr
    type(hecmwST_result_data) :: result

    sname = 'hecmwST_result_data'

    vname = 'ng_component'
    call hecmw_result_copy_c2f_set_if(sname, vname, result%ng_component, ierr)
    if(ierr /= 0) return

    if(result%ng_component > 0) then
      vname = 'ng_dof'
      allocate(result%ng_dof(result%ng_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%ng_dof, ierr)
      if(ierr /= 0) return

      vname = 'global_label'
      allocate(result%global_label(result%ng_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%global_label, ierr)
      if(ierr /= 0) return

      vname = 'global_val_item'
      allocate(result%global_val_item(sum(result%ng_dof)*n_global))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%global_val_item, ierr)
      if(ierr /= 0) return
    endif
  end subroutine get_global_component


  subroutine get_node_component(result, n_node, ierr)
    integer(kind=kint) :: n_node, ierr
    type(hecmwST_result_data) :: result

    sname = 'hecmwST_result_data'

    vname = 'nn_component'
    call hecmw_result_copy_c2f_set_if(sname, vname, result%nn_component, ierr)
    if(ierr /= 0) return

    if(result%nn_component > 0) then
      vname = 'nn_dof'
      allocate(result%nn_dof(result%nn_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%nn_dof, ierr)
      if(ierr /= 0) return

      vname = 'node_label'
      allocate(result%node_label(result%nn_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%node_label, ierr)
      if(ierr /= 0) return

      vname = 'node_val_item'
      allocate(result%node_val_item(sum(result%nn_dof)*n_node))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%node_val_item, ierr)
      if(ierr /= 0) return
    endif
  end subroutine get_node_component


  subroutine get_elem_component(result, n_elem, ierr)
    integer(kind=kint) :: n_elem, ierr
    type(hecmwST_result_data) :: result

    sname = 'hecmwST_result_data'

    vname = 'ne_component'
    call hecmw_result_copy_c2f_set_if(sname, vname, result%ne_component, ierr)
    if(ierr /= 0) return

    if(result%ne_component > 0) then
      vname = 'ne_dof'
      allocate(result%ne_dof(result%ne_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%ne_dof, ierr)
      if(ierr /= 0) return

      vname = 'elem_label'
      allocate(result%elem_label(result%ne_component))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%elem_label, ierr)
      if(ierr /= 0) return

      vname = 'elem_val_item'
      allocate(result%elem_val_item(sum(result%ne_dof)*n_elem))
      call hecmw_result_copy_c2f_set_if(sname, vname, result%elem_val_item, ierr)
      if(ierr /= 0) return
    endif
  end subroutine get_elem_component


  subroutine  hecmw_result_free( result_data )
    type(hecmwST_result_data), intent(inout) :: result_data
    integer(kind=kint)                       :: ierr

    ierr = 0

    if( associated( result_data%ng_dof ) )  then
      deallocate( result_data%ng_dof, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%global_label ) )  then
      deallocate( result_data%global_label, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%global_val_item ) )  then
      deallocate( result_data%global_val_item, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%nn_dof ) )  then
      deallocate( result_data%nn_dof, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%node_label ) )  then
      deallocate( result_data%node_label, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%node_val_item ) )  then
      deallocate( result_data%node_val_item, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%ne_dof ) )  then
      deallocate( result_data%ne_dof, stat=ierr )
      if ( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%elem_label ) )  then
      deallocate( result_data%elem_label, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif

    if( associated( result_data%elem_val_item ) )  then
      deallocate( result_data%elem_val_item, stat=ierr )
      if( ierr /= 0 )  then
        print *, "Error: Deallocation error"
        call  hecmw_abort( hecmw_comm_get_comm( ) )
      endif
    endif
  end subroutine  hecmw_result_free

end module hecmw_result
