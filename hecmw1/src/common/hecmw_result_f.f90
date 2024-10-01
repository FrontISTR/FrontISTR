!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module hecmw_result
  use hecmw_util
  use hecmw_etype
  implicit none

  private
  public :: hecmwST_result_data
  public :: HECMW_RESULT_DTYPE_NODE
  public :: HECMW_RESULT_DTYPE_ELEM
  public :: HECMW_RESULT_DTYPE_GLOBAL
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

  ! constants defined in hecmw_result_io.h
  integer(kind=kint), parameter :: HECMW_RESULT_DTYPE_NODE   = 1
  integer(kind=kint), parameter :: HECMW_RESULT_DTYPE_ELEM   = 2
  integer(kind=kint), parameter :: HECMW_RESULT_DTYPE_GLOBAL = 3

  character(len=HECMW_NAME_LEN) :: sname,vname

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
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN) :: comment

    nnode = hecMESH%n_node
    nelem = hecMESH%n_elem

    call hecmw_result_init_if(nnode, nelem, hecMESH%global_node_ID, hecMESH%global_elem_ID, &
        hecMESH%n_elem_type, hecMESH%elem_type_index, hecMESH%elem_type_item, &
        i_step, header, comment, ierr)

    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_result_init


  subroutine hecmw_result_add(dtype, n_dof, label, data)
    integer(kind=kint) :: dtype, n_dof, ierr
    character(len=HECMW_NAME_LEN) :: label
    real(kind=kreal) :: data(:)

    call hecmw_result_add_if(dtype, n_dof, label, data, ierr)

    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
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
