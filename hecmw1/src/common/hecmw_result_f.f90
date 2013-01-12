!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.3                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuaki Sakane (RIST)                          !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_result
    use hecmw_util
    implicit none
    private

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
    private :: put_node_component
    private :: put_elem_component
    private :: get_node_component
    private :: get_elem_component

    character(len=HECMW_NAME_LEN) :: sname,vname

    type hecmwST_result_data
        integer(kind=kint) :: nn_component
        integer(kind=kint) :: ne_component
        integer(kind=kint),pointer :: nn_dof(:)
        integer(kind=kint),pointer :: ne_dof(:)
        character(len=HECMW_NAME_LEN),pointer :: node_label(:)
        character(len=HECMW_NAME_LEN),pointer :: elem_label(:)
        real(kind=kreal),pointer :: node_val_item(:)
        real(kind=kreal),pointer :: elem_val_item(:)
    end type hecmwST_result_data

    contains

!C=============================================================================
!C nullify pointer
!C=============================================================================

    subroutine hecmw_nullify_result_data( P )
        type( hecmwST_result_data ) :: P
        nullify( P%nn_dof )
        nullify( P%ne_dof )
        nullify( P%node_label )
        nullify( P%elem_label )
        nullify( P%node_val_item )
        nullify( P%elem_val_item )
    end subroutine hecmw_nullify_result_data

!C=============================================================================
!C Write result data to file
!C=============================================================================

    subroutine hecmw_result_init(hecMESH, n_step, i_step, header)
        type(hecmwST_local_mesh):: hecMESH
        integer(kind=kint) :: nnode, nelem, n_step, i_step, ierr
        character(len=HECMW_HEADER_LEN) :: header

        nnode = hecMESH%n_node
        nelem = hecMESH%n_elem

        call hecmw_result_init_if(nnode, nelem, hecMESH%global_node_ID, hecMESH%global_elem_ID, n_step, i_step, header, ierr)
        if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_result_init


    subroutine hecmw_result_add(node_or_elem, n_dof, label, data)
        integer(kind=kint) :: node_or_elem, n_dof, ierr
        character(len=HECMW_NAME_LEN) :: label 
        real(kind=kreal) :: data(:)

        call hecmw_result_add_if(node_or_elem, n_dof, label, data, ierr)
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

        call  put_node_component( result_data, ierr )
        if( ierr /= 0 )  return
        call  put_elem_component( result_data, ierr )
        if( ierr /= 0 )  return
    end subroutine  hecmw_result_copy_f2c


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

    subroutine hecmw_result_read_by_name(name_ID, n_step, i_step, result)
        integer(kind=kint) :: n_node, n_elem, n_step, i_step, ierr
        character(len=HECMW_NAME_LEN) :: name_ID
        type(hecmwST_result_data) :: result

        call hecmw_result_read_by_name_if(name_ID, n_step, i_step, n_node, n_elem, ierr)
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_result_copy_c2f(result, n_node, n_elem, ierr)
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_result_read_finalize_if(ierr)
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_result_read_by_name


    subroutine hecmw_result_copy_c2f(result, n_node, n_elem, ierr)
        integer(kind=kint) :: n_node, n_elem, ierr
        type(hecmwST_result_data) :: result

        call get_node_component(result, n_node, ierr)
        if(ierr /= 0) return

        call get_elem_component(result, n_elem, ierr)
        if(ierr /= 0) return
    end subroutine hecmw_result_copy_c2f


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
