!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module  hecmw_visualizer

  use hecmw_util
  use hecmw_result
  use hecmw_dist_copy_f2c_f

  implicit none

  public  :: hecmw_visualize
  public  :: hecmw_visualize_init
  public  :: hecmw_visualize_finalize

  private
  character(len=100) :: sname, vname

contains

  subroutine  hecmw_visualize( mesh, result_data, step )
    implicit none
    type(hecmwST_local_mesh),  intent(in) :: mesh
    type(hecmwST_result_data), intent(in) :: result_data
    integer(kind=kint),        intent(in) :: step
    integer(kind=kint)                    :: ierr

    call  hecmw_visualize_init_if( mesh%n_node, mesh%n_elem, ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif

    call  hecmw_dist_copy_f2c( mesh, ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif

    call  hecmw_result_copy_f2c( result_data, ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif

    call  hecmw_visualize_if( step, ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif

    call  hecmw_visualize_finalize_if( ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif
  end subroutine  hecmw_visualize


  subroutine  hecmw_visualize_init( )
    implicit none
    integer(kind=kint) :: ierr

    call  hecmw_init_for_visual_if( ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif
  end subroutine  hecmw_visualize_init


  subroutine hecmw_visualize_finalize( )
    implicit none
    integer(kind=kint) :: ierr

    call  hecmw_finalize_for_visual_if( ierr )
    if( ierr /= 0 )  then
      call  hecmw_abort( hecmw_comm_get_comm( ) )
    endif
  end subroutine hecmw_visualize_finalize

end module  hecmw_visualizer
