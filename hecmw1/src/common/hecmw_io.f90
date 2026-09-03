!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief I/O and Utility

module hecmw_io
  use hecmw_util
  use hecmw_dist_copy_f2c_f
  use hecmw_dist_copy_c2f_f
  use hecmw_dist_free_f
  use hecmw_dist_print_f
  use hecmw_result
  use hecmw_restart
  use iso_c_binding
  implicit none

  public :: hecmw_get_mesh
  public :: hecmw_put_mesh

  interface

    subroutine hecmw_get_entire_mesh_init_if(filename,err,len) bind(C,name='hecmw_get_entire_mesh_init_if')
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: filename(*)
      integer(c_int) :: err
      integer(c_int), value :: len
    end subroutine

  end interface

contains

  !C====================================================================
  !C Get HEC-MW dist mesh from file
  !C====================================================================

  subroutine hecmw_get_mesh(name_ID, mesh)
    integer(kind=kint) :: ierr
    character(len=HECMW_NAME_LEN) :: name_ID
    type(hecmwST_local_mesh) :: mesh

    call hecmw_nullify_mesh(mesh)

    call hecmw_get_mesh_init_if(name_ID,ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_dist_copy_c2f(mesh, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_get_mesh_finalize_if(ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

  end subroutine hecmw_get_mesh


  !C====================================================================
  !C Put HEC-MW dist mesh to file
  !C====================================================================

  subroutine hecmw_put_mesh(name, mesh)
    integer(kind=kint) :: ierr
    character(len=HECMW_NAME_LEN) :: name
    type(hecmwST_local_mesh) :: mesh

    call hecmw_put_mesh_init_if(ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_dist_copy_f2c(mesh, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_put_mesh_if(name, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_put_mesh_finalize_if(ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
  end subroutine hecmw_put_mesh

  !C====================================================================
  !C Get HEC-MW entire mesh from msh file
  !C====================================================================

  subroutine hecmw_get_entire_mesh(msh_filename, mesh)
    use hecmw_util
    integer(kind=kint) :: ierr
    character(len=HECMW_FILENAME_LEN) :: msh_filename
    type(hecmwST_local_mesh) :: mesh

    call hecmw_nullify_mesh(mesh)

    call hecmw_get_entire_mesh_init_if(msh_filename,ierr,HECMW_FILENAME_LEN)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_dist_copy_c2f(mesh, ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_get_mesh_finalize_if(ierr)
    if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

  end subroutine hecmw_get_entire_mesh
  

end module hecmw_io

