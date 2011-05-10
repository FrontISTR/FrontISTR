!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Coupling Interface                                 !
!                                                                      !
!            Written by Shin'ichi Ezure (RIST)                         !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!


module hecmw_couple_init_f

  use hecmw_util
  use hecmw_dist_copy_f2c_f
  use hecmw_dist_free_f
  use hecmw_couple_define_f
  use hecmw_couple_info_f

  implicit none
  private
  public :: hecmw_couple_init

contains

subroutine hecmw_couple_init(boundary_id, mesh_unit1, mesh_unit2)

  character(len=HECMW_NAME_LEN), intent(in)    :: boundary_id
  type(hecmwST_local_mesh),      intent(inout) :: mesh_unit1
  type(hecmwST_local_mesh),      intent(inout) :: mesh_unit2
  integer(kind=kint)                           :: is_unit1_memb, is_unit2_memb
  integer(kind=kint)                           :: ierr

  is_unit1_memb = hecmw_couple_is_unit_member(boundary_id, HECMW_COUPLE_UNIT1)
  is_unit2_memb = hecmw_couple_is_unit_member(boundary_id, HECMW_COUPLE_UNIT2)
  if(is_unit1_memb < 0 .or. is_unit2_memb < 0) call hecmw_abort(hecmw_comm_get_comm())

  if(is_unit1_memb == 1) then
    call hecmw_couple_init_init_if(HECMW_COUPLE_UNIT1, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_dist_copy_f2c(mesh_unit1, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  endif

  if(is_unit2_memb == 1) then
    call hecmw_couple_init_init_if(HECMW_COUPLE_UNIT2, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

    call hecmw_dist_copy_f2c(mesh_unit2, ierr)
    if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())
  endif

  call hecmw_couple_init_if(boundary_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_init_final_if(ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

!  if(is_unit1_memb == 1) then
!    mesh_unit1%PETOT    = hecmw_intracomm_get_size(boundary_id, HECMW_COUPLE_UNIT1)
!    mesh_unit1%my_rank  = hecmw_intracomm_get_rank(boundary_id, HECMW_COUPLE_UNIT1)
!    mesh_unit1%MPI_COMM = hecmw_intracomm_get_comm(boundary_id, HECMW_COUPLE_UNIT1)
!    if(mesh_unit1%my_rank == 0) then
!      mesh_unit1%zero = 1
!    else
!      mesh_unit1%zero = 0
!    endif
!  endif
!
!  if(is_unit2_memb == 1) then
!    mesh_unit2%PETOT    = hecmw_intracomm_get_size(boundary_id, HECMW_COUPLE_UNIT2)
!    mesh_unit2%my_rank  = hecmw_intracomm_get_rank(boundary_id, HECMW_COUPLE_UNIT2)
!    mesh_unit2%MPI_COMM = hecmw_intracomm_get_comm(boundary_id, HECMW_COUPLE_UNIT2)
!    if(mesh_unit2%my_rank == 0) then
!      mesh_unit2%zero = 1
!    else
!      mesh_unit2%zero = 0
!    endif
!  endif

end subroutine hecmw_couple_init

end module hecmw_couple_init_f

