!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

subroutine hecmw_ML_get_nlocal(id, nlocal, nlocal_allcolumns, ierr)
  use hecmw_util
  use hecmw_mat_id
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: nlocal
  integer(kind=kint), intent(out) :: nlocal_allcolumns
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  nlocal = hecMAT%N * hecMAT%NDOF
  nlocal_allcolumns = hecMAT%NP * hecMAT%NDOF
  ierr = 0
end subroutine hecmw_ML_get_nlocal

subroutine hecmw_ML_get_coord(id, x, y, z, ierr)
  use hecmw_util
  use hecmw_mat_id
  implicit none
  integer(kind=kint), intent(in) :: id
  real(kind=kreal), intent(out) :: x(*), y(*), z(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  integer(kind=kint) :: offset, i
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  offset = 0
  do i = 1, hecMESH%nn_internal
    x(i) = hecMESH%node(offset+1)
    y(i) = hecMESH%node(offset+2)
    z(i) = hecMESH%node(offset+3)
    offset = offset + 3
  enddo
  ierr = 0
end subroutine hecmw_ML_get_coord

subroutine hecmw_ML_get_rbm(id, rbm, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_precond_rbm, only: hecmw_precond_rbm_from_mesh
  implicit none
  integer(kind=kint), intent(in) :: id
  real(kind=kreal), intent(out) :: rbm(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  ! shared rigid-body near-kernel (with SA-AMG)
  call hecmw_precond_rbm_from_mesh(hecMESH, hecMAT%NDOF, rbm)
  ierr = 0
end subroutine hecmw_ML_get_rbm

subroutine hecmw_ML_get_loglevel(id, level)
  use hecmw_util
  use hecmw_matrix_misc
  use hecmw_mat_id
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: level
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  level = hecmw_mat_get_loglevel(hecMAT)        ! independent LOGLEVEL ...
  if (level < 0) level = hecmw_mat_get_timelog(hecMAT)   ! ... unset -> fall back to TIMELOG
end subroutine hecmw_ML_get_loglevel

subroutine hecmw_ml_get_opt(id, opt, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_matrix_misc
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: opt(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  integer(kind=kint) :: iopt(10)
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_mat_get_solver_opt(hecMAT, iopt)
  opt(1:10) = iopt(1:10)
  ierr = 0
end subroutine hecmw_ml_get_opt

subroutine hecmw_ml_set_opt(id, opt, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_matrix_misc
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(in) :: opt(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  integer(kind=kint) :: iopt(10)
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  iopt(1:10) = opt(1:10)
  call hecmw_mat_set_solver_opt(hecMAT, iopt)
  ierr = 0
end subroutine hecmw_ml_set_opt
