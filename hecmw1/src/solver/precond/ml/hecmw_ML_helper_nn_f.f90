!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

subroutine hecmw_ML_getrow_nn(id, n_requested_rows, requested_rows, &
    allocated_space, cols, values, row_lengths, ierr)
  use hecmw_util
  use hecmw_mat_id
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(in) :: n_requested_rows
  integer(kind=kint), intent(in) :: requested_rows(n_requested_rows)
  integer(kind=kint), intent(in) :: allocated_space
  integer(kind=kint), intent(out) :: cols(allocated_space)
  real(kind=kreal), intent(out) :: values(allocated_space)
  integer(kind=kint), intent(out) :: row_lengths(n_requested_rows)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  integer(kind=kint) :: m, i, row, inod, idof, nl, nd, nu, js, je, j, jj, jdof, start, ndof
  ierr = 0
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  ndof = hecMAT%NDOF
  m = 1
  do i = 1, n_requested_rows
    row = requested_rows(i) + 1 ! '+1' for Fortran-numbering
    inod = (row-1)/ndof + 1
    idof = row - (inod-1)*ndof
    nl = (hecMAT%indexL(inod) - hecMAT%indexL(inod-1)) * ndof
    nd = ndof
    nu = (hecMAT%indexU(inod) - hecMAT%indexU(inod-1)) * ndof
    if (allocated_space < m + nl + nd + nu) return
    start = m
    js = hecMAT%indexL(inod-1)+1
    je = hecMAT%indexL(inod)
    do j = js, je
      jj = hecMAT%itemL(j)
      do jdof = 1, ndof
        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
        values(m) = hecMAT%AL((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
        m = m+1
      enddo
    enddo
    do jdof = 1, ndof
      cols(m) = (inod-1)*ndof + jdof - 1 ! '-1' for C-numbering
      values(m) = hecMAT%D((inod-1)*ndof*ndof + (idof-1)*ndof + jdof)
      m = m+1
    enddo
    js = hecMAT%indexU(inod-1)+1
    je = hecMAT%indexU(inod)
    do j = js, je
      jj = hecMAT%itemU(j)
      do jdof = 1, ndof
        cols(m) = (jj-1)*ndof + jdof - 1 ! '-1' for C-numbering
        values(m) = hecMAT%AU((j-1)*ndof*ndof + (idof-1)*ndof + jdof)
        m = m+1
      enddo
    enddo
    row_lengths(i) = m - start
  enddo
  ierr = 1
end subroutine hecmw_ML_getrow_nn

subroutine hecmw_ML_matvec_nn(id, in_length, p, out_length, ap, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_solver_las
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(in) :: in_length
  real(kind=kreal), intent(in) :: p(in_length)
  integer(kind=kint), intent(in) :: out_length
  real(kind=kreal), intent(out) :: ap(out_length)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  real(kind=kreal), allocatable :: w(:)
  integer(kind=kint) :: i
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  allocate(w(hecMAT%NP*hecMAT%NDOF))
  do i = 1, hecMAT%N*hecMAT%NDOF
    w(i) = p(i)
  enddo
  call hecmw_matvec(hecMESH, hecMAT, w, ap)
  deallocate(w)
  ierr = 0
end subroutine hecmw_ML_matvec_nn

subroutine hecmw_ML_comm_nn(id, x, ierr)
  use hecmw_util
  use hecmw_mat_id
  use m_hecmw_comm_f
  implicit none
  integer(kind=kint), intent(in) :: id
  real(kind=kreal), intent(inout) :: x(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_update_R (hecMESH, x, hecMAT%NP, hecMAT%NDOF)
  ierr = 0
end subroutine hecmw_ML_comm_nn

subroutine hecmw_ML_smoother_diag_setup_nn(id, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_precond_DIAG_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_precond_DIAG_nn_setup(hecMAT)
  ierr = 0
end subroutine hecmw_ML_smoother_diag_setup_nn

subroutine hecmw_ML_smoother_diag_apply_nn(id, x_length, x, rhs_length, rhs, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_matrix_misc
  use hecmw_solver_las
  use hecmw_precond_DIAG_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(in) :: x_length
  real(kind=kreal), intent(inout) :: x(x_length)
  integer(kind=kint), intent(in) :: rhs_length
  real(kind=kreal), intent(in) :: rhs(rhs_length)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH

  real(kind=kreal), allocatable :: resid(:)
  integer(kind=kint) :: i
  real(kind=kreal) :: COMMtime
  integer(kind=kint) :: num_sweeps, i_sweep
  integer(kind=kint) :: opt(10)

  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_mat_get_solver_opt(hecMAT, opt)
  num_sweeps = opt(6)
  allocate(resid(hecMAT%NP * hecMAT%NDOF))
  do i_sweep = 1, num_sweeps
    ! {resid} = {rhs} - [A] {x}
    call hecmw_matresid(hecMESH, hecMAT, x, rhs, resid, COMMtime)
    ! {delta_x} = [M]^-1 {resid}
    call hecmw_precond_DIAG_nn_apply(resid, hecMAT%NDOF)
    ! {x} = {x} + {delta_x}
    do i=1,x_length
      x(i) = x(i) + resid(i)
    enddo
  enddo
  deallocate(resid)
  ierr = 0
end subroutine hecmw_ML_smoother_diag_apply_nn

subroutine hecmw_ML_smoother_diag_clear_nn(id, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_precond_DIAG_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_precond_DIAG_nn_clear()
  ierr = 0
end subroutine hecmw_ML_smoother_diag_clear_nn

subroutine hecmw_ML_smoother_ssor_setup_nn(id, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_precond_SSOR_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_precond_SSOR_nn_setup(hecMAT)
  ierr = 0
end subroutine hecmw_ML_smoother_ssor_setup_nn

subroutine hecmw_ML_smoother_ssor_apply_nn(id, x_length, x, rhs_length, rhs, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_matrix_misc
  use hecmw_solver_las
  use hecmw_precond_SSOR_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(in) :: x_length
  real(kind=kreal), intent(inout) :: x(x_length)
  integer(kind=kint), intent(in) :: rhs_length
  real(kind=kreal), intent(in) :: rhs(rhs_length)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH

  real(kind=kreal), allocatable :: resid(:)
  integer(kind=kint) :: i
  real(kind=kreal) :: COMMtime
  integer(kind=kint) :: num_sweeps, i_sweep
  integer(kind=kint) :: opt(10)

  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_mat_get_solver_opt(hecMAT, opt)
  num_sweeps = opt(6)
  allocate(resid(hecMAT%NP * hecMAT%NDOF))
  do i_sweep = 1, num_sweeps
    ! {resid} = {rhs} - [A] {x}
    call hecmw_matresid(hecMESH, hecMAT, x, rhs, resid, COMMtime)
    ! {delta_x} = [M]^-1 {resid}
    call hecmw_precond_SSOR_nn_apply(resid, hecMAT%NDOF)
    ! {x} = {x} + {delta_x}
    do i=1,x_length
      x(i) = x(i) + resid(i)
    enddo
  enddo
  deallocate(resid)
  ierr = 0
end subroutine hecmw_ML_smoother_ssor_apply_nn

subroutine hecmw_ML_smoother_ssor_clear_nn(id, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_precond_SSOR_nn
  implicit none
  integer(kind=kint), intent(in) :: id
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  call hecmw_mat_id_get(id, hecMAT, hecMESH)
  call hecmw_precond_SSOR_nn_clear(hecMAT)
  ierr = 0
end subroutine hecmw_ML_smoother_ssor_clear_nn
