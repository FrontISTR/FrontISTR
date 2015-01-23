!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.7                                                !
!                                                                      !
!     Last Update : 2014/07/04                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  FrontISTR Forum,The University of Tokyo       !
!                                                                      !
!======================================================================!

  subroutine hecmw_ML_getrow(id, n_requested_rows, requested_rows, &
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
    integer(kind=kint) :: m, i, row, inod, idof, nl, nd, nu, js, je, j, jj, jdof, start
    ierr = 0
    call hecmw_mat_id_get(id, hecMAT, hecMESH)
    m = 1
    do i = 1, n_requested_rows
      row = requested_rows(i) + 1 ! '+1' for Fortran-numbering
      inod = (row-1)/3 + 1
      idof = row - (inod-1)*3
      nl = (hecMAT%indexL(inod) - hecMAT%indexL(inod-1)) * 3
      nd = 3
      nu = (hecMAT%indexU(inod) - hecMAT%indexU(inod-1)) * 3
      if (allocated_space < m + nl + nd + nu) return
      start = m
      js = hecMAT%indexL(inod-1)+1
      je = hecMAT%indexL(inod)
      do j = js, je
        jj = hecMAT%itemL(j)
        do jdof = 1, 3
          cols(m) = (jj-1)*3 + jdof - 1 ! '-1' for C-numbering
          values(m) = hecMAT%AL((j-1)*9 + (idof-1)*3 + jdof)
          m = m+1
        enddo
      enddo
      do jdof = 1, 3
        cols(m) = (inod-1)*3 + jdof - 1 ! '-1' for C-numbering
        values(m) = hecMAT%D((inod-1)*9 + (idof-1)*3 + jdof)
        m = m+1
      enddo
      js = hecMAT%indexU(inod-1)+1
      je = hecMAT%indexU(inod)
      do j = js, je
        jj = hecMAT%itemU(j)
        do jdof = 1, 3
          cols(m) = (jj-1)*3 + jdof - 1 ! '-1' for C-numbering
          values(m) = hecMAT%AU((j-1)*9 + (idof-1)*3 + jdof)
          m = m+1
        enddo
      enddo
      row_lengths(i) = m - start
    enddo
    ierr = 1
  end subroutine hecmw_ML_getrow

  subroutine hecmw_ML_matvec(id, in_length, p, out_length, ap, ierr)
    use hecmw_util
    use hecmw_mat_id
    use hecmw_solver_las_33
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
    call hecmw_matvec_33(hecMESH, hecMAT, w, ap)
    deallocate(w)
    ierr = 0
  end subroutine hecmw_ML_matvec

  subroutine hecmw_ML_comm(id, x, ierr)
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
    call hecmw_update_3_R (hecMESH, x, hecMESH%n_node)
    ierr = 0
  end subroutine hecmw_ML_comm

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
    level = hecmw_mat_get_timelog(hecMAT)
  end subroutine hecmw_ML_get_loglevel
