!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

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
end subroutine hecmw_ML_getrow

subroutine hecmw_ML_matvec(id, in_length, p, out_length, ap, ierr)
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
  call hecmw_update_m_R (hecMESH, x, hecMESH%n_node,hecMESH%n_dof)
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

subroutine hecmw_ML_get_rbm(id, rbm, ierr)
  use hecmw_util
  use hecmw_mat_id
  use hecmw_etype
  implicit none
  integer(kind=kint), intent(in) :: id
  real(kind=kreal), intent(out) :: rbm(*)
  integer(kind=kint), intent(out) :: ierr
  type(hecmwST_matrix), pointer :: hecMAT
  type(hecmwST_local_mesh), pointer :: hecMESH
  integer(kind=kint) :: Ndof, vec_leng, node, offset
  real(kind=kreal) :: x, y, z

  integer(kind=kint), allocatable :: mark(:)
  integer(kind=kint) :: itype, ic_type, nn, is, iE, icel, iiS, j, nod

  call hecmw_mat_id_get(id, hecMAT, hecMESH)

  ! Mark nodes for rotational DOF
  allocate(mark(hecMESH%n_node))
  mark = 0
  do itype = 1, hecMESH%n_elem_type
    ic_type = hecMESH%elem_type_item(itype)
    if (hecmw_is_etype_33struct(ic_type)) then
      nn = hecmw_get_max_node(ic_type)
      is = hecMESH%elem_type_index(itype-1)+1
      iE = hecMESH%elem_type_index(itype  )
      do icel = is, iE
        iiS = hecMESH%elem_node_index(icel-1)
        ! mark latter halves of the nodes
        do j = nn/2+1, nn
          nod = hecMESH%elem_node_item(iiS+j)
          mark(nod) = 1
        enddo
      enddo
    endif
  enddo
  Ndof = hecMAT%NDOF
  vec_leng = hecMESH%nn_internal * Ndof

  if (Ndof == 1) then

    rbm(1:vec_leng)=1.d0

  else if (Ndof == 2) then

    rbm(1:vec_leng*3)=0.0d0

    do node = 1, hecMESH%nn_internal
      x = hecMESH%node(3*node-2)
      y = hecMESH%node(3*node-1)

      ! translation x
      offset = (node-1)*Ndof
      rbm(offset+1)=1.d0
      rbm(offset+2)=0.d0

      ! translation y
      offset = offset + vec_leng
      rbm(offset+1)=0.d0
      rbm(offset+2)=1.d0

      ! rotation z
      offset = offset + vec_leng
      rbm(offset+1)= -y
      rbm(offset+2)=  x

    enddo

  else

    rbm(1:vec_leng*6)=0.0d0
    do node = 1, hecMESH%nn_internal
      if (mark(node) == 0) then
!!! translational DOF

        x = hecMESH%node(3*node-2)
        y = hecMESH%node(3*node-1)
        z = hecMESH%node(3*node  )

        ! translation x
        offset = (node-1)*Ndof
        rbm(offset+1)=1.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=0.d0

        ! translation y
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=1.d0
        rbm(offset+3)=0.d0

        ! translation z
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=1.d0

        ! rotation x
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)= -z
        rbm(offset+3)=  y

        ! rotation y
        offset = offset + vec_leng
        rbm(offset+1)=  z
        rbm(offset+2)=0.d0
        rbm(offset+3)= -x

        ! rotation z
        offset = offset + vec_leng
        rbm(offset+1)= -y
        rbm(offset+2)=  x
        rbm(offset+3)=0.d0

      else
!!! rotational DOF

        ! translation x
        offset = (node-1)*Ndof
        rbm(offset+1)=0.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=0.d0

        ! translation y
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=0.d0

        ! translation z
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=0.d0

        ! rotation x
        offset = offset + vec_leng
        rbm(offset+1)=1.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=0.d0

        ! rotation y
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=1.d0
        rbm(offset+3)=0.d0

        ! rotation z
        offset = offset + vec_leng
        rbm(offset+1)=0.d0
        rbm(offset+2)=0.d0
        rbm(offset+3)=1.d0
      endif
    enddo
  endif

  deallocate(mark)
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
  level = hecmw_mat_get_timelog(hecMAT)
end subroutine hecmw_ML_get_loglevel
