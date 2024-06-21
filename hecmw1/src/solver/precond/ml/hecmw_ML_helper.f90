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
