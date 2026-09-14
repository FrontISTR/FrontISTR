!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Essential boundary conditions kept as per-DOF marks and values so that
!> they can be imposed on the matrix after the MPC processing

module hecmw_ebc_defer
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private
  public :: hecmwST_ebc
  public :: hecmw_ebc_init
  public :: hecmw_ebc_set
  public :: hecmw_ebc_apply
  public :: hecmw_ebc_finalize

  type hecmwST_ebc
    integer(kind=kint) :: ndof = 0
    integer(kind=kint), pointer :: mark(:) => null()
    real(kind=kreal), pointer :: val(:) => null()
  end type hecmwST_ebc

contains

  !C
  !C***
  !C*** hecmw_ebc_init
  !C***
  !C
  subroutine hecmw_ebc_init(hecMAT, hecEBC)
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_ebc), intent(inout) :: hecEBC
    integer(kind=kint) :: npndof

    call hecmw_ebc_finalize(hecEBC)
    hecEBC%ndof = hecMAT%NDOF
    npndof = hecMAT%NP * hecMAT%NDOF
    allocate(hecEBC%mark(npndof))
    allocate(hecEBC%val(npndof))
    hecEBC%mark(:) = 0
    hecEBC%val(:) = 0.d0
  end subroutine hecmw_ebc_init

  !C
  !C***
  !C*** hecmw_ebc_set
  !C***
  !C
  subroutine hecmw_ebc_set(hecEBC, inode, idof, val)
    implicit none
    type (hecmwST_ebc), intent(inout) :: hecEBC
    integer(kind=kint), intent(in) :: inode, idof
    real(kind=kreal), intent(in) :: val
    integer(kind=kint) :: k

    if (idof > hecEBC%ndof) return
    k = hecEBC%ndof * (inode - 1) + idof
    if (hecEBC%mark(k) /= 0 .and. hecEBC%val(k) /= val) then
      write(*,'(a,i0,a,i0,a,i0,a,2(1pe14.6))') 'WARNING: rank ', hecmw_comm_get_rank(), &
          ': boundary value of local node ', inode, ' dof ', idof, ' overwritten:', hecEBC%val(k), val
    endif
    hecEBC%mark(k) = 1
    hecEBC%val(k) = val
  end subroutine hecmw_ebc_set

  !C
  !C***
  !C*** hecmw_ebc_apply
  !C***
  !C
  subroutine hecmw_ebc_apply(hecMESH, hecMAT, hecEBC, conMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_ebc), intent(inout) :: hecEBC
    type (hecmwST_matrix), intent(inout), optional :: conMAT

    call hecmw_ebc_extend(hecMESH, hecMAT, hecEBC)
    call hecmw_ebc_convert_slave(hecMESH, hecMAT, hecEBC)
    call hecmw_ebc_impose(hecMAT, hecEBC, 1.d0)
    if (present(conMAT)) call hecmw_ebc_impose(conMAT, hecEBC, 0.d0)
  end subroutine hecmw_ebc_apply

  !C
  !C***
  !C*** hecmw_ebc_extend
  !C***
  !C
  !> The matrix may have more nodes than the mesh the marks were set on (the MPC
  !> reduction imports external nodes), and a mark may be set only on a rank that
  !> does not own the node, so the marks are summed onto the owner and distributed
  !> back instead of just being updated from the owner.
  subroutine hecmw_ebc_extend(hecMESH, hecMAT, hecEBC)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_ebc), intent(inout) :: hecEBC
    integer(kind=kint), pointer :: mark(:)
    real(kind=kreal), pointer :: val(:)
    integer(kind=kint) :: ndof, npndof, npndof_old, i

    ndof = hecMAT%NDOF
    npndof = hecMAT%NP * ndof
    npndof_old = size(hecEBC%mark)

    if (npndof > npndof_old) then
      allocate(mark(npndof))
      allocate(val(npndof))
      mark(:) = 0
      val(:) = 0.d0
      do i = 1, npndof_old
        mark(i) = hecEBC%mark(i)
        val(i) = hecEBC%val(i)
      enddo
      deallocate(hecEBC%mark)
      deallocate(hecEBC%val)
      hecEBC%mark => mark
      hecEBC%val => val
    endif

    call hecmw_assemble_I(hecMESH, hecEBC%mark, hecMAT%NP, ndof)
    call hecmw_assemble_R(hecMESH, hecEBC%val, hecMAT%NP, ndof)
    call hecmw_update_I(hecMESH, hecEBC%mark, hecMAT%NP, ndof)
    call hecmw_update_R(hecMESH, hecEBC%val, hecMAT%NP, ndof)

    do i = 1, npndof
      if (hecEBC%mark(i) == 0) cycle
      hecEBC%val(i) = hecEBC%val(i) / hecEBC%mark(i)
      hecEBC%mark(i) = 1
    enddo
  end subroutine hecmw_ebc_extend

  !C
  !C***
  !C*** hecmw_ebc_convert_slave
  !C***
  !C
  !> A condition on the slave DOF of a single-master constraint a_s u_s + a_m u_m = c
  !> is imposed as u_m = (c - a_s u_s) / a_m, since the slave DOF is eliminated.
  !> A constraint with several masters cannot be expressed as a Dirichlet condition.
  !> The change is summed onto the owner rank because the constraint is held only
  !> by the ranks that have all of its nodes.
  subroutine hecmw_ebc_convert_slave(hecMESH, hecMAT, hecEBC)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_ebc), intent(inout) :: hecEBC
    integer(kind=kint), allocatable :: dmark(:)
    real(kind=kreal), allocatable :: dval(:)
    integer(kind=kint) :: ndof, npndof, i, j, k, km, ks, kk, nchange
    real(kind=kreal) :: um, tol

    ndof = hecMAT%NDOF
    npndof = hecMAT%NP * ndof
    allocate(dmark(npndof))
    allocate(dval(npndof))
    dmark(:) = 0
    dval(:) = 0.d0
    nchange = 0

    OUTER: do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      ks = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      if (hecEBC%mark(ks) == 0) cycle
      if (hecMESH%mpc%mpc_index(i) - hecMESH%mpc%mpc_index(i-1) /= 2) then
        write(*,'(a,i0,a,i0,a)') 'ERROR: a boundary condition is given to node ', &
            hecMESH%global_node_ID(hecMESH%mpc%mpc_item(k)), ' dof ', hecMESH%mpc%mpc_dof(k), &
            ', the slave of an !EQUATION with more than one master; it cannot be imposed as a boundary condition'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      km = k + 1
      kk = ndof * (hecMESH%mpc%mpc_item(km) - 1) + hecMESH%mpc%mpc_dof(km)
      um = (hecMESH%mpc%mpc_const(i) - hecMESH%mpc%mpc_val(k) * hecEBC%val(ks)) / hecMESH%mpc%mpc_val(km)
      tol = 1.d-10 * max(1.d0, abs(um))
      if ((hecEBC%mark(kk) /= 0 .and. abs(hecEBC%val(kk) - um) > tol) .or. &
          (dmark(kk) /= 0 .and. abs(dval(kk) - um) > tol)) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a)') 'ERROR: the boundary condition on node ', &
            hecMESH%global_node_ID(hecMESH%mpc%mpc_item(k)), ' dof ', hecMESH%mpc%mpc_dof(k), &
            ', the slave of an !EQUATION, conflicts with the one on its master node ', &
            hecMESH%global_node_ID(hecMESH%mpc%mpc_item(km)), ' dof ', hecMESH%mpc%mpc_dof(km), ''
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      if (hecEBC%mark(kk) == 0) then
        dmark(kk) = 1
        dval(kk) = um
      endif
      dmark(ks) = -1
      nchange = nchange + 1
    enddo OUTER

    call hecmw_allreduce_I1(hecMESH, nchange, hecmw_sum)
    if (nchange > 0) then
      call hecmw_assemble_I(hecMESH, dmark, hecMAT%NP, ndof)
      call hecmw_assemble_R(hecMESH, dval, hecMAT%NP, ndof)
      do i = 1, npndof
        if (dmark(i) > 0) then
          hecEBC%mark(i) = 1
          hecEBC%val(i) = dval(i) / dmark(i)
        else if (dmark(i) < 0) then
          hecEBC%mark(i) = 0
          hecEBC%val(i) = 0.d0
        endif
      enddo
      call hecmw_update_I(hecMESH, hecEBC%mark, hecMAT%NP, ndof)
      call hecmw_update_R(hecMESH, hecEBC%val, hecMAT%NP, ndof)
    endif

    deallocate(dmark)
    deallocate(dval)
  end subroutine hecmw_ebc_convert_slave

  !C
  !C***
  !C*** hecmw_ebc_impose
  !C***
  !C
  !> Same operations as hecmw_mat_ass_bc for every marked DOF, but the columns are
  !> found by sweeping the rows instead of through the row of the constrained node,
  !> whose entries have been moved to the owner rank when the node is external to
  !> the reduced matrix.
  subroutine hecmw_ebc_impose(hecMAT, hecEBC, diag)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_ebc), intent(in) :: hecEBC
    real(kind=kreal), intent(in) :: diag
    integer(kind=kint) :: ndof, ndof2, i, j, k, idof, jdof, ir, jc, idx

    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof

    do i = 1, hecMAT%NP
      do idof = 1, ndof
        ir = ndof * (i - 1) + idof
        do jdof = 1, ndof
          jc = ndof * (i - 1) + jdof
          idx = ndof2 * (i - 1) + ndof * (idof - 1) + jdof
          if (hecEBC%mark(ir) /= 0) then
            hecMAT%D(idx) = 0.d0
          else if (hecEBC%mark(jc) /= 0) then
            hecMAT%B(ir) = hecMAT%B(ir) - hecMAT%D(idx) * hecEBC%val(jc)
            hecMAT%D(idx) = 0.d0
          endif
        enddo
      enddo
      do k = hecMAT%indexL(i-1) + 1, hecMAT%indexL(i)
        j = hecMAT%itemL(k)
        do idof = 1, ndof
          ir = ndof * (i - 1) + idof
          do jdof = 1, ndof
            jc = ndof * (j - 1) + jdof
            idx = ndof2 * (k - 1) + ndof * (idof - 1) + jdof
            if (hecEBC%mark(ir) /= 0) then
              hecMAT%AL(idx) = 0.d0
            else if (hecEBC%mark(jc) /= 0) then
              hecMAT%B(ir) = hecMAT%B(ir) - hecMAT%AL(idx) * hecEBC%val(jc)
              hecMAT%AL(idx) = 0.d0
            endif
          enddo
        enddo
      enddo
      do k = hecMAT%indexU(i-1) + 1, hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        do idof = 1, ndof
          ir = ndof * (i - 1) + idof
          do jdof = 1, ndof
            jc = ndof * (j - 1) + jdof
            idx = ndof2 * (k - 1) + ndof * (idof - 1) + jdof
            if (hecEBC%mark(ir) /= 0) then
              hecMAT%AU(idx) = 0.d0
            else if (hecEBC%mark(jc) /= 0) then
              hecMAT%B(ir) = hecMAT%B(ir) - hecMAT%AU(idx) * hecEBC%val(jc)
              hecMAT%AU(idx) = 0.d0
            endif
          enddo
        enddo
      enddo
    enddo

    do i = 1, hecMAT%NP
      do idof = 1, ndof
        ir = ndof * (i - 1) + idof
        if (hecEBC%mark(ir) == 0) cycle
        hecMAT%D(ndof2 * (i - 1) + ndof * (idof - 1) + idof) = diag
        hecMAT%B(ir) = diag * hecEBC%val(ir)
      enddo
    enddo
  end subroutine hecmw_ebc_impose

  !C
  !C***
  !C*** hecmw_ebc_finalize
  !C***
  !C
  subroutine hecmw_ebc_finalize(hecEBC)
    implicit none
    type (hecmwST_ebc), intent(inout) :: hecEBC

    hecEBC%ndof = 0
    if (associated(hecEBC%mark)) deallocate(hecEBC%mark)
    if (associated(hecEBC%val)) deallocate(hecEBC%val)
    nullify(hecEBC%mark)
    nullify(hecEBC%val)
  end subroutine hecmw_ebc_finalize

end module hecmw_ebc_defer
