!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides interface of iteratie linear equation solver for
!! contact problems using Lagrange multiplier.
module m_solve_LINEQ_iter_contact
  use hecmw_util
  use hecmw_local_matrix
  use m_hecmw_contact_comm
  use hecmw_solver_iterative
  use hecmw_matrix_misc
  use m_hecmw_comm_f
  use hecmw_solver_misc
  use hecmw_solver_las

  implicit none

  private
  public :: solve_LINEQ_iter_contact_init
  public :: solve_LINEQ_iter_contact

  logical, save :: INITIALIZED = .false.
  integer, save :: SymType = 0

  integer, parameter :: DEBUG = 0  ! 0: no message, 1: some messages, 2: more messages, 3: vector output, 4: matrix output

contains

  subroutine solve_LINEQ_iter_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type hecmwST_matrix_lagrange
    logical, intent(in) :: is_sym

    if (INITIALIZED) then
      INITIALIZED = .false.
    endif

    hecMAT%Iarray(98) = 1
    hecMAT%Iarray(97) = 1

    if (is_sym) then
      SymType = 1
    else
      SymType = 0
    endif

    INITIALIZED = .true.
  end subroutine solve_LINEQ_iter_contact_init

  subroutine solve_LINEQ_iter_contact(hecMESH,hecMAT,hecLagMAT,istat,conMAT,is_contact_active)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint), intent(out) :: istat
    type (hecmwST_matrix), intent(in) :: conMAT
    logical,intent(in)                :: is_contact_active
    integer :: method_org, precond_org
    logical :: fg_amg
    integer(kind=kint) :: is_contact
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    hecMAT%Iarray(97) = 1

    fg_amg = .false.
    precond_org = hecmw_mat_get_precond(hecMAT)
    if (precond_org == 5) fg_amg = .true.

    is_contact = 0
    if( is_contact_active ) is_contact = 1
    call hecmw_allreduce_I1(hecMESH, is_contact, hecmw_max)

    if ( is_contact == 0 ) then
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: no contact'
      ! use CG because the matrix is symmetric
      method_org = hecmw_mat_get_method(hecMAT)
      call hecmw_mat_set_method(hecMAT, 1)
      ! solve
      call hecmw_solve_iterative(hecMESH,hecMAT)
      if (fg_amg .and. hecmw_mat_get_flag_diverged(hecMAT) /= 0) then
        ! avoid ML and retry when diverged
        call hecmw_mat_set_precond(hecMAT, 3) ! set diag-scaling
        hecMAT%Iarray(97:98) = 1
        hecMAT%X(:) = 0.d0
        call hecmw_solve_iterative(hecMESH,hecMAT)
        call hecmw_mat_set_precond(hecMAT, 5) ! restore amg
      endif
      ! restore solver setting
      call hecmw_mat_set_method(hecMAT, method_org)
    else
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: with contact'
      call solve_eliminate(hecMESH, hecMAT, hecLagMAT, conMAT)
    endif

    istat = hecmw_mat_get_flag_diverged(hecMAT)
  end subroutine solve_LINEQ_iter_contact

  !> \brief Solve with elimination of Lagrange-multipliers
  !>
  subroutine solve_eliminate(hecMESH,hecMAT,hecLagMAT,conMAT)
    type (hecmwST_local_mesh), intent(in), target :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    type (hecmwST_matrix), intent(in) :: conMAT
    integer(kind=kint) :: ndof
    integer(kind=kint), allocatable :: dof_type(:)         !< type of dof; >0: slave (number is slave ID), <=0: other
    integer(kind=kint), allocatable :: slaves_detected(:)  !< list of slave dofs detected in this subdomain
    real(kind=kreal), allocatable   :: BLs_inv(:)   !< inverse of diagonal BLs matrix
    real(kind=kreal), allocatable   :: BUs_inv(:)   !< inverse of diagonal BUs matrix
    type(hecmwST_local_matrix) :: BTmat        !< blocked T matrix
    type(hecmwST_local_matrix) :: BTtmat       !< blocked T^t matrix
    type(hecmwST_matrix)       :: hecTKT
    type(hecmwST_local_mesh)   :: hecMESHtmp
    real(kind=kreal) :: t1
    integer(kind=kint) :: i, ilag
    type(hecmwST_local_matrix) :: BKmat, BTtKmat, BTtKTmat
    integer(kind=kint), allocatable :: mark(:)
    type(hecmwST_contact_comm) :: conCOMM
    integer(kind=kint)              :: n_contact_dof    !< num of contact dofs including masters and slaves
    integer(kind=kint), allocatable :: contact_dofs(:)  !< list of contact dofs including masters and slaves
    integer(kind=kint)              :: n_slave    !< num of internal dofs that are in contact in some subdomain
    integer(kind=kint), allocatable :: slaves(:)  !< list of internal dofs that are in contact in some subdomain
    real(kind=kreal), allocatable :: Btot(:)
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    t1 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate start', hecmw_wtime()-t1

    ndof=hecMAT%NDOF
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: num_lagrange',hecLagMAT%num_lagrange

    ! copy hecMESH to hecMESHtmp
    call copy_mesh(hecMESH, hecMESHtmp)

    ! choose slave DOFs to be eliminated with Lag. DOFs
    allocate(dof_type(hecMAT%NP*ndof))
    allocate(slaves_detected(hecLagMAT%num_lagrange), BLs_inv(hecLagMAT%num_lagrange), &
      BUs_inv(hecLagMAT%num_lagrange))
    call choose_slaves(hecMAT, hecLagMAT, dof_type, slaves_detected, BLs_inv, BUs_inv)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: slave DOFs chosen', hecmw_wtime()-t1

    ! make transformation matrix and its transpose
    allocate(mark(hecMAT%NP*ndof))
    call make_transformation_matrices(hecMAT, hecLagMAT, hecMESH, hecMESHtmp, &
        dof_type, slaves_detected, BLs_inv, BUs_inv, mark, n_slave, slaves, BTmat, BTtmat)
    deallocate(dof_type)

    ! init BKmat and substitute conMAT
    call hecmw_localmat_init_with_hecmat(BKmat, conMAT, hecLagMAT%num_lagrange)
    if (DEBUG >= 4) call debug_write_matrix(BKmat, 'BKmat (conMAT local)')

    if ( hecmw_comm_get_size() > 1) then
      ! communicate and complete BKmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BKmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble K (conMAT) done', hecmw_wtime()-t1
      if (BKmat%nc /= BTtmat%nc) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with BKmat',BKmat%nc-BTtmat%nc
        BTmat%nc = BKmat%nc
        BTtmat%nc = BKmat%nc
      endif
      if (DEBUG >= 4) call debug_write_matrix(BKmat, 'BKmat (conMAT assembled)')
    endif

    ! add hecMAT to BKmat
    call hecmw_localmat_add_hecmat(BKmat, hecMAT)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: add hecMAT to K done', hecmw_wtime()-t1
    if (DEBUG >= 4) call debug_write_matrix(BKmat, 'BKmat (hecMAT added)')

    ! compute BTtKmat = BTtmat * BKmat (update hecMESHtmp)
    call hecmw_localmat_multmat(BTtmat, BKmat, hecMESHtmp, BTtKmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply Tt and K done', hecmw_wtime()-t1
    if (DEBUG >= 4) call debug_write_matrix(BTtKmat, 'BTtKmat')

    ! compute BTtKTmat = BTtKmat * BTmat (update hecMESHtmp)
    call hecmw_localmat_multmat(BTtKmat, BTmat, hecMESHtmp, BTtKTmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply TtK and T done', hecmw_wtime()-t1
    if (DEBUG >= 4) call debug_write_matrix(BTtKTmat, 'BTtKTmat')
    call hecmw_localmat_free(BTtKmat)

    ! shrink comm_table
    ! call hecmw_localmat_shrink_comm_table(BTtKTmat, hecMESHtmp)

    call place_num_on_diag_of_marked_dof(BTtKTmat, 1.0d0, mark)
    if (DEBUG >= 4) call debug_write_matrix(BTtKTmat, 'BTtKTmat (place 1.0 on slave diag)')
    call hecmw_mat_init(hecTKT)
    call hecmw_localmat_make_hecmat(hecMAT, BTtKTmat, hecTKT)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: convert TtKT to hecTKT done', hecmw_wtime()-t1
    call hecmw_localmat_free(BTtKTmat)
    !
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: calculated TtKT', hecmw_wtime()-t1

    if (DEBUG >= 3) then
      write(1000+myrank,*) '======================================================================='
      call debug_write_vector(hecMAT%B, 'RHS(original)', 'hecMAT%B', ndof, hecMAT%N, &
          hecMAT%NP, .false., hecLagMAT%num_lagrange, n_slave, slaves)
      call debug_write_vector(conMAT%B, 'RHS(conMAT)', 'conMAT%B', ndof, conMAT%N, &
          conMAT%NP, .true., hecLagMAT%num_lagrange, n_slave, slaves)
    endif

    allocate(Btot(hecMAT%NP*ndof+hecLagMAT%num_lagrange))
    call assemble_b(hecMESH, hecMAT, conMAT, hecLagMAT%num_lagrange, Btot)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble conMAT%B done', hecmw_wtime()-t1
    if (DEBUG >= 3) call debug_write_vector(Btot, 'RHS(conMAT assembled)', 'Btot', ndof, conMAT%N, &
        conMAT%NP, .false., hecLagMAT%num_lagrange, n_slave, slaves)

    do i=1,hecMAT%N*ndof
      Btot(i)=Btot(i)+hecMAT%B(i)
    enddo
    if (DEBUG >= 3) call debug_write_vector(Btot, 'RHS(total)', 'Btot', ndof, conMAT%N, &
        conMAT%NP, .false., hecLagMAT%num_lagrange, n_slave, slaves)

    allocate(hecTKT%B(hecTKT%NP*ndof + hecLagMAT%num_lagrange))
    allocate(hecTKT%X(hecTKT%NP*ndof + hecLagMAT%num_lagrange))
    do i=1, hecTKT%N*ndof
      hecTKT%B(i) = Btot(i)
    enddo
    do ilag=1, hecLagMAT%num_lagrange
      hecTKT%B(hecTKT%NP*ndof + ilag) = hecMAT%B(hecMAT%NP*ndof + ilag)
    enddo
    do i=1, hecTKT%N*ndof
      hecTKT%X(i) = hecMAT%X(i)
    enddo
    do ilag=1,hecLagMAT%num_lagrange
      hecTKT%X(slaves_detected(ilag)) = 0.d0
    enddo

    ! make contact dof list
    call make_contact_dof_list(hecMAT, hecLagMAT, n_contact_dof, contact_dofs)

    ! make comm_table for contact dof
    call hecmw_contact_comm_init(conCOMM, hecMESH, ndof, n_contact_dof, contact_dofs)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make contact comm_table done', hecmw_wtime()-t1

    ! make RHS in converted system
    call make_converted_b(hecMAT, Btot, hecMESHtmp, hecTKT, BTtmat, BKmat, &
        slaves_detected, BLs_inv, hecLagMAT%num_lagrange, conCOMM, hecTKT%B)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: converted RHS', hecmw_wtime()-t1
    if (DEBUG >= 3) call debug_write_vector(hecTKT%B, 'RHS(converted)', 'hecTKT%B', ndof, hecTKT%N)

    ! solve converted system
    call hecmw_solve_iterative(hecMESHtmp,hecTKT)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: linear solver finished', hecmw_wtime()-t1
    if (DEBUG >= 3) call debug_write_vector(hecTKT%X, 'Solution(converted)', 'hecTKT%X', ndof, hecTKT%N)

    ! recover solution in original system
    hecMAT%Iarray=hecTKT%Iarray
    hecMAT%Rarray=hecTKT%Rarray

    call comp_x_slave(hecMAT, Btot, hecMESHtmp, hecTKT, BTmat, &
        hecLagMAT%num_lagrange, slaves_detected, BLs_inv, conCOMM, n_slave, slaves)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: recovered slave disp', hecmw_wtime()-t1

    call comp_lag(hecMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
        n_slave, slaves, hecLagMAT%num_lagrange, slaves_detected, BUs_inv, conCOMM)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: calculated lag', hecmw_wtime()-t1

    if (DEBUG >= 3) call debug_write_vector(hecMAT%X, 'Solution(original)', 'hecMAT%X', ndof, hecMAT%N, &
        hecMAT%NP, .false., hecLagMAT%num_lagrange, n_slave, slaves)

    ! check solution in the original system
    if (DEBUG >= 2) call check_solution(hecMESH, hecMAT, hecLagMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
          conCOMM, n_slave, slaves)

    ! free matrices
    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(BTtmat)
    call hecmw_mat_finalize(hecTKT)
    call hecmw_localmat_free(BKmat)
    call hecmw_contact_comm_finalize(conCOMM)
    call free_mesh(hecMESHtmp)
    deallocate(slaves_detected)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate end', hecmw_wtime()-t1
  end subroutine solve_eliminate

  !> \brief Choose slave dofs and compute BL_s^(-1) and BU_s^(-1)
  !>
  subroutine choose_slaves(hecMAT, hecLagMAT, dof_type, slaves_detected, BLs_inv, BUs_inv)
    type (hecmwST_matrix    ), intent(in) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer, intent(out) :: dof_type(:), slaves_detected(:)
    real(kind=kreal), intent(out) :: BLs_inv(:)
    real(kind=kreal), intent(out) :: BUs_inv(:)
    integer :: ndof, n, i, j, idof, jdof, l, ls, le, idx, imax, iwmin
    real(kind=kreal) :: val, vmax
    integer, allocatable :: iw1L(:), iw1U(:)
    integer(kind=kint) :: n_slave_in,n_slave_out
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    dof_type=-1

    if (hecLagMAT%num_lagrange == 0) return

    ndof=hecMAT%NDOF
    slaves_detected=0

    n = hecMAT%NP

    allocate(iw1L(n*ndof))
    allocate(iw1U(n*ndof))
    iw1L=0
    iw1U=0

    ! Count how many times each dof appear in Lagrange matrix
    ! lower
    do i=1,hecLagMAT%num_lagrange
      ls=hecLagMAT%indexL_lagrange(i-1)+1
      le=hecLagMAT%indexL_lagrange(i)
      do l=ls,le
        j=hecLagMAT%itemL_lagrange(l)
        do jdof=1,ndof
          idx=(j-1)*ndof+jdof
          iw1L(idx)=iw1L(idx)+1
        enddo
      enddo
    enddo
    ! upper
    do i=1,n
      ls=hecLagMAT%indexU_lagrange(i-1)+1
      le=hecLagMAT%indexU_lagrange(i)
      do l=ls,le
        j=hecLagMAT%itemU_lagrange(l)
        do idof=1,ndof
          idx=(i-1)*ndof+idof
          iw1U(idx)=iw1U(idx)+1
        enddo
      enddo
    enddo
    !!$    write(0,*) 'iw1L, iw1U:'
    !!$    do i=1,n*ndof
    !!$      if (iw1L(i) > 0 .or. iw1U(i) > 0) write(0,*) i, iw1L(i), iw1U(i)
    !!$    enddo

    ! Choose dofs that
    ! - appear only once in both lower and upper Lag. and
    ! - has greatest coefficient among them (in lower Lag.)
    do i=1,hecLagMAT%num_lagrange
      ls=hecLagMAT%indexL_lagrange(i-1)+1
      le=hecLagMAT%indexL_lagrange(i)
      vmax = 0.d0
      imax = -1
      iwmin = n
      do l=ls,le
        j=hecLagMAT%itemL_lagrange(l)
        do jdof=1,ndof
          idx=(j-1)*ndof+jdof
          val=hecLagMAT%AL_lagrange((l-1)*ndof+jdof)
          if ( iw1L(idx) < iwmin .and. iw1U(idx) < iwmin ) then
            iwmin = min(iw1L(idx),iw1U(idx))
            vmax = 0.d0
          endif
          if ( iw1L(idx) == iwmin .and. iw1U(idx) == iwmin ) then
            if ( abs(val) > abs(vmax) ) then
              imax=idx
              vmax=val
            endif
          endif
        enddo
      enddo
      if (imax == -1) stop "ERROR: iterative solver for contact failed"
      dof_type(imax)=i
      slaves_detected(i)=imax; BLs_inv(i)=1.d0/vmax
    enddo
    !!$    write(0,*) 'dof_type:'
    !!$    do i=1,n*ndof
    !!$      if (dof_type(i) > 0) write(0,*) i, dof_type(i), iw1L(i), iw1U(i)
    !!$    enddo
    !!$    write(0,*) 'slaves_detected:'
    !!$    write(0,*) slaves_detected(:)
    n_slave_in = 0
    do i=1,hecMAT%N*ndof
      if (dof_type(i) > 0) n_slave_in = n_slave_in + 1
    enddo
    n_slave_out = 0
    do i=hecMAT%N*ndof,n*ndof
      if (dof_type(i) > 0) n_slave_out = n_slave_out + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: n_slave(in,out,tot)',n_slave_in,n_slave_out,n_slave_in+n_slave_out

    deallocate(iw1L, iw1U)

    call make_BUs_inv(hecLagMAT, n, ndof, dof_type, BUs_inv)
  end subroutine choose_slaves

  !> \brief Compute BU_s^(-1)
  !>
  subroutine make_BUs_inv(hecLagMAT, n, ndof, dof_type, BUs_inv)
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: dof_type(:)
    real(kind=kreal), intent(out) :: BUs_inv(:)
    integer(kind=kint) :: i, idof, idx, js, je, j, k

    if (hecLagMAT%num_lagrange == 0) return

    BUs_inv=0.d0
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        if (dof_type(idx) > 0) then
          js=hecLagMAT%indexU_lagrange(i-1)+1
          je=hecLagMAT%indexU_lagrange(i)
          do j=js,je
            k=hecLagMAT%itemU_lagrange(j)
            if (k==dof_type(idx)) then
              BUs_inv(dof_type(idx)) = 1.0d0/hecLagMAT%AU_lagrange((j-1)*ndof+idof)
            endif
          enddo
        endif
      enddo
    enddo
    !write(0,*) BUs_inv
  end subroutine make_BUs_inv

  !> \brief Make transformation matrices
  !>
  subroutine make_transformation_matrices(hecMAT, hecLagMAT, hecMESH, hecMESHtmp, &
      dof_type, slaves_detected, BLs_inv, BUs_inv, mark, n_slave, slaves, BTmat, BTtmat)
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_mesh), intent(inout) :: hecMESHtmp
    integer(kind=kint), intent(in) :: dof_type(:)
    integer(kind=kint), intent(in) :: slaves_detected(:)
    real(kind=kreal), intent(out) :: BLs_inv(:)
    real(kind=kreal), intent(out) :: BUs_inv(:)
    integer(kind=kint), intent(out) :: mark(:)
    integer(kind=kint), intent(out) :: n_slave
    integer(kind=kint), allocatable, intent(out) :: slaves(:)
    type (hecmwST_local_matrix), intent(out) :: BTmat
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    call make_BTmat(hecMAT, hecLagMAT, dof_type, BLs_inv, BTmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make T done'
    if (DEBUG >= 4) call debug_write_matrix(BTmat, 'BTmat (local)')

    call make_BTtmat(hecMAT, hecLagMAT, dof_type, slaves_detected, BUs_inv, BTtmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make Tt done'
    if (DEBUG >= 4) call debug_write_matrix(BTtmat, 'BTtmat (local)')

    if ( hecmw_comm_get_size() > 1 ) then
      ! communicate and complete BTmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BTmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble T done'
      if (BTmat%nc /= hecMESH%n_node) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with T',BTmat%nc-hecMESH%n_node
      endif
      if (DEBUG >= 4) call debug_write_matrix(BTmat, 'BTmat (assembled)')

      ! communicate and complete BTtmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BTtmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble Tt done'
      if (BTtmat%nc /= BTmat%nc) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with BTtmat',BTtmat%nc-BTmat%nc
        BTmat%nc = BTtmat%nc
      endif
      if (DEBUG >= 4) call debug_write_matrix(BTtmat, 'BTtmat (assembled)')
    endif

    ! place 1 on diag of non-slave dofs of BTmat and BTtmat
    call mark_slave_dof(BTmat, mark, n_slave, slaves)
    call place_one_on_diag_of_unmarked_dof(BTmat, mark)
    call place_one_on_diag_of_unmarked_dof(BTtmat, mark)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: place 1 on diag of T and Tt done'
    if (DEBUG >= 4) call debug_write_matrix(BTmat, 'BTmat (1s on non-slave diag)')
    if (DEBUG >= 4) call debug_write_matrix(BTtmat, 'BTtmat (1s on non-slave diag)')
  end subroutine make_transformation_matrices

  !> \brief Make 3x3-blocked T matrix
  !>
  subroutine make_BTmat(hecMAT, hecLagMAT, dof_type, BLs_inv, BTmat)
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer, intent(in) :: dof_type(:)
    real(kind=kreal), intent(in) :: BLs_inv(:)
    type (hecmwST_local_matrix), intent(out) :: BTmat
    type (hecmwST_local_matrix) :: Tmat
    integer :: ndof, n, i, nnz, l, js, je, j, k, jdof, kk, jj
    real(kind=kreal) :: factor

    ndof=hecMAT%NDOF
    n=hecMAT%NP
    Tmat%nr=n*ndof
    Tmat%nc=Tmat%nr
    Tmat%nnz=hecLagMAT%numL_lagrange*ndof-hecLagMAT%num_lagrange
    Tmat%ndof=1

    allocate(Tmat%index(0:Tmat%nr))
    allocate(Tmat%item(Tmat%nnz), Tmat%A(Tmat%nnz))
    ! index
    Tmat%index(0)=0
    do i=1,Tmat%nr
      if (dof_type(i) > 0) then
        nnz=ndof*(hecLagMAT%indexL_lagrange(dof_type(i))-hecLagMAT%indexL_lagrange(dof_type(i)-1))-1
      else
        nnz = 0
      endif
      Tmat%index(i)=Tmat%index(i-1)+nnz
    enddo
    if (Tmat%nnz /= Tmat%index(Tmat%nr)) then
      write(0,*) Tmat%nnz, Tmat%index(Tmat%nr)
      Tmat%nnz = Tmat%index(Tmat%nr)
      !stop 'ERROR: Tmat%nnz wrong'
    endif
    ! item and A
    do i=1,Tmat%nr
      l=Tmat%index(i-1)+1
      if (dof_type(i) > 0) then
        js=hecLagMAT%indexL_lagrange(dof_type(i)-1)+1
        je=hecLagMAT%indexL_lagrange(dof_type(i))
        factor=-BLs_inv(dof_type(i))
        do j=js,je
          k=hecLagMAT%itemL_lagrange(j)
          do jdof=1,ndof
            kk=(k-1)*ndof+jdof
            jj=(j-1)*ndof+jdof
            if (kk==i) cycle
            Tmat%item(l)=kk
            Tmat%A(l)=hecLagMAT%AL_lagrange(jj)*factor
            l=l+1
          enddo
        enddo
      endif
      if (l /= Tmat%index(i)+1) then
        write(0,*) l, Tmat%index(i)+1
        stop 'ERROR: Tmat%index wrong'
      endif
    enddo
    !call hecmw_localmat_write(Tmat, 0)
    ! make 3x3-block version of Tmat
    call hecmw_localmat_blocking(Tmat, ndof, BTmat)
    call hecmw_localmat_free(Tmat)
  end subroutine make_BTmat

  !> \brief Make 3x3-blocked T^t matrix
  !>
  subroutine make_BTtmat(hecMAT, hecLagMAT, dof_type, slaves_detected, BUs_inv, BTtmat)
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer, intent(in) :: dof_type(:), slaves_detected(:)
    real(kind=kreal), intent(in) :: BUs_inv(:)
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    type (hecmwST_local_matrix) :: Ttmat
    integer :: ndof, n, i, nnz, l, js, je, j, k, idof, idx

    ndof=hecMAT%NDOF
    n=hecMAT%NP
    Ttmat%nr=n*ndof
    Ttmat%nc=Ttmat%nr
    Ttmat%nnz=hecLagMAT%numU_lagrange*ndof-hecLagMAT%num_lagrange
    Ttmat%ndof=1

    allocate(Ttmat%index(0:Ttmat%nr))
    allocate(Ttmat%item(Ttmat%nnz), Ttmat%A(Ttmat%nnz))
    ! index
    Ttmat%index(0)=0
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        if (dof_type(idx) <= 0) then
          if (hecLagMAT%num_lagrange == 0) then
            nnz=0
          else
            nnz=hecLagMAT%indexU_lagrange(i)-hecLagMAT%indexU_lagrange(i-1)
          endif
        else
          nnz=0
        endif
        Ttmat%index(idx)=Ttmat%index(idx-1)+nnz
      enddo
    enddo
    if (Ttmat%nnz /= Ttmat%index(Ttmat%nr)) then
      write(0,*) Ttmat%nnz, Ttmat%index(Ttmat%nr)
      !stop 'ERROR: Ttmat%nnz wrong'
      Ttmat%nnz = Ttmat%index(Ttmat%nr)
    endif
    ! item and A
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        l=Ttmat%index(idx-1)+1
        if (dof_type(idx) <= 0) then
          if (hecLagMAT%num_lagrange > 0) then
            ! offdiagonal
            js=hecLagMAT%indexU_lagrange(i-1)+1
            je=hecLagMAT%indexU_lagrange(i)
            do j=js,je
              k=hecLagMAT%itemU_lagrange(j)
              Ttmat%item(l)=slaves_detected(k)
              Ttmat%A(l)=-hecLagMAT%AU_lagrange((j-1)*ndof+idof)*BUs_inv(k)
              l=l+1
            enddo
          endif
        else
          ! no element
        endif
        if (l /= Ttmat%index(idx)+1) then
          write(0,*) l, Ttmat%index(idx)+1
          stop 'ERROR: Ttmat%index wrong'
        endif
      enddo
    enddo
    !call hecmw_localmat_write(Ttmat, 0)
    ! make 3x3-block version of Tmat
    call hecmw_localmat_blocking(Ttmat, ndof, BTtmat)
    call hecmw_localmat_free(Ttmat)
  end subroutine make_BTtmat

  !> \brief Make list of contact dofs including masters and slaves
  !>
  subroutine make_contact_dof_list(hecMAT, hecLagMAT, n_contact_dof, contact_dofs)
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    integer(kind=kint), intent(out) :: n_contact_dof
    integer(kind=kint), allocatable, intent(out) :: contact_dofs(:)
    integer(kind=kint) :: ndof, icnt, i, ls, le, l, j, jdof, idx, k, idof
    integer(kind=kint), allocatable :: iw(:)
    real(kind=kreal) :: val
    if (hecLagMAT%num_lagrange == 0) then
      n_contact_dof = 0
      return
    endif
    ! lower
    ndof = hecMAT%NDOF
    allocate(iw(hecMAT%NP * ndof))
    icnt = 0
    do i = 1, hecLagMAT%num_lagrange
      ls = hecLagMAT%indexL_lagrange(i-1)+1
      le = hecLagMAT%indexL_lagrange(i)
      do l = ls, le
        j = hecLagMAT%itemL_lagrange(l)
        ljdof1: do jdof = 1, ndof
          val = hecLagMAT%AL_lagrange((l-1)*ndof+jdof)
          !write(0,*) 'j,jdof,val',j,jdof,val
          if (abs(val) < tiny(0.0d0)) cycle
          idx = (j-1)*ndof+jdof
          do k = 1, icnt
            if (iw(k) == idx) cycle ljdof1
          enddo
          icnt = icnt + 1
          iw(icnt) = idx
          !write(0,*) 'icnt,idx',icnt,idx
        enddo ljdof1
      enddo
    enddo
    ! upper
    do i = 1, hecMAT%NP
      ls = hecLagMAT%indexU_lagrange(i-1)+1
      le = hecLagMAT%indexU_lagrange(i)
      do l = ls, le
        j = hecLagMAT%itemU_lagrange(l)
        lidof1: do idof = 1, ndof
          val = hecLagMAT%AU_lagrange((l-1)*ndof+idof)
          !write(0,*) 'i,idof,val',i,idof,val
          if (abs(val) < tiny(0.0d0)) cycle
          idx = (i-1)*ndof+idof
          do k = 1, icnt
            if (iw(k) == idx) cycle lidof1
          enddo
          icnt = icnt + 1
          iw(icnt) = idx
          !write(0,*) 'icnt,idx',icnt,idx
        enddo lidof1
      enddo
    enddo
    call quick_sort(iw, 1, icnt)
    allocate(contact_dofs(icnt))
    contact_dofs(1:icnt) = iw(1:icnt)
    n_contact_dof = icnt
    deallocate(iw)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: contact_dofs',contact_dofs(:)
  end subroutine make_contact_dof_list

  !> \brief Assemble conMAT%B into Btot
  !>
  subroutine assemble_b(hecMESH, hecMAT, conMAT, num_lagrange, Btot)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT, conMAT
    integer(kind=kint), intent(in) :: num_lagrange
    real(kind=kreal), intent(out) :: Btot(:)
    integer(kind=kint) :: ndof, nndof, npndof, i
    ndof = hecMAT%NDOF
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    do i=1,npndof+num_lagrange
      Btot(i) = conMAT%B(i)
    enddo
    call hecmw_assemble_R(hecMESH, Btot, hecMAT%NP, hecMAT%NDOF)
  end subroutine assemble_b

  !> \brief Make converted RHS vector
  !>
  subroutine make_converted_b(hecMAT, Btot, hecMESHtmp, hecTKT, BTtmat, BKmat, &
       slaves_detected, BLs_inv, num_lagrange, conCOMM, Bnew)
    type(hecmwST_local_mesh), intent(in) :: hecMESHtmp
    type(hecmwST_matrix), intent(in) :: hecMAT, hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type(hecmwST_local_matrix), intent(in) :: BTtmat, BKmat
    integer(kind=kint), intent(in) :: slaves_detected(:)
    real(kind=kreal), intent(in) :: BLs_inv(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_contact_comm), intent(in) :: conCOMM
    real(kind=kreal), intent(out) :: Bnew(:)
    real(kind=kreal), allocatable :: Btmp(:)
    integer(kind=kint) :: npndof, nndof, npndof_new, i
    ! SIZE:
    ! Btot <=> hecMAT, hecMESH
    ! Btmp <=> BKmat, hecMESHtmp
    ! Bnew <=> hecTKT, hecMESHtmp
    npndof     = hecMAT%NP*hecMAT%NDOF
    nndof      = hecMAT%N *hecMAT%NDOF
    npndof_new = hecTKT%NP*hecTKT%NDOF
    allocate(Btmp(npndof_new))
    !
    !! BTtmat*(B+K*(-Bs^-1)*Blag)
    !
    ! B2=-Bs^-1*Blag
    Bnew=0.d0
    do i=1,num_lagrange
      Bnew(slaves_detected(i))=-BLs_inv(i)*Btot(npndof+i)
    enddo
    ! send external contact dof => recv internal contact dof
    call hecmw_contact_comm_reduce_r(conCOMM, Bnew, HECMW_SUM)
    ! Btmp=B+K*B2 (including update of Bnew)
    call hecmw_update_R(hecMESHtmp, Bnew, hecTKT%NP, hecMAT%NDOF)
    call hecmw_localmat_mulvec(BKmat, Bnew, Btmp)
    do i=1,nndof
      Btmp(i)=Btot(i)+Btmp(i)
    enddo
    ! B2=BTtmat*Btmp
    call hecmw_update_R(hecMESHtmp, Btmp, hecTKT%NP, hecMAT%NDOF)
    call hecmw_localmat_mulvec(BTtmat, Btmp, Bnew)
    deallocate(Btmp)
  end subroutine make_converted_b

  !> \brief Recover solution for slave dofs
  !>
  subroutine comp_x_slave(hecMAT, Btot, hecMESHtmp, hecTKT, BTmat, &
       num_lagrange, slaves_detected, BLs_inv, conCOMM, n_slave, slaves)
    type (hecmwST_local_mesh), intent(in) :: hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), intent(in) :: hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint), intent(in) :: slaves_detected(:)
    real(kind=kreal), intent(in) :: BLs_inv(:)
    type (hecmwST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint), intent(in) :: n_slave
    integer(kind=kint), intent(in) :: slaves(:)
    integer(kind=kint) :: ndof, ndof2, npndof, nndof, ilag, islave, i
    real(kind=kreal), allocatable :: Xtmp(:)
    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N  * ndof
    !!
    !! {X} = [BT] {Xp} - [-Bs^-1] {c}
    !!
    ! compute {X} = [BT] {Xp}
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, 3)
    call hecmw_localmat_mulvec(BTmat, hecTKT%X, hecMAT%X)
    !
    ! compute {Xtmp} = [-Bs^-1] {c}
    allocate(Xtmp(npndof))
    Xtmp(:) = 0.0d0
    do ilag = 1, num_lagrange
      islave = slaves_detected(ilag)
      Xtmp(islave) = -BLs_inv(ilag) * Btot(npndof + ilag)
    enddo
    !
    ! send external contact dof => recv internal contact dof
    call hecmw_contact_comm_reduce_r(conCOMM, Xtmp, HECMW_SUM)
    !
    ! {X} = {X} - {Xtmp}
    do i = 1, n_slave
      islave = slaves(i)
      hecMAT%X(islave) = hecMAT%X(islave) - Xtmp(islave)
    enddo
    deallocate(Xtmp)
  end subroutine comp_x_slave

  !> \brief Recover solution for lagrange multipliers
  !>
  subroutine comp_lag(hecMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
       n_slave, slaves, num_lagrange, slaves_detected, BUs_inv, conCOMM)
    type (hecmwST_local_mesh), intent(in) :: hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT, hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BKmat
    integer(kind=kint), intent(in) :: n_slave, num_lagrange
    integer(kind=kint), intent(in) :: slaves(:), slaves_detected(:)
    real(kind=kreal), intent(in) :: BUs_inv(:)
    type (hecmwST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint) :: ndof, npndof, nndof, npndof_new, i, ilag, islave
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: xlag(:)
    ndof=hecMAT%ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N * ndof
    npndof_new = hecTKT%NP * ndof
    !! <SUMMARY>
    !! {lag} = [Bs^-T] ( {fs} - [Ksp Kss] {u} )
    !!
    ! 1. {Btmp} = [BKmat] {X}
    hecTKT%X(1:nndof) = hecMAT%X(1:nndof)
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, 3)
    allocate(Btmp(npndof))
    call hecmw_localmat_mulvec(BKmat, hecTKT%X, Btmp)
    !
    ! 2. {Btmp_s} = {fs} - {Btmp_s}
    do i = 1, n_slave
      islave = slaves(i)
      Btmp(islave) = Btot(islave) - Btmp(islave)
    enddo
    !
    ! 3. send internal contact dof => recv external contact dof
    call hecmw_contact_comm_bcast_r(conCOMM, Btmp)
    !
    ! 4. {lag} = [Bs^-T] {Btmp_s}
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    do ilag = 1, num_lagrange
      islave = slaves_detected(ilag)
      xlag(ilag)=BUs_inv(ilag)*Btmp(islave)
    enddo
    deallocate(Btmp)
  end subroutine comp_lag

  !> \brief Check solution in original system
  !>
  subroutine check_solution(hecMESH, hecMAT, hecLagMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
       conCOMM, n_slave, slaves)
    type (hecmwST_local_mesh), intent(in) :: hecMESH, hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT, hecTKT
    type (hecmwST_matrix_lagrange) , intent(in) :: hecLagMAT
    real(kind=kreal), target, intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BKmat
    type (hecmwST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint), intent(in) :: n_slave
    integer(kind=kint), intent(in) :: slaves(:)
    integer(kind=kint) :: ndof, nndof, npndof, num_lagrange, i, ls, le, l, j, idof, jdof
    real(kind=kreal), allocatable, target :: r(:)
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: rlag(:), blag(:), xlag(:)
    real(kind=kreal) :: rnrm2, rlagnrm2
    real(kind=kreal) :: bnrm2, blagnrm2
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    ndof = hecMAT%NDOF
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    num_lagrange = hecLagMAT%num_lagrange
    !
    allocate(r(npndof + num_lagrange))
    r(:) = 0.0d0
    allocate(Btmp(npndof))
    !
    rlag => r(npndof+1:npndof+num_lagrange)
    blag => Btot(npndof+1:npndof+num_lagrange)
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    !
    !! {r}    = {b} - [K] {x} - [Bt] {lag}
    !! {rlag} = {c} - [B] {x}
    !
    ! {r} = {b} - [K] {x}
    hecTKT%X(1:nndof) = hecMAT%X(1:nndof)
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, ndof)
    call hecmw_localmat_mulvec(BKmat, hecTKT%X, Btmp)
    r(1:nndof) = Btot(1:nndof) - Btmp(1:nndof)
    !
    ! {r} = {r} - [Bt] {lag}
    Btmp(:) = 0.0d0
    if (hecLagMAT%num_lagrange > 0) then
      do i = 1, hecMAT%NP
        ls = hecLagMAT%indexU_lagrange(i-1)+1
        le = hecLagMAT%indexU_lagrange(i)
        do l = ls, le
          j = hecLagMAT%itemU_lagrange(l)
          do idof = 1, ndof
            Btmp(ndof*(i-1)+idof) = Btmp(ndof*(i-1)+idof) + hecLagMAT%AU_lagrange(ndof*(l-1)+idof) * xlag(j)
          enddo
        enddo
      enddo
    endif
    call hecmw_contact_comm_reduce_r(conCOMM, Btmp, HECMW_SUM)
    r(1:nndof) = r(1:nndof) - Btmp(1:nndof)
    !
    ! {rlag} = {c} - [B] {x}
    call hecmw_update_R(hecMESH, hecMAT%X, hecMAT%NP, ndof)
    do i = 1, num_lagrange
      rlag(i) = blag(i)
      ls = hecLagMAT%indexL_lagrange(i-1)+1
      le = hecLagMAT%indexL_lagrange(i)
      do l = ls, le
        j = hecLagMAT%itemL_lagrange(l)
        do jdof = 1, ndof
          rlag(i) = rlag(i) - hecLagMAT%AL_lagrange(ndof*(l-1)+jdof) * hecMAT%X(ndof*(j-1)+jdof)
        enddo
      enddo
    enddo
    !
    ! residual in original system
    if (DEBUG >= 3) call debug_write_vector(r, 'Residual', 'R', ndof, hecMAT%N, &
        hecMAT%NP, .false., hecLagMAT%num_lagrange, n_slave, slaves)
    !
    call hecmw_InnerProduct_R(hecMESH, NDOF, r, r, rnrm2)
    call hecmw_InnerProduct_R(hecMESH, NDOF, Btot, Btot, bnrm2)
    rlagnrm2 = 0.0d0
    do i = 1, num_lagrange
      rlagnrm2 = rlagnrm2 + rlag(i)*rlag(i)
    enddo
    call hecmw_allreduce_R1(hecMESH, rlagnrm2, HECMW_SUM)
    blagnrm2 = 0.0d0
    do i = 1, num_lagrange
      blagnrm2 = blagnrm2 + blag(i)*blag(i)
    enddo
    call hecmw_allreduce_R1(hecMESH, blagnrm2, HECMW_SUM)
    !
    if (myrank == 0) then
      write(0,*) 'INFO: resid(x,lag,tot)',sqrt(rnrm2),sqrt(rlagnrm2),sqrt(rnrm2+rlagnrm2)
      write(0,*) 'INFO: rhs  (x,lag,tot)',sqrt(bnrm2),sqrt(blagnrm2),sqrt(bnrm2+blagnrm2)
    endif
  end subroutine check_solution

  !> \brief Find internal dofs that are in contact in some subdomain
  !>
  subroutine mark_slave_dof(BTmat, mark, n_slave, slaves)
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(out) :: mark(:)
    integer(kind=kint), intent(out) :: n_slave
    integer(kind=kint), allocatable, intent(out) :: slaves(:)
    integer(kind=kint) :: ndof, ndof2, irow, js, je, j, jcol, idof, jdof, i
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    mark(:) = 0
    do irow = 1, BTmat%nr
      js = BTmat%index(irow-1)+1
      je = BTmat%index(irow)
      do j = js, je
        jcol = BTmat%item(j)
        do idof = 1, ndof
          if (mark(ndof*(irow-1)+idof) == 1) cycle
          do jdof = 1, ndof
            if (irow == jcol .and. idof == jdof) cycle
            if (abs(BTmat%A(ndof2*(j-1)+ndof*(idof-1)+jdof)) > tiny(0.0d0)) then
              mark(ndof*(irow-1)+idof) = 1
              exit
            endif
          enddo
        enddo
      enddo
    enddo
    n_slave = 0
    do i = 1, BTmat%nr * ndof
      if (mark(i) /= 0) n_slave = n_slave + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: n_slave',n_slave
    allocate(slaves(n_slave))
    n_slave = 0
    do i = 1, BTmat%nr * ndof
      if (mark(i) /= 0) then
        n_slave = n_slave + 1
        slaves(n_slave) = i
      endif
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: slaves',slaves(:)
  end subroutine mark_slave_dof

  !> \brief Place one on diagonal of unmarked dof of T or T^t matrix
  !>
  subroutine place_one_on_diag_of_unmarked_dof(BTmat, mark)
    type (hecmwST_local_matrix), intent(inout) :: BTmat
    integer(kind=kint), intent(in) :: mark(:)
    type (hecmwST_local_matrix) :: Imat, Wmat
    integer(kind=kint) :: ndof, ndof2, i, irow, idof, n_slave, n_other
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    ! Imat: unit matrix except for slave dofs
    Imat%nr = BTmat%nr
    Imat%nc = BTmat%nc
    Imat%nnz = Imat%nr
    Imat%ndof = ndof
    allocate(Imat%index(0:Imat%nr))
    allocate(Imat%item(Imat%nnz))
    Imat%index(0) = 0
    do i = 1, Imat%nr
      Imat%index(i) = i
      Imat%item(i) = i
    enddo
    allocate(Imat%A(ndof2 * Imat%nnz))
    Imat%A(:) = 0.0d0
    n_slave = 0
    n_other = 0
    do irow = 1, Imat%nr
      do idof = 1, ndof
        if (mark(ndof*(irow-1)+idof) == 0) then
          Imat%A(ndof2*(irow-1)+ndof*(idof-1)+idof) = 1.0d0
          n_other = n_other + 1
        else
          n_slave = n_slave + 1
        endif
      enddo
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: n_slave,n_other',n_slave,n_other
    call hecmw_localmat_add(BTmat, Imat, Wmat)
    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(Imat)
    BTmat%nr = Wmat%nr
    BTmat%nc = Wmat%nc
    BTmat%nnz = Wmat%nnz
    BTmat%ndof = Wmat%ndof
    BTmat%index => Wmat%index
    BTmat%item => Wmat%item
    BTmat%A => Wmat%A
  end subroutine place_one_on_diag_of_unmarked_dof

  !> \brief Place num on diagonal of marked dof of T^tKT matrix
  !>
  subroutine place_num_on_diag_of_marked_dof(BTtKTmat, num, mark)
    type (hecmwST_local_matrix), intent(inout) :: BTtKTmat
    real(kind=kreal), intent(in) :: num
    integer(kind=kint), intent(in) :: mark(:)
    integer(kind=kint) :: ndof, ndof2, irow, js, je, j, jcol, idof
    integer(kind=kint) :: n_error = 0
    ndof = BTtKTmat%ndof
    ndof2 = ndof*ndof
    do irow = 1, BTtKTmat%nr
      js = BTtKTmat%index(irow-1)+1
      je = BTtKTmat%index(irow)
      do j = js, je
        jcol = BTtKTmat%item(j)
        if (irow /= jcol) cycle
        do idof = 1, ndof
          if (mark(ndof*(irow-1)+idof) == 1) then
            if (abs(BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof)) > tiny(0.0d0)) &
                 stop 'ERROR: nonzero diag on slave dof of BTtKTmat'
            BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof) = num
          else
            if (abs(BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof)) < tiny(0.0d0)) then
              !write(0,*) 'irow,idof',irow,idof
              n_error = n_error+1
            endif
          endif
        enddo
      enddo
    enddo
    if (n_error > 0) then
      write(0,*) 'n_error',n_error
      stop 'ERROR: zero diag on non-slave dof of BTtKTmat'
    endif
  end subroutine place_num_on_diag_of_marked_dof

  !> \brief Copy mesh
  !>
  subroutine copy_mesh(src, dst)
    type (hecmwST_local_mesh), intent(in) :: src
    type (hecmwST_local_mesh), intent(out) :: dst
    dst%zero          = src%zero
    dst%MPI_COMM      = src%MPI_COMM
    dst%PETOT         = src%PETOT
    dst%PEsmpTOT      = src%PEsmpTOT
    dst%my_rank       = src%my_rank
    dst%n_subdomain   = src%n_subdomain
    dst%n_node        = src%n_node
    dst%nn_internal   = src%nn_internal
    dst%n_elem        = src%n_elem
    dst%ne_internal   = src%ne_internal
    dst%n_elem_type   = src%n_elem_type
    dst%n_dof         = src%n_dof
    dst%n_neighbor_pe = src%n_neighbor_pe
    if (src%n_neighbor_pe > 0) then
      allocate(dst%neighbor_pe(dst%n_neighbor_pe))
      dst%neighbor_pe(:) = src%neighbor_pe(:)
      allocate(dst%import_index(0:dst%n_neighbor_pe))
      allocate(dst%export_index(0:dst%n_neighbor_pe))
      dst%import_index(:)= src%import_index(:)
      dst%export_index(:)= src%export_index(:)
      allocate(dst%import_item(dst%import_index(dst%n_neighbor_pe)))
      dst%import_item(1:dst%import_index(dst%n_neighbor_pe)) = src%import_item(1:dst%import_index(dst%n_neighbor_pe))
      allocate(dst%export_item(dst%export_index(dst%n_neighbor_pe)))
      dst%export_item(1:dst%export_index(dst%n_neighbor_pe)) = src%export_item(1:dst%export_index(dst%n_neighbor_pe))
      allocate(dst%global_node_ID(dst%n_node))
      dst%global_node_ID(1:dst%n_node) = src%global_node_ID(1:dst%n_node)
    else
      dst%neighbor_pe => null()
      dst%import_index => null()
      dst%export_index => null()
      dst%import_item => null()
      dst%export_item => null()
      dst%global_node_ID => src%global_node_ID
    endif
    allocate(dst%node_ID(2*dst%n_node))
    dst%node_ID(1:2*dst%n_node) = src%node_ID(1:2*dst%n_node)
    allocate(dst%elem_type_item(dst%n_elem_type))
    dst%elem_type_item(:) = src%elem_type_item(:)
    !dst%mpc            = src%mpc
    ! MPC is already set outside of here
    dst%mpc%n_mpc = 0
    dst%node => src%node
  end subroutine copy_mesh

  !> \brief Free mesh
  !>
  subroutine free_mesh(hecMESH)
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    if (hecMESH%n_neighbor_pe > 0) then
      deallocate(hecMESH%neighbor_pe)
      deallocate(hecMESH%import_index)
      deallocate(hecMESH%export_index)
      deallocate(hecMESH%import_item)
      deallocate(hecMESH%export_item)
      deallocate(hecMESH%global_node_ID)
    endif
    deallocate(hecMESH%node_ID)
    deallocate(hecMESH%elem_type_item)
    !hecMESH%node => null()
  end subroutine free_mesh

  !> \brief Quick sort for integer array
  !>
  recursive subroutine quick_sort(array, id1, id2)
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: id1, id2
    integer(kind=kint) :: pivot, center, left, right, tmp
    if (id1 >= id2) return
    center = (id1 + id2) / 2
    pivot = array(center)
    left = id1
    right = id2
    do
      do while (array(left) < pivot)
        left = left + 1
      end do
      do while (pivot < array(right))
        right = right - 1
      end do
      if (left >= right) exit
      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp
      left = left + 1
      right = right - 1
    end do
    if (id1 < left-1) call quick_sort(array, id1, left-1)
    if (right+1 < id2) call quick_sort(array, right+1, id2)
    return
  end subroutine quick_sort

  !> \brief Debug write matrix
  !>
  subroutine debug_write_matrix(Mat, label)
    type (hecmwST_local_matrix), intent(in) :: Mat
    character(len=*), intent(in) :: label
    integer(kind=kint) :: myrank
    myrank = hecmw_comm_get_rank()
    write(1000+myrank,*) trim(label)
    call hecmw_localmat_write(Mat, 1000+myrank)
  end subroutine debug_write_matrix

  !> \brief Debug write vector
  !>
  subroutine debug_write_vector(Vec, label, name, ndof, N, &
      NP, write_ext, num_lagrange, n_slave, slaves)
    real(kind=kreal),   intent(in) :: Vec(:)
    character(len=*),   intent(in) :: label, name
    integer(kind=kint), intent(in) :: ndof, N
    integer(kind=kint), intent(in), optional :: NP
    logical,            intent(in), optional :: write_ext
    integer(kind=kint), intent(in), optional :: num_lagrange
    integer(kind=kint), intent(in), optional :: n_slave
    integer(kind=kint), intent(in), optional :: slaves(:)
    integer(kind=kint) :: myrank
    myrank = hecmw_comm_get_rank()
    write(1000+myrank,*) label,'------------------------------------------------------------'
    write(1000+myrank,*) 'size of ',trim(name),size(Vec)
    write(1000+myrank,*) trim(name),': 1-',N*ndof
    write(1000+myrank,*) Vec(1:N*ndof)
    if (present(write_ext) .and. present(NP)) then
      if (write_ext) then
        write(1000+myrank,*) trim(name),'(external): ',N*ndof+1,'-',NP*ndof
        write(1000+myrank,*) Vec(N*ndof+1:NP*ndof)
      endif
    endif
    if (present(num_lagrange) .and. present(NP)) then
      if (num_lagrange > 0) then
        write(1000+myrank,*) trim(name),'(lag):',NP*ndof+1,'-',NP*ndof+num_lagrange
        write(1000+myrank,*) Vec(NP*ndof+1:NP*ndof+num_lagrange)
      endif
    endif
    if (present(n_slave) .and. present(slaves)) then
      if (n_slave > 0) then
        write(1000+myrank,*) trim(name),'(slave):',slaves(:)
        write(1000+myrank,*) Vec(slaves(:))
      endif
    endif
  end subroutine debug_write_vector

end module m_solve_LINEQ_iter_contact
