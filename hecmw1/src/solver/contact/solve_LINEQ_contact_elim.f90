!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides interface of iteratie linear equation solver for
!! contact problems using Lagrange multiplier.
module m_solve_LINEQ_contact_elim
  use hecmw_util
  use hecmw_local_matrix
  use m_hecmw_contact_comm
  use hecmw_solver
  use hecmw_matrix_misc
  use m_hecmw_comm_f
  use hecmw_solver_misc
  use hecmw_solver_las
  use hecmw_mpc_prepost

  implicit none

  private
  public :: solve_LINEQ_contact_elim_init
  public :: solve_LINEQ_contact_elim

  logical, save :: INITIALIZED = .false.
  integer, save :: SymType = 0

  integer, parameter :: DEBUG = 0  ! 0: no message, 1: some messages, 2: more messages, 3: even more messages
  logical, parameter :: DEBUG_VECTOR = .false.
  logical, parameter :: DEBUG_MATRIX = .false.

contains

  subroutine solve_LINEQ_contact_elim_init(hecMESH, hecMAT, hecLagMAT, is_sym)
    type(hecmwST_local_mesh),      intent(in)    :: hecMESH   !< mesh
    type(hecmwST_matrix),          intent(inout) :: hecMAT    !< matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(in)    :: hecLagMAT !< matrix for lagrange multipliers
    logical,                       intent(in)    :: is_sym    !< if matrix is symmetric or not

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
  end subroutine solve_LINEQ_contact_elim_init

  subroutine solve_LINEQ_contact_elim(hecMESH, hecMAT, hecLagMAT, istat, conMAT, is_contact_active)
    type(hecmwST_local_mesh),      intent(inout) :: hecMESH           !< mesh
    type(hecmwST_matrix),          intent(inout) :: hecMAT            !< matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT         !< matrix for lagrange multipliers
    integer(kind=kint),            intent(out)   :: istat             !< status
    type(hecmwST_matrix),          intent(in)    :: conMAT            !< matrix for contact
    logical,                       intent(in)    :: is_contact_active !< if contact is active or not
    !
    integer(kind=kint) :: solver_type, method_org
    integer(kind=kint) :: is_contact
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    hecMAT%Iarray(97) = 1

    is_contact = 0
    if (is_contact_active) is_contact = 1
    call hecmw_allreduce_I1(hecMESH, is_contact, hecmw_max)

    if (is_contact == 0) then
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: no contact'
      solver_type = hecmw_mat_get_solver_type(hecMAT)
      if (solver_type == 1) then
        ! use CG because the matrix is symmetric
        method_org = hecmw_mat_get_method(hecMAT)
        call hecmw_mat_set_method(hecMAT, 1)
      endif
      ! solve
      call solve_with_MPC(hecMESH, hecMAT)
      if (solver_type == 1) then
        ! restore solver setting
        call hecmw_mat_set_method(hecMAT, method_org)
      endif
    else
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: with contact'
      call solve_eliminate(hecMESH, hecMAT, hecLagMAT, conMAT)
    endif

    istat = hecmw_mat_get_flag_diverged(hecMAT)
  end subroutine solve_LINEQ_contact_elim

  subroutine solve_with_MPC(hecMESH, hecMAT)
    type(hecmwST_local_mesh),      intent(inout) :: hecMESH           !< mesh
    type(hecmwST_matrix),          intent(inout) :: hecMAT            !< matrix excl. contact
    !
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: method
    logical            :: fg_cg, fg_amg

    fg_cg = (hecmw_mat_get_method(hecMAT) == 1)
    fg_amg = (hecmw_mat_get_precond(hecMAT) == 5)

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    call hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)
    call hecmw_solve(hecMESHmpc,hecMATmpc)
    if (fg_cg .and. fg_amg .and. hecmw_mat_get_flag_diverged(hecMATmpc) /= 0) then
      ! avoid ML and retry when diverged
      call hecmw_mat_set_precond(hecMATmpc, 3) ! set diag-scaling
      hecMATmpc%Iarray(97:98) = 1
      hecMATmpc%X(:) = 0.d0
      call hecmw_solve(hecMESHmpc,hecMATmpc)
      call hecmw_mat_set_precond(hecMATmpc, 5) ! restore amg
    endif
    call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)
    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
  end subroutine solve_with_MPC

  !> \brief Solve with elimination of Lagrange-multipliers
  !>
  subroutine solve_eliminate(hecMESH,hecMAT,hecLagMAT,conMAT)
    type(hecmwST_local_mesh),      intent(inout) :: hecMESH   !< original mesh
    type(hecmwST_matrix),          intent(inout) :: hecMAT    !< original matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< original matrix for lagrange multipliers
    type(hecmwST_matrix),          intent(in)    :: conMAT    !< original matrix for contact
    !
    type(hecmwST_local_mesh)        :: hecMESHtmp     !< temoprary copy of mesh for migrating nodes
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint), allocatable :: slaves4lag(:)  !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),   allocatable :: BLs_inv(:)     !< inverse of diagonal BLs matrix
    real(kind=kreal),   allocatable :: BUs_inv(:)     !< inverse of diagonal BUs matrix
    type(hecmwST_local_matrix)      :: Tmat           !< blocked T matrix
    type(hecmwST_local_matrix)      :: Ttmat          !< blocked T^t matrix
    integer(kind=kint), allocatable :: slaves(:)      !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    type(hecmwST_contact_comm)      :: conCOMM        !< contact comm table for optimized communication
    type(hecmwST_local_matrix)      :: Kmat           !< total K matrix assembled from hecMAT and conMAT
    real(kind=kreal),   allocatable :: Btot(:)        !< total RHS vector assembled from hecMAT%B and conMAT%B
    type(hecmwST_matrix)            :: hecTKT         !< converted system's matrix and RHS
    integer(kind=kint)              :: ndof
    integer(kind=kint)              :: myrank
    real(kind=kreal)                :: t0, t1, t2

    myrank = hecmw_comm_get_rank()

    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate start'
    t0 = hecmw_wtime()

    ndof=hecMAT%NDOF
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: num_lagrange',hecLagMAT%num_lagrange

    call copy_mesh(hecMESH, hecMESHtmp)

    allocate(slaves4lag(hecLagMAT%num_lagrange), BLs_inv(hecLagMAT%num_lagrange), &
      BUs_inv(hecLagMAT%num_lagrange))

    t1 = hecmw_wtime()
    call make_transformation_matrices(hecMESH, hecMESHtmp, hecMAT, hecLagMAT, &
        slaves4lag, BLs_inv, BUs_inv, slaves, Tmat, Ttmat)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: made trans matrices', t2-t1

    t1 = t2
    call make_contact_comm_table(hecMESH, hecMAT, hecLagMAT, conCOMM)
    t2 = hecmw_wtime()
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make contact comm_table done', hecmw_wtime()-t1

    t1 = t2
    allocate(Btot(hecMAT%NP*ndof+hecLagMAT%num_lagrange))
    call assemble_equation(hecMESH, hecMESHtmp, hecMAT, conMAT, hecLagMAT%num_lagrange, &
        slaves, Kmat, Btot)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: assembled equation ', t2-t1

    if (hecmw_comm_get_size() > 1) then
      if (Kmat%nc /= Ttmat%nc) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with Kmat',Kmat%nc-Ttmat%nc
        Tmat%nc = Kmat%nc
        Ttmat%nc = Kmat%nc
      endif
    endif

    t1 = t2
    call convert_equation(hecMESHtmp, hecMAT, Kmat, Tmat, Ttmat, Btot, slaves, &
        slaves4lag, BLs_inv, conCOMM, hecTKT)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: converted equation ', t2-t1

    t1 = t2
    call solve_with_MPC(hecMESHtmp, hecTKT)
    if (DEBUG_VECTOR) call debug_write_vector(hecTKT%X, 'Solution(converted)', 'hecTKT%X', ndof, hecTKT%N)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: linear solver done ', t2-t1

    t1 = t2
    call recover_solution(hecMESHtmp, hecMAT, hecTKT, Tmat, Kmat, Btot, &
        slaves4lag, BLs_inv, BUs_inv, conCOMM, slaves)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: recovered solution ', t2-t1

    if (DEBUG >= 1) call check_solution(hecMESH, hecMESHtmp, hecMAT, hecTKT, hecLagMAT, Kmat, Btot, &
       conCOMM, slaves)
    if (DEBUG >= 2) call check_solution2(hecMESH, hecMAT, conMAT, hecLagMAT, conCOMM, slaves)

    call hecmw_localmat_free(Tmat)
    call hecmw_localmat_free(Ttmat)
    call hecmw_mat_finalize(hecTKT)
    call hecmw_localmat_free(Kmat)
    call hecmw_contact_comm_finalize(conCOMM)
    call free_mesh(hecMESHtmp)
    deallocate(slaves4lag)
    t2 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate end', t2-t0
  end subroutine solve_eliminate

  !> \brief Copy mesh
  !>
  subroutine copy_mesh(src, dst)
    type(hecmwST_local_mesh), intent(in)  :: src !< original mesh
    type(hecmwST_local_mesh), intent(out) :: dst !< copy of the mesh

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
      dst%import_index(:)= src%import_index(:)
      allocate(dst%export_index(0:dst%n_neighbor_pe))
      dst%export_index(:)= src%export_index(:)
      allocate(dst%import_item(dst%import_index(dst%n_neighbor_pe)))
      dst%import_item(1:dst%import_index(dst%n_neighbor_pe)) = src%import_item(1:dst%import_index(dst%n_neighbor_pe))
      allocate(dst%export_item(dst%export_index(dst%n_neighbor_pe)))
      dst%export_item(1:dst%export_index(dst%n_neighbor_pe)) = src%export_item(1:dst%export_index(dst%n_neighbor_pe))
    else
      dst%neighbor_pe => null()
      dst%import_index => null()
      dst%export_index => null()
      dst%import_item => null()
      dst%export_item => null()
    endif
    allocate(dst%global_node_ID(dst%n_node))
    dst%global_node_ID(1:dst%n_node) = src%global_node_ID(1:dst%n_node)
    allocate(dst%node_ID(2*dst%n_node))
    dst%node_ID(1:2*dst%n_node) = src%node_ID(1:2*dst%n_node)
    allocate(dst%elem_type_item(dst%n_elem_type))
    dst%elem_type_item(:) = src%elem_type_item(:)
    !
    dst%mpc%n_mpc =  src%mpc%n_mpc
    dst%mpc%mpc_index => src%mpc%mpc_index
    dst%mpc%mpc_item => src%mpc%mpc_item
    dst%mpc%mpc_dof => src%mpc%mpc_dof
    dst%mpc%mpc_val => src%mpc%mpc_val
    dst%mpc%mpc_const => src%mpc%mpc_const
    !
    dst%node_group%n_grp = src%node_group%n_grp
    dst%node_group%n_bc = src%node_group%n_bc
    dst%node_group%grp_name => src%node_group%grp_name
    dst%node_group%grp_index => src%node_group%grp_index
    dst%node_group%grp_item => src%node_group%grp_item
    dst%node_group%bc_grp_ID => src%node_group%bc_grp_ID
    dst%node_group%bc_grp_type => src%node_group%bc_grp_type
    dst%node_group%bc_grp_index => src%node_group%bc_grp_index
    dst%node_group%bc_grp_dof => src%node_group%bc_grp_dof
    dst%node_group%bc_grp_val => src%node_group%bc_grp_val
    !
    dst%node => src%node
  end subroutine copy_mesh

  !> \brief Free mesh
  !>
  subroutine free_mesh(hecMESH)
    type(hecmwST_local_mesh), intent(inout) :: hecMESH  !< mesh

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

  !> \brief Make transformation matrices T and T^t
  !>
  subroutine make_transformation_matrices(hecMESH, hecMESHtmp, hecMAT, hecLagMAT, &
      slaves4lag, BLs_inv, BUs_inv, slaves, Tmat, Ttmat)
    type(hecmwST_local_mesh),        intent(in)    :: hecMESH       !< original mesh
    type(hecmwST_local_mesh),        intent(inout) :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),            intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix_lagrange),   intent(inout) :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),              intent(out)   :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),                intent(out)   :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    real(kind=kreal),                intent(out)   :: BUs_inv(:)    !< inverse of diagonal BUs matrix
    integer(kind=kint), allocatable, intent(out)   :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    type(hecmwST_local_matrix),      intent(out)   :: Tmat          !< blocked T matrix
    type(hecmwST_local_matrix),      intent(out)   :: Ttmat         !< blocked T^t matrix
    !
    integer(kind=kint) :: myrank, n

    myrank = hecmw_comm_get_rank()
    n = hecMAT%NP

    ! choose slave DOFs to be eliminated with Lag. DOFs
    call choose_slaves(hecMAT, hecLagMAT, n, slaves4lag)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: slave DOFs chosen'

    call make_BLs_inv(hecLagMAT, hecMAT%NDOF, slaves4lag, BLs_inv)
    call make_BUs_inv(hecLagMAT, hecMAT%NDOF, slaves4lag, BUs_inv)

    call add_C_to_Tmat(hecMAT, hecLagMAT, n, slaves4lag, BLs_inv, Tmat)
    if (DEBUG_MATRIX) call debug_write_matrix(Tmat, 'Tmat (local, C only)')
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: add C to Tmat done'

    call add_Ct_to_Ttmat(hecMAT, hecLagMAT, n, slaves4lag, BUs_inv, Ttmat)
    if (DEBUG_MATRIX) call debug_write_matrix(Ttmat, 'Ttmat (local, Ct only)')
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: add Ct to Tt done'

    if (hecmw_comm_get_size() > 1) then
      ! communicate and assemble Tmat (updating hecMESHtmp)
      call hecmw_localmat_assemble(Tmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) then
        write(0,*) '  DEBUG2[',myrank,']: assemble T done'
        if (Tmat%nc /= hecMESH%n_node) write(0,*) '  DEBUG2[',myrank,']: node migrated with T',Tmat%nc-hecMESH%n_node
      endif
      if (DEBUG_MATRIX) call debug_write_matrix(Tmat, 'Tmat (assembled, C only)')

      ! communicate and assemble Ttmat (updating hecMESHtmp)
      call hecmw_localmat_assemble(Ttmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) then
        write(0,*) '  DEBUG2[',myrank,']: assemble Tt done'
        if (Ttmat%nc /= Tmat%nc) write(0,*) '  DEBUG2[',myrank,']: node migrated with Ttmat',Ttmat%nc-Tmat%nc
        Tmat%nc = Ttmat%nc
      endif
      if (DEBUG_MATRIX) call debug_write_matrix(Ttmat, 'Ttmat (assembled, Ct only)')
    endif

    ! add Ip to Tmat and Ttmat
    call make_slave_list(hecMESHtmp, Tmat%ndof, slaves4lag, slaves)
    call add_Ip_to_Tmat(Tmat, slaves)
    if (DEBUG_MATRIX) call debug_write_matrix(Tmat, 'Tmat (final)')
    call add_Ip_to_Tmat(Ttmat, slaves)
    if (DEBUG_MATRIX) call debug_write_matrix(Ttmat, 'Ttmat (final)')
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: place 1 on diag of T and Tt done'
  end subroutine make_transformation_matrices

  !> \brief Choose slave dofs and compute BL_s^(-1) and BU_s^(-1)
  !>
  subroutine choose_slaves(hecMAT, hecLagMAT, n, slaves4lag)
    type(hecmwST_matrix),          intent(in)  :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(in)  :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),            intent(in)  :: n             !< num of nodes in this subdomain
    integer(kind=kint),            intent(out) :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    !
    integer(kind=kint) :: ndof, i, j, idof, jdof, l, ls, le, idx, imax, iwmin
    real(kind=kreal)   :: val, vmax
    integer(kind=kint), allocatable :: mark_slave4lag(:)
    integer(kind=kint), allocatable :: iw1L(:), iw1U(:)
    integer(kind=kint) :: n_slave_in, n_slave_out, ilag
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    ndof=hecMAT%NDOF

    allocate(mark_slave4lag(n*ndof), source=0)
    slaves4lag=0

    if (hecLagMAT%num_lagrange == 0) return

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
          if (iw1L(idx) < iwmin .and. iw1U(idx) < iwmin) then
            iwmin = min(iw1L(idx),iw1U(idx))
            vmax = 0.d0
          endif
          if (iw1L(idx) == iwmin .and. iw1U(idx) == iwmin) then
            if (abs(val) > abs(vmax)) then
              imax=idx
              vmax=val
            endif
          endif
        enddo
      enddo
      if (imax == -1) stop "ERROR: iterative solver for contact failed"
      mark_slave4lag(imax)=i
      slaves4lag(i)=imax
    enddo
    !!$    write(0,*) 'mark_slave4lag:'
    !!$    do i=1,n*ndof
    !!$      if (mark_slave4lag(i) > 0) write(0,*) i, mark_slave4lag(i), iw1L(i), iw1U(i)
    !!$    enddo
    !!$    write(0,*) 'slaves4lag:'
    !!$    write(0,*) slaves4lag(:)
    if (DEBUG >= 2) then
      n_slave_in = 0
      n_slave_out = 0
      do ilag=1,hecLagMAT%num_lagrange
        i = slaves4lag(ilag)
        if (0 < i .and. i <= hecMAT%N*ndof) then
          n_slave_in = n_slave_in + 1
        elseif (hecMAT%N*ndof < i .and. i <= hecMAT%NP*ndof) then
          n_slave_out = n_slave_out + 1
        endif
      enddo
      write(0,*) '  DEBUG2[',myrank,']: n_slave(in,out,tot)',n_slave_in,n_slave_out,hecLagMAT%num_lagrange
    endif

    deallocate(mark_slave4lag)
    deallocate(iw1L, iw1U)
  end subroutine choose_slaves

  !> \brief Compute BL_s^(-1)
  !>
  subroutine make_BLs_inv(hecLagMAT, ndof, slaves4lag, BLs_inv)
    type(hecmwST_matrix_lagrange), intent(in)  :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),            intent(in)  :: ndof          !< num of DOF per node
    integer(kind=kint),            intent(in)  :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),              intent(out) :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    !
    integer(kind=kint) :: ilag, ls, le, l, j, jdof, idx

    if (hecLagMAT%num_lagrange == 0) return

    BLs_inv=0.d0
    do ilag=1,hecLagMAT%num_lagrange
      ls=hecLagMAT%indexL_lagrange(ilag-1)+1
      le=hecLagMAT%indexL_lagrange(ilag)
      lloop: do l=ls,le
        j=hecLagMAT%itemL_lagrange(l)
        do jdof=1,ndof
          idx=(j-1)*ndof+jdof
          if (idx==slaves4lag(ilag)) then
            BLs_inv(ilag) = 1.0d0/hecLagMAT%AL_lagrange((l-1)*ndof+jdof)
            exit lloop
          endif
        enddo
      enddo lloop
    enddo
    !write(0,*) BLs_inv
  end subroutine make_BLs_inv

  !> \brief Compute BU_s^(-1)
  !>
  subroutine make_BUs_inv(hecLagMAT, ndof, slaves4lag, BUs_inv)
    type(hecmwST_matrix_lagrange), intent(in)  :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),            intent(in)  :: ndof          !< num of DOF per node
    integer(kind=kint),            intent(in)  :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),              intent(out) :: BUs_inv(:)    !< inverse of diagonal BUs matrix
    !
    integer(kind=kint) :: ilag, i, idof, js, je, j, k

    if (hecLagMAT%num_lagrange == 0) return

    BUs_inv=0.d0
    do ilag=1,size(slaves4lag)
      i=(slaves4lag(ilag)+ndof-1)/ndof
      idof=slaves4lag(ilag)-(i-1)*ndof
      js=hecLagMAT%indexU_lagrange(i-1)+1
      je=hecLagMAT%indexU_lagrange(i)
      do j=js,je
        k=hecLagMAT%itemU_lagrange(j)
        if (k==ilag) then
          BUs_inv(ilag) = 1.0d0/hecLagMAT%AU_lagrange((j-1)*ndof+idof)
          exit
        endif
      enddo
    enddo
    !write(0,*) BUs_inv
  end subroutine make_BUs_inv

  !> \brief Make 3x3-blocked T matrix
  !>
  subroutine add_C_to_Tmat(hecMAT, hecLagMAT, n, slaves4lag, BLs_inv, Tmat)
    type(hecmwST_matrix),          intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),            intent(in)    :: n             !< num of nodes in this subdomain
    integer(kind=kint),            intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),              intent(in)    :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    type(hecmwST_local_matrix),    intent(out)   :: Tmat          !< blocked T matrix
    !
    type(hecmwST_local_matrix) :: Tmat11
    integer(kind=kint), allocatable :: nz_cnt(:)
    integer(kind=kint) :: ndof, i, ilag, l, js, je, j, k, jdof, kk, jj
    real(kind=kreal)   :: factor

    ndof=hecMAT%NDOF
    Tmat11%nr=n*ndof
    Tmat11%nc=Tmat11%nr
    Tmat11%nnz=hecLagMAT%numL_lagrange*ndof-hecLagMAT%num_lagrange
    Tmat11%ndof=1

    allocate(Tmat11%index(0:Tmat11%nr))
    allocate(Tmat11%item(Tmat11%nnz), Tmat11%A(Tmat11%nnz))
    allocate(nz_cnt(Tmat11%nr), source=0)
    ! index
    do ilag=1,size(slaves4lag)
      nz_cnt(slaves4lag(ilag))=ndof*(hecLagMAT%indexL_lagrange(ilag)-hecLagMAT%indexL_lagrange(ilag-1))-1
    enddo
    Tmat11%index(0)=0
    do i=1,Tmat11%nr
      Tmat11%index(i)=Tmat11%index(i-1)+nz_cnt(i)
    enddo
    deallocate(nz_cnt)
    if (Tmat11%nnz /= Tmat11%index(Tmat11%nr)) then
      write(0,*) Tmat11%nnz, Tmat11%index(Tmat11%nr)
      Tmat11%nnz = Tmat11%index(Tmat11%nr)
      !stop 'ERROR: Tmat11%nnz wrong'
    endif
    ! item and A
    do ilag=1,size(slaves4lag)
      i=slaves4lag(ilag)
      l=Tmat11%index(i-1)+1
      js=hecLagMAT%indexL_lagrange(ilag-1)+1
      je=hecLagMAT%indexL_lagrange(ilag)
      factor=-BLs_inv(ilag)
      do j=js,je
        k=hecLagMAT%itemL_lagrange(j)
        do jdof=1,ndof
          kk=(k-1)*ndof+jdof
          jj=(j-1)*ndof+jdof
          if (kk==i) cycle
          Tmat11%item(l)=kk
          Tmat11%A(l)=hecLagMAT%AL_lagrange(jj)*factor
          l=l+1
        enddo
      enddo
      if (l /= Tmat11%index(i)+1) then
        write(0,*) l, Tmat11%index(i)+1
        stop 'ERROR: Tmat11%index wrong'
      endif
    enddo
    !call hecmw_localmat_write(Tmat11, 0)
    ! make 3x3-block version of Tmat
    call hecmw_localmat_blocking(Tmat11, ndof, Tmat)
    call hecmw_localmat_free(Tmat11)
  end subroutine add_C_to_Tmat

  !> \brief Make 3x3-blocked T^t matrix
  !>
  subroutine add_Ct_to_Ttmat(hecMAT, hecLagMAT, n, slaves4lag, BUs_inv, Ttmat)
    type(hecmwST_matrix),          intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT     !< original matrix for lagrange multipliers
    integer(kind=kint),            intent(in)    :: n             !< num of nodes in this subdomain
    integer(kind=kint),            intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),              intent(in)    :: BUs_inv(:)    !< inverse of diagonal BUs matrix
    type(hecmwST_local_matrix),    intent(out)   :: Ttmat         !< blocked T^t matrix
    !
    type(hecmwST_local_matrix) :: Ttmat11
    integer(kind=kint), allocatable :: nz_cnt(:)
    integer(kind=kint) :: ndof, i, idof, idx, ilag, l, js, je, j, k

    ndof=hecMAT%NDOF
    Ttmat11%nr=n*ndof
    Ttmat11%nc=Ttmat11%nr
    Ttmat11%nnz=hecLagMAT%numU_lagrange*ndof-hecLagMAT%num_lagrange
    Ttmat11%ndof=1

    allocate(Ttmat11%index(0:Ttmat11%nr))
    allocate(Ttmat11%item(Ttmat11%nnz), Ttmat11%A(Ttmat11%nnz))
    allocate(nz_cnt(Ttmat11%nr), source=0)
    ! index
    if (hecLagMAT%num_lagrange > 0) then
      do i=1,n
        do idof=1,ndof
          idx=(i-1)*ndof+idof
          nz_cnt(idx)=hecLagMAT%indexU_lagrange(i)-hecLagMAT%indexU_lagrange(i-1)
        enddo
      enddo
      do ilag=1,size(slaves4lag)
        nz_cnt(slaves4lag(ilag))=0
      enddo
    endif
    Ttmat11%index(0)=0
    do i=1,Ttmat11%nr
      Ttmat11%index(i)=Ttmat11%index(i-1)+nz_cnt(i)
    enddo
    if (Ttmat11%nnz /= Ttmat11%index(Ttmat11%nr)) then
      write(0,*) Ttmat11%nnz, Ttmat11%index(Ttmat11%nr)
      !stop 'ERROR: Ttmat11%nnz wrong'
      Ttmat11%nnz = Ttmat11%index(Ttmat11%nr)
    endif
    ! item and A
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        l=Ttmat11%index(idx-1)+1
        if (nz_cnt(idx) > 0) then
          ! offdiagonal
          js=hecLagMAT%indexU_lagrange(i-1)+1
          je=hecLagMAT%indexU_lagrange(i)
          do j=js,je
            k=hecLagMAT%itemU_lagrange(j)
            Ttmat11%item(l)=slaves4lag(k)
            Ttmat11%A(l)=-hecLagMAT%AU_lagrange((j-1)*ndof+idof)*BUs_inv(k)
            l=l+1
          enddo
        endif
        if (l /= Ttmat11%index(idx)+1) then
          write(0,*) l, Ttmat11%index(idx)+1
          stop 'ERROR: Ttmat11%index wrong'
        endif
      enddo
    enddo
    deallocate(nz_cnt)
    !call hecmw_localmat_write(Ttmat11, 0)
    ! make 3x3-block version of Ttmat
    call hecmw_localmat_blocking(Ttmat11, ndof, Ttmat)
    call hecmw_localmat_free(Ttmat11)
  end subroutine add_Ct_to_Ttmat

  !> \brief Make list of INTERNAL dofs that are in contact in WHOLE MODEL
  !>
  subroutine make_slave_list(hecMESHtmp, ndof, slaves4lag, slaves)
    type(hecmwST_local_mesh),        intent(in)  :: hecMESHtmp    !< mesh updated for the current contact state
    integer(kind=kint),              intent(in)  :: ndof          !< num of dof per node
    integer(kind=kint),              intent(in)  :: slaves4lag(:) !< LOCALLY detected slave dofs (incl. EXTERNAL)
    integer(kind=kint), allocatable, intent(out) :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint), allocatable :: mark_slave(:)
    integer(kind=kint) :: n_slave
    integer(kind=kint) :: ilag, i
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    allocate(mark_slave(hecMESHtmp%n_node*ndof), source=0)
    do ilag=1,size(slaves4lag)
      mark_slave(slaves4lag(ilag))=1
    enddo
    call hecmw_assemble_I(hecMESHtmp, mark_slave, hecMESHtmp%n_node, ndof)
    n_slave = 0
    do i = 1, hecMESHtmp%nn_internal * ndof
      if (mark_slave(i) /= 0) n_slave = n_slave + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: n_slave',n_slave
    allocate(slaves(n_slave))
    n_slave = 0
    do i = 1, hecMESHtmp%nn_internal * ndof
      if (mark_slave(i) /= 0) then
        n_slave = n_slave + 1
        slaves(n_slave) = i
      endif
    enddo
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: slaves',slaves(:)
    deallocate(mark_slave)
  end subroutine make_slave_list

  !> \brief Place one on diagonal of nonslave dof of T or T^t matrix
  !>
  subroutine add_Ip_to_Tmat(Tmat, slaves)
    type(hecmwST_local_matrix), intent(inout) :: Tmat      !< blocked T matrix
    integer(kind=kint),         intent(in)    :: slaves(:) !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    type(hecmwST_local_matrix) :: Imat, Wmat
    integer(kind=kint) :: ndof, ndof2, i, irow, idof

    ndof = Tmat%ndof
    ndof2 = ndof*ndof
    ! Imat: unit matrix except for slave dofs
    Imat%nr = Tmat%nr
    Imat%nc = Tmat%nc
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
    do irow = 1, Imat%nr
      do idof = 1, ndof
        Imat%A(ndof2*(irow-1)+ndof*(idof-1)+idof) = 1.0d0
      enddo
    enddo
    do i = 1, size(slaves)
      irow = (slaves(i)+ndof-1)/ndof
      idof = slaves(i)-ndof*(irow-1)
      Imat%A(ndof2*(irow-1)+ndof*(idof-1)+idof) = 0.0d0
    enddo
    call hecmw_localmat_add(Tmat, Imat, Wmat)
    call hecmw_localmat_free(Tmat)
    call hecmw_localmat_free(Imat)
    Tmat%nr = Wmat%nr
    Tmat%nc = Wmat%nc
    Tmat%nnz = Wmat%nnz
    Tmat%ndof = Wmat%ndof
    Tmat%index => Wmat%index
    Tmat%item => Wmat%item
    Tmat%A => Wmat%A
  end subroutine add_Ip_to_Tmat

  !> \brief Make comm table for contact dofs for optimized communication
  !>
  subroutine make_contact_comm_table(hecMESH, hecMAT, hecLagMAT, conCOMM)
    type(hecmwST_local_mesh),      intent(in)  :: hecMESH   !< original mesh
    type(hecmwST_matrix),          intent(in)  :: hecMAT    !< original matrix excl. contact
    type(hecmwST_matrix_lagrange), intent(in)  :: hecLagMAT !< original matrix for lagrange multipliers
    type(hecmwST_contact_comm),    intent(out) :: conCOMM   !< contact comm table for optimized communication
    !
    integer(kind=kint)              :: n_contact_dof
    integer(kind=kint), allocatable :: contact_dofs(:)

    ! make list of contact dofs including masters and slaves
    call make_contact_dof_list(hecMAT, hecLagMAT, n_contact_dof, contact_dofs)

    ! make comm_table for contact dofs
    call hecmw_contact_comm_init(conCOMM, hecMESH, hecMAT%ndof, n_contact_dof, contact_dofs)
  end subroutine make_contact_comm_table

  !> \brief Make list of contact dofs including all dofs of master and slave nodes
  !>
  subroutine make_contact_dof_list(hecMAT, hecLagMAT, n_contact_dof, contact_dofs)
    type(hecmwST_matrix),            intent(in)  :: hecMAT          !< original matrix excl. contact
    type(hecmwST_matrix_lagrange),   intent(in)  :: hecLagMAT       !< original matrix for lagrange multipliers
    integer(kind=kint),              intent(out) :: n_contact_dof   !< num of contact DOFs
    integer(kind=kint), allocatable, intent(out) :: contact_dofs(:) !< list of contact DOFs
    !
    integer(kind=kint) :: ndof, icnt, ilag, ls, le, l, jnode, k, inode, idof, i
    integer(kind=kint), allocatable :: iw(:)
    logical :: found
    integer(kind=kint) :: myrank

    if (hecLagMAT%num_lagrange == 0) then
      n_contact_dof = 0
      return
    endif
    ! lower
    ndof = hecMAT%NDOF
    allocate(iw(hecMAT%NP))
    icnt = 0
    do ilag = 1, hecLagMAT%num_lagrange
      ls = hecLagMAT%indexL_lagrange(ilag-1)+1
      le = hecLagMAT%indexL_lagrange(ilag)
      lloop1: do l = ls, le
        jnode = hecLagMAT%itemL_lagrange(l)
        do k = 1, icnt
          if (iw(k) == jnode) cycle lloop1
        enddo
        icnt = icnt + 1
        iw(icnt) = jnode
      enddo lloop1
    enddo
    ! upper
    do inode = 1, hecMAT%NP
      ls = hecLagMAT%indexU_lagrange(inode-1)+1
      le = hecLagMAT%indexU_lagrange(inode)
      if (ls <= le) then
        found = .false.
        do k = 1, icnt
          if (iw(k) == inode) found = .true.
        enddo
        if (.not. found) then
          icnt = icnt + 1
          iw(icnt) = inode
        endif
      endif
    enddo
    call quick_sort(iw, 1, icnt)
    allocate(contact_dofs(icnt*ndof))
    do i = 1, icnt
      do idof = 1, ndof
        contact_dofs((i-1)*ndof+idof) = (iw(i)-1)*ndof+idof
      enddo
    enddo
    n_contact_dof = icnt*ndof
    deallocate(iw)
    myrank = hecmw_comm_get_rank()
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: n_contact_dof',n_contact_dof
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: contact_dofs',contact_dofs(:)
  end subroutine make_contact_dof_list

  !> \brief Quick sort for integer array
  !>
  recursive subroutine quick_sort(array, id1, id2)
    integer(kind=kint), intent(inout) :: array(:)  !< integer array
    integer(kind=kint), intent(in)    :: id1, id2  !< index from id1 to id2 are sorted
    !
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

  !> \brief Assemble hecMAT and conMAT into Kmat and Btot
  !>
  subroutine assemble_equation(hecMESH, hecMESHtmp, hecMAT, conMAT, num_lagrange, &
      slaves, Kmat, Btot)
    type(hecmwST_local_mesh),   intent(in)    :: hecMESH      !< original mesh
    type(hecmwST_local_mesh),   intent(inout) :: hecMESHtmp   !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(in)    :: hecMAT       !< original matrix excl. contact
    type(hecmwST_matrix),       intent(in)    :: conMAT       !< original matrix for contact
    integer(kind=kint),         intent(in)    :: num_lagrange !< num of lagrange multipliers in this subdomain
    integer(kind=kint),         intent(in)    :: slaves(:)    !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    ! type(hecmwST_contact_comm), intent(in)    :: conCOMM      !< contact comm table for optimized communication
    type(hecmwST_local_matrix), intent(out)   :: Kmat         !< total K matrix assembled from hecMAT and conMAT
    real(kind=kreal),           intent(out)   :: Btot(:)      !< total RHS vector assembled from hecMAT%B and conMAT%B
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    call assemble_matrix(hecMESH, hecMESHtmp, hecMAT, conMAT, num_lagrange, Kmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: assemble matrix done'

    call assemble_rhs(hecMESH, hecMAT, conMAT, num_lagrange, slaves, Btot)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: assemble rhs done'

  end subroutine assemble_equation

  !> \brief Assemble hecMAT and conMAT into Kmat
  !>
  subroutine assemble_matrix(hecMESH, hecMESHtmp, hecMAT, conMAT, num_lagrange, Kmat)
    type(hecmwST_local_mesh),   intent(in)    :: hecMESH      !< original mesh
    type(hecmwST_local_mesh),   intent(inout) :: hecMESHtmp   !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(in)    :: hecMAT       !< original matrix excl. contact
    type(hecmwST_matrix),       intent(in)    :: conMAT       !< original matrix for contact
    integer(kind=kint),         intent(in)    :: num_lagrange !< num of lagrange multipliers in this subdomain
    type(hecmwST_local_matrix), intent(out)   :: Kmat         !< total K matrix assembled from hecMAT and conMAT
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    ! init Kmat and substitute conMAT
    call hecmw_localmat_init_with_hecmat(Kmat, conMAT, num_lagrange)
    if (DEBUG_MATRIX) call debug_write_matrix(Kmat, 'Kmat (conMAT local)')

    if (hecmw_comm_get_size() > 1) then
      ! communicate and assemble Kmat (updating hecMESHtmp)
      call hecmw_localmat_assemble(Kmat, hecMESH, hecMESHtmp)
      if (DEBUG_MATRIX) call debug_write_matrix(Kmat, 'Kmat (conMAT assembled)')
      if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: assemble K (conMAT) done'
    endif

    ! add hecMAT to Kmat
    call hecmw_localmat_add_hecmat(Kmat, hecMAT)
    if (DEBUG_MATRIX) call debug_write_matrix(Kmat, 'Kmat (hecMAT added)')
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: add hecMAT to K done'
  end subroutine assemble_matrix

  !> \brief Assemble hecMAT%B and conMAT%B into Btot
  !>
  subroutine assemble_rhs(hecMESH, hecMAT, conMAT, num_lagrange, slaves, Btot)
    type(hecmwST_local_mesh),   intent(in) :: hecMESH      !< original mesh
    type(hecmwST_matrix),       intent(in) :: hecMAT       !< original matrix excl. contact
    type(hecmwST_matrix),       intent(in) :: conMAT       !< original matrix for contact
    integer(kind=kint),         intent(in) :: num_lagrange !< num of lagrange multipliers in this subdomain
    integer(kind=kint),         intent(in) :: slaves(:)    !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    ! type(hecmwST_contact_comm), intent(in) :: conCOMM      !< contact comm table for optimized communication
    real(kind=kreal),           intent(out) :: Btot(:)     !< total RHS vector assembled from hecMAT%B and conMAT%B
    !
    integer(kind=kint) :: ndof, nndof, npndof, i, myrank

    myrank = hecmw_comm_get_rank()

    ndof = hecMAT%NDOF
    npndof = hecMAT%NP*ndof
    nndof  = hecMAT%N *ndof

    if (DEBUG_VECTOR) call debug_write_vector(hecMAT%B, 'RHS(hecMAT)', 'hecMAT%B', ndof, hecMAT%N, &
        hecMAT%NP, .false., num_lagrange, slaves)
    if (DEBUG_VECTOR) call debug_write_vector(conMAT%B, 'RHS(conMAT)', 'conMAT%B', ndof, conMAT%N, &
        conMAT%NP, .true., num_lagrange, slaves)

    do i=1,npndof+num_lagrange
      Btot(i) = conMAT%B(i)
    enddo

    if (hecmw_comm_get_size() > 1) then
      ! next line cannot be: call hecmw_contact_comm_reduce_r(conCOMM, Btot, HECMW_SUM) !!! alag_tied fails
      call hecmw_assemble_R(hecMESH, Btot, hecMAT%NP, ndof)
      if (DEBUG_VECTOR) call debug_write_vector(Btot, 'RHS(conMAT assembled)', 'Btot', ndof, conMAT%N, &
          conMAT%NP, .false., num_lagrange, slaves)
      if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: assemble RHS (conMAT%B) done'
    endif

    ! add hecMAT%B to Btot
    do i=1,nndof
      Btot(i)=Btot(i)+hecMAT%B(i)
    enddo
    if (DEBUG_VECTOR) call debug_write_vector(Btot, 'RHS(total)', 'Btot', ndof, conMAT%N, &
        conMAT%NP, .false., num_lagrange, slaves)
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: add hecMAT%B to RHS done'
  end subroutine assemble_rhs

  !> \brief Make converted equation
  !>
  subroutine convert_equation(hecMESHtmp, hecMAT, Kmat, Tmat, Ttmat, Btot, slaves, &
      slaves4lag, BLs_inv, conCOMM, hecTKT)
    type(hecmwST_local_mesh),   intent(inout) :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(in)    :: hecMAT        !< original matrix excl. contact
    type(hecmwST_local_matrix), intent(inout) :: Kmat          !< total K matrix assembled from hecMAT and conMAT
    type(hecmwST_local_matrix), intent(inout) :: Tmat          !< blocked T matrix
    type(hecmwST_local_matrix), intent(in)    :: Ttmat         !< blocked T^t matrix
    real(kind=kreal),           intent(in)    :: Btot(:)       !< total RHS vector assembled from hecMAT%B and conMAT%B
    integer(kind=kint),         intent(in)    :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    integer(kind=kint),         intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),           intent(in)    :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    type(hecmwST_contact_comm), intent(in)    :: conCOMM       !< contact comm table for optimized communication
    type(hecmwST_matrix),       intent(out)   :: hecTKT        !< converted system's matrix and RHS
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    call convert_matrix(hecMESHtmp, hecMAT, Ttmat, Kmat, Tmat, slaves, hecTKT)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: converted matrix'

    call convert_rhs(hecMESHtmp, hecMAT, hecTKT, Ttmat, Kmat, &
        slaves4lag, BLs_inv, Btot, conCOMM)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: converted RHS'
  end subroutine convert_equation

  !> \brief Make converted matrix
  !>
  subroutine convert_matrix(hecMESHtmp, hecMAT, Ttmat, Kmat, Tmat, slaves, hecTKT)
    type(hecmwST_local_mesh),   intent(inout) :: hecMESHtmp !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(in)    :: hecMAT     !< original matrix excl. contact
    type(hecmwST_local_matrix), intent(in)    :: Ttmat      !< blocked T^t matrix
    type(hecmwST_local_matrix), intent(inout) :: Kmat       !< total K matrix assembled from hecMAT and conMAT
    type(hecmwST_local_matrix), intent(inout) :: Tmat       !< blocked T matrix
    integer(kind=kint),         intent(in)    :: slaves(:)  !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    type(hecmwST_matrix),       intent(out)   :: hecTKT     !< converted system's matrix and RHS
    !
    type(hecmwST_local_matrix) :: TtKmat, TtKTmat
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    ! compute TtKmat = Ttmat * Kmat (updating hecMESHtmp)
    call hecmw_localmat_multmat(Ttmat, Kmat, hecMESHtmp, TtKmat)
    if (DEBUG_MATRIX) call debug_write_matrix(TtKmat, 'TtKmat')
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: multiply Tt and K done'

    ! compute TtKTmat = TtKmat * Tmat (updating hecMESHtmp)
    call hecmw_localmat_multmat(TtKmat, Tmat, hecMESHtmp, TtKTmat)
    if (DEBUG_MATRIX) call debug_write_matrix(TtKTmat, 'TtKTmat')
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: multiply TtK and T done'
    call hecmw_localmat_free(TtKmat)

    ! shrink comm_table
    ! call hecmw_localmat_shrink_comm_table(TtKTmat, hecMESHtmp)

    call place_one_on_diag_of_slave_dof(TtKTmat, slaves)
    if (DEBUG_MATRIX) call debug_write_matrix(TtKTmat, 'TtKTmat (place 1.0 on slave diag)')

    call hecmw_mat_init(hecTKT)
    call hecmw_localmat_make_hecmat(hecMAT, TtKTmat, hecTKT)
    if (DEBUG >= 3) write(0,*) '    DEBUG3[',myrank,']: convert TtKT to hecTKT done'
    call hecmw_localmat_free(TtKTmat)
  end subroutine convert_matrix

  !> \brief Place 1.0 on diagonal of slave dof of T^tKT matrix
  !>
  subroutine place_one_on_diag_of_slave_dof(TtKTmat, slaves)
    type(hecmwST_local_matrix), intent(inout) :: TtKTmat   !< converted matrix
    integer(kind=kint),         intent(in)    :: slaves(:) !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: ndof, ndof2, i, irow, idof, js, je, j, jcol

    ndof = TtKTmat%ndof
    ndof2 = ndof*ndof
    do i = 1, size(slaves)
      irow = (slaves(i)+ndof-1)/ndof
      idof = slaves(i)-ndof*(irow-1)
      js = TtKTmat%index(irow-1)+1
      je = TtKTmat%index(irow)
      do j = js, je
        jcol = TtKTmat%item(j)
        if (irow /= jcol) cycle
        if (abs(TtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof)) > tiny(0.0d0)) &
            stop 'ERROR: nonzero diag on slave dof of TtKTmat'
        TtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof) = 1.0d0
      enddo
    enddo
  end subroutine place_one_on_diag_of_slave_dof

  !> \brief Make converted RHS vector
  !>
  subroutine convert_rhs(hecMESHtmp, hecMAT, hecTKT, Ttmat, Kmat, &
       slaves4lag, BLs_inv, Btot, conCOMM)
    type(hecmwST_local_mesh),   intent(in)    :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(in)    :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix),       intent(inout) :: hecTKT        !< converted system's matrix and RHS
    type(hecmwST_local_matrix), intent(in)    :: Ttmat         !< blocked T^t matrix
    type(hecmwST_local_matrix), intent(in)    :: Kmat          !< total K matrix assembled from hecMAT and conMAT
    integer(kind=kint),         intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),           intent(in)    :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    real(kind=kreal), target,   intent(in)    :: Btot(:)       !< total RHS vector assembled from hecMAT%B and conMAT%B
    type(hecmwST_contact_comm), intent(in)    :: conCOMM       !< contact comm table for optimized communication
    !
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: Blag(:)
    integer(kind=kint) :: ndof, npndof, nndof, npndof_new, num_lagrange, i

    ! SIZE:
    ! Btot     <=> hecMAT, hecMESH
    ! Btmp     <=> Kmat, hecMESHtmp
    ! hecTKT%B <=> hecTKT, hecMESHtmp
    ndof = hecMAT%NDOF
    npndof     = hecMAT%NP*ndof
    nndof      = hecMAT%N *ndof
    npndof_new = hecTKT%NP*ndof
    num_lagrange = size(slaves4lag)

    allocate(hecTKT%B(npndof_new), source=0.d0)
    allocate(hecTKT%X(npndof_new), source=0.d0)
    allocate(Btmp(npndof_new))
    !
    !! Ttmat*(B+K*(-Bs^-1)*Blag)
    !
    ! B2=-Bs^-1*Blag
    Blag => Btot(npndof+1:npndof+num_lagrange)
    hecTKT%B(slaves4lag(:))=-BLs_inv(:)*Blag(:)
    ! send external contact dof => recv internal contact dof
    ! next line can be: call hecmw_assemble_R(hecMESHtmp, hecTKT%B, hecTKT%NP, ndof)
    call hecmw_contact_comm_reduce_r(conCOMM, hecTKT%B, HECMW_SUM)
    ! Btmp=B+K*B2 (including update of hecTKT%B)
    call hecmw_update_R(hecMESHtmp, hecTKT%B, hecTKT%NP, ndof)
    call hecmw_localmat_mulvec(Kmat, hecTKT%B, Btmp)
    do i=1,nndof
      Btmp(i)=Btot(i)+Btmp(i)
    enddo
    ! B2=Ttmat*Btmp
    call hecmw_update_R(hecMESHtmp, Btmp, hecTKT%NP, ndof)
    call hecmw_localmat_mulvec(Ttmat, Btmp, hecTKT%B)
    deallocate(Btmp)

    if (DEBUG_VECTOR) call debug_write_vector(hecTKT%B, 'RHS(converted)', 'hecTKT%B', ndof, hecTKT%N)
  end subroutine convert_rhs

  !> \brief Recover solution in original system
  !>
  subroutine recover_solution(hecMESHtmp, hecMAT, hecTKT, Tmat, Kmat, Btot, &
      slaves4lag, BLs_inv, BUs_inv, conCOMM, slaves)
    ! type(hecmwST_local_mesh),   intent(in)    :: hecMESH  !< original mesh
    type(hecmwST_local_mesh),   intent(in)    :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix),       intent(inout) :: hecTKT        !< converted system's matrix and RHS
    type(hecmwST_local_matrix), intent(in)    :: Tmat          !< blocked T matrix
    type(hecmwST_local_matrix), intent(in)    :: Kmat          !< total K matrix assembled from hecMAT and conMAT
    real(kind=kreal),           intent(in)    :: Btot(:)       !< total RHS vector assembled from hecMAT%B and conMAT%B
    integer(kind=kint),         intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),           intent(in)    :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    real(kind=kreal),           intent(in)    :: BUs_inv(:)    !< inverse of diagonal BUs matrix
    type(hecmwST_contact_comm), intent(in)    :: conCOMM       !< contact comm table for optimized communication
    integer(kind=kint),         intent(in)    :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()

    hecMAT%Iarray=hecTKT%Iarray
    hecMAT%Rarray=hecTKT%Rarray

    call comp_x_slave(hecMESHtmp, hecMAT, hecTKT, Tmat, Btot, &
        slaves4lag, BLs_inv, conCOMM, slaves)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: recovered slave disp'

    call comp_lag(hecMESHtmp, hecMAT, hecTKT, Kmat, Btot, &
        slaves4lag, BUs_inv, conCOMM, slaves)
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: recovered lag'

    if (DEBUG_VECTOR) call debug_write_vector(hecMAT%X, 'Solution(original)', 'hecMAT%X', hecMAT%NDOF, hecMAT%N, &
        hecMAT%NP, .false., size(slaves4lag), slaves)
  end subroutine recover_solution

  !> \brief Recover solution for slave dofs
  !>
  subroutine comp_x_slave(hecMESHtmp, hecMAT, hecTKT, Tmat, Btot, &
       slaves4lag, BLs_inv, conCOMM, slaves)
    ! type(hecmwST_local_mesh),   intent(in)    :: hecMESH  !< original mesh
    type(hecmwST_local_mesh),   intent(in)    :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix),       intent(in)    :: hecTKT        !< converted system's matrix and RHS
    type(hecmwST_local_matrix), intent(in)    :: Tmat          !< blocked T matrix
    real(kind=kreal), target,   intent(in)    :: Btot(:)       !< total RHS vector assembled from hecMAT%B and conMAT%B
    integer(kind=kint),         intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),           intent(in)    :: BLs_inv(:)    !< inverse of diagonal BLs matrix
    type(hecmwST_contact_comm), intent(in)    :: conCOMM       !< contact comm table for optimized communication
    integer(kind=kint),         intent(in)    :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: ndof, ndof2, npndof, nndof, num_lagrange
    real(kind=kreal), allocatable :: Xtmp(:)
    real(kind=kreal), pointer :: Blag(:)

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N  * ndof
    num_lagrange = size(slaves4lag)
    !!
    !! {X} = [T] {Xp} - [-Bs^-1] {c}
    !!
    ! compute {X} = [T] {Xp}
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, ndof)
    call hecmw_localmat_mulvec(Tmat, hecTKT%X, hecMAT%X)
    !
    ! compute {Xtmp} = [-Bs^-1] {c}
    allocate(Xtmp(npndof), source=0.0d0)
    Blag => Btot(npndof+1:npndof+num_lagrange)
    Xtmp(slaves4lag(:)) = -BLs_inv(:) * Blag(:)
    !
    ! send external contact dof => recv internal contact dof
    ! next line can be: call hecmw_assemble_R(hecMESH, Xtmp, hecMAT%NP, ndof)
    call hecmw_contact_comm_reduce_r(conCOMM, Xtmp, HECMW_SUM)
    !
    ! {X} = {X} - {Xtmp}
    hecMAT%X(slaves(:)) = hecMAT%X(slaves(:)) - Xtmp(slaves(:))
    deallocate(Xtmp)
  end subroutine comp_x_slave

  !> \brief Recover solution for lagrange multipliers
  !>
  subroutine comp_lag(hecMESHtmp, hecMAT, hecTKT, Kmat, Btot, &
       slaves4lag, BUs_inv, conCOMM, slaves)
    ! type(hecmwST_local_mesh),   intent(in)    :: hecMESH  !< original mesh
    type(hecmwST_local_mesh),   intent(in)    :: hecMESHtmp    !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),       intent(inout) :: hecMAT        !< original matrix excl. contact
    type(hecmwST_matrix),       intent(inout) :: hecTKT        !< converted system's matrix and RHS
    type(hecmwST_local_matrix), intent(in)    :: Kmat          !< total K matrix assembled from hecMAT and conMAT
    real(kind=kreal),           intent(in)    :: Btot(:)       !< total RHS vector assembled from hecMAT%B and conMAT%B
    integer(kind=kint),         intent(in)    :: slaves4lag(:) !< list of slave dofs chosed for EACH Lag. in THIS SUBDOMAIN
    real(kind=kreal),           intent(in)    :: BUs_inv(:)    !< inverse of diagonal BUs matrix
    type(hecmwST_contact_comm), intent(in)    :: conCOMM       !< contact comm table for optimized communication
    integer(kind=kint),         intent(in)    :: slaves(:)     !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: ndof, npndof, nndof, npndof_new, num_lagrange
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: xlag(:)

    ndof=hecMAT%ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N * ndof
    npndof_new = hecTKT%NP * ndof
    num_lagrange = size(slaves4lag)

    !! <SUMMARY>
    !! {lag} = [Bs^-T] ( {fs} - [Ksp Kss] {u} )
    !!
    ! 1. {Btmp} = [Kmat] {X}
    hecTKT%X(1:nndof) = hecMAT%X(1:nndof)
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, ndof)
    allocate(Btmp(npndof))
    call hecmw_localmat_mulvec(Kmat, hecTKT%X, Btmp)
    !
    ! 2. {Btmp_s} = {fs} - {Btmp_s}
    Btmp(slaves(:)) = Btot(slaves(:)) - Btmp(slaves(:))
    !
    ! 3. send internal contact dof => recv external contact dof
    ! next line can be: call hecmw_update_R(hecMESH, Btmp, hecMAT%NP, ndof)
    call hecmw_contact_comm_bcast_r(conCOMM, Btmp)
    !
    ! 4. {lag} = [Bs^-T] {Btmp_s}
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    xlag(:)=BUs_inv(:)*Btmp(slaves4lag(:))
    deallocate(Btmp)
  end subroutine comp_lag

  !> \brief Check solution in original system
  !>
  subroutine check_solution(hecMESH, hecMESHtmp, hecMAT, hecTKT, hecLagMAT, Kmat, Btot, &
       conCOMM, slaves)
    type(hecmwST_local_mesh),       intent(in)    :: hecMESH    !< original mesh
    type(hecmwST_local_mesh),       intent(in)    :: hecMESHtmp !< temoprary copy of mesh for migrating nodes
    type(hecmwST_matrix),           intent(inout) :: hecMAT     !< original matrix excl. contact
    type(hecmwST_matrix),           intent(inout) :: hecTKT     !< converted system's matrix and RHS
    type(hecmwST_matrix_lagrange) , intent(in)    :: hecLagMAT  !< original matrix for lagrange multipliers
    type(hecmwST_local_matrix),     intent(in)    :: Kmat       !< total K matrix assembled from hecMAT and conMAT
    real(kind=kreal), target,       intent(in)    :: Btot(:)    !< total RHS vector assembled from hecMAT%B and conMAT%B
    type(hecmwST_contact_comm),     intent(in)    :: conCOMM    !< contact comm table for optimized communication
    integer(kind=kint),             intent(in)    :: slaves(:)  !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
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
    allocate(r(npndof + num_lagrange), source=0.0d0)
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
    do i = 1, nndof
      hecTKT%X(i) = hecMAT%X(i)
    enddo
    call hecmw_update_R(hecMESHtmp, hecTKT%X, hecTKT%NP, ndof)
    call hecmw_localmat_mulvec(Kmat, hecTKT%X, Btmp)
    do i = 1, nndof
      r(i) = Btot(i) - Btmp(i)
    enddo
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
    ! next line can be: call hecmw_assemble_R(hecMESH, Btmp, hecMAT%NP, ndof)
    call hecmw_contact_comm_reduce_r(conCOMM, Btmp, HECMW_SUM)
    do i = 1, nndof
      r(i) = r(i) - Btmp(i)
    enddo
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
    if (DEBUG_VECTOR) call debug_write_vector(r, 'Residual', 'R', ndof, hecMAT%N, &
        hecMAT%NP, .false., hecLagMAT%num_lagrange, slaves)
    !
    call hecmw_InnerProduct_R(hecMESH, NDOF, r, r, rnrm2)
    call hecmw_InnerProduct_R(hecMESH, NDOF, Btot, Btot, bnrm2)
    rlagnrm2 = dot_product(rlag, rlag)
    call hecmw_allreduce_R1(hecMESH, rlagnrm2, HECMW_SUM)
    blagnrm2 = dot_product(blag, blag)
    call hecmw_allreduce_R1(hecMESH, blagnrm2, HECMW_SUM)
    !
    if (myrank == 0) then
      write(0,*) 'INFO: resid(x,lag,tot)',sqrt(rnrm2),sqrt(rlagnrm2),sqrt(rnrm2+rlagnrm2)
      write(0,*) 'INFO: rhs  (x,lag,tot)',sqrt(bnrm2),sqrt(blagnrm2),sqrt(bnrm2+blagnrm2)
    endif
  end subroutine check_solution

  !> \brief Check solution in original system (another version)
  !>
  subroutine check_solution2(hecMESH, hecMAT, conMAT, hecLagMAT, conCOMM, slaves)
    implicit none
    type(hecmwST_local_mesh),       intent(in)    :: hecMESH   !< original mesh
    type(hecmwST_matrix),           intent(inout) :: hecMAT    !< original matrix excl. contact
    type(hecmwST_matrix),           intent(in)    :: conMAT    !< original matrix for contact
    type(hecmwST_matrix_lagrange) , intent(in)    :: hecLagMAT !< original matrix for lagrange multipliers
    type(hecmwST_contact_comm),     intent(in)    :: conCOMM   !< contact comm table for optimized communication
    integer(kind=kint),             intent(in)    :: slaves(:) !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: ndof, ndof2, nndof, npndof, num_lagrange
    integer(kind=kint) :: i, idof, j, jdof, ls, le, l
    integer(kind=kint) :: irow, js, je, jcol
    real(kind=kreal), allocatable, target :: r(:)
    real(kind=kreal), allocatable :: r_con(:)
    real(kind=kreal), pointer :: rlag(:), blag(:), xlag(:)
    real(kind=kreal) :: rnrm2, rlagnrm2
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    num_lagrange = hecLagMAT%num_lagrange
    !
    allocate(r(npndof + num_lagrange))
    r(:) = 0.0d0
    allocate(r_con(npndof))
    r_con(:) = 0.0d0
    !
    rlag => r(npndof+1:npndof+num_lagrange)
    blag => conMAT%B(npndof+1:npndof+num_lagrange)
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    !
    !! <SUMMARY>
    !! {r}    = {b} - [K] {x} - [Bt] {lag}
    !! {rlag} = {c} - [B] {x}
    !
    ! 1. {r} = {r_org} + {r_con}
    ! {r_org} = {b_org} - [hecMAT] {x}
    ! {r_con} = {b_con} - [conMAT] {x} - [Bt] {lag}
    !
    ! 1.1 {r_org}
    ! {r} = {b} - [K] {x}
    call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r)
    !
    if (DEBUG_VECTOR) call debug_write_vector(r, 'Residual(original)', 'R', ndof, hecMAT%N, &
        hecMAT%NP, .false., num_lagrange, slaves)
    !
    ! 1.2 {r_con}
    ! 1.2.1 {r_con} = {b_con} - [conMAT]{x}
    ! Note: compute not only internal DOFs but also external DOFs, locally
    ! !call hecmw_update_3_R(hecMESH, hecMAT%X, hecMAT%NP)  ! X is already updated
    do i = 1, npndof
      r_con(i) = conMAT%B(i)
    enddo
    do irow = 1,hecMAT%NP
      ! lower
      js = conMAT%indexL(irow-1)+1
      je = conMAT%indexL(irow)
      do j = js, je
        jcol = conMAT%itemL(j)
        do idof = 1, ndof
          i = ndof*(irow-1)+idof
          do jdof = 1, ndof
            r_con(i) = r_con(i) - conMAT%AL(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
          enddo
        enddo
      enddo
      ! diag
      do idof = 1, ndof
        i = ndof*(irow-1)+idof
        do jdof = 1, ndof
          r_con(i) = r_con(i) - conMAT%D(ndof2*(irow-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(irow-1)+jdof)
        enddo
      enddo
      ! upper
      js = conMAT%indexU(irow-1)+1
      je = conMAT%indexU(irow)
      do j = js, je
        jcol = conMAT%itemU(j)
        do idof = 1, ndof
          i = ndof*(irow-1)+idof
          do jdof = 1, ndof
            r_con(i) = r_con(i) - conMAT%AU(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
          enddo
        enddo
      enddo
    end do
    !
    ! 1.2.2 {r_con} = {r_con} - [Bt] {lag}
    if (num_lagrange > 0) then
      do i = 1, hecMAT%NP
        ls = hecLagMAT%indexU_lagrange(i-1)+1
        le = hecLagMAT%indexU_lagrange(i)
        do l = ls, le
          j = hecLagMAT%itemU_lagrange(l)
          do idof = 1, ndof
            r_con(ndof*(i-1)+idof) = r_con(ndof*(i-1)+idof) - hecLagMAT%AU_lagrange(ndof*(l-1)+idof) * xlag(j)
          enddo
        enddo
      enddo
    endif
    !
    if (DEBUG_VECTOR) call debug_write_vector(r, 'Residual(contact,local)', 'R_con', ndof, hecMAT%N, &
        hecMAT%NP, .true., num_lagrange, slaves)
    !
    ! 1.2.3 send external part of {r_con}
    call hecmw_contact_comm_reduce_r(conCOMM, r_con, HECMW_SUM)
    !
    if (DEBUG_VECTOR) call debug_write_vector(r, 'Residual(contact,assembled)', 'R_con', ndof, hecMAT%N, &
        hecMAT%NP, .false., num_lagrange, slaves)
    !
    ! 1.3 {r} = {r_org} + {r_con}
    do i = 1,nndof
      r(i) = r(i) + r_con(i)
    enddo
    !
    if (DEBUG_VECTOR) call debug_write_vector(r, 'Residual(total)', 'R', ndof, hecMAT%N, &
        hecMAT%NP, .false., num_lagrange, slaves)
    !
    ! 2.  {rlag} = {c} - [B] {x}
    !call hecmw_update_3_R(hecMESH, hecMAT%X, hecMAT%NP)  ! X is already updated
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
    if (DEBUG_VECTOR) then
      write(1000+myrank,*) 'Residual(lagrange)-----------------------------------------------------'
      if (num_lagrange > 0) then
        write(1000+myrank,*) 'R(lag):',npndof+1,'-',npndof+num_lagrange
        write(1000+myrank,*) r(npndof+1:npndof+num_lagrange)
      endif
    endif
    !
    call hecmw_InnerProduct_R(hecMESH, NDOF, r, r, rnrm2)
    rlagnrm2 = dot_product(rlag, rlag)
    call hecmw_allreduce_R1(hecMESH, rlagnrm2, HECMW_SUM)
    !
    if (myrank == 0) write(0,*) 'INFO: resid(x,lag,tot)',sqrt(rnrm2),sqrt(rlagnrm2),sqrt(rnrm2+rlagnrm2)
  end subroutine check_solution2

  !> \brief Debug write matrix
  !>
  subroutine debug_write_matrix(Mat, label)
    type(hecmwST_local_matrix), intent(in) :: Mat   !< matrix
    character(len=*),           intent(in) :: label !< label for matrix
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    write(1000+myrank,*) trim(label),'============================================================'
    call hecmw_localmat_write(Mat, 1000+myrank)
  end subroutine debug_write_matrix

  !> \brief Debug write vector
  !>
  subroutine debug_write_vector(Vec, label, name, ndof, N, &
      NP, write_ext, num_lagrange, slaves)
    real(kind=kreal),   intent(in) :: Vec(:) !< vector
    character(len=*),   intent(in) :: label  !< label for vector
    character(len=*),   intent(in) :: name   !< name of vector
    integer(kind=kint), intent(in) :: ndof   !< num of DOF per node
    integer(kind=kint), intent(in) :: N      !< num of nodes excl. external nodes
    integer(kind=kint), intent(in), optional :: NP           !< num of nodes incl. external nodes
    logical,            intent(in), optional :: write_ext    !< whether to write external dofs or not
    integer(kind=kint), intent(in), optional :: num_lagrange !< num of lagrange multipliers in this subdomain
    integer(kind=kint), intent(in), optional :: slaves(:)    !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: myrank

    myrank = hecmw_comm_get_rank()
    write(1000+myrank,*) trim(label),'------------------------------------------------------------'
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
    if (present(slaves)) then
      if (size(slaves) > 0) then
        write(1000+myrank,*) trim(name),'(slave):',slaves(:)
        write(1000+myrank,*) Vec(slaves(:))
      endif
    endif
  end subroutine debug_write_vector

end module m_solve_LINEQ_contact_elim
