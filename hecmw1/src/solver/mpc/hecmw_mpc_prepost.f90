!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_mpc_prepost
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_matrix_ass
  use hecmw_local_matrix
  use hecmw_solver_las
  implicit none

  private
  public :: hecmw_mpc_mat_init
  public :: hecmw_mpc_mat_init_explicit
  public :: hecmw_mpc_mat_finalize
  public :: hecmw_mpc_mat_finalize_explicit
  public :: hecmw_mpc_mat_ass
  public :: hecmw_mpc_trans_rhs
  public :: hecmw_mpc_tback_sol
  public :: hecmw_mpc_trans_mass
  public :: hecmw_mpc_tback_eigvec
  public :: hecmw_mpc_mark_slave

contains

  !C
  !C***
  !C*** hecmw_mpc_mat_init
  !C***
  !C
  subroutine hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(inout), target :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD, SOLVER_TYPE

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      hecMESHmpc => hecMESH
      hecMATmpc => hecMAT
      return
    endif

    call hecmw_mpc_scale(hecMESH)

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
    if (MPC_METHOD < 1 .or. 3 < MPC_METHOD) then
      SOLVER_TYPE = hecmw_mat_get_solver_type(hecMAT)
      if (SOLVER_TYPE > 1) then  ! DIRECT SOLVER
        MPC_METHOD = 1           !   default: penalty
      else                       ! ITERATIVE SOLVER
        MPC_METHOD = 3           !   default: elimination
      endif
      call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD)
    endif

    if (MPC_METHOD == 2) then
      write(*,*) 'WARNING: MPCMETHOD=2 (MPCCG) is deprecated; may not work correctly'
      ! MPC_METHOD = 3
      ! call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD)
    endif

    select case (MPC_METHOD)
    case (1)  ! penalty
      hecMESHmpc => hecMESH
      hecMATmpc => hecMAT
    case (2)  ! MPCCG
      hecMESHmpc => hecMESH
      hecMATmpc => hecMAT
    case (3)  ! elimination
      allocate(hecMESHmpc)
      call hecmw_mpc_mesh_copy(hecMESH, hecMESHmpc)
      allocate(hecMATmpc)
      call hecmw_mat_init(hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_init

  !C
  !C***
  !C*** hecmw_mpc_mat_init_explicit
  !C***
  !C
  subroutine hecmw_mpc_mat_init_explicit(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(inout), target :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      hecMATmpc => hecMAT
      return
    endif

    call hecmw_mpc_scale(hecMESH)

    ! Force MPC_METHOD=3
    MPC_METHOD = 3
    call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD)

    allocate(hecMATmpc)
    call hecmw_mat_init(hecMATmpc)

    hecMATmpc%N = hecMAT%N
    hecMATmpc%NP = hecMAT%NP
    hecMATmpc%NDOF = hecMAT%NDOF
    allocate(hecMATmpc%B(size(hecMAT%B)))
    allocate(hecMATmpc%X(size(hecMAT%X)))
  end subroutine hecmw_mpc_mat_init_explicit

  !C
  !C***
  !C*** hecmw_mpc_mat_finalize
  !C***
  !C
  subroutine hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
      return
    endif

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
    case (2) ! MPCCG
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
    case (3) ! elimination
      call hecmw_mpc_mesh_free(hecMESHmpc)
      deallocate(hecMESHmpc)
      nullify(hecMESHmpc)
      call hecmw_mat_finalize(hecMATmpc)
      deallocate(hecMATmpc)
      nullify(hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_finalize

  !C
  !C***
  !C*** hecmw_mpc_mat_finalize_explicit
  !C***
  !C
  subroutine hecmw_mpc_mat_finalize_explicit(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      nullify(hecMATmpc)
      return
    endif

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      nullify(hecMATmpc)
    case (2) ! MPCCG
      nullify(hecMATmpc)
    case (3) ! elimination
      call hecmw_mat_finalize(hecMATmpc)
      deallocate(hecMATmpc)
      nullify(hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_finalize_explicit

  !C
  !C***
  !C*** hecmw_mpc_mat_ass
  !C***
  !C
  subroutine hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Penalty"
      call hecmw_mat_ass_equation ( hecMESH, hecMAT )
    case (2)  ! MPCCG
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: MPC-CG"
    case (3)  ! elimination
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Elimination"
      call hecmw_trimatmul_TtKT_mpc(hecMESHmpc, hecMAT, hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_ass


  !C
  !C***
  !C*** hecmw_mpc_trans_rhs
  !C***
  !C
  subroutine hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal) :: time_dumm
    integer(kind=kint) :: totalmpc, MPC_METHOD, i

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      call hecmw_mat_ass_equation_rhs ( hecMESH, hecMATmpc )
    case (2) ! MPCCG
      allocate(Btmp(hecMAT%NP*hecMAT%NDOF))
      do i = 1, hecMAT%NP*hecMAT%NDOF
        Btmp(i) = hecMAT%B(i)
      enddo
      call hecmw_trans_b(hecMESH, hecMAT, Btmp, hecMATmpc%B, time_dumm)
      deallocate(Btmp)
    case (3) ! elimination
      call hecmw_trans_b(hecMESH, hecMAT, hecMAT%B, hecMATmpc%B, time_dumm)
      hecMATmpc%Iarray=hecMAT%Iarray
      hecMATmpc%Rarray=hecMAT%Rarray
    end select

  end subroutine hecmw_mpc_trans_rhs

  !C
  !C***
  !C*** hecmw_mpc_tback_sol
  !C***
  !C
  subroutine hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    real(kind=kreal) :: time_dumm
    integer(kind=kint) :: totalmpc, MPC_METHOD, i

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      ! do nothing
    case (2)  ! MPCCG
      call hecmw_tback_x(hecMESH, hecMAT%NDOF, hecMAT%X, time_dumm)
    case (3)  ! elimination
      do i = 1, size(hecMAT%X)
        hecMAT%X(i) = hecMATmpc%X(i)
      enddo
      call hecmw_tback_x(hecMESH, hecMAT%NDOF, hecMAT%X, time_dumm)
      hecMAT%Iarray=hecMATmpc%Iarray
      hecMAT%Rarray=hecMATmpc%Rarray
    end select
  end subroutine hecmw_mpc_tback_sol

  !C
  !C***
  !C*** hecmw_mpc_trans_mass
  !C***
  !C
  subroutine hecmw_mpc_trans_mass(hecMESH, hecMAT, mass)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: mass(:)

    real(kind=kreal), allocatable :: Mtmp(:)
    real(kind=kreal) :: time_dumm
    integer(kind=kint) :: totalmpc, MPC_METHOD, i

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      ! do nothing
    case (2,3) ! MPCCG or elimination
      allocate(Mtmp(hecMAT%NP*hecMAT%NDOF))
      !C-- {Mt} = [T'] {w}
      call hecmw_Ttvec(hecMESH, hecMAT%NDOF, mass, Mtmp, time_dumm)
      do i = 1, hecMAT%NP*hecMAT%NDOF
        mass(i) = Mtmp(i)
      enddo
      deallocate(Mtmp)
    end select

  end subroutine hecmw_mpc_trans_mass

  !C
  !C***
  !C*** hecmw_mpc_tback_eigvec
  !C***
  !C
  subroutine hecmw_mpc_tback_eigvec(hecMESH, hecMAT, neig, eigvec)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: neig
    real(kind=kreal), intent(inout) :: eigvec(:,:)

    real(kind=kreal) :: time_dumm
    integer(kind=kint) :: totalmpc, MPC_METHOD, i

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      ! do nothing
    case (2,3)  ! MPCCG or elimination
      do i = 1, neig
        call hecmw_tback_x(hecMESH, hecMAT%NDOF, eigvec(:,i), time_dumm)
        !!! need normalization???
      enddo
    end select
  end subroutine hecmw_mpc_tback_eigvec

  !C
  !C***
  !C*** hecmw_mpc_mark_slave
  !C***
  !C
  subroutine hecmw_mpc_mark_slave(hecMESH, hecMAT, mark)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: mark(:)

    integer(kind=kint) :: ndof, i, j, k, kk

    ndof = hecMAT%NDOF
    mark(:) = 0
    OUTER: do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1)+1
      kk = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      mark(kk) = 1
    enddo OUTER
  end subroutine hecmw_mpc_mark_slave

  !C
  !C***
  !C*** hecmw_mpc_scale
  !C***
  !C
  subroutine hecmw_mpc_scale(hecMESH)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint) :: i, j, k
    real(kind=kreal) :: WVAL

    !$omp parallel default(none),private(i,j,k,WVAL),shared(hecMESH)
    !$omp do
    do i = 1, hecMESH%mpc%n_mpc
      k = hecMESH%mpc%mpc_index(i-1)+1
      WVAL = 1.d0 / hecMESH%mpc%mpc_val(k)
      hecMESH%mpc%mpc_val(k) = 1.d0
      do j = hecMESH%mpc%mpc_index(i-1)+2, hecMESH%mpc%mpc_index(i)
        hecMESH%mpc%mpc_val(j) = hecMESH%mpc%mpc_val(j) * WVAL
      enddo
      hecMESH%mpc%mpc_const(i) = hecMESH%mpc%mpc_const(i) * WVAL
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_mpc_scale


  !C
  !C***
  !C*** hecmw_trans_b
  !C***
  !C
  subroutine hecmw_trans_b(hecMESH, hecMAT, B, BT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), intent(in) :: B(:)
    real(kind=kreal), intent(out), target :: BT(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal), allocatable :: W(:)
    real(kind=kreal), pointer :: XG(:)
    integer(kind=kint) :: ndof, i, j, k, kk, flg_bak

    ndof = hecMAT%NDOF

    allocate(W(hecMESH%n_node * ndof))

    !C===
    !C +---------------------------+
    !C | {bt}= [T']({b} - [A]{xg}) |
    !C +---------------------------+
    !C===
    XG => BT
    do i = 1, hecMAT%N * ndof
      XG(i) = 0.d0
    enddo

    !C-- Generate {xg} from mpc_const
    !$omp parallel default(none),private(i,k,kk),shared(hecMESH,XG),firstprivate(ndof)
    !$omp do
    OUTER: do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      XG(kk) = hecMESH%mpc%mpc_const(i)
    enddo OUTER
    !$omp end do
    !$omp end parallel

    !C-- {w} = {b} - [A]{xg}
    flg_bak = hecmw_mat_get_flag_mpcmatvec(hecMAT)
    call hecmw_mat_set_flag_mpcmatvec(hecMAT, 0)
    call hecmw_matresid(hecMESH, hecMAT, XG, B, W, COMMtime)
    call hecmw_mat_set_flag_mpcmatvec(hecMAT, flg_bak)

    !C-- {bt} = [T'] {w}
    call hecmw_Ttvec(hecMESH, ndof, W, BT, COMMtime)

    deallocate(W)
  end subroutine hecmw_trans_b


  !C
  !C***
  !C*** hecmw_tback_x
  !C***
  !C
  subroutine hecmw_tback_x(hecMESH, ndof, X, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(inout) :: X(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal), allocatable :: W(:)
    integer(kind=kint) :: i, j, k, kk

    allocate(W(hecMESH%n_node * ndof))

    !C-- {tx} = [T]{x}
    call hecmw_Tvec(hecMESH, ndof, X, W, COMMtime)

    !C-- {x} = {tx} + {xg}
    !$omp parallel default(none),private(i,k,kk),shared(hecMESH,X,W),firstprivate(ndof)
    !$omp do
    do i= 1, hecMESH%nn_internal * ndof
      X(i)= W(i)
    enddo
    !$omp end do

    !$omp do
    OUTER: do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      X(kk) = X(kk) + hecMESH%mpc%mpc_const(i)
    enddo OUTER
    !$omp end do
    !$omp end parallel

    deallocate(W)

    if (ndof == 3) then
      call hecmw_update_3_R(hecMESH, X, hecMESH%n_node)
    else if (ndof == 2) then
      call hecmw_update_2_R(hecMESH, X, hecMESH%n_node)
    else
      call hecmw_update_m_R(hecMESH, X, hecMESH%n_node, ndof)
    endif
  end subroutine hecmw_tback_x

  subroutine hecmw_mpc_mesh_copy(src, dst)
    implicit none
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
      dst%import_index(:)= src%import_index(:)
      allocate(dst%export_index(0:dst%n_neighbor_pe))
      dst%export_index(:)= src%export_index(:)
      allocate(dst%import_item(dst%import_index(dst%n_neighbor_pe)))
      dst%import_item(:) = src%import_item(:)
      allocate(dst%export_item(dst%export_index(dst%n_neighbor_pe)))
      dst%export_item(:) = src%export_item(:)
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
    dst%node      => src%node
  end subroutine hecmw_mpc_mesh_copy

  subroutine hecmw_mpc_mesh_free(hecMESH)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    if (hecMESH%n_neighbor_pe > 1) then
      deallocate(hecMESH%neighbor_pe)
      deallocate(hecMESH%import_index)
      deallocate(hecMESH%export_index)
      deallocate(hecMESH%import_item)
      deallocate(hecMESH%export_item)
    endif
    deallocate(hecMESH%global_node_ID)
    deallocate(hecMESH%node_ID)
    deallocate(hecMESH%elem_type_item)
  end subroutine hecmw_mpc_mesh_free
end module hecmw_mpc_prepost
