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
  public :: hecmw_mpc_check_solution
  public :: hecmw_mpc_check_equation

  integer, parameter :: DEBUG = 0
  logical, parameter :: DEBUG_VECTOR = .false.

contains

  !C
  !C***
  !C*** hecmw_mpc_mat_init
  !C***
  !C
  subroutine hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, conMAT, conMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(inout), target :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    type (hecmwST_matrix), intent(in), target, optional :: conMAT
    type (hecmwST_matrix), pointer, optional :: conMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD, SOLVER_TYPE

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      hecMESHmpc => hecMESH
      hecMATmpc => hecMAT
      if (present(conMAT).and.present(conMATmpc)) conMATmpc => conMAT
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
      if (present(conMAT).and.present(conMATmpc)) conMATmpc => conMAT
    case (2)  ! MPCCG
      hecMESHmpc => hecMESH
      hecMATmpc => hecMAT
      if (present(conMAT).and.present(conMATmpc)) conMATmpc => conMAT
    case (3)  ! elimination
      allocate(hecMESHmpc)
      call hecmw_mpc_mesh_copy(hecMESH, hecMESHmpc)
      allocate(hecMATmpc)
      call hecmw_mat_init(hecMATmpc)
      if (present(conMAT).and.present(conMATmpc)) then
        allocate(conMATmpc)
        call hecMW_mat_init(conMATmpc)
      endif
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
  subroutine hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, conMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    type (hecmwST_matrix), pointer, optional :: conMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
      if (present(conMATmpc)) nullify(conMATmpc)
      return
    endif

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
      if (present(conMATmpc)) nullify(conMATmpc)
    case (2) ! MPCCG
      nullify(hecMESHmpc)
      nullify(hecMATmpc)
      if (present(conMATmpc)) nullify(conMATmpc)
    case (3) ! elimination
      call hecmw_mpc_mesh_free(hecMESHmpc)
      deallocate(hecMESHmpc)
      nullify(hecMESHmpc)
      call hecmw_mat_finalize(hecMATmpc)
      deallocate(hecMATmpc)
      nullify(hecMATmpc)
      if (present(conMATmpc)) then
        call hecmw_mat_finalize(conMATmpc)
        deallocate(conMATmpc)
        nullify(conMATmpc)
      endif
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
  subroutine hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, conMAT, conMATmpc, hecLagMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    type (hecmwST_matrix), intent(in), optional :: conMAT
    type (hecmwST_matrix), pointer, optional :: conMATmpc
    type (hecmwST_matrix_lagrange), intent(inout), optional :: hecLagMAT
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
      if (present(conMAT).and.present(conMATmpc).and.present(hecLagMAT)) then
        call hecmw_trimatmul_TtKT_mpc(hecMESHmpc, conMAT, conMATmpc)
        call resize_hecLagMAT(conMAT%NP, conMATmpc%NP, conMAT%NDOF, hecLagMAT)
      endif
    end select

  end subroutine hecmw_mpc_mat_ass


  subroutine resize_hecLagMAT(NP_orig, NP_new, ndof, hecLagMAT)
    integer(kind=kint), intent(in) :: NP_orig, NP_new, ndof
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT
    integer(kind=kint), pointer :: itemp(:)

    if (hecLagMAT%num_lagrange == 0) return

    allocate(itemp(0:NP_new))
    itemp(0:NP_orig) = hecLagMAT%indexU_lagrange(0:NP_orig)
    itemp(NP_orig+1:NP_new) = hecLagMAT%indexU_lagrange(NP_orig)

    deallocate(hecLagMAT%indexU_lagrange)
    hecLagMAT%indexU_lagrange => itemp

  end subroutine resize_hecLagMAT

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
    integer(kind=kint) :: npndof, npndof_mpc, num_lagrange

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
      npndof = hecMAT%NP * hecMAT%NDOF
      do i = 1, npndof
        hecMAT%X(i) = hecMATmpc%X(i)
      enddo
      call hecmw_tback_x(hecMESH, hecMAT%NDOF, hecMAT%X, time_dumm)
      num_lagrange = size(hecMAT%X) - npndof
      npndof_mpc = hecMATmpc%NP * hecMATmpc%NDOF
      do i = 1, num_lagrange
        hecMAT%X(npndof+i) = hecMATmpc%X(npndof_mpc+i)
      enddo
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

    call debug_write_vector(B, 'original RHS', 'B', ndof, hecMAT%N, hecMAT%NP, .true.)

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

    call debug_write_vector(BT, 'transformed RHS', 'BT', ndof, hecMAT%N, hecMAT%NP, .true.)
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

    call debug_write_vector(X, 'solution for transformed eqn', 'X', ndof, hecMESH%nn_internal, &
        hecMESH%n_node, .true.)

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

    call hecmw_update_R(hecMESH, X, hecMESH%n_node, ndof)
    call debug_write_vector(X, 'recovered solution', 'X', ndof, hecMESH%nn_internal, &
        hecMESH%n_node, .true.)
  end subroutine hecmw_tback_x

  subroutine hecmw_mpc_check_solution(hecMESH, hecMAT)
    use hecmw_solver_misc
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT

    real(kind=kreal), allocatable :: kx(:), ttkx(:), ttb(:), resid(:)
    integer(kind=kint) :: ndof, nndof, npndof, i
    real(kind=kreal) :: bnrm, rnrm, commtime

    ndof = hecMAT%NDOF
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    allocate(kx(npndof), ttkx(npndof), ttb(npndof), resid(npndof))
    call hecmw_matvec(hecMESH, hecMAT, hecMAT%X, kx)
    call hecmw_Ttvec(hecMESH, ndof, kx, ttkx, commtime)
    call hecmw_Ttvec(hecMESH, ndof, hecMAT%B, ttb, commtime)
    do i = 1, nndof
      resid(i) = ttb(i) - ttkx(i)
    enddo
    call hecmw_innerProduct_R(hecMESH, ndof, ttb, ttb, bnrm)
    call hecmw_innerProduct_R(hecMESH, ndof, resid, resid, rnrm)
    bnrm = sqrt(bnrm)
    rnrm = sqrt(rnrm)
    write(0,*) 'DEBUG: hecmw_mpc_check_solution: resid(abs), resid(rel):',rnrm, rnrm/bnrm
  end subroutine hecmw_mpc_check_solution

  subroutine hecmw_mpc_check_equation(hecMESH, hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT

    integer(kind=kint) :: ndof
    integer(kind=kint) :: i, j, kk
    real(kind=kreal) :: xmax, sum, resid, resid_max

    ndof = hecMAT%NDOF
    xmax = 0.d0
    do i = 1, hecMESH%nn_internal * ndof
      if (xmax < abs(hecMAT%X(i))) xmax = abs(hecMAT%X(i))
    enddo
    call hecmw_allreduce_R1 (hecMESH, xmax, hecmw_max)

    resid_max = 0.d0
    OUTER: do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      sum = 0.d0
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        kk = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        sum = sum + hecMAT%X(kk) * hecMESH%mpc%mpc_val(j)
      enddo
      resid = abs(hecMESH%mpc%mpc_const(i) - sum)
      if (resid_max < resid) resid_max = resid
    enddo OUTER
    call hecmw_allreduce_R1 (hecMESH, resid_max, hecmw_max)
    write(0,*) 'DEBUG: hecmw_mpc_check_equation: resid(abs), resid(rel):', resid_max, resid_max/xmax
  end subroutine hecmw_mpc_check_equation

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

  !> \brief Debug write vector
  !>
  subroutine debug_write_vector(Vec, label, name, ndof, N, &
      NP, write_ext, slaves)
    real(kind=kreal),   intent(in) :: Vec(:) !< vector
    character(len=*),   intent(in) :: label  !< label for vector
    character(len=*),   intent(in) :: name   !< name of vector
    integer(kind=kint), intent(in) :: ndof   !< num of DOF per node
    integer(kind=kint), intent(in) :: N      !< num of nodes excl. external nodes
    integer(kind=kint), intent(in), optional :: NP           !< num of nodes incl. external nodes
    logical,            intent(in), optional :: write_ext    !< whether to write external dofs or not
    integer(kind=kint), intent(in), optional :: slaves(:)    !< list of INTERNAL dofs that are in contact in WHOLE MODEL
    !
    integer(kind=kint) :: iunit
    character(len=128) :: fmt

    if (.not. DEBUG_VECTOR) return

    write(fmt,'(a,i0,a)') '(',ndof,'f12.3)'

    iunit = 1000 + hecmw_comm_get_rank()
    write(iunit,*) trim(label),'------------------------------------------------------------'
    write(iunit,*) 'size of ',trim(name),size(Vec)
    write(iunit,*) trim(name),': 1-',N*ndof
    write(iunit,fmt) Vec(1:N*ndof)
    if (present(write_ext) .and. present(NP)) then
      if (write_ext) then
        write(iunit,*) trim(name),'(external): ',N*ndof+1,'-',NP*ndof
        write(iunit,fmt) Vec(N*ndof+1:NP*ndof)
      endif
    endif
    if (present(slaves)) then
      if (size(slaves) > 0) then
        write(iunit,*) trim(name),'(slave):',slaves(:)
        write(iunit,fmt) Vec(slaves(:))
      endif
    endif
  end subroutine debug_write_vector
end module hecmw_mpc_prepost
