!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2016/10/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_mpc_prepost
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_matrix_ass
  use hecmw_local_matrix
  use hecmw_solver_las_11
  use hecmw_solver_las_22
  use hecmw_solver_las_33
  use hecmw_solver_las_44
  use hecmw_solver_las_66
  implicit none

  private
  public :: hecmw_mpc_mat_init
  public :: hecmw_mpc_mat_finalize
  public :: hecmw_mpc_mat_ass
  public :: hecmw_mpc_trans_rhs
  public :: hecmw_mpc_tback_sol

contains

  subroutine hecmw_mpc_mat_init(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD, SOLVER_TYPE

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) then
      hecMATmpc => hecMAT
      return
    endif

    call hecmw_mpc_scale(hecMESH)

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
    if (MPC_METHOD < 1 .or. 3 < MPC_METHOD) MPC_METHOD = 3

    SOLVER_TYPE = hecmw_mat_get_solver_type(hecMAT)
    if (SOLVER_TYPE > 1) then  ! DIRECT SOLVER
      MPC_METHOD = 1
      call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD)
      ! if (MPC_METHOD == 2) then
      !   MPC_METHOD = 3
      !   call hecmw_mat_set_mpc_method(hecMAT, MPC_METHOD)
      ! endif
    endif

    select case (MPC_METHOD)
    case (1)  ! penalty
      hecMATmpc => hecMAT
    case (2)  ! MPCCG
      hecMATmpc => hecMAT
    case (3)  ! elimination
      allocate(hecMATmpc)
      call hecmw_mat_init(hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_init

  subroutine hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMATmpc)
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
      nullify(hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_finalize

  !C===
  !C +-------------+
  !C | MPC Preproc |
  !C +-------------+
  !C===
  subroutine hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: totalmpc, MPC_METHOD, SOLVER_TYPE

    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    if (totalmpc == 0) return

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)

    select case (MPC_METHOD)
    case (1)  ! penalty
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Penalty"
      call hecmw_mat_ass_equation ( hecMESH, hecMATmpc )
    case (2)  ! MPCCG
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: MPC-CG"
    case (3)  ! elimination
      !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Elimination"
      call hecmw_trimatmul_TtKT_mpc(hecMESH, hecMAT, hecMATmpc)
    end select

  end subroutine hecmw_mpc_mat_ass


  subroutine hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
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
    end select

  end subroutine hecmw_mpc_trans_rhs

  !C===
  !C +--------------+
  !C | MPC Postproc |
  !C +--------------+
  !C===
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

    MPC_METHOD = hecmw_mat_get_mpc_method(hecMATmpc)

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


  subroutine hecmw_matvec_set_mpcmatvec_flg(ndof, flag)
    implicit none
    integer(kind=kint), intent(in) :: ndof
    logical, intent(in) :: flag
    select case (ndof)
    case(1)
      !call hecmw_matvec_11_set_mpcmatvec_flg(flag)
      stop 'ERROR: MPC not supported for DOF=1'
    case(2)
      !call hecmw_matvec_22_set_mpcmatvec_flg(flag)
      stop 'ERROR: MPC not supported for DOF=2'
    case(3)
      call hecmw_matvec_33_set_mpcmatvec_flg(flag)
    case(4)
      !call hecmw_matvec_44_set_mpcmatvec_flg(flag)
      stop 'ERROR: MPC not supported for DOF=4'
    case(6)
      !call hecmw_matvec_66_set_mpcmatvec_flg(flag)
      stop 'ERROR: MPC not supported for DOF=6'
    case default
      stop 'ERROR: hecmw_matvec_set_mpcmatvec_flg: unknown DOF'
    end select
  end subroutine hecmw_matvec_set_mpcmatvec_flg


  subroutine hecmw_trans_b(hecMESH, hecMAT, B, BT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), intent(in) :: B(:)
    real(kind=kreal), intent(out), target :: BT(:)
    real(kind=kreal), intent(inout) :: COMMtime
    select case (hecMAT%NDOF)
    case(1)
      !call hecmw_trans_b_11(hecMESH, hecMAT, B, BT, COMMtime)
      stop 'ERROR: MPC not supported for DOF=1'
    case(2)
      !call hecmw_trans_b_22(hecMESH, hecMAT, B, BT, COMMtime)
      stop 'ERROR: MPC not supported for DOF=2'
    case(3)
      call hecmw_trans_b_33(hecMESH, hecMAT, B, BT, COMMtime)
    case(4)
      !call hecmw_trans_b_44(hecMESH, hecMAT, B, BT, COMMtime)
      stop 'ERROR: MPC not supported for DOF=4'
    case(6)
      !call hecmw_trans_b_66(hecMESH, hecMAT, B, BT, COMMtime)
      stop 'ERROR: MPC not supported for DOF=6'
    case default
      stop 'ERROR: hecmw_trans_b: unknown DOF'
    end select
  end subroutine hecmw_trans_b


  subroutine hecmw_tback_x(hecMESH, ndof, X, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(inout) :: X(:)
    real(kind=kreal), intent(inout) :: COMMtime
    select case (ndof)
    case(1)
      !call hecmw_tback_x_11(hecMESH, X, COMMtime)
      stop 'ERROR: MPC not supported for DOF=1'
    case(2)
      !call hecmw_tback_x_22(hecMESH, X, COMMtime)
      stop 'ERROR: MPC not supported for DOF=2'
    case(3)
      call hecmw_tback_x_33(hecMESH, X, COMMtime)
    case(4)
      !call hecmw_tback_x_44(hecMESH, X, COMMtime)
      stop 'ERROR: MPC not supported for DOF=4'
    case(6)
      !call hecmw_tback_x_66(hecMESH, X, COMMtime)
      stop 'ERROR: MPC not supported for DOF=6'
    case default
      stop 'ERROR: hecmw_tback_x: unknown DOF'
    end select
  end subroutine hecmw_tback_x

end module hecmw_mpc_prepost
