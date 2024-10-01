!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_22

  private

  public :: hecmw_matvec_22
  public :: hecmw_matvec_22_set_async
  public :: hecmw_matvec_22_unset_async
  public :: hecmw_matresid_22
  public :: hecmw_rel_resid_L2_22
  public :: hecmw_Tvec_22
  public :: hecmw_Ttvec_22
  public :: hecmw_TtmatTvec_22
  public :: hecmw_mat_diag_sr_22

  ! ! for communication hiding in matvec
  ! integer(kind=kint), save, allocatable :: index_o(:), item_o(:)
  ! real(kind=kreal), save, allocatable :: A_o(:)
  logical, save :: async_matvec_flg = .false.
contains

  !C
  !C***
  !C*** hecmw_matvec_22
  !C***
  !C
  subroutine hecmw_matvec_22 (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_matrix_misc
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout), optional :: COMMtime

    real(kind=kreal) :: Tcomm
    real(kind=kreal), allocatable :: WK(:)

    Tcomm = 0.d0

    if (hecmw_mat_get_flag_mpcmatvec(hecMAT) /= 0) then
      allocate(WK(hecMAT%NP * hecMAT%NDOF))
      call hecmw_TtmatTvec_22(hecMESH, hecMAT, X, Y, WK, time_Ax, Tcomm)
      deallocate(WK)
    else
      call hecmw_matvec_22_inner(hecMESH, hecMAT, X, Y, time_Ax, Tcomm)
    endif

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
  end subroutine hecmw_matvec_22

  subroutine hecmw_matvec_22_set_async (hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
  end subroutine hecmw_matvec_22_set_async
  subroutine hecmw_matvec_22_unset_async
    implicit none
  end subroutine hecmw_matvec_22_unset_async

  !C
  !C***
  !C*** hecmw_matvec_22_inner ( private subroutine )
  !C***
  !C
  subroutine hecmw_matvec_22_inner (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    use hecmw_matrix_contact
    use hecmw_matrix_misc
    use hecmw_jad_type
    use hecmw_tuning_fx
    !$ use omp_lib

    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout), optional :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME, Tcomm
    integer(kind=kint) :: i, j, jS, jE, in
    real(kind=kreal) :: YV1, YV2, X1, X2

    integer(kind=kint) :: N, NP
    integer(kind=kint), pointer :: indexL(:), itemL(:), indexU(:), itemU(:)
    real(kind=kreal), pointer :: AL(:), AU(:), D(:)

    ! added for turning >>>
    integer, parameter :: numOfBlockPerThread = 100
    logical, save :: isFirst = .true.
    integer, save :: numOfThread = 1
    integer, save, allocatable :: startPos(:), endPos(:)
    integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1
    integer(kind=kint) :: threadNum, blockNum, numOfBlock
    integer(kind=kint) :: numOfElement, elementCount, blockIndex
    real(kind=kreal) :: numOfElementPerBlock
    ! <<< added for turning

    if (hecmw_JAD_IS_INITIALIZED().ne.0) then
      Tcomm = 0.d0
      START_TIME = hecmw_Wtime()
      call hecmw_JAD_MATVEC(hecMESH, hecMAT, X, Y, Tcomm)
      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME - Tcomm
      if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    else

      N = hecMAT%N
      NP = hecMAT%NP
      indexL => hecMAT%indexL
      indexU => hecMAT%indexU
      itemL => hecMAT%itemL
      itemU => hecMAT%itemU
      AL => hecMAT%AL
      AU => hecMAT%AU
      D => hecMAT%D

      ! added for turning >>>
      if (.not. isFirst) then
        numOfBlock = numOfThread * numOfBlockPerThread
        if (endPos(numOfBlock-1) .ne. N-1) then
          deallocate(startPos, endPos)
          isFirst = .true.
        endif
      endif
      if (isFirst) then
        !$ numOfThread = omp_get_max_threads()
        numOfBlock = numOfThread * numOfBlockPerThread
        allocate (startPos(0 : numOfBlock - 1), endPos(0 : numOfBlock - 1))
        numOfElement = N + indexL(N) + indexU(N)
        numOfElementPerBlock = dble(numOfElement) / numOfBlock
        blockNum = 0
        elementCount = 0
        startPos(blockNum) = 1
        do i= 1, N
          elementCount = elementCount + 1
          elementCount = elementCount + (indexL(i) - indexL(i-1))
          elementCount = elementCount + (indexU(i) - indexU(i-1))
          if (elementCount > (blockNum + 1) * numOfElementPerBlock) then
            endPos(blockNum) = i
            blockNum = blockNum + 1
            startPos(blockNum) = i + 1
            if (blockNum == (numOfBlock - 1)) exit
          endif
        enddo
        endPos(blockNum) = N
        ! for irregular data
        do i= blockNum+1, numOfBlock-1
          startPos(i) = N
          endPos(i) = N-1
        end do

        call hecmw_tuning_fx_calc_sector_cache(NP, 2, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
      ! <<< added for turning

      START_TIME= HECMW_WTIME()
      call hecmw_update_R (hecMESH, X, NP, 2)
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

      !call fapp_start("loopInMatvec33", 1, 0)
      !call start_collection("loopInMatvec33")

      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
      !OCL CACHE_SUBSECTOR_ASSIGN(X)

      !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP&PRIVATE(i,X1,X2,YV1,YV2,jS,jE,j,in,threadNum,blockNum,blockIndex) &
        !$OMP&SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread,N,async_matvec_flg)
      threadNum = 0
      !$ threadNum = omp_get_thread_num()
      do blockNum = 0 , numOfBlockPerThread - 1
        blockIndex = blockNum * numOfThread  + threadNum
        do i = startPos(blockIndex), endPos(blockIndex)
          X1= X(2*i-1)
          X2= X(2*i  )
          YV1= D(4*i-3)*X1 + D(4*i-2)*X2
          YV2= D(4*i-1)*X1 + D(4*i  )*X2

          jS= indexL(i-1) + 1
          jE= indexL(i  )
          do j= jS, jE
            in  = itemL(j)
            X1 = X(2*in-1)
            X2 = X(2*in  )
            YV1= YV1 + AL(4*j-3)*X1 + AL(4*j-2)*X2
            YV2= YV2 + AL(4*j-1)*X1 + AL(4*j  )*X2
          enddo
          jS= indexU(i-1) + 1
          jE= indexU(i  )
          do j= jS, jE
            in  = itemU(j)
            X1 = X(2*in-1)
            X2 = X(2*in  )
            YV1= YV1 + AU(4*j-3)*X1 + AU(4*j-2)*X2
            YV2= YV2 + AU(4*j-1)*X1 + AU(4*j  )*X2
          enddo
          Y(2*i-1)= YV1
          Y(2*i  )= YV2
        enddo
      enddo

      !$OMP END PARALLEL

      !OCL END_CACHE_SUBSECTOR
      !OCL END_CACHE_SECTOR_SIZE


      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME


    endif

    if (hecMAT%cmat%n_val > 0) then
      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, NP * hecMAT%NDOF )
    end if

  end subroutine hecmw_matvec_22_inner






  !C
  !C***
  !C*** hecmw_matresid_22
  !C***
  !C
  subroutine hecmw_matresid_22 (hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
    use hecmw_util

    implicit none
    real(kind=kreal) :: X(:), B(:), R(:)
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    real(kind=kreal) :: time_Ax
    real(kind=kreal), optional :: COMMtime

    integer(kind=kint) :: i
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_matvec_22 (hecMESH, hecMAT, X, R, time_Ax, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    !$omp parallel default(none),private(i),shared(hecMAT,R,B)
    !$omp do
    do i = 1, hecMAT%N * 2
      R(i) = B(i) - R(i)
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine hecmw_matresid_22

  !C
  !C***
  !C*** hecmw_rel_resid_L2_22
  !C***
  !C
  function hecmw_rel_resid_L2_22 (hecMESH, hecMAT, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_solver_misc

    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_22
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH
    type ( hecmwST_matrix     ), intent(in) :: hecMAT
    real(kind=kreal) :: time_Ax
    real(kind=kreal), optional :: COMMtime

    real(kind=kreal) :: r(hecMAT%NDOF*hecMAT%NP)
    real(kind=kreal) :: bnorm2, rnorm2
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, Tcomm)
    if (bnorm2 == 0.d0) then
      bnorm2 = 1.d0
    endif
    call hecmw_matresid_22(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, time_Ax, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    hecmw_rel_resid_L2_22 = sqrt(rnorm2 / bnorm2)

  end function hecmw_rel_resid_L2_22
  !C
  !C***
  !C*** hecmw_Tvec_22
  !C***
  !C
  subroutine hecmw_Tvec_22 (hecMESH, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i, j, jj, k, kk

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 2)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    !$omp parallel default(none),private(i,k,kk,j,jj),shared(hecMESH,X,Y)
    !$omp do
    do i= 1, hecMESH%nn_internal * hecMESH%n_dof
      Y(i)= X(i)
    enddo
    !$omp end do

    !$omp do
    OUTER: do i= 1, hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > 2) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 2 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 2 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Tvec_22

  !C
  !C***
  !C*** hecmw_Ttvec_22
  !C***
  !C
  subroutine hecmw_Ttvec_22 (hecMESH, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i, j, jj, k, kk

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 2)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    !$omp parallel default(none),private(i,k,kk,j,jj),shared(hecMESH,X,Y)
    !$omp do
    do i= 1, hecMESH%nn_internal * hecMESH%n_dof
      Y(i)= X(i)
    enddo
    !$omp end do

    !$omp do
    OUTER: do i= 1, hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > 2) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 2 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 2 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        !omp atomic
        Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Ttvec_22

  !C
  !C***
  !C*** hecmw_TtmatTvec_22
  !C***
  !C
  subroutine hecmw_TtmatTvec_22 (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout) :: COMMtime

    call hecmw_Tvec_22(hecMESH, X, Y, COMMtime)
    call hecmw_matvec_22_inner(hecMESH, hecMAT, Y, W, time_Ax, COMMtime)
    call hecmw_Ttvec_22(hecMESH, W, Y, COMMtime)

  end subroutine hecmw_TtmatTvec_22


  !C
  !C***
  !C*** hecmw_mat_diag_sr_22
  !C***
  !C
  subroutine hecmw_mat_diag_sr_22(hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout), target :: hecMAT
    real(kind=kreal), intent(inout), optional :: COMMtime
    real(kind=kreal), allocatable :: W(:,:)
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: ip
    real(kind=kreal) :: START_TIME, END_TIME
    allocate(W(2*hecMAT%NP,2))
    D => hecMAT%D
    do ip= 1, hecMAT%N
      W(2*ip-1,1)= D(4*ip-3); W(2*ip-1,2)= D(4*ip-2);
      W(2*ip-0,1)= D(4*ip-1); W(2*ip-0,2)= D(4*ip-0);
    enddo
    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, W(:,1), hecMAT%NP, 2)
    call hecmw_update_R (hecMESH, W(:,2), hecMAT%NP, 2)
    END_TIME= HECMW_WTIME()
    if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME
    do ip= hecMAT%N+1, hecMAT%NP
      D(4*ip-3)= W(2*ip-1,1); D(4*ip-2)= W(2*ip-1,2);
      D(4*ip-1)= W(2*ip-0,1); D(4*ip-0)= W(2*ip-0,2);
    enddo
    deallocate(W)
  end subroutine hecmw_mat_diag_sr_22
end module hecmw_solver_las_22
