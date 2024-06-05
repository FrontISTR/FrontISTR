!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_44
  use hecmw_util
  implicit none

  private

  public :: hecmw_matvec_44
  public :: hecmw_matvec_44_set_async
  public :: hecmw_matvec_44_unset_async
  public :: hecmw_matresid_44
  public :: hecmw_rel_resid_L2_44
  public :: hecmw_Tvec_44
  public :: hecmw_Ttvec_44
  public :: hecmw_TtmatTvec_44
  public :: hecmw_mat_diag_sr_44

  ! ! for communication hiding in matvec
  ! integer(kind=kint), save, allocatable :: index_o(:), item_o(:)
  ! real(kind=kreal), save, allocatable :: A_o(:)
  logical, save :: async_matvec_flg = .false.

contains

  !C
  !C***
  !C*** hecmw_matvec_44
  !C***
  !C
  subroutine hecmw_matvec_44 (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
      call hecmw_TtmatTvec_44(hecMESH, hecMAT, X, Y, WK, time_Ax, Tcomm)
      deallocate(WK)
    else
      call hecmw_matvec_44_inner(hecMESH, hecMAT, X, Y, time_Ax, Tcomm)
    endif

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
  end subroutine hecmw_matvec_44

  !C
  !C***
  !C*** hecmw_matvec_44_set_async
  !C***
  !C
  subroutine hecmw_matvec_44_set_async (hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT

  end subroutine hecmw_matvec_44_set_async

  !C
  !C***
  !C*** hecmw_matvec_44_unset_async
  !C***
  !C
  subroutine hecmw_matvec_44_unset_async
    implicit none

  end subroutine hecmw_matvec_44_unset_async

  !C
  !C***
  !C*** hecmw_matvec_44_inner ( private subroutine )
  !C***
  !C
  subroutine hecmw_matvec_44_inner (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
    real(kind=kreal) :: YV1, YV2, YV3, YV4, X1, X2, X3, X4

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
            ! write(9000+hecMESH%my_rank,*) mod(blockNum, numOfThread), &
              !      startPos(blockNum), endPos(blockNum)
            blockNum = blockNum + 1
            startPos(blockNum) = i + 1
            if (blockNum == (numOfBlock - 1)) exit
          endif
        enddo
        endPos(blockNum) = N
        ! write(9000+hecMESH%my_rank,*) mod(blockNum, numOfThread), &
          !      startPos(blockNum), endPos(blockNum)
        ! for irregular data
        do i= blockNum+1, numOfBlock-1
          startPos(i) = N
          endPos(i) = N-1
          ! write(9000+hecMESH%my_rank,*) mod(i, numOfThread), &
            !      startPos(i), endPos(i)
        end do

        call hecmw_tuning_fx_calc_sector_cache(NP, 4, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
      ! <<< added for turning

      START_TIME= HECMW_WTIME()

      call hecmw_update_R (hecMESH, X, NP, 4)

      ! endif
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

      !call fapp_start("loopInMatvec44", 1, 0)
      !call start_collection("loopInMatvec44")

      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
      !OCL CACHE_SUBSECTOR_ASSIGN(X)

      !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP&PRIVATE(i,X1,X2,X3,X4,YV1,YV2,YV3,YV4,jS,jE,j,in,threadNum,blockNum,blockIndex) &
        !$OMP&SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread,N,async_matvec_flg)
      threadNum = 0
      !$ threadNum = omp_get_thread_num()
      do blockNum = 0 , numOfBlockPerThread - 1
        blockIndex = blockNum * numOfThread  + threadNum
        do i = startPos(blockIndex), endPos(blockIndex)
          X1= X(4*i-3)
          X2= X(4*i-2)
          X3= X(4*i-1)
          X4= X(4*i  )
          YV1= D(16*i-15)*X1 + D(16*i-14)*X2 + D(16*i-13)*X3 + D(16*i-12)*X4
          YV2= D(16*i-11)*X1 + D(16*i-10)*X2 + D(16*i- 9)*X3 + D(16*i- 8)*X4
          YV3= D(16*i- 7)*X1 + D(16*i- 6)*X2 + D(16*i- 5)*X3 + D(16*i- 4)*X4
          YV4= D(16*i- 3)*X1 + D(16*i- 2)*X2 + D(16*i- 1)*X3 + D(16*i   )*X4

          jS= indexL(i-1) + 1
          jE= indexL(i  )
          do j= jS, jE
            in  = itemL(j)
            X1= X(4*in-3)
            X2= X(4*in-2)
            X3= X(4*in-1)
            X4= X(4*in  )
            YV1= YV1 + AL(16*j-15)*X1 + AL(16*j-14)*X2 + AL(16*j-13)*X3 + AL(16*j-12)*X4
            YV2= YV2 + AL(16*j-11)*X1 + AL(16*j-10)*X2 + AL(16*j- 9)*X3 + AL(16*j- 8)*X4
            YV3= YV3 + AL(16*j- 7)*X1 + AL(16*j- 6)*X2 + AL(16*j- 5)*X3 + AL(16*j- 4)*X4
            YV4= YV4 + AL(16*j- 3)*X1 + AL(16*j- 2)*X2 + AL(16*j- 1)*X3 + AL(16*j   )*X4
          enddo
          jS= indexU(i-1) + 1
          jE= indexU(i  )
          do j= jS, jE
            in  = itemU(j)
            ! if (async_matvec_flg .and. in > N) cycle
            X1= X(4*in-3)
            X2= X(4*in-2)
            X3= X(4*in-1)
            X4= X(4*in  )
            YV1= YV1 + AU(16*j-15)*X1 + AU(16*j-14)*X2 + AU(16*j-13)*X3 + AU(16*j-12)*X4
            YV2= YV2 + AU(16*j-11)*X1 + AU(16*j-10)*X2 + AU(16*j- 9)*X3 + AU(16*j- 8)*X4
            YV3= YV3 + AU(16*j- 7)*X1 + AU(16*j- 6)*X2 + AU(16*j- 5)*X3 + AU(16*j- 4)*X4
            YV4= YV4 + AU(16*j- 3)*X1 + AU(16*j- 2)*X2 + AU(16*j- 1)*X3 + AU(16*j   )*X4
          enddo
          Y(4*i-3)= YV1
          Y(4*i-2)= YV2
          Y(4*i-1)= YV3
          Y(4*i  )= YV4
        enddo
      enddo
      !$OMP END PARALLEL

      !OCL END_CACHE_SUBSECTOR
      !OCL END_CACHE_SECTOR_SIZE

      !call stop_collection("loopInMatvec44")
      !call fapp_stop("loopInMatvec44", 1, 0)

      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME

    endif

    if (hecMAT%cmat%n_val > 0) then
      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, NP * hecMAT%NDOF )
    end if

  end subroutine hecmw_matvec_44_inner

  !C
  !C***
  !C*** hecmw_matresid_44
  !C***
  !C
  subroutine hecmw_matresid_44 (hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:), B(:)
    real(kind=kreal), intent(out) :: R(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout), optional :: COMMtime

    integer(kind=kint) :: i
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_matvec_44 (hecMESH, hecMAT, X, R, time_Ax, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    !$omp parallel default(none),private(i),shared(hecMAT,R,B)
    !$omp do
    do i = 1, hecMAT%N * 4
      R(i) = B(i) - R(i)
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine hecmw_matresid_44

  !C
  !C***
  !C*** hecmw_rel_resid_L2_44
  !C***
  !C
  function hecmw_rel_resid_L2_44 (hecMESH, hecMAT, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_solver_misc
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_44
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH
    type ( hecmwST_matrix     ), intent(in) :: hecMAT
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout), optional :: COMMtime

    real(kind=kreal), allocatable :: r(:)
    real(kind=kreal) :: bnorm2, rnorm2
    real(kind=kreal) :: Tcomm

    allocate(r(hecMAT%NDOF*hecMAT%NP))

    Tcomm = 0.d0
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, &
      hecMAT%B, hecMAT%B, bnorm2, Tcomm)
    if (bnorm2 == 0.d0) then
      bnorm2 = 1.d0
    endif
    call hecmw_matresid_44(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, time_Ax, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    hecmw_rel_resid_L2_44 = sqrt(rnorm2 / bnorm2)

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    deallocate(r)
  end function hecmw_rel_resid_L2_44

  !C
  !C***
  !C*** hecmw_Tvec_44
  !C***
  !C
  subroutine hecmw_Tvec_44 (hecMESH, X, Y, COMMtime)
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
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 4)
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
        if (hecMESH%mpc%mpc_dof(j) > 4) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 4 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 4 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Tvec_44

  !C
  !C***
  !C*** hecmw_Ttvec_44
  !C***
  !C
  subroutine hecmw_Ttvec_44 (hecMESH, X, Y, COMMtime)
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
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 4)
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
        if (hecMESH%mpc%mpc_dof(j) > 4) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 4 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 4 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        !omp atomic
        Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Ttvec_44

  !C
  !C***
  !C*** hecmw_TtmatTvec_44
  !C***
  !C
  subroutine hecmw_TtmatTvec_44 (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout) :: COMMtime

    call hecmw_Tvec_44(hecMESH, X, Y, COMMtime)
    call hecmw_matvec_44_inner(hecMESH, hecMAT, Y, W, time_Ax, COMMtime)
    call hecmw_Ttvec_44(hecMESH, W, Y, COMMtime)

  end subroutine hecmw_TtmatTvec_44


  !C
  !C***
  !C*** hecmw_mat_diag_sr_44
  !C***
  !C
  subroutine hecmw_mat_diag_sr_44(hecMESH, hecMAT, COMMtime)
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
    allocate(W(4*hecMAT%NP,4))
    D => hecMAT%D
    do ip= 1, hecMAT%N
      W(4*ip-3,1)= D(16*ip-15); W(4*ip-3,2)= D(16*ip-14); W(4*ip-3,3)= D(16*ip-13); W(4*ip-3,4)= D(16*ip-12)
      W(4*ip-2,1)= D(16*ip-11); W(4*ip-2,2)= D(16*ip-10); W(4*ip-2,3)= D(16*ip- 9); W(4*ip-2,4)= D(16*ip- 8)
      W(4*ip-1,1)= D(16*ip- 7); W(4*ip-1,2)= D(16*ip- 6); W(4*ip-1,3)= D(16*ip- 5); W(4*ip-1,4)= D(16*ip- 4)
      W(4*ip  ,1)= D(16*ip- 3); W(4*ip  ,2)= D(16*ip- 2); W(4*ip  ,3)= D(16*ip- 1); W(4*ip  ,4)= D(16*ip   )
    enddo
    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, W(:,1), hecMAT%NP, 4)
    call hecmw_update_R (hecMESH, W(:,2), hecMAT%NP, 4)
    call hecmw_update_R (hecMESH, W(:,3), hecMAT%NP, 4)
    call hecmw_update_R (hecMESH, W(:,4), hecMAT%NP, 4)
    END_TIME= HECMW_WTIME()
    if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME
    do ip= hecMAT%N+1, hecMAT%NP
      D(16*ip-15)= W(4*ip-3,1); D(16*ip-14)= W(4*ip-3,2); D(16*ip-13)= W(4*ip-3,3); D(16*ip-12)= W(4*ip-3,4)
      D(16*ip-11)= W(4*ip-2,1); D(16*ip-10)= W(4*ip-2,2); D(16*ip- 9)= W(4*ip-2,3); D(16*ip- 8)= W(4*ip-2,4)
      D(16*ip- 7)= W(4*ip-1,1); D(16*ip- 6)= W(4*ip-1,2); D(16*ip- 5)= W(4*ip-1,3); D(16*ip- 4)= W(4*ip-1,4)
      D(16*ip- 3)= W(4*ip  ,1); D(16*ip- 2)= W(4*ip  ,2); D(16*ip- 1)= W(4*ip  ,3); D(16*ip   )= W(4*ip  ,4)
    enddo
    deallocate(W)
  end subroutine hecmw_mat_diag_sr_44

end module hecmw_solver_las_44
