!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_66
  use hecmw_util
  implicit none

  private

  public :: hecmw_matvec_66
  public :: hecmw_matresid_66
  public :: hecmw_rel_resid_L2_66
  public :: hecmw_Tvec_66
  public :: hecmw_Ttvec_66
  public :: hecmw_TtmatTvec_66

contains

  !C
  !C***
  !C*** hecmw_matvec_66
  !C***
  !C
  subroutine hecmw_matvec_66 (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
    real(kind=kreal) :: YV1, YV2, YV3, X1, X2, X3
    real(kind=kreal) :: YV4, YV5, YV6, X4, X5, X6

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

    if (hecmw_mat_get_usejad(hecMAT).ne.0) then
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
      if (isFirst .eqv. .true.) then
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

        call hecmw_tuning_fx_calc_sector_cache(NP, 6, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
      ! <<< added for turning

      START_TIME= HECMW_WTIME()
      call hecmw_update_R (hecMESH, X, NP, 6)
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

      !call fapp_start("loopInMatvec66", 1, 0)
      !call start_collection("loopInMatvec66")

      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
      !OCL CACHE_SUBSECTOR_ASSIGN(X)

      !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP&PRIVATE(i,X1,X2,X3,X4,X5,X6,YV1,YV2,YV3,YV4,YV5,YV6,jS,jE,j,in,threadNum,blockNum,blockIndex) &
        !$OMP&SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread)
      threadNum = 0
      !$ threadNum = omp_get_thread_num()
      do blockNum = 0 , numOfBlockPerThread - 1
        blockIndex = blockNum * numOfThread  + threadNum
        do i = startPos(blockIndex), endPos(blockIndex)
          X1= X(6*i-5)
          X2= X(6*i-4)
          X3= X(6*i-3)
          X4= X(6*i-2)
          X5= X(6*i-1)
          X6= X(6*i  )
          YV1= D(36*i-35)*X1 + D(36*i-34)*X2 + D(36*i-33)*X3 + D(36*i-32)*X4 + D(36*i-31)*X5 + D(36*i-30)*X6
          YV2= D(36*i-29)*X1 + D(36*i-28)*X2 + D(36*i-27)*X3 + D(36*i-26)*X4 + D(36*i-25)*X5 + D(36*i-24)*X6
          YV3= D(36*i-23)*X1 + D(36*i-22)*X2 + D(36*i-21)*X3 + D(36*i-20)*X4 + D(36*i-19)*X5 + D(36*i-18)*X6
          YV4= D(36*i-17)*X1 + D(36*i-16)*X2 + D(36*i-15)*X3 + D(36*i-14)*X4 + D(36*i-13)*X5 + D(36*i-12)*X6
          YV5= D(36*i-11)*X1 + D(36*i-10)*X2 + D(36*i-9 )*X3 + D(36*i-8 )*X4 + D(36*i-7 )*X5 + D(36*i-6 )*X6
          YV6= D(36*i-5 )*X1 + D(36*i-4 )*X2 + D(36*i-3 )*X3 + D(36*i-2 )*X4 + D(36*i-1 )*X5 + D(36*i   )*X6

          jS= indexL(i-1) + 1
          jE= indexL(i  )
          do j= jS, jE
            in  = itemL(j)
            X1= X(6*in-5)
            X2= X(6*in-4)
            X3= X(6*in-3)
            X4= X(6*in-2)
            X5= X(6*in-1)
            X6= X(6*in  )
            YV1= YV1 + AL(36*j-35)*X1 + AL(36*j-34)*X2 + AL(36*j-33)*X3 + AL(36*j-32)*X4 + AL(36*j-31)*X5 + AL(36*j-30)*X6
            YV2= YV2 + AL(36*j-29)*X1 + AL(36*j-28)*X2 + AL(36*j-27)*X3 + AL(36*j-26)*X4 + AL(36*j-25)*X5 + AL(36*j-24)*X6
            YV3= YV3 + AL(36*j-23)*X1 + AL(36*j-22)*X2 + AL(36*j-21)*X3 + AL(36*j-20)*X4 + AL(36*j-19)*X5 + AL(36*j-18)*X6
            YV4= YV4 + AL(36*j-17)*X1 + AL(36*j-16)*X2 + AL(36*j-15)*X3 + AL(36*j-14)*X4 + AL(36*j-13)*X5 + AL(36*j-12)*X6
            YV5= YV5 + AL(36*j-11)*X1 + AL(36*j-10)*X2 + AL(36*j-9 )*X3 + AL(36*j-8 )*X4 + AL(36*j-7 )*X5 + AL(36*j-6 )*X6
            YV6= YV6 + AL(36*j-5 )*X1 + AL(36*j-4 )*X2 + AL(36*j-3 )*X3 + AL(36*j-2 )*X4 + AL(36*j-1 )*X5 + AL(36*j   )*X6
          enddo
          jS= indexU(i-1) + 1
          jE= indexU(i  )
          do j= jS, jE
            in  = itemU(j)
            X1= X(6*in-5)
            X2= X(6*in-4)
            X3= X(6*in-3)
            X4= X(6*in-2)
            X5= X(6*in-1)
            X6= X(6*in  )
            YV1= YV1 + AU(36*j-35)*X1 + AU(36*j-34)*X2 + AU(36*j-33)*X3 + AU(36*j-32)*X4 + AU(36*j-31)*X5 + AU(36*j-30)*X6
            YV2= YV2 + AU(36*j-29)*X1 + AU(36*j-28)*X2 + AU(36*j-27)*X3 + AU(36*j-26)*X4 + AU(36*j-25)*X5 + AU(36*j-24)*X6
            YV3= YV3 + AU(36*j-23)*X1 + AU(36*j-22)*X2 + AU(36*j-21)*X3 + AU(36*j-20)*X4 + AU(36*j-19)*X5 + AU(36*j-18)*X6
            YV4= YV4 + AU(36*j-17)*X1 + AU(36*j-16)*X2 + AU(36*j-15)*X3 + AU(36*j-14)*X4 + AU(36*j-13)*X5 + AU(36*j-12)*X6
            YV5= YV5 + AU(36*j-11)*X1 + AU(36*j-10)*X2 + AU(36*j-9 )*X3 + AU(36*j-8 )*X4 + AU(36*j-7 )*X5 + AU(36*j-6 )*X6
            YV6= YV6 + AU(36*j-5 )*X1 + AU(36*j-4 )*X2 + AU(36*j-3 )*X3 + AU(36*j-2 )*X4 + AU(36*j-1 )*X5 + AU(36*j   )*X6
          enddo
          Y(6*i-5)= YV1
          Y(6*i-4)= YV2
          Y(6*i-3)= YV3
          Y(6*i-2)= YV4
          Y(6*i-1)= YV5
          Y(6*i  )= YV6
        enddo
      enddo
      !$OMP END PARALLEL

      !OCL END_CACHE_SUBSECTOR
      !OCL END_CACHE_SECTOR_SIZE

      !call stop_collection("loopInMatvec66")
      !call fapp_stop("loopInMatvec66", 1, 0)

      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME

    endif

    if (hecMAT%cmat%n_val > 0) then
      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, NP * hecMAT%NDOF )
    end if

  end subroutine hecmw_matvec_66

  !C
  !C***
  !C*** hecmw_matresid_66
  !C***
  !C
  subroutine hecmw_matresid_66 (hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
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
    call hecmw_matvec_66 (hecMESH, hecMAT, X, R, time_Ax, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    do i = 1, hecMAT%N * 6
      R(i) = B(i) - R(i)
    enddo

  end subroutine hecmw_matresid_66

  !C
  !C***
  !C*** hecmw_rel_resid_L2_66
  !C***
  !C
  function hecmw_rel_resid_L2_66 (hecMESH, hecMAT, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_solver_misc
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_66
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
    call hecmw_matresid_66(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, time_Ax, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    hecmw_rel_resid_L2_66 = sqrt(rnorm2 / bnorm2)

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    deallocate(r)
  end function hecmw_rel_resid_L2_66

  !C
  !C***
  !C*** hecmw_Tvec_66
  !C***
  !C
  subroutine hecmw_Tvec_66 (hecMESH, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 6)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, hecMESH%nn_internal * hecMESH%n_dof
      Y(i)= X(i)
    enddo

    !    do i= 1, hecMESH%mpc%n_mpc
    !      k = hecMESH%mpc%mpc_index(i-1) + 1
    !      kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
    !      Y(kk) = 0.d0
    !      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
    !        jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
    !        Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
    !      enddo
    !    enddo

  end subroutine hecmw_Tvec_66

  !C
  !C***
  !C*** hecmw_Ttvec_66
  !C***
  !C
  subroutine hecmw_Ttvec_66 (hecMESH, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, 6)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, hecMESH%nn_internal * hecMESH%n_dof
      Y(i)= X(i)
    enddo

    !    do i= 1, hecMESH%mpc%n_mpc
    !      k = hecMESH%mpc%mpc_index(i-1) + 1
    !      kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
    !      Y(kk) = 0.d0
    !      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
    !        jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
    !        Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
    !      enddo
    !    enddo

  end subroutine hecmw_Ttvec_66

  !C
  !C***
  !C*** hecmw_TtmatTvec_66
  !C***
  !C
  subroutine hecmw_TtmatTvec_66 (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout) :: COMMtime

    call hecmw_Tvec_66(hecMESH, X, Y, COMMtime)
    call hecmw_matvec_66(hecMESH, hecMAT, Y, W, time_Ax, COMMtime)
    call hecmw_Ttvec_66(hecMESH, W, Y, COMMtime)

  end subroutine hecmw_TtmatTvec_66


end module hecmw_solver_las_66
