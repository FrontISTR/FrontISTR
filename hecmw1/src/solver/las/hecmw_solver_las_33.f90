!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_33
  use hecmw_util
  implicit none

  private

  public :: hecmw_matvec_setup_33
  public :: hecmw_matvec_33
  public :: hecmw_matvec_33_set_async
  public :: hecmw_matvec_33_unset_async
  public :: hecmw_matresid_33
  public :: hecmw_rel_resid_L2_33
  public :: hecmw_Tvec_33
  public :: hecmw_Ttvec_33
  public :: hecmw_TtmatTvec_33
  public :: hecmw_mat_diag_sr_33

  integer(kind=kint), save :: N, NP, NZIN
  integer(kind=kint), allocatable, save :: index_i(:), item_i(:)
  real(kind=kreal), allocatable, save   :: A_i(:)
  ! ! for communication hiding in matvec
  !integer(kind=kint), save :: NZOUT
  ! integer(kind=kint), save, allocatable :: index_o(:), item_o(:)
  ! real(kind=kreal), save, allocatable :: A_o(:)
  logical, save :: async_matvec_flg = .false.

  integer, parameter :: numOfBlockPerThread = 100

contains

  subroutine hecmw_matvec_setup_33 (hecMESH, hecMAT)
    use hecmw_jad_type
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT

    integer(kind=kint) :: i, j, k

    if (hecmw_JAD_IS_INITIALIZED().ne.0) return

    N = hecMAT%N
    NP = hecMAT%NP
    if( allocated(index_i) ) deallocate(index_i)
    if( allocated(item_i) ) deallocate(item_i)
    if( allocated(A_i) ) deallocate(A_i)

    !count nonzero elements
    NZIN = hecMAT%indexL(N) + N + hecMAT%indexU(N)
    !NZOUT = 0
    !! upper
    !do i=1,N
    !  do k=hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
    !    j = hecMAT%itemU(k)
    !    if( j > N ) NZOUT = NZOUT + 1
    !  end do
    !end do
    !NZIN = NZIN - NZOUT

    allocate(index_i(0:N),item_i(NZIN),A_i(9*NZIN))
    !allocate(index_o(0:N),item_o(NZOUT),A_o(9*NZOUT))

    ! set value
    NZIN = 0
    !NZOUT = 0
    index_i(0) = 0
    !index_o(0) = 0
    do i=1,N
      ! lower
      do k=hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
        NZIN = NZIN + 1
        j = hecMAT%itemL(k)
        item_i(NZIN) = j
        A_i(9*NZIN-8:9*NZIN) = hecMAT%AL(9*k-8:9*k)
      end do
      ! diag
      NZIN = NZIN + 1
      item_i(NZIN) = i
      A_i(9*NZIN-8:9*NZIN) = hecMAT%D(9*i-8:9*i)
      ! upper
      do k=hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        !if( j > N ) then
        !  NZOUT = NZOUT + 1
        !  if( OUT_TYPE == 1 ) index_out(NZOUT) = i
        !  item_out(NZOUT) = j
        !  aval_out(9*NZOUT-8:9*NZOUT) = hecMAT%AU(9*k-8:9*k)
        !else
          NZIN = NZIN + 1
          item_i(NZIN) = j
          A_i(9*NZIN-8:9*NZIN) = hecMAT%AU(9*k-8:9*k)
        !end if
      end do

      index_i(i) = NZIN
      !index_o(i) = NZOUT
    end do

  end subroutine hecmw_matvec_setup_33

  !C
  !C***
  !C*** hecmw_matvec_33
  !C***
  !C
  subroutine hecmw_matvec_33 (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
      call hecmw_TtmatTvec_33(hecMESH, hecMAT, X, Y, WK, Tcomm)
      deallocate(WK)
    else
      call hecmw_matvec_33_inner(hecMESH, hecMAT, X, Y, time_Ax, Tcomm)
    endif

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
  end subroutine hecmw_matvec_33

  !C
  !C***
  !C*** hecmw_matvec_33_set_async
  !C***
  !C
  subroutine hecmw_matvec_33_set_async (hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    ! integer(kind=kint) :: i, j, jS, jE, idx, in

    ! allocate(index_o(0:hecMAT%N))
    ! index_o(0) = 0
    ! do i = 1, hecMAT%N
    !   jS= hecMAT%indexU(i-1) + 1
    !   jE= hecMAT%indexU(i  )
    !   idx = index_o(i-1)
    !   do j= jS, jE
    !     in  = hecMAT%itemU(j)
    !     if (in <= hecMAT%N) cycle
    !     idx = idx + 1
    !   enddo
    !   index_o(i) = idx
    ! enddo
    ! allocate(item_o(idx))
    ! allocate(A_o(idx*9))
    ! do i = 1, hecMAT%N
    !   jS= hecMAT%indexU(i-1) + 1
    !   jE= hecMAT%indexU(i  )
    !   idx = index_o(i-1)
    !   do j= jS, jE
    !     in  = hecMAT%itemU(j)
    !     if (in <= hecMAT%N) cycle
    !     idx = idx + 1
    !     item_o(idx) = hecMAT%itemU(j) - hecMAT%N
    !     A_o(9*idx-8:9*idx) = hecMAT%AU(9*j-8:9*j)
    !   enddo
    ! enddo
    ! async_matvec_flg = .true.
  end subroutine hecmw_matvec_33_set_async

  !C
  !C***
  !C*** hecmw_matvec_33_unset_async
  !C***
  !C
  subroutine hecmw_matvec_33_unset_async
    implicit none
    ! if (allocated(index_o)) deallocate(index_o)
    ! if (allocated(item_o)) deallocate(item_o)
    ! if (allocated(A_o)) deallocate(A_o)
    ! async_matvec_flg = .false.
  end subroutine hecmw_matvec_33_unset_async

  subroutine hecmw_matvec_33_core(N,NP,NNZ,X,Y,aval,item,idx,numOfBlock,numOfThread,startPos,endPos,async_matvec_flg)
    !$ use omp_lib
    integer(kind=kint), intent(in)  :: N
    integer(kind=kint), intent(in)  :: NP
    integer(kind=kint), intent(in)  :: NNZ
    real(kind=kreal), intent(in)    :: X(3*NP)
    real(kind=kreal), intent(inout) :: Y(3*N)
    real(kind=kreal), intent(in)    :: aval(9*NNZ)
    integer(kind=kint), intent(in)  :: item(NNZ)
    integer(kind=kint), intent(in)  :: idx(0:N)
    integer(kind=kint), intent(in)  :: numOfBlock
    integer(kind=kint), intent(in)  :: numOfThread
    integer(kind=kint), intent(in)  :: startPos(0:numOfBlock)
    integer(kind=kint), intent(in)  :: endPos(0:numOfBlock)
    logical, intent(in)             :: async_matvec_flg

    integer(kind=kint) :: i, j, jS, jE, in
    real(kind=kreal) :: YV1, YV2, YV3, X1, X2, X3
    integer(kind=kint) :: threadNum, blockNum, blockIndex

    !call fapp_start("loopInMatvec33", 1, 0)
    !call start_collection("loopInMatvec33")

    !OCL CACHE_SUBSECTOR_ASSIGN(X)

    if( numOfThread > 1 ) then

    !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP&PRIVATE(i,X1,X2,X3,YV1,YV2,YV3,j,in,threadNum,blockNum,blockIndex) &
      !$OMP&SHARED(aval,item,idx,X,Y,startPos,endPos,numOfThread,async_matvec_flg)
    threadNum = 0
    !$ threadNum = omp_get_thread_num()
    do blockNum = 0 , numOfBlockPerThread - 1
      blockIndex = blockNum * numOfThread  + threadNum
      do i = startPos(blockIndex), endPos(blockIndex)
        YV1= 0.d0
        YV2= 0.d0
        YV3= 0.d0

        do j= idx(i-1)+1, idx(i)
          in  = item(j)
          X1= X(3*in-2)
          X2= X(3*in-1)
          X3= X(3*in  )
          YV1= YV1 + aval(9*j-8)*X1 + aval(9*j-7)*X2 + aval(9*j-6)*X3
          YV2= YV2 + aval(9*j-5)*X1 + aval(9*j-4)*X2 + aval(9*j-3)*X3
          YV3= YV3 + aval(9*j-2)*X1 + aval(9*j-1)*X2 + aval(9*j  )*X3
        enddo

        Y(3*i-2)= YV1
        Y(3*i-1)= YV2
        Y(3*i  )= YV3
      enddo
    enddo
    !$OMP END PARALLEL

    else

      do i = 1, N
        YV1= 0.d0
        YV2= 0.d0
        YV3= 0.d0

        do j= idx(i-1)+1, idx(i)
          in  = item(j)
          X1= X(3*in-2)
          X2= X(3*in-1)
          X3= X(3*in  )
          YV1= YV1 + aval(9*j-8)*X1 + aval(9*j-7)*X2 + aval(9*j-6)*X3
          YV2= YV2 + aval(9*j-5)*X1 + aval(9*j-4)*X2 + aval(9*j-3)*X3
          YV3= YV3 + aval(9*j-2)*X1 + aval(9*j-1)*X2 + aval(9*j  )*X3
        enddo

        Y(3*i-2)= YV1
        Y(3*i-1)= YV2
        Y(3*i  )= YV3
      enddo

    end if


    !OCL END_CACHE_SUBSECTOR

    !call stop_collection("loopInMatvec33")
    !call fapp_stop("loopInMatvec33", 1, 0)

  end subroutine

  !C
  !C***
  !C*** hecmw_matvec_33_inner ( private subroutine )
  !C***
  !C
  subroutine hecmw_matvec_33_inner (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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

    integer(kind=kint) :: N, NP
    integer(kind=kint), pointer :: indexL(:), itemL(:), indexU(:), itemU(:)
    real(kind=kreal), pointer :: AL(:), AU(:), D(:)

    ! for communication hiding
    integer(kind=kint) :: ireq

    ! added for turning >>>
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

        call hecmw_tuning_fx_calc_sector_cache(NP, 3, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
      ! <<< added for turning

      START_TIME= HECMW_WTIME()
      ! if (async_matvec_flg) then
      !   call hecmw_update_3_R_async (hecMESH, X, NP, ireq)
      ! else
      call hecmw_update_3_R (hecMESH, X, NP)
      ! endif
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)

      call hecmw_matvec_33_core(N,NP,NZIN,X,Y,A_i,item_i,index_i, &
        numOfBlock,numOfThread,startPos,endPos,async_matvec_flg)

      !OCL END_CACHE_SECTOR_SIZE

      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME

      ! if (async_matvec_flg) then
      !   START_TIME= HECMW_WTIME()
      !   call hecmw_update_3_R_wait (hecMESH, ireq)
      !   END_TIME= HECMW_WTIME()
      !   if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      !   START_TIME = hecmw_Wtime()

      !   do i = 1, N
      !     jS= index_o(i-1) + 1
      !     jE= index_o(i  )
      !     if (jS > jE) cycle
      !     YV1= 0.d0
      !     YV2= 0.d0
      !     YV3= 0.d0
      !     do j=jS, jE
      !       in = item_o(j)
      !       X1= X(3*(N+in)-2)
      !       X2= X(3*(N+in)-1)
      !       X3= X(3*(N+in)  )
      !       YV1= YV1 + A_o(9*j-8)*X1 + A_o(9*j-7)*X2 + A_o(9*j-6)*X3
      !       YV2= YV2 + A_o(9*j-5)*X1 + A_o(9*j-4)*X2 + A_o(9*j-3)*X3
      !       YV3= YV3 + A_o(9*j-2)*X1 + A_o(9*j-1)*X2 + A_o(9*j  )*X3
      !     enddo
      !     Y(3*i-2)= Y(3*i-2)+YV1
      !     Y(3*i-1)= Y(3*i-1)+YV2
      !     Y(3*i  )= Y(3*i  )+YV3
      !   enddo

      !   END_TIME = hecmw_Wtime()
      !   time_Ax = time_Ax + END_TIME - START_TIME
      ! endif

    endif

    if (hecMAT%cmat%n_val > 0) then
      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, NP * hecMAT%NDOF )
    end if

  end subroutine hecmw_matvec_33_inner

  !C
  !C***
  !C*** hecmw_matresid_33
  !C***
  !C
  subroutine hecmw_matresid_33 (hecMESH, hecMAT, X, B, R, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:), B(:)
    real(kind=kreal), intent(out) :: R(:)
    real(kind=kreal), intent(inout), optional :: COMMtime

    integer(kind=kint) :: i
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_matvec_33 (hecMESH, hecMAT, X, R, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    !$omp parallel default(none),private(i),shared(hecMAT,R,B)
    !$omp do
    do i = 1, hecMAT%N * 3
      R(i) = B(i) - R(i)
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine hecmw_matresid_33

  !C
  !C***
  !C*** hecmw_rel_resid_L2_33
  !C***
  !C
  function hecmw_rel_resid_L2_33 (hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use hecmw_solver_misc
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_33
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH
    type ( hecmwST_matrix     ), intent(in) :: hecMAT
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
    call hecmw_matresid_33(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    hecmw_rel_resid_L2_33 = sqrt(rnorm2 / bnorm2)

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    deallocate(r)
  end function hecmw_rel_resid_L2_33

  !C
  !C***
  !C*** hecmw_Tvec_33
  !C***
  !C
  subroutine hecmw_Tvec_33 (hecMESH, X, Y, COMMtime)
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
    call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
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
        if (hecMESH%mpc%mpc_dof(j) > 3) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Tvec_33

  !C
  !C***
  !C*** hecmw_Ttvec_33
  !C***
  !C
  subroutine hecmw_Ttvec_33 (hecMESH, X, Y, COMMtime)
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
    call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
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
        if (hecMESH%mpc%mpc_dof(j) > 3) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Ttvec_33

  !C
  !C***
  !C*** hecmw_TtmatTvec_33
  !C***
  !C
  subroutine hecmw_TtmatTvec_33 (hecMESH, hecMAT, X, Y, W, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: COMMtime

    call hecmw_Tvec_33(hecMESH, X, Y, COMMtime)
    call hecmw_matvec_33_inner(hecMESH, hecMAT, Y, W, COMMtime)
    call hecmw_Ttvec_33(hecMESH, W, Y, COMMtime)

  end subroutine hecmw_TtmatTvec_33


  !C
  !C***
  !C*** hecmw_mat_diag_sr_33
  !C***
  !C
  subroutine hecmw_mat_diag_sr_33(hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout), target :: hecMAT
    real(kind=kreal), intent(inout), optional :: COMMtime
    real(kind=kreal), allocatable :: W(:,:)
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: ip
    real(kind=kreal) :: START_TIME, END_TIME, Tcomm
    allocate(W(3*hecMAT%NP,3))
    D => hecMAT%D
    do ip= 1, hecMAT%N
      W(3*ip-2,1)= D(9*ip-8); W(3*ip-2,2)= D(9*ip-7); W(3*ip-2,3)= D(9*ip-6)
      W(3*ip-1,1)= D(9*ip-5); W(3*ip-1,2)= D(9*ip-4); W(3*ip-1,3)= D(9*ip-3)
      W(3*ip  ,1)= D(9*ip-2); W(3*ip  ,2)= D(9*ip-1); W(3*ip  ,3)= D(9*ip  )
    enddo
    START_TIME= HECMW_WTIME()
    call hecmw_update_3_R (hecMESH, W(:,1), hecMAT%NP)
    call hecmw_update_3_R (hecMESH, W(:,2), hecMAT%NP)
    call hecmw_update_3_R (hecMESH, W(:,3), hecMAT%NP)
    END_TIME= HECMW_WTIME()
    if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME
    do ip= hecMAT%N+1, hecMAT%NP
      D(9*ip-8)= W(3*ip-2,1); D(9*ip-7)= W(3*ip-2,2); D(9*ip-6)= W(3*ip-2,3)
      D(9*ip-5)= W(3*ip-1,1); D(9*ip-4)= W(3*ip-1,2); D(9*ip-3)= W(3*ip-1,3)
      D(9*ip-2)= W(3*ip  ,1); D(9*ip-1)= W(3*ip  ,2); D(9*ip  )= W(3*ip  ,3)
    enddo
    deallocate(W)
  end subroutine hecmw_mat_diag_sr_33

end module hecmw_solver_las_33
