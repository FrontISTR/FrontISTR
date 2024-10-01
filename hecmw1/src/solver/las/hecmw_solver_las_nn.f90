!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_nn
  use hecmw_util
  implicit none

  private

  public :: hecmw_matvec_nn
  public :: hecmw_matvec_nn_set_async
  public :: hecmw_matvec_nn_unset_async
  public :: hecmw_matresid_nn
  public :: hecmw_rel_resid_L2_nn
  public :: hecmw_Tvec_nn
  public :: hecmw_Ttvec_nn
  public :: hecmw_TtmatTvec_nn
  public :: hecmw_mat_diag_sr_nn
  public :: hecmw_mat_add_nn
  public :: hecmw_mat_multiple_nn

  ! ! for communication hiding in matvec
  ! integer(kind=kint), save, allocatable :: index_o(:), item_o(:)
  ! real(kind=kreal), save, allocatable :: A_o(:)
  logical, save :: async_matvec_flg = .false.

contains

  !C
  !C***
  !C*** hecmw_matvec_nn
  !C***
  !C
  subroutine hecmw_matvec_nn (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
      call hecmw_TtmatTvec_nn(hecMESH, hecMAT, X, Y, WK, time_Ax, Tcomm)
      deallocate(WK)
    else
      call hecmw_matvec_nn_inner(hecMESH, hecMAT, X, Y, time_Ax, Tcomm)
    endif

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
  end subroutine hecmw_matvec_nn

  !C
  !C***
  !C*** hecmw_matvec_nn_set_async
  !C***
  !C
  subroutine hecmw_matvec_nn_set_async (hecMAT)
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
  end subroutine hecmw_matvec_nn_set_async

  !C
  !C***
  !C*** hecmw_matvec_nn_unset_async
  !C***
  !C
  subroutine hecmw_matvec_nn_unset_async
    implicit none
    ! if (allocated(index_o)) deallocate(index_o)
    ! if (allocated(item_o)) deallocate(item_o)
    ! if (allocated(A_o)) deallocate(A_o)
    ! async_matvec_flg = .false.
  end subroutine hecmw_matvec_nn_unset_async

  !C
  !C***
  !C*** hecmw_matvec_nn_inner ( private subroutine )
  !C***
  !C
  subroutine hecmw_matvec_nn_inner (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
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
    integer(kind=kint) :: i, j, k, l, jS, jE, in
    real(kind=kreal) :: YV(hecMAT%NDOF), XV(hecMAT%NDOF)

    integer(kind=kint) :: N, NP, NDOF, NDOF2
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
      NDOF =  hecMAT%NDOF
      NDOF2 = NDOF*NDOF

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

        call hecmw_tuning_fx_calc_sector_cache(NP, NDOF, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
      ! <<< added for turning

      START_TIME= HECMW_WTIME()
      ! if (async_matvec_flg) then
      !   call hecmw_update_3_R_async (hecMESH, X, NP, ireq)
      ! else
      call hecmw_update_R (hecMESH, X, NP, NDOF)
      ! endif
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

      !call fapp_start("loopInMatvec33", 1, 0)
      !call start_collection("loopInMatvec33")

      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
      !OCL CACHE_SUBSECTOR_ASSIGN(X)

      !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP&PRIVATE(i,XV,YV,jS,jE,j,k,l,in,threadNum,blockNum,blockIndex) &
        !$OMP&SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread,N,NDOF,NDOF2,async_matvec_flg)
      threadNum = 0
      !$ threadNum = omp_get_thread_num()
      do blockNum = 0 , numOfBlockPerThread - 1
        blockIndex = blockNum * numOfThread  + threadNum
        do i = startPos(blockIndex), endPos(blockIndex)
          do k=1,NDOF
            XV(k) = X(NDOF*(i-1)+k)
          end do
          YV(:)=0.0d0
          do k=1,NDOF
            do l=1,NDOF
              YV(k)=YV(k)+D(NDOF2*(i-1)+(k-1)*NDOF+l)*XV(l)
            end do
          end do
          jS= indexL(i-1) + 1
          jE= indexL(i  )
          do j= jS, jE
            in  = itemL(j)
            do k=1,NDOF
              XV(k) = X(NDOF*(in-1)+k)
            end do
            do k=1,NDOF
              do l=1,NDOF
                YV(k)=YV(k)+AL(NDOF2*(j-1)+(k-1)*NDOF+l)*XV(l)
              end do
            end do
          enddo
          jS= indexU(i-1) + 1
          jE= indexU(i  )
          do j= jS, jE
            in  = itemU(j)
            ! if (async_matvec_flg .and. in > N) cycle
            do k=1,NDOF
              XV(k) = X(NDOF*(in-1)+k)
            end do
            do k=1,NDOF
              do l=1,NDOF
                YV(k)=YV(k)+AU(NDOF2*(j-1)+(k-1)*NDOF+l)*XV(l)
              end do
            end do
          enddo
          do k=1,NDOF
            Y(NDOF*(i-1)+k) = YV(k)
          end do
        enddo
      enddo
      !$OMP END PARALLEL

      !OCL END_CACHE_SUBSECTOR
      !OCL END_CACHE_SECTOR_SIZE

      !call stop_collection("loopInMatvec33")
      !call fapp_stop("loopInMatvec33", 1, 0)

      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME



    endif

    if (hecMAT%cmat%n_val > 0) then
      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, NP * NDOF )
    end if
  end subroutine hecmw_matvec_nn_inner



  !C
  !C***
  !C*** hecmw_matresid_nn
  !C***
  !C
  subroutine hecmw_matresid_nn (hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
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
    call hecmw_matvec_nn (hecMESH, hecMAT, X, R, time_Ax, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    !$omp parallel default(none),private(i),shared(hecMAT,R,B)
    !$omp do
    do i = 1, hecMAT%N * hecMAT%NDOF
      R(i) = B(i) - R(i)
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine hecmw_matresid_nn

  !C
  !C***
  !C*** hecmw_rel_resid_L2_nn
  !C***
  !C
  function hecmw_rel_resid_L2_nn (hecMESH, hecMAT, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_solver_misc
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_nn
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
    call hecmw_matresid_nn(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, time_Ax, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    hecmw_rel_resid_L2_nn = sqrt(rnorm2 / bnorm2)

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    deallocate(r)
  end function hecmw_rel_resid_L2_nn

  !C
  !C***
  !C*** hecmw_Tvec_nn
  !C***
  !C
  subroutine hecmw_Tvec_nn (hecMESH, ndof, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i, j, jj, k, kk

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node, ndof)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    !$omp parallel default(none),private(i,k,kk,j,jj),shared(hecMESH,X,Y),firstprivate(ndof)
    !$omp do
    do i= 1, hecMESH%nn_internal * ndof
      Y(i)= X(i)
    enddo
    !$omp end do

    !$omp do
    OUTER: do i= 1, hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_Tvec_nn

  !C
  !C***
  !C*** hecmw_Ttvec_nn
  !C***
  !C
  subroutine hecmw_Ttvec_nn (hecMESH, ndof, X, Y, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i, j, jj, k, kk

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMESH%n_node,ndof)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    !$omp parallel default(none),private(i,k,kk,j,jj),shared(hecMESH,X,Y),firstprivate(ndof)
    !$omp do
    do i= 1, hecMESH%nn_internal * ndof
      Y(i)= X(i)
    enddo
    !$omp end do

    !$omp do
    OUTER: do i= 1, hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      k = hecMESH%mpc%mpc_index(i-1) + 1
      kk = ndof * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
      Y(kk) = 0.d0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        !$omp atomic
        Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
      enddo
    enddo OUTER
    !$omp end do
    !$omp end parallel


  end subroutine hecmw_Ttvec_nn

  !C
  !C***
  !C*** hecmw_TtmatTvec_nn
  !C***
  !C
  subroutine hecmw_TtmatTvec_nn (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: time_Ax
    real(kind=kreal), intent(inout) :: COMMtime

    call hecmw_Tvec_nn(hecMESH, hecMAT%NDOF, X, Y, COMMtime)
    call hecmw_matvec_nn_inner(hecMESH, hecMAT, Y, W, time_Ax, COMMtime)
    call hecmw_Ttvec_nn(hecMESH, hecMAT%NDOF, W, Y, COMMtime)

  end subroutine hecmw_TtmatTvec_nn

  !C
  !C***
  !C*** hecmw_mat_diag_sr_nn
  !C***
  !C
  subroutine hecmw_mat_diag_sr_nn(hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout), target :: hecMAT
    real(kind=kreal), intent(inout), optional :: COMMtime
    real(kind=kreal), allocatable :: W(:,:)
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: ip, NDOF, i, j
    real(kind=kreal) :: START_TIME, END_TIME
    NDOF = hecMAT%NDOF
    allocate(W(NDOF*hecMAT%NP,NDOF))
    D => hecMAT%D
    do ip= 1, hecMAT%N
      do i=1,NDOF
        do j=1,NDOF
          W(NDOF*(ip-1)+i,j) = D(NDOF*NDOF*(ip-1)+(i-1)*NDOF+j)
        end do
      end do
    enddo
    START_TIME= HECMW_WTIME()
    do i=1,NDOF
      call hecmw_update_R (hecMESH, W(:,i), hecMAT%NP, NDOF)
    end do
    END_TIME= HECMW_WTIME()
    if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME
    do ip= hecMAT%N+1, hecMAT%NP
      do i=1,NDOF
        do j=1,NDOF
          D(NDOF*NDOF*(ip-1)+(i-1)*NDOF+j) = W(NDOF*(ip-1)+i,j)
        end do
      end do
    enddo
    deallocate(W)
  end subroutine hecmw_mat_diag_sr_nn

  subroutine hecmw_mat_add_nn(hecMAT1, hecMAT2, hecMAT3)
    use hecmw_util
    implicit none
    type (hecmwST_matrix)     :: hecMAT1, hecMAT2, hecMAT3
    integer(kind=kint) :: i

    do i = 1, hecMAT1%NP*hecMAT1%NDOF*hecMAT1%NDOF
      hecMAT3%D(i) = hecMAT1%D(i) + hecMAT2%D(i)
    enddo

    do i = 1, hecMAT1%NPU*hecMAT1%NDOF*hecMAT1%NDOF
      hecMAT3%AU(i) = hecMAT1%AU(i) + hecMAT2%AU(i)
    enddo

    do i = 1, hecMAT1%NPL*hecMAT1%NDOF*hecMAT1%NDOF
      hecMAT3%AL(i) = hecMAT1%AL(i) + hecMAT2%AL(i)
    enddo
  end subroutine hecmw_mat_add_nn

  subroutine hecmw_mat_multiple_nn(hecMAT, alpha)
    use hecmw_util
    implicit none
    type (hecmwST_matrix)     :: hecMAT
    real(kind=kreal), intent(in) :: alpha
    integer(kind=kint) :: i

    do i = 1, hecMAT%NP*hecMAT%NDOF*hecMAT%NDOF
      hecMAT%D(i) = alpha*hecMAT%D(i)
    enddo

    do i = 1, hecMAT%NPU*hecMAT%NDOF*hecMAT%NDOF
      hecMAT%AU(i) = alpha*hecMAT%AU(i)
    enddo

    do i = 1, hecMAT%NPL*hecMAT%NDOF*hecMAT%NDOF
      hecMAT%AL(i) = alpha*hecMAT%AL(i)
    enddo
  end subroutine hecmw_mat_multiple_nn
end module hecmw_solver_las_nn
