!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las_11
  use hecmw_util
  implicit none

  private

  public :: hecmw_matvec_11
  public :: hecmw_matresid_11
  public :: hecmw_rel_resid_L2_11

contains

  !C
  !C***
  !C*** hecmw_matvec_11
  !C***
  !C
  subroutine hecmw_matvec_11 (hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f
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
    real(kind=kreal) :: YV

    integer(kind=kint) :: N, NP
    integer(kind=kint), pointer :: indexL(:), itemL(:), indexU(:), itemU(:), indexA(:), itemA(:)
    real(kind=kreal), pointer :: AL(:), AU(:), D(:), A(:)

    ! added for tuning >>>
    integer, parameter :: numOfBlockPerThread = 100
    logical, save :: isFirst = .true.
    integer, save :: numOfThread = 1
    integer, save, allocatable :: startPos(:), endPos(:)
    integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1
    integer(kind=kint) :: threadNum, blockNum, numOfBlock
    integer(kind=kint) :: numOfElement, elementCount, blockIndex
    real(kind=kreal) :: numOfElementPerBlock
    ! <<< added for tuning

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
      indexA => hecMAT%indexA
      itemL => hecMAT%itemL
      itemU => hecMAT%itemU
      itemA => hecMAT%itemA
      AL => hecMAT%AL
      AU => hecMAT%AU
      D => hecMAT%D
      A => hecMAT%A

      ! added for tuning >>>
#ifndef _OPENACC
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

        call hecmw_tuning_fx_calc_sector_cache(NP, 1, &
          sectorCacheSize0, sectorCacheSize1)

        isFirst = .false.
      endif
#endif
      ! <<< added for tuning

      START_TIME= HECMW_WTIME()
      call hecmw_update_R (hecMESH, X, NP, 1)
      END_TIME= HECMW_WTIME()
      if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

      START_TIME = hecmw_Wtime()

#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent
      do i = 1, N
        YV= 0
        jS= indexA(i-1) + 1
        jE= indexA(i)
        do j= jS, jE
          in= itemA(j)
          YV= YV + A(j) * X(in)
        enddo
        Y(i)= YV
      enddo
      !$acc end kernels
#else
      !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
      !OCL CACHE_SUBSECTOR_ASSIGN(X)

      !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP&PRIVATE(i,YV,jS,jE,j,in,threadNum,blockNum,blockIndex) &
        !$OMP&SHARED(D,AL,AU,indexL,itemL,indexU,itemU,X,Y,startPos,endPos,numOfThread,N)
      threadNum = 0
      !$ threadNum = omp_get_thread_num()
      do blockNum = 0 , numOfBlockPerThread - 1
        blockIndex = blockNum * numOfThread  + threadNum
        do i = startPos(blockIndex), endPos(blockIndex)
          YV= D(i) * X(i)

          jS= indexL(i-1) + 1
          jE= indexL(i  )
          do j= jS, jE
            in= itemL(j)
            YV= YV + AL(j) * X(in)
          enddo
          jS= indexU(i-1) + 1
          jE= indexU(i  )
          do j= jS, jE
            in= itemU(j)
            YV= YV + AU(j) * X(in)
          enddo
          Y(i)= YV
        enddo
      enddo
      !$OMP END PARALLEL

      !OCL END_CACHE_SUBSECTOR
      !OCL END_CACHE_SECTOR_SIZE
#endif

      END_TIME = hecmw_Wtime()
      time_Ax = time_Ax + END_TIME - START_TIME

    endif
  end subroutine hecmw_matvec_11

  !C
  !C***
  !C*** hecmw_matresid_11
  !C***
  !C
  subroutine hecmw_matresid_11 (hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
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
    call hecmw_matvec_11 (hecMESH, hecMAT, X, R, time_Ax, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
#ifdef _OPENACC
    !$acc kernels
    !$acc loop independent
#else
    !$omp parallel default(none),private(i),shared(hecMAT,R,B)
    !$omp do
#endif
    do i = 1, hecMAT%N
      R(i) = B(i) - R(i)
    enddo
#ifdef _OPENACC
    !$acc end kernels
#else
    !$omp end do
    !$omp end parallel
#endif
  end subroutine hecmw_matresid_11

  !C
  !C***
  !C*** hecmw_rel_resid_L2_11
  !C***
  !C
  function hecmw_rel_resid_L2_11 (hecMESH, hecMAT, time_Ax, COMMtime)
    use hecmw_util
    use hecmw_solver_misc
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_11
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
    call hecmw_matresid_11(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, time_Ax, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    hecmw_rel_resid_L2_11 = sqrt(rnorm2 / bnorm2)

    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    deallocate(r)
  end function hecmw_rel_resid_L2_11

end module hecmw_solver_las_11
