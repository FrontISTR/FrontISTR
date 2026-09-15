!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_SSOR_11
!C***
!C
module hecmw_precond_SSOR_11
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_matrix_ordering_CM
  use m_hecmw_matrix_ordering_MC
  use hecmw_matrix_reorder
#ifndef _OPENACC
  !$ use omp_lib
#endif

  private

  public:: hecmw_precond_SSOR_11_setup
  public:: hecmw_precond_SSOR_11_apply
  public:: hecmw_precond_SSOR_11_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: D(:) => null()
  real(kind=kreal), pointer :: AL(:) => null()
  real(kind=kreal), pointer :: AU(:) => null()
  integer(kind=kint), pointer :: indexL(:) => null()
  integer(kind=kint), pointer :: indexU(:) => null()
  integer(kind=kint), pointer :: itemL(:) => null()
  integer(kind=kint), pointer :: itemU(:) => null()
  real(kind=kreal), pointer :: ALU(:) => null()

  integer(kind=kint) :: NColor
  integer(kind=kint), pointer :: COLORindex(:) => null()
  integer(kind=kint), pointer :: perm(:) => null()
  integer(kind=kint), pointer :: iperm(:) => null()

  logical, save :: isFirst = .true.

  logical, save :: INITIALIZED = .false.

  ! relaxation parameter of
  !   M = 1/(2-OMEGA) (D/OMEGA + L) (D/OMEGA)^-1 (D/OMEGA + U)
  ! OMEGA = 1 is the classic (D+L) D^-1 (D+U).  _apply is handed only ZP, so setup
  ! leaves the value here for it.
  real(kind=kreal), save :: OMEGA = 1.d0

  ! for tuning
  integer(kind=kint), parameter :: numOfBlockPerThread = 100
  integer(kind=kint), save :: numOfThread = 1, numOfBlock
  integer(kind=kint), save, allocatable :: icToBlockIndex(:)
  integer(kind=kint), save, allocatable :: blockIndexToColorIndex(:)
  integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1

  integer(kind=kint), parameter :: DEBUG = 0

contains

  subroutine hecmw_precond_SSOR_11_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NPL, NPU
    integer(kind=kint ) :: NCOLOR_IN
    real   (kind=kreal) :: OMEGA_INV
    real   (kind=kreal) :: ALUtmp(1,1), PW(1)
    integer(kind=kint ) :: ii, i, j, k
    integer(kind=kint ) :: nthreads = 1
    integer(kind=kint ), allocatable :: perm_tmp(:)
    real   (kind=kreal) :: t0

    if (DEBUG >= 1) then
      t0 = hecmw_Wtime()
      write(*,*) 'DEBUG: SSOR setup start', hecmw_Wtime()-t0
    endif

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 1) then ! need symbolic and numerical setup
        call hecmw_precond_SSOR_11_clear(hecMAT)
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_SSOR_11_clear(hecMAT) ! TEMPORARY
      else
        return
      endif
    endif

#ifndef _OPENACC
    !$ nthreads = omp_get_max_threads()
#endif

    N = hecMAT%N
    ! N = hecMAT%NP
    NCOLOR_IN = hecmw_mat_get_ncolor_in(hecMAT)
    OMEGA = hecmw_mat_get_omega(hecMAT)
    OMEGA_INV = 1.d0 / OMEGA

#ifdef _OPENACC
    allocate(COLORindex(0:N), perm_tmp(N), perm(N), iperm(N))
    call hecmw_matrix_ordering_RCM(N, hecMAT%indexL, hecMAT%itemL, &
      hecMAT%indexU, hecMAT%itemU, perm_tmp, iperm)
    if (DEBUG >= 1) write(*,*) 'DEBUG: RCM ordering done', hecmw_Wtime()-t0
    call hecmw_matrix_ordering_MC(N, hecMAT%indexL, hecMAT%itemL, &
      hecMAT%indexU, hecMAT%itemU, perm_tmp, &
      NCOLOR_IN, NColor, COLORindex, perm, iperm)
    if (DEBUG >= 1) write(*,*) 'DEBUG: MC ordering done', hecmw_Wtime()-t0
    deallocate(perm_tmp)

#else
    if (nthreads == 1) then
      NColor = 1
      allocate(COLORindex(0:1), perm(N), iperm(N))
      COLORindex(0) = 0
      COLORindex(1) = N
      do i=1,N
        perm(i) = i
        iperm(i) = i
      end do
    else
      allocate(COLORindex(0:N), perm_tmp(N), perm(N), iperm(N))
      call hecmw_matrix_ordering_RCM(N, hecMAT%indexL, hecMAT%itemL, &
        hecMAT%indexU, hecMAT%itemU, perm_tmp, iperm)
      if (DEBUG >= 1) write(*,*) 'DEBUG: RCM ordering done', hecmw_Wtime()-t0
      call hecmw_matrix_ordering_MC(N, hecMAT%indexL, hecMAT%itemL, &
        hecMAT%indexU, hecMAT%itemU, perm_tmp, &
        NCOLOR_IN, NColor, COLORindex, perm, iperm)
      if (DEBUG >= 1) write(*,*) 'DEBUG: MC ordering done', hecmw_Wtime()-t0
      deallocate(perm_tmp)

    endif
#endif

    NPL = 0
    do i=1,N
      do j=hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
        if( hecMAT%itemU(j) > N ) exit
        NPL = NPL + 1
      enddo
    enddo
    NPL = max(hecMAT%indexL(N),NPL)
    NPU = hecMAT%indexU(N)
    allocate(indexL(0:N), indexU(0:N), itemL(NPL), itemU(NPU))
    call hecmw_matrix_reorder_profile(N, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      indexL, indexU, itemL, itemU)
    if (DEBUG >= 1) write(*,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0


    allocate(D(N), AL(NPL), AU(NPU))
    call hecmw_matrix_reorder_values(N, 1, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      hecMAT%AL, hecMAT%AU, hecMAT%D, &
      indexL, indexU, itemL, itemU, AL, AU, D)
    if (DEBUG >= 1) write(*,*) 'DEBUG: reordering values done', hecmw_Wtime()-t0

    call hecmw_matrix_reorder_renum_item(N, perm, indexL, itemL)
    call hecmw_matrix_reorder_renum_item(N, perm, indexU, itemU)

    allocate(ALU(N))
    ALU  = 0.d0

    do ii= 1, N
      ALU(ii) = D(ii)
    enddo

#ifdef _OPENACC
    !$acc kernels
    !$acc loop independent private(ALUtmp)
#else
    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,OMEGA_INV)
    !$omp do
#endif
    do ii= 1, N
      ALUtmp(1,1)= ALU(ii) * OMEGA_INV
      ALUtmp(1,1)= 1.d0/ALUtmp(1,1)
      ALU(ii)= ALUtmp(1,1)
    enddo
#ifdef _OPENACC
    !$acc end kernels
#else
    !$omp end do
    !$omp end parallel
#endif

    isFirst = .true.

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

    if (DEBUG >= 1) write(*,*) 'DEBUG: SSOR setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_SSOR_11_setup

  subroutine setup_tuning_parameters
    use hecmw_tuning_fx
    implicit none
    integer(kind=kint) :: blockIndex, elementCount, numOfElement, ii
    real(kind=kreal) :: numOfElementPerBlock
    integer(kind=kint) :: ic, i
    if (DEBUG >= 1) write(*,*) 'DEBUG: setting up tuning parameters for SSOR'
#ifndef _OPENACC
    !$ numOfThread = omp_get_max_threads()
#endif

    numOfBlock = numOfThread * numOfBlockPerThread
    if (allocated(icToBlockIndex)) deallocate(icToBlockIndex)
    if (allocated(blockIndexToColorIndex)) deallocate(blockIndexToColorIndex)
    allocate (icToBlockIndex(0:NColor), &
         blockIndexToColorIndex(0:numOfBlock + NColor))
    numOfElement = N + indexL(N) + indexU(N)
    numOfElementPerBlock = dble(numOfElement) / numOfBlock
    blockIndex = 0
    icToBlockIndex = -1
    icToBlockIndex(0) = 0
    blockIndexToColorIndex = -1
    blockIndexToColorIndex(0) = 0
    do ic = 1, NColor
      elementCount = 0
      ii = 1
      do i = COLORindex(ic-1)+1, COLORindex(ic)
        elementCount = elementCount + 1
        elementCount = elementCount + (indexL(i) - indexL(i-1))
        elementCount = elementCount + (indexU(i) - indexU(i-1))
        if (elementCount > ii * numOfElementPerBlock &
             .or. i == COLORindex(ic)) then
          ii = ii + 1
          blockIndex = blockIndex + 1
          blockIndexToColorIndex(blockIndex) = i
        endif
      enddo
      icToBlockIndex(ic) = blockIndex
    enddo
    numOfBlock = blockIndex

    call hecmw_tuning_fx_calc_sector_cache( N, 1, &
         sectorCacheSize0, sectorCacheSize1 )
  end subroutine setup_tuning_parameters

  subroutine hecmw_precond_SSOR_11_apply(ZP)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint) :: ic, i, iold, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, X1

    real(kind=kreal) :: OMEGA_FAC

    ! added for tuning >>>
    integer(kind=kint) :: blockIndex

#ifndef _OPENACC
    if (isFirst) then
      call setup_tuning_parameters
      isFirst = .false.
    endif
#endif
    ! <<< added for tuning

    OMEGA_FAC = 2.d0 - OMEGA

#ifndef _OPENACC
    !call start_collection("loopInPrecond11")

    !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
    !OCL CACHE_SUBSECTOR_ASSIGN(ZP)

    !$omp parallel default(none) &
      !$omp&shared(NColor,indexL,itemL,indexU,itemU,AL,AU,D,ALU,perm,&
      !$omp&       ZP,icToBlockIndex,blockIndexToColorIndex,OMEGA_FAC) &
      !$omp&private(SW1,X1,ic,i,iold,isL,ieL,isU,ieU,j,k,blockIndex)
#endif

    !C-- FORWARD
    do ic=1,NColor
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent
      do i = COLORindex(ic-1)+1, COLORindex(ic)
#else
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
        do i = blockIndexToColorIndex(blockIndex-1)+1, &
            blockIndexToColorIndex(blockIndex)
#endif
          iold = perm(i)
          SW1= OMEGA_FAC * ZP(iold)
          isL= indexL(i-1)+1
          ieL= indexL(i)
          do j= isL, ieL
            k= itemL(j)
            X1= ZP(k)
            SW1= SW1 - AL(j)*X1
          enddo ! j

          X1= ALU(i  )*  SW1
          ZP(iold)= X1
#ifdef _OPENACC
      enddo
      !$acc end kernels
#else
        enddo ! i
      enddo ! blockIndex
      !$omp end do
#endif
    enddo ! ic

    !C-- BACKWARD
    do ic=NColor, 1, -1
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent
      do i = COLORindex(ic-1)+1, COLORindex(ic)
#else
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic), icToBlockIndex(ic-1)+1, -1
        do i = blockIndexToColorIndex(blockIndex), &
            blockIndexToColorIndex(blockIndex-1)+1, -1
#endif
          ! do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
          !   do i = blockIndexToColorIndex(blockIndex-1)+1, &
            !        blockIndexToColorIndex(blockIndex)
          !   do i = endPos(threadNum, ic), startPos(threadNum, ic), -1
          SW1= 0.d0
          isU= indexU(i-1) + 1
          ieU= indexU(i)
          do j= ieU, isU, -1
            k= itemU(j)
            X1= ZP(k)
            SW1= SW1 + AU(j)*X1
          enddo ! j

          X1= ALU(i)*  SW1

          iold = perm(i)
          ZP(iold)=  ZP(iold) - X1
#ifdef _OPENACC
      enddo
      !$acc end kernels
#else
        enddo ! i
      enddo ! blockIndex
      !$omp end do
#endif
    enddo ! ic
#ifndef _OPENACC
    !$omp end parallel

    !OCL END_CACHE_SUBSECTOR
    !OCL END_CACHE_SECTOR_SIZE

    !call stop_collection("loopInPrecond11")
#endif

  end subroutine hecmw_precond_SSOR_11_apply

  subroutine hecmw_precond_SSOR_11_clear(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    if (associated(COLORindex)) deallocate(COLORindex)
    if (associated(perm)) deallocate(perm)
    if (associated(iperm)) deallocate(iperm)
    if (associated(ALU)) deallocate(ALU)
    if (associated(D)) deallocate(D)
    if (associated(AL)) deallocate(AL)
    if (associated(AU)) deallocate(AU)
    if (associated(indexL)) deallocate(indexL)
    if (associated(indexU)) deallocate(indexU)
    if (associated(itemL)) deallocate(itemL)
    if (associated(itemU)) deallocate(itemU)
    nullify(COLORindex)
    nullify(perm)
    nullify(iperm)
    nullify(ALU)
    nullify(D)
    nullify(AL)
    nullify(AU)
    nullify(indexL)
    nullify(indexU)
    nullify(itemL)
    nullify(itemU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_SSOR_11_clear

end module     hecmw_precond_SSOR_11
