!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_SSOR_66
!C***
!C
module hecmw_precond_SSOR_66
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_matrix_ordering_CM
  use m_hecmw_matrix_ordering_MC
  use hecmw_matrix_reorder
  use hecmw_matrix_contact
  !$ use omp_lib

  private

  public:: hecmw_precond_SSOR_66_setup
  public:: hecmw_precond_SSOR_66_apply
  public:: hecmw_precond_SSOR_66_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: D(:) => null()
  real(kind=kreal), pointer :: AL(:) => null()
  real(kind=kreal), pointer :: AU(:) => null()
  integer(kind=kint), pointer :: indexL(:) => null()
  integer(kind=kint), pointer :: indexU(:) => null()
  integer(kind=kint), pointer :: itemL(:) => null()
  integer(kind=kint), pointer :: itemU(:) => null()
  real(kind=kreal), pointer :: ALU(:) => null()

  integer(kind=kint) :: NContact = 0
  real(kind=kreal), pointer :: CAL(:) => null()
  real(kind=kreal), pointer :: CAU(:) => null()
  integer(kind=kint), pointer :: indexCL(:) => null()
  integer(kind=kint), pointer :: indexCU(:) => null()
  integer(kind=kint), pointer :: itemCL(:) => null()
  integer(kind=kint), pointer :: itemCU(:) => null()

  integer(kind=kint) :: NColor
  integer(kind=kint), pointer :: COLORindex(:) => null()
  integer(kind=kint), pointer :: perm(:) => null()
  integer(kind=kint), pointer :: iperm(:) => null()

  logical, save :: isFirst = .true.

contains

  subroutine hecmw_precond_SSOR_66_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint ) :: NPL, NPU, NPCL, NPCU
    real   (kind=kreal), allocatable :: CD(:)
    integer(kind=kint ) :: NCOLOR_IN
    real   (kind=kreal) :: SIGMA_DIAG
    real   (kind=kreal) :: ALUtmp(6,6), PW(6)
    integer(kind=kint ) :: ii, i, j, k
    integer(kind=kint ) :: nthreads = 1
    integer(kind=kint ), allocatable :: perm_tmp(:)
    !real   (kind=kreal) :: t0

    !t0 = hecmw_Wtime()
    !write(*,*) 'DEBUG: SSOR setup start', hecmw_Wtime()-t0

    !$ nthreads = omp_get_max_threads()

    N = hecMAT%N
    ! N = hecMAT%NP
    NCOLOR_IN = hecmw_mat_get_ncolor_in(hecMAT)
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)
    NContact = hecMAT%cmat%n_val

    if (NContact.gt.0) then
      call hecmw_cmat_LU( hecMAT )
    endif

    if (nthreads == 1) then
      NColor = 1
      allocate(COLORindex(0:1), perm(N), iperm(N))
      COLORindex(0) = 0
      COLORindex(1) = N
      do i=1,N
        perm(i) = i
        iperm(i) = i
      end do

      D => hecMAT%D
      AL => hecMAT%AL
      AU => hecMAT%AU
      indexL => hecMAT%indexL
      indexU => hecMAT%indexU
      itemL => hecMAT%itemL
      itemU => hecMAT%itemU
      if (NContact.gt.0) then
        CAL => hecMAT%CAL
        CAU => hecMAT%CAU
        indexCL => hecMAT%indexCL
        indexCU => hecMAT%indexCU
        itemCL => hecMAT%itemCL
        itemCU => hecMAT%itemCU
      end if
    else
      allocate(COLORindex(0:N), perm_tmp(N), perm(N), iperm(N))
      call hecmw_matrix_ordering_RCM(N, hecMAT%indexL, hecMAT%itemL, &
        hecMAT%indexU, hecMAT%itemU, perm_tmp, iperm)
      !write(*,*) 'DEBUG: RCM ordering done', hecmw_Wtime()-t0
      call hecmw_matrix_ordering_MC(N, hecMAT%indexL, hecMAT%itemL, &
        hecMAT%indexU, hecMAT%itemU, perm_tmp, &
        NCOLOR_IN, NColor, COLORindex, perm, iperm)
      !write(*,*) 'DEBUG: MC ordering done', hecmw_Wtime()-t0
      deallocate(perm_tmp)

      !call write_debug_info

      NPL = hecMAT%indexL(N)
      NPU = hecMAT%indexU(N)
      allocate(indexL(0:N), indexU(0:N), itemL(NPL), itemU(NPU))
      call hecmw_matrix_reorder_profile(N, perm, iperm, &
        hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
        indexL, indexU, itemL, itemU)
      !write(*,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0

      call check_ordering

      allocate(D(36*N), AL(36*NPL), AU(36*NPU))
      call hecmw_matrix_reorder_values(N, 6, perm, iperm, &
        hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
        hecMAT%AL, hecMAT%AU, hecMAT%D, &
        indexL, indexU, itemL, itemU, AL, AU, D)
      !write(*,*) 'DEBUG: reordering values done', hecmw_Wtime()-t0

      call hecmw_matrix_reorder_renum_item(N, perm, indexL, itemL)
      call hecmw_matrix_reorder_renum_item(N, perm, indexU, itemU)

      if (NContact.gt.0) then
        NPCL = hecMAT%indexCL(N)
        NPCU = hecMAT%indexCU(N)
        allocate(indexCL(0:N), indexCU(0:N), itemCL(NPCL), itemCU(NPCU))
        call hecmw_matrix_reorder_profile(N, perm, iperm, &
          hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
          indexCL, indexCU, itemCL, itemCU)

        allocate(CD(36*N), CAL(36*NPCL), CAU(36*NPCU))
        call hecmw_matrix_reorder_values(N, 6, perm, iperm, &
          hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
          hecMAT%CAL, hecMAT%CAU, hecMAT%D, &
          indexCL, indexCU, itemCL, itemCU, CAL, CAU, CD)
        deallocate(CD)

        call hecmw_matrix_reorder_renum_item(N, perm, indexCL, itemCL)
        call hecmw_matrix_reorder_renum_item(N, perm, indexCU, itemCU)
      end if

    end if

    allocate(ALU(36*N))
    ALU  = 0.d0

    do ii= 1, 36*N
      ALU(ii) = D(ii)
    enddo

    !    if (NContact.gt.0) then
    !      do k= 1, hecMAT%cmat%n_val
    !        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
    !        ii = iperm( hecMAT%cmat%pair(k)%i )
    !        ALU(9*ii-8) = ALU(9*ii-8) + hecMAT%cmat%pair(k)%val(1, 1)
    !        ALU(9*ii-7) = ALU(9*ii-7) + hecMAT%cmat%pair(k)%val(1, 2)
    !        ALU(9*ii-6) = ALU(9*ii-6) + hecMAT%cmat%pair(k)%val(1, 3)
    !        ALU(9*ii-5) = ALU(9*ii-5) + hecMAT%cmat%pair(k)%val(2, 1)
    !        ALU(9*ii-4) = ALU(9*ii-4) + hecMAT%cmat%pair(k)%val(2, 2)
    !        ALU(9*ii-3) = ALU(9*ii-3) + hecMAT%cmat%pair(k)%val(2, 3)
    !        ALU(9*ii-2) = ALU(9*ii-2) + hecMAT%cmat%pair(k)%val(3, 1)
    !        ALU(9*ii-1) = ALU(9*ii-1) + hecMAT%cmat%pair(k)%val(3, 2)
    !        ALU(9*ii  ) = ALU(9*ii  ) + hecMAT%cmat%pair(k)%val(3, 3)
    !      enddo
    !    endif

    do ii= 1, N
      ALUtmp(1,1)= ALU(36*ii-35) * SIGMA_DIAG
      ALUtmp(1,2)= ALU(36*ii-34)
      ALUtmp(1,3)= ALU(36*ii-33)
      ALUtmp(1,4)= ALU(36*ii-32)
      ALUtmp(1,5)= ALU(36*ii-31)
      ALUtmp(1,6)= ALU(36*ii-30)

      ALUtmp(2,1)= ALU(36*ii-29)
      ALUtmp(2,2)= ALU(36*ii-28) * SIGMA_DIAG
      ALUtmp(2,3)= ALU(36*ii-27)
      ALUtmp(2,4)= ALU(36*ii-26)
      ALUtmp(2,5)= ALU(36*ii-25)
      ALUtmp(2,6)= ALU(36*ii-24)

      ALUtmp(3,1)= ALU(36*ii-23)
      ALUtmp(3,2)= ALU(36*ii-22)
      ALUtmp(3,3)= ALU(36*ii-21) * SIGMA_DIAG
      ALUtmp(3,4)= ALU(36*ii-20)
      ALUtmp(3,5)= ALU(36*ii-19)
      ALUtmp(3,6)= ALU(36*ii-18)

      ALUtmp(4,1)= ALU(36*ii-17)
      ALUtmp(4,2)= ALU(36*ii-16)
      ALUtmp(4,3)= ALU(36*ii-15)
      ALUtmp(4,4)= ALU(36*ii-14) * SIGMA_DIAG
      ALUtmp(4,5)= ALU(36*ii-13)
      ALUtmp(4,6)= ALU(36*ii-12)

      ALUtmp(5,1)= ALU(36*ii-11)
      ALUtmp(5,2)= ALU(36*ii-10)
      ALUtmp(5,3)= ALU(36*ii-9 )
      ALUtmp(5,4)= ALU(36*ii-8 )
      ALUtmp(5,5)= ALU(36*ii-7 ) * SIGMA_DIAG
      ALUtmp(5,6)= ALU(36*ii-6 )

      ALUtmp(6,1)= ALU(36*ii-5 )
      ALUtmp(6,2)= ALU(36*ii-4 )
      ALUtmp(6,3)= ALU(36*ii-3 )
      ALUtmp(6,4)= ALU(36*ii-2 )
      ALUtmp(6,5)= ALU(36*ii-1 )
      ALUtmp(6,6)= ALU(36*ii   ) * SIGMA_DIAG

      do k= 1, 6
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, 6
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, 6
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, 6
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo

      ALU(36*ii-35)= ALUtmp(1,1)
      ALU(36*ii-34)= ALUtmp(1,2)
      ALU(36*ii-33)= ALUtmp(1,3)
      ALU(36*ii-32)= ALUtmp(1,4)
      ALU(36*ii-31)= ALUtmp(1,5)
      ALU(36*ii-30)= ALUtmp(1,6)
      ALU(36*ii-29)= ALUtmp(2,1)
      ALU(36*ii-28)= ALUtmp(2,2)
      ALU(36*ii-27)= ALUtmp(2,3)
      ALU(36*ii-26)= ALUtmp(2,4)
      ALU(36*ii-25)= ALUtmp(2,5)
      ALU(36*ii-24)= ALUtmp(2,6)
      ALU(36*ii-23)= ALUtmp(3,1)
      ALU(36*ii-22)= ALUtmp(3,2)
      ALU(36*ii-21)= ALUtmp(3,3)
      ALU(36*ii-20)= ALUtmp(3,4)
      ALU(36*ii-19)= ALUtmp(3,5)
      ALU(36*ii-18)= ALUtmp(3,6)
      ALU(36*ii-17)= ALUtmp(4,1)
      ALU(36*ii-16)= ALUtmp(4,2)
      ALU(36*ii-15)= ALUtmp(4,3)
      ALU(36*ii-14)= ALUtmp(4,4)
      ALU(36*ii-13)= ALUtmp(4,5)
      ALU(36*ii-12)= ALUtmp(4,6)
      ALU(36*ii-11)= ALUtmp(5,1)
      ALU(36*ii-10)= ALUtmp(5,2)
      ALU(36*ii-9 )= ALUtmp(5,3)
      ALU(36*ii-8 )= ALUtmp(5,4)
      ALU(36*ii-7 )= ALUtmp(5,5)
      ALU(36*ii-6 )= ALUtmp(5,6)
      ALU(36*ii-5 )= ALUtmp(6,1)
      ALU(36*ii-4 )= ALUtmp(6,2)
      ALU(36*ii-3 )= ALUtmp(6,3)
      ALU(36*ii-2 )= ALUtmp(6,4)
      ALU(36*ii-1 )= ALUtmp(6,5)
      ALU(36*ii   )= ALUtmp(6,6)
    enddo

    isFirst = .true.

    !write(*,*) 'DEBUG: SSOR setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_SSOR_66_setup

  subroutine hecmw_precond_SSOR_66_apply(ZP)
    use hecmw_tuning_fx
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint) :: ic, i, iold, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: X1, X2, X3, X4, X5, X6
    real(kind=kreal) :: SW1, SW2, SW3, SW4, SW5, SW6

    ! added for turning >>>
    integer(kind=kint), parameter :: numOfBlockPerThread = 100
    integer(kind=kint), save :: numOfThread = 1, numOfBlock
    integer(kind=kint), save, allocatable :: icToBlockIndex(:)
    integer(kind=kint), save, allocatable :: blockIndexToColorIndex(:)
    integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1
    integer(kind=kint) :: blockIndex, elementCount, numOfElement, ii
    real(kind=kreal) :: numOfElementPerBlock
    integer(kind=kint) :: my_rank

    if (isFirst) then
      !$ numOfThread = omp_get_max_threads()
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
      my_rank = hecmw_comm_get_rank()
      ! write(9000+my_rank,*) &
        !      '# numOfElementPerBlock =', numOfElementPerBlock
      ! write(9000+my_rank,*) &
        !      '# ic, blockIndex, colorIndex, elementCount'
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
            ! write(9000+my_rank,*) ic, blockIndex, &
              !      blockIndexToColorIndex(blockIndex), elementCount
          endif
        enddo
        icToBlockIndex(ic) = blockIndex
      enddo
      numOfBlock = blockIndex

      call hecmw_tuning_fx_calc_sector_cache( N, 3, &
        sectorCacheSize0, sectorCacheSize1 )

      isFirst = .false.
    endif
    ! <<< added for turning

    !call start_collection("loopInPrecond66")

    !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
    !OCL CACHE_SUBSECTOR_ASSIGN(ZP)

    !$omp parallel default(none) &
      !$omp&shared(NColor,indexL,itemL,indexU,itemU,AL,AU,D,ALU,perm,&
      !$omp&       NContact,indexCL,itemCL,indexCU,itemCU,CAL,CAU,&
      !$omp&       ZP,icToBlockIndex,blockIndexToColorIndex) &
      !$omp&private(SW1,SW2,SW3,SW4,SW5,SW6,X1,X2,X3,X4,X5,X6,ic,i,iold,isL,ieL,isU,ieU,j,k,blockIndex)

    !C-- FORWARD
    do ic=1,NColor
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
        do i = blockIndexToColorIndex(blockIndex-1)+1, &
            blockIndexToColorIndex(blockIndex)
          ! do i = startPos(threadNum, ic), endPos(threadNum, ic)
          iold = perm(i)
          SW1= ZP(6*iold-5)
          SW2= ZP(6*iold-4)
          SW3= ZP(6*iold-3)
          SW4= ZP(6*iold-2)
          SW5= ZP(6*iold-1)
          SW6= ZP(6*iold  )
          isL= indexL(i-1)+1
          ieL= indexL(i)
          do j= isL, ieL
            !k= perm(itemL(j))
            k= itemL(j)
            X1= ZP(6*k-5)
            X2= ZP(6*k-4)
            X3= ZP(6*k-3)
            X4= ZP(6*k-2)
            X5= ZP(6*k-1)
            X6= ZP(6*k  )
            SW1= SW1 -AL(36*j-35)*X1 -AL(36*j-34)*X2 -AL(36*j-33)*X3 -AL(36*j-32)*X4 -AL(36*j-31)*X5 -AL(36*j-30)*X6
            SW2= SW2 -AL(36*j-29)*X1 -AL(36*j-28)*X2 -AL(36*j-27)*X3 -AL(36*j-26)*X4 -AL(36*j-25)*X5 -AL(36*j-24)*X6
            SW3= SW3 -AL(36*j-23)*X1 -AL(36*j-22)*X2 -AL(36*j-21)*X3 -AL(36*j-20)*X4 -AL(36*j-19)*X5 -AL(36*j-18)*X6
            SW4= SW4 -AL(36*j-17)*X1 -AL(36*j-16)*X2 -AL(36*j-15)*X3 -AL(36*j-14)*X4 -AL(36*j-13)*X5 -AL(36*j-12)*X6
            SW5= SW5 -AL(36*j-11)*X1 -AL(36*j-10)*X2 -AL(36*j-9 )*X3 -AL(36*j-8 )*X4 -AL(36*j-7 )*X5 -AL(36*j-6 )*X6
            SW6= SW6 -AL(36*j-5 )*X1 -AL(36*j-4 )*X2 -AL(36*j-3 )*X3 -AL(36*j-2 )*X4 -AL(36*j-1 )*X5 -AL(36*j   )*X6
          enddo ! j

          if (NContact.ne.0) then
            isL= indexCL(i-1)+1
            ieL= indexCL(i)
            do j= isL, ieL
              !k= perm(itemCL(j))
              k= itemCL(j)
              X1= ZP(6*k-5)
              X2= ZP(6*k-4)
              X3= ZP(6*k-3)
              X4= ZP(6*k-2)
              X5= ZP(6*k-1)
              X6= ZP(6*k  )
              SW1= SW1 -CAL(36*j-35)*X1 -CAL(36*j-34)*X2 -CAL(36*j-33)*X3 -CAL(36*j-32)*X4 -CAL(36*j-31)*X5 -CAL(36*j-30)*X6
              SW2= SW2 -CAL(36*j-29)*X1 -CAL(36*j-28)*X2 -CAL(36*j-27)*X3 -CAL(36*j-26)*X4 -CAL(36*j-25)*X5 -CAL(36*j-24)*X6
              SW3= SW3 -CAL(36*j-23)*X1 -CAL(36*j-22)*X2 -CAL(36*j-21)*X3 -CAL(36*j-20)*X4 -CAL(36*j-19)*X5 -CAL(36*j-18)*X6
              SW4= SW4 -CAL(36*j-17)*X1 -CAL(36*j-16)*X2 -CAL(36*j-15)*X3 -CAL(36*j-14)*X4 -CAL(36*j-13)*X5 -CAL(36*j-12)*X6
              SW5= SW5 -CAL(36*j-11)*X1 -CAL(36*j-10)*X2 -CAL(36*j-9 )*X3 -CAL(36*j-8 )*X4 -CAL(36*j-7 )*X5 -CAL(36*j-6 )*X6
              SW6= SW6 -CAL(36*j-5 )*X1 -CAL(36*j-4 )*X2 -CAL(36*j-3 )*X3 -CAL(36*j-2 )*X4 -CAL(36*j-1 )*X5 -CAL(36*j   )*X6
            enddo ! j
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X4= SW4
          X5= SW5
          X6= SW6
          X2= X2 -ALU(36*i-29)*X1
          X3= X3 -ALU(36*i-23)*X1 -ALU(36*i-22)*X2
          X4= X4 -ALU(36*i-17)*X1 -ALU(36*i-16)*X2 -ALU(36*i-5)*X3
          X5= X5 -ALU(36*i-11)*X1 -ALU(36*i-10)*X2 -ALU(36*i-9)*X3 -ALU(36*i-8)*X4
          X6= X6 -ALU(36*i-5 )*X1 -ALU(36*i-4 )*X2 -ALU(36*i-3)*X3 -ALU(36*i-2)*X4 -ALU(36*i-1)*X5
          X6= ALU(36*i   )*  X6
          X5= ALU(36*i-7 )*( X5 -ALU(36*i-6)*X6 )
          X4= ALU(36*i-14)*( X4 -ALU(36*i-12)*X6 -ALU(36*i-13)*X5)
          X3= ALU(36*i-21)*( X3 -ALU(36*i-18)*X6 -ALU(36*i-19)*X5 -ALU(36*i-20)*X4)
          X2= ALU(36*i-28)*( X2 -ALU(36*i-24)*X6 -ALU(36*i-25)*X5 -ALU(36*i-26)*X4 -ALU(36*i-27)*X3)
          X1= ALU(36*i-35)*( X1 -ALU(36*i-30)*X6 -ALU(36*i-31)*X5 -ALU(36*i-32)*X4 -ALU(36*i-33)*X3 -ALU(36*i-34)*X2)
          ZP(6*iold-5)= X1
          ZP(6*iold-4)= X2
          ZP(6*iold-3)= X3
          ZP(6*iold-2)= X4
          ZP(6*iold-1)= X5
          ZP(6*iold  )= X6
        enddo ! i
      enddo ! blockIndex
      !$omp end do
    enddo ! ic

    !C-- BACKWARD
    do ic=NColor, 1, -1
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic), icToBlockIndex(ic-1)+1, -1
        do i = blockIndexToColorIndex(blockIndex), &
            blockIndexToColorIndex(blockIndex-1)+1, -1
          ! do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
          !   do i = blockIndexToColorIndex(blockIndex-1)+1, &
            !        blockIndexToColorIndex(blockIndex)
          !   do i = endPos(threadNum, ic), startPos(threadNum, ic), -1
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          SW4= 0.d0
          SW5= 0.d0
          SW6= 0.d0
          isU= indexU(i-1) + 1
          ieU= indexU(i)
          do j= ieU, isU, -1
            !k= perm(itemU(j))
            k= itemU(j)
            X1= ZP(6*k-5)
            X2= ZP(6*k-4)
            X3= ZP(6*k-3)
            X4= ZP(6*k-2)
            X5= ZP(6*k-1)
            X6= ZP(6*k  )
            SW1= SW1 +AU(36*j-35)*X1 +AU(36*j-34)*X2 +AU(36*j-33)*X3 +AU(36*j-32)*X4 +AU(36*j-31)*X5 +AU(36*j-30)*X6
            SW2= SW2 +AU(36*j-29)*X1 +AU(36*j-28)*X2 +AU(36*j-27)*X3 +AU(36*j-26)*X4 +AU(36*j-25)*X5 +AU(36*j-24)*X6
            SW3= SW3 +AU(36*j-23)*X1 +AU(36*j-22)*X2 +AU(36*j-21)*X3 +AU(36*j-20)*X4 +AU(36*j-19)*X5 +AU(36*j-18)*X6
            SW4= SW4 +AU(36*j-17)*X1 +AU(36*j-16)*X2 +AU(36*j-15)*X3 +AU(36*j-14)*X4 +AU(36*j-13)*X5 +AU(36*j-12)*X6
            SW5= SW5 +AU(36*j-11)*X1 +AU(36*j-10)*X2 +AU(36*j-9 )*X3 +AU(36*j-8 )*X4 +AU(36*j-7 )*X5 +AU(36*j-6 )*X6
            SW6= SW6 +AU(36*j-5 )*X1 +AU(36*j-4 )*X2 +AU(36*j-3 )*X3 +AU(36*j-2 )*X4 +AU(36*j-1 )*X5 +AU(36*j   )*X6
          enddo ! j

          if (NContact.gt.0) then
            isU= indexCU(i-1) + 1
            ieU= indexCU(i)
            do j= ieU, isU, -1
              !k= perm(itemCU(j))
              k= itemCU(j)
              X1= ZP(6*k-5)
              X2= ZP(6*k-4)
              X3= ZP(6*k-3)
              X4= ZP(6*k-2)
              X5= ZP(6*k-1)
              X6= ZP(6*k  )
              SW1= SW1 +CAU(36*j-35)*X1 +CAU(36*j-34)*X2 +CAU(36*j-33)*X3 +CAU(36*j-32)*X4 +CAU(36*j-31)*X5 +CAU(36*j-30)*X6
              SW2= SW2 +CAU(36*j-29)*X1 +CAU(36*j-28)*X2 +CAU(36*j-27)*X3 +CAU(36*j-26)*X4 +CAU(36*j-25)*X5 +CAU(36*j-24)*X6
              SW3= SW3 +CAU(36*j-23)*X1 +CAU(36*j-22)*X2 +CAU(36*j-21)*X3 +CAU(36*j-20)*X4 +CAU(36*j-19)*X5 +CAU(36*j-18)*X6
              SW4= SW4 +CAU(36*j-17)*X1 +CAU(36*j-16)*X2 +CAU(36*j-15)*X3 +CAU(36*j-14)*X4 +CAU(36*j-13)*X5 +CAU(36*j-12)*X6
              SW5= SW5 +CAU(36*j-11)*X1 +CAU(36*j-10)*X2 +CAU(36*j-9 )*X3 +CAU(36*j-8 )*X4 +CAU(36*j-7 )*X5 +CAU(36*j-6 )*X6
              SW6= SW6 +CAU(36*j-5 )*X1 +CAU(36*j-4 )*X2 +CAU(36*j-3 )*X3 +CAU(36*j-2 )*X4 +CAU(36*j-1 )*X5 +CAU(36*j   )*X6
            enddo ! j
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X4= SW4
          X5= SW5
          X6= SW6
          X2= X2 -ALU(36*i-29)*X1
          X3= X3 -ALU(36*i-23)*X1 -ALU(36*i-22)*X2
          X4= X4 -ALU(36*i-17)*X1 -ALU(36*i-16)*X2 -ALU(36*i-5)*X3
          X5= X5 -ALU(36*i-11)*X1 -ALU(36*i-10)*X2 -ALU(36*i-9)*X3 -ALU(36*i-8)*X4
          X6= X6 -ALU(36*i-5 )*X1 -ALU(36*i-4 )*X2 -ALU(36*i-3)*X3 -ALU(36*i-2)*X4 -ALU(36*i-1)*X5
          X6= ALU(36*i   )*  X6
          X5= ALU(36*i-7 )*( X5 -ALU(36*i-6)*X6 )
          X4= ALU(36*i-14)*( X4 -ALU(36*i-12)*X6 -ALU(36*i-13)*X5)
          X3= ALU(36*i-21)*( X3 -ALU(36*i-18)*X6 -ALU(36*i-19)*X5 -ALU(36*i-20)*X4)
          X2= ALU(36*i-28)*( X2 -ALU(36*i-24)*X6 -ALU(36*i-25)*X5 -ALU(36*i-26)*X4 -ALU(36*i-27)*X3)
          X1= ALU(36*i-35)*( X1 -ALU(36*i-30)*X6 -ALU(36*i-31)*X5 -ALU(36*i-32)*X4 -ALU(36*i-33)*X3 -ALU(36*i-34)*X2)
          iold = perm(i)
          ZP(6*iold-5)= ZP(6*iold-5) -X1
          ZP(6*iold-4)= ZP(6*iold-4) -X2
          ZP(6*iold-3)= ZP(6*iold-3) -X3
          ZP(6*iold-2)= ZP(6*iold-2) -X4
          ZP(6*iold-1)= ZP(6*iold-1) -X5
          ZP(6*iold  )= ZP(6*iold  ) -X6
        enddo ! i
      enddo ! blockIndex
      !$omp end do
    enddo ! ic
    !$omp end parallel

    !OCL END_CACHE_SUBSECTOR
    !OCL END_CACHE_SECTOR_SIZE

    !call stop_collection("loopInPrecond66")

  end subroutine hecmw_precond_SSOR_66_apply

  subroutine hecmw_precond_SSOR_66_clear(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: nthreads = 1
    !$ nthreads = omp_get_max_threads()
    if (associated(COLORindex)) deallocate(COLORindex)
    if (associated(perm)) deallocate(perm)
    if (associated(iperm)) deallocate(iperm)
    if (associated(ALU)) deallocate(ALU)
    if (nthreads >= 1) then
      if (associated(D)) deallocate(D)
      if (associated(AL)) deallocate(AL)
      if (associated(AU)) deallocate(AU)
      if (associated(indexL)) deallocate(indexL)
      if (associated(indexU)) deallocate(indexU)
      if (associated(itemL)) deallocate(itemL)
      if (associated(itemU)) deallocate(itemU)
      if (NContact.ne.0) then
        if (associated(CAL)) deallocate(CAL)
        if (associated(CAU)) deallocate(CAU)
        if (associated(indexCL)) deallocate(indexCL)
        if (associated(indexCU)) deallocate(indexCU)
        if (associated(itemCL)) deallocate(itemCL)
        if (associated(itemCU)) deallocate(itemCU)
      end if
    end if
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
    if (NContact.ne.0) then
      nullify(CAL)
      nullify(CAU)
      nullify(indexCL)
      nullify(indexCU)
      nullify(itemCL)
      nullify(itemCU)
      call hecmw_cmat_LU_free( hecMAT )
    endif
    NContact = 0
  end subroutine hecmw_precond_SSOR_66_clear

  subroutine write_debug_info
    implicit none
    integer(kind=kint) :: my_rank, ic, in
    my_rank = hecmw_comm_get_rank()
    !--------------------> debug: shizawa
    if (my_rank.eq.0) then
      write(*,*) 'DEBUG: Output fort.19000+myrank and fort.29000+myrank for coloring information'
    endif
    write(19000+my_rank,'(a)') '#NCOLORTot'
    write(19000+my_rank,*) NColor
    write(19000+my_rank,'(a)') '#ic  COLORindex(ic-1)+1  COLORindex(ic)'
    do ic=1,NColor
      write(19000+my_rank,*) ic, COLORindex(ic-1)+1,COLORindex(ic)
    enddo ! ic
    write(29000+my_rank,'(a)') '#n_node'
    write(29000+my_rank,*) N
    write(29000+my_rank,'(a)') '#in  OLDtoNEW(in)  NEWtoOLD(in)'
    do in=1,N
      write(29000+my_rank,*) in, iperm(in), perm(in)
      if (perm(iperm(in)) .ne. in) then
        write(29000+my_rank,*) '** WARNING **: NEWtoOLD and OLDtoNEW: ',in
      endif
    enddo
  end subroutine write_debug_info

  subroutine check_ordering
    implicit none
    integer(kind=kint) :: ic, i, j, k
    integer(kind=kint), allocatable :: iicolor(:)
    ! check color dependence of neighbouring nodes
    if (NColor.gt.1) then
      allocate(iicolor(N))
      do ic=1,NColor
        do i= COLORindex(ic-1)+1, COLORindex(ic)
          iicolor(i) = ic
        enddo ! i
      enddo ! ic
      ! FORWARD: L-part
      do ic=1,NColor
        do i= COLORindex(ic-1)+1, COLORindex(ic)
          do j= indexL(i-1)+1, indexL(i)
            k= itemL(j)
            if (iicolor(i).eq.iicolor(k)) then
              write(*,*) '** ERROR **: L-part: iicolor(i).eq.iicolor(k)',i,k,iicolor(i)
            endif
          enddo ! j
        enddo ! i
      enddo ! ic
      ! BACKWARD: U-part
      do ic=NColor, 1, -1
        do i= COLORindex(ic), COLORindex(ic-1)+1, -1
          do j= indexU(i-1)+1, indexU(i)
            k= itemU(j)
            if (iicolor(i).eq.iicolor(k)) then
              write(*,*) '** ERROR **: U-part: iicolor(i).eq.iicolor(k)',i,k,iicolor(i)
            endif
          enddo ! j
        enddo ! i
      enddo ! ic
      deallocate(iicolor)
    endif ! if (NColor.gt.1)
    !--------------------< debug: shizawa
  end subroutine check_ordering

end module     hecmw_precond_SSOR_66
