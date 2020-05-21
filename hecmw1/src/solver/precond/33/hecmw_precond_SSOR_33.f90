!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_SSOR_33
!C***
!C
module hecmw_precond_SSOR_33
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_matrix_ordering_CM
  use m_hecmw_matrix_ordering_MC
  use hecmw_matrix_reorder
  use hecmw_matrix_contact
  !$ use omp_lib

  private

  public:: hecmw_precond_SSOR_33_setup
  public:: hecmw_precond_SSOR_33_apply
  public:: hecmw_precond_SSOR_33_clear

  integer(kind=kint) :: N,NP
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

  logical, save :: INITIALIZED = .false.

  ! for tuning
  integer(kind=kint), parameter :: numOfBlockPerThread = 100
  integer(kind=kint), save :: numOfThread = 1, numOfBlock
  integer(kind=kint), save, allocatable :: icToBlockIndex(:)
  integer(kind=kint), save, allocatable :: blockIndexToColorIndex(:)
  integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1

  integer(kind=kint), parameter :: DEBUG = 0

contains

  subroutine hecmw_precond_SSOR_33_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NPL, NPU, NPCL, NPCU
    real   (kind=kreal), allocatable :: CD(:)
    integer(kind=kint ) :: NCOLOR_IN
    real   (kind=kreal) :: SIGMA_DIAG
    real   (kind=kreal) :: ALUtmp(3,3), PW(3), Atmp(9)
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
        call hecmw_precond_SSOR_33_clear(hecMAT)
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_SSOR_33_clear(hecMAT) ! TEMPORARY
      else
        return
      endif
    endif

    !$ nthreads = omp_get_max_threads()

    N = hecMAT%N
    NP = hecMAT%NP
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

      !call write_debug_info
    endif

    NPL = hecMAT%indexL(N)
    NPU = hecMAT%indexU(N)
    allocate(indexL(0:N), indexU(0:N), itemL(NPL), itemU(NPU))
    call hecmw_matrix_reorder_profile(N, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      indexL, indexU, itemL, itemU)
    if (DEBUG >= 1) write(*,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0

    !call check_ordering

    allocate(D(9*N), AL(9*NPL), AU(9*NPU))
    call hecmw_matrix_reorder_values(N, 3, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      hecMAT%AL, hecMAT%AU, hecMAT%D, &
      indexL, indexU, itemL, itemU, AL, AU, D)
    if (DEBUG >= 1) write(*,*) 'DEBUG: reordering values done', hecmw_Wtime()-t0

    call hecmw_matrix_reorder_renum_item(N, perm, indexL, itemL)
    call hecmw_matrix_reorder_renum_item(N, perm, indexU, itemU)

    if (NContact.gt.0) then
      NPCL = hecMAT%indexCL(N)
      NPCU = hecMAT%indexCU(N)
      allocate(indexCL(0:N), indexCU(0:N), itemCL(NPCL), itemCU(NPCU))
      call hecmw_matrix_reorder_profile(N, perm, iperm, &
        hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
        indexCL, indexCU, itemCL, itemCU)

      allocate(CD(9*N), CAL(9*NPCL), CAU(9*NPCU))
      call hecmw_matrix_reorder_values(N, 3, perm, iperm, &
        hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
        hecMAT%CAL, hecMAT%CAU, hecMAT%D, &
        indexCL, indexCU, itemCL, itemCU, CAL, CAU, CD)
      deallocate(CD)

      call hecmw_matrix_reorder_renum_item(N, perm, indexCL, itemCL)
      call hecmw_matrix_reorder_renum_item(N, perm, indexCU, itemCU)
    else
      NPCL = 1
      NPCU = 1
      allocate(indexCL(0:N), indexCU(0:N), itemCL(NPCL), itemCU(NPCU))
      indexCL(0:N) = 1
      indexCU(0:N) = 1
      allocate(CAL(9*NPCL), CAU(9*NPCU))
    end if

    allocate(ALU(9*N))
    ALU  = 0.d0

    do ii= 1, 9*N
      ALU(ii) = D(ii)
    enddo

    if (NContact.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = iperm( hecMAT%cmat%pair(k)%i )
        ALU(9*ii-8) = ALU(9*ii-8) + hecMAT%cmat%pair(k)%val(1, 1)
        ALU(9*ii-7) = ALU(9*ii-7) + hecMAT%cmat%pair(k)%val(1, 2)
        ALU(9*ii-6) = ALU(9*ii-6) + hecMAT%cmat%pair(k)%val(1, 3)
        ALU(9*ii-5) = ALU(9*ii-5) + hecMAT%cmat%pair(k)%val(2, 1)
        ALU(9*ii-4) = ALU(9*ii-4) + hecMAT%cmat%pair(k)%val(2, 2)
        ALU(9*ii-3) = ALU(9*ii-3) + hecMAT%cmat%pair(k)%val(2, 3)
        ALU(9*ii-2) = ALU(9*ii-2) + hecMAT%cmat%pair(k)%val(3, 1)
        ALU(9*ii-1) = ALU(9*ii-1) + hecMAT%cmat%pair(k)%val(3, 2)
        ALU(9*ii  ) = ALU(9*ii  ) + hecMAT%cmat%pair(k)%val(3, 3)
      enddo
    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW,Atmp),shared(N,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N
      ALUtmp(1,1)= ALU(9*ii-8) * SIGMA_DIAG
      ALUtmp(1,2)= ALU(9*ii-7)
      ALUtmp(1,3)= ALU(9*ii-6)
      ALUtmp(2,1)= ALU(9*ii-5)
      ALUtmp(2,2)= ALU(9*ii-4) * SIGMA_DIAG
      ALUtmp(2,3)= ALU(9*ii-3)
      ALUtmp(3,1)= ALU(9*ii-2)
      ALUtmp(3,2)= ALU(9*ii-1)
      ALUtmp(3,3)= ALU(9*ii  ) * SIGMA_DIAG
      do k= 1, 3
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, 3
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, 3
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, 3
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo
      Atmp(1)= ALUtmp(1,1)
      Atmp(2)= ALUtmp(1,2)
      Atmp(3)= ALUtmp(1,3)
      Atmp(4)= ALUtmp(2,1)
      Atmp(5)= ALUtmp(2,2)
      Atmp(6)= ALUtmp(2,3)
      Atmp(7)= ALUtmp(3,1)
      Atmp(8)= ALUtmp(3,2)
      Atmp(9)= ALUtmp(3,3)

      ALU(9*ii  )= Atmp(9)
      ALU(9*ii-1)= -Atmp(8)*Atmp(9)
      ALU(9*ii-2)= (Atmp(4)*Atmp(8)-Atmp(7))*Atmp(9)
      ALU(9*ii-3)= -Atmp(5)*Atmp(6)*ALU(9*ii)
      ALU(9*ii-4)= Atmp(5)*(1.d0-Atmp(6)*ALU(9*ii-1))
      ALU(9*ii-5)= -Atmp(5)*(Atmp(4)+Atmp(6)*ALU(9*ii-2))
      ALU(9*ii-6)= -Atmp(1)*(Atmp(3)*ALU(9*ii)+Atmp(2)*ALU(9*ii-3))
      ALU(9*ii-7)= -Atmp(1)*(Atmp(3)*ALU(9*ii-1)+Atmp(2)*ALU(9*ii-4))
      ALU(9*ii-8)= Atmp(1)*(1.d0-Atmp(3)*ALU(9*ii-2)-Atmp(2)*ALU(9*ii-5))

      !ALU(9*ii-8)= ALUtmp(1,1)
      !ALU(9*ii-7)= ALUtmp(1,2)
      !ALU(9*ii-6)= ALUtmp(1,3)
      !ALU(9*ii-5)= ALUtmp(2,1)
      !ALU(9*ii-4)= ALUtmp(2,2)
      !ALU(9*ii-3)= ALUtmp(2,3)
      !ALU(9*ii-2)= ALUtmp(3,1)
      !ALU(9*ii-1)= ALUtmp(3,2)
      !ALU(9*ii  )= ALUtmp(3,3)
    enddo
    !$omp end do
    !$omp end parallel

    isFirst = .true.

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

    if (DEBUG >= 1) write(*,*) 'DEBUG: SSOR setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_SSOR_33_setup

  subroutine setup_tuning_parameters
    use hecmw_tuning_fx
    implicit none
    integer(kind=kint) :: blockIndex, elementCount, numOfElement, ii
    real(kind=kreal) :: numOfElementPerBlock
    integer(kind=kint) :: my_rank
    integer(kind=kint) :: ic, i
    if (DEBUG >= 1) write(*,*) 'DEBUG: setting up tuning parameters for SSOR'
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
  end subroutine setup_tuning_parameters



  subroutine hecmw_precond_SSOR_33_apply(ZP)
    use hecmw_tuning_fx
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint) :: ic, i

    ! added for turning >>>
    integer(kind=kint) :: blockIndex

    if (isFirst) then
      call setup_tuning_parameters
      isFirst = .false.
    endif
    ! <<< added for turning

    call hecmw_precond_SSOR_33_apply_inner( &
    N,NP,ZP,AL,AU,D,ALU,itemL,itemU,indexL,indexU,perm,icToBlockIndex,blockIndexToColorIndex, &
    indexL(N),indexU(N),NColor,numOfBlock, &
    NContact,indexCL(N),indexCU(N),indexCL,itemCL,indexCU,itemCU,CAL,CAU)

  end subroutine hecmw_precond_SSOR_33_apply

  subroutine hecmw_precond_SSOR_33_apply_inner( &
    N,NP,ZP,AL,AU,D,ALU,itemL,itemU,indexL,indexU,perm,icToBlockIndex,blockIndexToColorIndex, &
    NPL,NPU,NColor,numOfBlock,NContact,NPCL,NPCU,indexCL,itemCL,indexCU,itemCU,CAL,CAU)
    integer(kind=kint), intent(in)  :: N
    integer(kind=kint), intent(in)  :: NP
    real(kind=kreal), intent(inout) :: ZP(3*NP)
    real(kind=kreal), intent(in)    :: AL(9*NPL)
    real(kind=kreal), intent(in)    :: AU(9*NPU)
    real(kind=kreal), intent(in)    :: D(9*N)
    real(kind=kreal), intent(in)    :: ALU(9*N)
    integer(kind=kint), intent(in)  :: itemL(NPL)
    integer(kind=kint), intent(in)  :: itemU(NPU)
    integer(kind=kint), intent(in)  :: indexL(0:N)
    integer(kind=kint), intent(in)  :: indexU(0:N)
    integer(kind=kint), intent(in)  :: perm(N)
    integer(kind=kint), intent(in)  :: icToBlockIndex(0:NColor)
    integer(kind=kint), intent(in)  :: blockIndexToColorIndex(0:numOfBlock+NColor)
    integer(kind=kint), intent(in)  :: NPL
    integer(kind=kint), intent(in)  :: NPU
    integer(kind=kint), intent(in)  :: NColor
    integer(kind=kint), intent(in)  :: numOfBlock
    integer(kind=kint), intent(in)  :: NContact
    integer(kind=kint), intent(in)  :: NPCL
    integer(kind=kint), intent(in)  :: NPCU
    integer(kind=kint), intent(in)  :: indexCL(0:NPCL)
    integer(kind=kint), intent(in)  :: itemCL(N)
    integer(kind=kint), intent(in)  :: indexCU(0:NPCU)
    integer(kind=kint), intent(in)  :: itemCU(N)
    real(kind=kreal), intent(in)    :: CAL(9*NPL)
    real(kind=kreal), intent(in)    :: CAU(9*NPU)

    integer(kind=kint) :: ic, i, iold, j, isL, ieL, isU, ieU, k
    integer(kind=kint) :: blockIndex
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3

    !call start_collection("loopInPrecond33")

    !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
    !OCL CACHE_SUBSECTOR_ASSIGN(ZP)

    if( NColor > 1 .or. NContact > 0 ) then

    !$omp parallel default(none) &
      !$omp&shared(NColor,indexL,itemL,indexU,itemU,AL,AU,D,ALU,perm,&
      !$omp&       NContact,indexCL,itemCL,indexCU,itemCU,CAL,CAU,&
      !$omp&       ZP,icToBlockIndex,blockIndexToColorIndex) &
      !$omp&private(SW1,SW2,SW3,X1,X2,X3,ic,i,iold,isL,ieL,isU,ieU,j,k,blockIndex)

    !C-- FORWARD
    do ic=1,NColor
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
        do i = blockIndexToColorIndex(blockIndex-1)+1, &
            blockIndexToColorIndex(blockIndex)
          ! do i = startPos(threadNum, ic), endPos(threadNum, ic)
          iold = perm(i)
          SW1= ZP(3*iold-2)
          SW2= ZP(3*iold-1)
          SW3= ZP(3*iold  )
          isL= indexL(i-1)+1
          ieL= indexL(i)
          do j= isL, ieL
            !k= perm(itemL(j))
            k= itemL(j)
            X1= ZP(3*k-2)
            X2= ZP(3*k-1)
            X3= ZP(3*k  )
            SW1= SW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            SW2= SW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            SW3= SW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
          enddo ! j

          if (NContact.ne.0) then
            isL= indexCL(i-1)+1
            ieL= indexCL(i)
            do j= isL, ieL
              !k= perm(itemCL(j))
              k= itemCL(j)
              X1= ZP(3*k-2)
              X2= ZP(3*k-1)
              X3= ZP(3*k  )
              SW1= SW1 - CAL(9*j-8)*X1 - CAL(9*j-7)*X2 - CAL(9*j-6)*X3
              SW2= SW2 - CAL(9*j-5)*X1 - CAL(9*j-4)*X2 - CAL(9*j-3)*X3
              SW3= SW3 - CAL(9*j-2)*X1 - CAL(9*j-1)*X2 - CAL(9*j  )*X3
            enddo ! j
          endif

          ZP(3*iold-2)= ALU(9*i-8)*SW1+ALU(9*i-7)*SW2+ALU(9*i-6)*SW3
          ZP(3*iold-1)= ALU(9*i-5)*SW1+ALU(9*i-4)*SW2+ALU(9*i-3)*SW3
          ZP(3*iold  )= ALU(9*i-2)*SW1+ALU(9*i-1)*SW2+ALU(9*i  )*SW3
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
          isU= indexU(i-1) + 1
          ieU= indexU(i)
          do j= ieU, isU, -1
            !k= perm(itemU(j))
            k= itemU(j)
            X1= ZP(3*k-2)
            X2= ZP(3*k-1)
            X3= ZP(3*k  )
            SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
          enddo ! j

          if (NContact.gt.0) then
            isU= indexCU(i-1) + 1
            ieU= indexCU(i)
            do j= ieU, isU, -1
              !k= perm(itemCU(j))
              k= itemCU(j)
              X1= ZP(3*k-2)
              X2= ZP(3*k-1)
              X3= ZP(3*k  )
              SW1= SW1 + CAU(9*j-8)*X1 + CAU(9*j-7)*X2 + CAU(9*j-6)*X3
              SW2= SW2 + CAU(9*j-5)*X1 + CAU(9*j-4)*X2 + CAU(9*j-3)*X3
              SW3= SW3 + CAU(9*j-2)*X1 + CAU(9*j-1)*X2 + CAU(9*j  )*X3
            enddo ! j
          endif

          iold = perm(i)
          ZP(3*iold-2)= ZP(3*iold-2) - ALU(9*i-8)*SW1 - ALU(9*i-7)*SW2 - ALU(9*i-6)*SW3
          ZP(3*iold-1)= ZP(3*iold-1) - ALU(9*i-5)*SW1 - ALU(9*i-4)*SW2 - ALU(9*i-3)*SW3
          ZP(3*iold  )= ZP(3*iold  ) - ALU(9*i-2)*SW1 - ALU(9*i-1)*SW2 - ALU(9*i  )*SW3
        enddo ! i
      enddo ! blockIndex
      !$omp end do
    enddo ! ic
    !$omp end parallel

    else

      !C-- FORWARD
      do i = 1, N
        iold = perm(i)
        SW1= ZP(3*iold-2)
        SW2= ZP(3*iold-1)
        SW3= ZP(3*iold  )
        isL= indexL(i-1)+1
        ieL= indexL(i)
        do j= isL, ieL
          !k= perm(itemL(j))
          k= itemL(j)
          X1= ZP(3*k-2)
          X2= ZP(3*k-1)
          X3= ZP(3*k  )
          SW1= SW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
          SW2= SW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
          SW3= SW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
        enddo ! j

        ZP(3*iold-2)= ALU(9*i-8)*SW1+ALU(9*i-7)*SW2+ALU(9*i-6)*SW3
        ZP(3*iold-1)= ALU(9*i-5)*SW1+ALU(9*i-4)*SW2+ALU(9*i-3)*SW3
        ZP(3*iold  )= ALU(9*i-2)*SW1+ALU(9*i-1)*SW2+ALU(9*i  )*SW3
      enddo ! i

      !C-- BACKWARD
      do i = N, 1, -1
        SW1= 0.d0
        SW2= 0.d0
        SW3= 0.d0
        isU= indexU(i-1) + 1
        ieU= indexU(i)
        do j= ieU, isU, -1
          !k= perm(itemU(j))
          k= itemU(j)
          X1= ZP(3*k-2)
          X2= ZP(3*k-1)
          X3= ZP(3*k  )
          SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
          SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
          SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
        enddo ! j

        iold = perm(i)
        ZP(3*iold-2)= ZP(3*iold-2) - ALU(9*i-8)*SW1 - ALU(9*i-7)*SW2 - ALU(9*i-6)*SW3
        ZP(3*iold-1)= ZP(3*iold-1) - ALU(9*i-5)*SW1 - ALU(9*i-4)*SW2 - ALU(9*i-3)*SW3
        ZP(3*iold  )= ZP(3*iold  ) - ALU(9*i-2)*SW1 - ALU(9*i-1)*SW2 - ALU(9*i  )*SW3
      enddo ! i

    end if


  end subroutine

  subroutine hecmw_precond_SSOR_33_clear(hecMAT)
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
    INITIALIZED = .false.
  end subroutine hecmw_precond_SSOR_33_clear

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

end module     hecmw_precond_SSOR_33
