!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_SSOR_22
!C***
!C
module hecmw_precond_SSOR_22
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_matrix_ordering_CM
  use m_hecmw_matrix_ordering_MC
  use hecmw_matrix_reorder
  use hecmw_matrix_contact
  !$ use omp_lib

  private

  public:: hecmw_precond_SSOR_22_setup
  public:: hecmw_precond_SSOR_22_apply
  public:: hecmw_precond_SSOR_22_clear

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

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_SSOR_22_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NPL, NPU, NPCL, NPCU
    real   (kind=kreal), allocatable :: CD(:)
    integer(kind=kint ) :: NCOLOR_IN
    real   (kind=kreal) :: SIGMA_DIAG
    real   (kind=kreal) :: ALUtmp(2,2), PW(2)
    integer(kind=kint ) :: ii, i, j, k
    integer(kind=kint ) :: nthreads = 1
    integer(kind=kint ), allocatable :: perm_tmp(:)
    !real   (kind=kreal) :: t0

    !t0 = hecmw_Wtime()
    !write(*,*) 'DEBUG: SSOR setup start', hecmw_Wtime()-t0

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 1) then ! need symbolic and numerical setup
        call hecmw_precond_SSOR_22_clear(hecMAT)
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_SSOR_22_clear(hecMAT) ! TEMPORARY
      else
        return
      endif
    endif

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
    endif

    NPL = hecMAT%indexL(N)
    NPU = hecMAT%indexU(N)
    allocate(indexL(0:N), indexU(0:N), itemL(NPL), itemU(NPU))
    call hecmw_matrix_reorder_profile(N, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      indexL, indexU, itemL, itemU)
    !write(*,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0

    !call check_ordering

    allocate(D(4*N), AL(4*NPL), AU(4*NPU))
    call hecmw_matrix_reorder_values(N, 2, perm, iperm, &
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

      allocate(CD(4*N), CAL(4*NPCL), CAU(4*NPCU))
      call hecmw_matrix_reorder_values(N, 2, perm, iperm, &
        hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
        hecMAT%CAL, hecMAT%CAU, hecMAT%D, &
        indexCL, indexCU, itemCL, itemCU, CAL, CAU, CD)
      deallocate(CD)

      call hecmw_matrix_reorder_renum_item(N, perm, indexCL, itemCL)
      call hecmw_matrix_reorder_renum_item(N, perm, indexCU, itemCU)
    end if

    allocate(ALU(4*N))
    ALU  = 0.d0

    do ii= 1, 4*N
      ALU(ii) = D(ii)
    enddo

    if (NContact.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = iperm( hecMAT%cmat%pair(k)%i )
        ALU(4*ii-3) = ALU(4*ii-3) + hecMAT%cmat%pair(k)%val(1, 1)
        ALU(4*ii-2) = ALU(4*ii-2) + hecMAT%cmat%pair(k)%val(1, 2)
        ALU(4*ii-1) = ALU(4*ii-1) + hecMAT%cmat%pair(k)%val(2, 1)
        ALU(4*ii-0) = ALU(4*ii-0) + hecMAT%cmat%pair(k)%val(2, 2)
      enddo
    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N
      ALUtmp(1,1)= ALU(4*ii-3) * SIGMA_DIAG
      ALUtmp(1,2)= ALU(4*ii-2)
      ALUtmp(2,1)= ALU(4*ii-1)
      ALUtmp(2,2)= ALU(4*ii-0) * SIGMA_DIAG
      do k= 1, 2
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, 2
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, 2
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, 2
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo
      ALU(4*ii-3)= ALUtmp(1,1)
      ALU(4*ii-2)= ALUtmp(1,2)
      ALU(4*ii-1)= ALUtmp(2,1)
      ALU(4*ii-0)= ALUtmp(2,2)
    enddo
    !$omp end do
    !$omp end parallel

    isFirst = .true.

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

    !write(*,*) 'DEBUG: SSOR setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_SSOR_22_setup

  subroutine hecmw_precond_SSOR_22_apply(ZP)
    use hecmw_tuning_fx
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint) :: ic, i, iold, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW(2), X(2)

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

      call hecmw_tuning_fx_calc_sector_cache( N, 2, &
        sectorCacheSize0, sectorCacheSize1 )

      isFirst = .false.
    endif
    ! <<< added for turning

    !call start_collection("loopInPrecond33")

    !OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
    !OCL CACHE_SUBSECTOR_ASSIGN(ZP)

    !$omp parallel default(none) &
      !$omp&shared(NColor,indexL,itemL,indexU,itemU,AL,AU,D,ALU,perm,&
      !$omp&       NContact,indexCL,itemCL,indexCU,itemCU,CAL,CAU,&
      !$omp&       ZP,icToBlockIndex,blockIndexToColorIndex) &
      !$omp&private(SW,X,ic,i,iold,isL,ieL,isU,ieU,j,k,blockIndex)

    !C-- FORWARD
    do ic=1,NColor
      !$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
        do i = blockIndexToColorIndex(blockIndex-1)+1, &
            blockIndexToColorIndex(blockIndex)
          iold = perm(i)
          SW(1)= ZP(2*iold-1)
          SW(2)= ZP(2*iold-0)
          isL= indexL(i-1)+1
          ieL= indexL(i)
          do j= isL, ieL
            k= itemL(j)
            X(1)= ZP(2*k-1)
            X(2)= ZP(2*k-0)
            SW(1)= SW(1) - AL(4*j-3)*X(1) - AL(4*j-2)*X(2)
            SW(2)= SW(2) - AL(4*j-1)*X(1) - AL(4*j-0)*X(2)
          enddo ! j

          if (NContact.ne.0) then
            isL= indexCL(i-1)+1
            ieL= indexCL(i)
            do j= isL, ieL
              k= itemCL(j)
              X(1)= ZP(2*k-1)
              X(2)= ZP(2*k-0)
              SW(1)= SW(1) - CAL(4*j-3)*X(1) - CAL(4*j-2)*X(2)
              SW(2)= SW(2) - CAL(4*j-1)*X(1) - CAL(4*j-0)*X(2)
            enddo ! j
          endif

          X = SW
          X(2)= X(2) - ALU(4*i-1)*X(1)
          X(2)= ALU(4*i  )*  X(2)
          X(1)= ALU(4*i-3)*( X(1) - ALU(4*i-2)*X(2))
          ZP(2*iold-1)= X(1)
          ZP(2*iold-0)= X(2)
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
          SW= 0.d0
          isU= indexU(i-1) + 1
          ieU= indexU(i)
          do j= ieU, isU, -1
            k= itemU(j)
            X(1)= ZP(2*k-1)
            X(2)= ZP(2*k-0)
            SW(1)= SW(1) + AU(4*j-3)*X(1) + AU(4*j-2)*X(2)
            SW(2)= SW(2) + AU(4*j-1)*X(1) + AU(4*j-0)*X(2)
          enddo ! j

          if (NContact.gt.0) then
            isU= indexCU(i-1) + 1
            ieU= indexCU(i)
            do j= ieU, isU, -1
              k= itemCU(j)
              X(1)= ZP(2*k-1)
              X(2)= ZP(2*k-0)
              SW(1)= SW(1) + CAU(4*j-3)*X(1) + CAU(4*j-2)*X(2)
              SW(2)= SW(2) + CAU(4*j-1)*X(1) + CAU(4*j-0)*X(2)
            enddo ! j
          endif

          X = SW
          X(2)= X(2) - ALU(4*i-1)*X(1)


          X(2)= ALU(4*i  )*  X(2)
          X(1)= ALU(4*i-3)*( X(1) - ALU(4*i-2)*X(2) )

          iold = perm(i)
          ZP(2*iold-1)=  ZP(2*iold-1) - X(1)
          ZP(2*iold  )=  ZP(2*iold  ) - X(2)
        enddo ! i
      enddo ! blockIndex
      !$omp end do
    enddo ! ic
    !$omp end parallel

    !OCL END_CACHE_SUBSECTOR
    !OCL END_CACHE_SECTOR_SIZE

    !call stop_collection("loopInPrecond33")

  end subroutine hecmw_precond_SSOR_22_apply

  subroutine hecmw_precond_SSOR_22_clear(hecMAT)
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
  end subroutine hecmw_precond_SSOR_22_clear

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

end module     hecmw_precond_SSOR_22
