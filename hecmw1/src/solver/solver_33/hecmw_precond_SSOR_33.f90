!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.5                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

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
  !$ use omp_lib

  private

  public:: hecmw_precond_SSOR_33_setup
  public:: hecmw_precond_SSOR_33_apply
  public:: hecmw_precond_SSOR_33_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()
  real(kind=kreal), pointer :: AU(:) => null()
  real(kind=kreal), pointer :: AL(:) => null()
  real(kind=kreal), pointer :: D(:) => null()
  integer(kind=kint), pointer :: indexL(:) => null()
  integer(kind=kint), pointer :: itemL(:) => null()
  integer(kind=kint), pointer :: indexU(:) => null()
  integer(kind=kint), pointer :: itemU(:) => null()

  integer(kind=kint) :: NContact = 0
  real(kind=kreal), pointer :: CAU(:) => null()
  real(kind=kreal), pointer :: CAL(:) => null()
  integer(kind=kint), pointer :: indexCL(:) => null()
  integer(kind=kint), pointer :: itemCL(:) => null()
  integer(kind=kint), pointer :: indexCU(:) => null()
  integer(kind=kint), pointer :: itemCU(:) => null()

  integer(kind=kint) :: NColor
  integer(kind=kint), pointer :: COLORindex(:) => null()
  integer(kind=kint), pointer :: perm(:) => null()
  integer(kind=kint), pointer :: iperm(:) => null()

contains

  subroutine hecmw_precond_SSOR_33_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint ) :: NP, NPL, NPU, NPCL, NPCU
    real   (kind=kreal), allocatable :: CD(:)
    integer(kind=kint ) :: NCOLOR_IN
    real   (kind=kreal) :: SIGMA_DIAG
    real   (kind=kreal) :: ALUtmp(3,3), PW(3)
    integer(kind=kint ) :: ii, i, j, k
    integer(kind=kint ) :: nthreads = 1
    integer(kind=kint ), allocatable :: perm_tmp(:)
    real   (kind=kreal) :: t0

    t0 = hecmw_Wtime()
    write(*,*) 'DEBUG: SSOR setup start', hecmw_Wtime()-t0

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
      allocate(COLORindex(0:1))
      COLORindex(0) = 0
      COLORindex(1) = N

      AU => hecMAT%AU
      AL => hecMAT%AL
      D  => hecMAT%D
      indexL => hecMAT%indexL
      itemL  => hecMAT%itemL
      indexU => hecMAT%indexU
      itemU  => hecMAT%itemU
      if (NContact.gt.0) then
        CAL => hecMAT%CAL
        CAU => hecMAT%CAU
        indexCL => hecMAT%indexCL
        itemCL => hecMAT%itemCL
        indexCU => hecMAT%indexCU
        itemCU => hecMAT%itemCU
      end if
    else
      allocate(COLORindex(0:NP), perm_tmp(NP), perm(NP), iperm(NP))
      call hecmw_matrix_ordering_RCM(N, hecMAT%indexL, hecMAT%itemL, &
           hecMAT%indexU, hecMAT%itemU, perm_tmp, iperm)
      write(*,*) 'DEBUG: RCM ordering done', hecmw_Wtime()-t0
      call hecmw_matrix_ordering_MC(N, hecMAT%indexL, hecMAT%itemL, &
           hecMAT%indexU, hecMAT%itemU, perm_tmp, &
           NCOLOR_IN, NColor, COLORindex, perm, iperm)
      write(*,*) 'DEBUG: MC ordering done', hecmw_Wtime()-t0
      deallocate(perm_tmp)

      call write_debug_info

      NPL = hecMAT%NPL
      NPU = hecMAT%NPU
      allocate(indexL(0:NP), indexU(0:NP), itemL(NPL), itemU(NPU))
      call hecmw_matrix_reorder_profile(N, perm, iperm, &
           hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
           indexL, indexU, itemL, itemU)
      write(*,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0

      call check_ordering

      allocate(D(9*NP), AL(9*NPL), AU(9*NPU))
      call hecmw_matrix_reorder_values(N, 3, perm, iperm, &
           hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
           hecMAT%AL, hecMAT%AU, hecMAT%D, &
           indexL, indexU, itemL, itemU, AL, AU, D)
      write(*,*) 'DEBUG: reordering values done', hecmw_Wtime()-t0

      if (NContact.gt.0) then
        NPCL = hecMAT%indexCL(NP)
        NPCU = hecMAT%indexCU(NP)
        allocate(indexCL(0:NP), indexCU(0:NP), itemCL(NPCL), itemCU(NPCU))
        call hecmw_matrix_reorder_profile(N, perm, iperm, &
             hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
             indexCL, indexCU, itemCL, itemCU)

        allocate(CD(9*NP), CAL(9*NPCL), CAU(9*NPCU))
        call hecmw_matrix_reorder_values(N, 3, perm, iperm, &
             hecMAT%indexCL, hecMAT%indexCU, hecMAT%itemCL, hecMAT%itemCU, &
             hecMAT%CAL, hecMAT%CAU, hecMAT%D, &
             indexCL, indexCU, itemCL, itemCU, CAL, CAU, CD)
        deallocate(CD)
      end if
    end if

    allocate(ALU(9*NP))
    ALU  = 0.d0

    do ii= 1, 9*N
      ALU(ii) = D(ii)
    enddo

    if (NContact.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = hecMAT%cmat%pair(k)%i
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
      ALU(9*ii-8)= ALUtmp(1,1)
      ALU(9*ii-7)= ALUtmp(1,2)
      ALU(9*ii-6)= ALUtmp(1,3)
      ALU(9*ii-5)= ALUtmp(2,1)
      ALU(9*ii-4)= ALUtmp(2,2)
      ALU(9*ii-3)= ALUtmp(2,3)
      ALU(9*ii-2)= ALUtmp(3,1)
      ALU(9*ii-1)= ALUtmp(3,2)
      ALU(9*ii  )= ALUtmp(3,3)
    enddo

    write(*,*) 'DEBUG: SSOR setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_SSOR_33_setup

  subroutine hecmw_precond_SSOR_33_apply(WW)
    implicit none
    real(kind=kreal), intent(inout), target :: WW(:)
    integer(kind=kint) :: ic, i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3
    real(kind=kreal), pointer :: ZP(:)

    ! added for turning >>>
    logical, save :: isFirst = .true.
    integer(kind=kint), parameter :: numOfBlockPerThread = 100
    integer(kind=kint), save :: numOfThread = 1, numOfBlock
    integer(kind=kint), save, allocatable :: icToBlockIndex(:)
    integer(kind=kint), save, allocatable :: blockIndexToColorIndex(:)
    integer(kind=kint), save :: sectorCacheSize0, sectorCacheSize1
    integer(kind=kint) :: blockIndex, elementCount, numOfElement, ii
    real(kind=kreal) :: numOfElementPerBlock
    integer(kind=kint) :: my_rank

    if (isFirst .eqv. .true.) then
      !$ numOfThread = omp_get_max_threads()
      numOfBlock = numOfThread * numOfBlockPerThread
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
      write(9000+my_rank,*) &
           '# numOfElementPerBlock =', numOfElementPerBlock
      write(9000+my_rank,*) &
           '# ic, blockIndex, colorIndex, elementCount'
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
            write(9000+my_rank,*) ic, blockIndex, &
                 blockIndexToColorIndex(blockIndex), elementCount
          endif
        enddo
        icToBlockIndex(ic) = blockIndex
      enddo
      numOfBlock = blockIndex

      ! calculate sector cache size
      sectorCacheSize1 = int((dble(N) * 3 * kreal / (4096 * 128)) + 0.999)
      if (sectorCacheSize1 > 6 ) sectorCacheSize1 = 6
      sectorCacheSize0 = 12 - sectorCacheSize1
      write(*,*) 'ZP size =', N * 3 * kreal, '[byte]  ', &
           'sectorCache0 =', sectorCacheSize0, '[way]  ', &
           'sectorCache1 =', sectorCacheSize1, '[way]'

      isFirst = .false.
    endif
    ! <<< added for turning

    if (numOfThread > 1) then
      allocate(ZP(N*3))
      call hecmw_matrix_reorder_vector(N, 3, perm, WW, ZP)
    else
      ZP => WW
    end if

!call start_collection("loopInPrecond33")

!xx!OCL CACHE_SECTOR_SIZE(sectorCacheSize0,sectorCacheSize1)
!xx!OCL CACHE_SUBSECTOR_ASSIGN(ZP)

!$omp parallel default(none) &
!$omp&shared(NColor,indexL,itemL,indexU,itemU,AL,AU,D,ALU,&
!$omp&       NContact,indexCL,itemCL,indexCU,itemCU,CAL,CAU,&
!$omp&       ZP,icToBlockIndex,blockIndexToColorIndex) &
!$omp&private(SW1,SW2,SW3,X1,X2,X3,ic,i,isL,ieL,isU,ieU,j,k,blockIndex)

    !C-- FORWARD
    do ic=1,NColor
!$omp do schedule (static, 1)
      do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
        do i = blockIndexToColorIndex(blockIndex-1)+1, &
             blockIndexToColorIndex(blockIndex)
          ! do i = startPos(threadNum, ic), endPos(threadNum, ic)
          SW1= ZP(3*i-2)
          SW2= ZP(3*i-1)
          SW3= ZP(3*i  )
          isL= indexL(i-1)+1
          ieL= indexL(i)
          do j= isL, ieL
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
              k= itemCL(j)
              X1= ZP(3*k-2)
              X2= ZP(3*k-1)
              X3= ZP(3*k  )
              SW1= SW1 - CAL(9*j-8)*X1 - CAL(9*j-7)*X2 - CAL(9*j-6)*X3
              SW2= SW2 - CAL(9*j-5)*X1 - CAL(9*j-4)*X2 - CAL(9*j-3)*X3
              SW3= SW3 - CAL(9*j-2)*X1 - CAL(9*j-1)*X2 - CAL(9*j  )*X3
            enddo ! j
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          ZP(3*i-2)= X1
          ZP(3*i-1)= X2
          ZP(3*i  )= X3
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
              k= itemCU(j)
              X1= ZP(3*k-2)
              X2= ZP(3*k-1)
              X3= ZP(3*k  )
              SW1= SW1 + CAU(9*j-8)*X1 + CAU(9*j-7)*X2 + CAU(9*j-6)*X3
              SW2= SW2 + CAU(9*j-5)*X1 + CAU(9*j-4)*X2 + CAU(9*j-3)*X3
              SW3= SW3 + CAU(9*j-2)*X1 + CAU(9*j-1)*X2 + CAU(9*j  )*X3
            enddo ! j
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          ZP(3*i-2)=  ZP(3*i-2) - X1
          ZP(3*i-1)=  ZP(3*i-1) - X2
          ZP(3*i  )=  ZP(3*i  ) - X3
        enddo ! i
      enddo ! blockIndex
!$omp end do
    enddo ! ic
!$omp end parallel

!xx!OCL END_CACHE_SUBSECTOR
!xx!OCL END_CACHE_SECTOR_SIZE

!call stop_collection("loopInPrecond33")

    if (numOfThread > 1) then
      call hecmw_matrix_reorder_back_vector(N, 3, perm, ZP, WW)
      deallocate(ZP)
    end if
  end subroutine hecmw_precond_SSOR_33_apply

  subroutine hecmw_precond_SSOR_33_clear(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: nthreads = 1
    !$ nthreads = omp_get_max_threads()
    if (associated(ALU)) deallocate(ALU)
    if (nthreads > 1) then
      if (associated(AL)) deallocate(AL)
      if (associated(AU)) deallocate(AU)
      if (associated(D)) deallocate(D)
      if (associated(indexL)) deallocate(indexL)
      if (associated(itemL)) deallocate(itemL)
      if (associated(indexU)) deallocate(indexU)
      if (associated(itemU)) deallocate(itemU)
      if (NContact.ne.0) then
        if (associated(CAL)) deallocate(CAL)
        if (associated(CAU)) deallocate(CAU)
        if (associated(indexCL)) deallocate(indexCL)
        if (associated(itemCL)) deallocate(itemCL)
        if (associated(indexCU)) deallocate(indexCU)
        if (associated(itemCU)) deallocate(itemCU)
      end if
    end if
    nullify(ALU)
    nullify(AU)
    nullify(AL)
    nullify(D)
    nullify(indexL)
    nullify(itemL)
    nullify(indexU)
    nullify(itemU)
    if (NContact.ne.0) then
      nullify(CAU)
      nullify(CAL)
      nullify(indexCL)
      nullify(itemCL)
      nullify(indexCU)
      nullify(itemCU)
      call hecmw_cmat_LU_free( hecMAT )
    endif
    NContact = 0
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
