!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_BILU_33
!C***
!C
module hecmw_precond_BILU_33
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_matrix_ordering_CM
  use m_hecmw_matrix_ordering_MC
  use hecmw_matrix_reorder
  use hecmw_matrix_contact
  !$ use omp_lib

  private
  integer(kind=kint), parameter :: DEBUG = 1

  public:: hecmw_precond_BILU_33_setup
  public:: hecmw_precond_BILU_33_apply
  public:: hecmw_precond_BILU_33_apply_opt1
  public:: hecmw_precond_BILU_33_apply_opt2
  public:: hecmw_precond_BILU_33_apply_opt3
  public:: hecmw_precond_BILU_33_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: Dlu0(:) => null()
  real(kind=kreal), pointer :: ALlu0(:) => null()
  real(kind=kreal), pointer :: AUlu0(:) => null()
  integer(kind=kint), pointer :: inumFI1L(:) => null()
  integer(kind=kint), pointer :: inumFI1U(:) => null()
  integer(kind=kint), pointer :: FI1L(:) => null()
  integer(kind=kint), pointer :: FI1U(:) => null()
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

  ! for tuning
  integer(kind=kint), parameter :: numOfBlockPerThread = 100
  integer(kind=kint), save :: numOfThread = 1, numOfBlock
  integer(kind=kint), save, allocatable :: icToBlockIndex(:)
  integer(kind=kint), save, allocatable :: blockIndexToColorIndex(:)

  integer(kind=kint) :: jmaxU, jmaxL
  real(kind=kreal), pointer :: SW1(:) => null()
  real(kind=kreal), pointer :: SW2(:) => null()
  real(kind=kreal), pointer :: SW3(:) => null()

  logical, save :: isFirst = .true.
  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_BILU_33_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP, NPU, NPL, NPCL, NPCU
    real   (kind=kreal), allocatable :: CD(:)
    integer(kind=kint ) :: PRECOND
    real   (kind=kreal) :: SIGMA, SIGMA_DIAG

    real(kind=kreal), pointer :: D(:)
    real(kind=kreal), pointer :: AL(:)
    real(kind=kreal), pointer :: AU(:)
    integer(kind=kint ) :: ii, i, j, k

    integer(kind=kint ), pointer :: INL(:), INU(:)
    integer(kind=kint ), pointer :: IAL(:)
    integer(kind=kint ), pointer :: IAU(:)
    integer(kind=kint ) :: NCOLOR_IN
    integer(kind=kint ), allocatable :: perm_tmp(:)
    real   (kind=kreal) :: t0
    integer(kind=kint) :: isL, ieL, isU, ieU

    if (DEBUG >= 1) then
      t0 = hecmw_Wtime()
      write(0,*) 'DEBUG: BILU setup start', hecmw_Wtime()-t0
    endif

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 1) then ! need symbolic and numerical setup
        call hecmw_precond_BILU_33_clear()
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_BILU_33_clear() ! TEMPORARY
      else
        return
      endif
    endif

    N = hecMAT%N
    NP = hecMAT%NP
    NPL = hecMAT%NPL
    NPU = hecMAT%NPU
    D => hecMAT%D
    AL => hecMAT%AL
    AU => hecMAT%AU
    INL => hecMAT%indexL
    INU => hecMAT%indexU
    IAL => hecMAT%itemL
    IAU => hecMAT%itemU
    PRECOND = hecmw_mat_get_precond(hecMAT)
    SIGMA = hecmw_mat_get_sigma(hecMAT)
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    if (PRECOND.eq.10) call FORM_ILU0_33 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.11) call FORM_ILU1_33 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.12) call FORM_ILU2_33 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)

    NCOLOR_IN = hecmw_mat_get_ncolor_in(hecMAT)
    NContact = hecMAT%cmat%n_val

    if (NContact.gt.0) then
      call hecmw_cmat_LU( hecMAT )
    endif

    allocate(COLORindex(0:N), perm_tmp(N), perm(N), iperm(N))
    call hecmw_matrix_ordering_RCM(N, hecMAT%indexL, hecMAT%itemL, &
      hecMAT%indexU, hecMAT%itemU, perm_tmp, iperm)
    if (DEBUG >= 1) write(0,*) 'DEBUG: RCM ordering done', hecmw_Wtime()-t0
    call hecmw_matrix_ordering_MC(N, hecMAT%indexL, hecMAT%itemL, &
      hecMAT%indexU, hecMAT%itemU, perm_tmp, &
      NCOLOR_IN, NColor, COLORindex, perm, iperm)
    if (DEBUG >= 1) write(0,*) 'DEBUG: MC ordering done', hecmw_Wtime()-t0
    deallocate(perm_tmp)
    NPL = hecMAT%indexL(N)
    NPU = hecMAT%indexU(N)
    allocate(indexL(0:N), indexU(0:N), itemL(NPL), itemU(NPU))
    call hecmw_matrix_reorder_profile(N, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      indexL, indexU, itemL, itemU)
    if (DEBUG >= 1) write(0,*) 'DEBUG: reordering profile done', hecmw_Wtime()-t0

    !call check_ordering

    allocate(D(9*N), AL(9*NPL), AU(9*NPU))
    call hecmw_matrix_reorder_values(N, 3, perm, iperm, &
      hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU, &
      hecMAT%AL, hecMAT%AU, hecMAT%D, &
      indexL, indexU, itemL, itemU, AL, AU, D)
    if (DEBUG >= 1) write(0,*) 'DEBUG: reordering values done', hecmw_Wtime()-t0

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

    jmaxU=0
    jmaxL=0
    do i=1, N
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      if(jmaxL .lt. ieL-isL) then
        jmaxL=ieL-isL+1
      endif
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      if(jmaxU .lt. ieU-isU) then
        jmaxU=ieU-isU+1
      endif
    enddo
    allocate(SW1(N))
    allocate(SW2(N))
    allocate(SW3(N))


    isFirst = .true.
    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

    if (DEBUG >= 1) write(0,*) "DEBUG: NCOLOR = ",NColor
    if (DEBUG >= 1) write(0,*) "DEBUG: NContact = ",NContact
    if (DEBUG >= 1) write(0,*) 'DEBUG: BILU setup done', hecmw_Wtime()-t0

  end subroutine hecmw_precond_BILU_33_setup

  subroutine setup_tuning_parameters
    implicit none
    integer(kind=kint) :: blockIndex, elementCount, numOfElement, ii
    real(kind=kreal) :: numOfElementPerBlock
    integer(kind=kint) :: my_rank
    integer(kind=kint) :: ic, i
    if (DEBUG >= 1) write(*,*) 'DEBUG: setting up tuning parameters for BILU'
#ifdef _OPENMP
    numOfThread = omp_get_max_threads()
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
  end subroutine setup_tuning_parameters

  subroutine hecmw_precond_BILU_33_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3
    !C
    !C-- FORWARD

    do i= 1, N
      SW1= WW(3*i-2)
      SW2= WW(3*i-1)
      SW3= WW(3*i  )
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      do j= isL, ieL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      enddo

      X1= SW1
      X2= SW2
      X3= SW3
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)= X1
      WW(3*i-1)= X2
      WW(3*i  )= X3
    enddo

    !C
    !C-- BACKWARD

    do i= N, 1, -1
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      SW1= 0.d0
      SW2= 0.d0
      SW3= 0.d0
      do j= ieU, isU, -1
        k= FI1U(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
        SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
        SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      enddo
      X1= SW1
      X2= SW2
      X3= SW3
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)=  WW(3*i-2) - X1
      WW(3*i-1)=  WW(3*i-1) - X2
      WW(3*i  )=  WW(3*i  ) - X3
    enddo
  end subroutine hecmw_precond_BILU_33_apply

  subroutine hecmw_precond_BILU_33_apply_opt1(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    integer(kind=kint) :: j_offset
    real(kind=kreal) ::  X1, X2, X3
    !C
    !C-- FORWARD

    do i= 1, N
      SW1(i)= WW(3*i-2)
      SW2(i)= WW(3*i-1)
      SW3(i)= WW(3*i  )
    enddo
    do j_offset=0, jmaxL
!NEC$ ivdep
      do i= 1, N
        isL= inumFI1L(i-1)+1
        ieL= inumFI1L(i)
        j=isL+j_offset
        if(j .le. ieL) then
          k= FI1L(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1(i)= SW1(i) - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
          SW2(i)= SW2(i) - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
          SW3(i)= SW3(i) - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        end if
      enddo
      X1= SW1(i)
      X2= SW2(i)
      X3= SW3(i)
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)= X1
      WW(3*i-1)= X2
      WW(3*i  )= X3
    enddo


    !C
    !C-- BACKWARD

    do i= N, 1, -1
      SW1(i)= 0.d0
      SW2(i)= 0.d0
      SW3(i)= 0.d0
    enddo
    do j_offset=0, jmaxU
!NEC$ ivdep
      do i= N, 1, -1
        isU= inumFI1U(i-1) + 1
        ieU= inumFI1U(i)
        j=ieU-j_offset
        if(j .ge. isU) then
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1(i)= SW1(i) + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2(i)= SW2(i) + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3(i)= SW3(i) + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
        end if
      enddo
      X1= SW1(i)
      X2= SW2(i)
      X3= SW3(i)
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)=  WW(3*i-2) - X1
      WW(3*i-1)=  WW(3*i-1) - X2
      WW(3*i  )=  WW(3*i  ) - X3
    enddo
  end subroutine hecmw_precond_BILU_33_apply_opt1

  subroutine hecmw_precond_BILU_33_apply_opt2(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k, j_offset
    integer(kind=kint) :: ic
    real(kind=kreal) ::  X1, X2, X3

    integer(kind=kint) :: blockIndex
    !C
    !C-- FORWARD

    if (isFirst) then
      call setup_tuning_parameters
      isFirst = .false.
    endif

    do i= 1, N
      SW1(i)= WW(3*i-2)
      SW2(i)= WW(3*i-1)
      SW3(i)= WW(3*i  )
    enddo
    do j_offset=0, jmaxL
      do ic=1,NColor
        do blockIndex = icToBlockIndex(ic-1)+1, icToBlockIndex(ic)
!NEC$ ivdep
          do i = blockIndexToColorIndex(blockIndex-1)+1, &
              blockIndexToColorIndex(blockIndex)
          isL= inumFI1L(i-1)+1
          ieL= inumFI1L(i)
          j=isL+j_offset
          if(j .le. ieL) then
            k= FI1L(j)
            X1= WW(3*k-2)
            X2= WW(3*k-1)
            X3= WW(3*k  )
            SW1(i)= SW1(i) - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
            SW2(i)= SW2(i) - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
            SW3(i)= SW3(i) - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
          end if
          X1= SW1(i)
          X2= SW2(i)
          X3= SW3(i)
          X2= X2 - Dlu0(9*i-5)*X1
          X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
          X3= Dlu0(9*i  )*  X3
          X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
          X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
          WW(3*i-2)= X1
          WW(3*i-1)= X2
          WW(3*i  )= X3
        enddo
      enddo
    enddo
    enddo

    !C
    !C-- BACKWARD

    do i= N, 1, -1
      SW1(i)= 0.d0
      SW2(i)= 0.d0
      SW3(i)= 0.d0
    enddo
    do j_offset=0, jmaxU
      do ic=NColor, 1, -1
        do blockIndex = icToBlockIndex(ic), icToBlockIndex(ic-1)+1, -1
!NEC$ ivdep
          do i = blockIndexToColorIndex(blockIndex), &
              blockIndexToColorIndex(blockIndex-1)+1, -1
            isU= inumFI1U(i-1) + 1
            ieU= inumFI1U(i)
            SW1(i)= 0.d0
            SW2(i)= 0.d0
            SW3(i)= 0.d0
            j=ieU-j_offset
            if(j .ge. isU) then
              k= FI1U(j)
              X1= WW(3*k-2)
              X2= WW(3*k-1)
              X3= WW(3*k  )
              SW1(i)= SW1(i) + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
              SW2(i)= SW2(i) + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
              SW3(i)= SW3(i) + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
            end if
          enddo
          X1= SW1(i)
          X2= SW2(i)
          X3= SW3(i)
          X2= X2 - Dlu0(9*i-5)*X1
          X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
          X3= Dlu0(9*i  )*  X3
          X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
          X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
          WW(3*i-2)=  WW(3*i-2) - X1
          WW(3*i-1)=  WW(3*i-1) - X2
          WW(3*i  )=  WW(3*i  ) - X3
        enddo
      enddo
    enddo
  end subroutine hecmw_precond_BILU_33_apply_opt2

  subroutine hecmw_precond_BILU_33_apply_opt3(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3
    !C
    !C-- FORWARD

!NEC$ ivdep
    do i= 1, N
      SW1= WW(3*i-2)
      SW2= WW(3*i-1)
      SW3= WW(3*i  )
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      if(ieL - iSL .eq. 1) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else if(ieL - iSL .eq.2) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else if(ieL - iSL .eq.3) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else if(ieL - iSL .eq. 4) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else if(ieL - iSL .eq. 5) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else if(ieL - iSL .eq. 6) then
        j=isL
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        j=j+1
        k= FI1L(j)
        X1= WW(3*k-2)
        X2= WW(3*k-1)
        X3= WW(3*k  )
        SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
        SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
        SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
      else
!NEC$ novector
        do j= isL, ieL
          k= FI1L(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 - ALlu0(9*j-8)*X1-ALlu0(9*j-7)*X2-ALlu0(9*j-6)*X3
          SW2= SW2 - ALlu0(9*j-5)*X1-ALlu0(9*j-4)*X2-ALlu0(9*j-3)*X3
          SW3= SW3 - ALlu0(9*j-2)*X1-ALlu0(9*j-1)*X2-ALlu0(9*j  )*X3
        enddo
      endif
      X1= SW1
      X2= SW2
      X3= SW3
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)= X1
      WW(3*i-1)= X2
      WW(3*i  )= X3
    enddo

    !C
    !C-- BACKWARD

!NEC$ ivdep
    do i= N, 1, -1
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      SW1= 0.d0
      SW2= 0.d0
      SW3= 0.d0
      if (ieU-isU .eq. 1) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 2) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 3) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 4) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 5) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 6) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 7) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 8) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 9) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 10) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 11) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 12) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 13) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 14) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else if(ieU-isU .eq. 15) then
          j=ieU
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
          j=j-1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
      else
!NEC$ novector
        do j= ieU, isU, -1
          k= FI1U(j)
          X1= WW(3*k-2)
          X2= WW(3*k-1)
          X3= WW(3*k  )
          SW1= SW1 + AUlu0(9*j-8)*X1+AUlu0(9*j-7)*X2+AUlu0(9*j-6)*X3
          SW2= SW2 + AUlu0(9*j-5)*X1+AUlu0(9*j-4)*X2+AUlu0(9*j-3)*X3
          SW3= SW3 + AUlu0(9*j-2)*X1+AUlu0(9*j-1)*X2+AUlu0(9*j  )*X3
        enddo
      end if
      X1= SW1
      X2= SW2
      X3= SW3
      X2= X2 - Dlu0(9*i-5)*X1
      X3= X3 - Dlu0(9*i-2)*X1 - Dlu0(9*i-1)*X2
      X3= Dlu0(9*i  )*  X3
      X2= Dlu0(9*i-4)*( X2 - Dlu0(9*i-3)*X3 )
      X1= Dlu0(9*i-8)*( X1 - Dlu0(9*i-6)*X3 - Dlu0(9*i-7)*X2)
      WW(3*i-2)=  WW(3*i-2) - X1
      WW(3*i-1)=  WW(3*i-1) - X2
      WW(3*i  )=  WW(3*i  ) - X3
    enddo
  end subroutine hecmw_precond_BILU_33_apply_opt3

  subroutine hecmw_precond_BILU_33_clear()
    implicit none
    if (associated(Dlu0)) deallocate(Dlu0)
    if (associated(ALlu0)) deallocate(ALlu0)
    if (associated(AUlu0)) deallocate(AUlu0)
    if (associated(inumFI1L)) deallocate(inumFI1L)
    if (associated(inumFI1U)) deallocate(inumFI1U)
    if (associated(FI1L)) deallocate(FI1L)
    if (associated(FI1U)) deallocate(FI1U)
    nullify(Dlu0)
    nullify(ALlu0)
    nullify(AUlu0)
    nullify(inumFI1L)
    nullify(inumFI1U)
    nullify(FI1L)
    nullify(FI1U)
    INITIALIZED = .false.
  end subroutine hecmw_precond_BILU_33_clear

  !C
  !C***
  !C*** FORM_ILU0_33
  !C***
  !C
  !C    form ILU(0) matrix
  !C
  subroutine FORM_ILU0_33                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(9*NPL), intent(in):: AL
    real(kind=kreal), dimension(9*NPU), intent(in):: AU
    real(kind=kreal), dimension(9*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    real (kind=kreal),  dimension(3,3) :: RHS_Aij, DkINV, Aik, Akj
    integer(kind=kint) :: i,jj,ij0,kk
    integer(kind=kint) :: j,k
    allocate (IW1(NP) , IW2(NP))
    allocate(Dlu0(9*NP), ALlu0(9*NPL), AUlu0(9*NPU))
    allocate(inumFI1L(0:NP), inumFI1U(0:NP), FI1L(NPL), FI1U(NPU))

    do i=1,9*NP
      Dlu0(i) = D(i)
    end do
    do i=1,9*NPL
      ALlu0(i) = AL(i)
    end do
    do i=1,9*NPU
      AUlu0(i) = AU(i)
    end do
    do i=0,NP
      inumFI1L(i) = INL(i)
      inumFI1U(i) = INU(i)
    end do
    do i=1,NPL
      FI1L(i) = IAL(i)
    end do
    do i=1,NPU
      FI1U(i) = IAU(i)
    end do

    !C
    !C +----------------------+
    !C | ILU(0) factorization |
    !C +----------------------+
    !C===
    do i=1,NP
      Dlu0(9*i-8)=Dlu0(9*i-8)*SIGMA_DIAG
      Dlu0(9*i-4)=Dlu0(9*i-4)*SIGMA_DIAG
      Dlu0(9*i  )=Dlu0(9*i  )*SIGMA_DIAG
    enddo

    i = 1
    call ILU1a33 (DkINV, &
      Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
      Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
      Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
    Dlu0(9*i-8)= DkINV(1,1)
    Dlu0(9*i-7)= DkINV(1,2)
    Dlu0(9*i-6)= DkINV(1,3)
    Dlu0(9*i-5)= DkINV(2,1)
    Dlu0(9*i-4)= DkINV(2,2)
    Dlu0(9*i-3)= DkINV(2,3)
    Dlu0(9*i-2)= DkINV(3,1)
    Dlu0(9*i-1)= DkINV(3,2)
    Dlu0(9*i  )= DkINV(3,3)

    do i= 2, NP
      IW1= 0
      IW2= 0

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        IW1(FI1L(k))= k
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        IW2(FI1U(k))= k
      enddo

      do kk= INL(i-1)+1, INL(i)
        k= IAL(kk)

        DkINV(1,1)= Dlu0(9*k-8)
        DkINV(1,2)= Dlu0(9*k-7)
        DkINV(1,3)= Dlu0(9*k-6)
        DkINV(2,1)= Dlu0(9*k-5)
        DkINV(2,2)= Dlu0(9*k-4)
        DkINV(2,3)= Dlu0(9*k-3)
        DkINV(3,1)= Dlu0(9*k-2)
        DkINV(3,2)= Dlu0(9*k-1)
        DkINV(3,3)= Dlu0(9*k  )

        Aik(1,1)= ALlu0(9*kk-8)
        Aik(1,2)= ALlu0(9*kk-7)
        Aik(1,3)= ALlu0(9*kk-6)
        Aik(2,1)= ALlu0(9*kk-5)
        Aik(2,2)= ALlu0(9*kk-4)
        Aik(2,3)= ALlu0(9*kk-3)
        Aik(3,1)= ALlu0(9*kk-2)
        Aik(3,2)= ALlu0(9*kk-1)
        Aik(3,3)= ALlu0(9*kk  )

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          if (IW1(j).eq.0.and.IW2(j).eq.0) cycle

          Akj(1,1)= AUlu0(9*jj-8)
          Akj(1,2)= AUlu0(9*jj-7)
          Akj(1,3)= AUlu0(9*jj-6)
          Akj(2,1)= AUlu0(9*jj-5)
          Akj(2,2)= AUlu0(9*jj-4)
          Akj(2,3)= AUlu0(9*jj-3)
          Akj(3,1)= AUlu0(9*jj-2)
          Akj(3,2)= AUlu0(9*jj-1)
          Akj(3,3)= AUlu0(9*jj  )

          call ILU1b33 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(9*i-8)= Dlu0(9*i-8) - RHS_Aij(1,1)
            Dlu0(9*i-7)= Dlu0(9*i-7) - RHS_Aij(1,2)
            Dlu0(9*i-6)= Dlu0(9*i-6) - RHS_Aij(1,3)
            Dlu0(9*i-5)= Dlu0(9*i-5) - RHS_Aij(2,1)
            Dlu0(9*i-4)= Dlu0(9*i-4) - RHS_Aij(2,2)
            Dlu0(9*i-3)= Dlu0(9*i-3) - RHS_Aij(2,3)
            Dlu0(9*i-2)= Dlu0(9*i-2) - RHS_Aij(3,1)
            Dlu0(9*i-1)= Dlu0(9*i-1) - RHS_Aij(3,2)
            Dlu0(9*i  )= Dlu0(9*i  ) - RHS_Aij(3,3)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(9*ij0-8)= ALlu0(9*ij0-8) - RHS_Aij(1,1)
            ALlu0(9*ij0-7)= ALlu0(9*ij0-7) - RHS_Aij(1,2)
            ALlu0(9*ij0-6)= ALlu0(9*ij0-6) - RHS_Aij(1,3)
            ALlu0(9*ij0-5)= ALlu0(9*ij0-5) - RHS_Aij(2,1)
            ALlu0(9*ij0-4)= ALlu0(9*ij0-4) - RHS_Aij(2,2)
            ALlu0(9*ij0-3)= ALlu0(9*ij0-3) - RHS_Aij(2,3)
            ALlu0(9*ij0-2)= ALlu0(9*ij0-2) - RHS_Aij(3,1)
            ALlu0(9*ij0-1)= ALlu0(9*ij0-1) - RHS_Aij(3,2)
            ALlu0(9*ij0  )= ALlu0(9*ij0  ) - RHS_Aij(3,3)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(9*ij0-8)= AUlu0(9*ij0-8) - RHS_Aij(1,1)
            AUlu0(9*ij0-7)= AUlu0(9*ij0-7) - RHS_Aij(1,2)
            AUlu0(9*ij0-6)= AUlu0(9*ij0-6) - RHS_Aij(1,3)
            AUlu0(9*ij0-5)= AUlu0(9*ij0-5) - RHS_Aij(2,1)
            AUlu0(9*ij0-4)= AUlu0(9*ij0-4) - RHS_Aij(2,2)
            AUlu0(9*ij0-3)= AUlu0(9*ij0-3) - RHS_Aij(2,3)
            AUlu0(9*ij0-2)= AUlu0(9*ij0-2) - RHS_Aij(3,1)
            AUlu0(9*ij0-1)= AUlu0(9*ij0-1) - RHS_Aij(3,2)
            AUlu0(9*ij0  )= AUlu0(9*ij0  ) - RHS_Aij(3,3)
          endif

        enddo
      enddo

      call ILU1a33 (DkINV, &
        Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
        Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
        Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
      Dlu0(9*i-8)= DkINV(1,1)
      Dlu0(9*i-7)= DkINV(1,2)
      Dlu0(9*i-6)= DkINV(1,3)
      Dlu0(9*i-5)= DkINV(2,1)
      Dlu0(9*i-4)= DkINV(2,2)
      Dlu0(9*i-3)= DkINV(2,3)
      Dlu0(9*i-2)= DkINV(3,1)
      Dlu0(9*i-1)= DkINV(3,2)
      Dlu0(9*i  )= DkINV(3,3)
    enddo

    deallocate (IW1, IW2)
  end subroutine FORM_ILU0_33

  !C
  !C***
  !C*** FORM_ILU1_33
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU1_33                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(9*NPL), intent(in):: AL
    real(kind=kreal), dimension(9*NPU), intent(in):: AU
    real(kind=kreal), dimension(9*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    integer(kind=kint), dimension(:), allocatable :: IWsL, IWsU
    real (kind=kreal),  dimension(3,3) :: RHS_Aij, DkINV, Aik, Akj
    integer(kind=kint) :: NPLf1,NPUf1
    integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
    integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3
    integer(kind=kint) :: j,k,iSL,iSU
    !C
    !C +--------------+
    !C | find fill-in |
    !C +--------------+
    !C===

    !C
    !C-- count fill-in
    allocate (IW1(NP) , IW2(NP))
    allocate (inumFI1L(0:NP), inumFI1U(0:NP))

    inumFI1L= 0
    inumFI1U= 0

    NPLf1= 0
    NPUf1= 0
    do i= 2, NP
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 1
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            inumFI1L(i)= inumFI1L(i)+1
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            inumFI1U(i)= inumFI1U(i)+1
            IW1(jj)= 1
          endif
        enddo
      enddo
      NPLf1= NPLf1 + inumFI1L(i)
      NPUf1= NPUf1 + inumFI1U(i)
    enddo

    !C
    !C-- specify fill-in
    allocate (IWsL(0:NP), IWsU(0:NP))
    allocate (FI1L (NPL+NPLf1), FI1U (NPU+NPUf1))
    allocate (ALlu0(9*(NPL+NPLf1)), AUlu0(9*(NPU+NPUf1)))

    FI1L= 0
    FI1U= 0

    IWsL= 0
    IWsU= 0
    do i= 1, NP
      IWsL(i)= INL(i)-INL(i-1) + inumFI1L(i) + IWsL(i-1)
      IWsU(i)= INU(i)-INU(i-1) + inumFI1U(i) + IWsU(i-1)
    enddo

    do i= 2, NP
      icouL= 0
      icouU= 0
      inumFI1L(i)= inumFI1L(i-1) + inumFI1L(i)
      inumFI1U(i)= inumFI1U(i-1) + inumFI1U(i)
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 1
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            icouL           = icouL + 1
            FI1L(icouL+IWsL(i-1)+INL(i)-INL(i-1))= jj
            IW1(jj)          = 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            icouU           = icouU + 1
            FI1U(icouU+IWsU(i-1)+INU(i)-INU(i-1))= jj
            IW1(jj)          = 1
          endif
        enddo
      enddo
    enddo
    !C===

    !C
    !C +-------------------------------------------------+
    !C | SORT and RECONSTRUCT matrix considering fill-in |
    !C +-------------------------------------------------+
    !C===
    ALlu0= 0.d0
    AUlu0= 0.d0
    iSL  = 0
    iSU  = 0
    do i= 1, NP
      icouL1=      INL(i) -      INL(i-1)
      icouL2= inumFI1L(i) - inumFI1L(i-1)
      icouL3= icouL1 + icouL2
      icouU1=      INU(i) -      INU(i-1)
      icouU2= inumFI1U(i) - inumFI1U(i-1)
      icouU3= icouU1 + icouU2
      !C
      !C-- LOWER part
      icou0= 0
      do k= INL(i-1)+1, INL(i)
        icou0 = icou0 + 1
        IW1(icou0)= IAL(k)
      enddo

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        icou0 = icou0 + 1
        IW1(icou0)= FI1L(icou0+IWsL(i-1))
      enddo

      do k= 1, icouL3
        IW2(k)= k
      enddo
      call fill_in_S33_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 9*( k+isL)
          kk2= 9*(ik+INL(i-1))
          ALlu0(kk1-8)= AL(kk2-8)
          ALlu0(kk1-7)= AL(kk2-7)
          ALlu0(kk1-6)= AL(kk2-6)
          ALlu0(kk1-5)= AL(kk2-5)
          ALlu0(kk1-4)= AL(kk2-4)
          ALlu0(kk1-3)= AL(kk2-3)
          ALlu0(kk1-2)= AL(kk2-2)
          ALlu0(kk1-1)= AL(kk2-1)
          ALlu0(kk1  )= AL(kk2  )
        endif
      enddo
      !C
      !C-- UPPER part
      icou0= 0
      do k= INU(i-1)+1, INU(i)
        icou0 = icou0 + 1
        IW1(icou0)= IAU(k)
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        icou0 = icou0 + 1
        IW1(icou0)= FI1U(icou0+IWsU(i-1))
      enddo

      do k= 1, icouU3
        IW2(k)= k
      enddo
      call fill_in_S33_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 9*( k+isU)
          kk2= 9*(ik+INU(i-1))
          AUlu0(kk1-8)= AU(kk2-8)
          AUlu0(kk1-7)= AU(kk2-7)
          AUlu0(kk1-6)= AU(kk2-6)
          AUlu0(kk1-5)= AU(kk2-5)
          AUlu0(kk1-4)= AU(kk2-4)
          AUlu0(kk1-3)= AU(kk2-3)
          AUlu0(kk1-2)= AU(kk2-2)
          AUlu0(kk1-1)= AU(kk2-1)
          AUlu0(kk1  )= AU(kk2  )
        endif
      enddo

      iSL= iSL + icouL3
      iSU= iSU + icouU3
    enddo

    !C===
    do i= 1, NP
      inumFI1L(i)= IWsL(i)
      inumFI1U(i)= IWsU(i)
    enddo
    deallocate (IWsL, IWsU)

    !C
    !C +----------------------+
    !C | ILU(1) factorization |
    !C +----------------------+
    !C===
    allocate (Dlu0(9*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(9*i-8)=Dlu0(9*i-8)*SIGMA_DIAG
      Dlu0(9*i-4)=Dlu0(9*i-4)*SIGMA_DIAG
      Dlu0(9*i  )=Dlu0(9*i  )*SIGMA_DIAG
    enddo

    i = 1
    call ILU1a33 (DkINV, &
      Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
      Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
      Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
    Dlu0(9*i-8)= DkINV(1,1)
    Dlu0(9*i-7)= DkINV(1,2)
    Dlu0(9*i-6)= DkINV(1,3)
    Dlu0(9*i-5)= DkINV(2,1)
    Dlu0(9*i-4)= DkINV(2,2)
    Dlu0(9*i-3)= DkINV(2,3)
    Dlu0(9*i-2)= DkINV(3,1)
    Dlu0(9*i-1)= DkINV(3,2)
    Dlu0(9*i  )= DkINV(3,3)

    do i= 2, NP
      IW1= 0
      IW2= 0

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        IW1(FI1L(k))= k
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        IW2(FI1U(k))= k
      enddo

      do kk= INL(i-1)+1, INL(i)
        k= IAL(kk)

        DkINV(1,1)= Dlu0(9*k-8)
        DkINV(1,2)= Dlu0(9*k-7)
        DkINV(1,3)= Dlu0(9*k-6)
        DkINV(2,1)= Dlu0(9*k-5)
        DkINV(2,2)= Dlu0(9*k-4)
        DkINV(2,3)= Dlu0(9*k-3)
        DkINV(3,1)= Dlu0(9*k-2)
        DkINV(3,2)= Dlu0(9*k-1)
        DkINV(3,3)= Dlu0(9*k  )

        do kk1= inumFI1L(i-1)+1, inumFI1L(i)
          if (k.eq.FI1L(kk1)) then
            Aik(1,1)= ALlu0(9*kk1-8)
            Aik(1,2)= ALlu0(9*kk1-7)
            Aik(1,3)= ALlu0(9*kk1-6)
            Aik(2,1)= ALlu0(9*kk1-5)
            Aik(2,2)= ALlu0(9*kk1-4)
            Aik(2,3)= ALlu0(9*kk1-3)
            Aik(3,1)= ALlu0(9*kk1-2)
            Aik(3,2)= ALlu0(9*kk1-1)
            Aik(3,3)= ALlu0(9*kk1  )
            exit
          endif
        enddo

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          do jj1= inumFI1U(k-1)+1, inumFI1U(k)
            if (j.eq.FI1U(jj1)) then
              Akj(1,1)= AUlu0(9*jj1-8)
              Akj(1,2)= AUlu0(9*jj1-7)
              Akj(1,3)= AUlu0(9*jj1-6)
              Akj(2,1)= AUlu0(9*jj1-5)
              Akj(2,2)= AUlu0(9*jj1-4)
              Akj(2,3)= AUlu0(9*jj1-3)
              Akj(3,1)= AUlu0(9*jj1-2)
              Akj(3,2)= AUlu0(9*jj1-1)
              Akj(3,3)= AUlu0(9*jj1  )
              exit
            endif
          enddo

          call ILU1b33 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(9*i-8)= Dlu0(9*i-8) - RHS_Aij(1,1)
            Dlu0(9*i-7)= Dlu0(9*i-7) - RHS_Aij(1,2)
            Dlu0(9*i-6)= Dlu0(9*i-6) - RHS_Aij(1,3)
            Dlu0(9*i-5)= Dlu0(9*i-5) - RHS_Aij(2,1)
            Dlu0(9*i-4)= Dlu0(9*i-4) - RHS_Aij(2,2)
            Dlu0(9*i-3)= Dlu0(9*i-3) - RHS_Aij(2,3)
            Dlu0(9*i-2)= Dlu0(9*i-2) - RHS_Aij(3,1)
            Dlu0(9*i-1)= Dlu0(9*i-1) - RHS_Aij(3,2)
            Dlu0(9*i  )= Dlu0(9*i  ) - RHS_Aij(3,3)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(9*ij0-8)= ALlu0(9*ij0-8) - RHS_Aij(1,1)
            ALlu0(9*ij0-7)= ALlu0(9*ij0-7) - RHS_Aij(1,2)
            ALlu0(9*ij0-6)= ALlu0(9*ij0-6) - RHS_Aij(1,3)
            ALlu0(9*ij0-5)= ALlu0(9*ij0-5) - RHS_Aij(2,1)
            ALlu0(9*ij0-4)= ALlu0(9*ij0-4) - RHS_Aij(2,2)
            ALlu0(9*ij0-3)= ALlu0(9*ij0-3) - RHS_Aij(2,3)
            ALlu0(9*ij0-2)= ALlu0(9*ij0-2) - RHS_Aij(3,1)
            ALlu0(9*ij0-1)= ALlu0(9*ij0-1) - RHS_Aij(3,2)
            ALlu0(9*ij0  )= ALlu0(9*ij0  ) - RHS_Aij(3,3)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(9*ij0-8)= AUlu0(9*ij0-8) - RHS_Aij(1,1)
            AUlu0(9*ij0-7)= AUlu0(9*ij0-7) - RHS_Aij(1,2)
            AUlu0(9*ij0-6)= AUlu0(9*ij0-6) - RHS_Aij(1,3)
            AUlu0(9*ij0-5)= AUlu0(9*ij0-5) - RHS_Aij(2,1)
            AUlu0(9*ij0-4)= AUlu0(9*ij0-4) - RHS_Aij(2,2)
            AUlu0(9*ij0-3)= AUlu0(9*ij0-3) - RHS_Aij(2,3)
            AUlu0(9*ij0-2)= AUlu0(9*ij0-2) - RHS_Aij(3,1)
            AUlu0(9*ij0-1)= AUlu0(9*ij0-1) - RHS_Aij(3,2)
            AUlu0(9*ij0  )= AUlu0(9*ij0  ) - RHS_Aij(3,3)
          endif

        enddo
      enddo

      call ILU1a33 (DkINV, &
        Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
        Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
        Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
      Dlu0(9*i-8)= DkINV(1,1)
      Dlu0(9*i-7)= DkINV(1,2)
      Dlu0(9*i-6)= DkINV(1,3)
      Dlu0(9*i-5)= DkINV(2,1)
      Dlu0(9*i-4)= DkINV(2,2)
      Dlu0(9*i-3)= DkINV(2,3)
      Dlu0(9*i-2)= DkINV(3,1)
      Dlu0(9*i-1)= DkINV(3,2)
      Dlu0(9*i  )= DkINV(3,3)
    enddo

    deallocate (IW1, IW2)
    !C===
  end subroutine FORM_ILU1_33

  !C
  !C***
  !C*** FORM_ILU2_33
  !C***
  !C
  !C    form ILU(2) matrix
  !C
  subroutine FORM_ILU2_33 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(9*NPL), intent(in):: AL
    real(kind=kreal), dimension(9*NPU), intent(in):: AU
    real(kind=kreal), dimension(9*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable:: IW1 , IW2
    integer(kind=kint), dimension(:), allocatable:: IWsL, IWsU
    integer(kind=kint), dimension(:), allocatable:: iconFI1L, iconFI1U
    integer(kind=kint), dimension(:), allocatable:: inumFI2L, inumFI2U
    integer(kind=kint), dimension(:), allocatable::     FI2L,     FI2U
    real (kind=kreal), dimension(3,3) :: RHS_Aij, DkINV, Aik, Akj
    integer(kind=kint) :: NPLf1,NPLf2,NPUf1,NPUf2,iAS,iconIK,iconKJ
    integer(kind=kint) :: i,jj,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
    integer(kind=kint) :: icou,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3
    integer(kind=kint) :: j,k,iSL,iSU

    !C
    !C +------------------+
    !C | find fill-in (1) |
    !C +------------------+
    !C===

    !C
    !C-- count fill-in
    allocate (IW1(NP) , IW2(NP))
    allocate (inumFI2L(0:NP), inumFI2U(0:NP))

    inumFI2L= 0
    inumFI2U= 0

    NPLf1= 0
    NPUf1= 0
    do i= 2, NP
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 1
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            inumFI2L(i)= inumFI2L(i)+1
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            inumFI2U(i)= inumFI2U(i)+1
            IW1(jj)= 1
          endif
        enddo
      enddo
      NPLf1= NPLf1 + inumFI2L(i)
      NPUf1= NPUf1 + inumFI2U(i)
    enddo

    !C
    !C-- specify fill-in
    allocate (IWsL(0:NP), IWsU(0:NP))
    allocate (FI2L (NPLf1), FI2U (NPUf1))

    FI2L= 0
    FI2U= 0

    do i= 2, NP
      icouL= 0
      icouU= 0
      inumFI2L(i)= inumFI2L(i-1) + inumFI2L(i)
      inumFI2U(i)= inumFI2U(i-1) + inumFI2U(i)
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 1
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            icouL = icouL + 1
            FI2L(icouL+inumFI2L(i-1))= jj
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            icouU = icouU + 1
            FI2U(icouU+inumFI2U(i-1))= jj
            IW1(jj)= 1
          endif
        enddo
      enddo
    enddo
    !C===

    !C
    !C +------------------+
    !C | find fill-in (2) |
    !C +------------------+
    !C===
    allocate (inumFI1L(0:NP), inumFI1U(0:NP))

    NPLf2= 0
    NPUf2= 0
    inumFI1L= 0
    inumFI1U= 0
    !C
    !C-- count fill-in
    do i= 2, NP
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 2
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 2
      enddo

      do L= inumFI2L(i-1)+1, inumFI2L(i)
        IW1(FI2L(L))= 1
      enddo

      do L= inumFI2U(i-1)+1, inumFI2U(i)
        IW1(FI2U(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= inumFI2U(kk-1) + 1
        iEj= inumFI2U(kk)
        do j= iSj, iEj
          jj= FI2U(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            inumFI1L(i)= inumFI1L(i) + 1
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            inumFI1U(i)= inumFI1U(i) + 1
            IW1(jj)= 1
          endif
        enddo
      enddo

      iSk= inumFI2L(i-1)+1
      iEk= inumFI2L(i)
      do k= iSk, iEk
        kk= FI2L(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            inumFI1L(i)= inumFI1L(i) + 1
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            inumFI1U(i)= inumFI1U(i) + 1
            IW1(jj)= 1
          endif
        enddo
      enddo
      NPLf2= NPLf2 + inumFI1L(i)
      NPUf2= NPUf2 + inumFI1U(i)
    enddo

    !C
    !C-- specify fill-in
    allocate (FI1L(NPL+NPLf1+NPLf2))
    allocate (FI1U(NPU+NPUf1+NPUf2))

    allocate (iconFI1L(NPL+NPLf1+NPLf2))
    allocate (iconFI1U(NPU+NPUf1+NPUf2))

    IWsL= 0
    IWsU= 0
    do i= 1, NP
      IWsL(i)= INL(i)-INL(i-1) + inumFI2L(i)-inumFI2L(i-1) +          &
        &                             inumFI1L(i) + IWsL(i-1)
      IWsU(i)= INU(i)-INU(i-1) + inumFI2U(i)-inumFI2U(i-1) +          &
        &                             inumFI1U(i) + IWsU(i-1)
    enddo

    do i= 2, NP
      icouL= 0
      icouU= 0
      inumFI1L(i)= inumFI1L(i-1) + inumFI1L(i)
      inumFI1U(i)= inumFI1U(i-1) + inumFI1U(i)
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= INL(i-1)+1, INL(i)
        IW1(IAL(L))= 1
      enddo
      do L= INU(i-1)+1, INU(i)
        IW1(IAU(L))= 1
      enddo

      do L= inumFI2L(i-1)+1, inumFI2L(i)
        IW1(FI2L(L))= 1
      enddo

      do L= inumFI2U(i-1)+1, inumFI2U(i)
        IW1(FI2U(L))= 1
      enddo

      iSk= INL(i-1) + 1
      iEk= INL(i)
      do k= iSk, iEk
        kk= IAL(k)
        iSj= inumFI2U(kk-1) + 1
        iEj= inumFI2U(kk  )
        do j= iSj, iEj
          jj= FI2U(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            iAS= INL(i)-INL(i-1)+inumFI2L(i)-inumFI2L(i-1)+IWsL(i-1)
            icouL     = icouL + 1
            FI1L(icouL+iAS)= jj
            IW1(jj)    = 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            iAS= INU(i)-INU(i-1)+inumFI2U(i)-inumFI2U(i-1)+IWsU(i-1)
            icouU     = icouU + 1
            FI1U(icouU+iAS)= jj
            IW1(jj)    = 1
          endif
        enddo
      enddo

      iSk= inumFI2L(i-1) + 1
      iEk= inumFI2L(i)
      do k= iSk, iEk
        kk= FI2L(k)
        iSj= INU(kk-1) + 1
        iEj= INU(kk  )
        do j= iSj, iEj
          jj= IAU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            iAS= INL(i)-INL(i-1)+inumFI2L(i)-inumFI2L(i-1)+IWsL(i-1)
            icouL     = icouL + 1
            FI1L(icouL+iAS)= jj
            IW1(jj)    = 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            iAS= INU(i)-INU(i-1)+inumFI2U(i)-inumFI2U(i-1)+IWsU(i-1)
            icouU     = icouU + 1
            FI1U(icouU+iAS)= jj
            IW1(jj)    = 1
          endif
        enddo
      enddo
    enddo
    !C===

    !C
    !C +-------------------------------------------------+
    !C | SORT and RECONSTRUCT matrix considering fill-in |
    !C +-------------------------------------------------+
    !C===
    allocate (ALlu0(9*(NPL+NPLf1+NPLf2)))
    allocate (AUlu0(9*(NPU+NPUf1+NPUf2)))

    ALlu0= 0.d0
    AUlu0= 0.d0
    iSL  = 0
    iSU  = 0

    iconFI1L= 0
    iconFI1U= 0

    do i= 1, NP
      icouL1=      INL(i) -      INL(i-1)
      icouL2= inumFI2L(i) - inumFI2L(i-1) + icouL1
      icouL3= inumFI1L(i) - inumFI1L(i-1) + icouL2

      icouU1=      INU(i) -      INU(i-1)
      icouU2= inumFI2U(i) - inumFI2U(i-1) + icouU1
      icouU3= inumFI1U(i) - inumFI1U(i-1) + icouU2

      !C
      !C-- LOWER part
      icou= 0
      do k= INL(i-1)+1, INL(i)
        icou = icou + 1
        IW1(icou)= IAL(k)
      enddo

      icou= 0
      do k= inumFI2L(i-1)+1, inumFI2L(i)
        icou        = icou + 1
        IW1(icou+icouL1)= FI2L(k)
      enddo

      icou= 0
      do k= inumFI1L(i-1)+1, inumFI1L(i)
        icou        = icou + 1
        IW1(icou+icouL2)= FI1L(icou+icouL2+iSL)
      enddo

      do k= 1, icouL3
        IW2(k)= k
      enddo

      call fill_in_S33_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 9*( k+isL)
          kk2= 9*(ik+INL(i-1))
          ALlu0(kk1-8)= AL(kk2-8)
          ALlu0(kk1-7)= AL(kk2-7)
          ALlu0(kk1-6)= AL(kk2-6)
          ALlu0(kk1-5)= AL(kk2-5)
          ALlu0(kk1-4)= AL(kk2-4)
          ALlu0(kk1-3)= AL(kk2-3)
          ALlu0(kk1-2)= AL(kk2-2)
          ALlu0(kk1-1)= AL(kk2-1)
          ALlu0(kk1  )= AL(kk2  )
        endif
      enddo

      icou= 0
      do k= INL(i-1)+1, INL(i)
        icou = icou + 1
        IW1(icou)= 0
      enddo

      icou= 0
      do k= inumFI2L(i-1)+1, inumFI2L(i)
        icou        = icou + 1
        IW1(icou+icouL1)= 1
      enddo

      icou= 0
      do k= inumFI1L(i-1)+1, inumFI1L(i)
        icou        = icou + 1
        IW1(icou+icouL2)= 2
      enddo

      do k= 1, icouL3
        iconFI1L(k+iSL)= IW1(IW2(k))
      enddo
      !C
      !C-- UPPER part
      icou= 0
      do k= INU(i-1)+1, INU(i)
        icou = icou + 1
        IW1(icou)= IAU(k)
      enddo

      icou= 0
      do k= inumFI2U(i-1)+1, inumFI2U(i)
        icou        = icou + 1
        IW1(icou+icouU1)= FI2U(k)
      enddo

      icou= 0
      do k= inumFI1U(i-1)+1, inumFI1U(i)
        icou        = icou + 1
        IW1(icou+icouU2)= FI1U(icou+icouU2+iSU)
      enddo

      do k= 1, icouU3
        IW2(k)= k
      enddo
      call fill_in_S33_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 9*( k+isU)
          kk2= 9*(ik+INU(i-1))
          AUlu0(kk1-8)= AU(kk2-8)
          AUlu0(kk1-7)= AU(kk2-7)
          AUlu0(kk1-6)= AU(kk2-6)
          AUlu0(kk1-5)= AU(kk2-5)
          AUlu0(kk1-4)= AU(kk2-4)
          AUlu0(kk1-3)= AU(kk2-3)
          AUlu0(kk1-2)= AU(kk2-2)
          AUlu0(kk1-1)= AU(kk2-1)
          AUlu0(kk1  )= AU(kk2  )
        endif
      enddo

      icou= 0
      do k= INU(i-1)+1, INU(i)
        icou = icou + 1
        IW1(icou)= 0
      enddo

      icou= 0
      do k= inumFI2U(i-1)+1, inumFI2U(i)
        icou        = icou + 1
        IW1(icou+icouU1)= 1
      enddo

      icou= 0
      do k= inumFI1U(i-1)+1, inumFI1U(i)
        icou        = icou + 1
        IW1(icou+icouU2)= 2
      enddo

      do k= 1, icouU3
        iconFI1U(k+iSU)= IW1(IW2(k))
      enddo

      iSL= iSL + icouL3
      iSU= iSU + icouU3
    enddo
    !C===
    do i= 1, NP
      inumFI1L(i)= IWsL(i)
      inumFI1U(i)= IWsU(i)
    enddo

    deallocate (IWsL, IWsU)
    deallocate (inumFI2L, inumFI2U)
    deallocate (    FI2L,     FI2U)

    !C
    !C +----------------------+
    !C | ILU(2) factorization |
    !C +----------------------+
    !C===
    allocate (Dlu0(9*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(9*i-8)=Dlu0(9*i-8)*SIGMA_DIAG
      Dlu0(9*i-4)=Dlu0(9*i-4)*SIGMA_DIAG
      Dlu0(9*i  )=Dlu0(9*i  )*SIGMA_DIAG
    enddo

    i = 1
    call ILU1a33 (DkINV, &
      Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
      Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
      Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
    Dlu0(9*i-8)= DkINV(1,1)
    Dlu0(9*i-7)= DkINV(1,2)
    Dlu0(9*i-6)= DkINV(1,3)
    Dlu0(9*i-5)= DkINV(2,1)
    Dlu0(9*i-4)= DkINV(2,2)
    Dlu0(9*i-3)= DkINV(2,3)
    Dlu0(9*i-2)= DkINV(3,1)
    Dlu0(9*i-1)= DkINV(3,2)
    Dlu0(9*i  )= DkINV(3,3)

    do i= 2, NP
      IW1= 0
      IW2= 0

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        IW1(FI1L(k))= k
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        IW2(FI1U(k))= k
      enddo

      do kk= inumFI1L(i-1)+1, inumFI1L(i)
        k= FI1L(kk)
        iconIK= iconFI1L(kk)

        DkINV(1,1)= Dlu0(9*k-8)
        DkINV(1,2)= Dlu0(9*k-7)
        DkINV(1,3)= Dlu0(9*k-6)
        DkINV(2,1)= Dlu0(9*k-5)
        DkINV(2,2)= Dlu0(9*k-4)
        DkINV(2,3)= Dlu0(9*k-3)
        DkINV(3,1)= Dlu0(9*k-2)
        DkINV(3,2)= Dlu0(9*k-1)
        DkINV(3,3)= Dlu0(9*k  )

        Aik(1,1)= ALlu0(9*kk-8)
        Aik(1,2)= ALlu0(9*kk-7)
        Aik(1,3)= ALlu0(9*kk-6)
        Aik(2,1)= ALlu0(9*kk-5)
        Aik(2,2)= ALlu0(9*kk-4)
        Aik(2,3)= ALlu0(9*kk-3)
        Aik(3,1)= ALlu0(9*kk-2)
        Aik(3,2)= ALlu0(9*kk-1)
        Aik(3,3)= ALlu0(9*kk  )

        do jj= inumFI1U(k-1)+1, inumFI1U(k)
          j= FI1U(jj)
          iconKJ= iconFI1U(jj)

          if ((iconIK+iconKJ).lt.2) then
            Akj(1,1)= AUlu0(9*jj-8)
            Akj(1,2)= AUlu0(9*jj-7)
            Akj(1,3)= AUlu0(9*jj-6)
            Akj(2,1)= AUlu0(9*jj-5)
            Akj(2,2)= AUlu0(9*jj-4)
            Akj(2,3)= AUlu0(9*jj-3)
            Akj(3,1)= AUlu0(9*jj-2)
            Akj(3,2)= AUlu0(9*jj-1)
            Akj(3,3)= AUlu0(9*jj  )

            call ILU1b33 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(9*i-8)= Dlu0(9*i-8) - RHS_Aij(1,1)
              Dlu0(9*i-7)= Dlu0(9*i-7) - RHS_Aij(1,2)
              Dlu0(9*i-6)= Dlu0(9*i-6) - RHS_Aij(1,3)
              Dlu0(9*i-5)= Dlu0(9*i-5) - RHS_Aij(2,1)
              Dlu0(9*i-4)= Dlu0(9*i-4) - RHS_Aij(2,2)
              Dlu0(9*i-3)= Dlu0(9*i-3) - RHS_Aij(2,3)
              Dlu0(9*i-2)= Dlu0(9*i-2) - RHS_Aij(3,1)
              Dlu0(9*i-1)= Dlu0(9*i-1) - RHS_Aij(3,2)
              Dlu0(9*i  )= Dlu0(9*i  ) - RHS_Aij(3,3)
            endif

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(9*ij0-8)= ALlu0(9*ij0-8) - RHS_Aij(1,1)
              ALlu0(9*ij0-7)= ALlu0(9*ij0-7) - RHS_Aij(1,2)
              ALlu0(9*ij0-6)= ALlu0(9*ij0-6) - RHS_Aij(1,3)
              ALlu0(9*ij0-5)= ALlu0(9*ij0-5) - RHS_Aij(2,1)
              ALlu0(9*ij0-4)= ALlu0(9*ij0-4) - RHS_Aij(2,2)
              ALlu0(9*ij0-3)= ALlu0(9*ij0-3) - RHS_Aij(2,3)
              ALlu0(9*ij0-2)= ALlu0(9*ij0-2) - RHS_Aij(3,1)
              ALlu0(9*ij0-1)= ALlu0(9*ij0-1) - RHS_Aij(3,2)
              ALlu0(9*ij0  )= ALlu0(9*ij0  ) - RHS_Aij(3,3)
            endif

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(9*ij0-8)= AUlu0(9*ij0-8) - RHS_Aij(1,1)
              AUlu0(9*ij0-7)= AUlu0(9*ij0-7) - RHS_Aij(1,2)
              AUlu0(9*ij0-6)= AUlu0(9*ij0-6) - RHS_Aij(1,3)
              AUlu0(9*ij0-5)= AUlu0(9*ij0-5) - RHS_Aij(2,1)
              AUlu0(9*ij0-4)= AUlu0(9*ij0-4) - RHS_Aij(2,2)
              AUlu0(9*ij0-3)= AUlu0(9*ij0-3) - RHS_Aij(2,3)
              AUlu0(9*ij0-2)= AUlu0(9*ij0-2) - RHS_Aij(3,1)
              AUlu0(9*ij0-1)= AUlu0(9*ij0-1) - RHS_Aij(3,2)
              AUlu0(9*ij0  )= AUlu0(9*ij0  ) - RHS_Aij(3,3)
            endif
          endif
        enddo
      enddo

      call ILU1a33 (DkINV, &
        Dlu0(9*i-8), Dlu0(9*i-7), Dlu0(9*i-6), &
        Dlu0(9*i-5), Dlu0(9*i-4), Dlu0(9*i-3), &
        Dlu0(9*i-2), Dlu0(9*i-1), Dlu0(9*i  ))
      Dlu0(9*i-8)= DkINV(1,1)
      Dlu0(9*i-7)= DkINV(1,2)
      Dlu0(9*i-6)= DkINV(1,3)
      Dlu0(9*i-5)= DkINV(2,1)
      Dlu0(9*i-4)= DkINV(2,2)
      Dlu0(9*i-3)= DkINV(2,3)
      Dlu0(9*i-2)= DkINV(3,1)
      Dlu0(9*i-1)= DkINV(3,2)
      Dlu0(9*i  )= DkINV(3,3)
    enddo

    deallocate (IW1, IW2)
    deallocate (iconFI1L, iconFI1U)
    !C===
  end subroutine FORM_ILU2_33


  !C
  !C***
  !C*** fill_in_S33_SORT
  !C***
  !C
  subroutine fill_in_S33_SORT (STEM, INUM, N, NP)
    use hecmw_util
    implicit none
    integer(kind=kint) :: N, NP
    integer(kind=kint), dimension(NP) :: STEM
    integer(kind=kint), dimension(NP) :: INUM
    integer(kind=kint), dimension(:), allocatable :: ISTACK
    integer(kind=kint) :: M,NSTACK,jstack,l,ir,ip,i,j,k,ss,ii,temp,it

    allocate (ISTACK(-NP:+NP))

    M     = 100
    NSTACK= NP

    jstack= 0
    l     = 1
    ir    = N

    ip= 0
    1   continue
    ip= ip + 1

    if (ir-l.lt.M) then
      do j= l+1, ir
        ss= STEM(j)
        ii= INUM(j)

        do i= j-1,1,-1
          if (STEM(i).le.ss) goto 2
          STEM(i+1)= STEM(i)
          INUM(i+1)= INUM(i)
        end do
        i= 0

        2       continue
        STEM(i+1)= ss
        INUM(i+1)= ii
      end do

      if (jstack.eq.0) then
        deallocate (ISTACK)
        return
      endif

      ir = ISTACK(jstack)
      l = ISTACK(jstack-1)
      jstack= jstack - 2
    else

      k= (l+ir) / 2
      temp = STEM(k)
      STEM(k)  = STEM(l+1)
      STEM(l+1)= temp

      it = INUM(k)
      INUM(k)  = INUM(l+1)
      INUM(l+1)= it

      if (STEM(l+1).gt.STEM(ir)) then
        temp = STEM(l+1)
        STEM(l+1)= STEM(ir)
        STEM(ir )= temp
        it = INUM(l+1)
        INUM(l+1)= INUM(ir)
        INUM(ir )= it
      endif

      if (STEM(l).gt.STEM(ir)) then
        temp = STEM(l)
        STEM(l )= STEM(ir)
        STEM(ir)= temp
        it = INUM(l)
        INUM(l )= INUM(ir)
        INUM(ir)= it
      endif

      if (STEM(l+1).gt.STEM(l)) then
        temp = STEM(l+1)
        STEM(l+1)= STEM(l)
        STEM(l  )= temp
        it = INUM(l+1)
        INUM(l+1)= INUM(l)
        INUM(l  )= it
      endif

      i= l + 1
      j= ir

      ss= STEM(l)
      ii= INUM(l)

      3     continue
      i= i + 1
      if (STEM(i).lt.ss) goto 3

      4     continue
      j= j - 1
      if (STEM(j).gt.ss) goto 4

      if (j.lt.i)        goto 5

      temp   = STEM(i)
      STEM(i)= STEM(j)
      STEM(j)= temp

      it     = INUM(i)
      INUM(i)= INUM(j)
      INUM(j)= it

      goto 3

      5     continue

      STEM(l)= STEM(j)
      STEM(j)= ss
      INUM(l)= INUM(j)
      INUM(j)= ii

      jstack= jstack + 2

      if (jstack.gt.NSTACK) then
        write (*,*) 'NSTACK overflow'
        stop
      endif

      if (ir-i+1.ge.j-1) then
        ISTACK(jstack  )= ir
        ISTACK(jstack-1)= i
        ir= j-1
      else
        ISTACK(jstack  )= j-1
        ISTACK(jstack-1)= l
        l= i
      endif

    endif

    goto 1

  end subroutine fill_in_S33_SORT

  !C
  !C***
  !C*** ILU1a33
  !C***
  !C
  !C    computes LU factorization of 3*3 Diagonal Block
  !C
  subroutine ILU1a33 (ALU, D11,D12,D13,D21,D22,D23,D31,D32,D33)
    use hecmw_util
    implicit none
    real(kind=kreal) :: ALU(3,3), PW(3)
    real(kind=kreal) :: D11,D12,D13,D21,D22,D23,D31,D32,D33
    integer(kind=kint) :: i,j,k

    ALU(1,1)= D11
    ALU(1,2)= D12
    ALU(1,3)= D13
    ALU(2,1)= D21
    ALU(2,2)= D22
    ALU(2,3)= D23
    ALU(3,1)= D31
    ALU(3,2)= D32
    ALU(3,3)= D33

    do k= 1, 3
      if (ALU(k,k) == 0.d0) then
        !write(*,*) ALU(1:3,1:3)
        stop 'ERROR: Divide by zero in ILU setup'
      endif
      ALU(k,k)= 1.d0/ALU(k,k)
      do i= k+1, 3
        ALU(i,k)= ALU(i,k) * ALU(k,k)
        do j= k+1, 3
          PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
        enddo
        do j= k+1, 3
          ALU(i,j)= PW(j)
        enddo
      enddo
    enddo

    return
  end subroutine ILU1a33

  !C
  !C***
  !C*** ILU1b33
  !C***
  !C
  !C    computes L_ik * D_k_INV * U_kj at ILU factorization
  !C    for 3*3 Block Type Matrix
  !C
  subroutine ILU1b33 (RHS_Aij, DkINV, Aik, Akj)
    use hecmw_util
    implicit none
    real(kind=kreal) :: RHS_Aij(3,3), DkINV(3,3), Aik(3,3), Akj(3,3)
    real(kind=kreal) :: X1,X2,X3

    !C
    !C-- 1st Col.
    X1= Akj(1,1)
    X2= Akj(2,1)
    X3= Akj(3,1)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

    X3= DkINV(3,3)*  X3
    X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,1)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
    RHS_Aij(2,1)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
    RHS_Aij(3,1)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

    !C
    !C-- 2nd Col.
    X1= Akj(1,2)
    X2= Akj(2,2)
    X3= Akj(3,2)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

    X3= DkINV(3,3)*  X3
    X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,2)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
    RHS_Aij(2,2)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
    RHS_Aij(3,2)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

    !C
    !C-- 3rd Col.
    X1= Akj(1,3)
    X2= Akj(2,3)
    X3= Akj(3,3)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2

    X3= DkINV(3,3)*  X3
    X2= DkINV(2,2)*( X2 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,3)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3
    RHS_Aij(2,3)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3
    RHS_Aij(3,3)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3

    return
  end subroutine ILU1b33

end module     hecmw_precond_BILU_33
