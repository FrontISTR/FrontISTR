!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_BILU_44
!C***
!C
module hecmw_precond_BILU_44
  use hecmw_util
  use hecmw_matrix_misc

  private

  public:: hecmw_precond_BILU_44_setup
  public:: hecmw_precond_BILU_44_apply
  public:: hecmw_precond_BILU_44_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: Dlu0(:) => null()
  real(kind=kreal), pointer :: ALlu0(:) => null()
  real(kind=kreal), pointer :: AUlu0(:) => null()
  integer(kind=kint), pointer :: inumFI1L(:) => null()
  integer(kind=kint), pointer :: inumFI1U(:) => null()
  integer(kind=kint), pointer :: FI1L(:) => null()
  integer(kind=kint), pointer :: FI1U(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_BILU_44_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP, NPU, NPL
    integer(kind=kint ) :: PRECOND
    real   (kind=kreal) :: SIGMA, SIGMA_DIAG

    real(kind=kreal), pointer :: D(:)
    real(kind=kreal), pointer :: AL(:)
    real(kind=kreal), pointer :: AU(:)

    integer(kind=kint ), pointer :: INL(:), INU(:)
    integer(kind=kint ), pointer :: IAL(:)
    integer(kind=kint ), pointer :: IAU(:)

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 1) then ! need symbolic and numerical setup
        call hecmw_precond_BILU_44_clear()
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_BILU_44_clear() ! TEMPORARY
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

    !if (PRECOND.eq.10) call FORM_ILU0_44 &
      call FORM_ILU0_44 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    !if (PRECOND.eq.11) call FORM_ILU1_44 &
      !     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      !     &    SIGMA, SIGMA_DIAG)
    !if (PRECOND.eq.12) call FORM_ILU2_44 &
      !     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      !     &    SIGMA, SIGMA_DIAG)

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_BILU_44_setup

  subroutine hecmw_precond_BILU_44_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, SW4, X1, X2, X3, X4
    !C
    !C-- FORWARD

    do i= 1, N
      SW1= WW(4*i-3)
      SW2= WW(4*i-2)
      SW3= WW(4*i-1)
      SW4= WW(4*i  )
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      do j= isL, ieL
        k= FI1L(j)
        X1= WW(4*k-3)
        X2= WW(4*k-2)
        X3= WW(4*k-1)
        X4= WW(4*k  )
        SW1= SW1 - ALlu0(16*j-15)*X1-ALlu0(16*j-14)*X2-ALlu0(16*j-13)*X3-ALlu0(16*j-12)*X4
        SW2= SW2 - ALlu0(16*j-11)*X1-ALlu0(16*j-10)*X2-ALlu0(16*j- 9)*X3-ALlu0(16*j- 8)*X4
        SW3= SW3 - ALlu0(16*j- 7)*X1-ALlu0(16*j- 6)*X2-ALlu0(16*j- 5)*X3-ALlu0(16*j- 4)*X4
        SW4= SW4 - ALlu0(16*j- 3)*X1-ALlu0(16*j- 2)*X2-ALlu0(16*j- 1)*X3-ALlu0(16*j   )*X4
      enddo

      X1= SW1
      X2= SW2
      X3= SW3
      X4= SW4
      X2= X2 - Dlu0(16*i-11)*X1
      X3= X3 - Dlu0(16*i- 7)*X1 - Dlu0(16*i-6)*X2
      X4= X4 - Dlu0(16*i- 3)*X1 - Dlu0(16*i-2)*X2 - Dlu0(16*i-1)*X3
      X4= Dlu0(16*i   )* X4
      X3= Dlu0(16*i- 5)*(X3 - Dlu0(16*i- 4)*X4)
      X2= Dlu0(16*i-10)*(X2 - Dlu0(16*i- 8)*X4 - Dlu0(16*i- 9)*X3 )
      X1= Dlu0(16*i-15)*(X1 - Dlu0(16*i-12)*X4 - Dlu0(16*i-13)*X3 - Dlu0(16*i-14)*X2)

      WW(4*i-3)= X1
      WW(4*i-2)= X2
      WW(4*i-1)= X3
      WW(4*i  )= X4
    enddo

    !C
    !C-- BACKWARD

    do i= N, 1, -1
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      SW1= 0.d0
      SW2= 0.d0
      SW3= 0.d0
      SW4= 0.d0
      do j= ieU, isU, -1
        k= FI1U(j)
        X1= WW(4*k-3)
        X2= WW(4*k-2)
        X3= WW(4*k-1)
        X4= WW(4*k  )
        SW1= SW1 + AUlu0(16*j-15)*X1+AUlu0(16*j-14)*X2+AUlu0(16*j-13)*X3+AUlu0(16*j-12)*X4
        SW2= SW2 + AUlu0(16*j-11)*X1+AUlu0(16*j-10)*X2+AUlu0(16*j- 9)*X3+AUlu0(16*j- 8)*X4
        SW3= SW3 + AUlu0(16*j- 7)*X1+AUlu0(16*j- 6)*X2+AUlu0(16*j- 5)*X3+AUlu0(16*j- 4)*X4
        SW4= SW4 + AUlu0(16*j- 3)*X1+AUlu0(16*j- 2)*X2+AUlu0(16*j- 1)*X3+AUlu0(16*j   )*X4
      enddo
      X1= SW1
      X2= SW2
      X3= SW3
      X4= SW4
      X2= X2 - Dlu0(16*i-11)*X1
      X3= X3 - Dlu0(16*i- 7)*X1 - Dlu0(16*i-6)*X2
      X4= X4 - Dlu0(16*i- 3)*X1 - Dlu0(16*i-2)*X2 - Dlu0(16*i-1)*X3
      X4= Dlu0(16*i   )*  X4
      X3= Dlu0(16*i- 5)*( X3 - Dlu0(16*i- 4)*X4 )
      X2= Dlu0(16*i-10)*( X2 - Dlu0(16*i- 8)*X4 - Dlu0(16*i- 9)*X3 )
      X1= Dlu0(16*i-15)*( X1 - Dlu0(16*i-12)*X4 - Dlu0(16*i-13)*X3 - Dlu0(16*i-14)*X2)
      WW(4*i-3)=  WW(4*i-3) - X1
      WW(4*i-2)=  WW(4*i-2) - X2
      WW(4*i-1)=  WW(4*i-1) - X3
      WW(4*i  )=  WW(4*i  ) - X4
    enddo
  end subroutine hecmw_precond_BILU_44_apply

  subroutine hecmw_precond_BILU_44_clear()
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
  end subroutine hecmw_precond_BILU_44_clear

  !C
  !C***
  !C*** FORM_ILU1_44
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU0_44                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(16*NPL), intent(in):: AL
    real(kind=kreal), dimension(16*NPU), intent(in):: AU
    real(kind=kreal), dimension(16*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    real (kind=kreal),  dimension(4,4) :: RHS_Aij, DkINV, Aik, Akj
    integer(kind=kint) :: i,jj,jj1,ij0,kk,kk1
    integer(kind=kint) :: j,k
    allocate (IW1(NP) , IW2(NP))
    allocate(Dlu0(9*NP), ALlu0(9*NPL), AUlu0(9*NPU))
    allocate(inumFI1L(0:NP), inumFI1U(0:NP), FI1L(NPL), FI1U(NPU))

    do i=1,16*NP
      Dlu0(i) = D(i)
    end do
    do i=1,16*NPL
      ALlu0(i) = AL(i)
    end do
    do i=1,16*NPU
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
      Dlu0(16*i-15)=Dlu0(16*i-15)*SIGMA_DIAG
      Dlu0(16*i-10)=Dlu0(16*i-10)*SIGMA_DIAG
      Dlu0(16*i- 5)=Dlu0(16*i- 5)*SIGMA_DIAG
      Dlu0(16*i   )=Dlu0(16*i   )*SIGMA_DIAG
    enddo

    i = 1
    call ILU1a44 (DkINV, &
      Dlu0(16*i-15), Dlu0(16*i-14), Dlu0(16*i-13), Dlu0(16*i-12), &
      Dlu0(16*i-11), Dlu0(16*i-10), Dlu0(16*i- 9), Dlu0(16*i- 8), &
      Dlu0(16*i- 7), Dlu0(16*i- 6), Dlu0(16*i- 5), Dlu0(16*i- 4), &
      Dlu0(16*i- 3), Dlu0(16*i- 2), Dlu0(16*i- 1), Dlu0(16*i   ) )
    Dlu0(16*i-15)= DkINV(1,1)
    Dlu0(16*i-14)= DkINV(1,2)
    Dlu0(16*i-13)= DkINV(1,3)
    Dlu0(16*i-12)= DkINV(1,4)
    Dlu0(16*i-11)= DkINV(2,1)
    Dlu0(16*i-10)= DkINV(2,2)
    Dlu0(16*i- 9)= DkINV(2,3)
    Dlu0(16*i- 8)= DkINV(2,4)
    Dlu0(16*i- 7)= DkINV(3,1)
    Dlu0(16*i- 6)= DkINV(3,2)
    Dlu0(16*i- 5)= DkINV(3,3)
    Dlu0(16*i- 4)= DkINV(3,4)
    Dlu0(16*i- 3)= DkINV(4,1)
    Dlu0(16*i- 2)= DkINV(4,2)
    Dlu0(16*i- 1)= DkINV(4,3)
    Dlu0(16*i   )= DkINV(4,4)

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

        DkINV(1,1)= Dlu0(16*k-15)
        DkINV(1,2)= Dlu0(16*k-14)
        DkINV(1,3)= Dlu0(16*k-13)
        DkINV(1,4)= Dlu0(16*k-12)
        DkINV(2,1)= Dlu0(16*k-11)
        DkINV(2,2)= Dlu0(16*k-10)
        DkINV(2,3)= Dlu0(16*k- 9)
        DkINV(2,4)= Dlu0(16*k- 8)
        DkINV(3,1)= Dlu0(16*k- 7)
        DkINV(3,2)= Dlu0(16*k- 6)
        DkINV(3,3)= Dlu0(16*k- 5)
        DkINV(3,4)= Dlu0(16*k- 4)
        DkINV(4,1)= Dlu0(16*k- 3)
        DkINV(4,2)= Dlu0(16*k- 2)
        DkINV(4,3)= Dlu0(16*k- 1)
        DkINV(4,4)= Dlu0(16*k   )

        Aik(1,1)= ALlu0(16*kk-15)
        Aik(1,2)= ALlu0(16*kk-14)
        Aik(1,3)= ALlu0(16*kk-13)
        Aik(1,4)= ALlu0(16*kk-12)
        Aik(2,1)= ALlu0(16*kk-11)
        Aik(2,2)= ALlu0(16*kk-10)
        Aik(2,3)= ALlu0(16*kk- 9)
        Aik(2,4)= ALlu0(16*kk- 8)
        Aik(3,1)= ALlu0(16*kk- 7)
        Aik(3,2)= ALlu0(16*kk- 6)
        Aik(3,3)= ALlu0(16*kk- 5)
        Aik(3,4)= ALlu0(16*kk- 4)
        Aik(4,1)= ALlu0(16*kk- 3)
        Aik(4,2)= ALlu0(16*kk- 2)
        Aik(4,3)= ALlu0(16*kk- 1)
        Aik(4,4)= ALlu0(16*kk   )

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          if (IW1(j).eq.0.and.IW2(j).eq.0) cycle

          Akj(1,1)= AUlu0(16*jj-15)
          Akj(1,2)= AUlu0(16*jj-14)
          Akj(1,3)= AUlu0(16*jj-13)
          Akj(1,4)= AUlu0(16*jj-12)
          Akj(2,1)= AUlu0(16*jj-11)
          Akj(2,2)= AUlu0(16*jj-10)
          Akj(2,3)= AUlu0(16*jj- 9)
          Akj(2,4)= AUlu0(16*jj- 8)
          Akj(3,1)= AUlu0(16*jj- 7)
          Akj(3,2)= AUlu0(16*jj- 6)
          Akj(3,3)= AUlu0(16*jj- 5)
          Akj(3,4)= AUlu0(16*jj- 4)
          Akj(4,1)= AUlu0(16*jj- 3)
          Akj(4,2)= AUlu0(16*jj- 2)
          Akj(4,3)= AUlu0(16*jj- 1)
          Akj(4,4)= AUlu0(16*jj   )

          call ILU1b44 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(16*i-15)= Dlu0(16*i-15) - RHS_Aij(1,1)
            Dlu0(16*i-14)= Dlu0(16*i-14) - RHS_Aij(1,2)
            Dlu0(16*i-13)= Dlu0(16*i-13) - RHS_Aij(1,3)
            Dlu0(16*i-12)= Dlu0(16*i-12) - RHS_Aij(1,4)
            Dlu0(16*i-11)= Dlu0(16*i-11) - RHS_Aij(2,1)
            Dlu0(16*i-10)= Dlu0(16*i-10) - RHS_Aij(2,2)
            Dlu0(16*i- 9)= Dlu0(16*i- 9) - RHS_Aij(2,3)
            Dlu0(16*i- 8)= Dlu0(16*i- 8) - RHS_Aij(2,4)
            Dlu0(16*i- 7)= Dlu0(16*i- 7) - RHS_Aij(3,1)
            Dlu0(16*i- 6)= Dlu0(16*i- 6) - RHS_Aij(3,2)
            Dlu0(16*i- 5)= Dlu0(16*i- 5) - RHS_Aij(3,3)
            Dlu0(16*i- 4)= Dlu0(16*i- 4) - RHS_Aij(3,4)
            Dlu0(16*i- 3)= Dlu0(16*i- 3) - RHS_Aij(4,1)
            Dlu0(16*i- 2)= Dlu0(16*i- 2) - RHS_Aij(4,2)
            Dlu0(16*i- 1)= Dlu0(16*i- 1) - RHS_Aij(4,3)
            Dlu0(16*i   )= Dlu0(16*i   ) - RHS_Aij(4,4)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(16*ij0-15)= ALlu0(16*ij0-15) - RHS_Aij(1,1)
            ALlu0(16*ij0-14)= ALlu0(16*ij0-14) - RHS_Aij(1,2)
            ALlu0(16*ij0-13)= ALlu0(16*ij0-13) - RHS_Aij(1,3)
            ALlu0(16*ij0-12)= ALlu0(16*ij0-12) - RHS_Aij(1,4)
            ALlu0(16*ij0-11)= ALlu0(16*ij0-11) - RHS_Aij(2,1)
            ALlu0(16*ij0-10)= ALlu0(16*ij0-10) - RHS_Aij(2,2)
            ALlu0(16*ij0- 9)= ALlu0(16*ij0- 9) - RHS_Aij(2,3)
            ALlu0(16*ij0- 8)= ALlu0(16*ij0- 8) - RHS_Aij(2,4)
            ALlu0(16*ij0- 7)= ALlu0(16*ij0- 7) - RHS_Aij(3,1)
            ALlu0(16*ij0- 6)= ALlu0(16*ij0- 6) - RHS_Aij(3,2)
            ALlu0(16*ij0- 5)= ALlu0(16*ij0- 5) - RHS_Aij(3,3)
            ALlu0(16*ij0- 4)= ALlu0(16*ij0- 4) - RHS_Aij(3,4)
            ALlu0(16*ij0- 3)= ALlu0(16*ij0- 3) - RHS_Aij(4,1)
            ALlu0(16*ij0- 2)= ALlu0(16*ij0- 2) - RHS_Aij(4,2)
            ALlu0(16*ij0- 1)= ALlu0(16*ij0- 1) - RHS_Aij(4,3)
            ALlu0(16*ij0   )= ALlu0(16*ij0   ) - RHS_Aij(4,4)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(16*ij0-15)= AUlu0(16*ij0-15) - RHS_Aij(1,1)
            AUlu0(16*ij0-14)= AUlu0(16*ij0-14) - RHS_Aij(1,2)
            AUlu0(16*ij0-13)= AUlu0(16*ij0-13) - RHS_Aij(1,3)
            AUlu0(16*ij0-12)= AUlu0(16*ij0-12) - RHS_Aij(1,4)
            AUlu0(16*ij0-11)= AUlu0(16*ij0-11) - RHS_Aij(2,1)
            AUlu0(16*ij0-10)= AUlu0(16*ij0-10) - RHS_Aij(2,2)
            AUlu0(16*ij0- 9)= AUlu0(16*ij0- 9) - RHS_Aij(2,3)
            AUlu0(16*ij0- 8)= AUlu0(16*ij0- 8) - RHS_Aij(2,4)
            AUlu0(16*ij0- 7)= AUlu0(16*ij0- 7) - RHS_Aij(3,1)
            AUlu0(16*ij0- 6)= AUlu0(16*ij0- 6) - RHS_Aij(3,2)
            AUlu0(16*ij0- 5)= AUlu0(16*ij0- 5) - RHS_Aij(3,3)
            AUlu0(16*ij0- 4)= AUlu0(16*ij0- 4) - RHS_Aij(3,4)
            AUlu0(16*ij0- 3)= AUlu0(16*ij0- 3) - RHS_Aij(4,1)
            AUlu0(16*ij0- 2)= AUlu0(16*ij0- 2) - RHS_Aij(4,2)
            AUlu0(16*ij0- 1)= AUlu0(16*ij0- 1) - RHS_Aij(4,3)
            AUlu0(16*ij0   )= AUlu0(16*ij0   ) - RHS_Aij(4,4)
          endif

        enddo
      enddo

      call ILU1a44 (DkINV, &
        Dlu0(16*i-15), Dlu0(16*i-14), Dlu0(16*i-13), Dlu0(16*i-12), &
        Dlu0(16*i-11), Dlu0(16*i-10), Dlu0(16*i- 9), Dlu0(16*i- 8), &
        Dlu0(16*i- 7), Dlu0(16*i- 6), Dlu0(16*i- 5), Dlu0(16*i- 4), &
        Dlu0(16*i- 3), Dlu0(16*i- 2), Dlu0(16*i- 1), Dlu0(16*i   ) )
      Dlu0(16*i-15)= DkINV(1,1)
      Dlu0(16*i-14)= DkINV(1,2)
      Dlu0(16*i-13)= DkINV(1,3)
      Dlu0(16*i-12)= DkINV(1,4)
      Dlu0(16*i-11)= DkINV(2,1)
      Dlu0(16*i-10)= DkINV(2,2)
      Dlu0(16*i- 9)= DkINV(2,3)
      Dlu0(16*i- 8)= DkINV(2,4)
      Dlu0(16*i- 7)= DkINV(3,1)
      Dlu0(16*i- 6)= DkINV(3,2)
      Dlu0(16*i- 5)= DkINV(3,3)
      Dlu0(16*i- 4)= DkINV(3,4)
      Dlu0(16*i- 3)= DkINV(4,1)
      Dlu0(16*i- 2)= DkINV(4,2)
      Dlu0(16*i- 1)= DkINV(4,3)
      Dlu0(16*i   )= DkINV(4,4)
    enddo

    deallocate (IW1, IW2)
  end subroutine FORM_ILU0_44

  !C
  !C***
  !C*** FORM_ILU1_44
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU1_44                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(16*NPL), intent(in):: AL
    real(kind=kreal), dimension(16*NPU), intent(in):: AU
    real(kind=kreal), dimension(16*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    integer(kind=kint), dimension(:), allocatable :: IWsL, IWsU
    real (kind=kreal),  dimension(4,4) :: RHS_Aij, DkINV, Aik, Akj
    real (kind=kreal)  :: D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44
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
    allocate (ALlu0(16*(NPL+NPLf1)), AUlu0(16*(NPU+NPUf1)))

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
      call fill_in_S44_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 16*( k+isL)
          kk2= 16*(ik+INL(i-1))
          ALlu0(kk1-15)= AL(kk2-15)
          ALlu0(kk1-14)= AL(kk2-14)
          ALlu0(kk1-13)= AL(kk2-13)
          ALlu0(kk1-12)= AL(kk2-12)
          ALlu0(kk1-11)= AL(kk2-11)
          ALlu0(kk1-10)= AL(kk2-10)
          ALlu0(kk1- 9)= AL(kk2- 9)
          ALlu0(kk1- 8)= AL(kk2- 8)
          ALlu0(kk1- 7)= AL(kk2- 7)
          ALlu0(kk1- 6)= AL(kk2- 6)
          ALlu0(kk1- 5)= AL(kk2- 5)
          ALlu0(kk1- 4)= AL(kk2- 4)
          ALlu0(kk1- 3)= AL(kk2- 3)
          ALlu0(kk1- 2)= AL(kk2- 2)
          ALlu0(kk1- 1)= AL(kk2- 1)
          ALlu0(kk1- 0)= AL(kk2- 0)
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
      call fill_in_S44_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 16*( k+isU)
          kk2= 16*(ik+INU(i-1))
          AUlu0(kk1-15)= AU(kk2-15)
          AUlu0(kk1-14)= AU(kk2-14)
          AUlu0(kk1-13)= AU(kk2-13)
          AUlu0(kk1-12)= AU(kk2-12)
          AUlu0(kk1-11)= AU(kk2-11)
          AUlu0(kk1-10)= AU(kk2-10)
          AUlu0(kk1- 9)= AU(kk2- 9)
          AUlu0(kk1- 8)= AU(kk2- 8)
          AUlu0(kk1- 7)= AU(kk2- 7)
          AUlu0(kk1- 6)= AU(kk2- 6)
          AUlu0(kk1- 5)= AU(kk2- 5)
          AUlu0(kk1- 4)= AU(kk2- 4)
          AUlu0(kk1- 3)= AU(kk2- 3)
          AUlu0(kk1- 2)= AU(kk2- 2)
          AUlu0(kk1- 1)= AU(kk2- 1)
          AUlu0(kk1- 0)= AU(kk2- 0)
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
    allocate (Dlu0(16*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(16*i-15)=Dlu0(16*i-15)*SIGMA_DIAG
      Dlu0(16*i-10)=Dlu0(16*i-10)*SIGMA_DIAG
      Dlu0(16*i- 5)=Dlu0(16*i- 5)*SIGMA_DIAG
      Dlu0(16*i- 0)=Dlu0(16*i- 0)*SIGMA_DIAG
    enddo

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
        D11= Dlu0(16*k-15)
        D12= Dlu0(16*k-14)
        D13= Dlu0(16*k-13)
        D14= Dlu0(16*k-12)
        D21= Dlu0(16*k-11)
        D22= Dlu0(16*k-10)
        D23= Dlu0(16*k- 9)
        D24= Dlu0(16*k- 8)
        D31= Dlu0(16*k- 7)
        D32= Dlu0(16*k- 6)
        D33= Dlu0(16*k- 5)
        D34= Dlu0(16*k- 4)
        D41= Dlu0(16*k- 3)
        D42= Dlu0(16*k- 2)
        D43= Dlu0(16*k- 1)
        D44= Dlu0(16*k- 0)

        call ILU1a44 (DkINV, D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44)

        do kk1= inumFI1L(i-1)+1, inumFI1L(i)
          if (k.eq.FI1L(kk1)) then
            Aik(1,1)= ALlu0(16*kk1-15)
            Aik(1,2)= ALlu0(16*kk1-14)
            Aik(1,3)= ALlu0(16*kk1-13)
            Aik(1,4)= ALlu0(16*kk1-12)
            Aik(2,1)= ALlu0(16*kk1-11)
            Aik(2,2)= ALlu0(16*kk1-10)
            Aik(2,3)= ALlu0(16*kk1- 9)
            Aik(2,4)= ALlu0(16*kk1- 8)
            Aik(3,1)= ALlu0(16*kk1- 7)
            Aik(3,2)= ALlu0(16*kk1- 6)
            Aik(3,3)= ALlu0(16*kk1- 5)
            Aik(3,4)= ALlu0(16*kk1- 4)
            Aik(4,1)= ALlu0(16*kk1- 3)
            Aik(4,2)= ALlu0(16*kk1- 2)
            Aik(4,3)= ALlu0(16*kk1- 1)
            Aik(4,4)= ALlu0(16*kk1- 0)
            exit
          endif
        enddo

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          do jj1= inumFI1U(k-1)+1, inumFI1U(k)
            if (j.eq.FI1U(jj1)) then
              Akj(1,1)= AUlu0(16*jj1-15)
              Akj(1,2)= AUlu0(16*jj1-14)
              Akj(1,3)= AUlu0(16*jj1-13)
              Akj(1,4)= AUlu0(16*jj1-12)
              Akj(2,1)= AUlu0(16*jj1-11)
              Akj(2,2)= AUlu0(16*jj1-10)
              Akj(2,3)= AUlu0(16*jj1- 9)
              Akj(2,4)= AUlu0(16*jj1- 8)
              Akj(3,1)= AUlu0(16*jj1- 7)
              Akj(3,2)= AUlu0(16*jj1- 6)
              Akj(3,3)= AUlu0(16*jj1- 5)
              Akj(3,4)= AUlu0(16*jj1- 4)
              Akj(4,1)= AUlu0(16*jj1- 3)
              Akj(4,2)= AUlu0(16*jj1- 2)
              Akj(4,3)= AUlu0(16*jj1- 1)
              Akj(4,4)= AUlu0(16*jj1- 0)
              exit
            endif
          enddo

          call ILU1b44 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(16*i-15)= Dlu0(16*i-15) - RHS_Aij(1,1)
            Dlu0(16*i-14)= Dlu0(16*i-14) - RHS_Aij(1,2)
            Dlu0(16*i-13)= Dlu0(16*i-13) - RHS_Aij(1,3)
            Dlu0(16*i-12)= Dlu0(16*i-12) - RHS_Aij(1,4)
            Dlu0(16*i-11)= Dlu0(16*i-11) - RHS_Aij(2,1)
            Dlu0(16*i-10)= Dlu0(16*i-10) - RHS_Aij(2,2)
            Dlu0(16*i- 9)= Dlu0(16*i- 9) - RHS_Aij(2,3)
            Dlu0(16*i- 8)= Dlu0(16*i- 8) - RHS_Aij(2,4)
            Dlu0(16*i- 7)= Dlu0(16*i- 7) - RHS_Aij(3,1)
            Dlu0(16*i- 6)= Dlu0(16*i- 6) - RHS_Aij(3,2)
            Dlu0(16*i- 5)= Dlu0(16*i- 5) - RHS_Aij(3,3)
            Dlu0(16*i- 4)= Dlu0(16*i- 4) - RHS_Aij(3,4)
            Dlu0(16*i- 3)= Dlu0(16*i- 3) - RHS_Aij(4,1)
            Dlu0(16*i- 2)= Dlu0(16*i- 2) - RHS_Aij(4,2)
            Dlu0(16*i- 1)= Dlu0(16*i- 1) - RHS_Aij(4,3)
            Dlu0(16*i- 0)= Dlu0(16*i- 0) - RHS_Aij(4,4)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(16*ij0-15)= ALlu0(16*ij0-15) - RHS_Aij(1,1)
            ALlu0(16*ij0-14)= ALlu0(16*ij0-14) - RHS_Aij(1,2)
            ALlu0(16*ij0-13)= ALlu0(16*ij0-13) - RHS_Aij(1,3)
            ALlu0(16*ij0-12)= ALlu0(16*ij0-12) - RHS_Aij(1,4)
            ALlu0(16*ij0-11)= ALlu0(16*ij0-11) - RHS_Aij(2,1)
            ALlu0(16*ij0-10)= ALlu0(16*ij0-10) - RHS_Aij(2,2)
            ALlu0(16*ij0- 9)= ALlu0(16*ij0- 9) - RHS_Aij(2,3)
            ALlu0(16*ij0- 8)= ALlu0(16*ij0- 8) - RHS_Aij(2,4)
            ALlu0(16*ij0- 7)= ALlu0(16*ij0- 7) - RHS_Aij(3,1)
            ALlu0(16*ij0- 6)= ALlu0(16*ij0- 6) - RHS_Aij(3,2)
            ALlu0(16*ij0- 5)= ALlu0(16*ij0- 5) - RHS_Aij(3,3)
            ALlu0(16*ij0- 4)= ALlu0(16*ij0- 4) - RHS_Aij(3,4)
            ALlu0(16*ij0- 3)= ALlu0(16*ij0- 3) - RHS_Aij(4,1)
            ALlu0(16*ij0- 2)= ALlu0(16*ij0- 2) - RHS_Aij(4,2)
            ALlu0(16*ij0- 1)= ALlu0(16*ij0- 1) - RHS_Aij(4,3)
            ALlu0(16*ij0- 0)= ALlu0(16*ij0- 0) - RHS_Aij(4,4)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(16*ij0-15)= AUlu0(16*ij0-15) - RHS_Aij(1,1)
            AUlu0(16*ij0-14)= AUlu0(16*ij0-14) - RHS_Aij(1,2)
            AUlu0(16*ij0-13)= AUlu0(16*ij0-13) - RHS_Aij(1,3)
            AUlu0(16*ij0-12)= AUlu0(16*ij0-12) - RHS_Aij(1,4)
            AUlu0(16*ij0-11)= AUlu0(16*ij0-11) - RHS_Aij(2,1)
            AUlu0(16*ij0-10)= AUlu0(16*ij0-10) - RHS_Aij(2,2)
            AUlu0(16*ij0- 9)= AUlu0(16*ij0- 9) - RHS_Aij(2,3)
            AUlu0(16*ij0- 8)= AUlu0(16*ij0- 8) - RHS_Aij(2,4)
            AUlu0(16*ij0- 7)= AUlu0(16*ij0- 7) - RHS_Aij(3,1)
            AUlu0(16*ij0- 6)= AUlu0(16*ij0- 6) - RHS_Aij(3,2)
            AUlu0(16*ij0- 5)= AUlu0(16*ij0- 5) - RHS_Aij(3,3)
            AUlu0(16*ij0- 4)= AUlu0(16*ij0- 4) - RHS_Aij(3,4)
            AUlu0(16*ij0- 3)= AUlu0(16*ij0- 3) - RHS_Aij(4,1)
            AUlu0(16*ij0- 2)= AUlu0(16*ij0- 2) - RHS_Aij(4,2)
            AUlu0(16*ij0- 1)= AUlu0(16*ij0- 1) - RHS_Aij(4,3)
            AUlu0(16*ij0- 0)= AUlu0(16*ij0- 0) - RHS_Aij(4,4)
          endif

        enddo
      enddo
    enddo

    deallocate (IW1, IW2)
    !C===
  end subroutine FORM_ILU1_44

  !C
  !C***
  !C*** FORM_ILU2_44
  !C***
  !C
  !C    form ILU(2) matrix
  !C
  subroutine FORM_ILU2_44 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(16*NPL), intent(in):: AL
    real(kind=kreal), dimension(16*NPU), intent(in):: AU
    real(kind=kreal), dimension(16*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable:: IW1 , IW2
    integer(kind=kint), dimension(:), allocatable:: IWsL, IWsU
    integer(kind=kint), dimension(:), allocatable:: iconFI1L, iconFI1U
    integer(kind=kint), dimension(:), allocatable:: inumFI2L, inumFI2U
    integer(kind=kint), dimension(:), allocatable::     FI2L,     FI2U
    real (kind=kreal), dimension(4,4) :: RHS_Aij, DkINV, Aik, Akj
    real (kind=kreal)  :: D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44
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
    allocate (ALlu0(16*(NPL+NPLf1+NPLf2)))
    allocate (AUlu0(16*(NPU+NPUf1+NPUf2)))

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

      call fill_in_S44_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 16*( k+isL)
          kk2= 16*(ik+INL(i-1))
          ALlu0(kk1-15)= AL(kk2-15)
          ALlu0(kk1-14)= AL(kk2-14)
          ALlu0(kk1-13)= AL(kk2-13)
          ALlu0(kk1-12)= AL(kk2-12)
          ALlu0(kk1-11)= AL(kk2-11)
          ALlu0(kk1-10)= AL(kk2-10)
          ALlu0(kk1- 9)= AL(kk2- 9)
          ALlu0(kk1- 8)= AL(kk2- 8)
          ALlu0(kk1- 7)= AL(kk2- 7)
          ALlu0(kk1- 6)= AL(kk2- 6)
          ALlu0(kk1- 5)= AL(kk2- 5)
          ALlu0(kk1- 4)= AL(kk2- 4)
          ALlu0(kk1- 3)= AL(kk2- 3)
          ALlu0(kk1- 2)= AL(kk2- 2)
          ALlu0(kk1- 1)= AL(kk2- 1)
          ALlu0(kk1- 0)= AL(kk2- 0)
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
      call fill_in_S44_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 16*( k+isU)
          kk2= 16*(ik+INU(i-1))
          AUlu0(kk1-15)= AU(kk2-15)
          AUlu0(kk1-14)= AU(kk2-14)
          AUlu0(kk1-13)= AU(kk2-13)
          AUlu0(kk1-12)= AU(kk2-12)
          AUlu0(kk1-11)= AU(kk2-11)
          AUlu0(kk1-10)= AU(kk2-10)
          AUlu0(kk1- 9)= AU(kk2- 9)
          AUlu0(kk1- 8)= AU(kk2- 8)
          AUlu0(kk1- 7)= AU(kk2- 7)
          AUlu0(kk1- 6)= AU(kk2- 6)
          AUlu0(kk1- 5)= AU(kk2- 5)
          AUlu0(kk1- 4)= AU(kk2- 4)
          AUlu0(kk1- 3)= AU(kk2- 3)
          AUlu0(kk1- 2)= AU(kk2- 2)
          AUlu0(kk1- 1)= AU(kk2- 1)
          AUlu0(kk1- 0)= AU(kk2- 0)
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
    allocate (Dlu0(16*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(16*i-15)=Dlu0(16*i-15)*SIGMA_DIAG
      Dlu0(16*i-10)=Dlu0(16*i-10)*SIGMA_DIAG
      Dlu0(16*i- 5)=Dlu0(16*i- 5)*SIGMA_DIAG
      Dlu0(16*i- 0)=Dlu0(16*i- 0)*SIGMA_DIAG
    enddo

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

        D11= Dlu0(16*k-15)
        D12= Dlu0(16*k-14)
        D13= Dlu0(16*k-13)
        D14= Dlu0(16*k-12)
        D21= Dlu0(16*k-11)
        D22= Dlu0(16*k-10)
        D23= Dlu0(16*k- 9)
        D24= Dlu0(16*k- 8)
        D31= Dlu0(16*k- 7)
        D32= Dlu0(16*k- 6)
        D33= Dlu0(16*k- 5)
        D34= Dlu0(16*k- 4)
        D41= Dlu0(16*k- 3)
        D42= Dlu0(16*k- 2)
        D43= Dlu0(16*k- 1)
        D44= Dlu0(16*k- 0)

        call ILU1a44 (DkINV, D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44)

        Aik(1,1)= ALlu0(16*kk1-15)
        Aik(1,2)= ALlu0(16*kk1-14)
        Aik(1,3)= ALlu0(16*kk1-13)
        Aik(1,4)= ALlu0(16*kk1-12)
        Aik(2,1)= ALlu0(16*kk1-11)
        Aik(2,2)= ALlu0(16*kk1-10)
        Aik(2,3)= ALlu0(16*kk1- 9)
        Aik(2,4)= ALlu0(16*kk1- 8)
        Aik(3,1)= ALlu0(16*kk1- 7)
        Aik(3,2)= ALlu0(16*kk1- 6)
        Aik(3,3)= ALlu0(16*kk1- 5)
        Aik(3,4)= ALlu0(16*kk1- 4)
        Aik(4,1)= ALlu0(16*kk1- 3)
        Aik(4,2)= ALlu0(16*kk1- 2)
        Aik(4,3)= ALlu0(16*kk1- 1)
        Aik(4,4)= ALlu0(16*kk1- 0)

        do jj= inumFI1U(k-1)+1, inumFI1U(k)
          j= FI1U(jj)
          iconKJ= iconFI1U(jj)

          if ((iconIK+iconKJ).lt.2) then
            Akj(1,1)= AUlu0(16*jj-15)
            Akj(1,2)= AUlu0(16*jj-14)
            Akj(1,3)= AUlu0(16*jj-13)
            Akj(1,4)= AUlu0(16*jj-12)
            Akj(2,1)= AUlu0(16*jj-11)
            Akj(2,2)= AUlu0(16*jj-10)
            Akj(2,3)= AUlu0(16*jj- 9)
            Akj(2,4)= AUlu0(16*jj- 8)
            Akj(3,1)= AUlu0(16*jj- 7)
            Akj(3,2)= AUlu0(16*jj- 6)
            Akj(3,3)= AUlu0(16*jj- 5)
            Akj(3,4)= AUlu0(16*jj- 4)
            Akj(4,1)= AUlu0(16*jj- 3)
            Akj(4,2)= AUlu0(16*jj- 2)
            Akj(4,3)= AUlu0(16*jj- 1)
            Akj(4,4)= AUlu0(16*jj- 0)

            call ILU1b44 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(16*i-15)= Dlu0(16*i-15) - RHS_Aij(1,1)
              Dlu0(16*i-14)= Dlu0(16*i-14) - RHS_Aij(1,2)
              Dlu0(16*i-13)= Dlu0(16*i-13) - RHS_Aij(1,3)
              Dlu0(16*i-12)= Dlu0(16*i-12) - RHS_Aij(1,4)
              Dlu0(16*i-11)= Dlu0(16*i-11) - RHS_Aij(2,1)
              Dlu0(16*i-10)= Dlu0(16*i-10) - RHS_Aij(2,2)
              Dlu0(16*i- 9)= Dlu0(16*i- 9) - RHS_Aij(2,3)
              Dlu0(16*i- 8)= Dlu0(16*i- 8) - RHS_Aij(2,4)
              Dlu0(16*i- 7)= Dlu0(16*i- 7) - RHS_Aij(3,1)
              Dlu0(16*i- 6)= Dlu0(16*i- 6) - RHS_Aij(3,2)
              Dlu0(16*i- 5)= Dlu0(16*i- 5) - RHS_Aij(3,3)
              Dlu0(16*i- 4)= Dlu0(16*i- 4) - RHS_Aij(3,4)
              Dlu0(16*i- 3)= Dlu0(16*i- 3) - RHS_Aij(4,1)
              Dlu0(16*i- 2)= Dlu0(16*i- 2) - RHS_Aij(4,2)
              Dlu0(16*i- 1)= Dlu0(16*i- 1) - RHS_Aij(4,3)
              Dlu0(16*i- 0)= Dlu0(16*i- 0) - RHS_Aij(4,4)
            endif

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(16*ij0-15)= ALlu0(16*ij0-15) - RHS_Aij(1,1)
              ALlu0(16*ij0-14)= ALlu0(16*ij0-14) - RHS_Aij(1,2)
              ALlu0(16*ij0-13)= ALlu0(16*ij0-13) - RHS_Aij(1,3)
              ALlu0(16*ij0-12)= ALlu0(16*ij0-12) - RHS_Aij(1,4)
              ALlu0(16*ij0-11)= ALlu0(16*ij0-11) - RHS_Aij(2,1)
              ALlu0(16*ij0-10)= ALlu0(16*ij0-10) - RHS_Aij(2,2)
              ALlu0(16*ij0- 9)= ALlu0(16*ij0- 9) - RHS_Aij(2,3)
              ALlu0(16*ij0- 8)= ALlu0(16*ij0- 8) - RHS_Aij(2,4)
              ALlu0(16*ij0- 7)= ALlu0(16*ij0- 7) - RHS_Aij(3,1)
              ALlu0(16*ij0- 6)= ALlu0(16*ij0- 6) - RHS_Aij(3,2)
              ALlu0(16*ij0- 5)= ALlu0(16*ij0- 5) - RHS_Aij(3,3)
              ALlu0(16*ij0- 4)= ALlu0(16*ij0- 4) - RHS_Aij(3,4)
              ALlu0(16*ij0- 3)= ALlu0(16*ij0- 3) - RHS_Aij(4,1)
              ALlu0(16*ij0- 2)= ALlu0(16*ij0- 2) - RHS_Aij(4,2)
              ALlu0(16*ij0- 1)= ALlu0(16*ij0- 1) - RHS_Aij(4,3)
              ALlu0(16*ij0- 0)= ALlu0(16*ij0- 0) - RHS_Aij(4,4)
            endif

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(16*ij0-15)= AUlu0(16*ij0-15) - RHS_Aij(1,1)
              AUlu0(16*ij0-14)= AUlu0(16*ij0-14) - RHS_Aij(1,2)
              AUlu0(16*ij0-13)= AUlu0(16*ij0-13) - RHS_Aij(1,3)
              AUlu0(16*ij0-12)= AUlu0(16*ij0-12) - RHS_Aij(1,4)
              AUlu0(16*ij0-11)= AUlu0(16*ij0-11) - RHS_Aij(2,1)
              AUlu0(16*ij0-10)= AUlu0(16*ij0-10) - RHS_Aij(2,2)
              AUlu0(16*ij0- 9)= AUlu0(16*ij0- 9) - RHS_Aij(2,3)
              AUlu0(16*ij0- 8)= AUlu0(16*ij0- 8) - RHS_Aij(2,4)
              AUlu0(16*ij0- 7)= AUlu0(16*ij0- 7) - RHS_Aij(3,1)
              AUlu0(16*ij0- 6)= AUlu0(16*ij0- 6) - RHS_Aij(3,2)
              AUlu0(16*ij0- 5)= AUlu0(16*ij0- 5) - RHS_Aij(3,3)
              AUlu0(16*ij0- 4)= AUlu0(16*ij0- 4) - RHS_Aij(3,4)
              AUlu0(16*ij0- 3)= AUlu0(16*ij0- 3) - RHS_Aij(4,1)
              AUlu0(16*ij0- 2)= AUlu0(16*ij0- 2) - RHS_Aij(4,2)
              AUlu0(16*ij0- 1)= AUlu0(16*ij0- 1) - RHS_Aij(4,3)
              AUlu0(16*ij0- 0)= AUlu0(16*ij0- 0) - RHS_Aij(4,4)
            endif
          endif
        enddo
      enddo
    enddo

    deallocate (IW1, IW2)
    deallocate (iconFI1L, iconFI1U)
    !C===
  end subroutine FORM_ILU2_44


  !C
  !C***
  !C*** fill_in_S44_SORT
  !C***
  !C
  subroutine fill_in_S44_SORT (STEM, INUM, N, NP)
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

  end subroutine fill_in_S44_SORT

  !C
  !C***
  !C*** ILU1a44
  !C***
  !C
  !C    computes LU factorization of 4*4 Diagonal Block
  !C
  subroutine ILU1a44 (ALU, D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44)
    use hecmw_util
    implicit none
    real(kind=kreal) :: ALU(4,4), PW(4)
    real(kind=kreal) :: D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44
    integer(kind=kint) :: i,j,k

    ALU(1,1)= D11
    ALU(1,2)= D12
    ALU(1,3)= D13
    ALU(1,4)= D14
    ALU(2,1)= D21
    ALU(2,2)= D22
    ALU(2,3)= D23
    ALU(2,4)= D24
    ALU(3,1)= D31
    ALU(3,2)= D32
    ALU(3,3)= D33
    ALU(3,4)= D34
    ALU(4,1)= D41
    ALU(4,2)= D42
    ALU(4,3)= D43
    ALU(4,4)= D44

    do k= 1, 4
      ALU(k,k)= 1.d0/ALU(k,k)
      do i= k+1, 4
        ALU(i,k)= ALU(i,k) * ALU(k,k)
        do j= k+1, 4
          PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
        enddo
        do j= k+1, 4
          ALU(i,j)= PW(j)
        enddo
      enddo
    enddo

    return
  end subroutine ILU1a44

  !C
  !C***
  !C*** ILU1b44
  !C***
  !C
  !C    computes L_ik * D_k_INV * U_kj at ILU factorization
  !C    for 4*4 Block Type Matrix
  !C
  subroutine ILU1b44 (RHS_Aij, DkINV, Aik, Akj)
    use hecmw_util
    implicit none
    real(kind=kreal) :: RHS_Aij(4,4), DkINV(4,4), Aik(4,4), Akj(4,4)
    real(kind=kreal) :: X1,X2,X3,X4

    !C
    !C-- 1st Col.
    X1= Akj(1,1)
    X2= Akj(2,1)
    X3= Akj(3,1)
    X4= Akj(4,1)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2
    X4= X4 - DkINV(4,1)*X1 - DkINV(4,2)*X2 - DkINV(4,3)*X3

    X4= DkINV(4,4)*  X4
    X3= DkINV(3,3)*( X3 - DkINV(3,4)*X4 )
    X2= DkINV(2,2)*( X2 - DkINV(2,4)*X4 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,4)*X4 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,1)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3 + Aik(1,4)*X4
    RHS_Aij(2,1)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3 + Aik(2,4)*X4
    RHS_Aij(3,1)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3 + Aik(3,4)*X4
    RHS_Aij(4,1)=  Aik(4,1)*X1 + Aik(4,2)*X2 + Aik(4,3)*X3 + Aik(4,4)*X4

    !C
    !C-- 2nd Col.
    X1= Akj(1,2)
    X2= Akj(2,2)
    X3= Akj(3,2)
    X4= Akj(4,2)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2
    X4= X4 - DkINV(4,1)*X1 - DkINV(4,2)*X2 - DkINV(4,3)*X3

    X4= DkINV(4,4)*  X4
    X3= DkINV(3,3)*( X3 - DkINV(3,4)*X4 )
    X2= DkINV(2,2)*( X2 - DkINV(2,4)*X4 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,4)*X4 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,2)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3 + Aik(1,4)*X4
    RHS_Aij(2,2)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3 + Aik(2,4)*X4
    RHS_Aij(3,2)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3 + Aik(3,4)*X4
    RHS_Aij(4,2)=  Aik(4,1)*X1 + Aik(4,2)*X2 + Aik(4,3)*X3 + Aik(4,4)*X4

    !C
    !C-- 3rd Col.
    X1= Akj(1,3)
    X2= Akj(2,3)
    X3= Akj(3,3)
    X4= Akj(4,3)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2
    X4= X4 - DkINV(4,1)*X1 - DkINV(4,2)*X2 - DkINV(4,3)*X3

    X4= DkINV(4,4)*  X4
    X3= DkINV(3,3)*( X3 - DkINV(3,4)*X4 )
    X2= DkINV(2,2)*( X2 - DkINV(2,4)*X4 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,4)*X4 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,3)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3 + Aik(1,4)*X4
    RHS_Aij(2,3)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3 + Aik(2,4)*X4
    RHS_Aij(3,3)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3 + Aik(3,4)*X4
    RHS_Aij(4,3)=  Aik(4,1)*X1 + Aik(4,2)*X2 + Aik(4,3)*X3 + Aik(4,4)*X4

    !C
    !C-- 4th Col.
    X1= Akj(1,4)
    X2= Akj(2,4)
    X3= Akj(3,4)
    X4= Akj(4,4)

    X2= X2 - DkINV(2,1)*X1
    X3= X3 - DkINV(3,1)*X1 - DkINV(3,2)*X2
    X4= X4 - DkINV(4,1)*X1 - DkINV(4,2)*X2 - DkINV(4,3)*X3

    X4= DkINV(4,4)*  X4
    X3= DkINV(3,3)*( X3 - DkINV(3,4)*X4 )
    X2= DkINV(2,2)*( X2 - DkINV(2,4)*X4 - DkINV(2,3)*X3 )
    X1= DkINV(1,1)*( X1 - DkINV(1,4)*X4 - DkINV(1,3)*X3 - DkINV(1,2)*X2)

    RHS_Aij(1,4)=  Aik(1,1)*X1 + Aik(1,2)*X2 + Aik(1,3)*X3 + Aik(1,4)*X4
    RHS_Aij(2,4)=  Aik(2,1)*X1 + Aik(2,2)*X2 + Aik(2,3)*X3 + Aik(2,4)*X4
    RHS_Aij(3,4)=  Aik(3,1)*X1 + Aik(3,2)*X2 + Aik(3,3)*X3 + Aik(3,4)*X4
    RHS_Aij(4,4)=  Aik(4,1)*X1 + Aik(4,2)*X2 + Aik(4,3)*X3 + Aik(4,4)*X4

    return
  end subroutine ILU1b44

end module     hecmw_precond_BILU_44

