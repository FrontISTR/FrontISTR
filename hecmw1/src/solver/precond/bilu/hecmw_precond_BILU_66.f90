!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_BILU_66
!C***
!C
module hecmw_precond_BILU_66
  use hecmw_util
  use hecmw_matrix_misc

  private

  public:: hecmw_precond_BILU_66_setup
  public:: hecmw_precond_BILU_66_apply
  public:: hecmw_precond_BILU_66_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: Dlu0(:) => null()
  real(kind=kreal), pointer :: ALlu0(:) => null()
  real(kind=kreal), pointer :: AUlu0(:) => null()
  integer(kind=kint), pointer :: inumFI1L(:) => null()
  integer(kind=kint), pointer :: inumFI1U(:) => null()
  integer(kind=kint), pointer :: FI1L(:) => null()
  integer(kind=kint), pointer :: FI1U(:) => null()
  real(kind=kreal), pointer :: ALU(:) => null()

contains

  subroutine hecmw_precond_BILU_66_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint ) :: NP, NPU, NPL
    integer(kind=kint ) :: PRECOND
    real   (kind=kreal) :: SIGMA, SIGMA_DIAG

    real(kind=kreal), pointer :: D(:)
    real(kind=kreal), pointer :: AL(:)
    real(kind=kreal), pointer :: AU(:)

    integer(kind=kint ), pointer :: INL(:), INU(:)
    integer(kind=kint ), pointer :: IAL(:)
    integer(kind=kint ), pointer :: IAU(:)

    real   (kind=kreal):: ALUtmp(6,6)
    integer(kind=kint ):: ip

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

    if (PRECOND.eq.10) call FORM_ILU0_66 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.11) call FORM_ILU1_66 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.12) call FORM_ILU2_66 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)

    allocate(ALU(36*NP))

    do ip= 1, N
      call ILU1a66 (ALUtmp,                                    &
        & Dlu0(36*ip-35), Dlu0(36*ip-34), Dlu0(36*ip-33), Dlu0(36*ip-32),&
        & Dlu0(36*ip-31), Dlu0(36*ip-30), Dlu0(36*ip-29), Dlu0(36*ip-28),&
        & Dlu0(36*ip-27), Dlu0(36*ip-26), Dlu0(36*ip-25), Dlu0(36*ip-24),&
        & Dlu0(36*ip-23), Dlu0(36*ip-22), Dlu0(36*ip-21), Dlu0(36*ip-20),&
        & Dlu0(36*ip-19), Dlu0(36*ip-18), Dlu0(36*ip-17), Dlu0(36*ip-16),&
        & Dlu0(36*ip-15), Dlu0(36*ip-14), Dlu0(36*ip-13), Dlu0(36*ip-12),&
        & Dlu0(36*ip-11), Dlu0(36*ip-10), Dlu0(36*ip-9 ), Dlu0(36*ip-8 ),&
        & Dlu0(36*ip-7), Dlu0(36*ip-6), Dlu0(36*ip-5), Dlu0(36*ip-4),&
        & Dlu0(36*ip-3), Dlu0(36*ip-2), Dlu0(36*ip-1), Dlu0(36*ip  ))
      ALU(36*ip-35)= ALUtmp(1,1)
      ALU(36*ip-34)= ALUtmp(1,2)
      ALU(36*ip-33)= ALUtmp(1,3)
      ALU(36*ip-32)= ALUtmp(1,4)
      ALU(36*ip-31)= ALUtmp(1,5)
      ALU(36*ip-30)= ALUtmp(1,6)
      ALU(36*ip-29)= ALUtmp(2,1)
      ALU(36*ip-28)= ALUtmp(2,2)
      ALU(36*ip-27)= ALUtmp(2,3)
      ALU(36*ip-26)= ALUtmp(2,4)
      ALU(36*ip-25)= ALUtmp(2,5)
      ALU(36*ip-24)= ALUtmp(2,6)
      ALU(36*ip-23)= ALUtmp(3,1)
      ALU(36*ip-22)= ALUtmp(3,2)
      ALU(36*ip-21)= ALUtmp(3,3)
      ALU(36*ip-20)= ALUtmp(3,4)
      ALU(36*ip-19)= ALUtmp(3,5)
      ALU(36*ip-18)= ALUtmp(3,6)
      ALU(36*ip-17)= ALUtmp(4,1)
      ALU(36*ip-16)= ALUtmp(4,2)
      ALU(36*ip-15)= ALUtmp(4,3)
      ALU(36*ip-14)= ALUtmp(4,4)
      ALU(36*ip-13)= ALUtmp(4,5)
      ALU(36*ip-12)= ALUtmp(4,6)
      ALU(36*ip-11)= ALUtmp(5,1)
      ALU(36*ip-10)= ALUtmp(5,2)
      ALU(36*ip-9 )= ALUtmp(5,3)
      ALU(36*ip-8 )= ALUtmp(5,4)
      ALU(36*ip-7 )= ALUtmp(5,5)
      ALU(36*ip-6 )= ALUtmp(5,6)
      ALU(36*ip-5 )= ALUtmp(6,1)
      ALU(36*ip-4 )= ALUtmp(6,2)
      ALU(36*ip-3 )= ALUtmp(6,3)
      ALU(36*ip-2 )= ALUtmp(6,4)
      ALU(36*ip-1 )= ALUtmp(6,5)
      ALU(36*ip   )= ALUtmp(6,6)
    enddo
  end subroutine hecmw_precond_BILU_66_setup

  subroutine hecmw_precond_BILU_66_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3
    real(kind=kreal) :: SW4, SW5, SW6, X4, X5, X6
    !C
    !C-- FORWARD

    do i= 1, N
      SW1= WW(6*i-5)
      SW2= WW(6*i-4)
      SW3= WW(6*i-3)
      SW4= WW(6*i-2)
      SW5= WW(6*i-1)
      SW6= WW(6*i  )
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      do j= isL, ieL
        k= FI1L(j)
        X1= WW(6*k-5)
        X2= WW(6*k-4)
        X3= WW(6*k-3)
        X4= WW(6*k-2)
        X5= WW(6*k-1)
        X6= WW(6*k  )
        SW1= SW1 -ALlu0(36*j-35)*X1-ALlu0(36*j-34)*X2-ALlu0(36*j-33)*X3-ALlu0(36*j-32)*X4-ALlu0(36*j-31)*X5-ALlu0(36*j-30)*X6
        SW2= SW2 -ALlu0(36*j-29)*X1-ALlu0(36*j-28)*X2-ALlu0(36*j-27)*X3-ALlu0(36*j-26)*X4-ALlu0(36*j-25)*X5-ALlu0(36*j-24)*X6
        SW3= SW3 -ALlu0(36*j-23)*X1-ALlu0(36*j-22)*X2-ALlu0(36*j-21)*X3-ALlu0(36*j-20)*X4-ALlu0(36*j-19)*X5-ALlu0(36*j-18)*X6
        SW4= SW4 -ALlu0(36*j-17)*X1-ALlu0(36*j-16)*X2-ALlu0(36*j-15)*X3-ALlu0(36*j-14)*X4-ALlu0(36*j-13)*X5-ALlu0(36*j-12)*X6
        SW5= SW5 -ALlu0(36*j-11)*X1-ALlu0(36*j-10)*X2-ALlu0(36*j-9 )*X3-ALlu0(36*j-8 )*X4-ALlu0(36*j-7 )*X5-ALlu0(36*j-6 )*X6
        SW6= SW6 -ALlu0(36*j-5 )*X1-ALlu0(36*j-4 )*X2-ALlu0(36*j-3 )*X3-ALlu0(36*j-2 )*X4-ALlu0(36*j-1 )*X5-ALlu0(36*j   )*X6
      enddo

      X1= SW1
      X2= SW2
      X3= SW3
      X4= SW4
      X5= SW5
      X6= SW6
      X2= X2 -ALU(36*i-29)*X1
      X3= X3 -ALU(36*i-23)*X1 -ALU(36*i-22)*X2
      X4= X4 -ALU(36*i-17)*X1 -ALU(36*i-16)*X2 -ALU(36*i-15)*X3
      X5= X5 -ALU(36*i-11)*X1 -ALU(36*i-10)*X2 -ALU(36*i-9 )*X3 -ALU(36*i-8)*X4
      X6= X6 -ALU(36*i-5 )*X1 -ALU(36*i-4 )*X2 -ALU(36*i-3 )*X3 -ALU(36*i-2)*X4 -ALU(36*i-1)*X5
      X6= ALU(36*i   )*  X6
      X5= ALU(36*i-7 )*( X5 -ALU(36*i-6 )*X6 )
      X4= ALU(36*i-14)*( X4 -ALU(36*i-12)*X6 -ALU(36*i-13)*X5)
      X3= ALU(36*i-21)*( X3 -ALU(36*i-18)*X6 -ALU(36*i-19)*X5 -ALU(36*i-20)*X4)
      X2= ALU(36*i-28)*( X2 -ALU(36*i-24)*X6 -ALU(36*i-25)*X5 -ALU(36*i-26)*X4 -ALU(36*i-27)*X3)
      X1= ALU(36*i-35)*( X1 -ALU(36*i-30)*X6 -ALU(36*i-31)*X5 -ALU(36*i-32)*X4 -ALU(36*i-33)*X3 -ALU(36*i-34)*X2)
      WW(6*i-5)= X1
      WW(6*i-4)= X2
      WW(6*i-3)= X3
      WW(6*i-2)= X4
      WW(6*i-1)= X5
      WW(6*i  )= X6
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
      SW5= 0.d0
      SW6= 0.d0
      do j= ieU, isU, -1
        k= FI1U(j)
        X1= WW(6*k-5)
        X2= WW(6*k-4)
        X3= WW(6*k-3)
        X4= WW(6*k-2)
        X5= WW(6*k-1)
        X6= WW(6*k  )
        SW1= SW1 +AUlu0(36*j-35)*X1+AUlu0(36*j-34)*X2+AUlu0(36*j-33)*X3+AUlu0(36*j-32)*X4+AUlu0(36*j-31)*X5+AUlu0(36*j-30)*X6
        SW2= SW2 +AUlu0(36*j-29)*X1+AUlu0(36*j-28)*X2+AUlu0(36*j-27)*X3+AUlu0(36*j-26)*X4+AUlu0(36*j-25)*X5+AUlu0(36*j-24)*X6
        SW3= SW3 +AUlu0(36*j-23)*X1+AUlu0(36*j-22)*X2+AUlu0(36*j-21)*X3+AUlu0(36*j-20)*X4+AUlu0(36*j-19)*X5+AUlu0(36*j-18)*X6
        SW4= SW4 +AUlu0(36*j-17)*X1+AUlu0(36*j-16)*X2+AUlu0(36*j-15)*X3+AUlu0(36*j-14)*X4+AUlu0(36*j-13)*X5+AUlu0(36*j-12)*X6
        SW5= SW5 +AUlu0(36*j-11)*X1+AUlu0(36*j-10)*X2+AUlu0(36*j-9 )*X3+AUlu0(36*j-8 )*X4+AUlu0(36*j-7 )*X5+AUlu0(36*j-6 )*X6
        SW6= SW6 +AUlu0(36*j-5 )*X1+AUlu0(36*j-4 )*X2+AUlu0(36*j-3 )*X3+AUlu0(36*j-2 )*X4+AUlu0(36*j-1 )*X5+AUlu0(36*j   )*X6
      enddo
      X1= SW1
      X2= SW2
      X3= SW3
      X4= SW4
      X5= SW5
      X6= SW6
      X2= X2 -ALU(36*i-29)*X1
      X3= X3 -ALU(36*i-23)*X1 -ALU(36*i-22)*X2
      X4= X4 -ALU(36*i-17)*X1 -ALU(36*i-16)*X2 -ALU(36*i-15)*X3
      X5= X5 -ALU(36*i-11)*X1 -ALU(36*i-10)*X2 -ALU(36*i-9 )*X3 -ALU(36*i-8)*X4
      X6= X6 -ALU(36*i-5 )*X1 -ALU(36*i-4 )*X2 -ALU(36*i-3 )*X3 -ALU(36*i-2)*X4 -ALU(36*i-1)*X5
      X6= ALU(36*i   )*  X6
      X5= ALU(36*i-7 )*( X5 -ALU(36*i-6 )*X6 )
      X4= ALU(36*i-14)*( X4 -ALU(36*i-12)*X6 -ALU(36*i-13)*X5)
      X3= ALU(36*i-21)*( X3 -ALU(36*i-18)*X6 -ALU(36*i-19)*X5 -ALU(36*i-20)*X4)
      X2= ALU(36*i-28)*( X2 -ALU(36*i-24)*X6 -ALU(36*i-25)*X5 -ALU(36*i-26)*X4 -ALU(36*i-27)*X3)
      X1= ALU(36*i-35)*( X1 -ALU(36*i-30)*X6 -ALU(36*i-31)*X5 -ALU(36*i-32)*X4 -ALU(36*i-33)*X3 -ALU(36*i-34)*X2)
      WW(6*i-5)=  WW(6*i-5) -X1
      WW(6*i-4)=  WW(6*i-4) -X2
      WW(6*i-3)=  WW(6*i-3) -X3
      WW(6*i-2)=  WW(6*i-2) -X4
      WW(6*i-1)=  WW(6*i-1) -X5
      WW(6*i  )=  WW(6*i  ) -X6
    enddo
  end subroutine hecmw_precond_BILU_66_apply

  subroutine hecmw_precond_BILU_66_clear()
    implicit none
    if (associated(Dlu0)) deallocate(Dlu0)
    if (associated(ALlu0)) deallocate(ALlu0)
    if (associated(AUlu0)) deallocate(AUlu0)
    if (associated(inumFI1L)) deallocate(inumFI1L)
    if (associated(inumFI1U)) deallocate(inumFI1U)
    if (associated(FI1L)) deallocate(FI1L)
    if (associated(FI1U)) deallocate(FI1U)
    if (associated(ALU)) deallocate(ALU)
    nullify(Dlu0)
    nullify(ALlu0)
    nullify(AUlu0)
    nullify(inumFI1L)
    nullify(inumFI1U)
    nullify(FI1L)
    nullify(FI1U)
    nullify(ALU)
  end subroutine hecmw_precond_BILU_66_clear

  !C
  !C***
  !C*** FORM_ILU1_66
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU0_66                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(36*NPL), intent(in):: AL
    real(kind=kreal), dimension(36*NPU), intent(in):: AU
    real(kind=kreal), dimension(36*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    real (kind=kreal),  dimension(6,6) :: RHS_Aij, DkINV, Aik, Akj
    real(kind=kreal) :: D11,D12,D13,D14,D15,D16
    real(kind=kreal) :: D21,D22,D23,D24,D25,D26
    real(kind=kreal) :: D31,D32,D33,D34,D35,D36
    real(kind=kreal) :: D41,D42,D43,D44,D45,D46
    real(kind=kreal) :: D51,D52,D53,D54,D55,D56
    real(kind=kreal) :: D61,D62,D63,D64,D65,D66
    integer(kind=kint) :: i,jj,ij0,kk
    integer(kind=kint) :: j,k
    allocate (IW1(NP) , IW2(NP))
    allocate(Dlu0(36*NP), ALlu0(36*NPL), AUlu0(36*NPU))
    allocate(inumFI1L(0:NP), inumFI1U(0:NP), FI1L(NPL), FI1U(NPU))

    do i=1,36*NP
      Dlu0(i) = D(i)
    end do
    do i=1,36*NPL
      ALlu0(i) = AL(i)
    end do
    do i=1,36*NPU
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
      Dlu0(36*i-35)=Dlu0(36*i-35)*SIGMA_DIAG
      Dlu0(36*i-28)=Dlu0(36*i-28)*SIGMA_DIAG
      Dlu0(36*i-21)=Dlu0(36*i-21)*SIGMA_DIAG
      Dlu0(36*i-14)=Dlu0(36*i-14)*SIGMA_DIAG
      Dlu0(36*i-7 )=Dlu0(36*i-7 )*SIGMA_DIAG
      Dlu0(36*i   )=Dlu0(36*i   )*SIGMA_DIAG
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
        D11= Dlu0(36*k-35)
        D12= Dlu0(36*k-34)
        D13= Dlu0(36*k-33)
        D14= Dlu0(36*k-32)
        D15= Dlu0(36*k-31)
        D16= Dlu0(36*k-30)
        D21= Dlu0(36*k-29)
        D22= Dlu0(36*k-28)
        D23= Dlu0(36*k-27)
        D24= Dlu0(36*k-26)
        D25= Dlu0(36*k-25)
        D26= Dlu0(36*k-24)
        D31= Dlu0(36*k-23)
        D32= Dlu0(36*k-22)
        D33= Dlu0(36*k-21)
        D34= Dlu0(36*k-20)
        D35= Dlu0(36*k-19)
        D36= Dlu0(36*k-18)
        D41= Dlu0(36*k-17)
        D42= Dlu0(36*k-16)
        D43= Dlu0(36*k-15)
        D44= Dlu0(36*k-14)
        D45= Dlu0(36*k-13)
        D46= Dlu0(36*k-12)
        D51= Dlu0(36*k-11)
        D52= Dlu0(36*k-10)
        D53= Dlu0(36*k-9 )
        D54= Dlu0(36*k-8 )
        D55= Dlu0(36*k-7 )
        D56= Dlu0(36*k-6 )
        D61= Dlu0(36*k-5 )
        D62= Dlu0(36*k-4 )
        D63= Dlu0(36*k-3 )
        D64= Dlu0(36*k-2 )
        D65= Dlu0(36*k-1 )
        D66= Dlu0(36*k   )

        call ILU1a66 (DkINV,D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26, &
          & D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56, &
          & D61,D62,D63,D64,D65,D66)

        Aik(1,1)= ALlu0(36*kk-35)
        Aik(1,2)= ALlu0(36*kk-34)
        Aik(1,3)= ALlu0(36*kk-33)
        Aik(1,4)= ALlu0(36*kk-32)
        Aik(1,5)= ALlu0(36*kk-31)
        Aik(1,6)= ALlu0(36*kk-30)
        Aik(2,1)= ALlu0(36*kk-29)
        Aik(2,2)= ALlu0(36*kk-28)
        Aik(2,3)= ALlu0(36*kk-27)
        Aik(2,4)= ALlu0(36*kk-26)
        Aik(2,5)= ALlu0(36*kk-25)
        Aik(2,6)= ALlu0(36*kk-24)
        Aik(3,1)= ALlu0(36*kk-23)
        Aik(3,2)= ALlu0(36*kk-22)
        Aik(3,3)= ALlu0(36*kk-21)
        Aik(3,4)= ALlu0(36*kk-20)
        Aik(3,5)= ALlu0(36*kk-19)
        Aik(3,6)= ALlu0(36*kk-18)
        Aik(4,1)= ALlu0(36*kk-17)
        Aik(4,2)= ALlu0(36*kk-16)
        Aik(4,3)= ALlu0(36*kk-15)
        Aik(4,4)= ALlu0(36*kk-14)
        Aik(4,5)= ALlu0(36*kk-13)
        Aik(4,6)= ALlu0(36*kk-12)
        Aik(5,1)= ALlu0(36*kk-11)
        Aik(5,2)= ALlu0(36*kk-10)
        Aik(5,3)= ALlu0(36*kk-9)
        Aik(5,4)= ALlu0(36*kk-8)
        Aik(5,5)= ALlu0(36*kk-7)
        Aik(5,6)= ALlu0(36*kk-6)
        Aik(6,1)= ALlu0(36*kk-5)
        Aik(6,2)= ALlu0(36*kk-4)
        Aik(6,3)= ALlu0(36*kk-3)
        Aik(6,4)= ALlu0(36*kk-2)
        Aik(6,5)= ALlu0(36*kk-1)
        Aik(6,6)= ALlu0(36*kk  )

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          if (IW1(j).eq.0.and.IW2(j).eq.0) cycle

          Akj(1,1)= AUlu0(36*jj-35)
          Akj(1,2)= AUlu0(36*jj-34)
          Akj(1,3)= AUlu0(36*jj-33)
          Akj(1,4)= AUlu0(36*jj-32)
          Akj(1,5)= AUlu0(36*jj-31)
          Akj(1,6)= AUlu0(36*jj-30)
          Akj(2,1)= AUlu0(36*jj-29)
          Akj(2,2)= AUlu0(36*jj-28)
          Akj(2,3)= AUlu0(36*jj-27)
          Akj(2,4)= AUlu0(36*jj-26)
          Akj(2,5)= AUlu0(36*jj-25)
          Akj(2,6)= AUlu0(36*jj-24)
          Akj(3,1)= AUlu0(36*jj-23)
          Akj(3,2)= AUlu0(36*jj-22)
          Akj(3,3)= AUlu0(36*jj-21)
          Akj(3,4)= AUlu0(36*jj-20)
          Akj(3,5)= AUlu0(36*jj-19)
          Akj(3,6)= AUlu0(36*jj-18)
          Akj(4,1)= AUlu0(36*jj-17)
          Akj(4,2)= AUlu0(36*jj-16)
          Akj(4,3)= AUlu0(36*jj-15)
          Akj(4,4)= AUlu0(36*jj-14)
          Akj(4,5)= AUlu0(36*jj-13)
          Akj(4,6)= AUlu0(36*jj-12)
          Akj(5,1)= AUlu0(36*jj-11)
          Akj(5,2)= AUlu0(36*jj-10)
          Akj(5,3)= AUlu0(36*jj-9)
          Akj(5,4)= AUlu0(36*jj-8)
          Akj(5,5)= AUlu0(36*jj-7)
          Akj(5,6)= AUlu0(36*jj-6)
          Akj(6,1)= AUlu0(36*jj-5)
          Akj(6,2)= AUlu0(36*jj-4)
          Akj(6,3)= AUlu0(36*jj-3)
          Akj(6,4)= AUlu0(36*jj-2)
          Akj(6,5)= AUlu0(36*jj-1)
          Akj(6,6)= AUlu0(36*jj  )

          call ILU1b66 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(36*i-35)= Dlu0(36*i-35) - RHS_Aij(1,1)
            Dlu0(36*i-34)= Dlu0(36*i-34) - RHS_Aij(1,2)
            Dlu0(36*i-33)= Dlu0(36*i-33) - RHS_Aij(1,3)
            Dlu0(36*i-32)= Dlu0(36*i-32) - RHS_Aij(1,4)
            Dlu0(36*i-31)= Dlu0(36*i-31) - RHS_Aij(1,5)
            Dlu0(36*i-30)= Dlu0(36*i-30) - RHS_Aij(1,6)
            Dlu0(36*i-29)= Dlu0(36*i-29) - RHS_Aij(2,1)
            Dlu0(36*i-28)= Dlu0(36*i-28) - RHS_Aij(2,2)
            Dlu0(36*i-27)= Dlu0(36*i-27) - RHS_Aij(2,3)
            Dlu0(36*i-26)= Dlu0(36*i-26) - RHS_Aij(2,4)
            Dlu0(36*i-25)= Dlu0(36*i-25) - RHS_Aij(2,5)
            Dlu0(36*i-24)= Dlu0(36*i-24) - RHS_Aij(2,6)
            Dlu0(36*i-23)= Dlu0(36*i-23) - RHS_Aij(3,1)
            Dlu0(36*i-22)= Dlu0(36*i-22) - RHS_Aij(3,2)
            Dlu0(36*i-21)= Dlu0(36*i-21) - RHS_Aij(3,3)
            Dlu0(36*i-20)= Dlu0(36*i-20) - RHS_Aij(3,4)
            Dlu0(36*i-19)= Dlu0(36*i-19) - RHS_Aij(3,5)
            Dlu0(36*i-18)= Dlu0(36*i-18) - RHS_Aij(3,6)
            Dlu0(36*i-17)= Dlu0(36*i-17) - RHS_Aij(4,1)
            Dlu0(36*i-16)= Dlu0(36*i-16) - RHS_Aij(4,2)
            Dlu0(36*i-15)= Dlu0(36*i-15) - RHS_Aij(4,3)
            Dlu0(36*i-14)= Dlu0(36*i-14) - RHS_Aij(4,4)
            Dlu0(36*i-13)= Dlu0(36*i-13) - RHS_Aij(4,5)
            Dlu0(36*i-12)= Dlu0(36*i-12) - RHS_Aij(4,6)
            Dlu0(36*i-11)= Dlu0(36*i-11) - RHS_Aij(5,1)
            Dlu0(36*i-10)= Dlu0(36*i-10) - RHS_Aij(5,2)
            Dlu0(36*i-9 )= Dlu0(36*i-9 ) - RHS_Aij(5,3)
            Dlu0(36*i-8 )= Dlu0(36*i-8 ) - RHS_Aij(5,4)
            Dlu0(36*i-7 )= Dlu0(36*i-7 ) - RHS_Aij(5,5)
            Dlu0(36*i-6 )= Dlu0(36*i-6 ) - RHS_Aij(5,6)
            Dlu0(36*i-5 )= Dlu0(36*i-5 ) - RHS_Aij(6,1)
            Dlu0(36*i-4 )= Dlu0(36*i-4 ) - RHS_Aij(6,2)
            Dlu0(36*i-3 )= Dlu0(36*i-3 ) - RHS_Aij(6,3)
            Dlu0(36*i-2 )= Dlu0(36*i-2 ) - RHS_Aij(6,4)
            Dlu0(36*i-1 )= Dlu0(36*i-1 ) - RHS_Aij(6,5)
            Dlu0(36*i   )= Dlu0(36*i   ) - RHS_Aij(6,6)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(36*ij0-35)= ALlu0(36*ij0-35) - RHS_Aij(1,1)
            ALlu0(36*ij0-34)= ALlu0(36*ij0-34) - RHS_Aij(1,2)
            ALlu0(36*ij0-33)= ALlu0(36*ij0-33) - RHS_Aij(1,3)
            ALlu0(36*ij0-32)= ALlu0(36*ij0-32) - RHS_Aij(1,4)
            ALlu0(36*ij0-31)= ALlu0(36*ij0-31) - RHS_Aij(1,5)
            ALlu0(36*ij0-30)= ALlu0(36*ij0-30) - RHS_Aij(1,6)
            ALlu0(36*ij0-29)= ALlu0(36*ij0-29) - RHS_Aij(2,1)
            ALlu0(36*ij0-28)= ALlu0(36*ij0-28) - RHS_Aij(2,2)
            ALlu0(36*ij0-27)= ALlu0(36*ij0-27) - RHS_Aij(2,3)
            ALlu0(36*ij0-26)= ALlu0(36*ij0-26) - RHS_Aij(2,4)
            ALlu0(36*ij0-25)= ALlu0(36*ij0-25) - RHS_Aij(2,5)
            ALlu0(36*ij0-24)= ALlu0(36*ij0-24) - RHS_Aij(2,6)
            ALlu0(36*ij0-23)= ALlu0(36*ij0-23) - RHS_Aij(3,1)
            ALlu0(36*ij0-22)= ALlu0(36*ij0-22) - RHS_Aij(3,2)
            ALlu0(36*ij0-21)= ALlu0(36*ij0-21) - RHS_Aij(3,3)
            ALlu0(36*ij0-20)= ALlu0(36*ij0-20) - RHS_Aij(3,4)
            ALlu0(36*ij0-19)= ALlu0(36*ij0-19) - RHS_Aij(3,5)
            ALlu0(36*ij0-18)= ALlu0(36*ij0-18) - RHS_Aij(3,6)
            ALlu0(36*ij0-17)= ALlu0(36*ij0-17) - RHS_Aij(4,1)
            ALlu0(36*ij0-16)= ALlu0(36*ij0-16) - RHS_Aij(4,2)
            ALlu0(36*ij0-15)= ALlu0(36*ij0-15) - RHS_Aij(4,3)
            ALlu0(36*ij0-14)= ALlu0(36*ij0-14) - RHS_Aij(4,4)
            ALlu0(36*ij0-13)= ALlu0(36*ij0-13) - RHS_Aij(4,5)
            ALlu0(36*ij0-12)= ALlu0(36*ij0-12) - RHS_Aij(4,6)
            ALlu0(36*ij0-11)= ALlu0(36*ij0-11) - RHS_Aij(5,1)
            ALlu0(36*ij0-10)= ALlu0(36*ij0-10) - RHS_Aij(5,2)
            ALlu0(36*ij0-9 )= ALlu0(36*ij0-9 ) - RHS_Aij(5,3)
            ALlu0(36*ij0-8 )= ALlu0(36*ij0-8 ) - RHS_Aij(5,4)
            ALlu0(36*ij0-7 )= ALlu0(36*ij0-7 ) - RHS_Aij(5,5)
            ALlu0(36*ij0-6 )= ALlu0(36*ij0-6 ) - RHS_Aij(5,6)
            ALlu0(36*ij0-5 )= ALlu0(36*ij0-5 ) - RHS_Aij(6,1)
            ALlu0(36*ij0-4 )= ALlu0(36*ij0-4 ) - RHS_Aij(6,2)
            ALlu0(36*ij0-3 )= ALlu0(36*ij0-3 ) - RHS_Aij(6,3)
            ALlu0(36*ij0-2 )= ALlu0(36*ij0-2 ) - RHS_Aij(6,4)
            ALlu0(36*ij0-1 )= ALlu0(36*ij0-1 ) - RHS_Aij(6,5)
            ALlu0(36*ij0   )= ALlu0(36*ij0   ) - RHS_Aij(6,6)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(36*ij0-35)= AUlu0(36*ij0-35) - RHS_Aij(1,1)
            AUlu0(36*ij0-34)= AUlu0(36*ij0-34) - RHS_Aij(1,2)
            AUlu0(36*ij0-33)= AUlu0(36*ij0-33) - RHS_Aij(1,3)
            AUlu0(36*ij0-32)= AUlu0(36*ij0-32) - RHS_Aij(1,4)
            AUlu0(36*ij0-31)= AUlu0(36*ij0-31) - RHS_Aij(1,5)
            AUlu0(36*ij0-30)= AUlu0(36*ij0-30) - RHS_Aij(1,6)
            AUlu0(36*ij0-29)= AUlu0(36*ij0-29) - RHS_Aij(2,1)
            AUlu0(36*ij0-28)= AUlu0(36*ij0-28) - RHS_Aij(2,2)
            AUlu0(36*ij0-27)= AUlu0(36*ij0-27) - RHS_Aij(2,3)
            AUlu0(36*ij0-26)= AUlu0(36*ij0-26) - RHS_Aij(2,4)
            AUlu0(36*ij0-25)= AUlu0(36*ij0-25) - RHS_Aij(2,5)
            AUlu0(36*ij0-24)= AUlu0(36*ij0-24) - RHS_Aij(2,6)
            AUlu0(36*ij0-23)= AUlu0(36*ij0-23) - RHS_Aij(3,1)
            AUlu0(36*ij0-22)= AUlu0(36*ij0-22) - RHS_Aij(3,2)
            AUlu0(36*ij0-21)= AUlu0(36*ij0-21) - RHS_Aij(3,3)
            AUlu0(36*ij0-20)= AUlu0(36*ij0-20) - RHS_Aij(3,4)
            AUlu0(36*ij0-19)= AUlu0(36*ij0-19) - RHS_Aij(3,5)
            AUlu0(36*ij0-18)= AUlu0(36*ij0-18) - RHS_Aij(3,6)
            AUlu0(36*ij0-17)= AUlu0(36*ij0-17) - RHS_Aij(4,1)
            AUlu0(36*ij0-16)= AUlu0(36*ij0-16) - RHS_Aij(4,2)
            AUlu0(36*ij0-15)= AUlu0(36*ij0-15) - RHS_Aij(4,3)
            AUlu0(36*ij0-14)= AUlu0(36*ij0-14) - RHS_Aij(4,4)
            AUlu0(36*ij0-13)= AUlu0(36*ij0-13) - RHS_Aij(4,5)
            AUlu0(36*ij0-12)= AUlu0(36*ij0-12) - RHS_Aij(4,6)
            AUlu0(36*ij0-11)= AUlu0(36*ij0-11) - RHS_Aij(5,1)
            AUlu0(36*ij0-10)= AUlu0(36*ij0-10) - RHS_Aij(5,2)
            AUlu0(36*ij0-9 )= AUlu0(36*ij0-9 ) - RHS_Aij(5,3)
            AUlu0(36*ij0-8 )= AUlu0(36*ij0-8 ) - RHS_Aij(5,4)
            AUlu0(36*ij0-7 )= AUlu0(36*ij0-7 ) - RHS_Aij(5,5)
            AUlu0(36*ij0-6 )= AUlu0(36*ij0-6 ) - RHS_Aij(5,6)
            AUlu0(36*ij0-5 )= AUlu0(36*ij0-5 ) - RHS_Aij(6,1)
            AUlu0(36*ij0-4 )= AUlu0(36*ij0-4 ) - RHS_Aij(6,2)
            AUlu0(36*ij0-3 )= AUlu0(36*ij0-3 ) - RHS_Aij(6,3)
            AUlu0(36*ij0-2 )= AUlu0(36*ij0-2 ) - RHS_Aij(6,4)
            AUlu0(36*ij0-1 )= AUlu0(36*ij0-1 ) - RHS_Aij(6,5)
            AUlu0(36*ij0   )= AUlu0(36*ij0   ) - RHS_Aij(6,6)
          endif

        enddo
      enddo
    enddo

    deallocate (IW1, IW2)
  end subroutine FORM_ILU0_66

  !C
  !C***
  !C*** FORM_ILU1_66
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU1_66                                   &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(36*NPL), intent(in):: AL
    real(kind=kreal), dimension(36*NPU), intent(in):: AU
    real(kind=kreal), dimension(36*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    integer(kind=kint), dimension(:), allocatable :: IWsL, IWsU
    real (kind=kreal),  dimension(6,6) :: RHS_Aij, DkINV, Aik, Akj
    real(kind=kreal) :: D11,D12,D13,D14,D15,D16
    real(kind=kreal) :: D21,D22,D23,D24,D25,D26
    real(kind=kreal) :: D31,D32,D33,D34,D35,D36
    real(kind=kreal) :: D41,D42,D43,D44,D45,D46
    real(kind=kreal) :: D51,D52,D53,D54,D55,D56
    real(kind=kreal) :: D61,D62,D63,D64,D65,D66
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
    allocate (ALlu0(36*(NPL+NPLf1)), AUlu0(36*(NPU+NPUf1)))

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
      call fill_in_S66_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 36*( k+isL)
          kk2= 36*(ik+INL(i-1))
          ALlu0(kk1-35)= AL(kk2-35)
          ALlu0(kk1-34)= AL(kk2-34)
          ALlu0(kk1-33)= AL(kk2-33)
          ALlu0(kk1-32)= AL(kk2-32)
          ALlu0(kk1-31)= AL(kk2-31)
          ALlu0(kk1-30)= AL(kk2-30)
          ALlu0(kk1-29)= AL(kk2-29)
          ALlu0(kk1-28)= AL(kk2-28)
          ALlu0(kk1-27)= AL(kk2-27)
          ALlu0(kk1-26)= AL(kk2-26)
          ALlu0(kk1-25)= AL(kk2-25)
          ALlu0(kk1-24)= AL(kk2-24)
          ALlu0(kk1-23)= AL(kk2-23)
          ALlu0(kk1-22)= AL(kk2-22)
          ALlu0(kk1-21)= AL(kk2-21)
          ALlu0(kk1-20)= AL(kk2-20)
          ALlu0(kk1-19)= AL(kk2-19)
          ALlu0(kk1-18)= AL(kk2-18)
          ALlu0(kk1-17)= AL(kk2-17)
          ALlu0(kk1-16)= AL(kk2-16)
          ALlu0(kk1-15)= AL(kk2-15)
          ALlu0(kk1-14)= AL(kk2-14)
          ALlu0(kk1-13)= AL(kk2-13)
          ALlu0(kk1-12)= AL(kk2-12)
          ALlu0(kk1-11)= AL(kk2-11)
          ALlu0(kk1-10)= AL(kk2-10)
          ALlu0(kk1-9 )= AL(kk2-9 )
          ALlu0(kk1-8 )= AL(kk2-8 )
          ALlu0(kk1-7 )= AL(kk2-7 )
          ALlu0(kk1-6 )= AL(kk2-6 )
          ALlu0(kk1-5 )= AL(kk2-5 )
          ALlu0(kk1-4 )= AL(kk2-4 )
          ALlu0(kk1-3 )= AL(kk2-3 )
          ALlu0(kk1-2 )= AL(kk2-2 )
          ALlu0(kk1-1 )= AL(kk2-1 )
          ALlu0(kk1   )= AL(kk2   )
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
      call fill_in_S66_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 36*( k+isU)
          kk2= 36*(ik+INU(i-1))
          AUlu0(kk1-35)= AU(kk2-35)
          AUlu0(kk1-34)= AU(kk2-34)
          AUlu0(kk1-33)= AU(kk2-33)
          AUlu0(kk1-32)= AU(kk2-32)
          AUlu0(kk1-31)= AU(kk2-31)
          AUlu0(kk1-30)= AU(kk2-30)
          AUlu0(kk1-29)= AU(kk2-29)
          AUlu0(kk1-28)= AU(kk2-28)
          AUlu0(kk1-27)= AU(kk2-27)
          AUlu0(kk1-26)= AU(kk2-26)
          AUlu0(kk1-25)= AU(kk2-25)
          AUlu0(kk1-24)= AU(kk2-24)
          AUlu0(kk1-23)= AU(kk2-23)
          AUlu0(kk1-22)= AU(kk2-22)
          AUlu0(kk1-21)= AU(kk2-21)
          AUlu0(kk1-20)= AU(kk2-20)
          AUlu0(kk1-19)= AU(kk2-19)
          AUlu0(kk1-18)= AU(kk2-18)
          AUlu0(kk1-17)= AU(kk2-17)
          AUlu0(kk1-16)= AU(kk2-16)
          AUlu0(kk1-15)= AU(kk2-15)
          AUlu0(kk1-14)= AU(kk2-14)
          AUlu0(kk1-13)= AU(kk2-13)
          AUlu0(kk1-12)= AU(kk2-12)
          AUlu0(kk1-11)= AU(kk2-11)
          AUlu0(kk1-10)= AU(kk2-10)
          AUlu0(kk1-9 )= AU(kk2-9 )
          AUlu0(kk1-8 )= AU(kk2-8 )
          AUlu0(kk1-7 )= AU(kk2-7 )
          AUlu0(kk1-6 )= AU(kk2-6 )
          AUlu0(kk1-5 )= AU(kk2-5 )
          AUlu0(kk1-4 )= AU(kk2-4 )
          AUlu0(kk1-3 )= AU(kk2-3 )
          AUlu0(kk1-2 )= AU(kk2-2 )
          AUlu0(kk1-1 )= AU(kk2-1 )
          AUlu0(kk1   )= AU(kk2   )
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
    allocate (Dlu0(36*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(36*i-35)=Dlu0(36*i-35)*SIGMA_DIAG
      Dlu0(36*i-28)=Dlu0(36*i-28)*SIGMA_DIAG
      Dlu0(36*i-21)=Dlu0(36*i-21)*SIGMA_DIAG
      Dlu0(36*i-14)=Dlu0(36*i-14)*SIGMA_DIAG
      Dlu0(36*i-7 )=Dlu0(36*i-7 )*SIGMA_DIAG
      Dlu0(36*i   )=Dlu0(36*i   )*SIGMA_DIAG
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
        D11= Dlu0(36*k-35)
        D12= Dlu0(36*k-34)
        D13= Dlu0(36*k-33)
        D14= Dlu0(36*k-32)
        D15= Dlu0(36*k-31)
        D16= Dlu0(36*k-30)
        D21= Dlu0(36*k-29)
        D22= Dlu0(36*k-28)
        D23= Dlu0(36*k-27)
        D24= Dlu0(36*k-26)
        D25= Dlu0(36*k-25)
        D26= Dlu0(36*k-24)
        D31= Dlu0(36*k-23)
        D32= Dlu0(36*k-22)
        D33= Dlu0(36*k-21)
        D34= Dlu0(36*k-20)
        D35= Dlu0(36*k-19)
        D36= Dlu0(36*k-18)
        D41= Dlu0(36*k-17)
        D42= Dlu0(36*k-16)
        D43= Dlu0(36*k-15)
        D44= Dlu0(36*k-14)
        D45= Dlu0(36*k-13)
        D46= Dlu0(36*k-12)
        D51= Dlu0(36*k-11)
        D52= Dlu0(36*k-10)
        D53= Dlu0(36*k-9 )
        D54= Dlu0(36*k-8 )
        D55= Dlu0(36*k-7 )
        D56= Dlu0(36*k-6 )
        D61= Dlu0(36*k-5 )
        D62= Dlu0(36*k-4 )
        D63= Dlu0(36*k-3 )
        D64= Dlu0(36*k-2 )
        D65= Dlu0(36*k-1 )
        D66= Dlu0(36*k   )

        call ILU1a66 (DkINV,D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26, &
          & D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56, &
          & D61,D62,D63,D64,D65,D66)

        do kk1= inumFI1L(i-1)+1, inumFI1L(i)
          if (k.eq.FI1L(kk1)) then
            Aik(1,1)= ALlu0(36*kk1-35)
            Aik(1,2)= ALlu0(36*kk1-34)
            Aik(1,3)= ALlu0(36*kk1-33)
            Aik(1,4)= ALlu0(36*kk1-32)
            Aik(1,5)= ALlu0(36*kk1-31)
            Aik(1,6)= ALlu0(36*kk1-30)
            Aik(2,1)= ALlu0(36*kk1-29)
            Aik(2,2)= ALlu0(36*kk1-28)
            Aik(2,3)= ALlu0(36*kk1-27)
            Aik(2,4)= ALlu0(36*kk1-26)
            Aik(2,5)= ALlu0(36*kk1-25)
            Aik(2,6)= ALlu0(36*kk1-24)
            Aik(3,1)= ALlu0(36*kk1-23)
            Aik(3,2)= ALlu0(36*kk1-22)
            Aik(3,3)= ALlu0(36*kk1-21)
            Aik(3,4)= ALlu0(36*kk1-20)
            Aik(3,5)= ALlu0(36*kk1-19)
            Aik(3,6)= ALlu0(36*kk1-18)
            Aik(4,1)= ALlu0(36*kk1-17)
            Aik(4,2)= ALlu0(36*kk1-16)
            Aik(4,3)= ALlu0(36*kk1-15)
            Aik(4,4)= ALlu0(36*kk1-14)
            Aik(4,5)= ALlu0(36*kk1-13)
            Aik(4,6)= ALlu0(36*kk1-12)
            Aik(5,1)= ALlu0(36*kk1-11)
            Aik(5,2)= ALlu0(36*kk1-10)
            Aik(5,3)= ALlu0(36*kk1-9)
            Aik(5,4)= ALlu0(36*kk1-8)
            Aik(5,5)= ALlu0(36*kk1-7)
            Aik(5,6)= ALlu0(36*kk1-6)
            Aik(6,1)= ALlu0(36*kk1-5)
            Aik(6,2)= ALlu0(36*kk1-4)
            Aik(6,3)= ALlu0(36*kk1-3)
            Aik(6,4)= ALlu0(36*kk1-2)
            Aik(6,5)= ALlu0(36*kk1-1)
            Aik(6,6)= ALlu0(36*kk1  )
            exit
          endif
        enddo

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          do jj1= inumFI1U(k-1)+1, inumFI1U(k)
            if (j.eq.FI1U(jj1)) then
              Akj(1,1)= AUlu0(36*jj1-35)
              Akj(1,2)= AUlu0(36*jj1-34)
              Akj(1,3)= AUlu0(36*jj1-33)
              Akj(1,4)= AUlu0(36*jj1-32)
              Akj(1,5)= AUlu0(36*jj1-31)
              Akj(1,6)= AUlu0(36*jj1-30)
              Akj(2,1)= AUlu0(36*jj1-29)
              Akj(2,2)= AUlu0(36*jj1-28)
              Akj(2,3)= AUlu0(36*jj1-27)
              Akj(2,4)= AUlu0(36*jj1-26)
              Akj(2,5)= AUlu0(36*jj1-25)
              Akj(2,6)= AUlu0(36*jj1-24)
              Akj(3,1)= AUlu0(36*jj1-23)
              Akj(3,2)= AUlu0(36*jj1-22)
              Akj(3,3)= AUlu0(36*jj1-21)
              Akj(3,4)= AUlu0(36*jj1-20)
              Akj(3,5)= AUlu0(36*jj1-19)
              Akj(3,6)= AUlu0(36*jj1-18)
              Akj(4,1)= AUlu0(36*jj1-17)
              Akj(4,2)= AUlu0(36*jj1-16)
              Akj(4,3)= AUlu0(36*jj1-15)
              Akj(4,4)= AUlu0(36*jj1-14)
              Akj(4,5)= AUlu0(36*jj1-13)
              Akj(4,6)= AUlu0(36*jj1-12)
              Akj(5,1)= AUlu0(36*jj1-11)
              Akj(5,2)= AUlu0(36*jj1-10)
              Akj(5,3)= AUlu0(36*jj1-9)
              Akj(5,4)= AUlu0(36*jj1-8)
              Akj(5,5)= AUlu0(36*jj1-7)
              Akj(5,6)= AUlu0(36*jj1-6)
              Akj(6,1)= AUlu0(36*jj1-5)
              Akj(6,2)= AUlu0(36*jj1-4)
              Akj(6,3)= AUlu0(36*jj1-3)
              Akj(6,4)= AUlu0(36*jj1-2)
              Akj(6,5)= AUlu0(36*jj1-1)
              Akj(6,6)= AUlu0(36*jj1  )
              exit
            endif
          enddo

          call ILU1b66 (RHS_Aij, DkINV, Aik, Akj)

          if (j.eq.i) then
            Dlu0(36*i-35)= Dlu0(36*i-35) - RHS_Aij(1,1)
            Dlu0(36*i-34)= Dlu0(36*i-34) - RHS_Aij(1,2)
            Dlu0(36*i-33)= Dlu0(36*i-33) - RHS_Aij(1,3)
            Dlu0(36*i-32)= Dlu0(36*i-32) - RHS_Aij(1,4)
            Dlu0(36*i-31)= Dlu0(36*i-31) - RHS_Aij(1,5)
            Dlu0(36*i-30)= Dlu0(36*i-30) - RHS_Aij(1,6)
            Dlu0(36*i-29)= Dlu0(36*i-29) - RHS_Aij(2,1)
            Dlu0(36*i-28)= Dlu0(36*i-28) - RHS_Aij(2,2)
            Dlu0(36*i-27)= Dlu0(36*i-27) - RHS_Aij(2,3)
            Dlu0(36*i-26)= Dlu0(36*i-26) - RHS_Aij(2,4)
            Dlu0(36*i-25)= Dlu0(36*i-25) - RHS_Aij(2,5)
            Dlu0(36*i-24)= Dlu0(36*i-24) - RHS_Aij(2,6)
            Dlu0(36*i-23)= Dlu0(36*i-23) - RHS_Aij(3,1)
            Dlu0(36*i-22)= Dlu0(36*i-22) - RHS_Aij(3,2)
            Dlu0(36*i-21)= Dlu0(36*i-21) - RHS_Aij(3,3)
            Dlu0(36*i-20)= Dlu0(36*i-20) - RHS_Aij(3,4)
            Dlu0(36*i-19)= Dlu0(36*i-19) - RHS_Aij(3,5)
            Dlu0(36*i-18)= Dlu0(36*i-18) - RHS_Aij(3,6)
            Dlu0(36*i-17)= Dlu0(36*i-17) - RHS_Aij(4,1)
            Dlu0(36*i-16)= Dlu0(36*i-16) - RHS_Aij(4,2)
            Dlu0(36*i-15)= Dlu0(36*i-15) - RHS_Aij(4,3)
            Dlu0(36*i-14)= Dlu0(36*i-14) - RHS_Aij(4,4)
            Dlu0(36*i-13)= Dlu0(36*i-13) - RHS_Aij(4,5)
            Dlu0(36*i-12)= Dlu0(36*i-12) - RHS_Aij(4,6)
            Dlu0(36*i-11)= Dlu0(36*i-11) - RHS_Aij(5,1)
            Dlu0(36*i-10)= Dlu0(36*i-10) - RHS_Aij(5,2)
            Dlu0(36*i-9 )= Dlu0(36*i-9 ) - RHS_Aij(5,3)
            Dlu0(36*i-8 )= Dlu0(36*i-8 ) - RHS_Aij(5,4)
            Dlu0(36*i-7 )= Dlu0(36*i-7 ) - RHS_Aij(5,5)
            Dlu0(36*i-6 )= Dlu0(36*i-6 ) - RHS_Aij(5,6)
            Dlu0(36*i-5 )= Dlu0(36*i-5 ) - RHS_Aij(6,1)
            Dlu0(36*i-4 )= Dlu0(36*i-4 ) - RHS_Aij(6,2)
            Dlu0(36*i-3 )= Dlu0(36*i-3 ) - RHS_Aij(6,3)
            Dlu0(36*i-2 )= Dlu0(36*i-2 ) - RHS_Aij(6,4)
            Dlu0(36*i-1 )= Dlu0(36*i-1 ) - RHS_Aij(6,5)
            Dlu0(36*i   )= Dlu0(36*i   ) - RHS_Aij(6,6)
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            ALlu0(36*ij0-35)= ALlu0(36*ij0-35) - RHS_Aij(1,1)
            ALlu0(36*ij0-34)= ALlu0(36*ij0-34) - RHS_Aij(1,2)
            ALlu0(36*ij0-33)= ALlu0(36*ij0-33) - RHS_Aij(1,3)
            ALlu0(36*ij0-32)= ALlu0(36*ij0-32) - RHS_Aij(1,4)
            ALlu0(36*ij0-31)= ALlu0(36*ij0-31) - RHS_Aij(1,5)
            ALlu0(36*ij0-30)= ALlu0(36*ij0-30) - RHS_Aij(1,6)
            ALlu0(36*ij0-29)= ALlu0(36*ij0-29) - RHS_Aij(2,1)
            ALlu0(36*ij0-28)= ALlu0(36*ij0-28) - RHS_Aij(2,2)
            ALlu0(36*ij0-27)= ALlu0(36*ij0-27) - RHS_Aij(2,3)
            ALlu0(36*ij0-26)= ALlu0(36*ij0-26) - RHS_Aij(2,4)
            ALlu0(36*ij0-25)= ALlu0(36*ij0-25) - RHS_Aij(2,5)
            ALlu0(36*ij0-24)= ALlu0(36*ij0-24) - RHS_Aij(2,6)
            ALlu0(36*ij0-23)= ALlu0(36*ij0-23) - RHS_Aij(3,1)
            ALlu0(36*ij0-22)= ALlu0(36*ij0-22) - RHS_Aij(3,2)
            ALlu0(36*ij0-21)= ALlu0(36*ij0-21) - RHS_Aij(3,3)
            ALlu0(36*ij0-20)= ALlu0(36*ij0-20) - RHS_Aij(3,4)
            ALlu0(36*ij0-19)= ALlu0(36*ij0-19) - RHS_Aij(3,5)
            ALlu0(36*ij0-18)= ALlu0(36*ij0-18) - RHS_Aij(3,6)
            ALlu0(36*ij0-17)= ALlu0(36*ij0-17) - RHS_Aij(4,1)
            ALlu0(36*ij0-16)= ALlu0(36*ij0-16) - RHS_Aij(4,2)
            ALlu0(36*ij0-15)= ALlu0(36*ij0-15) - RHS_Aij(4,3)
            ALlu0(36*ij0-14)= ALlu0(36*ij0-14) - RHS_Aij(4,4)
            ALlu0(36*ij0-13)= ALlu0(36*ij0-13) - RHS_Aij(4,5)
            ALlu0(36*ij0-12)= ALlu0(36*ij0-12) - RHS_Aij(4,6)
            ALlu0(36*ij0-11)= ALlu0(36*ij0-11) - RHS_Aij(5,1)
            ALlu0(36*ij0-10)= ALlu0(36*ij0-10) - RHS_Aij(5,2)
            ALlu0(36*ij0-9 )= ALlu0(36*ij0-9 ) - RHS_Aij(5,3)
            ALlu0(36*ij0-8 )= ALlu0(36*ij0-8 ) - RHS_Aij(5,4)
            ALlu0(36*ij0-7 )= ALlu0(36*ij0-7 ) - RHS_Aij(5,5)
            ALlu0(36*ij0-6 )= ALlu0(36*ij0-6 ) - RHS_Aij(5,6)
            ALlu0(36*ij0-5 )= ALlu0(36*ij0-5 ) - RHS_Aij(6,1)
            ALlu0(36*ij0-4 )= ALlu0(36*ij0-4 ) - RHS_Aij(6,2)
            ALlu0(36*ij0-3 )= ALlu0(36*ij0-3 ) - RHS_Aij(6,3)
            ALlu0(36*ij0-2 )= ALlu0(36*ij0-2 ) - RHS_Aij(6,4)
            ALlu0(36*ij0-1 )= ALlu0(36*ij0-1 ) - RHS_Aij(6,5)
            ALlu0(36*ij0   )= ALlu0(36*ij0   ) - RHS_Aij(6,6)
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            AUlu0(36*ij0-35)= AUlu0(36*ij0-35) - RHS_Aij(1,1)
            AUlu0(36*ij0-34)= AUlu0(36*ij0-34) - RHS_Aij(1,2)
            AUlu0(36*ij0-33)= AUlu0(36*ij0-33) - RHS_Aij(1,3)
            AUlu0(36*ij0-32)= AUlu0(36*ij0-32) - RHS_Aij(1,4)
            AUlu0(36*ij0-31)= AUlu0(36*ij0-31) - RHS_Aij(1,5)
            AUlu0(36*ij0-30)= AUlu0(36*ij0-30) - RHS_Aij(1,6)
            AUlu0(36*ij0-29)= AUlu0(36*ij0-29) - RHS_Aij(2,1)
            AUlu0(36*ij0-28)= AUlu0(36*ij0-28) - RHS_Aij(2,2)
            AUlu0(36*ij0-27)= AUlu0(36*ij0-27) - RHS_Aij(2,3)
            AUlu0(36*ij0-26)= AUlu0(36*ij0-26) - RHS_Aij(2,4)
            AUlu0(36*ij0-25)= AUlu0(36*ij0-25) - RHS_Aij(2,5)
            AUlu0(36*ij0-24)= AUlu0(36*ij0-24) - RHS_Aij(2,6)
            AUlu0(36*ij0-23)= AUlu0(36*ij0-23) - RHS_Aij(3,1)
            AUlu0(36*ij0-22)= AUlu0(36*ij0-22) - RHS_Aij(3,2)
            AUlu0(36*ij0-21)= AUlu0(36*ij0-21) - RHS_Aij(3,3)
            AUlu0(36*ij0-20)= AUlu0(36*ij0-20) - RHS_Aij(3,4)
            AUlu0(36*ij0-19)= AUlu0(36*ij0-19) - RHS_Aij(3,5)
            AUlu0(36*ij0-18)= AUlu0(36*ij0-18) - RHS_Aij(3,6)
            AUlu0(36*ij0-17)= AUlu0(36*ij0-17) - RHS_Aij(4,1)
            AUlu0(36*ij0-16)= AUlu0(36*ij0-16) - RHS_Aij(4,2)
            AUlu0(36*ij0-15)= AUlu0(36*ij0-15) - RHS_Aij(4,3)
            AUlu0(36*ij0-14)= AUlu0(36*ij0-14) - RHS_Aij(4,4)
            AUlu0(36*ij0-13)= AUlu0(36*ij0-13) - RHS_Aij(4,5)
            AUlu0(36*ij0-12)= AUlu0(36*ij0-12) - RHS_Aij(4,6)
            AUlu0(36*ij0-11)= AUlu0(36*ij0-11) - RHS_Aij(5,1)
            AUlu0(36*ij0-10)= AUlu0(36*ij0-10) - RHS_Aij(5,2)
            AUlu0(36*ij0-9 )= AUlu0(36*ij0-9 ) - RHS_Aij(5,3)
            AUlu0(36*ij0-8 )= AUlu0(36*ij0-8 ) - RHS_Aij(5,4)
            AUlu0(36*ij0-7 )= AUlu0(36*ij0-7 ) - RHS_Aij(5,5)
            AUlu0(36*ij0-6 )= AUlu0(36*ij0-6 ) - RHS_Aij(5,6)
            AUlu0(36*ij0-5 )= AUlu0(36*ij0-5 ) - RHS_Aij(6,1)
            AUlu0(36*ij0-4 )= AUlu0(36*ij0-4 ) - RHS_Aij(6,2)
            AUlu0(36*ij0-3 )= AUlu0(36*ij0-3 ) - RHS_Aij(6,3)
            AUlu0(36*ij0-2 )= AUlu0(36*ij0-2 ) - RHS_Aij(6,4)
            AUlu0(36*ij0-1 )= AUlu0(36*ij0-1 ) - RHS_Aij(6,5)
            AUlu0(36*ij0   )= AUlu0(36*ij0   ) - RHS_Aij(6,6)
          endif

        enddo
      enddo
    enddo

    deallocate (IW1, IW2)
    !C===
  end subroutine FORM_ILU1_66

  !C
  !C***
  !C*** FORM_ILU2_66
  !C***
  !C
  !C    form ILU(2) matrix
  !C
  subroutine FORM_ILU2_66 &
      &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(36*NPL), intent(in):: AL
    real(kind=kreal), dimension(36*NPU), intent(in):: AU
    real(kind=kreal), dimension(36*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable:: IW1 , IW2
    integer(kind=kint), dimension(:), allocatable:: IWsL, IWsU
    integer(kind=kint), dimension(:), allocatable:: iconFI1L, iconFI1U
    integer(kind=kint), dimension(:), allocatable:: inumFI2L, inumFI2U
    integer(kind=kint), dimension(:), allocatable::     FI2L,     FI2U
    real (kind=kreal), dimension(6,6) :: RHS_Aij, DkINV, Aik, Akj
    real(kind=kreal) :: D11,D12,D13,D14,D15,D16
    real(kind=kreal) :: D21,D22,D23,D24,D25,D26
    real(kind=kreal) :: D31,D32,D33,D34,D35,D36
    real(kind=kreal) :: D41,D42,D43,D44,D45,D46
    real(kind=kreal) :: D51,D52,D53,D54,D55,D56
    real(kind=kreal) :: D61,D62,D63,D64,D65,D66
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

      call fill_in_S66_SORT (IW1, IW2, icouL3, NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
        ik= IW2(k)
        if (ik.le.INL(i)-INL(i-1)) then
          kk1= 36*( k+isL)
          kk2= 36*(ik+INL(i-1))
          ALlu0(kk1-35)= AL(kk2-35)
          ALlu0(kk1-34)= AL(kk2-34)
          ALlu0(kk1-33)= AL(kk2-33)
          ALlu0(kk1-32)= AL(kk2-32)
          ALlu0(kk1-31)= AL(kk2-31)
          ALlu0(kk1-30)= AL(kk2-30)
          ALlu0(kk1-29)= AL(kk2-29)
          ALlu0(kk1-28)= AL(kk2-28)
          ALlu0(kk1-27)= AL(kk2-27)
          ALlu0(kk1-26)= AL(kk2-26)
          ALlu0(kk1-25)= AL(kk2-25)
          ALlu0(kk1-24)= AL(kk2-24)
          ALlu0(kk1-23)= AL(kk2-23)
          ALlu0(kk1-22)= AL(kk2-22)
          ALlu0(kk1-21)= AL(kk2-21)
          ALlu0(kk1-20)= AL(kk2-20)
          ALlu0(kk1-19)= AL(kk2-19)
          ALlu0(kk1-18)= AL(kk2-18)
          ALlu0(kk1-17)= AL(kk2-17)
          ALlu0(kk1-16)= AL(kk2-16)
          ALlu0(kk1-15)= AL(kk2-15)
          ALlu0(kk1-14)= AL(kk2-14)
          ALlu0(kk1-13)= AL(kk2-13)
          ALlu0(kk1-12)= AL(kk2-12)
          ALlu0(kk1-11)= AL(kk2-11)
          ALlu0(kk1-10)= AL(kk2-10)
          ALlu0(kk1-9 )= AL(kk2-9 )
          ALlu0(kk1-8 )= AL(kk2-8 )
          ALlu0(kk1-7 )= AL(kk2-7 )
          ALlu0(kk1-6 )= AL(kk2-6 )
          ALlu0(kk1-5 )= AL(kk2-5 )
          ALlu0(kk1-4 )= AL(kk2-4 )
          ALlu0(kk1-3 )= AL(kk2-3 )
          ALlu0(kk1-2 )= AL(kk2-2 )
          ALlu0(kk1-1 )= AL(kk2-1 )
          ALlu0(kk1   )= AL(kk2   )
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
      call fill_in_S66_SORT (IW1, IW2, icouU3, NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
        ik= IW2(k)
        if (ik.le.INU(i)-INU(i-1)) then
          kk1= 36*( k+isU)
          kk2= 36*(ik+INU(i-1))
          AUlu0(kk1-35)= AU(kk2-35)
          AUlu0(kk1-34)= AU(kk2-34)
          AUlu0(kk1-33)= AU(kk2-33)
          AUlu0(kk1-32)= AU(kk2-32)
          AUlu0(kk1-31)= AU(kk2-31)
          AUlu0(kk1-30)= AU(kk2-30)
          AUlu0(kk1-29)= AU(kk2-29)
          AUlu0(kk1-28)= AU(kk2-28)
          AUlu0(kk1-27)= AU(kk2-27)
          AUlu0(kk1-26)= AU(kk2-26)
          AUlu0(kk1-25)= AU(kk2-25)
          AUlu0(kk1-24)= AU(kk2-24)
          AUlu0(kk1-23)= AU(kk2-23)
          AUlu0(kk1-22)= AU(kk2-22)
          AUlu0(kk1-21)= AU(kk2-21)
          AUlu0(kk1-20)= AU(kk2-20)
          AUlu0(kk1-19)= AU(kk2-19)
          AUlu0(kk1-18)= AU(kk2-18)
          AUlu0(kk1-17)= AU(kk2-17)
          AUlu0(kk1-16)= AU(kk2-16)
          AUlu0(kk1-15)= AU(kk2-15)
          AUlu0(kk1-14)= AU(kk2-14)
          AUlu0(kk1-13)= AU(kk2-13)
          AUlu0(kk1-12)= AU(kk2-12)
          AUlu0(kk1-11)= AU(kk2-11)
          AUlu0(kk1-10)= AU(kk2-10)
          AUlu0(kk1-9 )= AU(kk2-9 )
          AUlu0(kk1-8 )= AU(kk2-8 )
          AUlu0(kk1-7 )= AU(kk2-7 )
          AUlu0(kk1-6 )= AU(kk2-6 )
          AUlu0(kk1-5 )= AU(kk2-5 )
          AUlu0(kk1-4 )= AU(kk2-4 )
          AUlu0(kk1-3 )= AU(kk2-3 )
          AUlu0(kk1-2 )= AU(kk2-2 )
          AUlu0(kk1-1 )= AU(kk2-1 )
          AUlu0(kk1   )= AU(kk2   )
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
    allocate (Dlu0(36*NP))
    Dlu0= D
    do i=1,NP
      Dlu0(36*i-35)=Dlu0(36*i-35)*SIGMA_DIAG
      Dlu0(36*i-28)=Dlu0(36*i-28)*SIGMA_DIAG
      Dlu0(36*i-21)=Dlu0(36*i-21)*SIGMA_DIAG
      Dlu0(36*i-14)=Dlu0(36*i-14)*SIGMA_DIAG
      Dlu0(36*i-7 )=Dlu0(36*i-7 )*SIGMA_DIAG
      Dlu0(36*i   )=Dlu0(36*i   )*SIGMA_DIAG
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

        D11= Dlu0(36*k-35)
        D12= Dlu0(36*k-34)
        D13= Dlu0(36*k-33)
        D14= Dlu0(36*k-32)
        D15= Dlu0(36*k-31)
        D16= Dlu0(36*k-30)
        D21= Dlu0(36*k-29)
        D22= Dlu0(36*k-28)
        D23= Dlu0(36*k-27)
        D24= Dlu0(36*k-26)
        D25= Dlu0(36*k-25)
        D26= Dlu0(36*k-24)
        D31= Dlu0(36*k-23)
        D32= Dlu0(36*k-22)
        D33= Dlu0(36*k-21)
        D34= Dlu0(36*k-20)
        D35= Dlu0(36*k-19)
        D36= Dlu0(36*k-18)
        D41= Dlu0(36*k-17)
        D42= Dlu0(36*k-16)
        D43= Dlu0(36*k-15)
        D44= Dlu0(36*k-14)
        D45= Dlu0(36*k-13)
        D46= Dlu0(36*k-12)
        D51= Dlu0(36*k-11)
        D52= Dlu0(36*k-10)
        D53= Dlu0(36*k-9 )
        D54= Dlu0(36*k-8 )
        D55= Dlu0(36*k-7 )
        D56= Dlu0(36*k-6 )
        D61= Dlu0(36*k-5 )
        D62= Dlu0(36*k-4 )
        D63= Dlu0(36*k-3 )
        D64= Dlu0(36*k-2 )
        D65= Dlu0(36*k-1 )
        D66= Dlu0(36*k   )

        call ILU1a66 (DkINV,D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26, &
          & D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56, &
          & D61,D62,D63,D64,D65,D66)

        Aik(1,1)= ALlu0(36*kk-35)
        Aik(1,2)= ALlu0(36*kk-34)
        Aik(1,3)= ALlu0(36*kk-33)
        Aik(1,4)= ALlu0(36*kk-32)
        Aik(1,5)= ALlu0(36*kk-31)
        Aik(1,6)= ALlu0(36*kk-30)
        Aik(2,1)= ALlu0(36*kk-29)
        Aik(2,2)= ALlu0(36*kk-28)
        Aik(2,3)= ALlu0(36*kk-27)
        Aik(2,4)= ALlu0(36*kk-26)
        Aik(2,5)= ALlu0(36*kk-25)
        Aik(2,6)= ALlu0(36*kk-24)
        Aik(3,1)= ALlu0(36*kk-23)
        Aik(3,2)= ALlu0(36*kk-22)
        Aik(3,3)= ALlu0(36*kk-21)
        Aik(3,4)= ALlu0(36*kk-20)
        Aik(3,5)= ALlu0(36*kk-19)
        Aik(3,6)= ALlu0(36*kk-18)
        Aik(4,1)= ALlu0(36*kk-17)
        Aik(4,2)= ALlu0(36*kk-16)
        Aik(4,3)= ALlu0(36*kk-15)
        Aik(4,4)= ALlu0(36*kk-14)
        Aik(4,5)= ALlu0(36*kk-13)
        Aik(4,6)= ALlu0(36*kk-12)
        Aik(5,1)= ALlu0(36*kk-11)
        Aik(5,2)= ALlu0(36*kk-10)
        Aik(5,3)= ALlu0(36*kk-9)
        Aik(5,4)= ALlu0(36*kk-8)
        Aik(5,5)= ALlu0(36*kk-7)
        Aik(5,6)= ALlu0(36*kk-6)
        Aik(6,1)= ALlu0(36*kk-5)
        Aik(6,2)= ALlu0(36*kk-4)
        Aik(6,3)= ALlu0(36*kk-3)
        Aik(6,4)= ALlu0(36*kk-2)
        Aik(6,5)= ALlu0(36*kk-1)
        Aik(6,6)= ALlu0(36*kk  )

        do jj= inumFI1U(k-1)+1, inumFI1U(k)
          j= FI1U(jj)
          iconKJ= iconFI1U(jj)

          if ((iconIK+iconKJ).lt.2) then
            Akj(1,1)= AUlu0(36*jj-35)
            Akj(1,2)= AUlu0(36*jj-34)
            Akj(1,3)= AUlu0(36*jj-33)
            Akj(1,4)= AUlu0(36*jj-32)
            Akj(1,5)= AUlu0(36*jj-31)
            Akj(1,6)= AUlu0(36*jj-30)
            Akj(2,1)= AUlu0(36*jj-29)
            Akj(2,2)= AUlu0(36*jj-28)
            Akj(2,3)= AUlu0(36*jj-27)
            Akj(2,4)= AUlu0(36*jj-26)
            Akj(2,5)= AUlu0(36*jj-25)
            Akj(2,6)= AUlu0(36*jj-24)
            Akj(3,1)= AUlu0(36*jj-23)
            Akj(3,2)= AUlu0(36*jj-22)
            Akj(3,3)= AUlu0(36*jj-21)
            Akj(3,4)= AUlu0(36*jj-20)
            Akj(3,5)= AUlu0(36*jj-19)
            Akj(3,6)= AUlu0(36*jj-18)
            Akj(4,1)= AUlu0(36*jj-17)
            Akj(4,2)= AUlu0(36*jj-16)
            Akj(4,3)= AUlu0(36*jj-15)
            Akj(4,4)= AUlu0(36*jj-14)
            Akj(4,5)= AUlu0(36*jj-13)
            Akj(4,6)= AUlu0(36*jj-12)
            Akj(5,1)= AUlu0(36*jj-11)
            Akj(5,2)= AUlu0(36*jj-10)
            Akj(5,3)= AUlu0(36*jj-9)
            Akj(5,4)= AUlu0(36*jj-8)
            Akj(5,5)= AUlu0(36*jj-7)
            Akj(5,6)= AUlu0(36*jj-6)
            Akj(6,1)= AUlu0(36*jj-5)
            Akj(6,2)= AUlu0(36*jj-4)
            Akj(6,3)= AUlu0(36*jj-3)
            Akj(6,4)= AUlu0(36*jj-2)
            Akj(6,5)= AUlu0(36*jj-1)
            Akj(6,6)= AUlu0(36*jj  )

            call ILU1b66 (RHS_Aij, DkINV, Aik, Akj)

            if (j.eq.i) then
              Dlu0(36*i-35)= Dlu0(36*i-35) - RHS_Aij(1,1)
              Dlu0(36*i-34)= Dlu0(36*i-34) - RHS_Aij(1,2)
              Dlu0(36*i-33)= Dlu0(36*i-33) - RHS_Aij(1,3)
              Dlu0(36*i-32)= Dlu0(36*i-32) - RHS_Aij(1,4)
              Dlu0(36*i-31)= Dlu0(36*i-31) - RHS_Aij(1,5)
              Dlu0(36*i-30)= Dlu0(36*i-30) - RHS_Aij(1,6)
              Dlu0(36*i-29)= Dlu0(36*i-29) - RHS_Aij(2,1)
              Dlu0(36*i-28)= Dlu0(36*i-28) - RHS_Aij(2,2)
              Dlu0(36*i-27)= Dlu0(36*i-27) - RHS_Aij(2,3)
              Dlu0(36*i-26)= Dlu0(36*i-26) - RHS_Aij(2,4)
              Dlu0(36*i-25)= Dlu0(36*i-25) - RHS_Aij(2,5)
              Dlu0(36*i-24)= Dlu0(36*i-24) - RHS_Aij(2,6)
              Dlu0(36*i-23)= Dlu0(36*i-23) - RHS_Aij(3,1)
              Dlu0(36*i-22)= Dlu0(36*i-22) - RHS_Aij(3,2)
              Dlu0(36*i-21)= Dlu0(36*i-21) - RHS_Aij(3,3)
              Dlu0(36*i-20)= Dlu0(36*i-20) - RHS_Aij(3,4)
              Dlu0(36*i-19)= Dlu0(36*i-19) - RHS_Aij(3,5)
              Dlu0(36*i-18)= Dlu0(36*i-18) - RHS_Aij(3,6)
              Dlu0(36*i-17)= Dlu0(36*i-17) - RHS_Aij(4,1)
              Dlu0(36*i-16)= Dlu0(36*i-16) - RHS_Aij(4,2)
              Dlu0(36*i-15)= Dlu0(36*i-15) - RHS_Aij(4,3)
              Dlu0(36*i-14)= Dlu0(36*i-14) - RHS_Aij(4,4)
              Dlu0(36*i-13)= Dlu0(36*i-13) - RHS_Aij(4,5)
              Dlu0(36*i-12)= Dlu0(36*i-12) - RHS_Aij(4,6)
              Dlu0(36*i-11)= Dlu0(36*i-11) - RHS_Aij(5,1)
              Dlu0(36*i-10)= Dlu0(36*i-10) - RHS_Aij(5,2)
              Dlu0(36*i-9 )= Dlu0(36*i-9 ) - RHS_Aij(5,3)
              Dlu0(36*i-8 )= Dlu0(36*i-8 ) - RHS_Aij(5,4)
              Dlu0(36*i-7 )= Dlu0(36*i-7 ) - RHS_Aij(5,5)
              Dlu0(36*i-6 )= Dlu0(36*i-6 ) - RHS_Aij(5,6)
              Dlu0(36*i-5 )= Dlu0(36*i-5 ) - RHS_Aij(6,1)
              Dlu0(36*i-4 )= Dlu0(36*i-4 ) - RHS_Aij(6,2)
              Dlu0(36*i-3 )= Dlu0(36*i-3 ) - RHS_Aij(6,3)
              Dlu0(36*i-2 )= Dlu0(36*i-2 ) - RHS_Aij(6,4)
              Dlu0(36*i-1 )= Dlu0(36*i-1 ) - RHS_Aij(6,5)
              Dlu0(36*i   )= Dlu0(36*i   ) - RHS_Aij(6,6)
            endif

            if (j.lt.i) then
              ij0= IW1(j)
              ALlu0(36*ij0-35)= ALlu0(36*ij0-35) - RHS_Aij(1,1)
              ALlu0(36*ij0-34)= ALlu0(36*ij0-34) - RHS_Aij(1,2)
              ALlu0(36*ij0-33)= ALlu0(36*ij0-33) - RHS_Aij(1,3)
              ALlu0(36*ij0-32)= ALlu0(36*ij0-32) - RHS_Aij(1,4)
              ALlu0(36*ij0-31)= ALlu0(36*ij0-31) - RHS_Aij(1,5)
              ALlu0(36*ij0-30)= ALlu0(36*ij0-30) - RHS_Aij(1,6)
              ALlu0(36*ij0-29)= ALlu0(36*ij0-29) - RHS_Aij(2,1)
              ALlu0(36*ij0-28)= ALlu0(36*ij0-28) - RHS_Aij(2,2)
              ALlu0(36*ij0-27)= ALlu0(36*ij0-27) - RHS_Aij(2,3)
              ALlu0(36*ij0-26)= ALlu0(36*ij0-26) - RHS_Aij(2,4)
              ALlu0(36*ij0-25)= ALlu0(36*ij0-25) - RHS_Aij(2,5)
              ALlu0(36*ij0-24)= ALlu0(36*ij0-24) - RHS_Aij(2,6)
              ALlu0(36*ij0-23)= ALlu0(36*ij0-23) - RHS_Aij(3,1)
              ALlu0(36*ij0-22)= ALlu0(36*ij0-22) - RHS_Aij(3,2)
              ALlu0(36*ij0-21)= ALlu0(36*ij0-21) - RHS_Aij(3,3)
              ALlu0(36*ij0-20)= ALlu0(36*ij0-20) - RHS_Aij(3,4)
              ALlu0(36*ij0-19)= ALlu0(36*ij0-19) - RHS_Aij(3,5)
              ALlu0(36*ij0-18)= ALlu0(36*ij0-18) - RHS_Aij(3,6)
              ALlu0(36*ij0-17)= ALlu0(36*ij0-17) - RHS_Aij(4,1)
              ALlu0(36*ij0-16)= ALlu0(36*ij0-16) - RHS_Aij(4,2)
              ALlu0(36*ij0-15)= ALlu0(36*ij0-15) - RHS_Aij(4,3)
              ALlu0(36*ij0-14)= ALlu0(36*ij0-14) - RHS_Aij(4,4)
              ALlu0(36*ij0-13)= ALlu0(36*ij0-13) - RHS_Aij(4,5)
              ALlu0(36*ij0-12)= ALlu0(36*ij0-12) - RHS_Aij(4,6)
              ALlu0(36*ij0-11)= ALlu0(36*ij0-11) - RHS_Aij(5,1)
              ALlu0(36*ij0-10)= ALlu0(36*ij0-10) - RHS_Aij(5,2)
              ALlu0(36*ij0-9 )= ALlu0(36*ij0-9 ) - RHS_Aij(5,3)
              ALlu0(36*ij0-8 )= ALlu0(36*ij0-8 ) - RHS_Aij(5,4)
              ALlu0(36*ij0-7 )= ALlu0(36*ij0-7 ) - RHS_Aij(5,5)
              ALlu0(36*ij0-6 )= ALlu0(36*ij0-6 ) - RHS_Aij(5,6)
              ALlu0(36*ij0-5 )= ALlu0(36*ij0-5 ) - RHS_Aij(6,1)
              ALlu0(36*ij0-4 )= ALlu0(36*ij0-4 ) - RHS_Aij(6,2)
              ALlu0(36*ij0-3 )= ALlu0(36*ij0-3 ) - RHS_Aij(6,3)
              ALlu0(36*ij0-2 )= ALlu0(36*ij0-2 ) - RHS_Aij(6,4)
              ALlu0(36*ij0-1 )= ALlu0(36*ij0-1 ) - RHS_Aij(6,5)
              ALlu0(36*ij0   )= ALlu0(36*ij0   ) - RHS_Aij(6,6)
            endif

            if (j.gt.i) then
              ij0= IW2(j)
              AUlu0(36*ij0-35)= AUlu0(36*ij0-35) - RHS_Aij(1,1)
              AUlu0(36*ij0-34)= AUlu0(36*ij0-34) - RHS_Aij(1,2)
              AUlu0(36*ij0-33)= AUlu0(36*ij0-33) - RHS_Aij(1,3)
              AUlu0(36*ij0-32)= AUlu0(36*ij0-32) - RHS_Aij(1,4)
              AUlu0(36*ij0-31)= AUlu0(36*ij0-31) - RHS_Aij(1,5)
              AUlu0(36*ij0-30)= AUlu0(36*ij0-30) - RHS_Aij(1,6)
              AUlu0(36*ij0-29)= AUlu0(36*ij0-29) - RHS_Aij(2,1)
              AUlu0(36*ij0-28)= AUlu0(36*ij0-28) - RHS_Aij(2,2)
              AUlu0(36*ij0-27)= AUlu0(36*ij0-27) - RHS_Aij(2,3)
              AUlu0(36*ij0-26)= AUlu0(36*ij0-26) - RHS_Aij(2,4)
              AUlu0(36*ij0-25)= AUlu0(36*ij0-25) - RHS_Aij(2,5)
              AUlu0(36*ij0-24)= AUlu0(36*ij0-24) - RHS_Aij(2,6)
              AUlu0(36*ij0-23)= AUlu0(36*ij0-23) - RHS_Aij(3,1)
              AUlu0(36*ij0-22)= AUlu0(36*ij0-22) - RHS_Aij(3,2)
              AUlu0(36*ij0-21)= AUlu0(36*ij0-21) - RHS_Aij(3,3)
              AUlu0(36*ij0-20)= AUlu0(36*ij0-20) - RHS_Aij(3,4)
              AUlu0(36*ij0-19)= AUlu0(36*ij0-19) - RHS_Aij(3,5)
              AUlu0(36*ij0-18)= AUlu0(36*ij0-18) - RHS_Aij(3,6)
              AUlu0(36*ij0-17)= AUlu0(36*ij0-17) - RHS_Aij(4,1)
              AUlu0(36*ij0-16)= AUlu0(36*ij0-16) - RHS_Aij(4,2)
              AUlu0(36*ij0-15)= AUlu0(36*ij0-15) - RHS_Aij(4,3)
              AUlu0(36*ij0-14)= AUlu0(36*ij0-14) - RHS_Aij(4,4)
              AUlu0(36*ij0-13)= AUlu0(36*ij0-13) - RHS_Aij(4,5)
              AUlu0(36*ij0-12)= AUlu0(36*ij0-12) - RHS_Aij(4,6)
              AUlu0(36*ij0-11)= AUlu0(36*ij0-11) - RHS_Aij(5,1)
              AUlu0(36*ij0-10)= AUlu0(36*ij0-10) - RHS_Aij(5,2)
              AUlu0(36*ij0-9 )= AUlu0(36*ij0-9 ) - RHS_Aij(5,3)
              AUlu0(36*ij0-8 )= AUlu0(36*ij0-8 ) - RHS_Aij(5,4)
              AUlu0(36*ij0-7 )= AUlu0(36*ij0-7 ) - RHS_Aij(5,5)
              AUlu0(36*ij0-6 )= AUlu0(36*ij0-6 ) - RHS_Aij(5,6)
              AUlu0(36*ij0-5 )= AUlu0(36*ij0-5 ) - RHS_Aij(6,1)
              AUlu0(36*ij0-4 )= AUlu0(36*ij0-4 ) - RHS_Aij(6,2)
              AUlu0(36*ij0-3 )= AUlu0(36*ij0-3 ) - RHS_Aij(6,3)
              AUlu0(36*ij0-2 )= AUlu0(36*ij0-2 ) - RHS_Aij(6,4)
              AUlu0(36*ij0-1 )= AUlu0(36*ij0-1 ) - RHS_Aij(6,5)
              AUlu0(36*ij0   )= AUlu0(36*ij0   ) - RHS_Aij(6,6)
            endif
          endif
        enddo
      enddo
    enddo

    deallocate (IW1, IW2)
    deallocate (iconFI1L, iconFI1U)
    !C===
  end subroutine FORM_ILU2_66


  !C
  !C***
  !C*** fill_in_S66_SORT
  !C***
  !C
  subroutine fill_in_S66_SORT (STEM, INUM, N, NP)
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

  end subroutine fill_in_S66_SORT

  !C
  !C***
  !C*** ILU1a66
  !C***
  !C
  !C    computes LU factorization of 3*3 Diagonal Block
  !C
  subroutine ILU1a66 (ALU, D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26, &
      & D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56, &
      & D61,D62,D63,D64,D65,D66)
    use hecmw_util
    implicit none
    real(kind=kreal) :: ALU(6,6), PW(6)
    real(kind=kreal) :: D11,D12,D13,D14,D15,D16
    real(kind=kreal) :: D21,D22,D23,D24,D25,D26
    real(kind=kreal) :: D31,D32,D33,D34,D35,D36
    real(kind=kreal) :: D41,D42,D43,D44,D45,D46
    real(kind=kreal) :: D51,D52,D53,D54,D55,D56
    real(kind=kreal) :: D61,D62,D63,D64,D65,D66
    integer(kind=kint) :: i,j,k

    ALU(1,1)= D11
    ALU(1,2)= D12
    ALU(1,3)= D13
    ALU(1,4)= D14
    ALU(1,5)= D15
    ALU(1,6)= D16
    ALU(2,1)= D21
    ALU(2,2)= D22
    ALU(2,3)= D23
    ALU(2,4)= D24
    ALU(2,5)= D25
    ALU(2,6)= D26
    ALU(3,1)= D31
    ALU(3,2)= D32
    ALU(3,3)= D33
    ALU(3,4)= D34
    ALU(3,5)= D35
    ALU(3,6)= D36
    ALU(4,1)= D41
    ALU(4,2)= D42
    ALU(4,3)= D43
    ALU(4,4)= D44
    ALU(4,5)= D45
    ALU(4,6)= D46
    ALU(5,1)= D51
    ALU(5,2)= D52
    ALU(5,3)= D53
    ALU(5,4)= D54
    ALU(5,5)= D55
    ALU(5,6)= D56
    ALU(6,1)= D61
    ALU(6,2)= D62
    ALU(6,3)= D63
    ALU(6,4)= D64
    ALU(6,5)= D65
    ALU(6,6)= D66

    do k= 1, 6
      ALU(k,k)= 1.d0/ALU(k,k)
      do i= k+1, 6
        ALU(i,k)= ALU(i,k) * ALU(k,k)
        do j= k+1, 6
          PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
        enddo
        do j= k+1, 6
          ALU(i,j)= PW(j)
        enddo
      enddo
    enddo

    return
  end subroutine ILU1a66

  !C
  !C***
  !C*** ILU1b66
  !C***
  !C
  !C    computes L_ik * D_k_INV * U_kj at ILU factorization
  !C    for 3*3 Block Type Matrix
  !C
  subroutine ILU1b66 (RHS_Aij, DkINV, Aik, Akj)
    use hecmw_util
    implicit none
    real(kind=kreal) :: RHS_Aij(6,6), DkINV(6,6), Aik(6,6), Akj(6,6)
    real(kind=kreal) :: X1,X2,X3,X4,X5,X6

    !C
    !C-- 1st Col.
    X1= Akj(1,1)
    X2= Akj(2,1)
    X3= Akj(3,1)
    X4= Akj(4,1)
    X5= Akj(5,1)
    X6= Akj(6,1)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,1)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,1)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,1)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,1)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,1)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,1)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    !C
    !C-- 2nd Col.
    X1= Akj(1,2)
    X2= Akj(2,2)
    X3= Akj(3,2)
    X4= Akj(4,2)
    X5= Akj(5,2)
    X6= Akj(6,2)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,2)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,2)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,2)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,2)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,2)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,2)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    !C
    !C-- 3rd Col.
    X1= Akj(1,3)
    X2= Akj(2,3)
    X3= Akj(3,3)
    X4= Akj(4,3)
    X5= Akj(5,3)
    X6= Akj(6,3)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,3)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,3)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,3)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,3)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,3)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,3)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    !C
    !C-- 4th Col.
    X1= Akj(1,4)
    X2= Akj(2,4)
    X3= Akj(3,4)
    X4= Akj(4,4)
    X5= Akj(5,4)
    X6= Akj(6,4)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,4)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,4)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,4)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,4)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,4)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,4)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    !C
    !C-- 5th Col.
    X1= Akj(1,5)
    X2= Akj(2,5)
    X3= Akj(3,5)
    X4= Akj(4,5)
    X5= Akj(5,5)
    X6= Akj(6,5)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,5)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,5)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,5)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,5)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,5)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,5)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    !C
    !C-- 6th Col.
    X1= Akj(1,6)
    X2= Akj(2,6)
    X3= Akj(3,6)
    X4= Akj(4,6)
    X5= Akj(5,6)
    X6= Akj(6,6)

    X2= X2 -DkINV(2,1)*X1
    X3= X3 -DkINV(3,1)*X1 -DkINV(3,2)*X2
    X4= X4 -DkINV(4,1)*X1 -DkINV(4,2)*X2 -DkINV(4,3)*X3
    X5= X5 -DkINV(5,1)*X1 -DkINV(5,2)*X2 -DkINV(5,3)*X3 -DkINV(5,4)*X4
    X6= X6 -DkINV(6,1)*X1 -DkINV(6,2)*X2 -DkINV(6,3)*X3 -DkINV(6,4)*X4 -DkINV(6,5)*X5

    X6= DkINV(6,6)*  X6
    X5= DkINV(5,5)*( X5 -DkINV(5,6)*X6 )
    X4= DkINV(4,4)*( X4 -DkINV(4,6)*X6 -DkINV(4,5)*X5)
    X3= DkINV(3,3)*( X3 -DkINV(3,6)*X6 -DkINV(3,5)*X5 -DkINV(3,4)*X4)
    X2= DkINV(2,2)*( X2 -DkINV(2,6)*X6 -DkINV(2,5)*X5 -DkINV(2,4)*X4 -DkINV(2,3)*X3)
    X1= DkINV(1,1)*( X1 -DkINV(1,6)*X6 -DkINV(1,5)*X5 -DkINV(1,4)*X4 -DkINV(1,3)*X3 -DkINV(1,2)*X2)

    RHS_Aij(1,6)=  Aik(1,1)*X1 +Aik(1,2)*X2 +Aik(1,3)*X3 +Aik(1,4)*X4 +Aik(1,5)*X5 +Aik(1,6)*X6
    RHS_Aij(2,6)=  Aik(2,1)*X1 +Aik(2,2)*X2 +Aik(2,3)*X3 +Aik(2,4)*X4 +Aik(2,5)*X5 +Aik(2,6)*X6
    RHS_Aij(3,6)=  Aik(3,1)*X1 +Aik(3,2)*X2 +Aik(3,3)*X3 +Aik(3,4)*X4 +Aik(3,5)*X5 +Aik(3,6)*X6
    RHS_Aij(4,6)=  Aik(4,1)*X1 +Aik(4,2)*X2 +Aik(4,3)*X3 +Aik(4,4)*X4 +Aik(4,5)*X5 +Aik(4,6)*X6
    RHS_Aij(5,6)=  Aik(5,1)*X1 +Aik(5,2)*X2 +Aik(5,3)*X3 +Aik(5,4)*X4 +Aik(5,5)*X5 +Aik(5,6)*X6
    RHS_Aij(6,6)=  Aik(6,1)*X1 +Aik(6,2)*X2 +Aik(6,3)*X3 +Aik(6,4)*X4 +Aik(6,5)*X5 +Aik(6,6)*X6

    return
  end subroutine ILU1b66

end module     hecmw_precond_BILU_66
