!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_BILU_nn
!C***
!C
module hecmw_precond_BILU_nn
  use hecmw_util
  use hecmw_matrix_misc

  private

  public:: hecmw_precond_BILU_nn_setup
  public:: hecmw_precond_BILU_nn_apply
  public:: hecmw_precond_BILU_nn_clear

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

  subroutine hecmw_precond_BILU_nn_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NDOF, NP, NPU, NPL
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
        call hecmw_precond_BILU_nn_clear()
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_BILU_nn_clear() ! TEMPORARY
      else
        return
      endif
    endif
    NDOF = hecMAT%NDOF
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

    if (PRECOND.eq.10) call FORM_ILU0_nn &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.11) call FORM_ILU1_nn &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    if (PRECOND.eq.12) call FORM_ILU2_nn &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_BILU_nn_setup

  subroutine hecmw_precond_BILU_nn_apply(WW,NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i, j, ii, ij, isL, ieL, isU, ieU, k, NDOF
    real(kind=kreal) :: SW(NDOF), X(NDOF)
    !C
    !C-- FORWARD

    do i= 1, N
      do ii = 1, NDOF
        SW(ii)= WW(NDOF*(i-1)+ii)
      end do
      isL= inumFI1L(i-1)+1
      ieL= inumFI1L(i)
      do j= isL, ieL
        k= FI1L(j)
        do ii = 1, NDOF
          X(ii)= WW(NDOF*(k-1)+ii)
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            SW(ii)= SW(ii) - ALlu0(NDOF*NDOF*(j-1)+NDOF*(ii-1)+ij)*X(ij)
          end do
        end do
      enddo

      X= SW
      do ii=2,NDOF
        do ij = 1,ii-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
      end do
      do ii=NDOF,1,-1
        do ij = NDOF,ii+1,-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
        X(ii)=Dlu0(NDOF*NDOF*(i-1)+(NDOF+1)*(ii-1)+1 )*X(ii)
      end do
      do ii = 1, NDOF
        WW(NDOF*(i-1)+ii)=X(ii)
      end do
    enddo

    !C
    !C-- BACKWARD

    do i= N, 1, -1
      isU= inumFI1U(i-1) + 1
      ieU= inumFI1U(i)
      SW= 0.d0

      do j= ieU, isU, -1
        k= FI1U(j)
        do ii = 1, NDOF
          X(ii)= WW(NDOF*(k-1)+ii)
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            SW(ii)= SW(ii) + AUlu0(NDOF*NDOF*(j-1)+NDOF*(ii-1)+ij)*X(ij)
          end do
        end do
      enddo
      X= SW
      do ii=2,NDOF
        do ij = 1,ii-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
      end do
      do ii=NDOF,1,-1
        do ij = NDOF,ii+1,-1
          X(ii)=X(ii)-Dlu0(NDOF*NDOF*(i-1)+NDOF*(ii-1)+ij )*X(ij)
        end do
        X(ii)=Dlu0(NDOF*NDOF*(i-1)+(NDOF+1)*(ii-1)+1 )*X(ii)
      end do
      do ii = 1, NDOF
        WW(NDOF*(i-1)+ii)= WW(NDOF*(i-1)+ii)-X(ii)
      end do
    enddo
  end subroutine hecmw_precond_BILU_nn_apply

  subroutine hecmw_precond_BILU_nn_clear()
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
  end subroutine hecmw_precond_BILU_nn_clear

  !C
  !C***
  !C*** FORM_ILU0_nn
  !C***
  !C
  !C    form ILU(0) matrix
  !C
  subroutine FORM_ILU0_nn                                   &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NDOF, NP, NPU, NPL
    real   (kind=kreal), intent(in):: SIGMA, SIGMA_DIAG

    real(kind=kreal), dimension(NDOF*NDOF*NPL), intent(in):: AL
    real(kind=kreal), dimension(NDOF*NDOF*NPU), intent(in):: AU
    real(kind=kreal), dimension(NDOF*NDOF*NP ), intent(in):: D

    integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
    integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
    integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

    integer(kind=kint), dimension(:), allocatable :: IW1, IW2
    real (kind=kreal),  dimension(NDOF,NDOF) :: RHS_Aij, DkINV, Aik, Akj
    integer(kind=kint) :: i,jj,ij0,kk,NDOF2
    integer(kind=kint) :: j,k,ii,ij
    NDOF2=NDOF*NDOF
    allocate (IW1(NP) , IW2(NP))
    allocate(Dlu0(NDOF2*NP), ALlu0(NDOF2*NPL), AUlu0(NDOF2*NPU))
    allocate(inumFI1L(0:NP), inumFI1U(0:NP), FI1L(NPL), FI1U(NPU))

    do i=1,NDOF2*NP
      Dlu0(i) = D(i)
    end do
    do i=1,NDOF2*NPL
      ALlu0(i) = AL(i)
    end do
    do i=1,NDOF2*NPU
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
      do ii=1,NDOF
        Dlu0(NDOF2*(i-1)+(ii-1)*NDOF+ii)=Dlu0(NDOF2*(i-1)+(ii-1)*NDOF+ii)*SIGMA_DIAG
      end do
    enddo

    i = 1
    call ILU1aNN (DkINV, Dlu0(NDOF2*(i-1)+1:NDOF2*NDOF2),NDOF)
    do ii=1,NDOF
      do ij=1,NDOF
        Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij)= DkINV(ii,ij)
      end do
    end do

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
        do ii=1,NDOF
          do ij=1,NDOF
            DkINV(ii,ij) = Dlu0(NDOF2*(k-1)+NDOF*(ii-1)+ij)
          end do
        end do
        do ii=1,NDOF
          do ij=1,NDOF
            Aik(ii,ij) = ALlu0(NDOF2*(kk-1)+NDOF*(ii-1)+ij)
          end do
        end do

        do jj= INU(k-1)+1, INU(k)
          j= IAU(jj)
          if (IW1(j).eq.0.and.IW2(j).eq.0) cycle
          do ii=1,NDOF
            do ij=1,NDOF
              Akj(ii,ij) = AUlu0(NDOF2*(jj-1)+NDOF*(ii-1)+ij)
            end do
          end do

          call ILU1bNN (RHS_Aij, DkINV, Aik, Akj,NDOF)

          if (j.eq.i) then
            do ii=1,NDOF
              do ij=1,NDOF
                Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) = Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

          if (j.lt.i) then
            ij0= IW1(j)
            do ii=1,NDOF
              do ij=1,NDOF
                ALlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) = ALlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

          if (j.gt.i) then
            ij0= IW2(j)
            do ii=1,NDOF
              do ij=1,NDOF
                AUlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) = AUlu0(NDOF2*(ij0-1)+NDOF*(ii-1)+ij) - RHS_Aij(ii,ij)
              end do
            end do
          endif

        enddo
      enddo
      call ILU1aNN (DkINV, Dlu0(NDOF2*(i-1)+1:NDOF2*NDOF2),NDOF)

      do ii=1,NDOF
        do ij=1,NDOF
          Dlu0(NDOF2*(i-1)+NDOF*(ii-1)+ij) = DkINV(ii,ij)
        end do
      end do
    enddo

    deallocate (IW1, IW2)
  end subroutine FORM_ILU0_nn

  !C
  !C***
  !C*** FORM_ILU1_nn
  !C***
  !C
  !C    form ILU(1) matrix
  !C
  subroutine FORM_ILU1_nn                                   &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NDOF, NP, NPU, NPL
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
    integer(kind=kint) :: NPLf1,NPUf1,NDOF2
    integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
    integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3
    integer(kind=kint) :: j,k,iSL,iSU
    !C
    !C +--------------+
    !C | find fill-in |
    !C +--------------+
    !C===
    NDOF2=NDOF*NDOF
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

          call ILU1bNN (RHS_Aij, DkINV, Aik, Akj,3)

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
  end subroutine FORM_ILU1_nn

  !C
  !C***
  !C*** FORM_ILU2_nn
  !C***
  !C
  !C    form ILU(2) matrix
  !C
  subroutine FORM_ILU2_nn &
      &   (N, NDOF, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU, &
      &    SIGMA, SIGMA_DIAG)
    implicit none
    integer(kind=kint ), intent(in):: N, NDOF, NP, NPU, NPL
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
    integer(kind=kint) :: NPLf1,NPLf2,NPUf1,NPUf2,iAS,iconIK,iconKJ,NDOF2
    integer(kind=kint) :: i,jj,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
    integer(kind=kint) :: icou,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3
    integer(kind=kint) :: j,k,iSL,iSU

    !C
    !C +------------------+
    !C | find fill-in (1) |
    !C +------------------+
    !C===
    NDOF2=NDOF*NDOF
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

            call ILU1bNN (RHS_Aij, DkINV, Aik, Akj,3)

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
  end subroutine FORM_ILU2_nn


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
  !C*** ILU1aNN
  !C***
  !C
  !C    computes LU factorization of N*N Diagonal Block
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
  subroutine ILU1aNN (ALU, D, NDOF)
    use hecmw_util
    implicit none
    real(kind=kreal) :: ALU(NDOF,NDOF), D(NDOF*NDOF), PW(NDOF)
    integer(kind=kint) :: NDOF, i,j,k

    do i = 1, NDOF
      do j = 1, NDOF
        ALU(i,j) = D(NDOF*(i-1)+j)
      end do
    end do

    do k= 1, NDOF
      if (ALU(k,k) == 0.d0) then
        stop 'ERROR: Divide by zero in ILU setup'
      endif
      ALU(k,k)= 1.d0/ALU(k,k)
      do i= k+1, NDOF
        ALU(i,k)= ALU(i,k) * ALU(k,k)
        do j= k+1, NDOF
          PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
        enddo
        do j= k+1, NDOF
          ALU(i,j)= PW(j)
        enddo
      enddo
    enddo

    return
  end subroutine ILU1aNN
  !C
  !C***
  !C*** ILU1bNN
  !C***
  !C
  !C    computes L_ik * D_k_INV * U_kj at ILU factorization
  !C    for 3*3 Block Type Matrix
  !C
  subroutine ILU1bNN (RHS_Aij, DkINV, Aik, Akj, NDOF)
    use hecmw_util
    implicit none
    real(kind=kreal) :: RHS_Aij(NDOF,NDOF), DkINV(NDOF,NDOF), Aik(NDOF,NDOF), Akj(NDOF,NDOF)
    real(kind=kreal) :: X(NDOF)
    integer(kind=kint) :: NDOF,i,j,k

    do k=1,NDOF
      X(1:NDOF)= Akj(1:NDOF,k)
      do i=2,NDOF
        do j = 1,i-1
          X(i)=X(i)-DkINV(i,j)*X(j)
        end do
      end do
      do i=NDOF,1,-1
        do j = NDOF,i+1,-1
          X(i)=X(i)-DkINV(i,j)*X(j)
        end do
        X(i)=DkINV(i,i)*X(i)
      end do
      RHS_Aij(:,k)=0.0d0
      do i=1,NDOF
        do j=1,NDOF
          RHS_Aij(i,k) = RHS_Aij(i,k)+Aik(i,j)*X(j)
        end do
      end do
    end do
    return
  end subroutine ILU1bNN

end module     hecmw_precond_BILU_nn
