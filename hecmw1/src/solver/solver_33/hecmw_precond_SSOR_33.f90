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

  private

  public:: hecmw_precond_SSOR_33_setup
  public:: hecmw_precond_SSOR_33_apply
  public:: hecmw_precond_SSOR_33_clear

  real(kind=kreal), pointer :: ALU(:) => null()
  real(kind=kreal), pointer :: AU(:) => null()
  real(kind=kreal), pointer :: AL(:) => null()
  real(kind=kreal), pointer :: D(:) => null()
  integer(kind=kint), pointer :: indexL(:) => null()
  integer(kind=kint), pointer :: itemL(:) => null()
  integer(kind=kint), pointer :: indexU(:) => null()
  integer(kind=kint), pointer :: itemU(:) => null()

  integer(kind=kint) :: NC = 0
  real(kind=kreal), pointer :: CAU(:) => null()
  real(kind=kreal), pointer :: CAL(:) => null()
  integer(kind=kint), pointer :: indexCL(:) => null()
  integer(kind=kint), pointer :: itemCL(:) => null()
  integer(kind=kint), pointer :: indexCU(:) => null()
  integer(kind=kint), pointer :: itemCU(:) => null()

contains

  subroutine hecmw_precond_SSOR_33_setup(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint ) :: N, NP
    real   (kind=kreal) :: SIGMA_DIAG
    real   (kind=kreal) :: ALUtmp(3,3), PW(3)
    integer(kind=kint ) :: ii, i, j, k

    N = hecMAT%N
    NP = hecMAT%NP
    AU => hecMAT%AU
    AL => hecMAT%AL
    D  => hecMAT%D
    indexL => hecMAT%indexL
    itemL  => hecMAT%itemL
    indexU => hecMAT%indexU
    itemU  => hecMAT%itemU
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(9*NP))
    ALU  = 0.d0

    do ii= 1, N
      ALU(9*ii-8) = D(9*ii-8)
      ALU(9*ii-7) = D(9*ii-7)
      ALU(9*ii-6) = D(9*ii-6)
      ALU(9*ii-5) = D(9*ii-5)
      ALU(9*ii-4) = D(9*ii-4)
      ALU(9*ii-3) = D(9*ii-3)
      ALU(9*ii-2) = D(9*ii-2)
      ALU(9*ii-1) = D(9*ii-1)
      ALU(9*ii  ) = D(9*ii  )
    enddo

    if (hecMAT%cmat%n_val.gt.0) then
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

      call hecmw_cmat_LU( hecMAT )

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

    NC = hecMAT%cmat%n_val
    CAL => hecMAT%CAL
    CAU => hecMAT%CAU
    indexCL => hecMAT%indexCL
    itemCL => hecMAT%itemCL
    indexCU => hecMAT%indexCU
    itemCU => hecMAT%itemCU

  end subroutine hecmw_precond_SSOR_33_setup

  subroutine hecmw_precond_SSOR_33_apply(N, ZP)
    implicit none
    integer(kind=kint) :: N
    real(kind=kreal) :: ZP(:)
    integer(kind=kint) :: i, j, isL, ieL, isU, ieU, k
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3

    !C-- FORWARD

    do i= 1, N
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
      enddo

      if (NC.ne.0) then
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
        enddo
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
    enddo

    !C-- BACKWARD

    do i= N, 1, -1
      isU= indexU(i-1) + 1
      ieU= indexU(i)
      SW1= 0.d0
      SW2= 0.d0
      SW3= 0.d0
      do j= ieU, isU, -1
        k= itemU(j)
        X1= ZP(3*k-2)
        X2= ZP(3*k-1)
        X3= ZP(3*k  )
        SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
        SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
        SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
      enddo

      if (NC.gt.0) then
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
        enddo
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
    enddo
  end subroutine hecmw_precond_SSOR_33_apply

  subroutine hecmw_precond_SSOR_33_clear(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    nullify(AU)
    nullify(AL)
    nullify(D)
    nullify(indexL)
    nullify(itemL)
    nullify(indexU)
    nullify(itemU)

    if (NC.ne.0) then
      nullify(CAU)
      nullify(CAL)
      nullify(indexCL)
      nullify(itemCL)
      nullify(indexCU)
      nullify(itemCU)
      call hecmw_cmat_LU_free( hecMAT )
    endif
    NC = 0
  end subroutine hecmw_precond_SSOR_33_clear

end module     hecmw_precond_SSOR_33
