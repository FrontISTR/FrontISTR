!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_33
!C***
!C
module hecmw_precond_DIAG_33
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_33_setup
  public:: hecmw_precond_DIAG_33_apply
  public:: hecmw_precond_DIAG_33_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_DIAG_33_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(3,3), PW(3)
    integer(kind=kint ):: ii, i, j, k

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 0 .and. hecMAT%Iarray(97) == 0) return
      call hecmw_precond_DIAG_33_clear()
    endif

    N = hecMAT%N
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(9*NP))
    ALU = 0.d0

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

      !call hecmw_cmat_LU( hecMAT )

    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,SIGMA_DIAG)
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
    !$omp end do
    !$omp end parallel

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_DIAG_33_setup

  subroutine hecmw_precond_DIAG_33_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2, X3

    !C
    !C== Block SCALING

    !$omp parallel default(none),private(i,X1,X2,X3),shared(N,WW,ALU)
    !$omp do
    do i= 1, N
      X1= WW(3*i-2)
      X2= WW(3*i-1)
      X3= WW(3*i  )
      X2= X2 - ALU(9*i-5)*X1
      X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
      X3= ALU(9*i  )*  X3
      X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
      X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
      WW(3*i-2)= X1
      WW(3*i-1)= X2
      WW(3*i  )= X3
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_33_apply

  subroutine hecmw_precond_DIAG_33_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_DIAG_33_clear

end module     hecmw_precond_DIAG_33
