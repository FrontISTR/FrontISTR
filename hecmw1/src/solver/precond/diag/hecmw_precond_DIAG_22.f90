!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_22
!C***
!C
module hecmw_precond_DIAG_22
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_22_setup
  public:: hecmw_precond_DIAG_22_apply
  public:: hecmw_precond_DIAG_22_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_DIAG_22_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(2,2), PW(2)
    integer(kind=kint ):: ii, i, j, k

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 0 .and. hecMAT%Iarray(97) == 0) return
      call hecmw_precond_DIAG_22_clear()
    endif

    N = hecMAT%N
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(4*NP))
    ALU = 0.d0

    do ii= 1, N
      ALU(4*ii-3) = D(4*ii-3)
      ALU(4*ii-2) = D(4*ii-2)
      ALU(4*ii-1) = D(4*ii-1)
      ALU(4*ii-0) = D(4*ii-0)
    enddo

    if (hecMAT%cmat%n_val.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = hecMAT%cmat%pair(k)%i
        ALU(4*ii-3) = ALU(4*ii-3) + hecMAT%cmat%pair(k)%val(1, 1)
        ALU(4*ii-2) = ALU(4*ii-2) + hecMAT%cmat%pair(k)%val(1, 2)
        ALU(4*ii-1) = ALU(4*ii-1) + hecMAT%cmat%pair(k)%val(2, 1)
        ALU(4*ii-0) = ALU(4*ii-0) + hecMAT%cmat%pair(k)%val(2, 2)
      enddo

      !call hecmw_cmat_LU( hecMAT )

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

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_DIAG_22_setup

  subroutine hecmw_precond_DIAG_22_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2

    !C
    !C== Block SCALING

    !$omp parallel default(none),private(i,X1,X2),shared(N,WW,ALU)
    !$omp do
    do i= 1, N
      X1= WW(2*i-1)
      X2= WW(2*i-0)
      X2= X2 - ALU(4*i-1)*X1
      X2= ALU(4*i  )*  X2
      X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2 )
      WW(2*i-1)= X1
      WW(2*i-0)= X2
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_22_apply

  subroutine hecmw_precond_DIAG_22_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_DIAG_22_clear

end module     hecmw_precond_DIAG_22
