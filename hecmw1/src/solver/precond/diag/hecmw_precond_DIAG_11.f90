!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_11
!C***
!C
module hecmw_precond_DIAG_11
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_11_setup
  public:: hecmw_precond_DIAG_11_apply
  public:: hecmw_precond_DIAG_11_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_DIAG_11_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(1,1), PW(1)
    integer(kind=kint ):: ii, i, j, k

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 0 .and. hecMAT%Iarray(97) == 0) return
      call hecmw_precond_DIAG_11_clear()
    endif

    N = hecMAT%N
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(NP))
    ALU = 0.d0

    do ii= 1, N
      ALU(ii) = D(ii)
    enddo

    if (hecMAT%cmat%n_val.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = hecMAT%cmat%pair(k)%i
        ALU(ii) = ALU(ii) + hecMAT%cmat%pair(k)%val(1, 1)
      enddo

      !call hecmw_cmat_LU( hecMAT )

    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N
      ALUtmp(1,1)= ALU(ii) * SIGMA_DIAG
      ALUtmp(1,1)= 1.d0/ALUtmp(1,1)
      ALU(ii)= ALUtmp(1,1)
    enddo
    !$omp end do
    !$omp end parallel

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_DIAG_11_setup

  subroutine hecmw_precond_DIAG_11_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i

    !C
    !C== Block SCALING

    !$omp parallel default(none),private(i),shared(N,WW,ALU)
    !$omp do
    do i= 1, N
      WW(i)= ALU(i)*WW(i)
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_11_apply

  subroutine hecmw_precond_DIAG_11_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_DIAG_11_clear

end module     hecmw_precond_DIAG_11
