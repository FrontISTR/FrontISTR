!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_nn
!C***
!C
module hecmw_precond_DIAG_nn
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_nn_setup
  public:: hecmw_precond_DIAG_nn_apply
  public:: hecmw_precond_DIAG_nn_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_DIAG_nn_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP,NDOF,NDOF2
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(hecMAT%NDOF,hecMAT%NDOF), PW(hecMAT%NDOF)
    integer(kind=kint ):: i, j, k, ii

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 0 .and. hecMAT%Iarray(97) == 0) return
      call hecmw_precond_DIAG_nn_clear()
    endif

    N = hecMAT%N
    NDOF = hecMAT%NDOF
    NDOF2 = NDOF*NDOF
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(NDOF2*NP))
    ALU = 0.d0

    do ii= 1, N
      do i = 1, NDOF2
        ALU(NDOF2*(ii-1)+i)=D(NDOF2*(ii-1)+i)
      end do
    enddo

    if (hecMAT%cmat%n_val.gt.0) then
      do k= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
        ii = hecMAT%cmat%pair(k)%i
        do i = 1,NDOF
          do j = 1,NDOF
            ALU(NDOF2*(ii-1)+(i-1)*NDOF+j) = ALU(NDOF2*(ii-1)+(i-1)*NDOF+j) + hecMAT%cmat%pair(k)%val(i, j)
          enddo
        enddo
      enddo
    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,NDOF,NDOF2,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N

      do i = 1, NDOF
        do j =  1, NDOF
          ALUtmp(i,j) = ALU(NDOF2*(ii-1)+(i-1)*NDOF+j)
          if (i==j) ALUtmp(i,j)=ALUtmp(i,j)*SIGMA_DIAG
        end do
      end do
      do k= 1, NDOF
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, NDOF
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, NDOF
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, NDOF
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo
      do i = 1, NDOF
        do j =  1, NDOF
          ALU(NDOF2*(ii-1)+(i-1)*NDOF+j)= ALUtmp(i,j)
        end do
      end do
    enddo
    !$omp end do
    !$omp end parallel
    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done

  end subroutine hecmw_precond_DIAG_nn_setup

  subroutine hecmw_precond_DIAG_nn_apply(WW,NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i,j,k,NDOF
    real(kind=kreal) :: X(NDOF)

    !C
    !C== Block SCALING
    !$omp parallel default(none),private(i,j,k,X),shared(N,WW,ALU,NDOF)
    !$omp do
    do i= 1, N
      do j=1,NDOF
        X(j)=WW(NDOF*(i-1)+j)
      end do
      do j=2,NDOF
        do k = 1,j-1
          X(j)=X(j)-ALU(NDOF*NDOF*(i-1)+NDOF*(j-1)+k )*X(k)
        end do
      end do
      do j=NDOF,1,-1
        do k = NDOF,j+1,-1
          X(j)=X(j)-ALU(NDOF*NDOF*(i-1)+NDOF*(j-1)+k )*X(k)
        end do
        X(j)=ALU(NDOF*NDOF*(i-1)+(NDOF+1)*(j-1)+1 )*X(j)
      end do
      do j=1,NDOF
        WW(NDOF*(i-1)+j)=X(j)
      end do

    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_nn_apply

  subroutine hecmw_precond_DIAG_nn_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_DIAG_nn_clear

end module     hecmw_precond_DIAG_nn
