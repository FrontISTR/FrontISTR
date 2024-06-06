!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_44
!C***
!C
module hecmw_precond_DIAG_44
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_44_setup
  public:: hecmw_precond_DIAG_44_apply
  public:: hecmw_precond_DIAG_44_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_DIAG_44_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint ) :: NP
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(4,4), PW(4)
    integer(kind=kint ):: ii, i, j, k

    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 0 .and. hecMAT%Iarray(97) == 0) return
      call hecmw_precond_DIAG_44_clear()
    endif

    N = hecMAT%N
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(16*NP))
    ALU = 0.d0

    do ii= 1, N
      ALU(16*ii-15) = D(16*ii-15)
      ALU(16*ii-14) = D(16*ii-14)
      ALU(16*ii-13) = D(16*ii-13)
      ALU(16*ii-12) = D(16*ii-12)
      ALU(16*ii-11) = D(16*ii-11)
      ALU(16*ii-10) = D(16*ii-10)
      ALU(16*ii-9 ) = D(16*ii-9 )
      ALU(16*ii-8 ) = D(16*ii-8 )
      ALU(16*ii-7 ) = D(16*ii-7 )
      ALU(16*ii-6 ) = D(16*ii-6 )
      ALU(16*ii-5 ) = D(16*ii-5 )
      ALU(16*ii-4 ) = D(16*ii-4 )
      ALU(16*ii-3 ) = D(16*ii-3 )
      ALU(16*ii-2 ) = D(16*ii-2 )
      ALU(16*ii-1 ) = D(16*ii-1 )
      ALU(16*ii   ) = D(16*ii   )
    enddo

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N
      ALUtmp(1,1)= ALU(16*ii-15) * SIGMA_DIAG
      ALUtmp(1,2)= ALU(16*ii-14)
      ALUtmp(1,3)= ALU(16*ii-13)
      ALUtmp(1,4)= ALU(16*ii-12)

      ALUtmp(2,1)= ALU(16*ii-11)
      ALUtmp(2,2)= ALU(16*ii-10) * SIGMA_DIAG
      ALUtmp(2,3)= ALU(16*ii- 9)
      ALUtmp(2,4)= ALU(16*ii- 8)

      ALUtmp(3,1)= ALU(16*ii- 7)
      ALUtmp(3,2)= ALU(16*ii- 6)
      ALUtmp(3,3)= ALU(16*ii- 5) * SIGMA_DIAG
      ALUtmp(3,4)= ALU(16*ii- 4)

      ALUtmp(4,1)= ALU(16*ii- 3)
      ALUtmp(4,2)= ALU(16*ii- 2)
      ALUtmp(4,3)= ALU(16*ii- 1)
      ALUtmp(4,4)= ALU(16*ii   ) * SIGMA_DIAG

      do k= 1, 4
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, 4
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, 4
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, 4
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo

      ALU(16*ii-15)= ALUtmp(1,1)
      ALU(16*ii-14)= ALUtmp(1,2)
      ALU(16*ii-13)= ALUtmp(1,3)
      ALU(16*ii-12)= ALUtmp(1,4)
      ALU(16*ii-11)= ALUtmp(2,1)
      ALU(16*ii-10)= ALUtmp(2,2)
      ALU(16*ii- 9)= ALUtmp(2,3)
      ALU(16*ii- 8)= ALUtmp(2,4)
      ALU(16*ii- 7)= ALUtmp(3,1)
      ALU(16*ii- 6)= ALUtmp(3,2)
      ALU(16*ii- 5)= ALUtmp(3,3)
      ALU(16*ii- 4)= ALUtmp(3,4)
      ALU(16*ii- 3)= ALUtmp(4,1)
      ALU(16*ii- 2)= ALUtmp(4,2)
      ALU(16*ii- 1)= ALUtmp(4,3)
      ALU(16*ii   )= ALUtmp(4,4)
    enddo
    !$omp end do
    !$omp end parallel

    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done
  end subroutine hecmw_precond_DIAG_44_setup

  subroutine hecmw_precond_DIAG_44_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2, X3, X4

    !C
    !C== Block SCALING
    !$omp parallel default(none),private(i,X1,X2,X3,X4),shared(N,WW,ALU)
    !$omp do
    do i= 1, N
      X1= WW(4*i-3)
      X2= WW(4*i-2)
      X3= WW(4*i-1)
      X4= WW(4*i  )
      X2= X2 -ALU(16*i-11)*X1
      X3= X3 -ALU(16*i- 7)*X1 -ALU(16*i- 6)*X2
      X4= X4 -ALU(16*i- 3)*X1 -ALU(16*i- 2)*X2 -ALU(16*i- 1)*X3
      X4= ALU(16*i   )* X4
      X3= ALU(16*i- 5)*(X3 -ALU(16*i- 4)*X4)
      X2= ALU(16*i-10)*(X2 -ALU(16*i- 8)*X4 -ALU(16*i- 9)*X3)
      X1= ALU(16*i-15)*(X1 -ALU(16*i-12)*X4 -ALU(16*i-13)*X3 -ALU(16*i-14)*X2)
      WW(4*i-3)= X1
      WW(4*i-2)= X2
      WW(4*i-1)= X3
      WW(4*i  )= X4
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_44_apply

  subroutine hecmw_precond_DIAG_44_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
    INITIALIZED = .false.
  end subroutine hecmw_precond_DIAG_44_clear

end module     hecmw_precond_DIAG_44
