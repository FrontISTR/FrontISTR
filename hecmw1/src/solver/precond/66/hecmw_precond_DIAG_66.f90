!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_DIAG_66
!C***
!C
module hecmw_precond_DIAG_66
  use  hecmw_util

  private

  public:: hecmw_precond_DIAG_66_setup
  public:: hecmw_precond_DIAG_66_apply
  public:: hecmw_precond_DIAG_66_clear

  integer(kind=kint) :: N
  real(kind=kreal), pointer :: ALU(:) => null()

contains

  subroutine hecmw_precond_DIAG_66_setup(hecMAT)
    use hecmw_matrix_misc
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint ) :: NP
    real   (kind=kreal) :: SIGMA_DIAG
    real(kind=kreal), pointer:: D(:)

    real   (kind=kreal):: ALUtmp(6,6), PW(6)
    integer(kind=kint ):: ii, i, j, k

    N = hecMAT%N
    NP = hecMAT%NP
    D => hecMAT%D
    SIGMA_DIAG = hecmw_mat_get_sigma_diag(hecMAT)

    allocate(ALU(36*NP))
    ALU = 0.d0

    do ii= 1, N
      ALU(36*ii-35) = D(36*ii-35)
      ALU(36*ii-34) = D(36*ii-34)
      ALU(36*ii-33) = D(36*ii-33)
      ALU(36*ii-32) = D(36*ii-32)
      ALU(36*ii-31) = D(36*ii-31)
      ALU(36*ii-30) = D(36*ii-30)
      ALU(36*ii-29) = D(36*ii-29)
      ALU(36*ii-28) = D(36*ii-28)
      ALU(36*ii-27) = D(36*ii-27)
      ALU(36*ii-26) = D(36*ii-26)
      ALU(36*ii-25) = D(36*ii-25)
      ALU(36*ii-24) = D(36*ii-24)
      ALU(36*ii-23) = D(36*ii-23)
      ALU(36*ii-22) = D(36*ii-22)
      ALU(36*ii-21) = D(36*ii-21)
      ALU(36*ii-20) = D(36*ii-20)
      ALU(36*ii-19) = D(36*ii-19)
      ALU(36*ii-18) = D(36*ii-18)
      ALU(36*ii-17) = D(36*ii-17)
      ALU(36*ii-16) = D(36*ii-16)
      ALU(36*ii-15) = D(36*ii-15)
      ALU(36*ii-14) = D(36*ii-14)
      ALU(36*ii-13) = D(36*ii-13)
      ALU(36*ii-12) = D(36*ii-12)
      ALU(36*ii-11) = D(36*ii-11)
      ALU(36*ii-10) = D(36*ii-10)
      ALU(36*ii-9 ) = D(36*ii-9 )
      ALU(36*ii-8 ) = D(36*ii-8 )
      ALU(36*ii-7 ) = D(36*ii-7 )
      ALU(36*ii-6 ) = D(36*ii-6 )
      ALU(36*ii-5 ) = D(36*ii-5 )
      ALU(36*ii-4 ) = D(36*ii-4 )
      ALU(36*ii-3 ) = D(36*ii-3 )
      ALU(36*ii-2 ) = D(36*ii-2 )
      ALU(36*ii-1 ) = D(36*ii-1 )
      ALU(36*ii   ) = D(36*ii   )
    enddo

    !    if (hecMAT%cmat%n_val.gt.0) then
    !      do k= 1, hecMAT%cmat%n_val
    !        if (hecMAT%cmat%pair(k)%i.ne.hecMAT%cmat%pair(k)%j) cycle
    !        ii = hecMAT%cmat%pair(k)%i
    !        ALU(9*ii-8) = ALU(9*ii-8) + hecMAT%cmat%pair(k)%val(1, 1)
    !        ALU(9*ii-7) = ALU(9*ii-7) + hecMAT%cmat%pair(k)%val(1, 2)
    !        ALU(9*ii-6) = ALU(9*ii-6) + hecMAT%cmat%pair(k)%val(1, 3)
    !        ALU(9*ii-5) = ALU(9*ii-5) + hecMAT%cmat%pair(k)%val(2, 1)
    !        ALU(9*ii-4) = ALU(9*ii-4) + hecMAT%cmat%pair(k)%val(2, 2)
    !        ALU(9*ii-3) = ALU(9*ii-3) + hecMAT%cmat%pair(k)%val(2, 3)
    !        ALU(9*ii-2) = ALU(9*ii-2) + hecMAT%cmat%pair(k)%val(3, 1)
    !        ALU(9*ii-1) = ALU(9*ii-1) + hecMAT%cmat%pair(k)%val(3, 2)
    !        ALU(9*ii  ) = ALU(9*ii  ) + hecMAT%cmat%pair(k)%val(3, 3)
    !      enddo
    !
    !      !call hecmw_cmat_LU( hecMAT )
    !
    !    endif

    !$omp parallel default(none),private(ii,ALUtmp,k,i,j,PW),shared(N,ALU,SIGMA_DIAG)
    !$omp do
    do ii= 1, N
      ALUtmp(1,1)= ALU(36*ii-35) * SIGMA_DIAG
      ALUtmp(1,2)= ALU(36*ii-34)
      ALUtmp(1,3)= ALU(36*ii-33)
      ALUtmp(1,4)= ALU(36*ii-32)
      ALUtmp(1,5)= ALU(36*ii-31)
      ALUtmp(1,6)= ALU(36*ii-30)

      ALUtmp(2,1)= ALU(36*ii-29)
      ALUtmp(2,2)= ALU(36*ii-28) * SIGMA_DIAG
      ALUtmp(2,3)= ALU(36*ii-27)
      ALUtmp(2,4)= ALU(36*ii-26)
      ALUtmp(2,5)= ALU(36*ii-25)
      ALUtmp(2,6)= ALU(36*ii-24)

      ALUtmp(3,1)= ALU(36*ii-23)
      ALUtmp(3,2)= ALU(36*ii-22)
      ALUtmp(3,3)= ALU(36*ii-21) * SIGMA_DIAG
      ALUtmp(3,4)= ALU(36*ii-20)
      ALUtmp(3,5)= ALU(36*ii-19)
      ALUtmp(3,6)= ALU(36*ii-18)

      ALUtmp(4,1)= ALU(36*ii-17)
      ALUtmp(4,2)= ALU(36*ii-16)
      ALUtmp(4,3)= ALU(36*ii-15)
      ALUtmp(4,4)= ALU(36*ii-14) * SIGMA_DIAG
      ALUtmp(4,5)= ALU(36*ii-13)
      ALUtmp(4,6)= ALU(36*ii-12)

      ALUtmp(5,1)= ALU(36*ii-11)
      ALUtmp(5,2)= ALU(36*ii-10)
      ALUtmp(5,3)= ALU(36*ii-9 )
      ALUtmp(5,4)= ALU(36*ii-8 )
      ALUtmp(5,5)= ALU(36*ii-7 ) * SIGMA_DIAG
      ALUtmp(5,6)= ALU(36*ii-6 )

      ALUtmp(6,1)= ALU(36*ii-5 )
      ALUtmp(6,2)= ALU(36*ii-4 )
      ALUtmp(6,3)= ALU(36*ii-3 )
      ALUtmp(6,4)= ALU(36*ii-2 )
      ALUtmp(6,5)= ALU(36*ii-1 )
      ALUtmp(6,6)= ALU(36*ii   ) * SIGMA_DIAG

      do k= 1, 6
        ALUtmp(k,k)= 1.d0/ALUtmp(k,k)
        do i= k+1, 6
          ALUtmp(i,k)= ALUtmp(i,k) * ALUtmp(k,k)
          do j= k+1, 6
            PW(j)= ALUtmp(i,j) - ALUtmp(i,k)*ALUtmp(k,j)
          enddo
          do j= k+1, 6
            ALUtmp(i,j)= PW(j)
          enddo
        enddo
      enddo

      ALU(36*ii-35)= ALUtmp(1,1)
      ALU(36*ii-34)= ALUtmp(1,2)
      ALU(36*ii-33)= ALUtmp(1,3)
      ALU(36*ii-32)= ALUtmp(1,4)
      ALU(36*ii-31)= ALUtmp(1,5)
      ALU(36*ii-30)= ALUtmp(1,6)
      ALU(36*ii-29)= ALUtmp(2,1)
      ALU(36*ii-28)= ALUtmp(2,2)
      ALU(36*ii-27)= ALUtmp(2,3)
      ALU(36*ii-26)= ALUtmp(2,4)
      ALU(36*ii-25)= ALUtmp(2,5)
      ALU(36*ii-24)= ALUtmp(2,6)
      ALU(36*ii-23)= ALUtmp(3,1)
      ALU(36*ii-22)= ALUtmp(3,2)
      ALU(36*ii-21)= ALUtmp(3,3)
      ALU(36*ii-20)= ALUtmp(3,4)
      ALU(36*ii-19)= ALUtmp(3,5)
      ALU(36*ii-18)= ALUtmp(3,6)
      ALU(36*ii-17)= ALUtmp(4,1)
      ALU(36*ii-16)= ALUtmp(4,2)
      ALU(36*ii-15)= ALUtmp(4,3)
      ALU(36*ii-14)= ALUtmp(4,4)
      ALU(36*ii-13)= ALUtmp(4,5)
      ALU(36*ii-12)= ALUtmp(4,6)
      ALU(36*ii-11)= ALUtmp(5,1)
      ALU(36*ii-10)= ALUtmp(5,2)
      ALU(36*ii-9 )= ALUtmp(5,3)
      ALU(36*ii-8 )= ALUtmp(5,4)
      ALU(36*ii-7 )= ALUtmp(5,5)
      ALU(36*ii-6 )= ALUtmp(5,6)
      ALU(36*ii-5 )= ALUtmp(6,1)
      ALU(36*ii-4 )= ALUtmp(6,2)
      ALU(36*ii-3 )= ALUtmp(6,3)
      ALU(36*ii-2 )= ALUtmp(6,4)
      ALU(36*ii-1 )= ALUtmp(6,5)
      ALU(36*ii   )= ALUtmp(6,6)
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_66_setup

  subroutine hecmw_precond_DIAG_66_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2, X3, X4, X5, X6

    !C
    !C== Block SCALING
    !$omp parallel default(none),private(i,X1,X2,X3,X4,X5,X6),shared(N,WW,ALU)
    !$omp do
    do i= 1, N
      X1= WW(6*i-5)
      X2= WW(6*i-4)
      X3= WW(6*i-3)
      X4= WW(6*i-2)
      X5= WW(6*i-1)
      X6= WW(6*i  )
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
    !$omp end do
    !$omp end parallel

  end subroutine hecmw_precond_DIAG_66_apply

  subroutine hecmw_precond_DIAG_66_clear()
    implicit none
    if (associated(ALU)) deallocate(ALU)
    nullify(ALU)
  end subroutine hecmw_precond_DIAG_66_clear

end module     hecmw_precond_DIAG_66
