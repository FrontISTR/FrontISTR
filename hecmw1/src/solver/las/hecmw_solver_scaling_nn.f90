!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_scaling_nn
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  implicit none

  private
  real(kind=kreal), private, allocatable :: scale(:)

  public :: hecmw_solver_scaling_fw_nn
  public :: hecmw_solver_scaling_bk_nn

contains

  subroutine hecmw_solver_scaling_fw_nn(hecMESH, hecMAT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint) :: N, NP, NDOF, NDOF2
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i,j,k,ii,ij, ip(hecMAT%NDOF),iq(hecMAT%NDOF)
    integer(kind=kint) :: isL, ieL, isU, ieU, inod
    real(kind=kreal) :: START_TIME, END_TIME

    if (hecmw_mat_get_scaling(hecMAT).eq.0) return

    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
    NDOF2 = NDOF*NDOF
    D => hecMAT%D
    AL => hecMAT%AL
    AU => hecMAT%AU
    INL => hecMAT%indexL
    IAL => hecMAT%itemL
    INU => hecMAT%indexU
    IAU => hecMAT%itemU
    B => hecMAT%B

    allocate(scale(NDOF*NP))

    do i= 1, N
      do k=1, NDOF
        scale (NDOF*(i-1)+k)= 1.d0/dsqrt(dabs(D(NDOF*NDOF*(i-1)+(k-1)*(NDOF+1)+1)))
      end do
    enddo

    START_TIME= HECMW_WTIME()
    call hecmw_update_m_R (hecMESH, scale, hecMESH%n_node, NDOF)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, NP
      do j = 1, NDOF
        ip(j)=NDOF*(i-1)+j
      end do
      do j = 1, NDOF
        do k = 1, NDOF
          D(NDOF2*(i-1)+NDOF*(j-1)+k)=D(NDOF2*(i-1)+NDOF*(j-1)+k)*scale(ip(j))*scale(ip(k))
        end do
      end do

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1)+ii
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            AL(NDOF2*(k-1)+NDOF*(ii-1)+ij)=AL(NDOF2*(k-1)+NDOF*(ii-1)+ij)*scale(ip(ii))*scale(iq(ij))
          end do
        end do
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1)+ii
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            AU(NDOF2*(k-1)+NDOF*(ii-1)+ij)=AU(NDOF2*(k-1)+NDOF*(ii-1)+ij)*scale(ip(ii))*scale(iq(ij))
          end do
        end do
      enddo
    enddo
    !*voption indep (B,SCALE)
    do i= 1, N
      do k = 1, NDOF
        B(NDOF*(i-1)+k)=B(NDOF*(i-1)+k)*scale(NDOF*(i-1)+k)
      end do
    enddo
  end subroutine hecmw_solver_scaling_fw_nn

  subroutine hecmw_solver_scaling_bk_nn(hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: N, NP, NDOF, NDOF2
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:), X(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i,j,k,ii,ij, ip(hecMAT%NDOF),iq(hecMAT%NDOF)
    integer(kind=kint) :: iq1, iq2, iq3, iq4, iq5, iq6
    integer(kind=kint) :: isL, ieL, isU, ieU, inod

    if (hecmw_mat_get_scaling(hecMAT).eq.0) return

    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
    NDOF2 = NDOF*NDOF
    D => hecMAT%D
    AL => hecMAT%AL
    AU => hecMAT%AU
    INL => hecMAT%indexL
    IAL => hecMAT%itemL
    INU => hecMAT%indexU
    IAU => hecMAT%itemU
    B => hecMAT%B
    X => hecMAT%X

    !*voption indep (X,B,SCALE)
    do i= 1, N
      do k=1,NDOF
        X(NDOF*(i-1)+k)=X(NDOF*(i-1)+k)*scale(NDOF*(i-1)+k)
        B(NDOF*(i-1)+k)=B(NDOF*(i-1)+k)/scale(NDOF*(i-1)+k)
      end do
    enddo

    do i= 1, NP
      do j = 1, NDOF
        ip(j)=NDOF*(i-1)+j
      end do
      do j = 1, NDOF
        do k = 1, NDOF
          D(NDOF2*(i-1)+NDOF*(j-1)+k)=D(NDOF2*(i-1)+NDOF*(j-1)+k)/(scale(ip(j))*scale(ip(k)))
        end do
      end do

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1)+ii
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            AL(NDOF2*(k-1)+NDOF*(ii-1)+ij)=AL(NDOF2*(k-1)+NDOF*(ii-1)+ij)/(scale(ip(ii))*scale(iq(ij)))
          end do
        end do
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        do ii = 1, NDOF
          iq(ii) = NDOF*(inod-1)+ii
        end do
        do ii = 1, NDOF
          do ij = 1, NDOF
            AU(NDOF2*(k-1)+NDOF*(ii-1)+ij)=AU(NDOF2*(k-1)+NDOF*(ii-1)+ij)/(scale(ip(ii))*scale(iq(ij)))
          end do
        end do
      enddo
    enddo

    deallocate(scale)
  end subroutine hecmw_solver_scaling_bk_nn

end module hecmw_solver_scaling_nn
