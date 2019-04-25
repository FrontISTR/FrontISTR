!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_scaling_33
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  implicit none

  private
  real(kind=kreal), private, allocatable :: scale(:)

  public :: hecmw_solver_scaling_fw_33
  public :: hecmw_solver_scaling_bk_33

contains

  subroutine hecmw_solver_scaling_fw_33(hecMESH, hecMAT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, iq1, iq2, iq3
    integer(kind=kint) :: isL, ieL, isU, ieU, inod
    real(kind=kreal) :: START_TIME, END_TIME

    if (hecmw_mat_get_scaling(hecMAT).eq.0) return

    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
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
      scale (3*i-2)= 1.d0/dsqrt(dabs(D(9*i-8)))
      scale (3*i-1)= 1.d0/dsqrt(dabs(D(9*i-4)))
      scale (3*i  )= 1.d0/dsqrt(dabs(D(9*i  )))
    enddo

    START_TIME= HECMW_WTIME()
    call hecmw_update_3_R (hecMESH, scale, hecMESH%n_node)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, NP
      ip1= 3*i-2
      ip2= 3*i-1
      ip3= 3*i
      D(9*i-8)= D(9*i-8)*scale(ip1)*scale(ip1)
      D(9*i-7)= D(9*i-7)*scale(ip1)*scale(ip2)
      D(9*i-6)= D(9*i-6)*scale(ip1)*scale(ip3)
      D(9*i-5)= D(9*i-5)*scale(ip2)*scale(ip1)
      D(9*i-4)= D(9*i-4)*scale(ip2)*scale(ip2)
      D(9*i-3)= D(9*i-3)*scale(ip2)*scale(ip3)
      D(9*i-2)= D(9*i-2)*scale(ip3)*scale(ip1)
      D(9*i-1)= D(9*i-1)*scale(ip3)*scale(ip2)
      D(9*i  )= D(9*i  )*scale(ip3)*scale(ip3)

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 3*inod - 2
        iq2= 3*inod - 1
        iq3= 3*inod
        AL(9*k-8)= AL(9*k-8)*scale(ip1)*scale(iq1)
        AL(9*k-7)= AL(9*k-7)*scale(ip1)*scale(iq2)
        AL(9*k-6)= AL(9*k-6)*scale(ip1)*scale(iq3)
        AL(9*k-5)= AL(9*k-5)*scale(ip2)*scale(iq1)
        AL(9*k-4)= AL(9*k-4)*scale(ip2)*scale(iq2)
        AL(9*k-3)= AL(9*k-3)*scale(ip2)*scale(iq3)
        AL(9*k-2)= AL(9*k-2)*scale(ip3)*scale(iq1)
        AL(9*k-1)= AL(9*k-1)*scale(ip3)*scale(iq2)
        AL(9*k  )= AL(9*k  )*scale(ip3)*scale(iq3)
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 3*inod - 2
        iq2= 3*inod - 1
        iq3= 3*inod
        AU(9*k-8)= AU(9*k-8)*scale(ip1)*scale(iq1)
        AU(9*k-7)= AU(9*k-7)*scale(ip1)*scale(iq2)
        AU(9*k-6)= AU(9*k-6)*scale(ip1)*scale(iq3)
        AU(9*k-5)= AU(9*k-5)*scale(ip2)*scale(iq1)
        AU(9*k-4)= AU(9*k-4)*scale(ip2)*scale(iq2)
        AU(9*k-3)= AU(9*k-3)*scale(ip2)*scale(iq3)
        AU(9*k-2)= AU(9*k-2)*scale(ip3)*scale(iq1)
        AU(9*k-1)= AU(9*k-1)*scale(ip3)*scale(iq2)
        AU(9*k  )= AU(9*k  )*scale(ip3)*scale(iq3)
      enddo
    enddo
    !*voption indep (B,SCALE)
    do i= 1, N
      B(3*i-2)= B(3*i-2) * scale(3*i-2)
      B(3*i-1)= B(3*i-1) * scale(3*i-1)
      B(3*i  )= B(3*i  ) * scale(3*i  )
    enddo
  end subroutine hecmw_solver_scaling_fw_33

  subroutine hecmw_solver_scaling_bk_33(hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:), X(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, iq1, iq2, iq3
    integer(kind=kint) :: isL, ieL, isU, ieU, inod

    if (hecmw_mat_get_scaling(hecMAT).eq.0) return

    N = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
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
      X(3*i-2)= X(3*i-2) * scale(3*i-2)
      X(3*i-1)= X(3*i-1) * scale(3*i-1)
      X(3*i  )= X(3*i  ) * scale(3*i  )
      B(3*i-2)= B(3*i-2) / scale(3*i-2)
      B(3*i-1)= B(3*i-1) / scale(3*i-1)
      B(3*i  )= B(3*i  ) / scale(3*i  )
    enddo

    do i= 1, NP
      ip1= 3*i-2
      ip2= 3*i-1
      ip3= 3*i
      D(9*i-8)= D(9*i-8)/(scale(ip1)*scale(ip1))
      D(9*i-7)= D(9*i-7)/(scale(ip1)*scale(ip2))
      D(9*i-6)= D(9*i-6)/(scale(ip1)*scale(ip3))
      D(9*i-5)= D(9*i-5)/(scale(ip2)*scale(ip1))
      D(9*i-4)= D(9*i-4)/(scale(ip2)*scale(ip2))
      D(9*i-3)= D(9*i-3)/(scale(ip2)*scale(ip3))
      D(9*i-2)= D(9*i-2)/(scale(ip3)*scale(ip1))
      D(9*i-1)= D(9*i-1)/(scale(ip3)*scale(ip2))
      D(9*i  )= D(9*i  )/(scale(ip3)*scale(ip3))

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 3*inod - 2
        iq2= 3*inod - 1
        iq3= 3*inod
        AL(9*k-8)= AL(9*k-8)/(scale(ip1)*scale(iq1))
        AL(9*k-7)= AL(9*k-7)/(scale(ip1)*scale(iq2))
        AL(9*k-6)= AL(9*k-6)/(scale(ip1)*scale(iq3))
        AL(9*k-5)= AL(9*k-5)/(scale(ip2)*scale(iq1))
        AL(9*k-4)= AL(9*k-4)/(scale(ip2)*scale(iq2))
        AL(9*k-3)= AL(9*k-3)/(scale(ip2)*scale(iq3))
        AL(9*k-2)= AL(9*k-2)/(scale(ip3)*scale(iq1))
        AL(9*k-1)= AL(9*k-1)/(scale(ip3)*scale(iq2))
        AL(9*k  )= AL(9*k  )/(scale(ip3)*scale(iq3))
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 3*inod - 2
        iq2= 3*inod - 1
        iq3= 3*inod
        AU(9*k-8)= AU(9*k-8)/(scale(ip1)*scale(iq1))
        AU(9*k-7)= AU(9*k-7)/(scale(ip1)*scale(iq2))
        AU(9*k-6)= AU(9*k-6)/(scale(ip1)*scale(iq3))
        AU(9*k-5)= AU(9*k-5)/(scale(ip2)*scale(iq1))
        AU(9*k-4)= AU(9*k-4)/(scale(ip2)*scale(iq2))
        AU(9*k-3)= AU(9*k-3)/(scale(ip2)*scale(iq3))
        AU(9*k-2)= AU(9*k-2)/(scale(ip3)*scale(iq1))
        AU(9*k-1)= AU(9*k-1)/(scale(ip3)*scale(iq2))
        AU(9*k  )= AU(9*k  )/(scale(ip3)*scale(iq3))
      enddo
    enddo

    deallocate(scale)
  end subroutine hecmw_solver_scaling_bk_33

end module hecmw_solver_scaling_33
