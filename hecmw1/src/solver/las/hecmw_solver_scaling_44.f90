!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_scaling_44
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  implicit none

  private
  real(kind=kreal), private, allocatable :: SCALE(:)

  public :: hecmw_solver_scaling_fw_44
  public :: hecmw_solver_scaling_bk_44

contains

  subroutine hecmw_solver_scaling_fw_44(hecMESH, hecMAT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, ip4, iq1, iq2, iq3, iq4
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

    allocate(SCALE(NDOF*NP))

    do i= 1, N
      SCALE (4*i-3)= 1.d0/dsqrt(dabs(D(16*i-15)))
      SCALE (4*i-2)= 1.d0/dsqrt(dabs(D(16*i-10)))
      SCALE (4*i-1)= 1.d0/dsqrt(dabs(D(16*i- 5)))
      SCALE (4*i  )= 1.d0/dsqrt(dabs(D(16*i   )))
    enddo

    START_TIME= HECMW_WTIME()
    call hecmw_update_4_R (hecMESH, SCALE, hecMESH%n_node)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, NP
      ip1= 4*i-3
      ip2= 4*i-2
      ip3= 4*i-1
      ip4= 4*i
      D(16*i-15)= D(16*i-15)*SCALE(ip1)*SCALE(ip1)
      D(16*i-14)= D(16*i-14)*SCALE(ip1)*SCALE(ip2)
      D(16*i-13)= D(16*i-13)*SCALE(ip1)*SCALE(ip3)
      D(16*i-12)= D(16*i-12)*SCALE(ip1)*SCALE(ip4)
      D(16*i-11)= D(16*i-11)*SCALE(ip2)*SCALE(ip1)
      D(16*i-10)= D(16*i-10)*SCALE(ip2)*SCALE(ip2)
      D(16*i- 9)= D(16*i- 9)*SCALE(ip2)*SCALE(ip3)
      D(16*i- 8)= D(16*i- 8)*SCALE(ip2)*SCALE(ip4)
      D(16*i- 7)= D(16*i- 7)*SCALE(ip3)*SCALE(ip1)
      D(16*i- 6)= D(16*i- 6)*SCALE(ip3)*SCALE(ip2)
      D(16*i- 5)= D(16*i- 5)*SCALE(ip3)*SCALE(ip3)
      D(16*i- 4)= D(16*i- 4)*SCALE(ip3)*SCALE(ip4)
      D(16*i- 3)= D(16*i- 3)*SCALE(ip4)*SCALE(ip1)
      D(16*i- 2)= D(16*i- 2)*SCALE(ip4)*SCALE(ip2)
      D(16*i- 1)= D(16*i- 1)*SCALE(ip4)*SCALE(ip3)
      D(16*i- 0)= D(16*i- 0)*SCALE(ip4)*SCALE(ip4)

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 4*inod - 3
        iq2= 4*inod - 2
        iq3= 4*inod - 1
        iq4= 4*inod
        AL(16*k-15)= AL(16*k-15)*SCALE(ip1)*SCALE(iq1)
        AL(16*k-14)= AL(16*k-14)*SCALE(ip1)*SCALE(iq2)
        AL(16*k-13)= AL(16*k-13)*SCALE(ip1)*SCALE(iq3)
        AL(16*k-12)= AL(16*k-12)*SCALE(ip1)*SCALE(iq4)
        AL(16*k-11)= AL(16*k-11)*SCALE(ip2)*SCALE(iq1)
        AL(16*k-10)= AL(16*k-10)*SCALE(ip2)*SCALE(iq2)
        AL(16*k- 9)= AL(16*k- 9)*SCALE(ip2)*SCALE(iq3)
        AL(16*k- 8)= AL(16*k- 8)*SCALE(ip2)*SCALE(iq4)
        AL(16*k- 7)= AL(16*k- 7)*SCALE(ip3)*SCALE(iq1)
        AL(16*k- 6)= AL(16*k- 6)*SCALE(ip3)*SCALE(iq2)
        AL(16*k- 5)= AL(16*k- 5)*SCALE(ip3)*SCALE(iq3)
        AL(16*k- 4)= AL(16*k- 4)*SCALE(ip3)*SCALE(iq4)
        AL(16*k- 3)= AL(16*k- 3)*SCALE(ip4)*SCALE(iq1)
        AL(16*k- 2)= AL(16*k- 2)*SCALE(ip4)*SCALE(iq2)
        AL(16*k- 1)= AL(16*k- 1)*SCALE(ip4)*SCALE(iq3)
        AL(16*k- 0)= AL(16*k- 0)*SCALE(ip4)*SCALE(iq4)
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 4*inod - 3
        iq2= 4*inod - 2
        iq3= 4*inod - 1
        iq4= 4*inod
        AU(16*k-15)= AU(16*k-15)*SCALE(ip1)*SCALE(iq1)
        AU(16*k-14)= AU(16*k-14)*SCALE(ip1)*SCALE(iq2)
        AU(16*k-13)= AU(16*k-13)*SCALE(ip1)*SCALE(iq3)
        AU(16*k-12)= AU(16*k-12)*SCALE(ip1)*SCALE(iq4)
        AU(16*k-11)= AU(16*k-11)*SCALE(ip2)*SCALE(iq1)
        AU(16*k-10)= AU(16*k-10)*SCALE(ip2)*SCALE(iq2)
        AU(16*k- 9)= AU(16*k- 9)*SCALE(ip2)*SCALE(iq3)
        AU(16*k- 8)= AU(16*k- 8)*SCALE(ip2)*SCALE(iq4)
        AU(16*k- 7)= AU(16*k- 7)*SCALE(ip3)*SCALE(iq1)
        AU(16*k- 6)= AU(16*k- 6)*SCALE(ip3)*SCALE(iq2)
        AU(16*k- 5)= AU(16*k- 5)*SCALE(ip3)*SCALE(iq3)
        AU(16*k- 4)= AU(16*k- 4)*SCALE(ip3)*SCALE(iq4)
        AU(16*k- 3)= AU(16*k- 3)*SCALE(ip4)*SCALE(iq1)
        AU(16*k- 2)= AU(16*k- 2)*SCALE(ip4)*SCALE(iq2)
        AU(16*k- 1)= AU(16*k- 1)*SCALE(ip4)*SCALE(iq3)
        AU(16*k- 0)= AU(16*k- 0)*SCALE(ip4)*SCALE(iq4)
      enddo
    enddo
    !*voption indep (B,SCALE)
    do i= 1, N
      B(4*i-3)= B(4*i-3) * SCALE(4*i-3)
      B(4*i-2)= B(4*i-2) * SCALE(4*i-2)
      B(4*i-1)= B(4*i-1) * SCALE(4*i-1)
      B(4*i  )= B(4*i  ) * SCALE(4*i  )
    enddo
  end subroutine hecmw_solver_scaling_fw_44

  subroutine hecmw_solver_scaling_bk_44(hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:), X(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, ip4, iq1, iq2, iq3, iq4
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
      X(4*i-3)= X(4*i-3) * SCALE(4*i-3)
      X(4*i-2)= X(4*i-2) * SCALE(4*i-2)
      X(4*i-1)= X(4*i-1) * SCALE(4*i-1)
      X(4*i  )= X(4*i  ) * SCALE(4*i  )
      B(4*i-3)= B(4*i-3) / SCALE(4*i-3)
      B(4*i-2)= B(4*i-2) / SCALE(4*i-2)
      B(4*i-1)= B(4*i-1) / SCALE(4*i-1)
      B(4*i  )= B(4*i  ) / SCALE(4*i  )
    enddo

    do i= 1, NP
      ip1= 4*i-3
      ip2= 4*i-2
      ip3= 4*i-1
      ip4= 4*i
      D(16*i-15)= D(16*i-15)/(SCALE(ip1)*SCALE(ip1))
      D(16*i-14)= D(16*i-14)/(SCALE(ip1)*SCALE(ip2))
      D(16*i-13)= D(16*i-13)/(SCALE(ip1)*SCALE(ip3))
      D(16*i-12)= D(16*i-12)/(SCALE(ip1)*SCALE(ip4))
      D(16*i-11)= D(16*i-11)/(SCALE(ip2)*SCALE(ip1))
      D(16*i-10)= D(16*i-10)/(SCALE(ip2)*SCALE(ip2))
      D(16*i- 9)= D(16*i- 9)/(SCALE(ip2)*SCALE(ip3))
      D(16*i- 8)= D(16*i- 8)/(SCALE(ip2)*SCALE(ip4))
      D(16*i- 7)= D(16*i- 7)/(SCALE(ip3)*SCALE(ip1))
      D(16*i- 6)= D(16*i- 6)/(SCALE(ip3)*SCALE(ip2))
      D(16*i- 5)= D(16*i- 5)/(SCALE(ip3)*SCALE(ip3))
      D(16*i- 4)= D(16*i- 4)/(SCALE(ip3)*SCALE(ip4))
      D(16*i- 3)= D(16*i- 3)/(SCALE(ip4)*SCALE(ip1))
      D(16*i- 2)= D(16*i- 2)/(SCALE(ip4)*SCALE(ip2))
      D(16*i- 1)= D(16*i- 1)/(SCALE(ip4)*SCALE(ip3))
      D(16*i- 0)= D(16*i- 0)/(SCALE(ip4)*SCALE(ip4))
      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 4*inod - 3
        iq2= 4*inod - 2
        iq3= 4*inod - 1
        iq4= 4*inod
        AL(16*k-15)= AL(16*k-15)/(SCALE(ip1)*SCALE(iq1))
        AL(16*k-14)= AL(16*k-14)/(SCALE(ip1)*SCALE(iq2))
        AL(16*k-13)= AL(16*k-13)/(SCALE(ip1)*SCALE(iq3))
        AL(16*k-12)= AL(16*k-12)/(SCALE(ip1)*SCALE(iq4))
        AL(16*k-11)= AL(16*k-11)/(SCALE(ip2)*SCALE(iq1))
        AL(16*k-10)= AL(16*k-10)/(SCALE(ip2)*SCALE(iq2))
        AL(16*k- 9)= AL(16*k- 9)/(SCALE(ip2)*SCALE(iq3))
        AL(16*k- 8)= AL(16*k- 8)/(SCALE(ip2)*SCALE(iq4))
        AL(16*k- 7)= AL(16*k- 7)/(SCALE(ip3)*SCALE(iq1))
        AL(16*k- 6)= AL(16*k- 6)/(SCALE(ip3)*SCALE(iq2))
        AL(16*k- 5)= AL(16*k- 5)/(SCALE(ip3)*SCALE(iq3))
        AL(16*k- 4)= AL(16*k- 4)/(SCALE(ip3)*SCALE(iq4))
        AL(16*k- 3)= AL(16*k- 3)/(SCALE(ip4)*SCALE(iq1))
        AL(16*k- 2)= AL(16*k- 2)/(SCALE(ip4)*SCALE(iq2))
        AL(16*k- 1)= AL(16*k- 1)/(SCALE(ip4)*SCALE(iq3))
        AL(16*k- 0)= AL(16*k- 0)/(SCALE(ip4)*SCALE(iq4))
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 4*inod - 3
        iq2= 4*inod - 2
        iq3= 4*inod - 1
        iq4= 4*inod
        AU(16*k-15)= AU(16*k-15)/(SCALE(ip1)*SCALE(iq1))
        AU(16*k-14)= AU(16*k-14)/(SCALE(ip1)*SCALE(iq2))
        AU(16*k-13)= AU(16*k-13)/(SCALE(ip1)*SCALE(iq3))
        AU(16*k-12)= AU(16*k-12)/(SCALE(ip1)*SCALE(iq4))
        AU(16*k-11)= AU(16*k-11)/(SCALE(ip2)*SCALE(iq1))
        AU(16*k-10)= AU(16*k-10)/(SCALE(ip2)*SCALE(iq2))
        AU(16*k- 9)= AU(16*k- 9)/(SCALE(ip2)*SCALE(iq3))
        AU(16*k- 8)= AU(16*k- 8)/(SCALE(ip2)*SCALE(iq4))
        AU(16*k- 7)= AU(16*k- 7)/(SCALE(ip3)*SCALE(iq1))
        AU(16*k- 6)= AU(16*k- 6)/(SCALE(ip3)*SCALE(iq2))
        AU(16*k- 5)= AU(16*k- 5)/(SCALE(ip3)*SCALE(iq3))
        AU(16*k- 4)= AU(16*k- 4)/(SCALE(ip3)*SCALE(iq4))
        AU(16*k- 3)= AU(16*k- 3)/(SCALE(ip4)*SCALE(iq1))
        AU(16*k- 2)= AU(16*k- 2)/(SCALE(ip4)*SCALE(iq2))
        AU(16*k- 1)= AU(16*k- 1)/(SCALE(ip4)*SCALE(iq3))
        AU(16*k- 0)= AU(16*k- 0)/(SCALE(ip4)*SCALE(iq4))
      enddo
    enddo

    deallocate(SCALE)
  end subroutine hecmw_solver_scaling_bk_44

end module hecmw_solver_scaling_44
