!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                       Naoki Morita (Univ. of Tokyo)                  !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
module hecmw_solver_scaling_66
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  implicit none

  private
  real(kind=kreal), private, allocatable :: SCALE(:)

  public :: hecmw_solver_scaling_fw_66
  public :: hecmw_solver_scaling_bk_66

contains

  subroutine hecmw_solver_scaling_fw_66(hecMESH, hecMAT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, ip4, ip5, ip6
    integer(kind=kint) :: iq1, iq2, iq3, iq4, iq5, iq6
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
      SCALE (6*i-5)= 1.d0/dsqrt(dabs(D(36*i-35)))
      SCALE (6*i-4)= 1.d0/dsqrt(dabs(D(36*i-28)))
      SCALE (6*i-3)= 1.d0/dsqrt(dabs(D(36*i-21)))
      SCALE (6*i-2)= 1.d0/dsqrt(dabs(D(36*i-14)))
      SCALE (6*i-1)= 1.d0/dsqrt(dabs(D(36*i-7)))
      SCALE (6*i  )= 1.d0/dsqrt(dabs(D(36*i  )))
    enddo

    START_TIME= HECMW_WTIME()
    call hecmw_update_3_R (hecMESH, SCALE, hecMESH%n_node)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, NP
      ip1= 6*i-5
      ip2= 6*i-4
      ip3= 6*i-3
      ip4= 6*i-2
      ip5= 6*i-1
      ip6= 6*i

      D(36*i-35)= D(36*i-35)*SCALE(ip1)*SCALE(ip1)
      D(36*i-34)= D(36*i-34)*SCALE(ip1)*SCALE(ip2)
      D(36*i-33)= D(36*i-33)*SCALE(ip1)*SCALE(ip3)
      D(36*i-32)= D(36*i-32)*SCALE(ip1)*SCALE(ip4)
      D(36*i-31)= D(36*i-31)*SCALE(ip1)*SCALE(ip5)
      D(36*i-30)= D(36*i-30)*SCALE(ip1)*SCALE(ip6)

      D(36*i-29)= D(36*i-29)*SCALE(ip2)*SCALE(ip1)
      D(36*i-28)= D(36*i-28)*SCALE(ip2)*SCALE(ip2)
      D(36*i-27)= D(36*i-27)*SCALE(ip2)*SCALE(ip3)
      D(36*i-26)= D(36*i-26)*SCALE(ip2)*SCALE(ip4)
      D(36*i-25)= D(36*i-25)*SCALE(ip2)*SCALE(ip5)
      D(36*i-24)= D(36*i-24)*SCALE(ip2)*SCALE(ip6)

      D(36*i-23)= D(36*i-23)*SCALE(ip3)*SCALE(ip1)
      D(36*i-22)= D(36*i-22)*SCALE(ip3)*SCALE(ip2)
      D(36*i-21)= D(36*i-21)*SCALE(ip3)*SCALE(ip3)
      D(36*i-20)= D(36*i-20)*SCALE(ip3)*SCALE(ip4)
      D(36*i-19)= D(36*i-19)*SCALE(ip3)*SCALE(ip5)
      D(36*i-18)= D(36*i-18)*SCALE(ip3)*SCALE(ip6)

      D(36*i-17)= D(36*i-17)*SCALE(ip4)*SCALE(ip1)
      D(36*i-16)= D(36*i-16)*SCALE(ip4)*SCALE(ip2)
      D(36*i-15)= D(36*i-15)*SCALE(ip4)*SCALE(ip3)
      D(36*i-14)= D(36*i-14)*SCALE(ip4)*SCALE(ip4)
      D(36*i-13)= D(36*i-13)*SCALE(ip4)*SCALE(ip5)
      D(36*i-12)= D(36*i-12)*SCALE(ip4)*SCALE(ip6)

      D(36*i-11)= D(36*i-11)*SCALE(ip5)*SCALE(ip1)
      D(36*i-10)= D(36*i-10)*SCALE(ip5)*SCALE(ip2)
      D(36*i-9 )= D(36*i-9 )*SCALE(ip5)*SCALE(ip3)
      D(36*i-8 )= D(36*i-8 )*SCALE(ip5)*SCALE(ip4)
      D(36*i-7 )= D(36*i-7 )*SCALE(ip5)*SCALE(ip5)
      D(36*i-6 )= D(36*i-6 )*SCALE(ip5)*SCALE(ip6)

      D(36*i-5 )= D(36*i-5 )*SCALE(ip6)*SCALE(ip1)
      D(36*i-4 )= D(36*i-4 )*SCALE(ip6)*SCALE(ip2)
      D(36*i-3 )= D(36*i-3 )*SCALE(ip6)*SCALE(ip3)
      D(36*i-2 )= D(36*i-2 )*SCALE(ip6)*SCALE(ip4)
      D(36*i-1 )= D(36*i-1 )*SCALE(ip6)*SCALE(ip5)
      D(36*i   )= D(36*i   )*SCALE(ip6)*SCALE(ip6)

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 6*inod - 5
        iq2= 6*inod - 4
        iq3= 6*inod - 3
        iq4= 6*inod - 2
        iq5= 6*inod - 1
        iq6= 6*inod

        AL(36*k-35)= AL(36*k-35)*SCALE(ip1)*SCALE(iq1)
        AL(36*k-34)= AL(36*k-34)*SCALE(ip1)*SCALE(iq2)
        AL(36*k-33)= AL(36*k-33)*SCALE(ip1)*SCALE(iq3)
        AL(36*k-32)= AL(36*k-32)*SCALE(ip1)*SCALE(iq4)
        AL(36*k-31)= AL(36*k-31)*SCALE(ip1)*SCALE(iq5)
        AL(36*k-30)= AL(36*k-30)*SCALE(ip1)*SCALE(iq6)

        AL(36*k-29)= AL(36*k-29)*SCALE(ip2)*SCALE(iq1)
        AL(36*k-28)= AL(36*k-28)*SCALE(ip2)*SCALE(iq2)
        AL(36*k-27)= AL(36*k-27)*SCALE(ip2)*SCALE(iq3)
        AL(36*k-26)= AL(36*k-26)*SCALE(ip2)*SCALE(iq4)
        AL(36*k-25)= AL(36*k-25)*SCALE(ip2)*SCALE(iq5)
        AL(36*k-24)= AL(36*k-24)*SCALE(ip2)*SCALE(iq6)

        AL(36*k-23)= AL(36*k-23)*SCALE(ip3)*SCALE(iq1)
        AL(36*k-22)= AL(36*k-22)*SCALE(ip3)*SCALE(iq2)
        AL(36*k-21)= AL(36*k-21)*SCALE(ip3)*SCALE(iq3)
        AL(36*k-20)= AL(36*k-20)*SCALE(ip3)*SCALE(iq4)
        AL(36*k-19)= AL(36*k-19)*SCALE(ip3)*SCALE(iq5)
        AL(36*k-18)= AL(36*k-18)*SCALE(ip3)*SCALE(iq6)

        AL(36*k-17)= AL(36*k-17)*SCALE(ip4)*SCALE(iq1)
        AL(36*k-16)= AL(36*k-16)*SCALE(ip4)*SCALE(iq2)
        AL(36*k-15)= AL(36*k-15)*SCALE(ip4)*SCALE(iq3)
        AL(36*k-14)= AL(36*k-14)*SCALE(ip4)*SCALE(iq4)
        AL(36*k-13)= AL(36*k-13)*SCALE(ip4)*SCALE(iq5)
        AL(36*k-12)= AL(36*k-12)*SCALE(ip4)*SCALE(iq6)

        AL(36*k-11)= AL(36*k-11)*SCALE(ip5)*SCALE(iq1)
        AL(36*k-10)= AL(36*k-10)*SCALE(ip5)*SCALE(iq2)
        AL(36*k-9 )= AL(36*k-9 )*SCALE(ip5)*SCALE(iq3)
        AL(36*k-8 )= AL(36*k-8 )*SCALE(ip5)*SCALE(iq4)
        AL(36*k-7 )= AL(36*k-7 )*SCALE(ip5)*SCALE(iq5)
        AL(36*k-6 )= AL(36*k-6 )*SCALE(ip5)*SCALE(iq6)

        AL(36*k-5 )= AL(36*k-5 )*SCALE(ip6)*SCALE(iq1)
        AL(36*k-4 )= AL(36*k-4 )*SCALE(ip6)*SCALE(iq2)
        AL(36*k-3 )= AL(36*k-3 )*SCALE(ip6)*SCALE(iq3)
        AL(36*k-2 )= AL(36*k-2 )*SCALE(ip6)*SCALE(iq4)
        AL(36*k-1 )= AL(36*k-1 )*SCALE(ip6)*SCALE(iq5)
        AL(36*k   )= AL(36*k   )*SCALE(ip6)*SCALE(iq6)
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 6*inod - 5
        iq2= 6*inod - 4
        iq3= 6*inod - 3
        iq4= 6*inod - 2
        iq5= 6*inod - 1
        iq6= 6*inod

        AU(36*k-35)= AU(36*k-35)*SCALE(ip1)*SCALE(iq1)
        AU(36*k-34)= AU(36*k-34)*SCALE(ip1)*SCALE(iq2)
        AU(36*k-33)= AU(36*k-33)*SCALE(ip1)*SCALE(iq3)
        AU(36*k-32)= AU(36*k-32)*SCALE(ip1)*SCALE(iq4)
        AU(36*k-31)= AU(36*k-31)*SCALE(ip1)*SCALE(iq5)
        AU(36*k-30)= AU(36*k-30)*SCALE(ip1)*SCALE(iq6)

        AU(36*k-29)= AU(36*k-29)*SCALE(ip2)*SCALE(iq1)
        AU(36*k-28)= AU(36*k-28)*SCALE(ip2)*SCALE(iq2)
        AU(36*k-27)= AU(36*k-27)*SCALE(ip2)*SCALE(iq3)
        AU(36*k-26)= AU(36*k-26)*SCALE(ip2)*SCALE(iq4)
        AU(36*k-25)= AU(36*k-25)*SCALE(ip2)*SCALE(iq5)
        AU(36*k-24)= AU(36*k-24)*SCALE(ip2)*SCALE(iq6)

        AU(36*k-23)= AU(36*k-23)*SCALE(ip3)*SCALE(iq1)
        AU(36*k-22)= AU(36*k-22)*SCALE(ip3)*SCALE(iq2)
        AU(36*k-21)= AU(36*k-21)*SCALE(ip3)*SCALE(iq3)
        AU(36*k-20)= AU(36*k-20)*SCALE(ip3)*SCALE(iq4)
        AU(36*k-19)= AU(36*k-19)*SCALE(ip3)*SCALE(iq5)
        AU(36*k-18)= AU(36*k-18)*SCALE(ip3)*SCALE(iq6)

        AU(36*k-17)= AU(36*k-17)*SCALE(ip4)*SCALE(iq1)
        AU(36*k-16)= AU(36*k-16)*SCALE(ip4)*SCALE(iq2)
        AU(36*k-15)= AU(36*k-15)*SCALE(ip4)*SCALE(iq3)
        AU(36*k-14)= AU(36*k-14)*SCALE(ip4)*SCALE(iq4)
        AU(36*k-13)= AU(36*k-13)*SCALE(ip4)*SCALE(iq5)
        AU(36*k-12)= AU(36*k-12)*SCALE(ip4)*SCALE(iq6)

        AU(36*k-11)= AU(36*k-11)*SCALE(ip5)*SCALE(iq1)
        AU(36*k-10)= AU(36*k-10)*SCALE(ip5)*SCALE(iq2)
        AU(36*k-9 )= AU(36*k-9 )*SCALE(ip5)*SCALE(iq3)
        AU(36*k-8 )= AU(36*k-8 )*SCALE(ip5)*SCALE(iq4)
        AU(36*k-7 )= AU(36*k-7 )*SCALE(ip5)*SCALE(iq5)
        AU(36*k-6 )= AU(36*k-6 )*SCALE(ip5)*SCALE(iq6)

        AU(36*k-5 )= AU(36*k-5 )*SCALE(ip6)*SCALE(iq1)
        AU(36*k-4 )= AU(36*k-4 )*SCALE(ip6)*SCALE(iq2)
        AU(36*k-3 )= AU(36*k-3 )*SCALE(ip6)*SCALE(iq3)
        AU(36*k-2 )= AU(36*k-2 )*SCALE(ip6)*SCALE(iq4)
        AU(36*k-1 )= AU(36*k-1 )*SCALE(ip6)*SCALE(iq5)
        AU(36*k   )= AU(36*k   )*SCALE(ip6)*SCALE(iq6)
      enddo
    enddo
    !*voption indep (B,SCALE)
    do i= 1, N
      B(6*i-5)= B(6*i-5) * SCALE(6*i-5)
      B(6*i-4)= B(6*i-4) * SCALE(6*i-4)
      B(6*i-3)= B(6*i-3) * SCALE(6*i-3)
      B(6*i-2)= B(6*i-2) * SCALE(6*i-2)
      B(6*i-1)= B(6*i-1) * SCALE(6*i-1)
      B(6*i  )= B(6*i  ) * SCALE(6*i  )
    enddo
  end subroutine hecmw_solver_scaling_fw_66

  subroutine hecmw_solver_scaling_bk_66(hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: N, NP, NDOF
    real(kind=kreal), pointer :: D(:), AL(:), AU(:), B(:), X(:)
    integer(kind=kint), pointer :: INL(:), IAL(:), INU(:), IAU(:)
    integer(kind=kint) :: i, k, ip1, ip2, ip3, ip4, ip5, ip6
    integer(kind=kint) :: iq1, iq2, iq3, iq4, iq5, iq6
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
      X(6*i-5)= X(6*i-5) * SCALE(6*i-5)
      X(6*i-4)= X(6*i-4) * SCALE(6*i-4)
      X(6*i-3)= X(6*i-3) * SCALE(6*i-3)
      X(6*i-2)= X(6*i-2) * SCALE(6*i-2)
      X(6*i-1)= X(6*i-1) * SCALE(6*i-1)
      X(6*i  )= X(6*i  ) * SCALE(6*i  )

      B(6*i-5)= B(6*i-5) / SCALE(6*i-5)
      B(6*i-4)= B(6*i-4) / SCALE(6*i-4)
      B(6*i-3)= B(6*i-3) / SCALE(6*i-3)
      B(6*i-2)= B(6*i-2) / SCALE(6*i-2)
      B(6*i-1)= B(6*i-1) / SCALE(6*i-1)
      B(6*i  )= B(6*i  ) / SCALE(6*i  )
    enddo

    do i= 1, NP
      ip1= 6*i-5
      ip2= 6*i-4
      ip3= 6*i-3
      ip4= 6*i-2
      ip5= 6*i-1
      ip6= 6*i

      D(36*i-35)= D(36*i-35)/(SCALE(ip1)*SCALE(ip1))
      D(36*i-34)= D(36*i-34)/(SCALE(ip1)*SCALE(ip2))
      D(36*i-33)= D(36*i-33)/(SCALE(ip1)*SCALE(ip3))
      D(36*i-32)= D(36*i-32)/(SCALE(ip1)*SCALE(ip4))
      D(36*i-31)= D(36*i-31)/(SCALE(ip1)*SCALE(ip5))
      D(36*i-30)= D(36*i-30)/(SCALE(ip1)*SCALE(ip6))

      D(36*i-29)= D(36*i-29)/(SCALE(ip2)*SCALE(ip1))
      D(36*i-28)= D(36*i-28)/(SCALE(ip2)*SCALE(ip2))
      D(36*i-27)= D(36*i-27)/(SCALE(ip2)*SCALE(ip3))
      D(36*i-26)= D(36*i-26)/(SCALE(ip2)*SCALE(ip4))
      D(36*i-25)= D(36*i-25)/(SCALE(ip2)*SCALE(ip5))
      D(36*i-24)= D(36*i-24)/(SCALE(ip2)*SCALE(ip6))

      D(36*i-23)= D(36*i-23)/(SCALE(ip3)*SCALE(ip1))
      D(36*i-22)= D(36*i-22)/(SCALE(ip3)*SCALE(ip2))
      D(36*i-21)= D(36*i-21)/(SCALE(ip3)*SCALE(ip3))
      D(36*i-20)= D(36*i-20)/(SCALE(ip3)*SCALE(ip4))
      D(36*i-19)= D(36*i-19)/(SCALE(ip3)*SCALE(ip5))
      D(36*i-18)= D(36*i-18)/(SCALE(ip3)*SCALE(ip6))

      D(36*i-17)= D(36*i-17)/(SCALE(ip4)*SCALE(ip1))
      D(36*i-16)= D(36*i-16)/(SCALE(ip4)*SCALE(ip2))
      D(36*i-15)= D(36*i-15)/(SCALE(ip4)*SCALE(ip3))
      D(36*i-14)= D(36*i-14)/(SCALE(ip4)*SCALE(ip4))
      D(36*i-13)= D(36*i-13)/(SCALE(ip4)*SCALE(ip5))
      D(36*i-12)= D(36*i-12)/(SCALE(ip4)*SCALE(ip6))

      D(36*i-11)= D(36*i-11)/(SCALE(ip5)*SCALE(ip1))
      D(36*i-10)= D(36*i-10)/(SCALE(ip5)*SCALE(ip2))
      D(36*i-9 )= D(36*i-9 )/(SCALE(ip5)*SCALE(ip3))
      D(36*i-8 )= D(36*i-8 )/(SCALE(ip5)*SCALE(ip4))
      D(36*i-7 )= D(36*i-7 )/(SCALE(ip5)*SCALE(ip5))
      D(36*i-6 )= D(36*i-6 )/(SCALE(ip5)*SCALE(ip6))

      D(36*i-5 )= D(36*i-5 )/(SCALE(ip6)*SCALE(ip1))
      D(36*i-4 )= D(36*i-4 )/(SCALE(ip6)*SCALE(ip2))
      D(36*i-3 )= D(36*i-3 )/(SCALE(ip6)*SCALE(ip3))
      D(36*i-2 )= D(36*i-2 )/(SCALE(ip6)*SCALE(ip4))
      D(36*i-1 )= D(36*i-1 )/(SCALE(ip6)*SCALE(ip5))
      D(36*i   )= D(36*i   )/(SCALE(ip6)*SCALE(ip6))

      isL= INL(i-1) + 1
      ieL= INL(i  )
      !*voption indep (IAL,AL,SCALE)
      do k= isL, ieL
        inod= IAL(k)
        iq1= 6*inod - 5
        iq2= 6*inod - 5
        iq3= 6*inod - 3
        iq4= 6*inod - 2
        iq5= 6*inod - 1
        iq6= 6*inod

        AL(36*k-35)= AL(36*k-35)/(SCALE(ip1)*SCALE(iq1))
        AL(36*k-34)= AL(36*k-34)/(SCALE(ip1)*SCALE(iq2))
        AL(36*k-33)= AL(36*k-33)/(SCALE(ip1)*SCALE(iq3))
        AL(36*k-32)= AL(36*k-32)/(SCALE(ip1)*SCALE(iq4))
        AL(36*k-31)= AL(36*k-31)/(SCALE(ip1)*SCALE(iq5))
        AL(36*k-30)= AL(36*k-30)/(SCALE(ip1)*SCALE(iq6))

        AL(36*k-29)= AL(36*k-29)/(SCALE(ip2)*SCALE(iq1))
        AL(36*k-28)= AL(36*k-28)/(SCALE(ip2)*SCALE(iq2))
        AL(36*k-27)= AL(36*k-27)/(SCALE(ip2)*SCALE(iq3))
        AL(36*k-26)= AL(36*k-26)/(SCALE(ip2)*SCALE(iq4))
        AL(36*k-25)= AL(36*k-25)/(SCALE(ip2)*SCALE(iq5))
        AL(36*k-24)= AL(36*k-24)/(SCALE(ip2)*SCALE(iq6))

        AL(36*k-23)= AL(36*k-23)/(SCALE(ip3)*SCALE(iq1))
        AL(36*k-22)= AL(36*k-22)/(SCALE(ip3)*SCALE(iq2))
        AL(36*k-21)= AL(36*k-21)/(SCALE(ip3)*SCALE(iq3))
        AL(36*k-20)= AL(36*k-20)/(SCALE(ip3)*SCALE(iq4))
        AL(36*k-19)= AL(36*k-19)/(SCALE(ip3)*SCALE(iq5))
        AL(36*k-18)= AL(36*k-18)/(SCALE(ip3)*SCALE(iq6))

        AL(36*k-17)= AL(36*k-17)/(SCALE(ip4)*SCALE(iq1))
        AL(36*k-16)= AL(36*k-16)/(SCALE(ip4)*SCALE(iq2))
        AL(36*k-15)= AL(36*k-15)/(SCALE(ip4)*SCALE(iq3))
        AL(36*k-14)= AL(36*k-14)/(SCALE(ip4)*SCALE(iq4))
        AL(36*k-13)= AL(36*k-13)/(SCALE(ip4)*SCALE(iq5))
        AL(36*k-12)= AL(36*k-12)/(SCALE(ip4)*SCALE(iq6))

        AL(36*k-11)= AL(36*k-11)/(SCALE(ip5)*SCALE(iq1))
        AL(36*k-10)= AL(36*k-10)/(SCALE(ip5)*SCALE(iq2))
        AL(36*k-9 )= AL(36*k-9 )/(SCALE(ip5)*SCALE(iq3))
        AL(36*k-8 )= AL(36*k-8 )/(SCALE(ip5)*SCALE(iq4))
        AL(36*k-7 )= AL(36*k-7 )/(SCALE(ip5)*SCALE(iq5))
        AL(36*k-6 )= AL(36*k-6 )/(SCALE(ip5)*SCALE(iq6))

        AL(36*k-5 )= AL(36*k-5 )/(SCALE(ip6)*SCALE(iq1))
        AL(36*k-4 )= AL(36*k-4 )/(SCALE(ip6)*SCALE(iq2))
        AL(36*k-3 )= AL(36*k-3 )/(SCALE(ip6)*SCALE(iq3))
        AL(36*k-2 )= AL(36*k-2 )/(SCALE(ip6)*SCALE(iq4))
        AL(36*k-1 )= AL(36*k-1 )/(SCALE(ip6)*SCALE(iq5))
        AL(36*k   )= AL(36*k   )/(SCALE(ip6)*SCALE(iq6))
      enddo

      isU= INU(i-1) + 1
      ieU= INU(i  )
      !*voption indep (IAU,AU,SCALE)
      do k= isU, ieU
        inod= IAU(k)
        iq1= 6*inod - 5
        iq2= 6*inod - 4
        iq3= 6*inod - 3
        iq4= 6*inod - 2
        iq5= 6*inod - 1
        iq6= 6*inod

        AU(36*k-35)= AU(36*k-35)/(SCALE(ip1)*SCALE(iq1))
        AU(36*k-34)= AU(36*k-34)/(SCALE(ip1)*SCALE(iq2))
        AU(36*k-33)= AU(36*k-33)/(SCALE(ip1)*SCALE(iq3))
        AU(36*k-32)= AU(36*k-32)/(SCALE(ip1)*SCALE(iq4))
        AU(36*k-31)= AU(36*k-31)/(SCALE(ip1)*SCALE(iq5))
        AU(36*k-30)= AU(36*k-30)/(SCALE(ip1)*SCALE(iq6))

        AU(36*k-29)= AU(36*k-29)/(SCALE(ip2)*SCALE(iq1))
        AU(36*k-28)= AU(36*k-28)/(SCALE(ip2)*SCALE(iq2))
        AU(36*k-27)= AU(36*k-27)/(SCALE(ip2)*SCALE(iq3))
        AU(36*k-26)= AU(36*k-26)/(SCALE(ip2)*SCALE(iq4))
        AU(36*k-25)= AU(36*k-25)/(SCALE(ip2)*SCALE(iq5))
        AU(36*k-24)= AU(36*k-24)/(SCALE(ip2)*SCALE(iq6))

        AU(36*k-23)= AU(36*k-23)/(SCALE(ip3)*SCALE(iq1))
        AU(36*k-22)= AU(36*k-22)/(SCALE(ip3)*SCALE(iq2))
        AU(36*k-21)= AU(36*k-21)/(SCALE(ip3)*SCALE(iq3))
        AU(36*k-20)= AU(36*k-20)/(SCALE(ip3)*SCALE(iq4))
        AU(36*k-19)= AU(36*k-19)/(SCALE(ip3)*SCALE(iq5))
        AU(36*k-18)= AU(36*k-18)/(SCALE(ip3)*SCALE(iq6))

        AU(36*k-17)= AU(36*k-17)/(SCALE(ip4)*SCALE(iq1))
        AU(36*k-16)= AU(36*k-16)/(SCALE(ip4)*SCALE(iq2))
        AU(36*k-15)= AU(36*k-15)/(SCALE(ip4)*SCALE(iq3))
        AU(36*k-14)= AU(36*k-14)/(SCALE(ip4)*SCALE(iq4))
        AU(36*k-13)= AU(36*k-13)/(SCALE(ip4)*SCALE(iq5))
        AU(36*k-12)= AU(36*k-12)/(SCALE(ip4)*SCALE(iq6))

        AU(36*k-11)= AU(36*k-11)/(SCALE(ip5)*SCALE(iq1))
        AU(36*k-10)= AU(36*k-10)/(SCALE(ip5)*SCALE(iq2))
        AU(36*k-9 )= AU(36*k-9 )/(SCALE(ip5)*SCALE(iq3))
        AU(36*k-8 )= AU(36*k-8 )/(SCALE(ip5)*SCALE(iq4))
        AU(36*k-7 )= AU(36*k-7 )/(SCALE(ip5)*SCALE(iq5))
        AU(36*k-6 )= AU(36*k-6 )/(SCALE(ip5)*SCALE(iq6))

        AU(36*k-5 )= AU(36*k-5 )/(SCALE(ip6)*SCALE(iq1))
        AU(36*k-4 )= AU(36*k-4 )/(SCALE(ip6)*SCALE(iq2))
        AU(36*k-3 )= AU(36*k-3 )/(SCALE(ip6)*SCALE(iq3))
        AU(36*k-2 )= AU(36*k-2 )/(SCALE(ip6)*SCALE(iq4))
        AU(36*k-1 )= AU(36*k-1 )/(SCALE(ip6)*SCALE(iq5))
        AU(36*k   )= AU(36*k   )/(SCALE(ip6)*SCALE(iq6))
      enddo
    enddo

    deallocate(SCALE)
  end subroutine hecmw_solver_scaling_bk_66

end module hecmw_solver_scaling_66
