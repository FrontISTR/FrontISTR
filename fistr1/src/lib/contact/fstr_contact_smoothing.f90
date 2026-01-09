!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact surface smoothing using Nagata patch interpolation
module m_fstr_contact_smoothing
  use hecmw, only: kint, kreal
  use elementInfo
  use m_utilities, only: calInverse
  implicit none

  ! Tikhonov regularization parameter
  real(kind=kreal), parameter, private :: NAGATA_EPSILON = 1.0d-8

  private
  public :: compute_Cab
  public :: compute_interpolation_matrix_P

contains

  !> Compute Cab correction matrices for edge (a,b)
  !! Kab = (Aab^T*Aab + eps*I)^-1 * Aab^T * D * Aab
  !! Cab_a = 0.5*I - Kab, Cab_b = 0.5*I + Kab
  subroutine compute_Cab(na_vec, nb_vec, Cab_a, Cab_b)
    implicit none
    real(kind=kreal), intent(in)  :: na_vec(3)   ! unit normal at vertex a
    real(kind=kreal), intent(in)  :: nb_vec(3)   ! unit normal at vertex b
    real(kind=kreal), intent(out) :: Cab_a(3,3)  ! correction matrix for vertex a
    real(kind=kreal), intent(out) :: Cab_b(3,3)  ! correction matrix for vertex b

    real(kind=kreal) :: Aab(2,3), D(2,2), AtA(3,3), AtA_inv(3,3)
    real(kind=kreal) :: Kab(3,3), temp(3,2), temp2(3,3)
    integer(kind=kint) :: i

    ! Build Aab = [na^T; nb^T]
    Aab(1,:) = na_vec(:)
    Aab(2,:) = nb_vec(:)

    ! D = diag(-1/4, +1/4)
    D = 0.0d0
    D(1,1) = -0.25d0
    D(2,2) =  0.25d0

    ! Compute AtA = Aab^T * Aab + eps*I
    AtA = matmul(transpose(Aab), Aab)
    do i = 1, 3
      AtA(i,i) = AtA(i,i) + NAGATA_EPSILON
    enddo

    ! Invert AtA using calInverse
    AtA_inv = AtA
    call calInverse(3, AtA_inv)

    ! Compute Kab = AtA_inv * Aab^T * D * Aab
    temp = matmul(transpose(Aab), D)     ! (3,2) = (3,2) * (2,2)
    temp2 = matmul(temp, Aab)            ! (3,3) = (3,2) * (2,3)
    Kab = matmul(AtA_inv, temp2)         ! (3,3) = (3,3) * (3,3)

    ! Cab_a = 0.5*I - Kab
    Cab_a = -Kab
    do i = 1, 3
      Cab_a(i,i) = Cab_a(i,i) + 0.5d0
    enddo

    ! Cab_b = 0.5*I + Kab
    Cab_b = Kab
    do i = 1, 3
      Cab_b(i,i) = Cab_b(i,i) + 0.5d0
    enddo

  end subroutine compute_Cab

  !> Compute interpolation matrix P for Nagata patch
  !! P_matrix(3, nnode*3) represents x = P * [X1; X2; X3; ...]
  subroutine compute_interpolation_matrix_P(etype, nnode, lpos, vertex_normals, P_matrix)
    implicit none
    integer(kind=kint), intent(in)  :: etype          ! element type
    integer(kind=kint), intent(in)  :: nnode          ! number of nodes
    real(kind=kreal), intent(in)    :: lpos(2)        ! local coordinates (xi, eta)
    real(kind=kreal), intent(in)    :: vertex_normals(3,nnode)  ! unit normals at vertices
    real(kind=kreal), intent(out)   :: P_matrix(3,nnode*3)      ! interpolation matrix

    real(kind=kreal) :: N(6)
    real(kind=kreal) :: Cab_12_1(3,3), Cab_12_2(3,3)
    real(kind=kreal) :: Cab_23_2(3,3), Cab_23_3(3,3)
    real(kind=kreal) :: Cab_31_3(3,3), Cab_31_1(3,3)
    real(kind=kreal) :: P1(3,3), P2(3,3), P3(3,3)
    real(kind=kreal) :: I3(3,3)
    integer(kind=kint) :: i, j

    P_matrix = 0.0d0

    ! Identity matrix
    I3 = 0.0d0
    I3(1,1) = 1.0d0
    I3(2,2) = 1.0d0
    I3(3,3) = 1.0d0

    select case(etype)
    case(fe_tri3n)
      ! Get quadratic shape functions for 6-node triangle
      ! N(1:3) are vertex functions, N(4:6) are edge functions
      call getShapeFunc(fe_tri6n, lpos, N)
      
      ! For Nagata patch, use hierarchical decomposition:
      ! Vertex contributions use linear parts, edge contributions use quadratic parts

      ! Compute Cab for 3 edges: (1,2), (2,3), (3,1)
      call compute_Cab(vertex_normals(:,1), vertex_normals(:,2), Cab_12_1, Cab_12_2)
      call compute_Cab(vertex_normals(:,2), vertex_normals(:,3), Cab_23_2, Cab_23_3)
      call compute_Cab(vertex_normals(:,3), vertex_normals(:,1), Cab_31_3, Cab_31_1)

      ! P1 = N1*I + N4*Cab_12_1 + N6*Cab_31_1
      P1 = N(1) * I3 + N(4) * Cab_12_1 + N(6) * Cab_31_1

      ! P2 = N2*I + N4*Cab_12_2 + N5*Cab_23_2
      P2 = N(2) * I3 + N(4) * Cab_12_2 + N(5) * Cab_23_2

      ! P3 = N3*I + N5*Cab_23_3 + N6*Cab_31_3
      P3 = N(3) * I3 + N(5) * Cab_23_3 + N(6) * Cab_31_3

      ! Assemble P = [P1, P2, P3]
      P_matrix(1:3, 1:3) = P1
      P_matrix(1:3, 4:6) = P2
      P_matrix(1:3, 7:9) = P3

    case default
      write(*,*) "Error: compute_interpolation_matrix_P - Unsupported element type for Nagata patch.",etype
      stop 
    end select

  end subroutine compute_interpolation_matrix_P

end module m_fstr_contact_smoothing
