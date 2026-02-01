!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact surface smoothing using Nagata patch interpolation
module m_fstr_contact_smoothing
  use hecmw, only: kint, kreal
  use hecmw_util, only: HECMW_NAME_LEN
  use elementInfo
  use m_utilities, only: calInverse
  use mSurfElement, only: tSurfElement
  implicit none

  ! Tikhonov regularization parameter
  real(kind=kreal), parameter, private :: NAGATA_EPSILON = 1.0d-4

  integer(kind=kint), parameter, private :: DEBUG = 1

  private
  public :: compute_Cab
  public :: compute_interpolation_matrix_P
  public :: update_surface_normal
  public :: create_intermediate_points
  public :: reorder_tri3n_to_tri6n

contains

  !> Reorder triangle vertices and midpoints from fe_tri3n to fe_tri6n order
  !! fe_tri3n: v1(xi), v2(et), v3(st), m1(xi-et), m2(et-st), m3(st-xi)
  !! fe_tri6n: v1(st), v2(xi), v3(et), m4(xi-st), m5(xi-et), m6(et-st)
  subroutine reorder_tri3n_to_tri6n(vertices_in, midpoints_in, nodes_out)
    real(kind=kreal), intent(in)  :: vertices_in(3,3)   ! fe_tri3n vertices
    real(kind=kreal), intent(in)  :: midpoints_in(3,3)  ! fe_tri3n midpoints
    real(kind=kreal), intent(out) :: nodes_out(3,6)     ! fe_tri6n nodes
    
    ! Reorder vertices: xi,et,st -> st,xi,et
    nodes_out(1:3, 1) = vertices_in(1:3, 3)  ! st
    nodes_out(1:3, 2) = vertices_in(1:3, 1)  ! xi
    nodes_out(1:3, 3) = vertices_in(1:3, 2)  ! et
    
    ! Reorder midpoints: xi-et,et-st,st-xi -> xi-st,xi-et,et-st
    nodes_out(1:3, 4) = midpoints_in(1:3, 3)  ! xi-st (was st-xi, same point)
    nodes_out(1:3, 5) = midpoints_in(1:3, 1)  ! xi-et
    nodes_out(1:3, 6) = midpoints_in(1:3, 2)  ! et-st
  end subroutine reorder_tri3n_to_tri6n

  !> Reorder vertex normals from fe_tri3n to fe_tri6n order
  subroutine reorder_normals_tri3n_to_tri6n(normals_in, normals_out)
    real(kind=kreal), intent(in)  :: normals_in(3,3)   ! fe_tri3n order
    real(kind=kreal), intent(out) :: normals_out(3,3)  ! fe_tri6n order
    
    ! xi,et,st -> st,xi,et
    normals_out(1:3, 1) = normals_in(1:3, 3)  ! st
    normals_out(1:3, 2) = normals_in(1:3, 1)  ! xi
    normals_out(1:3, 3) = normals_in(1:3, 2)  ! et
  end subroutine reorder_normals_tri3n_to_tri6n

  subroutine update_surface_normal( surf, currpos )
    use mSurfElement, only: calc_all_surf_vertex_normals
    type(tSurfElement), intent(inout) :: surf(:)    !< surface elements
    real(kind=kreal), intent(in) :: currpos(:)      !< current coordinate of all nodes

    integer(kind=kint) :: i, j, nn, gnode, max_node
    real(kind=kreal) :: normal(3)
    real(kind=kreal), allocatable :: vnormal(:,:)
    integer(kind=kint), allocatable :: vcount(:)
    
    if (size(surf) == 0) return
    
    ! Calculate geometric normals at vertices
    call calc_all_surf_vertex_normals(surf, currpos)
    
    ! Apply smoothing by averaging normals at shared vertices
    max_node = maxval([(maxval(surf(i)%nodes), i=1,size(surf))])
    allocate(vnormal(3, max_node), vcount(max_node))
    vnormal = 0.0d0
    vcount = 0
    
    do i = 1, size(surf)
      nn = size(surf(i)%nodes)
      do j = 1, nn
        gnode = surf(i)%nodes(j)
        vnormal(:, gnode) = vnormal(:, gnode) + surf(i)%vertex_normals(:, j)
        vcount(gnode) = vcount(gnode) + 1
      enddo
    enddo
    
    do i = 1, size(surf)
      nn = size(surf(i)%nodes)
      do j = 1, nn
        gnode = surf(i)%nodes(j)
        if (vcount(gnode) > 0) then
          normal = vnormal(:, gnode) / dble(vcount(gnode))
          normal = normal / dsqrt(dot_product(normal, normal))
          surf(i)%vertex_normals(:, j) = normal
        endif
      enddo
    enddo
    
    deallocate(vnormal, vcount)

  end subroutine update_surface_normal

  !> Compute Cab correction matrices for edge (a,b)
  !! Kab = (Aab^T*Aab + eps*I)^-1 * Aab^T * D * Aab
  !! Cab_a = 0.5*I - Kab, Cab_b = 0.5*I + Kab
  subroutine compute_Cab_old(na_vec, nb_vec, Cab_a, Cab_b)
    implicit none
    real(kind=kreal), intent(in)  :: na_vec(3)   ! unit normal at vertex a
    real(kind=kreal), intent(in)  :: nb_vec(3)   ! unit normal at vertex b
    real(kind=kreal), intent(out) :: Cab_a(3,3)  ! correction matrix for vertex a
    real(kind=kreal), intent(out) :: Cab_b(3,3)  ! correction matrix for vertex b

    real(kind=kreal) :: Aab(2,3), AtA(3,3), AtA_orig(3,3), AtA_inv(3,3)
    real(kind=kreal) :: P(3,3)
    real(kind=kreal) :: epsilon
    integer(kind=kint) :: i

    ! Build Aab = [na^T; nb^T]
    Aab(1,:) = na_vec(:)
    Aab(2,:) = nb_vec(:)

    ! Compute AtA = Aab^T * Aab
    AtA = matmul(transpose(Aab), Aab)
    AtA_orig = AtA

    ! Regularize: AtA_reg = AtA + eps*I
    epsilon = (AtA(1,1)+AtA(2,2)+AtA(3,3)) * NAGATA_EPSILON
    do i = 1, 3
      AtA(i,i) = AtA(i,i) + epsilon
    enddo

    ! Invert AtA_reg using calInverse
    AtA_inv = AtA
    call calInverse(3, AtA_inv)

    ! Compute P = AtA_inv * AtA_orig
    P = matmul(AtA_inv, AtA_orig)

    ! Cab_a = 0.5*I + 0.25*P
    Cab_a = 0.25d0 * P
    do i = 1, 3
      Cab_a(i,i) = Cab_a(i,i) + 0.5d0
    enddo

    ! Cab_b = 0.5*I - 0.25*P
    Cab_b = -0.25d0 * P
    do i = 1, 3
      Cab_b(i,i) = Cab_b(i,i) + 0.5d0
    enddo

  end subroutine compute_Cab_old

  subroutine compute_Cab(na_vec, nb_vec, Cab_a, Cab_b)
    implicit none
    real(kind=kreal), intent(in)  :: na_vec(3)   ! unit normal at vertex a (n1)
    real(kind=kreal), intent(in)  :: nb_vec(3)   ! unit normal at vertex b (n2)
    real(kind=kreal), intent(out) :: Cab_a(3,3)
    real(kind=kreal), intent(out) :: Cab_b(3,3)

    real(kind=kreal) :: I3(3,3)
    real(kind=kreal) :: m(3), mnorm
    real(kind=kreal) :: a1, a2
    real(kind=kreal) :: beta, den, K
    real(kind=kreal) :: r(3)      ! r = a2*nb - a1*na
    real(kind=kreal) :: Q(3,3)    ! outer(m, r)
    integer(kind=kint) :: i, j

    ! Identity
    I3 = 0.0d0
    do i=1,3
      I3(i,i) = 1.0d0
    enddo

    ! m = normalize(na + nb)
    m(:) = na_vec(:) + nb_vec(:)
    mnorm = sqrt( m(1)*m(1) + m(2)*m(2) + m(3)*m(3) )

    ! Fallback when na ~ -nb (average nearly zero)
    if (mnorm < 1.0d-12) then
      m(:) = na_vec(:)
      mnorm = 1.0d0
    else
      m(:) = m(:) / mnorm
    endif

    ! a1 = m·na, a2 = m·nb
    a1 = m(1)*na_vec(1) + m(2)*na_vec(2) + m(3)*na_vec(3)
    a2 = m(1)*nb_vec(1) + m(2)*nb_vec(2) + m(3)*nb_vec(3)

    ! Regularization for alpha (beta >= 0)
    ! Option 1: fixed dimensionless parameter
    beta = NAGATA_EPSILON   ! << define e.g. 1e-3 .. 1e-1 (see notes below)

    den = 16.0d0*(a1*a1 + a2*a2) + beta

    ! If den too small, do no curvature (midpoint)
    if (den < 1.0d-20) then
      Cab_a = 0.5d0 * I3
      Cab_b = 0.5d0 * I3
      return
    endif

    K = 4.0d0 / den

    ! r = a2*nb - a1*na  (3-vector)
    r(:) = a2*nb_vec(:) - a1*na_vec(:)

    ! Q = outer(m, r) = m * r^T  (3x3)
    do i=1,3
      do j=1,3
        Q(i,j) = m(i) * r(j)
      enddo
    enddo

    ! Cab_a = 0.5*I - K*Q
    ! Cab_b = 0.5*I + K*Q
    Cab_a = 0.5d0*I3 - K*Q
    Cab_b = 0.5d0*I3 + K*Q

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
    real(kind=kreal) :: normals_tri6n(3,3)
    integer(kind=kint) :: i, j

    P_matrix = 0.0d0

    ! Identity matrix
    I3 = 0.0d0
    I3(1,1) = 1.0d0
    I3(2,2) = 1.0d0
    I3(3,3) = 1.0d0

    select case(etype)
    case(fe_tri3n)
      call getShapeFunc(fe_tri6n, lpos, N)
      call reorder_normals_tri3n_to_tri6n(vertex_normals, normals_tri6n)
      
      ! Compute Cab for fe_tri6n edges: (st-xi), (xi-et), (et-st)
      call compute_Cab(normals_tri6n(:,1), normals_tri6n(:,2), Cab_31_3, Cab_31_1)
      call compute_Cab(normals_tri6n(:,2), normals_tri6n(:,3), Cab_12_1, Cab_12_2)
      call compute_Cab(normals_tri6n(:,3), normals_tri6n(:,1), Cab_23_2, Cab_23_3)

      ! P matrices for fe_tri6n: st, xi, et
      P1 = N(1) * I3 + N(4) * Cab_31_3 + N(6) * Cab_23_3
      P2 = N(2) * I3 + N(4) * Cab_31_1 + N(5) * Cab_12_1
      P3 = N(3) * I3 + N(5) * Cab_12_2 + N(6) * Cab_23_2

      ! Assemble P in fe_tri3n order: xi, et, st
      P_matrix(1:3, 1:3) = P2
      P_matrix(1:3, 4:6) = P3
      P_matrix(1:3, 7:9) = P1

    case default
      write(*,*) "Error: compute_interpolation_matrix_P - Unsupported element type for Nagata patch.",etype
      stop 
    end select

  end subroutine compute_interpolation_matrix_P

  subroutine create_intermediate_points( surf, currpos, contact_name )
    type(tSurfElement), intent(inout) :: surf(:)    !< surface elements
    real(kind=kreal), intent(in) :: currpos(:)      !< current coordinate of all nodes
    character(len=HECMW_NAME_LEN), intent(in) :: contact_name   !< name of contact surface

    integer(kind=kint) :: i
    
    if (size(surf) == 0) return

    !$omp parallel do private(i)
    do i = 1, size(surf)
      call create_intermediate_points_single(surf(i), currpos)
    enddo
    !$omp end parallel do

    if( DEBUG > 0 ) then
      call print_surface_elements_to_vtk(surf, currpos, contact_name)
    end if  

  end subroutine create_intermediate_points

  subroutine create_intermediate_points_single( surf, currpos )
    type(tSurfElement), intent(inout) :: surf    !< surface element
    real(kind=kreal), intent(in) :: currpos(:)      !< current coordinate of all nodes
    
    integer(kind=kint) :: i, etype
    real(kind=kreal) :: elem(3,3)
    real(kind=kreal), pointer :: vertex_normals(:,:)
    real(kind=kreal) :: Cab_12_1(3,3), Cab_12_2(3,3)
    real(kind=kreal) :: Cab_23_2(3,3), Cab_23_3(3,3)
    real(kind=kreal) :: Cab_31_3(3,3), Cab_31_1(3,3)

    etype = surf%etype
    vertex_normals => surf%vertex_normals

    select case(etype)
    case(fe_tri3n)

      if (.not. associated(surf%intermediate_points)) then
        allocate(surf%intermediate_points(3,3))
      end if

      do i=1,3
        elem(1:3,i) = currpos(3*(surf%nodes(i)-1)+1 : 3*surf%nodes(i))
      enddo

      ! Compute Cab for 3 edges: (1,2), (2,3), (3,1)
      call compute_Cab(vertex_normals(:,1), vertex_normals(:,2), Cab_12_1, Cab_12_2)
      call compute_Cab(vertex_normals(:,2), vertex_normals(:,3), Cab_23_2, Cab_23_3)
      call compute_Cab(vertex_normals(:,3), vertex_normals(:,1), Cab_31_3, Cab_31_1)

      ! intermediate point
      surf%intermediate_points(:,1) = matmul(Cab_12_1(1:3,1:3), elem(1:3,1)) + matmul(Cab_12_2(1:3,1:3), elem(1:3,2))
      surf%intermediate_points(:,2) = matmul(Cab_23_2(1:3,1:3), elem(1:3,2)) + matmul(Cab_23_3(1:3,1:3), elem(1:3,3))
      surf%intermediate_points(:,3) = matmul(Cab_31_3(1:3,1:3), elem(1:3,3)) + matmul(Cab_31_1(1:3,1:3), elem(1:3,1))
    case default
      write(*,*) "Error: create_intermediate_points - Unsupported element type for Nagata patch.",etype
      stop 
    end select

  end subroutine create_intermediate_points_single
  
  subroutine print_surface_elements_to_vtk( surf, currpos, contact_name )
    type(tSurfElement), intent(in) :: surf(:)    !< surface elements
    real(kind=kreal), intent(in) :: currpos(:)   !< current coordinate of all nodes
    character(len=HECMW_NAME_LEN), intent(in) :: contact_name   !< name of contact surface

    integer(kind=kint), parameter :: VTK_UNIT = 99
    integer(kind=kint) :: i, j, nn, n_tri, n_quad, n_elem, n_points, n_cells_size, node_id, inode
    real(kind=kreal) :: coord(3), normal(3)
    character(len=256) :: filename

    ! Count valid elements
    n_tri = 0
    n_quad = 0
    do i = 1, size(surf)
      if (surf(i)%etype == fe_tri3n) then
        n_tri = n_tri + 1
      else if (surf(i)%etype == fe_quad4n) then
        n_quad = n_quad + 1
      endif
    enddo

    n_elem = n_tri + n_quad
    if (n_elem == 0) return

    ! Calculate total points and cells size
    n_points = n_tri * 6 + n_quad * 8
    n_cells_size = n_tri * 7 + n_quad * 9  ! tri: 1+6, quad: 1+8

    ! Open VTK file
    filename = trim(contact_name) // '_nagata_debug.vtk'
    open(unit=VTK_UNIT, file=filename, status='replace', action='write', form='formatted')

    ! Write header
    write(VTK_UNIT, '(A)') '# vtk DataFile Version 3.0'
    write(VTK_UNIT, '(A,A)') 'Nagata patch debug: ', trim(contact_name)
    write(VTK_UNIT, '(A)') 'ASCII'
    write(VTK_UNIT, '(A)') 'DATASET UNSTRUCTURED_GRID'
    write(VTK_UNIT, *)

    ! Write POINTS
    write(VTK_UNIT, '(A,I0,A)') 'POINTS ', n_points, ' float'
    do i = 1, size(surf)
      if (surf(i)%etype == fe_tri3n .or. surf(i)%etype == fe_quad4n) then
        nn = size(surf(i)%nodes)
        ! Output vertex coordinates
        do j = 1, nn
          inode = surf(i)%nodes(j)
          coord(1:3) = currpos(3*(inode-1)+1 : 3*inode)
          write(VTK_UNIT, '(3(E15.7,1X))') coord(1), coord(2), coord(3)
        enddo
        ! Output intermediate point coordinates
        do j = 1, nn
          coord(1:3) = surf(i)%intermediate_points(1:3, j)
          write(VTK_UNIT, '(3(E15.7,1X))') coord(1), coord(2), coord(3)
        enddo
      endif
    enddo
    write(VTK_UNIT, *)

    ! Write CELLS
    write(VTK_UNIT, '(A,I0,1X,I0)') 'CELLS ', n_elem, n_cells_size
    node_id = 0
    do i = 1, size(surf)
      if (surf(i)%etype == fe_tri3n) then
        write(VTK_UNIT, '(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
          6, node_id, node_id+1, node_id+2, node_id+3, node_id+4, node_id+5
        node_id = node_id + 6
      else if (surf(i)%etype == fe_quad4n) then
        write(VTK_UNIT, '(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
          8, node_id, node_id+1, node_id+2, node_id+3, node_id+4, node_id+5, node_id+6, node_id+7
        node_id = node_id + 8
      endif
    enddo
    write(VTK_UNIT, *)

    ! Write CELL_TYPES
    write(VTK_UNIT, '(A,I0)') 'CELL_TYPES ', n_elem
    do i = 1, size(surf)
      if (surf(i)%etype == fe_tri3n) then
        write(VTK_UNIT, '(I0)') 22  ! VTK_QUADRATIC_TRIANGLE
      else if (surf(i)%etype == fe_quad4n) then
        write(VTK_UNIT, '(I0)') 23  ! VTK_QUADRATIC_QUAD
      endif
    enddo
    write(VTK_UNIT, *)

    ! Write POINT_DATA - VECTORS (vertex normals)
    write(VTK_UNIT, '(A,I0)') 'POINT_DATA ', n_points
    write(VTK_UNIT, '(A)') 'VECTORS vertex_normal float'
    do i = 1, size(surf)
      if (surf(i)%etype == fe_tri3n .or. surf(i)%etype == fe_quad4n) then
        nn = size(surf(i)%nodes)
        ! Output vertex normals
        do j = 1, nn
          normal(1:3) = surf(i)%vertex_normals(1:3, j)
          write(VTK_UNIT, '(3(E15.7,1X))') normal(1), normal(2), normal(3)
        enddo
        ! Output zero normals for intermediate points
        do j = 1, nn
          write(VTK_UNIT, '(A)') '0.0 0.0 0.0'
        enddo
      endif
    enddo

    close(VTK_UNIT)

    write(*,'(A,A)') 'VTK debug file written: ', trim(filename)

  end subroutine print_surface_elements_to_vtk


end module m_fstr_contact_smoothing
