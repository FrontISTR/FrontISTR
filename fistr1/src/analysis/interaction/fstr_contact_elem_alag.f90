!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Alag method implementations for contact element calculations
module m_fstr_contact_elem_alag
  use hecmw
  use elementInfo
  use mContactDef
  use mSurfElement
  use m_fstr_contact_geom
  use m_fstr_contact_elem_common
  implicit none

  public :: getContactStiffness_Alag
  public :: getContactNodalForce_Alag
  public :: getTiedStiffness_Alag
  public :: getTiedNodalForce_Alag
  public :: updateContactMultiplier_Alag

contains

  subroutine getContactStiffness_Alag(cstate, tSurf, ele, fcoeff, symm, stiff, force)

    type(tContactState), intent(inout) :: cstate       !< contact state (inout for projection info)
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: ele(:,:)        !< coord of surface element
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient
    logical, intent(in)             :: symm            !< symmetricalize
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force direction

    integer          :: i, j, nnode
    real(kind=kreal) :: Bn(size(tSurf%nodes)*3+3), Ht(2,size(tSurf%nodes)*3+3), Gt(2,size(tSurf%nodes)*3+3)
    real(kind=kreal) :: metric(2,2)
    real(kind=kreal) :: A(2,2)       !< 2D local tangent operator
    real(kind=kreal) :: norm_trial, alpha_proj, that(2)
    real(kind=kreal) :: K_fric(size(tSurf%nodes)*3+3,size(tSurf%nodes)*3+3)  !< friction stiffness
    real(kind=kreal) :: HtA(2,size(tSurf%nodes)*3+3), tmp_vec(2)
    real(kind=kreal) :: dummy_force(size(tSurf%nodes)*3+3)  !< dummy for computeFrictionForce_ALag
    real(kind=kreal) :: zero_disp(size(tSurf%nodes)*3+3)   !< zero displacement for trial friction

    nnode = size(tSurf%nodes)

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(cstate, tSurf, ele, Bn, metric, Ht, Gt)

    ! Normal stiffness: stiff = mu * Bn * Bn^T
    do j = 1, nnode*3+3
      do i = 1, nnode*3+3
        stiff(i,j) = mu * Bn(i) * Bn(j)
      enddo
    enddo
    force(1:nnode*3+3) = Bn(:)

    ! frictional component
    if( fcoeff /= 0.d0 ) then
      ! Compute trial friction info for consistent tangent
      ! Use zero displacement to get projection parameters from current multiplier state
      zero_disp = 0.0d0  ! zero displacement increment
      call computeFrictionForce_ALag(cstate, fcoeff, cstate%multiplier(1), metric, &
                                      Ht, Gt, zero_disp, nnode*3+3, dummy_force, &
                                      norm_trial=norm_trial, alpha=alpha_proj, that=that)

      ! Construct 2D local tangent operator A
      if( cstate%multiplier(1) <= 0.0d0 .or. alpha_proj <= 1.0d-20 ) then
        ! No friction: A = 0
        A = 0.0d0
      else if( alpha_proj >= 0.999d0 .or. norm_trial < 1.0d-20 ) then
        ! Stick: A = mut*I
        A(1,1) = mut
        A(2,2) = mut
        A(1,2) = 0.0d0
        A(2,1) = 0.0d0
      else
        ! Slip: A = alpha*mut*(I - that⊗that)
        A(1,1) = alpha_proj * mut * (1.0d0 - that(1)*that(1))
        A(1,2) = alpha_proj * mut * (-that(1)*that(2))
        A(2,1) = alpha_proj * mut * (-that(2)*that(1))
        A(2,2) = alpha_proj * mut * (1.0d0 - that(2)*that(2))
      endif

      ! Compute friction stiffness: K_fric = Ht^T * A * Gt
      ! First: HtA = A^T * Ht (2 x edof)
      HtA(1:2,1:nnode*3+3) = matmul(transpose(A), Ht(1:2,1:nnode*3+3))
      ! Then: K_fric = Ht^T * HtA = Ht^T * A^T * Ht (but we want Ht^T * A * Gt)
      ! Correct: K_fric(i,j) = Ht(:,i)^T * A * Gt(:,j)
      do j = 1, nnode*3+3
        tmp_vec = matmul(A, Gt(1:2,j))
        do i = 1, nnode*3+3
          K_fric(i,j) = dot_product(Ht(1:2,i), tmp_vec)
        enddo
      enddo

      stiff(1:nnode*3+3,1:nnode*3+3) = stiff(1:nnode*3+3,1:nnode*3+3) + K_fric(1:nnode*3+3,1:nnode*3+3)
    endif

  end subroutine getContactStiffness_Alag

  subroutine getContactNodalForce_Alag(ctState,tSurf,ndCoord,ndDu,tPenalty,fcoeff,lagrange,ctNForce,ctTForce,cflag)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint) :: nnode !< number of nodes of master segment
    integer(kind=kint) :: j
    real(kind=kreal)   :: fcoeff, tPenalty !< friction coefficient; tangential penalty (tPenalty not used, mut used instead)
    real(kind=kreal)   :: lagrange !< not used for ALagrange (kept for interface compatibility)
    real(kind=kreal)   :: ndCoord(:), ndDu(:) !< nodal coordinates (coord+disp+ddisp); nodal displacement increment (ddisp)
    real(kind=kreal)   :: ctNForce(:) !< contact normal force vector
    real(kind=kreal)   :: ctTForce(:) !< contact tangential force vector
    logical            :: cflag  !< not used for ALagrange (kept for interface compatibility)

    real(kind=kreal)   :: normal(3) !< normal vector at target point
    real(kind=kreal)   :: Bn(3*l_max_elem_node+3) !< normal distribution vector
    real(kind=kreal)   :: Ht(2,3*l_max_elem_node+3), Gt(2,3*l_max_elem_node+3) !< tangent and covariant maps
    real(kind=kreal)   :: elemcrd(3, l_max_elem_node) !< master node coords (coord+disp, for computeContactMaps_ALag)
    real(kind=kreal)   :: edisp(3*l_max_elem_node+3) !< displacement increment
    real(kind=kreal)   :: dgn, nrlforce !< normal gap; normal force
    real(kind=kreal)   :: metric(2,2)
    integer(kind=kint) :: edof  !< element vector size (nnode*3+3)

    nnode = size(tSurf%nodes)
    edof = nnode*3+3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    normal(1:3) = ctState%direction(1:3)

    ! Prepare elemcrd = ndCoord - ndDu (i.e., coord + disp) for computeContactMaps_ALag
    do j = 1, nnode
      elemcrd(1:3, j) = ndCoord(j*3+1:j*3+3) - ndDu(j*3+1:j*3+3)
    enddo

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(ctState, tSurf, elemcrd(:,1:nnode), &
                                  Bn, metric, Ht, Gt)

    ! Normal gap: dgn = Bn^T * ndCoord (using normal distribution vector)
    dgn = dot_product( Bn(1:edof), ndCoord(1:edof) )

    ! Normal force: multiplier + penalty * gap
    nrlforce = ctState%multiplier(1) + mu*dgn

    ! Distribute normal force using Bn: ctNForce = -nrlforce * Bn
    ctNForce(1:edof) = -nrlforce * Bn(1:edof)

    ! Lagrange row (not used in ALagrange, set to 0)
    ctNForce((nnode+1)*3+1) = 0.d0

    if( fcoeff == 0.d0 ) return

    ! --- Tangent component ---

    ! Prepare edisp from ndDu
    edisp(1:3) = ndDu(1:3)  ! slave
    do j = 1, nnode
      edisp(j*3+1:j*3+3) = ndDu(j*3+1:j*3+3)  ! master nodes
    enddo

    ! Compute friction force using common routine
    call computeFrictionForce_ALag(ctState, fcoeff, ctState%multiplier(1), metric, &
                                    Ht, Gt, edisp, edof, ctTForce)

    ! Lagrange row (not used in ALagrange, set to 0)
    ctTForce((nnode+1)*3+1) = 0.d0

  end subroutine getContactNodalForce_Alag

  subroutine updateContactMultiplier_Alag(ctState,ndLocal,coord,disp,ddisp,mu,mut,fcoeff,etype,lgnt,ctchanged,ctNForce,ctTForce)

    type(tContactState), intent(inout)   :: ctState             !< contact state
    integer(kind=kint), intent(in)       :: ndLocal(:)          !< global node numbers (slave + master)
    real(kind=kreal), intent(in)         :: coord(:)            !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)             !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)            !< disp till current substep
    real(kind=kreal), intent(in)         :: mu, mut             !< penalty parameters
    real(kind=kreal), intent(in)         :: fcoeff              !< friction coefficient
    integer(kind=kint), intent(in)       :: etype               !< element type
    real(kind=kreal), intent(inout)      :: lgnt(2)             !< convergence metrics
    logical, intent(inout)               :: ctchanged           !< contact state changed flag
    real(kind=kreal), intent(out)        :: ctNForce(:)         !< contact normal force vector
    real(kind=kreal), intent(out)        :: ctTForce(:)         !< contact tangential force vector

    integer(kind=kint)  :: nnode !< number of master nodes
    integer(kind=kint)  :: slave, j
    real(kind=kreal)    :: Bn(3*l_max_elem_node+3) !< normal distribution vector
    real(kind=kreal)    :: Ht(2,3*l_max_elem_node+3), Gt(2,3*l_max_elem_node+3) !< tangent and covariant maps
    real(kind=kreal)    :: elemcrd(3,l_max_elem_node) !< master node coords (coord+disp, for computeContactMaps_ALag)
    real(kind=kreal)    :: curpos(3*l_max_elem_node+3) !< current positions (coord+disp+ddisp)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3) !< displacement increment
    real(kind=kreal)    :: dgn, nrlforce !< normal gap; normal force
    real(kind=kreal)    :: metric(2,2)
    real(kind=kreal)    :: dxy(2)        !< for convergence check
    integer(kind=kint)  :: edof  !< element vector size (nnode*3+3)
    type(tSurfElement)  :: tSurf_tmp  !< temporary surface element for computeContactMaps_ALag

    nnode = size(ndLocal) - 1
    slave = ndLocal(1)
    edof = nnode*3+3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Prepare elemcrd (coord+disp) and current positions (coord+disp+ddisp)
    curpos(1:3) = coord(3*slave-2:3*slave) + disp(3*slave-2:3*slave) + ddisp(3*slave-2:3*slave)
    edisp(1:3) = ddisp(3*slave-2:3*slave)
    do j = 1, nnode
      elemcrd(1:3,j) = coord(3*ndLocal(j+1)-2:3*ndLocal(j+1)) + disp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
      curpos(j*3+1:j*3+3) = elemcrd(1:3,j) + ddisp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
      edisp(j*3+1:j*3+3) = ddisp(3*ndLocal(j+1)-2:3*ndLocal(j+1))
    enddo

    ! Build temporary tSurf for computeContactMaps_ALag
    allocate(tSurf_tmp%nodes(nnode))
    tSurf_tmp%etype = etype
    tSurf_tmp%nodes = ndLocal(2:nnode+1)

    ! Use common mapping routine to compute Bn, metric, Ht, Gt
    call computeContactMaps_ALag(ctState, tSurf_tmp, elemcrd(:,1:nnode), &
                                  Bn, metric, Ht, Gt)

    ! Normal gap: dgn = Bn^T * curpos (using normal distribution vector)
    dgn = dot_product( Bn(1:edof), curpos(1:edof) )

    ! Update multiplier and working distance
    ctState%wkdist = -dgn
    ctState%multiplier(1) = ctState%multiplier(1) - mu*ctState%wkdist
    ctState%distance = ctState%wkdist
    lgnt(1) = lgnt(1) - ctState%wkdist

    ! Normal force: use updated multiplier
    nrlforce = ctState%multiplier(1)

    ! Distribute normal force using Bn: ctNForce = -nrlforce * Bn
    ctNForce(1:edof) = -nrlforce * Bn(1:edof)

    if( fcoeff == 0.d0 ) return

    ! --- Tangent component ---

    ! Compute friction force with multiplier update and state check
    call computeFrictionForce_ALag(ctState, fcoeff, ctState%multiplier(1), metric, &
                                    Ht, Gt, edisp, edof, ctTForce, &
                                    update_multiplier=.true., slave_id=slave, ctchanged=ctchanged)

    ! Tangent displacement for convergence check: use Gt to project curpos directly
    dxy = matmul( Gt(:,1:edof), curpos(1:edof) )
    lgnt(2) = lgnt(2) + dsqrt( dxy(1)*dxy(1) + dxy(2)*dxy(2) )

    ! Clean up temporary surface element
    deallocate(tSurf_tmp%nodes)

  end subroutine updateContactMultiplier_Alag

  subroutine getTiedStiffness_Alag(cstate, etype, nnode, stiff, force)

    type(tContactState), intent(in) :: cstate          !< contact state
    integer, intent(in)             :: etype           !< type of contacting surface
    integer, intent(in)             :: nnode           !< number of elemental nodes
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force direction

    integer          :: i, j, k
    real(kind=kreal) :: shapefunc(nnode)
    real(kind=kreal) :: N(nnode*3+3)

    stiff = 0.d0

    call getShapeFunc( etype, cstate%lpos(1:2), shapefunc )
    N(1) = 1.d0
    N(2:nnode+1) = -shapefunc(1:nnode)

    do j = 1, nnode+1
      do k = 1, nnode+1
        do i = 1, 3
          stiff(3*k-3+i, 3*j-3+i) = mu * N(k) * N(j)
        enddo
      enddo
    enddo
    force(1:nnode*3+3) = N(:)

  end subroutine getTiedStiffness_Alag

  subroutine getTiedNodalForce_Alag(ctState,tSurf,ndu,lagrange,ctNForce,ctTForce)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    real(kind=kreal)   :: ndu(:) !< nodal total displacement (disp+ddisp)
    real(kind=kreal)   :: lagrange !< not used for ALagrange (kept for interface compatibility)
    real(kind=kreal)   :: ctNForce(:)  !< contact force vector
    real(kind=kreal)   :: ctTForce(:)  !< contact force vector (not used for tied)

    integer(kind=kint) :: j
    real(kind=kreal)   :: shapefunc(l_max_surface_node)
    real(kind=kreal)   :: edisp(3*l_max_elem_node+3) !< total displacement
    real(kind=kreal)   :: dg(3) !< gap vector
    real(kind=kreal)   :: nrlforce(3) !< nodal force (3 components)

    nnode = size(tSurf%nodes)

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Prepare total displacement: edisp = ndu
    ! edisp(1:3) = slave node
    edisp(1:3) = ndu(1:3)

    ! edisp(4:6, 7:9, ...) = master nodes
    do j = 1, nnode
      edisp(j*3+1:j*3+3) = ndu(j*3+1:j*3+3)
    enddo

    call getShapeFunc( tSurf%etype, ctState%lpos(1:2), shapefunc )

    ! Gap vector: slave displacement - interpolated master displacement
    ! dg = edisp(slave) - sum_j shapefunc(j) * edisp(master_j)
    dg(1:3) = edisp(1:3)
    do j = 1, nnode
      dg(1:3) = dg(1:3) - shapefunc(j) * edisp(j*3+1:j*3+3)
    enddo

    ! Force: multiplier + penalty * gap (3 components)
    nrlforce(1:3) = ctState%multiplier(1:3) + mu*dg(1:3)

    ! Slave node force: -nrlforce
    ctNForce(1:3) = -nrlforce(1:3)

    ! Master node forces: +shapefunc(j) * nrlforce
    do j = 1, nnode
      ctNForce(j*3+1:j*3+3) = shapefunc(j) * nrlforce(1:3)
    enddo

    ! Lagrange row (not used in ALagrange, set to 0)
    ctNForce((nnode+1)*3+1) = 0.d0
    ctTForce((nnode+1)*3+1) = 0.d0

  end subroutine getTiedNodalForce_Alag

end module m_fstr_contact_elem_alag
