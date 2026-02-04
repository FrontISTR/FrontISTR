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

contains

  subroutine getContactStiffness_Alag(cstate, etype, nnode, ele, fcoeff, symm, stiff, force)

    type(tContactState), intent(in) :: cstate          !< contact state
    integer, intent(in)             :: etype           !< type of contacting surface
    integer, intent(in)             :: nnode           !< number of elemental nodes
    real(kind=kreal), intent(in)    :: ele(3,nnode)    !< coord of surface element
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient
    logical, intent(in)             :: symm            !< symmetricalize
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force direction

    integer          :: i, j
    real(kind=kreal) :: shapefunc(nnode)
    real(kind=kreal) :: N(nnode*3+3), dispmat(2,nnode*3+3)
    real(kind=kreal) :: metric(2,2)
    real(kind=kreal) :: det, inverse(2,2), ff(2), cff(2)
    real(kind=kreal) :: dum11(nnode*3+3,nnode*3+3), dum12(nnode*3+3,nnode*3+3)
    real(kind=kreal) :: dum21(nnode*3+3,nnode*3+3), dum22(nnode*3+3,nnode*3+3)
    real(kind=kreal) :: tangent(3,2)

    call getShapeFunc( etype, cstate%lpos(1:2), shapefunc )
    N(1:3) = cstate%direction(1:3)
    do i = 1, nnode
      N(i*3+1:i*3+3) = -shapefunc(i) * cstate%direction(1:3)
    enddo
    do j = 1, nnode*3+3
      do i = 1, nnode*3+3
        stiff(i,j) = mu * N(i) * N(j)
      enddo
    enddo
    force(1:nnode*3+3) = N(:)

    if( fcoeff /= 0.d0 ) &
      call DispIncreMatrix( cstate%lpos(1:2), etype, nnode, ele, tangent, metric, dispmat )

    ! frictional component
    if( fcoeff /= 0.d0 ) then
      forall(i=1:nnode*3+3, j=1:nnode*3+3)
        dum11(i,j) = mut * dispmat(1,i) * dispmat(1,j)
        dum12(i,j) = mut * dispmat(1,i) * dispmat(2,j)
        dum21(i,j) = mut * dispmat(2,i) * dispmat(1,j)
        dum22(i,j) = mut * dispmat(2,i) * dispmat(2,j)
      end forall
      stiff(1:nnode*3+3,1:nnode*3+3) = stiff(1:nnode*3+3,1:nnode*3+3) &
        + metric(1,1)*dum11 + metric(1,2)*dum12 &
        + metric(2,1)*dum21 + metric(2,2)*dum22

      if( cstate%state == CONTACTSLIP ) then
        det = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
        if( det == 0.d0 ) stop "Math error in contact stiff calculation"
        inverse(1,1) = metric(2,2) / det
        inverse(2,1) = -metric(2,1) / det
        inverse(1,2) = -metric(1,2) / det
        inverse(2,2) = metric(1,1) / det
        ff(:) = cstate%multiplier(2:3)
        cff(:) = matmul( inverse, ff )
        ff(:) = ff(:) / dsqrt( ff(1)*ff(1) + ff(2)*ff(2) )
        cff(:) = cff(:) / dsqrt( cff(1)*cff(1) + cff(2)*cff(2) )
        stiff(1:nnode*3+3,1:nnode*3+3) = stiff(1:nnode*3+3,1:nnode*3+3) - &
          ( cff(1)*ff(1)*metric(1,1) + cff(2)*ff(1)*metric(1,2) )*dum11 - &
          ( cff(2)*ff(2)*metric(1,2) + cff(1)*ff(2)*metric(1,1) )*dum21 - &
          ( cff(1)*ff(1)*metric(1,2) + cff(2)*ff(1)*metric(2,2) )*dum12 - &
          ( cff(2)*ff(2)*metric(2,2) + cff(1)*ff(2)*metric(1,2) )*dum22
      endif
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
    real(kind=kreal)   :: shapefunc(l_max_surface_node)
    real(kind=kreal)   :: elemcrd(3, l_max_elem_node) !< master node coords (coord+disp, for DispIncreMatrix)
    real(kind=kreal)   :: edisp(3*l_max_elem_node+3) !< displacement increment
    real(kind=kreal)   :: elemg(3), dg(3) !< interpolated master position; gap vector
    real(kind=kreal)   :: dgn, nrlforce !< normal gap; normal force
    real(kind=kreal)   :: tangent(3,2), metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)   :: dxi(2), dxy(2), fric(2), f3(3*l_max_elem_node+3)

    nnode = size(tSurf%nodes)

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    normal(1:3) = ctState%direction(1:3)

    call getShapeFunc( tSurf%etype, ctState%lpos(1:2), shapefunc )

    ! Compute interpolated master position using current coordinates
    ! elemg = sum_j shapefunc(j) * ndCoord(master_j)
    elemg = 0.d0
    do j = 1, nnode
      elemg(1:3) = elemg(1:3) + shapefunc(j) * ndCoord(j*3+1:j*3+3)
    enddo

    ! Gap vector: slave current position - interpolated master position
    dg(1:3) = ndCoord(1:3) - elemg(1:3)

    ! Normal gap
    dgn = dot_product( normal, dg )

    ! Normal force: multiplier + penalty * gap
    nrlforce = ctState%multiplier(1) + mu*dgn

    ! Slave node force: -nrlforce * normal
    ctNForce(1:3) = -nrlforce * normal(1:3)

    ! Master node forces: +nrlforce * shapefunc(j) * normal
    do j = 1, nnode
      ctNForce(j*3+1:j*3+3) = nrlforce * shapefunc(j) * normal(1:3)
    enddo

    ! Lagrange row (not used in ALagrange, set to 0)
    ctNForce((nnode+1)*3+1) = 0.d0

    if( fcoeff == 0.d0 ) return

    ! --- Tangent component ---

    ! Prepare elemcrd = ndCoord - ndDu (i.e., coord + disp) for DispIncreMatrix
    do j = 1, nnode
      elemcrd(1:3, j) = ndCoord(j*3+1:j*3+3) - ndDu(j*3+1:j*3+3)
    enddo

    ! Prepare edisp from ndDu
    edisp(1:3) = ndDu(1:3)  ! slave
    do j = 1, nnode
      edisp(j*3+1:j*3+3) = ndDu(j*3+1:j*3+3)  ! master nodes
    enddo

    call DispIncreMatrix( ctState%lpos(1:2), tSurf%etype, nnode, elemcrd(:,1:nnode), tangent, metric, dispmat )

    dxi(1) = dot_product( dispmat(1,1:nnode*3+3), edisp(1:nnode*3+3) )
    dxi(2) = dot_product( dispmat(2,1:nnode*3+3), edisp(1:nnode*3+3) )
    dxy(1:2) = matmul( metric, dxi )

    fric(1:2) = ctState%multiplier(2:3) + mut*dxy(1:2)
    f3(1:nnode*3+3) = fric(1)*dispmat(1,1:nnode*3+3) + fric(2)*dispmat(2,1:nnode*3+3)

    if( ctState%state == CONTACTSLIP ) then
      dgn = dsqrt( f3(1)*f3(1) + f3(2)*f3(2) + f3(3)*f3(3) )
      f3(1:nnode*3+3) = f3(1:nnode*3+3) * fcoeff * ctState%multiplier(1) / dgn
    endif

    ! Slave node friction force
    ctTForce(1:3) = -f3(1:3)

    ! Master node friction forces
    do j = 1, nnode
      ctTForce(j*3+1:j*3+3) = -f3(j*3+1:j*3+3)
    enddo

    ! Lagrange row (not used in ALagrange, set to 0)
    ctTForce((nnode+1)*3+1) = 0.d0

  end subroutine getContactNodalForce_Alag

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
