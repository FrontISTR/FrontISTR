!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Common utilities for contact element calculations
module m_fstr_contact_elem_common
  use hecmw
  use elementInfo
  use mContactDef
  use m_fstr_contact_geom
  use m_fstr_contact_smoothing
  implicit none

  public :: computeTm_Tt
  public :: computeContactMaps_ALag
  public :: computeFrictionForce_ALag
  public :: getTrialFricForceANDcheckFricState
  public :: update_TangentForce

contains

  !> \brief Compute Tm (relative displacement mapping) and optionally Tt (tangential mapping)
  !! This subroutine constructs the mapping matrices based on contact_disc.md section 10:
  !!   Tm = [I_3; -N_1*I_3; -N_2*I_3; ...; -N_n*I_3]  (size: 3 x 3*(nnode+1))
  !!   Tt = Pt * Tm  where Pt = I_3 - n x n            (computed only if fcoeff /= 0)
  subroutine computeTm_Tt(ctState, tSurf, fcoeff, Tm, Tt, smoothing_type)
    implicit none

    ! Input arguments
    type(tContactState), intent(in) :: ctState         !< contact state (contains lpos and direction)
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient (if 0, Tt not computed)
    integer(kind=kint), optional, intent(in) :: smoothing_type  !< kcsNONE or kcsNAGATA

    ! Output arguments
    real(kind=kreal), intent(out)   :: Tm(3, 3*(l_max_surface_node+1)) !< relative displacement mapping matrix
    real(kind=kreal), intent(out)   :: Tt(3, 3*(l_max_surface_node+1)) !< tangential mapping matrix

    ! Local variables
    integer(kind=kint) :: i, j
    integer(kind=kint) :: nnode !< number of nodes of master segment
    integer(kind=kint) :: smoothing
    real(kind=kreal)   :: shapefunc(l_max_surface_node) !< shape functions [N_1, N_2, ..., N_n]
    real(kind=kreal)   :: P_matrix(3, 3*l_max_surface_node)  !< Nagata interpolation matrix
    real(kind=kreal)   :: normal(3)       !< normal vector (unit vector)
    real(kind=kreal)   :: Pt(3,3)         !< tangential projection operator Pt = I - n⊗n

    Tm = 0.0d0
    Tt = 0.0d0

    nnode = size(tSurf%nodes)
    normal(1:3) = ctState%direction(1:3)

    ! Determine smoothing type
    smoothing = kcsNONE
    if (present(smoothing_type)) smoothing = smoothing_type

    ! Construct Tm
    if ( smoothing == kcsNONE ) then
      ! Standard linear shape functions
      call getShapeFunc(tSurf%etype, ctState%lpos, shapefunc)

      ! First block (slave node): Identity matrix I_3
      Tm(1,1) = 1.0d0
      Tm(2,2) = 1.0d0
      Tm(3,3) = 1.0d0

      ! Remaining blocks (master nodes): -N_i * I_3
      do i = 1, nnode
        Tm(1, i*3+1) = -shapefunc(i)
        Tm(2, i*3+2) = -shapefunc(i)
        Tm(3, i*3+3) = -shapefunc(i)
      enddo

    else if ( smoothing == kcsNAGATA ) then
      ! Use Nagata patch interpolation
      call compute_interpolation_matrix_P(tSurf%etype, nnode, ctState%lpos, &
                                          tSurf%vertex_normals, P_matrix)

      ! Tm = [I_3; -P_matrix]
      Tm(1,1) = 1.0d0
      Tm(2,2) = 1.0d0
      Tm(3,3) = 1.0d0
      Tm(1:3, 4:3*(nnode+1)) = -P_matrix(1:3, 1:3*nnode)

    endif

    ! Compute Tt only if friction coefficient is non-zero
    if (fcoeff /= 0.0d0) then
      ! Construct tangential projection operator Pt = I_3 - n x n
      do i = 1, 3
        do j = 1, 3
          if (i == j) then
            Pt(i,j) = 1.0d0 - normal(i)*normal(j)
          else
            Pt(i,j) = -normal(i)*normal(j)
          endif
        enddo
      enddo

      ! Compute Tt = Pt * Tm using matrix multiplication
      Tt(1:3, 1:3*(l_max_surface_node+1)) = matmul(Pt(1:3,1:3), Tm(1:3, 1:3*(l_max_surface_node+1)))
    endif

  end subroutine computeTm_Tt

  !> \brief Compute contact maps for ALag method (normal distribution and tangential displacement maps)
  !! This subroutine constructs the mapping matrices for ALag contact using computeTm_Tt as the core:
  !!   Bn: normal force distribution vector [n; -N_1*n; -N_2*n; ...; -N_n*n] (size: nnode*3+3)
  !!   Ht: tangent displacement map from DispIncreMatrix (size: 2 x nnode*3+3)
  !!   Gt: covariant displacement map Gt = metric * Ht (size: 2 x nnode*3+3)
  !! Where: dxy = Gt * edisp, f3 = Ht^T * fric
  !! 
  !! This routine now uses computeTm_Tt to generate Bn consistently,
  !! enabling future smoothing extensions in one place.
  subroutine computeContactMaps_ALag(ctState, tSurf, ele, &
                                      Bn, metric, Ht, Gt, smoothing_type)
    implicit none

    ! Input arguments
    type(tContactState), intent(in) :: ctState         !< contact state (contains lpos and direction)
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: ele(:,:)        !< master node coords (3, nnode): coord+disp
    integer(kind=kint), optional, intent(in) :: smoothing_type  !< kcsNONE or kcsNAGATA

    ! Output arguments
    real(kind=kreal), intent(out)   :: Bn(:)           !< normal distribution vector (nnode*3+3)
    real(kind=kreal), intent(out)   :: metric(2,2)     !< metric tensor
    real(kind=kreal), intent(out)   :: Ht(:,:)         !< tangent displacement map (2, nnode*3+3)
    real(kind=kreal), intent(out)   :: Gt(:,:)         !< covariant displacement map (2, nnode*3+3)

    ! Local variables
    integer(kind=kint) :: nnode, ndof
    real(kind=kreal)   :: normal(3)       !< normal vector (unit vector)
    real(kind=kreal)   :: tangent(3,2)    !< tangent basis (temporary for DispIncreMatrix)
    real(kind=kreal)   :: Tm(3, 3*(l_max_surface_node+1)) !< relative displacement mapping (from computeTm_Tt)
    real(kind=kreal)   :: Tt(3, 3*(l_max_surface_node+1)) !< tangential mapping (unused for ALag)

    ! Initialize arrays to zero for safety
    Bn(:) = 0.0d0
    metric(:,:) = 0.0d0
    Ht(:,:) = 0.0d0
    Gt(:,:) = 0.0d0
    Tm(:,:) = 0.0d0
    Tt(:,:) = 0.0d0

    nnode = size(tSurf%nodes)
    ndof = nnode*3 + 3

    ! Call computeTm_Tt as the core routine (fcoeff=0 since we don't need Tt for ALag)
    call computeTm_Tt(ctState, tSurf, 0.0d0, Tm, Tt, smoothing_type)

    ! Get normal vector from contact state
    normal(1:3) = ctState%direction(1:3)

    ! Construct normal distribution vector Bn from Tm
    ! Bn = Tm^T * normal
    Bn(1:ndof) = matmul(transpose(Tm(1:3, 1:ndof)), normal(1:3))

    ! Call DispIncreMatrix to get tangent displacement map and metric
    ! DispIncreMatrix returns: tangent basis, metric tensor, and displacement increment matrix
    call DispIncreMatrix(ctState%lpos(1:2), tSurf%etype, nnode, ele(1:3,1:nnode), &
                         tangent, metric, Ht(1:2,1:ndof))

    ! Compute covariant displacement map: Gt = metric * Ht
    ! This gives: dxy = Gt * edisp = metric * (Ht * edisp)
    Gt(1:2,1:ndof) = matmul(metric(1:2,1:2), Ht(1:2,1:ndof))

  end subroutine computeContactMaps_ALag

  !> \brief Compute friction force for ALag method
  !! Given tangent maps Ht, Gt and displacement increment, compute the friction force vector.
  !! Handles both stick and slip states.
  !!
  !! Slip criterion uses 2D local norm: ||fric||_g = sqrt(fric^T * inv(metric) * fric)
  !! This ensures consistency between force calculation and stiffness matrix.
  !!
  !! Input:
  !!   ctState: contact state (contains multiplier(2:3) for tangent Lagrange multiplier)
  !!   fcoeff: friction coefficient
  !!   lambda_n: normal Lagrange multiplier (for slip scaling)
  !!   metric: metric tensor (2x2) for 2D local tangent space
  !!   Ht: tangent displacement map (2 x edof)
  !!   Gt: covariant displacement map (2 x edof)
  !!   edisp: displacement increment vector (edof)
  !!   edof: element DOF size
  !!   update_multiplier: (optional) if true, update ctState%multiplier(2:3) and check state change
  !!   slave_id: (optional) slave node ID for debug print
  !!
  !! Output:
  !!   ctTForce: friction force vector (edof)
  !!   ctState: (if update_multiplier=true) updated multiplier(2:3) and state
  !!   ctchanged: (optional) flag indicating if contact state changed
  !!   norm_trial: (optional) norm of fric_trial in metric space
  !!   alpha: (optional) projection ratio = min(1, r/norm_trial), r=fcoeff*lambda_n
  !!   that: (optional) normalized direction = fric_trial/norm_trial (2D)
  subroutine computeFrictionForce_ALag(ctState, fcoeff, lambda_n, metric, Ht, Gt, edisp, edof, ctTForce, &
                                        mut, update_multiplier, slave_id, ctchanged, &
                                        norm_trial, alpha, that, jump_ratio)
    implicit none

    type(tContactState), intent(inout) :: ctState       !< contact state (inout for multiplier update)
    real(kind=kreal), intent(in)       :: fcoeff        !< friction coefficient
    real(kind=kreal), intent(in)       :: lambda_n      !< normal Lagrange multiplier
    real(kind=kreal), intent(in)       :: metric(2,2)   !< metric tensor (2x2)
    real(kind=kreal), intent(in)       :: Ht(:,:)       !< tangent displacement map (2, edof)
    real(kind=kreal), intent(in)       :: Gt(:,:)       !< covariant displacement map (2, edof)
    real(kind=kreal), intent(in)       :: edisp(:)      !< displacement increment
    integer(kind=kint), intent(in)     :: edof          !< element DOF size
    real(kind=kreal), intent(out)      :: ctTForce(:)   !< friction force vector
    real(kind=kreal), intent(in)       :: mut           !< tangential penalty parameter
    logical, optional, intent(in)      :: update_multiplier !< if true, update multiplier and check state
    integer(kind=kint), optional, intent(in) :: slave_id     !< slave node ID for debug print
    logical, optional, intent(inout)   :: ctchanged     !< flag for state change
    real(kind=kreal), optional, intent(out) :: norm_trial !< norm of fric_trial (for stiffness)
    real(kind=kreal), optional, intent(out) :: alpha      !< projection ratio (for stiffness)
    real(kind=kreal), optional, intent(out) :: that(2)    !< normalized direction (for stiffness)
    real(kind=kreal), optional, intent(out) :: jump_ratio !< stick trial / slip limit force ratio

    real(kind=kreal) :: dxy(2)       !< tangent displacement (2D local)
    real(kind=kreal) :: fric(2)      !< trial friction force (2D local)
    real(kind=kreal) :: f3(edof)     !< friction force (3D nodal)
    real(kind=kreal) :: invmetric(2,2) !< inverse of metric tensor
    real(kind=kreal) :: det          !< determinant of metric
    real(kind=kreal) :: norm_fric    !< norm in 2D local space: sqrt(fric^T * inv(metric) * fric)
    real(kind=kreal) :: tmp(2)       !< temporary for matrix-vector product
    real(kind=kreal) :: alpha_local  !< local copy of projection ratio
    logical          :: do_update    !< internal flag

    ctTForce = 0.0d0
    do_update = .false.
    if(present(jump_ratio)) jump_ratio = 0.0d0
    if(present(update_multiplier)) do_update = update_multiplier

    ! Compute inverse of metric tensor (2x2)
    det = metric(1,1)*metric(2,2) - metric(1,2)*metric(2,1)
    if( abs(det) < 1.0d-20 ) then
      ! Degenerate metric: set friction to zero
      ctTForce = 0.0d0
      return
    endif
    invmetric(1,1) = metric(2,2) / det
    invmetric(2,2) = metric(1,1) / det
    invmetric(1,2) = -metric(1,2) / det
    invmetric(2,1) = -metric(2,1) / det

    ! Compute tangent displacement using Gt map: dxy = Gt * edisp
    dxy(1:2) = matmul( Gt(1:2,1:edof), edisp(1:edof) )

    ! Trial friction force: fric = λ_t + ρ_t * dxy
    fric(1:2) = ctState%multiplier(2:3) + mut*dxy(1:2)

    ! Compute 2D local norm: ||fric||_g = sqrt(fric^T * inv(metric) * fric)
    tmp(1:2) = matmul(invmetric(1:2,1:2), fric(1:2))
    norm_fric = dsqrt( dot_product(fric(1:2), tmp(1:2)) )

    ! Return projection information for stiffness (optional outputs)
    if(present(norm_trial)) norm_trial = norm_fric
    if(present(that) .and. norm_fric > 1.0d-20) then
      that(1:2) = fric(1:2) / norm_fric
    else if(present(that)) then
      that(1:2) = 0.0d0
    endif

    ! Compute projection ratio: alpha = min(1, r/norm_trial)
    if( lambda_n > 0.0d0 .and. norm_fric > 1.0d-20 ) then
      alpha_local = min(1.0d0, fcoeff*lambda_n / norm_fric)
    else
      alpha_local = 0.0d0
    endif
    if(present(alpha)) alpha = alpha_local

    ! Check stick/slip state using projection ratio
    if( do_update .and. lambda_n > 0.0d0 ) then
      ! Update mode: check state and update multiplier
      if( alpha_local < 1.0d0 ) then
        ! Slip state: projected to friction cone
        if( ctState%state == CONTACTSTICK ) then
          if(present(ctchanged)) ctchanged = .true.
          if(present(slave_id)) print *, "Node", slave_id, "to slip state", norm_fric, fcoeff*lambda_n
          if(present(jump_ratio) .and. lambda_n > 1.0e-20 .and. norm_fric > 1.0e-20) &
            &  jump_ratio = norm_fric / (fcoeff*lambda_n)
        endif
        ctState%state = CONTACTSLIP
        fric(1:2) = fric(1:2) * alpha_local
      else
        ! Stick state
        if( ctState%state == CONTACTSLIP ) then
          if(present(ctchanged)) ctchanged = .true.
          if(present(slave_id)) print *, "Node", slave_id, "to stick state", norm_fric, fcoeff*lambda_n
        endif
        ctState%state = CONTACTSTICK
      endif
      ! Update multiplier
      ctState%multiplier(2:3) = fric(1:2)
    else if( lambda_n <= 0.0d0 ) then
      ! No normal force: no friction
      fric(1:2) = 0.0d0
    else
      ! Read-only mode: just scale if alpha < 1
      if( alpha_local < 1.0d0 ) then
        fric(1:2) = fric(1:2) * alpha_local
      endif
    endif

    ! Compute friction force vector using Ht map: f3 = Ht^T * fric
    f3(1:edof) = matmul( transpose(Ht(1:2,1:edof)), fric(1:2) )

    ! Distribute friction force (negative for equilibrium)
    ctTForce(1:edof) = -f3(1:edof)

  end subroutine computeFrictionForce_ALag

  !> \brief This subroutine calculates trial friction force and checks friction state
  subroutine getTrialFricForceANDcheckFricState(ctstate,fcoeff,tPenalty,lagrange)

    type(tContactState) :: ctState !< type tContactState
    real(kind=kreal)   :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier

    integer(kind=kint) :: i
    real(kind=kreal)   :: tf_yield

    do i = 1, 3
      ctstate%tangentForce_trial(i) = ctstate%tangentForce1(i) + tPenalty*ctstate%reldisp(i)
    enddo

    tf_yield = fcoeff*dabs(lagrange)
    if(ctstate%state == contactSlip) tf_yield =0.99d0*tf_yield
    if( dsqrt(dot_product(ctstate%tangentForce_trial,ctstate%tangentForce_trial)) <= tf_yield ) then
      ctstate%state = contactStick
    else
      ctstate%state = contactSlip
    endif

  end subroutine getTrialFricForceANDcheckFricState

  !> This subroutine find the projection of a slave point onto master surface
  subroutine update_TangentForce(etype,nn,elemt0,elemt,cstate)
    integer, intent(in)                 :: etype         !< surface element type
    integer, intent(in)                 :: nn            !< number of elemental nodes
    real(kind=kreal),intent(in)         :: elemt0(3,nn)  !< nodes coordinates of surface element at t
    real(kind=kreal),intent(in)         :: elemt(3,nn)   !< nodes coordinates of surface element at t+dt
    type(tContactState), intent(inout)  :: cstate        !< Recorde of contact information

    integer           ::  i
    real(kind=kreal)  ::  tangent0(3,2), tangent(3,2)    ! base vectors in tangent space
    real(kind=kreal)  ::  coeff(2), norm, norm_tmp
    real(kind=kreal)  ::  tangentForce_tmp(3)

    call TangentBase( etype, nn, cstate%lpos(1:2), elemt0, tangent0 )
    call TangentBase( etype, nn, cstate%lpos(1:2), elemt, tangent )

    !project tangentforce to base vector tangent0
    do i=1,2
      coeff(i) = dot_product(cstate%tangentForce(1:3),tangent0(1:3,i))
      coeff(i) = coeff(i)/dot_product(tangent0(1:3,i),tangent0(1:3,i))
    enddo
    tangentForce_tmp(1:3) = coeff(1)*tangent0(1:3,1) + coeff(2)*tangent0(1:3,2)
    norm_tmp = dsqrt(dot_product(tangentForce_tmp,tangentForce_tmp))
    !adjust tangent force of slave point which moved over element boundary
    if( norm_tmp > 1.d-6 ) then
      norm = dsqrt(dot_product(cstate%tangentForce,cstate%tangentForce))
      coeff(1:2) = (norm/norm_tmp)*coeff(1:2)
    end if

    !set rotated tangentforce to tangentforce1
    cstate%tangentForce1(1:3) = coeff(1)*tangent(1:3,1) + coeff(2)*tangent(1:3,2)

  end subroutine update_TangentForce

end module m_fstr_contact_elem_common
