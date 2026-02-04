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
  implicit none

  public :: computeTm_Tt
  public :: getTrialFricForceANDcheckFricState
  public :: update_TangentForce

contains

  !> \brief Compute Tm (relative displacement mapping) and optionally Tt (tangential mapping)
  !! This subroutine constructs the mapping matrices based on contact_disc.md section 10:
  !!   Tm = [I_3; -N_1*I_3; -N_2*I_3; ...; -N_n*I_3]  (size: 3 x 3*(nnode+1))
  !!   Tt = Pt * Tm  where Pt = I_3 - n x n            (computed only if fcoeff /= 0)
  subroutine computeTm_Tt(ctState, tSurf, fcoeff, Tm, Tt)
    implicit none

    ! Input arguments
    type(tContactState), intent(in) :: ctState         !< contact state (contains lpos and direction)
    type(tSurfElement), intent(in)  :: tSurf           !< surface element structure
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient (if 0, Tt not computed)

    ! Output arguments
    real(kind=kreal), intent(out)   :: Tm(3, 3*(l_max_surface_node+1)) !< relative displacement mapping matrix
    real(kind=kreal), intent(out)   :: Tt(3, 3*(l_max_surface_node+1)) !< tangential mapping matrix

    ! Local variables
    integer(kind=kint) :: i, j
    integer(kind=kint) :: nnode !< number of nodes of master segment
    real(kind=kreal)   :: shapefunc(l_max_surface_node) !< shape functions [N_1, N_2, ..., N_n]
    real(kind=kreal)   :: normal(3)       !< normal vector (unit vector)
    real(kind=kreal)   :: Pt(3,3)         !< tangential projection operator Pt = I - n⊗n

    Tm = 0.0d0
    Tt = 0.0d0

    nnode = size(tSurf%nodes)

    call getShapeFunc(tSurf%etype, ctState%lpos(:), shapefunc)

    ! Get normal vector from contact state
    normal(1:3) = ctState%direction(1:3)

    ! Construct Tm = [I_3; -N_1*I_3; -N_2*I_3; ...; -N_n*I_3]
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
