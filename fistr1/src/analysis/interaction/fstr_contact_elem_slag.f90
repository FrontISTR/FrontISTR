!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Slag method implementations for contact element calculations
module m_fstr_contact_elem_slag
  use hecmw
  use elementInfo
  use mContactDef
  use mSurfElement
  use m_fstr_contact_geom
  use m_fstr_contact_elem_common
  implicit none

  public :: getContactStiffness_Slag
  public :: getContactNodalForce_Slag
  public :: getTiedStiffness_Slag
  public :: getTiedNodalForce_Slag

contains

  subroutine getContactStiffness_Slag(ctState,tSurf,iter,tPenalty,fcoeff,lagrange,stiffness)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint)  :: iter
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    integer(kind=kint)  :: i, j
    real(kind=kreal)    :: normal(3)
    real(kind=kreal)    :: shapefunc(l_max_surface_node) !< normal vector at target point; shape functions
    real(kind=kreal)    :: nTm((l_max_surface_node + 1)*3) !< vector
    real(kind=kreal)    :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)    :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)    :: tf_trial(3), length_tft
    real(kind=kreal)    :: tangent(3)
    real(kind=kreal)    :: stiffness(:,:) !< contact stiffness matrix
    real(kind=kreal)    :: Tm(3, 3*(l_max_surface_node+1)), Tt(3, 3*(l_max_surface_node+1)) !< mapping matrices
    real(kind=kreal)    :: q((l_max_surface_node+1)*3)  !< q = Tm^T * tangent (for slip)
    real(kind=kreal)    :: dot_tn  !< tangent · normal (for slip)

    nnode = size(tSurf%nodes)

    stiffness = 0.0d0

    ! Compute Tm and Tt matrices using standard shape functions
    call computeTm_Tt(ctState, tSurf, fcoeff, Tm, Tt)

    normal(1:3) = ctState%direction(1:3)

    ! Compute nTm = Tm^T * n for KKT constraint (Lagrange multiplier row/column)
    nTm(1:(nnode+1)*3) = matmul(transpose(Tm(1:3, 1:3*(nnode+1))), normal)

    i = (nnode+1)*3 + 1
    do j = 1, (nnode+1)*3
      stiffness(i,j) = nTm(j);  stiffness(j,i) = nTm(j)
    enddo

    if( fcoeff /= 0.0d0 ) then
      if( lagrange>0.0d0 .or. iter==1 ) then

        ! Stick friction: K_uu^stick = ρ * Tt^T * Tt
        do i = 1, (nnode+1)*3
          do j = 1, i
            stiffness(i,j) = stiffness(i,j) + tPenalty * dot_product(Tt(1:3,i), Tt(1:3,j))
            ! Fill symmetric part
            if (i /= j) then
              stiffness(j,i) = stiffness(i,j)
            endif
          enddo
        enddo

        if( ctstate%state == contactSlip ) then

          tf_trial(1:3) = ctstate%tangentForce_trial(1:3)
          length_tft = dsqrt(dot_product(tf_trial,tf_trial))
          tangent(1:3) = tf_trial(1:3)/length_tft

          ! Compute q = Tm^T * tangent
          q(1:(nnode+1)*3) = matmul(transpose(Tm(1:3, 1:3*(nnode+1))), tangent)

          ! Slip correction:
          dot_tn = dot_product(tangent, normal)

          do i = 1, (nnode+1)*3
            do j = 1, i
              stiffness(i,j) = stiffness(i,j) + tPenalty * (-q(i)*q(j) + q(i)*nTm(j)*dot_tn)
              if (i /= j) then
                stiffness(j,i) = stiffness(i,j)
              endif
            enddo
          enddo

          ! Scale (u,u) block by (μλ/|t_tr|)
          stiffness(1:(nnode+1)*3,1:(nnode+1)*3) = (fcoeff*lagrange/length_tft) * stiffness(1:(nnode+1)*3,1:(nnode+1)*3)

          ! (u,λ) cross term:)
          stiffness(1:(nnode+1)*3, (nnode+1)*3+1) = stiffness(1:(nnode+1)*3, (nnode+1)*3+1) + fcoeff * q(1:(nnode+1)*3)

        endif
      endif
    endif

  end subroutine getContactStiffness_Slag

  subroutine getContactNodalForce_Slag(ctState,tSurf,ndCoord,ndDu,tPenalty,fcoeff,lagrange,ctNForce,ctTForce,cflag)

    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint) :: nnode !< type of master segment; number of nodes of master segment
    integer(kind=kint) :: j
    real(kind=kreal)   :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)   :: ndCoord(:), ndDu(:) !< nodal coordinates; nodal displacement increment
    real(kind=kreal)   :: normal(3) !< normal vector at target point
    real(kind=kreal)   :: nTm((l_max_surface_node + 1)*3) !< vector
    real(kind=kreal)   :: tf_trial(3), length_tft, tangent(3), tf_final(3)
    real(kind=kreal)       :: ctNForce(:)                     !< contact force vector
    real(kind=kreal)       :: ctTForce(:)                     !< contact force vector
    logical            :: cflag  !< is necessary to update tangentForce_final

    real(kind=kreal)   :: Tm(3, 3*(l_max_surface_node+1)), Tt(3, 3*(l_max_surface_node+1)) !< mapping matrices
    real(kind=kreal)   :: relativeDisp(3) !< relative displacement

    ctState%multiplier(1) = lagrange

    nnode = size(tSurf%nodes)

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Compute Tm matrix (Tt not needed, pass fcoeff=0 to skip)
    call computeTm_Tt(ctState, tSurf, fcoeff, Tm, Tt)

    normal(1:3) = ctState%direction(1:3)

    ! Compute nTm = Tm^T * normal
    nTm(1:(nnode+1)*3) = -matmul(transpose(Tm(1:3, 1:3*(nnode+1))), normal)

    do j = 1, (nnode+1)*3
      ctNForce(j) = lagrange*nTm(j)
    enddo
    j = (nnode+1)*3 + 1
    ctNForce(j) = dot_product(nTm(1:(nnode+1)*3),ndCoord(1:(nnode+1)*3))

    if(fcoeff /= 0.0d0 .and. lagrange > 0.0d0)then

      if( cflag ) then !calc tf_final and set it to tangentForce_final
        ctstate%reldisp(1:3) = matmul(Tt(1:3,1:3*(nnode+1)), ndDu(1:3*(nnode+1)))

        call getTrialFricForceANDcheckFricState(ctstate,fcoeff,tPenalty,lagrange)

        if( ctstate%state == contactStick ) then
          tf_final(1:3) = ctstate%tangentForce_trial(1:3)
        elseif( ctstate%state == contactSlip ) then
          tf_trial(1:3) = ctstate%tangentForce_trial(1:3)
          length_tft = dsqrt(dot_product(tf_trial,tf_trial))
          tangent(1:3) = tf_trial(1:3)/length_tft
          tf_final(1:3) = fcoeff*dabs(lagrange)*tangent(1:3)
        endif
        ctstate%tangentForce_final(1:3) = tf_final(1:3)
      else ! just set tangentForce_final to tf_final (used for fstr_ass_load_contact)
        tf_final(1:3) = ctstate%tangentForce_final(1:3)
      endif

      ! Distribute friction force using Tm: ctTForce = -Tm^T * tf_final
      ctTForce(1:(nnode+1)*3) = -matmul(transpose(Tm(1:3, 1:3*(nnode+1))), tf_final)

    endif

  end subroutine getContactNodalForce_Slag

  subroutine getTiedStiffness_Slag(ctState,tSurf,idof,stiffness)
    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    integer(kind=kint)  :: idof
    real(kind=kreal)    :: stiffness(:,:) !< contact stiffness matrix

    integer(kind=kint)  :: i, j
    real(kind=kreal)    :: nTm((l_max_surface_node+1)*3)
    real(kind=kreal)    :: Tm(3, 3*(l_max_surface_node+1)), Tt(3, 3*(l_max_surface_node+1)) !< mapping matrices
    real(kind=kreal)    :: e_idof(3)  !< unit vector in idof direction

    nnode = size(tSurf%nodes)
    stiffness = 0.0d0

    ! Compute Tm matrix
    call computeTm_Tt(ctState, tSurf, 0.0d0, Tm, Tt)

    ! Create unit vector in idof direction
    e_idof = 0.0d0
    e_idof(idof) = 1.0d0

    ! Compute nTm = Tm^T * e_idof
    ! This gives constraint direction for tied contact
    nTm(1:(nnode+1)*3) = matmul(transpose(Tm(1:3, 1:3*(nnode+1))), e_idof)

    i = (nnode+1)*3 + 1
    do j = 1, (nnode+1)*3
      stiffness(i,j) = nTm(j);  stiffness(j,i) = nTm(j)
    enddo

  end subroutine getTiedStiffness_Slag

  subroutine getTiedNodalForce_Slag(ctState,tSurf,idof,ndu,lagrange,ctNForce,ctTForce)
    use mSurfElement
    type(tContactState) :: ctState !< type tContactState
    type(tSurfElement)  :: tSurf !< surface element structure
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    integer(kind=kint)  :: idof
    real(kind=kreal)   :: ndu(:) !< nodal displacement
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)       :: ctNForce(:)  !< tied contact force vector
    real(kind=kreal)       :: ctTForce(:)  !< tied contact force vector

    integer(kind=kint) :: j
    real(kind=kreal)   :: normal(3) !< normal vector at target point
    real(kind=kreal)   :: nTm((l_max_surface_node + 1)*3) !< normal vector
    real(kind=kreal)   :: tTm((l_max_surface_node + 1)*3) !< tangential vector
    real(kind=kreal)   :: Tm(3, 3*(l_max_surface_node+1)), Tt(3, 3*(l_max_surface_node+1)) !< mapping matrices
    real(kind=kreal)   :: e_idof(3), e_normal(3), e_tangent(3)  !< unit vectors

    nnode = size(tSurf%nodes)

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Compute Tm matrix
    call computeTm_Tt(ctState, tSurf, 0.0d0, Tm, Tt)

    normal(1:3) = ctState%direction(1:3)

    ! Create unit vector in idof direction
    e_idof = 0.0d0
    e_idof(idof) = 1.0d0

    ! Normal component: -normal(idof) * normal (OLD: nTm(1:3) = -normal(idof)*normal(1:3))
    e_normal = -normal(idof) * normal

    ! Tangential component: -e_idof - e_normal (OLD: tTm(idof)=-1.0; tTm(1:3)-=nTm(1:3))
    e_tangent = -e_idof - e_normal

    ! Compute nTm and tTm using Tm^T
    nTm(1:(nnode+1)*3) = matmul(transpose(Tm(1:3, 1:3*(nnode+1))), e_normal)
    tTm(1:(nnode+1)*3) = matmul(transpose(Tm(1:3, 1:3*(nnode+1))), e_tangent)

    do j = 1, (nnode+1)*3
      ctNForce(j) = lagrange*nTm(j)
      ctTForce(j) = lagrange*tTm(j)
    enddo
    j = (nnode+1)*3 + 1
    ctNForce(j) = dot_product(nTm(1:(nnode+1)*3),ndu(1:(nnode+1)*3))
    ctTForce(j) = dot_product(tTm(1:(nnode+1)*3),ndu(1:(nnode+1)*3))
  end subroutine getTiedNodalForce_Slag

end module m_fstr_contact_elem_slag
