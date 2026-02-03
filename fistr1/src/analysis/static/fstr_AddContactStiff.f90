!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions:
!!  1) obtain contact stiffness matrix of each contact pair and assemble
!!     it into global stiffness matrix.
!!  2) obtain contact nodal force vector of each contact pair and assemble
!!      it into right-hand side vector to update non-equilibrated nodal force vector.
!!  3) Modify Lagrange multiplier-related part of stiffness matrix and right-hand side
!!     vector for dealing with prescribed displacement boundary condition.

module m_addContactStiffness

  use m_fstr
  use elementInfo
  use mContact
  use m_contact_lib
  use fstr_matrix_con_contact
  use hecmw_matrix_ass
  use m_fstr_Residual

  implicit none

  private

  public :: fstr_AddContactStiffness
  public :: fstr_Update_NDForce_contact

contains

  !> \brief This subroutine obtains contact stiffness matrix of each contact pair
  !!  and assembles it into global stiffness matrix.
  subroutine fstr_AddContactStiffness(cstep,ctAlgo, iter,hecMESH,hecMAT,hecLagMAT,fstrSOLID)

    integer(kind=kint)                   :: cstep !< current loading step
    integer(kind=kint)                   :: ctAlgo !< contact algorithm type
    integer(kind=kint)                   :: iter
    type(hecmwST_local_mesh)             :: hecMESH !< type hecmwST_local_mesh
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint)                   :: i, grpid, id_lagrange

    if( associated(hecLagMAT%AL_lagrange) ) hecLagMAT%AL_lagrange = 0.0d0
    if( associated(hecLagMAT%AU_lagrange) ) hecLagMAT%AU_lagrange = 0.0d0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_stiffness_NodeSurf( ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), iter, hecLagMAT%Lagrange(:), &
        hecMAT, hecLagMAT )

    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_stiffness_NodeSurf( ctAlgo, fstrSOLID%embeds(i), hecMESH%node(:), iter, hecLagMAT%Lagrange(:), &
        hecMAT, hecLagMAT )

    enddo

  end subroutine fstr_AddContactStiffness

  !>\brief This subroutine calculates contact stiffness for each contact pair
  !! and assembles it into global stiffness matrix
  subroutine calcu_contact_stiffness_NodeSurf( ctAlgo, contact, coord, iter, lagrange_array, &
      hecMAT, hecLagMAT)
    use mContact, only: mu, mut
    integer(kind=kint), intent(in)             :: ctAlgo          !< contact analysis algorithm
    type(tContact), intent(inout)              :: contact         !< contact info
    real(kind=kreal), intent(in)               :: coord(:)        !< mesh coordinate
    integer(kind=kint), intent(in)             :: iter            !< iteration number
    real(kind=kreal), intent(in)               :: lagrange_array(:) !< Lagrange multiplier array
    type(hecmwST_matrix), intent(inout)        :: hecMAT          !< global stiffness matrix
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT     !< Lagrange matrix

    integer(kind=kint), parameter :: NDOF = 3
    integer(kind=kint) :: ctsurf, nnode, ndLocal(21), etype
    integer(kind=kint) :: j, k, m, nd, algtype, id_lagrange
    real(kind=kreal)   :: lagrange
    real(kind=kreal)   :: stiffness((l_max_surface_node+1)*3+1, (l_max_surface_node+1)*3+1)
    real(kind=kreal)   :: elecoord(3, l_max_surface_node)  !< master node coordinates
    real(kind=kreal)   :: force(l_max_surface_node*3+3)    !< contact force direction
    real(kind=kreal)   :: nrlforce                         !< normal force magnitude

    algtype = contact%algtype

    do j = 1, size(contact%slave)

      if( contact%states(j)%state == CONTACTFREE ) cycle

      ctsurf = contact%states(j)%surface
      etype = contact%master(ctsurf)%etype
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)

      ! Prepare master node coordinates for ALagrange
      do k = 1, nnode
        elecoord(1:3, k) = coord(3*ndLocal(k+1)-2:3*ndLocal(k+1))
      enddo

      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then

        if( ctAlgo == kcaSLagrange ) then
          id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
          id_lagrange = id_lagrange + 1
          lagrange = lagrange_array(id_lagrange)
          call getContactStiffness_Slag(contact%states(j), contact%master(ctsurf), iter, &
            contact%tPenalty, contact%fcoeff, lagrange, stiffness)

          ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
          call hecmw_mat_ass_contactlag(nnode, ndLocal, id_lagrange, contact%fcoeff, stiffness, hecMAT, hecLagMAT)

        else if( ctAlgo == kcaALagrange ) then
          call getContactStiffness_Alag(contact%states(j), etype, nnode, elecoord(:,1:nnode), &
            contact%fcoeff, contact%symmetric, stiffness, force)

          ! Assemble contact stiffness matrix into global stiffness matrix
          call hecmw_mat_ass_elem(hecMAT, nnode+1, ndLocal, stiffness)

        end if


      else if( algtype == CONTACTTIED ) then

        if( ctAlgo == kcaSLagrange ) then
          id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
          do k = 1, 3
            id_lagrange = id_lagrange + 1
            lagrange = lagrange_array(id_lagrange)

            call getTiedStiffness_Slag(contact%states(j), contact%master(ctsurf), k, stiffness)
            ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
            call hecmw_mat_ass_contactlag(nnode, ndLocal, id_lagrange, 0.d0, stiffness, hecMAT, hecLagMAT)
          enddo

        else if( ctAlgo == kcaALagrange ) then
          call getTiedStiffness_Alag(contact%states(j), etype, nnode, stiffness, force)

          ! Assemble contact stiffness matrix into global stiffness matrix
          call hecmw_mat_ass_elem(hecMAT, nnode+1, ndLocal, stiffness)

        end if

      endif

      ! Initial contact: add enforced displacement constraint (ALagrange only, iter==1)
      if( ctAlgo == kcaALagrange .and. iter == 1 ) then
        contact%states(j)%wkdist = contact%states(j)%distance
        nrlforce = -mu * contact%states(j)%distance
        force(1:nnode*NDOF+NDOF) = force(1:nnode*NDOF+NDOF) * nrlforce
        do m = 1, nnode+1
          nd = ndLocal(m)
          do k = 1, NDOF
            hecMAT%B(NDOF*(nd-1)+k) = hecMAT%B(NDOF*(nd-1)+k) - force((m-1)*NDOF+k)
          enddo
        enddo
      endif
      

    enddo

  end subroutine calcu_contact_stiffness_NodeSurf

  !> \brief This subroutine obtains contact stiffness matrix for ALagrange method
  !! Equivalent to contact2stiff but outputs only stiffness (force computed separately)
  subroutine getContactStiffness_Alag(cstate, etype, nnode, ele, fcoeff, symm, stiff, force)
    use mContact, only: mu, mut

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

  !> \brief This subroutine obtains tied contact stiffness matrix for ALagrange method
  !! Equivalent to tied2stiff but outputs only stiffness
  subroutine getTiedStiffness_Alag(cstate, etype, nnode, stiff, force)
    use mContact, only: mu

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

  !> \brief This subroutine obtains contact stiffness matrix of contact pair
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

  !> \brief This subroutine obtains contact stiffness matrix of contact pair
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

  !> \brief This subroutine obtains contact nodal force vector of each contact pair
  !! and assembles it into right-hand side vector to update non-equilibrated nodal force vector.
  subroutine fstr_Update_NDForce_contact(cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)

    integer(kind=kint), intent(in)       :: cstep   !< current calculation step
    integer(kind=kint), intent(in)       :: ctAlgo  !< contact analysis algorithm
    type(hecmwST_local_mesh)             :: hecMESH !< type hecmwST_local_mesh
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange) :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(hecmwST_matrix)                 :: conMAT !< type hecmwST_matrix for contact part only
    integer(kind=kint) :: ctsurf, nnode, ndLocal(21) !< contents of type tContact
    integer(kind=kint) :: i, j, k, nlag, id_lagrange, algtype
    real(kind=kreal)   :: ndCoord(21*3) !< nodal coordinates
    real(kind=kreal)   :: ndu(21*3), ndDu(21*3) !< nodal displacement and its increment
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)   :: ctNForce(21*3+1)     !< nodal normal contact force vector
    real(kind=kreal)   :: ctTForce(21*3+1)     !< nodal tangential contact force vector

    integer(kind=kint) :: grpid
    logical            :: if_flag
    real(kind=kreal)   :: ctime, etime
    integer(kind=kint) :: if_type

    if( associated(fstrSOLID%CONT_NFORCE) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
    if( associated(fstrSOLID%CONT_FRIC) ) fstrSOLID%CONT_FRIC(:) = 0.d0
    if( associated(fstrSOLID%EMBED_NFORCE) ) fstrSOLID%EMBED_NFORCE(:) = 0.d0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_ndforce_NodeSurf( ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:), &
        fstrSOLID%dunode(:), hecLagMAT%Lagrange(:), conMAT, &
        fstrSOLID%CONT_NFORCE, fstrSOLID%CONT_FRIC, hecLagMAT )

    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_ndforce_NodeSurf( ctAlgo, fstrSOLID%embeds(i), hecMESH%node(:), fstrSOLID%unode(:), &
        fstrSOLID%dunode(:), hecLagMAT%Lagrange(:), conMAT, &
        fstrSOLID%EMBED_NFORCE, fstrSOLID%EMBED_NFORCE, hecLagMAT )

    enddo

    !    Consider SPC condition
    call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, hecMAT%B)
    call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, conMAT%B)

  end subroutine fstr_Update_NDForce_contact

  !>\brief This subroutine calculates contact nodal force for each contact pair
  !! and assembles it into contact matrix and force arrays
  subroutine calcu_contact_ndforce_NodeSurf( ctAlgo, contact, coord, disp, ddisp, lagrange_array, &
      conMAT, CONT_NFORCE, CONT_FRIC, hecLagMAT )
    integer(kind=kint), intent(in)       :: ctAlgo          !< contact analysis algorithm
    type( tContact ), intent(inout)      :: contact         !< contact info
    real(kind=kreal), intent(in)         :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)         !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)        !< disp till current substep
    real(kind=kreal), intent(in)         :: lagrange_array(:) !< Lagrange multiplier array
    type(hecmwST_matrix), intent(inout)  :: conMAT          !< contact matrix
    real(kind=kreal), pointer            :: CONT_NFORCE(:)  !< contact normal force
    real(kind=kreal), pointer            :: CONT_FRIC(:)    !< contact friction force
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT  !< Lagrange matrix

    integer(kind=kint) :: ctsurf, nnode, ndLocal(21)
    integer(kind=kint) :: j, k, algtype, id_lagrange
    real(kind=kreal)   :: ndCoord(21*3)
    real(kind=kreal)   :: ndu(21*3), ndDu(21*3)
    real(kind=kreal)   :: lagrange
    real(kind=kreal)   :: ctNForce(21*3+1)
    real(kind=kreal)   :: ctTForce(21*3+1)
    logical            :: if_flag
    real(kind=kreal)   :: ctime, etime
    integer(kind=kint) :: if_type

    algtype = contact%algtype
    if_flag = (contact%if_type /= 0)
    if(if_flag)then
      ctime = contact%ctime
      etime = contact%if_etime
      if_type = contact%if_type
    end if

    do j = 1, size(contact%slave)

      if( contact%states(j)%state == CONTACTFREE ) cycle
      if(if_flag) call set_shrink_factor(ctime, contact%states(j), etime, if_type)

      ctsurf = contact%states(j)%surface
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(j)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)
      do k = 1, nnode+1
        ndDu((k-1)*3+1:(k-1)*3+3) = ddisp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
        ndu((k-1)*3+1:(k-1)*3+3) = disp((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
        ndCoord((k-1)*3+1:(k-1)*3+3) = coord((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
      enddo

      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        ! Obtain contact nodal force vector of contact pair
        if(if_flag) call get_shrink_elemact_surf(contact%states(j),ndCoord, nnode)

        if( ctAlgo == kcaSLagrange ) then
          id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
          id_lagrange = id_lagrange + 1
          lagrange = lagrange_array(id_lagrange)
          call getContactNodalForce_Slag(contact%states(j),contact%master(ctsurf),ndCoord,ndDu,    &
          contact%tPenalty,contact%fcoeff,lagrange,ctNForce,ctTForce,.true.)

        else if( ctAlgo == kcaALagrange ) then
          id_lagrange = 0
          lagrange = 0.d0
          call getContactNodalForce_Alag(contact%states(j),contact%master(ctsurf),ndCoord,ndDu,    &
          contact%tPenalty,contact%fcoeff,lagrange,ctNForce,ctTForce,.true.)

        end if

        ! Update non-eqilibrited force vector
        call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce,  &
          &  conMAT,CONT_NFORCE,CONT_FRIC)

      else if( algtype == CONTACTTIED ) then

        if( ctAlgo == kcaSLagrange ) then
          id_lagrange = hecLagMAT%lag_node_table(ndLocal(1)) - 1
          do k=1,3
            id_lagrange = id_lagrange + 1
            lagrange = lagrange_array(id_lagrange)
            contact%states(j)%multiplier(k) = lagrange

            call getTiedNodalForce_Slag(contact%states(j),contact%master(ctsurf),k,ndu, & 
              &  lagrange,ctNForce,ctTForce)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,1.d0,ctNForce,ctTForce,  &
              &  conMAT,CONT_NFORCE)
          end do

        else if( ctAlgo == kcaALagrange ) then
          id_lagrange = 0
          lagrange = 0.d0
          call getTiedNodalForce_Alag(contact%states(j),contact%master(ctsurf),ndu,    &
          lagrange,ctNForce,ctTForce)
          ! Update non-eqilibrited force vector
          call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce,  &
            &  conMAT,CONT_NFORCE)

        end if

      endif

    enddo

  end subroutine calcu_contact_ndforce_NodeSurf

  !> \brief This subroutine obtains contact nodal force vector of contact pair
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

  !> \brief This subroutine obtains contact nodal force vector of contact pair
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

  !> \brief This subroutine obtains contact nodal force vector of contact pair for ALagrange
  subroutine getContactNodalForce_Alag(ctState,tSurf,ndCoord,ndDu,tPenalty,fcoeff,lagrange,ctNForce,ctTForce,cflag)

    use mSurfElement
    use mContact, only: mu, mut
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


  !> \brief This subroutine obtains tied contact nodal force vector for ALagrange
  subroutine getTiedNodalForce_Alag(ctState,tSurf,ndu,lagrange,ctNForce,ctTForce)

    use mSurfElement
    use mContact, only: mu
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

  !> \brief This subroutine assembles contact nodal force vector into right-hand side vector
  !! to update non-equilibrated nodal force vector.
  subroutine update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce,hecMAT,cont_nforce,cont_fric)

    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    integer(kind=kint) :: nnode, ndLocal(nnode + 1) !< number of nodes of master segment
    !< global number of nodes of contact pair
    integer(kind=kint) :: id_lagrange !< number of Lagrange multiplier
    real(kind=kreal)                        :: lagrange                        !< value of Lagrange multiplier
    integer(kind=kint) :: np, ndof !< total number of nodes; degree of freedom
    integer (kind=kint)                     :: i, inod, idx
    real(kind=kreal)                        :: ctNForce((nnode+1)*3+1)         !< contact force vector
    real(kind=kreal)                        :: ctTForce((nnode+1)*3+1)         !< contact force vector
    real(kind=kreal), pointer               :: cont_nforce(:)         !< contact force vector
    real(kind=kreal), pointer, optional     :: cont_fric(:)         !< contact force vector

    np = hecMAT%NP; ndof = hecMAT%NDOF

    do i = 1, nnode + 1
      inod = ndLocal(i)
      idx = (inod-1)*3+1
      hecMAT%B(idx:idx+2) = hecMAT%B(idx:idx+2) + ctNForce((i-1)*3+1:(i-1)*3+3) + ctTForce((i-1)*3+1:(i-1)*3+3)
      cont_nforce(idx:idx+2) = cont_nforce(idx:idx+2) + ctNForce((i-1)*3+1:(i-1)*3+3)
      if( present(cont_fric) ) cont_fric(idx:idx+2) = cont_fric(idx:idx+2) + ctTForce((i-1)*3+1:(i-1)*3+3)
    enddo

    if( id_lagrange > 0 ) hecMAT%B(np*ndof+id_lagrange) = ctNForce((nnode+1)*3+1)+ctTForce((nnode+1)*3+1)

  end subroutine update_NDForce_contact

  subroutine get_shrink_elemact_surf(cstate, coords, nnode)
    use m_fstr_TimeInc
    integer(kind=kint)              :: nnode, i
    type(tContactState)             :: cstate 
    real(kind=kreal)                :: coords(:)
    real(kind=kreal)   :: d(3)

    d = - cstate%shrink_factor * cstate%direction(1:3)
    
    do i = 1, nnode
      coords(1+i*3:(i+1)*3) = coords(1+i*3:(i+1)*3) + d(1:3)
    enddo

  end subroutine get_shrink_elemact_surf


end module m_addContactStiffness



