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

  public :: fstr_AddContactStiffness
  public :: fstr_Update_NDForce_contact
  public :: update_NDForce_contact
  public :: fstr_ass_load_contact

  private :: computeTm_Tt
  private :: getContactStiffness
  private :: getContactNodalForce
  private :: getTrialFricForceANDcheckFricState

contains

  !> \brief This subroutine obtains contact stiffness matrix of each contact pair
  !!  and assembles it into global stiffness matrix.
  subroutine fstr_AddContactStiffness(cstep,iter,hecMAT,hecLagMAT,fstrSOLID)

    integer(kind=kint)                   :: cstep !< current loading step
    integer(kind=kint)                   :: iter
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint)                   :: ctsurf, nnode, ndLocal(21) !< contents of type tContact
    integer(kind=kint)                   :: i, j, k, nlag, id_lagrange, grpid, algtype
    real(kind=kreal)                     :: lagrange
    real(kind=kreal)                     :: stiffness(21*3 + 1, 21*3 + 1)

    hecLagMAT%AL_lagrange = 0.0d0
    hecLagMAT%AU_lagrange = 0.0d0

    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,nlag
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
            ! Obtain contact stiffness matrix of contact pair
            call getContactStiffness(fstrSOLID%contacts(i)%states(j),fstrSOLID%contacts(i)%master(ctsurf),iter,  &
              fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,stiffness)
          else if( algtype == CONTACTTIED ) then
            call getTiedStiffness(fstrSOLID%contacts(i)%states(j),fstrSOLID%contacts(i)%master(ctsurf),k,stiffness)
          endif 
    
          ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
          call hecmw_mat_ass_contactlag(nnode,ndLocal,id_lagrange,fstrSOLID%contacts(i)%fcoeff,stiffness,hecMAT,hecLagMAT)
        enddo
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      do j = 1, size(fstrSOLID%embeds(i)%slave)

        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%embeds(i)%states(j)%surface
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,3
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          call getTiedStiffness(fstrSOLID%embeds(i)%states(j),fstrSOLID%embeds(i)%master(ctsurf),k,stiffness)
    
          ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
          call hecmw_mat_ass_contactlag(nnode,ndLocal,id_lagrange,0.d0,stiffness,hecMAT,hecLagMAT)
        enddo
      enddo
    enddo

  end subroutine fstr_AddContactStiffness

  !> \brief This subroutine obtains contact stiffness matrix of contact pair
  subroutine getTiedStiffness(ctState,tSurf,idof,stiffness)
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

  end subroutine getTiedStiffness

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
  subroutine getContactStiffness(ctState,tSurf,iter,tPenalty,fcoeff,lagrange,stiffness)

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

  end subroutine getContactStiffness

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

    id_lagrange = 0
    if( associated(fstrSOLID%CONT_NFORCE) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
    if( associated(fstrSOLID%CONT_FRIC) ) fstrSOLID%CONT_FRIC(:) = 0.d0
    if( associated(fstrSOLID%EMBED_NFORCE) ) fstrSOLID%EMBED_NFORCE(:) = 0.d0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_ndforce_NodeSurf( ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:), &
        fstrSOLID%dunode(:), hecLagMAT%Lagrange(:), id_lagrange, conMAT, &
        fstrSOLID%CONT_NFORCE, fstrSOLID%CONT_FRIC )

    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      do j = 1, size(fstrSOLID%embeds(i)%slave)

        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%embeds(i)%states(j)%surface
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(1:nnode)
        do k = 1, nnode+1
          ndDu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%dunode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          ndu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
          ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
        enddo

        do k=1,3
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          fstrSOLID%embeds(i)%states(j)%multiplier(k) = hecLagMAT%Lagrange(id_lagrange)
    
          call getTiedNodalForce(fstrSOLID%embeds(i)%states(j),fstrSOLID%embeds(i)%master(ctsurf),k,ndu,lagrange,ctNForce,ctTForce)
          ! Update non-eqilibrited force vector
          call update_NDForce_contact(nnode,ndLocal,id_lagrange,-1.d0,ctNForce,ctTForce,  &
          &  conMAT,fstrSOLID%EMBED_NFORCE,fstrSOLID%EMBED_NFORCE)

        enddo
      enddo
    enddo

    !    Consider SPC condition
    call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, hecMAT%B)
    call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, conMAT%B)

  end subroutine fstr_Update_NDForce_contact

  !>\brief This subroutine calculates contact nodal force for each contact pair
  !! and assembles it into contact matrix and force arrays
  subroutine calcu_contact_ndforce_NodeSurf( ctAlgo, contact, coord, disp, ddisp, lagrange_array, &
      id_lagrange, conMAT, CONT_NFORCE, CONT_FRIC )
    integer(kind=kint), intent(in)       :: ctAlgo          !< contact analysis algorithm
    type( tContact ), intent(inout)      :: contact         !< contact info
    real(kind=kreal), intent(in)         :: coord(:)        !< mesh coordinate
    real(kind=kreal), intent(in)         :: disp(:)         !< disp till current step
    real(kind=kreal), intent(in)         :: ddisp(:)        !< disp till current substep
    real(kind=kreal), intent(in)         :: lagrange_array(:) !< Lagrange multiplier array
    integer(kind=kint), intent(inout)    :: id_lagrange     !< Lagrange multiplier index
    type(hecmwST_matrix), intent(inout)  :: conMAT          !< contact matrix
    real(kind=kreal), pointer            :: CONT_NFORCE(:)  !< contact normal force
    real(kind=kreal), pointer            :: CONT_FRIC(:)    !< contact friction force

    integer(kind=kint) :: ctsurf, nnode, ndLocal(21)
    integer(kind=kint) :: j, k, algtype
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

    if( ctAlgo == kcaALagrange) then
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call calcu_contact_force0( contact, coord(:), disp(:)  &
        , ddisp(:), contact%fcoeff, mu, mut, conMAT%B )
      else if( algtype == CONTACTTIED ) then
        call calcu_tied_force0( contact, disp(:), ddisp(:), mu, conMAT%B )
      endif
      return
    endif

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
          id_lagrange = id_lagrange + 1
          lagrange = lagrange_array(id_lagrange)
          contact%states(j)%multiplier(1) = lagrange

          call getContactNodalForce(contact%states(j),contact%master(ctsurf),ndCoord,ndDu,    &
          contact%tPenalty,contact%fcoeff,lagrange,ctNForce,ctTForce,.true.)
          ! Update non-eqilibrited force vector
          call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce,  &
            &  conMAT,CONT_NFORCE,CONT_FRIC)
          
        else if( ctAlgo == kcaALagrange ) then
          id_lagrange = 0
          lagrange = 0.d0
        end if

      else if( algtype == CONTACTTIED ) then

        if( ctAlgo == kcaSLagrange ) then
          do k=1,3
            id_lagrange = id_lagrange + 1
            lagrange = lagrange_array(id_lagrange)
            contact%states(j)%multiplier(k) = lagrange

            call getTiedNodalForce(contact%states(j),contact%master(ctsurf),k,ndu, & 
              &  lagrange,ctNForce,ctTForce)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,1.d0,ctNForce,ctTForce,  &
              &  conMAT,CONT_NFORCE,CONT_FRIC)
          end do

        else if( ctAlgo == kcaALagrange ) then
          id_lagrange = 0
          lagrange = 0.d0
        end if

      endif

    enddo

  end subroutine calcu_contact_ndforce_NodeSurf

  !>\brief This subroutine updates contact condition as follows:
  !!-# Contact force from multiplier and disp increment
  !!-# Update nodal force residual
  subroutine calcu_contact_force0( contact, coord, disp, ddisp, fcoeff, mu,     &
      mut, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: fcoeff         !< frictional coeff
    real(kind=kreal), intent(in)      :: mu, mut        !< penalty
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: nrlforce, elemdisp(3,l_max_elem_node), tangent(3,2)
    real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(3*l_max_elem_node+3), edisp(3*l_max_elem_node+3)

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemdisp(:,j) = ddisp(3*iSS-2:3*iSS)
        edisp(3*j+1:3*j+3) = ddisp(3*iSS-2:3*iSS)
        elemcrd(:,j) = coord(3*iSS-2:3*iSS)+disp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      elemg = 0.d0
      do j=1,nn
        elemg(:) = elemg(:)+shapefunc(j)*(elemcrd(:,j)+elemdisp(:,j))
      enddo
      dg(:) = coord(3*slave-2:3*slave)+disp(3*slave-2:3*slave)+ddisp(3*slave-2:3*slave) -  elemg(:)
      dgn = dot_product( contact%states(i)%direction, dg )
      nrlforce = contact%states(i)%multiplier(1) + mu*dgn
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(i)%direction

      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
          shapefunc(j)*contact%states(i)%direction
      enddo

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      dxi(1) = dot_product( dispmat(1,1:nn*3+3), edisp(1:nn*3+3) )
      dxi(2) = dot_product( dispmat(2,1:nn*3+3), edisp(1:nn*3+3) )
      dxy(:) = matmul( metric, dxi )
      fric(1:2) = contact%states(i)%multiplier(2:3) + mut*dxy(1:2)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)

      if(  contact%states(i)%state==CONTACTSLIP ) then
        dgn = dsqrt( f3(1)*f3(1)+f3(2)*f3(2)+f3(3)*f3(3) )
        f3(:) = f3(:)*fcoeff*contact%states(i)%multiplier(1)/dgn
      endif
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
      enddo
    enddo
  end subroutine calcu_contact_force0

  !>\brief This subroutine updates contact condition as follows:
  !!-# Contact force from multiplier and disp increment
  !!-# Update nodal force residual
  subroutine calcu_tied_force0( contact, disp, ddisp, mu, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: mu             !< penalty
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: nrlforce(3)
    real(kind=kreal)    :: dg(3)
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3)

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = disp(3*slave-2:3*slave)+ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        edisp(3*j+1:3*j+3) = disp(3*iSS-2:3*iSS)+ddisp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      dg(1:3) = edisp(1:3)
      do j=1,nn
        dg(1:3) = dg(1:3)-shapefunc(j)*edisp(3*j+1:3*j+3)
      enddo

      nrlforce(1:3) = contact%states(i)%multiplier(1:3)+mu*dg(1:3)

      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS) + shapefunc(j)*nrlforce(1:3)
      enddo
    enddo

  end subroutine 

  !> \brief This subroutine obtains contact nodal force vector of contact pair
  subroutine getTiedNodalForce(ctState,tSurf,idof,ndu,lagrange,ctNForce,ctTForce)
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
  end subroutine 

  !> \brief This subroutine obtains contact nodal force vector of contact pair
  subroutine getContactNodalForce(ctState,tSurf,ndCoord,ndDu,tPenalty,fcoeff,lagrange,ctNForce,ctTForce,cflag)

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

  end subroutine getContactNodalForce


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


  !> \brief This subroutine adds initial contact force vector to the right-hand side vector
  !> \at the beginning of each substep calculation
  subroutine fstr_ass_load_contact(cstep, hecMESH, hecMAT, fstrSOLID, hecLagMAT)

    type(hecmwST_local_mesh)             :: hecMESH !< type hecmwST_local_mesh
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint) :: cstep !< current step
    integer(kind=kint) :: np, ndof !< total number of nodes; degree of freedom
    integer(kind=kint) :: i, j, k, id_lagrange, grpid
    integer(kind=kint) :: ctsurf, nnode, ndLocal(21) !< contents of type tContact
    real(kind=kreal)   :: ndu(21*3), ndCoord(21*3), lagrange !< nodal coordinates; value of Lagrange mutiplier

    integer(kind=kint) :: algtype, nlag    
    real(kind=kreal)   :: ctNForce(21*3+1)     !< nodal normal contact force vector
    real(kind=kreal)   :: ctTForce(21*3+1)     !< nodal tangential contact force vector
    logical            :: if_flag
    real(kind=kreal)   :: ctime, etime
    integer(kind=kint) :: if_type

    np = hecMAT%NP ; ndof = hecMAT%NDOF

    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      if_flag = (fstrSOLID%contacts(i)%if_type /= 0)
      if(if_flag)then
        ctime = fstrSOLID%contacts(i)%ctime
        etime = fstrSOLID%contacts(i)%if_etime
        if_type = fstrSOLID%contacts(i)%if_type
      end if

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle
        if(if_flag) call set_shrink_factor(ctime, fstrSOLID%contacts(i)%states(j), etime, if_type)

        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)
        do k = 1, nnode+1
          ndu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
        enddo

        do k=1,nlag
          id_lagrange = id_lagrange + 1
          lagrange = fstrSOLID%contacts(i)%states(j)%multiplier(k)
          if(if_flag) call get_shrink_elemact_surf(fstrSOLID%contacts(i)%states(j),ndCoord, nnode)
          if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
            ! Obtain contact nodal force vector of contact pair
            call getContactNodalForce(fstrSOLID%contacts(i)%states(j),fstrSOLID%contacts(i)%master(ctsurf),ndCoord,ndu,    &
            fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,ctNForce,ctTForce,.false.)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce, & 
               &  hecMAT,fstrSOLID%CONT_NFORCE,fstrSOLID%CONT_FRIC)
          else if( algtype == CONTACTTIED ) then
            call getTiedNodalForce(fstrSOLID%contacts(i)%states(j),fstrSOLID%contacts(i)%master(ctsurf),k,ndu,&
              &  lagrange,ctNForce,ctTForce)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,-1.d0,ctNForce,ctTForce,hecMAT,fstrSOLID%CONT_NFORCE)
          endif 

        enddo
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle


      do j = 1, size(fstrSOLID%embeds(i)%slave)

        if( fstrSOLID%embeds(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%embeds(i)%states(j)%surface
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(1:nnode)
        do k = 1, nnode+1
          ndu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
        enddo

        do k=1,3
          id_lagrange = id_lagrange + 1
          lagrange = fstrSOLID%embeds(i)%states(j)%multiplier(k)
    
          call getTiedNodalForce(fstrSOLID%embeds(i)%states(j),fstrSOLID%embeds(i)%master(ctsurf),k,ndu,lagrange,ctNForce,ctTForce)
          ! Update non-eqilibrited force vector
          call update_NDForce_contact(nnode,ndLocal,id_lagrange,-1.d0,ctNForce,ctTForce,hecMAT,fstrSOLID%EMBED_NFORCE)

        enddo
      enddo
    enddo

  end subroutine fstr_ass_load_contact

  !> Update lagrangian multiplier
  subroutine fstr_ass_load_contactAlag( hecMESH, fstrSOLID, B )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH
    type(fstr_solid), intent(inout)        :: fstrSOLID
    real(kind=kreal), intent(inout)        :: B(:)        !< nodal force residual

    integer(kind=kint) :: i, algtype

    do i = 1, fstrSOLID%n_contacts
      algtype = fstrSOLID%contacts(i)%algtype
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call ass_contact_force( fstrSOLID%contacts(i), hecMESH%node, fstrSOLID%unode, B )
      else if( algtype == CONTACTTIED ) then
        call calcu_tied_force0( fstrSOLID%contacts(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
      endif
    enddo

    do i = 1, fstrSOLID%n_embeds
      call calcu_tied_force0( fstrSOLID%embeds(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
    enddo
  end subroutine

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



