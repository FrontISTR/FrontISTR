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
  use m_contact_lib
  use fstr_matrix_con_contact
  use hecmw_matrix_ass
  use m_fstr_Residual

  implicit none

  public :: fstr_AddContactStiffness
  public :: fstr_Update_NDForce_contact
  public :: update_NDForce_contact
  public :: fstr_ass_load_contact

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
    integer(kind=kint)                   :: ctsurf, etype, nnode, ndLocal(21) !< contents of type tContact
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
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,nlag
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
            ! Obtain contact stiffness matrix of contact pair
            call getContactStiffness(iter,etype,nnode,fstrSOLID%contacts(i)%states(j),  &
              fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,stiffness)
          else if( algtype == CONTACTTIED ) then
            call getTiedStiffness(etype,nnode,k,fstrSOLID%contacts(i)%states(j),lagrange,stiffness)
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
        etype = fstrSOLID%embeds(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(1:nnode)

        do k=1,3
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          call getTiedStiffness(etype,nnode,k,fstrSOLID%embeds(i)%states(j),lagrange,stiffness)
    
          ! Assemble contact stiffness matrix of contact pair into global stiffness matrix
          call hecmw_mat_ass_contactlag(nnode,ndLocal,id_lagrange,0.d0,stiffness,hecMAT,hecLagMAT)
        enddo
      enddo
    enddo

  end subroutine fstr_AddContactStiffness

  !> \brief This subroutine obtains contact stiffness matrix of contact pair
  subroutine getTiedStiffness(etype,nnode,idof,ctState,lagrange,stiffness)
    integer(kind=kint)  :: etype !< type of master segment
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    integer(kind=kint)  :: idof
    type(tContactState) :: ctState !< type tContactState
    real(kind=kreal)    :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)    :: stiffness(:,:) !< contact stiffness matrix

    integer(kind=kint)  :: i, j, k, l
    real(kind=kreal)    :: shapefunc(nnode) !< normal vector at target point; shape functions
    real(kind=kreal)    :: nTm((nnode+1)*3)

    stiffness = 0.0d0

    call getShapeFunc( etype, ctState%lpos(:), shapefunc )

    nTm = 0.d0
    nTm(idof) = 1.d0
    do i = 1, nnode
      nTm(i*3+idof) = -shapefunc(i)
    enddo

    i = (nnode+1)*3 + 1
    do j = 1, (nnode+1)*3
      stiffness(i,j) = nTm(j);  stiffness(j,i) = nTm(j)
    enddo

  end subroutine getTiedStiffness

  !> \brief This subroutine obtains contact stiffness matrix of contact pair
  subroutine getContactStiffness(iter,etype,nnode,ctState,tPenalty,fcoeff,lagrange,stiffness)

    type(tContactState) :: ctState !< type tContactState
    integer(kind=kint)  :: iter
    integer(kind=kint)  :: etype, nnode !< type of master segment; number of nodes of master segment
    integer(kind=kint)  :: i, j, k, l
    real(kind=kreal)    :: normal(3), shapefunc(nnode) !< normal vector at target point; shape functions
    real(kind=kreal)    :: nTm((nnode + 1)*3) !< vector
    real(kind=kreal)    :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)    :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)    :: tf_trial(3), length_tft
    real(kind=kreal)    :: tangent(3), tTm((nnode + 1)*3)
    real(kind=kreal)    :: stiffness(:,:) !< contact stiffness matrix

    stiffness = 0.0d0

    call getShapeFunc( etype, ctState%lpos(:), shapefunc )

    normal(1:3) = ctState%direction(1:3)

    nTm(1:3) = normal(1:3)
    do i = 1, nnode
      nTm(i*3+1:i*3+3) = -shapefunc(i)*normal(1:3)
    enddo

    i = (nnode+1)*3 + 1
    do j = 1, (nnode+1)*3
      stiffness(i,j) = nTm(j);  stiffness(j,i) = nTm(j)
    enddo


    if( fcoeff /= 0.0d0 ) then
      if( lagrange>0.0d0 .or. iter==1 ) then

        do i = 1, nnode+1
          do j = 1, i
            do k = 1, 3
              do l = 1, 3
                stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - tPenalty*nTm((i-1)*3+k)*nTm((j-1)*3+l)
                if( k==l ) then
                  if(i==1 .and. j==1)then
                    stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) + tPenalty
                  elseif(i>1 .and. j==1)then
                    stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - tPenalty*shapefunc(i-1)
                  elseif(i>1 .and. j>1)then
                    stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) + tPenalty*shapefunc(i-1)*shapefunc(j-1)
                  endif
                endif
                if(i==j) cycle
                stiffness((j-1)*3+l,(i-1)*3+k) = stiffness((i-1)*3+k,(j-1)*3+l)
              enddo
            enddo
          enddo
        enddo

        if( ctstate%state == contactSlip ) then

          tf_trial(1:3) = ctstate%tangentForce_trial(1:3)
          length_tft = dsqrt(dot_product(tf_trial,tf_trial))
          tangent(1:3) = tf_trial(1:3)/length_tft

          tTm(1:3) = -tangent(1:3)
          do i = 1, nnode
            tTm(i*3+1:i*3+3) = shapefunc(i)*tangent(1:3)
          enddo

          do i = 1, nnode+1
            do j = 1, nnode+1
              do k = 1, 3
                do l = 1, 3
                  stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l)       &
                    + tPenalty*(-tTm((i-1)*3+k)*tTm((j-1)*3+l)  &
                    +tTm((i-1)*3+k)*nTm((j-1)*3+l)*dot_product(tangent,normal))
                enddo
              enddo
            enddo
          enddo
          stiffness(1:(nnode+1)*3,1:(nnode+1)*3) = (fcoeff*lagrange/length_tft)*stiffness(1:(nnode+1)*3,1:(nnode+1)*3)

          !
          !        do i = 1, (nnode + 1)*3
          !          do j = 1, i
          !            do k = 1, 3
          !              do l = 1, k
          !                stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - mut*tTm((i-1)*3+k)*tTm((j-1)*3+l))
          !                if( k/=l )stiffness((j-1)*3+l,(i-1)*3+k) = stiffness((i-1)*3+k,(j-1)*3+l)
          !              enddo !- l
          !            enddo !- k
          !          enddo !- j
          !        enddo !- i
          !        stiffness(1:(nnode+1)*3,1:(nnode+1)*3) = (fcoeff*dabs(lagrange)/length_tft)*stiffness(1:(nnode+1)*3,1:(nnode+1)*3)

          !         j = (nnode+1)*3 + 1
          !         do i = 1, (nnode+1)*3
          !           stiffness(i,j) = stiffness(i,j) - fcoeff*tTm(i)
          !         enddo
          stiffness(1:(nnode+1)*3,(nnode+1)*3+1) = stiffness(1:(nnode+1)*3,(nnode+1)*3+1) - fcoeff*tTm(1:(nnode+1)*3)

        endif
      endif
    endif

  end subroutine getContactStiffness

  !> \brief This subroutine obtains contact nodal force vector of each contact pair
  !! and assembles it into right-hand side vector to update non-equilibrated nodal force vector.
  subroutine fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)

    type(hecmwST_local_mesh)             :: hecMESH !< type hecmwST_local_mesh
    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange) :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(hecmwST_matrix)                 :: conMAT !< type hecmwST_matrix for contact part only
    integer(kind=kint) :: ctsurf, etype, nnode, ndLocal(21) !< contents of type tContact
    integer(kind=kint) :: i, j, k, nlag, id_lagrange, algtype
    real(kind=kreal)   :: ndCoord(21*3) !< nodal coordinates
    real(kind=kreal)   :: ndu(21*3), ndDu(21*3) !< nodal displacement and its increment
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)   :: ctNForce(21*3+1)     !< nodal normal contact force vector
    real(kind=kreal)   :: ctTForce(21*3+1)     !< nodal tangential contact force vector

    integer(kind=kint) :: cstep !< current calculation step
    integer(kind=kint) :: grpid

    id_lagrange = 0
    if( associated(fstrSOLID%CONT_NFORCE) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
    if( associated(fstrSOLID%CONT_FRIC) ) fstrSOLID%CONT_FRIC(:) = 0.d0
    if( associated(fstrSOLID%EMBED_NFORCE) ) fstrSOLID%EMBED_NFORCE(:) = 0.d0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)
        do k = 1, nnode+1
          ndDu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%dunode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          ndu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndDu((k-1)*3+1:(k-1)*3+3)
          ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3) + ndu((k-1)*3+1:(k-1)*3+3)
        enddo

        do k=1,nlag
          id_lagrange = id_lagrange + 1
          lagrange = hecLagMAT%Lagrange(id_lagrange)
  
          fstrSOLID%contacts(i)%states(j)%multiplier(k) = hecLagMAT%Lagrange(id_lagrange)
    
          if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
            ! Obtain contact nodal force vector of contact pair
            call getContactNodalForce(etype,nnode,ndCoord,ndDu,fstrSOLID%contacts(i)%states(j),    &
            fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,ctNForce,ctTForce,.true.)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce,  &
              &  conMAT,fstrSOLID%CONT_NFORCE,fstrSOLID%CONT_FRIC)
          else if( algtype == CONTACTTIED ) then
            call getTiedNodalForce(etype,nnode,k,ndu,fstrSOLID%contacts(i)%states(j),lagrange,ctNForce,ctTForce)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,1.d0,ctNForce,ctTForce,  &
              &  conMAT,fstrSOLID%CONT_NFORCE,fstrSOLID%CONT_FRIC)
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
        etype = fstrSOLID%embeds(i)%master(ctsurf)%etype
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
    
          call getTiedNodalForce(etype,nnode,k,ndu,fstrSOLID%embeds(i)%states(j),lagrange,ctNForce,ctTForce)
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

  !> \brief This subroutine obtains contact nodal force vector of contact pair
  subroutine getTiedNodalForce(etype,nnode,idof,ndu,ctState,lagrange,ctNForce,ctTForce)
    integer(kind=kint)  :: etype !< type of master segment
    integer(kind=kint)  :: nnode !< number of nodes of master segment
    integer(kind=kint)  :: idof
    real(kind=kreal)   :: ndu((nnode + 1)*3) !< nodal displacement
    type(tContactState) :: ctState !< type tContactState
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)       :: ctNForce((nnode+1)*3+1)  !< tied contact force vector
    real(kind=kreal)       :: ctTForce((nnode+1)*3+1)  !< tied contact force vector

    integer(kind=kint) :: i, j
    real(kind=kreal)   :: normal(3), shapefunc(nnode) !< normal vector at target point; shape functions
    real(kind=kreal)   :: nTm((nnode + 1)*3) !< normal vector
    real(kind=kreal)   :: tTm((nnode + 1)*3) !< tangential vector

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    call getShapeFunc( etype, ctState%lpos, shapefunc )

    normal(1:3) = ctState%direction(1:3)

    nTm(1:3) = -normal(idof)*normal(1:3)
    do i = 1, nnode
      nTm(i*3+1:i*3+3) = -shapefunc(i)*nTm(1:3)
    enddo
    tTm = 0.d0
    tTm(idof) = -1.0
    tTm(1:3) = tTm(1:3)-nTm(1:3)
    do i = 1, nnode
      tTm(i*3+1:i*3+3) = -shapefunc(i)*tTm(1:3)
    enddo

    do j = 1, (nnode+1)*3
      ctNForce(j) = lagrange*nTm(j)
      ctTForce(j) = lagrange*tTm(j)
    enddo
    j = (nnode+1)*3 + 1
    ctNForce(j) = dot_product(nTm,ndu)
    ctTForce(j) = dot_product(tTm,ndu)
  end subroutine 

  !> \brief This subroutine obtains contact nodal force vector of contact pair
  subroutine getContactNodalForce(etype,nnode,ndCoord,ndDu,ctState,tPenalty,fcoeff,lagrange,ctNForce,ctTForce,cflag)

    type(tContactState) :: ctState !< type tContactState
    integer(kind=kint) :: etype, nnode !< type of master segment; number of nodes of master segment
    integer(kind=kint) :: i, j
    real(kind=kreal)   :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)   :: ndCoord((nnode + 1)*3), ndDu((nnode + 1)*3) !< nodal coordinates; nodal displacement increment
    real(kind=kreal)   :: normal(3), shapefunc(nnode) !< normal vector at target point; shape functions
    real(kind=kreal)   :: nTm((nnode + 1)*3) !< vector
    real(kind=kreal)   :: tf_trial(3), length_tft, tangent(3), tf_final(3)
    real(kind=kreal)       :: ctNForce((nnode+1)*3+1)                     !< contact force vector
    real(kind=kreal)       :: ctTForce((nnode+1)*3+1)                     !< contact force vector
    logical            :: cflag  !< is necessary to update tangentForce_final

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    call getShapeFunc( etype, ctState%lpos, shapefunc )

    normal(1:3) = ctState%direction(1:3)

    nTm(1:3) = -normal(1:3)
    do i = 1, nnode
      nTm(i*3+1:i*3+3) = shapefunc(i)*normal(1:3)
    enddo

    do j = 1, (nnode+1)*3
      ctNForce(j) = lagrange*nTm(j)
    enddo
    j = (nnode+1)*3 + 1
    ctNForce(j) = dot_product(nTm,ndCoord)

    if(fcoeff /= 0.0d0 .and. lagrange > 0.0d0)then

      if( cflag ) then !calc tf_final and set it to tangentForce_final
        call getTrialFricForceANDcheckFricState(nnode,tPenalty,fcoeff,lagrange,normal,shapefunc,nTm,ndDu,ctstate)

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

      ctTForce(1:3) = -tf_final(1:3)
      do j = 1, nnode
        ctTForce(j*3+1:j*3+3) = shapefunc(j)*tf_final(1:3)
      enddo

    endif

  end subroutine getContactNodalForce


  !> \brief This subroutine calculates trial friction force and checks friction state
  subroutine getTrialFricForceANDcheckFricState(nnode,tPenalty,fcoeff,lagrange,normal,shapefunc,nTm,ndDu,ctstate)

    type(tContactState) :: ctState !< type tContactState
    integer(kind=kint) :: nnode !< number of nodes of master segment
    integer(kind=kint) :: i, j
    real(kind=kreal)   :: fcoeff, tPenalty !< friction coefficient; tangential penalty
    real(kind=kreal)   :: lagrange !< value of Lagrange multiplier
    real(kind=kreal)   :: ndDu((nnode + 1)*3) !< nodal displacement increment
    real(kind=kreal)   :: normal(3), shapefunc(nnode) !< normal vector at target point; shape functions
    real(kind=kreal)   :: nTm((nnode + 1)*3) !< vector
    real(kind=kreal)   :: dotP
    real(kind=kreal)   :: relativeDisp(3) !< relative displacement
    real(kind=kreal)   :: tf_yield

    relativeDisp = 0.0d0

    dotP = dot_product(nTm,ndDu)
    do i = 1, 3
      relativeDisp(i) = - ndDu(i)
      do j = 1, nnode
        relativeDisp(i) = relativeDisp(i) + shapefunc(j)*ndDu(j*3+i)
      enddo
      relativeDisp(i) = relativeDisp(i) - dotP*normal(i)
      ctstate%reldisp(i) = -relativeDisp(i)
      ctstate%tangentForce_trial(i) = ctstate%tangentForce1(i) -tPenalty*relativeDisp(i)
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

    hecMAT%B(np*ndof+id_lagrange) = ctNForce((nnode+1)*3+1)+ctTForce((nnode+1)*3+1)

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
    integer(kind=kint) :: i, j, k, l, id_lagrange, lnod, grpid
    integer(kind=kint) :: ctsurf, etype, nnode, ndLocal(21) !< contents of type tContact
    real(kind=kreal)   :: ndu(21*3), ndCoord(21*3), lagrange !< nodal coordinates; value of Lagrange mutiplier
    real(kind=kreal)   :: normal(3), shapefunc(21) !< normal vector; shape functions
    real(kind=kreal)   :: nTm(10*3) !< vector
    real(kind=kreal)   :: tf_final(3) !< final friciton force vector
    real(kind=kreal)   :: ctForce(21*3 + 1) !< initial nodal contact force vector

    integer(kind=kint) :: algtype, nlag    
    real(kind=kreal)   :: ctNForce(21*3+1)     !< nodal normal contact force vector
    real(kind=kreal)   :: ctTForce(21*3+1)     !< nodal tangential contact force vector

    np = hecMAT%NP ; ndof = hecMAT%NDOF

    id_lagrange = 0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype
      nlag = fstr_get_num_lagrange_pernode(algtype)

      do j = 1, size(fstrSOLID%contacts(i)%slave)

        if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

        ctsurf = fstrSOLID%contacts(i)%states(j)%surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
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
    
          if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
            ! Obtain contact nodal force vector of contact pair
            call getContactNodalForce(etype,nnode,ndCoord,ndu,fstrSOLID%contacts(i)%states(j),    &
            fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,ctNForce,ctTForce,.false.)
            ! Update non-eqilibrited force vector
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,lagrange,ctNForce,ctTForce, & 
               &  hecMAT,fstrSOLID%CONT_NFORCE,fstrSOLID%CONT_FRIC)
          else if( algtype == CONTACTTIED ) then
            call getTiedNodalForce(etype,nnode,k,ndu,fstrSOLID%contacts(i)%states(j),lagrange,ctNForce,ctTForce)
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
        etype = fstrSOLID%embeds(i)%master(ctsurf)%etype
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
    
          call getTiedNodalForce(etype,nnode,k,ndu,fstrSOLID%embeds(i)%states(j),lagrange,ctNForce,ctTForce)
          ! Update non-eqilibrited force vector
          call update_NDForce_contact(nnode,ndLocal,id_lagrange,-1.d0,ctNForce,ctTForce,hecMAT,fstrSOLID%EMBED_NFORCE)

        enddo
      enddo
    enddo

  end subroutine fstr_ass_load_contact

end module m_addContactStiffness



