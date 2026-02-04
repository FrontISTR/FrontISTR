!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides functions of contact calculation
module m_fstr_contact_lib
  use hecmw
  use elementInfo
  use mContactDef
  implicit none

contains

  !> Transfer contact condition int mpc bundary conditions
  subroutine contact2mpcval( cstate, etype, nnode, mpcval )
    type(tContactState), intent(in) :: cstate !< contact state
    integer, intent(in)             :: etype !< type of contacting surface
    integer, intent(in)             :: nnode !< number of elemental nodes
    real(kind=kreal), intent(out)   :: mpcval(nnode*3 + 4) !< MPC constraint

    integer          :: i,j
    real(kind=kreal) :: shapefunc(nnode)

    call getShapeFunc( etype, cstate%lpos(1:2), shapefunc )
    mpcval(1:3) = cstate%direction(1:3)
    do i=1,nnode
      do j=1,3
        mpcval( i*3+j ) = -cstate%direction(j)*shapefunc(i)
      enddo
    enddo
    mpcval( 3*nnode+4 )=cstate%distance
  end subroutine

  !> This subroutine calculate the metric tensor of a elemental surface
  subroutine getMetricTensor( pos, etype, ele, tensor )
    real(kind=kreal), intent(in)  :: pos(2)        !< current position(local coordinate)
    integer, intent(in)           :: etype         !< surface element type
    real(kind=kreal), intent(in)  :: ele(:,:)      !< elemental coordinates
    real(kind=kreal), intent(out) :: tensor(2,2)   !< metric tensor

    integer          :: nn
    real(kind=kreal) :: tangent(3,2)
    nn= getNumberOfNodes(etype)
    call TangentBase( etype, nn, pos, ele, tangent )
    tensor(1,1)= dot_product( tangent(:,1), tangent(:,1) )
    tensor(1,2)= dot_product( tangent(:,1), tangent(:,2) )
    tensor(2,1)= dot_product( tangent(:,2), tangent(:,1) )
    tensor(2,2)= dot_product( tangent(:,2), tangent(:,2) )
  end subroutine

  !> This subroutine calculate the relation between global disp and displacement
  !> along natural coordinate of master surface supposing penetration is small
  subroutine DispIncreMatrix( pos, etype, nnode, ele, tangent, tensor, matrix )
    real(kind=kreal), intent(in)  :: pos(2)        !< current position(local coordinate)
    integer, intent(in)           :: etype         !< surface element type
    integer, intent(in)           :: nnode         !< number of nodes in surface
    real(kind=kreal), intent(in)  :: ele(3,nnode)  !< elemental coordinates
    real(kind=kreal), intent(out) :: tangent(3,2)  !< tangent basis
    real(kind=kreal), intent(out) :: tensor(2,2)   !< metric tensor
    real(kind=kreal), intent(out) :: matrix(2,nnode*3+3) !< relation between local and global disp increment

    integer          :: i,j
    real(kind=kreal) :: det
    real(kind=kreal) :: shapefunc(nnode), t1(nnode*3+3), t2(nnode*3+3)
    call TangentBase( etype, nnode, pos, ele, tangent )
    tensor(1,1)= dot_product( tangent(:,1), tangent(:,1) )
    tensor(1,2)= dot_product( tangent(:,1), tangent(:,2) )
    tensor(2,1)= dot_product( tangent(:,2), tangent(:,1) )
    tensor(2,2)= dot_product( tangent(:,2), tangent(:,2) )
    det = tensor(1,1)*tensor(2,2)-tensor(1,2)*tensor(2,1)
    if( det==0.d0 ) stop "Error in calculate DispIncreMatrix"
    !  inverse(1,1) = tensor(2,2)/det
    !  inverse(1,2) = -tensor(1,2)/det
    ! inverse(2,1) = -tensor(2,1)/det
    !  inverse(2,2) = tensor(1,1)/det

    call getShapeFunc( etype, pos(:), shapefunc )
    forall( j=1:3 )
      t1( j ) = tangent(j,1)
      t2( j ) = tangent(j,2)
    end forall
    forall( i=1:nnode, j=1:3 )
      t1( i*3+j ) = -tangent(j,1)*shapefunc(i)
      t2( i*3+j ) = -tangent(j,2)*shapefunc(i)
    end forall
    !matrix( 1:2,: ) = matmul( inverse(:,:), matrix )
    matrix(1,:) = (tensor(2,2)*t1(:)-tensor(1,2)*t2(:))/det
    matrix(2,:) = (tensor(1,1)*t2(:)-tensor(2,1)*t1(:))/det
    tangent(:,1) = tangent(:,1)/dsqrt(dot_product(tangent(:,1),tangent(:,1)))
    tangent(:,2) = tangent(:,2)/dsqrt(dot_product(tangent(:,2),tangent(:,2)))
  end subroutine

  !> This subroutine find the projection of a slave point onto master surface
  subroutine project_Point2Element(xyz,etype,nn,elemt,reflen,cstate,isin,distclr,ctpos,localclr)
    real(kind=kreal),intent(in)       :: xyz(3)        !< Coordinates of a spacial point, whose projecting point is to be computed
    integer, intent(in)               :: etype         !< surface element type
    integer, intent(in)               :: nn            !< number of elemental nodes
    real(kind=kreal),intent(in)       :: elemt(:,:)   !< nodes coordinates of surface element
    real(kind=kreal),intent(in)       :: reflen        !< reference length of surface element
    type(tContactState),intent(inout) :: cstate        !< Recorde of contact information
    logical, intent(out)              :: isin          !< in contact or not
    real(kind=kreal), intent(in)      :: distclr       !< clearance of contact distance
    real(kind=kreal), optional        :: ctpos(2)      !< curr contact position( natural coord )
    real(kind=kreal), optional        :: localclr      !< clearance of contact local coord

    integer          ::  count,order, initstate
    real(kind=kreal)  ::  determ, inverse(2,2)
    real(kind=kreal)  ::  sfunc(nn), curv(3,2,2)
    real(kind=kreal)  ::  r(2), dr(2), r_tmp(2)        ! natural coordinate
    real(kind=kreal)  ::  xyz_out(3)                   ! curr. projection position
    real(kind=kreal)  ::  dist_last,dist_now, dxyz(3)  ! dist between the point and its projection
    real(kind=kreal)  ::  tangent(3,2)                 ! base vectors in tangent space
    real(kind=kreal)  ::  dF(2),d2F(2,2),normal(3)
    real(kind=kreal),parameter :: eps = 1.0D-8
    real(kind=kreal)  ::  clr, tol, factor

    initstate = cstate%state
    clr = 1.d-4
    if( present( localclr ) ) clr=localclr
    if( present( ctpos ) ) then
      r(:)= ctpos
    else
      call getElementCenter( etype, r(:) )
    endif

    tol = 1.0D0
    do count=1,100
      call getShapeFunc( etype, r, sfunc )
      xyz_out = matmul( elemt(1:3,1:nn), sfunc )
      dxyz(1:3) = xyz_out(1:3) - xyz(1:3)
      dist_last = dot_product( dxyz, dxyz(:) )

      call TangentBase( etype, nn, r, elemt(1:3,1:nn), tangent )
      call Curvature( etype, nn, r, elemt(1:3,1:nn), curv )

      !     dF(1:2)
      dF(1:2) = -matmul( dxyz(:), tangent(:,:) )
      !     d2F(1:2,1:2)
      d2F(1,1)= dot_product( tangent(:,1), tangent(:,1) ) - dot_product( dxyz, curv(:,1,1) )
      d2F(1,2)= dot_product( tangent(:,1), tangent(:,2) ) - dot_product( dxyz, curv(:,1,2) )
      d2F(2,1)= dot_product( tangent(:,2), tangent(:,1) ) - dot_product( dxyz, curv(:,2,1) )
      d2F(2,2)= dot_product( tangent(:,2), tangent(:,2) ) - dot_product( dxyz, curv(:,2,2) )

      !     inverse of d2F
      determ = d2F(1,1)*d2F(2,2) - d2F(1,2)*d2f(2,1)
      if( determ==0.d0 ) stop "Math error in contact searching"
      inverse(1,1) = d2F(2,2) / determ
      inverse(2,2) = d2F(1,1) / determ
      inverse(1,2) = -d2F(1,2) / determ
      inverse(2,1) = -d2F(2,1) / determ
      dr=matmul(inverse,dF)

      tol=dot_product(dr,dr)
      if( dsqrt(tol)> 3.d0 ) then   ! too far away
        r= -100.d0; exit
      endif

      factor = 1.d0
      do order=1,10
        r_tmp(1:2) = r(1:2) + factor*dr(1:2)
        call getShapeFunc( etype, r_tmp, sfunc )
        xyz_out(1:3) = matmul( elemt(1:3,1:nn), sfunc(:) )
        dxyz(1:3) = xyz(1:3)-xyz_out(:)
        dist_now = dot_product( dxyz, dxyz )
        if(dist_now <= dist_last) exit
        factor = factor*0.7D0
      enddo
      r(1:2) = r_tmp(1:2)

      if( tol<eps ) exit
    enddo

    isin = .false.
    cstate%state = CONTACTFREE
    if( isInsideElement( etype, r, clr )>=0 ) then
      dxyz(:)=xyz_out(:)-xyz(:)
      normal(:) = SurfaceNormal( etype, nn, r, elemt(1:3,1:nn) )
      normal(:) = normal(:)/dsqrt( dot_product(normal, normal) )
      do count = 1,3
        if( dabs(normal(count))<1.D-10 ) normal(count) =0.d0
        if( dabs(1.d0-dabs(normal(count)))<1.D-10 ) normal(count) =sign(1.d0, normal(count))
      enddo
      cstate%distance = dot_product( dxyz, normal )
    
      if( cstate%interference_flag == C_IF_SLAVE) then
        if( cstate%init_pos == 0.d0 .and. cstate%distance < cstate%end_pos) then
          cstate%init_pos = cstate%distance
          cstate%time_factor = (cstate%end_pos - cstate%distance) / cstate%time_factor
          cstate%shrink_factor = cstate%distance
        end if
      end if

      if( cstate%interference_flag == 0)then ! not shrink-node
        if( cstate%distance < distclr*reflen .and. cstate%distance > -5.0d-01*reflen ) isin = .true.
      else
        if( cstate%distance < cstate%shrink_factor + distclr*reflen ) isin = .true.
      end if

      if( isin ) then
        if( initstate== CONTACTFREE ) then
          cstate%state = CONTACTSTICK
        else
          cstate%state = initstate
        endif
        cstate%gpos(:)=xyz_out(:)
        cstate%lpos(1:2)=r(:)
        cstate%direction(:) = normal(:)
        cstate%wkdist = cstate%distance
      endif
    endif
  end subroutine project_Point2Element

  !> This subroutine find the projection of a slave point onto master surface
  subroutine project_Point2SolidElement(xyz,etype,nn,elemt,reflen,cstate,isin,distclr,ctpos,localclr)
    use m_utilities
    real(kind=kreal),intent(in)       :: xyz(3)        !< Coordinates of a spacial point, whose projecting point is to be computed
    integer, intent(in)               :: etype         !< surface element type
    integer, intent(in)               :: nn            !< number of elemental nodes
    real(kind=kreal),intent(in)       :: elemt(:,:)   !< nodes coordinates of surface element
    real(kind=kreal),intent(in)       :: reflen        !< reference length of surface element
    type(tContactState),intent(inout) :: cstate        !< Recorde of contact information
    logical, intent(out)              :: isin          !< in contact or not
    real(kind=kreal), intent(in)      :: distclr       !< clearance of contact distance
    real(kind=kreal), optional        :: ctpos(3)      !< curr contact position( natural coord )
    real(kind=kreal), optional        :: localclr      !< clearance of contact local coord

    integer          ::  count,order, initstate
    real(kind=kreal)  ::  determ, inverse(3,3)
    real(kind=kreal)  ::  sfunc(nn), deriv(nn,3)
    real(kind=kreal)  ::  r(3), dr(3), r_tmp(3)        ! natural coordinate
    real(kind=kreal)  ::  xyz_out(3)                   ! curr. projection position
    real(kind=kreal)  ::  dist_last,dist_now, dxyz(3)  ! dist between the point and its projection
    real(kind=kreal)  ::  tangent(3,2)                 ! base vectors in tangent space
    real(kind=kreal)  ::  dF(2),d2F(2,2),normal(3)
    real(kind=kreal),parameter :: eps = 1.0D-8
    real(kind=kreal)  ::  clr, tol, factor

    initstate = cstate%state
    clr = 1.d-4
    if( present( localclr ) ) clr=localclr
    if( present( ctpos ) ) then
      r(:)= ctpos
    else
      call getElementCenter( etype, r(:) )
    endif

    tol = 1.0D0
    do count=1,100
      call getShapeFunc( etype, r, sfunc )
      xyz_out = matmul( elemt(1:3,1:nn), sfunc(1:nn) )
      dxyz(1:3) = xyz_out(1:3) - xyz(1:3)
      dist_last = dot_product( dxyz, dxyz(:) )

      call getShapeDeriv( etype, r, deriv )
      inverse(1:3,1:3) = matmul(elemt(1:3,1:nn),deriv(1:nn,1:3))
      call calInverse(3, inverse(1:3,1:3))
      dr(1:3) = -matmul(inverse(1:3,1:3),dxyz(1:3))

      tol=dot_product(dr,dr)
      if( count > 1 .and. dsqrt(tol)> 3.d0 ) then   ! too far away
        r= -100.d0; exit
      endif

      factor = 1.d0
      do order=1,10
        r_tmp(1:3) = r(1:3) + factor*dr(1:3)
        call getShapeFunc( etype, r_tmp, sfunc )
        xyz_out(1:3) = matmul( elemt(1:3,1:nn), sfunc(1:nn) )
        dxyz(1:3) = xyz(1:3)-xyz_out(1:3)
        dist_now = dot_product( dxyz, dxyz )
        if(dist_now <= dist_last) exit
        factor = factor*0.7D0
      enddo
      r(1:3) = r_tmp(1:3)

      if( tol<eps ) exit
    enddo

    isin = .false.
    cstate%state = CONTACTFREE
    if( isInside3DElement( etype, r, clr )>=0 ) then
      dxyz(:)=xyz_out(:)-xyz(:)
      cstate%distance = dsqrt(dot_product( dxyz, dxyz ))

      if( cstate%distance < distclr ) isin = .true.

      if( isin ) then
        if( initstate== CONTACTFREE ) then
          cstate%state = CONTACTSTICK
        else
          cstate%state = initstate
        endif
        cstate%gpos(:)=xyz_out(:)
        cstate%direction(:) = (/1.d0,0.d0,0.d0/)
        cstate%lpos(1:3)=r(:)
      endif
    endif
  end subroutine

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

  subroutine set_shrink_factor(ctime, cstate, etime, if_type)
    real(kind=kreal),intent(in)         :: ctime, etime
    type(tContactState), intent(inout)  :: cstate        !< Recorde of contact information
    integer, intent(in)                 :: if_type

    if (if_type == C_IF_SLAVE .and. cstate%init_pos == 0.d0) then
      cstate%shrink_factor = 0.0d0; return
    end if
    cstate%shrink_factor = cstate%time_factor*ctime + cstate%init_pos
    if(ctime >= etime) cstate%shrink_factor = cstate%end_pos

  end subroutine set_shrink_factor

  !> Wrapper for project_Point2Element that takes tSurfElement structure
  !! This subroutine handles element coordinate extraction from tSurfElement
  subroutine project_Point2SurfElement(xyz, surf, currpos, cstate, isin, distclr, ctpos, localclr)
    real(kind=kreal), intent(in)              :: xyz(3)        !< coordinates of slave point
    type(tSurfElement), intent(in)            :: surf          !< surface element structure
    real(kind=kreal), intent(in)              :: currpos(:)    !< current coordinate of all nodes
    type(tContactState), intent(inout)        :: cstate        !< contact state
    logical, intent(out)                      :: isin          !< in contact or not
    real(kind=kreal), intent(in)              :: distclr       !< clearance of contact distance
    real(kind=kreal), optional, intent(in)    :: ctpos(2)      !< current contact position (natural coord)
    real(kind=kreal), optional, intent(in)    :: localclr      !< clearance of contact local coord

    integer(kind=kint) :: nn, j, iSS
    real(kind=kreal)   :: elem(3, l_max_elem_node)

    ! Extract element coordinates from surf structure
    nn = size(surf%nodes)
    do j = 1, nn
      iSS = surf%nodes(j)
      elem(1:3, j) = currpos(3*iSS-2:3*iSS)
    enddo

    call project_Point2Element(xyz, surf%etype, nn, elem, surf%reflen, cstate, &
                               isin, distclr, ctpos, localclr)

  end subroutine project_Point2SurfElement

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_contact_state( flag_ctAlgo, contact, currpos, currdisp, ndforce, infoCTChange, &
      nodeID, elemID, is_init, active, mu, B )
    character(len=9), intent(in)                     :: flag_ctAlgo  !< contact analysis algorithm flag
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: ndforce(:)   !< nodal force
    integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                              :: is_init      !< whether initial scan or not
    logical, intent(out)                             :: active       !< if any in contact
    real(kind=kreal), intent(in)                     :: mu           !< penalty
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3)
    real(kind=kreal)    :: nlforce, slforce(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_cand, is_present_B
    real(kind=kreal), pointer :: Bp(:)

    if( is_init ) then
      distclr = contact%cparam%DISTCLR_INIT
    else
      distclr = contact%cparam%DISTCLR_FREE
      if( contact%algtype == CONTACTTIED ) then
        active = .false.
        do i= 1, size(contact%slave)
          if( contact%states(i)%state==CONTACTSTICK ) then
            active = .true.
            exit
          endif
        enddo
      endif
    endif

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(contact%slave)))
    contact_surf(:) = huge(0)
    do i = 1, size(contact%slave)
      states_prev(i) = contact%states(i)%state
    enddo

    call update_surface_box_info( contact%master, currpos )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )

    ! for gfortran-10: optional parameter seems not allowed within omp parallel
    is_present_B = present(B)
    if( is_present_B ) Bp => B

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,is_cand,idm,etype,isin, &
      !$omp&         bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax,is_present_B) &
      !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,mu,nodeID,elemID,Bp,distclr,contact_surf,is_init) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      if(contact%if_type /= 0) call set_shrink_factor(contact%ctime, contact%states(i), contact%if_etime, contact%if_type)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
        slforce(1:3)=ndforce(3*slave-2:3*slave)
        id = contact%states(i)%surface
        nlforce = contact%states(i)%multiplier(1)

        ! update direction of TIED contact
        if( contact%algtype == CONTACTTIED ) then
          call update_direction( i, contact, currpos )
          if (.not.is_init) cycle
        endif

        if( nlforce < contact%cparam%TENSILE_FORCE ) then
          contact%states(i)%state = CONTACTFREE
          contact%states(i)%multiplier(:) = 0.d0
          write(*,'(A,i10,A,i10,A,e12.3)') "Node",nodeID(slave)," free from contact with element", &
            elemID(contact%master(id)%eid), " with tensile force ", nlforce
          cycle
        endif
        if( contact%algtype /= CONTACTFSLID .or. (.not. is_present_B) ) then   ! small slide problem
          contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
        else
          call track_contact_position( flag_ctAlgo, i, contact, currpos, currdisp, mu, infoCTChange, nodeID, elemID, Bp )
          if( contact%states(i)%state /= CONTACTFREE ) then
            id = contact%states(i)%surface
            contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
          endif
        endif

      else if( contact%states(i)%state==CONTACTFREE ) then
        if( contact%algtype == CONTACTTIED .and. .not. is_init ) cycle 
        coord(:) = currpos(3*slave-2:3*slave)

        ! get master candidates from bucketDB
        bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
        nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
        if (nCand == 0) cycle
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)

        nMasterMax = nCand
        allocate(indexMaster(nMasterMax))
        nMaster = 0

        ! narrow down candidates
        do idm= 1, nCand
          id = indexCand(idm)
          if (.not. is_in_surface_box( contact%master(id), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          nMaster = nMaster + 1
          indexMaster(nMaster) = id
        enddo
        deallocate(indexCand)

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif

        do idm = 1,nMaster
          id = indexMaster(idm)
          call project_Point2SurfElement( coord, contact%master(id), currpos, &
            contact%states(i), isin, distclr, localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          etype = contact%master(id)%etype
          iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
          if( iSS>0 ) &
            call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%lpos(1:2), &
              contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = elemID(contact%master(id)%eid)
          write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3,A,i6)') "Node",nodeID(slave)," contact with element", &
            elemID(contact%master(id)%eid),       &
            " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(1:2), &
            " along direction ", contact%states(i)%direction," rank=",hecmw_comm_get_rank()
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    if( contact%algtype == CONTACTTIED .and. .not. is_init ) then
      deallocate(contact_surf)
      deallocate(states_prev)
      return
    endif

    call hecmw_contact_comm_allreduce_i(contact%comm, contact_surf, HECMW_MIN)
    nactive = 0
    do i = 1, size(contact%slave)
      if (contact%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
        id = contact%states(i)%surface
        if (abs(contact_surf(contact%slave(i))) /= elemID(contact%master(id)%eid)) then ! that is in contact with other surface
          contact%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i10,A,i6,A,i6,A)') "Node",nodeID(contact%slave(i))," contact with element", &
            &  elemID(contact%master(id)%eid), " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else
          nactive = nactive + 1
        endif
      endif
      if (states_prev(i) == CONTACTFREE .and. contact%states(i)%state /= CONTACTFREE) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (states_prev(i) /= CONTACTFREE .and. contact%states(i)%state == CONTACTFREE) then
        infoCTChange%contact2free = infoCTChange%contact2free + 1
      endif
    enddo
    active = (nactive > 0)
    deallocate(contact_surf)
    deallocate(states_prev)
  end subroutine scan_contact_state

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_embed_state( flag_ctAlgo, embed, currpos, currdisp, ndforce, infoCTChange, &
    nodeID, elemID, is_init, active, mu, B )
    character(len=9), intent(in)                     :: flag_ctAlgo  !< contact analysis algorithm flag
    type( tContact ), intent(inout)                  :: embed      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: ndforce(:)   !< nodal force
    integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                              :: is_init      !< whether initial scan or not
    logical, intent(out)                             :: active       !< if any in contact
    real(kind=kreal), intent(in)                     :: mu           !< penalty
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
    real(kind=kreal)    :: nlforce, slforce(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_cand, is_present_B
    real(kind=kreal), pointer :: Bp(:)

    if( is_init ) then
      distclr = embed%cparam%DISTCLR_INIT
    else
      distclr = embed%cparam%DISTCLR_FREE
      active = .false.
      do i= 1, size(embed%slave)
        if( embed%states(i)%state==CONTACTSTICK ) then
          active = .true.
          exit
        endif
      enddo
      return
    endif

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(embed%slave)))
    contact_surf(:) = huge(0)
    do i = 1, size(embed%slave)
      states_prev(i) = embed%states(i)%state
    enddo

    call update_surface_box_info( embed%master, currpos )
    call update_surface_bucket_info( embed%master, embed%master_bktDB )

    ! for gfortran-10: optional parameter seems not allowed within omp parallel
    is_present_B = present(B)
    if( is_present_B ) Bp => B

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,elem,is_cand,idm,etype,isin, &
      !$omp&         bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax,is_present_B) &
      !$omp& shared(embed,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,mu,nodeID,elemID,Bp,distclr,contact_surf) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(embed%slave)
      slave = embed%slave(i)
      if( embed%states(i)%state==CONTACTFREE ) then
        coord(:) = currpos(3*slave-2:3*slave)

        ! get master candidates from bucketDB
        bktID = bucketDB_getBucketID(embed%master_bktDB, coord)
        nCand = bucketDB_getNumCand(embed%master_bktDB, bktID)
        if (nCand == 0) cycle
        allocate(indexCand(nCand))
        call bucketDB_getCand(embed%master_bktDB, bktID, nCand, indexCand)

        nMasterMax = nCand
        allocate(indexMaster(nMasterMax))
        nMaster = 0

        ! narrow down candidates
        do idm= 1, nCand
          id = indexCand(idm)
          if (.not. is_in_surface_box( embed%master(id), coord(1:3), embed%cparam%BOX_EXP_RATE )) cycle
          nMaster = nMaster + 1
          indexMaster(nMaster) = id
        enddo
        deallocate(indexCand)

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif

        do idm = 1,nMaster
          id = indexMaster(idm)
          etype = embed%master(id)%etype
          if( mod(etype,10) == 2 ) etype = etype - 1 !search by 1st-order shape function
          nn = getNumberOfNodes(etype)
          do j=1,nn
            iSS = embed%master(id)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2SolidElement( coord,etype,nn,elem,embed%master(id)%reflen,embed%states(i), &
            isin,distclr,localclr=embed%cparam%CLEARANCE )
          if( .not. isin ) cycle
          embed%states(i)%surface = id
          embed%states(i)%multiplier(:) = 0.d0
          contact_surf(embed%slave(i)) = elemID(embed%master(id)%eid)
          write(*,'(A,i10,A,i10,A,3f7.3,A,i6)') "Node",nodeID(slave)," embeded to element", &
            elemID(embed%master(id)%eid), " at ",embed%states(i)%lpos(:)," rank=",hecmw_comm_get_rank()
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    call hecmw_contact_comm_allreduce_i(embed%comm, contact_surf, HECMW_MIN)
    nactive = 0
    do i = 1, size(embed%slave)
      if (embed%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
        id = embed%states(i)%surface
        if (abs(contact_surf(embed%slave(i))) /= elemID(embed%master(id)%eid)) then ! that is in contact with other surface
          embed%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i10,A,i6,A,i6,A)') "Node",nodeID(embed%slave(i))," contact with element", &
            &  elemID(embed%master(id)%eid), " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else
          nactive = nactive + 1
        endif
      endif
      if (states_prev(i) == CONTACTFREE .and. embed%states(i)%state /= CONTACTFREE) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (states_prev(i) /= CONTACTFREE .and. embed%states(i)%state == CONTACTFREE) then
        infoCTChange%contact2free = infoCTChange%contact2free + 1
      endif
    enddo
    active = (nactive > 0)
    deallocate(contact_surf)
    deallocate(states_prev)
  end subroutine scan_embed_state


  !> Calculate averaged nodal normal
  subroutine cal_node_normal( csurf, isin, surf, currpos, lpos, normal )
    use elementInfo, only:getVertexCoord, SurfaceNormal
    integer, intent(in)            :: csurf       !< current surface element
    integer, intent(in)            :: isin        !< return value from isInsideElement()
    type(tSurfElement), intent(in) :: surf(:)     !< surface elements
    real(kind=kreal), intent(in)   :: currpos(:)  !< current coordinate of each nodes
    real(kind=kreal), intent(in)   :: lpos(:)     !< local coordinate of contact position
    real(kind=kreal), intent(out)  :: normal(3)   !< averaged node nomral
    integer(kind=kint) :: cnode, i, j, cnt, nd1, gn, etype, iSS, nn,cgn
    real(kind=kreal) :: cnpos(2), elem(3, l_max_elem_node)
    integer(kind=kint) :: cnode1, cnode2, gn1, gn2, nsurf, cgn1, cgn2, isin_n
    real(kind=kreal) :: x=0, normal_n(3), lpos_n(2)

    if( 1 <= isin .and. isin <= 4 ) then  ! corner
      cnode = isin
      gn = surf(csurf)%nodes(cnode)
      etype = surf(csurf)%etype
      call getVertexCoord( etype, cnode, cnpos )
      nn = size( surf(csurf)%nodes )
      do j=1,nn
        iSS = surf(csurf)%nodes(j)
        elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      enddo
      normal = SurfaceNormal( etype, nn, cnpos, elem )
      cnt = 1
      do i=1,surf(csurf)%n_neighbor
        nd1 = surf(csurf)%neighbor(i)
        nn = size( surf(nd1)%nodes )
        etype = surf(nd1)%etype
        cgn = 0
        do j=1,nn
          iSS = surf(nd1)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          if( iSS==gn ) cgn=j
        enddo
        if( cgn>0 ) then
          call getVertexCoord( etype, cgn, cnpos )
          !normal = normal+SurfaceNormal( etype, nn, cnpos, elem )
          normal_n = SurfaceNormal( etype, nn, cnpos, elem )
          normal = normal+normal_n
          cnt = cnt+1
        endif
      enddo
      !normal = normal/cnt                                        !!-???
    elseif( 12 <= isin .and. isin <= 41 ) then  ! edge
      cnode1 = isin / 10
      cnode2 = mod(isin, 10)
      gn1 = surf(csurf)%nodes(cnode1)
      gn2 = surf(csurf)%nodes(cnode2)
      etype = surf(csurf)%etype
      nn = size( surf(csurf)%nodes )
      do j=1,nn
        iSS = surf(csurf)%nodes(j)
        elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      enddo
      normal = SurfaceNormal( etype, nn, lpos, elem )
      select case (etype)
      case (fe_tri3n, fe_tri6n, fe_tri6nc)
        if    ( isin==12 ) then; x=lpos(2)-lpos(1)
        elseif( isin==23 ) then; x=1.d0-2.d0*lpos(2)
        elseif( isin==31 ) then; x=2.d0*lpos(1)-1.d0
        else; stop "Error: cal_node_normal: invalid isin"
        endif
      case (fe_quad4n, fe_quad8n)
        if    ( isin==12 ) then; x=lpos(1)
        elseif( isin==23 ) then; x=lpos(2)
        elseif( isin==34 ) then; x=-lpos(1)
        elseif( isin==41 ) then; x=-lpos(2)
        else; stop "Error: cal_node_normal: invalid isin"
        endif
      end select
      ! find neighbor surf that includes cnode1 and cnode2
      nsurf = 0
      NEIB_LOOP: do i=1, surf(csurf)%n_neighbor
        nd1 = surf(csurf)%neighbor(i)
        nn = size( surf(nd1)%nodes )
        etype = surf(nd1)%etype
        cgn1 = 0
        cgn2 = 0
        do j=1,nn
          iSS = surf(nd1)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          if( iSS==gn1 ) cgn1=j
          if( iSS==gn2 ) cgn2=j
        enddo
        if( cgn1>0 .and. cgn2>0 ) then
          nsurf = nd1
          isin_n = 10*cgn2 + cgn1
          x = -x
          select case (etype)
          case (fe_tri3n, fe_tri6n, fe_tri6nc)
            if    ( isin_n==12 ) then; lpos_n(1)=0.5d0*(1.d0-x); lpos_n(2)=0.5d0*(1.d0+x)
            elseif( isin_n==23 ) then; lpos_n(1)=0.d0;           lpos_n(2)=0.5d0*(1.d0-x)
            elseif( isin_n==31 ) then; lpos_n(1)=0.5d0*(1.d0+x); lpos_n(2)=0.d0
            else; stop "Error: cal_node_normal: invalid isin_n"
            endif
          case (fe_quad4n, fe_quad8n)
            if    ( isin_n==12 ) then; lpos_n(1)= x;    lpos_n(2)=-1.d0
            elseif( isin_n==23 ) then; lpos_n(1)= 1.d0; lpos_n(2)= x
            elseif( isin_n==34 ) then; lpos_n(1)=-x;    lpos_n(2)= 1.d0
            elseif( isin_n==41 ) then; lpos_n(1)=-1.d0; lpos_n(2)=-x
            else; stop "Error: cal_node_normal: invalid isin_n"
            endif
          end select
          !normal = normal + SurfaceNormal( etype, nn, lpos_n, elem )
          normal_n = SurfaceNormal( etype, nn, lpos_n, elem )
          normal = normal+normal_n
          exit NEIB_LOOP
        endif
      enddo NEIB_LOOP
      !if( nsurf==0 ) write(0,*) "Warning: cal_node_normal: neighbor surf not found"
      !normal = normal/2
    endif
    normal = normal/ dsqrt( dot_product( normal, normal ) )
  end subroutine cal_node_normal

  !> This subroutine tracks down next contact position after a finite slide
  subroutine track_contact_position( flag_ctAlgo, nslave, contact, currpos, currdisp, mu, infoCTChange, nodeID, elemID, B )
    character(len=9), intent(in)                    :: flag_ctAlgo  !< contact analysis algorithm flag
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)    !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: mu           !< penalty
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    real(kind=kreal), intent(inout)                  :: B(:)         !< nodal force residual

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node), elem0(3, l_max_elem_node)
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos(1:2), contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( .not. isin ) then   ! such case is considered to rarely or never occur
      write(*,*) 'Warning: contact moved beyond neighbor elements'
      ! get master candidates from bucketDB
      bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
      nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
      if (nCand > 0) then
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)
        do idm= 1, nCand
          sid = indexCand(idm)
          if( sid==sid0 ) cycle
          if( associated(contact%master(sid0)%neighbor) ) then
            if( any(sid==contact%master(sid0)%neighbor(:)) ) cycle
          endif
          if (.not. is_in_surface_box( contact%master(sid), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          call project_Point2SurfElement( coord, contact%master(sid), currpos, &
            contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
            localclr=contact%cparam%CLEARANCE )
          if( isin ) then
            contact%states(nslave)%surface = sid
            exit
          endif
        enddo
        deallocate(indexCand)
      endif
    endif

    if( isin ) then
      if( contact%states(nslave)%surface==sid0 ) then
        if(any(dabs(contact%states(nslave)%lpos(1:2)-opos(:)) >= contact%cparam%CLR_DIFFLPOS))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Node",nodeID(slave)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          contact%states(nslave)%distance," at ",contact%states(nslave)%lpos(1:2)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
        if( flag_ctAlgo=='ALagrange' )  &
          call reset_contact_force( contact, currpos, nslave, sid0, opos, odirec, B )
      endif
      if( flag_ctAlgo=='SLagrange' ) then
        ! Setup elem array for update_TangentForce
        etype = contact%master(contact%states(nslave)%surface)%etype
        nn = size(contact%master(contact%states(nslave)%surface)%nodes)
        do j=1,nn
          iSS = contact%master(contact%states(nslave)%surface)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call update_TangentForce(etype,nn,elem0,elem,contact%states(nslave))
      endif
      iSS = isInsideElement( etype, contact%states(nslave)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( contact%states(nslave)%surface, iSS, contact%master, currpos, &
          contact%states(nslave)%lpos(1:2), contact%states(nslave)%direction(:) )
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      contact%states(nslave)%state = CONTACTFREE
      contact%states(nslave)%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position

  !> This subroutine tracks down next contact position after a finite slide
  subroutine update_direction( nslave, contact, currpos )
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes

    integer(kind=kint) :: slave, sid0, etype, iSS
    real(kind=kreal)    :: coord(3)
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    type(tContactState) :: cstate_tmp  !< Recorde of contact information

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction

    cstate_tmp%state = contact%states(nslave)%state
    cstate_tmp%surface = sid0
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      cstate_tmp, isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos, contact%cparam%CLR_SAME_ELEM )

    if( isin ) then
      etype = contact%master(sid0)%etype
      iSS = isInsideElement( etype, cstate_tmp%lpos, contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( cstate_tmp%surface, iSS, contact%master, currpos, &
        cstate_tmp%lpos, cstate_tmp%direction(:) )
    endif

    contact%states(nslave)%direction = cstate_tmp%direction

  end subroutine 


  !>\brief This subroutine update contact force in case that contacting element
  !> is changed
  subroutine reset_contact_force( contact, currpos, lslave, omaster, opos, odirec, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: currpos(:)     !< current coordinate of each nodes
    integer, intent(in)               :: lslave         !< slave node
    integer, intent(in)               :: omaster        !< former master element
    real(kind=kreal), intent(in)      :: opos(2)        !< former contact pos
    real(kind=kreal), intent(in)      :: odirec(3)      !< former contact direction
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, j, iSS
    real(kind=kreal)    :: nrlforce, fcoeff, tangent(3,2)
    real(kind=kreal)    :: elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)
    integer(kind=kint)  :: i, idx0

    slave = contact%slave(lslave)
    fcoeff = contact%fcoeff
    ! clear contact force in former contacted element
    nrlforce = contact%states(lslave)%multiplier(1)
    B(3*slave-2:3*slave) = B(3*slave-2:3*slave)+nrlforce*odirec
    nn = size( contact%master(omaster)%nodes )
    etype = contact%master(omaster)%etype
    call getShapeFunc( etype, opos(:), shapefunc )
    do j=1,nn
      iSS = contact%master(omaster)%nodes(j)
      ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-nrlforce*shapefunc(j)*odirec
      idx0 = 3*(iSS-1)
      do i=1,3
        !$omp atomic
        B(idx0+i) = B(idx0+i)-nrlforce*shapefunc(j)*odirec(i)
      enddo
    enddo
    if( fcoeff/=0.d0 ) then
      do j=1,nn
        iSS = contact%master(omaster)%nodes(j)
        elemcrd(:,j) = currpos(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( opos(:), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      fric(1:2) = contact%states(lslave)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)+f3(1:3)
      do j=1,nn
        iSS = contact%master(omaster)%nodes(j)
        ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+f3(3*j+1:3*j+3)
        idx0 = 3*(iSS-1)
        do i=1,3
          !$omp atomic
          B(idx0+i) = B(idx0+i)+f3(3*j+i)
        enddo
      enddo
    endif

    ! reset contact force in new contacting element
    master = contact%states(lslave)%surface
    nn = size( contact%master(master)%nodes )
    etype = contact%master(master)%etype
    call getShapeFunc( etype, contact%states(lslave)%lpos(1:2), shapefunc )
    B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(lslave)%direction
    do j=1,nn
      iSS = contact%master(master)%nodes(j)
      ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
        !            shapefunc(j)*contact%states(lslave)%direction
      idx0 = 3*(iSS-1)
      do i=1,3
        !$omp atomic
        B(idx0+i) = B(idx0+i)+nrlforce*        &
          shapefunc(j)*contact%states(lslave)%direction(i)
      enddo
    enddo
    if( fcoeff/=0.d0 ) then
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemcrd(:,j) = currpos(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( contact%states(lslave)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      fric(1:2) = contact%states(lslave)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
        idx0 = 3*(iSS-1)
        do i=1,3
          !$omp atomic
          B(idx0+i) = B(idx0+i)-f3(3*j+i)
        enddo
      enddo
    endif

  end subroutine reset_contact_force

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_contact_multiplier( contact, coord, disp, ddisp, fcoeff, mu, mut,   &
      gnt, ctchanged )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: fcoeff         !< frictional coeff
    real(kind=kreal), intent(in)      :: mu, mut        !< penalty
    real(kind=kreal), intent(out)     :: gnt(2)         !< convergency information
    logical, intent(inout)            :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS, cnt
    real(kind=kreal)    :: elemdisp(3,l_max_elem_node), tangent(3,2)
    real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: lgnt(2), fric(2), f3(3*l_max_elem_node+3), edisp(3*l_max_elem_node+3)

    cnt =0; lgnt(:)=0.d0
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
      contact%states(i)%wkdist = -dgn
      contact%states(i)%multiplier(1) = contact%states(i)%multiplier(1) - mu*contact%states(i)%wkdist
      contact%states(i)%distance = contact%states(i)%wkdist
      cnt = cnt+1
      lgnt(1) = lgnt(1)- contact%states(i)%wkdist

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      dxi(1) = dot_product( dispmat(1,1:nn*3+3), edisp(1:nn*3+3) )
      dxi(2) = dot_product( dispmat(2,1:nn*3+3), edisp(1:nn*3+3) )
      dxy(:) = matmul( metric, dxi )
      fric(1:2) = contact%states(i)%multiplier(2:3) + mut*dxy(1:2)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      dgn = dsqrt( f3(1)*f3(1)+f3(2)*f3(2)+f3(3)*f3(3) )
      if( contact%states(i)%multiplier(1)>0.d0 ) then
        if(  dgn > fcoeff*contact%states(i)%multiplier(1) ) then
          if( contact%states(i)%state==CONTACTSTICK ) then
            ctchanged= .true.
            print *, "Node", slave, "to slip state",dgn, fcoeff*contact%states(i)%multiplier(1)
          endif
          contact%states(i)%state = CONTACTSLIP
          fric(:) = fric(:)*fcoeff*contact%states(i)%multiplier(1)/dgn
        else
          if( contact%states(i)%state==CONTACTSLIP ) then
            ctchanged= .true.
            print *, "Node", slave, "to stick state",dgn, fcoeff*contact%states(i)%multiplier(1)
          endif
          contact%states(i)%state = CONTACTSTICK
        endif
      endif
      contact%states(i)%multiplier(2:3) = fric(:)
      dxy(:) = matmul( dg, tangent )
      lgnt(2) = lgnt(2)+dsqrt( dxy(1)*dxy(1)+dxy(2)*dxy(2) )
    enddo
    if(cnt>0) lgnt(:) = lgnt(:)/cnt
    gnt = gnt + lgnt
  end subroutine update_contact_multiplier

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_tied_multiplier( contact, disp, ddisp, mu, ctchanged )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: mu             !< penalty
    logical, intent(inout)            :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave, etype, master
    integer(kind=kint)  :: nn, i, j, iSS, cnt
    real(kind=kreal)    :: dg(3), dgmax
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

      contact%states(i)%multiplier(1:3) = contact%states(i)%multiplier(1:3) + mu*dg(1:3)

      ! check if tied constraint converged
      dgmax = 0.d0
      do j=1,(nn+1)*3
        dgmax = dgmax + dabs(edisp(j))
      enddo
      dgmax = dgmax/dble((nn+1)*3)
      do j=1,3
        if( dabs(dg(j))/dmax1(1.d0,dgmax) > 1.d-3 ) ctchanged = .true.
      enddo

    enddo
  end subroutine 

  !>\brief This subroutine assemble contact force into contacing nodes
  subroutine ass_contact_force( contact, coord, disp, B )
    type( tContact ), intent(in)      :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current now
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: fcoeff, nrlforce, tangent(3,2)
    real(kind=kreal)    :: elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)
    fcoeff = contact%fcoeff
    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      nrlforce = contact%states(i)%multiplier(1)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(i)%direction

      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
          shapefunc(j)*contact%states(i)%direction
      enddo

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemcrd(:,j) = coord(3*iSS-2:3*iSS)+disp(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )

      fric(1:2) = contact%states(i)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
      enddo
    enddo

  end subroutine ass_contact_force

  !>\brief This subroutine setup contact output nodal vectors
  subroutine set_contact_state_vector( contact, dt, relvel_vec, state_vec )
      type( tContact ), intent(in)      :: contact        !< contact info
      real(kind=kreal), intent(in)      :: dt
      real(kind=kreal), intent(inout)   :: relvel_vec(:)       !< mesh coordinate
      real(kind=kreal), intent(inout)   :: state_vec(:)        !< disp till current now

      integer(kind=kint)  :: i, slave

      do i= 1, size(contact%slave)
        slave = contact%slave(i)
        if( state_vec(slave) < 0.1d0 .or. contact%states(i)%state > 0 ) &
          &  state_vec(slave) = dble(contact%states(i)%state)

        if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
        if( dt < 1.d-16 ) cycle ! too small delta t
        relvel_vec(3*slave-2:3*slave) = contact%states(i)%reldisp(1:3)/dt
      enddo

  end subroutine set_contact_state_vector

  subroutine update_contact_TangentForce( contact )
    type( tContact ), intent(inout)   :: contact        !< contact info

    integer(kind=kint)  :: i

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) then
        contact%states(i)%tangentForce(1:3) = 0.d0
        contact%states(i)%tangentForce_trial(1:3) = 0.d0
        contact%states(i)%tangentForce_final(1:3) = 0.d0
      else
        contact%states(i)%tangentForce(1:3) = contact%states(i)%tangentForce_final(1:3)
      end if
      contact%states(i)%tangentForce1(1:3) = contact%states(i)%tangentForce(1:3)
    enddo
  end subroutine update_contact_TangentForce

  !> This subroutine tracks down next contact position after a finite slide
  subroutine track_contact_position_exp( nslave, contact, currpos, currdisp, infoCTChange, nodeID, elemID )
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                 :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)  :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                    :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                    :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos(1:2), contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( .not. isin ) then   ! such case is considered to rarely or never occur
      write(*,*) 'Warning: contact moved beyond neighbor elements'
      ! get master candidates from bucketDB
      bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
      nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
      if (nCand > 0) then
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)
        do idm= 1, nCand
          sid = indexCand(idm)
          if( sid==sid0 ) cycle
          if( any(sid==contact%master(sid0)%neighbor(:)) ) cycle
          if (.not. is_in_surface_box( contact%master(sid), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          call project_Point2SurfElement( coord, contact%master(sid), currpos, &
            contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
            localclr=contact%cparam%CLEARANCE )
          if( isin ) then
            contact%states(nslave)%surface = sid
            exit
          endif
        enddo
        deallocate(indexCand)
      endif
    endif

    if( isin ) then
      if( contact%states(nslave)%surface==sid0 ) then
        if(any(dabs(contact%states(nslave)%lpos(1:2)-opos(:)) >= contact%cparam%CLR_DIFFLPOS))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Node",nodeID(slave)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          contact%states(nslave)%distance," at ",contact%states(nslave)%lpos(1:2)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
      endif
      iSS = isInsideElement( etype, contact%states(nslave)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) then
        call cal_node_normal( contact%states(nslave)%surface, iSS, contact%master, currpos, &
          contact%states(nslave)%lpos(1:2), contact%states(nslave)%direction(:) )
      endif
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      contact%states(nslave)%state = CONTACTFREE
      contact%states(nslave)%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position_exp

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_contact_state_exp( contact, currpos, currdisp, infoCTChange, &
      nodeID, elemID, is_init, active )
    type( tContact ), intent(inout)                 :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)  :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                    :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                    :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                             :: is_init      !< whether initial scan or not
    logical, intent(out)                            :: active       !< if any in contact

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3)
    real(kind=kreal)    :: nlforce
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_cand

    if( is_init ) then
      distclr = contact%cparam%DISTCLR_INIT
    else
      distclr = contact%cparam%DISTCLR_FREE
      if( contact%algtype == CONTACTTIED ) then
        active = .false.
        do i= 1, size(contact%slave)
          if( contact%states(i)%state==CONTACTSTICK ) then
            active = .true.
            exit
          endif
        enddo
        return
      endif
    endif

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(contact%slave)))
    contact_surf(:) = size(elemID)+1
    do i = 1, size(contact%slave)
      states_prev(i) = contact%states(i)%state
    enddo

    call update_surface_box_info( contact%master, currpos )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,is_cand,idm,etype,isin, &
      !$omp&         bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax) &
      !$omp& shared(contact,infoCTChange,currpos,currdisp,nodeID,elemID,distclr,contact_surf) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
        call track_contact_position_exp( i, contact, currpos, currdisp, infoCTChange, nodeID, elemID )
        if( contact%states(i)%state /= CONTACTFREE ) then
          contact_surf(contact%slave(i)) = -contact%states(i)%surface
        endif
      else if( contact%states(i)%state==CONTACTFREE ) then
        coord(:) = currpos(3*slave-2:3*slave)

        ! get master candidates from bucketDB
        bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
        nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
        if (nCand == 0) cycle
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)

        nMasterMax = nCand
        allocate(indexMaster(nMasterMax))
        nMaster = 0

        ! narrow down candidates
        do idm= 1, nCand
          id = indexCand(idm)
          if (.not. is_in_surface_box( contact%master(id), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          nMaster = nMaster + 1
          indexMaster(nMaster) = id
        enddo
        deallocate(indexCand)

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif

        do idm = 1,nMaster
          id = indexMaster(idm)
          call project_Point2SurfElement( coord, contact%master(id), currpos, &
            contact%states(i), isin, distclr, localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          etype = contact%master(id)%etype
          iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
          if( iSS>0 ) &
            call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%lpos(1:2), &
              contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = id
          write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3)') "Node",nodeID(slave)," contact with element", &
            elemID(contact%master(id)%eid),       &
            " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(1:2), &
            " along direction ", contact%states(i)%direction
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    call hecmw_contact_comm_allreduce_i(contact%comm, contact_surf, HECMW_MIN)
    nactive = 0
    do i = 1, size(contact%slave)
      if (contact%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
        if (abs(contact_surf(contact%slave(i))) /= contact%states(i)%surface) then ! that is in contact with other surface
          contact%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i6,A,i6,A)') "Node",nodeID(contact%slave(i)), &
            " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else
          nactive = nactive + 1
        endif
      endif
      if (states_prev(i) == CONTACTFREE .and. contact%states(i)%state /= CONTACTFREE) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (states_prev(i) /= CONTACTFREE .and. contact%states(i)%state == CONTACTFREE) then
        infoCTChange%contact2free = infoCTChange%contact2free + 1
      endif
    enddo
    active = (nactive > 0)
    deallocate(contact_surf)
    deallocate(states_prev)
  end subroutine scan_contact_state_exp

end module m_fstr_contact_lib
