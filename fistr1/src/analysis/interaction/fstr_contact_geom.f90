!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides geometric calculations for contact
module m_fstr_contact_geom
  use hecmw
  use elementInfo
  use mContactDef
  use m_fstr_contact_smoothing
  implicit none

  public :: cal_node_normal

contains

  ! TODO: Move from fstr_contact_lib.f90
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
    real(kind=kreal)  ::  inverse(3,3)
    real(kind=kreal)  ::  sfunc(nn), deriv(nn,3)
    real(kind=kreal)  ::  r(3), dr(3), r_tmp(3)        ! natural coordinate
    real(kind=kreal)  ::  xyz_out(3)                   ! curr. projection position
    real(kind=kreal)  ::  dist_last,dist_now, dxyz(3)  ! dist between the point and its projection
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

  !> Wrapper for project_Point2Element that takes tSurfElement structure
  !! This subroutine handles element coordinate extraction from tSurfElement
  subroutine project_Point2SurfElement(xyz, surf, currpos, cstate, isin, distclr, ctpos, localclr, smoothing)
    real(kind=kreal), intent(in)              :: xyz(3)        !< coordinates of slave point
    type(tSurfElement), intent(in)            :: surf          !< surface element structure
    real(kind=kreal), intent(in)              :: currpos(:)    !< current coordinate of all nodes
    type(tContactState), intent(inout)        :: cstate        !< contact state
    logical, intent(out)                      :: isin          !< in contact or not
    real(kind=kreal), intent(in)              :: distclr       !< clearance of contact distance
    real(kind=kreal), optional, intent(in)    :: ctpos(2)      !< current contact position (natural coord)
    real(kind=kreal), optional, intent(in)    :: localclr      !< clearance of contact local coord
    integer(kind=kint), intent(in)            :: smoothing     !< kcsNONE or kcsNAGATA

    integer(kind=kint) :: nn, j, iSS, etype_use, nn_use
    real(kind=kreal)   :: elem(3, l_max_elem_node)
    real(kind=kreal)   :: elem_tri3n(3, 3)

    ! Extract element coordinates from surf structure
    nn = size(surf%nodes)
    do j = 1, nn
      iSS = surf%nodes(j)
      elem(1:3, j) = currpos(3*iSS-2:3*iSS)
    enddo

    if (smoothing == kcsNAGATA) then
      if (surf%etype == fe_tri3n) then
        elem_tri3n = elem(1:3, 1:3)
        call reorder_tri3n_to_tri6n(elem_tri3n, surf%intermediate_points, elem(1:3,1:6))
        etype_use = fe_tri6n
        nn_use = 6
      else if (surf%etype == fe_quad4n) then
        do j = 1, nn
          elem(1:3, nn+j) = surf%intermediate_points(1:3, j)
        enddo
        etype_use = fe_quad8n
        nn_use = 8
      else
        etype_use = surf%etype
        nn_use = nn
      endif
    else
      etype_use = surf%etype
      nn_use = nn
    endif

    call project_Point2Element(xyz, etype_use, nn_use, elem, surf%reflen, cstate, &
      isin, distclr, ctpos, localclr)

  end subroutine project_Point2SurfElement

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
      contact%states(nslave)%lpos, contact%cparam%CLR_SAME_ELEM, &
      smoothing=contact%smoothing )

    if( isin ) then
      etype = contact%master(sid0)%etype
      iSS = isInsideElement( etype, cstate_tmp%lpos, contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( cstate_tmp%surface, iSS, contact%master, currpos, &
        cstate_tmp%lpos, cstate_tmp%direction(:) )
    endif

    contact%states(nslave)%direction = cstate_tmp%direction

  end subroutine

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

end module m_fstr_contact_geom
