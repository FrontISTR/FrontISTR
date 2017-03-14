!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provide functions of contact stiffness calculation
module m_contact_lib
use elementInfo
implicit none

integer, parameter, private    :: kreal = kind(0.0d0)
integer, parameter, private    :: l_max_surface_node =20
integer, parameter, private    :: l_max_elem_node = 100

integer, parameter :: CONTACTUNKNOWN = -1
!> contact state definition
integer, parameter :: CONTACTFREE = -1
integer, parameter :: CONTACTSTICK = 1
integer, parameter :: CONTACTSLIP = 2

!> contact type or algorithm definition
integer, parameter :: CONTACTTIED = 1
integer, parameter :: CONTACTGLUED = 2
integer, parameter :: CONTACTSSLID = 3
integer, parameter :: CONTACTFSLID = 4

      !> This structure records contact status
      type tContactState
        integer             :: state        !< -1:free, 1:in contact, or other needed
        integer             :: surface      !< contacting surface number
        real(kind=kreal)    :: distance      !< penetration value
        real(kind=kreal)    :: wkdist        !< copy of penetration value
        real(kind=kreal)    :: lpos(2)       !< contact position(local coordinate)
        real(kind=kreal)    :: gpos(3)       !< contact position(global coordinate)
        real(kind=kreal)    :: direction(3)  !< contact direction
        real(kind=kreal)    :: multiplier(3) !< Lagrangian multiplier or contact force
		                                     !< 1: normal 2:tangent component
	    real(kind=kreal)    :: tangentForce(3)                       !< friction force
	    real(kind=kreal)    :: tangentForce_trial(3)                 !< trial friction force
	    real(kind=kreal)    :: tangentForce_final(3)                 !< final friction force
      end type

contains

  !> Initializer
  subroutine contact_state_init(cstate)
    type(tContactState), intent(inout) :: cstate   !< contact state
    cstate%state=-1
    cstate%surface=-1
  end subroutine

  !> Print out contact state
  subroutine print_contact_state(fnum, cstate)
    integer, intent(in)             :: fnum        !< file number
    type(tContactState), intent(in) :: cstate      !< contact state
    write(fnum, *) "--Contact state=",cstate%state
    write(fnum, *) cstate%surface, cstate%distance
    write(fnum, *) cstate%lpos
    write(fnum, *) cstate%direction
    write(fnum, *) cstate%multiplier
  end subroutine

  !> Transfer contact condition int mpc bundary conditions
  subroutine contact2mpcval( cstate, etype, nnode, mpcval )
    type(tContactState), intent(in) :: cstate              !< contact state
    integer, intent(in)             :: etype               !< type of contacting surface
    integer, intent(in)             :: nnode               !< number of elemental nodes
    real(kind=kreal), intent(out)   :: mpcval(nnode*3+4)   !< MPC constraint

    integer          :: i,j
    real(kind=kreal) :: shapefunc(nnode)

    call getShapeFunc( etype, cstate%lpos(:), shapefunc )
    mpcval(1:3) = cstate%direction(1:3)
    do i=1,nnode
    do j=1,3
      mpcval( i*3+j ) = -cstate%direction(j)*shapefunc(i)
    enddo
    enddo
    mpcval( 3*nnode+4 )=cstate%distance
  end subroutine

  !> This subroutine calculate contact stiff matrix and contact force
  subroutine contact2stiff( flag, cstate, etype, nnode, ele, mu, mut,  &
                         fcoeff, symm, stiff, force )
    integer, intent(in)             :: flag            !< small slid or finite slide
    type(tContactState), intent(in) :: cstate          !< contact state
    integer, intent(in)             :: etype           !< type of contacting surface
    integer, intent(in)             :: nnode           !< number of elemental nodes
    real(kind=kreal), intent(in)    :: ele(3,nnode)    !< coord of surface element
    real(kind=kreal), intent(in)    :: mu              !< penalty
    real(kind=kreal), intent(in)    :: mut             !< penalty along tangent
    real(kind=kreal), intent(in)    :: fcoeff          !< friction coefficient
    logical, intent(in)             :: symm            !< symmtricalize
    real(kind=kreal), intent(out)   :: stiff(:,:)      !< contact stiffness
    real(kind=kreal), intent(out)   :: force(:)        !< contact force

    integer          :: i,j
    real(kind=kreal) :: shapefunc(nnode), shapederiv(nnode,2)
    real(kind=kreal) :: N(nnode*3+3), NN(nnode*3+3,2), dispmat(2,nnode*3+3)
    real(kind=kreal) :: curv(2,2), metric(2,2), l2ndderiv(3,2,2)
    real(kind=kreal) :: det, inverse(2,2), ff(2), cff(2)
    real(kind=kreal) :: dum11(nnode*3+3,nnode*3+3), dum12(nnode*3+3,nnode*3+3)
    real(kind=kreal) :: dum21(nnode*3+3,nnode*3+3), dum22(nnode*3+3,nnode*3+3)
    real(kind=kreal) :: fstiff(nnode*3+3,nnode*3+3), tangent(3,2)

    call getShapeFunc( etype, cstate%lpos(:), shapefunc )
    N(1:3) = cstate%direction(1:3)
    do i=1,nnode
      N( i*3+1:i*3+3 ) = -shapefunc(i)*cstate%direction(1:3)
    enddo
    forall( i=1:nnode*3+3, j=1:nnode*3+3 )
      stiff(i,j) = mu* N(i)*N(j)
    end forall
    force(1:nnode*3+3) = N(:)

    if( fcoeff/=0.d0 .or. flag==CONTACTFSLID ) &
      call DispIncreMatrix( cstate%lpos, etype, nnode, ele, tangent, metric, dispmat )

    ! frictional component
    if( fcoeff/=0.d0 ) then
      forall(i=1:nnode*3+3, j=1:nnode*3+3)
        dum11(i,j) = mut*dispmat(1,i)*dispmat(1,j)
        dum12(i,j) = mut*dispmat(1,i)*dispmat(2,j)
        dum21(i,j) = mut*dispmat(2,i)*dispmat(1,j)
        dum22(i,j) = mut*dispmat(2,i)*dispmat(2,j)
      end forall
      stiff(1:nnode*3+3,1:nnode*3+3)                   &
           = stiff(1:nnode*3+3,1:nnode*3+3)            &
           + metric(1,1)*dum11 + metric(1,2)*dum12     &
           + metric(2,1)*dum21 + metric(2,2)*dum22

      if( cstate%state == CONTACTSLIP ) then
        det = metric(1,1)*metric(2,2)-metric(1,2)*metric(2,1)
        if( det==0.d0 ) stop "Math error in contact stiff calculation"
        inverse(1,1) = metric(2,2)/det
        inverse(2,1) = -metric(2,1)/det
        inverse(1,2) = -metric(1,2)/det
        inverse(2,2) = metric(1,1)/det
        ff(:) = cstate%multiplier(2:3)
        cff(:) = matmul( inverse, ff )
        ff(:) = ff(:)/dsqrt( ff(1)*ff(1)+ff(2)*ff(2) )
        cff(:) = cff(:)/dsqrt( cff(1)*cff(1)+cff(2)*cff(2) )
        stiff(1:nnode*3+3,1:nnode*3+3) = stiff(1:nnode*3+3,1:nnode*3+3) -         &
          ( cff(1)*ff(1)*metric(1,1)+ cff(2)*ff(1)*metric(1,2) )*dum11 -   &
          ( cff(2)*ff(2)*metric(1,2)+ cff(1)*ff(2)*metric(1,1) )*dum21 -   &
          ( cff(1)*ff(1)*metric(1,2)+ cff(2)*ff(1)*metric(2,2) )*dum12 -   &
          ( cff(2)*ff(2)*metric(2,2)+ cff(1)*ff(2)*metric(1,2) )*dum22
      endif
    endif

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
    real(kind=kreal) :: det, inverse(2,2)
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
  subroutine project_Point2Element(xyz,etype,nn,elemt,cstate,isin,distclr,ctpos,localclr)
     real(kind=kreal),intent(in)       :: xyz(3)        !< Coordinates of a spacial point, whose projecting point is to be computed
     integer, intent(in)               :: etype         !< surface element type
     integer, intent(in)               :: nn            !< number of elemental nodes
     real(kind=kreal),intent(in)       :: elemt(3,nn)   !< nodes coordinates of surface element
     type(tContactState), intent(out)  :: cstate        !< Recorde of contact information
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
       xyz_out = matmul( elemt(:,:), sfunc )
       dxyz(1:3) = xyz_out(1:3) - xyz(1:3)
       dist_last = dot_product( dxyz, dxyz(:) )

       call TangentBase( etype, nn, r, elemt, tangent )
       call Curvature( etype, nn, r, elemt, curv )

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
            xyz_out(1:3) = matmul( elemt(:,:), sfunc(:) )
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
       normal(:) = SurfaceNormal( etype, nn, r, elemt )
       normal(:) = normal(:)/dsqrt( dot_product(normal, normal) )
       do count = 1,3
         if( dabs(normal(count))<1.D-10 ) normal(count) =0.d0
         if( dabs(1.d0-dabs(normal(count)))<1.D-10 ) normal(count) =sign(1.d0, normal(count))
       enddo
       cstate%distance = dot_product( dxyz, normal )

       if( cstate%distance < distclr .and. cstate%distance > -5.0d-01 ) isin = .true.

       if( isin ) then
         if( initstate== CONTACTFREE ) then
            cstate%state = CONTACTSTICK
         else
            cstate%state = initstate
         endif
         cstate%gpos(:)=xyz_out(:)
         cstate%lpos(:)=r(:)
         cstate%direction(:) = normal(:)
         cstate%wkdist = cstate%distance
       endif
     endif
  end subroutine project_Point2Element

end module m_contact_lib


