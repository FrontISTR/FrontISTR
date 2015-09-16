!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.6                                   !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Dynamic Transit Analysis                           !
!                                                                      !
!                    Written by Noburu Imai (Univ. of Tokyo)           !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module contains functions relates to coupling analysis
module m_dynamic_mat_ass_couple
use m_fstr
private :: area_of_triangle
private :: area_of_triangle2
private :: area_of_squre
private :: node_on_surface

integer(kind=kint), private :: icall = 0
contains

  subroutine DYNAMIC_MAT_ASS_COUPLE( hecMESH, hecMAT, fstrSOLID, fstrCPL )
    use m_fstr
    implicit none
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid        ) :: fstrSOLID
    type (fstr_couple       ) :: fstrCPL
    ! local
    integer( kind=kint ) :: ig0,ig,is,ie,ik
    integer( kind=kint ) :: i,j,count
    integer( kind=kint ) :: eid, sid, etype
    integer( kind=kint ) :: node(20)
    integer( kind=kint ) :: node_n
    real(kind=kreal) :: px, py, pz ! traction on surface
    real(kind=kreal) :: vx, vy, vz ! force on vertex

! ================== modified by K. Tagami ============= 2010/02/25 ====
!!!!      real(kind=kreal) :: xx(4), yy(4), zz(4)
    integer :: MaxNumNodesOnSurf
    Parameter( MaxNumNodesOnSurf = 8 )
    real(kind=kreal) :: xx( MaxNumNodesOnSurf )
    real(kind=kreal) :: yy( MaxNumNodesOnSurf )
    real(kind=kreal) :: zz( MaxNumNodesOnSurf )
! ================================================================

    real(kind=kreal) :: area, wg
    integer( kind=kint ) :: ierr

! ============= added by K. Tagami =========for debug ====== 2010/03/02 ====
    integer( kind=kint ) :: node_global
! ==========================================================

    icall = icall + 1
    ierr = 0

    do ig0= 1, fstrSOLID%COUPLE_ngrp_tot
       ig= fstrSOLID%COUPLE_ngrp_ID(ig0)
       is= hecMESH%surf_group%grp_index(ig-1) + 1
       ie= hecMESH%surf_group%grp_index(ig  )
       do ik= is, ie
          sid   = hecMESH%surf_group%grp_item(2*ik)
          eid   = hecMESH%surf_group%grp_item(2*ik-1)
          etype = hecMESH%elem_type(eid)
          call node_on_surface( hecMESH, etype, eid, sid, node, node_n )

         !-------------------------------------------------------------
          ! traction on surface
          !-------------------------------------------------------------
          px = 0
          py = 0
          pz = 0
          count=0
          do i=1, node_n
             j=3*fstrCPL%index( node(i) )

! ============= added by K. Tagami =========for debug ====== 2010/03/02 ====
!             node_global = hecMESH%global_node_ID( node(i) )
! ======================================================================
!
             if( icall == 1 ) then
                if( j<=0 ) then
                           write(IDBG,'(a,i0,a)') "dynamic_mat_ass_couple: traction for node ", &
                                &   hecMESH%global_node_ID(node(i)), " on coupling surface not found"
! ============ added by K. Tagami ============== for debug ==== 2010/03/02
!                           write(IDBG,*) "local : ", node(i), " Global :", node_global, " not found"
! ==============================================================
                   ierr = 1
                   cycle
                else
                   write(IDBG,'(a,i0,a)') "dynamic_mat_ass_couple: traction for node ", &
                        & hecMESH%global_node_ID(node(i)), " on coupling surface OK"
! ============ added by K. Tagami ============== for debug ==== 2010/03/02
!                           write(IDBG,*) "local : ", node(i), " Global :", node_global, " OK"
! ==============================================================
                endif
             endif
             px = px + fstrCPL%trac(j-2)
             py = py + fstrCPL%trac(j-1)
             pz = pz + fstrCPL%trac(j  )
             count = count +1
          end do
          if( count == 0 ) cycle
          px = px / count
          py = py / count
          pz = pz / count
          !-------------------------------------------------------------
          ! force on vertex
          !-------------------------------------------------------------
          do i=1, node_n
             j=3*node(i)
             xx(i)=hecMESH%node(j-2)
             yy(i)=hecMESH%node(j-1)
             zz(i)=hecMESH%node(j  )
          end do
          if( node_n == 3 ) then
             area = area_of_triangle( xx,yy,zz )
          else if( node_n == 4 ) then
             area = area_of_squre( xx,yy,zz )
          else if( node_n == 6 ) then
             area = area_of_triangle2(xx,yy,zz)
          else
             write(*,*) "#Error : in FSTR_MAT_ASS_COUPLE "
             call hecmw_abort( hecmw_comm_get_comm())
          end if
          wg = area / node_n

          vx = px * wg
          vy = py * wg
          vz = pz * wg
          !-------------------------------------------------------------
          ! add in B
          !-------------------------------------------------------------
          do i=1, node_n
             j=3*node(i)

             hecMAT%B(j-2)=hecMAT%B(j-2) + vx
             hecMAT%B(j-1)=hecMAT%B(j-1) + vy
             hecMAT%B(j  )=hecMAT%B(j  ) + vz
          end do
       end do
    end do
    call hecmw_barrier(hecMESH)
    if( ierr == 1 ) then
       write(*,*) "#Error : in FSTR_MAT_ASS_COUPLE"
       call hecmw_abort( hecmw_comm_get_comm())
    endif
  end subroutine DYNAMIC_MAT_ASS_COUPLE

  subroutine fstr_get_trac(inpX,inpY,inpVx,inpVy,inpP,idx,idy,cx,cy,cz,normal,trac)
      use elementInfo
      implicit none

      integer(kind = kint) :: i, j, k, idx, idy
      integer(kind = kint) :: ix1, ix2, iy1, iy2
      real(kind = kreal)   :: inpX(:), inpY(:)
      real(kind = kreal)   :: trac(3), normal(3)
      real(kind = kreal)   :: cx, cy, cz, p
      real(kind = kreal)   :: x1, x2, y1, y2
      real(kind = kreal)   :: p11, p12, p21, p22
      real(kind = kreal)   :: x, r, theta
      real(kind = kreal)   :: xx, yy, weight(4), coord(2)
      real(kind = kreal)   :: inpVx(:,:), inpVy(:,:), inpP(:,:)

      x = cx
      r = dsqrt(cy*cy + cz*cz)
      theta = datan(cy/cz)

      do i=1,idx-1
        if(inpX(i) <= x &
          .and. x < inpX(i+1))then
          x1  = inpX(i)
          x2  = inpX(i+1)
          ix1 = i
          ix2 = i+1
          exit
        endif
      enddo

      do i=1,idy-1
        if(inpY(i) <= r &
          .and. r < inpY(i+1))then
          y1  = inpY(i)
          y2  = inpY(i+1)
          iy1 = i
          iy2 = i+1
          exit
        endif
      enddo

      coord(1) = (x-x1)/(x2-x1)
      coord(2) = (r-y1)/(y2-y1)
      call ShapeFunc_quad4n(coord,weight)

      p = inpP(ix1,iy1)*weight(1) &
        + inpP(ix2,iy1)*weight(2) &
        + inpP(ix2,iy2)*weight(3) &
        + inpP(ix1,iy2)*weight(4)

      trac(1) = p*normal(1)
      trac(2) = p*normal(2)
      trac(3) = p*normal(3)

  end subroutine

    subroutine DYNAMIC_MAT_ASS_COUPLE_INPUT(hecMESH, hecMAT, fstrSOLID)
      USE elementInfo
      IMPLICIT NONE

!--------------------------------------------------------------------

      TYPE(hecmwST_local_mesh) :: hecMESH
      TYPE(hecmwST_matrix)     :: hecMAT
      TYPE(fstr_solid)         :: fstrSOLID

!--------------------------------------------------------------------

      INTEGER(KIND = kint) :: n_nodes
      INTEGER(KIND = kint) :: NumOfSurface
      INTEGER(KIND = kint) :: i, j, k, idx, idy
      INTEGER(KIND = kint) :: nid, eid, fid, etype
      INTEGER(KIND = kint) :: node(27)
      INTEGER(KIND = kint) :: node_n
      INTEGER(KIND = kint), allocatable :: table_element(:)
      INTEGER(KIND = kint), allocatable :: table_face(:)
      INTEGER(KIND = kint) :: naa, fetype, nn
      INTEGER(KIND = kint) :: TotalTimeStep, TimeStep
      INTEGER(KIND = kint) :: TimeStep_Period
      INTEGER(KIND = kint) :: number

      REAL(KIND = kreal) :: trac(3), normal(3)  ! traction on surface
      REAL(KIND = kreal) :: px, py, pz
      REAL(KIND = kreal) :: vx, vy, vz
      REAL(KIND = kreal) :: cx, cy, cz
      REAL(KIND = kreal) :: xx(9), yy(9), zz(9)
      REAL(KIND = kreal) :: coord(3,4), local(2)
      REAL(KIND = kreal) :: area, wg
      REAL(KIND = kreal) :: TimeData
      REAL(KIND = kreal), allocatable :: inpX(:), inpY(:)
      REAL(KIND = kreal), allocatable :: inpVx(:,:), inpVy(:,:), inpP(:,:)

      CHARACTER(1)  :: ch
      CHARACTER(72) :: filename
      CHARACTER(72) :: dataname(2)
      CHARACTER(72) :: header(2)
      CHARACTER(72) :: testcase

!--------------------------------------------------------------------
      ! input data
      filename = "fluent.dat"

      open(9, file=filename, status="old")
        read(9,*)idx,idy
        read(9,*)ch

        allocate(inpX(idx))
        allocate(inpY(idy))
        allocate(inpVx(idx,idy))
        allocate(inpVy(idx,idy))
        allocate(inpP (idx,idy))
        inpX = 0.0d0
        inpY = 0.0d0
        inpVx= 0.0d0
        inpVy= 0.0d0
        inpP = 0.0d0

        do i=1,idx
          do j=1,idy
            read(9,*)inpX(i),inpY(j),inpVx(i,j),inpVy(i,j),inpP(i,j)
          enddo
        enddo
      close(9)


      filename = "surface.dat"

      open(9, file=filename, status="old")
        read(9,*)NumOfSurface

        allocate( table_element(NumOfSurface) )
        allocate( table_face(NumOfSurface   ) )
        table_element = 0
        table_face = 0

        do i=1,NumOfSurface
            read(9,*)table_element(i),table_face(i)
        enddo
      close(9)

!--------------------------------------------------------------------
      !allocate( table_element(NumOfSurface) )
      !allocate( table_face(NumOfSurface   ) )
      !call fstr_get_surface(hecMESH,fstrSOLID,NumOfSurface,table_element,table_face)

      local(1) = 0.5d0
      local(2) = 0.5d0

      DO j = 1, NumOfSurface

        fid = table_face(j)
        eid = table_element(j)

        call node_on_surface(hecMESH, etype, eid, fid, node, node_n)

        !--------------------------------------------------------
        DO naa =1, node_n
          nid = node(naa)
          xx(naa) = hecMESH%node( 3*(nid-1)+1 )
          yy(naa) = hecMESH%node( 3*(nid-1)+2 )
          zz(naa) = hecMESH%node( 3*(nid-1)+3 )
          coord(1,naa) = xx(naa)
          coord(2,naa) = yy(naa)
          coord(3,naa) = zz(naa)
        END DO

        fetype = 741
        nn     = 4
        normal = SurfaceNormal(fetype,nn,local,coord)

        !--------------------------------------------------------
        ! traction on surface

        px = 0.0D0
        py = 0.0D0
        pz = 0.0D0

        DO naa = 1, node_n
          cx = xx(naa)
          cy = yy(naa)
          cz = zz(naa)

          call fstr_get_trac(inpX,inpY,inpVx,inpVy,inpP,idx,idy,cx,cy,cz,normal,trac)

          px = px + trac(1)
          py = py + trac(2)
          pz = pz + trac(3)
        END DO

        ! Average in an element
        px = px/node_n
        py = py/node_n
        pz = pz/node_n

        !--------------------------------------------------------
        if( node_n == 3 ) THEN
          area = area_of_triangle(xx, yy, zz)
        elseif( node_n == 4 ) then
          area = area_of_squre(xx, yy, zz)
        elseif( node_n == 6 ) then
          area = area_of_triangle2(xx, yy, zz)
        else
          WRITE(*, *) "#Error : in FSTR_MAT_ASS_COUPLE "
          CALL hecmw_abort( hecmw_comm_get_comm() )
        endif

        !--------------------------------------------------------
        ! force on vertex
        wg = area/node_n

        vx = px*wg
        vy = py*wg
        vz = pz*wg

        !--------------------------------------------------------
        ! add in B
        DO naa = 1, node_n
          nid = node(naa)
          hecMAT%B( 3*(nid-1)+1 ) = hecMAT%B( 3*(nid-1)+1 ) + vx
          hecMAT%B( 3*(nid-1)+2 ) = hecMAT%B( 3*(nid-1)+2 ) + vy
          hecMAT%B( 3*(nid-1)+3 ) = hecMAT%B( 3*(nid-1)+3 ) + vz
        END DO

      END DO

!--------------------------------------------------------------------
      deallocate( table_element )
      deallocate( table_face )
      deallocate( inpX )
      deallocate( inpY )
      deallocate( inpVx )
      deallocate( inpVy )
      deallocate( inpP )

!--------------------------------------------------------------------
      RETURN

      END SUBROUTINE DYNAMIC_MAT_ASS_COUPLE_INPUT

!==============================================================================
! CALC AREA
!==============================================================================

function area_of_triangle( XX,YY,ZZ )
      implicit none
      real(kind=kreal) XX(*),YY(*),ZZ(*)
      real(kind=kreal) V1X,V1Y,V1Z
      real(kind=kreal) V2X,V2Y,V2Z
      real(kind=kreal) V3X,V3Y,V3Z
      real(kind=kreal) area_of_triangle
      V1X=XX(2)-XX(1)
      V1Y=YY(2)-YY(1)
      V1Z=ZZ(2)-ZZ(1)
      V2X=XX(3)-XX(1)
      V2Y=YY(3)-YY(1)
      V2Z=ZZ(3)-ZZ(1)
      V3X= V1Y*V2Z-V1Z*V2Y
      V3Y=-V1X*V2Z+V1Z*V2X
      V3Z= V1X*V2Y-V1Y*V2X
      area_of_triangle = SQRT( V3X*V3X + V3Y*V3Y + V3Z*V3Z )*0.5
end function area_of_triangle

function area_of_triangle2( XX,YY,ZZ )
      implicit none
      real(kind=kreal) XX(*),YY(*),ZZ(*)
      real(kind=kreal) area_of_triangle2
      real(kind=kreal) x(3),y(3),z(3)
      x(1)=XX(1)
      x(2)=XX(3)
      x(3)=XX(5)
      y(1)=YY(1)
      y(2)=YY(3)
      y(3)=YY(5)
      z(1)=ZZ(1)
      z(2)=ZZ(3)
      z(3)=ZZ(5)
      area_of_triangle2 = area_of_triangle(x,y,z)
end function area_of_triangle2

function area_of_triangle2A( XX,YY,ZZ )
      implicit none
      real(kind=kreal) XX(*),YY(*),ZZ(*)
      real(kind=kreal) area_of_triangle2A
      real(kind=kreal) x(3),y(3),z(3)
      x(1)=XX(1)
      x(2)=XX(2)
      x(3)=XX(3)
      y(1)=YY(1)
      y(2)=YY(2)
      y(3)=YY(3)
      z(1)=ZZ(1)
      z(2)=ZZ(2)
      z(3)=ZZ(3)
      area_of_triangle2A = area_of_triangle(x,y,z)
end function area_of_triangle2A


function area_of_squre ( XX,YY,ZZ )
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*)
      REAL(kind=kreal) area_of_squre
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=8,NG=2)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) IG1,IG2,LX,LY,LZ,I
      REAL(kind=kreal) VX,VY,VZ,XCOD,YCOD,ZCOD
      REAL(kind=kreal) AX,AY,AZ,RX,RY,RZ,HX,HY,HZ,VAL
      REAL(kind=kreal) PHX,PHY,PHZ
      REAL(kind=kreal) G1X,G1Y,G1Z
      REAL(kind=kreal) G2X,G2Y,G2Z
      REAL(kind=kreal) G3X,G3Y,G3Z
      REAL(kind=kreal) XSUM,COEFX,COEFY,COEFZ
      REAL(kind=kreal) area, XG(2),WGT(2)
      DATA WGT/1.0,1.0/
      DATA XG/-0.5773502691896, 0.5773502691896/

      area = 0.0
! INTEGRATION OVER SURFACE
        DO IG2=1,NG
          SI=XG(IG2)
          DO IG1=1,NG
            RI=XG(IG1)
            H(1)=0.25*(1.0-RI)*(1.0-SI)
            H(2)=0.25*(1.0+RI)*(1.0-SI)
            H(3)=0.25*(1.0+RI)*(1.0+SI)
            H(4)=0.25*(1.0-RI)*(1.0+SI)
            HR(1)=-.25*(1.0-SI)
            HR(2)= .25*(1.0-SI)
            HR(3)= .25*(1.0+SI)
            HR(4)=-.25*(1.0+SI)
            HS(1)=-.25*(1.0-RI)
            HS(2)=-.25*(1.0+RI)
            HS(3)= .25*(1.0+RI)
            HS(4)= .25*(1.0-RI)
            G1X=0.0
            G1Y=0.0
            G1Z=0.0
            G2X=0.0
            G2Y=0.0
            G2Z=0.0
            DO I=1,NN
              G1X=G1X+HR(I)*XX(I)
              G1Y=G1Y+HR(I)*YY(I)
              G1Z=G1Z+HR(I)*ZZ(I)
              G2X=G2X+HS(I)*XX(I)
              G2Y=G2Y+HS(I)*YY(I)
              G2Z=G2Z+HS(I)*ZZ(I)
            ENDDO
            G3X=G1Y*G2Z-G1Z*G2Y
            G3Y=G1Z*G2X-G1X*G2Z
            G3Z=G1X*G2Y-G1Y*G2X
            XSUM=DSQRT(G3X**2+G3Y**2+G3Z**2)
            G3X=G3X/XSUM
            G3Y=G3Y/XSUM
            G3Z=G3Z/XSUM
!JACOBI MATRIX
            XJ11=G1X
            XJ12=G1Y
            XJ13=G1Z
            XJ21=G2X
            XJ22=G2Y
            XJ23=G2Z
            XJ31=G3X
            XJ32=G3Y
            XJ33=G3Z
!DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
            WG=WGT(IG1)*WGT(IG2)*DET
            do i = 1, NN
              area = area + H(i)*WG
            enddo
          ENDDO
        ENDDO
      area_of_squre = area;
end function area_of_squre

!==============================================================================
! GET NODES ON SURFACE
!==============================================================================

subroutine node_on_surface( hecMESH, etype, eid, sid, node, node_n )
      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer( kind=kint ) :: eid
      integer( kind=kint ) :: sid
      integer( kind=kint ) :: node(*)
      integer( kind=kint ) :: node_n
      ! local parameters
      integer( kind=kint ) :: etype
      integer( kind=kint ) :: is
      integer( kind=kint ) :: tbl341(3,4)
      integer( kind=kint ) :: tbl342(6,4)
      integer( kind=kint ) :: tbl361(4,6)
      !! vertex id tables by definition of fstr
      data tbl341 / 1,2,3,  1,2,4,  2,3,4,  3,1,4 /
      data tbl342 / 1,5,2,6,3,7, 1,5,2,9,4,8, 2,6,3,10,4,9, 3,7,1,10,4,8 /
      data tbl361 / 1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8 /

! =============== added by K. Tagami ======== experimental ==== 2010/02/25 ==
      integer( kind=kint ) :: tbl351(4,5)
      data tbl351 / 1,2,3,1, 4,5,6,4, 1,2,5,4, 2,3,6,5, 3,1,4,6 /

      integer( kind=kint ) :: tbl362(8,6)
      data tbl362 / 1, 9, 2,10, 3,11, 4,12, &
           &        5,13, 6,14, 7,15, 8,16, &
           &        1, 9, 2,18, 6,13, 5,17, &
           &        2,10, 3,19, 7,14, 6,18, &
           &        3,11, 4,20, 8,15, 7,19, &
           &        4,12, 1,17, 5,16, 8,20 /
! =================================================================

      etype = hecMESH%elem_type( eid )
      is = hecMESH%elem_node_index( eid-1 )
      select case( etype )
      case ( 341 )
            node_n = 3
            node(1) = hecMESH%elem_node_item (is + tbl341(1,sid) )
            node(2) = hecMESH%elem_node_item (is + tbl341(2,sid) )
            node(3) = hecMESH%elem_node_item (is + tbl341(3,sid) )
      case ( 342 )
            node_n = 6
            node(1) = hecMESH%elem_node_item (is + tbl342(1,sid) )
            node(2) = hecMESH%elem_node_item (is + tbl342(2,sid) )
            node(3) = hecMESH%elem_node_item (is + tbl342(3,sid) )
            node(4) = hecMESH%elem_node_item (is + tbl342(4,sid) )
            node(5) = hecMESH%elem_node_item (is + tbl342(5,sid) )
            node(6) = hecMESH%elem_node_item (is + tbl342(6,sid) )
      case ( 361 )
            node_n = 4
            node(1) = hecMESH%elem_node_item (is + tbl361(1,sid) )
            node(2) = hecMESH%elem_node_item (is + tbl361(2,sid) )
            node(3) = hecMESH%elem_node_item (is + tbl361(3,sid) )
            node(4) = hecMESH%elem_node_item (is + tbl361(4,sid) )

! ==================== Added by K. Tagami ========== experimental == 2010/02/25
      case ( 351 )
            node_n = 4
            node(1) = hecMESH%elem_node_item (is + tbl351(1,sid) )
            node(2) = hecMESH%elem_node_item (is + tbl351(2,sid) )
            node(3) = hecMESH%elem_node_item (is + tbl351(3,sid) )
            node(4) = hecMESH%elem_node_item (is + tbl351(4,sid) )
            if ( node(1) == node(4) ) then
               node_n = 3
            endif
! --------------- K. Tagami -------------- under construction ---
!      case ( 362 )
!            node_n = 8
!            node(1) = hecMESH%elem_node_item (is + tbl362(1,sid) )
!            node(2) = hecMESH%elem_node_item (is + tbl362(2,sid) )
!            node(3) = hecMESH%elem_node_item (is + tbl362(3,sid) )
!            node(4) = hecMESH%elem_node_item (is + tbl362(4,sid) )
!            node(5) = hecMESH%elem_node_item (is + tbl362(5,sid) )
!            node(6) = hecMESH%elem_node_item (is + tbl362(6,sid) )
!            node(7) = hecMESH%elem_node_item (is + tbl362(7,sid) )
!            node(8) = hecMESH%elem_node_item (is + tbl362(8,sid) )
! --------------------------------------------------------------
! ================================================================
      case default
            if( hecMESH%my_rank==0 ) then
                  write(*,*) '##Error: not supported element type for coupling analysis ', etype
! =================== modified by K. Tagami ===== 2010/02/25
!                  write(*,*) '         This version supports element type 341,342 and 361 only.'
                  write(*,*) '         This version supports element type 341,342, 351, and 361 only.'
! ----------------------------------------------------------
            endif
            call hecmw_abort( hecmw_comm_get_comm() )
      end select

end subroutine node_on_surface

end module m_dynamic_mat_ass_couple



