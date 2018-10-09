!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
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
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_couple)        :: fstrCPL
    ! local
    integer(kind=kint) :: ig0, ig, is, ie, ik
    integer(kind=kint) :: i, j, count
    integer(kind=kint) :: eid, sid, etype
    integer(kind=kint) :: node(20)
    integer(kind=kint) :: node_n
    real(kind=kreal)   :: px, py, pz ! traction on surface
    real(kind=kreal)   :: vx, vy, vz ! force on vertex

    ! ================== modified by K. Tagami ============= 2010/02/25 ====
    !!!!      real(kind=kreal) :: xx(4), yy(4), zz(4)
    integer, parameter :: MaxNumNodesOnSurf = 8
    !parameter(MaxNumNodesOnSurf=8)
    real(kind=kreal) :: xx(MaxNumNodesOnSurf)
    real(kind=kreal) :: yy(MaxNumNodesOnSurf)
    real(kind=kreal) :: zz(MaxNumNodesOnSurf)
    ! ================================================================

    real(kind=kreal)   :: area, wg
    integer(kind=kint) :: ierr

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
        area = 0.0d0
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

  !==============================================================================
  ! CALC AREA
  !==============================================================================

  function area_of_triangle( XX,YY,ZZ )
    implicit none
    real(kind=kreal) :: XX(*), YY(*), ZZ(*)
    real(kind=kreal) :: V1X, V1Y, V1Z
    real(kind=kreal) :: V2X, V2Y, V2Z
    real(kind=kreal) :: V3X, V3Y, V3Z
    real(kind=kreal) :: area_of_triangle

    V1X=XX(2)-XX(1)
    V1Y=YY(2)-YY(1)
    V1Z=ZZ(2)-ZZ(1)
    V2X=XX(3)-XX(1)
    V2Y=YY(3)-YY(1)
    V2Z=ZZ(3)-ZZ(1)
    V3X= V1Y*V2Z-V1Z*V2Y
    V3Y=-V1X*V2Z+V1Z*V2X
    V3Z= V1X*V2Y-V1Y*V2X
    area_of_triangle = sqrt( V3X*V3X + V3Y*V3Y + V3Z*V3Z )*0.5
  end function area_of_triangle

  function area_of_triangle2( XX,YY,ZZ )
    implicit none
    real(kind=kreal) :: XX(*), YY(*), ZZ(*)
    real(kind=kreal) :: area_of_triangle2
    real(kind=kreal) :: x(3), y(3), z(3)

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
    implicit none
    ! I/F VARIABLES
    real(kind=kreal) :: XX(*), YY(*), ZZ(*)
    real(kind=kreal) :: area_of_squre
    ! LOCAL VARIABLES
    integer(kind=kint), parameter :: NN = 8
    integer(kind=kint), parameter :: NG = 2
    !parameter(NN=8, NG=2)
    real(kind=kreal)   :: H(NN), HR(NN), HS(NN), HT(NN)
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33, DET, WG
    integer(kind=kint) :: IG1, IG2, LX, LY, LZ, I
    real(kind=kreal)   :: VX, VY, VZ, XCOD, YCOD, ZCOD
    real(kind=kreal)   :: AX, AY, AZ, RX, RY, RZ, HX, HY, HZ, val
    real(kind=kreal)   :: PHX, PHY, PHZ
    real(kind=kreal)   :: G1X, G1Y, G1Z
    real(kind=kreal)   :: G2X, G2Y, G2Z
    real(kind=kreal)   :: G3X, G3Y, G3Z
    real(kind=kreal)   :: XSUM, COEFX, COEFY, COEFZ
    real(kind=kreal)   :: area, XG(2), WGT(2)
    data WGT/1.0, 1.0/
    data XG/-0.5773502691896, 0.5773502691896/

    area = 0.0
    ! INTEGRATION OVER SURFACE
    do IG2=1,NG
      SI=XG(IG2)
      do IG1=1,NG
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
        do I=1,NN
          G1X=G1X+HR(I)*XX(I)
          G1Y=G1Y+HR(I)*YY(I)
          G1Z=G1Z+HR(I)*ZZ(I)
          G2X=G2X+HS(I)*XX(I)
          G2Y=G2Y+HS(I)*YY(I)
          G2Z=G2Z+HS(I)*ZZ(I)
        enddo
        G3X=G1Y*G2Z-G1Z*G2Y
        G3Y=G1Z*G2X-G1X*G2Z
        G3Z=G1X*G2Y-G1Y*G2X
        XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
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
        DET=XJ11*XJ22*XJ33 &
          +XJ12*XJ23*XJ31 &
          +XJ13*XJ21*XJ32 &
          -XJ13*XJ22*XJ31 &
          -XJ12*XJ21*XJ33 &
          -XJ11*XJ23*XJ32
        WG=WGT(IG1)*WGT(IG2)*DET
        do i = 1, NN
          area = area + H(i)*WG
        enddo
      enddo
    enddo
    area_of_squre = area;
  end function area_of_squre

  !==============================================================================
  ! GET NODES ON SURFACE
  !==============================================================================

  subroutine node_on_surface( hecMESH, etype, eid, sid, node, node_n )
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: eid
    integer(kind=kint) :: sid
    integer(kind=kint) :: node(*)
    integer(kind=kint) :: node_n
    ! local parameters
    integer(kind=kint) :: etype
    integer(kind=kint) :: is
    integer(kind=kint) :: tbl341(3, 4)
    integer(kind=kint) :: tbl342(6, 4)
    integer(kind=kint) :: tbl361(4, 6)
    !! vertex id tables by definition of fstr
    data tbl341 / 1,2,3,  1,2,4,  2,3,4,  3,1,4 /
    data tbl342 / 1,5,2,6,3,7, 1,5,2,9,4,8, 2,6,3,10,4,9, 3,7,1,10,4,8 /
    data tbl361 / 1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8 /

    ! =============== added by K. Tagami ======== experimental ==== 2010/02/25 ==
    integer(kind=kint) :: tbl351(4,5)
    data tbl351 / 1,2,3,1, 4,5,6,4, 1,2,5,4, 2,3,6,5, 3,1,4,6 /

    integer(kind=kint) :: tbl362(8,6)
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



