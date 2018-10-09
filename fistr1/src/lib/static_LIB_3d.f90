!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>   This module provide common functions of Solid elements
module m_static_LIB_3d

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

  real(kind=kreal), parameter, private :: I3(3,3) = reshape((/ &
      &  1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/3,3/))

contains
  !----------------------------------------------------------------------*
  subroutine GEOMAT_C3(stress, mat)
    !----------------------------------------------------------------------*

    real(kind=kreal), intent(in)  :: stress(6) !> stress
    real(kind=kreal), intent(out) :: mat(6, 6) !> geometric stiff matrix

    !---------------------------------------------------------------------

    mat(1, 1) = 2.0D0*stress(1);               mat(1, 2) = 0.0D0;           mat(1, 3) = 0.0D0
    mat(1, 4) = stress(4);                     mat(1, 5) = 0.0D0;           mat(1, 6) = stress(6)
    mat(2, 1) = mat(1, 2);                     mat(2, 2) = 2.0D0*stress(2); mat(2, 3) = 0.0D0
    mat(2, 4) = stress(4);                     mat(2, 5) = stress(5);       mat(2, 6) = 0.0D0
    mat(3, 1) = mat(1, 3);                     mat(3, 2) = mat(2, 3);       mat(3, 3) = 2.0D0*stress(3)
    mat(3, 4) = 0.0D0;                         mat(3, 5) = stress(5);       mat(3, 6) = stress(6)

    mat(4, 1) = mat(1, 4);                     mat(4, 2) = mat(2, 4);                     mat(4, 3) = mat(3, 4)
    mat(4, 4) = 0.5D0*( stress(1)+stress(2) ); mat(4, 5) = 0.5D0*stress(6);               mat(4, 6) = 0.5D0*stress(5)
    mat(5, 1) = mat(1, 5);                     mat(5, 2) = mat(2, 5);                     mat(5, 3) = mat(3, 5)
    mat(5, 4) = mat(4, 5);                     mat(5, 5) = 0.5d0*( stress(3)+stress(2) ); mat(5, 6) = 0.5D0*stress(4)
    mat(6, 1) = mat(1, 6);                     mat(6, 2) = mat(2, 6);                     mat(6, 3) = mat(3, 6)
    mat(6, 4) = mat(4, 6);                     mat(6, 5) = mat(5, 6);                     mat(6, 6) = 0.5D0*( stress(1)+stress(3) );

  end subroutine


  !=====================================================================*
  !>  This subroutine calculate stiff matrix of general solid elements
  !
  !>  \author     X. YUAN, K. SATO (AdavanceSoft)
  !>  \date       2009/08/03
  !>  \version    0.00
  !----------------------------------------------------------------------*
  subroutine STF_C3                                                &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, u ,temperature)
    !----------------------------------------------------------------------*

    use mMechGauss
    use m_MatMatrix
    use m_common_struct

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)  :: etype                  !< element type
    integer(kind=kint), intent(in)  :: nn                     !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3,3)            !< variables to define matreial coordinate system
    real(kind=kreal), intent(in)    :: time                   !< current time
    real(kind=kreal), intent(in)    :: tincr                  !< time increment
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    real(kind=kreal), intent(in), optional :: u(:,:)          !< nodal displacemwent

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6), B(6, NDOF*nn), DB(6, NDOF*nn)
    real(kind=kreal) :: gderiv(nn, 3), stress(6), mat(6, 6)
    real(kind=kreal) :: det, wg
    integer(kind=kint) :: i, j, LX, serr
    real(kind=kreal) :: temp, naturalCoord(3)
    real(kind=kreal) :: spfunc(nn), gdispderiv(3, 3)
    real(kind=kreal) :: B1(6, NDOF*nn), coordsys(3, 3)
    real(kind=kreal) :: Smat(9, 9), elem(3, nn)
    real(kind=kreal) :: BN(9, NDOF*nn), SBN(9, NDOF*nn)

    !---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITE    ! enforce to infinite deformation analysis
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)

    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      if( present(temperature) ) then
        call getShapeFunc(etype, naturalcoord, spfunc)
        temp = dot_product(temperature, spfunc)
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )
      else
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys )
      end if

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
        B(1, 3*j-2)=gderiv(j, 1)
        B(2, 3*j-1)=gderiv(j, 2)
        B(3, 3*j  )=gderiv(j, 3)
        B(4, 3*j-2)=gderiv(j, 2)
        B(4, 3*j-1)=gderiv(j, 1)
        B(5, 3*j-1)=gderiv(j, 3)
        B(5, 3*j  )=gderiv(j, 2)
        B(6, 3*j-2)=gderiv(j, 3)
        B(6, 3*j  )=gderiv(j, 1)
      enddo

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG ) then
        ! ---dudx(i, j) ==> gdispderiv(i, j)
        gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*NDOF)=0.0D0
        do j=1, nn
          B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
          B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
          B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
          B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
          B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
          B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
          B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
          B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
          B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
          B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
          B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
          B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
          B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
          B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
          B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
          B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
          B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
          B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
        end do
        ! ---BL = BL0 + BL1
        do j=1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do
      end if

      DB(1:6, 1:nn*ndof) = matmul( D, B(1:6, 1:nn*ndof) )
      forall( i=1:nn*ndof,  j=1:nn*ndof )
        stiff(i, j) = stiff(i, j)+dot_product( B(:, i),  DB(:, j) )*WG
      end forall

      ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG .OR. flag==UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        BN(1:9, 1:nn*ndof) = 0.0D0
        do j = 1,  nn
          BN(1, 3*j-2) = gderiv(j, 1)
          BN(2, 3*j-1) = gderiv(j, 1)
          BN(3, 3*j  ) = gderiv(j, 1)
          BN(4, 3*j-2) = gderiv(j, 2)
          BN(5, 3*j-1) = gderiv(j, 2)
          BN(6, 3*j  ) = gderiv(j, 2)
          BN(7, 3*j-2) = gderiv(j, 3)
          BN(8, 3*j-1) = gderiv(j, 3)
          BN(9, 3*j  ) = gderiv(j, 3)
        end do
        Smat(:, :) = 0.0D0
        do j = 1, 3
          Smat(j  , j  ) = stress(1)
          Smat(j  , j+3) = stress(4)
          Smat(j  , j+6) = stress(6)
          Smat(j+3, j  ) = stress(4)
          Smat(j+3, j+3) = stress(2)
          Smat(j+3, j+6) = stress(5)
          Smat(j+6, j  ) = stress(6)
          Smat(j+6, j+3) = stress(5)
          Smat(j+6, j+6) = stress(3)
        end do
        SBN(1:9, 1:nn*ndof) = matmul( Smat(1:9, 1:9),  BN(1:9, 1:nn*ndof) )
        forall( i=1:nn*ndof,  j=1:nn*ndof )
          stiff(i, j) = stiff(i, j)+dot_product( BN(:, i),  SBN(:, j) )*WG
        end forall

      end if

    enddo ! gauss roop

  end subroutine STF_C3


  !> Distrubuted external load
  !----------------------------------------------------------------------*
  subroutine DL_C3(etype, nn, XX, YY, ZZ, RHO, LTYPE, PARAMS, VECT, NSIZE)
    !----------------------------------------------------------------------*
    !**
    !**  SET DLOAD
    !**
    !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
    !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
    !   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
    !   GRAV LTYPE=4  :GRAVITY FORCE
    !   CENT LTYPE=5  :CENTRIFUGAL LOAD
    !   P1   LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR FACE-1
    !   P2   LTYPE=20 :TRACTION IN NORMAL-DIRECTION FOR FACE-2
    !   P3   LTYPE=30 :TRACTION IN NORMAL-DIRECTION FOR FACE-3
    !   P4   LTYPE=40 :TRACTION IN NORMAL-DIRECTION FOR FACE-4
    !   P5   LTYPE=50 :TRACTION IN NORMAL-DIRECTION FOR FACE-5
    !   P6   LTYPE=60 :TRACTION IN NORMAL-DIRECTION FOR FACE-6
    ! I/F VARIABLES
    integer(kind=kint), intent(in)  :: etype, nn
    real(kind=kreal), intent(in)    :: XX(:),YY(:),ZZ(:)
    real(kind=kreal), intent(in)    :: PARAMS(0:6)
    real(kind=kreal), intent(inout) :: VECT(:)
    real(kind=kreal) RHO
    integer(kind=kint) LTYPE,NSIZE

    ! LOCAL VARIABLES
    integer(kind=kint) NDOF
    parameter(NDOF=3)
    real(kind=kreal) H(nn)
    real(kind=kreal) PLX(nn), PLY(nn), PLZ(nn)
    real(kind=kreal) XJ(3, 3), DET, WG
    integer(kind=kint) IVOL, ISUF
    integer(kind=kint) NOD(nn)
    integer(kind=kint) IG2, LX, I , SURTYPE, NSUR
    real(kind=kreal) VX, VY, VZ, XCOD, YCOD, ZCOD
    real(kind=kreal) AX, AY, AZ, RX, RY, RZ, HX, HY, HZ, val
    real(kind=kreal) PHX, PHY, PHZ
    real(kind=kreal) COEFX, COEFY, COEFZ
    real(kind=kreal) normal(3), localcoord(3), elecoord(3, nn), deriv(nn, 3)

    AX = 0.0d0; AY = 0.0d0; AZ = 0.0d0; RX = 0.0d0; RY = 0.0d0; RZ = 0.0d0;
    !
    ! SET VALUE
    !
    val = PARAMS(0)
    !
    ! SELCTION OF LOAD TYPE
    !
    IVOL=0
    ISUF=0
    if( LTYPE.LT.10 ) then
      IVOL=1
      if( LTYPE.EQ.5 ) then
        AX=PARAMS(1)
        AY=PARAMS(2)
        AZ=PARAMS(3)
        RX=PARAMS(4)
        RY=PARAMS(5)
        RZ=PARAMS(6)
      endif
    else if( LTYPE.GE.10 ) then
      ISUF=1
      call getSubFace( ETYPE, LTYPE/10, SURTYPE, NOD )
      NSUR = getNumberOfNodes( SURTYPE )
    endif
    ! CLEAR VECT
    NSIZE=nn*NDOF
    VECT(1:NSIZE)=0.0D0
    !** SURFACE LOAD
    if( ISUF==1 ) then
      ! INTEGRATION OVER SURFACE
      do I=1,NSUR
        elecoord(1,i)=XX(NOD(I))
        elecoord(2,i)=YY(NOD(i))
        elecoord(3,i)=ZZ(NOD(i))
      enddo
      do IG2=1,NumOfQuadPoints( SURTYPE )
        call getQuadPoint( SURTYPE, IG2, localcoord(1:2) )
        call getShapeFunc( SURTYPE, localcoord(1:2), H(1:NSUR) )

        WG=getWeight( SURTYPE, IG2 )
        normal=SurfaceNormal( SURTYPE, NSUR, localcoord(1:2), elecoord(:,1:NSUR) )
        do I=1,NSUR
          VECT(3*NOD(I)-2)=VECT(3*NOD(I)-2)+val*WG*H(I)*normal(1)
          VECT(3*NOD(I)-1)=VECT(3*NOD(I)-1)+val*WG*H(I)*normal(2)
          VECT(3*NOD(I)  )=VECT(3*NOD(I)  )+val*WG*H(I)*normal(3)
        enddo
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL==1 ) then
      PLX(:)=0.0D0
      PLY(:)=0.0D0
      PLZ(:)=0.0D0
      ! LOOP FOR INTEGRATION POINTS
      do  LX=1,NumOfQuadPoints( ETYPE )
        call getQuadPoint( ETYPE, LX, localcoord )
        call getShapeFunc( ETYPE, localcoord, H(1:nn) )
        call getShapeDeriv( ETYPE, localcoord, deriv )
        !  JACOBI MATRIX
        XJ(1,1:3)= matmul( xx(1:nn), deriv(1:nn,1:3) )
        XJ(2,1:3)= matmul( yy(1:nn), deriv(1:nn,1:3) )
        XJ(3,1:3)= matmul( zz(1:nn), deriv(1:nn,1:3) )
        !DETERMINANT OF JACOBIAN
        DET=XJ(1,1)*XJ(2,2)*XJ(3,3)                                                 &
          +XJ(2,1)*XJ(3,2)*XJ(1,3)                                                 &
          +XJ(3,1)*XJ(1,2)*XJ(2,3)                                                 &
          -XJ(3,1)*XJ(2,2)*XJ(1,3)                                                 &
          -XJ(2,1)*XJ(1,2)*XJ(3,3)                                                 &
          -XJ(1,1)*XJ(3,2)*XJ(2,3)

        COEFX=1.0
        COEFY=1.0
        COEFZ=1.0
        ! CENTRIFUGAL LOAD
        if( LTYPE==5 ) then
          XCOD=dot_product( H(1:nn),XX(1:nn) )
          YCOD=dot_product( H(1:nn),YY(1:nn) )
          ZCOD=dot_product( H(1:nn),ZZ(1:nn) )
          HX=AX+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RX
          HY=AY+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RY
          HZ=AZ+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RZ
          PHX=XCOD-HX
          PHY=YCOD-HY
          PHZ=ZCOD-HZ
          COEFX=RHO*val*val*PHX
          COEFY=RHO*val*val*PHY
          COEFZ=RHO*val*val*PHZ
        end if

        WG=getWeight( etype, LX )*DET
        do I=1,nn
          PLX(I)=PLX(I)+H(I)*WG*COEFX
          PLY(I)=PLY(I)+H(I)*WG*COEFY
          PLZ(I)=PLZ(I)+H(I)*WG*COEFZ
        enddo
      enddo
      if( LTYPE.EQ.1) then
        do I=1,nn
          VECT(3*I-2)=val*PLX(I)
        enddo
      else if( LTYPE.EQ.2 ) then
        do I=1,nn
          VECT(3*I-1)=val*PLY(I)
        enddo
      else if( LTYPE.EQ.3 ) then
        do I=1,nn
          VECT(3*I  )=val*PLZ(I)
        enddo
      else if( LTYPE.EQ.4 ) then
        VX=PARAMS(1)
        VY=PARAMS(2)
        VZ=PARAMS(3)
        VX=VX/sqrt(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
        VY=VY/sqrt(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
        VZ=VZ/sqrt(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
        do I=1,nn
          VECT(3*I-2)=val*PLX(I)*RHO*VX
          VECT(3*I-1)=val*PLY(I)*RHO*VY
          VECT(3*I  )=val*PLZ(I)*RHO*VZ
        enddo
      else if( LTYPE.EQ.5 ) then
        do I=1,nn
          VECT(3*I-2)=PLX(I)
          VECT(3*I-1)=PLY(I)
          VECT(3*I  )=PLZ(I)
        enddo
      end if
    endif

  end subroutine DL_C3

  !> This subroutien calculate thermal loading
  !----------------------------------------------------------------------*
  subroutine TLOAD_C3                                 &
      (etype, nn, XX, YY, ZZ, TT, T0, gausses, &
      VECT, cdsys_ID, coords)
    !----------------------------------------------------------------------*

    use m_fstr
    use mMechGauss
    use m_MatMatrix
    use m_utilities

    !---------------------------------------------------------------------

    integer(kind=kint), parameter   :: ndof = 3
    integer(kind=kint), intent(in)  :: etype, nn
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal), intent(in)    :: XX(nn), YY(nn), ZZ(nn)
    real(kind=kreal), intent(in)    :: TT(nn),T0(nn)
    real(kind=kreal), intent(out)   :: VECT(nn*NDOF)
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)           !< variables to define matreial coordinate system

    !---------------------------------------------------------------------

    real(kind=kreal) :: ALP, ALP0, D(6, 6), B(6, ndof*nn)
    real(kind=kreal) :: det, ecoord(3, nn)
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal) :: EPSTH(6),SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    real(kind=kreal) :: wg, outa(1), ina(1)
    real(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    matlaniso = .FALSE.

    if( cdsys_ID > 0 ) then   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .TRUE.
    end if

    VECT(:) = 0.0D0

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)

    ! LOOP FOR INTEGRATION POINTS
    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc( ETYPE, naturalcoord, H(1:nn) )
      call getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv )

      if( matlaniso ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      ! WEIGHT VALUE AT GAUSSIAN POINT
      wg = getWeight(etype, LX)*det
      B(1:6,1:nn*NDOF)=0.0D0
      do J=1,nn
        B(1,3*J-2)=gderiv(j,1)
        B(2,3*J-1)=gderiv(j,2)
        B(3,3*J  )=gderiv(j,3)
        B(4,3*J-2)=gderiv(j,2)
        B(4,3*J-1)=gderiv(j,1)
        B(5,3*J-1)=gderiv(j,3)
        B(5,3*J  )=gderiv(j,2)
        B(6,3*J-2)=gderiv(j,3)
        B(6,3*J  )=gderiv(j,1)
      enddo

      TEMPC = dot_product( H(1:nn), TT(1:nn) )
      TEMP0 = dot_product( H(1:nn), T0(1:nn) )

      call MatlMatrix( gausses(LX), D3, D, 1.d0, 0.0D0, coordsys, tempc )

      ina(1) = TEMPC
      if( matlaniso ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      end if
      ina(1) = TEMP0
      if( matlaniso  ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      end if

      !**
      !** THERMAL strain
      !**
      if( matlaniso ) then
        do j=1,3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        end do
        EPSTH(4:6) = 0.0D0
        call transformation(coordsys, tm)
        EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      else
        THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      end if

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = matmul( D(:, :), EPSTH(:) )
      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:nn*NDOF) = VECT(1:nn*NDOF)+matmul( SGM(:),B(:, :) )*wg

    end do

  end subroutine TLOAD_C3

  subroutine Cal_Thermal_expansion_C3( tt0, ttc, material, coordsys, matlaniso, EPSTH )
    use m_fstr
    use m_utilities
    real(kind=kreal),   INTENT(IN)     :: tt0
    real(kind=kreal),   INTENT(IN)     :: ttc
    type( tMaterial ),  INTENT(IN)     :: material
    real(kind=kreal),   INTENT(IN)     :: coordsys(3,3)
    logical,            INTENT(IN)     :: matlaniso
    real(kind=kreal),   INTENT(OUT)    :: EPSTH(6)

    integer(kind=kint) :: j
    real(kind=kreal)   :: ina(1), outa(1)
    logical            :: ierr
    real(kind=kreal)   :: alp, alp0, alpo(3), alpo0(3), tm(6,6)

    ina(1) = ttc
    if( matlaniso ) then
      call fetch_TableData( MC_ORTHOEXP, material%dict, alpo(:), ierr, ina )
      if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
    else
      call fetch_TableData( MC_THEMOEXP, material%dict, outa(:), ierr, ina )
      if( ierr ) outa(1) = material%variables(M_EXAPNSION)
      alp = outa(1)
    end if
    ina(1) = tt0
    if( matlaniso  ) then
      call fetch_TableData( MC_ORTHOEXP, material%dict, alpo0(:), ierr, ina )
      if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
    else
      call fetch_TableData( MC_THEMOEXP, material%dict, outa(:), ierr, ina )
      if( ierr ) outa(1) = material%variables(M_EXAPNSION)
      alp0 = outa(1)
    end if
    if( matlaniso ) then
      do j=1,3
        EPSTH(j) = ALPO(j)*(ttc-ref_temp)-alpo0(j)*(tt0-ref_temp)
      end do
      call transformation( coordsys(:, :), tm)
      EPSTH(:) = matmul( EPSTH(:), tm  ) ! to global coord
      EPSTH(4:6) = EPSTH(4:6)*2.0D0
    else
      EPSTH(1:3)=alp*(ttc-ref_temp)-alp0*(tt0-ref_temp)
    end if
  end subroutine

  !Hughes,  T. J. R., and J. Winget, ?gFinite Rotation Effects in Numerical Integration of
  ! Rate Constitutive Equations Arising in Large Deformation Analysis,?h
  ! International Journal for Numerical Methods in Engineering, vol. 15, pp. 1862-1867, 1980.
  subroutine Hughes_Winget_rotation_3D( rot, stress_in, stress_out )
    use m_utilities
    real(kind=kreal), intent(in)     :: rot(3,3)
    real(kind=kreal), intent(in)     :: stress_in(3,3)
    real(kind=kreal), intent(out)    :: stress_out(3,3)

    real(kind=kreal) :: dr(3,3)

    !calc dR=(I-0.5*dW)^-1*(I+0.5*dW)
    dr = I3-0.5d0*rot
    call calInverse(3, dr)
    dr = matmul(dr,I3+0.5d0*rot)

    stress_out = matmul(dr,stress_in)
    stress_out = matmul(stress_out,transpose(dr))
  end subroutine

  subroutine Update_Stress3D( flag, gauss, rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )
    use m_fstr
    use m_MatMatrix
    use mMechGauss
    use m_utilities

    type(tGaussStatus), intent(inout)       :: gauss
    integer(kind=kint), intent(in)          :: flag
    real(kind=kreal), intent(in)            :: rot(3,3)
    real(kind=kreal), intent(in)            :: dstrain(6)
    real(kind=kreal), intent(in)            :: F(3,3)        !deformation gradient (used for ss_out)
    real(kind=kreal), intent(in)            :: coordsys(3,3)
    real(kind=kreal), intent(in)            :: time
    real(kind=kreal), intent(in)            :: tincr
    real(kind=kreal), intent(in), optional  :: ttc
    real(kind=kreal), intent(in), optional  :: tt0
    real(kind=kreal), intent(in), optional  :: ttn

    integer(kind=kint) :: mtype, i, j, k
    integer(kind=kint) :: isEp
    real(kind=kreal)   :: D(6,6), dstress(6), dumstress(3,3), dum(3,3), trD, det
    real(kind=kreal)   :: tensor(6)     !< tensor
    real(kind=kreal)   :: eigval(3)     !< vector containing the eigvalches
    real(kind=kreal)   :: princ(3,3), norm !< matrix containing the three principal column vectors

    mtype = gauss%pMaterial%mtype

    if( isElastoplastic(mtype) .OR. mtype == NORTON )then
      isEp = 1
    else
      isEp = 0
    endif

    if( present(ttc) .AND. present(ttn) ) then
      call MatlMatrix( gauss, D3, D, time, tincr, coordsys, ttc, isEp )
    else
      call MatlMatrix( gauss, D3, D, time, tincr, coordsys, isEp=isEp)
    end if

    if( flag == INFINITE ) then

      gauss%stress(1:6) = matmul( D(1:6, 1:6), dstrain(1:6) )
      if( isViscoelastic(mtype) .AND. tincr /= 0.0D0 ) then
        if( present(ttc) .AND. present(ttn) ) then
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr, ttc, ttn )
        else
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr )
        end if
        gauss%stress = real(gauss%stress)
      end if

    else if( flag == TOTALLAG ) then

      if( mtype == NEOHOOKE .OR. mtype == MOONEYRIVLIN .OR.  mtype == ARRUDABOYCE  .OR.   &
          mtype==USERELASTIC .OR. mtype==USERHYPERELASTIC .OR. mtype==USERMATERIAL ) then
        call StressUpdate( gauss, D3, dstrain, gauss%stress )
      else if( ( isViscoelastic(mtype) .OR. mtype == NORTON ) .AND. tincr /= 0.0D0 ) then
        gauss%pMaterial%mtype=mtype
        if( present(ttc) .AND. present(ttn) ) then
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr, ttc, ttn )
        else
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr )
        end if
      else
        gauss%stress(1:6) = matmul( D(1:6, 1:6), dstrain(1:6) )
      end if

    else if( flag == UPDATELAG ) then
      !  CALL GEOMAT_C3( gauss%stress, mat )
      !  D(:, :) = D(:, :)+mat(:, :)

      if( isViscoelastic(mtype) .AND. tincr /= 0.0D0 ) then
        if( present(ttc) .AND. present(ttn) ) then
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr, ttc, tt0 )
        else
          call StressUpdate( gauss, D3, dstrain, gauss%stress, time, tincr )
        end if
      else
        dstress = real( matmul( D(1:6,1:6), dstrain(1:6) ) )
        dumstress(1,1) = gauss%stress_bak(1)
        dumstress(2,2) = gauss%stress_bak(2)
        dumstress(3,3) = gauss%stress_bak(3)
        dumstress(1,2) = gauss%stress_bak(4);  dumstress(2,1)=dumstress(1,2)
        dumstress(2,3) = gauss%stress_bak(5);  dumstress(3,2)=dumstress(2,3)
        dumstress(3,1) = gauss%stress_bak(6);  dumstress(1,3)=dumstress(3,1)

        !stress integration
        trD = dstrain(1)+dstrain(2)+dstrain(3)
        dum(:,:) = dumstress + matmul( rot,dumstress ) - matmul( dumstress, rot ) + dumstress*trD
        !call Hughes_Winget_rotation_3D( rot, dumstress, dum )

        gauss%stress(1) = dum(1,1) + dstress(1)
        gauss%stress(2) = dum(2,2) + dstress(2)
        gauss%stress(3) = dum(3,3) + dstress(3)
        gauss%stress(4) = dum(1,2) + dstress(4)
        gauss%stress(5) = dum(2,3) + dstress(5)
        gauss%stress(6) = dum(3,1) + dstress(6)

        if( mtype == USERMATERIAL ) then
          call StressUpdate( gauss, D3, dstrain, gauss%stress )
        else if( mtype == NORTON ) then
          !gauss%pMaterial%mtype = mtype
          if( tincr /= 0.0D0 .AND. any( gauss%stress /= 0.0D0 ) ) then
            !gauss%pMaterial%mtype = mtype
            if( present(ttc) .AND. present(ttn) ) then
              call StressUpdate( gauss, D3, gauss%strain, gauss%stress, time, tincr, ttc, ttn )
            else
              call StressUpdate( gauss, D3, gauss%strain, gauss%stress, time, tincr )
            end if
          end if
        end if
      end if

    end if

    if( isElastoplastic(mtype) ) then
      if( present(ttc) ) then
        call BackwardEuler( gauss%pMaterial, gauss%stress, gauss%plstrain, &
          gauss%istatus(1), gauss%fstatus, ttc )
      else
        call BackwardEuler( gauss%pMaterial, gauss%stress, gauss%plstrain, &
          gauss%istatus(1), gauss%fstatus )
      end if
    end if

    !convert stress/strain measure for output
    if( OPSSTYPE == kOPSS_SOLUTION ) then

      if( flag == INFINITE ) then !linear
        gauss%stress_out(1:6) = gauss%stress(1:6)
        gauss%strain_out(1:6) = gauss%strain(1:6)
      else !nonlinear
        !convert stress
        if( flag == TOTALLAG ) then
          dumstress(1,1) = gauss%stress(1)
          dumstress(2,2) = gauss%stress(2)
          dumstress(3,3) = gauss%stress(3)
          dumstress(1,2) = gauss%stress(4);  dumstress(2,1)=dumstress(1,2)
          dumstress(2,3) = gauss%stress(5);  dumstress(3,2)=dumstress(2,3)
          dumstress(3,1) = gauss%stress(6);  dumstress(1,3)=dumstress(3,1)

          det = Determinant33(F)
          if( det == 0.d0 ) stop "Fail to convert stress: detF=0"
          ! cauchy stress = (1/detF)*F*(2ndPK stress)*F^T
          dumstress(1:3,1:3) = matmul(dumstress(1:3,1:3),transpose(F(1:3,1:3)))
          dumstress(1:3,1:3) = (1.d0/det)*matmul(F(1:3,1:3),dumstress(1:3,1:3))

          gauss%stress_out(1) = dumstress(1,1)
          gauss%stress_out(2) = dumstress(2,2)
          gauss%stress_out(3) = dumstress(3,3)
          gauss%stress_out(4) = dumstress(1,2)
          gauss%stress_out(5) = dumstress(2,3)
          gauss%stress_out(6) = dumstress(3,1)
        else if( flag == UPDATELAG ) then
          gauss%stress_out(1:6) = gauss%stress(1:6)
        endif

        !calc logarithmic strain
        dum(1:3,1:3) = matmul(F(1:3,1:3),transpose(F(1:3,1:3)))
        tensor(1) = dum(1,1)
        tensor(2) = dum(2,2)
        tensor(3) = dum(3,3)
        tensor(4) = dum(1,2)
        tensor(5) = dum(2,3)
        tensor(6) = dum(3,1)
        call get_principal(tensor, eigval, princ)

        do k=1,3
          if( eigval(k) <= 0.d0 ) stop "Fail to calc log strain: stretch<0"
          eigval(k) = 0.5d0*dlog(eigval(k)) !log(sqrt(lambda))
          norm = dsqrt(dot_product(princ(1:3,k),princ(1:3,k)))
          if( norm <= 0.d0 ) stop "Fail to calc log strain: stretch direction vector=0"
          princ(1:3,k) = princ(1:3,k)/norm
        end do
        do i=1,3
          do j=1,3
            dum(i,j) = 0.d0
            do k=1,3
              dum(i,j) = dum(i,j) + eigval(k)*princ(i,k)*princ(j,k)
            end do
          end do
        end do
        gauss%strain_out(1) = dum(1,1)
        gauss%strain_out(2) = dum(2,2)
        gauss%strain_out(3) = dum(3,3)
        gauss%strain_out(4) = 2.d0*dum(1,2)
        gauss%strain_out(5) = 2.d0*dum(2,3)
        gauss%strain_out(6) = 2.d0*dum(3,1)
      endif

    else
      gauss%stress_out(1:6) = gauss%stress(1:6)
      gauss%strain_out(1:6) = gauss%strain(1:6)
    end if

  end subroutine

  !> Update strain and stress inside element
  !---------------------------------------------------------------------*
  subroutine UPDATE_C3                                       &
      (etype,nn,ecoord, u, ddu, cdsys_ID, coords, qf, &
      gausses, iter, time, tincr, TT, T0, TN)
    !---------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities

    integer(kind=kint), intent(in)    :: etype         !< \param [in] element type
    integer(kind=kint), intent(in)    :: nn            !< \param [in] number of elemental nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: u(3, nn)      !< \param [in] nodal dislplacements
    real(kind=kreal), intent(in)      :: ddu(3, nn)    !< \param [in] nodal displacement
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define matreial coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    integer, intent(in)               :: iter
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(in), optional :: TT(nn)   !< current temperature
    real(kind=kreal), intent(in), optional :: T0(nn)   !< reference temperature
    real(kind=kreal), intent(in), optional :: TN(nn)   !< reference temperature

    ! LCOAL VARIAVLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: D(6,6), B(6,ndof*nn), B1(6,ndof*nn), spfunc(nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gderiv1(nn,3), gdispderiv(3,3), F(3,3), det, det1, WG, ttc,tt0, ttn,outa(1)
    integer(kind=kint) :: i, j, k, LX, serr
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), mat(6,6), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn), elem(3,nn), elem1(3,nn), coordsys(3,3)
    real(kind=kreal)   :: dstrain(6)
    real(kind=kreal)   :: alpo(3)
    logical            :: ierr, matlaniso

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:,:) = ecoord(:,:)
    totaldisp(:,:) = u(:,:)+ ddu(:,:)
    if( flag == UPDATELAG ) then
      elem(:,:) = (0.5D0*ddu(:,:)+u(:,:) ) +ecoord(:,:)
      elem1(:,:) = (ddu(:,:)+u(:,:) ) +ecoord(:,:)
      ! elem = elem1
      totaldisp(:,:) = ddu(:,:)
    end if

    matlaniso = .FALSE.
    if( present(TT) .and. cdsys_ID > 0 ) then
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .true.
    end if

    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      ! ========================================================
      ! UPDATE STRAIN and STRESS
      ! ========================================================

      ! Thermal Strain
      EPSTH = 0.0D0
      if( present(tt) .AND. present(t0) ) then
        call getShapeFunc(etype, naturalcoord, spfunc)
        ttc = dot_product(TT, spfunc)
        tt0 = dot_product(T0, spfunc)
        ttn = dot_product(TN, spfunc)
        call Cal_Thermal_expansion_C3( tt0, ttc, gausses(LX)%pMaterial, coordsys, matlaniso, EPSTH )
      end if

      ! Update strain
      ! Small strain
      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)   ! allright?

      F(1:3,1:3) = 0.d0; F(1,1)=1.d0; F(2,2)=1.d0; F(3,3)=1.d0; !deformation gradient
      if( flag == INFINITE ) then
        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(1:6)

      else if( flag == TOTALLAG ) then
        ! Green-Lagrange strain
        dstrain(1) = dstrain(1)+0.5d0*dot_product( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5d0*dot_product( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5d0*dot_product( gdispderiv(:, 3), gdispderiv(:, 3) )
        dstrain(4) = dstrain(4)+( gdispderiv(1, 1)*gdispderiv(1, 2)                                     &
          +gdispderiv(2, 1)*gdispderiv(2, 2)+gdispderiv(3, 1)*gdispderiv(3, 2) )
        dstrain(5) = dstrain(5)+( gdispderiv(1, 2)*gdispderiv(1, 3)                                     &
          +gdispderiv(2, 2)*gdispderiv(2, 3)+gdispderiv(3, 2)*gdispderiv(3, 3) )
        dstrain(6) = dstrain(6)+( gdispderiv(1, 1)*gdispderiv(1, 3)                                     &
          +gdispderiv(2, 1)*gdispderiv(2, 3)+gdispderiv(3, 1)*gdispderiv(3, 3) )

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
        F(1:3,1:3) = F(1:3,1:3) + gdispderiv(1:3,1:3)

      else if( flag == UPDATELAG ) then
        rot = 0.0D0
        rot(1, 2)= 0.5d0*(gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5d0*(gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5d0*(gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)

        call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det1, gderiv1)
        F(1:3,1:3) = F(1:3,1:3) + matmul( u(1:ndof, 1:nn)+ddu(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )

      end if

      ! Update stress
      if( present(tt) .AND. present(t0) ) then
        call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )
      else
        call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr )
      end if

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      do J=1,nn
        B(1,3*j-2) = gderiv(j, 1)
        B(2,3*j-1) = gderiv(j, 2)
        B(3,3*j  ) = gderiv(j, 3)
        B(4,3*j-2) = gderiv(j, 2)
        B(4,3*j-1) = gderiv(j, 1)
        B(5,3*j-1) = gderiv(j, 3)
        B(5,3*j  ) = gderiv(j, 2)
        B(6,3*j-2) = gderiv(j, 3)
        B(6,3*j  ) = gderiv(j, 1)
      end do

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == INFINITE ) then

      else if( flag == TOTALLAG ) then

        gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*ndof)=0.0D0
        do j = 1,nn
          B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
          B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
          B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
          B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
          B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
          B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
          B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
          B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
          B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
          B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
          B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
          B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
          B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
          B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
          B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
          B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
          B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
          B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
        end do
        ! BL = BL0 + BL1
        do j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        end do

      else if( flag == UPDATELAG ) then

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
        B(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B(1, 3*J-2) = gderiv(j, 1)
          B(2, 3*J-1) = gderiv(j, 2)
          B(3, 3*J  ) = gderiv(j, 3)
          B(4, 3*J-2) = gderiv(j, 2)
          B(4, 3*J-1) = gderiv(j, 1)
          B(5, 3*J-1) = gderiv(j, 3)
          B(5, 3*J  ) = gderiv(j, 2)
          B(6, 3*J-2) = gderiv(j, 3)
          B(6, 3*J  ) = gderiv(j, 1)
        end do

      end if

      ! calculate the Internal Force
      WG=getWeight( etype, LX )*DET
      qf(1:nn*ndof)                                                          &
        = qf(1:nn*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG

    end do

  end subroutine UPDATE_C3
  !
  !----------------------------------------------------------------------*
  subroutine NodalStress_C3(etype, nn, gausses, ndstrain, ndstress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in) :: etype, nn
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: ndstrain(nn,6)
    real(kind=kreal), intent(out)  :: ndstress(nn,6)

    !---------------------------------------------------------------------

    integer :: i, ic
    real(kind=kreal) :: TEMP(12)

    !---------------------------------------------------------------------

    TEMP(:) = 0.0D0

    IC = NumOfQuadPoints(etype)

    do i = 1, IC
      TEMP(1:6)  = TEMP(1:6) +gausses(i)%strain_out(1:6)
      TEMP(7:12) = TEMP(7:12)+gausses(i)%stress_out(1:6)
    end do

    TEMP(1:12) = TEMP(1:12)/IC

    forall( i=1:nn )
      ndstrain(i, 1:6) = TEMP(1:6)
      ndstress(i, 1:6) = TEMP(7:12)
    end forall

  end subroutine NodalStress_C3


  !----------------------------------------------------------------------*
  subroutine ElementStress_C3(etype, gausses, strain, stress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in) :: etype
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: strain(6)
    real(kind=kreal), intent(out)  :: stress(6)

    !---------------------------------------------------------------------

    integer :: i, ic

    !---------------------------------------------------------------------

    strain(:) = 0.0D0; stress(:) = 0.0D0

    IC = NumOfQuadPoints(etype)

    do i = 1, IC
      strain(:) = strain(:)+gausses(i)%strain_out(1:6)
      stress(:) = stress(:)+gausses(i)%stress_out(1:6)
    enddo

    strain(:) = strain(:)/IC
    stress(:) = stress(:)/IC

  end subroutine ElementStress_C3


  !> Volume of element
  !----------------------------------------------------------------------*
  real(kind=kreal) function VOLUME_C3(etype, nn, XX, YY, ZZ)
    !----------------------------------------------------------------------*

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in)   :: XX(:), YY(:), ZZ(:)

    !---------------------------------------------------------------------

    real(kind=kreal) :: XJ(3, 3), det, wg
    integer(kind=kint) :: LX, i
    real(kind=kreal) :: localcoord(3), deriv(nn, 3)

    !---------------------------------------------------------------------

    VOLUME_C3 = 0.0D0

    ! LOOP FOR INTEGRATION POINTS
    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint(etype, LX, localcoord)
      call getShapeDeriv(etype, localcoord, deriv)

      ! JACOBI MATRIX
      XJ(1, 1:3)= matmul( xx(1:nn), deriv(1:nn,1:3) )
      XJ(2, 1:3)= matmul( yy(1:nn), deriv(1:nn,1:3) )
      XJ(3, 1:3)= matmul( zz(1:nn), deriv(1:nn,1:3) )

      ! DETERMINANT OF JACOBIAN
      det = XJ(1, 1)*XJ(2, 2)*XJ(3, 3) &
        +XJ(2, 1)*XJ(3, 2)*XJ(1, 3) &
        +XJ(3, 1)*XJ(1, 2)*XJ(2, 3) &
        -XJ(3, 1)*XJ(2, 2)*XJ(1, 3) &
        -XJ(2, 1)*XJ(1, 2)*XJ(3, 3) &
        -XJ(1, 1)*XJ(3, 2)*XJ(2, 3)

      VOLUME_C3 = VOLUME_C3+getWeight(etype, LX)*det

    end do

  end function VOLUME_C3


end module m_static_LIB_3d
