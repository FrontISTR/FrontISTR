!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.7                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN, K. SATO (AdavanceSoft)                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!>   This module provide common functions of Solid elements
!>
!>  \author                date                  version
!>  X.Yuan(Advancesoft)    2009/08/03        original
!>  X.Yuan                 2013/03/18        consider anisotropic materials
!======================================================================!
MODULE m_static_LIB_3d

   USE hecmw, only : kint, kreal
   USE elementInfo

   IMPLICIT NONE

   CONTAINS


!----------------------------------------------------------------------*
   SUBROUTINE GEOMAT_C3(stress, mat)
!----------------------------------------------------------------------*

    REAL(kind=kreal), INTENT(IN)  :: stress(6) !> stress
    REAL(kind=kreal), INTENT(OUT) :: mat(6, 6) !> geometric stiff matrix

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

   END SUBROUTINE


!=====================================================================*
!>  This subroutine calculate stiff matrix of general solid elements
!
!>  \author     X. YUAN, K. SATO (AdavanceSoft)
!>  \date       2009/08/03
!>  \version    0.00
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3                                                &
              (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
               tincr, u ,temperature)
!----------------------------------------------------------------------*

    USE mMechGauss
    USE m_MatMatrix
    USE m_common_struct

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN)  :: etype                  !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                     !< number of elemental nodes
    REAL(kind=kreal),   INTENT(IN)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)             !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)             !< stiff matrix
    INTEGER(kind=kint), INTENT(IN)  :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)            !< variables to define matreial coordinate system
    REAL(kind=kreal), INTENT(IN)    :: tincr                  !< time increment
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: temperature(nn) !< temperature
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: u(:,:)          !< nodal displacemwent

!---------------------------------------------------------------------

    INTEGER(kind=kint) :: flag
    INTEGER(kind=kint), PARAMETER :: ndof = 3
    REAL(kind=kreal) :: D(6, 6), B(6, NDOF*nn), DB(6, NDOF*nn)
    REAL(kind=kreal) :: gderiv(nn, 3), stress(6), mat(6, 6)
    REAL(kind=kreal) :: det, wg
    INTEGER(kind=kint) :: i, j, LX, serr
    REAL(kind=kreal) :: temp, naturalCoord(3)
    REAL(kind=kreal) :: spfunc(nn), gdispderiv(3, 3)
    REAL(kind=kreal) :: B1(6, NDOF*nn), coordsys(3, 3)
    REAL(kind=kreal) :: Smat(9, 9), elem(3, nn)
    REAL(kind=kreal) :: BN(9, NDOF*nn), SBN(9, NDOF*nn)

!---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    IF( .NOT. PRESENT(u) ) flag = INFINITE    ! enforce to infinite deformation analysis
    elem(:, :) = ecoord(:, :)
    IF( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)

    DO LX = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      IF( PRESENT(temperature) ) THEN
        CALL getShapeFunc(etype, naturalcoord, spfunc)
        temp = DOT_PRODUCT(temperature, spfunc)
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, temp )
      ELSE
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys )
      END IF

      IF( flag == UPDATELAG ) then
        CALL GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      ENDIF

      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      DO j = 1, nn
        B(1, 3*j-2)=gderiv(j, 1)
        B(2, 3*j-1)=gderiv(j, 2)
        B(3, 3*j  )=gderiv(j, 3)
        B(4, 3*j-2)=gderiv(j, 2)
        B(4, 3*j-1)=gderiv(j, 1)
        B(5, 3*j-1)=gderiv(j, 3)
        B(5, 3*j  )=gderiv(j, 2)
        B(6, 3*j-2)=gderiv(j, 3)
        B(6, 3*j  )=gderiv(j, 1)
      ENDDO

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      IF( flag == TOTALLAG ) THEN
        ! ---dudx(i, j) ==> gdispderiv(i, j)
        gdispderiv(1:ndof, 1:ndof) = MATMUL( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*NDOF)=0.0D0
        DO j=1, nn
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
        END DO
        ! ---BL = BL0 + BL1
        DO j=1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        END DO
      END IF

      DB(1:6, 1:nn*ndof) = MATMUL( D, B(1:6, 1:nn*ndof) )
      forall( i=1:nn*ndof,  j=1:nn*ndof )
        stiff(i, j) = stiff(i, j)+DOT_PRODUCT( B(:, i),  DB(:, j) )*WG
      end forall

      ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      IF( flag == TOTALLAG .OR. flag==UPDATELAG ) THEN
        stress(1:6) = gausses(LX)%stress
        BN(1:9, 1:nn*ndof) = 0.0D0
        DO j = 1,  nn
          BN(1, 3*j-2) = gderiv(j, 1)
          BN(2, 3*j-1) = gderiv(j, 1)
          BN(3, 3*j  ) = gderiv(j, 1)
          BN(4, 3*j-2) = gderiv(j, 2)
          BN(5, 3*j-1) = gderiv(j, 2)
          BN(6, 3*j  ) = gderiv(j, 2)
          BN(7, 3*j-2) = gderiv(j, 3)
          BN(8, 3*j-1) = gderiv(j, 3)
          BN(9, 3*j  ) = gderiv(j, 3)
        END DO
        Smat(:, :) = 0.0D0
        DO j = 1, 3
          Smat(j  , j  ) = stress(1)
          Smat(j  , j+3) = stress(4)
          Smat(j  , j+6) = stress(6)
          Smat(j+3, j  ) = stress(4)
          Smat(j+3, j+3) = stress(2)
          Smat(j+3, j+6) = stress(5)
          Smat(j+6, j  ) = stress(6)
          Smat(j+6, j+3) = stress(5)
          Smat(j+6, j+6) = stress(3)
        END DO
        SBN(1:9, 1:nn*ndof) = MATMUL( Smat(1:9, 1:9),  BN(1:9, 1:nn*ndof) )
        forall( i=1:nn*ndof,  j=1:nn*ndof )
          stiff(i, j) = stiff(i, j)+DOT_PRODUCT( BN(:, i),  SBN(:, j) )*WG
        end forall

      END IF

    ENDDO ! gauss roop

   END SUBROUTINE STF_C3


!> Distrubuted external load
!----------------------------------------------------------------------*
   SUBROUTINE DL_C3(etype, nn, XX, YY, ZZ, RHO, LTYPE, PARAMS, VECT, NSIZE)
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
      INTEGER(kind=kint), INTENT(IN)  :: etype, nn
      REAL(kind=kreal), INTENT(IN)    :: XX(:),YY(:),ZZ(:)
      REAL(kind=kreal), INTENT(IN)    :: PARAMS(0:6)
      REAL(kind=kreal), INTENT(INOUT) :: VECT(:)
      REAL(kind=kreal) RHO
      INTEGER(kind=kint) LTYPE,NSIZE

! LOCAL VARIABLES
      INTEGER(kind=kint) NDOF
      PARAMETER(NDOF=3)
      REAL(kind=kreal) H(nn)
      REAL(kind=kreal) PLX(nn), PLY(nn), PLZ(nn)
      REAL(kind=kreal) XJ(3, 3), DET, WG
      INTEGER(kind=kint) IVOL, ISUF
      INTEGER(kind=kint) NOD(nn)
      INTEGER(kind=kint) IG2, LX, I , SURTYPE, NSUR
      REAL(kind=kreal) VX, VY, VZ, XCOD, YCOD, ZCOD
      REAL(kind=kreal) AX, AY, AZ, RX, RY, RZ, HX, HY, HZ, VAL
      REAL(kind=kreal) PHX, PHY, PHZ
      REAL(kind=kreal) COEFX, COEFY, COEFZ
      REAL(kind=kreal) normal(3), localcoord(3), elecoord(3, nn), deriv(nn, 3)
!
! SET VALUE
!
      VAL = PARAMS(0)
!
! SELCTION OF LOAD TYPE
!
      IVOL=0
      ISUF=0
      IF( LTYPE.LT.10 ) THEN
        IVOL=1
        IF( LTYPE.EQ.5 ) THEN
          AX=PARAMS(1)
          AY=PARAMS(2)
          AZ=PARAMS(3)
          RX=PARAMS(4)
          RY=PARAMS(5)
          RZ=PARAMS(6)
        ENDIF
      ELSE IF( LTYPE.GE.10 ) THEN
        ISUF=1
        CALL getSubFace( ETYPE, LTYPE/10, SURTYPE, NOD )
        NSUR = getNumberOfNodes( SURTYPE )
      ENDIF
! CLEAR VECT
      NSIZE=nn*NDOF
      VECT(1:NSIZE)=0.0D0
!** SURFACE LOAD
      IF( ISUF==1 ) THEN
! INTEGRATION OVER SURFACE
        DO I=1,NSUR
          elecoord(1,i)=XX(NOD(I))
          elecoord(2,i)=YY(NOD(i))
          elecoord(3,i)=ZZ(NOD(i))
        ENDDO
        DO IG2=1,NumOfQuadPoints( SURTYPE )
            CALL getQuadPoint( SURTYPE, IG2, localcoord(1:2) )
            CALL getShapeFunc( SURTYPE, localcoord(1:2), H(1:NSUR) )

            WG=getWeight( SURTYPE, IG2 )
            normal=SurfaceNormal( SURTYPE, NSUR, localcoord(1:2), elecoord(:,1:NSUR) )
            DO I=1,NSUR
              VECT(3*NOD(I)-2)=VECT(3*NOD(I)-2)+VAL*WG*H(I)*normal(1)
              VECT(3*NOD(I)-1)=VECT(3*NOD(I)-1)+VAL*WG*H(I)*normal(2)
              VECT(3*NOD(I)  )=VECT(3*NOD(I)  )+VAL*WG*H(I)*normal(3)
            ENDDO
        ENDDO
      ENDIF
!** VOLUME LOAD
      IF( IVOL==1 ) THEN
        PLX(:)=0.0D0
        PLY(:)=0.0D0
        PLZ(:)=0.0D0
! LOOP FOR INTEGRATION POINTS
        DO  LX=1,NumOfQuadPoints( ETYPE )
              CALL getQuadPoint( ETYPE, LX, localcoord )
              CALL getShapeFunc( ETYPE, localcoord, H(1:nn) )
              CALL getShapeDeriv( ETYPE, localcoord, deriv )
!  JACOBI MATRIX
              XJ(1,1:3)= MATMUL( xx(1:nn), deriv(1:nn,1:3) )
              XJ(2,1:3)= MATMUL( yy(1:nn), deriv(1:nn,1:3) )
              XJ(3,1:3)= MATMUL( zz(1:nn), deriv(1:nn,1:3) )
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
              IF( LTYPE==5 ) THEN
                XCOD=DOT_PRODUCT( H(1:nn),XX(1:nn) )
                YCOD=DOT_PRODUCT( H(1:nn),YY(1:nn) )
                ZCOD=DOT_PRODUCT( H(1:nn),ZZ(1:nn) )
                HX=AX+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RX
                HY=AY+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RY
                HZ=AZ+((XCOD-AX)*RX+(YCOD-AY)*RY+(ZCOD-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RZ
                PHX=XCOD-HX
                PHY=YCOD-HY
                PHZ=ZCOD-HZ
                COEFX=RHO*VAL*VAL*PHX
                COEFY=RHO*VAL*VAL*PHY
                COEFZ=RHO*VAL*VAL*PHZ
              END IF

              WG=getWeight( etype, LX )*DET
              DO I=1,nn
                PLX(I)=PLX(I)+H(I)*WG*COEFX
                PLY(I)=PLY(I)+H(I)*WG*COEFY
                PLZ(I)=PLZ(I)+H(I)*WG*COEFZ
              ENDDO
        ENDDO
        IF( LTYPE.EQ.1) THEN
          DO I=1,nn
            VECT(3*I-2)=VAL*PLX(I)
          ENDDO
        ELSE IF( LTYPE.EQ.2 ) THEN
          DO I=1,nn
            VECT(3*I-1)=VAL*PLY(I)
          ENDDO
        ELSE IF( LTYPE.EQ.3 ) THEN
          DO I=1,nn
            VECT(3*I  )=VAL*PLZ(I)
          ENDDO
        ELSE IF( LTYPE.EQ.4 ) THEN
          VX=PARAMS(1)
          VY=PARAMS(2)
          VZ=PARAMS(3)
          VX=VX/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          VY=VY/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          VZ=VZ/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          DO I=1,nn
            VECT(3*I-2)=VAL*PLX(I)*RHO*VX
            VECT(3*I-1)=VAL*PLY(I)*RHO*VY
            VECT(3*I  )=VAL*PLZ(I)*RHO*VZ
          ENDDO
        ELSE IF( LTYPE.EQ.5 ) THEN
          DO I=1,nn
            VECT(3*I-2)=PLX(I)
            VECT(3*I-1)=PLY(I)
            VECT(3*I  )=PLZ(I)
          ENDDO
        END IF
      ENDIF

    end subroutine DL_C3

!> This subroutien calculate thermal loading
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C3                                 &
              (etype, nn, XX, YY, ZZ, TT, T0, gausses, &
               VECT, cdsys_ID, coords)
!----------------------------------------------------------------------*

    USE m_fstr
    USE mMechGauss
    USE m_MatMatrix
    USE m_utilities

!---------------------------------------------------------------------

    INTEGER(kind=kint), PARAMETER   :: ndof = 3
    INTEGER(kind=kint), INTENT(IN)  :: etype, nn
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)             !< status of qudrature points
    REAL(kind=kreal), INTENT(IN)    :: XX(nn), YY(nn), ZZ(nn)
    REAL(kind=kreal), INTENT(IN)    :: TT(nn),T0(nn)
    REAL(kind=kreal), INTENT(OUT)   :: VECT(nn*NDOF)
    INTEGER(kind=kint), INTENT(IN)  :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT) :: coords(3, 3)           !< variables to define matreial coordinate system

!---------------------------------------------------------------------

    REAL(kind=kreal) :: ALP, ALP0, D(6, 6), B(6, ndof*nn)
    REAL(kind=kreal) :: det, ecoord(3, nn)
    INTEGER(kind=kint) :: j, LX, serr
    REAL(kind=kreal) :: EPSTH(6),SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3)
    REAL(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    REAL(kind=kreal) :: wg, outa(1), ina(1)
    REAL(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    matlaniso = .FALSE.

    IF( cdsys_ID > 0 ) THEN   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .TRUE.
    END IF

    VECT(:) = 0.0D0

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)

    ! LOOP FOR INTEGRATION POINTS
    DO LX = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getShapeFunc( ETYPE, naturalcoord, H(1:nn) )
      CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv )

      IF( matlaniso ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      ! WEIGHT VALUE AT GAUSSIAN POINT
      wg = getWeight(etype, LX)*det
      B(1:6,1:nn*NDOF)=0.0D0
      DO J=1,nn
        B(1,3*J-2)=gderiv(j,1)
        B(2,3*J-1)=gderiv(j,2)
        B(3,3*J  )=gderiv(j,3)
        B(4,3*J-2)=gderiv(j,2)
        B(4,3*J-1)=gderiv(j,1)
        B(5,3*J-1)=gderiv(j,3)
        B(5,3*J  )=gderiv(j,2)
        B(6,3*J-2)=gderiv(j,3)
        B(6,3*J  )=gderiv(j,1)
      ENDDO

      TEMPC = DOT_PRODUCT( H(1:nn), TT(1:nn) )
      TEMP0 = DOT_PRODUCT( H(1:nn), T0(1:nn) )

      CALL MatlMatrix( gausses(LX), D3, D, 0.0D0, coordsys, tempc )

      ina(1) = TEMPC
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      END IF
      ina(1) = TEMP0
      IF( matlaniso  ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      END IF

      !**
      !** THERMAL strain
      !**
      IF( matlaniso ) THEN
        DO j=1,3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        END DO
        EPSTH(4:6) = 0.0D0
        CALL transformation(coordsys, tm)
        EPSTH(:) = MATMUL( EPSTH(:), tm  )      ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      ELSE
        THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      END IF

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = MATMUL( D(:, :), EPSTH(:) )
      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:nn*NDOF) = VECT(1:nn*NDOF)+MATMUL( SGM(:),B(:, :) )*wg

    END DO

   END SUBROUTINE TLOAD_C3


!> Update strain and stress inside element
!---------------------------------------------------------------------*
   SUBROUTINE UPDATE_C3                                       &
              (etype,nn,ecoord, u, ddu, cdsys_ID, coords, qf, &
               gausses, iter, tincr, TT, T0, TN)
!---------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities

    integer(kind=kint), INTENT(IN)    :: etype         !< \param [in] element type
    integer(kind=kint), INTENT(IN)    :: nn            !< \param [in] number of elemental nodes
    real(kind=kreal), INTENT(IN)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), INTENT(IN)      :: u(3, nn)      !< \param [in] nodal dislplacements
    real(kind=kreal), INTENT(IN)      :: ddu(3, nn)    !< \param [in] nodal displacement
    INTEGER(kind=kint), INTENT(IN)    :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT)   :: coords(3, 3)  !< variables to define matreial coordinate system
    real(kind=kreal), INTENT(OUT)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), INTENT(INOUT) :: gausses(:)    !< \param [out] status of qudrature points
    integer, intent(in)               :: iter
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    REAL(kind=kreal), INTENT(IN), optional :: TT(nn)   !< current temperature
    REAL(kind=kreal), INTENT(IN), optional :: T0(nn)   !< reference temperature
    REAL(kind=kreal), INTENT(IN), optional :: TN(nn)   !< reference temperature

    ! LCOAL VARIAVLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: D(6,6), B(6,ndof*nn), B1(6,ndof*nn), spfunc(nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gdispderiv(3,3), det, WG, ttc,tt0, ttn,outa(1)
    integer(kind=kint) :: i, j, k, LX, mtype, serr
    integer(kind=kint) :: isEp
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), mat(6,6), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn), elem(3,nn), elem1(3,nn), coordsys(3,3), tm(6,6)
    real(kind=kreal)   :: dstrain(6),dstress(6),dumstress(3,3),dum(3,3)
    real(kind=kreal)   :: alp, alp0, alpo(3),alpo0(3)
    logical            :: ierr, matlaniso

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:,:) = ecoord(:,:)
    totaldisp(:,:) = u(:,:)+ ddu(:,:)
    IF( flag == UPDATELAG ) THEN
      elem(:,:) = (0.5D0*ddu(:,:)+u(:,:) ) +ecoord(:,:)
      elem1(:,:) = (ddu(:,:)+u(:,:) ) +ecoord(:,:)
      ! elem = elem1
      totaldisp(:,:) = ddu(:,:)
    END IF

    matlaniso = .FALSE.
    IF( cdsys_ID > 0 ) THEN
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .true.
    END IF

    DO LX = 1, NumOfQuadPoints(etype)

      mtype = gausses(LX)%pMaterial%mtype

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      ! ========================================================
      ! UPDATE STRAIN and STRESS
      ! ========================================================

      if( isElastoplastic(mtype) .OR. mtype == NORTON )then
        isEp = 1
      else
        isEp = 0
      endif
      !gausses(LX)%pMaterial%mtype = ELASTIC

      EPSTH = 0.0D0
      IF( PRESENT(tt) .AND. PRESENT(t0) ) THEN
        CALL getShapeFunc(etype, naturalcoord, spfunc)
        ttc = DOT_PRODUCT(TT, spfunc)
        tt0 = DOT_PRODUCT(T0, spfunc)
        ttn = DOT_PRODUCT(TN, spfunc)
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, ttc, isEp )

        ina(1) = ttc
        IF( matlaniso ) THEN
          CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
          IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
        ELSE
          CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
          IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
          alp = outa(1)
        END IF
        ina(1) = tt0
        IF( matlaniso  ) THEN
          CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
          IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
        ELSE
          CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
          IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
          alp0 = outa(1)
        END IF
        IF( matlaniso ) THEN
          DO j=1,3
            EPSTH(j) = ALPO(j)*(ttc-ref_temp)-alpo0(j)*(tt0-ref_temp)
          END DO
          CALL transformation( coordsys(:, :), tm)
          EPSTH(:) = MATMUL( EPSTH(:), tm  ) ! to global coord
          EPSTH(4:6) = EPSTH(4:6)*2.0D0
        ELSE
          EPSTH(1:3)=ALP*(ttc-ref_temp)-alp0*(tt0-ref_temp)
        END IF

      ELSE

        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, isEp=isEp)

      END IF

      ! Small strain
      gdispderiv(1:ndof, 1:ndof) = MATMUL( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)   ! allright?

      IF( flag == INFINITE ) THEN

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
        gausses(LX)%stress(1:6) = MATMUL( D(1:6, 1:6), dstrain(1:6) )
        IF( isViscoelastic(mtype) .AND. tincr /= 0.0D0 ) THEN
          IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
          ELSE
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
          END IF
          gausses(LX)%stress = real(gausses(LX)%stress)
        END IF

      ELSE IF( flag == TOTALLAG ) THEN

        ! Green-Lagrange strain
        dstrain(1) = dstrain(1)+0.5d0*DOT_PRODUCT( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5d0*DOT_PRODUCT( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5d0*DOT_PRODUCT( gdispderiv(:, 3), gdispderiv(:, 3) )
        dstrain(4) = dstrain(4)+( gdispderiv(1, 1)*gdispderiv(1, 2)                                     &
                                 +gdispderiv(2, 1)*gdispderiv(2, 2)+gdispderiv(3, 1)*gdispderiv(3, 2) )
        dstrain(5) = dstrain(5)+( gdispderiv(1, 2)*gdispderiv(1, 3)                                     &
                                 +gdispderiv(2, 2)*gdispderiv(2, 3)+gdispderiv(3, 2)*gdispderiv(3, 3) )
        dstrain(6) = dstrain(6)+( gdispderiv(1, 1)*gdispderiv(1, 3)                                     &
                                 +gdispderiv(2, 1)*gdispderiv(2, 3)+gdispderiv(3, 1)*gdispderiv(3, 3) )

        IF( mtype == NEOHOOKE .OR. mtype == MOONEYRIVLIN .OR.  mtype == ARRUDABOYCE  .OR.   &
            mtype==USERELASTIC .OR. mtype==USERHYPERELASTIC .OR. mtype==USERMATERIAL ) THEN
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
        ELSE IF( ( isViscoelastic(mtype) .OR. mtype == NORTON ) .AND. tincr /= 0.0D0 ) THEN
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%pMaterial%mtype=mtype
          IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
          ELSE
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
          END IF
        ELSE
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%stress(1:6) = MATMUL( D(1:6, 1:6), dstrain(1:6) )
        END IF

      ELSE IF( flag == UPDATELAG ) THEN

        !  CALL GEOMAT_C3( gausses(LX)%stress, mat )
        !  D(:, :) = D(:, :)+mat(:, :)
        rot = 0.0D0
        rot(1, 2)= 0.5d0*(gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5d0*(gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5d0*(gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+ dstrain(1:6)+EPSTH(:)

        IF( isViscoelastic(mtype) .AND. tincr /= 0.0D0 ) THEN
          !(LX)%pMaterial%mtype = mtype
          IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, tt0 )
          ELSE
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
          END IF
        ELSE

          dstress = real( MATMUL( D(1:6,1:6), dstrain(1:6) ) )
          dumstress(1,1) = gausses(LX)%stress_bak(1)
          dumstress(2,2) = gausses(LX)%stress_bak(2)
          dumstress(3,3) = gausses(LX)%stress_bak(3)
          dumstress(1,2) = gausses(LX)%stress_bak(4);  dumstress(2,1)=dumstress(1,2)
          dumstress(2,3) = gausses(LX)%stress_bak(5);  dumstress(3,2)=dumstress(2,3)
          dumstress(3,1) = gausses(LX)%stress_bak(6);  dumstress(1,3)=dumstress(3,1)

          dum(:,:) = MATMUL( rot,dumstress ) -MATMUL( dumstress, rot )
          gausses(LX)%stress(1) = gausses(LX)%stress_bak(1)+dstress(1)+ dum(1,1)
          gausses(LX)%stress(2) = gausses(LX)%stress_bak(2)+dstress(2)+ dum(2,2)
          gausses(LX)%stress(3) = gausses(LX)%stress_bak(3)+dstress(3)+ dum(3,3)
          gausses(LX)%stress(4) = gausses(LX)%stress_bak(4)+dstress(4)+ dum(1,2)
          gausses(LX)%stress(5) = gausses(LX)%stress_bak(5)+dstress(5)+ dum(2,3)
          gausses(LX)%stress(6) = gausses(LX)%stress_bak(6)+dstress(6)+ dum(3,1)

          IF( mtype == USERMATERIAL ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
          ELSE IF( mtype == NORTON ) THEN
            !gausses(LX)%pMaterial%mtype = mtype
            IF( tincr /= 0.0D0 .AND. ANY( gausses(LX)%stress /= 0.0D0 ) ) THEN
              !gausses(LX)%pMaterial%mtype = mtype
              IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
                CALL StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr, ttc, ttn )
              ELSE
                CALL StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr )
              END IF
            END IF
          END IF
        END IF

      END IF

      !gausses(LX)%pMaterial%mtype = mtype

      IF( isElastoplastic(mtype) ) THEN
        IF( PRESENT(tt) ) THEN
          CALL BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain, &
                              gausses(LX)%istatus(1), gausses(LX)%fstatus, ttc )
        ELSE
          CALL BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain, &
                              gausses(LX)%istatus(1), gausses(LX)%fstatus )
        END IF
      END IF

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      DO J=1,nn
        B(1,3*j-2) = gderiv(j, 1)
        B(2,3*j-1) = gderiv(j, 2)
        B(3,3*j  ) = gderiv(j, 3)
        B(4,3*j-2) = gderiv(j, 2)
        B(4,3*j-1) = gderiv(j, 1)
        B(5,3*j-1) = gderiv(j, 3)
        B(5,3*j  ) = gderiv(j, 2)
        B(6,3*j-2) = gderiv(j, 3)
        B(6,3*j  ) = gderiv(j, 1)
      END DO

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      IF( flag == INFINITE ) THEN

      ELSE IF( flag == TOTALLAG ) THEN

        gdispderiv(1:ndof, 1:ndof) = MATMUL( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*ndof)=0.0D0
        DO j = 1,nn
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
        END DO
        ! BL = BL0 + BL1
        DO j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        END DO

      ELSE IF( flag == UPDATELAG ) THEN

        CALL getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
        B(1:6, 1:nn*ndof) = 0.0D0
        DO j = 1, nn
          B(1, 3*J-2) = gderiv(j, 1)
          B(2, 3*J-1) = gderiv(j, 2)
          B(3, 3*J  ) = gderiv(j, 3)
          B(4, 3*J-2) = gderiv(j, 2)
          B(4, 3*J-1) = gderiv(j, 1)
          B(5, 3*J-1) = gderiv(j, 3)
          B(5, 3*J  ) = gderiv(j, 2)
          B(6, 3*J-2) = gderiv(j, 3)
          B(6, 3*J  ) = gderiv(j, 1)
        END DO

      END IF

      ! calculate the Internal Force
      WG=getWeight( etype, LX )*DET
      qf(1:nn*ndof)                                                          &
      = qf(1:nn*ndof)+MATMUL( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG

    END DO

   END SUBROUTINE UPDATE_C3


!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C3                            &
              (etype, nn, XX, YY, ZZ, TT, T0, edisp, &
               gausses, cdsys_ID, coords)
!----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    USE m_fstr
    USE mMechGauss
    USE m_MatMatrix
    USE m_utilities

!---------------------------------------------------------------------

    INTEGER(kind=kint), PARAMETER :: ndof = 3
    ! I/F VARIABLES
    INTEGER(kind=kint), INTENT(IN)    :: etype, nn
    REAL(kind=kreal), INTENT(IN)      :: XX(nn), YY(nn), ZZ(nn)
    REAL(kind=kreal), INTENT(IN)      :: TT(nn), T0(nn), edisp(nn*ndof)
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)
    INTEGER(kind=kint), INTENT(IN)    :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT)   :: coords(3, 3)

!---------------------------------------------------------------------

    ! LOCAL VARIABLES
    REAL(kind=kreal) :: ALP, ALP0,D(6, 6), B(6, ndof*nn)
    REAL(kind=kreal) :: det, ecoord(3, nn)
    INTEGER(kind=kint) :: j, k, IC, serr
    REAL(kind=kreal) :: EPSA(6), EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3)
    REAL(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    REAL(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6),outa(1),ina(1)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    matlaniso = .FALSE.

    IF( cdsys_ID > 0 ) THEN   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .TRUE.
    END IF

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)

    ! LOOP FOR INTEGRATION POINTS
    DO IC = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, IC, naturalCoord(:) )
      CALL getShapeFunc( etype, naturalcoord, H(1:nn) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      IF( matlaniso ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      B(1:6, 1:nn*ndof) = 0.0D0
      DO j = 1, nn
        B(1,3*j-2) = gderiv(j,1)
        B(2,3*j-1) = gderiv(j,2)
        B(3,3*j  ) = gderiv(j,3)
        B(4,3*j-2) = gderiv(j,2)
        B(4,3*j-1) = gderiv(j,1)
        B(5,3*j-1) = gderiv(j,3)
        B(5,3*j  ) = gderiv(j,2)
        B(6,3*j-2) = gderiv(j,3)
        B(6,3*j  ) = gderiv(j,1)
      END DO

      TEMPC = DOT_PRODUCT( H(1:nn), TT(1:nn) )
      TEMP0 = DOT_PRODUCT( H(1:nn), T0(1:nn) )

      CALL MatlMatrix( gausses(IC), D3, D, 0.0D0, coordsys, TEMPC )

      ina(1) = TEMPC
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      END IF
      ina(1) = TEMP0
      IF( matlaniso  ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      END IF

      !**
      !** THERMAL strain
      !**
      IF( matlaniso ) THEN
        DO j = 1,3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        END DO
        EPSTH(4:6) = 0.0D0
        CALL transformation(coordsys, tm)
        EPSTH(:) = MATMUL( EPSTH(:), tm  )      ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      ELSE
        THERMAL_EPS = ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      END IF
      !**
      !** SET EPS  {e}=[B]{u}
      !**
      EPSA(1:6) = MATMUL( B(1:6, :), edisp )
      !**
      !** SET SGM  {S}=[D]{e}
      !**
      DO j = 1, 6
        SGM(j) = 0.0D0
        DO k = 1, 6
          SGM(j) = SGM(j)+D(j, k)*( EPSA(k)-EPSTH(k) )
        END DO
      END DO
      !**
      !** Adding stress in each gauss points
      !**
      gausses(IC)%strain(1:6) = EPSA(1:6)
      gausses(IC)%stress(1:6) = SGM(1:6)

    END DO

   END SUBROUTINE UpdateST_C3


!----------------------------------------------------------------------*
   SUBROUTINE NodalStress_C3(etype, nn, gausses, ndstrain, ndstress)
!----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    USE mMechGauss

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN) :: etype, nn
    TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
    REAL(kind=kreal), INTENT(OUT)  :: ndstrain(nn,6)
    REAL(kind=kreal), INTENT(OUT)  :: ndstress(nn,6)

!---------------------------------------------------------------------

    INTEGER :: i, ic
    REAL(kind=kreal) :: TEMP(12)

!---------------------------------------------------------------------

    TEMP(:) = 0.0D0

    IC = NumOfQuadPoints(etype)

    DO i = 1, IC
      TEMP(1:6)  = TEMP(1:6) +gausses(i)%strain(1:6)
      TEMP(7:12) = TEMP(7:12)+gausses(i)%stress(1:6)
    END DO

    TEMP(1:12) = TEMP(1:12)/IC

    FORALL( i=1:nn )
      ndstrain(i, 1:6) = TEMP(1:6)
      ndstress(i, 1:6) = TEMP(7:12)
    END FORALL

   END SUBROUTINE NodalStress_C3


!----------------------------------------------------------------------*
   SUBROUTINE ElementStress_C3(etype, gausses, strain, stress)
!----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    USE mMechGauss

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN) :: etype
    TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
    REAL(kind=kreal), INTENT(OUT)  :: strain(6)
    REAL(kind=kreal), INTENT(OUT)  :: stress(6)

!---------------------------------------------------------------------

    INTEGER :: i, ic

!---------------------------------------------------------------------

    strain(:) = 0.0D0; stress(:) = 0.0D0

    IC = NumOfQuadPoints(etype)

    DO i = 1, IC
      strain(:) = strain(:)+gausses(i)%strain(1:6)
      stress(:) = stress(:)+gausses(i)%stress(1:6)
    ENDDO

    strain(:) = strain(:)/IC
    stress(:) = stress(:)/IC

   END SUBROUTINE ElementStress_C3


!> Volume of element
!----------------------------------------------------------------------*
   REAL(kind=kreal) FUNCTION VOLUME_C3(etype, nn, XX, YY, ZZ)
!----------------------------------------------------------------------*

    INTEGER(kind=kint), INTENT(IN) :: etype, nn
    REAL(kind=kreal), INTENT(IN)   :: XX(:), YY(:), ZZ(:)

!---------------------------------------------------------------------

    REAL(kind=kreal) :: XJ(3, 3), det, wg
    INTEGER(kind=kint) :: LX, i
    REAL(kind=kreal) :: localcoord(3), deriv(nn, 3)

!---------------------------------------------------------------------

    VOLUME_C3 = 0.0D0

    ! LOOP FOR INTEGRATION POINTS
    DO LX = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint(etype, LX, localcoord)
      CALL getShapeDeriv(etype, localcoord, deriv)

      ! JACOBI MATRIX
      XJ(1, 1:3)= MATMUL( xx(1:nn), deriv(1:nn,1:3) )
      XJ(2, 1:3)= MATMUL( yy(1:nn), deriv(1:nn,1:3) )
      XJ(3, 1:3)= MATMUL( zz(1:nn), deriv(1:nn,1:3) )

      ! DETERMINANT OF JACOBIAN
      det = XJ(1, 1)*XJ(2, 2)*XJ(3, 3) &
           +XJ(2, 1)*XJ(3, 2)*XJ(1, 3) &
           +XJ(3, 1)*XJ(1, 2)*XJ(2, 3) &
           -XJ(3, 1)*XJ(2, 2)*XJ(1, 3) &
           -XJ(2, 1)*XJ(1, 2)*XJ(3, 3) &
           -XJ(1, 1)*XJ(3, 2)*XJ(2, 3)

      VOLUME_C3 = VOLUME_C3+getWeight(etype, LX)*det

    END DO

   END FUNCTION VOLUME_C3


END MODULE m_static_LIB_3d
