!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
module m_static_LIB_3d
   use hecmw, only : kint, kreal
   use elementInfo
   implicit none

   contains

!=====================================================================*
!>  This subroutine calculate stiff matrix of general solid elements
!
!>  \author     X. YUAN, K. SATO (AdavanceSoft)
!>  \date       2009/08/03
!>  \version    0.00

  SUBROUTINE STF_C3( etype,nn,ecoord,gausses,stiff, tincr, coords,u ,temperature)
    USE mMechGauss
    use m_MatMatrix
    use m_common_struct
    INTEGER(kind=kint), INTENT(IN)  :: etype               !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                  !< number of elemental nodes
    REAL(kind=kreal),   INTENT(IN)  :: ecoord(3,nn)        !< coordinates of elemental nodes
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)          !< stiff matrix
    real(kind=kreal),    intent(in) :: tincr               !< time increment
    REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)     !< variables to define matreial coordinate system
    REAL(kind=kreal),   INTENT(IN), optional :: temperature(nn)     !< temperature
    REAL(kind=kreal),   INTENT(IN), optional :: u(:,:)     !< nodal displacemwent
!
    INTEGER             :: flag
    INTEGER(kind=kint), PARAMETER :: NDOF=3
    REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
    REAL(kind=kreal) gderiv(NN,3),stress(6),mat(6,6)
    REAL(kind=kreal) DET,WG
    INTEGER(kind=kint) I,J,LX, cdsys_ID, serr
    REAL(kind=kreal) temp, naturalCoord(3)
    REAL(kind=kreal) spfunc(nn), gdispderiv(3,3)
    REAL(kind=kreal) B1(6,NDOF*NN), coordsys(3,3)
    REAL(kind=kreal) Smat(9,9), elem(3,nn)
    REAL(kind=kreal) BN(9,NDOF*NN), SBN(9,NDOF*NN)

    stiff(:,:) = 0.d0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag=INFINITE    ! enforce to infinite deformation analysis
    elem(:,:) = ecoord(:,:)
    if( flag == UPDATELAG ) elem(:,:) = ecoord(:,:) + u(:,:)
	
    cdsys_ID = gausses(1)%pMaterial%cdsys_ID

    DO LX=1,NumOfQuadPoints(etype)
      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )
	  
      if( cdsys_ID>0 ) then
        call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
        if( serr==-1 ) stop "Fail to setup local coordinate"
        if( serr==-2 ) then
           write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        endif
      endif
	  
      if( present(temperature) ) then
        CALL getShapeFunc( etype, naturalcoord, spfunc )
        temp = dot_product( temperature, spfunc )
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, temp )
      else
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys )
      endif

      WG=getWeight( etype, LX )*DET
      B(1:6,1:NN*NDOF)=0.d0
      DO J=1,NN
        B(1,3*J-2)=gderiv(J,1)
        B(2,3*J-1)=gderiv(J,2)
        B(3,3*J  )=gderiv(J,3)
        B(4,3*J-2)=gderiv(J,2)
        B(4,3*J-1)=gderiv(J,1)
        B(5,3*J-1)=gderiv(J,3)
        B(5,3*J  )=gderiv(J,2)
        B(6,3*J-2)=gderiv(J,3)
        B(6,3*J  )=gderiv(J,1)
      ENDDO
!
!    calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag==TOTALLAG ) then
!      ---dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof,1:ndof) = matmul( u(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
        B1(1:6,1:NN*NDOF)=0.d0
        do j=1,nn
          B1(1,3*J-2)=gdispderiv(1,1)*gderiv(J,1)
          B1(1,3*J-1)=gdispderiv(2,1)*gderiv(J,1)
          B1(1,3*J  )=gdispderiv(3,1)*gderiv(J,1)
          B1(2,3*J-2)=gdispderiv(1,2)*gderiv(J,2)
          B1(2,3*J-1)=gdispderiv(2,2)*gderiv(J,2)
          B1(2,3*J  )=gdispderiv(3,2)*gderiv(J,2)
          B1(3,3*J-2)=gdispderiv(1,3)*gderiv(J,3)
          B1(3,3*J-1)=gdispderiv(2,3)*gderiv(J,3)
          B1(3,3*J  )=gdispderiv(3,3)*gderiv(J,3)
          B1(4,3*J-2)=gdispderiv(1,2)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,2)
          B1(4,3*J-1)=gdispderiv(2,2)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,2)
          B1(4,3*J  )=gdispderiv(3,2)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,2)
          B1(5,3*J-2)=gdispderiv(1,2)*gderiv(J,3)+gdispderiv(1,3)*gderiv(J,2)
          B1(5,3*J-1)=gdispderiv(2,2)*gderiv(J,3)+gdispderiv(2,3)*gderiv(J,2)
          B1(5,3*J  )=gdispderiv(3,2)*gderiv(J,3)+gdispderiv(3,3)*gderiv(J,2)
          B1(6,3*J-2)=gdispderiv(1,3)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,3)
          B1(6,3*J-1)=gdispderiv(2,3)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,3)
          B1(6,3*J  )=gdispderiv(3,3)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,3)
        enddo
!    ---BL = BL0 + BL1
        do j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        enddo
      endif

      DB(1:6,1:NN*NDOF) = matmul( D,B(1:6,1:NN*NDOF) )
      forall( i=1:nn*ndof, j=1:nn*ndof )
        stiff(i,j) = stiff(i,j)+dot_product( B(:,i), DB(:,j) )*WG
      end forall

!    calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag==TOTALLAG .or. flag==UPDATELAG ) then
        stress(1:6)=gausses(LX)%stress
        BN(1:9,1:nn*ndof) = 0.d0
        do j=1,nn
          BN(1,3*J-2) = gderiv(J,1)
          BN(2,3*J-1) = gderiv(J,1)
          BN(3,3*J  ) = gderiv(J,1)
          BN(4,3*J-2) = gderiv(J,2)
          BN(5,3*J-1) = gderiv(J,2)
          BN(6,3*J  ) = gderiv(J,2)
          BN(7,3*J-2) = gderiv(J,3)
          BN(8,3*J-1) = gderiv(J,3)
          BN(9,3*J  ) = gderiv(J,3)
        enddo
        Smat(:,:) = 0.d0
        do j=1,3
          Smat(j  ,j  ) = stress(1)
          Smat(j  ,j+3) = stress(4)
          Smat(j  ,j+6) = stress(6)
          Smat(j+3,j  ) = stress(4)
          Smat(j+3,j+3) = stress(2)
          Smat(j+3,j+6) = stress(5)
          Smat(j+6,j  ) = stress(6)
          Smat(j+6,j+3) = stress(5)
          Smat(j+6,j+6) = stress(3)
        enddo
        SBN(1:9,1:nn*ndof) = matmul( Smat(1:9,1:9), BN(1:9,1:nn*ndof) )
        forall( i=1:nn*ndof, j=1:nn*ndof )
          stiff(i,j) = stiff(i,j)+dot_product( BN(:,i), SBN(:,j) )*WG
        end forall

      endif
    ENDDO      ! gauss roop

  end subroutine STF_C3
!
!
!> Distrubuted external load
!----------------------------------------------------------------------*
   SUBROUTINE DL_C3(ETYPE,NN,XX,YY,ZZ,RHO,LTYPE,PARAMS,VECT,NSIZE)
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
      INTEGER(kind=kint), INTENT(IN)  :: ETYPE, NN
      REAL(kind=kreal), INTENT(IN)    :: XX(:),YY(:),ZZ(:)
      REAL(kind=kreal), INTENT(IN)    :: PARAMS(0:6)
      REAL(kind=kreal), INTENT(INOUT) :: VECT(:)
      REAL(kind=kreal) RHO
      INTEGER(kind=kint) LTYPE,NSIZE
! LOCAL VARIABLES
      INTEGER(kind=kint) NDOF
      PARAMETER(NDOF=3)
      REAL(kind=kreal) H(NN)
      REAL(kind=kreal) PLX(NN),PLY(NN),PLZ(NN)
      REAL(kind=kreal) XJ(3,3),DET,WG
      INTEGER(kind=kint) IVOL,ISUF
      INTEGER(kind=kint) NOD(NN)
      INTEGER(kind=kint) IG2,LX,I ,SURTYPE,NSUR
      REAL(kind=kreal) VX,VY,VZ,XCOD,YCOD,ZCOD
      REAL(kind=kreal) AX,AY,AZ,RX,RY,RZ,HX,HY,HZ,VAL
      REAL(kind=kreal) PHX,PHY,PHZ
      REAL(kind=kreal) COEFX,COEFY,COEFZ
      REAL(kind=kreal) normal(3),localcoord(3),elecoord(3,NN),deriv(NN,3)
!
! SET VALUE
!
      VAL=PARAMS(0)
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
      NSIZE=NN*NDOF
      VECT(1:NSIZE)=0.d0
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
        PLX(:)=0.d0
        PLY(:)=0.d0
        PLZ(:)=0.d0
! LOOP FOR INTEGRATION POINTS
        DO  LX=1,NumOfQuadPoints( ETYPE )
              CALL getQuadPoint( ETYPE, LX, localcoord )
              CALL getShapeFunc( ETYPE, localcoord, H(1:NN) )
              CALL getShapeDeriv( ETYPE, localcoord, deriv )
!  JACOBI MATRIX
              XJ(1,1:3)= matmul( xx(1:NN), deriv(1:NN,1:3) )
              XJ(2,1:3)= matmul( yy(1:NN), deriv(1:NN,1:3) )
              XJ(3,1:3)= matmul( zz(1:NN), deriv(1:NN,1:3) )
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
                XCOD=DOT_PRODUCT( H(1:NN),XX(1:NN) )
                YCOD=DOT_PRODUCT( H(1:NN),YY(1:NN) )
                ZCOD=DOT_PRODUCT( H(1:NN),ZZ(1:NN) )
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
              DO I=1,NN
                PLX(I)=PLX(I)+H(I)*WG*COEFX
                PLY(I)=PLY(I)+H(I)*WG*COEFY
                PLZ(I)=PLZ(I)+H(I)*WG*COEFZ
              ENDDO
        ENDDO
        IF( LTYPE.EQ.1) THEN
          DO I=1,NN
            VECT(3*I-2)=VAL*PLX(I)
          ENDDO
        ELSE IF( LTYPE.EQ.2 ) THEN
          DO I=1,NN
            VECT(3*I-1)=VAL*PLY(I)
          ENDDO
        ELSE IF( LTYPE.EQ.3 ) THEN
          DO I=1,NN
            VECT(3*I  )=VAL*PLZ(I)
          ENDDO
        ELSE IF( LTYPE.EQ.4 ) THEN
          VX=PARAMS(1)
          VY=PARAMS(2)
          VZ=PARAMS(3)
          VX=VX/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          VY=VY/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          VZ=VZ/SQRT(PARAMS(1)**2+PARAMS(2)**2+PARAMS(3)**2)
          DO I=1,NN
            VECT(3*I-2)=VAL*PLX(I)*RHO*VX
            VECT(3*I-1)=VAL*PLY(I)*RHO*VY
            VECT(3*I  )=VAL*PLZ(I)*RHO*VZ
          ENDDO
        ELSE IF( LTYPE.EQ.5 ) THEN
          DO I=1,NN
            VECT(3*I-2)=PLX(I)
            VECT(3*I-1)=PLY(I)
            VECT(3*I  )=PLZ(I)
          ENDDO
        END IF
      ENDIF

   end subroutine DL_C3

!> This subroutien calculate thermal loading
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C3(ETYPE,NN,XX,YY,ZZ,TT,T0,gausses,VECT,coords)
!----------------------------------------------------------------------*
      use m_fstr
      USE mMechGauss
      use m_MatMatrix
      use m_utilities
      INTEGER(kind=kint), PARAMETER   :: NDOF=3
      INTEGER(kind=kint), INTENT(IN)  :: ETYPE,NN
      TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
      REAL(kind=kreal), INTENT(IN)    :: XX(NN),YY(NN),ZZ(NN),TT(NN),T0(NN)
      REAL(kind=kreal), INTENT(OUT)   :: VECT(NN*NDOF)
      REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)       !< variables to define matreial coordinate system

      REAL(kind=kreal) ALP,ALP0,D(6,6),B(6,NDOF*NN)
      REAL(kind=kreal) DET,ecoord(3,NN)
      INTEGER(kind=kint) J,LX, cdsys_ID, serr
      REAL(kind=kreal) EPSTH(6),SGM(6),H(NN),alpo(3),alpo0(3), coordsys(3,3)
      REAL(kind=kreal) naturalcoord(3),gderiv(NN,3)
      REAL(kind=kreal) WG, outa(1), ina(1)
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS, tm(6,6)
      logical   :: ierr, matlaniso
	  
        matlaniso = .false.
        cdsys_ID = gausses(1)%pMaterial%cdsys_ID
        if( cdsys_ID>0 ) then   ! cannot define aniso exapansion when no local coord defined
          ina = TT(1)
          call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
          if( .not. ierr ) matlaniso = .true.
        endif

      VECT(:)=0.d0
	  
      ecoord(1,:)=XX(:)
      ecoord(2,:)=YY(:)
      ecoord(3,:)=ZZ(:)
! LOOP FOR INTEGRATION POINTS
      DO LX=1,NumOfQuadPoints(etype)  
        CALL getQuadPoint( etype, LX, naturalCoord(:) )
        CALL getShapeFunc( ETYPE, naturalcoord, H(1:NN) )
        CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv )
		
        if( matlaniso ) then
          call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
          if( serr==-1 ) stop "Fail to setup local coordinate"
          if( serr==-2 ) then
            write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
          endif
        endif

!  WEIGHT VALUE AT GAUSSIAN POINT
            WG=getWeight( etype, LX )*DET
            B(1:6,1:NN*NDOF)=0.d0
            DO J=1,NN
              B(1,3*J-2)=gderiv(J,1)
              B(2,3*J-1)=gderiv(J,2)
              B(3,3*J  )=gderiv(J,3)
              B(4,3*J-2)=gderiv(J,2)
              B(4,3*J-1)=gderiv(J,1)
              B(5,3*J-1)=gderiv(J,3)
              B(5,3*J  )=gderiv(J,2)
              B(6,3*J-2)=gderiv(J,3)
              B(6,3*J  )=gderiv(J,1)
            ENDDO

            TEMPC=DOT_PRODUCT( H(1:NN),TT(1:NN) )
            TEMP0=DOT_PRODUCT( H(1:NN),T0(1:NN) )
			
        CALL MatlMatrix( gausses(LX), D3, D, 0.d0, coordsys, tempc )
		
            ina(1) = TEMPC
            if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
               call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
            endif
            ina(1) = TEMP0
            if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
              call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
              alp0 = outa(1)
            endif

!**
!** THERMAL strain
!**
            if( matlaniso ) then
              do j=1,3
                EPSTH(j)=ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
              enddo
              EPSTH(4:6)=0.d0
              call transformation(coordsys, tm) 
              EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord  
              EPSTH(4:6)=EPSTH(4:6)*2.d0 
            else
              THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
              EPSTH(1:3)=THERMAL_EPS
              EPSTH(4:6)=0.d0
            endif
			
!**
!** SET SGM  {s}=[D]{e}
!**
            SGM(:)=MATMUL( D(:,:),EPSTH(:) )
!**
!** CALCULATE LOAD {F}=[B]T{e}
!**
            VECT(1:NN*NDOF)=VECT(1:NN*NDOF)+MATMUL( SGM(:),B(:,:))*WG
      ENDDO
   end subroutine TLOAD_C3
   
!
!> Update strain and stress inside element
!---------------------------------------------------------------------*
  SUBROUTINE UPDATE_C3( etype,nn,ecoord, u, ddu, coords, qf ,gausses, iter, tincr, TT,T0, TN  )
!---------------------------------------------------------------------*
    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: ecoord(3,nn)    !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: u(3,nn)         !< \param [in] nodal dislplacements 
    real(kind=kreal),   INTENT(IN)     :: ddu(3,nn)       !< \param [in] nodal displacement
    REAL(kind=kreal), INTENT(INOUT)    :: coords(3,3)     !< variables to define matreial coordinate system
    real(kind=kreal),   INTENT(OUT)    :: qf(nn*3)        !< \param [out] Internal Force    
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points
    integer, intent(in)                :: iter
	real(kind=kreal),    intent(in)    :: tincr           !< time increment
    REAL(kind=kreal),   INTENT(IN), optional :: TT(nn)    !< current temperature
    REAL(kind=kreal),   INTENT(IN), optional :: T0(nn)    !< reference temperature
    REAL(kind=kreal),   INTENT(IN), optional :: TN(nn)    !< reference temperature
! LCOAL VARIAVLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof=3
    real(kind=kreal)   :: D(6,6), B(6,ndof*nn), B1(6,ndof*nn), spfunc(nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gdispderiv(3,3), det, WG, ttc,tt0, ttn,outa(1)
    integer(kind=kint) :: i, j, k, LX, mtype, cdsys_ID, serr
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), mat(6,6), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn), elem(3,nn), elem1(3,nn), coordsys(3,3), tm(6,6)
    real(kind=kreal)   :: dstrain(6),dstress(6),dumstress(3,3),dum(3,3)
    real(kind=kreal)   :: alp, alp0, alpo(3),alpo0(3)
    logical            :: ierr, matlaniso

    qf(:)              = 0.d0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:,:) = ecoord(:,:)
    totaldisp(:,:) = u(:,:)+ ddu(:,:)
    if( flag == UPDATELAG ) then
      elem(:,:) = (0.5D0*ddu(:,:)+u(:,:) ) +ecoord(:,:)
      elem1(:,:) = (ddu(:,:)+u(:,:) ) +ecoord(:,:) 
    !  elem = elem1
      totaldisp(:,:) = ddu(:,:)
    endif
	
    cdsys_ID = gausses(1)%pMaterial%cdsys_ID
    matlaniso = .false.
    if( cdsys_ID>0 ) then
       ina = TT(1)
       call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
       if( .not. ierr ) matlaniso = .true.
    endif

    do LX=1, NumOfQuadPoints(etype)
      mtype = gausses(LX)%pMaterial%mtype
      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )
	  
      if( cdsys_ID>0 ) then
        call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
        if( serr==-1 ) stop "Fail to setup local coordinate"
        if( serr==-2 ) then
           write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        endif
      endif

!
! ========================================================
!     UPDATE STRAIN and STRESS
! ========================================================
      

        if( isElastoplastic(gausses(LX)%pMaterial%mtype)    &
         .or. gausses(LX)%pMaterial%mtype==NORTON )         &   
           gausses(LX)%pMaterial%mtype = ELASTIC

      EPSTH = 0.d0		   
      if( present(tt) .and. present(t0) ) then
        CALL getShapeFunc( etype, naturalcoord, spfunc )
        ttc = dot_product( TT, spfunc )
        tt0 = dot_product( T0, spfunc )
        ttn = dot_product( TN, spfunc )
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, ttc )
         
          ina(1) = ttc
          if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
          else
               call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
          endif
          ina(1) = tt0
          if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
          else
              call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
              if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
              alp0 = outa(1)
          endif
          if( matlaniso ) then
              do j=1,3
                EPSTH(j)=ALPO(j)*(ttc-ref_temp)-alpo0(j)*(tt0-ref_temp)
              enddo
              call transformation(coordsys(:,:), tm)  
              EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord 
              EPSTH(4:6)=EPSTH(4:6)*2.d0 
          else
              EPSTH(1:3)=ALP*(ttc-ref_temp)-alp0*(tt0-ref_temp)
          endif

      else
        CALL MatlMatrix( gausses(LX), D3, D, tincr , coordsys)
      endif

!       Small strain
        gdispderiv(1:ndof,1:ndof) = matmul( totaldisp(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
        dstrain(1) = gdispderiv(1,1) 
        dstrain(2) = gdispderiv(2,2) 
        dstrain(3) = gdispderiv(3,3)
        dstrain(4) = ( gdispderiv(1,2)+gdispderiv(2,1) )
        dstrain(5) = ( gdispderiv(2,3)+gdispderiv(3,2) )
        dstrain(6) = ( gdispderiv(3,1)+gdispderiv(1,3) )
        dstrain(:) = dstrain(:)-EPSTH(:)   ! allright?
		
        if( flag==INFINITE ) then
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%stress(1:6) = matmul( D(1:6,1:6), dstrain(1:6) )
          if( isViscoelastic(mtype) .and. tincr/=0.d0 ) then
              if( present(TT) .and. present(T0) ) then
                call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
              else
                call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
              endif
              gausses(LX)%stress = real(gausses(LX)%stress)
          endif	
        else if( flag==TOTALLAG ) then
!       Green-Lagrange strain
          dstrain(1) = dstrain(1) +0.5d0*dot_product( gdispderiv(:,1), gdispderiv(:,1) )
          dstrain(2) = dstrain(2) +0.5d0*dot_product( gdispderiv(:,2), gdispderiv(:,2) )
          dstrain(3) = dstrain(3) +0.5d0*dot_product( gdispderiv(:,3), gdispderiv(:,3) )
          dstrain(4) = dstrain(4)+ (gdispderiv(1,1)*gdispderiv(1,2)    &
                      +gdispderiv(2,1)*gdispderiv(2,2)+gdispderiv(3,1)*gdispderiv(3,2) )
          dstrain(5) = dstrain(5)+ (gdispderiv(1,2)*gdispderiv(1,3)    &
                      +gdispderiv(2,2)*gdispderiv(2,3)+gdispderiv(3,2)*gdispderiv(3,3) )
          dstrain(6) = dstrain(6)+ (gdispderiv(1,1)*gdispderiv(1,3)    &
                      +gdispderiv(2,1)*gdispderiv(2,3)+gdispderiv(3,1)*gdispderiv(3,3) )
          			  
          if( mtype == NEOHOOKE .or. mtype == MOONEYRIVLIN .or.  mtype == ARRUDABOYCE  .or.      &
            mtype==USERELASTIC .or. mtype==USERHYPERELASTIC .or. mtype==USERMATERIAL ) then
             gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
          else if( ( isViscoelastic(mtype) &
              .or. mtype==NORTON) .and. tincr/=0.d0 ) then
            gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            gausses(LX)%pMaterial%mtype=mtype
            if( present(TT) .and. present(T0) ) then
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
            else
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
            endif
          else
            gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            gausses(LX)%stress(1:6) = matmul( D(1:6,1:6), dstrain(1:6) )
          endif
        else if( flag==UPDATELAG ) then
        !  call GEOMAT_C3( gausses(LX)%stress, mat )
        !  D(:,:) = D(:,:)+mat(:,:)
          rot =0.d0
          rot(1,2)= 0.5d0*(gdispderiv(1,2)-gdispderiv(2,1) );  rot(2,1) = -rot(1,2)
          rot(2,3)= 0.5d0*(gdispderiv(2,3)-gdispderiv(3,2) );  rot(3,2) = -rot(2,3)
          rot(1,3)= 0.5d0*(gdispderiv(1,3)-gdispderiv(3,1) );  rot(3,1) = -rot(1,3)

          gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+ dstrain(1:6)+EPSTH(:)
		  
          if( isViscoelastic(mtype) .and. tincr/=0.d0 ) then
            gausses(LX)%pMaterial%mtype=mtype
            if( present(TT) .and. present(T0) ) then
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, tt0 )
            else
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
            endif
          else	

          dstress = real( matmul( D(1:6,1:6), dstrain(1:6) ) )
          dumstress(1,1) = gausses(LX)%stress_bak(1)
          dumstress(2,2) = gausses(LX)%stress_bak(2)
          dumstress(3,3) = gausses(LX)%stress_bak(3)
          dumstress(1,2) = gausses(LX)%stress_bak(4);  dumstress(2,1)=dumstress(1,2)
          dumstress(2,3) = gausses(LX)%stress_bak(5);  dumstress(3,2)=dumstress(2,3)
          dumstress(3,1) = gausses(LX)%stress_bak(6);  dumstress(1,3)=dumstress(3,1)
    
          dum(:,:) = matmul( rot,dumstress ) -matmul( dumstress, rot )
          gausses(LX)%stress(1) = gausses(LX)%stress_bak(1)+dstress(1)+ dum(1,1)
          gausses(LX)%stress(2) = gausses(LX)%stress_bak(2)+dstress(2)+ dum(2,2)
          gausses(LX)%stress(3) = gausses(LX)%stress_bak(3)+dstress(3)+ dum(3,3)
          gausses(LX)%stress(4) = gausses(LX)%stress_bak(4)+dstress(4)+ dum(1,2)
          gausses(LX)%stress(5) = gausses(LX)%stress_bak(5)+dstress(5)+ dum(2,3)
          gausses(LX)%stress(6) = gausses(LX)%stress_bak(6)+dstress(6)+ dum(3,1)
		  
          if( mtype==USERMATERIAL ) then
            call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
          else if( mtype==NORTON ) then
            gausses(LX)%pMaterial%mtype=mtype
		    if( tincr/=0.d0 .and. any(gausses(LX)%stress/=0.d0) ) then
              gausses(LX)%pMaterial%mtype=mtype
              if( present(TT) .and. present(T0) ) then
                call StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr, ttc, ttn )
              else
                call StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr )
              endif
            endif
          endif
          endif
        endif
        gausses(LX)%pMaterial%mtype=mtype
        if( isElastoplastic(mtype) ) then
          if( present(tt) ) then
            call BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain  &
             , gausses(LX)%istatus(1), gausses(LX)%fstatus, ttc )
          else
            call BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain  &
             , gausses(LX)%istatus(1), gausses(LX)%fstatus )
          endif
        endif
!
!
! ========================================================
!     calculate the internal force ( equivalent nodal force )
! ========================================================
!       Small strain
      B(1:6, 1:nn*ndof) = 0.0d0
      DO J=1,NN
          B(1,3*J-2) = gderiv(J,1)
          B(2,3*J-1) = gderiv(J,2)
          B(3,3*J  ) = gderiv(J,3)
          B(4,3*J-2) = gderiv(J,2)
          B(4,3*J-1) = gderiv(J,1)
          B(5,3*J-1) = gderiv(J,3)
          B(5,3*J  ) = gderiv(J,2)
          B(6,3*J-2) = gderiv(J,3)
          B(6,3*J  ) = gderiv(J,1)
      ENDDO
!     calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag==INFINITE ) then
      else if( flag==TOTALLAG ) then
        gdispderiv(1:ndof,1:ndof) = matmul( totaldisp(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
        B1(1:6,1:NN*NDOF)=0.d0
        do j=1,nn
          B1(1,3*J-2) = gdispderiv(1,1)*gderiv(J,1)
          B1(1,3*J-1) = gdispderiv(2,1)*gderiv(J,1)
          B1(1,3*J  ) = gdispderiv(3,1)*gderiv(J,1)
          B1(2,3*J-2) = gdispderiv(1,2)*gderiv(J,2)
          B1(2,3*J-1) = gdispderiv(2,2)*gderiv(J,2)
          B1(2,3*J  ) = gdispderiv(3,2)*gderiv(J,2)
          B1(3,3*J-2) = gdispderiv(1,3)*gderiv(J,3)
          B1(3,3*J-1) = gdispderiv(2,3)*gderiv(J,3)
          B1(3,3*J  ) = gdispderiv(3,3)*gderiv(J,3)
          B1(4,3*J-2) = gdispderiv(1,2)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,2)
          B1(4,3*J-1) = gdispderiv(2,2)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,2)
          B1(4,3*J  ) = gdispderiv(3,2)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,2)
          B1(5,3*J-2) = gdispderiv(1,2)*gderiv(J,3)+gdispderiv(1,3)*gderiv(J,2)
          B1(5,3*J-1) = gdispderiv(2,2)*gderiv(J,3)+gdispderiv(2,3)*gderiv(J,2)
          B1(5,3*J  ) = gdispderiv(3,2)*gderiv(J,3)+gdispderiv(3,3)*gderiv(J,2)
          B1(6,3*J-2) = gdispderiv(1,3)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,3)
          B1(6,3*J-1) = gdispderiv(2,3)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,3)
          B1(6,3*J  ) = gdispderiv(3,3)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,3)
        enddo
!       -BL = BL0 + BL1
        do j=1,nn*ndof
            B(:,j) = B(:,j)+B1(:,j)
        enddo
!
      else if( flag == UPDATELAG ) then
        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv )
        B(1:6, 1:nn*ndof) = 0.0d0
        DO J=1,NN
          B(1,3*J-2) = gderiv(J,1)
          B(2,3*J-1) = gderiv(J,2)
          B(3,3*J  ) = gderiv(J,3)
          B(4,3*J-2) = gderiv(J,2)
          B(4,3*J-1) = gderiv(J,1)
          B(5,3*J-1) = gderiv(J,3)
          B(5,3*J  ) = gderiv(J,2)
          B(6,3*J-2) = gderiv(J,3)
          B(6,3*J  ) = gderiv(J,1)
        ENDDO
      endif

!!  calculate the Internal Force 
      WG=getWeight( etype, LX )*DET
      qf(1:nn*ndof) = qf(1:nn*ndof) +                                   &
          matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG
    enddo

  end subroutine UPDATE_C3
!
!
!
!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C3(ETYPE,NN,XX,YY,ZZ,TT,T0,EDISP,gausses,coords)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
      use m_fstr
      use mMechGauss
      use m_MatMatrix
      use m_utilities
!
      INTEGER(kind=kint), PARAMETER :: NDOF=3
! I/F VARIABLES      
      INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
      REAL(kind=kreal), INTENT(IN)   :: XX(NN),YY(NN),ZZ(NN)
      REAL(kind=kreal), INTENT(IN)   :: TT(NN),T0(NN),EDISP(NN*NDOF)
      TYPE(tGaussStatus), INTENT(INOUT)  :: gausses(:)
      REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)
! LOCAL VARIABLES
      REAL(kind=kreal) ALP,ALP0,D(6,6),B(6,NDOF*NN)
      REAL(kind=kreal) DET,ecoord(3,NN)
      INTEGER(kind=kint) J,K,IC,cdsys_ID,serr
      REAL(kind=kreal) EPSA(6),EPSTH(6),SGM(6),H(NN),alpo(3),alpo0(3),coordsys(3,3)
      REAL(kind=kreal) naturalcoord(3),gderiv(NN,3)
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS,tm(6,6),outa(1),ina(1)
      logical   :: ierr, matlaniso

        matlaniso = .false.
        cdsys_ID = gausses(1)%pMaterial%cdsys_ID
        if( cdsys_ID>0 ) then   ! cannot define aniso exapansion when no local coord defined
          ina = TT(1)
          call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
          if( .not. ierr ) matlaniso = .true.
        endif

      ecoord(1,:)=XX(:)
      ecoord(2,:)=YY(:)
      ecoord(3,:)=ZZ(:)
! LOOP FOR INTEGRATION POINTS
      DO IC=1,NumOfQuadPoints(etype)
        CALL getQuadPoint( etype, IC, naturalCoord(:) )
        CALL getShapeFunc( etype, naturalcoord, H(1:NN) )
        CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv )

        if( matlaniso ) then
          call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
          if( serr==-1 ) stop "Fail to setup local coordinate"
          if( serr==-2 ) then
            write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
          endif
        endif

            B(1:6,1:NN*NDOF)=0.d0
            DO J=1,NN
              B(1,3*J-2)=gderiv(J,1)
              B(2,3*J-1)=gderiv(J,2)
              B(3,3*J  )=gderiv(J,3)
              B(4,3*J-2)=gderiv(J,2)
              B(4,3*J-1)=gderiv(J,1)
              B(5,3*J-1)=gderiv(J,3)
              B(5,3*J  )=gderiv(J,2)
              B(6,3*J-2)=gderiv(J,3)
              B(6,3*J  )=gderiv(J,1)
            ENDDO

            TEMPC=DOT_PRODUCT( H(1:NN),TT(1:NN) )
            TEMP0=DOT_PRODUCT( H(1:NN),T0(1:NN) )

        CALL MatlMatrix( gausses(IC), D3, D, 0.d0, coordsys, TEMPC )

            ina(1) = TEMPC
            if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
               call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
            endif
            ina(1) = TEMP0
            if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
              call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
              alp0 = outa(1)
            endif

!**
!** THERMAL strain
!**
            if( matlaniso ) then
              do j=1,3
                EPSTH(j)=ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
              enddo
              EPSTH(4:6)=0.d0
              call transformation(coordsys, tm) 
              EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord  
              EPSTH(4:6)=EPSTH(4:6)*2.d0 
            else
              THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
              EPSTH(1:3)=THERMAL_EPS
              EPSTH(4:6)=0.d0
            endif
!**
!** SET EPS  {e}=[B]{u}
!**
            EPSA(1:6) = MATMUL( B(1:6,:), EDISP )
!**
!** SET SGM  {S}=[D]{e}
!**
            DO J=1,6
              SGM(J)=0.d0
              DO K=1,6
                SGM(J)=SGM(J)+D(J,K)*(EPSA(K)-EPSTH(K))
              ENDDO
            ENDDO
!**
!** Adding stress in each gauss points
!**
            gausses(IC)%strain(1:6)=EPSA(1:6)
            gausses(IC)%stress(1:6)=SGM(1:6)

     ENDDO

   end subroutine UpdateST_C3
   

!----------------------------------------------------------------------*
   SUBROUTINE NodalStress_C3(ETYPE,NN,gausses,ndstrain,ndstress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: ndstrain(NN,6)
   REAL(kind=kreal), INTENT(OUT)  :: ndstress(NN,6)
!
   INTEGER          :: i,ic
   REAL(kind=kreal) :: TEMP(12)
!
   TEMP(:)=0.d0
   IC = NumOfQuadPoints(etype)
   DO i=1,IC
     TEMP(1:6) = TEMP(1:6) + gausses(i)%strain(1:6)
     TEMP(7:12) = TEMP(7:12) + gausses(i)%stress(1:6)
   ENDDO
   TEMP(1:12) = TEMP(1:12)/IC
   FORALL( i=1:NN )
     ndstrain(i,1:6) = TEMP(1:6)
     ndstress(i,1:6) = TEMP(7:12)
   END FORALL
!
   END SUBROUTINE
!
!
!
!----------------------------------------------------------------------*
   SUBROUTINE ElementStress_C3(ETYPE,gausses,strain,stress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: strain(6)
   REAL(kind=kreal), INTENT(OUT)  :: stress(6)

   INTEGER          :: i,ic

   strain(:)=0.d0; stress(:)=0.d0
   IC = NumOfQuadPoints(etype)
   DO i=1,IC
     strain(:) = strain(:) + gausses(i)%strain(1:6)
     stress(:) = stress(:) + gausses(i)%stress(1:6)
   ENDDO
   strain(:) = strain(:)/IC
   stress(:) = stress(:)/IC

   END SUBROUTINE
   
!> Volume of element
!----------------------------------------------------------------------*
   real(kind=kreal) function VOLUME_C3(ETYPE,NN,XX,YY,ZZ)
!----------------------------------------------------------------------*
      INTEGER(kind=kint), INTENT(IN)  :: ETYPE, NN
      REAL(kind=kreal), INTENT(IN)    :: XX(:),YY(:),ZZ(:)
      
      REAL(kind=kreal) XJ(3,3),DET,WG
      INTEGER(kind=kint) LX,I 
      REAL(kind=kreal) localcoord(3),deriv(NN,3)

      VOLUME_C3 = 0.d0
! LOOP FOR INTEGRATION POINTS
        DO  LX=1,NumOfQuadPoints( ETYPE )
              CALL getQuadPoint( ETYPE, LX, localcoord )
              CALL getShapeDeriv( ETYPE, localcoord, deriv )
!  JACOBI MATRIX
              XJ(1,1:3)= matmul( xx(1:NN), deriv(1:NN,1:3) )
              XJ(2,1:3)= matmul( yy(1:NN), deriv(1:NN,1:3) )
              XJ(3,1:3)= matmul( zz(1:NN), deriv(1:NN,1:3) )
!DETERMINANT OF JACOBIAN
              DET=XJ(1,1)*XJ(2,2)*XJ(3,3)                                                 &
                 +XJ(2,1)*XJ(3,2)*XJ(1,3)                                                 &
                 +XJ(3,1)*XJ(1,2)*XJ(2,3)                                                 &
                 -XJ(3,1)*XJ(2,2)*XJ(1,3)                                                 &
                 -XJ(2,1)*XJ(1,2)*XJ(3,3)                                                 &
                 -XJ(1,1)*XJ(3,2)*XJ(2,3)

              VOLUME_C3 = VOLUME_C3+getWeight( etype, LX )*DET
        ENDDO
       
   end function VOLUME_C3
   
end module m_static_LIB_3d
