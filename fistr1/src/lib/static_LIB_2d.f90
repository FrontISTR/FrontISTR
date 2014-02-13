!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN (AdavanceSoft)                         !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!>   This module provide common functions of Plane deformation elements
!>
!!
!>  \author     X. YUAN (AdavanceSoft)
!>  \date       2009/06/12
!>  \version    0.00
!!
!======================================================================!
module m_static_LIB_2d
   use hecmw, only : kint, kreal
   use elementInfo
   implicit none

   contains
!***********************************************************************
!  2D Element:
!  STF_C2   :Calculate stiff matrix of 2d elements
!  DL_C2    :Deal with DLOAD conditions
!  TLOAD_C2 :Deal with thermal expansion force
!  UpdateST_C2  : Update Strain and stress
!***********************************************************************
!----------------------------------------------------------------------*
   SUBROUTINE STF_C2( ETYPE,NN,ecoord,gausses,PARAM1,stiff         &
                      ,ISET, u)
!----------------------------------------------------------------------*
      USE mMechGauss
      use m_MatMatrix
      INTEGER(kind=kint), INTENT(IN) :: ETYPE
      INTEGER(kind=kint), INTENT(IN) :: NN
      TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
      REAL(kind=kreal), INTENT(IN)   :: ecoord(2,NN),PARAM1
      REAL(kind=kreal), INTENT(OUT)  :: stiff(:,:)
      INTEGER(kind=kint),INTENT(IN)  :: ISET
      REAL(kind=kreal), INTENT(IN), optional :: u(:,:)

! LOCAL VARIABLES
      integer(kind=kint) :: flag
      INTEGER(kind=kint) NDOF
      PARAMETER(NDOF=2)
      REAL(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
      REAL(kind=kreal) H(NN),stress(4)
      REAL(kind=kreal) THICK,PAI
      REAL(kind=kreal) DET,RR,WG
      REAL(kind=kreal) localcoord(2)
      REAL(kind=kreal) gderiv(NN,2), cdsys(3,3)
      INTEGER(kind=kint) J,LX
      real(kind=kreal) gdispderiv(2,2)
      real(kind=kreal) B1(4,nn*ndof)
      real(kind=kreal) Smat(4,4)
      real(kind=kreal) BN(4,nn*ndof), SBN(4,nn*ndof)

      stiff =0.d0
      flag = gausses(1)%pMaterial%nlgeom_flag
! THICKNESS
      THICK=PARAM1
!  FOR AX-SYM. ANALYSIS
      IF(ISET==2) THEN
        THICK=1.d0
        PAI=4.d0*ATAN(1.d0)
      ENDIF
!* LOOP OVER ALL INTEGRATION POINTS
      DO LX=1,NumOfQuadPoints( ETYPE )
        CALL MatlMatrix( gausses(LX), ISET, D,1.d0,cdsys )
        if( .not. present(u) ) flag=INFINITE    ! enforce to infinite deformation analysis

        if( flag==1 .and. ISET == 2 ) then
          write(*,'(a)') '    PROGRAM STOP : non-TL element for axixsymmetric element'
          stop
        endif

        CALL getQuadPoint( ETYPE, LX, localcoord )
        CALL getGlobalDeriv( etype, nn, localcoord, ecoord, det, gderiv )
!
        IF(ISET==2) THEN
          CALL getShapeFunc( ETYPE, localcoord, H(:) )
          RR=DOT_PRODUCT( H(1:NN), ecoord(1,1:NN) )
          WG=getWeight( ETYPE, LX )*DET*RR*2.d0*PAI
        ELSE
          RR=THICK
          H(:)=0.d0
          WG=getWeight( ETYPE, LX )*DET*RR
        END IF
        DO J=1,NN
          B(1,2*J-1)=gderiv(J,1)
          B(2,2*J-1)=0.d0
          B(3,2*J-1)=gderiv(J,2)
          B(1,2*J  )=0.d0
          B(2,2*J  )=gderiv(J,2)
          B(3,2*J  )=gderiv(J,1)
          B(4,2*J-1)=H(J)/RR
          B(4,2*J  )=0.d0
        ENDDO
!
! ----- calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
        if( flag==1 ) then
!       ----- dudx(i,j) ==> gdispderiv(i,j)
          gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
          B1(1:4,1:nn*ndof) = 0.d0
          do j=1,nn
            B1(1,2*j-1) = gdispderiv(1,1)*gderiv(j,1)
            B1(2,2*j-1) = gdispderiv(1,2)*gderiv(j,2)
            B1(3,2*j-1) = gdispderiv(1,2)*gderiv(j,1)+gdispderiv(1,1)*gderiv(j,2)
            B1(1,2*j  ) = gdispderiv(2,1)*gderiv(j,1)
            B1(2,2*j  ) = gdispderiv(2,2)*gderiv(j,2)
            B1(3,2*j  ) = gdispderiv(2,2)*gderiv(j,1)+gdispderiv(2,1)*gderiv(j,2)
            B1(4,2*j-1) = 0.d0
            B1(4,2*j  ) = 0.d0
          enddo
          do j=1,nn*ndof
            B(:,j) = B(:,j)+B1(:,j)
          enddo
        endif
!
        DB(1:4,1:nn*ndof) = 0.d0
        DB(1:4,1:nn*ndof) = MATMUL( D(1:4,1:4), B(1:4,1:nn*ndof) )
        stiff(1:nn*ndof,1:nn*ndof) = stiff(1:nn*ndof,1:nn*ndof) +             &
           matmul( transpose( B(1:4,1:nn*ndof) ), DB(1:4,1:nn*ndof) )*WG
!
!    ----- calculate the stress matrix ( TOTAL LAGRANGE METHOD )
        if( flag==1 ) then
          stress(1:4)=gausses(LX)%stress(1:4)
          BN(1:4,1:nn*ndof) = 0.d0
          do j=1,nn
            BN(1,2*j-1) = gderiv(j,1)
            BN(2,2*j  ) = gderiv(j,1)
            BN(3,2*j-1) = gderiv(j,2)
            BN(4,2*j  ) = gderiv(j,2)
          enddo
          Smat(:,:) = 0.d0
          do j=1,2
            Smat(j  ,j  ) = stress(1)
            Smat(j  ,j+2) = stress(3)
            Smat(j+2,j  ) = stress(3)
            Smat(j+2,j+2) = stress(2)
          enddo
          SBN(1:4,1:nn*ndof) = matmul( Smat(1:4,1:4), BN(1:4,1:nn*ndof) )
          stiff(1:nn*ndof,1:nn*ndof) = stiff(1:nn*ndof,1:nn*ndof) +           &
             matmul( transpose( BN(1:4,1:nn*ndof) ), SBN(1:4,1:nn*ndof) )*WG
        endif
!
      ENDDO

   end subroutine STF_C2


!----------------------------------------------------------------------*
   SUBROUTINE DL_C2(ETYPE,NN,XX,YY,RHO,PARAM1,LTYPE,PARAMS,VECT,NSIZE,ISET)
!----------------------------------------------------------------------*
!**
!**  Deal with DLOAD conditions
!** 
!   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
!   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
!   GRAV LTYPE=4  :GRAVITY FORCE
!   CENT LTYPE=5  :CENTRIFUGAL LOAD
!   P1   LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR FACE-1
!   P2   LTYPE=20 :TRACTION IN NORMAL-DIRECTION FOR FACE-2
!   P3   LTYPE=30 :TRACTION IN NORMAL-DIRECTION FOR FACE-3
!   P4   LTYPE=40 :TRACTION IN NORMAL-DIRECTION FOR FACE-4
! I/F VARIABLES
      INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
      REAL(kind=kreal), INTENT(IN)   :: XX(NN),YY(NN),PARAMS(0:6)
      REAL(kind=kreal), INTENT(OUT)  :: VECT(NN*2)
      REAL(kind=kreal) RHO,PARAM1
      INTEGER(kind=kint) LTYPE,NSIZE,ISET
! LOCAL VARIABLES
      INTEGER(kind=kint), PARAMETER :: NDOF=2
      REAL(kind=kreal) H(NN)
      REAL(kind=kreal) PLX(NN),PLY(NN)
      REAL(kind=kreal) XJ(2,2),DET,RR,WG
      REAL(kind=kreal) PAI,VAL,THICK
      INTEGER(kind=kint) IVOL,ISUF
      REAL(kind=kreal) AX,AY,RX,RY
      REAL(kind=kreal) normal(2)
      REAL(kind=kreal) COEFX,COEFY,XCOD,YCOD,HX,HY,PHX,PHY
      INTEGER(kind=kint) NOD(NN)
      INTEGER(kind=kint) LX,I,SURTYPE,NSUR
      REAL(kind=kreal) VX,VY,localcoord(2),deriv(NN,2),elecoord(2,NN)
!  CONSTANT
      PAI=4.d0*ATAN(1.d0)
! SET VALUE
      VAL=PARAMS(0)
! THICKNESS
      THICK=PARAM1
! CLEAR VECT
      NSIZE=NN*NDOF
      VECT(1:NSIZE)=0.d0
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
          RX=PARAMS(4)
          RY=PARAMS(5)
        ENDIF
      ELSE IF( LTYPE.GE.10 ) THEN 
        ISUF=1
        CALL getSubFace( ETYPE, LTYPE/10, SURTYPE, NOD )
        NSUR = getNumberOfNodes( SURTYPE )
      ENDIF
!** SURFACE LOAD
      IF( ISUF==1 ) THEN
        DO I=1,NSUR
          elecoord(1,i)=XX(NOD(I))
          elecoord(2,i)=YY(NOD(i))
        ENDDO
!** INTEGRATION OVER SURFACE
        DO LX=1,NumOfQuadPoints( SURTYPE )
          CALL getQuadPoint( SURTYPE, LX, localcoord(1:1) )
          CALL getShapeFunc( SURTYPE, localcoord(1:1), H(1:NSUR) )
          normal=EdgeNormal( SURTYPE, NSUR, localcoord(1:1), elecoord(:,1:NSUR) )
          WG = getWeight( SURTYPE, LX )
          IF( ISET==2 ) THEN
            RR=0.d0
            DO I=1,NSUR
              RR=RR+H(I)*XX(NOD(I))
            ENDDO
            WG=WG*RR*2.d0*PAI
          ELSE
            WG=WG*THICK
          ENDIF
          DO I=1,NSUR
            VECT(2*NOD(I)-1)=VECT(2*NOD(I)-1)+VAL*WG*H(I)*normal(1)
            VECT(2*NOD(I)  )=VECT(2*NOD(I)  )+VAL*WG*H(I)*normal(2)
          ENDDO
        ENDDO
      ENDIF
!** VOLUME LOAD
      IF( IVOL==1 ) THEN
        PLX(:)=0.d0
        PLY(:)=0.d0
        DO LX=1,NumOfQuadPoints( ETYPE )
          CALL getQuadPoint( ETYPE, LX, localcoord(1:2) )
          CALL getShapeDeriv( ETYPE, localcoord(1:2), deriv )
          CALL getShapeFunc( ETYPE, localcoord(1:2), H(1:NN) )
          XJ(1,1:2)=MATMUL( XX(1:NN), deriv(1:NN,1:2) )
          XJ(2,1:2)=MATMUL( YY(1:NN), deriv(1:NN,1:2) )

          DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)

            WG = getWeight( ETYPE, LX )
            IF(ISET==2) THEN
              RR=DOT_PRODUCT( H(1:NN),XX(1:NN) )
              WG=WG*DET*RR*2.d0*PAI
            ELSE
              RR=THICK
              WG=WG*DET*RR
            END IF
            COEFX=1.d0
            COEFY=1.d0
! CENTRIFUGAL LOAD
            IF( LTYPE==5 ) THEN
              XCOD=DOT_PRODUCT( H(1:NN),XX(1:NN) )
              YCOD=DOT_PRODUCT( H(1:NN),YY(1:NN) )
              HX=AX + ( (XCOD-AX)*RX+(YCOD-AY)*RY )/(RX**2+RY**2)*RX
              HY=AY + ( (XCOD-AX)*RX+(YCOD-AY)*RY )/(RX**2+RY**2)*RY
              PHX=XCOD-HX
              PHY=YCOD-HY
              COEFX=RHO*VAL*VAL*PHX
              COEFY=RHO*VAL*VAL*PHY
            END IF
            DO I=1,NN
              PLX(I)=PLX(I)+H(I)*WG*COEFX
              PLY(I)=PLY(I)+H(I)*WG*COEFY
            ENDDO
        ENDDO
        IF( LTYPE.EQ.1) THEN
          DO I=1,NN
            VECT(2*I-1)=VAL*PLX(I)
            VECT(2*I  )=0.d0
          ENDDO
        ELSE IF( LTYPE.EQ.2 ) THEN
          DO I=1,NN
            VECT(2*I-1)=0.d0
            VECT(2*I  )=VAL*PLY(I)
          ENDDO
        ELSE IF( LTYPE.EQ.4 ) THEN
          VX=PARAMS(1)
          VY=PARAMS(2)
          VX=VX/SQRT( PARAMS(1)**2 + PARAMS(2)**2 )
          VY=VY/SQRT( PARAMS(1)**2 + PARAMS(2)**2 )
          DO I=1,NN
            VECT(2*I-1)=VAL*PLX(I)*RHO*VX
            VECT(2*I  )=VAL*PLY(I)*RHO*VY
          ENDDO
        ELSE IF( LTYPE==5 ) THEN
          DO I=1,NN
            VECT(2*I-1)=PLX(I)
            VECT(2*I  )=PLY(I)
          ENDDO  
        END IF
      ENDIF
    end subroutine DL_C2
   
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C2(ETYPE,NN,XX,YY,TT,T0,gausses,PARAM1,ISET,VECT)
!----------------------------------------------------------------------*
!
! THERMAL LOAD OF PLANE ELEMENT
!
      use mMaterial
      use m_MatMatrix
      use m_ElasticLinear
! I/F VARIABLES
      INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
      REAL(kind=kreal),INTENT(IN)    :: XX(NN),YY(NN),TT(NN),T0(NN),PARAM1
      TYPE(tGaussStatus), INTENT(IN) :: gausses(:)          !< status of qudrature points
      REAL(kind=kreal),INTENT(OUT)   :: VECT(NN*2)
      INTEGER(kind=kint) ISET
! LOCAL VARIABLES
      INTEGER(kind=kint) NDOF
      PARAMETER(NDOF=2)
      REAL(kind=kreal) ALP,PP,D(4,4),B(4,NDOF*NN)
      REAL(kind=kreal) H(NN)
      REAL(kind=kreal) EPS(4),SGM(4),localcoord(2)
      REAL(kind=kreal) deriv(NN,2),DNDE(2,NN)
      REAL(kind=kreal) XJ(2,2),DET,RR,DUM
      REAL(kind=kreal) XJI(2,2)
      REAL(kind=kreal) PAI,THICK,WG
      INTEGER(kind=kint) J,LX
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS
!*************************
!  CONSTANT
!*************************
      PAI=4.d0*ATAN(1.d0)
! CLEAR VECT
      VECT(1:NN*NDOF)=0.d0
! THICKNESS
      THICK=PARAM1
!  FOR AX-SYM. ANALYSIS
      IF(ISET==2) THICK=1.d0
! We suppose material properties doesn't varies inside element
      CALL calElasticMatrix( gausses(1)%pMaterial, ISET, D  )
	  ALP = gausses(1)%pMaterial%variables(M_EXAPNSION)
      pp = gausses(1)%pMaterial%variables(M_POISSON)
!* LOOP OVER ALL INTEGRATION POINTS
      DO LX=1,NumOfQuadPoints( ETYPE )    
          CALL getQuadPoint( ETYPE, LX, localcoord )
          CALL getShapeFunc( ETYPE, localcoord, H(1:NN) )
          CALL getShapeDeriv( ETYPE, localcoord, deriv(:,:) )
          XJ(1,1:2)=MATMUL( XX(1:NN), deriv(1:NN,1:2) )
          XJ(2,1:2)=MATMUL( YY(1:NN), deriv(1:NN,1:2) )

          DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)

          TEMPC=DOT_PRODUCT(H(1:NN),TT(1:NN))
          TEMP0=DOT_PRODUCT(H(1:NN),T0(1:NN))
          IF(ISET==2) THEN
            RR=DOT_PRODUCT(H(1:NN),XX(1:NN))
            WG=getWeight( ETYPE, LX )*DET*RR*2.d0*PAI
          ELSE
            RR=THICK
            WG=getWeight( ETYPE, LX )*DET*RR
          END IF
          DUM=1.d0/DET
          XJI(1,1)= XJ(2,2)*DUM
          XJI(1,2)=-XJ(2,1)*DUM
          XJI(2,1)=-XJ(1,2)*DUM
          XJI(2,2)= XJ(1,1)*DUM

          DNDE=MATMUL( XJI, TRANSPOSE(deriv) )
          DO J=1,NN
            B(1,2*J-1)=DNDE(1,J)
            B(2,2*J-1)=0.d0
            B(3,2*J-1)=DNDE(2,J)
            B(1,2*J  )=0.d0
            B(2,2*J  )=DNDE(2,J)
            B(3,2*J  )=DNDE(1,J)
            B(4,2*J-1)=0.d0
            B(4,2*J  )=0.d0
            IF( ISET==2 ) THEN
              B(4,2*J-1)=H(J)/RR
            ENDIF
          ENDDO
!**
!** THERMAL EPS
!**
          THERMAL_EPS=ALP*(TEMPC-TEMP0)
          IF( ISET .EQ. 2 ) THEN
            EPS(1)=THERMAL_EPS
            EPS(2)=THERMAL_EPS
            EPS(3)=0.d0
            EPS(4)=THERMAL_EPS
          ELSE IF( ISET.EQ.0 ) THEN
            EPS(1)=THERMAL_EPS*(1.d0+PP)
            EPS(2)=THERMAL_EPS*(1.d0+PP)
            EPS(3)=0.d0
            EPS(4)=0.d0
          ELSE IF( ISET.EQ.1 ) THEN
            EPS(1)=THERMAL_EPS
            EPS(2)=THERMAL_EPS
            EPS(3)=0.d0
            EPS(4)=0.d0
          ENDIF
!**
!** SET SGM  {s}=[D]{e}
!**
          SGM = MATMUL( EPS(1:4), D(1:4,1:4) )
!**
!** CALCULATE LOAD {F}=[B]T{e}
!**
          VECT(1:NN*NDOF)=VECT(1:NN*NDOF)+MATMUL( SGM(1:4), B(1:4,1:NN*NDOF) )*WG
      ENDDO

   end subroutine TLOAD_C2
!
!
!---------------------------------------------------------------------*
  SUBROUTINE UPDATE_C2( etype,nn,ecoord,gausses,PARAM1,iset,         &
                        u, ddu, qf  )
!---------------------------------------------------------------------*
    use mMechGauss
    use m_MatMatrix
! I/F VARIAVLES
    integer(kind=kint), INTENT(IN)     :: etype, nn
    real(kind=kreal),   INTENT(IN)     :: ecoord(3,nn),PARAM1
    INTEGER(kind=kint), INTENT(IN)     :: ISET
    TYPE(tGaussStatus), INTENT(INOUT)  :: gausses(:)
    real(kind=kreal),   INTENT(IN)     :: u(2,nn)
    real(kind=kreal),   INTENT(IN)     :: ddu(2,nn)
    real(kind=kreal),   INTENT(OUT)    :: qf(:)
    

    integer(kind=kint), parameter :: ndof=2
    real(kind=kreal)   :: D(4,4), B(4,ndof*nn), B1(4,ndof*nn)
    real(kind=kreal)   :: H(nn)
    real(kind=kreal)   :: thick,pai,rr
    real(kind=kreal)   :: det, WG
    real(kind=kreal)   :: localCoord(2)
    real(kind=kreal)   :: gderiv(nn,2), dstrain(4)
    real(kind=kreal)   :: gdispderiv(2,2), cdsys(3,3)
    integer(kind=kint) :: j, LX
    real(kind=kreal)   :: totaldisp(2*nn)
!
    qf(:)              = 0.d0
    do j=1,nn
      totaldisp((j-1)*2+1) = u(1,j) + ddu(1,j)
      totaldisp(j*2)       = u(2,j) + ddu(2,j)
    enddo
!
! THICKNESS
    THICK=PARAM1
!  FOR AX-SYM. ANALYSIS
    IF(ISET==2) THEN
      THICK=1.d0
      PAI=4.d0*ATAN(1.d0)
    ENDIF
!* LOOP OVER ALL INTEGRATION POINTS
    do LX=1, NumOfQuadPoints(etype)
      call MatlMatrix( gausses(LX), ISET, D, 1.d0, cdsys  )

      call getQuadPoint( etype, LX, localCoord(:) )
      call getGlobalDeriv( etype, nn, localcoord, ecoord, det, gderiv )
!
      WG=getWeight( etype, LX )*DET
      IF(ISET==2) THEN
        CALL getShapeFunc( ETYPE, localcoord, H(:) )
        RR=DOT_PRODUCT( H(1:NN), ecoord(1,1:NN) )
      ELSE
        RR=THICK
        H(:)=0.d0
      END IF
      DO J=1,NN
        B(1,2*J-1)=gderiv(J,1)
        B(2,2*J-1)=0.d0
        B(3,2*J-1)=gderiv(J,2)
        B(1,2*J  )=0.d0
        B(2,2*J  )=gderiv(J,2)
        B(3,2*J  )=gderiv(J,1)
        B(4,2*J-1)=H(J)/RR
        B(4,2*J  )=0.d0
      ENDDO
	  
      gausses(LX)%strain(1:4) = matmul( B(:,:), totaldisp )
	  gausses(LX)%stress(1:4) = matmul( D(1:4, 1:4), gausses(LX)%strain(1:4) )
	  
!
!    ----- calculate the Internal Force
      qf(1:nn*ndof) = qf(1:nn*ndof) +                                     &
          matmul( gausses(LX)%stress(1:4), B(1:4,1:nn*ndof) )*WG
!
    enddo
!
end subroutine UPDATE_C2
!

!
!
!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C2(ETYPE,NN,XX,YY,TT,T0,PARAM1,ISET,EDISP,gauss)
!----------------------------------------------------------------------*
!
! Update strain and stress of PLANE ELEMENT
!
      use mMechGauss
      use m_MatMatrix
! I/F VARIABLES
      INTEGER(kind=kint), INTENT(IN) :: ETYPE, NN
      REAL(kind=kreal), INTENT(IN)   :: XX(NN),YY(NN),TT(NN),T0(NN),EDISP(NN*2)
      TYPE(tGaussStatus), INTENT(INOUT)  :: gauss(:)
      REAL(kind=kreal) ALP,PARAM1
      INTEGER(kind=kint) ISET
! LOCAL VARIABLES
      INTEGER(kind=kint) NDOF
      PARAMETER(NDOF=2)
      REAL(kind=kreal) D(4,4),B(4,NDOF*NN)
      REAL(kind=kreal) H(NN)
      REAL(kind=kreal) EPS(4),SGM(4),EPSTH(4),deriv(NN,2)
      REAL(kind=kreal) XJ(2,2),DNDE(2,NN),DET,RR,DUM
      REAL(kind=kreal) XJI(2,2),localcoord(2)
      REAL(kind=kreal) PAI,THICK,EE,PP, cdsys(3,3)
      INTEGER(kind=kint) J,K
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS
      REAL(kind=kreal) ARRAY_TEMP(8)
      INTEGER(kind=kint) IC
!*************************
!  CONSTANT
!*************************
      PAI=4.d0*ATAN(1.d0)
      ARRAY_TEMP=0.d0
! THICKNESS
      THICK=PARAM1
!  FOR AX-SYM. ANALYSIS
      IF(ISET==2) THICK=1.d0
!* LOOP OVER INTEGRATION POINT
      DO IC=1,NumOfQuadPoints( ETYPE )
          EE = gauss(IC)%pMaterial%variables(M_YOUNGS)
          PP = gauss(IC)%pMaterial%variables(M_POISSON)
          ALP = gauss(IC)%pMaterial%variables(M_EXAPNSION)
          CALL MatlMatrix( gauss(IC), ISET, D, 1.d0, cdsys  )
          CALL getQuadPoint( ETYPE, IC, localcoord )
          CALL getShapeDeriv( ETYPE, localcoord, deriv(:,:) )
          CALL getShapeFunc( ETYPE, localcoord, H(1:NN) )

          XJ(1,1:2)=MATMUL( XX(1:NN), deriv(1:NN,1:2) )
          XJ(2,1:2)=MATMUL( YY(1:NN), deriv(1:NN,1:2) )
          DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)
          DUM=1.d0/DET
          XJI(1,1)= XJ(2,2)*DUM
          XJI(1,2)=-XJ(2,1)*DUM
          XJI(2,1)=-XJ(1,2)*DUM
          XJI(2,2)= XJ(1,1)*DUM

          DNDE=MATMUL( XJI, TRANSPOSE(deriv) )

          TEMPC=DOT_PRODUCT( H(1:NN),TT(1:NN) )
          TEMP0=DOT_PRODUCT( H(1:NN),T0(1:NN) )
          IF(ISET==2) THEN
            RR=DOT_PRODUCT( H(1:NN),XX(1:NN) )
          ELSE
            RR=THICK
          END IF
          DO J=1,NN
            B(1,2*J-1)=DNDE(1,J)
            B(2,2*J-1)=0.d0
            B(3,2*J-1)=DNDE(2,J)
            B(1,2*J  )=0.d0
            B(2,2*J  )=DNDE(2,J)
            B(3,2*J  )=DNDE(1,J)
            B(4,2*J-1)=0.d0
            B(4,2*J  )=0.d0
            IF( ISET==2 ) THEN
              B(4,2*J-1)=H(J)/RR
            ENDIF
          ENDDO
!**
!** THERMAL EPS
!**
          THERMAL_EPS=ALP*(TEMPC-TEMP0)
          IF( ISET .EQ. 2 ) THEN
            EPSTH(1)=THERMAL_EPS
            EPSTH(2)=THERMAL_EPS
            EPSTH(3)=0.d0
            EPSTH(4)=THERMAL_EPS
          ELSE IF( ISET.EQ.0 ) THEN
            EPSTH(1)=THERMAL_EPS*(1.0+PP)
            EPSTH(2)=THERMAL_EPS*(1.0+PP)
            EPSTH(3)=0.d0
            EPSTH(4)=0.d0
          ELSE IF( ISET.EQ.1 ) THEN
            EPSTH(1)=THERMAL_EPS
            EPSTH(2)=THERMAL_EPS
            EPSTH(3)=0.d0
            EPSTH(4)=0.d0
          ENDIF
!**
!** SET EPS  {e}=[B]{u}
!**
          EPS(1:4)=MATMUL( B(1:4,1:NN*NDOF), EDISP(1:NN*NDOF) )
!**
!** SET SGM  {S}=[D]{e}
!**
          DO J=1,4
            SGM(J)=0.d0
            DO K=1,4
              SGM(J)=SGM(J)+D(J,K)*(EPS(K)-EPSTH(K))
            ENDDO
          ENDDO
!**
!** FOR PLANE STRAIN
!**
          IF( ISET==0 ) THEN
            SGM(4)=PP*(SGM(1)+SGM(2))-EE*THERMAL_EPS
          ENDIF
!**
!** Elemental average
!**
          gauss(IC)%strain(1:4)=EPS(1:4)
          gauss(IC)%stress(1:4)=SGM(1:4)
       ENDDO

   end subroutine UpdateST_C2
   
!----------------------------------------------------------------------*
   SUBROUTINE NodalStress_C2(ETYPE,NN,gausses,ndstrain,ndstress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss

   INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: ndstrain(NN,4)
   REAL(kind=kreal), INTENT(OUT)  :: ndstress(NN,4)

   INTEGER          :: i,ic
   REAL(kind=kreal) :: TEMP(8)
   
   TEMP(:)=0.d0
   IC = NumOfQuadPoints(etype)
   DO i=1,IC
     TEMP(1:4) = TEMP(1:4) + gausses(i)%strain(1:4)
     TEMP(5:8) = TEMP(5:8) + gausses(i)%stress(1:4)
   ENDDO
   TEMP(1:8) = TEMP(1:8)/IC
   FORALL( i=1:NN )
     ndstrain(i,1:4) = TEMP(1:4)
     ndstress(i,1:4) = TEMP(5:8)
   END FORALL

   END SUBROUTINE
   
!----------------------------------------------------------------------*
   SUBROUTINE ElementStress_C2(ETYPE,gausses,strain,stress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: strain(4)
   REAL(kind=kreal), INTENT(OUT)  :: stress(4)

   INTEGER          :: i,ic

   strain(:)=0.d0; stress(:)=0.d0
   IC = NumOfQuadPoints(etype)
   DO i=1,IC
     strain(:) = strain(:) + gausses(i)%strain(1:4)
     stress(:) = stress(:) + gausses(i)%stress(1:4)
   ENDDO
   strain(:) = strain(:)/IC
   stress(:) = stress(:)/IC

   END SUBROUTINE

end module m_static_LIB_2d
