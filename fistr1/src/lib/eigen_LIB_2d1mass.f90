!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!C================================================================C
!>  This module contains subroutines used in 2d eigen analysis for
!!  element 231, 241
!C================================================================C

module m_eigen_LIB_2d1mass
contains

!C*--------------------------------------------------------------------*
      SUBROUTINE MASS_C2D4(XX,YY,EE,PP,RHO,PARAM1,SS,ISET,myEIG)
!C*--------------------------------------------------------------------*
!C*
!C* CALCULATION 2D 4 NODE PLANE ELEMENT
!C*
      use hecmw
      use gauss_integration
      use lczparm
 
      IMPLICIT NONE
!C* I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),SS(*),EE,PP,PARAM1
      INTEGER(kind=kint) ISET
!C* LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=4,NDOF=2,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=3)
      REAL(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN)
      REAL(kind=kreal) RHO,THICK,COEF1,COEF2,PAI
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,RR,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI12,XJI22
      INTEGER(kind=kint) I,J,J2,K,LX,LY,NUM
      INTEGER(kind=kint) ind1, ind2 
      REAL(kind=kreal) totdiag, totmass
      TYPE(lczparam) :: myEIG

      totdiag = 0.0
      totmass = 0.0
!C**************************
!C*  CONSTANT
!C**************************
      PAI=4.0*ATAN(1.0)
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
!C*THICKNESS 
      THICK=PARAM1
!C*FOR AX-SYM. ANALYSIS
      IF(ISET.EQ.2) THEN
        THICK=1.0
      END IF
!C*CLEAR [D] MATRIX
      DO J=1,4
        DO I=1,4
          D(I,J)=0.0
        ENDDO
      ENDDO
!C*LOOP OVER ALL INTEGRATION POINTS
      DO LX=1,NG
        RI=XG(NG,LX)
        DO LY=1,NG
          SI=XG(NG,LY)
          RP=1.0+RI
          SP=1.0+SI
          RM=1.0-RI
          SM=1.0-SI
          H(1)=0.25*RM*SM
          H(2)=0.25*RP*SM
          H(3)=0.25*RP*SP
          H(4)=0.25*RM*SP
          HR(1)=-.25*SM
          HR(2)= .25*SM
          HR(3)= .25*SP
          HR(4)=-.25*SP
          HS(1)=-.25*RM
          HS(2)=-.25*RP
          HS(3)= .25*RP
          HS(4)= .25*RM
          XJ11=0.0
          XJ21=0.0
          XJ12=0.0
          XJ22=0.0
          DO I=1,NN
            XJ11=XJ11+HR(I)*XX(I)
            XJ21=XJ21+HS(I)*XX(I)
            XJ12=XJ12+HR(I)*YY(I)
            XJ22=XJ22+HS(I)*YY(I)
          ENDDO
          DET=XJ11*XJ22-XJ21*XJ12
          IF(ISET.EQ.2) THEN
            RR=0.0
            DO I=1,NN
              RR=RR+H(I)*XX(I)
            ENDDO
            WG=WGT(NG,LX)*WGT(NG,LY)*DET*RR*2.0*PAI
          ELSE
            RR=THICK
            WG=WGT(NG,LX)*WGT(NG,LY)*DET*RR
          END IF
          DUM=1.0/DET
          XJI11= XJ22*DUM
          XJI12=-XJ12*DUM
          XJI21=-XJ21*DUM
          XJI22= XJ11*DUM

          ind1 = 1
          DO J=1,NN*NDOF
            ind2 = 1
            DO J2=1,J
              NUM=(J-1)*J/2+J2
              IF(J.EQ.J2) THEN
                SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
              ENDIF
              IF(MOD(J2,ndof).EQ.MOD(J,ndof)) THEN
                totmass=totmass+H(ind1)*H(ind2)*WG*RHO
              ENDIF
              IF(MOD(J2,NDOF).EQ.0) ind2 = ind2 + 1
            ENDDO
            IF(MOD(J,NDOF).EQ.0) ind1 = ind1 + 1
          ENDDO
!C
        ENDDO
      ENDDO
!C
      DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
      END DO
      IF(2*totmass-totdiag .LT. 0) STOP "Mass error 241!"

      RETURN
      END SUBROUTINE MASS_C2D4
!C*--------------------------------------------------------------------*
      SUBROUTINE MASS_C2D3(XX,YY,EE,PP,RHO,PARAM1,SS,ISET,myEIG)
!C*--------------------------------------------------------------------*
!C*
!C* CALCULATION 2D 3 NODE PLANE ELEMENT
!C*
      use hecmw
      use gauss_integration
      use lczparm
      IMPLICIT NONE
!C*I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),SS(*),EE,PP,PARAM1
      INTEGER(kind=kint) ISET
!C*LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=3,NDOF=2,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=2)
      REAL(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN)
      INTEGER(kind=kint) L1,L2
      REAL(kind=kreal) X1,X2,X3,XL1,XL2
      REAL(kind=kreal) DNDX,DNDY
      REAL(kind=kreal) RHO,THICK,COEF1,COEF2,PAI
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,RR,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI12,XJI22
      INTEGER(kind=kint) I,J,J2,K,LX,LY,NUM
      INTEGER(kind=kint) ind1,ind2
      REAL(kind=kreal) totdiag, totmass
      TYPE(lczparam) :: myEIG
      totdiag = 0.0
      totmass = 0.0
!C**************************
!C*  CONSTANT
!C**************************
      PAI=4.0*ATAN(1.0)
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
!C*THICKNESS 
      THICK=PARAM1
!C*FOR AX-SYM. ANALYSIS
      IF(ISET.EQ.2) THEN
        THICK=1.0
      END IF
!C*CLEAR [D] MATRIX
      DO J=1,4
        DO I=1,4
          D(I,J)=0.0
        ENDDO
      ENDDO
!C*LOOP OVER ALL INTEGRATION POINTS
      DO L2=1,NG
        XL2=XG(NG,L2)
        X2 =(XL2+1.0)*0.5
        DO L1=1,NG
          XL1=XG(NG,L1)
          X1=0.5*(1.0-X2)*(XL1+1.0)
!C*INTERPOLATION FUNCTION
          X3=1.0-X1-X2
          H(1)=X1
          H(2)=X2
          H(3)=X3
!C*DERIVATIVE OF INTERPOLATION FUNCTION
!C*FOR L1-COORDINATE
          HL1(1)=1.0
          HL1(2)=0.0
          HL1(3)=0.0
!C*FOR L2-COORDINATE
          HL2(1)=0.0
          HL2(2)=1.0
          HL2(3)=0.0
!C*FOR L3-COORDINATE
          HL3(1)=0.0
          HL3(2)=0.0
          HL3(3)=1.0
!C*JACOBI MATRIX
          XJ11=0.0
          XJ21=0.0
          XJ12=0.0
          XJ22=0.0
          DO I=1,NN
            XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
            XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
            XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
            XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
          ENDDO
          DET=XJ11*XJ22-XJ21*XJ12
          IF(ISET.EQ.2) THEN
            RR=0.0
          DO I=1,NN
            RR=RR+H(I)*XX(I)
          ENDDO
            WG=WGT(NG,L1)*WGT(NG,L2)*DET*RR*2.0*PAI*(1.0-X2)*0.25
          ELSE
            RR=THICK
            WG=WGT(NG,L1)*WGT(NG,L2)*DET*RR*(1.0-X2)*0.25
          END IF
!C
          ind1 = 1
          DO J=1,NN*NDOF
            ind2 = 1
            DO J2=1,J
              NUM=(J-1)*J/2+J2
!C
              IF(J.EQ.J2) THEN
                SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
              ENDIF
              IF(MOD(J2,ndof).EQ.MOD(J,ndof)) THEN
                totmass=totmass+H(ind1)*H(ind2)*WG*RHO
              ENDIF
              IF(MOD(J2,NDOF).EQ.0) ind2 = ind2 + 1
            ENDDO
            IF(MOD(J,NDOF).EQ.0) ind1 = ind1 + 1
          ENDDO
!C
        ENDDO
      ENDDO
!C
      DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
      END DO
      IF(2*totmass-totdiag .LT. 0) STOP "Mass error 231!"
!C
      RETURN
      END SUBROUTINE MASS_C2D3

end module m_eigen_LIB_2d1mass
