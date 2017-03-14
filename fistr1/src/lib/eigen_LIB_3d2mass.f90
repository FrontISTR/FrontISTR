!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  This module contains subroutines used in 3d eigen analysis for
!!   elements 342, 352, 362
module m_eigen_LIB_3d2mass

contains



!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D20(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 20-NODE SOLID ELEMENT
!
      use hecmw
      USE gauss_integration
      USE lczparm
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),SS(*),EE,PP,RHO
! LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=20,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=3)
      REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,J,L,LX,LY,LZ,NUM,J2
      REAL(kind=kreal) DNDX,DNDY,DNDZ
!*EHM CONSISTENT MASS MATRIX 18Apr2004
      REAL(kind=kreal) totdiag, totmass
      INTEGER(kind=kint) ind1, ind2
      TYPE(lczparam) :: myEIG

      totdiag = 0.0
      totmass = 0.0
! ZERO CLEAR MATRIX S()
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
! LOOP FOR INTEGRATION POINTS
      DO LX=1,NG
        RI=XG(NG,LX)
        DO LY=1,NG
          SI=XG(NG,LY)
          DO LZ=1,NG
            TI=XG(NG,LZ)
            RP=1.0+RI
            SP=1.0+SI
            TP=1.0+TI
            RM=1.0-RI
            SM=1.0-SI
            TM=1.0-TI
!  INTERPOLATION FUNCTION
            H(1)=-0.125*RM*SM*TM*(2.0+RI+SI+TI)
            H(2)=-0.125*RP*SM*TM*(2.0-RI+SI+TI)
            H(3)=-0.125*RP*SP*TM*(2.0-RI-SI+TI)
            H(4)=-0.125*RM*SP*TM*(2.0+RI-SI+TI)
            H(5)=-0.125*RM*SM*TP*(2.0+RI+SI-TI)
            H(6)=-0.125*RP*SM*TP*(2.0-RI+SI-TI)
            H(7)=-0.125*RP*SP*TP*(2.0-RI-SI-TI)
            H(8)=-0.125*RM*SP*TP*(2.0+RI-SI-TI)
            H(9)=0.25*(1.0-RI**2)*SM*TM
            H(10)=0.25*RP*(1.0-SI**2)*TM
            H(11)=0.25*(1.0-RI**2)*SP*TM
            H(12)=0.25*RM*(1.0-SI**2)*TM
            H(13)=0.25*(1.0-RI**2)*SM*TP
            H(14)=0.25*RP*(1.0-SI**2)*TP
            H(15)=0.25*(1.0-RI**2)*SP*TP
            H(16)=0.25*RM*(1.0-SI**2)*TP
            H(17)=0.25*RM*SM*(1.0-TI**2)
            H(18)=0.25*RP*SM*(1.0-TI**2)
            H(19)=0.25*RP*SP*(1.0-TI**2)
            H(20)=0.25*RM*SP*(1.0-TI**2)
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR R-COORDINATE
            HR(1)=-0.125*RM*SM*TM+0.125*SM*TM*(2.0+RI+SI+TI)
            HR(2)=+0.125*RP*SM*TM-0.125*SM*TM*(2.0-RI+SI+TI)
            HR(3)=+0.125*RP*SP*TM-0.125*SP*TM*(2.0-RI-SI+TI)
            HR(4)=-0.125*RM*SP*TM+0.125*SP*TM*(2.0+RI-SI+TI)
            HR(5)=-0.125*RM*SM*TP+0.125*SM*TP*(2.0+RI+SI-TI)
            HR(6)=+0.125*RP*SM*TP-0.125*SM*TP*(2.0-RI+SI-TI)
            HR(7)=+0.125*RP*SP*TP-0.125*SP*TP*(2.0-RI-SI-TI)
            HR(8)=-0.125*RM*SP*TP+0.125*SP*TP*(2.0+RI-SI-TI)
            HR(9 )=-0.50*RI*SM*TM
            HR(10)=+0.25*(1.0-SI**2)*TM
            HR(11)=-0.50*RI*SP*TM
            HR(12)=-0.25*(1.0-SI**2)*TM
            HR(13)=-0.50*RI*SM*TP
            HR(14)=+0.25*(1.0-SI**2)*TP
            HR(15)=-0.50*RI*SP*TP
            HR(16)=-0.25*(1.0-SI**2)*TP
            HR(17)=-0.25*SM*(1.0-TI**2)
            HR(18)=+0.25*SM*(1.0-TI**2)
            HR(19)=+0.25*SP*(1.0-TI**2)
            HR(20)=-0.25*SP*(1.0-TI**2)
!  FOR S-COORDINATE
            HS(1)=-0.125*RM*SM*TM+0.125*RM*TM*(2.0+RI+SI+TI)
            HS(2)=-0.125*RP*SM*TM+0.125*RP*TM*(2.0-RI+SI+TI)
            HS(3)=+0.125*RP*SP*TM-0.125*RP*TM*(2.0-RI-SI+TI)
            HS(4)=+0.125*RM*SP*TM-0.125*RM*TM*(2.0+RI-SI+TI)
            HS(5)=-0.125*RM*SM*TP+0.125*RM*TP*(2.0+RI+SI-TI)
            HS(6)=-0.125*RP*SM*TP+0.125*RP*TP*(2.0-RI+SI-TI)
            HS(7)=+0.125*RP*SP*TP-0.125*RP*TP*(2.0-RI-SI-TI)
            HS(8)=+0.125*RM*SP*TP-0.125*RM*TP*(2.0+RI-SI-TI)
            HS(9)=-0.25*(1.0-RI**2)*TM
            HS(10)=-0.50*RP*SI*TM
            HS(11)=+0.25*(1.0-RI**2)*TM
            HS(12)=-0.50*RM*SI*TM
            HS(13)=-0.25*(1.0-RI**2)*TP
            HS(14)=-0.50*RP*SI*TP
            HS(15)=+0.25*(1.0-RI**2)*TP
            HS(16)=-0.50*RM*SI*TP
            HS(17)=-0.25*RM*(1.0-TI**2)
            HS(18)=-0.25*RP*(1.0-TI**2)
            HS(19)=+0.25*RP*(1.0-TI**2)
            HS(20)=+0.25*RM*(1.0-TI**2)
!  FOR T-COORDINATE
            HT(1)=-0.125*RM*SM*TM+0.125*RM*SM*(2.0+RI+SI+TI)
            HT(2)=-0.125*RP*SM*TM+0.125*RP*SM*(2.0-RI+SI+TI)
            HT(3)=-0.125*RP*SP*TM+0.125*RP*SP*(2.0-RI-SI+TI)
            HT(4)=-0.125*RM*SP*TM+0.125*RM*SP*(2.0+RI-SI+TI)
            HT(5)=+0.125*RM*SM*TP-0.125*RM*SM*(2.0+RI+SI-TI)
            HT(6)=+0.125*RP*SM*TP-0.125*RP*SM*(2.0-RI+SI-TI)
            HT(7)=+0.125*RP*SP*TP-0.125*RP*SP*(2.0-RI-SI-TI)
            HT(8)=+0.125*RM*SP*TP-0.125*RM*SP*(2.0+RI-SI-TI)
            HT(9)=-0.25*(1.0-RI**2)*SM
            HT(10)=-0.25*RP*(1.0-SI**2)
            HT(11)=-0.25*(1.0-RI**2)*SP
            HT(12)=-0.25*RM*(1.0-SI**2)
            HT(13)=0.25*(1.0-RI**2)*SM
            HT(14)=0.25*RP*(1.0-SI**2)
            HT(15)=0.25*(1.0-RI**2)*SP
            HT(16)=0.25*RM*(1.0-SI**2)
            HT(17)=-0.5*RM*SM*TI
            HT(18)=-0.5*RP*SM*TI
            HT(19)=-0.5*RP*SP*TI
            HT(20)=-0.5*RM*SP*TI
!  JACOBI MATRIX
            XJ11=0.0
            XJ21=0.0
            XJ31=0.0
            XJ12=0.0
            XJ22=0.0
            XJ32=0.0
            XJ13=0.0
            XJ23=0.0
            XJ33=0.0
            DO I=1,NN
              XJ11=XJ11+HR(I)*XX(I)
              XJ21=XJ21+HS(I)*XX(I)
              XJ31=XJ31+HT(I)*XX(I)
              XJ12=XJ12+HR(I)*YY(I)
              XJ22=XJ22+HS(I)*YY(I)
              XJ32=XJ32+HT(I)*YY(I)
              XJ13=XJ13+HR(I)*ZZ(I)
              XJ23=XJ23+HS(I)*ZZ(I)
              XJ33=XJ33+HT(I)*ZZ(I)
            ENDDO
!DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33   &
            +XJ12*XJ23*XJ31      &
            +XJ13*XJ21*XJ32      &
            -XJ13*XJ22*XJ31      &
            -XJ12*XJ21*XJ33      &
            -XJ11*XJ23*XJ32
! INVERSION OF JACOBIAN
            DUM=1.0/DET
            XJI11=DUM*( XJ22*XJ33-XJ23*XJ32)
            XJI21=DUM*(-XJ21*XJ33+XJ23*XJ31)
            XJI31=DUM*( XJ21*XJ32-XJ22*XJ31)
            XJI12=DUM*(-XJ12*XJ33+XJ13*XJ32)
            XJI22=DUM*( XJ11*XJ33-XJ13*XJ31)
            XJI32=DUM*(-XJ11*XJ32+XJ12*XJ31)
            XJI13=DUM*( XJ12*XJ23-XJ13*XJ22)
            XJI23=DUM*(-XJ11*XJ23+XJ13*XJ21)
            XJI33=DUM*( XJ11*XJ22-XJ12*XJ21)
!  WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,LX)*WGT(NG,LY)*WGT(NG,LZ)*DET

            ind1 = 1
            DO J=1,NN*NDOF
              ind2 = 1
              DO J2=1,J
               NUM=(J-1)*J/2+J2
!DEBUG
!               IF(myEIG%mastyp.EQ.2) THEN
!                IF(MOD(j,NDOF).EQ.MOD(j2,NDOF)) THEN
!                  SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
!                ENDIF
!               ELSE IF(myEIG%mastyp.EQ.1) THEN
                IF(j.EQ.j2) THEN
                 SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                 totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
                ENDIF
                IF(MOD(j2,ndof).EQ.MOD(j,ndof)) THEN
                 totmass=totmass+H(ind1)*H(ind2)*WG*RHO
                ENDIF
!               ENDIF
               IF(MOD(j2,NDOF).EQ.0) ind2 = ind2 + 1
              ENDDO
               IF(MOD(j,NDOF).EQ.0) ind1 = ind1 + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!DEBUG
!      IF(myEIG%mastyp.EQ.1) THEN
       DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
       END DO
       IF(2*totmass-totdiag .LT. 0) STOP "Mass error 362!"
!      ENDIF
      RETURN
      END SUBROUTINE MASS_C3D20
!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D15(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 15-NODE SOLID ELEMENT
!
      use hecmw
      USE gauss_integration
      USE lczparm
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),SS(*),EE,PP,RHO
! LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=15,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=3)
      REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HZ(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,J,L,L1,L2,LZ,NUM,J2
      REAL(kind=kreal) X1,X2,X3,XL1,XL2,ZI
      REAL(kind=kreal) DNDX,DNDY,DNDZ
!*EHM CONSISTENT MASS MATRIX 18 Apr 2004
      REAL(kind=kreal) totdiag, totmass
      INTEGER(kind=kint) ind1, ind2
      TYPE(lczparam) :: myEIG

      totdiag = 0.0
      totmass = 0.0
! ZERO CLEAR MATRIX S()
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
! LOOP FOR INTEGRATION POINTS
      DO LZ=1,NG
        ZI=XG(NG,LZ)
        DO L2=1,NG
          XL2=XG(NG,L2)
          X2 =(XL2+1.0)*0.5
          DO L1=1,NG
            XL1=XG(NG,L1)
            X1=0.5*(1.0-X2)*(XL1+1.0)
!  INTERPOLATION FUNCTION
            X3=1.0-X1-X2
            H(1)= 0.5*X1*(2.0*X1-2.-ZI)*(1.0-ZI)
            H(2)= 0.5*X2*(2.0*X2-2.-ZI)*(1.0-ZI)
            H(3)= 0.5*X3*(2.0*X3-2.-ZI)*(1.0-ZI)
            H(4)= 0.5*X1*(2.0*X1-2.+ZI)*(1.0+ZI)
            H(5)= 0.5*X2*(2.0*X2-2.+ZI)*(1.0+ZI)
            H(6)= 0.5*X3*(2.0*X3-2.+ZI)*(1.0+ZI)
            H(7)= 2.0*X1*X2*(1.0-ZI)
            H(8)= 2.0*X2*X3*(1.0-ZI)
            H(9)= 2.0*X1*X3*(1.0-ZI)
            H(10)=2.0*X1*X2*(1.0+ZI)
            H(11)=2.0*X2*X3*(1.0+ZI)
            H(12)=2.0*X1*X3*(1.0+ZI)
            H(13)=X1*(1.0-ZI**2)
            H(14)=X2*(1.0-ZI**2)
            H(15)=X3*(1.0-ZI**2)
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR L1-COORDINATE
            HL1(1)= 0.5*(4.0*X1-2.-ZI)*(1.0-ZI)
            HL1(2)= 0.0
            HL1(3)= 0.0
            HL1(4)= 0.5*(4.0*X1-2.+ZI)*(1.0+ZI)
            HL1(5)= 0.0
            HL1(6)= 0.0
            HL1(7)= 2.0*X2*(1.0-ZI)
            HL1(8)= 0.0
            HL1(9)= 2.0*X3*(1.0-ZI)
            HL1(10)=2.0*X2*(1.0+ZI)
            HL1(11)=0.0
            HL1(12)=2.0*X3*(1.0+ZI)
            HL1(13)=1.0-ZI**2
            HL1(14)=0.0
            HL1(15)=0.0
!  FOR L2-COORDINATE
            HL2(1)= 0.0
            HL2(2)= 0.5*(4.0*X2-2.-ZI)*(1.0-ZI)
            HL2(3)= 0.0
            HL2(4)= 0.0
            HL2(5)= 0.5*(4.0*X2-2.+ZI)*(1.0+ZI)
            HL2(6)= 0.0
            HL2(7)= 2.0*X1*(1.0-ZI)
            HL2(8)= 2.0*X3*(1.0-ZI)
            HL2(9)= 0.0
            HL2(10)=2.0*X1*(1.0+ZI)
            HL2(11)=2.0*X3*(1.0+ZI)
            HL2(12)=0.0
            HL2(13)=0.0
            HL2(14)=1.0-ZI**2
            HL2(15)=0.0
!  FOR L3-COORDINATE
            HL3(1)= 0.0
            HL3(2)= 0.0
            HL3(3)= 0.5*(4.0*X3-2.-ZI)*(1.0-ZI)
            HL3(4)= 0.0
            HL3(5)= 0.0
            HL3(6)= 0.5*(4.0*X3-2.+ZI)*(1.0+ZI)
            HL3(7)= 0.0
            HL3(8)= 2.0*X2*(1.0-ZI)
            HL3(9)= 2.0*X1*(1.0-ZI)
            HL3(10)=0.0
            HL3(11)=2.0*X2*(1.0+ZI)
            HL3(12)=2.0*X1*(1.0+ZI)
            HL3(13)=0.0
            HL3(14)=0.0
            HL3(15)=1.0-ZI**2
!  FOR Z-COORDINATE
            HZ(1)= 0.5*X1*(-2.0*X1+1.0+2.0*ZI)
            HZ(2)= 0.5*X2*(-2.0*X2+1.0+2.0*ZI)
            HZ(3)= 0.5*X3*(-2.0*X3+1.0+2.0*ZI)
            HZ(4)= 0.5*X1*( 2.0*X1-1.0+2.0*ZI)
            HZ(5)= 0.5*X2*( 2.0*X2-1.0+2.0*ZI)
            HZ(6)= 0.5*X3*( 2.0*X3-1.0+2.0*ZI)
            HZ(7)= -2.0*X1*X2
            HZ(8)= -2.0*X2*X3
            HZ(9)= -2.0*X1*X3
            HZ(10)= 2.0*X1*X2
            HZ(11)= 2.0*X2*X3
            HZ(12)= 2.0*X1*X3
            HZ(13)=-2.0*X1*ZI
            HZ(14)=-2.0*X2*ZI
            HZ(15)=-2.0*X3*ZI
!  JACOBI MATRIX
            XJ11=0.0
            XJ21=0.0
            XJ31=0.0
            XJ12=0.0
            XJ22=0.0
            XJ32=0.0
            XJ13=0.0
            XJ23=0.0
            XJ33=0.0
            DO I=1,NN
              XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
              XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
              XJ31=XJ31+HZ(I)*XX(I)
              XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
              XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
              XJ32=XJ32+HZ(I)*YY(I)
              XJ13=XJ13+(HL1(I)-HL3(I))*ZZ(I)
              XJ23=XJ23+(HL2(I)-HL3(I))*ZZ(I)
              XJ33=XJ33+HZ(I)*ZZ(I)
            ENDDO
!DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33    &
               +XJ12*XJ23*XJ31    &
               +XJ13*XJ21*XJ32    &
               -XJ13*XJ22*XJ31    &
               -XJ12*XJ21*XJ33    &
               -XJ11*XJ23*XJ32
! INVERSION OF JACOBIAN
            DUM=1.0/DET
            XJI11=DUM*( XJ22*XJ33-XJ23*XJ32)
            XJI21=DUM*(-XJ21*XJ33+XJ23*XJ31)
            XJI31=DUM*( XJ21*XJ32-XJ22*XJ31)
            XJI12=DUM*(-XJ12*XJ33+XJ13*XJ32)
            XJI22=DUM*( XJ11*XJ33-XJ13*XJ31)
            XJI32=DUM*(-XJ11*XJ32+XJ12*XJ31)
            XJI13=DUM*( XJ12*XJ23-XJ13*XJ22)
            XJI23=DUM*(-XJ11*XJ23+XJ13*XJ21)
            XJI33=DUM*( XJ11*XJ22-XJ12*XJ21)
!  WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,LZ)*DET*(1.0-X2)*0.25

            ind1 = 1
            DO J=1,NN*NDOF
              ind2 = 1
              DO J2=1,J
               NUM=(J-1)*J/2+J2
!DEBUG
!               IF(myEIG%mastyp.EQ.2) THEN
!                IF(MOD(j,NDOF).EQ.MOD(j2,NDOF)) THEN
!                  SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
!                ENDIF
!               ELSE IF(myEIG%mastyp.EQ.1) THEN
                IF(j.EQ.j2) THEN
                 SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                 totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
                ENDIF
                IF(MOD(j2,ndof).EQ.MOD(j,ndof)) THEN
                 totmass=totmass+H(ind1)*H(ind2)*WG*RHO
                ENDIF
!               ENDIF
               IF(MOD(j2,NDOF).EQ.0) ind2 = ind2 + 1
              ENDDO
               IF(MOD(j,NDOF).EQ.0) ind1 = ind1 + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!DEBUG
!      IF(myEIG%mastyp.EQ.1) THEN
       DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
       END DO
       IF(2*totmass-totdiag .LT. 0) STOP "Mass error 352!"
!      ENDIF
      RETURN
      END SUBROUTINE MASS_C3D15
!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D10(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 10-NODE SOLID ELEMENT
!
      use hecmw
      USE gauss_integration
      USE lczparm
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),SS(*),EE,PP,RHO
! LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=10,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=3)
      REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HL4(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,J,L,L1,L2,L3,NUM,J2
      REAL(kind=kreal) XL1,XL2,XL3
      REAL(kind=kreal) X1,X2,X3,X4
      REAL(kind=kreal) DNDX,DNDY,DNDZ
!*EHM CONSISTENT MASS MATRIX 18 Apr 2004
      REAL(kind=kreal) totdiag, totmass
      INTEGER(kind=kint) ind1, ind2
      TYPE(lczparam) :: myEIG

      totdiag = 0.0
      totmass = 0.0
! ZERO CLEAR MATRIX S()
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
! LOOP FOR INTEGRATION POINTS
      DO L3=1,NG
        XL3=XG(NG,L3)
        X3 =(XL3+1.0)*0.5
        DO L2=1,NG
          XL2=XG(NG,L2)
          X2 =(1.0-X3)*(XL2+1.0)*0.5
          DO L1=1,NG
            XL1=XG(NG,L1)
            X1=(1.0-X2-X3)*(XL1+1.0)*0.5
!  INTERPOLATION FUNCTION
            X4=1.0-X1-X2-X3
            H(1)=X1*(2.0*X1-1.0)
            H(2)=X2*(2.0*X2-1.0)
            H(3)=X3*(2.0*X3-1.0)
            H(4)=X4*(2.0*X4-1.0)
            H(5)=4.0*X1*X2
            H(6)=4.0*X2*X3
            H(7)=4.0*X1*X3
            H(8)=4.0*X1*X4
            H(9)=4.0*X2*X4
            H(10)=4.0*X3*X4
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR L1-COORDINATE
            HL1(1)= 4.0*X1-1.0
            HL1(2)= 0.0
            HL1(3)= 0.0
            HL1(4)= 0.0
            HL1(5)= 4.0*X2
            HL1(6)= 0.0
            HL1(7)= 4.0*X3
            HL1(8)= 4.0*X4
            HL1(9)= 0.0
            HL1(10)= 0.0
!  FOR L2-COORDINATE
            HL2(1)= 0.0
            HL2(2)= 4.0*X2-1.0
            HL2(3)= 0.0
            HL2(4)= 0.0
            HL2(5)= 4.0*X1
            HL2(6)= 4.0*X3
            HL2(7)= 0.0
            HL2(8)= 0.0
            HL2(9)= 4.0*X4
            HL2(10)= 0.0
!  FOR L3-COORDINATE
            HL3(1)= 0.0
            HL3(2)= 0.0
            HL3(3)= 4.0*X3-1.0
            HL3(4)= 0.0
            HL3(5)= 0.0
            HL3(6)= 4.0*X2
            HL3(7)= 4.0*X1
            HL3(8)= 0.0
            HL3(9)= 0.0
            HL3(10)= 4.0*X4
!  FOR L4-COORDINATE
            HL4(1)= 0.0
            HL4(2)= 0.0
            HL4(3)= 0.0
            HL4(4)= 4.0*X4-1.0
            HL4(5)= 0.0
            HL4(6)= 0.0
            HL4(7)= 0.0
            HL4(8)= 4.0*X1
            HL4(9)= 4.0*X2
            HL4(10)= 4.0*X3
!  JACOBI MATRIX
            XJ11=0.0
            XJ21=0.0
            XJ31=0.0
            XJ12=0.0
            XJ22=0.0
            XJ32=0.0
            XJ13=0.0
            XJ23=0.0
            XJ33=0.0
            DO I=1,NN
              XJ11=XJ11+(HL1(I)-HL4(I))*XX(I)
              XJ21=XJ21+(HL2(I)-HL4(I))*XX(I)
              XJ31=XJ31+(HL3(I)-HL4(I))*XX(I)
              XJ12=XJ12+(HL1(I)-HL4(I))*YY(I)
              XJ22=XJ22+(HL2(I)-HL4(I))*YY(I)
              XJ32=XJ32+(HL3(I)-HL4(I))*YY(I)
              XJ13=XJ13+(HL1(I)-HL4(I))*ZZ(I)
              XJ23=XJ23+(HL2(I)-HL4(I))*ZZ(I)
              XJ33=XJ33+(HL3(I)-HL4(I))*ZZ(I)
            ENDDO
!DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33    &
               +XJ12*XJ23*XJ31    &
               +XJ13*XJ21*XJ32    &
               -XJ13*XJ22*XJ31    &
               -XJ12*XJ21*XJ33    &
               -XJ11*XJ23*XJ32
! INVERSION OF JACOBIAN
            DUM=1.0/DET
            XJI11=DUM*( XJ22*XJ33-XJ23*XJ32)
            XJI21=DUM*(-XJ21*XJ33+XJ23*XJ31)
            XJI31=DUM*( XJ21*XJ32-XJ22*XJ31)
            XJI12=DUM*(-XJ12*XJ33+XJ13*XJ32)
            XJI22=DUM*( XJ11*XJ33-XJ13*XJ31)
            XJI32=DUM*(-XJ11*XJ32+XJ12*XJ31)
            XJI13=DUM*( XJ12*XJ23-XJ13*XJ22)
            XJI23=DUM*(-XJ11*XJ23+XJ13*XJ21)
            XJI33=DUM*( XJ11*XJ22-XJ12*XJ21)
!
            DET=-DET
!  WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125

            ind1 = 1
            DO J=1,NN*NDOF
              ind2 = 1
              DO J2=1,J
               NUM=(J-1)*J/2+J2
!DEBUG
!               IF(myEIG%mastyp.EQ.2) THEN
!                IF(MOD(j,NDOF).EQ.MOD(j2,NDOF)) THEN
!                 SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
!                ENDIF
!               ELSE IF(myEIG%mastyp.EQ.1) THEN
                IF(j.EQ.j2) THEN
                 SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                 totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
                ENDIF
                IF(MOD(j2,ndof).EQ.MOD(j,ndof)) THEN
                 totmass=totmass+H(ind1)*H(ind2)*WG*RHO
                ENDIF
!               ENDIF
               IF(MOD(j2,NDOF).EQ.0) ind2 = ind2 + 1
              ENDDO
               IF(MOD(j,NDOF).EQ.0) ind1 = ind1 + 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!DEBUG
!      IF(myEIG%mastyp.EQ.1) THEN
       DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
       END DO
       IF(2*totmass-totdiag .LT. 0) STOP "Mass error 342!"
!      ENDIF
      RETURN
      END SUBROUTINE MASS_C3D10

end module m_eigen_LIB_3d2mass
