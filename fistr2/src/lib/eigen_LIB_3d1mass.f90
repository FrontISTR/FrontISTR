!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
!>  This module contains subroutines used in 3d eigen analysis for
!!  elements 341, 351, 361
!C================================================================C

!***********************************************************************
!  3D SOLID Element:
!  MASS_C3D8(XX,YY,ZZ,EE,PP,SS,IFLG)
!  MASS_C3D6(XX,YY,ZZ,EE,PP,SS)
!  MASS_C3D4(XX,YY,ZZ,EE,PP,SS)
!***********************************************************************

module m_eigen_LIB_3d1mass

contains


!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D8(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 8-NODE SOLID ELEMENT
!
      use hecmw
      USE gauss_integration
      USE lczparm
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),SS(*),EE,PP,RHO
      INTEGER(kind=kint) IFLG
! LOCAL VARIABLES
      INTEGER(kind=kint) NN,NDOF,ISIZE
      INTEGER(kind=kint) NG
      PARAMETER(NN=8,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=2)
      REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,J,L,LX,LY,LZ,NUM,J2,K1
      REAL(kind=kreal) DNDX,DNDY,DNDZ
      REAL(kind=kreal) DRDX,DRDY,DRDZ
      REAL(kind=kreal) DSDX,DSDY,DSDZ
      REAL(kind=kreal) DTDX,DTDY,DTDZ
!**FOR CORRECTION 
      REAL(kind=kreal) C(6,9),DC(6,9),CC(9,9),V(9,24),CCV(9,24)
!*EHM CONSISTENT MASS MATRIX 17Apr04
      REAL(kind=kreal) totdiag, totmass
      INTEGER(kind=kint) ind1, ind2
      REAL(kind=kreal) ej2, ej, sstest(NN*NDOF,NN*NDOF)
      TYPE(lczparam) :: myEIG

      !print *,'RHO:',RHO
      !pause

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
            H(1)=0.125*RM*SM*TM
            H(2)=0.125*RP*SM*TM
            H(3)=0.125*RP*SP*TM
            H(4)=0.125*RM*SP*TM
            H(5)=0.125*RM*SM*TP
            H(6)=0.125*RP*SM*TP
            H(7)=0.125*RP*SP*TP
            H(8)=0.125*RM*SP*TP
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR R-COORDINATE
            HR(1)=-.125*SM*TM
            HR(2)= .125*SM*TM
            HR(3)= .125*SP*TM
            HR(4)=-.125*SP*TM
            HR(5)=-.125*SM*TP
            HR(6)= .125*SM*TP
            HR(7)= .125*SP*TP
            HR(8)=-.125*SP*TP
!  FOR S-COORDINATE
            HS(1)=-.125*RM*TM
            HS(2)=-.125*RP*TM
            HS(3)= .125*RP*TM
            HS(4)= .125*RM*TM
            HS(5)=-.125*RM*TP
            HS(6)=-.125*RP*TP
            HS(7)= .125*RP*TP
            HS(8)= .125*RM*TP
!  FOR T-COORDINATE
            HT(1)=-.125*RM*SM
            HT(2)=-.125*RP*SM
            HT(3)=-.125*RP*SP
            HT(4)=-.125*RM*SP
            HT(5)= .125*RM*SM
            HT(6)= .125*RP*SM
            HT(7)= .125*RP*SP
            HT(8)= .125*RM*SP
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
            DET=XJ11*XJ22*XJ33     &
               +XJ12*XJ23*XJ31     &
               +XJ13*XJ21*XJ32     &
               -XJ13*XJ22*XJ31     &
               -XJ12*XJ21*XJ33     &
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
!                IF(MOD(j,ndof).EQ.MOD(j2,ndof)) THEN
!                 SS(NUM)=SS(NUM)+H(ind2)*H(ind1)*WG*RHO
!!*EHM CONSISTENT MASS MATRIX 22 Apr 2004: Uncomment for debug
!!                 sstest(j,j2) = SS(NUM)
!                 ENDIF
!               ELSE IF(myEIG%mastyp.EQ.1) THEN
                IF(j.EQ.j2) THEN
                 SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
                 totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
                ENDIF
                IF(MOD(j2,ndof).EQ.MOD(j,ndof)) THEN
                 totmass=totmass+H(ind1)*H(ind2)*WG*RHO
                ENDIF
!               ENDIF
                IF(MOD(J2,NDOF).EQ.0) ind2 = ind2 + 1 
              ENDDO
                IF(MOD(J,NDOF).EQ.0) ind1 = ind1 + 1 
            ENDDO

!*EHM CONSISTENT MASS MATRIX 15 Apr 2004
         END DO
        END DO
       END DO

!DEBUG
!      IF(myEIG%mastyp.EQ.1) THEN
       DO j = 1,nn*ndof
        num = j*(j+1)/2
        ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
       END DO
       IF(2*totmass-totdiag .LT. 0) then
          write(*,*) "Mass error 361!"
          call hecmw_abort( hecmw_comm_get_comm())
       endif
!      ENDIF
!DEBUG

!*EHM DEBUG CONSISTENT MASS MATRIX 18Apr04: Uncomment for debug
!	OPEN(55,FILE="sstest.dat",STATUS="UNKNOWN",ACCESS="APPEND")
!        WRITE(55,*)'*--------------------------------------------*'
!	DO j=1,NN*NDOF
!	 WRITE(55,10) (sstest(j,j2),j2=1,NN*NDOF)
!	END DO
! 10     FORMAT(24F7.5)
!        CLOSE(55)
!	!PAUSE 'Check SSTEST'
      RETURN
      END SUBROUTINE MASS_C3D8
!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D6(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 6-NODE SOLID ELEMENT
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
      PARAMETER(NN=6,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=2)
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
            H(1)=X1*(1.0-ZI)*0.5
            H(2)=X2*(1.0-ZI)*0.5
            H(3)=X3*(1.0-ZI)*0.5
            H(4)=X1*(1.0+ZI)*0.5
            H(5)=X2*(1.0+ZI)*0.5
            H(6)=X3*(1.0+ZI)*0.5
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR L1-COORDINATE
            HL1(1)=(1.0-ZI)*0.5
            HL1(2)=0.0
            HL1(3)=0.0
            HL1(4)=(1.0+ZI)*0.5
            HL1(5)=0.0
            HL1(6)=0.0
!  FOR L2-COORDINATE
            HL2(1)=0.0
            HL2(2)=(1.0-ZI)*0.5
            HL2(3)=0.0
            HL2(4)=0.0
            HL2(5)=(1.0+ZI)*0.5
            HL2(6)=0.0
!  FOR L3-COORDINATE
            HL3(1)=0.0
            HL3(2)=0.0
            HL3(3)=(1.0-ZI)*0.5
            HL3(4)=0.0
            HL3(5)=0.0
            HL3(6)=(1.0+ZI)*0.5
!  FOR Z-COORDINATE
            HZ(1)=-X1*0.5
            HZ(2)=-X2*0.5
            HZ(3)=-X3*0.5
            HZ(4)= X1*0.5
            HZ(5)= X2*0.5
            HZ(6)= X3*0.5
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
            DET=XJ11*XJ22*XJ33     &
               +XJ12*XJ23*XJ31     &
               +XJ13*XJ21*XJ32     &
               -XJ13*XJ22*XJ31     &
               -XJ12*XJ21*XJ33     &
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
                IF(MOD(j2,NDOF).EQ.0) ind2 = ind2+1
              ENDDO
                IF(MOD(j,NDOF).EQ.0) ind1 = ind1+1
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
       IF(2*totmass-totdiag .LT. 0) STOP "Mass error 351!"
!      ENDIF

      RETURN
      END SUBROUTINE MASS_C3D6
!----------------------------------------------------------------------*
      SUBROUTINE MASS_C3D4(XX,YY,ZZ,EE,PP,RHO,SS,myEIG)
!----------------------------------------------------------------------*
!
! CALCULATION 3D 4-NODE SOLID ELEMENT
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
      PARAMETER(NN=4,NDOF=3,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=2)
      REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HL4(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG,DUM
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,J,L,L1,L2,L3,NUM,J2
      REAL(kind=kreal) XL1,XL2,XL3
      REAL(kind=kreal) X1,X2,X3,X4
      REAL(kind=kreal) DNDX,DNDY,DNDZ
!*EHM CONSISTENT MASS MATRIX 18 Apr 2004
      REAL(kind=kreal) totdiag,totmass
      INTEGER(kind=kint) ind1, ind2
      TYPE(lczparam) :: myEIG

      totdiag = 0.0
      totmass = 0.0
! ZERO CLEAR MATRIX S()
      DO I=1,ISIZE
        SS(I)=0.0
      ENDDO
      totdiag = 0.
      totmass = 0.
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
            H(1)=X1
            H(2)=X2
            H(3)=X3
            H(4)=X4
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR L1-COORDINATE
            HL1(1)= 1.0
            HL1(2)= 0.0
            HL1(3)= 0.0
            HL1(4)= 0.0
!  FOR L2-COORDINATE
            HL2(1)= 0.0 
            HL2(2)= 1.0 
            HL2(3)= 0.0
            HL2(4)= 0.0
!  FOR L3-COORDINATE
            HL3(1)= 0.0
            HL3(2)= 0.0
            HL3(3)= 1.0 
            HL3(4)= 0.0
!  FOR L4-COORDINATE
            HL4(1)= 0.0
            HL4(2)= 0.0
            HL4(3)= 0.0 
            HL4(4)= 1.0
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
            DET=XJ11*XJ22*XJ33     &
                 +XJ12*XJ23*XJ31   &
                 +XJ13*XJ21*XJ32   &
                 -XJ13*XJ22*XJ31   &
                 -XJ12*XJ21*XJ33   &
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
            DET=-DET
!  WEIGT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
            ind1 = 1
            DO J=1,NN*NDOF
              ind2 = 1
              DO J2=1,J
               NUM=(J-1)*J/2+J2
!DEBUG
!               IF(myEIG%mastyp.EQ.2) THEN
!                IF(MOD(j,NDOF).EQ.MOD(j2,NDOF))THEN
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
       IF(2*totmass-totdiag .LT. 0) STOP "Mass error 341!"
!      ENDIF
      RETURN
      END SUBROUTINE MASS_C3D4
!--------------------------------------------------------------------*

end module m_eigen_LIB_3d1mass
