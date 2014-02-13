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
!>   This module provide common functions of Shell elements

module m_static_LIB_shell
   contains
!***********************************************************************
!  LIBshellelem.f
!***********************************************************************
!  SHELL Element:
!  STF_S4(XX,YY,ZZ,EE,PP,PARAM1,SS)
!  STF_S3(XX,YY,ZZ,EE,PP,PARAM1,SS)
!  DL_S4(XX,YY,ZZ,RHO,PARAM1,LTYPE,PARAMS,VECT,NN)
!  DL_S3(XX,YY,ZZ,RHO,PARAM1,LTYPE,PARAMS,VECT,NN)
!  RCV_S4(XX,YY,ZZ,EE,PP,PARAM1,EDISP,RI,SI,TI,ARRAY)
!  RCV_S3(XX,YY,ZZ,EE,PP,PARAM1,EDISP,TI,ARRAY)
!***********************************************************************
!----------------------------------------------------------------------*
   SUBROUTINE STF_S4(XX,YY,ZZ,EE,PP,THICK,STIFF)
!----------------------------------------------------------------------*
! 
! CALCULATION OF 4-NODES ISOPARAMETRIC DEGENERATED SHELL ELEMENT
!
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(4),YY(4),ZZ(4),STIFF(300)
! LOCAL VARIABLES
      DIMENSION COD(3,4)
      DIMENSION INDX(4,15)
      DIMENSION G1(3),G2(3),G3(3),REF(3)
      DIMENSION G1A(3),G2A(3),G3A(3)
      DIMENSION G1B(3),G2B(3),G3B(3)
      DIMENSION G1C(3),G2C(3),G3C(3)
      DIMENSION G1D(3),G2D(3),G3D(3)
      DIMENSION PHA12(4),PHA13(4),PHA23(4)
      DIMENSION B1(3,24),B2(3,24),B3(3,24)
      DIMENSION B1A(3,24),B2A(3,24),B3A(3,24)
      DIMENSION B1B(3,24),B2B(3,24),B3B(3,24)
      DIMENSION B1C(3,24),B2C(3,24),B3C(3,24)
      DIMENSION B1D(3,24),B2D(3,24),B3D(3,24)
      DIMENSION B11(24),B22(24),B12(24),B13(24),B23(24)
      DIMENSION EN(3,4),THE(3,3),AMAT(3,3)
      DIMENSION H(4),HR(4),HS(4),XG(2),WGT(2)
      DIMENSION ELAS(9),XELAS(15)
!
! TUNING PARAMETER
!
      DATA SCF/1.2/
      DATA ILOC/1/
!
! SET CORD
!
      DO I=1,4
        COD(1,I)=XX(I)
        COD(2,I)=YY(I)
        COD(3,I)=ZZ(I)
      ENDDO
!
! SET INDEX
!
      CALL S4INDX(INDX)
!
! SET GAUUUSIAN INTEGRATION PARAMETER
!
      XG(1) =-0.5773502691896258
      XG(2) =-XG(1)
      WGT(1)=1.0
      WGT(2)=WGT(1)
!
! ZERO CLEAR ARRAY STIFF()
!
      DO I=1,300
        STIFF(I)=0.0
      ENDDO
!
! SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM
!
      CALL S4REFV(COD,REF,ILOC)
!
! CALCULATE [PHAI] AT EACH NODE
!
      CALL S4PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!
! LOOP OVER LAMINATES
!
! ELASTICITY MATRIX
      ELAS(1)=EE/(1.0-PP*PP)
      ELAS(2)=EE/(1.0-PP*PP)*PP
      ELAS(3)=EE/(1.0-PP*PP)
      ELAS(4)=0.0
      ELAS(5)=0.0
      ELAS(6)=EE/(1.0+PP)/2.0
      ELAS(7)=EE/(1.0+PP)/2.0/SCF
      ELAS(8)=0.0
      ELAS(9)=EE/(1.0+PP)/2.0/SCF
!
!   LOOP FOR GAUSS INTEGRATION POINT
!
      DO IG3=1,2
        TI=XG(IG3)
!
! METRICS AT A,B,C,D
!
        A0= 0.0
        AP= 1.0
        AM=-1.0
        CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    A0,AP, TI,G1A,G2A,G3A,B1A,B2A,B3A)
        CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    AM,A0, TI,G1B,G2B,G3B,B1B,B2B,B3B)
        CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    A0,AM, TI,G1C,G2C,G3C,B1C,B2C,B3C)
        CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    AP,A0, TI,G1D,G2D,G3D,B1D,B2D,B3D)
        DO IG2=1,2
          SI=XG(IG2)
          DO IG1=1,2
            RI=XG(IG1)
            CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,RI,SI,TI,G1,G2,G3,B1,B2,B3)
! JACOBI MATRIX 
            CALL S4JAC(G1,G2,G3,DET,AMAT)
! SET INTEGRATION WEIGHT
            WEIGHT=WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
! SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
            CALL S4THE(G3,REF,THE)
! SET STRAIN IN NATURAL COORDINATE
!
            DO I=1,24
              B11(I)=G1(1)*B1(1,I)+G1(2)*B1(2,I)+G1(3)*B1(3,I)
              B22(I)=G2(1)*B2(1,I)+G2(2)*B2(2,I)+G2(3)*B2(3,I)
              B12(I)=(G2(1)*B1(1,I)+G2(2)*B1(2,I)+G2(3)*B1(3,I)                &
                    +G1(1)*B2(1,I)+G1(2)*B2(2,I)+G1(3)*B2(3,I))
              B13(I)=0.5*(1.0+SI)*                                             &
                    (G3A(1)*B1A(1,I)+G3A(2)*B1A(2,I)+G3A(3)*B1A(3,I)           &
                    +G1A(1)*B3A(1,I)+G1A(2)*B3A(2,I)+G1A(3)*B3A(3,I))          &
                    +0.5*(1.0-SI)*                                             &
                    (G3C(1)*B1C(1,I)+G3C(2)*B1C(2,I)+G3C(3)*B1C(3,I)           &
                    +G1C(1)*B3C(1,I)+G1C(2)*B3C(2,I)+G1C(3)*B3C(3,I))
              B23(I)=0.5*(1.0+RI)*                                             &
                     (G3D(1)*B2D(1,I)+G3D(2)*B2D(2,I)+G3D(3)*B2D(3,I)          &
                     +G2D(1)*B3D(1,I)+G2D(2)*B3D(2,I)+G2D(3)*B3D(3,I))         &
                     +0.5*(1.0-RI)*                                            &
                     (G3B(1)*B2B(1,I)+G3B(2)*B2B(2,I)+G3B(3)*B2B(3,I)          &
                     +G2B(1)*B3B(1,I)+G2B(2)*B3B(2,I)+G2B(3)*B3B(3,I)) 
            enddo 
!
! CONSTITUTIVE LAW
            CALL S4CONS(THE,G1,G2,G3,DET,ELAS,INDX,XELAS)
! ELEMENT STIFFNESS
!
            DO JJ=1,24
              DO II=1,JJ
                NPOIN=JJ*(JJ-1)/2+II
                VAL=                                                           &
                  +XELAS( 1)*B11(II)*B11(JJ)                                   &
                  +XELAS( 2)*(B11(II)*B22(JJ)+B22(II)*B11(JJ))                 &
                  +XELAS( 3)*B22(II)*B22(JJ)                                   &
                  +XELAS( 4)*(B11(II)*B12(JJ)+B12(II)*B11(JJ))                 &
                  +XELAS( 5)*(B22(II)*B12(JJ)+B12(II)*B22(JJ))                 &
                  +XELAS( 6)* B12(II)*B12(JJ)                                  &
                  +XELAS( 7)*(B11(II)*B13(JJ)+B13(II)*B11(JJ))                 &
                  +XELAS( 8)*(B22(II)*B13(JJ)+B13(II)*B22(JJ))                 &
                  +XELAS( 9)*(B12(II)*B13(JJ)+B13(II)*B12(JJ))                 &
                  +XELAS(10)* B13(II)*B13(JJ)                                  &
                  +XELAS(11)*(B11(II)*B23(JJ)+B23(II)*B11(JJ))                 &
                  +XELAS(12)*(B22(II)*B23(JJ)+B23(II)*B22(JJ))                 &
                  +XELAS(13)*(B12(II)*B23(JJ)+B23(II)*B12(JJ))                 &
                  +XELAS(14)*(B13(II)*B23(JJ)+B23(II)*B13(JJ))                 &
                  +XELAS(15)* B23(II)*B23(JJ)
                  STIFF(NPOIN)=STIFF(NPOIN)+VAL*WEIGHT
              enddo
            enddo
!
!  ADD TORSIONAL STIFFNESS 
!
            CALL S4PSET(RI,SI,H,HR,HS)
            CALL S4TORT(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,STIFF,ELAS(6),WEIGHT)
!
          enddo
        enddo
      enddo
      RETURN

   end subroutine STF_S4
!----------------------------------------------------------------------*
   SUBROUTINE S4CONS(THE,G1,G2,G3,DET,ELAS,INDX,XELAS)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION THE(3,3),G1(3),G2(3),G3(3),ELAS(9),INDX(4,15),XELAS(15)
      DIMENSION GG1(3),GG2(3),GG3(3),A(3,3)
! CONTRAVARIANT BASE VECTOR
      GG1(1)=(G2(2)*G3(3)-G2(3)*G3(2))/DET
      GG1(2)=(G2(3)*G3(1)-G2(1)*G3(3))/DET
      GG1(3)=(G2(1)*G3(2)-G2(2)*G3(1))/DET
      GG2(1)=(G3(2)*G1(3)-G3(3)*G1(2))/DET
      GG2(2)=(G3(3)*G1(1)-G3(1)*G1(3))/DET
      GG2(3)=(G3(1)*G1(2)-G3(2)*G1(1))/DET
      GG3(1)=(G1(2)*G2(3)-G1(3)*G2(2))/DET
      GG3(2)=(G1(3)*G2(1)-G1(1)*G2(3))/DET
      GG3(3)=(G1(1)*G2(2)-G1(2)*G2(1))/DET
      DO I=1,3
        A(1,I)=GG1(1)*THE(1,I)+GG1(2)*THE(2,I)+GG1(3)*THE(3,I)
        A(2,I)=GG2(1)*THE(1,I)+GG2(2)*THE(2,I)+GG2(3)*THE(3,I)
        A(3,I)=GG3(1)*THE(1,I)+GG3(2)*THE(2,I)+GG3(3)*THE(3,I)
      enddo
      DO I=1,15
        I1=INDX(1,I)
        I2=INDX(2,I)
        I3=INDX(3,I)
        I4=INDX(4,I)
        XELAS(I)=A(I1,1)*A(I2,1)*A(I3,1)*A(I4,1)*ELAS(1)                       &
                +A(I1,1)*A(I2,1)*A(I3,2)*A(I4,2)*ELAS(2)                       &
                +A(I1,2)*A(I2,2)*A(I3,1)*A(I4,1)*ELAS(2)                       &
                +A(I1,2)*A(I2,2)*A(I3,2)*A(I4,2)*ELAS(3)                       &
                +A(I1,1)*A(I2,1)*A(I3,1)*A(I4,2)*ELAS(4)                       &
                +A(I1,1)*A(I2,2)*A(I3,1)*A(I4,1)*ELAS(4)                       &
                +A(I1,2)*A(I2,1)*A(I3,1)*A(I4,1)*ELAS(4)                       &
                +A(I1,1)*A(I2,1)*A(I3,2)*A(I4,1)*ELAS(4)                       &
                +A(I1,2)*A(I2,2)*A(I3,1)*A(I4,2)*ELAS(5)                       &
                +A(I1,2)*A(I2,2)*A(I3,2)*A(I4,1)*ELAS(5)                       &
                +A(I1,1)*A(I2,2)*A(I3,2)*A(I4,2)*ELAS(5)                       &
                +A(I1,2)*A(I2,1)*A(I3,2)*A(I4,2)*ELAS(5)                       &
                +A(I1,1)*A(I2,2)*A(I3,1)*A(I4,2)*ELAS(6)                       &
                +A(I1,2)*A(I2,1)*A(I3,2)*A(I4,1)*ELAS(6)                       &
                +A(I1,1)*A(I2,2)*A(I3,2)*A(I4,1)*ELAS(6)                       &
                +A(I1,2)*A(I2,1)*A(I3,1)*A(I4,2)*ELAS(6)
        XELAS(I)=XELAS(I)                                                      &
                +A(I1,1)*A(I2,3)*A(I3,1)*A(I4,3)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,3)*A(I4,1)*ELAS(7)                       &
                +A(I1,1)*A(I2,3)*A(I3,3)*A(I4,1)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,1)*A(I4,3)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,3)*A(I4,2)*ELAS(8)                       &
                +A(I1,3)*A(I2,1)*A(I3,2)*A(I4,3)*ELAS(8)                       &
                +A(I1,1)*A(I2,3)*A(I3,3)*A(I4,2)*ELAS(8)                       &
                +A(I1,1)*A(I2,3)*A(I3,2)*A(I4,3)*ELAS(8)                       &
                +A(I1,3)*A(I2,2)*A(I3,3)*A(I4,1)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,3)*A(I4,1)*ELAS(8)                       &
                +A(I1,3)*A(I2,2)*A(I3,1)*A(I4,3)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,1)*A(I4,3)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,2)*A(I4,3)*ELAS(9)                       &
                +A(I1,3)*A(I2,2)*A(I3,3)*A(I4,2)*ELAS(9)                       &
                +A(I1,2)*A(I2,3)*A(I3,3)*A(I4,2)*ELAS(9)                       &
                +A(I1,3)*A(I2,2)*A(I3,2)*A(I4,3)*ELAS(9)
      enddo
      RETURN

   end subroutine S4CONS
!----------------------------------------------------------------------*
   SUBROUTINE S4COV(COD,EN,H,HR,HS,TI,THICK,G1,G2,G3)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,4),EN(3,4),H(4),HR(4),HS(4),G1(3),G2(3),G3(3)
      DO I=1,3
        G1(I)=0.0
        G2(I)=0.0
        G3(I)=0.0
        DO INOD=1,4
          VAR=COD(I,INOD)+THICK*0.5*TI*EN(I,INOD)
          G1(I)=G1(I)+HR(INOD)*VAR
          G2(I)=G2(I)+HS(INOD)*VAR
          G3(I)=G3(I)+THICK*0.5*H(INOD)*EN(I,INOD)
        enddo
      enddo
      RETURN

   end subroutine S4COV
!----------------------------------------------------------------------*
      SUBROUTINE S4INDX(INDX)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION INDX(4,15)
      INDX(1,1)=1
      INDX(2,1)=1
      INDX(3,1)=1
      INDX(4,1)=1
      INDX(1,2)=1
      INDX(2,2)=1
      INDX(3,2)=2
      INDX(4,2)=2
      INDX(1,3)=2
      INDX(2,3)=2
      INDX(3,3)=2
      INDX(4,3)=2
      INDX(1,4)=1
      INDX(2,4)=1
      INDX(3,4)=1
      INDX(4,4)=2
      INDX(1,5)=2
      INDX(2,5)=2
      INDX(3,5)=1
      INDX(4,5)=2
      INDX(1,6)=1
      INDX(2,6)=2
      INDX(3,6)=1
      INDX(4,6)=2
      INDX(1,7)=1
      INDX(2,7)=1
      INDX(3,7)=1
      INDX(4,7)=3
      INDX(1,8)=2
      INDX(2,8)=2
      INDX(3,8)=1
      INDX(4,8)=3
      INDX(1,9)=1
      INDX(2,9)=2
      INDX(3,9)=1
      INDX(4,9)=3
      INDX(1,10)=1
      INDX(2,10)=3
      INDX(3,10)=1
      INDX(4,10)=3
      INDX(1,11)=1
      INDX(2,11)=1
      INDX(3,11)=2
      INDX(4,11)=3
      INDX(1,12)=2
      INDX(2,12)=2
      INDX(3,12)=2
      INDX(4,12)=3
      INDX(1,13)=1
      INDX(2,13)=2
      INDX(3,13)=2
      INDX(4,13)=3
      INDX(1,14)=1
      INDX(2,14)=3
      INDX(3,14)=2
      INDX(4,14)=3
      INDX(1,15)=2
      INDX(2,15)=3
      INDX(3,15)=2
      INDX(4,15)=3
      RETURN

   end subroutine S4INDX
!----------------------------------------------------------------------*
   SUBROUTINE S4JAC(G1,G2,G3,DET,AMAT)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION G1(3),G2(3),G3(3),AMAT(3,3)
      XJ11=G1(1)
      XJ12=G1(2)  
      XJ13=G1(3)  
      XJ21=G2(1)
      XJ22=G2(2)  
      XJ23=G2(3)  
      XJ31=G3(1)
      XJ32=G3(2)  
      XJ33=G3(3)  
! DETERMINANT OF JACOBIAN
      DET=XJ11*XJ22*XJ33                                                       &
         +XJ12*XJ23*XJ31                                                       &
         +XJ13*XJ21*XJ32                                                       &
         -XJ13*XJ22*XJ31                                                       &
         -XJ12*XJ21*XJ33                                                       &
         -XJ11*XJ23*XJ32
! INVERSION OF JACOBIAN
      DUM=1.0/DET
      AMAT(1,1)=DUM*( XJ22*XJ33-XJ23*XJ32)
      AMAT(2,1)=DUM*(-XJ21*XJ33+XJ23*XJ31)
      AMAT(3,1)=DUM*( XJ21*XJ32-XJ22*XJ31)
      AMAT(1,2)=DUM*(-XJ12*XJ33+XJ13*XJ32)
      AMAT(2,2)=DUM*( XJ11*XJ33-XJ13*XJ31)
      AMAT(3,2)=DUM*(-XJ11*XJ32+XJ12*XJ31)
      AMAT(1,3)=DUM*( XJ12*XJ23-XJ13*XJ22)
      AMAT(2,3)=DUM*(-XJ11*XJ23+XJ13*XJ21)
      AMAT(3,3)=DUM*( XJ11*XJ22-XJ12*XJ21)
      RETURN

   end subroutine S4JAC
!----------------------------------------------------------------------*
      SUBROUTINE S4PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,4),REF(3),PHA12(4),PHA13(4),PHA23(4),EN(3,4)
      DIMENSION G1(3),G2(3),G3(3),E1(3),E2(3),E3(3)
      DO I=1,4
        IF(I.EQ.1) THEN
          G1(1)=-0.5*COD(1,1)+0.5*COD(1,2)
          G1(2)=-0.5*COD(2,1)+0.5*COD(2,2)
          G1(3)=-0.5*COD(3,1)+0.5*COD(3,2)
          G2(1)=-0.5*COD(1,1)+0.5*COD(1,4)
          G2(2)=-0.5*COD(2,1)+0.5*COD(2,4)
          G2(3)=-0.5*COD(3,1)+0.5*COD(3,4)
        ELSE IF(I.EQ.2) THEN
          G1(1)=-0.5*COD(1,1)+0.5*COD(1,2)
          G1(2)=-0.5*COD(2,1)+0.5*COD(2,2)
          G1(3)=-0.5*COD(3,1)+0.5*COD(3,2)
          G2(1)=-0.5*COD(1,2)+0.5*COD(1,3)
          G2(2)=-0.5*COD(2,2)+0.5*COD(2,3)
          G2(3)=-0.5*COD(3,2)+0.5*COD(3,3)
        ELSE IF(I.EQ.3) THEN
          G1(1)= 0.5*COD(1,3)-0.5*COD(1,4)
          G1(2)= 0.5*COD(2,3)-0.5*COD(2,4)
          G1(3)= 0.5*COD(3,3)-0.5*COD(3,4)
          G2(1)=-0.5*COD(1,2)+0.5*COD(1,3)
          G2(2)=-0.5*COD(2,2)+0.5*COD(2,3)
          G2(3)=-0.5*COD(3,2)+0.5*COD(3,3)
        ELSE IF(I.EQ.4) THEN
          G1(1)= 0.5*COD(1,3)-0.5*COD(1,4)
          G1(2)= 0.5*COD(2,3)-0.5*COD(2,4)
          G1(3)= 0.5*COD(3,3)-0.5*COD(3,4)
          G2(1)=-0.5*COD(1,1)+0.5*COD(1,4)
          G2(2)=-0.5*COD(2,1)+0.5*COD(2,4)
          G2(3)=-0.5*COD(3,1)+0.5*COD(3,4)
        END IF
! G3()=G1() X G2()
        G3(1)=G1(2)*G2(3)-G1(3)*G2(2)
        G3(2)=G1(3)*G2(1)-G1(1)*G2(3)
        G3(3)=G1(1)*G2(2)-G1(2)*G2(1)
!  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
        XSUM=DSQRT(G3(1)**2+G3(2)**2+G3(3)**2)
        E3(1)=G3(1)/XSUM
        E3(2)=G3(2)/XSUM
        E3(3)=G3(3)/XSUM
        EN(1,I)=E3(1)
        EN(2,I)=E3(2)
        EN(3,I)=E3(3)
        E2(1)=-REF(2)*E3(3)+REF(3)*E3(2)
        E2(2)=-REF(3)*E3(1)+REF(1)*E3(3)
        E2(3)=-REF(1)*E3(2)+REF(2)*E3(1)
        XSUM=DSQRT(E2(1)**2+E2(2)**2+E2(3)**2)
        IF (XSUM.GT.1.E-15) THEN
          E2(1)=E2(1)/XSUM
          E2(2)=E2(2)/XSUM
          E2(3)=E2(3)/XSUM
          E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
          E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
          E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
          XSUM=DSQRT(E1(1)**2+E1(2)**2+E1(3)**2)
          E1(1)=E1(1)/XSUM
          E1(2)=E1(2)/XSUM
          E1(3)=E1(3)/XSUM
        ELSE
          E1(1)=  0.0
          E1(2)=  0.0
          E1(3)=- 1.0
          E2(1) = 0.0
          E2(2) = 1.0
          E2(3) = 0.0
        END IF
        PHA12(I)=-E2(1)*E1(2)+E1(1)*E2(2)
        PHA13(I)=-E2(1)*E1(3)+E1(1)*E2(3)
        PHA23(I)=-E2(2)*E1(3)+E1(2)*E2(3)
      enddo
      RETURN

   end subroutine S4PHAI
!----------------------------------------------------------------------*
   SUBROUTINE S4PSET(RI,SI,H,HR,HS)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION H(4),HR(4),HS(4)
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
      RETURN
   
   end subroutine S4PSET
!----------------------------------------------------------------------*
   SUBROUTINE S4REFV(COD,REF,ILOC)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,4),REF(3)
      IF(ILOC.EQ.0) THEN
        REF(1)=1.0  
        REF(2)=0.0 
        REF(3)=0.0  
      ELSE 
        REF(1)=0.25*(COD(1,2)+COD(1,3)-COD(1,1)-COD(1,4))
        REF(2)=0.25*(COD(2,2)+COD(2,3)-COD(2,1)-COD(2,4))
        REF(3)=0.25*(COD(3,2)+COD(3,3)-COD(3,1)-COD(3,4))
      END IF
      RETURN

   end subroutine S4REFV
!----------------------------------------------------------------------*
   SUBROUTINE S4THE(G3,REF,THE)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION G3(3),REF(3),THE(3,3)
      DIMENSION E1(3),E2(3),E3(3)
      XSUM=DSQRT(G3(1)**2+G3(2)**2+G3(3)**2)
      E3(1)=G3(1)/XSUM
      E3(2)=G3(2)/XSUM
      E3(3)=G3(3)/XSUM
      E2(1)=-REF(2)*E3(3)+REF(3)*E3(2)
      E2(2)=-REF(3)*E3(1)+REF(1)*E3(3)
      E2(3)=-REF(1)*E3(2)+REF(2)*E3(1)
      E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
      E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
      E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
      XSUM=DSQRT(E2(1)**2+E2(2)**2+E2(3)**2)
      IF (XSUM.GT.1.E-15) THEN
        E2(1)=E2(1)/XSUM
        E2(2)=E2(2)/XSUM
        E2(3)=E2(3)/XSUM
        E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
        E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
        E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
        XSUM=DSQRT(E1(1)**2+E1(2)**2+E1(3)**2)
        E1(1)=E1(1)/XSUM
        E1(2)=E1(2)/XSUM
        E1(3)=E1(3)/XSUM
      ELSE
        E1(1)=  0.0
        E1(2)=  0.0
        E1(3)= -1.0
        E2(1) = 0.0
        E2(2) = 1.0
        E2(3) = 0.0
      END IF
      THE(1,1)=E1(1)
      THE(1,2)=E2(1)
      THE(1,3)=E3(1)
      THE(2,1)=E1(2)
      THE(2,2)=E2(2)
      THE(2,3)=E3(2)
      THE(3,1)=E1(3)
      THE(3,2)=E2(3)
      THE(3,3)=E3(3)
      RETURN

   end subroutine S4THE
!----------------------------------------------------------------------*
   SUBROUTINE S4TORT(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,STIFF,GVAL,WEIGHT)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
!*****
!* TORTIONAL STIFFNESS FOR S4 ELEMENT
!*****
      DIMENSION H(4),HR(4),HS(4),PHA12(4),PHA13(4),PHA23(4)
      DIMENSION AMAT(3,3),THE(3,3),STIFF(300)
      DIMENSION B11(24),B21(24),B31(24)
      DIMENSION B12(24),B22(24),B32(24)
      DIMENSION B13(24),B23(24),B33(24)
      DIMENSION WK(3,24),CC(24),BMAT(3,3),WORK1(3,3),WORK2(3,3)
!
      DO I=1,24
        B11(I)=0.0
        B21(I)=0.0
        B31(I)=0.0
        B12(I)=0.0
        B22(I)=0.0
        B32(I)=0.0
        B13(I)=0.0
        B23(I)=0.0
        B33(I)=0.0
        WK(1,I)=0.0
        WK(2,I)=0.0
        WK(3,I)=0.0
      enddo

      DO INOD=1,4
        B11(6*INOD-5)= HR(INOD)
        B11(6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B11(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B21(6*INOD-4)= HR(INOD)
        B21(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B21(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B31(6*INOD-3)= HR(INOD)
        B31(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B31(6*INOD-1)=-THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B12(6*INOD-5)= HS(INOD)
        B12(6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B12(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B22(6*INOD-4)= HS(INOD)
        B22(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B22(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B32(6*INOD-3)= HS(INOD)
        B32(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B32(6*INOD-1)=-THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B13(6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B13(6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B23(6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B23(6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B33(6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B33(6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
        WK(1,6*INOD-2)=H(INOD)
        WK(2,6*INOD-1)=H(INOD)
        WK(3,6*INOD  )=H(INOD)
      enddo

      DO I=1,24
        BMAT(1,1)=B11(I)
        BMAT(1,2)=B21(I)
        BMAT(1,3)=B31(I)
        BMAT(2,1)=B12(I)
        BMAT(2,2)=B22(I)
        BMAT(2,3)=B32(I)
        BMAT(3,1)=B13(I)
        BMAT(3,2)=B23(I)
        BMAT(3,3)=B33(I)

        DO II=1,3
          DO JJ=1,3
            WORK1(II,JJ)=0.0
            DO KK=1,3
              WORK1(II,JJ)=WORK1(II,JJ)+AMAT(II,KK)*BMAT(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            WORK2(II,JJ)=0.0
            DO KK=1,3
              WORK2(II,JJ)=WORK2(II,JJ)+WORK1(II,KK)*THE(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            BMAT(II,JJ)=0.0
            DO KK=1,3
              BMAT(II,JJ)=BMAT(II,JJ)+THE(KK,II)*WORK2(KK,JJ)
            enddo
          enddo
        enddo

        B21(I)=BMAT(1,2)
        B12(I)=BMAT(2,1)
        CC(I)=THE(1,3)*WK(1,I)+THE(2,3)*WK(2,I)+THE(3,3)*WK(3,I)-0.5*(B21(I)-B12(I))

      enddo
!
!****ADD TORTIONAL STIFFNESS 
      DO JJ=1,24
        DO II=1,JJ
          NPOIN=JJ*(JJ-1)/2+II
          VAL=GVAL*CC(II)*CC(JJ)
          STIFF(NPOIN)=STIFF(NPOIN)+VAL*WEIGHT
        enddo
      enddo
      RETURN

   end subroutine S4TORT
!----------------------------------------------------------------------*
   SUBROUTINE S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,RI,SI,TI,G1,G2,G3,B1,B2,B3)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,4),EN(3,4),PHA12(4),PHA13(4),PHA23(4)
      DIMENSION G1(3),G2(3),G3(3),B1(3,24),B2(3,24),B3(3,24)
      DIMENSION H(4),HR(4),HS(4)
! SET INTERPOLATION FUNCTION
      CALL S4PSET(RI,SI,H,HR,HS)
! SET G1,G2,G3
      CALL S4COV(COD,EN,H,HR,HS,TI,THICK,G1,G2,G3)
! CREATE B1,B2,B3
      DO INOD=1,4
! D{U}/DR
        B1(1,6*INOD-5)= HR(INOD)
        B1(1,6*INOD-4)= 0.0
        B1(1,6*INOD-3)= 0.0
        B1(1,6*INOD-2)= 0.0
        B1(1,6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B1(1,6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B1(2,6*INOD-5)= 0.0
        B1(2,6*INOD-4)= HR(INOD)
        B1(2,6*INOD-3)= 0.0
        B1(2,6*INOD-2)= -THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B1(2,6*INOD-1)= 0.0
        B1(2,6*INOD  )=  THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B1(3,6*INOD-5)= 0.0
        B1(3,6*INOD-4)= 0.0
        B1(3,6*INOD-3)= HR(INOD)
        B1(3,6*INOD-2)= -THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B1(3,6*INOD-1)= -THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B1(3,6*INOD  )= 0.0
!* D{U}/DS
        B2(1,6*INOD-5)= HS(INOD)
        B2(1,6*INOD-4)= 0.0
        B2(1,6*INOD-3)= 0.0
        B2(1,6*INOD-2)= 0.0
        B2(1,6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B2(1,6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B2(2,6*INOD-5)= 0.0
        B2(2,6*INOD-4)= HS(INOD)
        B2(2,6*INOD-3)= 0.0
        B2(2,6*INOD-2)= -THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B2(2,6*INOD-1)= 0.0
        B2(2,6*INOD  )=  THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B2(3,6*INOD-5)= 0.0
        B2(3,6*INOD-4)= 0.0
        B2(3,6*INOD-3)= HS(INOD)
        B2(3,6*INOD-2)= -THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B2(3,6*INOD-1)= -THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B2(3,6*INOD   )= 0.0
!* D{U}/DT
        B3(1,6*INOD-5)= 0.0
        B3(1,6*INOD-4)= 0.0
        B3(1,6*INOD-3)= 0.0
        B3(1,6*INOD-2)= 0.0
        B3(1,6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B3(1,6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B3(2,6*INOD-5)= 0.0
        B3(2,6*INOD-4)= 0.0
        B3(2,6*INOD-3)= 0.0
        B3(2,6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B3(2,6*INOD-1)= 0.0
        B3(2,6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B3(3,6*INOD-5)= 0.0
        B3(3,6*INOD-4)= 0.0
        B3(3,6*INOD-3)= 0.0
        B3(3,6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B3(3,6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
        B3(3,6*INOD  )= 0.0
      enddo
      RETURN

   end subroutine S4BMT1
!----------------------------------------------------------------------*
   SUBROUTINE STF_S3(XX,YY,ZZ,EE,PP,THICK,STIFF)
!----------------------------------------------------------------------*
! 
! CALCULATION OF 3-NODES ISOPARAMETRIC DEGENERATED SHELL ELEMENT
!
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(3),YY(3),ZZ(3),STIFF(171)
! LOCAL VARIABLES
      DIMENSION COD(3,3)
      DIMENSION INDX(4,15)
      DIMENSION G1(3),G2(3),G3(3),REF(3)
      DIMENSION G1A(3),G2A(3),G3A(3)
      DIMENSION G1B(3),G2B(3),G3B(3)
      DIMENSION G1C(3),G2C(3),G3C(3)
      DIMENSION G1D(3),G2D(3),G3D(3)
      DIMENSION PHA12(3),PHA13(3),PHA23(3)
      DIMENSION B1(3,18),B2(3,18),B3(3,18)
      DIMENSION B1A(3,18),B2A(3,18),B3A(3,18)
      DIMENSION B1B(3,18),B2B(3,18),B3B(3,18)
      DIMENSION B1C(3,18),B2C(3,18),B3C(3,18)
      DIMENSION B1D(3,18),B2D(3,18),B3D(3,18)
      DIMENSION B11(18),B22(18),B12(18),B13(18),B23(18)
      DIMENSION EN(3,3),THE(3,3),AMAT(3,3)
      DIMENSION H(3),HR(3),HS(3),XG(2),WGT(2)
      DIMENSION ELAS(9),XELAS(15)
!
! TUNING PARAMETER
      DATA SCF/1.2/
      DATA ILOC/1/
!
! SET CORD
      DO I=1,3
        COD(1,I)=XX(I)
        COD(2,I)=YY(I)
        COD(3,I)=ZZ(I)
      ENDDO
!
! SET INDEX
      CALL S3INDX(INDX)
!
! SET GAUUUSIAN INTEGRATION PARAMETER
      XG(1) =-0.5773502691896258
      XG(2) =-XG(1)
      WGT(1)=1.0
      WGT(2)=WGT(1)
!
! ZERO CLEAR ARRAY STIFF()
      DO I=1,171
        STIFF(I)=0.0
      ENDDO
!
! SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM
      CALL S3REFV(COD,REF,ILOC)
!
! CALCULATE [PHAI] AT EACH NODE
      CALL S3PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!
! LOOP OVER LAMINATES
!
! ELASTICITY MATRIX
      ELAS(1)=EE/(1.0-PP*PP)
      ELAS(2)=EE/(1.0-PP*PP)*PP
      ELAS(3)=EE/(1.0-PP*PP)
      ELAS(4)=0.0
      ELAS(5)=0.0
      ELAS(6)=EE/(1.0+PP)/2.0
      ELAS(7)=EE/(1.0+PP)/2.0/SCF
      ELAS(8)=0.0
      ELAS(9)=EE/(1.0+PP)/2.0/SCF
!
!   LOOP FOR GAUSS INTEGRATION POINT
!
      DO IG3=1,2
        TI=XG(IG3)
!
! METRICS AT A,B,C,D
        A0= 0.0
        AP= 1.0
        AM=-1.0
        CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    A0,AP, TI,G1A,G2A,G3A,B1A,B2A,B3A)
        CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    AM,A0, TI,G1B,G2B,G3B,B1B,B2B,B3B)
        CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    A0,AM, TI,G1C,G2C,G3C,B1C,B2C,B3C)
        CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,                            &
                    AP,A0, TI,G1D,G2D,G3D,B1D,B2D,B3D)
        DO IG2=1,2
          SI=XG(IG2)
          DO IG1=1,2
            RI=XG(IG1)
            CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,RI,SI,TI,G1,G2,G3,B1,B2,B3)
!
! JACOBI MATRIX 
            CALL S3JAC(G1,G2,G3,DET,AMAT)
! SET INTEGRATION WEIGHT
            WEIGHT=WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
! SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
            CALL S3THE(G3,REF,THE)
! SET STRAIN IN NATURAL COORDINATE
!
            DO I=1,18
              B11(I)=G1(1)*B1(1,I)+G1(2)*B1(2,I)+G1(3)*B1(3,I)
              B22(I)=G2(1)*B2(1,I)+G2(2)*B2(2,I)+G2(3)*B2(3,I)
              B12(I)=(G2(1)*B1(1,I)+G2(2)*B1(2,I)+G2(3)*B1(3,I)                &
                    +G1(1)*B2(1,I)+G1(2)*B2(2,I)+G1(3)*B2(3,I))
              B13(I)=0.5*(1.0+SI)*                                             &
                    (G3A(1)*B1A(1,I)+G3A(2)*B1A(2,I)+G3A(3)*B1A(3,I)           &
                    +G1A(1)*B3A(1,I)+G1A(2)*B3A(2,I)+G1A(3)*B3A(3,I))          &
                    +0.5*(1.0-SI)*                                             &
                    (G3C(1)*B1C(1,I)+G3C(2)*B1C(2,I)+G3C(3)*B1C(3,I)           &
                    +G1C(1)*B3C(1,I)+G1C(2)*B3C(2,I)+G1C(3)*B3C(3,I))
              B23(I)=0.5*(1.0+RI)*                                             &
                     (G3D(1)*B2D(1,I)+G3D(2)*B2D(2,I)+G3D(3)*B2D(3,I)          &
                     +G2D(1)*B3D(1,I)+G2D(2)*B3D(2,I)+G2D(3)*B3D(3,I))         &
                     +0.5*(1.0-RI)*                                            &
                     (G3B(1)*B2B(1,I)+G3B(2)*B2B(2,I)+G3B(3)*B2B(3,I)          &
                     +G2B(1)*B3B(1,I)+G2B(2)*B3B(2,I)+G2B(3)*B3B(3,I))
            enddo
!
! CONSTITUTIVE LAW
            CALL S3CONS(THE,G1,G2,G3,DET,ELAS,INDX,XELAS)
!
! ELEMENT STIFFNESS
!
            DO JJ=1,18
              DO II=1,JJ
                NPOIN=JJ*(JJ-1)/2+II
                VAL=                                                           &
                  +XELAS( 1)*B11(II)*B11(JJ)                                   &
                  +XELAS( 2)*(B11(II)*B22(JJ)+B22(II)*B11(JJ))                 &
                  +XELAS( 3)*B22(II)*B22(JJ)                                   &
                  +XELAS( 4)*(B11(II)*B12(JJ)+B12(II)*B11(JJ))                 &
                  +XELAS( 5)*(B22(II)*B12(JJ)+B12(II)*B22(JJ))                 &
                  +XELAS( 6)* B12(II)*B12(JJ)                                  &
                  +XELAS( 7)*(B11(II)*B13(JJ)+B13(II)*B11(JJ))                 &
                  +XELAS( 8)*(B22(II)*B13(JJ)+B13(II)*B22(JJ))                 &
                  +XELAS( 9)*(B12(II)*B13(JJ)+B13(II)*B12(JJ))                 &
                  +XELAS(10)* B13(II)*B13(JJ)                                  &
                  +XELAS(11)*(B11(II)*B23(JJ)+B23(II)*B11(JJ))                 &
                  +XELAS(12)*(B22(II)*B23(JJ)+B23(II)*B22(JJ))                 &
                  +XELAS(13)*(B12(II)*B23(JJ)+B23(II)*B12(JJ))                 &
                  +XELAS(14)*(B13(II)*B23(JJ)+B23(II)*B13(JJ))                 &
                  +XELAS(15)* B23(II)*B23(JJ)
                  STIFF(NPOIN)=STIFF(NPOIN)+VAL*WEIGHT
              enddo
            enddo
!
!  ADD TORSIONAL STIFFNESS 
            CALL S3PSET(RI,SI,H,HR,HS)
            CALL S3TORT(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,STIFF,ELAS(6),WEIGHT)
          enddo
        enddo
      enddo
      RETURN

   end subroutine STF_S3
!----------------------------------------------------------------------*
   SUBROUTINE S3CONS(THE,G1,G2,G3,DET,ELAS,INDX,XELAS)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION THE(3,3),G1(3),G2(3),G3(3),ELAS(9),INDX(4,15),XELAS(15)
      DIMENSION GG1(3),GG2(3),GG3(3),A(3,3)
! CONTRAVARIANT BASE VECTOR
      GG1(1)=(G2(2)*G3(3)-G2(3)*G3(2))/DET
      GG1(2)=(G2(3)*G3(1)-G2(1)*G3(3))/DET
      GG1(3)=(G2(1)*G3(2)-G2(2)*G3(1))/DET
      GG2(1)=(G3(2)*G1(3)-G3(3)*G1(2))/DET
      GG2(2)=(G3(3)*G1(1)-G3(1)*G1(3))/DET
      GG2(3)=(G3(1)*G1(2)-G3(2)*G1(1))/DET
      GG3(1)=(G1(2)*G2(3)-G1(3)*G2(2))/DET
      GG3(2)=(G1(3)*G2(1)-G1(1)*G2(3))/DET
      GG3(3)=(G1(1)*G2(2)-G1(2)*G2(1))/DET
      DO I=1,3
        A(1,I)=GG1(1)*THE(1,I)+GG1(2)*THE(2,I)+GG1(3)*THE(3,I)
        A(2,I)=GG2(1)*THE(1,I)+GG2(2)*THE(2,I)+GG2(3)*THE(3,I)
        A(3,I)=GG3(1)*THE(1,I)+GG3(2)*THE(2,I)+GG3(3)*THE(3,I)
      enddo

      DO I=1,15
        I1=INDX(1,I)
        I2=INDX(2,I)
        I3=INDX(3,I)
        I4=INDX(4,I)
        XELAS(I)=A(I1,1)*A(I2,1)*A(I3,1)*A(I4,1)*ELAS(1)                       &
                +A(I1,1)*A(I2,1)*A(I3,2)*A(I4,2)*ELAS(2)                       &
                +A(I1,2)*A(I2,2)*A(I3,1)*A(I4,1)*ELAS(2)                       &
                +A(I1,2)*A(I2,2)*A(I3,2)*A(I4,2)*ELAS(3)                       &
                +A(I1,1)*A(I2,1)*A(I3,1)*A(I4,2)*ELAS(4)                       &
                +A(I1,1)*A(I2,2)*A(I3,1)*A(I4,1)*ELAS(4)                       &
                +A(I1,2)*A(I2,1)*A(I3,1)*A(I4,1)*ELAS(4)                       &
                +A(I1,1)*A(I2,1)*A(I3,2)*A(I4,1)*ELAS(4)                       &
                +A(I1,2)*A(I2,2)*A(I3,1)*A(I4,2)*ELAS(5)                       &
                +A(I1,2)*A(I2,2)*A(I3,2)*A(I4,1)*ELAS(5)                       &
                +A(I1,1)*A(I2,2)*A(I3,2)*A(I4,2)*ELAS(5)                       &
                +A(I1,2)*A(I2,1)*A(I3,2)*A(I4,2)*ELAS(5)                       &
                +A(I1,1)*A(I2,2)*A(I3,1)*A(I4,2)*ELAS(6)                       &
                +A(I1,2)*A(I2,1)*A(I3,2)*A(I4,1)*ELAS(6)                       &
                +A(I1,1)*A(I2,2)*A(I3,2)*A(I4,1)*ELAS(6)                       &
                +A(I1,2)*A(I2,1)*A(I3,1)*A(I4,2)*ELAS(6)
        XELAS(I)=XELAS(I)                                                      &
                +A(I1,1)*A(I2,3)*A(I3,1)*A(I4,3)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,3)*A(I4,1)*ELAS(7)                       &
                +A(I1,1)*A(I2,3)*A(I3,3)*A(I4,1)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,1)*A(I4,3)*ELAS(7)                       &
                +A(I1,3)*A(I2,1)*A(I3,3)*A(I4,2)*ELAS(8)                       &
                +A(I1,3)*A(I2,1)*A(I3,2)*A(I4,3)*ELAS(8)                       &
                +A(I1,1)*A(I2,3)*A(I3,3)*A(I4,2)*ELAS(8)                       &
                +A(I1,1)*A(I2,3)*A(I3,2)*A(I4,3)*ELAS(8)                       &
                +A(I1,3)*A(I2,2)*A(I3,3)*A(I4,1)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,3)*A(I4,1)*ELAS(8)                       &
                +A(I1,3)*A(I2,2)*A(I3,1)*A(I4,3)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,1)*A(I4,3)*ELAS(8)                       &
                +A(I1,2)*A(I2,3)*A(I3,2)*A(I4,3)*ELAS(9)                       &
                +A(I1,3)*A(I2,2)*A(I3,3)*A(I4,2)*ELAS(9)                       &
                +A(I1,2)*A(I2,3)*A(I3,3)*A(I4,2)*ELAS(9)                       &
                +A(I1,3)*A(I2,2)*A(I3,2)*A(I4,3)*ELAS(9)
      enddo
      RETURN

   end subroutine S3CONS
!----------------------------------------------------------------------*
   SUBROUTINE S3COV(COD,EN,H,HR,HS,TI,THICK,G1,G2,G3)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,3),EN(3,3),H(3),HR(3),HS(3),G1(3),G2(3),G3(3)
      DO I=1,3
        G1(I)=0.0
        G2(I)=0.0
        G3(I)=0.0
        DO INOD=1,3
          VAR=COD(I,INOD)+THICK*0.5*TI*EN(I,INOD)
          G1(I)=G1(I)+HR(INOD)*VAR
          G2(I)=G2(I)+HS(INOD)*VAR
          G3(I)=G3(I)+THICK*0.5*H(INOD)*EN(I,INOD)
        enddo
      enddo
      RETURN

   end subroutine S3COV
!----------------------------------------------------------------------*
      SUBROUTINE S3INDX(INDX)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION INDX(4,15)
      INDX(1,1)=1
      INDX(2,1)=1
      INDX(3,1)=1
      INDX(4,1)=1
      INDX(1,2)=1
      INDX(2,2)=1
      INDX(3,2)=2
      INDX(4,2)=2
      INDX(1,3)=2
      INDX(2,3)=2
      INDX(3,3)=2
      INDX(4,3)=2
      INDX(1,4)=1
      INDX(2,4)=1
      INDX(3,4)=1
      INDX(4,4)=2
      INDX(1,5)=2
      INDX(2,5)=2
      INDX(3,5)=1
      INDX(4,5)=2
      INDX(1,6)=1
      INDX(2,6)=2
      INDX(3,6)=1
      INDX(4,6)=2
      INDX(1,7)=1
      INDX(2,7)=1
      INDX(3,7)=1
      INDX(4,7)=3
      INDX(1,8)=2
      INDX(2,8)=2
      INDX(3,8)=1
      INDX(4,8)=3
      INDX(1,9)=1
      INDX(2,9)=2
      INDX(3,9)=1
      INDX(4,9)=3
      INDX(1,10)=1
      INDX(2,10)=3
      INDX(3,10)=1
      INDX(4,10)=3
      INDX(1,11)=1
      INDX(2,11)=1
      INDX(3,11)=2
      INDX(4,11)=3
      INDX(1,12)=2
      INDX(2,12)=2
      INDX(3,12)=2
      INDX(4,12)=3
      INDX(1,13)=1
      INDX(2,13)=2
      INDX(3,13)=2
      INDX(4,13)=3
      INDX(1,14)=1
      INDX(2,14)=3
      INDX(3,14)=2
      INDX(4,14)=3
      INDX(1,15)=2
      INDX(2,15)=3
      INDX(3,15)=2
      INDX(4,15)=3
      RETURN

   end subroutine S3INDX
!----------------------------------------------------------------------*
   SUBROUTINE S3JAC(G1,G2,G3,DET,AMAT)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION G1(3),G2(3),G3(3),AMAT(3,3)
      XJ11=G1(1)
      XJ12=G1(2)  
      XJ13=G1(3)  
      XJ21=G2(1)
      XJ22=G2(2)  
      XJ23=G2(3)  
      XJ31=G3(1)
      XJ32=G3(2)  
      XJ33=G3(3)  
! DETERMINANT OF JACOBIAN
      DET=XJ11*XJ22*XJ33                                                       &
         +XJ12*XJ23*XJ31                                                       &
         +XJ13*XJ21*XJ32                                                       &
         -XJ13*XJ22*XJ31                                                       &
         -XJ12*XJ21*XJ33                                                       &
         -XJ11*XJ23*XJ32
! INVERSION OF JACOBIAN
      DUM=1.0/DET
      AMAT(1,1)=DUM*( XJ22*XJ33-XJ23*XJ32)
      AMAT(2,1)=DUM*(-XJ21*XJ33+XJ23*XJ31)
      AMAT(3,1)=DUM*( XJ21*XJ32-XJ22*XJ31)
      AMAT(1,2)=DUM*(-XJ12*XJ33+XJ13*XJ32)
      AMAT(2,2)=DUM*( XJ11*XJ33-XJ13*XJ31)
      AMAT(3,2)=DUM*(-XJ11*XJ32+XJ12*XJ31)
      AMAT(1,3)=DUM*( XJ12*XJ23-XJ13*XJ22)
      AMAT(2,3)=DUM*(-XJ11*XJ23+XJ13*XJ21)
      AMAT(3,3)=DUM*( XJ11*XJ22-XJ12*XJ21)
      RETURN

   end subroutine S3JAC
!----------------------------------------------------------------------*
   SUBROUTINE S3PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,3),REF(3),PHA12(3),PHA13(3),PHA23(3),EN(3,3)
      DIMENSION G1(3),G2(3),G3(3),E1(3),E2(3),E3(3)
! IN-PLANE VECTOR
      G1(1)=-0.5*COD(1,1)+0.5*COD(1,2)
      G1(2)=-0.5*COD(2,1)+0.5*COD(2,2)
      G1(3)=-0.5*COD(3,1)+0.5*COD(3,2)
      G2(1)=-0.5*COD(1,1)+0.5*COD(1,3)
      G2(2)=-0.5*COD(2,1)+0.5*COD(2,3)
      G2(3)=-0.5*COD(3,1)+0.5*COD(3,3)
! G3()=G1() X G2()
      G3(1)=G1(2)*G2(3)-G1(3)*G2(2)
      G3(2)=G1(3)*G2(1)-G1(1)*G2(3)
      G3(3)=G1(1)*G2(2)-G1(2)*G2(1)
!  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
      XSUM=DSQRT(G3(1)**2+G3(2)**2+G3(3)**2)
      E3(1)=G3(1)/XSUM
      E3(2)=G3(2)/XSUM
      E3(3)=G3(3)/XSUM
      E2(1)=-REF(2)*E3(3)+REF(3)*E3(2)
      E2(2)=-REF(3)*E3(1)+REF(1)*E3(3)
      E2(3)=-REF(1)*E3(2)+REF(2)*E3(1)
      XSUM=DSQRT(E2(1)**2+E2(2)**2+E2(3)**2)
      E2(1)=E2(1)/XSUM
      E2(2)=E2(2)/XSUM
      E2(3)=E2(3)/XSUM
      E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
      E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
      E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
      XSUM=DSQRT(E1(1)**2+E1(2)**2+E1(3)**2)
      E1(1)=E1(1)/XSUM
      E1(2)=E1(2)/XSUM
      E1(3)=E1(3)/XSUM
!** EN & E1,E2,E3 IS CONSTANT FOR 3-NODE ELEMENT
      DO I=1,3
        EN(1,I)=E3(1)
        EN(2,I)=E3(2)
        EN(3,I)=E3(3)
        PHA12(I)=-E2(1)*E1(2)+E1(1)*E2(2)
        PHA13(I)=-E2(1)*E1(3)+E1(1)*E2(3)
        PHA23(I)=-E2(2)*E1(3)+E1(2)*E2(3)
      enddo
      RETURN

   end subroutine S3PHAI
!----------------------------------------------------------------------*
      SUBROUTINE S3PSET(RI,SI,H,HR,HS)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION H(3),HR(3),HS(3)
      RP=1.0+RI
      SP=1.0+SI
      RM=1.0-RI
      SM=1.0-SI
      H(1)=0.25*RM*SM
      H(2)=0.25*RP*SM
      H(3)=0.5*SP
      HR(1)=-.25*SM
      HR(2)= .25*SM
      HR(3)= 0.0
      HS(1)=-.25*RM
      HS(2)=-.25*RP
      HS(3)= 0.5
      RETURN

   end subroutine S3PSET
!----------------------------------------------------------------------*
   SUBROUTINE S3REFV(COD,REF,ILOC)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,3),REF(3)
      IF(ILOC.EQ.0) THEN
        REF(1)=1.0  
        REF(2)=0.0 
        REF(3)=0.0  
      ELSE 
        REF(1)=(COD(1,2)-COD(1,1))
        REF(2)=(COD(2,2)-COD(2,1))
        REF(3)=(COD(3,2)-COD(3,1))
      END IF
      RETURN

   end subroutine S3REFV
!----------------------------------------------------------------------*
      SUBROUTINE S3THE(G3,REF,THE)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION G3(3),REF(3),THE(3,3)
      DIMENSION E1(3),E2(3),E3(3)
      XSUM=DSQRT(G3(1)**2+G3(2)**2+G3(3)**2)
      E3(1)=G3(1)/XSUM
      E3(2)=G3(2)/XSUM
      E3(3)=G3(3)/XSUM
      E2(1)=-REF(2)*E3(3)+REF(3)*E3(2)
      E2(2)=-REF(3)*E3(1)+REF(1)*E3(3)
      E2(3)=-REF(1)*E3(2)+REF(2)*E3(1)
      E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
      E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
      E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
      XSUM=DSQRT(E2(1)**2+E2(2)**2+E2(3)**2)
      IF (XSUM.GT.1.E-15) THEN
        E2(1)=E2(1)/XSUM
        E2(2)=E2(2)/XSUM
        E2(3)=E2(3)/XSUM
        E1(1)=-E3(2)*E2(3)+E3(3)*E2(2)
        E1(2)=-E3(3)*E2(1)+E3(1)*E2(3)
        E1(3)=-E3(1)*E2(2)+E3(2)*E2(1)
        XSUM=DSQRT(E1(1)**2+E1(2)**2+E1(3)**2)
        E1(1)=E1(1)/XSUM
        E1(2)=E1(2)/XSUM
        E1(3)=E1(3)/XSUM
      ELSE
        E1(1)=  0.0
        E1(2)=  0.0
        E1(3)=- 1.0
        E2(1) = 0.0
        E2(2) = 1.0
        E2(3) = 0.0
      END IF
      THE(1,1)=E1(1)
      THE(1,2)=E2(1)
      THE(1,3)=E3(1)
      THE(2,1)=E1(2)
      THE(2,2)=E2(2)
      THE(2,3)=E3(2)
      THE(3,1)=E1(3)
      THE(3,2)=E2(3)
      THE(3,3)=E3(3)
      RETURN

   end subroutine S3THE
!----------------------------------------------------------------------*
   SUBROUTINE S3TORT(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,STIFF,GVAL,WEIGHT)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
!*****
!* TORTIONAL STIFFNESS FOR S3 ELEMENT
!*****
      DIMENSION H(3),HR(3),HS(3),PHA12(3),PHA13(3),PHA23(3)          
      DIMENSION AMAT(3,3),THE(3,3),STIFF(171)        
      DIMENSION B11(18),B21(18),B31(18) 
      DIMENSION B12(18),B22(18),B32(18)
      DIMENSION B13(18),B23(18),B33(18)
      DIMENSION WK(3,18),CC(18),BMAT(3,3),WORK1(3,3),WORK2(3,3)
!
      DO I=1,18
        B11(I)=0.0
        B21(I)=0.0
        B31(I)=0.0
        B12(I)=0.0
        B22(I)=0.0
        B32(I)=0.0
        B13(I)=0.0
        B23(I)=0.0
        B33(I)=0.0
        WK(1,I)=0.0
        WK(2,I)=0.0
        WK(3,I)=0.0
      enddo

      DO INOD=1,3
        B11(6*INOD-5)= HR(INOD)
        B11(6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B11(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B21(6*INOD-4)= HR(INOD)
        B21(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B21(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B31(6*INOD-3)= HR(INOD)
        B31(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B31(6*INOD-1)=-THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B12(6*INOD-5)= HS(INOD)
        B12(6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B12(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B22(6*INOD-4)= HS(INOD)
        B22(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B22(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B32(6*INOD-3)= HS(INOD)
        B32(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B32(6*INOD-1)=-THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B13(6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B13(6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B23(6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B23(6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B33(6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B33(6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
        WK(1,6*INOD-2)=H(INOD)
        WK(2,6*INOD-1)=H(INOD)
        WK(3,6*INOD  )=H(INOD)
      enddo

      DO I=1,18
        BMAT(1,1)=B11(I)
        BMAT(1,2)=B21(I)
        BMAT(1,3)=B31(I)
        BMAT(2,1)=B12(I)
        BMAT(2,2)=B22(I)
        BMAT(2,3)=B32(I)
        BMAT(3,1)=B13(I)
        BMAT(3,2)=B23(I)
        BMAT(3,3)=B33(I)
!
        DO II=1,3
          DO  JJ=1,3
            WORK1(II,JJ)=0.0
            DO  KK=1,3
              WORK1(II,JJ)=WORK1(II,JJ)+AMAT(II,KK)*BMAT(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            WORK2(II,JJ)=0.0
            DO KK=1,3
              WORK2(II,JJ)=WORK2(II,JJ)+WORK1(II,KK)*THE(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            BMAT(II,JJ)=0.0
            DO KK=1,3
              BMAT(II,JJ)=BMAT(II,JJ)+THE(KK,II)*WORK2(KK,JJ)
            enddo
          enddo
        enddo

        B21(I)=BMAT(1,2)
        B12(I)=BMAT(2,1)
        CC(I)=THE(1,3)*WK(1,I)+THE(2,3)*WK(2,I)+THE(3,3)*WK(3,I)-0.5*(B21(I)-B12(I))

      enddo
!****ADD TORTIONAL STIFFNESS 
      DO JJ=1,18
        DO II=1,JJ
          NPOIN=JJ*(JJ-1)/2+II
          VAL=GVAL*CC(II)*CC(JJ)
          STIFF(NPOIN)=STIFF(NPOIN)+VAL*WEIGHT
        enddo
      enddo
      RETURN

   end subroutine S3TORT
!----------------------------------------------------------------------*
   SUBROUTINE S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,RI,SI,TI,G1,G2,G3,B1,B2,B3)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION COD(3,3),EN(3,3),PHA12(3),PHA13(3),PHA23(3)
      DIMENSION G1(3),G2(3),G3(3),B1(3,18),B2(3,18),B3(3,18)
      DIMENSION H(3),HR(3),HS(3)
! SET INTERPOLATION FUNCTION
      CALL S3PSET(RI,SI,H,HR,HS)
! SET G1,G2,G3
      CALL S3COV(COD,EN,H,HR,HS,TI,THICK,G1,G2,G3)
! CREATE B1,B2,B3
      DO INOD=1,3
! D{U}/DR
        B1(1,6*INOD-5)= HR(INOD)
        B1(1,6*INOD-4)= 0.0
        B1(1,6*INOD-3)= 0.0
        B1(1,6*INOD-2)= 0.0
        B1(1,6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B1(1,6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B1(2,6*INOD-5)= 0.0
        B1(2,6*INOD-4)= HR(INOD)
        B1(2,6*INOD-3)= 0.0
        B1(2,6*INOD-2)= -THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B1(2,6*INOD-1)= 0.0
        B1(2,6*INOD  )=  THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B1(3,6*INOD-5)= 0.0
        B1(3,6*INOD-4)= 0.0
        B1(3,6*INOD-3)= HR(INOD)
        B1(3,6*INOD-2)= -THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B1(3,6*INOD-1)= -THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B1(3,6*INOD  )= 0.0
!* D{U}/DS
        B2(1,6*INOD-5)= HS(INOD)
        B2(1,6*INOD-4)= 0.0
        B2(1,6*INOD-3)= 0.0
        B2(1,6*INOD-2)= 0.0
        B2(1,6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B2(1,6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B2(2,6*INOD-5)= 0.0
        B2(2,6*INOD-4)= HS(INOD)
        B2(2,6*INOD-3)= 0.0
        B2(2,6*INOD-2)= -THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B2(2,6*INOD-1)= 0.0
        B2(2,6*INOD  )=  THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B2(3,6*INOD-5)= 0.0
        B2(3,6*INOD-4)= 0.0
        B2(3,6*INOD-3)= HS(INOD)
        B2(3,6*INOD-2)= -THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B2(3,6*INOD-1)= -THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B2(3,6*INOD   )= 0.0
!* D{U}/DT
        B3(1,6*INOD-5)= 0.0
        B3(1,6*INOD-4)= 0.0
        B3(1,6*INOD-3)= 0.0
        B3(1,6*INOD-2)= 0.0
        B3(1,6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B3(1,6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B3(2,6*INOD-5)= 0.0
        B3(2,6*INOD-4)= 0.0
        B3(2,6*INOD-3)= 0.0
        B3(2,6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B3(2,6*INOD-1)= 0.0
        B3(2,6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B3(3,6*INOD-5)= 0.0
        B3(3,6*INOD-4)= 0.0
        B3(3,6*INOD-3)= 0.0
        B3(3,6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B3(3,6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
        B3(3,6*INOD  )= 0.0
      enddo
      RETURN

   end subroutine S3BMT1
!----------------------------------------------------------------------*
   SUBROUTINE DL_S4(XX,YY,ZZ,RHO,THICK,LTYPE,PARAMS,VECT,NSIZE)
!----------------------------------------------------------------------*
!**
!**  SET S4 DLOAD
!** 
      use hecmw
      use m_utilities
!
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(*),YY(*),ZZ(*),PARAMS(*),VECT(*)
!   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
!   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
!   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
!   CRAV LTYPE=4  :GRAVITY FORCE
!   CENT LTYPE=5  :CENTRIFUGAL LOAD
!   P    LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR SHELL SURFACE
! LOCAL VARIABLES
      DIMENSION COD(3,4),G1(3),G2(3),G3(3),REF(3)
      DIMENSION PHA12(4),PHA13(4),PHA23(4)
      DIMENSION EN(3,4),AMAT(3,3)
      DIMENSION H(4),HR(4),HS(4),XG(2),WGT(2),UVW(3,24)
      DATA ILOC/1/
!*************************
!  GAUSS INTEGRATION POINT
!*************************
      DATA XG/-0.5773502691896,0.5773502691896/
      DATA WGT/1.0,1.0/
      INTEGER(kind=kint) NSIZE,NN,NDOF
!
! INITIALIZE
!
      NN = 4     
      NDOF = 6
      NSIZE = NN*NDOF
      DO I = 1, NN
        COD(1,I) = XX(I)
        COD(2,I) = YY(I)
        COD(3,I) = ZZ(I)
      ENDDO
      VAL = PARAMS(1)
      AX  = PARAMS(2)
      AY  = PARAMS(3)
      AZ  = PARAMS(4)
      RX  = PARAMS(5)
      RY  = PARAMS(6)
      RZ  = PARAMS(7)
!
! LOCAL LOAD VECTOR
!
      DO I=1,NSIZE
        VECT(I)=0.0
      ENDDO
              
! SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM
      CALL S4REFV(COD,REF,ILOC)
! CALCULATE [PHAI] AT EACH NODE
      CALL S4PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!
! SELCTION OF LOAD TYPE
!
      IVOL=0
      ISUF=0
      IF( LTYPE.LT.10 ) THEN 
        IVOL=1
      ELSE IF( LTYPE.GE.10 ) THEN 
        ISUF=1
      ENDIF
!** SURFACE LOAD
      IF( ISUF.EQ.1 ) THEN
! INTEGRATION OVER VOLUME 
        DO IG2=1,2
          SI=XG(IG2)
          DO IG1=1,2
            RI=XG(IG1)
! INTERPOLATION  FUNCTIONS
            CALL S4PSET(RI,SI,H,HR,HS)
! COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
            zero=0.0
            CALL S4COV(COD,EN,H,HR,HS,zero,THICK,G1,G2,G3)
! JACOBI MATRIX
            CALL S4JAC(G1,G2,G3,DET,AMAT)
            XSUM=DSQRT(G3(1)**2+G3(2)**2+G3(3)**2)
            G3(1)=G3(1)/XSUM
            G3(2)=G3(2)/XSUM
            G3(3)=G3(3)/XSUM
! SET INTEGRATION WEIGHT
            WEIGHT=WGT(IG1)*WGT(IG2)*DET*2.0
! CREATE [UVW]:3*24 MATRIX
            UVW(:,:)=0.d0
            DO INOD=1,NN
              UVW(1,6*INOD-5)=H(INOD)
              UVW(2,6*INOD-4)=H(INOD)
              UVW(3,6*INOD-3)=H(INOD)
            ENDDO
            DO I=1,NSIZE
              VECT(I)=VECT(I)+(UVW(1,I)*G3(1)+UVW(2,I)*G3(2)+UVW(3,I)*G3(3))*VAL*WEIGHT 
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!** VOLUME LOAD
      IF( IVOL.EQ.1 ) THEN
! INTEGRATION OVER VOLUME
        DO IG3=1,2
          TI=XG(IG3)
          DO IG2=1,2
            SI=XG(IG2)
            DO IG1=1,2
              RI=XG(IG1)
! INTERPOLATION  FUNCTIONS
              CALL S4PSET(RI,SI,H,HR,HS)
! COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
              CALL S4COV(COD,EN,H,HR,HS,TI,THICK,G1,G2,G3)
! JACOBI MATRIX
              CALL S4JAC(G1,G2,G3,DET,AMAT)
! SET INTEGRATION WEIGHT
              WEIGHT=WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
! CREATE [UVW]:3*24 MATRIX
              UVW(:,:)=0.d0
              DO INOD=1,NN
                UVW(1,6*INOD-5)=H(INOD)
                UVW(1,6*INOD-1)=0.5*THICK*TI*H(INOD)*PHA12(INOD)
                UVW(1,6*INOD  )=0.5*THICK*TI*H(INOD)*PHA13(INOD)
                UVW(2,6*INOD-4)=H(INOD)
                UVW(2,6*INOD-2)=-0.5*THICK*TI*H(INOD)*PHA12(INOD)
                UVW(2,6*INOD  )= 0.5*THICK*TI*H(INOD)*PHA23(INOD)
                UVW(3,6*INOD-3)=H(INOD)
                UVW(3,6*INOD-2)=-0.5*THICK*TI*H(INOD)*PHA13(INOD)
                UVW(3,6*INOD-1)=-0.5*THICK*TI*H(INOD)*PHA23(INOD)
              ENDDO
              IF( LTYPE.EQ.1) THEN
                DO I=1,NSIZE
                  VECT(I)=VECT(I)+VAL*UVW(1,I)*WEIGHT
                ENDDO
              ELSE IF( LTYPE.EQ.2 ) THEN
                DO I=1,NSIZE
                  VECT(I)=VECT(I)+VAL*UVW(2,I)*WEIGHT
                ENDDO
              ELSE IF( LTYPE.EQ.3 ) THEN
                DO I=1,NSIZE
                  VECT(I)=VECT(I)+VAL*UVW(3,I)*WEIGHT
                ENDDO
              ELSE IF( LTYPE.EQ.4 ) THEN
                DO I=1,NSIZE
                  VECT(I)=VECT(I)+VAL*RHO*AX*UVW(1,I)*WEIGHT
                  VECT(I)=VECT(I)+VAL*RHO*AY*UVW(2,I)*WEIGHT
                  VECT(I)=VECT(I)+VAL*RHO*AZ*UVW(3,I)*WEIGHT
                ENDDO
              ELSE IF( LTYPE.EQ.5 ) THEN
                XXX = 0.0
                YYY = 0.0
                ZZZ = 0.0
                DO INOD = 1, NN
                  XXX = XXX + H(INOD)*XX(INOD)
                  YYY = YYY + H(INOD)*YY(INOD)
                  ZZZ = ZZZ + H(INOD)*ZZ(INOD)
                ENDDO
                HX=AX+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RX
                 
                HY=AY+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RY
                 
                HZ=AZ+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RZ
                 
                PHX=XXX-HX
                PHY=YYY-HY
                PHZ=ZZZ-HZ
!                COEFX=PHX*RHO*VAL**2
!                COEFY=PHY*RHO*VAL**2
!                COEFZ=PHZ*RHO*VAL**2
                COEFX=PHX*VAL*RHO*VAL
                COEFY=PHY*VAL*RHO*VAL
                COEFZ=PHZ*VAL*RHO*VAL

                DO I=1,NSIZE
                  VECT(I)=VECT(I)+(UVW(1,I)*COEFX+UVW(2,I)*COEFY+UVW(3,I)*COEFZ)*WEIGHT
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      RETURN
   end subroutine DL_S4
!----------------------------------------------------------------------*
   SUBROUTINE DL_S3(XX,YY,ZZ,RHO,THICK,LTYPE,PARAMS,VECT,NSIZE)
!----------------------------------------------------------------------*
!**
!**  SET S3 DLOAD
!** 
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(*),YY(*),ZZ(*),PARAMS(*),VECT(*)
!   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
!   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
!   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
!   GRAV LTYPE=4  :GRAVITY FORCE
!   CENT LTYPE=5  :CENTRIFUGAL LOAD
!   P    LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR SHELL SURFACE
! LOCAL VARIABLES
      DIMENSION XG(2),WGT(2),H(3)
      DIMENSION PLX(3),PLY(3),PLZ(3)
!*************************
!  GAUSS INTEGRATION POINT
!*************************
      DATA XG/-0.5773502691896,0.5773502691896/
      DATA WGT/1.0,1.0/
      INTEGER(kind=kint) NSIZE,NN,NDOF
!
! SELCTION OF LOAD TYPE
!
      NN=3
      NDOF=6
      NSIZE=NN*NDOF

      IVOL=0
      ISUF=0
      IF( LTYPE.LT.10 ) THEN 
        IVOL=1
      ELSE IF( LTYPE.GE.10 ) THEN 
        ISUF=1
      ENDIF
!
! LOCAL LOAD VECTOR
!
      DO I=1,NSIZE
        VECT(I)=0.0
      ENDDO
      VAL = PARAMS(1)
      AX  = PARAMS(2)
      AY  = PARAMS(3)
      AZ  = PARAMS(4)
      RX  = PARAMS(5)
      RY  = PARAMS(6)
      RZ  = PARAMS(7)
!** SURFACE LOAD
      IF( ISUF.EQ.1 ) THEN
!** FACE 1-2-3
        V1X=XX(2)-XX(1)
        V1Y=YY(2)-YY(1)
        V1Z=ZZ(2)-ZZ(1)
        V2X=XX(3)-XX(1)
        V2Y=YY(3)-YY(1)
        V2Z=ZZ(3)-ZZ(1)
        V3X= V1Y*V2Z-V1Z*V2Y
        V3Y=-V1X*V2Z+V1Z*V2X
        V3Z= V1X*V2Y-V1Y*V2X
        VECT(1)=VAL*V3X/6.0
        VECT(2)=VAL*V3Y/6.0
        VECT(3)=VAL*V3Z/6.0
        VECT(7)=VAL*V3X/6.0
        VECT(8)=VAL*V3Y/6.0
        VECT(9)=VAL*V3Z/6.0
        VECT(13)=VAL*V3X/6.0
        VECT(14)=VAL*V3Y/6.0
        VECT(15)=VAL*V3Z/6.0
      ENDIF
!** VOLUME LOAD
      IF( IVOL.EQ.1 ) THEN
        DO I=1,NN
          PLX(I)=0.0
          PLY(I)=0.0
          PLZ(I)=0.0
        ENDDO
!** AREA OF SURFACE
        V1X=XX(2)-XX(1)
        V1Y=YY(2)-YY(1)
        V1Z=ZZ(2)-ZZ(1)
        V2X=XX(3)-XX(1)
        V2Y=YY(3)-YY(1)
        V2Z=ZZ(3)-ZZ(1)
        V3X= V1Y*V2Z-V1Z*V2Y
        V3Y=-V1X*V2Z+V1Z*V2X
        V3Z= V1X*V2Y-V1Y*V2X
        AREA=0.5*SQRT(V3X**2+V3Y**2+V3Z**2)
        DO L2=1,2
          XL2=XG(L2)
          X2 =0.5*(XL2+1.0)
          DO L1=1,2
            XL1=XG(L1)
            X1=0.5*(XL1+1.0)
!  INTERPOLATION FUNCTION
            H(1)=1.0-X1
            H(2)=X1*(1.0-X2)
            H(3)=X1*X2
!  WEIGT VALUE AT GAUSSIAN POINT
            WG=0.5*X1*WGT(L1)*WGT(L2)*AREA*THICK
            COEFX=1.0
            COEFY=1.0
            COEFZ=1.0
            IF( LTYPE.EQ.5 ) THEN
              XXX=H(1)*XX(1)+H(2)*XX(2)+H(3)*XX(3)
              YYY=H(1)*YY(1)+H(2)*YY(2)+H(3)*YY(3)
              ZZZ=H(1)*ZZ(1)+H(2)*ZZ(2)+H(3)*ZZ(3)
              HX=AX+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RX
              HY=AY+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RY
              HZ=AZ+((XXX-AX)*RX+(YYY-AY)*RY+(ZZZ-AZ)*RZ)/(RX**2+RY**2+RZ**2)*RZ
              COEFX=(XXX-HX)*VAL*RHO*VAL
              COEFY=(YYY-HY)*VAL*RHO*VAL
              COEFZ=(ZZZ-HZ)*VAL*RHO*VAL
            ENDIF
            DO I=1,NN
              PLX(I)=PLX(I)+H(I)*WG*COEFX
              PLY(I)=PLY(I)+H(I)*WG*COEFY
              PLZ(I)=PLZ(I)+H(I)*WG*COEFZ
            ENDDO
          ENDDO
        ENDDO
        
        IF( LTYPE.EQ.1) THEN
          DO I=1,NN
            VECT(6*I-5)=VAL*PLX(I)
          ENDDO
        ELSE IF( LTYPE.EQ.2 ) THEN
          DO I=1,NN
            VECT(6*I-4)=VAL*PLY(I)
          ENDDO
        ELSE IF( LTYPE.EQ.3 ) THEN
          DO I=1,NN
            VECT(6*I-3)=VAL*PLZ(I)
          ENDDO
        ELSE IF( LTYPE.EQ.4 ) THEN
          DO I=1,NN
            VECT(6*I-5)=VAL*RHO*AX*PLX(I)
            VECT(6*I-4)=VAL*RHO*AY*PLY(I)
            VECT(6*I-3)=VAL*RHO*AZ*PLZ(I)
          ENDDO
        ELSEIF( LTYPE.EQ.5 ) THEN
          DO I=1,NN
            VECT(6*I-5)=PLX(I)
            VECT(6*I-4)=PLY(I)
            VECT(6*I-3)=PLZ(I)
          ENDDO
        END IF
      ENDIF
      RETURN

   end subroutine DL_S3
!----------------------------------------------------------------------*
   SUBROUTINE RCV_S4(XX,YY,ZZ,EE,PP,THICK,EDISP,RI,SI,TI,ARRAY)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(4),YY(4),ZZ(4),EDISP(24),ARRAY(11)
! LOCAL VARIABLES
      DIMENSION COD(3,4)
      DIMENSION G1(3),G2(3),G3(3),REF(3)
      DIMENSION PHA12(4),PHA13(4),PHA23(4)
      DIMENSION B1(3,24),B2(3,24),B3(3,24)
      DIMENSION B11(24),B12(24),B13(24)
      DIMENSION B21(24),B22(24),B23(24)
      DIMENSION B31(24),B32(24),B33(24)
      DIMENSION EN(3,4),THE(3,3),AMAT(3,3)
      DIMENSION H(4),HR(4),HS(4)
      DIMENSION ELAS(9)
      DIMENSION SGM(5),EPS(5)
!
! TUNING PARAMETER
!
      DATA SCF/1.2/
      DATA ILOC/1/
!
! SET CORD
!
      DO I=1,4
        COD(1,I)=XX(I)
        COD(2,I)=YY(I)
        COD(3,I)=ZZ(I)
      ENDDO
!
! SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM
      CALL S4REFV(COD,REF,ILOC)
!
! CALCULATE [PHAI] AT EACH NODE
      CALL S4PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!
! ELASTICITY MATRIX
      ELAS(1)=EE/(1.0-PP*PP)
      ELAS(2)=EE/(1.0-PP*PP)*PP
      ELAS(3)=EE/(1.0-PP*PP)
      ELAS(4)=0.0
      ELAS(5)=0.0
      ELAS(6)=EE/(1.0+PP)/2.0
      ELAS(7)=EE/(1.0+PP)/2.0/SCF
      ELAS(8)=0.0
      ELAS(9)=EE/(1.0+PP)/2.0/SCF
!
! ECALUATION POINT (RI,SI,TI)
      CALL S4BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,RI,SI,TI,G1,G2,G3,B1,B2,B3)
!
! JACOBI MATRIX 
      CALL S4JAC(G1,G2,G3,DET,AMAT)
!
! SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
      CALL S4THE(G3,REF,THE)
!
! SET H,HR,HS
      CALL S4PSET(RI,SI,H,HR,HS)
!
! SET COMPONETS OF [B] MATRIX
      CALL S4BM0(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,                  &
                                  B11,B21,B31,B12,B22,B32,B13,B23,B33)
!**
!** SET EPS  {e}=[B]{u}
!**
      EPS(1)=0.0
      EPS(2)=0.0
      EPS(3)=0.0
      EPS(4)=0.0
      EPS(5)=0.0
      DO I=1,24
        EPS(1)=EPS(1)+B11(I)*EDISP(I)
        EPS(2)=EPS(2)+B22(I)*EDISP(I)
        EPS(3)=EPS(3)+(B12(I)+B21(I))*EDISP(I)
        EPS(4)=EPS(4)+(B13(I)+B31(I))*EDISP(I)
        EPS(5)=EPS(5)+(B23(I)+B32(I))*EDISP(I)
      enddo
!**
!** SET SGM  {S}=[D]{e}
!**
      SGM(1)=ELAS(1)*EPS(1)+ELAS(2)*EPS(2)+ELAS(4)*EPS(3)
      SGM(2)=ELAS(2)*EPS(1)+ELAS(3)*EPS(2)+ELAS(5)*EPS(3)
      SGM(3)=ELAS(4)*EPS(1)+ELAS(5)*EPS(2)+ELAS(6)*EPS(3)
      SGM(4)=ELAS(7)*EPS(4)+ELAS(8)*EPS(5)
      SGM(5)=ELAS(8)*EPS(4)+ELAS(9)*EPS(5)
!
!**  DISPLAY EPS & SGM
      PS = (SGM(1)+SGM(2))/3.0
      SMISES = DSQRT(3.0*( 0.5*( (SGM(1)-PS)**2                                &
                                +(SGM(2)-PS)**2                                &
                                +(SGM(4)-PS)**2                                &
                                )                                              &
                           )    +SGM(3)**2                                     &
                    )
      ARRAY(1)=EPS(1)
      ARRAY(2)=EPS(2)
      ARRAY(3)=EPS(3)
      ARRAY(4)=EPS(4)
      ARRAY(5)=EPS(5)
      ARRAY(6)=SGM(1)
      ARRAY(7)=SGM(2)
      ARRAY(8)=SGM(3)
      ARRAY(9)=SGM(4)
      ARRAY(10)=SGM(5)
      ARRAY(11)=SMISES
      RETURN

   end subroutine RCV_S4
!----------------------------------------------------------------------*
   SUBROUTINE S4BM0(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,               &
                                         B11,B21,B31,B12,B22,B32,B13,B23,B33)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION H(4),HR(4),HS(4),PHA12(4),PHA13(4),PHA23(4)
      DIMENSION AMAT(3,3),THE(3,3)
      DIMENSION B11(24),B21(24),B31(24)
      DIMENSION B12(24),B22(24),B32(24)
      DIMENSION B13(24),B23(24),B33(24)
      DIMENSION BMAT(3,3),WORK1(3,3),WORK2(3,3)
!** CLEAR MATRIX
      DO I=1,24
        B11(I)=0.0
        B21(I)=0.0
        B31(I)=0.0
        B12(I)=0.0
        B22(I)=0.0
        B32(I)=0.0
        B13(I)=0.0
        B23(I)=0.0
        B33(I)=0.0
      enddo
      DO INOD=1,4
        B11(6*INOD-5)= HR(INOD)
        B11(6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B11(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B21(6*INOD-4)= HR(INOD)
        B21(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B21(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B31(6*INOD-3)= HR(INOD)
        B31(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B31(6*INOD-1)=-THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B12(6*INOD-5)= HS(INOD)
        B12(6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B12(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B22(6*INOD-4)= HS(INOD)
        B22(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B22(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B32(6*INOD-3)= HS(INOD)
        B32(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B32(6*INOD-1)=-THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B13(6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B13(6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B23(6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B23(6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B33(6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B33(6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
      enddo
      DO I=1,24
        BMAT(1,1)=B11(I)
        BMAT(1,2)=B21(I)
        BMAT(1,3)=B31(I)
        BMAT(2,1)=B12(I)
        BMAT(2,2)=B22(I)
        BMAT(2,3)=B32(I)
        BMAT(3,1)=B13(I)
        BMAT(3,2)=B23(I)
        BMAT(3,3)=B33(I)

        DO II=1,3
          DO JJ=1,3
            WORK1(II,JJ)=0.0
            DO KK=1,3
              WORK1(II,JJ)=WORK1(II,JJ)+AMAT(II,KK)*BMAT(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            WORK2(II,JJ)=0.0
            DO KK=1,3
              WORK2(II,JJ)=WORK2(II,JJ)+WORK1(II,KK)*THE(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            BMAT(II,JJ)=0.0
            DO KK=1,3
              BMAT(II,JJ)=BMAT(II,JJ)+THE(KK,II)*WORK2(KK,JJ)
            enddo
          enddo
        enddo

        B11(I)=BMAT(1,1)
        B21(I)=BMAT(1,2)
        B31(I)=BMAT(1,3)
        B12(I)=BMAT(2,1)
        B22(I)=BMAT(2,2)
        B32(I)=BMAT(2,3)
        B13(I)=BMAT(3,1)
        B23(I)=BMAT(3,2)
        B33(I)=BMAT(3,3)

      enddo
      RETURN

   end subroutine S4BM0
!----------------------------------------------------------------------*
   SUBROUTINE RCV_S3(XX,YY,ZZ,EE,PP,THICK,EDISP,TI,ARRAY)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(3),YY(3),ZZ(3),EDISP(24),ARRAY(11)
! LOCAL VARIABLES
      DIMENSION COD(3,3)
      DIMENSION G1(3),G2(3),G3(3),REF(3)
      DIMENSION PHA12(3),PHA13(3),PHA23(3)
      DIMENSION B1(3,24),B2(3,24),B3(3,24)
      DIMENSION B11(18),B12(18),B13(18)
      DIMENSION B21(18),B22(18),B23(18)
      DIMENSION B31(18),B32(18),B33(18)
      DIMENSION EN(3,3),THE(3,3),AMAT(3,3)
      DIMENSION H(3),HR(3),HS(3)
      DIMENSION ELAS(9)
      DIMENSION SGM(5),EPS(5)
!
! TUNING PARAMETER
!
      DATA SCF/1.2/
      DATA ILOC/1/
!
! SET CORD
!
      DO I=1,3
        COD(1,I)=XX(I)
        COD(2,I)=YY(I)
        COD(3,I)=ZZ(I)
      ENDDO
!
! SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM
      CALL S3REFV(COD,REF,ILOC)
!
! CALCULATE [PHAI] AT EACH NODE
      CALL S3PHAI(COD,REF,PHA12,PHA13,PHA23,EN)
!
! ELASTICITY MATRIX
      ELAS(1)=EE/(1.0-PP*PP)
      ELAS(2)=EE/(1.0-PP*PP)*PP
      ELAS(3)=EE/(1.0-PP*PP)
      ELAS(4)=0.0
      ELAS(5)=0.0
      ELAS(6)=EE/(1.0+PP)/2.0
      ELAS(7)=EE/(1.0+PP)/2.0/SCF
      ELAS(8)=0.0
      ELAS(9)=EE/(1.0+PP)/2.0/SCF
!
! EVALUATION POINT (0.0,0.0,TI)
      ZR=0.0
      CALL S3BMT1(COD,EN,PHA12,PHA13,PHA23,THICK,ZR,ZR,TI,G1,G2,G3,B1,B2,B3)
!
! JACOBI MATRIX 
      CALL S3JAC(G1,G2,G3,DET,AMAT)
!
! SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
      CALL S3THE(G3,REF,THE)
!
! SET H,HR,HS
      CALL S3PSET(ZR,ZR,H,HR,HS)
!
! SET COMPONETS OF [B] MATRIX
      CALL S3BM0(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,                  &
                                  B11,B21,B31,B12,B22,B32,B13,B23,B33)
!**
!** SET EPS  {e}=[B]{u}
!**
      EPS(1)=0.0
      EPS(2)=0.0
      EPS(3)=0.0
      EPS(4)=0.0
      EPS(5)=0.0
      DO I=1,18
      EPS(1)=EPS(1)+B11(I)*EDISP(I)
        EPS(2)=EPS(2)+B22(I)*EDISP(I)
        EPS(3)=EPS(3)+(B12(I)+B21(I))*EDISP(I)
        EPS(4)=EPS(4)+(B13(I)+B31(I))*EDISP(I)
        EPS(5)=EPS(5)+(B23(I)+B32(I))*EDISP(I)
      enddo
!**
!** SET SGM  {S}=[D]{e}
!**
      SGM(1)=ELAS(1)*EPS(1)+ELAS(2)*EPS(2)+ELAS(4)*EPS(3)
      SGM(2)=ELAS(2)*EPS(1)+ELAS(3)*EPS(2)+ELAS(5)*EPS(3)
      SGM(3)=ELAS(4)*EPS(1)+ELAS(5)*EPS(2)+ELAS(6)*EPS(3)
      SGM(4)=ELAS(7)*EPS(4)+ELAS(8)*EPS(5)
      SGM(5)=ELAS(8)*EPS(4)+ELAS(9)*EPS(5)
!
!**  DISPLAY EPS & SGM
      PS = (SGM(1)+SGM(2))/3.0
      SMISES = DSQRT(3.0*(0.5*((SGM(1)-PS)**2+(SGM(2)-PS)**2+(SGM(4)-PS)**2    &
                              )                                                &
                         ) + SGM(3)**2                                         &
                    )                                               
!
      ARRAY(1)=EPS(1)
      ARRAY(2)=EPS(2)
      ARRAY(3)=EPS(3)
      ARRAY(4)=EPS(4)
      ARRAY(5)=EPS(5)
      ARRAY(6)=SGM(1)
      ARRAY(7)=SGM(2)
      ARRAY(8)=SGM(3)
      ARRAY(9)=SGM(4)
      ARRAY(10)=SGM(5)
      ARRAY(11)=SMISES
      RETURN

   end subroutine RCV_S3 
!----------------------------------------------------------------------*
   SUBROUTINE S3BM0(H,HR,HS,TI,THICK,PHA12,PHA13,PHA23,AMAT,THE,               &
                                     B11,B21,B31,B12,B22,B32,B13,B23,B33)
!----------------------------------------------------------------------*
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION H(3),HR(3),HS(3),PHA12(3),PHA13(3),PHA23(3)
      DIMENSION AMAT(3,3),THE(3,3)
      DIMENSION B11(18),B21(18),B31(18)
      DIMENSION B12(18),B22(18),B32(18)
      DIMENSION B13(18),B23(18),B33(18)
      DIMENSION BMAT(3,3),WORK1(3,3),WORK2(3,3)
!** CLEAR MATRIX
      DO I=1,18
        B11(I)=0.0
        B21(I)=0.0
        B31(I)=0.0
        B12(I)=0.0
        B22(I)=0.0
        B32(I)=0.0
        B13(I)=0.0
        B23(I)=0.0
        B33(I)=0.0
      enddo
!
      DO INOD=1,3
        B11(6*INOD-5)= HR(INOD)
        B11(6*INOD-1)= THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B11(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B21(6*INOD-4)= HR(INOD)
        B21(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA12(INOD)
        B21(6*INOD  )= THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B31(6*INOD-3)= HR(INOD)
        B31(6*INOD-2)=-THICK*TI*0.5*HR(INOD)*PHA13(INOD)
        B31(6*INOD-1)=-THICK*TI*0.5*HR(INOD)*PHA23(INOD)
        B12(6*INOD-5)= HS(INOD)
        B12(6*INOD-1)= THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B12(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B22(6*INOD-4)= HS(INOD)
        B22(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA12(INOD)
        B22(6*INOD  )= THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B32(6*INOD-3)= HS(INOD)
        B32(6*INOD-2)=-THICK*TI*0.5*HS(INOD)*PHA13(INOD)
        B32(6*INOD-1)=-THICK*TI*0.5*HS(INOD)*PHA23(INOD)
        B13(6*INOD-1)= THICK*0.5*H(INOD)*PHA12(INOD)
        B13(6*INOD  )= THICK*0.5*H(INOD)*PHA13(INOD)
        B23(6*INOD-2)=-THICK*0.5*H(INOD)*PHA12(INOD)
        B23(6*INOD  )= THICK*0.5*H(INOD)*PHA23(INOD)
        B33(6*INOD-2)=-THICK*0.5*H(INOD)*PHA13(INOD)
        B33(6*INOD-1)=-THICK*0.5*H(INOD)*PHA23(INOD)
      enddo
!
      DO I=1,18
        BMAT(1,1)=B11(I)
        BMAT(1,2)=B21(I)
        BMAT(1,3)=B31(I)
        BMAT(2,1)=B12(I)
        BMAT(2,2)=B22(I)
        BMAT(2,3)=B32(I)
        BMAT(3,1)=B13(I)
        BMAT(3,2)=B23(I)
        BMAT(3,3)=B33(I)
        DO II=1,3
          DO JJ=1,3
            WORK1(II,JJ)=0.0
            DO KK=1,3
              WORK1(II,JJ)=WORK1(II,JJ)+AMAT(II,KK)*BMAT(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            WORK2(II,JJ)=0.0
            DO KK=1,3
              WORK2(II,JJ)=WORK2(II,JJ)+WORK1(II,KK)*THE(KK,JJ)
            enddo
          enddo
        enddo
        DO II=1,3
          DO JJ=1,3
            BMAT(II,JJ)=0.0
            DO KK=1,3
              BMAT(II,JJ)=BMAT(II,JJ)+THE(KK,II)*WORK2(KK,JJ)
            enddo
          enddo
        enddo

        B11(I)=BMAT(1,1)
        B21(I)=BMAT(1,2)
        B31(I)=BMAT(1,3)
        B12(I)=BMAT(2,1)
        B22(I)=BMAT(2,2)
        B32(I)=BMAT(2,3)
        B13(I)=BMAT(3,1)
        B23(I)=BMAT(3,2)
        B33(I)=BMAT(3,3)
 
      enddo
      RETURN

   end subroutine S3BM0
end module m_static_LIB_shell
