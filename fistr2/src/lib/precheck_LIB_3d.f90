!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provides function to check input data of 3d static analysis
module m_precheck_LIB_3d
   contains
!***********************************************************************
!  3D SOLID Element: PreCheck
!  PRE_341(XX,YY,ZZ,vol,almax,almin)
!  PRE_351(XX,YY,ZZ,vol,almax,almin)
!  PRE_361(XX,YY,ZZ,vol,almax,almin)
!  PRE_342(XX,YY,ZZ,vol,almax,almin)
!  PRE_352(XX,YY,ZZ,vol,almax,almin)
!  PRE_362(XX,YY,ZZ,vol,almax,almin)
!----------------------------------------------------------------------*
   subroutine PRE_341( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 4-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=4,NG=2)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HL4(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) I,L,L1,L2,L3
      REAL(kind=kreal) XL1,XL2,XL3
      REAL(kind=kreal) X1,X2,X3,X4
      REAL(kind=kreal) a1,a2,a3,a4,a5,a6
!C
      vol = 0.0
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
! INTERPOLATION FUNCTION
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
! FOR L2-COORDINATE
            HL2(1)= 0.0 
            HL2(2)= 1.0 
            HL2(3)= 0.0
            HL2(4)= 0.0
! FOR L3-COORDINATE
            HL3(1)= 0.0
            HL3(2)= 0.0
            HL3(3)= 1.0 
            HL3(4)= 0.0
! FOR L4-COORDINATE
            HL4(1)= 0.0
            HL4(2)= 0.0
            HL4(3)= 0.0 
            HL4(4)= 1.0
! JACOBI MATRIX 
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
              XJ11=XJ11+(HL4(I)-HL1(I))*XX(I)
              XJ21=XJ21+(HL4(I)-HL2(I))*XX(I)
              XJ31=XJ31+(HL4(I)-HL3(I))*XX(I)
              XJ12=XJ12+(HL4(I)-HL1(I))*YY(I)
              XJ22=XJ22+(HL4(I)-HL2(I))*YY(I)
              XJ32=XJ32+(HL4(I)-HL3(I))*YY(I)
              XJ13=XJ13+(HL4(I)-HL1(I))*ZZ(I)
              XJ23=XJ23+(HL4(I)-HL2(I))*ZZ(I)
              XJ33=XJ33+(HL4(I)-HL3(I))*ZZ(I)
            ENDDO
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
            DO I = 1, NN
              vol = vol + H(I)*WG
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2+(zz(3)-zz(2))**2 )
      a3 = SQRT( (xx(1)-xx(3))**2+(yy(1)-yy(3))**2+(zz(1)-zz(3))**2 )
      a4 = SQRT( (xx(4)-xx(1))**2+(yy(4)-yy(1))**2+(zz(4)-zz(1))**2 )
      a5 = SQRT( (xx(4)-xx(2))**2+(yy(4)-yy(2))**2+(zz(4)-zz(2))**2 )
      a6 = SQRT( (xx(4)-xx(3))**2+(yy(4)-yy(3))**2+(zz(4)-zz(3))**2 )
      almax = DMAX1( a1,a2,a3,a4,a5,a6 )
      almin = DMIN1( a1,a2,a3,a4,a5,a6 )

   end subroutine PRE_341
!----------------------------------------------------------------------*
   subroutine PRE_351( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 6-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      INTEGER(kind=kint) nline
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,tline,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=6,NG=2)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HZ(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      REAL(kind=kreal) XJI11,XJI21,XJI31,XJI12,XJI22,XJI32,XJI13,XJI23,XJI33
      INTEGER(kind=kint) I,L1,L2,LZ
      REAL(kind=kreal) X1,X2,X3,XL1,XL2,ZI
      REAL(kind=kreal) a1,a2,a3,a4,a5,a6,a7,a8,a9
!C
      vol = 0.0
! LOOP FOR INTEGRATION POINTS
      DO LZ=1,NG
        ZI=XG(NG,LZ)
        DO L2=1,NG
          XL2=XG(NG,L2)
          X2 =(XL2+1.0)*0.5
          DO L1=1,NG
            XL1=XG(NG,L1)
            X1=0.5*(1.0-X2)*(XL1+1.0)
! INTERPOLATION FUNCTION
            X3=1.0-X1-X2
            H(1)=X1*(1.0-ZI)*0.5
            H(2)=X2*(1.0-ZI)*0.5
            H(3)=X3*(1.0-ZI)*0.5
            H(4)=X1*(1.0+ZI)*0.5
            H(5)=X2*(1.0+ZI)*0.5
            H(6)=X3*(1.0+ZI)*0.5
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR L1-COORDINATE
            HL1(1)=(1.0-ZI)*0.5
            HL1(2)=0.0
            HL1(3)=0.0
            HL1(4)=(1.0+ZI)*0.5
            HL1(5)=0.0
            HL1(6)=0.0
! FOR L2-COORDINATE
            HL2(1)=0.0
            HL2(2)=(1.0-ZI)*0.5
            HL2(3)=0.0
            HL2(4)=0.0
            HL2(5)=(1.0+ZI)*0.5
            HL2(6)=0.0
! FOR L3-COORDINATE
            HL3(1)=0.0
            HL3(2)=0.0
            HL3(3)=(1.0-ZI)*0.5
            HL3(4)=0.0
            HL3(5)=0.0
            HL3(6)=(1.0+ZI)*0.5
! FOR Z-COORDINATE
            HZ(1)=-X1*0.5
            HZ(2)=-X2*0.5
            HZ(3)=-X3*0.5
            HZ(4)= X1*0.5
            HZ(5)= X2*0.5
            HZ(6)= X3*0.5
! JACOBI MATRIX 
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
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,LZ)*DET*(1.0-X2)*0.25
            DO I = 1, NN
              vol = vol + H(I)*WG 
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2+(zz(3)-zz(2))**2 )
      a3 = SQRT( (xx(1)-xx(3))**2+(yy(1)-yy(3))**2+(zz(1)-zz(3))**2 )
      a4 = SQRT( (xx(5)-xx(4))**2+(yy(5)-yy(4))**2+(zz(5)-zz(4))**2 )
      a5 = SQRT( (xx(6)-xx(5))**2+(yy(6)-yy(5))**2+(zz(6)-zz(5))**2 )
      a6 = SQRT( (xx(4)-xx(6))**2+(yy(4)-yy(6))**2+(zz(4)-zz(6))**2 )
      a7 = SQRT( (xx(4)-xx(1))**2+(yy(4)-yy(1))**2+(zz(4)-zz(1))**2 )
      a8 = SQRT( (xx(5)-xx(2))**2+(yy(5)-yy(2))**2+(zz(5)-zz(2))**2 )
      a9 = SQRT( (xx(6)-xx(3))**2+(yy(6)-yy(3))**2+(zz(6)-zz(3))**2 )
      almax = DMAX1( a1,a2,a3,a4,a5,a6,a7,a8,a9 )
      almin = DMIN1( a1,a2,a3,a4,a5,a6,a7,a8,a9 )
      
   end subroutine PRE_351
!----------------------------------------------------------------------*
   subroutine PRE_361( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 8-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=8,NG=2)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) I,LX,LY,LZ
      REAL(kind=kreal) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
!C
      vol = 0.0
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
! INTERPOLATION FUNCTION
            H(1)=0.125*RM*SM*TM
            H(2)=0.125*RP*SM*TM
            H(3)=0.125*RP*SP*TM
            H(4)=0.125*RM*SP*TM
            H(5)=0.125*RM*SM*TP
            H(6)=0.125*RP*SM*TP
            H(7)=0.125*RP*SP*TP
            H(8)=0.125*RM*SP*TP
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR R-COORDINATE
            HR(1)=-.125*SM*TM
            HR(2)= .125*SM*TM
            HR(3)= .125*SP*TM
            HR(4)=-.125*SP*TM
            HR(5)=-.125*SM*TP
            HR(6)= .125*SM*TP
            HR(7)= .125*SP*TP
            HR(8)=-.125*SP*TP
! FOR S-COORDINATE
            HS(1)=-.125*RM*TM
            HS(2)=-.125*RP*TM
            HS(3)= .125*RP*TM
            HS(4)= .125*RM*TM
            HS(5)=-.125*RM*TP
            HS(6)=-.125*RP*TP
            HS(7)= .125*RP*TP
            HS(8)= .125*RM*TP
! FOR T-COORDINATE
            HT(1)=-.125*RM*SM
            HT(2)=-.125*RP*SM
            HT(3)=-.125*RP*SP
            HT(4)=-.125*RM*SP
            HT(5)= .125*RM*SM
            HT(6)= .125*RP*SM
            HT(7)= .125*RP*SP
            HT(8)= .125*RM*SP
! JACOBI MATRIX 
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
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,LX)*WGT(NG,LY)*WGT(NG,LZ)*DET
            DO I=1,NN
              vol = vol + H(i)*WG
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2+(zz(3)-zz(2))**2 )
      a3 = SQRT( (xx(4)-xx(3))**2+(yy(4)-yy(3))**2+(zz(4)-zz(3))**2 )
      a4 = SQRT( (xx(1)-xx(4))**2+(yy(1)-yy(4))**2+(zz(1)-zz(4))**2 )

      a5 = SQRT( (xx(6)-xx(5))**2+(yy(6)-yy(5))**2+(zz(6)-zz(5))**2 )
      a6 = SQRT( (xx(7)-xx(6))**2+(yy(7)-yy(6))**2+(zz(7)-zz(6))**2 )
      a7 = SQRT( (xx(8)-xx(7))**2+(yy(8)-yy(7))**2+(zz(8)-zz(7))**2 )
      a8 = SQRT( (xx(5)-xx(8))**2+(yy(5)-yy(8))**2+(zz(5)-zz(8))**2 )

      a9 = SQRT( (xx(5)-xx(1))**2+(yy(5)-yy(1))**2+(zz(5)-zz(1))**2 )
      a10 = SQRT( (xx(6)-xx(2))**2+(yy(6)-yy(2))**2+(zz(6)-zz(2))**2 )
      a11 = SQRT( (xx(7)-xx(3))**2+(yy(7)-yy(3))**2+(zz(7)-zz(3))**2 )
      a12 = SQRT( (xx(8)-xx(4))**2+(yy(8)-yy(4))**2+(zz(8)-zz(4))**2 )
      almax = DMAX1( a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12 )
      almin = DMIN1( a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12 )

   end subroutine PRE_361
!----------------------------------------------------------------------*
   subroutine PRE_342( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 10-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration 
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=10,NG=3)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HL4(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) I,L1,L2,L3
      REAL(kind=kreal) XL1,XL2,XL3
      REAL(kind=kreal) X1,X2,X3,X4
      REAL(kind=kreal)  a1,a2,AL1,AL2,AL3,AL4,AL5,AL6
!
   VOL = 0.0
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
! INTERPOLATION FUNCTION
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
! FOR L2-COORDINATE
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
! FOR L3-COORDINATE
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
! FOR L4-COORDINATE
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
! JACOBI MATRIX 
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
              XJ11=XJ11+(HL4(I)-HL1(I))*XX(I)
              XJ21=XJ21+(HL4(I)-HL2(I))*XX(I)
              XJ31=XJ31+(HL4(I)-HL3(I))*XX(I)
              XJ12=XJ12+(HL4(I)-HL1(I))*YY(I)
              XJ22=XJ22+(HL4(I)-HL2(I))*YY(I)
              XJ32=XJ32+(HL4(I)-HL3(I))*YY(I)
              XJ13=XJ13+(HL4(I)-HL1(I))*ZZ(I)
              XJ23=XJ23+(HL4(I)-HL2(I))*ZZ(I)
              XJ33=XJ33+(HL4(I)-HL3(I))*ZZ(I)
            ENDDO
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
            DO I = 1, NN
              vol = vol + H(I)*WG
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      a1 = SQRT( (xx(5)-xx(1))**2+(yy(5)-yy(1))**2+(zz(5)-zz(1))**2 )
      a2 = SQRT( (xx(2)-xx(5))**2+(yy(2)-yy(5))**2+(zz(2)-zz(5))**2 )
      AL1 = a1 + a2
      a1 = SQRT( (xx(6)-xx(2))**2+(yy(6)-yy(2))**2+(zz(6)-zz(2))**2 )
      a2 = SQRT( (xx(3)-xx(6))**2+(yy(3)-yy(6))**2+(zz(3)-zz(6))**2 )
      AL2 = a1 + a2
      a1 = SQRT( (xx(7)-xx(3))**2+(yy(7)-yy(3))**2+(zz(7)-zz(3))**2 )
      a2 = SQRT( (xx(1)-xx(7))**2+(yy(1)-yy(7))**2+(zz(1)-zz(7))**2 )
      AL3 = a1 + a2
      a1 = SQRT( (xx(8)-xx(1))**2+(yy(8)-yy(1))**2+(zz(8)-zz(1))**2 )
      a2 = SQRT( (xx(4)-xx(8))**2+(yy(4)-yy(8))**2+(zz(4)-zz(8))**2 )
      AL4 = a1 + a2
      a1 = SQRT( (xx(9)-xx(2))**2+(yy(9)-yy(2))**2+(zz(9)-zz(2))**2 )
      a2 = SQRT( (xx(4)-xx(9))**2+(yy(4)-yy(9))**2+(zz(4)-zz(9))**2 )
      AL5 = a1 + a2
      a1 = SQRT( (xx(10)-xx(3))**2+(yy(10)-yy(3))**2+(zz(10)-zz(3))**2 )
      a2 = SQRT( (xx(4)-xx(10))**2+(yy(4)-yy(10))**2+(zz(4)-zz(10))**2 )
      AL6 = a1 + a2
      almax = DMAX1( AL1,AL2,AL3,AL4,AL5,AL6 )
      almin = DMIN1( AL1,AL2,AL3,AL4,AL5,AL6 )

   end subroutine PRE_342
!----------------------------------------------------------------------*
   subroutine PRE_352( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 15-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration 
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,tline,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=15,NG=3)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN),HZ(NN)
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) I,L1,L2,LZ
      REAL(kind=kreal) X1,X2,X3,XL1,XL2,ZI
      REAL(kind=kreal) a1,a2,AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9
!C
      vol = 0.0
! LOOP FOR INTEGRATION POINTS      
      DO LZ=1,NG
        ZI=XG(NG,LZ)
        DO L2=1,NG
          XL2=XG(NG,L2)
          X2 =(XL2+1.0)*0.5
          DO L1=1,NG
            XL1=XG(NG,L1)
            X1=0.5*(1.0-X2)*(XL1+1.0)
! INTERPOLATION FUNCTION
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
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR L1-COORDINATE
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
! FOR L2-COORDINATE
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
! FOR L3-COORDINATE
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
! FOR Z-COORDINATE
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
! JACOBI MATRIX 
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
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,L1)*WGT(NG,L2)*WGT(NG,LZ)*DET*(1.0-X2)*0.25
            DO I = 1, NN
              vol = vol + H(I)*WG
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      a1 = SQRT( (xx(7)-xx(1))**2+(yy(7)-yy(1))**2+(zz(7)-zz(1))**2 )
      a2 = SQRT( (xx(2)-xx(7))**2+(yy(2)-yy(7))**2+(zz(2)-zz(7))**2 )
      AL1 = a1 + a2
      a1 = SQRT( (xx(8)-xx(2))**2+(yy(8)-yy(2))**2+(zz(8)-zz(2))**2 )
      a2 = SQRT( (xx(3)-xx(8))**2+(yy(3)-yy(8))**2+(zz(3)-zz(8))**2 )
      AL2 = a1 + a2
      a1 = SQRT( (xx(9)-xx(3))**2+(yy(9)-yy(3))**2+(zz(9)-zz(3))**2 )
      a2 = SQRT( (xx(1)-xx(9))**2+(yy(1)-yy(9))**2+(zz(1)-zz(9))**2 )
      AL3 = a1 + a2

      a1 = SQRT( (xx(10)-xx(4))**2+(yy(10)-yy(4))**2+(zz(10)-zz(4))**2 )
      a2 = SQRT( (xx(5)-xx(10))**2+(yy(5)-yy(10))**2+(zz(5)-zz(10))**2 )
      AL4 = a1 + a2
      a1 = SQRT( (xx(11)-xx(5))**2+(yy(11)-yy(5))**2+(zz(11)-zz(5))**2 )
      a2 = SQRT( (xx(6)-xx(11))**2+(yy(6)-yy(11))**2+(zz(6)-zz(11))**2 )
      AL5 = a1 + a2
      a1 = SQRT( (xx(12)-xx(6))**2+(yy(12)-yy(6))**2+(zz(12)-zz(6))**2 )
      a2 = SQRT( (xx(4)-xx(12))**2+(yy(4)-yy(12))**2+(zz(4)-zz(12))**2 )
      AL6 = a1 + a2

      a1 = SQRT( (xx(13)-xx(1))**2+(yy(13)-yy(1))**2+(zz(13)-zz(1))**2 )
      a2 = SQRT( (xx(4)-xx(13))**2+(yy(4)-yy(13))**2+(zz(4)-zz(13))**2 )
      AL7 = a1 + a2
      a1 = SQRT( (xx(14)-xx(2))**2+(yy(14)-yy(2))**2+(zz(14)-zz(2))**2 )
      a2 = SQRT( (xx(5)-xx(14))**2+(yy(5)-yy(14))**2+(zz(5)-zz(14))**2 )
      AL8 = a1 + a2
      a1 = SQRT( (xx(15)-xx(3))**2+(yy(15)-yy(3))**2+(zz(15)-zz(3))**2 )
      a2 = SQRT( (xx(6)-xx(15))**2+(yy(6)-yy(15))**2+(zz(6)-zz(15))**2 )
      AL9 = a1 + a2

      almax = DMAX1( AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9 )
      almin = DMIN1( AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9 )
      
   end subroutine PRE_352
!----------------------------------------------------------------------*
   subroutine PRE_362( XX,YY,ZZ,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 3D 20-NODE SOLID ELEMENT
!
      use hecmw
      use gauss_integration 
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=20,NG=3)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) I,LX,LY,LZ
      REAL(kind=kreal)  a1,a2,AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9,AL10,AL11,AL12
!C
      vol = 0.0
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
! INTERPOLATION FUNCTION
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
! DERIVATIVE OF INTERPOLATION FUNCTION
! FOR R-COORDINATE
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
! FOR S-COORDINATE
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
! FOR T-COORDINATE
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
! JACOBI MATRIX 
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
! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
! WEIGHT VALUE AT GAUSSIAN POINT
            WG=WGT(NG,LX)*WGT(NG,LY)*WGT(NG,LZ)*DET
            DO I=1,NN
              vol = vol + H(I)*WG
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
      a1 = SQRT( (xx(9)-xx(1))**2+(yy(9)-yy(1))**2+(zz(9)-zz(1))**2 )
      a2 = SQRT( (xx(2)-xx(9))**2+(yy(2)-yy(9))**2+(zz(2)-zz(9))**2 )
      AL1 = a1 + a2
      a1 = SQRT( (xx(10)-xx(2))**2+(yy(10)-yy(2))**2+(zz(10)-zz(2))**2 )
      a2 = SQRT( (xx(3)-xx(10))**2+(yy(3)-yy(10))**2+(zz(3)-zz(10))**2 )
      AL2 = a1 + a2
      a1 = SQRT( (xx(11)-xx(3))**2+(yy(11)-yy(3))**2+(zz(11)-zz(3))**2 )
      a2 = SQRT( (xx(4)-xx(11))**2+(yy(4)-yy(11))**2+(zz(4)-zz(11))**2 )
      AL3 = a1 + a2
      a1 = SQRT( (xx(12)-xx(4))**2+(yy(12)-yy(4))**2+(zz(12)-zz(4))**2 )
      a2 = SQRT( (xx(1)-xx(12))**2+(yy(1)-yy(12))**2+(zz(1)-zz(12))**2 )
      AL4 = a1 + a2
      
      a1 = SQRT( (xx(13)-xx(5))**2+(yy(13)-yy(5))**2+(zz(13)-zz(5))**2 )
      a2 = SQRT( (xx(6)-xx(13))**2+(yy(6)-yy(13))**2+(zz(6)-zz(13))**2 )
      AL5 = a1 + a2
      a1 = SQRT( (xx(14)-xx(6))**2+(yy(14)-yy(6))**2+(zz(14)-zz(6))**2 )
      a2 = SQRT( (xx(7)-xx(14))**2+(yy(7)-yy(14))**2+(zz(7)-zz(14))**2 )
      AL6 = a1 + a2
      a1 = SQRT( (xx(15)-xx(7))**2+(yy(15)-yy(7))**2+(zz(15)-zz(7))**2 )
      a2 = SQRT( (xx(8)-xx(15))**2+(yy(8)-yy(15))**2+(zz(8)-zz(15))**2 )
      AL7 = a1 + a2
      a1 = SQRT( (xx(16)-xx(8))**2+(yy(16)-yy(8))**2+(zz(16)-zz(8))**2 )
      a2 = SQRT( (xx(5)-xx(16))**2+(yy(5)-yy(16))**2+(zz(5)-zz(16))**2 )
      AL8 = a1 + a2
      
      a1 = SQRT( (xx(17)-xx(1))**2+(yy(17)-yy(1))**2+(zz(17)-zz(1))**2 )
      a2 = SQRT( (xx(5)-xx(17))**2+(yy(5)-yy(17))**2+(zz(5)-zz(17))**2 )
      AL9 = a1 + a2
      a1 = SQRT( (xx(18)-xx(2))**2+(yy(18)-yy(2))**2+(zz(18)-zz(2))**2 )
      a2 = SQRT( (xx(6)-xx(18))**2+(yy(6)-yy(18))**2+(zz(6)-zz(18))**2 )
      AL10 = a1 + a2
      a1 = SQRT( (xx(19)-xx(3))**2+(yy(19)-yy(3))**2+(zz(19)-zz(3))**2 )
      a2 = SQRT( (xx(7)-xx(19))**2+(yy(7)-yy(19))**2+(zz(7)-zz(19))**2 )
      AL11 = a1 + a2
      a1 = SQRT( (xx(20)-xx(4))**2+(yy(20)-yy(4))**2+(zz(20)-zz(4))**2 )
      a2 = SQRT( (xx(8)-xx(20))**2+(yy(8)-yy(20))**2+(zz(8)-zz(20))**2 )
      AL12 = a1 + a2

      almax = DMAX1( AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9,AL10,AL11,AL12 )
      almin = DMIN1( AL1,AL2,AL3,AL4,AL5,AL6,AL7,AL8,AL9 ,AL10,AL11,AL12)

   end subroutine PRE_362
end module m_precheck_LIB_3d
