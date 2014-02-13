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
!> This module provides function to check input data of 2d static analysis
module m_precheck_LIB_2d
   contains
!***********************************************************************
!  2D Element: PreCheck
!    PRE_231( xx,yy,thick,vol,almax,almin )
!    PRE_241( xx,yy,thick,vol,almax,almin )
!    PRE_232( xx,yy,thick,vol,almax,almin )
!    PRE_242( xx,yy,thick,vol,almax,almin )
!***********************************************************************
!----------------------------------------------------------------------*
   subroutine PRE_231( XX,YY,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 2D 3 NODE PLANE ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),thick,vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=3,NG=2)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN)
      INTEGER(kind=kint) L1,L2,I
      REAL(kind=kreal) X1,X2,X3,XL1,XL2
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,WG,area
      REAL(kind=kreal) a1,a2,a3
!C
      area = 0.0
      vol  = 0.0
!* LOOP OVER ALL INTEGRATION POINTS
      DO L2=1,NG
        XL2=XG(NG,L2)
        X2 =(XL2+1.0)*0.5
        DO L1=1,NG
          XL1=XG(NG,L1)
          X1=0.5*(1.0-X2)*(XL1+1.0)
!  INTERPOLATION FUNCTION
          X3=1.0-X1-X2
          H(1)=X1
          H(2)=X2
          H(3)=X3
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR L1-COORDINATE
          HL1(1)=1.0
          HL1(2)=0.0
          HL1(3)=0.0
!  FOR L2-COORDINATE
          HL2(1)=0.0
          HL2(2)=1.0
          HL2(3)=0.0
!  FOR L3-COORDINATE
          HL3(1)=0.0
          HL3(2)=0.0
          HL3(3)=1.0
!  JACOBI MATRIX
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
          WG=WGT(NG,L1)*WGT(NG,L2)*DET*(1.0-X2)*0.25
          DO I = 1, NN
            area = area + H(I)*WG
          ENDDO
        ENDDO
      ENDDO

      vol = area * thick
      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2 )
      a3 = SQRT( (xx(1)-xx(3))**2+(yy(1)-yy(3))**2 )
      almax = DMAX1( a1,a2,a3 )
      almin = DMIN1( a1,a2,a3 )

   end subroutine PRE_231
!----------------------------------------------------------------------*
   subroutine PRE_241( XX,YY,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 2D 4 NODE PLANE ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),thick,vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=4,NG=2)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN)
      INTEGER(kind=kint) LX,LY,I
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,WG,area
      REAL(kind=kreal) a1,a2,a3,a4
!C
      area = 0.0
      vol  = 0.0
!* LOOP OVER ALL INTEGRATION POINTS
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
          WG=WGT(NG,LX)*WGT(NG,LY)*DET
          DO I = 1, NN
            area = area + H(I)*WG
          ENDDO
        ENDDO
      ENDDO

      vol = area * thick
      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2 )
      a3 = SQRT( (xx(4)-xx(3))**2+(yy(4)-yy(3))**2 )
      a4 = SQRT( (xx(1)-xx(4))**2+(yy(1)-yy(4))**2 )
      almax = DMAX1( a1,a2,a3,a4 )
      almin = DMIN1( a1,a2,a3,a4 )

   end subroutine PRE_241
!----------------------------------------------------------------------*
   subroutine PRE_232( XX,YY,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 2D 6 NODE PLANE ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),thick,vol,tline,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=6,NG=3)
      REAL(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN)
      INTEGER(kind=kint) L1,L2,I
      REAL(kind=kreal) X1,X2,X3,XL1,XL2
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,WG,area
      REAL(kind=kreal)  a1,a2,AL1,AL2,AL3

!*************************
      area = 0.0
      vol  = 0.0
!* LOOP OVER ALL INTEGRATION POINTS
      DO L2=1,NG
        XL2=XG(NG,L2)
        X2 =(XL2+1.0)*0.5
        DO L1=1,NG
          XL1=XG(NG,L1)
          X1=0.5*(1.0-X2)*(XL1+1.0)
!  INTERPOLATION FUNCTION
          X3=1.0-X1-X2
          H(1)= X1*(2.0*X1-1.)
          H(2)= X2*(2.0*X2-1.)
          H(3)= X3*(2.0*X3-1.)
          H(4)= 4.0*X1*X2
          H(5)= 4.0*X2*X3
          H(6)= 4.0*X1*X3
!  DERIVATIVE OF INTERPOLATION FUNCTION
!  FOR L1-COORDINATE
          HL1(1)=4.0*X1-1.0
          HL1(2)= 0.0
          HL1(3)= 0.0
          HL1(4)= 4.0*X2
          HL1(5)= 0.0
          HL1(6)= 4.0*X3
!  FOR L2-COORDINATE
          HL2(1)= 0.0
          HL2(2)= 4.0*X2-1.0
          HL2(3)= 0.0
          HL2(4)= 4.0*X1
          HL2(5)= 4.0*X3
          HL2(6)= 0.0
!  FOR L3-COORDINATE
          HL3(1)= 0.0
          HL3(2)= 0.0
          HL3(3)= 4.0*X3-1.0
          HL3(4)= 0.0
          HL3(5)= 4.0*X2
          HL3(6)= 4.0*X1
!  JACOBI MATRIX
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
          WG=WGT(NG,L1)*WGT(NG,L2)*DET*(1.0-X2)*0.25
          DO I = 1, NN
            area = area + H(I)*WG
          ENDDO
        ENDDO
      ENDDO

      vol = area * thick
      a1 = SQRT( (xx(4)-xx(1))**2+(yy(4)-yy(1))**2 )
      a2 = SQRT( (xx(2)-xx(4))**2+(yy(2)-yy(4))**2 )
      AL1 = a1 + a2
      a1 = SQRT( (xx(5)-xx(2))**2+(yy(5)-yy(2))**2 )
      a2 = SQRT( (xx(3)-xx(5))**2+(yy(3)-yy(5))**2 )
      AL2 = a1 + a2
      a1 = SQRT( (xx(6)-xx(3))**2+(yy(6)-yy(3))**2 )
      a2 = SQRT( (xx(1)-xx(6))**2+(yy(1)-yy(6))**2 )
      AL3 = a1 + a2
      almax = DMAX1( AL1,AL2,AL3 )
      almin = DMIN1( AL1,AL2,AL3 )
   
   end subroutine  PRE_232
!----------------------------------------------------------------------*
   subroutine PRE_242( XX,YY,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!
! CALCULATION 2D 8 NODE PLANE ELEMENT
!
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),thick,vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=8,NG=3)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN)
      INTEGER(kind=kint) LX,LY,I
      REAL(kind=kreal) RI,SI,RP,SP,RM,SM
      REAL(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,WG,area
      REAL(kind=kreal) RX(8),SX(8)
      REAL(kind=kreal)  a1,a2,AL1,AL2,AL3,AL4
      DATA RX/-1., 1.,1.,-1., 0., 1., 0., -1./
      DATA SX/-1.,-1.,1., 1.,-1., 0., 1.,  0./
!*************************
      area = 0.0
      vol  = 0.0
!* LOOP OVER ALL INTEGRATION POINTS
      DO LX=1,NG
        RI=XG(NG,LX)
        DO LY=1,NG
          SI=XG(NG,LY)
          RP=1.0+RI
          SP=1.0+SI
          RM=1.0-RI
          SM=1.0-SI
          H(1)=0.25*RM*SM*(-1.0-RI-SI)
          H(2)=0.25*RP*SM*(-1.0+RI-SI)
          H(3)=0.25*RP*SP*(-1.0+RI+SI)
          H(4)=0.25*RM*SP*(-1.0-RI+SI)
          H(5)=0.5*(1.0-RI*RI)*(1.0-SI)
          H(6)=0.5*(1.0-SI*SI)*(1.0+RI)
          H(7)=0.5*(1.0-RI*RI)*(1.0+SI)
          H(8)=0.5*(1.0-SI*SI)*(1.0-RI)
          HR(1)=-.25*SM*(-1.0-RI-SI)-0.25*RM*SM
          HR(2)= .25*SM*(-1.0+RI-SI)+0.25*RP*SM
          HR(3)= .25*SP*(-1.0+RI+SI)+0.25*RP*SP
          HR(4)=-.25*SP*(-1.0-RI+SI)-0.25*RM*SP
          HR(5)=-RI*(1.0-SI)
          HR(6)= 0.5*(1.0-SI*SI)
          HR(7)=-RI*(1.0+SI)
          HR(8)=-0.5*(1.0-SI*SI)
          HS(1)=-.25*RM*(-1.0-RI-SI)-0.25*RM*SM
          HS(2)=-.25*RP*(-1.0+RI-SI)-0.25*RP*SM
          HS(3)= .25*RP*(-1.0+RI+SI)+0.25*RP*SP
          HS(4)= .25*RM*(-1.0-RI+SI)+0.25*RM*SP
          HS(5)=-0.5*(1.0-RI*RI)
          HS(6)=-SI *(1.0+RI)
          HS(7)= 0.5*(1.0-RI*RI)
          HS(8)= -SI*(1.0-RI)
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
          WG=WGT(NG,LX)*WGT(NG,LY)*DET
          DO I = 1, NN
            area = area + H(I)*WG
          ENDDO
        ENDDO
      ENDDO

      vol = area * thick
      a1 = SQRT( (xx(5)-xx(1))**2+(yy(5)-yy(1))**2 )
      a2 = SQRT( (xx(2)-xx(5))**2+(yy(2)-yy(5))**2 )
      AL1 = a1 + a2
      a1 = SQRT( (xx(6)-xx(2))**2+(yy(6)-yy(2))**2 )
      a2 = SQRT( (xx(3)-xx(6))**2+(yy(3)-yy(6))**2 )
      AL2 = a1 + a2
      a1 = SQRT( (xx(7)-xx(3))**2+(yy(7)-yy(3))**2 )
      a2 = SQRT( (xx(4)-xx(7))**2+(yy(4)-yy(7))**2 )
      AL3 = a1 + a2
      a1 = SQRT( (xx(8)-xx(4))**2+(yy(8)-yy(4))**2 )
      a2 = SQRT( (xx(1)-xx(8))**2+(yy(1)-yy(8))**2 )
      AL4 = a1 + a2
      almax = DMAX1( AL1,AL2,AL3,AL4 )
      almin = DMIN1( AL1,AL2,AL3,AL4 )
      
   end subroutine PRE_242
end module m_precheck_LIB_2d
