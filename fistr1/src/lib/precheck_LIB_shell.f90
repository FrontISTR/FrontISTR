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

!> This module provides function to check input data of shell elements
module m_precheck_LIB_shell
   contains
!***********************************************************************
!  SHELL Element:
!  PRE_731( XX,YY,ZZ,thick,vol,almax,almin )
!  PRE_741( XX,YY,ZZ,thick,vol,almax,almin )
!----------------------------------------------------------------------*
   subroutine PRE_731( XX,YY,ZZ,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!**
!**  Precheck for 3nodes SHELL
!** 
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),thick,vol,almax,almin
! LOCAL VARIABLES
      REAL(kind=kreal) V1X,V1Y,V1Z
      REAL(kind=kreal) V2X,V2Y,V2Z
      REAL(kind=kreal) V3X,V3Y,V3Z
      REAL(kind=kreal) area,a1,a2,a3
      area = 0.0
      vol = 0.0
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
      area=SQRT( V3X*V3X + V3Y*V3Y + V3Z*V3Z )*0.5
      vol = area * thick
      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2+(zz(3)-zz(2))**2 )
      a3 = SQRT( (xx(1)-xx(3))**2+(yy(1)-yy(3))**2+(zz(1)-zz(3))**2 )
      almax = DMAX1( a1,a2,a3 )
      almin = DMIN1( a1,a2,a3 )
      
   end subroutine PRE_731 
!----------------------------------------------------------------------*
   subroutine PRE_741( XX,YY,ZZ,thick,vol,almax,almin )
!----------------------------------------------------------------------*
!**
!**  Precheck for 3nodes SHELL
!** 
      use hecmw
      use gauss_integration
      IMPLICIT NONE
! I/F VARIABLES
      REAL(kind=kreal) XX(*),YY(*),ZZ(*),thick,vol,almax,almin
! LOCAL VARIABLES
      INTEGER(kind=kint) NN
      INTEGER(kind=kint) NG
      PARAMETER(NN=8,NG=2)
      REAL(kind=kreal) H(NN),HR(NN),HS(NN),HT(NN)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      REAL(kind=kreal) XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,DET,WG
      INTEGER(kind=kint) IG1,IG2,LX,LY,LZ,I
      REAL(kind=kreal) VX,VY,VZ,XCOD,YCOD,ZCOD
      REAL(kind=kreal) AX,AY,AZ,RX,RY,RZ,HX,HY,HZ,VAL
      REAL(kind=kreal) PHX,PHY,PHZ
      REAL(kind=kreal) G1X,G1Y,G1Z
      REAL(kind=kreal) G2X,G2Y,G2Z
      REAL(kind=kreal) G3X,G3Y,G3Z
      REAL(kind=kreal) XSUM,COEFX,COEFY,COEFZ
      REAL(kind=kreal) area,a1,a2,a3,a4
!
      area = 0.0
      vol  = 0.0
! INTEGRATION OVER SURFACE
        DO IG2=1,NG
          SI=XG(NG,IG2)
          DO IG1=1,NG
            RI=XG(NG,IG1)
            H(1)=0.25*(1.0-RI)*(1.0-SI)
            H(2)=0.25*(1.0+RI)*(1.0-SI)
            H(3)=0.25*(1.0+RI)*(1.0+SI)
            H(4)=0.25*(1.0-RI)*(1.0+SI)
            HR(1)=-.25*(1.0-SI)
            HR(2)= .25*(1.0-SI)
            HR(3)= .25*(1.0+SI)
            HR(4)=-.25*(1.0+SI)
            HS(1)=-.25*(1.0-RI)
            HS(2)=-.25*(1.0+RI)
            HS(3)= .25*(1.0+RI)
            HS(4)= .25*(1.0-RI)
            G1X=0.0
            G1Y=0.0
            G1Z=0.0
            G2X=0.0
            G2Y=0.0
            G2Z=0.0
            DO I=1,NN
              G1X=G1X+HR(I)*XX(I)
              G1Y=G1Y+HR(I)*YY(I)
              G1Z=G1Z+HR(I)*ZZ(I)
              G2X=G2X+HS(I)*XX(I)
              G2Y=G2Y+HS(I)*YY(I)
              G2Z=G2Z+HS(I)*ZZ(I)
            ENDDO
            G3X=G1Y*G2Z-G1Z*G2Y
            G3Y=G1Z*G2X-G1X*G2Z
            G3Z=G1X*G2Y-G1Y*G2X
            XSUM=DSQRT(G3X**2+G3Y**2+G3Z**2)
            G3X=G3X/XSUM
            G3Y=G3Y/XSUM
            G3Z=G3Z/XSUM
!JACOBI MATRIX
            XJ11=G1X
            XJ12=G1Y
            XJ13=G1Z
            XJ21=G2X
            XJ22=G2Y
            XJ23=G2Z
            XJ31=G3X
            XJ32=G3Y
            XJ33=G3Z
!DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33                                                 &
               +XJ12*XJ23*XJ31                                                 &
               +XJ13*XJ21*XJ32                                                 &
               -XJ13*XJ22*XJ31                                                 &
               -XJ12*XJ21*XJ33                                                 &
               -XJ11*XJ23*XJ32
            WG=WGT(NG,IG1)*WGT(NG,IG2)*DET
            do i = 1, NN
              area = area + H(i)*WG
            enddo
          ENDDO
        ENDDO

      vol = area*thick
      a1 = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
      a2 = SQRT( (xx(3)-xx(2))**2+(yy(3)-yy(2))**2+(zz(3)-zz(2))**2 )
      a3 = SQRT( (xx(4)-xx(3))**2+(yy(4)-yy(3))**2+(zz(4)-zz(3))**2 )
      a4 = SQRT( (xx(1)-xx(4))**2+(yy(1)-yy(4))**2+(zz(1)-zz(4))**2 )
      almax = DMAX1( a1,a2,a3,a4 )
      almin = DMIN1( a1,a2,a3,a4 )

   end subroutine PRE_741
end module m_precheck_LIB_shell
