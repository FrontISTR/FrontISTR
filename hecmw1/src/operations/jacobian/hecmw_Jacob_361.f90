!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2008/03/13                                         !
!        Category : Jacobian calculation                               !
!                                                                      !
!            Written by Satoshi Ito (Univ. of Tokyo)                   !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_Jacob361
contains

   SUBROUTINE hecmw_Jacob_361 ( hecMESH, iElem, DET, W, N, NX, NY, NZ)
      use hecmw_util

      implicit none

      type(hecmwST_local_mesh):: hecMESH
      integer(kind=kint)::       iElem
      real(kind=kreal)::         DET
      real(kind=kreal)::         W(8), N(8,8),NX(8,8),NY(8,8),NZ(8,8)

      integer(kind=kint):: iLocal, j, jj
      integer(kind=kint):: LX, LY, LZ
      real(kind=kreal)::   DUM
      real(kind=kreal)::   XX(8), YY(8), ZZ(8)
      real(kind=kreal)::   XG(2), WGT(2), HR(8), HS(8), HT(8)
      real(kind=kreal)::   H(2,2,2,8),BX(2,2,2,8),BY(2,2,2,8),BZ(2,2,2,8)
      real(kind=kreal)::   RI, SI, TI, RP, SP, TP, RM, SM, TM
      real(kind=kreal)::   XJ11,XJ12,XJ13,XJ21,XJ22,XJ23,XJ31,XJ32,XJ33
      real(kind=kreal)::   XJI11,XJI12,XJI13,XJI21,XJI22,XJI23,XJI31,XJI32,XJI33

      DATA WGT/1.0D0,1.0D0/
      DATA XG/-0.5773502691896D0, 0.5773502691896D0/
!C

      do j = 1, 8
         jj = 8 - j
         iLocal = hecMESH%elem_node_item  ( 8*iElem -jj )
         XX(j)  = hecMESH%node( iLocal*3 -2 )
         YY(j)  = hecMESH%node( iLocal*3 -1 )
         ZZ(j)  = hecMESH%node( iLocal*3    )
         W(j)   = 1.0D0
      end do

      DO LX=1,2
        RI=XG(LX)
        DO LY=1,2
          SI=XG(LY)
          DO LZ=1,2
            TI=XG(LZ)
!C            
            RP=1.0+RI
            SP=1.0+SI
            TP=1.0+TI
            RM=1.0-RI
            SM=1.0-SI
            TM=1.0-TI
!C
!C*INTERPOLATION FUNCTION
            H(LX,LY,LZ,1)=0.125*RM*SM*TM
            H(LX,LY,LZ,2)=0.125*RP*SM*TM
            H(LX,LY,LZ,3)=0.125*RP*SP*TM
            H(LX,LY,LZ,4)=0.125*RM*SP*TM
            H(LX,LY,LZ,5)=0.125*RM*SM*TP
            H(LX,LY,LZ,6)=0.125*RP*SM*TP
            H(LX,LY,LZ,7)=0.125*RP*SP*TP
            H(LX,LY,LZ,8)=0.125*RM*SP*TP
!C
!C*DERIVATIVE OF INTERPOLATION FUNCTION
!C*  FOR R-COORDINATE
            HR(1)=-.125*SM*TM
            HR(2)= .125*SM*TM
            HR(3)= .125*SP*TM
            HR(4)=-.125*SP*TM
            HR(5)=-.125*SM*TP
            HR(6)= .125*SM*TP
            HR(7)= .125*SP*TP
            HR(8)=-.125*SP*TP
!C*  FOR S-COORDINATE
            HS(1)=-.125*RM*TM
            HS(2)=-.125*RP*TM
            HS(3)= .125*RP*TM
            HS(4)= .125*RM*TM
            HS(5)=-.125*RM*TP
            HS(6)=-.125*RP*TP
            HS(7)= .125*RP*TP
            HS(8)= .125*RM*TP
!C*  FOR T-COORDINATE
            HT(1)=-.125*RM*SM
            HT(2)=-.125*RP*SM
            HT(3)=-.125*RP*SP
            HT(4)=-.125*RM*SP
            HT(5)= .125*RM*SM
            HT(6)= .125*RP*SM
            HT(7)= .125*RP*SP
            HT(8)= .125*RM*SP
!C
!C*JACOBI MATRIX 
            XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4) &
     &          +HR(5)*XX(5)+HR(6)*XX(6)+HR(7)*XX(7)+HR(8)*XX(8)
            XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4) &
     &          +HS(5)*XX(5)+HS(6)*XX(6)+HS(7)*XX(7)+HS(8)*XX(8)
            XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4) &
     &          +HT(5)*XX(5)+HT(6)*XX(6)+HT(7)*XX(7)+HT(8)*XX(8)
!C
            XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4) &
     &          +HR(5)*YY(5)+HR(6)*YY(6)+HR(7)*YY(7)+HR(8)*YY(8)
            XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4) &
     &          +HS(5)*YY(5)+HS(6)*YY(6)+HS(7)*YY(7)+HS(8)*YY(8)
            XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4) &
     &          +HT(5)*YY(5)+HT(6)*YY(6)+HT(7)*YY(7)+HT(8)*YY(8)
!C
            XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4) &
     &          +HR(5)*ZZ(5)+HR(6)*ZZ(6)+HR(7)*ZZ(7)+HR(8)*ZZ(8)
            XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4) &
     &          +HS(5)*ZZ(5)+HS(6)*ZZ(6)+HS(7)*ZZ(7)+HS(8)*ZZ(8)
            XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4) &
     &          +HT(5)*ZZ(5)+HT(6)*ZZ(6)+HT(7)*ZZ(7)+HT(8)*ZZ(8)
!C  
!C*DETERMINANT OF JACOBIAN
!C 
            DET=XJ11*XJ22*XJ33 &
     &         +XJ12*XJ23*XJ31 &
     &         +XJ13*XJ21*XJ32 &
     &         -XJ13*XJ22*XJ31 &
     &         -XJ12*XJ21*XJ33 &
     &         -XJ11*XJ23*XJ32
!C
!C* INVERSION OF JACOBIAN
!C
            DUM= -1.0/DET
!C
            XJI11=DUM*( XJ22*XJ33-XJ23*XJ32)
            XJI21=DUM*(-XJ21*XJ33+XJ23*XJ31)
            XJI31=DUM*( XJ21*XJ32-XJ22*XJ31)
            XJI12=DUM*(-XJ12*XJ33+XJ13*XJ32)
            XJI22=DUM*( XJ11*XJ33-XJ13*XJ31)
            XJI32=DUM*(-XJ11*XJ32+XJ12*XJ31)
            XJI13=DUM*( XJ12*XJ23-XJ13*XJ22)
            XJI23=DUM*(-XJ11*XJ23+XJ13*XJ21)
            XJI33=DUM*( XJ11*XJ22-XJ12*XJ21)
!C
            DO J=1, 8
              BX(LX,LY,LZ,J)=XJI11*HR(J)+XJI12*HS(J)+XJI13*HT(J)
              BY(LX,LY,LZ,J)=XJI21*HR(J)+XJI22*HS(J)+XJI23*HT(J)
              BZ(LX,LY,LZ,J)=XJI31*HR(J)+XJI32*HS(J)+XJI33*HT(J)
           end DO
!C
!C
        end DO
     end DO
  end DO

  J = 1
  do LX = 1, 2
     do LY = 1, 2
        do LZ = 1, 2

           N(J,1) = H(LX,LY,LZ,1)
           N(J,2) = H(LX,LY,LZ,2)
           N(J,3) = H(LX,LY,LZ,3)
           N(J,4) = H(LX,LY,LZ,4)
           N(J,5) = H(LX,LY,LZ,5)
           N(J,6) = H(LX,LY,LZ,6)
           N(J,7) = H(LX,LY,LZ,7)
           N(J,8) = H(LX,LY,LZ,8)

           NX(J,1) = BX(LX,LY,LZ,1)
           NX(J,2) = BX(LX,LY,LZ,2)
           NX(J,3) = BX(LX,LY,LZ,3)
           NX(J,4) = BX(LX,LY,LZ,4)
           NX(J,5) = BX(LX,LY,LZ,5)
           NX(J,6) = BX(LX,LY,LZ,6)
           NX(J,7) = BX(LX,LY,LZ,7)
           NX(J,8) = BX(LX,LY,LZ,8)

           NY(J,1) = BY(LX,LY,LZ,1)
           NY(J,2) = BY(LX,LY,LZ,2)
           NY(J,3) = BY(LX,LY,LZ,3)
           NY(J,4) = BY(LX,LY,LZ,4)
           NY(J,5) = BY(LX,LY,LZ,5)
           NY(J,6) = BY(LX,LY,LZ,6)
           NY(J,7) = BY(LX,LY,LZ,7)
           NY(J,8) = BY(LX,LY,LZ,8)

           NZ(J,1) = BZ(LX,LY,LZ,1)
           NZ(J,2) = BZ(LX,LY,LZ,2)
           NZ(J,3) = BZ(LX,LY,LZ,3)
           NZ(J,4) = BZ(LX,LY,LZ,4)
           NZ(J,5) = BZ(LX,LY,LZ,5)
           NZ(J,6) = BZ(LX,LY,LZ,6)
           NZ(J,7) = BZ(LX,LY,LZ,7)
           NZ(J,8) = BZ(LX,LY,LZ,8)

           J = J + 1
        end do
     end do
  end do

!C
!C
end SUBROUTINE hecmw_Jacob_361

!*----------------------------------------------------------------------*
SUBROUTINE hecmw_Jacob_361_surface( hecMESH, iElem, iSurface, NOD, DET, H )
!*----------------------------------------------------------------------*
   use hecmw_util
   IMPLICIT NONE

   type(hecmwST_local_mesh):: hecMESH
   integer(kind=kint)::       iElem, iSurface, NOD(4)
   real(kind=kreal)::         DET

   real(kind=kreal):: XX(4),YY(4),ZZ(4)
   real(kind=kreal):: H(2,2,4),HR(4),HS(4),PL(4)
   real(kind=kreal):: XG(2),WGT(2)
   real(kind=kreal):: RI,SI,TI,RP,SP,TP,RM,SM,TM
   real(kind=kreal):: XJ11,XJ21,XJ31,XJ12,XJ22,XJ32,XJ13,XJ23,XJ33,WG
   integer(kind=kint):: iNode
   integer(kind=kint):: IG1,IG2,LX,LY,LZ,I
   real(kind=kreal):: VX,VY,VZ
   real(kind=kreal):: G1X,G1Y,G1Z
   real(kind=kreal):: G2X,G2Y,G2Z
   real(kind=kreal):: G3X,G3Y,G3Z

!**************************
!*  GAUSS INTEGRATION POINT
!**************************
   DATA XG/-0.5773502691896D0,0.5773502691896D0/
   DATA WGT/1.0D0,1.0D0/
!*
!* Set node numbers for selected surface
!*
   IF( iSurface.EQ.1 ) THEN
      NOD(1) = 1
      NOD(2) = 2
      NOD(3) = 3
      NOD(4) = 4
   ELSE IF( iSurface.EQ.2 ) THEN
      NOD(1) = 8
      NOD(2) = 7
      NOD(3) = 6
      NOD(4) = 5
   ELSE IF( iSurface.EQ.3 ) THEN
      NOD(1) = 5
      NOD(2) = 6
      NOD(3) = 2
      NOD(4) = 1
   ELSE IF( iSurface.EQ.4 ) THEN
      NOD(1) = 6
      NOD(2) = 7
      NOD(3) = 3
      NOD(4) = 2
   ELSE IF( iSurface.EQ.5 ) THEN
      NOD(1) = 7
      NOD(2) = 8
      NOD(3) = 4
      NOD(4) = 3
   ELSE IF( iSurface.EQ.6 ) THEN
      NOD(1) = 8
      NOD(2) = 5
      NOD(3) = 1
      NOD(4) = 4
   ENDIF

   do i = 1, 4
      iNode = hecMESH%elem_node_item( 8*(iElem-1) + NOD(i) )
      XX(i) = hecMESH%node( 3*iNode -2 )
      YY(i) = hecMESH%node( 3*iNode -1 )
      ZZ(i) = hecMESH%node( 3*iNode    )
   end do

!*
!*** SURFACE LOAD
!*
!* INTEGRATION OVER SURFACE
   DO IG2=1,2
      SI=XG(IG2)
      DO IG1=1,2
         RI=XG(IG1)
         H(IG1,IG2,1)=0.25*(1.0-RI)*(1.0-SI)
         H(IG1,IG2,2)=0.25*(1.0+RI)*(1.0-SI)
         H(IG1,IG2,3)=0.25*(1.0+RI)*(1.0+SI)
         H(IG1,IG2,4)=0.25*(1.0-RI)*(1.0+SI)
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
         DO I=1,4
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
         DET=DSQRT(G3X**2+G3Y**2+G3Z**2)
      ENDDO
   ENDDO

END SUBROUTINE hecmw_Jacob_361_surface

end module hecmw_Jacob361
