!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines for calculating distributed heat flux
!! for various elements
module m_heat_LIB_DFLUX
contains
  !C************************************************************************
  !C*  DFLUX_111(NN,XX,YY,ZZ,ASECT,LTYPE,VAL,VECT)
  !C*  DFLUX_231(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*  DFLUX_232(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*  DFLUX_241(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*  DFLUX_242(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*  DFLUX_341(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_342(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_351(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_352(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_361(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_362(NN,XX,YY,ZZ,LTYPE,VAL,VECT)
  !C*  DFLUX_731(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*  DFLUX_741(NN,XX,YY,ZZ,THICK,LTYPE,VAL,VECT)
  !C*--------------------------------------------------------------------*
  subroutine heat_DFLUX_111(NN,XX,YY,ZZ,ASECT,LTYPE,val,VECT)
    !C*--------------------------------------------------------------------*
    !C***
    !C***  SET 111 DFLUX
    !C***
    !C*   LTYPE=0  : BODY FLUX
    use hecmw
    implicit none
    !C* I/F VARIABLES
    integer(kind=kint) NN, LTYPE
    real(kind=kreal) XX(NN),YY(NN),ZZ(NN),ASECT,val,VECT(NN)
    !C* LOCAL VARIABLES
    integer(kind=kint) I
    real(kind=kreal) DX, DY, DZ, AL, VV
    !C*
    !C*
    if( LTYPE.EQ.0 ) then
      ASECT = 1.0d0 * ASECT
    else
      ASECT = 0.0d0
    endif

    do I=1,NN
      VECT(I)=0.0
    enddo
    !C*
    DX = XX(2) - XX(1)
    DY = YY(2) - YY(1)
    DZ = ZZ(2) - ZZ(1)
    AL = dsqrt( DX*DX + DY*DY + DZ*DZ )
    VV = ASECT * AL
    !
    VECT(1)= -val*VV/2.0d0
    VECT(2)= -val*VV/2.0d0
    !
    return

  end subroutine heat_DFLUX_111
  !C*--------------------------------------------------------------------*
  subroutine heat_DFLUX_231(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !C*--------------------------------------------------------------------*
    !C***
    !C***  SET 231 DFLUX
    !C***
    !C*   LTYPE=0  : BODY FLUX
    use hecmw
    implicit none
    !C* I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    !C* LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, LX
    real(kind=kreal)   :: RI, GX, GY, XSUM
    real(kind=kreal)   :: V1X, V1Y, V1Z, V2X, V2Y, V2Z, V3X, V3Y, V3Z, AA, VV
    real(kind=kreal)   :: XG(2), WGT(2), H(4), HR(4), HS(4)
    integer(kind=kint) :: NOD(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !C*
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else if( LTYPE.GE.1 ) then
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
      else if(LTYPE.EQ.2) then
        NOD(1)=2
        NOD(2)=3
      else if(LTYPE.EQ.3) then
        NOD(1)=3
        NOD(2)=1
      endif
    endif
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      do I=1,NN
        VECT(I)=0.0
      enddo
      !** INTEGRATION OVER SURFACE
      do LX=1,2
        RI=XG(LX)
        H(1)=0.5*(1.0-RI)
        H(2)=0.5*(1.0+RI)
        HR(1)=-0.5
        HR(2)= 0.5
        GX=0.0
        GY=0.0
        do I=1,2
          GX=GX+HR(I)*XX(NOD(I))
          GY=GY+HR(I)*YY(NOD(I))
        enddo
        XSUM=dsqrt(GX*GX+GY*GY)*THICK
        do I=1,2
          VECT(NOD(I))=VECT(NOD(I))-XSUM*WGT(LX)*H(I)*val
        enddo
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then
      do I=1,NN
        VECT(I)=0.0
      enddo

      V1X=XX(2)-XX(1)
      V1Y=YY(2)-YY(1)
      V1Z=ZZ(2)-ZZ(1)
      V2X=XX(3)-XX(1)
      V2Y=YY(3)-YY(1)
      V2Z=ZZ(3)-ZZ(1)
      V3X= V1Y*V2Z-V1Z*V2Y
      V3Y= V1Z*V2X-V1X*V2Z
      V3Z= V1X*V2Y-V1Y*V2X

      AA=0.5*dsqrt( V3X*v3X + V3Y*V3Y + V3Z*V3Z )
      VV=AA*THICK
      VECT(1)= -val*VV/3.0
      VECT(2)= -val*VV/3.0
      VECT(3)= -val*VV/3.0

    endif
    return

  end subroutine heat_DFLUX_231
  !C*--------------------------------------------------------------------*
  subroutine heat_DFLUX_232(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !C*--------------------------------------------------------------------*
    !C***
    !C***  SET 232 DFLUX
    !C***
    !C*   LTYPE=0  : BODY FLUX
    use hecmw
    implicit none
    !C* I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    !C* LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, L1, L2, LX
    real(kind=kreal)   :: RI, GX, GY, XSUM, X1, X2, X3, XL1, XL2
    real(kind=kreal)   :: XJ11, XJ12, XJ21, XJ22, DET, WG
    real(kind=kreal)   :: XG(3), WGT(3), H(6), HR(3)
    real(kind=kreal)   :: HL1(6), HL2(6), HL3(6)
    integer(kind=kint) :: NOD(3)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C*
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else if( LTYPE.GE.1 ) then
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
        NOD(3)=4
      else if(LTYPE.EQ.2) then
        NOD(1)=2
        NOD(2)=3
        NOD(3)=5
      else if(LTYPE.EQ.3) then
        NOD(1)=3
        NOD(2)=1
        NOD(3)=6
      endif
    endif
    do I=1,NN
      VECT(I)=0.0
    enddo
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      !** INTEGRATION OVER SURFACE
      do LX=1,3
        RI=XG(LX)
        H(1)=-0.5*(1.0-RI)*RI
        H(2)= 0.5*(1.0+RI)*RI
        H(3)=1.0-RI*RI
        HR(1)= RI-0.5
        HR(2)= RI+0.5
        HR(3)=-2.0*RI
        GX=0.0
        GY=0.0
        do I=1,3
          GX=GX+HR(I)*XX(NOD(I))
          GY=GY+HR(I)*YY(NOD(I))
        enddo
        XSUM=dsqrt(GX*GX+GY*GY)
        WG=WGT(LX)*val*THICK
        VECT(NOD(1))=VECT(NOD(1))-H(1)*XSUM*WG
        VECT(NOD(2))=VECT(NOD(2))-H(2)*XSUM*WG
        VECT(NOD(3))=VECT(NOD(3))-H(3)*XSUM*WG
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then

      !* LOOP OVER ALL INTEGRATION POINTS
      do L2=1,3
        XL2=XG(L2)
        X2 =(XL2+1.0)*0.5
        do L1=1,3
          XL1=XG(L1)
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
          do I=1,NN
            XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
            XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
            XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
            XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
          enddo
          DET=XJ11*XJ22-XJ21*XJ12
          WG=WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
          do I = 1, NN
            VECT(I) = VECT(I) - H(I)*WG*val
          enddo
        enddo
      enddo
    endif
    return

  end subroutine heat_DFLUX_232
  !C*--------------------------------------------------------------------*
  subroutine heat_DFLUX_241(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !C*--------------------------------------------------------------------*
    !C***
    !C***  SET 241 DFLUX
    !C***
    !C*   LTYPE=0  : BODY FLUX
    use hecmw
    implicit none
    !C* I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    !C* LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, LX, LY
    real(kind=kreal)   :: GX, GY, XSUM
    real(kind=kreal)   :: RI, SI, RP, SP, RM, SM
    real(kind=kreal)   :: XJ11, XJ12, XJ21, XJ22, DET, WG
    real(kind=kreal)   :: XG(2), WGT(2), H(4), HR(4), HS(4)
    integer(kind=kint) :: NOD(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !C*
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else if( LTYPE.GE.1 ) then
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
      else if(LTYPE.EQ.2) then
        NOD(1)=2
        NOD(2)=3
      else if(LTYPE.EQ.3) then
        NOD(1)=3
        NOD(2)=4
      else if(LTYPE.EQ.4) then
        NOD(1)=4
        NOD(2)=1
      endif
    endif
    do I=1,NN
      VECT(I)=0.0
    enddo
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      !** INTEGRATION OVER SURFACE
      do LX=1,2
        RI=XG(LX)
        H(1)=0.5*(1.0-RI)
        H(2)=0.5*(1.0+RI)
        HR(1)=-0.5
        HR(2)= 0.5
        GX=0.0
        GY=0.0
        do I=1,2
          GX=GX+HR(I)*XX(NOD(I))
          GY=GY+HR(I)*YY(NOD(I))
        enddo
        XSUM=dsqrt(GX*GX+GY*GY)*THICK
        do I=1,2
          VECT(NOD(I))=VECT(NOD(I))-XSUM*WGT(LX)*H(I)*val
        enddo
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then
      do LX=1,2
        RI=XG(LX)
        do LY=1,2
          SI=XG(LY)
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
          do I=1,NN
            XJ11=XJ11+HR(I)*XX(I)
            XJ21=XJ21+HS(I)*XX(I)
            XJ12=XJ12+HR(I)*YY(I)
            XJ22=XJ22+HS(I)*YY(I)
          enddo
          DET=XJ11*XJ22-XJ21*XJ12
          WG=WGT(LX)*WGT(LY)*DET*THICK
          do I =1,NN
            VECT(I)=VECT(I)-WG*H(I)*val
          enddo
        enddo
      enddo
    endif
    return

  end subroutine heat_DFLUX_241
  !C*--------------------------------------------------------------------*
  subroutine heat_DFLUX_242(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !C*--------------------------------------------------------------------*
    !C***
    !C***  SET 242 DFLUX
    !C***
    !C*   LTYPE=0  : BODY FLUX
    use hecmw
    implicit none
    !C* I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    !C* LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, LX, LY
    real(kind=kreal)   :: GX, GY, XSUM
    real(kind=kreal)   :: RI, SI, RP, SP, RM, SM
    real(kind=kreal)   :: XJ11, XJ12, XJ21, XJ22, DET, WG
    real(kind=kreal)   :: XG(3), WGT(3), H(8), HR(8), HS(8)
    integer(kind=kint) :: NOD(3)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C*
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else if( LTYPE.GE.1 ) then
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
        NOD(3)=5
      else if(LTYPE.EQ.2) then
        NOD(1)=2
        NOD(2)=3
        NOD(3)=6
      else if(LTYPE.EQ.3) then
        NOD(1)=3
        NOD(2)=4
        NOD(3)=7
      else if(LTYPE.EQ.4) then
        NOD(1)=4
        NOD(2)=1
        NOD(3)=8
      endif
    endif
    do I=1,NN
      VECT(I)=0.0
    enddo
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      !** INTEGRATION OVER SURFACE
      do LX=1,3
        RI=XG(LX)
        H(1)=-0.5*(1.0-RI)*RI
        H(2)= 0.5*(1.0+RI)*RI
        H(3)=1.0-RI*RI
        HR(1)= RI-0.5
        HR(2)= RI+0.5
        HR(3)=-2.0*RI
        GX=0.0
        GY=0.0
        do I=1,3
          GX=GX+HR(I)*XX(NOD(I))
          GY=GY+HR(I)*YY(NOD(I))
        enddo
        XSUM=dsqrt(GX*GX+GY*GY)*THICK
        VECT(NOD(1))=VECT(NOD(1))-H(1)*WGT(LX)*val*XSUM
        VECT(NOD(2))=VECT(NOD(2))-H(2)*WGT(LX)*val*XSUM
        VECT(NOD(3))=VECT(NOD(3))-H(3)*WGT(LX)*val*XSUM
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then
      do LX=1,3
        RI=XG(LX)
        do LY=1,3
          SI=XG(LY)
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
          do I=1,NN
            XJ11=XJ11+HR(I)*XX(I)
            XJ21=XJ21+HS(I)*XX(I)
            XJ12=XJ12+HR(I)*YY(I)
            XJ22=XJ22+HS(I)*YY(I)
          enddo
          DET=XJ11*XJ22-XJ21*XJ12
          WG=WGT(LX)*WGT(LY)*DET*THICK
          do I =1,NN
            VECT(I)=VECT(I)-WG*H(I)*val
          enddo
        enddo
      enddo
    endif
    return

  end subroutine heat_DFLUX_242
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_341(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 341 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), VECT(NN)
    ! LOCAL VARIABLES
    real(kind=kreal)   :: H(4), HR(4), HS(4), HT(4)
    real(kind=kreal)   :: XG(2), WGT(2)
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33, DET, WG
    integer(kind=kint) :: IVOL, ISUF
    integer(kind=kint) :: NOD(4)
    integer(kind=kint) :: L1, L2, L3, I
    real(kind=kreal)   :: val, AA
    real(kind=kreal)   :: V1X, V1Y, V1Z
    real(kind=kreal)   :: V2X, V2Y, V2Z
    real(kind=kreal)   :: V3X, V3Y, V3Z
    real(kind=kreal)   :: XL1, XL2, XL3
    real(kind=kreal)   :: X1, X2, X3
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !
    IVOL=0
    ISUF=0
    NOD = 0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
        NOD(3)=3
      else if(LTYPE.EQ.2) then
        NOD(1)=4
        NOD(2)=2
        NOD(3)=1
      else if(LTYPE.EQ.3) then
        NOD(1)=4
        NOD(2)=3
        NOD(3)=2
      else if(LTYPE.EQ.4) then
        NOD(1)=4
        NOD(2)=1
        NOD(3)=3
      endif
    endif
    ! CLEAR VECT
    do I=1,NN
      VECT(I)=0.0
    enddo
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      V1X=XX(NOD(2))-XX(NOD(1))
      V1Y=YY(NOD(2))-YY(NOD(1))
      V1Z=ZZ(NOD(2))-ZZ(NOD(1))
      V2X=XX(NOD(3))-XX(NOD(1))
      V2Y=YY(NOD(3))-YY(NOD(1))
      V2Z=ZZ(NOD(3))-ZZ(NOD(1))
      V3X= V1Y*V2Z-V1Z*V2Y
      V3Y= V1Z*V2X-V1X*V2Z
      V3Z= V1X*V2Y-V1Y*V2X
      AA=0.5*dsqrt( V3X*v3X + V3Y*V3Y + V3Z*V3Z )
      VECT(NOD(1))= -val*AA/3.0
      VECT(NOD(2))= -val*AA/3.0
      VECT(NOD(3))= -val*AA/3.0
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then
      ! LOOP FOR INTEGRATION POINTS
      do L3=1,2
        XL3=XG(L3)
        X3 =(XL3+1.0)*0.5
        do L2=1,2
          XL2=XG(L2)
          X2 =(1.0-X3)*(XL2+1.0)*0.5
          do L1=1,2
            XL1=XG(L1)
            X1=(1.0-X2-X3)*(XL1+1.0)*0.5
            !  INTERPOLATION FUNCTION
            H(1)=X1
            H(2)=X2
            H(3)=X3
            H(4)=1.0-X1-X2-X3
            ! DERIVATIVE OF INTERPOLATION FUNCTION
            ! FOR L1-COORDINATE
            HR(1)= 1.0
            HR(2)= 0.0
            HR(3)= 0.0
            HR(4)=-1.0
            !  FOR L2-COORDINATE
            HS(1)= 0.0
            HS(2)= 1.0
            HS(3)= 0.0
            HS(4)=-1.0
            !  FOR ZETA-COORDINATE
            HT(1)= 0.0
            HT(2)= 0.0
            HT(3)= 1.0
            HT(4)=-1.0
            !  JACOBI MATRIX
            XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
            XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
            XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4)
            XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
            XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
            XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4)
            XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
            XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)
            XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4)
            ! DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            !
            do I = 1, NN
              VECT(I) = VECT(I) + val*H(I)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
            enddo

          enddo
        enddo
      enddo
    endif
    !
    return

  end subroutine heat_DFLUX_341
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_342(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 342 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), VECT(NN)
    ! LOCAL VARIABLES
    real(kind=kreal)   :: H(10)
    real(kind=kreal)   :: HL1(10), HL2(10), HL3(10), HL4(10)
    real(kind=kreal)   :: XG(3), WGT(3)
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: DET, WG
    integer(kind=kint) :: IVOL, ISUF
    integer(kind=kint) :: NOD(10)
    integer(kind=kint) :: LX, LY, LZ, I, IG1, IG2, L1, L2, L3
    real(kind=kreal)   :: val, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z
    real(kind=kreal)   :: G2X, G2Y, G2Z
    real(kind=kreal)   :: G3X, G3Y, G3Z
    real(kind=kreal)   :: XL1, XL2, XL3
    real(kind=kreal)   :: X1, X2, X3, X4
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if(LTYPE.EQ.1) then
        NOD(1)=1
        NOD(2)=2
        NOD(3)=3
        NOD(4)=5
        NOD(5)=6
        NOD(6)=7
      else if(LTYPE.EQ.2) then
        NOD(1)=4
        NOD(2)=2
        NOD(3)=1
        NOD(4)=9
        NOD(5)=5
        NOD(6)=8
      else if(LTYPE.EQ.3) then
        NOD(1)=4
        NOD(2)=3
        NOD(3)=2
        NOD(4)=10
        NOD(5)=6
        NOD(6)=9
      else if(LTYPE.EQ.4) then
        NOD(1)=4
        NOD(2)=1
        NOD(3)=3
        NOD(4)=8
        NOD(5)=7
        NOD(6)=10
      endif
    endif
    ! CLEAR VECT
    do I=1,10
      VECT(I)=0.0
    enddo
    !** SURFACE LOAD
    if( ISUF.EQ.1 ) then
      ! INTEGRATION OVER SURFACE
      do L2=1,3
        XL2=XG(L2)
        X2 =(XL2+1.0)*0.5
        do L1=1,3
          XL1=XG(L1)
          X1=0.5*(1.0-X2)*(XL1+1.0)
          ! INTERPOLATION FUNCTION
          X3=1.0-X1-X2
          H(1)= X1*(2.0*X1-1.)
          H(2)= X2*(2.0*X2-1.)
          H(3)= X3*(2.0*X3-1.)
          H(4)= 4.0*X1*X2
          H(5)= 4.0*X2*X3
          H(6)= 4.0*X1*X3
          ! DERIVATIVE OF INTERPOLATION FUNCTION
          ! FOR L1-COORDINATE
          HL1(1)=4.0*X1-1.0
          HL1(2)= 0.0
          HL1(3)= 0.0
          HL1(4)= 4.0*X2
          HL1(5)= 0.0
          HL1(6)= 4.0*X3
          ! FOR L2-COORDINATE
          HL2(1)= 0.0
          HL2(2)= 4.0*X2-1.0
          HL2(3)= 0.0
          HL2(4)= 4.0*X1
          HL2(5)= 4.0*X3
          HL2(6)= 0.0
          ! FOR L3-COORDINATE
          HL3(1)= 0.0
          HL3(2)= 0.0
          HL3(3)= 4.0*X3-1.0
          HL3(4)= 0.0
          HL3(5)= 4.0*X2
          HL3(6)= 4.0*X1
          ! JACOBI MATRIX
          G1X=0.0
          G1Y=0.0
          G1Z=0.0
          G2X=0.0
          G2Y=0.0
          G2Z=0.0
          do I=1,6
            G1X=G1X+(HL1(I)-HL3(I))*XX(NOD(I))
            G1Y=G1Y+(HL1(I)-HL3(I))*YY(NOD(I))
            G1Z=G1Z+(HL1(I)-HL3(I))*ZZ(NOD(I))
            G2X=G2X+(HL2(I)-HL3(I))*XX(NOD(I))
            G2Y=G2Y+(HL2(I)-HL3(I))*YY(NOD(I))
            G2Z=G2Z+(HL2(I)-HL3(I))*ZZ(NOD(I))
          enddo
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          G3X=G3X/XSUM
          G3Y=G3Y/XSUM
          G3Z=G3Z/XSUM
          ! JACOBI MATRIX
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
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          WG=WGT(L1)*WGT(L2)*DET*(1.0-X2)*0.25
          do I=1,6
            VECT(NOD(I))=VECT(NOD(I))-val*WG*H(I)
          enddo
          !
        enddo
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL.EQ.1 ) then
      ! LOOP FOR INTEGRATION POINTS
      !DETERMINANT OF JACOBIAN
      do L3=1,3
        XL3=XG(L3)
        X3 =(XL3+1.0)*0.5
        do L2=1,3
          XL2=XG(L2)
          X2 =(1.0-X3)*(XL2+1.0)*0.5
          do L1=1,3
            XL1=XG(L1)
            X1=(1.0-X2-X3)*(XL1+1.0)*0.5
            !C* INTERPOLATION FUNCTION
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
            !C* DERIVATIVE OF INTERPOLATION FUNCTION
            !C* FOR L1-COORDINATE
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
            !C* FOR L2-COORDINATE
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
            !C* FOR L3-COORDINATE
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
            !C* FOR L4-COORDINATE
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
            !C* JACOBI MATRIX
            XJ11=0.0
            XJ21=0.0
            XJ31=0.0
            XJ12=0.0
            XJ22=0.0
            XJ32=0.0
            XJ13=0.0
            XJ23=0.0
            XJ33=0.0
            do I=1,NN
              XJ11=XJ11+(HL1(I)-HL4(I))*XX(I)
              XJ21=XJ21+(HL2(I)-HL4(I))*XX(I)
              XJ31=XJ31+(HL3(I)-HL4(I))*XX(I)
              XJ12=XJ12+(HL1(I)-HL4(I))*YY(I)
              XJ22=XJ22+(HL2(I)-HL4(I))*YY(I)
              XJ32=XJ32+(HL3(I)-HL4(I))*YY(I)
              XJ13=XJ13+(HL1(I)-HL4(I))*ZZ(I)
              XJ23=XJ23+(HL2(I)-HL4(I))*ZZ(I)
              XJ33=XJ33+(HL3(I)-HL4(I))*ZZ(I)
            enddo
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            DET=-DET
            WG=DET*WGT(L1)*WGT(L2)*WGT(L3)*(1.-X3)*(1.0-X2-X3)*0.125
            do I = 1, NN
              VECT(I) = VECT(I) - val*WG*H(I)
            enddo
            !
          enddo
        enddo
      enddo
    endif
    !
    return

  end subroutine heat_DFLUX_342
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_351(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 351 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    !   S5   LTYPE=5  :FLUX IN NORMAL-DIRECTION FOR FACE-5
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), val, VECT(NN)
    ! LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, IG1, IG2, L12, LZ
    real(kind=kreal)   :: RI, SI, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: V1X, V1Y, V1Z, V2X, V2Y, V2Z, V3X, V3Y, V3Z, AA
    real(kind=kreal)   :: ZI, X1, X2
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: DET, WG
    real(kind=kreal)   :: H(6), HR(6), HS(6), HT(6), PL(6)
    real(kind=kreal)   :: XG(2), WGT(2), XG1(3), XG2(3), WGT1(3)
    integer(kind=kint) :: NOD(4)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG  / -0.5773502691896,0.5773502691896 /
    data WGT /  1.0,            1.0 /
    data XG1 /  0.6666666667,   0.16666666667, 0.16666666667 /
    data XG2 /  0.1666666667,   0.66666666667, 0.16666666667 /
    data WGT1/  0.1666666667,   0.16666666667, 0.16666666667 /
    !
    ! SELECTION OF LOAD TYPE
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if( LTYPE.EQ.1 ) then
        NOD(1) = 1
        NOD(2) = 2
        NOD(3) = 3
        ISUF=2
      else if( LTYPE.EQ.2 ) then
        NOD(1) = 6
        NOD(2) = 5
        NOD(3) = 4
        ISUF=2
      else if( LTYPE.EQ.3 ) then
        NOD(1) = 4
        NOD(2) = 5
        NOD(3) = 2
        NOD(4) = 1
      else if( LTYPE.EQ.4 ) then
        NOD(1) = 5
        NOD(2) = 6
        NOD(3) = 3
        NOD(4) = 2
      else if( LTYPE.EQ.5 ) then
        NOD(1) = 6
        NOD(2) = 4
        NOD(3) = 1
        NOD(4) = 3
      endif
    endif
    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    !** SURFACE LOAD
    !
    if( ISUF.EQ.1 ) then
      do IG2=1,2
        SI=XG(IG2)
        do IG1=1,2
          RI=XG(IG1)
          !
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
          !
          G1X=0.0
          G1Y=0.0
          G1Z=0.0
          G2X=0.0
          G2Y=0.0
          G2Z=0.0
          !
          do I=1,4
            G1X=G1X+HR(I)*XX(NOD(I))
            G1Y=G1Y+HR(I)*YY(NOD(I))
            G1Z=G1Z+HR(I)*ZZ(NOD(I))
            G2X=G2X+HS(I)*XX(NOD(I))
            G2Y=G2Y+HS(I)*YY(NOD(I))
            G2Z=G2Z+HS(I)*ZZ(NOD(I))
          enddo
          !
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          do I=1,4
            VECT(NOD(I))=VECT(NOD(I))-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)
          enddo
          !
        enddo
      enddo
    elseif( ISUF.EQ.2 ) then
      V1X=XX(NOD(2))-XX(NOD(1))
      V1Y=YY(NOD(2))-YY(NOD(1))
      V1Z=ZZ(NOD(2))-ZZ(NOD(1))
      V2X=XX(NOD(3))-XX(NOD(1))
      V2Y=YY(NOD(3))-YY(NOD(1))
      V2Z=ZZ(NOD(3))-ZZ(NOD(1))
      V3X= V1Y*V2Z-V1Z*V2Y
      V3Y= V1Z*V2X-V1X*V2Z
      V3Z= V1X*V2Y-V1Y*V2X
      AA=0.5*dsqrt( V3X*v3X + V3Y*V3Y + V3Z*V3Z )
      VECT(NOD(1))= -val*AA/3.0
      VECT(NOD(2))= -val*AA/3.0
      VECT(NOD(3))= -val*AA/3.0
    endif
    !
    !** VOLUME LOAD
    !
    if( IVOL.EQ.1 ) then
      do I=1,6
        VECT(I)=0.0
      enddo
      !
      do LZ=1,2
        ZI=XG(LZ)
        do L12=1,3
          X1=XG1(L12)
          X2=XG2(L12)

          H(1)=X1*(1.0-ZI)*0.5
          H(2)=X2*(1.0-ZI)*0.5
          H(3)=(1.0-X1-X2)*(1.0-ZI)*0.5
          H(4)=X1*(1.0+ZI)*0.5
          H(5)=X2*(1.0+ZI)*0.5
          H(6)=(1.0-X1-X2)*(1.0+ZI)*0.5
          HR(1)= (1.0-ZI)*0.5
          HR(2)= 0.0
          HR(3)=-(1.0-ZI)*0.5
          HR(4)= (1.0+ZI)*0.5
          HR(5)= 0.0
          HR(6)=-(1.0+ZI)*0.5
          HS(1)=  0.0
          HS(2)=  (1.0-ZI)*0.5
          HS(3)= -(1.0-ZI)*0.5
          HS(4)= 0.0
          HS(5)=  (1.0+ZI)*0.5
          HS(6)= -(1.0+ZI)*0.5
          HT(1)= -0.5*X1
          HT(2)= -0.5*X2
          HT(3)= -0.5*(1.0-X1-X2)
          HT(4)=  0.5*X1
          HT(5)=  0.5*X2
          HT(6)=  0.5*(1.0-X1-X2)
          !
          XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4) &
            +HR(5)*XX(5)+HR(6)*XX(6)
          XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4) &
            +HS(5)*XX(5)+HS(6)*XX(6)
          XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4) &
            +HT(5)*XX(5)+HT(6)*XX(6)
          XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4) &
            +HR(5)*YY(5)+HR(6)*YY(6)
          XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4) &
            +HS(5)*YY(5)+HS(6)*YY(6)
          XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4) &
            +HT(5)*YY(5)+HT(6)*YY(6)
          XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4) &
            +HR(5)*ZZ(5)+HR(6)*ZZ(6)
          XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4) &
            +HS(5)*ZZ(5)+HS(6)*ZZ(6)
          XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4) &
            +HT(5)*ZZ(5)+HT(6)*ZZ(6)
          !*DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          WG=WGT1(L12)*WGT(LZ)*DET
          do I=1,NN
            VECT(I)=VECT(I)-val*H(I)*WG
          enddo
          !
        enddo
      enddo
    endif
    !
    return

  end subroutine heat_DFLUX_351
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_352(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 352 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    !   S5   LTYPE=5  :FLUX IN NORMAL-DIRECTION FOR FACE-5
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), val, VECT(NN)
    ! LOCAL VARIABLES
    integer(kind=kint) :: IVOL, ISUF, I, IG1, IG2, LX, LY, LZ
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: DET, WG
    real(kind=kreal)   :: H(15), HR(15), HS(15), HT(15), PL(15)
    real(kind=kreal)   :: XG(3), WGT(3)
    integer(kind=kint) :: NOD(8)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !
    ! SELECTION OF LOAD TYPE
    !
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if( LTYPE.EQ.1 ) then
        NOD(1) = 1
        NOD(2) = 2
        NOD(3) = 3
        NOD(4) = 7
        NOD(5) = 8
        NOD(6) = 9
        ISUF=2
      else if( LTYPE.EQ.2 ) then
        NOD(1) = 6
        NOD(2) = 5
        NOD(3) = 4
        NOD(4) = 11
        NOD(5) = 10
        NOD(6) = 12
        ISUF=2
      else if( LTYPE.EQ.3 ) then
        NOD(1) = 4
        NOD(2) = 5
        NOD(3) = 2
        NOD(4) = 1
        NOD(5) = 10
        NOD(6) = 14
        NOD(7) = 7
        NOD(8) = 13
      else if( LTYPE.EQ.4 ) then
        NOD(1) = 5
        NOD(2) = 6
        NOD(3) = 3
        NOD(4) = 2
        NOD(5) = 11
        NOD(6) = 15
        NOD(7) = 8
        NOD(8) = 14
      else if( LTYPE.EQ.5 ) then
        NOD(1) = 6
        NOD(2) = 4
        NOD(3) = 1
        NOD(4) = 3
        NOD(5) = 12
        NOD(6) = 13
        NOD(7) = 9
        NOD(8) = 15
      endif
    endif
    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    !** SURFACE LOAD
    !
    if( ISUF.EQ.1 ) then
      do IG2=1,3
        SI=XG(IG2)
        do IG1=1,3
          RI=XG(IG1)
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
          HR(6)= .5*(1.0-SI*SI)
          HR(7)=-RI*(1.0+SI)
          HR(8)=-.5*(1.0-SI*SI)
          HS(1)=-.25*RM*(-1.0-RI-SI)-0.25*RM*SM
          HS(2)=-.25*RP*(-1.0+RI-SI)-0.25*RP*SM
          HS(3)= .25*RP*(-1.0+RI+SI)+0.25*RP*SP
          HS(4)= .25*RM*(-1.0-RI+SI)+0.25*RM*SP
          HS(5)=-.5*(1.0-RI*RI)
          HS(6)=-SI*(1.0+RI)
          HS(7)= .5*(1.0-RI*RI)
          HS(8)=-SI*(1.0-RI)
          G1X=0.0
          G1Y=0.0
          G1Z=0.0
          G2X=0.0
          G2Y=0.0
          G2Z=0.0
          G3X=0.0
          G3Y=0.0
          G3Z=0.0
          do I=1,8
            G1X=G1X+HR(I)*XX(NOD(I))
            G1Y=G1Y+HR(I)*YY(NOD(I))
            G1Z=G1Z+HR(I)*ZZ(NOD(I))
            G2X=G2X+HS(I)*XX(NOD(I))
            G2Y=G2Y+HS(I)*YY(NOD(I))
            G2Z=G2Z+HS(I)*ZZ(NOD(I))
            G3X=G3X+HT(I)*XX(NOD(I))
            G3Y=G3Y+HT(I)*YY(NOD(I))
            G3Z=G3Z+HT(I)*ZZ(NOD(I))
          enddo
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          do I=1,8
            VECT(NOD(I))=VECT(NOD(I))-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)
          enddo
        enddo
      enddo
    elseif( ISUF.EQ.2 ) then
      do IG2=1,3
        SI=XG(IG2)
        do IG1=1,3
          RI=XG(IG1)
          RP=1.0+RI
          SP=1.0+SI
          RM=1.0-RI
          SM=1.0-SI
          H(1)=0.25*RM*SM*(-1.0-RI-SI)
          H(2)=0.25*RP*SM*(-1.0+RI-SI)
          H(3)=0.5*SP*SI
          H(4)=0.5*(1.0-RI*RI)*(1.0-SI)
          H(5)=0.5*(1.0-SI*SI)*(1.0+RI)
          H(6)=0.5*(1.0-SI*SI)*(1.0-RI)
          HR(1)=-.25*SM*(-1.0-RI-SI)-0.25*RM*SM
          HR(2)= .25*SM*(-1.0+RI-SI)+0.25*RP*SM
          HR(3)= 0.0
          HR(4)=-RI*(1.0-SI)
          HR(5)= 0.5*(1.0-SI*SI)
          HR(6)=-0.5*(1.0-SI*SI)
          HS(1)=-.25*RM*(-1.0-RI-SI)-0.25*RM*SM
          HS(2)=-.25*RP*(-1.0+RI-SI)-0.25*RP*SM
          HS(3)=0.5*(1.0+2.0*SI)
          HS(4)=-0.5*(1.0-RI*RI)
          HS(5)=-SI *(1.0+RI)
          HS(6)=-SI*(1.0-RI)
          G1X=0.0
          G1Y=0.0
          G1Z=0.0
          G2X=0.0
          G2Y=0.0
          G2Z=0.0
          G3X=0.0
          G3Y=0.0
          G3Z=0.0
          do I=1, 6
            G1X=G1X+HR(I)*XX(NOD(I))
            G1Y=G1Y+HR(I)*YY(NOD(I))
            G1Z=G1Z+HR(I)*ZZ(NOD(I))
            G2X=G2X+HS(I)*XX(NOD(I))
            G2Y=G2Y+HS(I)*YY(NOD(I))
            G2Z=G2Z+HS(I)*ZZ(NOD(I))
            G3X=G3X+HT(I)*XX(NOD(I))
            G3Y=G3Y+HT(I)*YY(NOD(I))
            G3Z=G3Z+HT(I)*ZZ(NOD(I))
          enddo
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          do I=1,6
            VECT(NOD(I))=VECT(NOD(I))-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)
          enddo
        enddo
      enddo
    endif
    !
    !** VOLUME LOAD
    !
    if( IVOL.EQ.1 ) then
      do LX=1,3
        RI=XG(LX)
        do LY=1,3
          SI=XG(LY)
          do LZ=1,3
            TI=XG(LZ)
            RP=1.0+RI
            SP=1.0+SI
            TP=1.0+TI
            RM=1.0-RI
            SM=1.0-SI
            TM=1.0-TI
            !  INTERPOLATION FUNCTION
            H(1)=-0.125*RM*SM*TM*(2.0+RI+SI+TI)
            H(2)=-0.125*RP*SM*TM*(2.0-RI+SI+TI)
            H(3)=-0.25*SP*TM*(1.0-SI+TI)
            H(4)=-0.125*RM*SM*TP*(2.0+RI+SI-TI)
            H(5)=-0.125*RP*SM*TP*(2.0-RI+SI-TI)
            H(6)=0.25*SP*TP*(SI+TI-1.0)
            H(7)=0.25*(1.0-RI**2)*SM*TM
            H(8)=0.25*RP*(1.0-SI**2)*TM
            H(9)=0.25*RM*(1.0-SI**2)*TM
            H(10)=0.25*(1.0-RI**2)*SM*TP
            H(11)=0.25*RP*(1.0-SI**2)*TP
            H(12)=0.25*RM*(1.0-SI**2)*TP
            H(13)=0.25*RM*SM*(1.0-TI**2)
            H(14)=0.25*RP*SM*(1.0-TI**2)
            H(15)=0.5*SP*(1.0-TI**2)
            !  DERIVATIVE OF INTERPOLATION FUNCTION
            !  FOR R-COORDINATE
            HR(1)=-0.125*RM*SM*TM+0.125*SM*TM*(2.0+RI+SI+TI)
            HR(2)=+0.125*RP*SM*TM-0.125*SM*TM*(2.0-RI+SI+TI)
            HR(3)=0.0
            HR(4)=-0.125*RM*SM*TP+0.125*SM*TP*(2.0+RI+SI-TI)
            HR(5)=+0.125*RP*SM*TP-0.125*SM*TP*(2.0-RI+SI-TI)
            HR(6)=0.0
            HR(7)=-0.50*RI*SM*TM
            HR(8)=+0.25*(1.0-SI**2)*TM
            HR(9)=-0.25*(1.0-SI**2)*TM
            HR(10)=-0.50*RI*SM*TP
            HR(11)=+0.25*(1.0-SI**2)*TP
            HR(12)=-0.25*(1.0-SI**2)*TP
            HR(13)=-0.25*SM*(1.0-TI**2)
            HR(14)=+0.25*SM*(1.0-TI**2)
            HR(15)=0.0
            !  FOR S-COORDINATE
            HS(1)=-0.125*RM*SM*TM+0.125*RM*TM*(2.0+RI+SI+TI)
            HS(2)=-0.125*RP*SM*TM+0.125*RP*TM*(2.0-RI+SI+TI)
            HS(3)=+0.25*TM*(2.0*SI-TI)
            HS(4)=-0.125*RM*SM*TP+0.125*RM*TP*(2.0+RI+SI-TI)
            HS(5)=-0.125*RP*SM*TP+0.125*RP*TP*(2.0-RI+SI-TI)
            HS(6)=+0.25*TP*(2.0*SI+TI)
            HS(7)=-0.25*(1.0-RI**2)*TM
            HS(8)=-0.50*RP*SI*TM
            HS(9)=-0.50*RM*SI*TM
            HS(10)=-0.25*(1.0-RI**2)*TP
            HS(11)=-0.50*RP*SI*TP
            HS(12)=-0.50*RM*SI*TP
            HS(13)=-0.25*RM*(1.0-TI**2)
            HS(14)=-0.25*RP*(1.0-TI**2)
            HS(15)=+0.50*(1.0-TI**2)
            !  FOR T-COORDINATE
            HT(1)=-0.125*RM*SM*TM+0.125*RM*SM*(2.0+RI+SI+TI)
            HT(2)=-0.125*RP*SM*TM+0.125*RP*SM*(2.0-RI+SI+TI)
            HT(3)=+0.25*SP*(2.0*TI-SI)
            HT(4)=+0.125*RM*SM*TP-0.125*RM*SM*(2.0+RI+SI-TI)
            HT(5)=+0.125*RP*SM*TP-0.125*RP*SM*(2.0-RI+SI-TI)
            HT(6)=+0.25*SP*(SI+2.0*TI)
            HT(7)=-0.25*(1.0-RI**2)*SM
            HT(8)=-0.25*RP*(1.0-SI**2)
            HT(9)=-0.25*RM*(1.0-SI**2)
            HT(10)=0.25*(1.0-RI**2)*SM
            HT(11)=0.25*RP*(1.0-SI**2)
            HT(12)=0.25*RM*(1.0-SI**2)
            HT(13)=-0.5*RM*SM*TI
            HT(14)=-0.5*RP*SM*TI
            HT(15)=-SP*TI
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
            do I=1,NN
              XJ11=XJ11+HR(I)*XX(I)
              XJ21=XJ21+HS(I)*XX(I)
              XJ31=XJ31+HT(I)*XX(I)
              XJ12=XJ12+HR(I)*YY(I)
              XJ22=XJ22+HS(I)*YY(I)
              XJ32=XJ32+HT(I)*YY(I)
              XJ13=XJ13+HR(I)*ZZ(I)
              XJ23=XJ23+HS(I)*ZZ(I)
              XJ33=XJ33+HT(I)*ZZ(I)
            enddo
            !DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            WG=DET*WGT(LX)*WGT(LY)*WGT(LZ)
            do I=1,NN
              VECT(I)=VECT(I) - val*WG*H(I)
            enddo
          enddo
        enddo
      enddo
    endif
    return

  end subroutine heat_DFLUX_352
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_361(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 361 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    !   S5   LTYPE=5  :FLUX IN NORMAL-DIRECTION FOR FACE-5
    !   S6   LTYPE=6  :FLUX IN NORMAL-DIRECTION FOR FACE-6
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), val, VECT(NN)
    ! LOCAL VARIABLES
    real(kind=kreal)   :: H(8), HR(8), HS(8), HT(8), PL(8)
    real(kind=kreal)   :: XG(2), WGT(2)
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33, DET, WG
    integer(kind=kint) :: IVOL, ISUF
    integer(kind=kint) :: NOD(4)
    integer(kind=kint) :: IG1, IG2, LX, LY, LZ, I
    real(kind=kreal)   :: VX, VY, VZ
    real(kind=kreal)   :: G1X, G1Y, G1Z
    real(kind=kreal)   :: G2X, G2Y, G2Z
    real(kind=kreal)   :: G3X, G3Y, G3Z
    real(kind=kreal)   :: XSUM
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !
    ! SELECTION OF LOAD TYPE
    !
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if( LTYPE.EQ.1 ) then
        NOD(1) = 1
        NOD(2) = 2
        NOD(3) = 3
        NOD(4) = 4
      else if( LTYPE.EQ.2 ) then
        NOD(1) = 8
        NOD(2) = 7
        NOD(3) = 6
        NOD(4) = 5
      else if( LTYPE.EQ.3 ) then
        NOD(1) = 5
        NOD(2) = 6
        NOD(3) = 2
        NOD(4) = 1
      else if( LTYPE.EQ.4 ) then
        NOD(1) = 6
        NOD(2) = 7
        NOD(3) = 3
        NOD(4) = 2
      else if( LTYPE.EQ.5 ) then
        NOD(1) = 7
        NOD(2) = 8
        NOD(3) = 4
        NOD(4) = 3
      else if( LTYPE.EQ.6 ) then
        NOD(1) = 8
        NOD(2) = 5
        NOD(3) = 1
        NOD(4) = 4
      endif
    endif
    ! CLEAR VECT
    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    !** SURFACE LOAD
    !
    if( ISUF.EQ.1 ) then
      ! INTEGRATION OVER SURFACE
      do IG2=1,2
        SI=XG(IG2)
        do IG1=1,2
          RI=XG(IG1)
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
          do I=1,4
            G1X=G1X+HR(I)*XX(NOD(I))
            G1Y=G1Y+HR(I)*YY(NOD(I))
            G1Z=G1Z+HR(I)*ZZ(NOD(I))
            G2X=G2X+HS(I)*XX(NOD(I))
            G2Y=G2Y+HS(I)*YY(NOD(I))
            G2Z=G2Z+HS(I)*ZZ(NOD(I))
          enddo
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          do I=1,4
            VECT(NOD(I))=VECT(NOD(I))-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)
          enddo
        enddo
      enddo
    endif
    !
    !** VOLUME LOAD
    !
    if( IVOL.EQ.1 ) then
      do I=1,NN
        PL(I)=0.0
      enddo
      ! LOOP FOR INTEGRATION POINTS
      do  LX=1,2
        RI=XG(LX)
        do LY=1,2
          SI=XG(LY)
          do LZ=1,2
            TI=XG(LZ)
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
            do I=1,NN
              XJ11=XJ11+HR(I)*XX(I)
              XJ21=XJ21+HS(I)*XX(I)
              XJ31=XJ31+HT(I)*XX(I)
              XJ12=XJ12+HR(I)*YY(I)
              XJ22=XJ22+HS(I)*YY(I)
              XJ32=XJ32+HT(I)*YY(I)
              XJ13=XJ13+HR(I)*ZZ(I)
              XJ23=XJ23+HS(I)*ZZ(I)
              XJ33=XJ33+HT(I)*ZZ(I)
            enddo
            !DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            WG=WGT(LX)*WGT(LY)*WGT(LZ)*DET
            do I=1,NN
              VECT(I)=VECT(I)-H(I)*WG*val
            enddo
          enddo
        enddo
      enddo
    endif
    !*
    return

  end subroutine heat_DFLUX_361
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_362(NN,XX,YY,ZZ,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET 362 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :FLUX IN NORMAL-DIRECTION FOR FACE-1
    !   S2   LTYPE=2  :FLUX IN NORMAL-DIRECTION FOR FACE-2
    !   S3   LTYPE=3  :FLUX IN NORMAL-DIRECTION FOR FACE-3
    !   S4   LTYPE=4  :FLUX IN NORMAL-DIRECTION FOR FACE-4
    !   S5   LTYPE=5  :FLUX IN NORMAL-DIRECTION FOR FACE-5
    !   S6   LTYPE=6  :FLUX IN NORMAL-DIRECTION FOR FACE-6
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), val, VECT(NN)
    ! LOCAL VARIABLES
    real(kind=kreal)   :: H(20), HR(20), HS(20), HT(20)
    real(kind=kreal)   :: XG(3), WGT(3)
    real(kind=kreal)   :: RI, SI, TI, RP, SP, TP, RM, SM, TM
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: DET, WG
    integer(kind=kint) :: IVOL, ISUF
    integer(kind=kint) :: NOD(8)
    integer(kind=kint) :: IG1, IG2, LX, LY, LZ, I
    real(kind=kreal)   :: VX, VY, VZ
    real(kind=kreal)   :: G1X, G1Y, G1Z
    real(kind=kreal)   :: G2X, G2Y, G2Z
    real(kind=kreal)   :: G3X, G3Y, G3Z
    real(kind=kreal)   :: XSUM
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !
    ! SELECTION OF LOAD TYPE
    !
    IVOL=0
    ISUF=0
    if( LTYPE.EQ.0 ) then
      IVOL=1
    else
      ISUF=1
      if( LTYPE.EQ.1 ) then
        NOD(1) = 1
        NOD(2) = 2
        NOD(3) = 3
        NOD(4) = 4
        NOD(5) = 9
        NOD(6) = 10
        NOD(7) = 11
        NOD(8) = 12
      else if( LTYPE.EQ.2 ) then
        NOD(1) = 8
        NOD(2) = 7
        NOD(3) = 6
        NOD(4) = 5
        NOD(5) = 15
        NOD(6) = 14
        NOD(7) = 13
        NOD(8) = 16
      else if( LTYPE.EQ.3 ) then
        NOD(1) = 5
        NOD(2) = 6
        NOD(3) = 2
        NOD(4) = 1
        NOD(5) = 13
        NOD(6) = 18
        NOD(7) = 9
        NOD(8) = 17
      else if( LTYPE.EQ.4 ) then
        NOD(1) = 6
        NOD(2) = 7
        NOD(3) = 3
        NOD(4) = 2
        NOD(5) = 14
        NOD(6) = 19
        NOD(7) = 10
        NOD(8) = 18
      else if( LTYPE.EQ.5 ) then
        NOD(1) = 7
        NOD(2) = 8
        NOD(3) = 4
        NOD(4) = 3
        NOD(5) = 15
        NOD(6) = 20
        NOD(7) = 11
        NOD(8) = 19
      else if( LTYPE.EQ.6 ) then
        NOD(1) = 8
        NOD(2) = 5
        NOD(3) = 1
        NOD(4) = 4
        NOD(5) = 16
        NOD(6) = 17
        NOD(7) = 12
        NOD(8) = 20
      endif
    endif
    ! CLEAR VECT
    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    !** SURFACE LOAD
    !
    if( ISUF.EQ.1 ) then
      ! INTEGRATION OVER SURFACE
      do IG2=1,3
        SI=XG(IG2)
        do IG1=1,3
          RI=XG(IG1)
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
          HR(6)= .5*(1.0-SI*SI)
          HR(7)=-RI*(1.0+SI)
          HR(8)=-.5*(1.0-SI*SI)
          HS(1)=-.25*RM*(-1.0-RI-SI)-0.25*RM*SM
          HS(2)=-.25*RP*(-1.0+RI-SI)-0.25*RP*SM
          HS(3)= .25*RP*(-1.0+RI+SI)+0.25*RP*SP
          HS(4)= .25*RM*(-1.0-RI+SI)+0.25*RM*SP
          HS(5)=-.5*(1.0-RI*RI)
          HS(6)=-SI*(1.0+RI)
          HS(7)= .5*(1.0-RI*RI)
          HS(8)=-SI*(1.0-RI)
          G1X=0.0
          G1Y=0.0
          G1Z=0.0
          G2X=0.0
          G2Y=0.0
          G2Z=0.0
          do I=1,8
            G1X=G1X+HR(I)*XX(NOD(I))
            G1Y=G1Y+HR(I)*YY(NOD(I))
            G1Z=G1Z+HR(I)*ZZ(NOD(I))
            G2X=G2X+HS(I)*XX(NOD(I))
            G2Y=G2Y+HS(I)*YY(NOD(I))
            G2Z=G2Z+HS(I)*ZZ(NOD(I))
          enddo
          G3X=G1Y*G2Z-G1Z*G2Y
          G3Y=G1Z*G2X-G1X*G2Z
          G3Z=G1X*G2Y-G1Y*G2X
          XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
          do I=1,8
            VECT(NOD(I))=VECT(NOD(I))-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)
          enddo
        enddo
      enddo
    endif
    !
    !** VOLUME LOAD
    !
    if( IVOL.EQ.1 ) then
      ! LOOP FOR INTEGRATION POINTS
      do  LX=1,3
        RI=XG(LX)
        do LY=1,3
          SI=XG(LY)
          do LZ=1,3
            TI=XG(LZ)
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
            do I=1,NN
              XJ11=XJ11+HR(I)*XX(I)
              XJ21=XJ21+HS(I)*XX(I)
              XJ31=XJ31+HT(I)*XX(I)
              XJ12=XJ12+HR(I)*YY(I)
              XJ22=XJ22+HS(I)*YY(I)
              XJ32=XJ32+HT(I)*YY(I)
              XJ13=XJ13+HR(I)*ZZ(I)
              XJ23=XJ23+HS(I)*ZZ(I)
              XJ33=XJ33+HT(I)*ZZ(I)
            enddo
            !DETERMINANT OF JACOBIAN
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            WG=WGT(LX)*WGT(LY)*WGT(LZ)*DET
            do I=1,NN
              VECT(I)=VECT(I)-H(I)*WG*val
            enddo
          enddo
        enddo
      enddo
    endif
    !*
    return

  end subroutine heat_DFLUX_362
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_731(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET S3 DFLUX
    !**
    !   BF   LTYPE=0  :BODY FLUX
    !   S1   LTYPE=1  :SURFACE FLUX
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I
    real(kind=kreal)   :: V1X, V1Y, V1Z, V2X, V2Y, V2Z, V3X, V3Y, V3Z, AA
    !
    if( LTYPE.EQ.0 ) then
      THICK = THICK
    elseif( LTYPE.EQ.1 ) then
      THICK = 1.0d0
    else
      THICK = 0.0d0
    endif

    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    V1X=XX(2)-XX(1)
    V1Y=YY(2)-YY(1)
    V1Z=ZZ(2)-ZZ(1)
    V2X=XX(3)-XX(1)
    V2Y=YY(3)-YY(1)
    V2Z=ZZ(3)-ZZ(1)
    V3X= V1Y*V2Z-V1Z*V2Y
    V3Y= V1Z*V2X-V1X*V2Z
    V3Z= V1X*V2Y-V1Y*V2X
    AA=0.5*dsqrt( V3X*v3X + V3Y*V3Y + V3Z*V3Z )
    VECT(1)= -val*AA*THICK/3.0
    VECT(2)= -val*AA*THICK/3.0
    VECT(3)= -val*AA*THICK/3.0
    !
    return

  end subroutine heat_DFLUX_731
  !----------------------------------------------------------------------*
  subroutine heat_DFLUX_741(NN,XX,YY,ZZ,THICK,LTYPE,val,VECT)
    !----------------------------------------------------------------------*
    !**
    !**  SET S4 DFLUX
    !**
    !   BF   LTYPE=0  : BODY FLUX
    !   S1   LTYPE=1  : SURFACE FLUX
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), THICK, val, VECT(NN)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IG1, IG2
    real(kind=kreal)   :: RI, SI, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: H(4), HR(4), HS(4), HT(4), PL(4)
    real(kind=kreal)   :: XG(2), WGT(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG  / -0.5773502691896,0.5773502691896 /
    data WGT /  1.0, 1.0 /
    !
    ! SELECTION OF LOAD TYPE
    !
    if    ( LTYPE.EQ.0 ) then
      THICK = 1.0d0 * THICK
    elseif( LTYPE.EQ.1 ) then
      THICK = 1.0d0
    else
      THICK = 0.0d0
    endif
    !
    do I=1,NN
      VECT(I)=0.0
    enddo
    !
    do IG2=1,2
      SI=XG(IG2)
      do IG1=1,2
        RI=XG(IG1)
        !
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
        !
        G1X=0.0
        G1Y=0.0
        G1Z=0.0
        G2X=0.0
        G2Y=0.0
        G2Z=0.0
        do I=1,NN
          G1X=G1X+HR(I)*XX(I)
          G1Y=G1Y+HR(I)*YY(I)
          G1Z=G1Z+HR(I)*ZZ(I)
          G2X=G2X+HS(I)*XX(I)
          G2Y=G2Y+HS(I)*YY(I)
          G2Z=G2Z+HS(I)*ZZ(I)
        enddo
        !
        G3X=G1Y*G2Z-G1Z*G2Y
        G3Y=G1Z*G2X-G1X*G2Z
        G3Z=G1X*G2Y-G1Y*G2X
        XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
        do I=1,NN
          VECT(I)=VECT(I)-XSUM*val*WGT(IG1)*WGT(IG2)*H(I)*THICK
        enddo
        !
      enddo
    enddo
    return

  end subroutine heat_DFLUX_741
end module m_heat_LIB_DFLUX
