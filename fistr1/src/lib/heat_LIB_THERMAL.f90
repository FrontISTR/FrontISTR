!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines for calculating heat conductive
!! matrix for various elements

module m_heat_lib_thermal
contains

  !> GET CONDUCTIVITY
  subroutine heat_GET_CONDUCTIVITY ( Tpoi,imat,COND,ntab,temp,funcA,funcB )
    !C***
    !C*** GET CONDUCTIVITY
    !C***
    use hecmw

    implicit real(kind=kreal) (A-H,O-Z)
    dimension COND(3),temp(ntab),funcA(ntab+1),funcB(ntab+1)

    COND(1)= funcB(1)
    COND(2)= COND(1)
    COND(3)= COND(1)

    if( ntab .GT. 1 ) then
      if( Tpoi.LE.temp(1) ) then
        itab= 1
      else
        itab= ntab + 1
        do ikk= 1, ntab - 1
          if( Tpoi.GT.temp(ikk).AND.Tpoi.LE.temp(ikk+1) ) then
            itab= ikk + 1
            exit
          endif
        enddo
      endif
      COND(1)= funcA(itab)*Tpoi+ funcB(itab)
      COND(2)= COND(1)
      COND(3)= COND(1)
    endif
    return
  end subroutine heat_GET_CONDUCTIVITY

  !C************************************************************************
  !C*  THERMAL_111 ( NN,XX,YY,ZZ,TT,IMAT,ASECT,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_231 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_241 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_341 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_351 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_361 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )
  !c*  THERMAL_531 ( NN,XXX,YYY,ZZZ,TEMP,TZERO,CROSS,HH,RR1,RR2,SS )
  !C*  THERMAL_541 ( NN,XXX,YYY,ZZZ,TEMP,TZERO,CROSS,HH,RR1,RR2,SS )
  !C*  THERMAL_731 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_741 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_232 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_242 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_342 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_352 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )
  !C*  THERMAL_362 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

  !> CALCULATION 1D 2 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_111 (NN,XX,YY,ZZ,TT,IMAT,ASECT,SS &
      ,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 1D 2 NODE CONDUCTANCE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !
    DX = XX(2) - XX(1)
    DY = YY(2) - YY(1)
    DZ = ZZ(2) - ZZ(1)

    AL = dsqrt( DX*DX + DY*DY + DZ*DZ )

    !*CONDUCTIVITY AT CURRENT TEMPERATURE

    CTEMP=0.0
    do I = 1, NN
      CTEMP = CTEMP + TT(I) * 0.5d0
    enddo

    call heat_GET_CONDUCTIVITY( CTEMP, IMAT, CC, ntab, temp, funcA, funcB )
    !
    SS(1) =  CC(1) * ASECT * AL
    SS(2) = -SS(1)
    SS(3) = -SS(1)
    SS(4) =  SS(1)
    !
    return
  end subroutine heat_THERMAL_111

  !> CALCULATION 2D 3 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_231 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !C*
    !C*CALCULATION 2D 3 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(2), WGT(2), H(3), HL1(3), HL2(3), HL3(3)
    dimension BX(3), BY(3), BZ(3), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    data WGT/1.0D0,1.0D0/
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
    do L2=1,2
      XL2=XG(L2)
      X2 =(XL2+1.0)*0.5
      do L1=1,2
        XL1=XG(L1)
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
        do I=1,NN
          XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
          XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
          XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
          XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
        enddo
        DET=XJ11*XJ22-XJ21*XJ12
        DUM=1.0/DET
        XJI11= XJ22*DUM
        XJI12=-XJ12*DUM
        XJI21=-XJ21*DUM
        XJI22= XJ11*DUM
        do J=1, NN
          BX(J)=XJI11*(HL1(J)-HL3(J))+XJI12*(HL2(J)-HL3(J))
          BY(J)=XJI21*(HL1(J)-HL3(J))+XJI22*(HL2(J)-HL3(J))
        enddo
        !C*CONDUCTIVITY AT CURRENT TEMPERATURE
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
        !C*WEIGT VALUE AT GAUSSIAN POINT
        WGX=CC(1)*WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        WGY=CC(2)*WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        !C
        IJ = 0
        do I = 1, NN
          do J = 1, NN
            IJ = IJ + 1
            SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX+BY(I)*BY(J)*WGY
          enddo
        enddo
        !C
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_231

  !> CALCULATION 2D 4 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_241 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !C*
    !C*CALCULATION 2D 4 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(2), WGT(2), H(4), HR(4), HS(4), HT(4)
    dimension BX(4), BY(4), BZ(4), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    data XG/-0.5773502691896D0, 0.5773502691896D0/
    data WGT/1.0D0,1.0D0/
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
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
        DUM=1.0/DET
        XJI11= XJ22*DUM
        XJI12=-XJ12*DUM
        XJI21=-XJ21*DUM
        XJI22= XJ11*DUM
        do J=1, NN
          BX(J)=XJI11*HR(J)+XJI12*HS(J)
          BY(J)=XJI21*HR(J)+XJI22*HS(J)
        enddo
        !*CONDUCTIVITY AT CURRENT TEMPERATURE
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
        !*WEIGT VALUE AT GAUSSIAN POINT
        WGX=CC(1)*WGT(LX)*WGT(LY)*DET*THICK
        WGY=CC(2)*WGT(LX)*WGT(LY)*DET*THICK
        !C
        IJ = 0
        do I = 1, NN
          do J = 1, NN
            IJ = IJ + 1
            SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX+BY(I)*BY(J)*WGY
          enddo
        enddo
        !C
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_241

  !> CALCULATION 3D 4 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_341 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 3D 4 NODE CONDUCTANCE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(2), WGT(2), H(4), HR(4), HS(4), HT(4)
    dimension BX(4), BY(4), BZ(4), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    !*      DATA XG/-1.0,1.0/
    data WGT/1.0D0,1.0D0/
    !
    !* LOOP FOR INTEGRATION POINTS
    do L3=1,2
      XL3=XG(L3)
      X3 =(XL3+1.0)*0.5
      do L2=1,2
        XL2=XG(L2)
        X2 =(1.0-X3)*(XL2+1.0)*0.5
        do L1=1,2
          XL1=XG(L1)
          X1=(1.0-X2-X3)*(XL1+1.0)*0.5
          !INTERPOLATION FUNCTION
          H(1)=X1
          H(2)=X2
          H(3)=X3
          H(4)=1.0-X1-X2-X3
          !DERIVATIVE OF INTERPOLATION FUNCTION
          !  FOR L1-COORDINATE
          HR(1)=1.0
          HR(2)=0.0
          HR(3)=0.0
          HR(4)=-1.0
          !  FOR L2-COORDINATE
          HS(1)=0.0
          HS(2)=1.0
          HS(3)=0.0
          HS(4)=-1.0
          !  FOR ZETA-COORDINATE
          HT(1)=0.0
          HT(2)=0.0
          HT(3)=1.0
          HT(4)=-1.0

          !JACOBI MATRIX
          XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
          XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
          XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4)
          XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
          XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
          XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4)
          XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
          XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)
          XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4)
          !
          !DETERMINANT OF JACOBIAN
          !
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32

          !CONDUCTIVITY AT CURRENT TEMPERATURE

          CTEMP=0.0
          do I=1,4
            CTEMP=CTEMP+H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)

          ! WEIGT VALUE AT GAUSSIAN POINT

          DET = -DET
          WGX=CC(1)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          WGY=CC(2)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          WGZ=CC(3)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125

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
          do J=1,4
            BX(J)=XJI11*HR(J)+XJI12*HS(J)+XJI13*HT(J)
            BY(J)=XJI21*HR(J)+XJI22*HS(J)+XJI23*HT(J)
            BZ(J)=XJI31*HR(J)+XJI32*HS(J)+XJI33*HT(J)
          enddo
          !
          SS( 1)=SS( 1)+BX(1)*BX(1)*WGX &
            +BY(1)*BY(1)*WGY &
            +BZ(1)*BZ(1)*WGZ
          SS( 2)=SS( 2)+BX(1)*BX(2)*WGX &
            +BY(1)*BY(2)*WGY &
            +BZ(1)*BZ(2)*WGZ
          SS( 6)=SS( 6)+BX(2)*BX(2)*WGX &
            +BY(2)*BY(2)*WGY &
            +BZ(2)*BZ(2)*WGZ
          SS( 3)=SS( 3)+BX(1)*BX(3)*WGX &
            +BY(1)*BY(3)*WGY &
            +BZ(1)*BZ(3)*WGZ
          SS( 7)=SS( 7)+BX(2)*BX(3)*WGX &
            +BY(2)*BY(3)*WGY &
            +BZ(2)*BZ(3)*WGZ
          SS(11)=SS(11)+BX(3)*BX(3)*WGX &
            +BY(3)*BY(3)*WGY &
            +BZ(3)*BZ(3)*WGZ
          SS( 4)=SS( 4)+BX(1)*BX(4)*WGX &
            +BY(1)*BY(4)*WGY &
            +BZ(1)*BZ(4)*WGZ
          SS( 8)=SS( 8)+BX(2)*BX(4)*WGX &
            +BY(2)*BY(4)*WGY &
            +BZ(2)*BZ(4)*WGZ
          SS(12)=SS(12)+BX(3)*BX(4)*WGX &
            +BY(3)*BY(4)*WGY &
            +BZ(3)*BZ(4)*WGZ
          SS(16)=SS(16)+BX(4)*BX(4)*WGX &
            +BY(4)*BY(4)*WGY &
            +BZ(4)*BZ(4)*WGZ
          !
        enddo
      enddo
    enddo
    !
    SS( 5) = SS( 2)
    SS( 9) = SS( 3)
    SS(10) = SS( 7)
    SS(13) = SS( 4)
    SS(14) = SS( 8)
    SS(15) = SS(12)
    !
    return
  end subroutine heat_THERMAL_341

  !> CALCULATION 3D 6 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_351 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 3D 6 NODE CONDUCTANCE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(2), WGT(2), XG1(3), XG2(3), WGT1(3) &
      , H(6), HR(6), HS(6), HT(6) &
      , BX(6), BY(6), BZ(6), CC(3)

    dimension temp(*),funcA(*),funcB(*)

    data XG   /-0.5773502691896D0,0.5773502691896D0/
    data WGT  /1.0D0,1.0D0/
    data XG1  /0.6666666667D0,0.16666666667D0,0.16666666667D0/
    data XG2  /0.1666666667D0,0.66666666667D0,0.16666666667D0/
    data WGT1 /0.1666666667D0,0.16666666667D0,0.16666666667D0/


    !*INTEGRATION FOR Z-COORDINATE
    do LZ=1,2
      ZI = XG(LZ)

      !*INTEGRAION FOR L1,L2-COORDINATE
      do L12=1,3
        X1 = XG1(L12)
        X2 = XG2(L12)

        !*INTERPOLATION FUNCTION
        H(1) = X1*(1.0-ZI)*0.5
        H(2) = X2*(1.0-ZI)*0.5
        H(3) = (1.0-X1-X2)*(1.0-ZI)*0.5
        H(4) = X1*(1.0+ZI)*0.5
        H(5) = X2*(1.0+ZI)*0.5
        H(6) = (1.0-X1-X2)*(1.0+ZI)*0.5

        !*DERIVATIVE OF INTERPOLATION FUNCTION
        !*FOR L1-COORDINATE
        HR(1) = (1.0-ZI)*0.5
        HR(2) = 0.0
        HR(3) =-(1.0-ZI)*0.5
        HR(4) = (1.0+ZI)*0.5
        HR(5) = 0.0
        HR(6) =-(1.0+ZI)*0.5
        !*FOR L2-COORDINATE
        HS(1) =  0.0
        HS(2) =  (1.0-ZI)*0.5
        HS(3) = -(1.0-ZI)*0.5
        HS(4) = 0.0
        HS(5) =  (1.0+ZI)*0.5
        HS(6) = -(1.0+ZI)*0.5
        !*FOR ZETA-COORDINATE
        HT(1) = -0.5*X1
        HT(2) = -0.5*X2
        HT(3) = -0.5*(1.0-X1-X2)
        HT(4) =  0.5*X1
        HT(5) =  0.5*X2
        HT(6) =  0.5*(1.0-X1-X2)

        !*JACOBI MATRIX
        XJ11 = HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4) &
          +HR(5)*XX(5)+HR(6)*XX(6)
        XJ21 = HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4) &
          +HS(5)*XX(5)+HS(6)*XX(6)
        XJ31 = HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4) &
          +HT(5)*XX(5)+HT(6)*XX(6)

        XJ12 = HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4) &
          +HR(5)*YY(5)+HR(6)*YY(6)
        XJ22 = HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4) &
          +HS(5)*YY(5)+HS(6)*YY(6)
        XJ32 = HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4) &
          +HT(5)*YY(5)+HT(6)*YY(6)

        XJ13 = HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4) &
          +HR(5)*ZZ(5)+HR(6)*ZZ(6)
        XJ23 = HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4) &
          +HS(5)*ZZ(5)+HS(6)*ZZ(6)
        XJ33 = HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4) &
          +HT(5)*ZZ(5)+HT(6)*ZZ(6)

        !*DETERMINANT OF JACOBIAN
        DET = XJ11*XJ22*XJ33 &
          +XJ12*XJ23*XJ31 &
          +XJ13*XJ21*XJ32 &
          -XJ13*XJ22*XJ31 &
          -XJ12*XJ21*XJ33 &
          -XJ11*XJ23*XJ32

        !*CONDUCTIVITY AT CURRENT TEMPERATURE
        CTEMP=0.0
        do I=1,6
          CTEMP = CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CONDUCTIVITY (CTEMP,IMAT,CC,ntab,temp,funcA,funcB)

        !* WEIGT VALUE AT GAUSSIAN POINT
        WGX = CC(1)*WGT1(L12)*WGT(LZ)*DET
        WGY = CC(2)*WGT1(L12)*WGT(LZ)*DET
        WGZ = CC(3)*WGT1(L12)*WGT(LZ)*DET
        !
        !* INVERSION OF JACOBIAN
        DUM   = 1.0/DET
        XJI11 = DUM*( XJ22*XJ33-XJ23*XJ32)
        XJI21 = DUM*(-XJ21*XJ33+XJ23*XJ31)
        XJI31 = DUM*( XJ21*XJ32-XJ22*XJ31)
        XJI12 = DUM*(-XJ12*XJ33+XJ13*XJ32)
        XJI22 = DUM*( XJ11*XJ33-XJ13*XJ31)
        XJI32 = DUM*(-XJ11*XJ32+XJ12*XJ31)
        XJI13 = DUM*( XJ12*XJ23-XJ13*XJ22)
        XJI23 = DUM*(-XJ11*XJ23+XJ13*XJ21)
        XJI33 = DUM*( XJ11*XJ22-XJ12*XJ21)
        !
        do J=1,6
          BX(J) = XJI11*HR(J)+XJI12*HS(J)+XJI13*HT(J)
          BY(J) = XJI21*HR(J)+XJI22*HS(J)+XJI23*HT(J)
          BZ(J) = XJI31*HR(J)+XJI32*HS(J)+XJI33*HT(J)
        enddo
        !
        SS( 1) = SS( 1) + BX(1)*BX(1)*WGX &
          + BY(1)*BY(1)*WGY &
          + BZ(1)*BZ(1)*WGZ
        SS( 2) = SS( 2) + BX(1)*BX(2)*WGX &
          + BY(1)*BY(2)*WGY &
          + BZ(1)*BZ(2)*WGZ
        SS( 3) = SS( 3) + BX(1)*BX(3)*WGX &
          + BY(1)*BY(3)*WGY &
          + BZ(1)*BZ(3)*WGZ
        SS( 4) = SS( 4) + BX(1)*BX(4)*WGX &
          + BY(1)*BY(4)*WGY &
          + BZ(1)*BZ(4)*WGZ
        SS( 5) = SS( 5) + BX(1)*BX(5)*WGX &
          + BY(1)*BY(5)*WGY &
          + BZ(1)*BZ(5)*WGZ
        SS( 6) = SS( 6) + BX(1)*BX(6)*WGX &
          + BY(1)*BY(6)*WGY &
          + BZ(1)*BZ(6)*WGZ
        SS( 8) = SS( 8) + BX(2)*BX(2)*WGX &
          + BY(2)*BY(2)*WGY &
          + BZ(2)*BZ(2)*WGZ
        SS( 9) = SS( 9) + BX(2)*BX(3)*WGX &
          + BY(2)*BY(3)*WGY &
          + BZ(2)*BZ(3)*WGZ
        SS(10) = SS(10) + BX(2)*BX(4)*WGX &
          + BY(2)*BY(4)*WGY &
          + BZ(2)*BZ(4)*WGZ
        SS(11) = SS(11) + BX(2)*BX(5)*WGX &
          + BY(2)*BY(5)*WGY &
          + BZ(2)*BZ(5)*WGZ
        SS(12) = SS(12) + BX(2)*BX(6)*WGX &
          + BY(2)*BY(6)*WGY &
          + BZ(2)*BZ(6)*WGZ
        SS(15) = SS(15) + BX(3)*BX(3)*WGX &
          + BY(3)*BY(3)*WGY &
          + BZ(3)*BZ(3)*WGZ
        SS(16) = SS(16) + BX(3)*BX(4)*WGX &
          + BY(3)*BY(4)*WGY &
          + BZ(3)*BZ(4)*WGZ
        SS(17) = SS(17) + BX(3)*BX(5)*WGX &
          + BY(3)*BY(5)*WGY &
          + BZ(3)*BZ(5)*WGZ
        SS(18) = SS(18) + BX(3)*BX(6)*WGX &
          + BY(3)*BY(6)*WGY &
          + BZ(3)*BZ(6)*WGZ
        SS(22) = SS(22) + BX(4)*BX(4)*WGX &
          + BY(4)*BY(4)*WGY &
          + BZ(4)*BZ(4)*WGZ
        SS(23) = SS(23) + BX(4)*BX(5)*WGX &
          + BY(4)*BY(5)*WGY &
          + BZ(4)*BZ(5)*WGZ
        SS(24) = SS(24) + BX(4)*BX(6)*WGX &
          + BY(4)*BY(6)*WGY &
          + BZ(4)*BZ(6)*WGZ
        SS(29) = SS(29) + BX(5)*BX(5)*WGX &
          + BY(5)*BY(5)*WGY &
          + BZ(5)*BZ(5)*WGZ
        SS(30) = SS(30) + BX(5)*BX(6)*WGX &
          + BY(5)*BY(6)*WGY &
          + BZ(5)*BZ(6)*WGZ
        SS(36) = SS(36) + BX(6)*BX(6)*WGX &
          + BY(6)*BY(6)*WGY &
          + BZ(6)*BZ(6)*WGZ

      enddo
    enddo

    SS( 7) = SS( 2)
    SS(13) = SS( 3)
    SS(14) = SS( 9)
    SS(19) = SS( 4)
    SS(20) = SS(10)
    SS(21) = SS(16)
    SS(25) = SS( 5)
    SS(26) = SS(11)
    SS(27) = SS(17)
    SS(28) = SS(23)
    SS(31) = SS( 6)
    SS(32) = SS(12)
    SS(33) = SS(18)
    SS(34) = SS(24)
    SS(35) = SS(30)

    return
  end subroutine heat_THERMAL_351

  !> CALCULATION 3D 8 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_361 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 3D 8 NODE CONDUCTANCE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(2), WGT(2), H(8), HR(8), HS(8), HT(8)
    dimension BX(8), BY(8), BZ(8), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !*      DATA XG/-1.0D0,1.0D0/
    data WGT/1.0D0,1.0D0/
    data XG/-0.5773502691896D0, 0.5773502691896D0/
    !
    do LX=1,2
      RI=XG(LX)
      do LY=1,2
        SI=XG(LY)
        do LZ=1,2
          TI=XG(LZ)
          !
          RP=1.0+RI
          SP=1.0+SI
          TP=1.0+TI
          RM=1.0-RI
          SM=1.0-SI
          TM=1.0-TI
          !
          !*INTERPOLATION FUNCTION
          H(1)=0.125*RM*SM*TM
          H(2)=0.125*RP*SM*TM
          H(3)=0.125*RP*SP*TM
          H(4)=0.125*RM*SP*TM
          H(5)=0.125*RM*SM*TP
          H(6)=0.125*RP*SM*TP
          H(7)=0.125*RP*SP*TP
          H(8)=0.125*RM*SP*TP
          !
          !*DERIVATIVE OF INTERPOLATION FUNCTION
          !*  FOR R-COORDINATE
          HR(1)=-.125*SM*TM
          HR(2)= .125*SM*TM
          HR(3)= .125*SP*TM
          HR(4)=-.125*SP*TM
          HR(5)=-.125*SM*TP
          HR(6)= .125*SM*TP
          HR(7)= .125*SP*TP
          HR(8)=-.125*SP*TP
          !*  FOR S-COORDINATE
          HS(1)=-.125*RM*TM
          HS(2)=-.125*RP*TM
          HS(3)= .125*RP*TM
          HS(4)= .125*RM*TM
          HS(5)=-.125*RM*TP
          HS(6)=-.125*RP*TP
          HS(7)= .125*RP*TP
          HS(8)= .125*RM*TP
          !*  FOR T-COORDINATE
          HT(1)=-.125*RM*SM
          HT(2)=-.125*RP*SM
          HT(3)=-.125*RP*SP
          HT(4)=-.125*RM*SP
          HT(5)= .125*RM*SM
          HT(6)= .125*RP*SM
          HT(7)= .125*RP*SP
          HT(8)= .125*RM*SP
          !
          !*JACOBI MATRIX
          XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4) &
            +HR(5)*XX(5)+HR(6)*XX(6)+HR(7)*XX(7)+HR(8)*XX(8)
          XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4) &
            +HS(5)*XX(5)+HS(6)*XX(6)+HS(7)*XX(7)+HS(8)*XX(8)
          XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4) &
            +HT(5)*XX(5)+HT(6)*XX(6)+HT(7)*XX(7)+HT(8)*XX(8)
          !
          XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4) &
            +HR(5)*YY(5)+HR(6)*YY(6)+HR(7)*YY(7)+HR(8)*YY(8)
          XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4) &
            +HS(5)*YY(5)+HS(6)*YY(6)+HS(7)*YY(7)+HS(8)*YY(8)
          XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4) &
            +HT(5)*YY(5)+HT(6)*YY(6)+HT(7)*YY(7)+HT(8)*YY(8)
          !
          XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4) &
            +HR(5)*ZZ(5)+HR(6)*ZZ(6)+HR(7)*ZZ(7)+HR(8)*ZZ(8)
          XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4) &
            +HS(5)*ZZ(5)+HS(6)*ZZ(6)+HS(7)*ZZ(7)+HS(8)*ZZ(8)
          XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4) &
            +HT(5)*ZZ(5)+HT(6)*ZZ(6)+HT(7)*ZZ(7)+HT(8)*ZZ(8)
          !
          !*DETERMINANT OF JACOBIAN
          !
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !
          !* INVERSION OF JACOBIAN
          !
          DET = -DET
          DUM=1.0/DET
          !
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
          do J=1, NN
            BX(J)=XJI11*HR(J)+XJI12*HS(J)+XJI13*HT(J)
            BY(J)=XJI21*HR(J)+XJI22*HS(J)+XJI23*HT(J)
            BZ(J)=XJI31*HR(J)+XJI32*HS(J)+XJI33*HT(J)
          enddo
          !
          !*CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
          !
          !* WEIGT VALUE AT GAUSSIAN POINT
          WGX=-CC(1)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          WGY=-CC(2)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          WGZ=-CC(3)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          !
          SS( 1)=SS( 1)+BX(1)*BX(1)*WGX &
            +BY(1)*BY(1)*WGY &
            +BZ(1)*BZ(1)*WGZ
          SS( 2)=SS( 2)+BX(1)*BX(2)*WGX &
            +BY(1)*BY(2)*WGY &
            +BZ(1)*BZ(2)*WGZ
          SS( 3)=SS( 3)+BX(1)*BX(3)*WGX &
            +BY(1)*BY(3)*WGY &
            +BZ(1)*BZ(3)*WGZ
          SS( 4)=SS( 4)+BX(1)*BX(4)*WGX &
            +BY(1)*BY(4)*WGY &
            +BZ(1)*BZ(4)*WGZ
          SS( 5)=SS( 5)+BX(1)*BX(5)*WGX &
            +BY(1)*BY(5)*WGY &
            +BZ(1)*BZ(5)*WGZ
          SS( 6)=SS( 6)+BX(1)*BX(6)*WGX &
            +BY(1)*BY(6)*WGY &
            +BZ(1)*BZ(6)*WGZ
          SS( 7)=SS( 7)+BX(1)*BX(7)*WGX &
            +BY(1)*BY(7)*WGY &
            +BZ(1)*BZ(7)*WGZ
          SS( 8)=SS( 8)+BX(1)*BX(8)*WGX &
            +BY(1)*BY(8)*WGY &
            +BZ(1)*BZ(8)*WGZ
          SS(10)=SS(10)+BX(2)*BX(2)*WGX &
            +BY(2)*BY(2)*WGY &
            +BZ(2)*BZ(2)*WGZ
          SS(11)=SS(11)+BX(2)*BX(3)*WGX &
            +BY(2)*BY(3)*WGY &
            +BZ(2)*BZ(3)*WGZ
          SS(12)=SS(12)+BX(2)*BX(4)*WGX &
            +BY(2)*BY(4)*WGY &
            +BZ(2)*BZ(4)*WGZ
          SS(13)=SS(13)+BX(2)*BX(5)*WGX &
            +BY(2)*BY(5)*WGY &
            +BZ(2)*BZ(5)*WGZ
          SS(14)=SS(14)+BX(2)*BX(6)*WGX &
            +BY(2)*BY(6)*WGY &
            +BZ(2)*BZ(6)*WGZ
          SS(15)=SS(15)+BX(2)*BX(7)*WGX &
            +BY(2)*BY(7)*WGY &
            +BZ(2)*BZ(7)*WGZ
          SS(16)=SS(16)+BX(2)*BX(8)*WGX &
            +BY(2)*BY(8)*WGY &
            +BZ(2)*BZ(8)*WGZ
          SS(19)=SS(19)+BX(3)*BX(3)*WGX &
            +BY(3)*BY(3)*WGY &
            +BZ(3)*BZ(3)*WGZ
          SS(20)=SS(20)+BX(3)*BX(4)*WGX &
            +BY(3)*BY(4)*WGY &
            +BZ(3)*BZ(4)*WGZ
          SS(21)=SS(21)+BX(3)*BX(5)*WGX &
            +BY(3)*BY(5)*WGY &
            +BZ(3)*BZ(5)*WGZ
          SS(22)=SS(22)+BX(3)*BX(6)*WGX &
            +BY(3)*BY(6)*WGY &
            +BZ(3)*BZ(6)*WGZ
          SS(23)=SS(23)+BX(3)*BX(7)*WGX &
            +BY(3)*BY(7)*WGY &
            +BZ(3)*BZ(7)*WGZ
          SS(24)=SS(24)+BX(3)*BX(8)*WGX &
            +BY(3)*BY(8)*WGY &
            +BZ(3)*BZ(8)*WGZ
          SS(28)=SS(28)+BX(4)*BX(4)*WGX &
            +BY(4)*BY(4)*WGY &
            +BZ(4)*BZ(4)*WGZ
          SS(29)=SS(29)+BX(4)*BX(5)*WGX &
            +BY(4)*BY(5)*WGY &
            +BZ(4)*BZ(5)*WGZ
          SS(30)=SS(30)+BX(4)*BX(6)*WGX &
            +BY(4)*BY(6)*WGY &
            +BZ(4)*BZ(6)*WGZ
          SS(31)=SS(31)+BX(4)*BX(7)*WGX &
            +BY(4)*BY(7)*WGY &
            +BZ(4)*BZ(7)*WGZ
          SS(32)=SS(32)+BX(4)*BX(8)*WGX &
            +BY(4)*BY(8)*WGY &
            +BZ(4)*BZ(8)*WGZ
          SS(37)=SS(37)+BX(5)*BX(5)*WGX &
            +BY(5)*BY(5)*WGY &
            +BZ(5)*BZ(5)*WGZ
          SS(38)=SS(38)+BX(5)*BX(6)*WGX &
            +BY(5)*BY(6)*WGY &
            +BZ(5)*BZ(6)*WGZ
          SS(39)=SS(39)+BX(5)*BX(7)*WGX &
            +BY(5)*BY(7)*WGY &
            +BZ(5)*BZ(7)*WGZ
          SS(40)=SS(40)+BX(5)*BX(8)*WGX &
            +BY(5)*BY(8)*WGY &
            +BZ(5)*BZ(8)*WGZ
          SS(46)=SS(46)+BX(6)*BX(6)*WGX &
            +BY(6)*BY(6)*WGY &
            +BZ(6)*BZ(6)*WGZ
          SS(47)=SS(47)+BX(6)*BX(7)*WGX &
            +BY(6)*BY(7)*WGY &
            +BZ(6)*BZ(7)*WGZ
          SS(48)=SS(48)+BX(6)*BX(8)*WGX &
            +BY(6)*BY(8)*WGY &
            +BZ(6)*BZ(8)*WGZ
          SS(55)=SS(55)+BX(7)*BX(7)*WGX &
            +BY(7)*BY(7)*WGY &
            +BZ(7)*BZ(7)*WGZ
          SS(56)=SS(56)+BX(7)*BX(8)*WGX &
            +BY(7)*BY(8)*WGY &
            +BZ(7)*BZ(8)*WGZ
          SS(64)=SS(64)+BX(8)*BX(8)*WGX &
            +BY(8)*BY(8)*WGY &
            +BZ(8)*BZ(8)*WGZ
          !
        enddo
      enddo
    enddo
    !
    SS( 9) = SS( 2)
    SS(17) = SS( 3)
    SS(18) = SS(11)
    SS(25) = SS( 4)
    SS(26) = SS(12)
    SS(27) = SS(20)
    SS(33) = SS( 5)
    SS(34) = SS(13)
    SS(35) = SS(21)
    SS(36) = SS(29)
    SS(41) = SS( 6)
    SS(42) = SS(14)
    SS(43) = SS(22)
    SS(44) = SS(30)
    SS(45) = SS(38)
    SS(49) = SS( 7)
    SS(50) = SS(15)
    SS(51) = SS(23)
    SS(52) = SS(31)
    SS(53) = SS(39)
    SS(54) = SS(47)
    SS(57) = SS( 8)
    SS(58) = SS(16)
    SS(59) = SS(24)
    SS(60) = SS(32)
    SS(61) = SS(40)
    SS(62) = SS(48)
    SS(63) = SS(56)
    !
    return
  end subroutine heat_THERMAL_361

  !> CALCULATION 3D 6 NODE INTERFACE ELEMENT
  subroutine heat_THERMAL_531 ( NN,XXX,YYY,ZZZ,TEMP,TZERO,THICK,HH,RR1,RR2,SS )
    !*
    !* CALCULATION 3D 6 NODE INTERFACE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XXX(NN), YYY(NN), ZZZ(NN), TEMP(NN), SS(NN*NN)
    dimension XX(3), YY(3), ZZ(3)

    write(*,*) 'TYPE=531 not yet available..'
    stop
  end subroutine heat_THERMAL_531

  !> CALCULATION 3D 8 NODE INTERFACE ELEMENT
  subroutine heat_THERMAL_541 ( NN,XXX,YYY,ZZZ,TEMP,TZERO,THICK,HH,RR1,RR2,SS )
    !*
    !* CALCULATION 3D 8 NODE INTERFACE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XXX(NN), YYY(NN), ZZZ(NN), TEMP(NN), SS(NN*NN)
    dimension XX(4), YY(4), ZZ(4)

    XX(1)=XXX(1)
    XX(2)=XXX(2)
    XX(3)=XXX(3)
    XX(4)=XXX(4)
    YY(1)=YYY(1)
    YY(2)=YYY(2)
    YY(3)=YYY(3)
    YY(4)=YYY(4)
    ZZ(1)=ZZZ(1)
    ZZ(2)=ZZZ(2)
    ZZ(3)=ZZZ(3)
    ZZ(4)=ZZZ(4)
    call heat_get_area ( XX,YY,ZZ,SA )

    XX(1)=XXX(5)
    XX(2)=XXX(6)
    XX(3)=XXX(7)
    XX(4)=XXX(8)
    YY(1)=YYY(5)
    YY(2)=YYY(6)
    YY(3)=YYY(7)
    YY(4)=YYY(8)
    ZZ(1)=ZZZ(5)
    ZZ(2)=ZZZ(6)
    ZZ(3)=ZZZ(7)
    ZZ(4)=ZZZ(8)
    call heat_get_area ( XX,YY,ZZ,SB )

    T1Z=TEMP(1)-TZERO
    T2Z=TEMP(2)-TZERO
    T3Z=TEMP(3)-TZERO
    T4Z=TEMP(4)-TZERO
    T5Z=TEMP(5)-TZERO
    T6Z=TEMP(6)-TZERO
    T7Z=TEMP(7)-TZERO
    T8Z=TEMP(8)-TZERO

    RRR1 = RR1**0.25
    RRR2 = RR2**0.25

    HA1= (( RRR1 * T1Z )**2 + ( RRR2 * T5Z )**2 ) &
      &    * ( RRR1 * T1Z      +   RRR2 * T5Z )    * RRR1
    HA2= (( RRR1 * T2Z )**2 + ( RRR2 * T6Z )**2 ) &
      &    * ( RRR1 * T2Z      +   RRR2 * T6Z )    * RRR1
    HA3= (( RRR1 * T3Z )**2 + ( RRR2 * T7Z )**2 ) &
      &    * ( RRR1 * T3Z      +   RRR2 * T7Z )    * RRR1
    HA4= (( RRR1 * T4Z )**2 + ( RRR2 * T8Z )**2 ) &
      &    * ( RRR1 * T4Z      +   RRR2 * T8Z )    * RRR1

    HB1= (( RRR1 * T1Z )**2 + ( RRR2 * T5Z )**2 ) &
      &    * ( RRR1 * T1Z      +   RRR2 * T5Z )    * RRR2
    HB2= (( RRR1 * T2Z )**2 + ( RRR2 * T6Z )**2 ) &
      &    * ( RRR1 * T2Z      +   RRR2 * T6Z )    * RRR2
    HB3= (( RRR1 * T3Z )**2 + ( RRR2 * T7Z )**2 ) &
      &    * ( RRR1 * T3Z      +   RRR2 * T7Z )    * RRR2
    HB4= (( RRR1 * T4Z )**2 + ( RRR2 * T8Z )**2 ) &
      &    * ( RRR1 * T4Z      +   RRR2 * T8Z )    * RRR2

    HHH = HH / THICK

    SS( 1) = ( HHH + HA1 ) * SA * 0.25
    SS(10) = ( HHH + HA2 ) * SA * 0.25
    SS(19) = ( HHH + HA3 ) * SA * 0.25
    SS(28) = ( HHH + HA4 ) * SA * 0.25

    SS(37) = ( HHH + HB1 ) * SB * 0.25
    SS(46) = ( HHH + HB2 ) * SB * 0.25
    SS(55) = ( HHH + HB3 ) * SB * 0.25
    SS(64) = ( HHH + HB4 ) * SB * 0.25

    SM = ( SA + SB ) * 0.5
    HH1 = ( HA1 + HB1 ) * 0.5
    HH2 = ( HA2 + HB2 ) * 0.5
    HH3 = ( HA3 + HB3 ) * 0.5
    HH4 = ( HA4 + HB4 ) * 0.5

    SS(33) = -( HHH + HH1 ) * SM * 0.25
    SS(42) = -( HHH + HH2 ) * SM * 0.25
    SS(51) = -( HHH + HH3 ) * SM * 0.25
    SS(60) = -( HHH + HH4 ) * SM * 0.25

    SS( 5) = SS(33)
    SS(14) = SS(42)
    SS(23) = SS(51)
    SS(32) = SS(60)

    !C    1                              ( 1- 1)    2  3  4  5  6  7  8
    !C    2  3                          9   (10- 3)   11 12 13 14 15 16
    !C    4  5  6                      17 18   (19- 6)   20 21 22 23 24
    !C    7  8  9 10          --->     25 26 27   (28-10)   29 30 31 32
    !C   11 12 13 14 15       --->     33 34 35 36   (37-15)   38 39 40
    !C   16 17 18 19 20 21             41 42 43 44 45   (46-21)   47 48
    !C   22 23 24 25 26 27 28          49 50 51 52 53 54   (55-28)   56
    !C   29 30 31 32 33 34 35 36       57 58 59 60 61 62 63   (64-36)

    return
  end subroutine heat_THERMAL_541

  !> CALCULATION SURFACE AREA FOR 4 POINTS
  subroutine heat_get_area ( XX,YY,ZZ,AA )
    !*
    !*  CALCULATION SURFACE AREA FOR 4 POINTS
    !*
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(4),YY(4),ZZ(4)
    dimension XG(2),H(4),HR(4),HS(4)
    !
    PI=4.0*atan(1.0)
    !
    XG(1) =-0.5773502691896258D0
    XG(2) =-XG(1)
    AA=0.0
    !
    do LX=1,2
      RI=XG(LX)
      do LY=1,2
        SI=XG(LY)
        RP=1.0+RI
        SP=1.0+SI
        RM=1.0-RI
        SM=1.0-SI
        !*INTERPOLATION FUNCTION
        H(1)=0.25*RP*SP
        H(2)=0.25*RM*SP
        H(3)=0.25*RM*SM
        H(4)=0.25*RP*SM
        !*DERIVATIVE OF INTERPOLATION FUNCTION
        !*  FOR R-COORDINATE
        HR(1)= .25*SP
        HR(2)=-.25*SP
        HR(3)=-.25*SM
        HR(4)= .25*SM
        !*  FOR S-COORDINATE
        HS(1)= .25*RP
        HS(2)= .25*RM
        HS(3)=-.25*RM
        HS(4)=-.25*RP

        !*JACOBI MATRIX
        XR=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
        XS=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
        YR=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
        YS=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
        ZR=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
        ZS=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)

        DET=(YR*ZS-ZR*YS)**2+(ZR*XS-XR*ZS)**2+(XR*YS-YR*XS)**2
        DET=sqrt(DET)
        !
        AA=AA+DET
        !
      enddo
    enddo

    return
  end subroutine heat_get_area

  !> CALCULATION 4 NODE SHELL ELEMENT
  subroutine heat_THERMAL_731 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 4 NODE SHELL ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !
    dimension XG(2), WGT(2), H(4), HR(4), HS(4)
    dimension COD(3, 4)
    dimension G1(3), G2(3), G3(3), E1(3), E2(3), E3(3), REF(3)
    dimension EN(3, 4), THE(3, 3), AMAT(3, 3), BV(3), WK(3)
    dimension DTDX(4), DTDY(4)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    data WGT/1.0D0,1.0D0/
    !
    !* SET COORDINATES
    !
    do I = 1, NN
      COD(1,I) = XX(I)
      COD(2,I) = YY(I)
      COD(3,I) = ZZ(I)
    enddo
    !
    !* SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM

    REF(1) = 0.25*( COD(1,2) + COD(1,3) - COD(1,1) - COD(1,4) )
    REF(2) = 0.25*( COD(2,2) + COD(2,3) - COD(2,1) - COD(2,4) )
    REF(3) = 0.25*( COD(3,2) + COD(3,3) - COD(3,1) - COD(3,4) )


    G1(1) = COD(1,1) - COD(1,2)
    G1(2) = COD(2,1) - COD(2,2)
    G1(3) = COD(3,1) - COD(3,2)
    G2(1) = COD(1,2) - COD(1,3)
    G2(2) = COD(2,2) - COD(2,3)
    G2(3) = COD(3,2) - COD(3,3)

    !* G3()=G1() X G2()

    G3(1) = G1(2)*G2(3) - G1(3)*G2(2)
    G3(2) = G1(3)*G2(1) - G1(1)*G2(3)
    G3(3) = G1(1)*G2(2) - G1(2)*G2(1)

    !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE

    XSUM = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )

    do I = 1, NN
      EN(1,I) = G3(1) / XSUM
      EN(2,I) = G3(2) / XSUM
      EN(3,I) = G3(3) / XSUM
    enddo

    !
    !*   LOOP FOR GAUSS INTEGRATION POINT
    !
    do IG3 = 1, 2
      TI = XG(IG3)

      do IG2 = 1, 2
        SI = XG(IG2)

        do IG1 = 1, 2
          RI = XG(IG1)

          RP = 1.0 + RI
          SP = 1.0 + SI
          RM = 1.0 - RI
          SM = 1.0 - SI

          H(1)  =  0.25*RM*SM
          H(2)  =  0.25*RP*SM
          H(3)  =  0.25*RP*SP
          H(4)  =  0.25*RM*SP

          HR(1) = -0.25*SM
          HR(2) =  0.25*SM
          HR(3) =  0.25*SP
          HR(4) = -0.25*SP

          HS(1) = -0.25*RM
          HS(2) = -0.25*RP
          HS(3) =  0.25*RP
          HS(4) =  0.25*RM
          !
          !* COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
          !
          do I = 1, 3

            G1(I) = 0.0
            G2(I) = 0.0
            G3(I) = 0.0

            do INOD = 1, NN
              VAR   = COD(I,INOD) + THICK*0.5*TI*EN(I,INOD)

              G1(I) = G1(I) + HR(INOD)*VAR
              G2(I) = G2(I) + HS(INOD)*VAR
              G3(I) = G3(I) + THICK*0.5*H(INOD)*EN(I,INOD)
            enddo

          enddo

          !*JACOBI MATRIX
          !
          XJ11 = G1(1)
          XJ12 = G1(2)
          XJ13 = G1(3)
          XJ21 = G2(1)
          XJ22 = G2(2)
          XJ23 = G2(3)
          XJ31 = G3(1)
          XJ32 = G3(2)
          XJ33 = G3(3)
          !
          !*DETERMINANT OF JACOBIAN
          !
          DET = XJ11*XJ22*XJ33 &
            + XJ12*XJ23*XJ31 &
            + XJ13*XJ21*XJ32 &
            - XJ13*XJ22*XJ31 &
            - XJ12*XJ21*XJ33 &
            - XJ11*XJ23*XJ32
          !
          !* INVERSION OF JACOBIAN
          !
          DUM = 1.0 / DET
          AMAT(1,1) = DUM*(  XJ22*XJ33 - XJ23*XJ32 )
          AMAT(2,1) = DUM*( -XJ21*XJ33 + XJ23*XJ31 )
          AMAT(3,1) = DUM*(  XJ21*XJ32 - XJ22*XJ31 )
          AMAT(1,2) = DUM*( -XJ12*XJ33 + XJ13*XJ32 )
          AMAT(2,2) = DUM*(  XJ11*XJ33 - XJ13*XJ31 )
          AMAT(3,2) = DUM*( -XJ11*XJ32 + XJ12*XJ31 )
          AMAT(1,3) = DUM*(  XJ12*XJ23 - XJ13*XJ22 )
          AMAT(2,3) = DUM*( -XJ11*XJ23 + XJ13*XJ21 )
          AMAT(3,3) = DUM*(  XJ11*XJ22 - XJ12*XJ21 )
          !
          !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
          !
          XSUM  = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )
          E3(1) = G3(1) / XSUM
          E3(2) = G3(2) / XSUM
          E3(3) = G3(3) / XSUM
          !
          E2(1) = -REF(2)*E3(3) + REF(3)*E3(2)
          E2(2) = -REF(3)*E3(1) + REF(1)*E3(3)
          E2(3) = -REF(1)*E3(2) + REF(2)*E3(1)
          !
          E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
          E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
          E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
          !
          XSUM = dsqrt(E2(1)**2 + E2(2)**2 + E2(3)**2)

          if ( XSUM .GT. 1.E-15 ) then

            E2(1) = E2(1) / XSUM
            E2(2) = E2(2) / XSUM
            E2(3) = E2(3) / XSUM
            !
            E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
            E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
            E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
            !
            XSUM = dsqrt( E1(1)**2 + E1(2)**2 + E1(3)**2 )

            E1(1) = E1(1) / XSUM
            E1(2) = E1(2) / XSUM
            E1(3) = E1(3) / XSUM

          else

            E1(1) =  0.D0
            E1(2) =  0.D0
            E1(3) = -1.D0
            E2(1) =  0.D0
            E2(2) =  1.D0
            E2(3) =  0.D0

          end if
          !
          THE(1,1) = E1(1)
          THE(1,2) = E1(2)
          THE(1,3) = E1(3)
          THE(2,1) = E2(1)
          THE(2,2) = E2(2)
          THE(2,3) = E2(3)
          THE(3,1) = E3(1)
          THE(3,2) = E3(2)
          THE(3,3) = E3(3)
          !
          do I = 1, NN

            BV(1) = HR(I)
            BV(2) = HS(I)
            BV(3) = 0.0

            WK(1) = AMAT(1,1)*BV(1) &
              + AMAT(1,2)*BV(2) &
              + AMAT(1,3)*BV(3)
            WK(2) = AMAT(2,1)*BV(1) &
              + AMAT(2,2)*BV(2) &
              + AMAT(2,3)*BV(3)
            WK(3) = AMAT(3,1)*BV(1) &
              + AMAT(3,2)*BV(2) &
              + AMAT(3,3)*BV(3)
            BV(1) = THE(1,1)*WK(1) &
              + THE(1,2)*WK(2) &
              + THE(1,3)*WK(3)
            BV(2) = THE(2,1)*WK(1) &
              + THE(2,2)*WK(2) &
              + THE(2,3)*WK(3)
            BV(3) = THE(3,1)*WK(1) &
              + THE(3,2)*WK(2) &
              + THE(3,3)*WK(3)

            DTDX(I) = BV(1)
            DTDY(I) = BV(2)

          enddo

          !
          !*CONDUCTIVITY AT CURRENT TEMPERATURE

          CTEMP=0.0
          do I = 1, NN
            CTEMP = CTEMP + H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY( CTEMP,IMAT,CC,ntab,temp,funcA,funcB )
          !
          !* SET INTEGRATION WEIGHT
          !
          VALX = CC(1)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          VALY = CC(2)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          !
          SS( 1) = SS( 1) + DTDX(1)*DTDX(1)*VALX  &
            + DTDY(1)*DTDY(1)*VALY
          SS( 2) = SS( 2) + DTDX(1)*DTDX(2)*VALX  &
            + DTDY(1)*DTDY(2)*VALY
          SS( 3) = SS( 3) + DTDX(1)*DTDX(3)*VALX  &
            + DTDY(1)*DTDY(3)*VALY
          SS( 4) = SS( 4) + DTDX(1)*DTDX(4)*VALX  &
            + DTDY(1)*DTDY(4)*VALY
          SS( 6) = SS( 6) + DTDX(2)*DTDX(2)*VALX  &
            + DTDY(2)*DTDY(2)*VALY
          SS( 7) = SS( 7) + DTDX(2)*DTDX(3)*VALX  &
            + DTDY(2)*DTDY(3)*VALY
          SS( 8) = SS( 8) + DTDX(2)*DTDX(4)*VALX  &
            + DTDY(2)*DTDY(4)*VALY
          SS(11) = SS(11) + DTDX(3)*DTDX(3)*VALX  &
            + DTDY(3)*DTDY(3)*VALY
          SS(12) = SS(12) + DTDX(3)*DTDX(4)*VALX  &
            + DTDY(3)*DTDY(4)*VALY
          SS(16) = SS(16) + DTDX(4)*DTDX(4)*VALX  &
            + DTDY(4)*DTDY(4)*VALY

        enddo
      enddo
    enddo
    !
    SS( 3) = SS( 3) + SS( 4)
    SS( 7) = SS( 7) + SS( 8)
    SS(11) = SS(11) + SS(16) + 2.0*SS(12)
    SS( 4) = 0.0D0
    SS( 8) = 0.0D0
    SS(12) = 0.0D0
    SS(16) = 0.0D0

    SS( 5) = SS( 2)
    SS( 9) = SS( 3)
    SS(10) = SS( 7)
    SS(13) = SS( 4)
    SS(14) = SS( 8)
    SS(15) = SS(12)
    !
    return
  end subroutine heat_THERMAL_731

  !> CALCULATION 4 NODE SHELL ELEMENT
  subroutine heat_THERMAL_741 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 4 NODE SHELL ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !
    dimension XG(2), WGT(2), H(4), HR(4), HS(4)
    dimension COD(3, 4)
    dimension G1(3), G2(3), G3(3), E1(3), E2(3), E3(3), REF(3)
    dimension EN(3, 4), THE(3, 3), AMAT(3, 3), BV(3), WK(3)
    dimension DTDX(4), DTDY(4)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    data WGT/1.0D0,1.0D0/

    COD = 0.0d0
    !
    !* SET COORDINATES
    !
    do I = 1, NN
      COD(1,I) = XX(I)
      COD(2,I) = YY(I)
      COD(3,I) = ZZ(I)
    enddo
    !
    !* SET REFFRENSE VECTOR TO DETERMINE LOCAL COORDINATE SYSTEM

    REF(1) = 0.25*( COD(1,2) + COD(1,3) - COD(1,1) - COD(1,4) )
    REF(2) = 0.25*( COD(2,2) + COD(2,3) - COD(2,1) - COD(2,4) )
    REF(3) = 0.25*( COD(3,2) + COD(3,3) - COD(3,1) - COD(3,4) )


    do I = 1, NN

      if     ( I .EQ. 1 ) then
        G1(1) = -0.5*COD(1,1) + 0.5*COD(1,2)
        G1(2) = -0.5*COD(2,1) + 0.5*COD(2,2)
        G1(3) = -0.5*COD(3,1) + 0.5*COD(3,2)
        G2(1) = -0.5*COD(1,1) + 0.5*COD(1,4)
        G2(2) = -0.5*COD(2,1) + 0.5*COD(2,4)
        G2(3) = -0.5*COD(3,1) + 0.5*COD(3,4)
      else if( I .EQ. 2 ) then
        G1(1) = -0.5*COD(1,1) + 0.5*COD(1,2)
        G1(2) = -0.5*COD(2,1) + 0.5*COD(2,2)
        G1(3) = -0.5*COD(3,1) + 0.5*COD(3,2)
        G2(1) = -0.5*COD(1,2) + 0.5*COD(1,3)
        G2(2) = -0.5*COD(2,2) + 0.5*COD(2,3)
        G2(3) = -0.5*COD(3,2) + 0.5*COD(3,3)
      else if( I .EQ. 3 ) then
        G1(1) =  0.5*COD(1,3) - 0.5*COD(1,4)
        G1(2) =  0.5*COD(2,3) - 0.5*COD(2,4)
        G1(3) =  0.5*COD(3,3) - 0.5*COD(3,4)
        G2(1) = -0.5*COD(1,2) + 0.5*COD(1,3)
        G2(2) = -0.5*COD(2,2) + 0.5*COD(2,3)
        G2(3) = -0.5*COD(3,2) + 0.5*COD(3,3)
      else if( I .EQ. 4 ) then
        G1(1) =  0.5*COD(1,3) - 0.5*COD(1,4)
        G1(2) =  0.5*COD(2,3) - 0.5*COD(2,4)
        G1(3) =  0.5*COD(3,3) - 0.5*COD(3,4)
        G2(1) = -0.5*COD(1,1) + 0.5*COD(1,4)
        G2(2) = -0.5*COD(2,1) + 0.5*COD(2,4)
        G2(3) = -0.5*COD(3,1) + 0.5*COD(3,4)
      end if

      !* G3()=G1() X G2()

      G3(1) = G1(2)*G2(3) - G1(3)*G2(2)
      G3(2) = G1(3)*G2(1) - G1(1)*G2(3)
      G3(3) = G1(1)*G2(2) - G1(2)*G2(1)

      !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE

      XSUM = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )
      E3(1)   = G3(1) / XSUM
      E3(2)   = G3(2) / XSUM
      E3(3)   = G3(3) / XSUM
      EN(1,I) = E3(1)
      EN(2,I) = E3(2)
      EN(3,I) = E3(3)

    enddo

    !
    !*   LOOP FOR GAUSS INTEGRATION POINT
    !
    do IG3 = 1, 2
      TI = XG(IG3)

      do IG2 = 1, 2
        SI = XG(IG2)

        do IG1 = 1, 2
          RI = XG(IG1)

          RP = 1.0 + RI
          SP = 1.0 + SI
          RM = 1.0 - RI
          SM = 1.0 - SI

          H(1)  =  0.25*RM*SM
          H(2)  =  0.25*RP*SM
          H(3)  =  0.25*RP*SP
          H(4)  =  0.25*RM*SP

          HR(1) = -0.25*SM
          HR(2) =  0.25*SM
          HR(3) =  0.25*SP
          HR(4) = -0.25*SP

          HS(1) = -0.25*RM
          HS(2) = -0.25*RP
          HS(3) =  0.25*RP
          HS(4) =  0.25*RM
          !
          !* COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
          !
          do I = 1, 3
            G1(I) = 0.0
            G2(I) = 0.0
            G3(I) = 0.0

            do INOD = 1, NN
              VAR   = COD(I,INOD) + THICK*0.5*TI*EN(I,INOD)
              G1(I) = G1(I) + HR(INOD)*VAR
              G2(I) = G2(I) + HS(INOD)*VAR
              G3(I) = G3(I) + THICK*0.5*H(INOD)*EN(I,INOD)
            enddo

          enddo

          !
          !*JACOBI MATRIX
          !
          XJ11 = G1(1)
          XJ12 = G1(2)
          XJ13 = G1(3)
          XJ21 = G2(1)
          XJ22 = G2(2)
          XJ23 = G2(3)
          XJ31 = G3(1)
          XJ32 = G3(2)
          XJ33 = G3(3)
          !
          !*DETERMINANT OF JACOBIAN
          !
          DET = XJ11*XJ22*XJ33 &
            + XJ12*XJ23*XJ31 &
            + XJ13*XJ21*XJ32 &
            - XJ13*XJ22*XJ31 &
            - XJ12*XJ21*XJ33 &
            - XJ11*XJ23*XJ32
          !
          !* INVERSION OF JACOBIAN
          !
          DUM = 1.0 / DET
          AMAT(1,1) = DUM*(  XJ22*XJ33 - XJ23*XJ32 )
          AMAT(2,1) = DUM*( -XJ21*XJ33 + XJ23*XJ31 )
          AMAT(3,1) = DUM*(  XJ21*XJ32 - XJ22*XJ31 )
          AMAT(1,2) = DUM*( -XJ12*XJ33 + XJ13*XJ32 )
          AMAT(2,2) = DUM*(  XJ11*XJ33 - XJ13*XJ31 )
          AMAT(3,2) = DUM*( -XJ11*XJ32 + XJ12*XJ31 )
          AMAT(1,3) = DUM*(  XJ12*XJ23 - XJ13*XJ22 )
          AMAT(2,3) = DUM*( -XJ11*XJ23 + XJ13*XJ21 )
          AMAT(3,3) = DUM*(  XJ11*XJ22 - XJ12*XJ21 )
          !
          !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
          !
          XSUM  = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )
          E3(1) = G3(1) / XSUM
          E3(2) = G3(2) / XSUM
          E3(3) = G3(3) / XSUM
          !
          E2(1) = -REF(2)*E3(3) + REF(3)*E3(2)
          E2(2) = -REF(3)*E3(1) + REF(1)*E3(3)
          E2(3) = -REF(1)*E3(2) + REF(2)*E3(1)
          !
          E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
          E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
          E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
          !
          XSUM = dsqrt( E2(1)**2 + E2(2)**2 + E2(3)**2 )

          if ( XSUM .GT. 1.E-15 ) then

            E2(1) = E2(1) / XSUM
            E2(2) = E2(2) / XSUM
            E2(3) = E2(3) / XSUM
            !
            E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
            E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
            E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
            !
            XSUM = dsqrt( E1(1)**2 + E1(2)**2 + E1(3)**2 )

            E1(1) = E1(1) / XSUM
            E1(2) = E1(2) / XSUM
            E1(3) = E1(3) / XSUM

          else

            E1(1) =  0.D0
            E1(2) =  0.D0
            E1(3) = -1.D0
            E2(1) =  0.D0
            E2(2) =  1.D0
            E2(3) =  0.D0

          end if
          !
          THE(1,1) = E1(1)
          THE(1,2) = E1(2)
          THE(1,3) = E1(3)
          THE(2,1) = E2(1)
          THE(2,2) = E2(2)
          THE(2,3) = E2(3)
          THE(3,1) = E3(1)
          THE(3,2) = E3(2)
          THE(3,3) = E3(3)
          !
          do I = 1, NN

            BV(1) = HR(I)
            BV(2) = HS(I)
            BV(3) = 0.0

            WK(1) = AMAT(1,1)*BV(1) &
              + AMAT(1,2)*BV(2) &
              + AMAT(1,3)*BV(3)
            WK(2) = AMAT(2,1)*BV(1) &
              + AMAT(2,2)*BV(2) &
              + AMAT(2,3)*BV(3)
            WK(3) = AMAT(3,1)*BV(1) &
              + AMAT(3,2)*BV(2) &
              + AMAT(3,3)*BV(3)
            BV(1) = THE(1,1)*WK(1) &
              + THE(1,2)*WK(2) &
              + THE(1,3)*WK(3)
            BV(2) = THE(2,1)*WK(1) &
              + THE(2,2)*WK(2) &
              + THE(2,3)*WK(3)
            BV(3) = THE(3,1)*WK(1) &
              + THE(3,2)*WK(2) &
              + THE(3,3)*WK(3)

            DTDX(I) = BV(1)
            DTDY(I) = BV(2)

          enddo


          !
          !*CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I = 1, NN
            CTEMP = CTEMP + H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY( CTEMP,IMAT,CC,ntab,temp,funcA,funcB )
          !
          !* SET INTEGRATION WEIGHT
          !
          VALX = CC(1)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          VALY = CC(2)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          !
          SS( 1) = SS( 1) + DTDX(1)*DTDX(1)*VALX  &
            + DTDY(1)*DTDY(1)*VALY
          SS( 2) = SS( 2) + DTDX(1)*DTDX(2)*VALX  &
            + DTDY(1)*DTDY(2)*VALY
          SS( 3) = SS( 3) + DTDX(1)*DTDX(3)*VALX  &
            + DTDY(1)*DTDY(3)*VALY
          SS( 4) = SS( 4) + DTDX(1)*DTDX(4)*VALX  &
            + DTDY(1)*DTDY(4)*VALY
          SS( 6) = SS( 6) + DTDX(2)*DTDX(2)*VALX  &
            + DTDY(2)*DTDY(2)*VALY
          SS( 7) = SS( 7) + DTDX(2)*DTDX(3)*VALX  &
            + DTDY(2)*DTDY(3)*VALY
          SS( 8) = SS( 8) + DTDX(2)*DTDX(4)*VALX  &
            + DTDY(2)*DTDY(4)*VALY
          SS(11) = SS(11) + DTDX(3)*DTDX(3)*VALX  &
            + DTDY(3)*DTDY(3)*VALY
          SS(12) = SS(12) + DTDX(3)*DTDX(4)*VALX  &
            + DTDY(3)*DTDY(4)*VALY
          SS(16) = SS(16) + DTDX(4)*DTDX(4)*VALX  &
            + DTDY(4)*DTDY(4)*VALY

        enddo
      enddo
    enddo
    !
    SS( 5) = SS( 2)
    SS( 9) = SS( 3)
    SS(10) = SS( 7)
    SS(13) = SS( 4)
    SS(14) = SS( 8)
    SS(15) = SS(12)
    !
    return
  end subroutine heat_THERMAL_741

  !> CALCULATION 2D 6 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_232 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !C*--------------------------------------------------------------------*
    !C*
    !C*CALCULATION 2D 6 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(3), WGT(3), H(6), HL1(6), HL2(6), HL3(6)
    dimension BX(6), BY(6), BZ(6), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
    do L2=1,3
      XL2=XG(L2)
      X2 =(XL2+1.0)*0.5
      do L1=1,3
        XL1=XG(L1)
        X1=0.5*(1.0-X2)*(XL1+1.0)
        !C*INTERPOLATION FUNCTION
        X3=1.0-X1-X2
        H(1)= X1*(2.0*X1-1.)
        H(2)= X2*(2.0*X2-1.)
        H(3)= X3*(2.0*X3-1.)
        H(4)= 4.0*X1*X2
        H(5)= 4.0*X2*X3
        H(6)= 4.0*X1*X3
        !C*DERIVATIVE OF INTERPOLATION FUNCTION
        !C*FOR L1-COORDINATE
        HL1(1)=4.0*X1-1.0
        HL1(2)= 0.0
        HL1(3)= 0.0
        HL1(4)= 4.0*X2
        HL1(5)= 0.0
        HL1(6)= 4.0*X3
        !C*FOR L2-COORDINATE
        HL2(1)= 0.0
        HL2(2)= 4.0*X2-1.0
        HL2(3)= 0.0
        HL2(4)= 4.0*X1
        HL2(5)= 4.0*X3
        HL2(6)= 0.0
        !C*FOR L3-COORDINATE
        HL3(1)= 0.0
        HL3(2)= 0.0
        HL3(3)= 4.0*X3-1.0
        HL3(4)= 0.0
        HL3(5)= 4.0*X2
        HL3(6)= 4.0*X1
        !C*JACOBI MATRIX
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
        DUM=1.0/DET
        XJI11= XJ22*DUM
        XJI12=-XJ12*DUM
        XJI21=-XJ21*DUM
        XJI22= XJ11*DUM
        do J=1, NN
          BX(J)=XJI11*(HL1(J)-HL3(J))+XJI12*(HL2(J)-HL3(J))
          BY(J)=XJI21*(HL1(J)-HL3(J))+XJI22*(HL2(J)-HL3(J))
        enddo
        !C*CONDUCTIVITY AT CURRENT TEMPERATURE
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
        !C*WEIGT VALUE AT GAUSSIAN POINT
        WGX=CC(1)*WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        WGY=CC(2)*WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        !C
        IJ = 0
        do I = 1, NN
          do J = 1, NN
            IJ = IJ + 1
            SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX &
              +BY(I)*BY(J)*WGY
          enddo
        enddo
        !C
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_232

  !> CALCULATION 2D 8 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_242 (NN,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB)
    !C*
    !C*CALCULATION 2D 8 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension XG(3), WGT(3), H(8), HR(8), HS(8), HT(8)
    dimension BX(8), BY(8), BZ(8), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
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
        DUM=1.0/DET
        XJI11= XJ22*DUM
        XJI12=-XJ12*DUM
        XJI21=-XJ21*DUM
        XJI22= XJ11*DUM
        do J=1, NN
          BX(J)=XJI11*HR(J)+XJI12*HS(J)
          BY(J)=XJI21*HR(J)+XJI22*HS(J)
        enddo
        !*CONDUCTIVITY AT CURRENT TEMPERATURE
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
        !*WEIGT VALUE AT GAUSSIAN POINT
        WGX=CC(1)*WGT(LX)*WGT(LY)*DET*THICK
        WGY=CC(2)*WGT(LX)*WGT(LY)*DET*THICK
        !C
        IJ = 0
        do I = 1, NN
          do J = 1, NN
            IJ = IJ + 1
            SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX &
              +BY(I)*BY(J)*WGY
          enddo
        enddo
        !C
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_242

  !> CALCULATION 3D 10 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_342 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !C*
    !C* CALCULATION 3D 10 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension H(10), HL1(10), HL2(10), HL3(10), HL4(10)
    dimension BX(10), BY(10), BZ(10), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    real(kind=kreal) XG(3), WGT(3)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C* LOOP FOR INTEGRATION POINTS
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
          !C
          !C* DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C* INVERSION OF JACOBIAN
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
          !C
          !C*CONDUCTIVITY AT CURRENT TEMPERATURE
          !C
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY( CTEMP,IMAT,CC,ntab,temp,funcA,funcB )
          !C
          !C* WEIGT VALUE AT GAUSSIAN POINT
          DET = -DET
          WGX=CC(1)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          WGY=CC(2)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          WGZ=CC(3)*WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          !C
          do J=1, NN
            BX(J)=XJI11*(HL1(J)-HL4(J)) &
              +XJI12*(HL2(J)-HL4(J)) &
              +XJI13*(HL3(J)-HL4(J))
            BY(J)=XJI21*(HL1(J)-HL4(J)) &
              +XJI22*(HL2(J)-HL4(J)) &
              +XJI23*(HL3(J)-HL4(J))
            BZ(J)=XJI31*(HL1(J)-HL4(J)) &
              +XJI32*(HL2(J)-HL4(J)) &
              +XJI33*(HL3(J)-HL4(J))
          enddo
          !C
          IJ = 0
          do I = 1, NN
            do J = 1, NN
              IJ = IJ + 1
              SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX &
                +BY(I)*BY(J)*WGY &
                +BZ(I)*BZ(J)*WGZ
            enddo
          enddo
          !C
        enddo
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_342

  !> CALCULATION 3D 15 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_352 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !*
    !* CALCULATION 3D 15 NODE CONDUCTANCE ELEMENT
    !*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension H(15), HL1(15), HL2(15), HL3(15), HZ(15), BX(15), BY(15), BZ(15), CC(3)
    dimension temp(*), funcA(*), funcB(*)
    real(kind=kreal) XG(3), WGT(3)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C* LOOP FOR INTEGRATION POINTS
    do LZ=1,3
      ZI=XG(LZ)
      do L2=1,3
        XL2=XG(L2)
        X2 =(XL2+1.0)*0.5
        do L1=1,3
          XL1=XG(L1)
          X1=0.5*(1.0-X2)*(XL1+1.0)
          !C* INTERPOLATION FUNCTION
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
          !C* DERIVATIVE OF INTERPOLATION FUNCTION
          !C* FOR L1-COORDINATE
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
          !C* FOR L2-COORDINATE
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
          !C* FOR L3-COORDINATE
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
          !C* FOR Z-COORDINATE
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
          !C
          !C*JACOBI MATRIX
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
            XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
            XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
            XJ31=XJ31+HZ(I)*XX(I)
            XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
            XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
            XJ32=XJ32+HZ(I)*YY(I)
            XJ13=XJ13+(HL1(I)-HL3(I))*ZZ(I)
            XJ23=XJ23+(HL2(I)-HL3(I))*ZZ(I)
            XJ33=XJ33+HZ(I)*ZZ(I)
          enddo
          !C
          !C* DETERMINANT OF JACOBIAN
          DET = XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C
          !C* CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I=1,NN
            CTEMP = CTEMP+H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
          !C* WEIGT VALUE AT GAUSSIAN POINT
          WGX = CC(1)*WGT(L1)*WGT(L2)*WGT(LZ)*DET*(1.0-X2)*0.25
          WGY = CC(2)*WGT(L1)*WGT(L2)*WGT(LZ)*DET*(1.0-X2)*0.25
          WGZ = CC(3)*WGT(L1)*WGT(L2)*WGT(LZ)*DET*(1.0-X2)*0.25
          !C
          !C* INVERSION OF JACOBIAN
          DUM   = 1.0/DET
          XJI11 = DUM*( XJ22*XJ33-XJ23*XJ32)
          XJI21 = DUM*(-XJ21*XJ33+XJ23*XJ31)
          XJI31 = DUM*( XJ21*XJ32-XJ22*XJ31)
          XJI12 = DUM*(-XJ12*XJ33+XJ13*XJ32)
          XJI22 = DUM*( XJ11*XJ33-XJ13*XJ31)
          XJI32 = DUM*(-XJ11*XJ32+XJ12*XJ31)
          XJI13 = DUM*( XJ12*XJ23-XJ13*XJ22)
          XJI23 = DUM*(-XJ11*XJ23+XJ13*XJ21)
          XJI33 = DUM*( XJ11*XJ22-XJ12*XJ21)
          !C
          do J=1,NN
            BX(J) = XJI11*(HL1(J)-HL3(J)) &
              +XJI12*(HL2(J)-HL3(J))+XJI13*HZ(J)
            BY(J) = XJI21*(HL1(J)-HL3(J)) &
              +XJI22*(HL2(J)-HL3(J))+XJI23*HZ(J)
            BZ(J) = XJI31*(HL1(J)-HL3(J)) &
              +XJI32*(HL2(J)-HL3(J))+XJI33*HZ(J)
          enddo
          !C
          IJ = 0
          do I = 1, NN
            do J = 1, NN
              IJ = IJ + 1
              SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX &
                +BY(I)*BY(J)*WGY &
                +BZ(I)*BZ(J)*WGZ
            enddo
          enddo
          !C
        enddo
      enddo
    enddo
    !C
    return
  end subroutine heat_THERMAL_352

  !> CALCULATION 3D 20 NODE CONDUCTANCE ELEMENT
  subroutine heat_THERMAL_362 (NN,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB)
    !C*
    !C* CALCULATION 3D 20 NODE CONDUCTANCE ELEMENT
    !C*
    use hecmw
    implicit real(kind=kreal) (A - H, O - Z)
    dimension XX(NN), YY(NN), ZZ(NN), TT(NN), SS(NN*NN)
    dimension H(20), HR(20), HS(20), HT(20)
    dimension BX(20), BY(20), BZ(20), CC(3)
    dimension temp(ntab), funcA(ntab + 1), funcB(ntab + 1)
    real(kind=kreal) XG(3), WGT(3)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C* LOOP FOR INTEGRATION POINTS
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
          !C* INTERPOLATION FUNCTION
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
          !C* DERIVATIVE OF INTERPOLATION FUNCTION
          !C* FOR R-COORDINATE
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
          !C* FOR S-COORDINATE
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
          !C* FOR T-COORDINATE
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
          !C*  JACOBI MATRIX
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
          !C* DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33 &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C
          !C* INVERSION OF JACOBIAN
          !C
          DET = -DET
          DUM=1.0/DET
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
          do J=1, NN
            BX(J)=XJI11*HR(J)+XJI12*HS(J)+XJI13*HT(J)
            BY(J)=XJI21*HR(J)+XJI22*HS(J)+XJI23*HT(J)
            BZ(J)=XJI31*HR(J)+XJI32*HS(J)+XJI33*HT(J)
          enddo
          !C
          !C* CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+H(I)*TT(I)
          enddo
          call heat_GET_CONDUCTIVITY(CTEMP,IMAT,CC,ntab,temp,funcA,funcB)
          !C
          !C* WEIGT VALUE AT GAUSSIAN POINT
          WGX=-CC(1)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          WGY=-CC(2)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          WGZ=-CC(3)*WGT(LX)*WGT(LY)*WGT(LZ)*DET
          !C
          IJ = 0
          do I = 1, NN
            do J = 1, NN
              IJ = IJ + 1
              SS(IJ)=SS(IJ)+BX(I)*BX(J)*WGX &
                +BY(I)*BY(J)*WGY &
                +BZ(I)*BZ(J)*WGZ
            enddo
          enddo
          !C
        enddo
      enddo
    enddo
    !
    return
  end subroutine heat_THERMAL_362
end module m_heat_lib_thermal
