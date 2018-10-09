!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines for calculating heat capacity matrix
!! for various elements
module m_heat_LIB_CAPACITY
contains
  !***********************************************************************
  !  CAPACITY_111 ( NN,XX,YY,ZZ,TT,IMAT,ASECT,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_231 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_241 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_341 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_351 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_361 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_731 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_741 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_232 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_242 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS
  !                                       ,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_342 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_352 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !  CAPACITY_362 ( NN,XX,YY,ZZ,TT,IMAT,SS,ntab1,temp1,funcA1,funcB1
  !                                       ,ntab2,temp2,funcA2,funcB2 )
  !***********************************************************************

  !> Build up the heat capacity matrix of line element 111
  subroutine heat_CAPACITY_111 ( NN,XX,YY,ZZ,TT,IMAT,ASECT,SS &
      ,ntab1,temp1,funcA1,funcB1  &
      ,ntab2,temp2,funcA2,funcB2 )
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    !
    DX = XX(2) - XX(1)
    DY = YY(2) - YY(1)
    DZ = ZZ(2) - ZZ(1)
    AL = dsqrt( DX*DX + DY*DY + DZ*DZ )
    VV = ASECT * AL

    TEMP = ( TT(1) + TT(2) ) * 0.5d0
    call heat_GET_CP     (TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1)
    call heat_GET_DENSITY(TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2)

    SS(1) = VV*SPE(1)*DEN(1)*0.5
    SS(2) = SS(1)
    write(*,*) 'VV,SPE,DEN,SS,'
    write(*,*) VV,SPE(1),DEN(1),SS(1)
    !
    return

  end subroutine heat_CAPACITY_111

  !> Build up the heat capacity matrix of element 231
  subroutine heat_CAPACITY_231 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !C*--------------------------------------------------------------------*
    !C*
    !C* HEAT CAPACITY FOR 231
    !C*
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension XG(2),WGT(2),H(3),HL1(3),HL2(3),HL3(3)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
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
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        WG=WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        do I=1,NN
          SS(I)=SS(I)+WG*SPE(1)*DEN(1)*H(I)
        enddo
        !C
      enddo
    enddo
    !C
    return
    !C
  end subroutine heat_CAPACITY_231

  !> Build up the heat capacity matrix of element 241
  subroutine heat_CAPACITY_241 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !C*
    !C* HEAT CAPACITY FOR 241
    !C*
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension XG(2),WGT(2),H(4),HR(4),HS(4),HT(4)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    data XG/-0.5773502691896D0,0.5773502691896D0/
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
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        !*WEIGT VALUE AT GAUSSIAN POINT
        WG=WGT(LX)*WGT(LY)*DET*THICK
        do I=1,NN
          SS(I)=SS(I)+WG*SPE(1)*DEN(1)*H(I)
        enddo
      enddo
    enddo

    return

  end subroutine heat_CAPACITY_241

  !> Build up the heat capacity matrix of element 341
  subroutine heat_CAPACITY_341 ( NN,XX,YY,ZZ,TT,IMAT,SS     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !*
    !
    ! HEAT CAPACITY FOR 341
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension XG(2),WGT(2),H(4),HR(4),HS(4),HT(4)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    !
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !
    VV = 0.0
    !
    do L3=1,2
      XL3=XG(L3)
      X3=(XL3+1.0)*0.5
      do L2=1,2
        XL2=XG(L2)
        X2=(1.0-X3)*(XL2+1.0)*0.5
        do L1 =1,2
          XL1=XG(L1)
          X1=(1.0-X2-X3 )*(XL1+1.0)*0.5

          !*INTERPOLATION FUNCTION

          H(1)= X1
          H(2)= X2
          H(3)= X3
          H(4)= 1.0 - X1 - X2 - X3
          !
          !*DERIVATIVE OF INTERPOLATION FUNCTION

          !*  FOR L1-COORDINATE
          HR(1)= 1.0
          HR(2)= 0.0
          HR(3)= 0.0
          HR(4)=-1.0

          !*  FOR L2-COORDINATE
          HS(1)= 0.0
          HS(2)= 1.0
          HS(3)= 0.0
          HS(4)=-1.0

          !*  FOR ZETA-COORDINATE
          HT(1)= 0.0
          HT(2)= 0.0
          HT(3)= 1.0
          HT(4)=-1.0
          !
          !*JACOBI MATRIX
          !
          XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
          XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
          XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4)
          XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
          XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
          XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4)
          XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
          XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)
          XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4)
          !*DETERMINANT OF JACOBIAN
          DET = XJ11*XJ22*XJ33   &
            + XJ12*XJ23*XJ31  &
            + XJ13*XJ21*XJ32  &
            - XJ13*XJ22*XJ31  &
            - XJ12*XJ21*XJ33  &
            - XJ11*XJ23*XJ32
          DET = -DET
          !C
          !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
          !C
          TEMP=0.0
          do I=1,NN
            TEMP = TT(I)
            call heat_GET_CP     ( TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
            call heat_GET_DENSITY( TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
          enddo
          WG=WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X2-X3)*(1.0-X3)*0.125
          do I=1,NN
            SS(I)=SS(I)+WG*SPE(1)*DEN(1)*H(I)
          enddo
          !
        enddo
      enddo
    enddo
    !
    return
    !
  end subroutine heat_CAPACITY_341

  !> Build up the heat capacity matrix of element 351
  subroutine heat_CAPACITY_351 ( NN,XX,YY,ZZ,TT,IMAT,S0     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 351
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension XG(2),WGT(2),XG1(3),XG2(3),WGT1(3),H(6),HR(6),HS(6),HT(6)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    data XG   /-0.5773502691896D0,0.5773502691896D0/
    data WGT  /1.0D0,1.0D0/
    data XG1  /0.6666666667D0,0.16666666667D0,0.16666666667D0/
    data XG2  /0.1666666667D0,0.66666666667D0,0.16666666667D0/
    data WGT1 /0.1666666667D0,0.16666666667D0,0.16666666667D0/
    !
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
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        TEMP=0.0
        do I=1,NN
          TEMP = TEMP+TT(I)*H(I)
        enddo

        call heat_GET_CP     ( TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        WG=WGT(LZ)*WGT1(L12)*DET
        do I=1,NN
          S0(I) = S0(I)+WG*SPE(1)*DEN(1)*H(I)
        enddo
        !
      enddo
    enddo
    !
    return

  end subroutine heat_CAPACITY_351

  !> Build up the heat capacity matrix of element 361
  subroutine heat_CAPACITY_361 ( NN,XX,YY,ZZ,TT,IMAT,S0     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 361
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension XG(2),WGT(2),H(8),HR(8),HS(8),HT(8)
    dimension BX(8),BY(8),BZ(8),CC(3)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    !
    ! SET GAUUUSIAN INTEGRATION PARAMETER
    !
    XG(1) =-0.5773502691896258
    XG(2) =-XG(1)
    WGT(1)= 1.0D0
    WGT(2)= WGT(1)
    !
    do LX=1,2
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
          !INTERPOLATION FUNCTION
          H(1)=0.125*RM*SM*TM
          H(2)=0.125*RP*SM*TM
          H(3)=0.125*RP*SP*TM
          H(4)=0.125*RM*SP*TM
          H(5)=0.125*RM*SM*TP
          H(6)=0.125*RP*SM*TP
          H(7)=0.125*RP*SP*TP
          H(8)=0.125*RM*SP*TP
          !DERIVATIVE OF INTERPOLATION FUNCTION
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
          !JACOBI MATRIX
          XJ11=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4) &
            +HR(5)*XX(5)+HR(6)*XX(6)+HR(7)*XX(7)+HR(8)*XX(8)
          XJ21=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4) &
            +HS(5)*XX(5)+HS(6)*XX(6)+HS(7)*XX(7)+HS(8)*XX(8)
          XJ31=HT(1)*XX(1)+HT(2)*XX(2)+HT(3)*XX(3)+HT(4)*XX(4) &
            +HT(5)*XX(5)+HT(6)*XX(6)+HT(7)*XX(7)+HT(8)*XX(8)
          XJ12=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4) &
            +HR(5)*YY(5)+HR(6)*YY(6)+HR(7)*YY(7)+HR(8)*YY(8)
          XJ22=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4) &
            +HS(5)*YY(5)+HS(6)*YY(6)+HS(7)*YY(7)+HS(8)*YY(8)
          XJ32=HT(1)*YY(1)+HT(2)*YY(2)+HT(3)*YY(3)+HT(4)*YY(4) &
            +HT(5)*YY(5)+HT(6)*YY(6)+HT(7)*YY(7)+HT(8)*YY(8)
          XJ13=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4) &
            +HR(5)*ZZ(5)+HR(6)*ZZ(6)+HR(7)*ZZ(7)+HR(8)*ZZ(8)
          XJ23=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4) &
            +HS(5)*ZZ(5)+HS(6)*ZZ(6)+HS(7)*ZZ(7)+HS(8)*ZZ(8)
          XJ33=HT(1)*ZZ(1)+HT(2)*ZZ(2)+HT(3)*ZZ(3)+HT(4)*ZZ(4) &
            +HT(5)*ZZ(5)+HT(6)*ZZ(6)+HT(7)*ZZ(7)+HT(8)*ZZ(8)
          !DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33   &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C
          !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
          !C
          TEMP= 0.0
          do I=1, NN
            TEMP =TEMP+TT(I)*H(I)
          enddo
          call heat_GET_CP     (TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1)
          call heat_GET_DENSITY(TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2)
          WG=WGT(LX)*WGT(LY)*WGT(LZ)*DET
          do I=1, NN
            S0(I)=S0(I)+WG*SPE(1)*DEN(1)*H(I)
          enddo
          !
        enddo
      enddo
    enddo
    !
    return

  end subroutine heat_CAPACITY_361


  !> Build up the heat capacity matrix of element 731
  subroutine heat_CAPACITY_731 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 731
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    !C
    A1 = ( XX(2)-XX(1) )**2+( YY(2)-YY(1) )**2+( ZZ(2)-ZZ(1) )**2
    A2 = ( XX(1)-XX(3) )*( XX(2)-XX(1) ) &
      + ( YY(1)-YY(3) )*( YY(2)-YY(1) ) &
      + ( ZZ(1)-ZZ(3) )*( ZZ(2)-ZZ(1) )
    A3 = ( XX(3)-XX(1) )**2+( YY(3)-YY(1) )**2+( ZZ(3)-ZZ(1) )**2
    AA = 0.5*sqrt( A1*A3-A2**2 )
    !C
    do I=1, NN
      TEMP=TT(I)
      call heat_GET_CP     ( TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
      call heat_GET_DENSITY( TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
      SS(I)=SS(I)+AA*THICK*SPE(1)*DEN(1)/NN
    enddo
    !C
    return
    !C
  end subroutine heat_CAPACITY_731

  !> Build up the heat capacity matrix of element 741
  subroutine heat_CAPACITY_741 ( NN,XX,YY,ZZ,TT,IMAT,THICK,SS &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 741
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),SS(NN)
    dimension XG(2),WGT(2),H(4),HR(4),HS(4)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    !C
    XG(1) =-0.5773502691896258D0
    XG(2) =-XG(1)
    WGT(1)= 1.0D0
    WGT(2)= WGT(1)
    AA=0.0
    do LX=1,2
      RI=XG(LX)
      do LY=1,2
        SI=XG(LY)
        RP=1.0+RI
        SP=1.0+SI
        RM=1.0-RI
        SM=1.0-SI
        !C*INTERPOLATION FUNCTION
        H(1)=0.25*RP*SP
        H(2)=0.25*RM*SP
        H(3)=0.25*RM*SM
        H(4)=0.25*RP*SM
        !C*DERIVATIVE OF INTERPOLATION FUNCTION
        !C*  FOR R-COORDINATE
        HR(1)= .25*SP
        HR(2)=-.25*SP
        HR(3)=-.25*SM
        HR(4)= .25*SM
        !C*  FOR S-COORDINATE
        HS(1)= .25*RP
        HS(2)= .25*RM
        HS(3)=-.25*RM
        HS(4)=-.25*RP
        !C*JACOBI MATRIX
        XR=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
        XS=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
        YR=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
        YS=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
        ZR=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
        ZS=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)
        DET=(YR*ZS-ZR*YS)**2+(ZR*XS-XR*ZS)**2+(XR*YS-YR*XS)**2
        DET=sqrt(DET)
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        TEMP=0.0
        do I=1, NN
          TEMP =TEMP+TT(I)*H(I)
        enddo
        call heat_GET_CP     ( TEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( TEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        WG=WGT(LX)*WGT(LY)*DET*THICK
        do I=1, NN
          SS(I)=SS(I)+WG*SPE(1)*DEN(1)*H(I)
        enddo
        !C
      enddo
    enddo
    !C
    return
    !C
  end subroutine heat_CAPACITY_741

  !> Build up the heat capacity matrix of element 232
  !C*--------------------------------------------------------------------*
  subroutine heat_CAPACITY_232 ( NN,XX,YY,ZZ,TT,IMAT,THICK,S0 &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !C*--------------------------------------------------------------------*
    !C*
    !C* HEAT CAPACITY FOR 232
    !C*
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension XG(3),WGT(3),H(6),HL1(6),HL2(6),HL3(6)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    dimension SS(21)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
    do I=1,21
      SS(I)=0.0
    enddo
    TOTD=0.0
    TOTM=0.0
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
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+H(I)*TT(I)
        enddo
        call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        WG=WGT(L1)*WGT(L2)*DET*THICK*(1.0-X2)*0.25
        !C
        ind1=1
        do J=1,NN
          ind2=1
          do J2=1,J
            NUM=(J-1)*J/2+J2
            if(J.EQ.J2) then
              SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              TOTD=TOTD+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
            endif
            TOTM=TOTM+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
            ind2=ind2+1
          enddo
          ind1=ind1+1
        enddo
        !C
      enddo
    enddo
    !C
    do J = 1,NN
      num=J*(J+1)/2
      S0(J)=SS(num)*(2.*TOTM-TOTD)/TOTD
    end do
    !C
    return
    !C
  end subroutine heat_CAPACITY_232

  !> Build up the heat capacity matrix of element 242
  !C*--------------------------------------------------------------------*
  subroutine heat_CAPACITY_242 ( NN,XX,YY,ZZ,TT,IMAT,THICK,S0 &
      ,ntab1,temp1,funcA1,funcB1   &
      ,ntab2,temp2,funcA2,funcB2 )
    !C*--------------------------------------------------------------------*
    !C*
    !C* HEAT CAPACITY FOR 242
    !C*
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension XG(3),WGT(3),H(8),HR(8),HS(8),HT(8)
    dimension SPE(3),DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    dimension SS(36)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
    do I=1,36
      SS(I)=0.0
    enddo
    TOTD=0.0
    TOTM=0.0
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
        !C
        !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
        !C
        CTEMP=0.0
        do I=1,NN
          CTEMP=CTEMP+TT(I)*H(I)
        enddo
        call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
        call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
        WG=WGT(LX)*WGT(LY)*DET*THICK
        !C
        ind1=1
        do J=1,NN
          ind2=1
          do J2=1,J
            NUM=(J-1)*J/2+J2
            if(J.EQ.J2) then
              SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              TOTD=TOTD+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
            endif
            TOTM=TOTM+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
            ind2=ind2+1
          enddo
          ind1=ind1+1
        enddo
        !C
      enddo
    enddo
    !C
    do J = 1,NN
      num=J*(J+1)/2
      S0(J)=SS(num)*(2.*TOTM-TOTD)/TOTD
    end do
    !C
    return
    !C
  end subroutine heat_CAPACITY_242

  !> Build up the heat capacity matrix of element 342
  subroutine heat_CAPACITY_342 ( NN,XX,YY,ZZ,TT,IMAT,S0     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !*
    !
    ! HEAT CAPACITY FOR 342
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension H(10),HL1(10),HL2(10),HL3(10),HL4(10),XG(3),WGT(3)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    dimension SS(55)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C*LOOP OVER ALL INTEGRATION POINTS
    do I=1,55
      SS(I)=0.0
    enddo
    TOTD=0.0
    TOTM=0.0
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
          DET=XJ11*XJ22*XJ33   &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C
          !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
          !C
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+TT(I)*H(I)
          enddo
          call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
          call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
          DET=-DET
          WG=WGT(L1)*WGT(L2)*WGT(L3)*DET*(1.0-X3)*(1.0-X2-X3)*0.125
          ind1=1
          do J=1,NN
            ind2=1
            do J2=1,J
              NUM=(J-1)*J/2+J2
              if(J.EQ.J2) then
                SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
                TOTD=TOTD+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              endif
              TOTM=TOTM+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              ind2=ind2+1
            enddo
            ind1=ind1+1
          enddo
          !C
        enddo
      enddo
    enddo
    !C
    do J = 1,NN
      num=J*(J+1)/2
      S0(J)=SS(num)*(2.*TOTM-TOTD)/TOTD
    end do
    !C
    return
    !C
  end subroutine heat_CAPACITY_342

  !> Build up the heat capacity matrix of element 352
  subroutine heat_CAPACITY_352 ( NN,XX,YY,ZZ,TT,IMAT,S0     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 352
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension H(15),HL1(15),HL2(15),HL3(15),HZ(15),XG(3),WGT(3)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    dimension SS(120)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    ! LOOP FOR INTEGRATION POINTS
    do I=1,120
      SS(I)=0.0
    enddo
    TOTD=0.0
    TOTM=0.0
    do LZ=1,3
      ZI=XG(LZ)
      do L2=1,3
        XL2=XG(L2)
        X2 =(XL2+1.0)*0.5
        do L1=1,3
          XL1=XG(L1)
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
          !DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33   &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !C
          !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
          !C
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+TT(I)*H(I)
          enddo
          call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
          call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
          WG=WGT(L1)*WGT(L2)*WGT(LZ)*DET*(1.0-X2)*0.25
          ind1=1
          do J=1,NN
            ind2=1
            do J2=1,J
              NUM=(J-1)*J/2+J2
              if(J.EQ.J2) then
                SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
                TOTD=TOTD+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              endif
              TOTM=TOTM+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              ind2=ind2+1
            enddo
            ind1=ind1+1
          enddo
          !C
        enddo
      enddo
    enddo
    !C
    do J = 1,NN
      num=J*(J+1)/2
      S0(J)=SS(num)*(2.*TOTM-TOTD)/TOTD
    end do
    !C
    return
    !C
  end subroutine heat_CAPACITY_352

  !> Build up the heat capacity matrix of element 362
  subroutine heat_CAPACITY_362 ( NN,XX,YY,ZZ,TT,IMAT,S0     &
      ,ntab1,temp1,funcA1,funcB1 &
      ,ntab2,temp2,funcA2,funcB2 )
    !
    ! HEAT CAPACITY FOR 362
    !
    use hecmw
    implicit real(kind=kreal)(A-H,O-Z)
    dimension XX(NN),YY(NN),ZZ(NN),TT(NN),S0(NN)
    dimension XG(3),WGT(3),H(20),HR(20),HS(20),HT(20)
    dimension SPE(3), DEN(3)
    dimension temp1(ntab1),funcA1(ntab1+1),funcB1(ntab1+1)
    dimension temp2(ntab2),funcA2(ntab2+1),funcB2(ntab2+1)
    dimension SS(210)
    !C
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    !C
    !C LOOP FOR INTEGRATION POINTS
    do I=1,210
      SS(I)=0.0
    enddo
    TOTD=0.0
    TOTM=0.0
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
          !C INTERPOLATION FUNCTION
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
          !C DERIVATIVE OF INTERPOLATION FUNCTION
          !C FOR R-COORDINATE
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
          !C FOR S-COORDINATE
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
          !C FOR T-COORDINATE
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
          !C JACOBI MATRIX
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
          !C DETERMINANT OF JACOBIAN
          DET=XJ11*XJ22*XJ33   &
            +XJ12*XJ23*XJ31 &
            +XJ13*XJ21*XJ32 &
            -XJ13*XJ22*XJ31 &
            -XJ12*XJ21*XJ33 &
            -XJ11*XJ23*XJ32
          !  WEIGHT VALUE AT GAUSSIAN POINT
          WG=WGT(LX)*WGT(LY)*WGT(LZ)*DET
          !C
          !C    THERMAL CAPACITY AT CURRENT TEMPERATURE
          !C
          CTEMP=0.0
          do I=1,NN
            CTEMP=CTEMP+TT(I)*H(I)
          enddo
          call heat_GET_CP     ( CTEMP,IMAT,SPE,ntab1,temp1,funcA1,funcB1 )
          call heat_GET_DENSITY( CTEMP,IMAT,DEN,ntab2,temp2,funcA2,funcB2 )
          ind1=1
          do J=1,NN
            ind2=1
            do J2=1,J
              NUM=(J-1)*J/2+J2
              if(J.EQ.J2) then
                SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
                TOTD=TOTD+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              endif
              TOTM=TOTM+H(ind1)*H(ind2)*WG*SPE(1)*DEN(1)
              ind2=ind2+1
            enddo
            ind1=ind1+1
          enddo
          !C
        enddo
      enddo
    enddo
    !C
    do J = 1,NN
      num=J*(J+1)/2
      S0(J)=SS(num)*(2.*TOTM-TOTD)/TOTD
    end do
    !C
    return
    !C
  end subroutine heat_CAPACITY_362
  !C***
  !C*** GET SPECIFIC HEAT
  !C***
  subroutine heat_GET_CP ( Tpoi,imat,CP,ntab,temp,funcA,funcB )

    use hecmw

    implicit real(kind=kreal) (A-H,O-Z)
    dimension CP(3),temp(ntab),funcA(ntab+1),funcB(ntab+1)

    itab = 0
    if( Tpoi.LT.temp(1) ) then
      itab= 1
    elseif ( Tpoi.GE.temp(ntab) ) then
      itab= ntab + 1
    else
      do ikk= 1, ntab - 1
        if ( Tpoi.GE.temp(ikk) .AND. Tpoi.LT.temp(ikk+1) ) then
          itab= ikk + 1
          exit
        endif
      enddo
    endif

    CP(1)= funcA(itab)*Tpoi+ funcB(itab)
    CP(2)= CP(1)
    CP(3)= CP(1)
    return

  end subroutine heat_GET_CP
  !C***
  !C*** GET DENSITY
  !C***
  subroutine heat_GET_DENSITY(Tpoi,imat,DENS,ntab,temp,funcA,funcB)

    use hecmw
    implicit real(kind=kreal) (A-H,O-Z)
    dimension DENS(3),temp(ntab),funcA(ntab+1),funcB(ntab+1)

    itab = 0
    if( Tpoi.LT.temp(1) ) then
      itab= 1
    elseif ( Tpoi.GE.temp(ntab) ) then
      itab= ntab + 1
    else
      do ikk= 1, ntab - 1
        if ( Tpoi.GE.temp(ikk) .AND. Tpoi.LT.temp(ikk+1) ) then
          itab= ikk + 1
          exit
        endif
      enddo
    endif

    DENS(1)= funcA(itab)*Tpoi+ funcB(itab)
    DENS(2)= DENS(1)
    DENS(3)= DENS(1)
    return

  end subroutine heat_GET_DENSITY
end module m_heat_LIB_CAPACITY
