!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to generate heat transfer
!! boundary
module m_heat_LIB_FILM
contains
  !***********************************************************************
  !  FILM_231(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_232(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_241(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_242(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_341(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_342(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_351(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_352(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_361(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_362(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
  !  FILM_731(NN,XX,YY,ZZ,LTYPE,HH,SINK,TERM1,TERM2)
  !  FILM_741(NN,XX,YY,ZZ,LTYPE,HH,SINK,TERM1,TERM2)
  !----------------------------------------------------------------------*
  subroutine heat_FILM_231(NN,XX,YY,ZZ,THICK,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 231 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) NN, LTYPE, MM
    integer(kind=kint) NOD(MM)
    real(kind=kreal)   THICK, HH, SINK
    real(kind=kreal)   XX(NN),YY(NN),ZZ(NN),TERM1(MM*MM),TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) I, IC, IP, JP, LX
    real(kind=kreal)   RI, GX, GY, XSUM
    real(kind=kreal)   XG(2),WGT(2),H(2),HR(2)
    data WGT/1.0,1.0/
    data XG/-0.5773502691896, 0.5773502691896/
    !
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
    !
    IC = 0
    do IP = 1, 2
      TERM2(IP) = 0.0
      do JP = 1, 2
        IC = IC + 1
        TERM1(IC) = 0.0
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
          XSUM = dsqrt( GX*GX+GY*GY )*THICK
          TERM1(IC) = TERM1(IC) - XSUM*WGT(LX)*H(IP)*H(JP)*HH
          if( IP.EQ.JP )  then
            TERM2(IP) = TERM2(IP) - XSUM*WGT(LX)*H(JP)*HH*SINK
          endif
        enddo
      enddo
    enddo
    !
    return

  end subroutine heat_FILM_231
  !----------------------------------------------------------------------*
  subroutine heat_FILM_232(NN,XX,YY,ZZ,THICK,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 232 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: THICK, HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, LX
    real(kind=kreal)   :: RI, GX, GY, XSUM
    real(kind=kreal)   :: XG(3), WGT(3), H(3), HR(3)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555

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

    IC = 0
    do IP = 1, 3
      TERM2(IP) = 0.0
      do JP = 1, 3
        IC = IC + 1
        TERM1(IC) = 0.0
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
          XSUM = dsqrt( GX*GX+GY*GY )*THICK
          TERM1(IC) = TERM1(IC) - XSUM*WGT(LX)*H(IP)*H(JP)*HH
          if( IP.EQ.JP ) then
            TERM2(IP) = TERM2(IP) - XSUM*WGT(LX)*H(JP)*HH*SINK
          endif
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_232
  !----------------------------------------------------------------------*
  subroutine heat_FILM_241(NN,XX,YY,ZZ,THICK,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 241 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: THICK, HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, LX
    real(kind=kreal)   :: RI, GX, GY, XSUM
    real(kind=kreal)   :: XG(2), WGT(2), H(2), HR(2)
    data WGT/1.0, 1.0/
    data XG/-0.5773502691896, 0.5773502691896/

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

    IC = 0
    do IP = 1, 2
      TERM2(IP) = 0.0
      do JP = 1, 2
        IC = IC + 1
        TERM1(IC) = 0.0
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
          XSUM = dsqrt( GX*GX+GY*GY )*THICK
          TERM1(IC) = TERM1(IC) - XSUM*WGT(LX)*H(IP)*H(JP)*HH
          if( IP.EQ.JP ) then
            TERM2(IP) = TERM2(IP) - XSUM*WGT(LX)*H(JP)*HH*SINK
          endif
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_241
  !----------------------------------------------------------------------*
  subroutine heat_FILM_242(NN,XX,YY,ZZ,THICK,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 242 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: THICK, HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, LX
    real(kind=kreal)   :: RI, GX, GY, XSUM
    real(kind=kreal)   :: XG(3), WGT(3), H(3), HR(3)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555

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

    IC = 0
    do IP = 1, 3
      TERM2(IP) = 0.0
      do JP = 1, 3
        IC = IC + 1
        TERM1(IC) = 0.0
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
          XSUM = dsqrt( GX*GX+GY*GY )*THICK
          TERM1(IC) = TERM1(IC) - XSUM*WGT(LX)*H(IP)*H(JP)*HH
          if( IP.EQ.JP ) then
            TERM2(IP) = TERM2(IP) - XSUM*WGT(LX)*H(JP)*HH*SINK
          endif
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_242
  !----------------------------------------------------------------------*
  subroutine heat_FILM_341(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 341 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: IC, IP, JP
    real(kind=kreal)   :: AX, AY, AZ, BX, BY, BZ, AA

    !
    if     ( LTYPE.EQ.1 ) then
      NOD(1) = 1
      NOD(2) = 2
      NOD(3) = 3
    else if( LTYPE.EQ.2 ) then
      NOD(1) = 4
      NOD(2) = 2
      NOD(3) = 1
    else if( LTYPE.EQ.3 ) then
      NOD(1) = 4
      NOD(2) = 3
      NOD(3) = 2
    else if( LTYPE.EQ.4 ) then
      NOD(1) = 4
      NOD(2) = 1
      NOD(3) = 3
    endif
    !
    AX = XX( NOD(2) ) - XX( NOD(1) )
    AY = YY( NOD(2) ) - YY( NOD(1) )
    AZ = ZZ( NOD(2) ) - ZZ( NOD(1) )
    BX = XX( NOD(3) ) - XX( NOD(1) )
    BY = YY( NOD(3) ) - YY( NOD(1) )
    BZ = ZZ( NOD(3) ) - ZZ( NOD(1) )
    !
    AA = dsqrt( (AY*BZ-AZ*BY)**2+(AZ*BX-AX*BZ)**2+(AX*BY-AY*BX)**2 )/6.0
    !
    IC = 0
    do IP = 1, 3
      TERM2(IP) = - AA*HH*SINK
      do JP = 1, 3
        IC = IC + 1
        TERM1(IC) = - AA*HH/3.0
      enddo
    enddo
    !
    return

  end subroutine heat_FILM_341
  !----------------------------------------------------------------------*
  subroutine heat_FILM_342(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 342 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, L1, L2
    real(kind=kreal)   :: X1, X2, X3, XL1, XL2
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: XSUM, DET, WG
    real(kind=kreal)   :: XG(3), WGT(3), H(6), HL1(6), HL2(6), HL3(6)
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
    !
    IC = 0
    do IP = 1, 6
      TERM2(IP) = 0.0
      do JP = 1, 6
        IC = IC + 1
        TERM1(IC) = 0.0
        ! INTEGRATION OVER SURFACE
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
            !
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
            DET=XJ11*XJ22*XJ33 &
              +XJ12*XJ23*XJ31 &
              +XJ13*XJ21*XJ32 &
              -XJ13*XJ22*XJ31 &
              -XJ12*XJ21*XJ33 &
              -XJ11*XJ23*XJ32
            !
            WG=WGT(L1)*WGT(L2)*DET*(1.0-X2)*0.25*HH
            TERM1(IC) = TERM1(IC)-WG*H(IP)*H(JP)
            if( IP.EQ.JP ) TERM2(IP) = TERM2(IP)-WG*H(JP)*SINK
            !
          enddo
        enddo
      enddo
    enddo
    !
    return

  end subroutine heat_FILM_342
  !----------------------------------------------------------------------*
  subroutine heat_FILM_351(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 351 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    !   F5   LTYPE=5  : FILM IN NORMAL-DIRECTION FOR FACE-5
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: ISUF, I, IC, IP, JP, IG1, IG2
    real(kind=kreal)   :: RI, SI, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: AX, AY, AZ, BX, BY, BZ, AA
    real(kind=kreal)   :: H(6), HR(6), HS(6), HT(6)
    real(kind=kreal)   :: XG(2), WGT(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/

    ISUF = 0
    !
    if( LTYPE.EQ.1 ) then
      NOD(1) = 1
      NOD(2) = 2
      NOD(3) = 3
      ISUF = 1
    else if( LTYPE.EQ.2 ) then
      NOD(1) = 6
      NOD(2) = 5
      NOD(3) = 4
      ISUF = 1
    else if( LTYPE.EQ.3 ) then
      NOD(1) = 4
      NOD(2) = 5
      NOD(3) = 2
      NOD(4) = 1
      ISUF = 2
    else if( LTYPE.EQ.4 ) then
      NOD(1) = 5
      NOD(2) = 6
      NOD(3) = 3
      NOD(4) = 2
      ISUF = 2
    else if( LTYPE.EQ.5 ) then
      NOD(1) = 6
      NOD(2) = 4
      NOD(3) = 1
      NOD(4) = 3
      ISUF = 2
    endif
    !
    if( ISUF.EQ.2 ) then

      IC = 0
      do IP = 1, 4
        TERM2(IP) = 0.0
        do JP = 1, 4
          IC = IC + 1
          TERM1(IC) = 0.0
          !
          do  IG2=1,2
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
              do I = 1, 4
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
              !
              TERM1(IC)=TERM1(IC)-XSUM*WGT(IG1)*WGT(IG2)*H(IP)*H(JP)*HH
              if( IP.EQ.JP ) then
                TERM2(IP)=TERM2(IP)-XSUM*WGT(IG1)*WGT(IG2)*H(JP)*HH*SINK
              endif
            enddo
          enddo
        enddo
      enddo
    elseif( ISUF.EQ.1 ) then
      AX = XX( NOD(2) ) - XX( NOD(1) )
      AY = YY( NOD(2) ) - YY( NOD(1) )
      AZ = ZZ( NOD(2) ) - ZZ( NOD(1) )
      BX = XX( NOD(3) ) - XX( NOD(1) )
      BY = YY( NOD(3) ) - YY( NOD(1) )
      BZ = ZZ( NOD(3) ) - ZZ( NOD(1) )
      AA = dsqrt( (AY*BZ-AZ*BY)**2+(AZ*BX-AX*BZ)**2+(AX*BY-AY*BX)**2 )/6.0
      !
      IC = 0
      do IP = 1, 3
        TERM2(IP) = - AA*HH*SINK
        do JP = 1, 3
          IC = IC + 1
          TERM1(IC) = - AA*HH/3.0
        enddo
      enddo

    endif

    return

  end subroutine heat_FILM_351
  !----------------------------------------------------------------------*
  subroutine heat_FILM_352(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 352 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    !   F5   LTYPE=5  : FILM IN NORMAL-DIRECTION FOR FACE-5
    !   F6   LTYPE=6  : FILM IN NORMAL-DIRECTION FOR FACE-6
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: ISUF, I, IC, IP, JP, IG1, IG2, L1, L2
    real(kind=kreal)   :: RI, SI, RP, SP, RM, SM
    real(kind=kreal)   :: XSUM, DET, WG
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: X1, X2, X3, XL1, XL2
    real(kind=kreal)   :: XJ11, XJ21, XJ31, XJ12, XJ22, XJ32, XJ13, XJ23, XJ33
    real(kind=kreal)   :: H(8), HR(8), HS(8), HT(8), HL1(6), HL2(6), HL3(6)
    real(kind=kreal)   :: XG(3), WGT(3)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    XG(1) = -0.7745966692
    XG(2) =  0.0
    XG(3) =  0.7745966692
    WGT(1) = 0.5555555555
    WGT(2) = 0.8888888888
    WGT(3) = 0.5555555555
    ISUF = 0
    !
    if( LTYPE.EQ.1 ) then
      NOD(1) = 1
      NOD(2) = 2
      NOD(3) = 3
      NOD(4) = 7
      NOD(5) = 8
      NOD(6) = 9
      ISUF = 1
    else if( LTYPE.EQ.2 ) then
      NOD(1) = 6
      NOD(2) = 5
      NOD(3) = 4
      NOD(4) = 11
      NOD(5) = 10
      NOD(6) = 12
      ISUF = 1
    else if( LTYPE.EQ.3 ) then
      NOD(1) = 4
      NOD(2) = 5
      NOD(3) = 2
      NOD(4) = 1
      NOD(5) = 10
      NOD(6) = 14
      NOD(7) = 7
      NOD(8) = 13
      ISUF = 2
    else if( LTYPE.EQ.4 ) then
      NOD(1) = 5
      NOD(2) = 6
      NOD(3) = 3
      NOD(4) = 2
      NOD(5) = 11
      NOD(6) = 15
      NOD(7) = 8
      NOD(8) = 14
      ISUF = 2
    else if( LTYPE.EQ.5 ) then
      NOD(1) = 6
      NOD(2) = 4
      NOD(3) = 1
      NOD(4) = 3
      NOD(5) = 12
      NOD(6) = 13
      NOD(7) = 9
      NOD(8) = 15
      ISUF = 2
    endif
    !
    if( ISUF.EQ.2 ) then
      IC = 0
      do IP = 1,8
        TERM2(IP) = 0.0
        do JP = 1,8
          IC = IC + 1
          TERM1(IC) = 0.0
          do  IG2=1,3
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
              do I = 1,8
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
              WG=XSUM*WGT(IG1)*WGT(IG2)*HH
              TERM1(IC)=TERM1(IC)-WG*H(IP)*H(JP)
              if( IP.EQ.JP ) TERM2(IP)=TERM2(IP)-WG*H(JP)*SINK
            enddo
          enddo
        enddo
      enddo
      !
    elseif( ISUF.EQ.1 ) then
      IC = 0
      do IP = 1, 6
        TERM2(IP) = 0.0
        do JP = 1, 6
          IC = IC + 1
          TERM1(IC) = 0.0
          ! INTEGRATION OVER SURFACE
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
              !
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
              DET=XJ11*XJ22*XJ33   &
                +XJ12*XJ23*XJ31   &
                +XJ13*XJ21*XJ32   &
                -XJ13*XJ22*XJ31   &
                -XJ12*XJ21*XJ33   &
                -XJ11*XJ23*XJ32
              WG=WGT(L1)*WGT(L2)*DET*(1.0-X2)*0.25*HH
              TERM1(IC) = TERM1(IC)-WG*H(IP)*H(JP)
              if( IP.EQ.JP ) TERM2(IP) = TERM2(IP)-WG*H(JP)*SINK
              !
            enddo
          enddo
        enddo
      enddo
    endif
    !
    return

  end subroutine heat_FILM_352
  !----------------------------------------------------------------------*
  subroutine heat_FILM_361(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 361 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    !   F5   LTYPE=5  : FILM IN NORMAL-DIRECTION FOR FACE-5
    !   F6   LTYPE=6  : FILM IN NORMAL-DIRECTION FOR FACE-6
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, IG1, IG2
    real(kind=kreal)   :: RI, SI, XSUM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: H(4), HR(4), HS(4), HT(4)
    real(kind=kreal)   :: XG(2), WGT(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !
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
    !
    IC = 0
    do IP = 1, 4
      TERM2(IP) = 0.0
      do JP = 1, 4
        IC = IC + 1
        TERM1(IC) = 0.0
        do  IG2=1,2
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
            do I = 1, 4
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
            !
            TERM1(IC)=TERM1(IC)-XSUM*WGT(IG1)*WGT(IG2)*H(IP)*H(JP)*HH
            if( IP.EQ.JP ) then
              TERM2(IP)=TERM2(IP)-XSUM*WGT(IG1)*WGT(IG2)*H(JP)*HH*SINK
            endif
            !
          enddo
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_361
  !----------------------------------------------------------------------*
  subroutine heat_FILM_362(NN,XX,YY,ZZ,LTYPE,HH,SINK,MM,TERM1,TERM2,NOD)
    !----------------------------------------------------------------------*
    !**
    !**  SET 362 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    !   F2   LTYPE=2  : FILM IN NORMAL-DIRECTION FOR FACE-2
    !   F3   LTYPE=3  : FILM IN NORMAL-DIRECTION FOR FACE-3
    !   F4   LTYPE=4  : FILM IN NORMAL-DIRECTION FOR FACE-4
    !   F5   LTYPE=5  : FILM IN NORMAL-DIRECTION FOR FACE-5
    !   F6   LTYPE=6  : FILM IN NORMAL-DIRECTION FOR FACE-6
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN, LTYPE, MM
    integer(kind=kint) :: NOD(MM)
    real(kind=kreal)   :: HH, SINK
    real(kind=kreal)   :: XX(NN), YY(NN), ZZ(NN), TERM1(MM*MM), TERM2(MM)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, IG1, IG2
    real(kind=kreal)   :: RI, SI, RP, SP, RM, SM
    real(kind=kreal)   :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal)   :: XSUM
    real(kind=kreal)   :: H(8), HR(8), HS(8), HT(8)
    real(kind=kreal)   :: XG(3), WGT(3)
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

    IC = 0
    do IP = 1,8
      TERM2(IP) = 0.0
      do JP = 1,8
        IC = IC + 1
        TERM1(IC) = 0.0
        !
        do  IG2=1,3
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

            do I = 1,8
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
            !
            TERM1(IC)=TERM1(IC)-XSUM*WGT(IG1)*WGT(IG2)*H(IP)*H(JP)*HH
            if( IP.EQ.JP ) then
              TERM2(IP) = TERM2(IP)- XSUM*WGT(IG1)*WGT(IG2)*H(JP)*HH*SINK
            endif
            !
          enddo
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_362
  !----------------------------------------------------------------------*
  subroutine heat_FILM_731( NN,XX,YY,ZZ,HH,SINK,TERM1,TERM2 )
    !----------------------------------------------------------------------*
    !**
    !**  SET 731 FILM
    !**
    !   F1   LTYPE=1  : SURFACE FILM
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN
    real(kind=kreal) :: XX(NN), YY(NN), ZZ(NN), HH, SINK, TERM1(NN*NN), TERM2(NN)
    ! LOCAL VARIABLES
    real(kind=kreal) :: AX, AY, AZ, BX, BY, BZ, AA
    integer(kind=kint) :: IC, IP, JP
    !
    AX = XX(2) - XX(1)
    AY = YY(2) - YY(1)
    AZ = ZZ(2) - ZZ(1)
    BX = XX(3) - XX(1)
    BY = YY(3) - YY(1)
    BZ = ZZ(3) - ZZ(1)
    AA = dsqrt( (AY*BZ-AZ*BY)**2+(AZ*BX-AX*BZ)**2+(AX*BY-AY*BX)**2 )/6.0
    !
    IC = 0
    do IP = 1, NN
      TERM2(IP) = - AA*HH*SINK
      do JP = 1, NN
        IC = IC + 1
        TERM1(IC) = - AA*HH/3.0
      enddo
    enddo

    return

  end subroutine heat_FILM_731
  !----------------------------------------------------------------------*
  subroutine heat_FILM_741( NN,XX,YY,ZZ,HH,SINK,TERM1,TERM2 )
    !----------------------------------------------------------------------*
    !**
    !**  SET 741 FILM
    !**
    !   F1   LTYPE=1  : FILM IN NORMAL-DIRECTION FOR FACE-1
    use hecmw
    implicit none
    ! I/F VARIABLES
    integer(kind=kint) :: NN
    real(kind=kreal) :: XX(NN), YY(NN), ZZ(NN), HH, SINK, TERM1(NN*NN), TERM2(NN)
    ! LOCAL VARIABLES
    integer(kind=kint) :: I, IC, IP, JP, IG1, IG2
    real(kind=kreal) :: RI, SI, XSUM
    real(kind=kreal) :: G1X, G1Y, G1Z, G2X, G2Y, G2Z, G3X, G3Y, G3Z
    real(kind=kreal) :: H(4), HR(4), HS(4), HT(4)
    real(kind=kreal) :: XG(2), WGT(2)
    !*************************
    !  GAUSS INTEGRATION POINT
    !*************************
    data XG/-0.5773502691896,0.5773502691896/
    data WGT/1.0,1.0/
    !
    IC = 0
    do IP = 1, NN
      TERM2(IP) = 0.0
      do JP = 1,  NN
        IC = IC + 1
        TERM1(IC) = 0.0
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
            do I = 1, NN
              G1X=G1X+HR(I)*XX(I)
              G1Y=G1Y+HR(I)*YY(I)
              G1Z=G1Z+HR(I)*ZZ(I)
              G2X=G2X+HS(I)*XX(I)
              G2Y=G2Y+HS(I)*YY(I)
              G2Z=G2Z+HS(I)*ZZ(I)
            enddo
            G3X=G1Y*G2Z-G1Z*G2Y
            G3Y=G1Z*G2X-G1X*G2Z
            G3Z=G1X*G2Y-G1Y*G2X
            XSUM=dsqrt(G3X**2+G3Y**2+G3Z**2)
            !
            TERM1(IC)=TERM1(IC)-XSUM*WGT(IG1)*WGT(IG2)*H(IP)*H(JP)*HH
            if( IP.EQ.JP ) then
              TERM2(IP)=TERM2(IP)-XSUM*WGT(IG1)*WGT(IG2)*H(JP)*HH*SINK
            endif
            !
          enddo
        enddo
      enddo
    enddo

    return

  end subroutine heat_FILM_741
end module m_heat_LIB_FILM
