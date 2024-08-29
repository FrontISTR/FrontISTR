!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_heat_lib_conductivity
contains

  subroutine heat_GET_coefficient(Tpoi, imat, coef, ntab, temp, funcA, funcB)
    use hecmw
    implicit none
    real(kind=kreal) :: coef(3)
    integer(kind=kint) :: i, in, imat, itab, ntab
    real(kind=kreal) :: Tpoi, temp(ntab), funcA(ntab+1), funcB(ntab+1)

    itab = 0
    if(Tpoi < temp(1)) then
      itab = 1
    elseif(temp(ntab) <= Tpoi)then
      itab = ntab + 1
    else
      do in = 1, ntab - 1
        if(temp(in) <= Tpoi .and. Tpoi < temp(in+1))then
          itab = in + 1
          exit
        endif
      enddo
    endif

    coef(1) = funcA(itab)*Tpoi + funcB(itab)
    coef(2) = coef(1)
    coef(3) = coef(1)
  end subroutine heat_GET_coefficient

  subroutine heat_conductivity_C1(etype, nn, ecoord, temperature, IMAT, surf, stiff, &
      ntab, temp, funcA, funcB)
    use hecmw
    implicit none
    integer(kind=kint), intent(in) :: etype !< element type
    integer(kind=kint), intent(in) :: nn !< number of elemental nodes
    real(kind=kreal), intent(in) :: temperature(nn) !< temperature
    real(kind=kreal), intent(out) :: stiff(:,:) !< stiff matrix
    integer(kind=kint), parameter :: ndof = 1
    integer(kind=kint) :: i, j, IMAT
    integer(kind=kint) :: ntab
    real(kind=kreal) :: ecoord(3, nn)
    real(kind=kreal) :: dx, dy, dz, surf, val, temp_i, length
    real(kind=kreal) :: CC(3)
    real(kind=kreal) :: temp(ntab), funcA(ntab+1), funcB(ntab+1)

    if(etype == 611)then
      ecoord(1,2) = ecoord(1,3); ecoord(2,2) = ecoord(2,3); ecoord(3,2) = ecoord(3,3)
    endif

    dx = ecoord(1,2) - ecoord(1,1)
    dy = ecoord(2,2) - ecoord(2,1)
    dz = ecoord(3,2) - ecoord(3,1)
    length = dsqrt(dx*dx + dy*dy + dz*dz)
    val = surf*length

    temp_i = (temperature(1) + temperature(2))*0.5d0
    call heat_GET_coefficient(temp_i, IMAT, CC, ntab, temp, funcA, funcB)

    stiff(1,1) =  val*CC(1)
    stiff(2,1) = -val*CC(1)
    stiff(1,2) = -val*CC(1)
    stiff(2,2) =  val*CC(1)
  end subroutine heat_conductivity_C1

  subroutine heat_conductivity_C2(etype, nn, ecoord, temperature, IMAT, thick, stiff, &
      ntab, temp, funcA, funcB)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(2,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: stiff(:,:)              !< stiff matrix
    real(kind=kreal), intent(in) :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab
    real(kind=kreal) :: naturalCoord(2)
    real(kind=kreal) :: func(nn), thick, temp_i
    real(kind=kreal) :: det, wg, rho, diag_stiff, total_stiff
    real(kind=kreal) :: D(2,2), N(2,nn), DN(2,nn)
    real(kind=kreal) :: gderiv(nn,2)
    real(kind=kreal) :: CC(3)
    real(kind=kreal) :: temp(ntab), funcA(ntab+1), funcB(ntab+1)

    stiff = 0.0d0
    !matl => gausses(1)%pMaterial

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      temp_i = dot_product(func, temperature)
      call heat_GET_coefficient(temp_i, IMAT, CC, ntab, temp, funcA, funcB)

      D = 0.0d0
      D(1,1) = CC(1)*thick
      D(2,2) = CC(1)*thick
      wg = getWeight(etype, LX)*det

      DN = matmul(D, transpose(gderiv))
      do j = 1,nn
        do i = 1,nn
          stiff(i,j) = stiff(i,j) + dot_product(gderiv(i,:), DN(:,j))*wg
        enddo
      enddo
    enddo
  end subroutine heat_conductivity_C2

  subroutine heat_conductivity_C3(etype, nn, ecoord, temperature, IMAT, stiff, &
      ntab, temp, funcA, funcB)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: stiff(:,:)              !< stiff matrix
    real(kind=kreal), intent(in) :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: func(nn), temp_i
    real(kind=kreal) :: det, wg, rho, diag_stiff, total_stiff
    real(kind=kreal) :: D(3, 3), N(3, nn), DN(3, nn)
    real(kind=kreal) :: gderiv(nn, 3)
    real(kind=kreal) :: SPE(3)
    real(kind=kreal) :: temp(ntab), funcA(ntab+1), funcB(ntab+1)

    stiff = 0.0d0
    !matl => gausses(1)%pMaterial

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      temp_i = dot_product(func, temperature)
      call heat_GET_coefficient(temp_i, IMAT, SPE, ntab, temp, funcA, funcB)

      D = 0.0d0
      D(1,1) = SPE(1)
      D(2,2) = SPE(1)
      D(3,3) = SPE(1)
      wg = getWeight(etype, LX)*det

      DN = matmul(D, transpose(gderiv))

      do j = 1,nn
        do i = 1,nn
          stiff(i,j) = stiff(i,j) + dot_product(gderiv(i,:), DN(:,j))*wg
        enddo
      enddo
    enddo
  end subroutine heat_conductivity_C3

  !> CALCULATION 4 NODE SHELL ELEMENT
  subroutine heat_conductivity_shell_731(etype, nn, ecoord, TT, IMAT, thick, SS, stiff, &
      ntab, temp, funcA, funcB)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    use m_dynamic_mass
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: stiff(:,:)              !< stiff matrix
    real(kind=kreal), intent(inout) :: SS(:)              !< stiff matrix
    real(kind=kreal), intent(inout) :: TT(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab
    integer(kind=kint) :: IG1, IG2, IG3, INOD
    real(kind=kreal) :: surf, thick, temp_i
    real(kind=kreal) :: RI, RM, RP, SI, SM, SP, TI, VALX, VALY, VALZ, VAR
    real(kind=kreal) :: XJ11, XJ12, XJ13, XJ21, XJ22, XJ23, XJ31, XJ32, XJ33, XSUM
    real(kind=kreal) :: det, wg, rho, diag_stiff, total_stiff
    real(kind=kreal) :: CC(3), CTEMP, DUM
    real(kind=kreal) :: temp(ntab), funcA(ntab+1), funcB(ntab+1)
    real(kind=kreal) ::  XG(2), WGT(2), H(4), HR(4), HS(4)
    real(kind=kreal) ::  COD(3, 4)
    real(kind=kreal) ::  G1(3), G2(3), G3(3), E1(3), E2(3), E3(3), REF(3)
    real(kind=kreal) ::  EN(3, 4), THE(3, 3), AMAT(3, 3), BV(3), WK(3)
    real(kind=kreal) ::  DTDX(4), DTDY(4)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    data WGT/1.0D0,1.0D0/

    do I = 1, NN
      COD(1,I) = ecoord(1,i)
      COD(2,I) = ecoord(2,i)
      COD(3,I) = ecoord(3,i)
    enddo

    TT(nn) = TT(nn-1)
    COD(1,nn) = COD(1,nn-1)
    COD(2,nn) = COD(2,nn-1)
    COD(3,nn) = COD(3,nn-1)

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

    !*   LOOP FOR GAUSS INTEGRATION POINT
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

          !* COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
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
          XJ11 = G1(1)
          XJ12 = G1(2)
          XJ13 = G1(3)
          XJ21 = G2(1)
          XJ22 = G2(2)
          XJ23 = G2(3)
          XJ31 = G3(1)
          XJ32 = G3(2)
          XJ33 = G3(3)

          !*DETERMINANT OF JACOBIAN
          DET = XJ11*XJ22*XJ33 &
            + XJ12*XJ23*XJ31 &
            + XJ13*XJ21*XJ32 &
            - XJ13*XJ22*XJ31 &
            - XJ12*XJ21*XJ33 &
            - XJ11*XJ23*XJ32

          !* INVERSION OF JACOBIAN
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

          !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
          XSUM  = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )
          E3(1) = G3(1) / XSUM
          E3(2) = G3(2) / XSUM
          E3(3) = G3(3) / XSUM
          E2(1) = -REF(2)*E3(3) + REF(3)*E3(2)
          E2(2) = -REF(3)*E3(1) + REF(1)*E3(3)
          E2(3) = -REF(1)*E3(2) + REF(2)*E3(1)
          E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
          E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
          E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
          XSUM = dsqrt(E2(1)**2 + E2(2)**2 + E2(3)**2)

          if ( XSUM .GT. 1.E-15 ) then
            E2(1) = E2(1) / XSUM
            E2(2) = E2(2) / XSUM
            E2(3) = E2(3) / XSUM
            E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
            E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
            E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
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

          THE(1,1) = E1(1)
          THE(1,2) = E1(2)
          THE(1,3) = E1(3)
          THE(2,1) = E2(1)
          THE(2,2) = E2(2)
          THE(2,3) = E2(3)
          THE(3,1) = E3(1)
          THE(3,2) = E3(2)
          THE(3,3) = E3(3)

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

          !*CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I = 1, NN
            CTEMP = CTEMP + H(I)*TT(I)
          enddo
          call heat_GET_coefficient( CTEMP,IMAT,CC,ntab,temp,funcA,funcB )

          !* SET INTEGRATION WEIGHT
          VALX = CC(1)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          VALY = CC(2)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
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

    SS( 3) = SS( 3) + SS( 4)
    SS( 7) = SS( 7) + SS( 8)
    SS(11) = SS(11) + SS(16) + 2.0*SS(12)
    SS( 4) = 0.0d0
    SS( 8) = 0.0d0
    SS(12) = 0.0d0
    SS(16) = 0.0d0

    SS( 5) = SS( 2)
    SS( 9) = SS( 3)
    SS(10) = SS( 7)
    SS(13) = SS( 4)
    SS(14) = SS( 8)
    SS(15) = SS(12)
  end subroutine heat_conductivity_shell_731

  subroutine heat_conductivity_shell_741(etype, nn, ecoord, TT, IMAT, thick, SS, stiff, &
      ntab, temp, funcA, funcB)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: stiff(:,:)              !< stiff matrix
    real(kind=kreal), intent(inout) :: SS(:)              !< stiff matrix
    real(kind=kreal), intent(in) :: TT(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, LY, IMAT, ntab
    integer(kind=kint) :: IG1, IG2, IG3, INOD
    real(kind=kreal) :: TI, VALX, VALY, VALZ, VAR
    real(kind=kreal) :: XJ11, XJ12, XJ13, XJ21, XJ22, XJ23, XJ31, XJ32, XJ33
    real(kind=kreal) :: naturalCoord(2), XSUM
    real(kind=kreal) :: func(nn), thick, temp_i
    real(kind=kreal) :: det, wg, rho, diag_stiff, total_stiff
    real(kind=kreal) :: D(1,1), N(1, nn), DN(1, nn)
    real(kind=kreal) :: gderiv(nn,2)
    real(kind=kreal) :: CC(3)
    real(kind=kreal) :: temp(ntab), funcA(ntab+1), funcB(ntab+1)
    real(kind=kreal) :: XG(2), RI, SI, RP, SP, RM, SM, HR(4), HS(4)
    real(kind=kreal) :: XR, XS, YR, YS, ZR, ZS
    real(kind=kreal) :: H(4), X(4), Y(4), Z(4)
    real(kind=kreal) ::  WGT(2), CTEMP, DUM
    real(kind=kreal) ::  COD(3, 4)
    real(kind=kreal) ::  G1(3), G2(3), G3(3), E1(3), E2(3), E3(3), REF(3)
    real(kind=kreal) ::  EN(3, 4), THE(3, 3), AMAT(3, 3), BV(3), WK(3)
    real(kind=kreal) ::  DTDX(4), DTDY(4)
    data XG/-0.5773502691896D0,0.5773502691896D0/
    data WGT/1.0D0,1.0D0/

    do I = 1, NN
      COD(1,I) = ecoord(1,i)
      COD(2,I) = ecoord(2,i)
      COD(3,I) = ecoord(3,i)
    enddo

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

    !*   LOOP FOR GAUSS INTEGRATION POINT
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

          !* COVARIANT BASE VECTOR AT A GAUSS INTEGRARION POINT
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
          XJ11 = G1(1)
          XJ12 = G1(2)
          XJ13 = G1(3)
          XJ21 = G2(1)
          XJ22 = G2(2)
          XJ23 = G2(3)
          XJ31 = G3(1)
          XJ32 = G3(2)
          XJ33 = G3(3)

          !*DETERMINANT OF JACOBIAN
          DET = XJ11*XJ22*XJ33 &
            + XJ12*XJ23*XJ31 &
            + XJ13*XJ21*XJ32 &
            - XJ13*XJ22*XJ31 &
            - XJ12*XJ21*XJ33 &
            - XJ11*XJ23*XJ32

          !* INVERSION OF JACOBIAN
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

          !*  SET BASE VECTOR IN LOCAL CARTESIAN COORDINATE
          XSUM  = dsqrt( G3(1)**2 + G3(2)**2 + G3(3)**2 )
          E3(1) = G3(1) / XSUM
          E3(2) = G3(2) / XSUM
          E3(3) = G3(3) / XSUM
          E2(1) = -REF(2)*E3(3) + REF(3)*E3(2)
          E2(2) = -REF(3)*E3(1) + REF(1)*E3(3)
          E2(3) = -REF(1)*E3(2) + REF(2)*E3(1)
          E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
          E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
          E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
          XSUM = dsqrt( E2(1)**2 + E2(2)**2 + E2(3)**2 )

          if ( XSUM .GT. 1.E-15 ) then
            E2(1) = E2(1) / XSUM
            E2(2) = E2(2) / XSUM
            E2(3) = E2(3) / XSUM
            E1(1) = -E3(2)*E2(3) + E3(3)*E2(2)
            E1(2) = -E3(3)*E2(1) + E3(1)*E2(3)
            E1(3) = -E3(1)*E2(2) + E3(2)*E2(1)
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
          endif
          THE(1,1) = E1(1)
          THE(1,2) = E1(2)
          THE(1,3) = E1(3)
          THE(2,1) = E2(1)
          THE(2,2) = E2(2)
          THE(2,3) = E2(3)
          THE(3,1) = E3(1)
          THE(3,2) = E3(2)
          THE(3,3) = E3(3)
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

          !*CONDUCTIVITY AT CURRENT TEMPERATURE
          CTEMP=0.0
          do I = 1, NN
            CTEMP = CTEMP + H(I)*TT(I)
          enddo
          call heat_GET_coefficient( CTEMP,IMAT,CC,ntab,temp,funcA,funcB )

          !* SET INTEGRATION WEIGHT
          VALX = CC(1)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
          VALY = CC(2)*WGT(IG1)*WGT(IG2)*WGT(IG3)*DET
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
    SS( 5) = SS( 2)
    SS( 9) = SS( 3)
    SS(10) = SS( 7)
    SS(13) = SS( 4)
    SS(14) = SS( 8)
    SS(15) = SS(12)
  end subroutine heat_conductivity_shell_741

  subroutine heat_conductivity_541(NN,ecoord,TEMP,TZERO,THICK,HH,RR1,RR2,SS,stiff )
    use hecmw
    implicit none
    integer(kind=kint), intent(in) :: NN
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: stiff(:,:)              !< stiff matrix
    real(kind=kreal), intent(in)  :: TEMP(NN), TZERO, THICK, HH, RR1, RR2
    real(kind=kreal), intent(inout) :: SS(NN*NN)
    real(kind=kreal) :: SA, SB, SM
    real(kind=kreal) :: T1Z, T2Z, T3Z, T4Z, T5Z, T6Z, T7Z, T8Z
    real(kind=kreal) :: RRR1, RRR2, HHH
    real(kind=kreal) :: HA1, HA2, HA3, HA4
    real(kind=kreal) :: HB1, HB2, HB3, HB4
    real(kind=kreal) :: HH1, HH2, HH3, HH4
    real(kind=kreal) :: XXX(NN), YYY(NN), ZZZ(NN)
    real(kind=kreal) :: XX(4), YY(4), ZZ(4)

    XX(1)=ecoord(1,1)
    XX(2)=ecoord(1,2)
    XX(3)=ecoord(1,3)
    XX(4)=ecoord(1,4)
    YY(1)=ecoord(2,1)
    YY(2)=ecoord(2,2)
    YY(3)=ecoord(2,3)
    YY(4)=ecoord(2,4)
    ZZ(1)=ecoord(3,1)
    ZZ(2)=ecoord(3,2)
    ZZ(3)=ecoord(3,3)
    ZZ(4)=ecoord(3,4)
    call heat_get_area(XX, YY, ZZ, SA)

    XX(1)=ecoord(1,5)
    XX(2)=ecoord(1,6)
    XX(3)=ecoord(1,7)
    XX(4)=ecoord(1,8)
    YY(1)=ecoord(2,5)
    YY(2)=ecoord(2,6)
    YY(3)=ecoord(2,7)
    YY(4)=ecoord(2,8)
    ZZ(1)=ecoord(3,5)
    ZZ(2)=ecoord(3,6)
    ZZ(3)=ecoord(3,7)
    ZZ(4)=ecoord(3,8)
    call heat_get_area(XX, YY, ZZ, SB)

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
  end subroutine heat_conductivity_541

  subroutine heat_get_area ( XX,YY,ZZ,AA )
    use hecmw
    implicit none
    real(kind=kreal), intent(inout) :: AA
    real(kind=kreal), intent(in) :: XX(4),YY(4),ZZ(4)
    integer(kind=kint) :: LX, LY
    real(kind=kreal) :: PI, RI, SI, RP, SP, RM, SM, XR, XS, YR, YS, ZR, ZS, DET
    real(kind=kreal) :: XG(2),H(4),HR(4),HS(4)

    PI=4.0*atan(1.0)

    XG(1) =-0.5773502691896258D0
    XG(2) =-XG(1)
    AA=0.0

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
        AA=AA+DET
      enddo
    enddo
  end subroutine heat_get_area

end module m_heat_lib_conductivity
