!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_heat_LIB_CAPACITY
  use m_dynamic_mass

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

  subroutine heat_capacity_C1(etype, nn, ecoord, temperature, IMAT, surf, lumped, mass, &
      ntab1, temp1, funcA1, funcB1, ntab2, temp2, funcA2, funcB2)
    use hecmw
    implicit none
    integer(kind=kint), intent(in) :: etype !< element type
    integer(kind=kint), intent(in) :: nn !< number of elemental nodes
    real(kind=kreal), intent(in) :: temperature(nn) !< temperature
    real(kind=kreal), intent(out) :: mass(:,:) !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:) !< mass matrix
    integer(kind=kint), parameter :: ndof = 1
    integer(kind=kint) :: i, j, IMAT
    integer(kind=kint) :: ntab1, ntab2
    real(kind=kreal) :: ecoord(3, nn)
    real(kind=kreal) :: dx, dy, dz, surf, val, temp, length
    real(kind=kreal) :: SPE(3), DEN(3)
    real(kind=kreal) :: temp1(ntab1), funcA1(ntab1+1), funcB1(ntab1+1)
    real(kind=kreal) :: temp2(ntab2), funcA2(ntab2+1), funcB2(ntab2+1)

    if(etype == 611)then
      ecoord(1,2) = ecoord(1,3); ecoord(2,2) = ecoord(2,3); ecoord(3,2) = ecoord(3,3)
    endif

    dx = ecoord(1,2) - ecoord(1,1)
    dy = ecoord(2,2) - ecoord(2,1)
    dz = ecoord(3,2) - ecoord(3,1)
    length = dsqrt(dx*dx + dy*dy + dz*dz)
    val = surf*length

    temp = (temperature(1) + temperature(2))*0.5d0
    call heat_GET_coefficient(temp, IMAT, SPE, ntab1, temp1, funcA1, funcB1)
    call heat_GET_coefficient(temp, IMAT, DEN, ntab2, temp2, funcA2, funcB2)

    lumped(1) = val*SPE(1)*DEN(1)*0.5d0
    lumped(2) = val*SPE(1)*DEN(1)*0.5d0
  end subroutine heat_capacity_C1

  subroutine heat_capacity_C2(etype, nn, ecoord, temperature, IMAT, thick, lumped, mass, &
      ntab1, temp1, funcA1, funcB1, ntab2, temp2, funcA2, funcB2)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(2,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in) :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab1, ntab2
    real(kind=kreal) :: naturalCoord(2)
    real(kind=kreal) :: func(nn), thick, temp
    real(kind=kreal) :: det, wg, rho, diag_mass, total_mass
    real(kind=kreal) :: D(1,1), N(1,nn), DN(1,nn)
    real(kind=kreal) :: gderiv(nn,2)
    real(kind=kreal) :: SPE(3), DEN(3)
    real(kind=kreal) :: temp1(ntab1), funcA1(ntab1+1), funcB1(ntab1+1)
    real(kind=kreal) :: temp2(ntab2), funcA2(ntab2+1), funcB2(ntab2+1)
    logical :: is_lumped

    mass = 0.0d0
    lumped = 0.0d0
    !matl => gausses(1)%pMaterial

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      temp = dot_product(func, temperature)
      call heat_GET_coefficient(temp, IMAT, SPE, ntab1, temp1, funcA1, funcB1)
      call heat_GET_coefficient(temp, IMAT, DEN, ntab2, temp2, funcA2, funcB2)

      D(1,1) = SPE(1)*DEN(1)*thick
      wg = getWeight(etype, LX)*det

      N = 0.0d0
      do i = 1, nn
        N(1,i) = func(i)
      enddo

      DN = matmul(D, N)
      do  j = 1,nn
        do i = 1,nn
          mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
        enddo
      enddo 
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, 1, mass, lumped)
  end subroutine heat_capacity_C2

  subroutine heat_capacity_C3(etype, nn, ecoord, temperature, IMAT, lumped, mass, &
      ntab1, temp1, funcA1, funcB1, ntab2, temp2, funcA2, funcB2)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab1, ntab2
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: func(nn), temp
    real(kind=kreal) :: det, wg, rho, diag_mass, total_mass
    real(kind=kreal) :: D(1, 1), N(1, nn), DN(1, nn)
    real(kind=kreal) :: gderiv(nn, 3)
    real(kind=kreal) :: SPE(3), DEN(3)
    real(kind=kreal) :: temp1(ntab1), funcA1(ntab1+1), funcB1(ntab1+1)
    real(kind=kreal) :: temp2(ntab2), funcA2(ntab2+1), funcB2(ntab2+1)
    logical :: is_lumped

    mass = 0.0d0
    lumped = 0.0d0
    !matl => gausses(1)%pMaterial

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      temp = dot_product(func, temperature)
      call heat_GET_coefficient(temp, IMAT, SPE, ntab1, temp1, funcA1, funcB1)
      call heat_GET_coefficient(temp, IMAT, DEN, ntab2, temp2, funcA2, funcB2)

      D = 0.0d0
      D(1,1) = SPE(1)*DEN(1)
      wg = getWeight(etype, LX)*det

      N = 0.0d0
      do i = 1, nn
        N(1,i) = func(i)
      enddo

      DN = matmul(D, N)
      do j = 1,nn
        do i = 1,nn
          mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
        enddo
      enddo
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, 1, mass, lumped)
  end subroutine heat_capacity_C3

  subroutine heat_capacity_shell_731(etype, nn, ecoord, temperature, IMAT, thick, lumped, mass, &
      ntab1, temp1, funcA1, funcB1, ntab2, temp2, funcA2, funcB2)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    use m_dynamic_mass
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, IMAT, ntab1, ntab2
    real(kind=kreal) :: surf, thick, temp
    real(kind=kreal) :: det, wg, rho, diag_mass, total_mass
    real(kind=kreal) :: SPE(3), DEN(3)
    real(kind=kreal) :: temp1(ntab1), funcA1(ntab1+1), funcB1(ntab1+1)
    real(kind=kreal) :: temp2(ntab2), funcA2(ntab2+1), funcB2(ntab2+1)
    logical :: is_lumped

    mass = 0.0d0
    lumped = 0.0d0
    !matl => gausses(1)%pMaterial

    surf = get_face3(ecoord)

    do i = 1, nn
      temp = temperature(i)
      call heat_GET_coefficient(temp, IMAT, SPE, ntab1, temp1, funcA1, funcB1)
      call heat_GET_coefficient(temp, IMAT, DEN, ntab2, temp2, funcA2, funcB2)

      mass(i,i) = mass(i,i) + surf*thick*SPE(1)*DEN(1)/3.0d0
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, 1, mass, lumped)
  end subroutine heat_capacity_shell_731

  subroutine heat_capacity_shell_741(etype, nn, ecoord, temperature, IMAT, thick, lumped, mass, &
      ntab1, temp1, funcA1, funcB1, ntab2, temp2, funcA2, funcB2)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    !type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint) :: i, j, LX, LY, IMAT, ntab1, ntab2
    real(kind=kreal) :: naturalCoord(2)
    real(kind=kreal) :: func(nn), thick, temp
    real(kind=kreal) :: det, wg, rho, diag_mass, total_mass
    real(kind=kreal) :: D(1,1), N(1, nn), DN(1, nn)
    real(kind=kreal) :: gderiv(nn,2)
    real(kind=kreal) :: SPE(3), DEN(3)
    real(kind=kreal) :: temp1(ntab1), funcA1(ntab1+1), funcB1(ntab1+1)
    real(kind=kreal) :: temp2(ntab2), funcA2(ntab2+1), funcB2(ntab2+1)
    real(kind=kreal) :: XG(2), RI, SI, RP, SP, RM, SM, HR(4), HS(4)
    real(kind=kreal) :: XR, XS, YR, YS, ZR, ZS
    real(kind=kreal) :: H(4), X(4), Y(4), Z(4)
    logical :: is_lumped

    mass = 0.0d0
    lumped = 0.0d0
    !matl => gausses(1)%pMaterial

    X(1) = ecoord(1,1); Y(1) = ecoord(2,1); Z(1) = ecoord(3,1)
    X(2) = ecoord(1,2); Y(2) = ecoord(2,2); Z(2) = ecoord(3,2)
    X(3) = ecoord(1,3); Y(3) = ecoord(2,3); Z(3) = ecoord(3,3)
    X(4) = ecoord(1,4); Y(4) = ecoord(2,4); Z(4) = ecoord(3,4)

    XG(1) = -0.5773502691896258D0
    XG(2) = -XG(1)

    do LX = 1, 2
      RI = XG(LX)
      do LY = 1, 2
        SI = XG(LY)
        RP = 1.0d0 + RI
        SP = 1.0d0 + SI
        RM = 1.0d0 - RI
        SM = 1.0d0 - SI

        !C*INTERPOLATION FUNCTION
        H(1) = 0.25d0*RP*SP
        H(2) = 0.25d0*RM*SP
        H(3) = 0.25d0*RM*SM
        H(4) = 0.25d0*RP*SM

        !C*  FOR R-COORDINATE
        HR(1) =  0.25d0*SP
        HR(2) = -0.25d0*SP
        HR(3) = -0.25d0*SM
        HR(4) =  0.25d0*SM

        !C*  FOR S-COORDINATE
        HS(1) =  0.25d0*RP
        HS(2) =  0.25d0*RM
        HS(3) = -0.25d0*RM
        HS(4) = -0.25d0*RP

        !C*JACOBI MATRIX
        XR = HR(1)*X(1) + HR(2)*X(2) + HR(3)*X(3) + HR(4)*X(4)
        XS = HS(1)*X(1) + HS(2)*X(2) + HS(3)*X(3) + HS(4)*X(4)
        YR = HR(1)*Y(1) + HR(2)*Y(2) + HR(3)*Y(3) + HR(4)*Y(4)
        YS = HS(1)*Y(1) + HS(2)*Y(2) + HS(3)*Y(3) + HS(4)*Y(4)
        ZR = HR(1)*Z(1) + HR(2)*Z(2) + HR(3)*Z(3) + HR(4)*Z(4)
        ZS = HS(1)*Z(1) + HS(2)*Z(2) + HS(3)*Z(3) + HS(4)*Z(4)

        det = (YR*ZS - ZR*YS)**2 + (ZR*XS - XR*ZS)**2 + (XR*YS - YR*XS)**2
        det = dsqrt(det)

        temp = 0.0d0
        do i = 1, 4
          temp = temp + temperature(i)*H(i)
        enddo
        call heat_GET_coefficient(temp, IMAT, SPE, ntab1, temp1, funcA1, funcB1)
        call heat_GET_coefficient(temp, IMAT, DEN, ntab2, temp2, funcA2, funcB2)

        do i = 1, 4
          mass(i,i) = mass(i,i) + SPE(1)*DEN(1)*H(i)*det*thick
        enddo
      enddo
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, 1, mass, lumped)
  end subroutine heat_capacity_shell_741
end module m_heat_LIB_CAPACITY
