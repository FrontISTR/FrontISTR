!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  This module contains subroutines used in 3d eigen analysis for
module m_dynamic_mass

contains

  subroutine mass_C2(etype, nn, ecoord, gausses, sec_opt, thick, mass, lumped, temperature)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(2,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint), parameter :: ndof = 2
    integer(kind=kint) :: i, j, LX, sec_opt
    real(kind=kreal) :: naturalCoord(2)
    real(kind=kreal) :: func(nn), thick
    real(kind=kreal) :: det, wg, rho
    real(kind=kreal) :: D(2,2), N(2, nn*ndof), DN(2, nn*ndof)
    real(kind=kreal) :: gderiv(nn,2)
    logical :: is_lumped

    mass(:,:) = 0.0d0
    lumped = 0.0d0
    matl => gausses(1)%pMaterial

    if(sec_opt == 2) thick = 1.0d0

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      if(present(temperature))then
        !ina(1) = temperature
        !call fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr, ina)
        !if(ierr)then
          rho = matl%variables(M_DENSITY)
        !else
        !  rho = outa(1)
        !endif
      else
        !call fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr)
        !if(ierr)then
          rho = matl%variables(M_DENSITY)
        !else
        !  rho = outa(1)
        !endif
      endif

      D = 0.0d0
      D(1,1) = rho*thick
      D(2,2) = rho*thick

      wg = getWeight(etype, LX)*det

      N = 0.0d0
      do i = 1, nn
        N(1,2*i-1) = func(i)
        N(2,2*i  ) = func(i)
      enddo

      DN(1:2, 1:nn*ndof) = matmul(D, N(1:2, 1:nn*ndof))
      do  j = 1,nn*ndof
        do i = 1,nn*ndof
          mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
        enddo
      enddo
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, ndof, mass, lumped)
  end subroutine mass_C2

  subroutine mass_C3(etype, nn, ecoord, gausses, mass, lumped, temperature)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    implicit none
    type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint), parameter :: ndof = 3
    integer(kind=kint) :: i, j, LX
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: func(nn)
    real(kind=kreal) :: det, wg, rho
    real(kind=kreal) :: D(3, 3), N(3, nn*ndof), DN(3, nn*ndof)
    real(kind=kreal) :: gderiv(nn, 3)
    logical :: is_lumped

    mass(:,:) = 0.0d0
    lumped = 0.0d0
    matl => gausses(1)%pMaterial

    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint(etype, LX, naturalCoord)
      call getShapeFunc(etype, naturalCoord, func)
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      if(present(temperature))then
        !ina(1) = temperature
        !call fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr, ina)
        !if(ierr)then
          rho = matl%variables(M_DENSITY)
        !else
        !  rho = outa(1)
        !endif
      else
        !call fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr)
        !if(ierr)then
          rho = matl%variables(M_DENSITY)
        !else
        !  rho = outa(1)
        !endif
      endif

      D = 0.0d0
      D(1,1) = rho
      D(2,2) = rho
      D(3,3) = rho

      wg = getWeight(etype,LX)*det

      N = 0.0d0
      do i = 1, nn
        N(1,3*i-2) = func(i)
        N(2,3*i-1) = func(i)
        N(3,3*i  ) = func(i)
      enddo

      DN(1:3, 1:nn*ndof) = matmul(D, N(1:3, 1:nn*ndof))
      do  j = 1,nn*ndof
        do i = 1,nn*ndof
          mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
        enddo
      enddo
    enddo

    is_lumped = .true.
    if(is_lumped) call get_lumped_mass(nn, ndof, mass, lumped)
  end subroutine mass_C3

  subroutine get_lumped_mass(nn, ndof, mass, lumped)
    use hecmw
    implicit none
    integer(kind=kint) :: i, j, nn, ndof
    real(kind=kreal) :: lumped(:), mass(:,:)
    real(kind=kreal) :: diag_mass, total_mass

    total_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      do j = 1, nn*ndof, ndof
        total_mass = total_mass + mass(j,i)
      enddo
    enddo

    diag_mass = 0.0d0
    do i = 1, nn*ndof, ndof
      diag_mass = diag_mass + mass(i,i)
    enddo

    diag_mass = 1.0d0/diag_mass
    do i = 1, nn*ndof
      lumped(i) = lumped(i) + mass(i,i)*total_mass*diag_mass
    enddo

    mass = 0.0d0
    do i = 1, nn*ndof
      mass(i,i) = lumped(i)
    enddo
  end subroutine get_lumped_mass

  function get_length(ecoord)
    use hecmw
    implicit none
    real(kind=kreal) :: get_length, ecoord(3,20)

    get_length = dsqrt( &
      (ecoord(1,2) - ecoord(1,1))**2 + &
      (ecoord(2,2) - ecoord(2,1))**2 + &
      (ecoord(3,2) - ecoord(3,1))**2 )
  end function get_length

  function get_face3(ecoord)
    use hecmw
    implicit none
    real(kind=kreal) :: get_face3, ecoord(3,20)
    real(kind=kreal) :: a1, a2, a3
    real(kind=kreal) :: X(3), Y(3), Z(3)

    X(1) = ecoord(1,1); Y(1) = ecoord(2,1); Z(1) = ecoord(3,1)
    X(2) = ecoord(1,2); Y(2) = ecoord(2,2); Z(2) = ecoord(3,2)
    X(3) = ecoord(1,3); Y(3) = ecoord(2,3); Z(3) = ecoord(3,3)

    a1 = (X(2) - X(1))**2 + (Y(2) - Y(1))**2 + (Z(2) - Z(1))**2
    a2 = (X(1) - X(3))*(X(2) - X(1)) &
     & + (Y(1) - Y(3))*(Y(2) - Y(1)) &
     & + (Z(1) - Z(3))*(Z(2) - Z(1))
    a3 = (X(3) - X(1))**2 + (Y(3) - Y(1))**2 + (Z(3) - Z(1))**2

    get_face3 = 0.5d0*dsqrt(a1*a3 - a2*a2)
  end function get_face3

  function get_face4(ecoord)
    use hecmw
    implicit none
    integer(kind=kint) :: LX, LY
    real(kind=kreal) :: get_face4, ecoord(3,20)
    real(kind=kreal) :: XG(2), RI, SI, RP, SP, RM, SM, HR(4), HS(4)
    real(kind=kreal) :: XR, XS, YR, YS, ZR, ZS
    real(kind=kreal) :: X(4), Y(4), Z(4), det

    X(1) = ecoord(1,1); Y(1) = ecoord(2,1); Z(1) = ecoord(3,1)
    X(2) = ecoord(1,2); Y(2) = ecoord(2,2); Z(2) = ecoord(3,2)
    X(3) = ecoord(1,3); Y(3) = ecoord(2,3); Z(3) = ecoord(3,3)
    X(4) = ecoord(1,4); Y(4) = ecoord(2,4); Z(4) = ecoord(3,4)

    XG(1) = -0.5773502691896258D0
    XG(2) = -XG(1)
    get_face4 = 0.0d0

    do LX = 1, 2
      RI = XG(LX)
      do LY = 1, 2
        SI = XG(LY)
        RP = 1.0d0 + RI
        SP = 1.0d0 + SI
        RM = 1.0d0 - RI
        SM = 1.0d0 - SI

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

        get_face4 = get_face4 + det
      enddo
    enddo
  end function get_face4

end module m_dynamic_mass
