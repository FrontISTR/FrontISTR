!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  This module contains subroutines used in 3d eigen analysis for
module m_dynamic_mass
contains

  subroutine mass_C2(etype, nn, ecoord, gausses, mass, lumped, temperature)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(2,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint), parameter :: ndof = 2
    integer(kind=kint) :: i, j, LX, serr
    real(kind=kreal) :: naturalCoord(2)
    real(kind=kreal) :: func(nn)
    real(kind=kreal) :: det, wg, rho
    real(kind=kreal) :: D(2, 2), N(2, nn*ndof), DN(2, nn*ndof)
    real(kind=kreal) :: gderiv(nn, 2)
    real(kind=kreal) :: coordsys(2, 2)
    logical :: is_lumped

    mass(:,:) = 0.0d0
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

      wg = getWeight(etype,LX)*det

      N = 0.0d0
      do i = 1, nn
        N(1,2*i-1) = func(i)
        N(2,2*i  ) = func(i)
      enddo

      DN(1:2, 1:nn*ndof) = matmul(D, N(1:2, 1:nn*ndof))
      forall(i = 1:nn*ndof,  j = 1:nn*ndof)
        mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
      end forall
    enddo

    is_lumped = .true.
    if(is_lumped)then
      lumped = 0.0d0
      do i = 1, nn*ndof
        do j = 1, nn*ndof
          lumped(i) = lumped(i) + mass(j,i)
        enddo
      enddo
      mass = 0.0d0
      do i = 1, nn*ndof
        mass(i,i) = lumped(i)
      enddo
    endif
  end subroutine mass_C2

  subroutine mass_C3(etype, nn, ecoord, gausses, mass, lumped, temperature)
    use mMechGauss
    use m_MatMatrix
    use elementInfo
    type(tGaussStatus), intent(in) :: gausses(:)             !< status of qudrature points
    integer(kind=kint), intent(in) :: etype                  !< element type
    integer(kind=kint), intent(in) :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(out) :: mass(:,:)              !< mass matrix
    real(kind=kreal), intent(out) :: lumped(:)              !< mass matrix
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    type(tMaterial), pointer :: matl !< material information
    integer(kind=kint), parameter :: ndof = 3
    integer(kind=kint) :: i, j, LX, serr
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: func(nn)
    real(kind=kreal) :: det, wg, rho
    real(kind=kreal) :: D(3, 3), N(3, nn*ndof), DN(3, nn*ndof)
    real(kind=kreal) :: gderiv(nn, 3)
    real(kind=kreal) :: coordsys(3, 3)
    logical :: is_lumped

    mass(:,:) = 0.0d0
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
      forall(i = 1:nn*ndof,  j = 1:nn*ndof)
        mass(i,j) = mass(i,j) + dot_product(N(:,i), DN(:,j))*wg
      end forall
    enddo

    is_lumped = .true.
    if(is_lumped)then
      lumped = 0.0d0
      do i = 1, nn*ndof
        do j = 1, nn*ndof
          lumped(i) = lumped(i) + mass(j,i)
        enddo
      enddo
      mass = 0.0d0
      do i = 1, nn*ndof
        mass(i,i) = lumped(i)
      enddo
    endif
  end subroutine mass_C3
end module m_dynamic_mass
