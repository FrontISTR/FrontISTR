!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Uniform strain, reduced integration of the 8 node hexahedron (FORM361=RI).
!>
!> STF_C3D8RI, UPDATE_C3D8RI and TLOAD_C3D8RI are written as STF_C3, UPDATE_C3 and
!> TLOAD_C3 on the single point of the rule. They differ in two places only:
!> getMeanDeriv replaces getGlobalDeriv, and the hourglass control is added after
!> the loop of the points.
module m_static_LIB_C3D8RI

  use hecmw, only : kint, kreal
  use elementInfo
  use m_static_LIB_3d, only: GEOMAT_C3, Update_Stress3D, Cal_Thermal_expansion_C3
  use m_static_LIB_hourglass

  implicit none

  private
  public :: STF_C3D8RI, UPDATE_C3D8RI, TLOAD_C3D8RI

contains

  !> Counterpart of getGlobalDeriv for the uniform strain: the gradient of the shape
  !> functions on the configuration elem and the Jacobian determinant of the point.
  !> The mean gradient g0 over ecoord is mapped with F = I + sum_I (elem_I - ecoord_I) g0_I,
  !> and the weight of the point times det is the volume V0 det F.
  subroutine getMeanDeriv( etype, nn, ecoord, elem, det, gderiv )
    use Quadrature, only: weight3d1
    use m_utilities, only: Determinant33, calInverse

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in)   :: ecoord(3,nn), elem(3,nn)
    real(kind=kreal), intent(out)  :: det, gderiv(nn,3)

    real(kind=kreal)   :: V0, g0(nn,3), F(3,3), detF
    integer(kind=kint) :: i

    call mean_gradient( etype, nn, 3, ecoord, V0, g0 )
    F = matmul( elem-ecoord, g0 )
    do i = 1, 3
      F(i,i) = F(i,i) + 1.d0
    enddo
    detF = Determinant33( F )
    if( detF == 0.d0 ) stop "Error in C3D8RI: the deformation gradient is singular"
    call calInverse( 3, F )
    gderiv = matmul( g0, F )
    det = V0*detF/weight3d1(1)
  end subroutine getMeanDeriv

  !> Tangent stiffness: STF_C3 on the single point, with the hourglass stiffness.
  subroutine STF_C3D8RI                                            &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, u ,temperature, hg_shear)
    !----------------------------------------------------------------------*

    use m_fstr, only: tSolidQuad, fstr_get_solid_quad
    use mMechGauss
    use m_MatMatrix
    use m_common_struct

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)  :: etype                  !< element type
    integer(kind=kint), intent(in)  :: nn                     !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3,3)            !< variables to define material coordinate system
    real(kind=kreal), intent(in)    :: time                   !< current time
    real(kind=kreal), intent(in)    :: tincr                  !< time increment
    real(kind=kreal), intent(in)    :: temperature(nn) !< temperature
    real(kind=kreal), intent(in), optional :: u(:,:)          !< nodal displacemwent
    real(kind=kreal), intent(in)    :: hg_shear               !< initial shear modulus

    !---------------------------------------------------------------------

    type(tSolidQuad) :: quad
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6), B(6, NDOF*nn), DB(6, NDOF*nn)
    real(kind=kreal) :: gderiv(nn, 3), stress(6), mat(6, 6)
    real(kind=kreal) :: det, wg
    integer(kind=kint) :: i, j, LX, serr
    real(kind=kreal) :: temp, naturalCoord(3)
    real(kind=kreal) :: spfunc(nn), gdispderiv(3, 3)
    real(kind=kreal) :: B1(6, NDOF*nn), coordsys(3, 3)
    real(kind=kreal) :: Smat(9, 9), elem(3, nn)
    real(kind=kreal) :: BN(9, NDOF*nn), SBN(9, NDOF*nn)

    !---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)

    call fstr_get_solid_quad( etype, size(gausses, kind=kint), quad )
    do LX = 1, quad%npoints

      naturalCoord(:) = quad%points(:, LX)
      call getMeanDeriv(etype, nn, ecoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          !$omp critical
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
          !$omp end critical
        end if
      end if

      call getShapeFunc(etype, naturalcoord, spfunc)
      temp = dot_product(temperature, spfunc)
      call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      wg = quad%weights(LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
        B(1, 3*j-2)=gderiv(j, 1)
        B(2, 3*j-1)=gderiv(j, 2)
        B(3, 3*j  )=gderiv(j, 3)
        B(4, 3*j-2)=gderiv(j, 2)
        B(4, 3*j-1)=gderiv(j, 1)
        B(5, 3*j-1)=gderiv(j, 3)
        B(5, 3*j  )=gderiv(j, 2)
        B(6, 3*j-2)=gderiv(j, 3)
        B(6, 3*j  )=gderiv(j, 1)
      enddo

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG ) then
        ! ---dudx(i, j) ==> gdispderiv(i, j)
        gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*NDOF)=0.0D0
        do j=1, nn
          B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
          B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
          B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
          B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
          B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
          B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
          B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
          B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
          B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
          B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
          B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
          B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
          B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
          B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
          B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
          B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
          B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
          B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
        end do
        ! ---BL = BL0 + BL1
        do j=1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do
      end if

      DB(1:6, 1:nn*ndof) = matmul( D, B(1:6, 1:nn*ndof) )
      do j=1,nn*ndof
        do i=1,nn*ndof
          stiff(i, j) = stiff(i, j)+dot_product( B(:, i),  DB(:, j) )*WG
        enddo
      enddo

      ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG .OR. flag==UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        BN(1:9, 1:nn*ndof) = 0.0D0
        do j = 1,  nn
          BN(1, 3*j-2) = gderiv(j, 1)
          BN(2, 3*j-1) = gderiv(j, 1)
          BN(3, 3*j  ) = gderiv(j, 1)
          BN(4, 3*j-2) = gderiv(j, 2)
          BN(5, 3*j-1) = gderiv(j, 2)
          BN(6, 3*j  ) = gderiv(j, 2)
          BN(7, 3*j-2) = gderiv(j, 3)
          BN(8, 3*j-1) = gderiv(j, 3)
          BN(9, 3*j  ) = gderiv(j, 3)
        end do
        Smat(:, :) = 0.0D0
        do j = 1, 3
          Smat(j  , j  ) = stress(1)
          Smat(j  , j+3) = stress(4)
          Smat(j  , j+6) = stress(6)
          Smat(j+3, j  ) = stress(4)
          Smat(j+3, j+3) = stress(2)
          Smat(j+3, j+6) = stress(5)
          Smat(j+6, j  ) = stress(6)
          Smat(j+6, j+3) = stress(5)
          Smat(j+6, j+6) = stress(3)
        end do
        SBN(1:9, 1:nn*ndof) = matmul( Smat(1:9, 1:9),  BN(1:9, 1:nn*ndof) )
        do j=1,nn*ndof
          do i=1,nn*ndof
            stiff(i, j) = stiff(i, j)+dot_product( BN(:, i),  SBN(:, j) )*WG
          enddo
        enddo

      end if

    enddo ! gauss roop

    call add_hourglass_stiffness( etype, nn, ndof, ecoord, 1.d0, hg_shear, stiff )

  end subroutine STF_C3D8RI

  !> Update of the state and the internal force: UPDATE_C3 on the single point, with the hourglass force.
  subroutine UPDATE_C3D8RI                                   &
      (etype,nn,ecoord, u, ddu, cdsys_ID, coords, qf, &
      gausses, time, tincr, TT, T0, TN, reftemp, hg_shear)
    !---------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities

    integer(kind=kint), intent(in)    :: etype         !< \param [in] element type
    integer(kind=kint), intent(in)    :: nn            !< \param [in] number of elemental nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: u(3, nn)      !< \param [in] nodal dislplacements
    real(kind=kreal), intent(in)      :: ddu(3, nn)    !< \param [in] nodal displacement
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define material coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(in)      :: TT(nn)   !< current temperature
    real(kind=kreal), intent(in)      :: T0(nn)   !< reference temperature
    real(kind=kreal), intent(in)      :: TN(nn)   !< reference temperature
    real(kind=kreal), intent(in)      :: reftemp  !< reference temperature
    real(kind=kreal), intent(in)      :: hg_shear !< initial shear modulus

    ! LOCAL VARIABLES
    type(tSolidQuad)   :: quad
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: B(6,ndof*nn), B1(6,ndof*nn), spfunc(nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gderiv1(nn,3), gdispderiv(3,3), F(3,3), det, det1, WG, ttc, tt0, ttn
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn), elem(3,nn), elem1(3,nn), coordsys(3,3)
    real(kind=kreal)   :: dstrain(6)
    real(kind=kreal)   :: alpo(3)
    logical            :: ierr, matlaniso

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:,:) = ecoord(:,:)
    totaldisp(:,:) = u(:,:)+ ddu(:,:)
    if( flag == UPDATELAG ) then
      elem(:,:) = (0.5D0*ddu(:,:)+u(:,:) ) +ecoord(:,:)
      elem1(:,:) = (ddu(:,:)+u(:,:) ) +ecoord(:,:)
      ! elem = elem1
      totaldisp(:,:) = ddu(:,:)
    end if

    matlaniso = .FALSE.
    ina = TT(1)
    call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
    if( .not. ierr ) matlaniso = .true.

    call fstr_get_solid_quad( etype, size(gausses, kind=kint), quad )
    do LX = 1, quad%npoints

      naturalCoord(:) = quad%points(:, LX)
      call getMeanDeriv(etype, nn, ecoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          !$omp critical
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
          !$omp end critical
        end if
      end if

      ! ========================================================
      ! UPDATE STRAIN and STRESS
      ! ========================================================

      ! Thermal Strain
      EPSTH = 0.0D0
      call getShapeFunc(etype, naturalcoord, spfunc)
      ttc = dot_product(TT, spfunc)
      tt0 = dot_product(T0, spfunc)
      ttn = dot_product(TN, spfunc)
      call Cal_Thermal_expansion_C3( tt0, ttc, reftemp, gausses(LX)%pMaterial, coordsys, matlaniso, EPSTH )

      ! Update strain
      ! Small strain
      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)   ! alright?

      F(1:3,1:3) = 0.d0; F(1,1)=1.d0; F(2,2)=1.d0; F(3,3)=1.d0; !deformation gradient
      if( flag == INFINITESIMAL ) then
        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(1:6)

      else if( flag == TOTALLAG ) then
        ! Green-Lagrange strain
        dstrain(1) = dstrain(1)+0.5d0*dot_product( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5d0*dot_product( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5d0*dot_product( gdispderiv(:, 3), gdispderiv(:, 3) )
        dstrain(4) = dstrain(4)+( gdispderiv(1, 1)*gdispderiv(1, 2)                                     &
          +gdispderiv(2, 1)*gdispderiv(2, 2)+gdispderiv(3, 1)*gdispderiv(3, 2) )
        dstrain(5) = dstrain(5)+( gdispderiv(1, 2)*gdispderiv(1, 3)                                     &
          +gdispderiv(2, 2)*gdispderiv(2, 3)+gdispderiv(3, 2)*gdispderiv(3, 3) )
        dstrain(6) = dstrain(6)+( gdispderiv(1, 1)*gdispderiv(1, 3)                                     &
          +gdispderiv(2, 1)*gdispderiv(2, 3)+gdispderiv(3, 1)*gdispderiv(3, 3) )

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
        F(1:3,1:3) = F(1:3,1:3) + gdispderiv(1:3,1:3)

      else if( flag == UPDATELAG ) then
        rot = 0.0D0
        rot(1, 2)= 0.5d0*(gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5d0*(gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5d0*(gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)

        call getMeanDeriv(etype, nn, ecoord, ecoord, det1, gderiv1)
        F(1:3,1:3) = F(1:3,1:3) + matmul( u(1:ndof, 1:nn)+ddu(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )

      end if

      ! Update stress
      call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      do J=1,nn
        B(1,3*j-2) = gderiv(j, 1)
        B(2,3*j-1) = gderiv(j, 2)
        B(3,3*j  ) = gderiv(j, 3)
        B(4,3*j-2) = gderiv(j, 2)
        B(4,3*j-1) = gderiv(j, 1)
        B(5,3*j-1) = gderiv(j, 3)
        B(5,3*j  ) = gderiv(j, 2)
        B(6,3*j-2) = gderiv(j, 3)
        B(6,3*j  ) = gderiv(j, 1)
      end do

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == INFINITESIMAL ) then

      else if( flag == TOTALLAG ) then

        gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*ndof)=0.0D0
        do j = 1,nn
          B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
          B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
          B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
          B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
          B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
          B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
          B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
          B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
          B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
          B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
          B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
          B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
          B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
          B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
          B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
          B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
          B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
          B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
        end do
        ! BL = BL0 + BL1
        do j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        end do

      else if( flag == UPDATELAG ) then

        call getMeanDeriv(etype, nn, ecoord, elem1, det, gderiv)
        B(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B(1, 3*J-2) = gderiv(j, 1)
          B(2, 3*J-1) = gderiv(j, 2)
          B(3, 3*J  ) = gderiv(j, 3)
          B(4, 3*J-2) = gderiv(j, 2)
          B(4, 3*J-1) = gderiv(j, 1)
          B(5, 3*J-1) = gderiv(j, 3)
          B(5, 3*J  ) = gderiv(j, 2)
          B(6, 3*J-2) = gderiv(j, 3)
          B(6, 3*J  ) = gderiv(j, 1)
        end do

      end if

      ! calculate the Internal Force
      WG=quad%weights(LX)*DET
      qf(1:nn*ndof)                                                          &
        = qf(1:nn*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG

      ! integrate strain energy
      gausses(LX)%strain_energy = gausses(LX)%strain_energy*WG
    end do

    call add_hourglass_force( etype, nn, ndof, ecoord, 1.d0, hg_shear, u+ddu, qf )

  end subroutine UPDATE_C3D8RI

  !> Equivalent nodal force of the thermal strain: TLOAD_C3 on the single point.
  subroutine TLOAD_C3D8RI                             &
      (etype, nn, XX, YY, ZZ, TT, T0, reftemp, gausses, &
      VECT, cdsys_ID, coords)
    !----------------------------------------------------------------------*

    use m_fstr
    use mMechGauss
    use m_MatMatrix
    use m_utilities

    !---------------------------------------------------------------------

    integer(kind=kint), parameter   :: ndof = 3
    integer(kind=kint), intent(in)  :: etype, nn
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal), intent(in)    :: XX(nn), YY(nn), ZZ(nn)
    real(kind=kreal), intent(in)    :: TT(nn),T0(nn)
    real(kind=kreal), intent(in)    :: reftemp
    real(kind=kreal), intent(out)   :: VECT(nn*NDOF)
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)           !< variables to define material coordinate system

    !---------------------------------------------------------------------

    type(tSolidQuad) :: quad
    real(kind=kreal) :: ALP, ALP0, D(6, 6), B(6, ndof*nn)
    real(kind=kreal) :: det, ecoord(3, nn)
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal) :: EPSTH(6),SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    real(kind=kreal) :: wg, outa(1), ina(1)
    real(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    matlaniso = .FALSE.

    if( cdsys_ID > 0 ) then   ! cannot define aniso expansion when no local coord defined
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .TRUE.
    end if

    VECT(:) = 0.0D0

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)

    ! LOOP FOR INTEGRATION POINTS
    call fstr_get_solid_quad( etype, size(gausses, kind=kint), quad )
    do LX = 1, quad%npoints

      naturalCoord(:) = quad%points(:, LX)
      call getShapeFunc( ETYPE, naturalcoord, H(1:nn) )
      call getMeanDeriv( etype, nn, ecoord, ecoord, det, gderiv )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      ! WEIGHT VALUE AT GAUSSIAN POINT
      wg = quad%weights(LX)*det
      B(1:6,1:nn*NDOF)=0.0D0
      do J=1,nn
        B(1,3*J-2)=gderiv(j,1)
        B(2,3*J-1)=gderiv(j,2)
        B(3,3*J  )=gderiv(j,3)
        B(4,3*J-2)=gderiv(j,2)
        B(4,3*J-1)=gderiv(j,1)
        B(5,3*J-1)=gderiv(j,3)
        B(5,3*J  )=gderiv(j,2)
        B(6,3*J-2)=gderiv(j,3)
        B(6,3*J  )=gderiv(j,1)
      enddo

      TEMPC = dot_product( H(1:nn), TT(1:nn) )
      TEMP0 = dot_product( H(1:nn), T0(1:nn) )

      call MatlMatrix( gausses(LX), D3, D, 1.d0, 0.0D0, coordsys, tempc )

      ina(1) = TEMPC
      if( matlaniso ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXPANSION)
        alp = outa(1)
      end if
      ina(1) = TEMP0
      if( matlaniso  ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXPANSION)
        alp0 = outa(1)
      end if

      !**
      !** THERMAL strain
      !**
      if( matlaniso ) then
        do j=1,3
          EPSTH(j) = ALPO(j)*(TEMPC-reftemp)-alpo0(j)*(TEMP0-reftemp)
        end do
        EPSTH(4:6) = 0.0D0
        call transformation(coordsys, tm)
        EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      else
        THERMAL_EPS=ALP*(TEMPC-reftemp)-alp0*(TEMP0-reftemp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      end if

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = matmul( D(:, :), EPSTH(:) )
      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:nn*NDOF) = VECT(1:nn*NDOF)+matmul( SGM(:),B(:, :) )*wg

    end do

  end subroutine TLOAD_C3D8RI

end module m_static_LIB_C3D8RI
