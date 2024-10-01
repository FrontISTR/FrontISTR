!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains several strategy to free locking problem
!> in Eight-node hexagonal element
module m_static_LIB_C3D8

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

contains


  !>  This subroutine calculate stiff matrix using b-bar method
  !>
  !> \see Hughes, T. J. "Generalization of Selective Integration Procedures
  !>  to Anisotropic and Nonlinear Media", Intl. J. Numer. Methods Engng, 15,
  !>  pp1413-1418,1980
  !----------------------------------------------------------------------*
  subroutine STF_C3D8Bbar &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, u, temperature)
    !----------------------------------------------------------------------*

    use mMechGauss
    use m_MatMatrix
    use m_common_struct
    use m_static_LIB_3d, only: GEOMAT_C3

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)  :: etype                  !< element type
    integer(kind=kint), intent(in)  :: nn                     !< number of elemental nodes
    real(kind=kreal), intent(in)    :: ecoord(3, nn)          !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal), intent(out)   :: stiff(:,:)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)           !< variables to define material coordinate system
    real(kind=kreal), intent(in)    :: time                   !< current time
    real(kind=kreal), intent(in)    :: tincr                  !< time increment
    real(kind=kreal), intent(in), optional :: u(:, :)         !< nodal displacemwent
    real(kind=kreal), intent(in)    :: temperature(nn)        !< temperature

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6),B(6, ndof*nn),DB(6, ndof*nn)
    real(kind=kreal) :: gderiv(nn, 3),stress(6),mat(6, 6)
    real(kind=kreal) :: det, wg, temp, spfunc(nn)
    integer(kind=kint) :: i, j, LX, serr
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: gdispderiv(3, 3)
    real(kind=kreal) :: B1(6, ndof*nn), Bbar(nn, 3)
    real(kind=kreal) :: Smat(9, 9), elem(3, nn)
    real(kind=kreal) :: BN(9, ndof*nn), SBN(9, ndof*nn)
    real(kind=kreal) :: B4, B6, B8, coordsys(3, 3)

    !---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)

    ! dilatation component at centroid
    naturalCoord = 0.0D0
    call getGlobalDeriv(etype, nn, naturalcoord, elem, det, Bbar)

    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      call getShapeFunc( etype, naturalcoord, spfunc )
      temp = dot_product( temperature, spfunc )
      call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
        B4 = ( Bbar(j, 1)-gderiv(j, 1) )/3.0D0
        B6 = ( Bbar(j, 2)-gderiv(j, 2) )/3.0D0
        B8 = ( Bbar(j, 3)-gderiv(j, 3) )/3.0D0
        B(1, 3*j-2) = gderiv(j, 1)+B4
        B(1, 3*j-1) = B6
        B(1, 3*j  ) = B8

        B(2, 3*j-2) = B4
        B(2, 3*j-1) = gderiv(j, 2)+B6
        B(2, 3*j  ) = B8

        B(3, 3*j-2) = B4
        B(3, 3*j-1) = B6
        B(3, 3*j  ) = gderiv(j, 3)+B8

        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      end do

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG ) then
        ! ---dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B1(1, 3*J-2) = gdispderiv(1, 1)*gderiv(J, 1)
          B1(1, 3*J-1) = gdispderiv(2, 1)*gderiv(J, 1)
          B1(1, 3*J  ) = gdispderiv(3, 1)*gderiv(J, 1)
          B1(2, 3*J-2) = gdispderiv(1, 2)*gderiv(J, 2)
          B1(2, 3*J-1) = gdispderiv(2, 2)*gderiv(J, 2)
          B1(2, 3*J  ) = gdispderiv(3, 2)*gderiv(J, 2)
          B1(3, 3*J-2) = gdispderiv(1, 3)*gderiv(J, 3)
          B1(3, 3*J-1) = gdispderiv(2, 3)*gderiv(J, 3)
          B1(3, 3*J  ) = gdispderiv(3, 3)*gderiv(J, 3)
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
        do j = 1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do
      end if

      DB(1:6, 1:nn*ndof) = matmul( D, B(1:6, 1:nn*ndof) )
      do j=1,nn*ndof 
        do i=1,nn*ndof
          stiff(i, j) = stiff(i, j)+dot_product( B(:, i), DB(:, j) )*wg
        end do
      end do

      ! calculate the initial stress matrix
      if( flag == TOTALLAG .or. flag == UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        BN(1:9, 1:nn*ndof) = 0.0D0
        do j = 1, nn
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
        do j= 1, 3
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
        SBN(1:9, 1:nn*ndof) = matmul( Smat(1:9, 1:9), BN(1:9, 1:nn*ndof) )
        do j=1,nn*ndof 
          do i=1,nn*ndof
            stiff(i, j) = stiff(i, j)+dot_product( BN(:, i), SBN(:, j) )*wg
          end do
        end do
      end if

    end do ! gauss roop

  end subroutine STF_C3D8Bbar


  !>  Update Strain stress of this element
  !----------------------------------------------------------------------*
  subroutine Update_C3D8Bbar                              &
      (etype, nn, ecoord, u, du, cdsys_ID, coords, &
      qf ,gausses, iter, time, tincr, TT,T0, TN  )
    !----------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use mHyperElastic
    use m_utilities
    use m_static_LIB_3d

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)    :: etype         !< \param [in] element type
    integer(kind=kint), intent(in)    :: nn            !< \param [in] number of elemental nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: u(3, nn)      !< \param [in] nodal dislplacements
    real(kind=kreal), intent(in)      :: du(3, nn)     !< \param [in] nodal displacement increment
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define material coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    integer, intent(in) :: iter
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(in)      :: TT(nn)        !< current temperature
    real(kind=kreal), intent(in)      :: T0(nn)        !< reference temperature
    real(kind=kreal), intent(in)      :: TN(nn)        !< reference temperature

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: B(6, ndof*nn), B1(6, ndof*nn)
    real(kind=kreal) :: gderiv(nn, 3), gderiv1(nn,3), gdispderiv(3, 3), F(3,3), det, det1, wg
    integer(kind=kint) :: j, LX, mtype, serr
    real(kind=kreal) :: naturalCoord(3), rot(3, 3), spfunc(nn)
    real(kind=kreal) :: totaldisp(3, nn), elem(3, nn), elem1(3, nn), coordsys(3, 3)
    real(kind=kreal) :: dstrain(6)
    real(kind=kreal) :: dvol, vol0, Bbar(nn, 3), derivdum(1:ndof, 1:ndof), BBar2(nn, 3)
    real(kind=kreal) :: B4, B6, B8, ttc, tt0, ttn, alpo(3), ina(1), EPSTH(6)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:, :) = ecoord(:, :)
    totaldisp(:, :) = u(:, :)+du(:, :)
    if( flag == UPDATELAG ) then
      elem(:, :) = ( 0.5D0*du(:, :)+u(:, :) )+ecoord(:, :)
      elem1(:, :) = ( du(:, :)+u(:, :) )+ecoord(:, :)
      !  elem = elem1
      totaldisp(:, :) = du(:, :)
    end if

    matlaniso = .FALSE.
    ina = TT(1)
    call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
    if( .not. ierr ) matlaniso = .TRUE.

    ! dilatation at centroid
    naturalCoord = 0.0D0
    call getGlobalDeriv(etype, nn, naturalcoord, elem, det, Bbar)
    derivdum = matmul( totaldisp(1:ndof, 1:nn), Bbar(1:nn, 1:ndof) )
    vol0 = ( derivdum(1, 1)+derivdum(2, 2)+derivdum(3, 3) )/3.0D0
    if( flag == UPDATELAG ) call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, Bbar2)

    do LX = 1, NumOfQuadPoints(etype)

      mtype = gausses(LX)%pMaterial%mtype

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dvol = vol0-( gdispderiv(1, 1)+gdispderiv(2, 2)+gdispderiv(3, 3) )/3.0D0
      !gdispderiv(1, 1) = gdispderiv(1, 1)+dvol
      !gdispderiv(2, 2) = gdispderiv(2, 2)+dvol
      !gdispderiv(3, 3) = gdispderiv(3, 3)+dvol

      ! ========================================================
      !     UPDATE STRAIN and STRESS
      ! ========================================================

      ! Thermal Strain
      call getShapeFunc(etype, naturalcoord, spfunc)
      ttc = dot_product(TT, spfunc)
      tt0 = dot_product(T0, spfunc)
      ttn = dot_product(TN, spfunc)
      call Cal_Thermal_expansion_C3( tt0, ttc, gausses(LX)%pMaterial, coordsys, matlaniso, EPSTH )

      ! Update strain
      ! Small strain
      dstrain(1) = gdispderiv(1, 1)+dvol
      dstrain(2) = gdispderiv(2, 2)+dvol
      dstrain(3) = gdispderiv(3, 3)+dvol
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)

      F(1:3,1:3) = 0.d0; F(1,1)=1.d0; F(2,2)=1.d0; F(3,3)=1.d0; !deformation gradient
      if( flag == INFINITESIMAL ) then
        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)

      else if( flag == TOTALLAG ) then
        ! Green-Lagrange strain
        dstrain(1) = dstrain(1)+0.5D0*dot_product( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5D0*dot_product( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5D0*dot_product( gdispderiv(:, 3), gdispderiv(:, 3) )
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
        rot(1, 2)= 0.5D0*( gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5D0*( gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5D0*( gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)
        call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det1, gderiv1)
        F(1:3,1:3) = F(1:3,1:3) + matmul( u(1:ndof, 1:nn)+du(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )

      end if

      ! Update stress
      call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
        B4 = ( Bbar(j, 1)-gderiv(j, 1) )/3.0D0
        B6 = ( Bbar(j, 2)-gderiv(j, 2) )/3.0D0
        B8 = ( Bbar(j, 3)-gderiv(j, 3) )/3.0D0
        B(1, 3*j-2) = gderiv(j, 1)+B4
        B(1, 3*j-1) = B6
        B(1, 3*j  ) = B8

        B(2, 3*j-2) = B4
        B(2, 3*j-1) = gderiv(j, 2)+B6
        B(2, 3*j  ) = B8

        B(3, 3*j-2) = B4
        B(3, 3*j-1) = B6
        B(3, 3*j  ) = gderiv(j, 3)+B8

        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      end do

      if( flag == INFINITESIMAL ) then

      else if( flag == TOTALLAG ) then

        B1(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B1(1, 3*J-2) = gdispderiv(1, 1)*gderiv(J, 1)
          B1(1, 3*J-1) = gdispderiv(2, 1)*gderiv(J, 1)
          B1(1, 3*J  ) = gdispderiv(3, 1)*gderiv(J, 1)
          B1(2, 3*J-2) = gdispderiv(1, 2)*gderiv(J, 2)
          B1(2, 3*J-1) = gdispderiv(2, 2)*gderiv(J, 2)
          B1(2, 3*J  ) = gdispderiv(3, 2)*gderiv(J, 2)
          B1(3, 3*J-2) = gdispderiv(1, 3)*gderiv(J, 3)
          B1(3, 3*J-1) = gdispderiv(2, 3)*gderiv(J, 3)
          B1(3, 3*J  ) = gdispderiv(3, 3)*gderiv(J, 3)
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
        do j = 1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do

      else if( flag == UPDATELAG ) then

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
        B(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B4 = ( Bbar2(j, 1)-gderiv(j, 1) )/3.0D0
          B6 = ( Bbar2(j, 2)-gderiv(j, 2) )/3.0D0
          B8 = ( Bbar2(j, 3)-gderiv(j, 3) )/3.0D0
          B(1, 3*j-2) = gderiv(j, 1)+B4
          B(1, 3*j-1) = B6
          B(1, 3*j  ) = B8

          B(2, 3*j-2) = B4
          B(2, 3*j-1) = gderiv(j, 2)+B6
          B(2, 3*j  ) = B8

          B(3, 3*j-2) = B4
          B(3, 3*j-1) = B6
          B(3, 3*j  ) = gderiv(j, 3)+B8

          B(4, 3*j-2) = gderiv(j, 2)
          B(4, 3*j-1) = gderiv(j, 1)
          B(5, 3*j-1) = gderiv(j, 3)
          B(5, 3*j  ) = gderiv(j, 2)
          B(6, 3*j-2) = gderiv(j, 3)
          B(6, 3*j  ) = gderiv(j, 1)
        end do
      end if

      !!  calculate the Internal Force
      wg = getWeight(etype, LX)*det
      qf(1:nn*ndof) = qf(1:nn*ndof)                                           &
        +matmul( gausses(LX)%stress(1:6), B(1:6, 1:nn*ndof) )*wg

    end do

  end subroutine Update_C3D8Bbar

  !> This subroutien calculate thermal loading
  !----------------------------------------------------------------------*
  subroutine TLOAD_C3D8Bbar                    &
      (etype, nn, XX, YY, ZZ, TT, T0,   &
      gausses, VECT, cdsys_ID, coords)
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
    real(kind=kreal), intent(in)    :: TT(nn), T0(nn)
    real(kind=kreal), intent(out)   :: VECT(nn*ndof)
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)           !< variables to define material coordinate system

    !---------------------------------------------------------------------

    real(kind=kreal) :: ALP, alp0, D(6, 6), B(6, ndof*nn)
    real(kind=kreal) :: B4, B6, B8, det, ecoord(3, nn)
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal) :: estrain(6), SGM(6), H(nn)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    real(kind=kreal) :: wg, outa(1), ina(1), Bbar(nn, 3), alpo(3), alpo0(3), coordsys(3, 3)
    real(kind=kreal) :: TEMPC, TEMP0, TEMPC0, TEMP00, THERMAL_EPS, tm(6,6)
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

    naturalCoord = 0.0D0
    call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, Bbar)
    call getShapeFunc( etype, naturalcoord, H(1:nn) )
    TEMPC0 = dot_product( H(1:nn), TT(1:nn) )
    TEMP00 = dot_product( H(1:nn), T0(1:nn) )

    ! LOOP FOR INTEGRATION POINTS
    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc( etype, naturalcoord, H(1:nn) )
      call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      !  WEIGHT VALUE AT GAUSSIAN POINT
      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
        B4 = ( Bbar(j, 1)-gderiv(j, 1) )/3.0D0
        B6 = ( Bbar(j, 2)-gderiv(j, 2) )/3.0D0
        B8 = ( Bbar(j, 3)-gderiv(j, 3) )/3.0D0
        B(1, 3*j-2) = gderiv(j, 1)+B4
        B(1, 3*j-1) = B6
        B(1, 3*j  ) = B8

        B(2, 3*j-2) = B4
        B(2, 3*j-1) = gderiv(j, 2)+B6
        B(2, 3*j  ) = B8

        B(3, 3*j-2) = B4
        B(3, 3*j-1) = B6
        B(3, 3*j  ) = gderiv(j, 3)+B8

        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      end do

      TEMPC = dot_product( H(1:nn),TT(1:nn) )
      TEMP0 = dot_product( H(1:nn),T0(1:nn) )

      call MatlMatrix( gausses(LX), D3, D, 1.d0, 0.0D0, coordsys, TEMPC )

      ina(1) = TEMPC
      if( matlaniso ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      end if
      ina(1) = TEMP0
      if( matlaniso ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      end if

      !**
      !** THERMAL strain
      !**
      if( matlaniso ) then
        do j = 1, 3
          estrain(j) = ALPO(j)*(TEMPC0-ref_temp)-alpo0(j)*(TEMP00-ref_temp)
        end do
        estrain(4:6) = 0.0D0
        call transformation(coordsys, tm)
        estrain(:) = matmul( estrain(:), tm  )      ! to global coord
        estrain(4:6) = estrain(4:6)*2.0D0
      else
        THERMAL_EPS = ALP*(TEMPC0-ref_temp)-alp0*(TEMP00-ref_temp)
        estrain(1:3) = THERMAL_EPS
        estrain(4:6) = 0.0D0
      end if

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = matmul( D(:, :), estrain(:) )

      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:nn*ndof) = VECT(1:nn*ndof)+matmul( SGM(:), B(:, :) )*wg

    end do

  end subroutine TLOAD_C3D8Bbar


end module m_static_LIB_C3D8
