!-------------------------------------------------------------------------------
! Copyright (c) 2021 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_static_LIB_C3D8RI

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

contains


  !----------------------------------------------------------------------*
  subroutine STF_C3D8RI &
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
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3,3)            !< variables to define matreial coordinate system
    real(kind=kreal), intent(in)    :: time                   !< current time
    real(kind=kreal), intent(in)    :: tincr                  !< time increment
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature
    real(kind=kreal), intent(in), optional :: u(:,:)          !< nodal displacemwent

    !---------------------------------------------------------------------

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

    !do LX = 1, NumOfQuadPoints(etype)
    do LX = 1, 1

      !call getQuadPoint( etype, LX, naturalCoord(:) )
      naturalCoord(:) = 0.d0
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      if( present(temperature) ) then
        call getShapeFunc(etype, naturalcoord, spfunc)
        temp = dot_product(temperature, spfunc)
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )
      else
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys )
      end if

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      !wg = getWeight(etype, LX)*det
      wg = 1.d0*det
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
      forall( i=1:nn*ndof,  j=1:nn*ndof )
        stiff(i, j) = stiff(i, j)+dot_product( B(:, i),  DB(:, j) )*WG
      end forall

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
        forall( i=1:nn*ndof,  j=1:nn*ndof )
          stiff(i, j) = stiff(i, j)+dot_product( BN(:, i),  SBN(:, j) )*WG
        end forall

      end if

    enddo ! gauss roop

  end subroutine STF_C3D8RI


  !>  Update Strain stress of this element
  !----------------------------------------------------------------------*
  subroutine Update_C3D8RI                             &
      (etype,nn,ecoord, u, ddu, cdsys_ID, coords, qf, &
      gausses, iter, time, tincr, TT, T0, TN)
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
    real(kind=kreal), intent(in)      :: ddu(3, nn)    !< \param [in] nodal displacement
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define matreial coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    integer, intent(in)               :: iter
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(in), optional :: TT(nn)   !< current temperature
    real(kind=kreal), intent(in), optional :: T0(nn)   !< reference temperature
    real(kind=kreal), intent(in), optional :: TN(nn)   !< reference temperature

    ! LCOAL VARIAVLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: B(6,ndof*nn), B1(6,ndof*nn), spfunc(nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gderiv1(nn,3), gdispderiv(3,3), F(3,3), det, det1, WG, ttc,tt0, ttn
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
    if( present(TT) .and. cdsys_ID > 0 ) then
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .true.
    end if

    !do LX = 1, NumOfQuadPoints(etype)
    do LX = 1, 1

      !call getQuadPoint( etype, LX, naturalCoord(:) )
      naturalCoord(:) = 0.d0
      call getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      ! ========================================================
      ! UPDATE STRAIN and STRESS
      ! ========================================================

      ! Thermal Strain
      EPSTH = 0.0D0
      if( present(tt) .AND. present(t0) ) then
        call getShapeFunc(etype, naturalcoord, spfunc)
        ttc = dot_product(TT, spfunc)
        tt0 = dot_product(T0, spfunc)
        ttn = dot_product(TN, spfunc)
        call Cal_Thermal_expansion_C3( tt0, ttc, gausses(LX)%pMaterial, coordsys, matlaniso, EPSTH )
      end if

      ! Update strain
      ! Small strain
      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)   ! allright?

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

        call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det1, gderiv1)
        F(1:3,1:3) = F(1:3,1:3) + matmul( u(1:ndof, 1:nn)+ddu(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )

      end if

      ! Update stress
      if( present(tt) .AND. present(t0) ) then
        call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )
      else
        call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr )
      end if

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

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
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
      !WG=getWeight( etype, LX )*DET
      WG = 1.d0*DET
      qf(1:nn*ndof)                                                          &
        = qf(1:nn*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG

    end do

  end subroutine Update_C3D8RI

end module m_static_LIB_C3D8RI
