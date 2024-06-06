!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains several strategy to free locking problem
!> in Eight-node hexagonal element
module m_static_LIB_Fbar

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

  real(kind=kreal), private, parameter :: I33(3,3) = &
    &  reshape( (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/3,3/) )

contains


  !>  This subroutine calculate stiff matrix using F-bar method
  !----------------------------------------------------------------------*
  subroutine STF_C3D8Fbar &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, u, temperature)
    !----------------------------------------------------------------------*

    use mMechGauss
    use m_MatMatrix
    use m_common_struct
    use m_static_LIB_3d, only: GEOMAT_C3
    use m_utilities

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
    integer(kind=kint) :: i, j, p, q, LX, serr
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: gdispderiv(3, 3)
    real(kind=kreal) :: B1(6, ndof*nn)
    real(kind=kreal) :: Smat(9, 9), elem(3, nn)
    real(kind=kreal) :: BN(9, ndof*nn), SBN(9, ndof*nn)
    real(kind=kreal) :: coordsys(3, 3)

    real(kind=kreal) :: elem0(3,nn), elem1(3, nn), gderiv1(nn,ndof), B2(6, ndof*nn), Z1(3,2)
    real(kind=kreal) :: V0, jacob, jacob_ave, gderiv1_ave(nn,ndof)
    real(kind=kreal) :: gderiv2_ave(ndof*nn,ndof*nn)
    real(kind=kreal) :: Fbar(3,3), Jratio(8), coeff, sff, dstrain(6), ddFS, FS(3,3), GFS(3,2)

    !---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis
    elem(:, :) = ecoord(:, :)
    elem0(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)
    elem1(:, :) = ecoord(:, :)+u(:, :)

    !cal volumetric average of J=detF and dN/dx
    V0 = 0.d0
    jacob_ave = 0.d0
    gderiv1_ave(1:nn,1:ndof) = 0.d0
    gderiv2_ave(1:ndof*nn,1:ndof*nn) = 0.d0
    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem0, det, gderiv)
      wg = getWeight(etype, LX)*det
      if( flag == INFINITESIMAL ) then
        jacob = 1.d0
        gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof) + jacob*wg*gderiv(1:nn, 1:ndof)
      else
        gdispderiv(1:3, 1:3) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        jacob = Determinant33( I33(1:ndof,1:ndof) + gdispderiv(1:ndof, 1:ndof) )
        Jratio(LX) = dsign(dabs(jacob)**(-1.d0/3.d0),jacob)
        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv1)
        gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof) + jacob*wg*gderiv1(1:nn, 1:ndof)
        do p=1,nn
          do q=1,nn
            do i=1,ndof
              do j=1,ndof
                gderiv2_ave(3*p-3+i,3*q-3+j) = gderiv2_ave(3*p-3+i,3*q-3+j)  &
                  & + jacob*wg*(gderiv1(p,i)*gderiv1(q,j)-gderiv1(q,i)*gderiv1(p,j))
              end do
            end do
          end do
        end do
      endif
      V0 = V0 + wg
      jacob_ave = jacob_ave + jacob*wg
    enddo
    if( dabs(V0) > 1.d-12 ) then
      if( dabs(jacob_ave) < 1.d-12 ) stop 'Error in Update_C3D8Fbar: too small average jacob'
      jacob_ave = jacob_ave/V0
      Jratio(1:8) = (dsign(dabs(jacob_ave)**(1.d0/3.d0),jacob_ave))*Jratio(1:8) !Jratio(1:8) = (jacob_ave**(1.d0/3.d0))*Jratio(1:8)
      gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof)/(V0*jacob_ave)
      gderiv2_ave(1:ndof*nn,1:ndof*nn) = gderiv2_ave(1:ndof*nn,1:ndof*nn)/(V0*jacob_ave)
    else
      stop 'Error in Update_C3D8Fbar: too small element volume'
    end if

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
        B(1, 3*j-2) = gderiv(j, 1)
        B(2, 3*j-1) = gderiv(j, 2)
        B(3, 3*j  ) = gderiv(j, 3)
        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      end do

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag == INFINITESIMAL ) then
        B2(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3,1) = (gderiv1_ave(j,1:3)-gderiv(j,1:3))/3.d0
          B2(1,3*j-2:3*j) = Z1(1:3,1)
          B2(2,3*j-2:3*j) = Z1(1:3,1)
          B2(3,3*j-2:3*j) = Z1(1:3,1)
        end do

        ! BL = Jratio*(BL0 + BL1)+BL2
        do j = 1, nn*ndof
          B(1:3, j) = B(1:3,j)+B2(1:3,j)
        end do

      else if( flag == TOTALLAG ) then
        gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        Fbar(1:ndof, 1:ndof) = Jratio(LX)*(I33(1:ndof,1:ndof) + gdispderiv(1:ndof, 1:ndof))
        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv1)

        ! ---dudx(i,j) ==> gdispderiv(i,j)
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

        dstrain(1) = 0.5d0*(dot_product(Fbar(1:3,1),Fbar(1:3,1))-1.d0)
        dstrain(2) = 0.5d0*(dot_product(Fbar(1:3,2),Fbar(1:3,2))-1.d0)
        dstrain(3) = 0.5d0*(dot_product(Fbar(1:3,3),Fbar(1:3,3))-1.d0)
        dstrain(4) = dot_product(Fbar(1:3,1),Fbar(1:3,2))
        dstrain(5) = dot_product(Fbar(1:3,2),Fbar(1:3,3))
        dstrain(6) = dot_product(Fbar(1:3,3),Fbar(1:3,1))

        B2(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3,1) = (gderiv1_ave(j,1:3)-gderiv1(j,1:3))/3.d0
          B2(1,3*j-2:3*j) = (2.d0*dstrain(1)+1.d0)*Z1(1:3,1)
          B2(2,3*j-2:3*j) = (2.d0*dstrain(2)+1.d0)*Z1(1:3,1)
          B2(3,3*j-2:3*j) = (2.d0*dstrain(3)+1.d0)*Z1(1:3,1)
          B2(4,3*j-2:3*j) = 2.d0*dstrain(4)*Z1(1:3,1)
          B2(5,3*j-2:3*j) = 2.d0*dstrain(5)*Z1(1:3,1)
          B2(6,3*j-2:3*j) = 2.d0*dstrain(6)*Z1(1:3,1)
        end do

        ! BL = Jratio*(BL0 + BL1)+BL2
        do j = 1, nn*ndof
          B(:, j) = Jratio(LX)*Jratio(LX)*(B(:,j)+B1(:,j))+B2(:,j)
        end do

      else if( flag == UPDATELAG ) then
        wg = (Jratio(LX)**3.d0)*getWeight(etype, LX)*det

        B2(1:3, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3,1) = (gderiv1_ave(j,1:3)-gderiv(j,1:3))/3.d0
          B2(1, 3*j-2:3*j) = Z1(1:3,1)
          B2(2, 3*j-2:3*j) = Z1(1:3,1)
          B2(3, 3*j-2:3*j) = Z1(1:3,1)
        end do

        do j = 1, nn*ndof
          B(1:3, j) = B(1:3,j)+B2(1:3,j)
        end do

      end if

      DB(1:6, 1:nn*ndof) = matmul( D, B(1:6, 1:nn*ndof) )
      do j=1,nn*ndof 
        do i=1,nn*ndof
          stiff(i, j) = stiff(i, j)+dot_product( B(:, i), DB(:, j) )*wg
        end do
      end do

      ! calculate the initial stress matrix(1): dFbar*dFbar*Stress
      if( flag == TOTALLAG .or. flag == UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        if( flag == TOTALLAG ) then
          coeff = Jratio(LX)
          sff = dot_product(stress(1:6),dstrain(1:6))
        else if( flag == UPDATELAG ) then
          coeff = 1.d0
          sff = stress(1)+stress(2)+stress(3)
          gderiv1 = gderiv
          Fbar(1:3,1:3) = I33(1:3,1:3)
        end if

        BN(1:9, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          BN(1, 3*j-2) = coeff*gderiv(j, 1)
          BN(2, 3*j-1) = coeff*gderiv(j, 1)
          BN(3, 3*j  ) = coeff*gderiv(j, 1)
          BN(4, 3*j-2) = coeff*gderiv(j, 2)
          BN(5, 3*j-1) = coeff*gderiv(j, 2)
          BN(6, 3*j  ) = coeff*gderiv(j, 2)
          BN(7, 3*j-2) = coeff*gderiv(j, 3)
          BN(8, 3*j-1) = coeff*gderiv(j, 3)
          BN(9, 3*j  ) = coeff*gderiv(j, 3)
          Z1(1:3,1) = (gderiv1_ave(j,1:3)-gderiv1(j,1:3))/3.d0
          BN(1, 3*j-2:3*j) = BN(1, 3*j-2:3*j) + Fbar(1,1)*Z1(1:3,1)
          BN(2, 3*j-2:3*j) = BN(2, 3*j-2:3*j) + Fbar(2,1)*Z1(1:3,1)
          BN(3, 3*j-2:3*j) = BN(3, 3*j-2:3*j) + Fbar(3,1)*Z1(1:3,1)
          BN(4, 3*j-2:3*j) = BN(4, 3*j-2:3*j) + Fbar(1,2)*Z1(1:3,1)
          BN(5, 3*j-2:3*j) = BN(5, 3*j-2:3*j) + Fbar(2,2)*Z1(1:3,1)
          BN(6, 3*j-2:3*j) = BN(6, 3*j-2:3*j) + Fbar(3,2)*Z1(1:3,1)
          BN(7, 3*j-2:3*j) = BN(7, 3*j-2:3*j) + Fbar(1,3)*Z1(1:3,1)
          BN(8, 3*j-2:3*j) = BN(8, 3*j-2:3*j) + Fbar(2,3)*Z1(1:3,1)
          BN(9, 3*j-2:3*j) = BN(9, 3*j-2:3*j) + Fbar(3,3)*Z1(1:3,1)
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

        ! calculate the initial stress matrix(2): d(dFbar)*Stress
        FS(1,1) = Fbar(1,1)*stress(1)+Fbar(1,2)*stress(4)+Fbar(1,3)*stress(6)
        FS(1,2) = Fbar(1,1)*stress(4)+Fbar(1,2)*stress(2)+Fbar(1,3)*stress(5)
        FS(1,3) = Fbar(1,1)*stress(6)+Fbar(1,2)*stress(5)+Fbar(1,3)*stress(3)
        FS(2,1) = Fbar(2,1)*stress(1)+Fbar(2,2)*stress(4)+Fbar(2,3)*stress(6)
        FS(2,2) = Fbar(2,1)*stress(4)+Fbar(2,2)*stress(2)+Fbar(2,3)*stress(5)
        FS(2,3) = Fbar(2,1)*stress(6)+Fbar(2,2)*stress(5)+Fbar(2,3)*stress(3)
        FS(3,1) = Fbar(3,1)*stress(1)+Fbar(3,2)*stress(4)+Fbar(3,3)*stress(6)
        FS(3,2) = Fbar(3,1)*stress(4)+Fbar(3,2)*stress(2)+Fbar(3,3)*stress(5)
        FS(3,3) = Fbar(3,1)*stress(6)+Fbar(3,2)*stress(5)+Fbar(3,3)*stress(3)
        do i=1,nn
          Z1(1:3,1) = (gderiv1_ave(i,1:3)-gderiv1(i,1:3))/3.d0
          GFS(1:3,1) = coeff*matmul(FS,gderiv(i,1:3))
          do j=1,nn
            Z1(1:3,2) = (gderiv1_ave(j,1:3)-gderiv1(j,1:3))/3.d0
            GFS(1:3,2) = coeff*matmul(FS,gderiv(j,1:3))
            do p=1,ndof
              do q=1,ndof
                ddFS = Z1(p,1)*Z1(q,2)
                ddFS = ddFS + (gderiv2_ave(3*i-3+p,3*j-3+q)-gderiv1_ave(i,p)*gderiv1_ave(j,q))/3.d0
                ddFS = ddFS + gderiv1(i,q)*gderiv1(j,p)/3.d0
                ddFS = sff*ddFS + Z1(p,1)*GFS(q,2)+Z1(q,2)*GFS(p,1)
                stiff(3*i-3+p, 3*j-3+q) = stiff(3*i-3+p, 3*j-3+q) + ddFS*wg
              end do
            end do
          end do
        end do
      end if


    end do ! gauss roop

  end subroutine STF_C3D8Fbar


  !>  Update Strain stress of this element
  !----------------------------------------------------------------------*
  subroutine Update_C3D8Fbar                              &
      (etype, nn, ecoord, u, du, cdsys_ID, coords, &
      qf, gausses, iter, time, tincr, TT, T0, TN  )
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
    real(kind=kreal) :: gderiv(nn, 3), gdispderiv(3, 3), det, wg
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal) :: naturalCoord(3), rot(3, 3), spfunc(nn)
    real(kind=kreal) :: totaldisp(3, nn), elem(3, nn), elem1(3, nn), coordsys(3, 3)
    real(kind=kreal) :: dstrain(6)
    real(kind=kreal) :: dvol
    real(kind=kreal) :: ttc, tt0, ttn, alpo(3), ina(1), EPSTH(6)
    logical :: ierr, matlaniso

    real(kind=kreal) :: elem0(3,nn), gderiv1(nn,ndof), B2(6, ndof*nn), Z1(3)
    real(kind=kreal) :: V0, jacob, jacob_ave, gderiv1_ave(nn,ndof)
    real(kind=kreal) :: Fbar(3,3), Jratio(8)
    real(kind=kreal) :: jacob_ave05, gderiv05_ave(nn,ndof)

    !---------------------------------------------------------------------

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem0(1:ndof,1:nn) = ecoord(1:ndof,1:nn)
    totaldisp(:, :) = u(:, :)+du(:, :)
    if( flag == INFINITESIMAL ) then
      elem(:, :) = ecoord(:, :)
      elem1(:, :) = ecoord(:, :)
    else if( flag == TOTALLAG ) then
      elem(:, :) = ecoord(:, :)
      elem1(:, :) = ( du(:, :)+u(:, :) )+ecoord(:, :)
    else if( flag == UPDATELAG ) then
      elem(:, :) = ( 0.5D0*du(:, :)+u(:, :) )+ecoord(:, :)
      elem1(:, :) = ( du(:, :)+u(:, :) )+ecoord(:, :)
      totaldisp(:, :) = du(:, :)
    end if

    matlaniso = .FALSE.
    ina = TT(1)
    call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
    if( .not. ierr ) matlaniso = .TRUE.

    !cal volumetric average of J=detF and dN/dx
    V0 = 0.d0
    jacob_ave = 0.d0
    gderiv1_ave(1:nn,1:ndof) = 0.d0
    if( flag == UPDATELAG ) then
      jacob_ave05 = 0.d0
      gderiv05_ave(1:nn,1:ndof) = 0.d0
    endif
    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem0, det, gderiv)
      wg = getWeight(etype, LX)*det
      if( flag == INFINITESIMAL ) then
        jacob = 1.d0
        gderiv1(1:nn, 1:ndof) = gderiv(1:nn, 1:ndof)
      else
        gdispderiv(1:3, 1:3) = matmul( du(1:ndof, 1:nn)+u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        jacob = Determinant33( I33(1:ndof,1:ndof) + gdispderiv(1:ndof, 1:ndof) )
        Jratio(LX) = dsign(dabs(jacob)**(-1.d0/3.d0),jacob) ! Jratio(LX) = jacob**(-1.d0/3.d0)

        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv1)
      endif
      V0 = V0 + wg
      jacob_ave = jacob_ave + jacob*wg
      gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof) + jacob*wg*gderiv1(1:nn, 1:ndof)
      if( flag == UPDATELAG ) then
        call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv1)
        wg = getWeight(etype, LX)*det
        jacob_ave05 = jacob_ave05 + wg
        gderiv05_ave(1:nn,1:ndof) = gderiv05_ave(1:nn,1:ndof) + wg*gderiv1(1:nn, 1:ndof)
      endif
    enddo
    if( dabs(V0) > 1.d-12 ) then
      if( dabs(jacob_ave) < 1.d-12 ) stop 'Error in Update_C3D8Fbar: too small average jacob'
      jacob_ave = jacob_ave/V0
      Jratio(1:8) = (dsign(dabs(jacob_ave)**(1.d0/3.d0),jacob_ave))*Jratio(1:8) !Jratio(1:8) = (jacob_ave**(1.d0/3.d0))*Jratio(1:8)
      gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof)/(V0*jacob_ave)
      if( flag == UPDATELAG ) gderiv05_ave(1:nn,1:ndof) = gderiv05_ave(1:nn,1:ndof)/jacob_ave05
    else
      stop 'Error in Update_C3D8Fbar: too small element volume'
    end if

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

      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )

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
      if( flag == INFINITESIMAL ) then
        dvol = dot_product( totaldisp(1,1:nn), gderiv1_ave(1:nn,1) ) !du1/dx1
        dvol = dvol + dot_product( totaldisp(2,1:nn), gderiv1_ave(1:nn,2) ) !du2/dx2
        dvol = dvol + dot_product( totaldisp(3,1:nn), gderiv1_ave(1:nn,3) ) !du3/dx3
        dvol = (dvol-(gdispderiv(1,1)+gdispderiv(2,2)+gdispderiv(3,3)))/3.d0
        dstrain(1) = gdispderiv(1, 1) + dvol
        dstrain(2) = gdispderiv(2, 2) + dvol
        dstrain(3) = gdispderiv(3, 3) + dvol
        dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
        dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
        dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
        dstrain(:) = dstrain(:)-EPSTH(:)

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)

      else if( flag == TOTALLAG ) then
        Fbar(1:ndof, 1:ndof) = Jratio(LX)*(I33(1:ndof,1:ndof) + gdispderiv(1:ndof, 1:ndof))

        ! Green-Lagrange strain
        dstrain(1) = 0.5d0*(dot_product(Fbar(1:3,1),Fbar(1:3,1))-1.d0)
        dstrain(2) = 0.5d0*(dot_product(Fbar(1:3,2),Fbar(1:3,2))-1.d0)
        dstrain(3) = 0.5d0*(dot_product(Fbar(1:3,3),Fbar(1:3,3))-1.d0)
        dstrain(4) = dot_product(Fbar(1:3,1),Fbar(1:3,2))
        dstrain(5) = dot_product(Fbar(1:3,2),Fbar(1:3,3))
        dstrain(6) = dot_product(Fbar(1:3,3),Fbar(1:3,1))
        dstrain(:) = dstrain(:)-EPSTH(:)

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)

      else if( flag == UPDATELAG ) then
        dvol = dot_product( totaldisp(1,1:nn), gderiv05_ave(1:nn,1) ) !du1/dx1
        dvol = dvol + dot_product( totaldisp(2,1:nn), gderiv05_ave(1:nn,2) ) !du2/dx2
        dvol = dvol + dot_product( totaldisp(3,1:nn), gderiv05_ave(1:nn,3) ) !du3/dx3
        dvol = (dvol-(gdispderiv(1,1)+gdispderiv(2,2)+gdispderiv(3,3)))/3.d0
        dstrain(1) = gdispderiv(1, 1) + dvol
        dstrain(2) = gdispderiv(2, 2) + dvol
        dstrain(3) = gdispderiv(3, 3) + dvol
        dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
        dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
        dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
        dstrain(:) = dstrain(:)-EPSTH(:)

        rot = 0.0D0
        rot(1, 2)= 0.5D0*( gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5D0*( gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5D0*( gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)

        call getGlobalDeriv( etype, nn, naturalcoord, elem0, det, gderiv1)
        gdispderiv(1:ndof, 1:ndof) = matmul( du(1:ndof, 1:nn)+u(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )
        Fbar(1:ndof, 1:ndof) = Jratio(LX)*(I33(1:ndof,1:ndof) + gdispderiv(1:ndof, 1:ndof))

      end if

      ! Update stress
      call Update_Stress3D( flag, gausses(LX), rot, dstrain, Fbar, coordsys, time, tincr, ttc, tt0, ttn )

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      do j = 1, nn
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

      WG=getWeight( etype, LX )*DET
      if( flag == INFINITESIMAL ) then
        gderiv1(1:nn, 1:ndof) = gderiv(1:nn, 1:ndof)

        B2(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3) = (gderiv1_ave(j,1:3)-gderiv1(j,1:3))/3.d0
          B2(1,3*j-2:3*j) = Z1(1:3)
          B2(2,3*j-2:3*j) = Z1(1:3)
          B2(3,3*j-2:3*j) = Z1(1:3)
        end do

        ! BL = BL0 + BL2
        do j = 1, nn*ndof
          B(:, j) = B(:,j)+B2(:,j)
        end do

      else if( flag == TOTALLAG ) then

        B1(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
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

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv1)

        B2(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3) = (gderiv1_ave(j,1:3)-gderiv1(j,1:3))/3.d0
          B2(1,3*j-2:3*j) = (2.d0*dstrain(1)+1.d0)*Z1(1:3)
          B2(2,3*j-2:3*j) = (2.d0*dstrain(2)+1.d0)*Z1(1:3)
          B2(3,3*j-2:3*j) = (2.d0*dstrain(3)+1.d0)*Z1(1:3)
          B2(4,3*j-2:3*j) = 2.d0*dstrain(4)*Z1(1:3)
          B2(5,3*j-2:3*j) = 2.d0*dstrain(5)*Z1(1:3)
          B2(6,3*j-2:3*j) = 2.d0*dstrain(6)*Z1(1:3)
        end do

        ! BL = R^2*(BL0 + BL1)+BL2
        do j = 1, nn*ndof
          B(:, j) = Jratio(LX)*Jratio(LX)*(B(:,j)+B1(:,j))+B2(:,j)
        end do

      else if( flag == UPDATELAG ) then

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
        wg = (Jratio(LX)**3.d0)*getWeight(etype, LX)*det
        B(1:6, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          B(1, 3*j-2) = gderiv(j, 1)
          B(2, 3*j-1) = gderiv(j, 2)
          B(3, 3*j  ) = gderiv(j, 3)
          B(4, 3*j-2) = gderiv(j, 2)
          B(4, 3*j-1) = gderiv(j, 1)
          B(5, 3*j-1) = gderiv(j, 3)
          B(5, 3*j  ) = gderiv(j, 2)
          B(6, 3*j-2) = gderiv(j, 3)
          B(6, 3*j  ) = gderiv(j, 1)
        end do

        B2(1:3, 1:nn*ndof) = 0.0D0
        do j = 1, nn
          Z1(1:3) = (gderiv1_ave(j,1:3)-gderiv(j,1:3))/3.d0
          B2(1, 3*j-2:3*j) = Z1(1:3)
          B2(2, 3*j-2:3*j) = Z1(1:3)
          B2(3, 3*j-2:3*j) = Z1(1:3)
        end do

        do j = 1, nn*ndof
          B(1:3, j) = B(1:3,j)+B2(1:3,j)
        end do

      end if

      qf(1:nn*ndof)                                                          &
        = qf(1:nn*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG

    end do

  end subroutine Update_C3D8Fbar

  !> This subroutien calculate thermal loading
  !----------------------------------------------------------------------*
  subroutine TLOAD_C3D8Fbar                    &
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
    real(kind=kreal) :: Z1(3), det, ecoord(3, nn)
    integer(kind=kint) :: j, LX, serr
    real(kind=kreal) :: estrain(6), SGM(6), H(nn)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    real(kind=kreal) :: wg, outa(1), ina(1), gderiv1_ave(nn, 3), alpo(3), alpo0(3), coordsys(3, 3)
    real(kind=kreal) :: TEMPC, TEMP0, V0, jacob_ave, THERMAL_EPS, tm(6,6)
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

    !cal volumetric average of J=detF and dN/dx
    V0 = 0.d0
    jacob_ave = 0.d0
    gderiv1_ave(1:nn,1:ndof) = 0.d0
    do LX = 1, NumOfQuadPoints(etype)
      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv)
      wg = getWeight(etype, LX)*det
      V0 = V0 + wg
      jacob_ave = jacob_ave + wg
      gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof) + wg*gderiv(1:nn, 1:ndof)
    enddo
    if( dabs(V0) > 1.d-12 ) then
      if( dabs(jacob_ave) < 1.d-12 ) stop 'Error in TLOAD_C3D8Fbar: too small average jacob'
      jacob_ave = jacob_ave/V0
      gderiv1_ave(1:nn,1:ndof) = gderiv1_ave(1:nn,1:ndof)/(V0*jacob_ave)
    else
      stop 'Error in TLOAD_C3D8Fbar: too small element volume'
    end if

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

      do j = 1, nn
        Z1(1:3) = (gderiv1_ave(j,1:3)-gderiv(j,1:3))/3.d0
        B(1,3*j-2:3*j) = B(1,3*j-2:3*j)+Z1(1:3)
        B(2,3*j-2:3*j) = B(2,3*j-2:3*j)+Z1(1:3)
        B(3,3*j-2:3*j) = B(3,3*j-2:3*j)+Z1(1:3)
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
          estrain(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        end do
        estrain(4:6) = 0.0D0
        call transformation(coordsys, tm)
        estrain(:) = matmul( estrain(:), tm  )      ! to global coord
        estrain(4:6) = estrain(4:6)*2.0D0
      else
        THERMAL_EPS = ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
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

  end subroutine TLOAD_C3D8Fbar


end module m_static_LIB_Fbar
