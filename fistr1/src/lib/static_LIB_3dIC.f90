!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Eight-node hexagonal element with incompatible mode
!> \see R.L.Taylor,P.J.Bereford, and E.L.Wilson, "A Nonconforming element
!>  for Stress Analysis", Intl. J. Numer. Methods Engng, 10(6), pp1211-1219
!>  ,1976
module m_static_LIB_3dIC

  use hecmw, only : kint, kreal
  use m_utilities
  use elementInfo

  implicit none

contains

  !>  CALCULATION STIFF Matrix for C3D8IC ELEMENT
  !----------------------------------------------------------------------*
  subroutine STF_C3D8IC &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, u, aux, temperature)
    !----------------------------------------------------------------------*

    use mMechGauss
    use m_MatMatrix
    use m_common_struct
    use m_static_LIB_3d, only: GEOMAT_C3

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)  :: etype                !< element type, not used here
    integer(kind=kint), intent(in)  :: nn                   !< number of elements nodes
    real(kind=kreal), intent(in)    :: ecoord(3, nn)        !< nodal coord of curr element
    type(tGaussStatus), intent(in)  :: gausses(:)           !< Info of qudrature points
    real(kind=kreal), intent(out)   :: stiff(:, :)          !< stiffness matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)         !< variables to define material coordinate system
    real(kind=kreal), intent(in)    :: time                 !< current time
    real(kind=kreal), intent(in)    :: tincr                !< time increment
    real(kind=kreal), intent(in), optional :: u(3, nn)      !< nodal displacemwent
    real(kind=kreal), intent(in), optional :: aux(3, 3)     !< enhanced disp of bending mode
    real(kind=kreal), intent(in)    :: temperature(nn)      !< temperature

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    real(kind=kreal) :: gderiv(nn+3, 3), stress(6)
    real(kind=kreal) :: xj(9, 9), jacobian(3, 3), inverse(3, 3)
    real(kind=kreal) :: tmpstiff((nn+3)*3, (nn+3)*3), tmpk(nn*3, 9)
    real(kind=kreal) :: det, wg, elem(3, nn), mat(6, 6)
    integer(kind=kint) :: i, j, LX
    integer(kind=kint) :: serr
    real(kind=kreal) :: naturalCoord(3), unode(3, nn+3)
    real(kind=kreal) :: gdispderiv(3, 3), coordsys(3, 3)
    real(kind=kreal) :: B1(6, ndof*(nn+3))
    real(kind=kreal) :: Smat(9, 9)
    real(kind=kreal) :: BN(9, ndof*(nn+3)), SBN(9, ndof*(nn+3))
    real(kind=kreal) :: spfunc(nn), temp

    !---------------------------------------------------------------------

    if( present(u) .AND. present(aux) ) then
      unode(:, 1:nn)      = u(:, :)
      unode(:, nn+1:nn+3) = aux(:, :)
    end if

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+unode(:, 1:nn)

    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(etype, nn, naturalcoord, elem, det, jacobian, inverse)
    inverse(:, :)= inverse(:, :)*det
    ! ---- We now calculate stiff matrix include incompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    tmpstiff(:, :) = 0.0D0
    B1(1:6, 1:(nn+3)*ndof) = 0.0D0
    BN(1:9, 1:(nn+3)*ndof) = 0.0D0

    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint(etype, LX, naturalCoord)
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv(1:nn, 1:3) )

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

      ! -- Derivative of shape function of incompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(etype, LX)*det
      B(1:6, 1:(nn+3)*ndof)=0.0D0
      do j = 1, nn+3
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
      if( flag == TOTALLAG ) then
        ! ---dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof,1:ndof) = matmul( unode(1:ndof, 1:nn+3), gderiv(1:nn+3, 1:ndof) )
        do j = 1, nn+3

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
        do j = 1, (nn+3)*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do

      end if

      DB(1:6, 1:(nn+3)*ndof) = matmul( D, B(1:6, 1:(nn+3)*ndof) )
      do  j=1,(nn+3)*ndof
        do i=1,(nn+3)*ndof
          tmpstiff(i, j) = tmpstiff(i, j)+dot_product( B(:, i), DB(:, j) )*wg
        enddo
      enddo

      ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG .OR. flag == UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        do j = 1, nn+3
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
        SBN(1:9, 1:(nn+3)*ndof) = matmul( Smat(1:9, 1:9), BN(1:9, 1:(nn+3)*ndof) )
        do j=1,(nn+3)*ndof 
          do i=1,(nn+3)*ndof
            tmpstiff(i, j) = tmpstiff(i, j)+dot_product( BN(:, i), SBN(:, j) )*wg
          enddo
        enddo
      end if

    end do

    ! -----Condense tmpstiff to stiff
    xj(1:9, 1:9)= tmpstiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    call calInverse(9, xj)
    tmpk = matmul( tmpstiff( 1:nn*ndof, nn*ndof+1:(nn+3)*ndof ), xj )
    stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)-matmul( tmpk, tmpstiff(nn*ndof+1:(nn+3)*ndof, 1:nn*ndof)  )

  end subroutine STF_C3D8IC


  !> Update strain and stress inside element
  !---------------------------------------------------------------------*
  subroutine UPDATE_C3D8IC                                   &
      (etype,nn,ecoord, u, du, ddu, cdsys_ID, coords, &
      qf, gausses, iter, time, tincr, aux, ddaux, TT, T0, TN )
    !---------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities
    use m_static_LIB_3d

    integer(kind=kint), intent(in)    :: etype         !< \param [in] element type
    integer(kind=kint), intent(in)    :: nn            !< \param [in] number of elemental nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: u(3, nn)      !< \param [in] nodal displacements
    real(kind=kreal), intent(in)      :: du(3, nn)     !< \param [in] increment of nodal displacements
    real(kind=kreal), intent(in)      :: ddu(3, nn)    !< \param [in] correction of nodal displacements
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define material coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    integer, intent(in)               :: iter
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(inout)   :: aux(3, 3)     !< \param [in] incompatible dof
    real(kind=kreal), intent(out)     :: ddaux(3, 3)   !< \param [in] increment of incompatible dof
    real(kind=kreal), intent(in)      :: TT(nn)        !< current temperature
    real(kind=kreal), intent(in)      :: T0(nn)        !< reference temperature
    real(kind=kreal), intent(in)      :: TN(nn)        !< reference temperature

    ! LOCAL VARIABLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: D(6,6), B(6,ndof*(nn+3)), B1(6,ndof*(nn+3)), spfunc(nn), ina(1)
    real(kind=kreal)   :: BN(9,ndof*(nn+3)), DB(6, ndof*(nn+3)), stress(6), Smat(9, 9), SBN(9, ndof*(nn+3))
    real(kind=kreal)   :: gderiv(nn+3,3), gderiv0(nn+3,3), gdispderiv(3,3), F(3,3), det, det0, WG, ttc, tt0, ttn
    integer(kind=kint) :: i, j, LX, mtype, serr
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), mat(6,6), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn+3), elem(3,nn), elem1(3,nn), coordsys(3,3)
    real(kind=kreal)   :: dstrain(6)
    real(kind=kreal)   :: alpo(3)
    logical            :: ierr, matlaniso
    real(kind=kreal)   :: stiff_ad(9, nn*3), stiff_aa(9, 9), qf_a(9)
    real(kind=kreal)   :: xj(9, 9)
    real(kind=kreal)   :: tmpforce(9), tmpdisp(9), tmp_d(nn*3), tmp_a(9)
    real(kind=kreal)   :: jacobian(3, 3), inverse(3, 3), inverse1(3, 3), inverse0(3, 3)

    !---------------------------------------------------------------------

    totaldisp(:, 1:nn)      = u(:, :) + du(:, :) - ddu(:, :)
    totaldisp(:, nn+1:nn+3) = aux(:, :)

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+totaldisp(:, 1:nn)

    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(etype, nn, naturalcoord, elem, det, jacobian, inverse)
    inverse(:, :)= inverse(:, :)*det
    ! ---- We now calculate stiff matrix include incompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff_ad(:, :) = 0.0D0
    stiff_aa(:, :) = 0.0D0
    B1(1:6, 1:(nn+3)*ndof) = 0.0D0
    BN(1:9, 1:(nn+3)*ndof) = 0.0D0

    do LX = 1, NumOfQuadPoints(etype)

      call getQuadPoint(etype, LX, naturalCoord)
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv(1:nn, 1:3) )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      call getShapeFunc( etype, naturalcoord, spfunc )
      ttc = dot_product( tt, spfunc )
      call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, ttc )

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      ! -- Derivative of shape function of incompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(etype, LX)*det
      B(1:6, 1:(nn+3)*ndof)=0.0D0
      do j = 1, nn+3
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
      if( flag == TOTALLAG ) then
        ! ---dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof,1:ndof) = matmul( totaldisp(1:ndof, 1:nn+3), gderiv(1:nn+3, 1:ndof) )
        do j = 1, nn+3

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
        do j = 1, (nn+3)*ndof
          B(:, j) = B(:, j)+B1(:, j)
        end do

      end if

      DB(1:6, 1:(nn+3)*ndof) = matmul( D, B(1:6, 1:(nn+3)*ndof) )
      do j=1,nn*ndof 
        do i=1,3*ndof
          stiff_ad(i, j) = stiff_ad(i, j)+dot_product( B(:, nn*ndof+i), DB(:, j) )*wg
        enddo
      enddo

      do j=1,3*ndof 
        do i=1,3*ndof
          stiff_aa(i, j) = stiff_aa(i, j)+dot_product( B(:, nn*ndof+i), DB(:, nn*ndof+j) )*wg
        enddo
      enddo

      ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag == TOTALLAG .OR. flag == UPDATELAG ) then
        stress(1:6) = gausses(LX)%stress
        do j = 1, nn+3
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
        SBN(1:9, 1:(nn+3)*ndof) = matmul( Smat(1:9, 1:9), BN(1:9, 1:(nn+3)*ndof) )

        do j=1,nn*ndof 
          do i=1,3*ndof
            stiff_ad(i, j) = stiff_ad(i, j)+dot_product( BN(:, nn*ndof+i), SBN(:, j) )*wg
          enddo
        enddo
        do j=1,3*ndof 
          do i=1,3*ndof
            stiff_aa(i, j) = stiff_aa(i, j)+dot_product( BN(:, nn*ndof+i), SBN(:, nn*ndof+j) )*wg
          enddo
        enddo
      end if

    end do

    ! -----recover incompatible dof of nodal displacement
    xj(1:3*ndof, 1:3*ndof)= stiff_aa(1:3*ndof, 1:3*ndof)
    call calInverse(3*ndof, xj)
    ! ---  [Kad] * ddu
    do i=1,nn
      tmp_d((i-1)*ndof+1:i*ndof) = ddu(1:ndof, i)
    enddo
    tmpforce(1:3*ndof) = matmul( stiff_ad(1:3*ndof, 1:nn*ndof), tmp_d(1:nn*ndof) )
    ! ---  ddaux = -[Kaa]-1 * ([Kad] * ddu)
    tmp_a(1:3*ndof) = -matmul( xj(1:3*ndof, 1:3*ndof), tmpforce(1:3*ndof) )
    do i=1,3
      ddaux(1:ndof, i) = tmp_a((i-1)*ndof+1:i*ndof)
    enddo


    !---------------------------------------------------------------------

    totaldisp(:,1:nn) = u(:,:)+ du(:,:)
    totaldisp(:,nn+1:nn+3) = aux(:,:) + ddaux(:,:)

    !---------------------------------------------------------------------

    qf(:) = 0.0D0
    qf_a(:) = 0.0D0

    stiff_ad(:, :) = 0.0D0
    stiff_aa(:, :) = 0.0D0

    !---------------------------------------------------------------------

    if( flag == UPDATELAG ) then
      elem(:,:) = (0.5D0*du(:,:)+u(:,:) ) +ecoord(:,:)
      elem1(:,:) = (du(:,:)+u(:,:) ) +ecoord(:,:)
      ! elem = elem1
      totaldisp(:,1:nn) = du(:,:)
      if( iter == 1 ) aux(:,:) = 0.0D0   !--- is this correct???
      totaldisp(:,nn+1:nn+3) = aux(:,:)
    end if

    matlaniso = .FALSE.
    ina = TT(1)
    call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
    if( .not. ierr ) matlaniso = .true.

    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(etype, nn, naturalcoord, elem, det, jacobian, inverse)
    inverse(:, :)= inverse(:, :)*det
    if( flag == UPDATELAG ) then
      call getJacobian(etype, nn, naturalcoord, elem1, det, jacobian, inverse1)
      inverse1(:, :)= inverse1(:, :)*det
      call getJacobian(etype, nn, naturalcoord, ecoord, det, jacobian, inverse0)
      inverse0(:, :)= inverse0(:, :)*det !for deformation gradient F
    endif

    do LX = 1, NumOfQuadPoints(etype)

      mtype = gausses(LX)%pMaterial%mtype

      call getQuadPoint( etype, LX, naturalCoord(:) )
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
      call getShapeFunc(etype, naturalcoord, spfunc)
      ttc = dot_product(TT, spfunc)
      tt0 = dot_product(T0, spfunc)
      ttn = dot_product(TN, spfunc)
      call Cal_Thermal_expansion_C3( tt0, ttc, gausses(LX)%pMaterial, coordsys, matlaniso, EPSTH )

      ! -- Derivative of shape function of incompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      ! Small strain
      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn+3), gderiv(1:nn+3, 1:ndof) )
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)   ! alright?

      F(1:3,1:3) = 0.d0; F(1,1)=1.d0; F(2,2)=1.d0; F(3,3)=1.d0; !deformation gradient
      if( flag == INFINITESIMAL ) then
        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
        F(1:3,1:3) = F(1:3,1:3) + gdispderiv(1:3,1:3)

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

      else if( flag == UPDATELAG ) then
        !  CALL GEOMAT_C3( gausses(LX)%stress, mat )
        !  D(:, :) = D(:, :)+mat(:, :)
        rot = 0.0D0
        rot(1, 2)= 0.5d0*(gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5d0*(gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5d0*(gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+ dstrain(1:6)+EPSTH(:)
        call getGlobalDeriv(etype, nn, naturalcoord, ecoord, det0, gderiv0)
        gderiv0(nn+1, :) = -2.0D0*naturalcoord(1)*inverse0(1, :)/det0
        gderiv0(nn+2, :) = -2.0D0*naturalcoord(2)*inverse0(2, :)/det0
        gderiv0(nn+3, :) = -2.0D0*naturalcoord(3)*inverse0(3, :)/det0
        F(1:3,1:3) = F(1:3,1:3) + matmul( totaldisp(1:ndof, 1:nn+3), gderiv0(1:nn+3, 1:ndof) )

      end if

      ! Update stress
      call Update_Stress3D( flag, gausses(LX), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn )

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:(nn+3)*ndof) = 0.0D0
      do J=1,nn+3
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

        gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn+3), gderiv(1:nn+3, 1:ndof) )
        B1(1:6, 1:(nn+3)*ndof)=0.0D0
        do j = 1,nn+3
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
        do j=1,(nn+3)*ndof
          B(:,j) = B(:,j)+B1(:,j)
        end do

      else if( flag == UPDATELAG ) then

        call getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv(1:nn,1:3))

        ! -- Derivative of shape function of incompatible mode --
        !     [ -2*a   0,   0   ]
        !     [   0,  -2*b, 0   ]
        !     [   0,   0,  -2*c ]
        ! we don't call getShapeDeriv but use above value directly for computation efficiency
        gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse1(1, :)/det
        gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse1(2, :)/det
        gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse1(3, :)/det

        B(1:6, 1:(nn+3)*ndof) = 0.0D0
        do j = 1, nn+3
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

      WG=getWeight( etype, LX )*DET

      DB(1:6, 1:(nn+3)*ndof) = matmul( D, B(1:6, 1:(nn+3)*ndof) )
      
      do j=1,nn*ndof 
        do i=1,3*ndof
          stiff_ad(i, j) = stiff_ad(i, j)+dot_product( B(:, nn*ndof+i), DB(:, j) )*wg
        enddo
      enddo
      do j=1,3*ndof 
        do i=1,3*ndof
          stiff_aa(i, j) = stiff_aa(i, j)+dot_product( B(:, nn*ndof+i), DB(:, nn*ndof+j) )*wg
        enddo
      enddo

      ! calculate the Internal Force
      qf(1:nn*ndof)                                                          &
        = qf(1:nn*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG
      qf_a(1:3*ndof)                                                         &
        = qf_a(1:3*ndof)+matmul( gausses(LX)%stress(1:6), B(1:6,nn*ndof+1:(nn+3)*ndof) )*WG
    end do

    ! condence ( qf - [Kda] * [Kaa]-1 * qf_a )
    xj(1:9, 1:9)= stiff_aa(1:3*ndof, 1:3*ndof)
    call calInverse(9, xj)
    tmpdisp(:) = matmul( xj(:,:), qf_a(:) )
    do i=1,nn*ndof
      qf(i) = qf(i) - dot_product( stiff_ad(:,i), tmpdisp(:) )
    enddo
  end subroutine UPDATE_C3D8IC


  !>  This subroutine calculates thermal loading
  !----------------------------------------------------------------------*
  subroutine TLOAD_C3D8IC &
      (etype, nn, xx, yy, zz, tt, t0, &
      gausses, VECT, cdsys_ID, coords)
    !----------------------------------------------------------------------*

    use m_fstr
    use mMechGauss
    use m_MatMatrix
    use m_common_struct

    !---------------------------------------------------------------------

    integer(kind=kint), parameter     :: ndof=3
    integer(kind=kint), intent(in)    :: etype                   !< element type, not used here
    integer(kind=kint), intent(in)    :: nn                      !< number of element nodes
    real(kind=kreal), intent(in)      :: xx(nn), yy(nn), zz(nn)  !< nodes coordinate of element
    real(kind=kreal), intent(in)      :: tt(nn), t0(nn)          !< current and ref temperature
    type(tGaussStatus), intent(inout) :: gausses(:)              !< info about qudrature points
    real(kind=kreal), intent(out)     :: VECT(nn*ndof)           !< load vector
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)            !< variables to define material coordinate system

    !---------------------------------------------------------------------

    real(kind=kreal) :: ALP, ALP0
    real(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB_a(6, ndof*3)
    real(kind=kreal) :: det, wg, ecoord(3, nn)
    integer(kind=kint) :: i, j, IC, serr
    real(kind=kreal) :: EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3), xj(9, 9)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn+3, 3)
    real(kind=kreal) :: jacobian(3, 3),inverse(3, 3)
    real(kind=kreal) :: stiff_da(nn*3, 3*3), stiff_aa(3*3, 3*3)
    real(kind=kreal) :: VECT_a(3*3)
    real(kind=kreal) :: icdisp(9)
    real(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6), outa(1), ina(1)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    matlaniso = .FALSE.
    if( cdsys_ID > 0 ) then   ! cannot define aniso expansion when no local coord defined
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .true.
    end if

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)
    ! ---- Calculate enhanced displacement at first
    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(etype, nn, naturalcoord, ecoord, det, jacobian, inverse)
    inverse(:, :) = inverse(:, :)*det
    ! ---- We now calculate stiff matrix include incompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff_da(:, :) = 0.0D0
    stiff_aa(:, :) = 0.0D0
    B(1:6, 1:(nn+3)*ndof) = 0.0D0
    VECT(1:nn*ndof) = 0.0D0
    VECT_a(1:3*ndof) = 0.0D0

    do IC = 1, NumOfQuadPoints(etype)

      call getQuadPoint(etype, IC, naturalCoord)
      call getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv(1:nn, 1:3) )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      call getShapeFunc( etype, naturalcoord, H(1:nn) )
      TEMPC = dot_product( H(1:nn), TT(1:nn) )
      TEMP0 = dot_product( H(1:nn), T0(1:nn) )
      call MatlMatrix( gausses(IC), D3, D, 1.d0, 1.0D0, coordsys, TEMPC )

      ! -- Derivative of shape function of incompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(etype, IC)*det
      do j = 1, nn+3
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

      DB_a(1:6, 1:3*ndof) = matmul( D, B(1:6, nn*ndof+1:(nn+3)*ndof) )
      do j=1,3*ndof 
        do i=1,nn*ndof
          stiff_da(i, j) = stiff_da(i, j) + dot_product( B(1:6, i), DB_a(1:6, j) ) * wg
        enddo
      enddo
      do j=1,3*ndof 
        do i=1,3*ndof
          stiff_aa(i, j) = stiff_aa(i, j) + dot_product( B(1:6, nn*ndof+i), DB_a(1:6, j) ) * wg
        enddo
      enddo

      ina(1) = TEMPC
      if( matlaniso ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      end if
      ina(1) = TEMP0
      if( matlaniso  ) then
        call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
        if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
      else
        call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      end if

      !**
      !** THERMAL strain
      !**
      if( matlaniso ) then
        do j = 1, 3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        end do
        EPSTH(4:6) = 0.0D0
        call transformation(coordsys, tm)
        EPSTH(:) = matmul( EPSTH(:), tm  ) ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      else
        THERMAL_EPS = ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
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
      VECT(1:nn*ndof) = VECT(1:nn*ndof)+matmul( SGM(1:6), B(1:6, 1:nn*ndof) )*wg
      VECT_a(1:3*ndof) = VECT_a(1:3*ndof)+matmul( SGM(1:6), B(1:6, nn*ndof+1:(nn+3)*ndof) )*wg

    end do

    ! --- condense vect
    xj(1:9,1:9)= stiff_aa(1:9, 1:9)
    call calInverse(9, xj)
    ! ---  [Kaa]-1 * fa
    icdisp(1:9) = matmul( xj(:, :), VECT_a(1:3*ndof) )
    ! ---  fd - [Kda] * [Kaa]-1 * fa
    VECT(1:nn*ndof) = VECT(1:nn*ndof) - matmul( stiff_da(1:nn*ndof, 1:3*ndof), icdisp(1:3*ndof) )

  end subroutine TLOAD_C3D8IC


end module m_static_LIB_3dIC
