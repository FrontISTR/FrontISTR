!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Eight-node hexagonal element with imcompatible mode
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
      time, tincr, nddisp, ehdisp, temperature)
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
    real(kind=kreal), intent(inout) :: coords(3, 3)         !< variables to define matreial coordinate system
    real(kind=kreal), intent(in)    :: time                 !< current time
    real(kind=kreal), intent(in)    :: tincr                !< time increment
    real(kind=kreal), intent(in), optional :: nddisp(3, nn) !< nodal displacemwent
    real(kind=kreal), intent(in), optional :: ehdisp(3, 3)  !< enhanced disp of bending mode
    real(kind=kreal), intent(in), optional :: temperature(nn) !< temperature

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    real(kind=kreal) :: gderiv(nn+3, 3), stress(6)
    real(kind=kreal) :: xj(9, 9), jacobian(3, 3), inverse(3, 3)
    real(kind=kreal) :: tmpstiff((nn+3)*3, (nn+3)*3), tmpk(nn*3, 9)
    real(kind=kreal) :: det, wg, elem(3, nn), mat(6, 6)
    integer(kind=kint) :: i, j, LX, fetype
    integer(kind=kint) :: serr
    real(kind=kreal) :: naturalCoord(3), unode(3, nn+3)
    real(kind=kreal) :: gdispderiv(3, 3), coordsys(3, 3)
    real(kind=kreal) :: B1(6, ndof*(nn+3))
    real(kind=kreal) :: Smat(9, 9)
    real(kind=kreal) :: BN(9, ndof*(nn+3)), SBN(9, ndof*(nn+3))
    real(kind=kreal) :: spfunc(nn), temp

    !---------------------------------------------------------------------

    fetype = fe_hex8n
    if( present(nddisp) .AND. present(ehdisp) ) then
      unode(:, 1:nn)      = nddisp(:, :)
      unode(:, nn+1:nn+3) = ehdisp(:, :)
    end if

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(nddisp) ) flag = INFINITE    ! enforce to infinite deformation analysis
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+unode(:, 1:nn)

    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(fetype, nn, naturalcoord, elem, det, jacobian, inverse)
    inverse(:, :)= inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    tmpstiff(:, :) = 0.0D0
    B1(1:6, 1:(nn+3)*ndof) = 0.0D0
    BN(1:9, 1:(nn+3)*ndof) = 0.0D0

    do LX = 1, NumOfQuadPoints(fetype)

      call getQuadPoint(fetype, LX, naturalCoord)
      call getGlobalDeriv( fetype, nn, naturalcoord, elem, det, gderiv(1:nn, 1:3) )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      if( present(temperature) ) then
        call getShapeFunc( etype, naturalcoord, spfunc )
        temp = dot_product( temperature, spfunc )
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )
      else
        call MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys )
      end if

      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      endif

      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(fetype, LX)*det
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
      forall( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
        tmpstiff(i, j) = tmpstiff(i, j)+dot_product( B(:, i), DB(:, j) )*wg
      end forall
    end do

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
      forall( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
        tmpstiff(i, j) = tmpstiff(i, j)+dot_product( BN(:, i), SBN(:, j) )*wg
      end forall
    end if

    ! -----Condense tmpstiff to stiff
    xj(1:9, 1:9)= tmpstiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    call calInverse(9, xj)
    tmpk = matmul( tmpstiff( 1:nn*ndof, nn*ndof+1:(nn+3)*ndof ), xj )
    stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)-matmul( tmpk, tmpstiff(nn*ndof+1:(nn+3)*ndof, 1:nn*ndof)  )

  end subroutine STF_C3D8IC


  !>  Update Strain stress of this element
  !----------------------------------------------------------------------*
  subroutine UpdateST_C3D8IC &
      (etype, nn, xx, yy, zz, edisp, &
      gausses, cdsys_ID, coords, qf, tt, t0)
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
    real(kind=kreal), intent(in)      :: edisp(nn*ndof)          !< nodal displacement
    type(tGaussStatus), intent(inout) :: gausses(:)              !< info about qudrature points
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)            !< variables to define matreial coordinate system
    real(kind=kreal), intent(in), optional  :: tt(nn), t0(nn)          !< current and ref temprature
    real(kind=kreal), intent(out), optional :: qf(nn*ndof)

    !---------------------------------------------------------------------

    real(kind=kreal) :: ALP, ALP0
    real(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    real(kind=kreal) :: det, wg, ecoord(3, nn)
    integer(kind=kint) :: j, k, IC, serr, fetype
    real(kind=kreal) :: EPSA(6), EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3), xj(9, 9)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn+3, 3)
    real(kind=kreal) :: jacobian(3, 3),inverse(3, 3)
    real(kind=kreal) :: stiff((nn+3)*3, (nn+3)*3)
    real(kind=kreal) :: tmpforce(9), cdisp((nn+3)*3), vect((nn+3)*3)
    real(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6), outa(1), ina(1)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    matlaniso = .FALSE.
    if( present(tt) .and. present(t0) .and. cdsys_ID > 0 ) then   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .true.
    end if

    fetype = fe_hex8n
    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)
    ! ---- Calculate enhanced displacement at first
    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(fetype, nn, naturalcoord, ecoord, det, jacobian, inverse)
    inverse(:, :) = inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff(:, :) = 0.0D0
    B(1:6, 1:(nn+3)*ndof) = 0.0D0

    do IC = 1, NumOfQuadPoints(fetype)

      call getQuadPoint(fetype, IC, naturalCoord)
      call getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn, 1:3) )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      if( present(tt) .and. present(t0) ) then
        call getShapeFunc( fetype, naturalcoord, H(1:nn) )
        TEMPC = dot_product( H(1:nn), TT(1:nn) )
        call MatlMatrix( gausses(IC), D3, D, 1.d0, 1.0D0, coordsys, TEMPC )
      else
        call MatlMatrix( gausses(IC), D3, D, 1.d0, 1.0D0, coordsys )
      endif

      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(fetype, IC)*det
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

      DB(1:6, 1:(nn+3)*ndof) = matmul( D, B(1:6, 1:(nn+3)*ndof) )
      stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof) &
        = stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof) &
        +matmul( transpose(B(1:6, 1:(nn+3)*ndof)), DB(1:6, 1:(nn+3)*ndof) )*wg
    end do
    xj(1:9,1:9)= stiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    call calInverse(9, xj)
    ! ---  [Kda]*edisp
    tmpforce(:) = matmul( stiff(nn*3+1:(nn+3)*ndof, 1:nn*ndof), edisp )
    cdisp(1:nn*3) = edisp(:)
    ! ---  -[Kaa]-1 * [Kda] * edisp
    cdisp(nn*3+1:(nn+3)*3) = -matmul( xj(:, :), tmpforce(:) )

    if( present(qf) ) then
      qf(1:nn*ndof) = matmul( stiff(1:nn*ndof, 1:(nn+3)*ndof), cdisp(1:(nn+3)*3))
      if( present(tt) .and. present(t0) ) then
        call TLOAD_C3D8IC(etype, nn, xx, yy, zz, tt, t0, gausses, vect, cdsys_ID, coords)
        qf(1:nn*ndof) = qf(1:nn*ndof) - vect(1:nn*ndof)
      end if
    endif

    ! ---- Now strain and stress calculation
    do IC = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, IC, naturalCoord(:) )
      call getShapeFunc( etype, naturalcoord, H(1:nn) )
      call getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      if( matlaniso ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      B(1:6, 1:(nn+3)*ndof) = 0.0D0
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

      if( present(tt) .and. present(t0) ) then

        TEMPC = dot_product( H(1:nn), TT(1:nn) )
        TEMP0 = dot_product( H(1:nn), T0(1:nn) )

        if( cdsys_ID > 0 ) then
          call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
          if( serr == -1 ) stop "Fail to setup local coordinate"
          if( serr == -2 ) then
            write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
          end if
        end if
        call MatlMatrix( gausses(IC), D3, D, 1.d0, 0.0D0, coordsys, TEMPC )

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
      else
        EPSTH = 0.d0
        call MatlMatrix( gausses(IC), D3, D, 1.d0, 0.0D0, coordsys )
      endif
      !**
      !** SET EPS  {e}=[B]{u}
      !**
      EPSA(1:6) = matmul( B(1:6,:), CDISP(:) )
      !**
      !** SET SGM  {S}=[D]{e}
      !**
      do j = 1, 6
        SGM(j) = 0.0D0
        do k = 1, 6
          SGM(j) = SGM(j)+D(j, k)*( EPSA(k)-EPSTH(k) )
        end do
      end do
      !**
      !** Adding stress in each gauss points
      !**
      gausses(IC)%strain(1:6) = EPSA(1:6)
      gausses(IC)%stress(1:6) = SGM(1:6)

    end do

  end subroutine UpdateST_C3D8IC


  !>  This subroutine calculatess thermal loading
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
    real(kind=kreal), intent(in)      :: tt(nn), t0(nn)          !< current and ref temprature
    type(tGaussStatus), intent(inout) :: gausses(:)              !< info about qudrature points
    real(kind=kreal), intent(out)     :: VECT(nn*ndof)           !< load vector
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)            !< variables to define matreial coordinate system

    !---------------------------------------------------------------------

    real(kind=kreal) :: ALP, ALP0
    real(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    real(kind=kreal) :: det, wg, ecoord(3, nn)
    integer(kind=kint) :: j, IC, serr, fetype
    real(kind=kreal) :: EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3), xj(9, 9)
    real(kind=kreal) :: naturalcoord(3), gderiv(nn+3, 3)
    real(kind=kreal) :: jacobian(3, 3),inverse(3, 3)
    real(kind=kreal) :: stiff((nn+3)*3, (nn+3)*3)
    real(kind=kreal) :: tmpvect((nn+3)*3)
    real(kind=kreal) :: tmpforce(nn*ndof), icdisp(9)
    real(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6), outa(1), ina(1)
    logical :: ierr, matlaniso

    !---------------------------------------------------------------------

    matlaniso = .FALSE.
    if( cdsys_ID > 0 ) then   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      if( .not. ierr ) matlaniso = .true.
    end if

    fetype = fe_hex8n
    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)
    ! ---- Calculate enhanced displacement at first
    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    call getJacobian(fetype, nn, naturalcoord, ecoord, det, jacobian, inverse)
    inverse(:, :) = inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff(:, :) = 0.0D0
    B(1:6, 1:(nn+3)*ndof) = 0.0D0
    tmpvect(1:(nn+3)*ndof) = 0.0D0

    do IC = 1, NumOfQuadPoints(fetype)

      call getQuadPoint(fetype, IC, naturalCoord)
      call getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn, 1:3) )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) then
          write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        end if
      end if

      call getShapeFunc( fetype, naturalcoord, H(1:nn) )
      TEMPC = dot_product( H(1:nn), TT(1:nn) )
      TEMP0 = dot_product( H(1:nn), T0(1:nn) )
      call MatlMatrix( gausses(IC), D3, D, 1.d0, 1.0D0, coordsys, TEMPC )

      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(fetype, IC)*det
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

      DB(1:6, 1:(nn+3)*ndof) = matmul( D, B(1:6, 1:(nn+3)*ndof) )
      stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof) &
        = stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof) &
        +matmul( transpose(B(1:6, 1:(nn+3)*ndof)), DB(1:6, 1:(nn+3)*ndof) )*wg

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
      TMPVECT(1:(nn+3)*ndof) = TMPVECT(1:(nn+3)*ndof)+matmul( SGM(1:6), B(1:6, 1:(nn+3)*ndof) )*wg

    end do

    ! --- condense tmpvect to vect
    xj(1:9,1:9)= stiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    call calInverse(9, xj)
    ! ---  [Kaa]-1 * fa
    icdisp(1:9) = matmul( xj(:, :), tmpvect(nn*ndof+1:(nn+3)*ndof) )
    ! ---  [Kda] * [Kaa]-1 * fa
    tmpforce(1:nn*ndof) = matmul( stiff(1:nn*ndof, nn*3+1:(nn+3)*ndof), icdisp(1:9) )
    ! ---  fd - [Kda] * [Kaa]-1 * fa
    vect(1:nn*ndof) = tmpvect(1:nn*ndof) - tmpforce(1:nn*ndof)

  end subroutine TLOAD_C3D8IC


end module m_static_LIB_3dIC
