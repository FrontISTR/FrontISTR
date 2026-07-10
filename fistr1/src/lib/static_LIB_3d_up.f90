!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions of u-p mixed (UP) solid elements
module m_static_LIB_3d_up

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

contains

  !>  This subroutine calculate stiff matrix of u-p mixed solid elements
  subroutine STF_C3_up                                        &
      (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords,   &
      time, tincr, nn_p, lambda, u, temperature)

    use mMechGauss
    use m_MatMatrix
    use m_common_struct
    use m_static_LIB_3d, only: GEOMAT_C3

    integer(kind=kint), intent(in)  :: etype                   !< element type
    integer(kind=kint), intent(in)  :: nn                      !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3, nn)           !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)              !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:, :)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3, 3)            !< variables to define material coordinate system
    real(kind=kreal),   intent(in)  :: time                    !< current time
    real(kind=kreal),   intent(in)  :: tincr                   !< time increment
    integer(kind=kint), intent(in)  :: nn_p                    !< number of elemental pressure nodes
    real(kind=kreal),   intent(in)  :: lambda(nn_p)            !< Lagrange multiplier
    real(kind=kreal),   intent(in), optional :: u(:, :)        !< nodal displacemwent
    real(kind=kreal),   intent(in), optional :: temperature(nn) !< temperature

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: d(6, 6), b(6, ndof*nn), db(6, ndof*nn)
    real(kind=kreal) :: gderiv(nn, 3), stress(6), mat(6, 6)
    real(kind=kreal) :: det, wg, temp, spfunc(nn)
    integer(kind=kint) :: i, j, lx, serr
    real(kind=kreal) :: naturalcoord(3), coordsys(3, 3)
    real(kind=kreal) :: gdispderiv(3, 3)
    real(kind=kreal) :: b1(6, ndof*nn)
    real(kind=kreal) :: smat(9, 9), elem(3, nn)
    real(kind=kreal) :: bn(9, ndof*nn), sbn(9, ndof*nn)
    integer(kind=kint) :: na, nb
    integer(kind=kint) :: na_p, nb_p
    integer(kind=kint) :: isize, jsize
    integer(kind=kint) :: jsize1, jsize2, jsize3
    real(kind=kreal) :: stiff_up(3*nn, nn_p)
    real(kind=kreal) :: stiff_pp(nn_p, nn_p), stiff_pp_inv(nn_p, nn_p)
    real(kind=kreal) :: stiff_up_stiff_pp_inv(3*nn, nn_p)
    real(kind=kreal) :: stiff_up_stiff_pp_inv_stiff_up(3*nn, 3*nn)
    real(kind=kreal) :: alpha_inv
    real(kind=kreal) :: g
    real(kind=kreal) :: d2(6)
    real(kind=kreal) :: bd2(3*nn)
    real(kind=kreal) :: lambda_lx
    real(kind=kreal) :: spfunc_p(nn_p)

    stiff(:, :) = 0.0D0
    stiff_up(:, :) = 0.0D0
    stiff_pp(:, :) = 0.0D0

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis
    elem(:, :) = ecoord(:, :)
    ! Updated Lagrangian : evaluate on the current (deformed) configuration
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)
    if( flag == INFINITESIMAL ) gdispderiv(:, :) = 0.0D0

    coordsys = 0.0D0
    coordsys(1, 1) = 1.0D0; coordsys(2, 2) = 1.0D0; coordsys(3, 3) = 1.0D0

    do lx = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, lx, naturalcoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
      end if

      ! Displacement gradient (B1 nonlinear term : Total Lagrangian only)
      if( flag == TOTALLAG ) then
        gdispderiv(1:3, 1:3) = matmul( u(1:3, 1:nn), gderiv(1:nn, 1:3) )
      else
        gdispderiv(1:3, 1:3) = 0.0D0
      end if

      if( nn_p == 1 ) spfunc_p(1) = 1.0D0

      ! Lagrange multiplier
      lambda_lx = 0.0D0
      do na_p = 1, nn_p
        lambda_lx = lambda_lx+spfunc_p(na_p)*lambda(na_p)
      end do

      if( present( temperature ) ) then
        call getShapeFunc( etype, naturalcoord, spfunc )
        temp = dot_product( temperature, spfunc )
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys, temp )
      else
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys )
      end if

      ! Updated Lagrangian : material tangent correction with the Cauchy stress
      if( flag == UPDATELAG ) then
        call GEOMAT_C3( gausses(lx)%stress, mat )
        d(:, :) = d(:, :)-mat
      end if

      wg = getWeight( etype, lx )*det

      do nb = 1, nn
        jsize1 = 3*(nb-1)+1
        jsize2 = 3*(nb-1)+2
        jsize3 = 3*(nb-1)+3
        b(1, jsize1) = gderiv(nb, 1)
        b(2, jsize1) = 0.0D0
        b(3, jsize1) = 0.0D0
        b(4, jsize1) = gderiv(nb, 2)
        b(5, jsize1) = 0.0D0
        b(6, jsize1) = gderiv(nb, 3)
        b(1, jsize2) = 0.0D0
        b(2, jsize2) = gderiv(nb, 2)
        b(3, jsize2) = 0.0D0
        b(4, jsize2) = gderiv(nb, 1)
        b(5, jsize2) = gderiv(nb, 3)
        b(6, jsize2) = 0.0D0
        b(1, jsize3) = 0.0D0
        b(2, jsize3) = 0.0D0
        b(3, jsize3) = gderiv(nb, 3)
        b(4, jsize3) = 0.0D0
        b(5, jsize3) = gderiv(nb, 2)
        b(6, jsize3) = gderiv(nb, 1)
      end do

      do nb = 1, nn
        jsize1 = 3*(nb-1)+1
        jsize2 = 3*(nb-1)+2
        jsize3 = 3*(nb-1)+3
        b1(1, jsize1) = gdispderiv(1, 1)*gderiv(nb, 1)
        b1(2, jsize1) = gdispderiv(1, 2)*gderiv(nb, 2)
        b1(3, jsize1) = gdispderiv(1, 3)*gderiv(nb, 3)
        b1(4, jsize1) = gdispderiv(1, 2)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 2)
        b1(5, jsize1) = gdispderiv(1, 2)*gderiv(nb, 3)+gdispderiv(1, 3)*gderiv(nb, 2)
        b1(6, jsize1) = gdispderiv(1, 3)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 3)
        b1(1, jsize2) = gdispderiv(2, 1)*gderiv(nb, 1)
        b1(2, jsize2) = gdispderiv(2, 2)*gderiv(nb, 2)
        b1(3, jsize2) = gdispderiv(2, 3)*gderiv(nb, 3)
        b1(4, jsize2) = gdispderiv(2, 2)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 2)
        b1(5, jsize2) = gdispderiv(2, 2)*gderiv(nb, 3)+gdispderiv(2, 3)*gderiv(nb, 2)
        b1(6, jsize2) = gdispderiv(2, 3)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 3)
        b1(1, jsize3) = gdispderiv(3, 1)*gderiv(nb, 1)
        b1(2, jsize3) = gdispderiv(3, 2)*gderiv(nb, 2)
        b1(3, jsize3) = gdispderiv(3, 3)*gderiv(nb, 3)
        b1(4, jsize3) = gdispderiv(3, 2)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 2)
        b1(5, jsize3) = gdispderiv(3, 2)*gderiv(nb, 3)+gdispderiv(3, 3)*gderiv(nb, 2)
        b1(6, jsize3) = gdispderiv(3, 3)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 3)
      end do

      do jsize = 1, 3*nn
        b(:, jsize) = b(:, jsize)+b1(:, jsize)
      end do

      db(1:6, 1:3*nn) = matmul( d(1:6, 1:6), b(1:6, 1:3*nn) )

      ! [ Kuu ] matrix
      forall( isize = 1:3*nn, jsize = 1:3*nn )
        stiff(isize, jsize) = stiff(isize, jsize)+wg*dot_product( b(:, isize), db(:, jsize) )
      end forall

      ! calculate the initial stress matrix (Total/Updated Lagrangian)
      if( flag == TOTALLAG .or. flag == UPDATELAG ) then
        stress(1:6) = gausses(lx)%stress

        do nb = 1, nn
          jsize1 = 3*(nb-1)+1
          jsize2 = 3*(nb-1)+2
          jsize3 = 3*(nb-1)+3
          bn(1, jsize1) = gderiv(nb, 1)
          bn(2, jsize1) = 0.0D0
          bn(3, jsize1) = 0.0D0
          bn(4, jsize1) = gderiv(nb, 2)
          bn(5, jsize1) = 0.0D0
          bn(6, jsize1) = 0.0D0
          bn(7, jsize1) = gderiv(nb, 3)
          bn(8, jsize1) = 0.0D0
          bn(9, jsize1) = 0.0D0
          bn(1, jsize2) = 0.0D0
          bn(2, jsize2) = gderiv(nb, 1)
          bn(3, jsize2) = 0.0D0
          bn(4, jsize2) = 0.0D0
          bn(5, jsize2) = gderiv(nb, 2)
          bn(6, jsize2) = 0.0D0
          bn(7, jsize2) = 0.0D0
          bn(8, jsize2) = gderiv(nb, 3)
          bn(9, jsize2) = 0.0D0
          bn(1, jsize3) = 0.0D0
          bn(2, jsize3) = 0.0D0
          bn(3, jsize3) = gderiv(nb, 1)
          bn(4, jsize3) = 0.0D0
          bn(5, jsize3) = 0.0D0
          bn(6, jsize3) = gderiv(nb, 2)
          bn(7, jsize3) = 0.0D0
          bn(8, jsize3) = 0.0D0
          bn(9, jsize3) = gderiv(nb, 3)
        end do

        smat(:, :) = 0.0D0
        do j = 1, 3
          smat(j  , j  ) = stress(1)
          smat(j  , j+3) = stress(4)
          smat(j  , j+6) = stress(6)
          smat(j+3, j  ) = stress(4)
          smat(j+3, j+3) = stress(2)
          smat(j+3, j+6) = stress(5)
          smat(j+6, j  ) = stress(6)
          smat(j+6, j+3) = stress(5)
          smat(j+6, j+6) = stress(3)
        end do

        sbn(1:9, 1:3*nn) = matmul( smat(1:9, 1:9), bn(1:9, 1:3*nn) )

        forall( isize = 1:3*nn, jsize = 1:3*nn )
          stiff(isize, jsize) = stiff(isize, jsize)+wg*dot_product( bn(:, isize), sbn(:, jsize) )
        end forall
      end if

      do isize = 1, 3*nn
        bd2(isize) = dot_product( b(:, isize), d2(:) )
      end do

      ! [ Kup ] matrix
      forall( isize = 1:3*nn, nb_p = 1:nn_p )
        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)
      end forall

      ! [ Kpp ] matrix
      forall( na_p = 1:nn_p, nb_p = 1:nn_p )
        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)
      end forall

    end do ! gauss loop

    ! [ Kpp ]^-1 matrix (constant pressure: nn_p = 1)
    if( nn_p == 1 ) then
      stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)
    else
      write(6, *) 'Error: nn_p should be equal to 1.'
      return
    end if

    stiff_up_stiff_pp_inv(1:3*nn, 1:nn_p) = matmul( stiff_up(1:3*nn, 1:nn_p), stiff_pp_inv(1:nn_p, 1:nn_p) )

    ! [ Kup ] [ Kpp ]^-1 [ Kup ]^T matrix
    forall( isize = 1:3*nn, jsize = 1:3*nn )
      stiff_up_stiff_pp_inv_stiff_up(isize, jsize) = dot_product( stiff_up_stiff_pp_inv(isize, :), stiff_up(jsize, :) )
    end forall

    ! static condensation : [ Kuu ] - [ Kup ] [ Kpp ]^-1 [ Kup ]^T
    forall( isize = 1:3*nn, jsize = 1:3*nn )
      stiff(isize, jsize) = stiff(isize, jsize)-stiff_up_stiff_pp_inv_stiff_up(isize, jsize)
    end forall

  end subroutine STF_C3_up


  !>  Update strain and stress inside u-p mixed (UP) solid element
  subroutine Update_C3_up                                     &
      (etype, nn, ecoord, u, du, ddu, cdsys_ID, coords, qf,   &
      gausses, iter, time, tincr,                             &
      nn_p, lambda, ddlambda, tt, t0)

    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use mHyperElastic
    use m_common_struct
    use m_utilities

    integer(kind=kint), intent(in)    :: etype          !< \param [in] element type
    integer(kind=kint), intent(in)    :: nn             !< \param [in] number of elemental nodes
    real(kind=kreal),   intent(in)    :: ecoord(3, nn)  !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   intent(in)    :: u(3, nn)       !< \param [in] nodal dislplacements
    real(kind=kreal),   intent(in)    :: du(3, nn)      !< \param [in] nodal displacement increment of the step
    real(kind=kreal),   intent(in)    :: ddu(3, nn)     !< \param [in] nodal displacement increment of the iteration
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)   !< variables to define material coordinate system
    real(kind=kreal),   intent(out)   :: qf(nn*3)       !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)     !< \param [out] status of qudrature points
    integer, intent(in)               :: iter
    real(kind=kreal),   intent(in)    :: time           !< current time
    real(kind=kreal),   intent(in)    :: tincr          !< time increment
    integer(kind=kint), intent(in)    :: nn_p
    real(kind=kreal),   intent(in)    :: lambda(nn_p)   !< Lagrange multiplier
    real(kind=kreal),   intent(inout) :: ddlambda(nn_p)
    real(kind=kreal),   intent(in), optional :: tt(nn)  !< current temperature
    real(kind=kreal),   intent(in), optional :: t0(nn)  !< reference temperature

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: d(6, 6), b(6, ndof*nn), b1(6, ndof*nn)
    real(kind=kreal) :: gderiv(nn, 3), gdispderiv(3, 3), det, wg
    integer(kind=kint) :: i, j, lx, mtype, serr
    real(kind=kreal) :: naturalcoord(3), rot(3, 3), spfunc(nn), coordsys(3, 3)
    real(kind=kreal) :: totaldisp(3, nn), elem(3, nn), elem1(3, nn)
    real(kind=kreal) :: dstrain(6), dstress(6), dumstress(3, 3), dum(3, 3)
    real(kind=kreal) :: trd, p_bak
    real(kind=kreal) :: ttc, tt0, outa(1), ina(1), epsth(6)
    logical :: ierr
    integer(kind=kint) :: na, nb
    integer(kind=kint) :: na_p, nb_p
    integer(kind=kint) :: isize, jsize
    integer(kind=kint) :: jsize1, jsize2, jsize3
    real(kind=kreal) :: totallambda(nn_p)
    real(kind=kreal) :: stiff_up(3*nn, nn_p)
    real(kind=kreal) :: stiff_pp(nn_p, nn_p), stiff_pp_inv(nn_p, nn_p)
    real(kind=kreal) :: stiff_up_stiff_pp_inv(3*nn, nn_p)
    real(kind=kreal) :: stiff_up_stiff_pp_inv_stiff_up(3*nn, 3*nn)
    real(kind=kreal) :: stiff_up_stiff_pp_inv_qf_p(3*nn)
    real(kind=kreal) :: alpha_inv
    real(kind=kreal) :: g
    real(kind=kreal) :: qf_p(nn_p)
    real(kind=kreal) :: stiff_up_ddu(nn_p)
    real(kind=kreal) :: d2(6)
    real(kind=kreal) :: bd2(3*nn)
    real(kind=kreal) :: lambda_lx
    real(kind=kreal) :: spfunc_p(nn_p)

    stiff_up(:, :) = 0.0D0
    stiff_pp(:, :) = 0.0D0
    qf_p(:) = 0.0D0

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:, :) = ecoord(:, :)
    totaldisp(:, :) = u(:, :)+( du(:, :)-ddu(:, :) )
    ! Updated Lagrangian : 1st pass evaluated on the previous-iterate configuration
    if( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)+( du(:, :)-ddu(:, :) )
    if( flag == INFINITESIMAL ) gdispderiv(:, :) = 0.0D0

    coordsys = 0.0D0
    coordsys(1, 1) = 1.0D0; coordsys(2, 2) = 1.0D0; coordsys(3, 3) = 1.0D0

    !========================================================
    ! 1st pass : evaluate [ Kup ], [ Kpp ] and the pressure
    !            residual to condense the increment ddlambda
    !========================================================
    do lx = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, lx, naturalcoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
      end if

      ! Displacement gradient (B1 nonlinear term : Total Lagrangian only)
      if( flag == TOTALLAG ) gdispderiv(1:3, 1:3) = matmul( totaldisp(1:3, 1:nn), gderiv(1:nn, 1:3) )

      if( nn_p == 1 ) spfunc_p(1) = 1.0D0

      ! Lagrange multiplier
      lambda_lx = 0.0D0
      do na_p = 1, nn_p
        lambda_lx = lambda_lx+spfunc_p(na_p)*lambda(na_p)
      end do

      mtype = gausses(lx)%pMaterial%mtype

      wg = getWeight( etype, lx )*det

      ! Update strain and stress : plasticity/creep are treated as elastic in UP
      if( isElastoplastic(gausses(lx)%pMaterial%mtype) .or. &
        ( gausses(lx)%pMaterial%mtype == NORTON ) ) gausses(lx)%pMaterial%mtype = ELASTIC

      epsth = 0.0D0

      if( present( tt ) .and. present( t0 ) ) then
        call getShapeFunc( etype, naturalcoord, spfunc )
        ttc = dot_product( tt, spfunc )
        tt0 = dot_product( t0, spfunc )
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys, ttc )
        if( ( iter <= 1 ) .or. ( flag == TOTALLAG ) ) then
          ina(1) = ttc
          call fetch_TableData( MC_THEMOEXP, gausses(lx)%pMaterial%dict, outa, ierr, ina )
          if( ierr ) outa(1) = gausses(lx)%pMaterial%variables(M_EXAPNSION)
          epsth(1:3) = outa(1)*( ttc-tt0 )
        end if
      else
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys )
      end if

      do nb = 1, nn
        jsize1 = 3*(nb-1)+1
        jsize2 = 3*(nb-1)+2
        jsize3 = 3*(nb-1)+3
        b(1, jsize1) = gderiv(nb, 1)
        b(2, jsize1) = 0.0D0
        b(3, jsize1) = 0.0D0
        b(4, jsize1) = gderiv(nb, 2)
        b(5, jsize1) = 0.0D0
        b(6, jsize1) = gderiv(nb, 3)
        b(1, jsize2) = 0.0D0
        b(2, jsize2) = gderiv(nb, 2)
        b(3, jsize2) = 0.0D0
        b(4, jsize2) = gderiv(nb, 1)
        b(5, jsize2) = gderiv(nb, 3)
        b(6, jsize2) = 0.0D0
        b(1, jsize3) = 0.0D0
        b(2, jsize3) = 0.0D0
        b(3, jsize3) = gderiv(nb, 3)
        b(4, jsize3) = 0.0D0
        b(5, jsize3) = gderiv(nb, 2)
        b(6, jsize3) = gderiv(nb, 1)
      end do

      ! B1 nonlinear term (Total Lagrangian only)
      if( flag == TOTALLAG ) then
        do nb = 1, nn
          jsize1 = 3*(nb-1)+1
          jsize2 = 3*(nb-1)+2
          jsize3 = 3*(nb-1)+3
          b1(1, jsize1) = gdispderiv(1, 1)*gderiv(nb, 1)
          b1(2, jsize1) = gdispderiv(1, 2)*gderiv(nb, 2)
          b1(3, jsize1) = gdispderiv(1, 3)*gderiv(nb, 3)
          b1(4, jsize1) = gdispderiv(1, 2)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 2)
          b1(5, jsize1) = gdispderiv(1, 2)*gderiv(nb, 3)+gdispderiv(1, 3)*gderiv(nb, 2)
          b1(6, jsize1) = gdispderiv(1, 3)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 3)
          b1(1, jsize2) = gdispderiv(2, 1)*gderiv(nb, 1)
          b1(2, jsize2) = gdispderiv(2, 2)*gderiv(nb, 2)
          b1(3, jsize2) = gdispderiv(2, 3)*gderiv(nb, 3)
          b1(4, jsize2) = gdispderiv(2, 2)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 2)
          b1(5, jsize2) = gdispderiv(2, 2)*gderiv(nb, 3)+gdispderiv(2, 3)*gderiv(nb, 2)
          b1(6, jsize2) = gdispderiv(2, 3)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 3)
          b1(1, jsize3) = gdispderiv(3, 1)*gderiv(nb, 1)
          b1(2, jsize3) = gdispderiv(3, 2)*gderiv(nb, 2)
          b1(3, jsize3) = gdispderiv(3, 3)*gderiv(nb, 3)
          b1(4, jsize3) = gdispderiv(3, 2)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 2)
          b1(5, jsize3) = gdispderiv(3, 2)*gderiv(nb, 3)+gdispderiv(3, 3)*gderiv(nb, 2)
          b1(6, jsize3) = gdispderiv(3, 3)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 3)
        end do
        do jsize = 1, 3*nn
          b(:, jsize) = b(:, jsize)+b1(:, jsize)
        end do
      end if

      do isize = 1, 3*nn
        bd2(isize) = dot_product( b(:, isize), d2 )
      end do

      ! [ Kup ] matrix
      forall( isize = 1:3*nn, nb_p = 1:nn_p )
        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)
      end forall

      ! [ Kpp ] matrix
      forall( na_p = 1:nn_p, nb_p = 1:nn_p )
        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)
      end forall

      qf_p(1:nn_p) = qf_p(1:nn_p)+wg*spfunc_p(1:nn_p)*( g-alpha_inv*lambda_lx )

    end do ! gauss loop

    ! [ Kpp ]^-1 matrix (constant pressure: nn_p = 1)
    if( nn_p == 1 ) then
      stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)
    else
      write(6, *) 'Error: nn_p should be equal to 1.'
      return
    end if

    do na_p = 1, nn_p
      stiff_up_ddu(na_p) = 0.0D0
      do nb = 1, nn
        do i = 1, 3
          jsize = 3*(nb-1)+i
          stiff_up_ddu(na_p) = stiff_up_ddu(na_p)+stiff_up(jsize, na_p)*ddu(i, nb)
        end do
      end do
    end do

    do na_p = 1, nn_p
      ddlambda(na_p) = dot_product( stiff_pp_inv(na_p, :), -qf_p-stiff_up_ddu )
    end do

    totaldisp(:, :) = u(:, :)+du(:, :)
    totallambda(:) = lambda(:)+ddlambda(:)

    stiff_up(:, :) = 0.0D0
    stiff_pp(:, :) = 0.0D0
    qf(:) = 0.0D0
    qf_p(:) = 0.0D0

    ! Updated Lagrangian : strain increment on the mid configuration, internal
    ! force assembled on the end (current) configuration
    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) then
      elem(:, :)  = ecoord(:, :)+u(:, :)+0.5D0*du(:, :)
      elem1(:, :) = ecoord(:, :)+u(:, :)+du(:, :)
    end if

    !========================================================
    ! 2nd pass : update strain/stress and assemble the
    !            condensed internal force at the new state
    !========================================================
    do lx = 1, NumOfQuadPoints(etype)

      call getQuadPoint( etype, lx, naturalcoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

      if( cdsys_ID > 0 ) then
        call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        if( serr == -1 ) stop "Fail to setup local coordinate"
        if( serr == -2 ) write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
      end if

      ! displacement gradient : step increment (UPDATELAG) or total (others)
      if( flag == UPDATELAG ) then
        gdispderiv(1:3, 1:3) = matmul( du(1:3, 1:nn), gderiv(1:nn, 1:3) )
      else
        gdispderiv(1:3, 1:3) = matmul( totaldisp(1:3, 1:nn), gderiv(1:nn, 1:3) )
      end if

      if( nn_p == 1 ) spfunc_p(1) = 1.0D0

      ! Lagrange multiplier
      lambda_lx = 0.0D0
      do na_p = 1, nn_p
        lambda_lx = lambda_lx+spfunc_p(na_p)*totallambda(na_p)
      end do

      mtype = gausses(lx)%pMaterial%mtype

      wg = getWeight( etype, lx )*det

      ! Small (linear) strain
      dstrain(1) = gdispderiv(1, 1)
      dstrain(2) = gdispderiv(2, 2)
      dstrain(3) = gdispderiv(3, 3)
      dstrain(4) = gdispderiv(1, 2)+gdispderiv(2, 1)
      dstrain(5) = gdispderiv(2, 3)+gdispderiv(3, 2)
      dstrain(6) = gdispderiv(3, 1)+gdispderiv(1, 3)
      dstrain(:) = dstrain(:)-epsth(:)

      ! Green-Lagrange strain (Total Lagrangian only)
      if( flag == TOTALLAG ) then
        dstrain(1) = dstrain(1)+0.5D0*dot_product( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5D0*dot_product( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5D0*dot_product( gdispderiv(:, 3), gdispderiv(:, 3) )
        dstrain(4) = dstrain(4)+dot_product( gdispderiv(:, 1), gdispderiv(:, 2) )
        dstrain(5) = dstrain(5)+dot_product( gdispderiv(:, 2), gdispderiv(:, 3) )
        dstrain(6) = dstrain(6)+dot_product( gdispderiv(:, 1), gdispderiv(:, 3) )
      end if

      if( flag == UPDATELAG ) then
        ! spin (skew part) for the objective stress rate
        rot = 0.0D0
        rot(1, 2) = 0.5D0*( gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3) = 0.5D0*( gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3) = 0.5D0*( gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)
        gausses(lx)%strain(1:6) = gausses(lx)%strain_bak(1:6)+dstrain(1:6)+epsth(:)
      else
        gausses(lx)%strain(1:6) = dstrain(1:6)+epsth(:)
      end if

      if( flag == UPDATELAG ) then
        ! deviatoric tangent (K=0) of the increment
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys )
        ! objective (Jaumann/Hughes-Winget) update of the deviatoric stress.
        ! pressure is imposed by the Lagrange multiplier (UP), not by the rate.
        p_bak = ( gausses(lx)%stress_bak(1)+gausses(lx)%stress_bak(2)+gausses(lx)%stress_bak(3) )/3.0D0
        dumstress(1, 1) = gausses(lx)%stress_bak(1)-p_bak
        dumstress(2, 2) = gausses(lx)%stress_bak(2)-p_bak
        dumstress(3, 3) = gausses(lx)%stress_bak(3)-p_bak
        dumstress(1, 2) = gausses(lx)%stress_bak(4);  dumstress(2, 1) = dumstress(1, 2)
        dumstress(2, 3) = gausses(lx)%stress_bak(5);  dumstress(3, 2) = dumstress(2, 3)
        dumstress(3, 1) = gausses(lx)%stress_bak(6);  dumstress(1, 3) = dumstress(3, 1)
        trd = dstrain(1)+dstrain(2)+dstrain(3)
        dum(:, :) = dumstress+matmul( rot, dumstress )-matmul( dumstress, rot )-dumstress*trd
        dstress(1:6) = matmul( d(1:6, 1:6), dstrain(1:6) )
        gausses(lx)%stress(1) = dum(1, 1)+dstress(1)+lambda_lx
        gausses(lx)%stress(2) = dum(2, 2)+dstress(2)+lambda_lx
        gausses(lx)%stress(3) = dum(3, 3)+dstress(3)+lambda_lx
        gausses(lx)%stress(4) = dum(1, 2)+dstress(4)
        gausses(lx)%stress(5) = dum(2, 3)+dstress(5)
        gausses(lx)%stress(6) = dum(3, 1)+dstress(6)
      else
        call StressUpdate_up( gausses(lx), D3, dstrain, gausses(lx)%stress, lambda_lx, g, coordsys )
      end if

      gausses(lx)%stress_out(1:6) = gausses(lx)%stress(1:6)
      gausses(lx)%strain_out(1:6) = gausses(lx)%strain(1:6)

      epsth = 0.0D0

      if( present( tt ) .and. present( t0 ) ) then
        call getShapeFunc( etype, naturalcoord, spfunc )
        ttc = dot_product( tt, spfunc )
        tt0 = dot_product( t0, spfunc )
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys, ttc )
        if( ( iter <= 1 ) .or. ( flag == TOTALLAG ) ) then
          ina(1) = ttc
          call fetch_TableData( MC_THEMOEXP, gausses(lx)%pMaterial%dict, outa, ierr, ina )
          if( ierr ) outa(1) = gausses(lx)%pMaterial%variables(M_EXAPNSION)
          epsth(1:3) = outa(1)*( ttc-tt0 )
        end if
      else
        call MatlMatrix_up( gausses(lx), D3, d, lambda_lx, d2, alpha_inv, g, time, tincr, coordsys )
      end if

      ! Updated Lagrangian : switch to the end (current) configuration for the
      ! B-matrix and the internal force / coupling assembly
      if( flag == UPDATELAG ) then
        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv )
        wg = getWeight( etype, lx )*det
      end if

      do nb = 1, nn
        jsize1 = 3*(nb-1)+1
        jsize2 = 3*(nb-1)+2
        jsize3 = 3*(nb-1)+3
        b(1, jsize1) = gderiv(nb, 1)
        b(2, jsize1) = 0.0D0
        b(3, jsize1) = 0.0D0
        b(4, jsize1) = gderiv(nb, 2)
        b(5, jsize1) = 0.0D0
        b(6, jsize1) = gderiv(nb, 3)
        b(1, jsize2) = 0.0D0
        b(2, jsize2) = gderiv(nb, 2)
        b(3, jsize2) = 0.0D0
        b(4, jsize2) = gderiv(nb, 1)
        b(5, jsize2) = gderiv(nb, 3)
        b(6, jsize2) = 0.0D0
        b(1, jsize3) = 0.0D0
        b(2, jsize3) = 0.0D0
        b(3, jsize3) = gderiv(nb, 3)
        b(4, jsize3) = 0.0D0
        b(5, jsize3) = gderiv(nb, 2)
        b(6, jsize3) = gderiv(nb, 1)
      end do

      ! B1 nonlinear term (Total Lagrangian only)
      if( flag == TOTALLAG ) then
        do nb = 1, nn
          jsize1 = 3*(nb-1)+1
          jsize2 = 3*(nb-1)+2
          jsize3 = 3*(nb-1)+3
          b1(1, jsize1) = gdispderiv(1, 1)*gderiv(nb, 1)
          b1(2, jsize1) = gdispderiv(1, 2)*gderiv(nb, 2)
          b1(3, jsize1) = gdispderiv(1, 3)*gderiv(nb, 3)
          b1(4, jsize1) = gdispderiv(1, 2)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 2)
          b1(5, jsize1) = gdispderiv(1, 2)*gderiv(nb, 3)+gdispderiv(1, 3)*gderiv(nb, 2)
          b1(6, jsize1) = gdispderiv(1, 3)*gderiv(nb, 1)+gdispderiv(1, 1)*gderiv(nb, 3)
          b1(1, jsize2) = gdispderiv(2, 1)*gderiv(nb, 1)
          b1(2, jsize2) = gdispderiv(2, 2)*gderiv(nb, 2)
          b1(3, jsize2) = gdispderiv(2, 3)*gderiv(nb, 3)
          b1(4, jsize2) = gdispderiv(2, 2)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 2)
          b1(5, jsize2) = gdispderiv(2, 2)*gderiv(nb, 3)+gdispderiv(2, 3)*gderiv(nb, 2)
          b1(6, jsize2) = gdispderiv(2, 3)*gderiv(nb, 1)+gdispderiv(2, 1)*gderiv(nb, 3)
          b1(1, jsize3) = gdispderiv(3, 1)*gderiv(nb, 1)
          b1(2, jsize3) = gdispderiv(3, 2)*gderiv(nb, 2)
          b1(3, jsize3) = gdispderiv(3, 3)*gderiv(nb, 3)
          b1(4, jsize3) = gdispderiv(3, 2)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 2)
          b1(5, jsize3) = gdispderiv(3, 2)*gderiv(nb, 3)+gdispderiv(3, 3)*gderiv(nb, 2)
          b1(6, jsize3) = gdispderiv(3, 3)*gderiv(nb, 1)+gdispderiv(3, 1)*gderiv(nb, 3)
        end do
        do jsize = 1, 3*nn
          b(:, jsize) = b(:, jsize)+b1(:, jsize)
        end do
      end if

      do isize = 1, 3*nn
        bd2(isize) = dot_product( b(:, isize), d2 )
      end do

      ! [ Kup ] matrix
      forall( isize = 1:3*nn, nb_p = 1:nn_p )
        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)
      end forall

      ! [ Kpp ] matrix
      forall( na_p = 1:nn_p, nb_p = 1:nn_p )
        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)
      end forall

      ! calculate the Internal Force
      qf(1:3*nn) = qf(1:3*nn)+wg*matmul( gausses(lx)%stress(1:6), b(1:6, 1:3*nn) )

      qf_p(1:nn_p) = qf_p(1:nn_p)+wg*spfunc_p(1:nn_p)*( g-alpha_inv*lambda_lx )

    end do ! gauss loop

    ! [ Kpp ]^-1 matrix (constant pressure: nn_p = 1)
    if( nn_p == 1 ) then
      stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)
    else
      write(6, *) 'Error: nn_p should be equal to 1.'
      return
    end if

    stiff_up_stiff_pp_inv(1:3*nn, 1:nn_p) = matmul( stiff_up(1:3*nn, 1:nn_p), stiff_pp_inv(1:nn_p, 1:nn_p) )

    do isize = 1, 3*nn
      stiff_up_stiff_pp_inv_qf_p(isize) = dot_product( stiff_up_stiff_pp_inv(isize, :), qf_p )
    end do

    ! condensed internal force : qf - [ Kup ] [ Kpp ]^-1 qf_p
    do isize = 1, 3*nn
      qf(isize) = qf(isize)-stiff_up_stiff_pp_inv_qf_p(isize)
    end do

  end subroutine Update_C3_up

end module m_static_LIB_3d_up
