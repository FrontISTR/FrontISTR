      MODULE m_static_LIB_3d_up
!####################################################################

      USE hecmw, ONLY : kint, kreal
      USE elementInfo

!--------------------------------------------------------------------

      IMPLICIT NONE

!--------------------------------------------------------------------

      CONTAINS


!####################################################################
      SUBROUTINE STF_C3_up                                   &
                 ( etype, nn, ecoord, gausses, stiff, tincr, &
                   nn_p, lambda, u, temperature )
!####################################################################

      USE mMechGauss
      USE m_MatMatrix
      USE m_common_struct

!--------------------------------------------------------------------

      INTEGER( KIND = kint ), INTENT(IN) :: etype                   !< element type
      INTEGER( KIND = kint ), INTENT(IN) :: nn                      !< number of elemental nodes
      REAL( KIND = kreal ), INTENT(IN) :: ecoord(3, nn)             !< coordinates of elemental nodes
      TYPE( tGaussStatus ), INTENT(IN) :: gausses(:)                !< status of qudrature points
      REAL( KIND = kreal ), INTENT(OUT) :: stiff(:, :)              !< stiff matrix
      REAL( KIND = kreal ), INTENT(IN) :: tincr                     !< time increment
      REAL( KIND = kreal ), INTENT(IN), OPTIONAL :: u(:, :)         !< nodal displacemwent
      REAL( KIND = kreal ), INTENT(IN), OPTIONAL :: temperature(nn) !< temperature

!--------------------------------------------------------------------

      INTEGER( KIND = kint ), INTENT(IN) :: nn_p                    !< number of elemental pressure nodes
      REAL( KIND = kreal ), INTENT(IN) :: lambda(nn_p)              !< Lagrange multiplier

!--------------------------------------------------------------------

      INTEGER :: flag
      INTEGER( KIND = kint ), PARAMETER :: ndof = 3
      REAL( KIND = kreal ) d(6, 6), d_p(6, 6), b(6, ndof*nn), db(6, ndof*nn)
      REAL( KIND = kreal ) gderiv(nn, 3), stress(6), mat(6, 6)
      REAL( KIND = kreal ) det, wg, temp, spfunc(nn)
      INTEGER( KIND = kint ) i, j, lx
      REAL( KIND = kreal ) naturalcoord(3)
      REAL( KIND = kreal ) gdispderiv(3, 3)
      REAL( KIND = kreal ) b1(6, ndof*nn)
      REAL( KIND = kreal ) smat(9, 9), elem(3, nn)
      REAL( KIND = kreal ) bn(9, ndof*nn), sbn(9, ndof*nn)
      REAL( KIND = kreal ) vol

!--------------------------------------------------------------------

      INTEGER( KIND = kint ) :: npoints

      INTEGER( KIND = kint ) :: na, nb
      INTEGER( KIND = kint ) :: na_p, nb_p
      INTEGER( KIND = kint ) :: isize, jsize
      INTEGER( KIND = kint ) :: jsize1, jsize2, jsize3
      INTEGER( KIND = kint ) :: ip, jp, kp

      REAL( KIND = kreal ) :: stiff_up(3*nn, nn_p)
      REAL( KIND = kreal ) :: stiff_pp(nn_p, nn_p), stiff_pp_inv(nn_p, nn_p)
      REAL( KIND = kreal ) :: stiff_up_stiff_pp_inv(3*nn, nn_p)
      REAL( KIND = kreal ) :: stiff_up_stiff_pp_inv_stiff_up(3*nn, 3*nn)
      REAL( KIND = kreal ) :: alpha_inv
      REAL( KIND = kreal ) :: g
      REAL( KIND = kreal ) :: d2(6)
      REAL( KIND = kreal ) :: bd2(3*nn)
      REAL( KIND = kreal ) :: lambda_lx
      REAL( KIND = kreal ) :: lambda_lx_2
      REAL( KIND = kreal ) :: spfunc_p(nn_p)

!--------------------------------------------------------------------

      stiff(:, :) = 0.0D0
      stiff_up(:, :) = 0.0D0
      stiff_pp(:, :) = 0.0D0

      ! we suppose the same material type in the element
      flag = gausses(1)%pMaterial%nlgeom_flag

      IF( .NOT. PRESENT( u ) ) flag = INFINITE    ! enforce to infinite deformation analysis

      elem(:, :) = ecoord(:, :)

      IF( flag .NE. 1 ) THEN

       WRITE(6, *) 'Error: this formulation should be the total Lagrangian.'

       RETURN

      END IF

!--------------------------------------------------------------------

      DO lx = 1, NumOfQuadPoints(etype)

       !--------------------------------------------------------

       CALL getQuadPoint( etype, lx, naturalcoord(:) )

       !--------------------------------------------------------

       CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

       !--------------------------------------------------------

       ! Displacement gradient
       gdispderiv(1:3, 1:3) = MATMUL( u(1:3, 1:nn), gderiv(1:nn, 1:3) )

       !--------------------------------------------------------

       IF( nn_p .EQ. 1 ) THEN

        spfunc_p(1) = 1.0D0

       END IF

       !--------------------------------------------------------

       ! Lagrange multiplier
       lambda_lx = 0.0D0

       DO na_p = 1, nn_p

        lambda_lx = lambda_lx+spfunc_p(na_p)*lambda(na_p)

       END DO

       !--------------------------------------------------------

       IF( PRESENT( temperature ) ) THEN

        CALL getShapeFunc( etype, naturalcoord, spfunc )

        temp = DOT_PRODUCT( temperature, spfunc )

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr, temp )

       ELSE

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr )

       END IF

       !--------------------------------------------------------

       wg = getWeight( etype, lx )*det

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       DO jsize = 1, 3*nn

        b(:, jsize) = b(:, jsize)+b1(:, jsize)

       END DO

       !--------------------------------------------------------

       db(1:6, 1:3*nn) = MATMUL( d(1:6, 1:6), b(1:6, 1:3*nn) )

       ! [ Kuu ] matrix
       FORALL( isize = 1:3*nn, jsize = 1:3*nn )

        stiff(isize, jsize) = stiff(isize, jsize)+wg*DOT_PRODUCT( b(:, isize), db(:, jsize) )

       END FORALL

       !--------------------------------------------------------

       ! calculate the initial stress matrix
       stress(1:6) = gausses(lx)%stress

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       smat(:, :) = 0.0D0

       DO j = 1, 3

        smat(j  , j  ) = stress(1) ! Sxx
        smat(j  , j+3) = stress(4) ! Sxy
        smat(j  , j+6) = stress(6) ! Sxz

        smat(j+3, j  ) = stress(4) ! Syx
        smat(j+3, j+3) = stress(2) ! Syy
        smat(j+3, j+6) = stress(5) ! Syz

        smat(j+6, j  ) = stress(6) ! Szx
        smat(j+6, j+3) = stress(5) ! Szy
        smat(j+6, j+6) = stress(3) ! Szz

       END DO

       !--------------------------------------------------------

       sbn(1:9, 1:3*nn) = MATMUL( smat(1:9, 1:9), bn(1:9, 1:3*nn) )

       ! [ Kuu ] matrix
       FORALL( isize = 1:3*nn, jsize = 1:3*nn )

        stiff(isize, jsize) = stiff(isize, jsize)+wg*DOT_PRODUCT( bn(:, isize), sbn(:, jsize) )

       END FORALL

       !--------------------------------------------------------

       DO isize = 1, 3*nn

        bd2(isize) = DOT_PRODUCT( b(:, isize), d2(:) )

       END DO

       ! [ Kup ] matrix
       FORALL( isize = 1:3*nn, nb_p = 1:nn_p )

        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)

       END FORALL

       !--------------------------------------------------------

       ! [ Kpp ] matrix
       FORALL( na_p = 1:nn_p, nb_p = 1:nn_p )

        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)

       END FORALL

       !--------------------------------------------------------

      END DO ! gauss roop

!--------------------------------------------------------------------

      ! [ Kpp ]^-1 matrix

      ! Constant
      IF( nn_p .EQ. 1 ) THEN

       stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)

      ELSE

       WRITE(6, *) 'Error: nn_p should be equal to 1.'

       RETURN

      END IF

!--------------------------------------------------------------------

      stiff_up_stiff_pp_inv(1:3*nn, 1:nn_p) = MATMUL( stiff_up(1:3*nn, 1:nn_p), stiff_pp_inv(1:nn_p, 1:nn_p) )

      ! [ Kup ] [ Kpp ]^-1 [ Kup ]^T matrix
      FORALL( isize = 1:3*nn, jsize = 1:3*nn )

       stiff_up_stiff_pp_inv_stiff_up(isize, jsize) = DOT_PRODUCT( stiff_up_stiff_pp_inv(isize, :), stiff_up(jsize, :) )

      END FORALL

      !---------------------------------------------------------

      ! [ Kuu ]
      FORALL( isize = 1:3*nn, jsize = 1:3*nn )

       stiff(isize, jsize) = stiff(isize, jsize)-stiff_up_stiff_pp_inv_stiff_up(isize, jsize)

      END FORALL

!--------------------------------------------------------------------

      RETURN

!####################################################################
      END SUBROUTINE STF_C3_up
!####################################################################


!####################################################################
      SUBROUTINE Update_C3_up                         &
                 ( etype, nn, ecoord, u, du, ddu, qf, &
                   gausses, iter, tincr,              &
                   nn_p, lambda, ddlambda, tt, t0 )
!####################################################################

      USE mMaterial
      USE mMechGauss
      USE m_MatMatrix
      USE m_ElastoPlastic
      USE mHyperElastic
      USE m_utilities

!--------------------------------------------------------------------

      INTEGER( KIND = kint ), INTENT(IN) :: etype           !< \param [in] element type
      INTEGER( KIND = kint ), INTENT(IN) :: nn              !< \param [in] number of elemental nodes
      REAL( KIND = kreal ), INTENT(IN) :: ecoord(3, nn)     !< \param [in] coordinates of elemental nodes
      REAL( KIND = kreal ), INTENT(IN) :: u(3, nn)          !< \param [in] nodal dislplacements
      REAL( KIND = kreal ), INTENT(IN) :: du(3, nn)         !< \param [in] nodal displacement ( solutions of solver )
      REAL( KIND = kreal ), INTENT(IN) :: ddu(3, nn)        !< \param [in] nodal displacement ( solutions of solver )
      REAL( KIND = kreal ), INTENT(OUT) :: qf(nn*3)         !< \param [out] Internal Force
      TYPE( tGaussStatus ), INTENT(INOUT) :: gausses(:)     !< \param [out] status of qudrature points

      INTEGER, INTENT(IN) :: iter
      REAL( KIND = kreal ), INTENT(IN) :: tincr             !< time increment
      REAL( KIND = kreal ), INTENT(IN), OPTIONAL :: tt(nn)  !< current temperature
      REAL( KIND = kreal ), INTENT(IN), OPTIONAL :: t0(nn)  !< reference temperature

!--------------------------------------------------------------------

      INTEGER( KIND = kint ), INTENT(IN) :: nn_p
      REAL( KIND = kreal ), INTENT(IN) :: lambda(nn_p)      !< Lagrange multiplier
      REAL( KIND = kreal ), INTENT(INOUT) :: ddlambda(nn_p)

!--------------------------------------------------------------------

      INTEGER( KIND = kint ) :: flag
      INTEGER( KIND = kint ), PARAMETER :: ndof = 3
      REAL( KIND = kreal ) :: d(6,6), b(6,ndof*nn), b1(6,ndof*nn)
      REAL( KIND = kreal ) :: gderiv(nn,3), gdispderiv(3,3), det, wg, dt
      INTEGER( KIND = kint ) :: i, j, k, lx, mtype
      REAL( KIND = kreal ) :: naturalCoord(3), rot(3,3), R(3,3), spfunc(nn)
      REAL( KIND = kreal ) :: totaldisp(3,nn), elem(3,nn), elem1(3,nn)
      REAL( KIND = kreal ) :: dstrain(6), dstress(6), dumstress(3,3), dum(3,3)
      REAL( KIND = kreal ) :: derivdum(1:ndof,1:ndof)
      REAL( KIND = kreal ) :: ttc, tt0, outa(1), ina(1), epsth(6)
      LOGICAL  :: ierr

!--------------------------------------------------------------------

      INTEGER( KIND = kint ) :: npoints

      INTEGER( KIND = kint ) :: na, nb
      INTEGER( KIND = kint ) :: na_p, nb_p
      INTEGER( KIND = kint ) :: isize, jsize
      INTEGER( KIND = kint ) :: jsize1, jsize2, jsize3

      REAL( KIND = kreal ) :: totallambda(nn_p)
      REAL( KIND = kreal ) :: stiff_up(3*nn, nn_p)
      REAL( KIND = kreal ) :: stiff_pp(nn_p, nn_p), stiff_pp_inv(nn_p, nn_p)
      REAL( KIND = kreal ) :: stiff_up_stiff_pp_inv(3*nn, nn_p)
      REAL( KIND = kreal ) :: stiff_up_stiff_pp_inv_stiff_up(3*nn, 3*nn)
      REAL( KIND = kreal ) :: stiff_up_stiff_pp_inv_qf_p(3*nn)
      REAL( KIND = kreal ) :: alpha_inv
      REAL( KIND = kreal ) :: g
      REAL( KIND = kreal ) :: qf_p(nn_p)
      REAL( KIND = kreal ) :: stiff_up_ddu(nn_p)
      REAL( KIND = kreal ) :: d2(6)
      REAL( KIND = kreal ) :: bd2(3*nn)
      REAL( KIND = kreal ) :: lambda_lx
      REAL( KIND = kreal ) :: spfunc_p(nn_p)

!--------------------------------------------------------------------

      stiff_up(:, :) = 0.0D0
      stiff_pp(:, :) = 0.0D0
      qf_p(:) = 0.0D0

      ! we suppose the same material type in the element
      flag = gausses(1)%pMaterial%nlgeom_flag

      elem(:, :) = ecoord(:, :)

      totaldisp(:, :) = u(:, :)+( du(:, :)-ddu(:, :) )

      IF( flag .NE. 1 ) THEN

       WRITE(6, *) 'Error: this formulation should be the total Lagrangian.'

       RETURN

      END IF

!--------------------------------------------------------------------

      DO lx = 1, NumOfQuadPoints(etype)

       !--------------------------------------------------------

       CALL getQuadPoint( etype, lx, naturalcoord(:) )

       !--------------------------------------------------------

       CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

       !--------------------------------------------------------

       ! Displacement gradient
       gdispderiv(1:3, 1:3) = MATMUL( totaldisp(1:3, 1:nn), gderiv(1:nn, 1:3) )
       !do isize = 1, 3
       ! do nb_p = 1, 3
       !  write(6, *) isize, nb_p, gdispderiv(isize, nb_p)
       ! end do
       !end do

       !--------------------------------------------------------

       IF( nn_p .EQ. 1 ) THEN

        spfunc_p(1) = 1.0D0

       END IF

       !--------------------------------------------------------

       ! Lagrange multiplier

       lambda_lx = 0.0D0

       DO na_p = 1, nn_p

        lambda_lx = lambda_lx+spfunc_p(na_p)*lambda(na_p)

       END DO

       !--------------------------------------------------------

       mtype = gausses(lx)%pMaterial%mtype

       !--------------------------------------------------------

       wg = getWeight( etype, lx )*det

       !--------------------------------------------------------

       ! ========================================================
       !     UPDATE STRAIN and STRESS
       ! ========================================================

       IF( isElastoplastic(gausses(lx)%pMaterial%mtype) .OR. &
       & ( gausses(lx)%pMaterial%mtype .EQ. NORTON ) ) gausses(lx)%pMaterial%mtype = ELASTIC

       epsth = 0.0D0

       IF( PRESENT( tt ) .AND. PRESENT( t0 ) ) THEN

        CALL getShapeFunc( etype, naturalcoord, spfunc )

        ttc = DOT_PRODUCT( tt, spfunc )
        tt0 = DOT_PRODUCT( t0, spfunc )

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr, ttc )

        IF( ( iter .LE. 1 ) .OR. ( flag .EQ. TOTALLAG ) ) THEN

         ina(1) = ttc

         CALL fetch_TableData( MC_THEMOEXP, gausses(lx)%pMaterial%dict, outa, ierr, ina )

         IF( ierr ) outa(1) = gausses(lx)%pMaterial%variables(M_EXAPNSION)

         epsth(1:3) = outa(1)*( ttc-tt0 )

        END IF

       ELSE

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr )

       END IF

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       DO jsize = 1, 3*nn

        b(:, jsize) = b(:, jsize)+b1(:, jsize)

       END DO

       !--------------------------------------------------------

       DO isize = 1, 3*nn

        bd2(isize) = DOT_PRODUCT( b(:, isize), d2 )

       END DO

       ! [ Kup ] matrix
       FORALL( isize = 1:3*nn, nb_p = 1:nn_p )

        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)

       END FORALL
       !do isize = 1, 3
       ! do nb_p = 1, nn_p
       !  write(6, *) isize, nb_p, stiff_up(isize, nb_p)
       ! end do
       !end do
       !--------------------------------------------------------

       ! [ Kpp ] matrix
       FORALL( na_p = 1:nn_p, nb_p = 1:nn_p )

        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)

       END FORALL

       !--------------------------------------------------------

       qf_p(1:nn_p) = qf_p(1:nn_p)+wg*spfunc_p(1:nn_p)*( g-alpha_inv*lambda_lx )

       !--------------------------------------------------------

      END DO ! gauss roop

!--------------------------------------------------------------------

      ! [ Kpp ]^-1 matrix

      ! Constant
      IF( nn_p .EQ. 1 ) THEN

       stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)

      ELSE

       WRITE(6, *) 'Error: nn_p should be equal to 1.'

       RETURN

      END IF

!--------------------------------------------------------------------

      DO na_p = 1, nn_p

       stiff_up_ddu(na_p) = 0.0D0

       DO nb = 1, nn

        DO i = 1, 3

         jsize = 3*(nb-1)+i

         stiff_up_ddu(na_p) = stiff_up_ddu(na_p)+stiff_up(jsize, na_p)*ddu(i, nb)

        END DO

       END DO

      END DO

      DO na_p = 1, nn_p

       ddlambda(na_p) = DOT_PRODUCT( stiff_pp_inv(na_p,:), -qf_p-stiff_up_ddu )

      END DO

!--------------------------------------------------------------------

      totaldisp(:, :) = u(:, :)+du(:, :)

      totallambda(:) = lambda(:)+ddlambda(:)

!--------------------------------------------------------------------

      stiff_up(:, :) = 0.0D0
      stiff_pp(:, :) = 0.0D0

      qf(:) = 0.0D0
      qf_p(:) = 0.0D0

!--------------------------------------------------------------------

      DO lx = 1, NumOfQuadPoints(etype)

       !--------------------------------------------------------

       CALL getQuadPoint( etype, lx, naturalcoord(:) )

       !--------------------------------------------------------

       CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )

       !--------------------------------------------------------

       gdispderiv(1:3, 1:3) = MATMUL( totaldisp(1:3, 1:nn), gderiv(1:nn, 1:3) )

       !--------------------------------------------------------

       IF( nn_p .EQ. 1 ) THEN

        spfunc_p(1) = 1.0D0

       END IF

       !--------------------------------------------------------

       ! Lagrange multiplier

       lambda_lx = 0.0D0

       DO na_p = 1, nn_p

        lambda_lx = lambda_lx+spfunc_p(na_p)*totallambda(na_p)

       END DO

       !--------------------------------------------------------

       mtype = gausses(lx)%pMaterial%mtype

       !--------------------------------------------------------

       wg = getWeight( etype, lx )*det

       !--------------------------------------------------------

       ! Small strain
       dstrain(1) = gdispderiv(1, 1)
       dstrain(2) = gdispderiv(2, 2)
       dstrain(3) = gdispderiv(3, 3)
       dstrain(4) = gdispderiv(1, 2)+gdispderiv(2, 1)
       dstrain(5) = gdispderiv(2, 3)+gdispderiv(3, 2)
       dstrain(6) = gdispderiv(3, 1)+gdispderiv(1, 3)
       dstrain(:) = dstrain(:)-epsth(:)

       ! Green-Lagrange strain
       dstrain(1) = dstrain(1)+0.5D0*DOT_PRODUCT( gdispderiv(:, 1), gdispderiv(:, 1) )
       dstrain(2) = dstrain(2)+0.5D0*DOT_PRODUCT( gdispderiv(:, 2), gdispderiv(:, 2) )
       dstrain(3) = dstrain(3)+0.5D0*DOT_PRODUCT( gdispderiv(:, 3), gdispderiv(:, 3) )
       dstrain(4) = dstrain(4)+DOT_PRODUCT( gdispderiv(:, 1), gdispderiv(:, 2) )
       dstrain(5) = dstrain(5)+DOT_PRODUCT( gdispderiv(:, 2), gdispderiv(:, 3) )
       dstrain(6) = dstrain(6)+DOT_PRODUCT( gdispderiv(:, 1), gdispderiv(:, 3) )

       !--------------------------------------------------------

       gausses(lx)%strain(1:6) = dstrain(1:6)+epsth(:)

       CALL StressUpdate_up( gausses(lx), d3, dstrain, gausses(lx)%stress, lambda_lx, g )

       !--------------------------------------------------------

       epsth = 0.0D0

       IF( PRESENT( tt ) .AND. PRESENT( t0 ) ) THEN

        CALL getShapeFunc( etype, naturalcoord, spfunc )

        ttc = DOT_PRODUCT( tt, spfunc )
        tt0 = DOT_PRODUCT( t0, spfunc )

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr, ttc )

        IF( ( iter .LE. 1 ) .OR. ( flag .EQ. TOTALLAG ) ) THEN

         ina(1) = ttc

         CALL fetch_TableData( MC_THEMOEXP, gausses(lx)%pMaterial%dict, outa, ierr, ina )

         IF( ierr ) outa(1) = gausses(lx)%pMaterial%variables(M_EXAPNSION)

         epsth(1:3) = outa(1)*( ttc-tt0 )

        END IF

       ELSE

        CALL MatlMatrix_up( gausses(lx), d3, d, lambda_lx, d2, alpha_inv, g, tincr )

       END IF

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       DO nb = 1, nn

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

       END DO

       !--------------------------------------------------------

       DO jsize = 1, 3*nn

        b(:, jsize) = b(:, jsize)+b1(:, jsize)

       END DO

       !--------------------------------------------------------

       DO isize = 1, 3*nn

        bd2(isize) = DOT_PRODUCT( b(:, isize), d2 )

       END DO

       ! [ Kup ] matrix
       FORALL( isize = 1:3*nn, nb_p = 1:nn_p )

        stiff_up(isize, nb_p) = stiff_up(isize, nb_p)+wg*bd2(isize)*spfunc_p(nb_p)

       END FORALL

       !--------------------------------------------------------

       ! [ Kpp ] matrix
       FORALL( na_p = 1:nn_p, nb_p = 1:nn_p )

        stiff_pp(na_p, nb_p) = stiff_pp(na_p, nb_p)-wg*alpha_inv*spfunc_p(na_p)*spfunc_p(nb_p)

       END FORALL

       !--------------------------------------------------------

       ! calculate the Internal Force
       qf(1:3*nn) = qf(1:3*nn)+wg*MATMUL( gausses(lx)%stress(1:6), b(1:6,1:3*nn) )

       qf_p(1:nn_p) = qf_p(1:nn_p)+wg*spfunc_p(1:nn_p)*( g-alpha_inv*lambda_lx )

       !--------------------------------------------

      END DO ! gauss roop

!--------------------------------------------------------------------

      ! [ Kpp ]^-1 matrix

      ! Constant
      IF( nn_p .EQ. 1 ) THEN

       stiff_pp_inv(1, 1) = 1.0D0/stiff_pp(1, 1)

      ELSE

       WRITE(6, *) 'Error: nn_p should be equal to 1.'

       RETURN

      END IF

!--------------------------------------------------------------------

      stiff_up_stiff_pp_inv(1:3*nn, 1:nn_p) = MATMUL( stiff_up(1:3*nn, 1:nn_p), stiff_pp_inv(1:nn_p, 1:nn_p) )

      DO isize = 1, 3*nn

       stiff_up_stiff_pp_inv_qf_p(isize) = DOT_PRODUCT( stiff_up_stiff_pp_inv(isize, :), qf_p )

      END DO

      !---------------------------------------------------------

      DO isize = 1, 3*nn

       qf(isize) = qf(isize)-stiff_up_stiff_pp_inv_qf_p(isize)

      END DO

!--------------------------------------------------------------------

      RETURN

!####################################################################
      END SUBROUTINE Update_C3_up
!####################################################################


!####################################################################
      END MODULE m_static_LIB_3d_up