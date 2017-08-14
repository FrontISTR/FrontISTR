!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains several strategy to free locking problem
!> in Eight-node hexagonal element
MODULE m_static_LIB_C3D8

   USE hecmw, only : kint, kreal
   USE elementInfo

   IMPLICIT NONE

   CONTAINS


!>  This subroutine calculate stiff matrix using b-bar method
!>
!> \see Hughes, T. J. "Generalization of Selective Integration Procedures
!>  to Anisotropic and Nonlinear Media", Intl. J. Numer. Methods Engng, 15,
!>  pp1413-1418,1980
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3D8Bbar                                          &
              (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
               time, tincr, u, temperature)
!----------------------------------------------------------------------*

    USE mMechGauss
    USE m_MatMatrix
    USE m_common_struct
    USE m_static_LIB_3d, ONLY: GEOMAT_C3

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN)  :: etype                  !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                     !< number of elemental nodes
    REAL(kind=kreal), INTENT(IN)    :: ecoord(3, nn)          !< coordinates of elemental nodes
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)             !< status of qudrature points
    REAL(kind=kreal), INTENT(OUT)   :: stiff(:,:)             !< stiff matrix
    INTEGER(kind=kint), INTENT(IN)  :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT) :: coords(3, 3)           !< variables to define matreial coordinate system
    REAL(kind=kreal), INTENT(IN)    :: time                   !< current time
    REAL(kind=kreal), INTENT(IN)    :: tincr                  !< time increment
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: u(:, :)         !< nodal displacemwent
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: temperature(nn) !< temperature

!---------------------------------------------------------------------

    INTEGER(kind=kint) :: flag
    INTEGER(kind=kint), PARAMETER :: ndof = 3
    REAL(kind=kreal) :: D(6, 6),B(6, ndof*nn),DB(6, ndof*nn)
    REAL(kind=kreal) :: gderiv(nn, 3),stress(6),mat(6, 6)
    REAL(kind=kreal) :: det, wg, temp, spfunc(nn)
    INTEGER(kind=kint) :: i, j, LX, serr
    REAL(kind=kreal) :: naturalCoord(3)
    REAL(kind=kreal) :: gdispderiv(3, 3)
    REAL(kind=kreal) :: B1(6, ndof*nn), Bbar(nn, 3)
    REAL(kind=kreal) :: Smat(9, 9), elem(3, nn)
    REAL(kind=kreal) :: BN(9, ndof*nn), SBN(9, ndof*nn)
    REAL(kind=kreal) :: B4, B6, B8, Vol, coordsys(3, 3)

!---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    IF( .NOT. PRESENT(u) ) flag = INFINITE    ! enforce to infinite deformation analysis
    elem(:, :) = ecoord(:, :)
    IF( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+u(:, :)

    ! dilatation component at centroid
    naturalCoord = 0.0D0
    CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, Bbar)

    DO LX = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      IF( PRESENT(temperature) ) THEN
        CALL getShapeFunc( etype, naturalcoord, spfunc )
        temp = DOT_PRODUCT( temperature, spfunc )
        CALL MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, temp )
      ELSE
        CALL MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys )
      END IF

      IF( flag == UPDATELAG ) then
        CALL GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      ENDIF

      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      DO j = 1, nn
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
      END DO

      ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      IF( flag == TOTALLAG ) THEN
        ! ---dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof, 1:ndof) = MATMUL( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
        B1(1:6, 1:nn*ndof) = 0.0D0
        DO j = 1, nn
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
        END DO
        ! ---BL = BL0 + BL1
        DO j = 1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        END DO
      END IF

      DB(1:6, 1:nn*ndof) = MATMUL( D, B(1:6, 1:nn*ndof) )
      FORALL( i=1:nn*ndof, j=1:nn*ndof )
        stiff(i, j) = stiff(i, j)+DOT_PRODUCT( B(:, i), DB(:, j) )*wg
      END FORALL

      ! calculate the initial stress matrix
      IF( flag == TOTALLAG .or. flag == UPDATELAG ) THEN
        stress(1:6) = gausses(LX)%stress
        BN(1:9, 1:nn*ndof) = 0.0D0
        DO j = 1, nn
          BN(1, 3*j-2) = gderiv(j, 1)
          BN(2, 3*j-1) = gderiv(j, 1)
          BN(3, 3*j  ) = gderiv(j, 1)
          BN(4, 3*j-2) = gderiv(j, 2)
          BN(5, 3*j-1) = gderiv(j, 2)
          BN(6, 3*j  ) = gderiv(j, 2)
          BN(7, 3*j-2) = gderiv(j, 3)
          BN(8, 3*j-1) = gderiv(j, 3)
          BN(9, 3*j  ) = gderiv(j, 3)
        END DO
        Smat(:, :) = 0.0D0
        DO j= 1, 3
          Smat(j  , j  ) = stress(1)
          Smat(j  , j+3) = stress(4)
          Smat(j  , j+6) = stress(6)
          Smat(j+3, j  ) = stress(4)
          Smat(j+3, j+3) = stress(2)
          Smat(j+3, j+6) = stress(5)
          Smat(j+6, j  ) = stress(6)
          Smat(j+6, j+3) = stress(5)
          Smat(j+6, j+6) = stress(3)
        END DO
        SBN(1:9, 1:nn*ndof) = MATMUL( Smat(1:9, 1:9), BN(1:9, 1:nn*ndof) )
        FORALL( i=1:nn*ndof, j=1:nn*ndof )
          stiff(i, j) = stiff(i, j)+DOT_PRODUCT( BN(:, i), SBN(:, j) )*wg
        END FORALL
      END IF

    END DO ! gauss roop

   END SUBROUTINE STF_C3D8Bbar


!>  Update Strain stress of this element
!----------------------------------------------------------------------*
   SUBROUTINE Update_C3D8Bbar                              &
              (etype, nn, ecoord, u, du, cdsys_ID, coords, &
               qf ,gausses, iter, time, tincr, TT,T0, TN  )
!----------------------------------------------------------------------*

    USE m_fstr
    USE mMaterial
    USE mMechGauss
    USE m_MatMatrix
    USE m_ElastoPlastic
    USE mHyperElastic
    USE m_utilities

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN)    :: etype         !< \param [in] element type
    INTEGER(kind=kint), INTENT(IN)    :: nn            !< \param [in] number of elemental nodes
    REAL(kind=kreal), INTENT(IN)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    REAL(kind=kreal), INTENT(IN)      :: u(3, nn)      !< \param [in] nodal dislplacements
    REAL(kind=kreal), INTENT(IN)      :: du(3, nn)     !< \param [in] nodal displacement increment
    INTEGER(kind=kint), INTENT(IN)    :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT)   :: coords(3, 3)  !< variables to define matreial coordinate system
    REAL(kind=kreal), INTENT(OUT)     :: qf(nn*3)      !< \param [out] Internal Force
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)    !< \param [out] status of qudrature points
    INTEGER, INTENT(IN) :: iter
    REAL(kind=kreal), INTENT(IN)      :: time          !< current time
    REAL(kind=kreal), INTENT(IN)      :: tincr         !< time increment
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: TT(nn)   !< current temperature
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: T0(nn)   !< reference temperature
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: TN(nn)   !< reference temperature

!---------------------------------------------------------------------

    INTEGER(kind=kint) :: flag
    INTEGER(kind=kint), PARAMETER :: ndof = 3
    REAL(kind=kreal) :: D(6, 6), B(6, ndof*nn), B1(6, ndof*nn)
    REAL(kind=kreal) :: gderiv(nn, 3), gdispderiv(3, 3), det, wg
    INTEGER(kind=kint) :: i, j, k, LX, mtype, serr
    integer(kind=kint) :: isEp
    REAL(kind=kreal) :: naturalCoord(3), rot(3, 3), R(3, 3), spfunc(nn)
    REAL(kind=kreal) :: totaldisp(3, nn), elem(3, nn), elem1(3, nn), coordsys(3, 3), tm(6, 6)
    REAL(kind=kreal) :: dstrain(6), dstress(6), dumstress(3, 3), dum(3, 3)
    REAL(kind=kreal) :: dvol, vol0, Bbar(nn, 3), derivdum(1:ndof, 1:ndof), BBar2(nn, 3)
    REAL(kind=kreal) :: B4, B6, B8, ttc, tt0, ttn, alp, alp0, alpo(3), alpo0(3), outa(1), ina(1), EPSTH(6)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:, :) = ecoord(:, :)
    totaldisp(:, :) = u(:, :)+du(:, :)
    IF( flag == UPDATELAG ) THEN
      elem(:, :) = ( 0.5D0*du(:, :)+u(:, :) )+ecoord(:, :)
      elem1(:, :) = ( du(:, :)+u(:, :) )+ecoord(:, :)
      !  elem = elem1
      totaldisp(:, :) = du(:, :)
    END IF

    matlaniso = .FALSE.
    IF( cdsys_ID > 0 .AND. PRESENT(TT) ) THEN
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .TRUE.
    END IF

    ! dilatation at centroid
    naturalCoord = 0.0D0
    CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, Bbar)
    derivdum = MATMUL( totaldisp(1:ndof, 1:nn), Bbar(1:nn, 1:ndof) )
    vol0 = ( derivdum(1, 1)+derivdum(2, 2)+derivdum(3, 3) )/3.0D0
    IF( flag == UPDATELAG ) CALL getGlobalDeriv(etype, nn, naturalcoord, elem1, det, Bbar2)

    DO LX = 1, NumOfQuadPoints(etype)

      mtype = gausses(LX)%pMaterial%mtype

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, elem, det, gderiv)

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      gdispderiv(1:ndof, 1:ndof) = MATMUL( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      dvol = vol0-( gdispderiv(1, 1)+gdispderiv(2, 2)+gdispderiv(3, 3) )/3.0D0
      !gdispderiv(1, 1) = gdispderiv(1, 1)+dvol
      !gdispderiv(2, 2) = gdispderiv(2, 2)+dvol
      !gdispderiv(3, 3) = gdispderiv(3, 3)+dvol

      ! ========================================================
      !     UPDATE STRAIN and STRESS
      ! ========================================================

      if( isElastoplastic(mtype) .OR. mtype == NORTON )then
        isEp = 1
      else
        isEp = 0
      endif
      !gausses(LX)%pMaterial%mtype = ELASTIC

      EPSTH = 0.0D0
      IF( PRESENT(tt) .AND. PRESENT(t0) ) THEN
        CALL getShapeFunc(etype, naturalcoord, spfunc)
        ttc = DOT_PRODUCT(TT, spfunc)
        tt0 = DOT_PRODUCT(T0, spfunc)
        ttn = DOT_PRODUCT(TN, spfunc)
        CALL MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, ttc, isEp )

        ina(1) = ttc
        IF( matlaniso ) THEN
          CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
          IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
        ELSE
          CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
          IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
          alp = outa(1)
        END IF
        ina(1) = tt0
        IF( matlaniso ) THEN
          CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
          IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
        ELSE
          CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
          IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
          alp0 = outa(1)
        END IF
        IF( matlaniso ) THEN
          DO j=1,3
            EPSTH(j) = ALPO(j)*( ttc-ref_temp )-alpo0(j)*( tt0-ref_temp )
          END DO
          CALL transformation( coordsys(:,:), tm )
          EPSTH(:) = MATMUL( EPSTH(:), tm  ) ! to global coord
          EPSTH(4:6) = EPSTH(4:6)*2.0D0
        ELSE
          EPSTH(1:3) = ALP*( ttc-ref_temp )-alp0*( tt0-ref_temp )
        END IF

      ELSE

        CALL MatlMatrix( gausses(LX), D3, D, time, tincr, coordsys, isEp=isEp )

      END IF

      ! Small strain
      dstrain(1) = gdispderiv(1, 1)+dvol
      dstrain(2) = gdispderiv(2, 2)+dvol
      dstrain(3) = gdispderiv(3, 3)+dvol
      dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
      dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
      dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
      dstrain(:) = dstrain(:)-EPSTH(:)

      IF( flag == INFINITE ) THEN

        gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
        gausses(LX)%stress(1:6) = matmul( D(1:6, 1:6), dstrain(1:6) )
        IF( isViscoelastic(mtype) .AND. tincr /= 0.0D0 ) THEN
          IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, time, tincr, ttc, ttn )
          ELSE
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, time, tincr )
          END IF
          gausses(LX)%stress = REAL( gausses(LX)%stress )
        END IF

      ELSE IF( flag == TOTALLAG ) THEN

        ! Green-Lagrange strain
        dstrain(1) = dstrain(1)+0.5D0*DOT_PRODUCT( gdispderiv(:, 1), gdispderiv(:, 1) )
        dstrain(2) = dstrain(2)+0.5D0*DOT_PRODUCT( gdispderiv(:, 2), gdispderiv(:, 2) )
        dstrain(3) = dstrain(3)+0.5D0*DOT_PRODUCT( gdispderiv(:, 3), gdispderiv(:, 3) )
        dstrain(4) = dstrain(4)+( gdispderiv(1, 1)*gdispderiv(1, 2)                                     &
                                 +gdispderiv(2, 1)*gdispderiv(2, 2)+gdispderiv(3, 1)*gdispderiv(3, 2) )
        dstrain(5) = dstrain(5)+( gdispderiv(1, 2)*gdispderiv(1, 3)                                     &
                                 +gdispderiv(2, 2)*gdispderiv(2, 3)+gdispderiv(3, 2)*gdispderiv(3, 3) )
        dstrain(6) = dstrain(6)+( gdispderiv(1, 1)*gdispderiv(1, 3)                                     &
                                 +gdispderiv(2, 1)*gdispderiv(2, 3)+gdispderiv(3, 1)*gdispderiv(3, 3) )

        IF( mtype == NEOHOOKE .OR. mtype == MOONEYRIVLIN .OR.  mtype == ARRUDABOYCE  .OR.         &
            mtype == USERELASTIC .OR. mtype == USERHYPERELASTIC .OR. mtype == USERMATERIAL ) THEN
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
        ELSE IF( ( isViscoelastic(mtype) .OR. mtype == NORTON ) .AND. tincr /= 0.0D0  ) THEN
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%stress(1:6) = MATMUL( D(1:6, 1:6), dstrain(1:6) )
          !gausses(LX)%pMaterial%mtype = mtype
          IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, time, tincr, ttc, ttn )
          ELSE
            CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, time, tincr )
          END IF
        ELSE
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%stress(1:6) = matmul( D(1:6, 1:6), dstrain(1:6) )
        END IF

      ELSE IF( flag == UPDATELAG ) THEN

        rot = 0.0D0
        rot(1, 2)= 0.5D0*( gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
        rot(2, 3)= 0.5D0*( gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
        rot(1, 3)= 0.5D0*( gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

        gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)

        dstress = ( MATMUL( D(1:6, 1:6), dstrain(1:6) ) )
        dumstress(1, 1) = gausses(LX)%stress_bak(1)
        dumstress(2, 2) = gausses(LX)%stress_bak(2)
        dumstress(3, 3) = gausses(LX)%stress_bak(3)
        dumstress(1, 2) = gausses(LX)%stress_bak(4);  dumstress(2, 1) = dumstress(1, 2)
        dumstress(2, 3) = gausses(LX)%stress_bak(5);  dumstress(3, 2) = dumstress(2, 3)
        dumstress(3, 1) = gausses(LX)%stress_bak(6);  dumstress(1, 3) = dumstress(3, 1)

        dum(:, :) = MATMUL(rot, dumstress)-MATMUL(dumstress, rot)

        gausses(LX)%stress(1) = gausses(LX)%stress_bak(1)+dstress(1)+dum(1,1)-gausses(LX)%stress_bak(1)*3.0D0*vol0
        gausses(LX)%stress(2) = gausses(LX)%stress_bak(2)+dstress(2)+dum(2,2)-gausses(LX)%stress_bak(2)*3.0D0*vol0
        gausses(LX)%stress(3) = gausses(LX)%stress_bak(3)+dstress(3)+dum(3,3)-gausses(LX)%stress_bak(3)*3.0D0*vol0
        gausses(LX)%stress(4) = gausses(LX)%stress_bak(4)+dstress(4)+dum(1,2)-gausses(LX)%stress_bak(4)*3.0D0*vol0
        gausses(LX)%stress(5) = gausses(LX)%stress_bak(5)+dstress(5)+dum(2,3)-gausses(LX)%stress_bak(5)*3.0D0*vol0
        gausses(LX)%stress(6) = gausses(LX)%stress_bak(6)+dstress(6)+dum(3,1)-gausses(LX)%stress_bak(6)*3.0D0*vol0

        IF( mtype == USERMATERIAL ) THEN
          CALL StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
        ELSEIF( mtype == NORTON ) THEN
          !gausses(LX)%pMaterial%mtype = mtype
          IF( tincr /= 0.0D0 .AND. any(gausses(LX)%stress /= 0.0D0) ) THEN
            IF( PRESENT(TT) .AND. PRESENT(T0) ) THEN
              CALL StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, time, tincr, ttc, ttn )
            ELSE
              CALL StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, time, tincr )
            END IF
          END IF
        END IF

      END IF

      IF( isElastoplastic(mtype) ) THEN
        !gausses(LX)%pMaterial%mtype = mtype
        IF( PRESENT(tt) ) THEN
          CALL BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain, &
                              gausses(LX)%istatus(1), gausses(LX)%fstatus, ttc )
        ELSE
          CALL BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain, &
                              gausses(LX)%istatus(1), gausses(LX)%fstatus )
        END IF
      END IF

      ! ========================================================
      ! calculate the internal force ( equivalent nodal force )
      ! ========================================================
      ! Small strain
      B(1:6, 1:nn*ndof) = 0.0D0
      DO j = 1, nn
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
      END DO

      IF( flag == INFINITE ) THEN

      ELSE IF( flag == TOTALLAG ) THEN

        B1(1:6, 1:nn*ndof) = 0.0D0
        DO j = 1, nn
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
        END DO
        ! BL = BL0 + BL1
        DO j = 1, nn*ndof
          B(:, j) = B(:, j)+B1(:, j)
        END DO

      ELSE IF( flag == UPDATELAG ) THEN

        CALL getGlobalDeriv(etype, nn, naturalcoord, elem1, det, gderiv)
        B(1:6, 1:nn*ndof) = 0.0D0
        DO j = 1, nn
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
        END DO
      END IF

      !!  calculate the Internal Force
      wg = getWeight(etype, LX)*det
      qf(1:nn*ndof) = qf(1:nn*ndof)                                           &
                     +MATMUL( gausses(LX)%stress(1:6), B(1:6, 1:nn*ndof) )*wg

    END DO

   END SUBROUTINE Update_C3D8Bbar

   !> This subroutien calculate thermal loading
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C3D8Bbar                    &
              (etype, nn, XX, YY, ZZ, TT, T0,   &
               gausses, VECT, cdsys_ID, coords)
!----------------------------------------------------------------------*

    use m_fstr
    USE mMechGauss
    use m_MatMatrix
    use m_utilities

!---------------------------------------------------------------------

    INTEGER(kind=kint), PARAMETER   :: ndof = 3
    INTEGER(kind=kint), INTENT(IN)  :: etype, nn
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)             !< status of qudrature points
    REAL(kind=kreal), INTENT(IN)    :: XX(nn), YY(nn), ZZ(nn)
    REAL(kind=kreal), INTENT(IN)    :: TT(nn), T0(nn)
    REAL(kind=kreal), INTENT(OUT)   :: VECT(nn*ndof)
    INTEGER(kind=kint), INTENT(IN)  :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT) :: coords(3, 3)           !< variables to define matreial coordinate system

!---------------------------------------------------------------------

    REAL(kind=kreal) :: ALP, alp0, D(6, 6), B(6, ndof*nn)
    REAL(kind=kreal) :: B4, B6, B8, det, ecoord(3, nn)
    INTEGER(kind=kint) :: j, LX, serr
    REAL(kind=kreal) :: estrain(6), SGM(6), H(nn)
    REAL(kind=kreal) :: naturalcoord(3), gderiv(nn, 3)
    REAL(kind=kreal) :: wg, outa(1), ina(1), Bbar(nn, 3), alpo(3), alpo0(3), coordsys(3, 3)
    REAL(kind=kreal) :: TEMPC, TEMP0, TEMPC0, TEMP00, THERMAL_EPS, tm(6,6)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    matlaniso = .FALSE.

    IF( cdsys_ID > 0 ) THEN   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .TRUE.
    END IF

    VECT(:) = 0.0D0

    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)

    naturalCoord = 0.0D0
    CALL getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, Bbar)
    CALL getShapeFunc( etype, naturalcoord, H(1:nn) )
    TEMPC0 = DOT_PRODUCT( H(1:nn), TT(1:nn) )
    TEMP00 = DOT_PRODUCT( H(1:nn), T0(1:nn) )

    ! LOOP FOR INTEGRATION POINTS
    DO LX = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getShapeFunc( etype, naturalcoord, H(1:nn) )
      CALL getGlobalDeriv(etype, nn, naturalcoord, ecoord, det, gderiv)

      IF( matlaniso ) THEN
        CALL set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      !  WEIGHT VALUE AT GAUSSIAN POINT
      wg = getWeight(etype, LX)*det
      B(1:6, 1:nn*ndof) = 0.0D0
      DO j = 1, nn
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
      END DO

      TEMPC = DOT_PRODUCT( H(1:nn),TT(1:nn) )
      TEMP0 = DOT_PRODUCT( H(1:nn),T0(1:nn) )

      CALL MatlMatrix( gausses(LX), D3, D, 1.d0, 0.0D0, coordsys, TEMPC )

      ina(1) = TEMPC
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      END IF
      ina(1) = TEMP0
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      END IF

      !**
      !** THERMAL strain
      !**
      IF( matlaniso ) THEN
        DO j = 1, 3
          estrain(j) = ALPO(j)*(TEMPC0-ref_temp)-alpo0(j)*(TEMP00-ref_temp)
        END DO
        estrain(4:6) = 0.0D0
        CALL transformation(coordsys, tm)
        estrain(:) = matmul( estrain(:), tm  )      ! to global coord
        estrain(4:6) = estrain(4:6)*2.0D0
      ELSE
        THERMAL_EPS = ALP*(TEMPC0-ref_temp)-alp0*(TEMP00-ref_temp)
        estrain(1:3) = THERMAL_EPS
        estrain(4:6) = 0.0D0
      END IF

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = MATMUL( D(:, :), estrain(:) )

      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:nn*ndof) = VECT(1:nn*ndof)+MATMUL( SGM(:), B(:, :) )*wg

    END DO

  END SUBROUTINE TLOAD_C3D8Bbar


END MODULE m_static_LIB_C3D8
