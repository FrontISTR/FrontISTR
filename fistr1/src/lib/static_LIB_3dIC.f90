!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.7                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN (AdavanceSoft)                         !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  Eight-node hexagonal element with imcompatible mode
!>
!> \see R.L.Taylor,P.J.Bereford, and E.L.Wilson, "A Nonconforming element
!>  for Stress Analysis", Intl. J. Numer. Methods Engng, 10(6), pp1211-1219
!>  ,1976
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/08/12
!>  \version    0.00
!                                                                    !
!======================================================================!
MODULE m_static_LIB_3dIC

   USE hecmw, only : kint, kreal
   USE m_utilities
   USE elementInfo

   IMPLICIT NONE

   CONTAINS

!>  CALCULATION STIFF Matrix for C3D8IC ELEMENT
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3D8IC                                            &
              (etype, nn, ecoord, gausses, stiff, cdsys_ID, coords, &
               tincr, nddisp, ehdisp, temperature)
!----------------------------------------------------------------------*

    USE mMechGauss
    USE m_MatMatrix
    USE m_common_struct
    USE m_static_LIB_3d, ONLY: GEOMAT_C3

!---------------------------------------------------------------------

    INTEGER(kind=kint), INTENT(IN)  :: etype                !< element type, not used here
    INTEGER(kind=kint), INTENT(IN)  :: nn                   !< number of elements nodes
    REAL(kind=kreal), INTENT(IN)    :: ecoord(3, nn)        !< nodal coord of curr element
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)           !< Info of qudrature points
    REAL(kind=kreal), INTENT(OUT)   :: stiff(:, :)          !< stiffness matrix
    INTEGER(kind=kint), INTENT(IN)  :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT) :: coords(3, 3)         !< variables to define matreial coordinate system
    REAL(kind=kreal), INTENT(IN)    :: tincr                !< time increment
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: nddisp(3, nn) !< nodal displacemwent
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: ehdisp(3, 3)  !< enhanced disp of bending mode
    REAL(kind=kreal), INTENT(IN), OPTIONAL :: temperature(nn) !< temperature

!---------------------------------------------------------------------

    INTEGER(kind=kint) :: flag
    INTEGER(kind=kint), PARAMETER :: ndof = 3
    REAL(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    REAL(kind=kreal) :: gderiv(nn+3, 3), stress(6)
    REAL(kind=kreal) :: xj(9, 9), jacobian(3, 3), inverse(3, 3)
    REAL(kind=kreal) :: tmpstiff((nn+3)*3, (nn+3)*3), tmpk(nn*3, 9)
    REAL(kind=kreal) :: det, wg, elem(3, nn), mat(6, 6)
    INTEGER(kind=kint) :: i, j, LX, fetype
    INTEGER(kind=kint) :: serr
    REAL(kind=kreal) :: naturalCoord(3), unode(3, nn+3)
    REAL(kind=kreal) :: gdispderiv(3, 3), coordsys(3, 3)
    REAL(kind=kreal) :: B1(6, ndof*(nn+3))
    REAL(kind=kreal) :: Smat(9, 9)
    REAL(kind=kreal) :: BN(9, ndof*(nn+3)), SBN(9, ndof*(nn+3))
    REAL(kind=kreal) :: spfunc(nn), temp

!---------------------------------------------------------------------

    fetype = fe_hex8n
    IF( PRESENT(nddisp) .AND. PRESENT(ehdisp) ) THEN
      unode(:, 1:nn)      = nddisp(:, :)
      unode(:, nn+1:nn+3) = ehdisp(:, :)
    END IF

    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    IF( .NOT. PRESENT(nddisp) ) flag = INFINITE    ! enforce to infinite deformation analysis
    elem(:, :) = ecoord(:, :)
    IF( flag == UPDATELAG ) elem(:, :) = ecoord(:, :)+unode(:, 1:nn)

    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    CALL getJacobian(fetype, nn, naturalcoord, elem, det, jacobian, inverse)
    inverse(:, :)= inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    tmpstiff(:, :) = 0.0D0
    B1(1:6, 1:(nn+3)*ndof) = 0.0D0
    BN(1:9, 1:(nn+3)*ndof) = 0.0D0

    DO LX = 1, NumOfQuadPoints(fetype)

      CALL getQuadPoint(fetype, LX, naturalCoord)
      CALL getGlobalDeriv( fetype, nn, naturalcoord, elem, det, gderiv(1:nn, 1:3) )

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
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, temp )
      ELSE
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys )
      END IF

      IF( flag == UPDATELAG ) then
        CALL GEOMAT_C3( gausses(LX)%stress, mat )
        D(:, :) = D(:, :)-mat
      ENDIF

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
      DO j = 1, nn+3
        B(1, 3*j-2) = gderiv(j, 1)
        B(2, 3*j-1) = gderiv(j, 2)
        B(3, 3*j  ) = gderiv(j, 3)
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
        gdispderiv(1:ndof,1:ndof) = MATMUL( unode(1:ndof, 1:nn+3), gderiv(1:nn+3, 1:ndof) )
        DO j = 1, nn+3

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

        END DO
        ! ---BL = BL0 + BL1
        DO j = 1, (nn+3)*ndof
          B(:, j) = B(:, j)+B1(:, j)
        END DO

      END IF

      DB(1:6, 1:(nn+3)*ndof) = MATMUL( D, B(1:6, 1:(nn+3)*ndof) )
      FORALL( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
        tmpstiff(i, j) = tmpstiff(i, j)+dot_product( B(:, i), DB(:, j) )*wg
      END FORALL
    END DO

    ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
    IF( flag == TOTALLAG .OR. flag == UPDATELAG ) THEN
      stress(1:6) = gausses(LX)%stress
      DO j = 1, nn+3
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
      DO j = 1, 3
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
      SBN(1:9, 1:(nn+3)*ndof) = MATMUL( Smat(1:9, 1:9), BN(1:9, 1:(nn+3)*ndof) )
      FORALL( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
        tmpstiff(i, j) = tmpstiff(i, j)+DOT_PRODUCT( BN(:, i), SBN(:, j) )*wg
      END FORALL
    END IF

    ! -----Condense tmpstiff to stiff
    xj(1:9, 1:9)= tmpstiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    CALL calInverse(9, xj)
    tmpk = MATMUL( tmpstiff( 1:nn*ndof, nn*ndof+1:(nn+3)*ndof ), xj )
    stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)-MATMUL( tmpk, tmpstiff(nn*ndof+1:(nn+3)*ndof, 1:nn*ndof)  )

   END SUBROUTINE STF_C3D8IC


!>  Update Strain stress of this element
!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C3D8IC                        &
              (etype, nn, xx, yy, zz, tt, t0, edisp, &
               gausses, cdsys_ID, coords)
!----------------------------------------------------------------------*

    USE m_fstr
    USE mMechGauss
    USE m_MatMatrix
    USE m_common_struct

!---------------------------------------------------------------------

    INTEGER(kind=kint), PARAMETER     :: ndof=3
    INTEGER(kind=kint), INTENT(IN)    :: etype                   !< element type, not used here
    INTEGER(kind=kint), INTENT(IN)    :: nn                      !< number of element nodes
    REAL(kind=kreal), INTENT(IN)      :: xx(nn), yy(nn), zz(nn)  !< nodes coordinate of element
    REAL(kind=kreal), INTENT(IN)      :: tt(nn), t0(nn)          !< current and ref temprature
    REAL(kind=kreal), INTENT(IN)      :: edisp(nn*ndof)          !< nodal displacement
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)              !< info about qudrature points
    INTEGER(kind=kint), INTENT(IN)    :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT)   :: coords(3, 3)            !< variables to define matreial coordinate system

!---------------------------------------------------------------------

    REAL(kind=kreal) :: ALP, ALP0
    REAL(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    REAL(kind=kreal) :: det, wg, ecoord(3, nn)
    INTEGER(kind=kint) :: j, k, IC, serr, fetype
    REAL(kind=kreal) :: EPSA(6), EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3), xj(9, 9)
    REAL(kind=kreal) :: naturalcoord(3), gderiv(nn+3, 3)
    REAL(kind=kreal) :: jacobian(3, 3),inverse(3, 3)
    REAL(kind=kreal) :: stiff((nn+3)*3, (nn+3)*3)
    REAL(kind=kreal) :: tmpforce(9), cdisp((nn+3)*3)
    REAL(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6), outa(1), ina(1)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    matlaniso = .FALSE.
    IF( cdsys_ID > 0 ) THEN   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .true.
    END IF

    fetype = fe_hex8n
    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)
    ! ---- Calculate enhanced displacement at first
    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    CALL getJacobian(fetype, nn, naturalcoord, ecoord, det, jacobian, inverse)
    inverse(:, :) = inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff(:, :) = 0.0D0
    B(1:6, 1:(nn+3)*ndof) = 0.0D0

    DO IC = 1, NumOfQuadPoints(fetype)

      CALL getQuadPoint(fetype, IC, naturalCoord)
      CALL getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn, 1:3) )

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      CALL getShapeFunc( fetype, naturalcoord, H(1:nn) )
      TEMPC = DOT_PRODUCT( H(1:nn), TT(1:nn) )
      CALL MatlMatrix( gausses(IC), D3, D, 1.0D0, coordsys, TEMPC )

      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(fetype, IC)*det
      DO j = 1, nn+3
        B(1, 3*j-2) = gderiv(j, 1)
        B(2, 3*j-1) = gderiv(j, 2)
        B(3, 3*j  ) = gderiv(j, 3)
        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      END DO

      DB(1:6, 1:(nn+3)*ndof) = MATMUL( D, B(1:6, 1:(nn+3)*ndof) )
      stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof)                                     &
      = stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof)                                   &
       +MATMUL( TRANSPOSE(B(1:6, 1:(nn+3)*ndof)), DB(1:6, 1:(nn+3)*ndof) )*wg
    END DO
    xj(1:9,1:9)= stiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    CALL calInverse(9, xj)
    ! ---  [Kda]*edisp
    tmpforce(:) = MATMUL( stiff(nn*3+1:(nn+3)*ndof, 1:nn*ndof), edisp )
    cdisp(1:nn*3) = edisp(:)
    ! ---  -[Kaa]-1 * [Kda] * edisp
    cdisp(nn*3+1:(nn+3)*3) = -MATMUL( xj(:, :), tmpforce(:) )

    ! ---- Now strain and stress calculation
    DO IC = 1, NumOfQuadPoints(etype)

      CALL getQuadPoint( etype, IC, naturalCoord(:) )
      CALL getShapeFunc( etype, naturalcoord, H(1:nn) )
      CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      IF( matlaniso ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys, serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      B(1:6, 1:(nn+3)*ndof) = 0.0D0
      DO j = 1, nn+3
        B(1, 3*j-2) = gderiv(j, 1)
        B(2, 3*j-1) = gderiv(j, 2)
        B(3, 3*j  ) = gderiv(j, 3)
        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      END DO

      TEMPC = DOT_PRODUCT( H(1:nn), TT(1:nn) )
      TEMP0 = DOT_PRODUCT( H(1:nn), T0(1:nn) )

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
           WRITE(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF
      CALL MatlMatrix( gausses(IC), D3, D, 0.0D0, coordsys, TEMPC )

      ina(1) = TEMPC
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      END IF
      ina(1) = TEMP0
      IF( matlaniso  ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      END IF

      !**
      !** THERMAL strain
      !**
      IF( matlaniso ) THEN
        DO j = 1, 3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        END DO
        EPSTH(4:6) = 0.0D0
        CALL transformation(coordsys, tm)
        EPSTH(:) = MATMUL( EPSTH(:), tm  ) ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      ELSE
        THERMAL_EPS = ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      END IF
      !**
      !** SET EPS  {e}=[B]{u}
      !**
      EPSA(1:6) = MATMUL( B(1:6,:), CDISP(:) )
      !**
      !** SET SGM  {S}=[D]{e}
      !**
      DO j = 1, 6
        SGM(j) = 0.0D0
        DO k = 1, 6
          SGM(j) = SGM(j)+D(j, k)*( EPSA(k)-EPSTH(k) )
        END DO
      END DO
      !**
      !** Adding stress in each gauss points
      !**
      gausses(IC)%strain(1:6) = EPSA(1:6)
      gausses(IC)%stress(1:6) = SGM(1:6)

    END DO

   END SUBROUTINE UpdateST_C3D8IC


!>  This subroutine calculatess thermal loading
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C3D8IC                      &
              (etype, nn, xx, yy, zz, tt, t0,   &
               gausses, VECT, cdsys_ID, coords)
!----------------------------------------------------------------------*

    USE m_fstr
    USE mMechGauss
    USE m_MatMatrix
    USE m_common_struct

!---------------------------------------------------------------------

    INTEGER(kind=kint), PARAMETER     :: ndof=3
    INTEGER(kind=kint), INTENT(IN)    :: etype                   !< element type, not used here
    INTEGER(kind=kint), INTENT(IN)    :: nn                      !< number of element nodes
    REAL(kind=kreal), INTENT(IN)      :: xx(nn), yy(nn), zz(nn)  !< nodes coordinate of element
    REAL(kind=kreal), INTENT(IN)      :: tt(nn), t0(nn)          !< current and ref temprature
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)              !< info about qudrature points
    REAL(kind=kreal), INTENT(OUT)     :: VECT(nn*ndof)           !< load vector
    INTEGER(kind=kint), INTENT(IN)    :: cdsys_ID
    REAL(kind=kreal), INTENT(INOUT)   :: coords(3, 3)            !< variables to define matreial coordinate system

!---------------------------------------------------------------------

    REAL(kind=kreal) :: ALP, ALP0
    REAL(kind=kreal) :: D(6, 6), B(6, ndof*(nn+3)), DB(6, ndof*(nn+3))
    REAL(kind=kreal) :: det, wg, ecoord(3, nn)
    INTEGER(kind=kint) :: j, IC, serr, fetype
    REAL(kind=kreal) :: EPSTH(6), SGM(6), H(nn), alpo(3), alpo0(3), coordsys(3, 3), xj(9, 9)
    REAL(kind=kreal) :: naturalcoord(3), gderiv(nn+3, 3)
    REAL(kind=kreal) :: jacobian(3, 3),inverse(3, 3)
    REAL(kind=kreal) :: stiff((nn+3)*3, (nn+3)*3)
    REAL(kind=kreal) :: tmpvect((nn+3)*3)
    REAL(kind=kreal) :: tmpforce(nn*ndof), icdisp(9)
    REAL(kind=kreal) :: TEMPC, TEMP0, THERMAL_EPS, tm(6, 6), outa(1), ina(1)
    LOGICAL :: ierr, matlaniso

!---------------------------------------------------------------------

    matlaniso = .FALSE.
    IF( cdsys_ID > 0 ) THEN   ! cannot define aniso exapansion when no local coord defined
      ina = TT(1)
      CALL fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
      IF( .NOT. ierr ) matlaniso = .true.
    END IF

    fetype = fe_hex8n
    ecoord(1, :) = XX(:)
    ecoord(2, :) = YY(:)
    ecoord(3, :) = ZZ(:)
    ! ---- Calculate enhanced displacement at first
    ! --- Inverse of Jacobian at elemental center
    naturalcoord(:) = 0.0D0
    CALL getJacobian(fetype, nn, naturalcoord, ecoord, det, jacobian, inverse)
    inverse(:, :) = inverse(:, :)*det
    ! ---- We now calculate stiff matrix include imcompatible mode
    !    [   Kdd   Kda ]
    !    [   Kad   Kaa ]
    stiff(:, :) = 0.0D0
    B(1:6, 1:(nn+3)*ndof) = 0.0D0
    tmpvect(1:(nn+3)*ndof) = 0.0D0

    DO IC = 1, NumOfQuadPoints(fetype)

      CALL getQuadPoint(fetype, IC, naturalCoord)
      CALL getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn, 1:3) )

      IF( cdsys_ID > 0 ) THEN
        CALL set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
        IF( serr == -1 ) STOP "Fail to setup local coordinate"
        IF( serr == -2 ) THEN
          WRITE(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        END IF
      END IF

      CALL getShapeFunc( fetype, naturalcoord, H(1:nn) )
      TEMPC = DOT_PRODUCT( H(1:nn), TT(1:nn) )
      TEMP0 = DOT_PRODUCT( H(1:nn), T0(1:nn) )
      CALL MatlMatrix( gausses(IC), D3, D, 1.0D0, coordsys, TEMPC )

      ! -- Derivative of shape function of imcompatible mode --
      !     [ -2*a   0,   0   ]
      !     [   0,  -2*b, 0   ]
      !     [   0,   0,  -2*c ]
      ! we don't call getShapeDeriv but use above value directly for computation efficiency
      gderiv(nn+1, :) = -2.0D0*naturalcoord(1)*inverse(1, :)/det
      gderiv(nn+2, :) = -2.0D0*naturalcoord(2)*inverse(2, :)/det
      gderiv(nn+3, :) = -2.0D0*naturalcoord(3)*inverse(3, :)/det

      wg = getWeight(fetype, IC)*det
      DO j = 1, nn+3
        B(1, 3*j-2) = gderiv(j, 1)
        B(2, 3*j-1) = gderiv(j, 2)
        B(3, 3*j  ) = gderiv(j, 3)
        B(4, 3*j-2) = gderiv(j, 2)
        B(4, 3*j-1) = gderiv(j, 1)
        B(5, 3*j-1) = gderiv(j, 3)
        B(5, 3*j  ) = gderiv(j, 2)
        B(6, 3*j-2) = gderiv(j, 3)
        B(6, 3*j  ) = gderiv(j, 1)
      END DO

      DB(1:6, 1:(nn+3)*ndof) = MATMUL( D, B(1:6, 1:(nn+3)*ndof) )
      stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof)                                     &
      = stiff(1:(nn+3)*ndof, 1:(nn+3)*ndof)                                   &
       +MATMUL( TRANSPOSE(B(1:6, 1:(nn+3)*ndof)), DB(1:6, 1:(nn+3)*ndof) )*wg

      ina(1) = TEMPC
      IF( matlaniso ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
      END IF
      ina(1) = TEMP0
      IF( matlaniso  ) THEN
        CALL fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
        IF( ierr ) STOP "Fails in fetching orthotropic expansion coefficient!"
      ELSE
        CALL fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
        IF( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
      END IF

      !**
      !** THERMAL strain
      !**
      IF( matlaniso ) THEN
        DO j = 1, 3
          EPSTH(j) = ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
        END DO
        EPSTH(4:6) = 0.0D0
        CALL transformation(coordsys, tm)
        EPSTH(:) = MATMUL( EPSTH(:), tm  ) ! to global coord
        EPSTH(4:6) = EPSTH(4:6)*2.0D0
      ELSE
        THERMAL_EPS = ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
        EPSTH(1:3) = THERMAL_EPS
        EPSTH(4:6) = 0.0D0
      END IF

      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM(:) = MATMUL( D(:, :), EPSTH(:) )

      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      TMPVECT(1:(nn+3)*ndof) = TMPVECT(1:(nn+3)*ndof)+MATMUL( SGM(1:6), B(1:6, 1:(nn+3)*ndof) )*wg

    END DO

    ! --- condense tmpvect to vect
    xj(1:9,1:9)= stiff(nn*ndof+1:(nn+3)*ndof, nn*ndof+1:(nn+3)*ndof)
    CALL calInverse(9, xj)
    ! ---  [Kaa]-1 * fa
    icdisp(1:9) = MATMUL( xj(:, :), tmpvect(nn*ndof+1:(nn+3)*ndof) )
    ! ---  [Kda] * [Kaa]-1 * fa
    tmpforce(1:nn*ndof) = MATMUL( stiff(1:nn*ndof, nn*3+1:(nn+3)*ndof), icdisp(1:9) )
    ! ---  fd - [Kda] * [Kaa]-1 * fa
    vect(1:nn*ndof) = tmpvect(1:nn*ndof) - tmpforce(1:nn*ndof)

   END SUBROUTINE TLOAD_C3D8IC


END MODULE m_static_LIB_3dIC
