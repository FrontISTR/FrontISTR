! Fluid (2016/09/08) <
MODULE m_static_LIB_3d_vp
    
    USE hecmw, ONLY : kint, kreal
    USE elementInfo
    
    IMPLICIT NONE
    
    CONTAINS
    
    
!--------------------------------------------------------------------
    SUBROUTINE STF_C3_vp                                  &
               (etype, nn, ecoord, gausses, stiff, tincr, &
                v, temperature)                           
!--------------------------------------------------------------------
    
    USE mMechGauss
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint), INTENT(IN) :: etype                     !< element type
    INTEGER(kind=kint), INTENT(IN) :: nn                        !< the number of elemental nodes
    REAL(kind=kreal),   INTENT(IN) :: ecoord(3, nn)             !< coordinates of elemental nodes
    TYPE(tGaussStatus), INTENT(IN) :: gausses(:)                !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)               !< stiff matrix
    REAL(kind=kreal),   INTENT(IN) :: tincr                     !< time increment
    REAL(kind=kreal),   INTENT(IN), OPTIONAL :: v(:, :)         !< nodal velocity
    REAL(kind=kreal),   INTENT(IN), OPTIONAL :: temperature(nn) !< temperature
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint) :: i, j
    INTEGER(kind=kint) :: na, nb
    INTEGER(kind=kint) :: isize, jsize
    INTEGER(kind=kint) :: LX
    
    REAL(kind=kreal) :: MM(nn, nn), AA(nn, nn), DD(nn, nn, 3, 3), &
                        trD(nn, nn), BB(nn, nn), CC(nn, nn, 3),   &
                        MS(nn, nn), AS(nn, nn), CS(nn, nn, 3),    &
                        MP(nn, nn, 3), AP(nn, nn, 3), CP(nn, nn)  
    REAL(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    REAL(kind=kreal) :: elem(3, nn)
    REAL(kind=kreal) :: naturalCoord(3)
    REAL(kind=kreal) :: dndx(nn, 3)
    REAL(kind=kreal) :: tincr_inv
    REAL(kind=kreal) :: volume, volume_inv
    REAL(kind=kreal) :: mu
    REAL(kind=kreal) :: rho, rho_inv
    REAL(kind=kreal) :: vx, vy, vz
    REAL(kind=kreal) :: t1, t2, t3
    REAL(kind=kreal) :: v_dot_v
    REAL(kind=kreal) :: d
    REAL(kind=kreal) :: det, wg
    REAL(kind=kreal) :: tau
    REAL(kind=kreal), PARAMETER :: gamma = 0.5D0
    
!--------------------------------------------------------------------
    
    tincr_inv = 1.0D0/tincr
    
!--------------------------------------------------------------------
    
    elem(:, :) = ecoord(:, :)
    
!--------------------------------------------------------------------
    
    t1 = 2.0D0*tincr_inv
    
    !---------------------------------------------------------------
    
    volume = 0.0D0
    
    loopVolume: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(:) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     volume = volume+wg
     
     !----------------------------------------------------------
     
    END DO loopVolume
    
    volume_inv = 1.0D0/volume
    
    !---------------------------------------------------------------
    
    naturalCoord(1) = 0.25D0
    naturalCoord(2) = 0.25D0
    naturalCoord(3) = 0.25D0
    
    CALL getShapeFunc(etype, naturalCoord, spfunc)
    
    vx = 0.0D0
    vy = 0.0D0
    vz = 0.0D0
    
    DO na = 1, nn
     
     vx = vx+spfunc(na)*v(1, na)
     vy = vy+spfunc(na)*v(2, na)
     vz = vz+spfunc(na)*v(3, na)
     
    END DO
    
    v_dot_v = vx*vx+vy*vy+vz*vz 
    
    !---------------------------------------------------------------
    
    mu  = 0.0D0
    rho = 0.0D0
    
    DO na = 1, nn
     
     dndx(na, 1) = 0.0D0
     dndx(na, 2) = 0.0D0
     dndx(na, 3) = 0.0D0
     
    END DO
    
    loopGlobalDeriv: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(1:3) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     mu = mu+wg*gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     rho = rho+wg*gausses(LX)%pMaterial%variables(M_DENSITY)
     
     !----------------------------------------------------------
     
     DO na = 1, nn
      
      dndx(na, 1) = dndx(na, 1)+wg*gderiv(na, 1)
      dndx(na, 2) = dndx(na, 2)+wg*gderiv(na, 2)
      dndx(na, 3) = dndx(na, 3)+wg*gderiv(na, 3)
      
     END DO
     
     !----------------------------------------------------------
     
    END DO loopGlobalDeriv
    
    mu  = volume_inv*mu
    rho = volume_inv*rho
    
    DO na = 1, nn
     
     dndx(na, 1) = volume_inv*dndx(na, 1)
     dndx(na, 2) = volume_inv*dndx(na, 2)
     dndx(na, 3) = volume_inv*dndx(na, 3)
     
    END DO
    
    !---------------------------------------------------------------
    
    d = 0.0D0
    
    DO na = 1, nn
     
     d = d+DABS( vx*dndx(na, 1)+vy*dndx(na, 2)+vz*dndx(na, 3) )
     
    END DO
    
    ! h_es3d = 2.0D0/( d/DSQRT( v_dot_v ) )
    
    !---------------------------------------------------------------
    
    ! t2 = 2.0D0*DSQRT( v_dot_v )/h_es3d
    t2 = d
    
    !----------------------------------------------------------------
    
    IF( v_dot_v .LT. 1.0D-15 ) THEN
     
     t3 = 4.0D0*mu/( rho*volume**(2.0D0/3.0D0) )
     
    ELSE
     
     t3 = mu*d*d/( rho*v_dot_v )
     
    END IF
    
    !----------------------------------------------------------------
    
    tau = 1.0D0/DSQRT( t1*t1+t2*t2+t3*t3 )
    
!--------------------------------------------------------------------
    
    stiff(:, :) = 0.0D0
    
    loopGauss: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     mu = gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     rho = gausses(LX)%pMaterial%variables(M_DENSITY)
     rho_inv = 1.0D0/rho
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(1:3) )
     CALL getShapeFunc( etype, naturalCoord(:), spfunc(:) )
     CALL getGlobalDeriv( etype, nn, naturalCoord(:), elem(:,:), &
                          det, gderiv(:,:) )                     
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     vx = 0.0D0
     vy = 0.0D0
     vz = 0.0D0
     
     DO na = 1, nn
      
      vx = vx+spfunc(na)*v(1, na)
      vy = vy+spfunc(na)*v(2, na)
      vz = vz+spfunc(na)*v(3, na)
      
     END DO
     
     !----------------------------------------------------------
     
     FORALL(na = 1:nn, nb = 1:nn)
      
      MM(na, nb) = spfunc(na)*spfunc(nb)
      AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
                  +vy*( spfunc(na)*gderiv(nb, 2) ) &
                  +vz*( spfunc(na)*gderiv(nb, 3) ) 
      DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
      DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
      DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
      DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
      DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
      DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
      DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
      DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
      DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
      trD(na, nb) = DD(na, nb, 1, 1) &
                   +DD(na, nb, 2, 2) &
                   +DD(na, nb, 3, 3) 
      BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
                  +( vx*vy )*DD(na, nb, 1, 2) &
                  +( vx*vz )*DD(na, nb, 1, 3) &
                  +( vy*vx )*DD(na, nb, 2, 1) &
                  +( vy*vy )*DD(na, nb, 2, 2) &
                  +( vy*vz )*DD(na, nb, 2, 3) &
                  +( vz*vx )*DD(na, nb, 3, 1) &
                  +( vz*vy )*DD(na, nb, 3, 2) &
                  +( vz*vz )*DD(na, nb, 3, 3) 
      CC(na, nb, 1) = gderiv(na, 1)*spfunc(nb)
      CC(na, nb, 2) = gderiv(na, 2)*spfunc(nb)
      CC(na, nb, 3) = gderiv(na, 3)*spfunc(nb)
      
      MS(nb, na) = AA(na, nb)
      AS(na, nb) = BB(na, nb)
      CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
                     +vy*DD(na, nb, 2, 1) &
                     +vz*DD(na, nb, 3, 1) 
      CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
                     +vy*DD(na, nb, 2, 2) &
                     +vz*DD(na, nb, 3, 2) 
      CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
                     +vy*DD(na, nb, 2, 3) &
                     +vz*DD(na, nb, 3, 3) 
      MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
      MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
      MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
      AP(nb, na, 1) = CS(na, nb, 1)
      AP(nb, na, 2) = CS(na, nb, 2)
      AP(nb, na, 3) = CS(na, nb, 3) 
      CP(na, nb) = trD(na, nb)
      
     END FORALL
     
     !----------------------------------------------------------
     
     DO nb = 1, nn
      
      DO na = 1, nn
       
       i = 1
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                              &
       = stiff(isize, jsize)                            &
        +wg                                             &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )        
       
       i = 1
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 2, 1) ) 
       
       i = 1
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 3, 1) ) 
       
       i = 1
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 1)+tau*CS(na, nb, 1) ) 
       
       i = 2
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 1, 2) ) 
       
       i = 2
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                              &
       = stiff(isize, jsize)                            &
        +wg                                             &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) ) 
       
       i = 2
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 3, 2) ) 
       
       i = 2
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 2)+tau*CS(na, nb, 2) ) 
        
       i = 3
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 1, 3) ) 
       
       i = 3
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 2, 3) ) 
       
       i = 3
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                               &
       = stiff(isize, jsize)                             &
        +wg                                              &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) )  &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )      &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) ) 
      
       i = 3
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 3)+tau*CS(na, nb, 3) ) 
       
       i = 4
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)            &
       = stiff(isize, jsize)          &
        +wg                           &
         *( rho_inv*tau*trD(na, nb) ) 
       
      END DO
      
     END DO
     
     !----------------------------------------------------------
     
    END DO loopGauss
    
!--------------------------------------------------------------------    
    END SUBROUTINE STF_C3_vp
!--------------------------------------------------------------------
    
    
!--------------------------------------------------------------------
    SUBROUTINE UPDATE_C3_vp                        &
               (etype, nn, ecoord, v, dv, gausses) 
!--------------------------------------------------------------------
    
    USE mMechGauss
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint), INTENT(IN) :: etype            !< element type
    INTEGER(kind=kint), INTENT(IN) :: nn               !< the number of elemental nodes
    REAL(kind=kreal), INTENT(IN) :: ecoord(3, nn)      !< coordinates of elemental nodes
    REAL(kind=kreal), INTENT(IN) :: v(4, nn)           !< nodal velcoity
    REAL(kind=kreal), INTENT(IN) :: dv(4, nn)          !< nodal velocity increment
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)    !< status of qudrature points
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint) :: LX
    
    REAL(kind=kreal) :: elem(3, nn)
    REAL(kind=kreal) :: totalvelo(4, nn)
    REAL(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    REAL(kind=kreal) :: gveloderiv(3, 3)
    REAL(kind=kreal) :: naturalCoord(3)
    REAL(kind=kreal) :: det
    REAL(kind=kreal) :: mu
    REAL(kind=kreal) :: p
    
!--------------------------------------------------------------------
    
    elem(:, :) = ecoord(:, :)
    
    totalvelo(:, :) = v(:, :)+dv(:, :)
    
!--------------------------------------------------------------------
    
    loopMatrix: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     mu = gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(:) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     ! Deformation rate tensor
     gveloderiv(1:3, 1:3) = MATMUL( totalvelo(1:3, 1:nn), gderiv(1:nn, 1:3) )
     gausses(LX)%strain(1) = gveloderiv(1, 1)
     gausses(LX)%strain(2) = gveloderiv(2, 2)
     gausses(LX)%strain(3) = gveloderiv(3, 3)
     gausses(LX)%strain(4) = 0.5D0*( gveloderiv(1, 2)+gveloderiv(2, 1) )
     gausses(LX)%strain(5) = 0.5D0*( gveloderiv(2, 3)+gveloderiv(3, 2) )
     gausses(LX)%strain(6) = 0.5D0*( gveloderiv(3, 1)+gveloderiv(1, 3) )
     
     !----------------------------------------------------------
     
     ! Pressure
     p = DOT_PRODUCT(totalvelo(4, 1:nn), spfunc(1:nn))
     
     ! Cauchy stress tensor
     gausses(LX)%stress(1) = -p+2.0D0*mu*gausses(LX)%strain(1)
     gausses(LX)%stress(2) = -p+2.0D0*mu*gausses(LX)%strain(2)
     gausses(LX)%stress(3) = -p+2.0D0*mu*gausses(LX)%strain(3)
     gausses(LX)%stress(4) = 2.0D0*mu*gausses(LX)%strain(4)
     gausses(LX)%stress(5) = 2.0D0*mu*gausses(LX)%strain(5)
     gausses(LX)%stress(6) = 2.0D0*mu*gausses(LX)%strain(6)
     
     !----------------------------------------------------------
     
    END DO loopMatrix
    
    !----------------------------------------------------------------
    
!--------------------------------------------------------------------
    END SUBROUTINE UPDATE_C3_vp
!--------------------------------------------------------------------
    
    
!--------------------------------------------------------------------
    SUBROUTINE LOAD_C3_vp                                    &
               (etype, nn, ecoord, v, dv, r, gausses, tincr) 
!--------------------------------------------------------------------
    
    USE mMechGauss
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint), INTENT(IN)    :: etype         !< element type
    INTEGER(kind=kint), INTENT(IN)    :: nn            !< the number of elemental nodes
    REAL(kind=kreal), INTENT(IN)      :: ecoord(3, nn) !< coordinates of elemental nodes
    REAL(kind=kreal), INTENT(IN)      :: v(4, nn)      !< nodal dislplacements
    REAL(kind=kreal), INTENT(IN)      :: dv(4, nn)     !< nodal velocity increment
    REAL(kind=kreal), INTENT(OUT)     :: r(4*nn)       !< elemental residual
    TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)    !< status of qudrature points
    REAL(kind=kreal), intent(in)      :: tincr         !< time increment
    
!--------------------------------------------------------------------
    
    INTEGER(kind=kint) :: i, j, k
    INTEGER(kind=kint) :: na, nb
    INTEGER(kind=kint) :: isize, jsize
    INTEGER(kind=kint) :: isize1, isize2, isize3, isize4
    INTEGER(kind=kint) :: LX
    
    REAL(kind=kreal) :: elem(3, nn)
    REAL(kind=kreal) :: totalvelo(4, nn)
    REAL(kind=kreal) :: velo_new(4*nn)
    REAL(kind=kreal) :: stiff(4*nn, 4*nn)
    REAL(kind=kreal) :: b(4*nn)
    REAL(kind=kreal) :: MM(nn, nn), AA(nn, nn), DD(nn, nn, 3, 3), &
                        trD(nn, nn), BB(nn, nn), CC(nn, nn, 3),   &
                        MS(nn, nn), AS(nn, nn), CS(nn, nn, 3),    &
                        MP(nn, nn, 3), AP(nn, nn, 3), CP(nn, nn)  
    REAL(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    REAL(kind=kreal) :: gveloderiv(3, 3)
    REAL(kind=kreal) :: naturalCoord(3)
    REAL(kind=kreal) :: dndx(nn, 3)
    REAL(kind=kreal) :: tincr_inv
    REAL(kind=kreal) :: volume, volume_inv
    REAL(kind=kreal) :: mu
    REAL(kind=kreal) :: rho, rho_inv
    REAL(kind=kreal) :: vx, vy, vz
    REAL(kind=kreal) :: t1, t2, t3
    REAL(kind=kreal) :: v_a_dot_v_a
    REAL(kind=kreal) :: d
    REAL(kind=kreal) :: det, wg
    REAL(kind=kreal) :: tau
    REAL(kind=kreal) :: m_v(3), a_v(3), d_v(3, 3, 3),        &
                        ms_v(3), as_v(3), mp_dot_v, ap_dot_v 
    REAL(kind=kreal) :: stiff_velo
    REAL(kind=kreal), PARAMETER :: gamma = 0.5D0
    
!--------------------------------------------------------------------
    
    tincr_inv = 1.0D0/tincr
    
!--------------------------------------------------------------------
    
    elem(:, :) = ecoord(:, :)
    
    FORALL(na = 1:nn, i = 1:4)
     
     velo_new( 4*(na-1)+i ) = v(i, na)+dv(i, na)
     
    END FORALL
    
!--------------------------------------------------------------------
    
    t1 = 2.0D0*tincr_inv
    
    !---------------------------------------------------------------
    
    volume = 0.0D0
    
    loopVolume: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(:) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     volume = volume+wg
     
     !----------------------------------------------------------
     
    END DO loopVolume
    
    volume_inv = 1.0D0/volume
    
    !---------------------------------------------------------------
    
    naturalCoord(1) = 0.25D0
    naturalCoord(2) = 0.25D0
    naturalCoord(3) = 0.25D0
    
    CALL getShapeFunc(etype, naturalCoord, spfunc)
    
    vx = 0.0D0
    vy = 0.0D0
    vz = 0.0D0
    
    DO na = 1, nn
     
     vx = vx+spfunc(na)*v(1, na)
     vy = vy+spfunc(na)*v(2, na)
     vz = vz+spfunc(na)*v(3, na)
     
    END DO
    
    v_a_dot_v_a = vx*vx+vy*vy+vz*vz 
    
    !---------------------------------------------------------------
    
    mu = 0.0D0
    
    DO na = 1, nn
     
     dndx(na, 1) = 0.0D0
     dndx(na, 2) = 0.0D0
     dndx(na, 3) = 0.0D0
     
    END DO
    
    loopGlobalDeriv: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(1:3) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     mu = mu+wg*gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     !----------------------------------------------------------
     
     DO na = 1, nn
      
      dndx(na, 1) = dndx(na, 1)+wg*gderiv(na, 1)
      dndx(na, 2) = dndx(na, 2)+wg*gderiv(na, 2)
      dndx(na, 3) = dndx(na, 3)+wg*gderiv(na, 3)
      
     END DO
     
     !----------------------------------------------------------
     
    END DO loopGlobalDeriv
    
    mu = volume_inv*mu
    
    DO na = 1, nn
     
     dndx(na, 1) = volume_inv*dndx(na, 1)
     dndx(na, 2) = volume_inv*dndx(na, 2)
     dndx(na, 3) = volume_inv*dndx(na, 3)
     
    END DO
    
    !---------------------------------------------------------------
    
    d = 0.0D0
    
    DO na = 1, nn
     
     d = d+DABS( vx*dndx(na, 1)+vy*dndx(na, 2)+vz*dndx(na, 3) )
     
    END DO
    
    ! h_es3d = 2.0D0/( d/DSQRT( v_dot_v ) )
    
    !---------------------------------------------------------------
    
    ! t2 = 2.0D0*DSQRT( v_dot_v )/h_es3d
    t2 = d
    
    !----------------------------------------------------------------
    
    IF( v_a_dot_v_a .LT. 1.0D-15 ) THEN
     
     t3 = 4.0D0*mu/( rho*volume**(2.0D0/3.0D0) )
     
    ELSE
     
     t3 = mu*d*d/( rho*v_a_dot_v_a )
     
    END IF
    
    !----------------------------------------------------------------
    
    tau = 1.0D0/DSQRT( t1*t1+t2*t2+t3*t3 )
    
!--------------------------------------------------------------------
    
    stiff(:, :) = 0.0D0
    
    loopMatrix: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     mu = gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     rho = gausses(LX)%pMaterial%variables(M_DENSITY)
     rho_inv = 1.0D0/rho
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(:) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     vx = 0.0D0
     vy = 0.0D0
     vz = 0.0D0
     
     DO na = 1, nn
      
      vx = vx+spfunc(na)*v(1, na)
      vy = vy+spfunc(na)*v(2, na)
      vz = vz+spfunc(na)*v(3, na)
      
     END DO
     
     !----------------------------------------------------------
     
     FORALL(na = 1:nn, nb = 1:nn)
      
      MM(na, nb) = spfunc(na)*spfunc(nb)
      AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
                  +vy*( spfunc(na)*gderiv(nb, 2) ) &
                  +vz*( spfunc(na)*gderiv(nb, 3) ) 
      DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
      DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
      DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
      DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
      DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
      DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
      DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
      DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
      DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
      trD(na, nb) = DD(na, nb, 1, 1) &
                   +DD(na, nb, 2, 2) &
                   +DD(na, nb, 3, 3) 
      BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
                  +( vx*vy )*DD(na, nb, 1, 2) &
                  +( vx*vz )*DD(na, nb, 1, 3) &
                  +( vy*vx )*DD(na, nb, 2, 1) &
                  +( vy*vy )*DD(na, nb, 2, 2) &
                  +( vy*vz )*DD(na, nb, 2, 3) &
                  +( vz*vx )*DD(na, nb, 3, 1) &
                  +( vz*vy )*DD(na, nb, 3, 2) &
                  +( vz*vz )*DD(na, nb, 3, 3) 
      CC(na, nb, 1) = gderiv(na, 1)*spfunc(nb)
      CC(na, nb, 2) = gderiv(na, 2)*spfunc(nb)
      CC(na, nb, 3) = gderiv(na, 3)*spfunc(nb)
      
      MS(nb, na) = AA(na, nb)
      AS(na, nb) = BB(na, nb)
      CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
                     +vy*DD(na, nb, 2, 1) &
                     +vz*DD(na, nb, 3, 1) 
      CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
                     +vy*DD(na, nb, 2, 2) &
                     +vz*DD(na, nb, 3, 2) 
      CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
                     +vy*DD(na, nb, 2, 3) &
                     +vz*DD(na, nb, 3, 3) 
      MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
      MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
      MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
      AP(nb, na, 1) = CS(na, nb, 1)
      AP(nb, na, 2) = CS(na, nb, 2)
      AP(nb, na, 3) = CS(na, nb, 3) 
      CP(na, nb) = trD(na, nb)
      
     END FORALL
     
     !----------------------------------------------------------
     
     DO nb = 1, nn
      
      DO na = 1, nn
       
       i = 1
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                              &
       = stiff(isize, jsize)                            &
        +wg                                             &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) ) 
       
       i = 1
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 2, 1) ) 
       
       i = 1
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 3, 1) ) 
       
       i = 1
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 1)+tau*CS(na, nb, 1) ) 
       
       i = 2
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 1, 2) ) 
       
       i = 2
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                              &
       = stiff(isize, jsize)                            &
        +wg                                             &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) ) 
       
       i = 2
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 3, 2) ) 
       
       i = 2
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 2)+tau*CS(na, nb, 2) ) 
       
       i = 3
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 1, 3) ) 
       
       i = 3
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( gamma*mu*DD(na, nb, 2, 3) ) 
       
       i = 3
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                              &
       = stiff(isize, jsize)                            &
        +wg                                             &
         *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
           +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
           +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) ) 
       
       i = 3
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)                     &
       = stiff(isize, jsize)                   &
        +wg                                    &
         *( -CC(na, nb, 3)+tau*CS(na, nb, 3) ) 
       
       i = 4
       j = 1
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 2
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 3
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)              &
       = stiff(isize, jsize)            &
        +wg                             &
         *( CC(nb, na, j)               &
           +tincr_inv*tau*MP(na, nb, j) &
           +gamma*tau*AP(na, nb, j) )   
       
       i = 4
       j = 4
       isize = 4*(na-1)+i
       jsize = 4*(nb-1)+j
       
       stiff(isize, jsize)            &
       = stiff(isize, jsize)          &
        +wg                           &
         *( rho_inv*tau*trD(na, nb) ) 
       
      END DO
      
     END DO
     
     !----------------------------------------------------------
     
    END DO loopMatrix
    
!--------------------------------------------------------------------
    
    b(:) = 0.0D0
    
    loopVector: DO LX = 1, NumOfQuadPoints(etype)
     
     !----------------------------------------------------------
     
     mu = gausses(LX)%pMaterial%variables(M_YOUNGS)
     
     rho = gausses(LX)%pMaterial%variables(M_DENSITY)
     rho_inv = 1.0D0/rho
     
     !----------------------------------------------------------
     
     CALL getQuadPoint( etype, LX, naturalCoord(:) )
     CALL getShapeFunc(etype, naturalCoord, spfunc)
     CALL getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)
     
     !----------------------------------------------------------
     
     wg = getWeight(etype, LX)*det
     
     !----------------------------------------------------------
     
     vx = 0.0D0
     vy = 0.0D0
     vz = 0.0D0
     
     DO na = 1, nn
      
      vx = vx+spfunc(na)*v(1, na)
      vy = vy+spfunc(na)*v(2, na)
      vz = vz+spfunc(na)*v(3, na)
      
     END DO
     
     !----------------------------------------------------------
     
     FORALL(na = 1:nn, nb = 1:nn)
      
      MM(na, nb) = spfunc(na)*spfunc(nb)
      AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
                  +vy*( spfunc(na)*gderiv(nb, 2) ) &
                  +vz*( spfunc(na)*gderiv(nb, 3) ) 
      DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
      DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
      DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
      DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
      DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
      DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
      DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
      DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
      DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
      BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
                  +( vx*vy )*DD(na, nb, 1, 2) &
                  +( vx*vz )*DD(na, nb, 1, 3) &
                  +( vy*vx )*DD(na, nb, 2, 1) &
                  +( vy*vy )*DD(na, nb, 2, 2) &
                  +( vy*vz )*DD(na, nb, 2, 3) &
                  +( vz*vx )*DD(na, nb, 3, 1) &
                  +( vz*vy )*DD(na, nb, 3, 2) &
                  +( vz*vz )*DD(na, nb, 3, 3) 
      
      MS(nb, na) = AA(na, nb)
      AS(na, nb) = BB(na, nb)
      CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
                     +vy*DD(na, nb, 2, 1) &
                     +vz*DD(na, nb, 3, 1) 
      CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
                     +vy*DD(na, nb, 2, 2) &
                     +vz*DD(na, nb, 3, 2) 
      CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
                     +vy*DD(na, nb, 2, 3) &
                     +vz*DD(na, nb, 3, 3) 
      MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
      MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
      MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
      AP(nb, na, 1) = CS(na, nb, 1)
      AP(nb, na, 2) = CS(na, nb, 2)
      AP(nb, na, 3) = CS(na, nb, 3) 
      
     END FORALL
     
     !----------------------------------------------------------
     
     DO na = 1, nn
      
      !----------------------------------------------------
      
      DO i = 1, 3
       
       m_v(i) = 0.0D0
       a_v(i) = 0.0D0
       DO j = 1, 3
        DO k = 1, 3
         d_v(j, k, i) = 0.0D0
        END DO
       END DO
       ms_v(i) = 0.0D0
       as_v(i) = 0.0D0
       mp_dot_v = 0.0D0
       ap_dot_v = 0.0D0
       
       DO nb = 1, nn
        
        ! Unsteady term
        m_v(i) = m_v(i)+MM(na, nb)*v(i, nb)
        ! Advection term
        a_v(i) = a_v(i)+AA(na, nb)*v(i, nb)
        ! Diffusion term
        DO j = 1, 3
         DO k = 1, 3
          d_v(j, k, i) = d_v(j, k, i)+DD(na, nb, j, k)*v(i, nb)
         END DO
        END DO
        ! Unsteady term (SUPG) 
        ms_v(i) = ms_v(i)+MS(na, nb)*v(i, nb)
        ! Advection term (SUPG)
        as_v(i) = as_v(i)+AS(na, nb)*v(i, nb)
        ! Unsteady term (PSPG) 
        mp_dot_v = mp_dot_v+( MP(na, nb, 1)*v(1, nb)   &
                             +MP(na, nb, 2)*v(2, nb)   &
                             +MP(na, nb, 3)*v(3, nb) ) 
        ! Advection term (PSPG) 
        ap_dot_v = ap_dot_v+( AP(na, nb, 1)*v(1, nb)   &
                             +AP(na, nb, 2)*v(2, nb)   &
                             +AP(na, nb, 3)*v(3, nb) ) 
       END DO
       
      END DO
      
      !----------------------------------------------------
      
      DO i = 1, 3
       
       isize = 4*(na-1)+i
       
       b(isize)                                                  &
       = b(isize)                                                &
        +wg                                                      &
         *( tincr_inv*rho*( m_v(i)+tau*ms_v(i) )                 &
           -( 1.0D0-gamma )*rho*( a_v(i)+tau*as_v(i) )           &
           -( 1.0D0-gamma )                                      &
            *mu*( ( d_v(1, 1, i)+d_v(2, 2, i)+d_v(3, 3, i) )     &
                 +( d_v(1, i, 1)+d_v(2, i, 2)+d_v(3, i, 3) ) ) ) 
       
      END DO
      
      i = 4
      isize = 4*(na-1)+i
      
      b(isize)                                &
      = b(isize)                              &
       +wg                                    &
        *( tincr_inv*tau*( mp_dot_v )         &
          -( 1.0D0-gamma )*tau*( ap_dot_v ) ) 
      
      !----------------------------------------------------
      
     END DO
     
     !----------------------------------------------------------
     
    END DO loopVector
    
    !----------------------------------------------------------------
    
    DO isize = 1, 4*nn
     
     stiff_velo = 0.0D0
     
     DO jsize = 1, 4*nn
      
      stiff_velo = stiff_velo+stiff(isize, jsize)*velo_new(jsize)
      
     END DO
     
     r(isize) = b(isize)-stiff_velo
     
    END DO
    
!--------------------------------------------------------------------
    END SUBROUTINE LOAD_C3_vp
!--------------------------------------------------------------------
    
    
END MODULE m_static_LIB_3d_vp
! > Fluid (2016/09/08) 