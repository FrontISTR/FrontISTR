!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Gaku Hashimoto (Univ. of Tokyo)        !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

      MODULE m_static_LIB_shell
!####################################################################
      
      USE hecmw, ONLY : kint, kreal
      USE elementInfo
      
!--------------------------------------------------------------------
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      !--------------------------------------------------------------
      !
      ! (Programmer) 
      ! Gaku Hashimoto
      ! Department of Human and Engineered Environmental Studies
      ! Graduate School of Frontier Sciences, The University of Tokyo
      ! 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-8563 JAPAN
      !
      ! (Ref.)
      ! [1] Noguchi, H. and Hisada, T., 
      !     "Sensitivity analysis in post-buckling problems of shell 
      !     structures," 
      !     Computers & Structures, Vol.47, No.4, pp.699-710, (1993).
      ! [2] Dvorkin, E.N. and Bathe, K.J., 
      !     "A Continuum Mechanics Based Four-node Shell Element for 
      !     General Non-linear Analysis,"
      !     Engineering Computations, Vol.1, pp.77-88, (1984).
      ! [3] Bucalem, M.L. and Bathe, K.J., 
      !     "Higher-order MITC general shell element," 
      !     International Journal for Numerical Methods in 
      !     Engineering, Vol.36, pp.3729-3754, (1993).
      ! [4] Lee, P.S. and Bathe, K.J., 
      !     "Development of MITC Isotropic Triangular Shell Finite 
      !     Elements," 
      !     Computers & Structures, Vol.82, pp.945-962, (2004).
      !
      !--------------------------------------------------------------
      
!--------------------------------------------------------------------
      
      CONTAINS
      
      
!####################################################################
      SUBROUTINE STF_Shell_MITC                                           &
                 (etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp) 
!####################################################################
      
      USE mMechGauss
      USE gauss_integration
      USE m_MatMatrix

!-------------------------------------------------------------------- 
      
      INTEGER(KIND = kint), INTENT(IN) :: etype
      INTEGER(KIND = kint), INTENT(IN) :: nn, mixflag
      INTEGER(KIND = kint), INTENT(IN) :: ndof
      REAL(KIND = kreal), INTENT(IN)   :: ecoord(3, nn)
      TYPE(tGaussStatus), INTENT(IN)   :: gausses(:)
      REAL(KIND = kreal), INTENT(OUT)  :: stiff(:, :)
      REAL(KIND = kreal), INTENT(IN)   :: thick
      
      REAL(KIND = kreal), INTENT(IN), OPTIONAL :: nddisp(3, nn)
	  
!--------------------------------------------------------------------
      
      INTEGER :: flag, flag_dof, shellmatl
      INTEGER :: i, j, m
      INTEGER :: lx, ly
      INTEGER :: fetype
      INTEGER :: ny
      INTEGER :: ntying
      INTEGER :: npoints_tying(3)
      INTEGER :: it, ip
      INTEGER :: na, nb
      INTEGER :: isize, jsize
      INTEGER :: jsize1, jsize2, jsize3, &
                 jsize4, jsize5, jsize6  
      INTEGER :: n_layer,n_total_layer, sstable(24)
      
      REAL(KIND = kreal) :: D(5, 5), B(5, ndof*nn), DB(5, ndof*nn)
      REAL(KIND = kreal) :: tmpstiff(ndof*nn, ndof*nn)
      REAL(KIND = kreal) :: elem(3, nn)
      REAL(KIND = kreal) :: unode(3, nn)
      REAL(KIND = kreal) :: xi_lx, eta_lx, zeta_ly
      REAL(KIND = kreal) :: w_w_lx, w_ly
      REAL(KIND = kreal) :: B_di(5, ndof*nn, 6, 3, 7)
      REAL(KIND = kreal) :: B1(3, ndof*nn), B2(3, ndof*nn), &
                            B3(3, ndof*nn)                  
      REAL(KIND = kreal) :: naturalcoord(2)
      REAL(KIND = kreal) :: tpcoord(6, 2, 3)
      REAL(KIND = kreal) :: nncoord(nn, 2)
      REAL(KIND = kreal) :: shapefunc(nn)
      REAL(KIND = kreal) :: shapederiv(nn, 2)
      REAL(KIND = kreal) :: aa1(3), aa2(3), aa3(3)
      REAL(KIND = kreal) :: bb1(3), bb2(3), bb3(3)
      REAL(KIND = kreal) :: cc1(3), cc2(3)
      REAL(KIND = kreal) :: alpha
      REAL(KIND = kreal) :: xxi_lx, eeta_lx
      REAL(KIND = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
      REAL(KIND = kreal) :: h(nn, 3)
      REAL(KIND = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
      REAL(KIND = kreal) :: v1_i(3), v2_i(3), v3_i(3)
      REAL(KIND = kreal) :: v1_abs, v2_abs, v3_abs
      REAL(KIND = kreal) :: a_over_2_v3(3, nn)
      REAL(KIND = kreal) :: u_rot(3, nn)
      REAL(KIND = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
                            dudzeta_rot(3, nn)                   
      REAL(KIND = kreal) :: g1(3), g2(3), g3(3)
      REAL(KIND = kreal) :: g1_abs, g2_abs, g3_abs
      REAL(KIND = kreal) :: e_0(3)
      REAL(KIND = kreal) :: cg1(3), cg2(3), cg3(3)
      REAL(KIND = kreal) :: det
      REAL(KIND = kreal) :: det_cg3(3)
      REAL(KIND = kreal) :: det_inv
      REAL(KIND = kreal) :: det_cg3_abs
      REAL(KIND = kreal) :: w_w_w_det
      REAL(KIND = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
      REAL(KIND = kreal) :: e1_hat_abs, e2_hat_abs
      REAL(KIND = kreal) :: Cv12(ndof*nn), Cv13(ndof*nn), &
                            Cv21(ndof*nn), Cv23(ndof*nn), &
                            Cv31(ndof*nn), Cv32(ndof*nn)  
      REAL(KIND = kreal) :: Cv_theta(ndof*nn), Cv_w(ndof*nn)
      REAL(KIND = kreal) :: Cv(ndof*nn)
      REAL(KIND = kreal) :: sigma_layer, thick_layer 

		sstable = 0
		flag_dof = 0
!--------------------------------------------------------------------
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN
       
       fetype = fe_mitc4_shell
       
       ny = 2
       
       ntying = 1
       npoints_tying(1)= 4
       
      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       fetype = fe_mitc9_shell
       
       ny = 3
       
       ntying = 3
       npoints_tying(1)= 6
       npoints_tying(2)= 6
       npoints_tying(3)= 4
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       fetype = fe_mitc3_shell
       
       ny = 2
       
       ntying = 1
       npoints_tying(1)= 3
       
      END IF
      
!--------------------------------------------------------------------
      
      IF( PRESENT( nddisp ) ) THEN
       
       unode(:, 1:nn) = nddisp(:, :)
       
      END IF
          
!--------------------------------------------------------------------
      
      flag = gausses(1)%pMaterial%nlgeom_flag
      
      IF( .NOT. PRESENT( nddisp ) ) flag = INFINITE
      
!--------------------------------------------------------------------
      
      elem(:, :) = ecoord(:, :)
      
!--------------------------------------------------------------------
      
      tmpstiff(:, :) = 0.0D0
      
!-------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      ! eta-coordinate at a node in a local element
      CALL getNodalNaturalCoord(fetype, nncoord)
      
!-------------------------------------------------------------------
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN
       
       !--------------------------------------------------------
       
       ! xi-coordinate at a tying pont in a local element
       tpcoord(1, 1, 1) =  0.0D0
       tpcoord(2, 1, 1) =  1.0D0
       tpcoord(3, 1, 1) =  0.0D0
       tpcoord(4, 1, 1) = -1.0D0
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = -1.0D0
       tpcoord(2, 2, 1) =  0.0D0
       tpcoord(3, 2, 1) =  1.0D0
       tpcoord(4, 2, 1) =  0.0D0
       
       !--------------------------------------------------------
       
      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       !--------------------------------------------------------
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(5, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(6, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(2, 2, 1) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(3, 2, 1) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(4, 2, 1) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(5, 2, 1) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(6, 2, 1) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 2) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(2, 1, 2) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(3, 1, 2) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(4, 1, 2) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(5, 1, 2) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(6, 1, 2) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(5, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(6, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 1, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 1, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 1, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 2, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 2, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 2, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       
       !--------------------------------------------------------
       
       ! Xi-coordinate at a tying point in a local element
       xxi_di(1, 1)  = -1.0D0
       xxi_di(2, 1)  =  1.0D0
       xxi_di(3, 1)  =  1.0D0
       xxi_di(4, 1)  = -1.0D0
       xxi_di(5, 1)  =  1.0D0
       xxi_di(6, 1)  = -1.0D0
       ! Eta-coordinate at a tying point in a local element
       eeta_di(1, 1) = -1.0D0
       eeta_di(2, 1) = -1.0D0
       eeta_di(3, 1) =  1.0D0
       eeta_di(4, 1) =  1.0D0
       eeta_di(5, 1) =  0.0D0
       eeta_di(6, 1) =  0.0D0
       
       ! Xi-coordinate at a tying point in a local element
       xxi_di(1, 2)  = -1.0D0
       xxi_di(2, 2)  =  0.0D0
       xxi_di(3, 2)  =  1.0D0
       xxi_di(4, 2)  =  1.0D0
       xxi_di(5, 2)  =  0.0D0
       xxi_di(6, 2)  = -1.0D0
       ! Eta-coordinate at a tying point in a local element
       eeta_di(1, 2) = -1.0D0
       eeta_di(2, 2) = -1.0D0
       eeta_di(3, 2) = -1.0D0
       eeta_di(4, 2) =  1.0D0
       eeta_di(5, 2) =  1.0D0
       eeta_di(6, 2) =  1.0D0
       
       !--------------------------------------------------------
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       !--------------------------------------------------------
        
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 1)  = 0.5D0
       tpcoord(2, 1, 1)  = 0.0D0
       tpcoord(3, 1, 1)  = 0.5D0
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = 0.0D0
       tpcoord(2, 2, 1) = 0.5D0
       tpcoord(3, 2, 1) = 0.5D0
       
       !--------------------------------------------------------
       
      END IF
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at the center point in a local element
      ! eta-coordinate at the center point in a local element
      naturalcoord(1) = 0.0D0
      naturalcoord(2) = 0.0D0
      
      CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
      
      !--------------------------------------------------------------
      
      ! Covariant basis vector
      DO i = 1, 3
       
       g1(i) = 0.0D0
       
       DO na = 1, nn
        
        g1(i) = g1(i)+shapederiv(na, 1) &
                      *elem(i, na)      
        
       END DO
       
      END DO
      
      e_0(1) = g1(1)
      e_0(2) = g1(2)
      e_0(3) = g1(3)
      
      !--------------------------------------------------------------
      
      DO nb = 1, nn
       
       !--------------------------------------------------------
       
       naturalcoord(1) = nncoord(nb, 1) 
       naturalcoord(2) = nncoord(nb, 2) 
       
       CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
       
       !--------------------------------------------------------
       
       ! Covariant basis vector
       DO i = 1, 3
        
        g1(i) = 0.0D0
        g2(i) = 0.0D0
        
        DO na = 1, nn
         
         g1(i) = g1(i)+shapederiv(na, 1) &
                       *elem(i, na)      
         g2(i) = g2(i)+shapederiv(na, 2) &
                       *elem(i, na)      
         
        END DO
        
       END DO
       
       !------------------------------------------
       
       det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
       det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
       det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)
       
       det_cg3_abs = DSQRT( det_cg3(1)*det_cg3(1)   &
                           +det_cg3(2)*det_cg3(2)   &
                           +det_cg3(3)*det_cg3(3) ) 
       
       v3(1, nb) = det_cg3(1)/det_cg3_abs
       v3(2, nb) = det_cg3(2)/det_cg3_abs
       v3(3, nb) = det_cg3(3)/det_cg3_abs
       
       !--------------------------------------------------------
       
       v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
       v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
       v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)
       
       v2_abs = DSQRT( v2(1, nb)*v2(1, nb)   &
                      +v2(2, nb)*v2(2, nb)   &
                      +v2(3, nb)*v2(3, nb) ) 
       
       IF( v2_abs .GT. 1.0D-15 ) THEN 
        
        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs
        
        v1(1, nb) = v2(2, nb)*v3(3, nb) &
                   -v2(3, nb)*v3(2, nb) 
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
                   -v2(1, nb)*v3(3, nb) 
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
                   -v2(2, nb)*v3(1, nb) 
        
        v1_abs = DSQRT( v1(1, nb)*v1(1, nb)   &
                       +v1(2, nb)*v1(2, nb)   &
                       +v1(3, nb)*v1(3, nb) ) 
        
        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs
        
       ELSE    ! YX: impossible
        
        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0
        
        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0
        
       END IF
       
       !---------------------------------------------------
       
       v3(1, nb) = v1(2, nb)*v2(3, nb) &
                  -v1(3, nb)*v2(2, nb) 
       v3(2, nb) = v1(3, nb)*v2(1, nb) &
                  -v1(1, nb)*v2(3, nb) 
       v3(3, nb) = v1(1, nb)*v2(2, nb) &
                  -v1(2, nb)*v2(1, nb) 
       
       v3_abs = DSQRT( v3(1, nb)*v3(1, nb)   &
                      +v3(2, nb)*v3(2, nb)   &
                      +v3(3, nb)*v3(3, nb) ) 
       
       v3(1, nb) = v3(1, nb)/v3_abs
       v3(2, nb) = v3(2, nb)/v3_abs
       v3(3, nb) = v3(3, nb)/v3_abs
       
       !--------------------------------------------------------
       
       a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
       a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
       a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)
       
       !--------------------------------------------------------
       
      END DO
      
!--------------------------------------------------------------------
!		MODIFIED to LAMINATED SHELL ANALYSIS
!--------------------------------------------------------------------

   n_total_layer =  int(gausses(1)%pMaterial%variables(M_TOTAL_LAYER))
   DO n_layer=1,n_total_layer
      DO ly = 1, ny
       
       !--------------------------------------------------------
       
       ! MITC4
       IF( etype .EQ. fe_mitc4_shell ) THEN
        
        zeta_ly = 0.0D0
        
       ! MITC9
       ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
        
        zeta_ly = xg(ny, ly)
        
       ! MITC3
       ELSE IF( etype .EQ. fe_mitc3_shell )THEN
        
        zeta_ly = 0.0D0
        
       END IF
       
       !--------------------------------------------------------
       
       DO it = 1, ntying
        
        DO ip = 1, npoints_tying(it)
         
         !-------------------------------------------------
         
         naturalcoord(1) = tpcoord(ip, 1, it)
         naturalcoord(2) = tpcoord(ip, 2, it) 
         
         CALL getShapeFunc(fetype, naturalcoord, shapefunc)
         
         CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
         
         !-------------------------------------------------
         
         DO na = 1, nn
          
          DO i = 1, 3
           
           u_rot(i, na)                      &
           = shapefunc(na)                   &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           
           dudxi_rot(i, na)                  &
           = shapederiv(na, 1)               &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           dudeta_rot(i, na)                 &
           = shapederiv(na, 2)               &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           dudzeta_rot(i, na)                &
           = shapefunc(na)                   &
             *( a_over_2_v3(i, na) )         
           
          END DO
          
         END DO
         
         !-------------------------------------------------
         
         ! Covariant basis vector
         DO i = 1, 3
          
          g1(i) = 0.0D0
          g2(i) = 0.0D0
          g3(i) = 0.0D0
          
          DO na = 1, nn
           
           g1(i) = g1(i)+shapederiv(na, 1)  &
                         *elem(i, na)       &
                        +dudxi_rot(i, na)   
           g2(i) = g2(i)+shapederiv(na, 2)  &
                         *elem(i, na)       &
                        +dudeta_rot(i, na)  
           g3(i) = g3(i)+dudzeta_rot(i, na) 
           
          END DO
          
         END DO
         
         !-------------------------------------------------
         
         ! [ B L ] matrix
         DO nb = 1, nn
          
          jsize1 = ndof*(nb-1)+1
          jsize2 = ndof*(nb-1)+2
          jsize3 = ndof*(nb-1)+3
          jsize4 = ndof*(nb-1)+4
          jsize5 = ndof*(nb-1)+5
          jsize6 = ndof*(nb-1)+6
          
          aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
          aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
          aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)
          
          aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
          aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
          aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)
          
          aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
          aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
          aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)
          
          bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
          bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
          bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)
          
          bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
          bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
          bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)
          
          bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
          bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
          bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)
          
          cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
          cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
          cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)
          
          cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
          cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
          cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)
          
          B_di(1, jsize1, ip, it, ly) = shapederiv(nb, 1)*g1(1)
          B_di(2, jsize1, ip, it, ly) = shapederiv(nb, 2)*g2(1)
          B_di(3, jsize1, ip, it, ly) = shapederiv(nb, 1)*g2(1) &
                                       +shapederiv(nb, 2)*g1(1) 
          B_di(4, jsize1, ip, it, ly) = shapederiv(nb, 2)*g3(1)
          B_di(5, jsize1, ip, it, ly) = shapederiv(nb, 1)*g3(1)
          
          B_di(1, jsize2, ip, it, ly) = shapederiv(nb, 1)*g1(2)
          B_di(2, jsize2, ip, it, ly) = shapederiv(nb, 2)*g2(2)
          B_di(3, jsize2, ip, it, ly) = shapederiv(nb, 1)*g2(2) &
                                       +shapederiv(nb, 2)*g1(2) 
          B_di(4, jsize2, ip, it, ly) = shapederiv(nb, 2)*g3(2)
          B_di(5, jsize2, ip, it, ly) = shapederiv(nb, 1)*g3(2)
          
          B_di(1, jsize3, ip, it, ly) = shapederiv(nb, 1)*g1(3)
          B_di(2, jsize3, ip, it, ly) = shapederiv(nb, 2)*g2(3)
          B_di(3, jsize3, ip, it, ly) = shapederiv(nb, 1)*g2(3) &
                                       +shapederiv(nb, 2)*g1(3) 
          B_di(4, jsize3, ip, it, ly) = shapederiv(nb, 2)*g3(3)
          B_di(5, jsize3, ip, it, ly) = shapederiv(nb, 1)*g3(3)
          
          B_di(1, jsize4, ip, it, ly) = aa1(1)
          B_di(2, jsize4, ip, it, ly) = bb2(1)
          B_di(3, jsize4, ip, it, ly) = aa2(1)+bb1(1)
          B_di(4, jsize4, ip, it, ly) = bb3(1)+cc2(1)
          B_di(5, jsize4, ip, it, ly) = aa3(1)+cc1(1)
          
          B_di(1, jsize5, ip, it, ly) = aa1(2)
          B_di(2, jsize5, ip, it, ly) = bb2(2)
          B_di(3, jsize5, ip, it, ly) = aa2(2)+bb1(2)
          B_di(4, jsize5, ip, it, ly) = bb3(2)+cc2(2)
          B_di(5, jsize5, ip, it, ly) = aa3(2)+cc1(2)
          
          B_di(1, jsize6, ip, it, ly) = aa1(3)
          B_di(2, jsize6, ip, it, ly) = bb2(3)
          B_di(3, jsize6, ip, it, ly) = aa2(3)+bb1(3)
          B_di(4, jsize6, ip, it, ly) = bb3(3)+cc2(3)
          B_di(5, jsize6, ip, it, ly) = aa3(3)+cc1(3)
          
         END DO
         
         !-------------------------------------------------
         
        END DO
        
       END DO
       
       !--------------------------------------------------------
       
        IF ( n_total_layer.ge. 2 ) THEN
		     shellmatl = int(gausses(1)%pMaterial%variables(M_SHELL_MATLTYPE))
	        if (shellmatl == 0)then
			   sigma_layer = 0.0D0
		       DO m = 1, n_layer
			      sigma_layer = sigma_layer + 2 * gausses(1)%pMaterial%variables(100+3*m) / thick
		       END DO
			   zeta_ly = -1 + sigma_layer - gausses(1)%pMaterial%variables(100+3*n_layer) / thick * (1-xg(ny, ly))
	        elseif (shellmatl == 1)then
			   sigma_layer = 0.0D0
		       DO m = 1, n_layer
			      sigma_layer = sigma_layer + 2 * gausses(1)%pMaterial%variables(100+8*m-5) / thick
		       END DO
               zeta_ly = -1 + sigma_layer - gausses(1)%pMaterial%variables(100+8*n_layer-5) / thick * (1-xg(ny, ly))			   
	        else
			   write(*,*)"ERROR : shellmatl isnot correct"; stop
            endif
        ELSE
			zeta_ly = xg(ny, ly)
        ENDIF
	   
		w_ly    = wgt(ny, ly)

       !--------------------------------------------------------
       
       DO lx = 1, NumOfQuadPoints(fetype)
        
        !--------------------------------------------------
        
        CALL getQuadPoint(fetype, lx, naturalcoord)
        
        xi_lx  = naturalcoord(1)
        eta_lx = naturalcoord(2)
        
        w_w_lx = getWeight(fetype, lx)
        
        CALL getShapeFunc(fetype, naturalcoord, shapefunc)
        
        CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
        
        !--------------------------------------------------
        
        DO i = 1, 3
          
         v1_i(i) = 0.0D0
         v2_i(i) = 0.0D0
         v3_i(i) = 0.0D0
        
         DO na = 1, nn
          
          v1_i(i) = v1_i(i)+shapefunc(na)*v1(i, na)
          v2_i(i) = v2_i(i)+shapefunc(na)*v2(i, na)
          v3_i(i) = v3_i(i)+shapefunc(na)*v3(i, na)
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        DO na = 1, nn
         
         DO i = 1, 3
          
          u_rot(i, na)                      &
          = shapefunc(na)                   &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          
          dudxi_rot(i, na)                  &
          = shapederiv(na, 1)               &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          dudeta_rot(i, na)                 &
          = shapederiv(na, 2)               &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          dudzeta_rot(i, na)                &
          = shapefunc(na)                   &
            *( a_over_2_v3(i, na) )         
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        ! Covariant basis vector
        DO i = 1, 3
         
         g1(i) = 0.0D0
         g2(i) = 0.0D0
         g3(i) = 0.0D0
         
         DO na = 1, nn
          
          g1(i) = g1(i)+shapederiv(na, 1)  &
                        *elem(i, na)       &
                       +dudxi_rot(i, na)   
          g2(i) = g2(i)+shapederiv(na, 2)  &
                        *elem(i, na)       &
                       +dudeta_rot(i, na)  
          g3(i) = g3(i)+dudzeta_rot(i, na) 
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        ! Jacobian
        det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
             +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
             +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) ) 
        
        det_inv = 1.0D0/det
        
        !--------------------------------------------------
        
        ! Contravariant basis vector
        cg1(1) = det_inv                      &
                 *( g2(2)*g3(3)-g2(3)*g3(2) ) 
        cg1(2) = det_inv                      &
                 *( g2(3)*g3(1)-g2(1)*g3(3) ) 
        cg1(3) = det_inv                      &
                 *( g2(1)*g3(2)-g2(2)*g3(1) ) 
        cg2(1) = det_inv                      &
                 *( g3(2)*g1(3)-g3(3)*g1(2) ) 
        cg2(2) = det_inv                      &
                 *( g3(3)*g1(1)-g3(1)*g1(3) ) 
        cg2(3) = det_inv                      &
                 *( g3(1)*g1(2)-g3(2)*g1(1) ) 
        cg3(1) = det_inv                      &
                 *( g1(2)*g2(3)-g1(3)*g2(2) ) 
        cg3(2) = det_inv                      &
                 *( g1(3)*g2(1)-g1(1)*g2(3) ) 
        cg3(3) = det_inv                      &
                 *( g1(1)*g2(2)-g1(2)*g2(1) ) 
        
        !--------------------------------------------------
        
        g3_abs = DSQRT( g3(1)*g3(1)   &
                       +g3(2)*g3(2)   &
                       +g3(3)*g3(3) ) 
        
        !--------------------------------------------------
        
        ! Orthonormal vectors
        
        e3_hat(1) = g3(1)/g3_abs
        e3_hat(2) = g3(2)/g3_abs
        e3_hat(3) = g3(3)/g3_abs
        
        e1_hat(1) = g2(2)*e3_hat(3) &
                   -g2(3)*e3_hat(2) 
        e1_hat(2) = g2(3)*e3_hat(1) &
                   -g2(1)*e3_hat(3) 
        e1_hat(3) = g2(1)*e3_hat(2) &
                   -g2(2)*e3_hat(1) 
        e1_hat_abs = DSQRT( e1_hat(1)*e1_hat(1)   &
                           +e1_hat(2)*e1_hat(2)   &
                           +e1_hat(3)*e1_hat(3) ) 
        e1_hat(1) = e1_hat(1)/e1_hat_abs
        e1_hat(2) = e1_hat(2)/e1_hat_abs
        e1_hat(3) = e1_hat(3)/e1_hat_abs
        
        e2_hat(1) = e3_hat(2)*e1_hat(3) &
                   -e3_hat(3)*e1_hat(2) 
        e2_hat(2) = e3_hat(3)*e1_hat(1) &
                   -e3_hat(1)*e1_hat(3) 
        e2_hat(3) = e3_hat(1)*e1_hat(2) &
                   -e3_hat(2)*e1_hat(1)  
        e2_hat_abs = DSQRT( e2_hat(1)*e2_hat(1)   &
                           +e2_hat(2)*e2_hat(2)   &
                           +e2_hat(3)*e2_hat(3) ) 
        e2_hat(1) = e2_hat(1)/e2_hat_abs
        e2_hat(2) = e2_hat(2)/e2_hat_abs
        e2_hat(3) = e2_hat(3)/e2_hat_abs
        
        !--------------------------------------------------
    !    write(*,*)  'Matlmatrix_layer', n_layer !下のgaussが非常によくわからない
        CALL MatlMatrix_Shell                        &
             (gausses(lx), Shell, D,                 &
              e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
              alpha, n_layer)                                 
        
        !--------------------------------------------------
        
        ! [ B L ] matrix
        DO nb = 1, nn
         
         jsize1 = ndof*(nb-1)+1
         jsize2 = ndof*(nb-1)+2
         jsize3 = ndof*(nb-1)+3
         jsize4 = ndof*(nb-1)+4
         jsize5 = ndof*(nb-1)+5
         jsize6 = ndof*(nb-1)+6
         
         aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
         aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
         aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)
         
         aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
         aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
         aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)
         
         aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
         aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
         aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)
         
         bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
         bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
         bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)
         
         bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
         bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
         bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)
         
         bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
         bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
         bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)
         
         cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
         cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
         cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)
         
         cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
         cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
         cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)
         
         B(1, jsize1) = shapederiv(nb, 1)*g1(1)
         B(2, jsize1) = shapederiv(nb, 2)*g2(1)
         B(3, jsize1) = shapederiv(nb, 1)*g2(1) &
                       +shapederiv(nb, 2)*g1(1) 
         B(4, jsize1) = shapederiv(nb, 2)*g3(1)
         B(5, jsize1) = shapederiv(nb, 1)*g3(1)
         
         B(1, jsize2) = shapederiv(nb, 1)*g1(2)
         B(2, jsize2) = shapederiv(nb, 2)*g2(2)
         B(3, jsize2) = shapederiv(nb, 1)*g2(2) &
                       +shapederiv(nb, 2)*g1(2) 
         B(4, jsize2) = shapederiv(nb, 2)*g3(2)
         B(5, jsize2) = shapederiv(nb, 1)*g3(2)
         
         B(1, jsize3) = shapederiv(nb, 1)*g1(3)
         B(2, jsize3) = shapederiv(nb, 2)*g2(3)
         B(3, jsize3) = shapederiv(nb, 1)*g2(3) &
                       +shapederiv(nb, 2)*g1(3) 
         B(4, jsize3) = shapederiv(nb, 2)*g3(3)
         B(5, jsize3) = shapederiv(nb, 1)*g3(3)
         
         B(1, jsize4) = aa1(1)
         B(2, jsize4) = bb2(1)
         B(3, jsize4) = aa2(1)+bb1(1)
         B(4, jsize4) = bb3(1)+cc2(1)
         B(5, jsize4) = aa3(1)+cc1(1)
         
         B(1, jsize5) = aa1(2)
         B(2, jsize5) = bb2(2)
         B(3, jsize5) = aa2(2)+bb1(2)
         B(4, jsize5) = bb3(2)+cc2(2)
         B(5, jsize5) = aa3(2)+cc1(2)
         
         B(1, jsize6) = aa1(3)
         B(2, jsize6) = bb2(3)
         B(3, jsize6) = aa2(3)+bb1(3)
         B(4, jsize6) = bb3(3)+cc2(3)
         B(5, jsize6) = aa3(3)+cc1(3)
         
        END DO
        
        !--------------------------------------------------
        
        ! MITC4
        IF( etype .EQ. fe_mitc4_shell ) THEN
                 
         DO jsize = 1, ndof*nn
          
          B(4, jsize) = 0.0D0
          B(5, jsize) = 0.0D0
          
          ! B_as(4, jsize)
          B(4, jsize)                                       &
          = 0.5D0*( 1.0D0-xi_lx  )*B_di(4, jsize, 4, 1, ly) &
           +0.5D0*( 1.0D0+xi_lx  )*B_di(4, jsize, 2, 1, ly) 
          ! B_as(5, jsize)
          B(5, jsize)                                       &
          = 0.5D0*( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
           +0.5D0*( 1.0D0+eta_lx )*B_di(5, jsize, 3, 1, ly) 
          
         END DO
         
        ! MITC9
        ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
         
         xxi_lx  = xi_lx /DSQRT( 1.0D0/3.0D0 )
         eeta_lx = eta_lx/DSQRT( 3.0D0/5.0D0 )
         
         DO ip = 1, npoints_tying(1)
          
          h(ip, 1)                                     &
          = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
            *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
               *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
              +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
               *( 1.0D0-eeta_lx*eeta_lx ) )            
          
         END DO
         
         xxi_lx  = xi_lx /DSQRT( 3.0D0/5.0D0 )
         eeta_lx = eta_lx/DSQRT( 1.0D0/3.0D0 )
         
         DO ip = 1, npoints_tying(2)
          
          h(ip, 2)                                      &
          = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
              *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
             +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
              *( 1.0D0-xxi_lx*xxi_lx ) )                &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) ) 
          
         END DO
         
         xxi_lx  = xi_lx /DSQRT( 1.0D0/3.0D0 )
         eeta_lx = eta_lx/DSQRT( 1.0D0/3.0D0 )
         
         DO ip = 1, npoints_tying(3)
          
          h(ip, 3)                                      &
          = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) ) 
          
         END DO
         
         DO jsize = 1, ndof*nn
          
          B(1, jsize) = 0.0D0
          B(2, jsize) = 0.0D0
          B(3, jsize) = 0.0D0
          B(4, jsize) = 0.0D0
          B(5, jsize) = 0.0D0
          
          DO ip = 1, npoints_tying(1)
           
           ! B_as(1, jsize)
           B(1, jsize)                                      &
           = B(1, jsize)+h(ip, 1)*B_di(1, jsize, ip, 1, ly) 
           ! B_as(5, jsize)
           B(5, jsize)                                      &
           = B(5, jsize)+h(ip, 1)*B_di(5, jsize, ip, 1, ly) 
           
          END DO
          
          DO ip = 1, npoints_tying(2)
           
           ! B_as(2, jsize)
           B(2, jsize)                                      &
           = B(2, jsize)+h(ip, 2)*B_di(2, jsize, ip, 2, ly) 
           ! B_as(4, jsize)
           B(4, jsize)                                      &
           = B(4, jsize)+h(ip, 2)*B_di(4, jsize, ip, 2, ly) 
           
          END DO
          
          DO ip = 1, npoints_tying(3)
           
           ! B_as(3, jsize)
           B(3, jsize)                                      &
           = B(3, jsize)+h(ip, 3)*B_di(3, jsize, ip, 3, ly) 
           
          END DO
          
         END DO
         
        ! MITC3
        ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
         
         DO jsize = 1, ndof*nn
          
          B(4, jsize) = 0.0D0
          B(5, jsize) = 0.0D0
          
          ! B_as(4, jsize)
          B(4, jsize)                                 &
          = ( 1.0D0-xi_lx  )*B_di(4, jsize, 2, 1, ly) &
           +xi_lx *B_di(5, jsize, 1, 1, ly)           &
           +xi_lx *( B_di(4, jsize, 3, 1, ly)         &
                    -B_di(5, jsize, 3, 1, ly)  )      
          
          ! B_as(5, jsize)
          B(5, jsize)                                 &
          = eta_lx*B_di(4, jsize, 2, 1, ly)           &
           +( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
           -eta_lx*( B_di(4, jsize, 3, 1, ly)         &
                    -B_di(5, jsize, 3, 1, ly)  )      
          
         END DO
         
        END IF
        
        !--------------------------------------------------
        
        w_w_w_det = w_w_lx*w_ly*det
        
        !--------------------------------------------------
        
        DB(1:5, 1:ndof*nn) = MATMUL( D, B(1:5, 1:ndof*nn ) )
        
        !--------------------------------------------------
        
          shellmatl = int(gausses(1)%pMaterial%variables(M_SHELL_MATLTYPE))
	        if (shellmatl == 0)then
			    FORALL( isize=1:ndof*nn, jsize=1:ndof*nn )
				   tmpstiff(isize, jsize)                               &
				   = tmpstiff(isize, jsize)                             &
				   +w_w_w_det*gausses(1)%pMaterial%variables(100+3*n_layer)/thick*DOT_PRODUCT( B(:, isize), DB(:, jsize) ) 
		       END FORALL
		    elseif (shellmatl == 1)then
			    FORALL( isize=1:ndof*nn, jsize=1:ndof*nn )
				   tmpstiff(isize, jsize)                               &
				   = tmpstiff(isize, jsize)                             &
				   +w_w_w_det*gausses(1)%pMaterial%variables(100+8*n_layer-5)/thick*DOT_PRODUCT( B(:, isize), DB(:, jsize) ) 
	           END FORALL
		    else
			   write(*,*)"ERROR : shellmatl isnot correct"; stop
            endif
        
        !--------------------------------------------------
        
        ! [ B_{i} ] matrix
        DO nb = 1, nn
         
         jsize1 = ndof*(nb-1)+1
         jsize2 = ndof*(nb-1)+2
         jsize3 = ndof*(nb-1)+3
         jsize4 = ndof*(nb-1)+4
         jsize5 = ndof*(nb-1)+5
         jsize6 = ndof*(nb-1)+6
         
         B1(1, jsize1) =  shapederiv(nb, 1)
         B1(2, jsize1) =  0.0D0
         B1(3, jsize1) =  0.0D0
         B1(1, jsize2) =  0.0D0
         B1(2, jsize2) =  shapederiv(nb, 1)
         B1(3, jsize2) =  0.0D0
         B1(1, jsize3) =  0.0D0
         B1(2, jsize3) =  0.0D0
         B1(3, jsize3) =  shapederiv(nb, 1)
         B1(1, jsize4) =  0.0D0
         B1(2, jsize4) = -dudxi_rot(3, nb)
         B1(3, jsize4) =  dudxi_rot(2, nb)
         B1(1, jsize5) =  dudxi_rot(3, nb)
         B1(2, jsize5) =  0.0D0
         B1(3, jsize5) = -dudxi_rot(1, nb)
         B1(1, jsize6) = -dudxi_rot(2, nb)
         B1(2, jsize6) =  dudxi_rot(1, nb)
         B1(3, jsize6) =  0.0D0
         
         B2(1, jsize1) =  shapederiv(nb, 2)
         B2(2, jsize1) =  0.0D0
         B2(3, jsize1) =  0.0D0
         B2(1, jsize2) =  0.0D0
         B2(2, jsize2) =  shapederiv(nb, 2)
         B2(3, jsize2) =  0.0D0
         B2(1, jsize3) =  0.0D0
         B2(2, jsize3) =  0.0D0
         B2(3, jsize3) =  shapederiv(nb, 2)
         B2(1, jsize4) =  0.0D0
         B2(2, jsize4) = -dudeta_rot(3, nb)
         B2(3, jsize4) =  dudeta_rot(2, nb)
         B2(1, jsize5) =  dudeta_rot(3, nb)
         B2(2, jsize5) =  0.0D0
         B2(3, jsize5) = -dudeta_rot(1, nb)
         B2(1, jsize6) = -dudeta_rot(2, nb)
         B2(2, jsize6) =  dudeta_rot(1, nb)
         B2(3, jsize6) =  0.0D0
         
         B3(1, jsize1) =  0.0D0
         B3(2, jsize1) =  0.0D0
         B3(3, jsize1) =  0.0D0
         B3(1, jsize2) =  0.0D0
         B3(2, jsize2) =  0.0D0
         B3(3, jsize2) =  0.0D0
         B3(1, jsize3) =  0.0D0
         B3(2, jsize3) =  0.0D0
         B3(3, jsize3) =  0.0D0
         B3(1, jsize4) =  0.0D0
         B3(2, jsize4) = -dudzeta_rot(3, nb)
         B3(3, jsize4) =  dudzeta_rot(2, nb)
         B3(1, jsize5) =  dudzeta_rot(3, nb)
         B3(2, jsize5) =  0.0D0
         B3(3, jsize5) = -dudzeta_rot(1, nb)
         B3(1, jsize6) = -dudzeta_rot(2, nb)
         B3(2, jsize6) =  dudzeta_rot(1, nb)
         B3(3, jsize6) =  0.0D0
         
        END DO
        
        !--------------------------------------------------
        
        ! { C_{ij} } vector
        DO jsize = 1, ndof*nn
         
         Cv12(jsize) = ( cg1(1)*B1(2, jsize)   &
                        +cg2(1)*B2(2, jsize)   &
                        +cg3(1)*B3(2, jsize) ) &
                      -( cg1(2)*B1(1, jsize)   &
                        +cg2(2)*B2(1, jsize)   &
                        +cg3(2)*B3(1, jsize) ) 
         Cv13(jsize) = ( cg1(1)*B1(3, jsize)   &
                        +cg2(1)*B2(3, jsize)   &
                        +cg3(1)*B3(3, jsize) ) &
                      -( cg1(3)*B1(1, jsize)   &
                        +cg2(3)*B2(1, jsize)   &
                        +cg3(3)*B3(1, jsize) ) 
         Cv21(jsize) = ( cg1(2)*B1(1, jsize)   &
                        +cg2(2)*B2(1, jsize)   &
                        +cg3(2)*B3(1, jsize) ) &
                      -( cg1(1)*B1(2, jsize)   &
                        +cg2(1)*B2(2, jsize)   &
                        +cg3(1)*B3(2, jsize) ) 
         Cv23(jsize) = ( cg1(2)*B1(3, jsize)   &
                        +cg2(2)*B2(3, jsize)   &
                        +cg3(2)*B3(3, jsize) ) &
                      -( cg1(3)*B1(2, jsize)   &
                        +cg2(3)*B2(2, jsize)   &
                        +cg3(3)*B3(2, jsize) ) 
         Cv31(jsize) = ( cg1(3)*B1(1, jsize)   &
                        +cg2(3)*B2(1, jsize)   &
                        +cg3(3)*B3(1, jsize) ) &
                      -( cg1(1)*B1(3, jsize)   &
                        +cg2(1)*B2(3, jsize)   &
                        +cg3(1)*B3(3, jsize) ) 
         Cv32(jsize) = ( cg1(3)*B1(2, jsize)   &
                        +cg2(3)*B2(2, jsize)   &
                        +cg3(3)*B3(2, jsize) ) &
                      -( cg1(2)*B1(3, jsize)   &
                        +cg2(2)*B2(3, jsize)   &
                        +cg3(2)*B3(3, jsize) ) 
						
        END DO
        
        !--------------------------------------------------
        
        ! { Cw } vector
        DO nb = 1, nn
         
         DO j = 1, ndof
          
          jsize = ndof*(nb-1)+j
          
          Cv_w(jsize)                       &
          = v1_i(1)*Cv12(jsize)*v2_i(2) &
           +v1_i(1)*Cv13(jsize)*v2_i(3) &
           +v1_i(2)*Cv21(jsize)*v2_i(1) &
           +v1_i(2)*Cv23(jsize)*v2_i(3) &
           +v1_i(3)*Cv31(jsize)*v2_i(1) &
           +v1_i(3)*Cv32(jsize)*v2_i(2) 
          
         END DO
         
        END DO
        
        ! { Ctheta } vector
        DO nb = 1, nn
         
         jsize1 = ndof*(nb-1)+1
         jsize2 = ndof*(nb-1)+2
         jsize3 = ndof*(nb-1)+3
         jsize4 = ndof*(nb-1)+4
         jsize5 = ndof*(nb-1)+5
         jsize6 = ndof*(nb-1)+6
         
         Cv_theta(jsize1) = 0.0D0
         Cv_theta(jsize2) = 0.0D0
         Cv_theta(jsize3) = 0.0D0
         Cv_theta(jsize4) = v3_i(1)*shapefunc(nb)
         Cv_theta(jsize5) = v3_i(2)*shapefunc(nb)
         Cv_theta(jsize6) = v3_i(3)*shapefunc(nb)
         
        END DO
        
        ! { C } vector
        DO jsize = 1, ndof*nn
         
         Cv(jsize) = Cv_theta(jsize)-0.5D0*Cv_w(jsize)
         
        END DO
        
        !--------------------------------------------------
        
        ! [ K L ] matrix
		  shellmatl = int(gausses(1)%pMaterial%variables(M_SHELL_MATLTYPE))
	        if (shellmatl == 0)then
			   DO jsize = 1, ndof*nn
              DO isize = 1, ndof*nn
          
          tmpstiff(isize, jsize)                &
          = tmpstiff(isize, jsize)              &
           +w_w_w_det*gausses(1)%pMaterial%variables(100+3*n_layer)/thick*alpha*Cv(isize)*Cv(jsize) 
          
              END DO
              END DO
	        elseif (shellmatl == 1)then
			   DO jsize = 1, ndof*nn
			   DO isize = 1, ndof*nn
			  
			  tmpstiff(isize, jsize)                &
			  = tmpstiff(isize, jsize)              &
			   +w_w_w_det*gausses(1)%pMaterial%variables(100+8*n_layer-5)/thick*alpha*Cv(isize)*Cv(jsize) 
			  
			   END DO
			   END DO
	        else
			   write(*,*)"ERROR : shellmatl isnot correct"; stop
            endif
        
        !--------------------------------------------------
        
       END DO
       
       !--------------------------------------------------------
       
      END DO
      
      !--------------------------------------------------------------
      
      stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)

!--------------------------------------------------------------------
		END DO     !< LAMINATED SHELL ANALYSIS
		
!		write(*,"(24E25.16)") stiff
		
!******************** Shell-Solid mixed analysis ********************  
!		mixglaf = 0 < natural shell (6 dof)
!		mixglaf = 1 < mixed 361 (3*2 dof)*8 nod
!		mixglaf = 2 < mixed 351 (3*2 dof)*6 nod
       if( mixflag == 1 )then
	   
!       write(*,*) 'convert for shell-solid mixed analysis'
		sstable(1) = 1
		sstable(2) = 2
		sstable(3) = 3
		sstable(4) = 7
		sstable(5) = 8
		sstable(6) = 9
		sstable(7) = 13
		sstable(8) = 14
		sstable(9) = 15
		sstable(10)= 19
		sstable(11)= 20
		sstable(12)= 21
		sstable(13)= 4
		sstable(14)= 5
		sstable(15)= 6
		sstable(16)= 10
		sstable(17)= 11
		sstable(18)= 12
		sstable(19)= 16
		sstable(20)= 17
		sstable(21)= 18
		sstable(22)= 22
		sstable(23)= 23
		sstable(24)= 24
		
		tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)
		
       DO i = 1, nn*ndof
       DO j = 1, nn*ndof   
			stiff(i,j) = tmpstiff(sstable(i),sstable(j))
       ENDDO
       ENDDO
	   
	   elseif( mixflag == 2 )then
!		write(*,*) 'convert for shell-solid mixed analysis 351'
		sstable(1) = 1
		sstable(2) = 2
		sstable(3) = 3
		sstable(4) = 7
		sstable(5) = 8
		sstable(6) = 9
		sstable(7) = 13
		sstable(8) = 14
		sstable(9) = 15
		sstable(10)= 4
		sstable(11)= 5
		sstable(12)= 6
		sstable(13)= 10
		sstable(14)= 11
		sstable(15)= 12
		sstable(16)= 16
		sstable(17)= 17
		sstable(18)= 18
		
		tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)
		
       DO i = 1, nn*ndof
       DO j = 1, nn*ndof   
			stiff(i,j) = tmpstiff(sstable(i),sstable(j))
       ENDDO
       ENDDO
	   
	   endif
	   
      RETURN
      
!####################################################################
      END SUBROUTINE STF_Shell_MITC
!####################################################################
      
      
!####################################################################
      SUBROUTINE ElementStress_Shell_MITC                  &
                 (etype, nn, ndof, ecoord, gausses, edisp, &
                  strain, stress, thick, zeta, n_layer, n_total_layer)             
!####################################################################
      
      USE mMechGauss
      USE gauss_integration
      USE m_MatMatrix
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype
      INTEGER(KIND = kint), INTENT(IN) :: nn
      INTEGER(KIND = kint), INTENT(IN) :: ndof
      REAL(KIND = kreal), INTENT(IN)   :: ecoord(3, nn)
      TYPE(tGaussStatus), INTENT(IN)   :: gausses(:)
      REAL(KIND = kreal), INTENT(IN)   :: edisp(6, nn)
      REAL(KIND = kreal), INTENT(OUT)  :: strain(:,:)
      REAL(KIND = kreal), INTENT(OUT)  :: stress(:,:)
      REAL(KIND = kreal), INTENT(IN)   :: thick
      REAL(KIND = kreal), INTENT(IN)   :: zeta
      
!--------------------------------------------------------------------
      
      INTEGER :: i, j, m
      INTEGER :: lx
      INTEGER :: fetype
      INTEGER :: ntying
      INTEGER :: npoints_tying(3)
      INTEGER :: it, ip
      INTEGER :: na, nb
      INTEGER :: isize, jsize
      INTEGER :: jsize1, jsize2, jsize3, &
                 jsize4, jsize5, jsize6
      INTEGER :: n_layer, n_total_layer, shellmatl
      
      REAL(KIND = kreal) :: D(5, 5)
      REAL(KIND = kreal) :: elem(3, nn)
      REAL(KIND = kreal) :: xi_lx, eta_lx, zeta_ly
      REAL(KIND = kreal) :: B_di(5, ndof*nn, 6, 3)
      REAL(KIND = kreal) :: B1(3, ndof*nn), B2(3, ndof*nn), &
                            B3(3, ndof*nn)                  
      REAL(KIND = kreal) :: naturalcoord(2)
      REAL(KIND = kreal) :: tpcoord(6, 2, 3)
      REAL(KIND = kreal) :: nncoord(nn, 2)
      REAL(KIND = kreal) :: shapefunc(nn)
      REAL(KIND = kreal) :: shapederiv(nn, 2)
      REAL(KIND = kreal) :: alpha
      REAL(KIND = kreal) :: xxi_lx, eeta_lx
      REAL(KIND = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
      REAL(KIND = kreal) :: h(nn, 3)
      REAL(KIND = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
      REAL(KIND = kreal) :: v1_abs, v2_abs, v3_abs
      REAL(KIND = kreal) :: a_over_2_v3(3, nn)
      REAL(KIND = kreal) :: a_over_2_theta_cross_v3(3, nn)
      REAL(KIND = kreal) :: u_rot(3, nn)
      REAL(KIND = kreal) :: theta(3, nn)
      REAL(KIND = kreal) :: dudxi(3), dudeta(3), dudzeta(3)
      REAL(KIND = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
                            dudzeta_rot(3, nn)                   
      REAL(KIND = kreal) :: g1(3), g2(3), g3(3)
      REAL(KIND = kreal) :: g1_abs, g2_abs, g3_abs
      REAL(KIND = kreal) :: e_0(3)
      REAL(KIND = kreal) :: cg1(3), cg2(3), cg3(3)
      REAL(KIND = kreal) :: det
      REAL(KIND = kreal) :: det_cg3(3)
      REAL(KIND = kreal) :: det_inv
      REAL(KIND = kreal) :: det_cg3_abs
      REAL(KIND = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
      REAL(KIND = kreal) :: e1_hat_abs, e2_hat_abs
      REAL(KIND = kreal) :: e11, e22, e12_2, e23_2, e31_2
      REAL(KIND = kreal) :: e11_di(6, 3), e22_di(6, 3),     &
                            e12_di_2(6, 3), e23_di_2(6, 3), &
                            e31_di_2(6, 3)                  
      REAL(KIND = kreal) :: E(3, 3), Ev(5)
      REAL(KIND = kreal) :: S(3, 3), Sv(5)
      REAL(KIND = kreal) :: sigma_layer
      
!--------------------------------------------------------------------

      ! for lamina stress 
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN
       
       fetype = fe_mitc4_shell
       
       ntying = 1
       npoints_tying(1)= 4

      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       fetype = fe_mitc9_shell
       
       ntying = 3
       npoints_tying(1)= 6
       npoints_tying(2)= 6
       npoints_tying(3)= 4
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       fetype = fe_mitc3_shell
       
       ntying = 1
       npoints_tying(1)= 3
       
      END IF
      
!--------------------------------------------------------------------
      
      elem(:, :) = ecoord(:, :)
	  
!--------------------------------------------------------------------
      
      DO na = 1, nn
       
       theta(1, na) = edisp(4, na)
       theta(2, na) = edisp(5, na)
       theta(3, na) = edisp(6, na)
       
      END DO
	  
!-------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      ! eta-coordinate at a node in a local element
      CALL getNodalNaturalCoord(fetype, nncoord)
      
!-------------------------------------------------------------------
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN
       
       !--------------------------------------------------------
       
       ! xi-coordinate at a tying pont in a local element
       tpcoord(1, 1, 1) =  0.0D0
       tpcoord(2, 1, 1) =  1.0D0
       tpcoord(3, 1, 1) =  0.0D0
       tpcoord(4, 1, 1) = -1.0D0
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = -1.0D0
       tpcoord(2, 2, 1) =  0.0D0
       tpcoord(3, 2, 1) =  1.0D0
       tpcoord(4, 2, 1) =  0.0D0
       
       !--------------------------------------------------------
       
      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       !--------------------------------------------------------
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(5, 1, 1) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(6, 1, 1) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(2, 2, 1) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(3, 2, 1) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(4, 2, 1) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(5, 2, 1) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(6, 2, 1) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 2) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(2, 1, 2) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(3, 1, 2) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(4, 1, 2) =  1.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(5, 1, 2) =  0.0D0*DSQRT( 3.0D0/5.0D0 )
       tpcoord(6, 1, 2) = -1.0D0*DSQRT( 3.0D0/5.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 2, 2) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(5, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(6, 2, 2) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 1, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 1, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 1, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(2, 2, 3) = -1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(3, 2, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       tpcoord(4, 2, 3) =  1.0D0*DSQRT( 1.0D0/3.0D0 )
       
       !--------------------------------------------------------
       
       ! Xi-coordinate at a tying point in a local element
       xxi_di(1, 1)  = -1.0D0
       xxi_di(2, 1)  =  1.0D0
       xxi_di(3, 1)  =  1.0D0
       xxi_di(4, 1)  = -1.0D0
       xxi_di(5, 1)  =  1.0D0
       xxi_di(6, 1)  = -1.0D0
       ! Eta-coordinate at a tying point in a local element
       eeta_di(1, 1) = -1.0D0
       eeta_di(2, 1) = -1.0D0
       eeta_di(3, 1) =  1.0D0
       eeta_di(4, 1) =  1.0D0
       eeta_di(5, 1) =  0.0D0
       eeta_di(6, 1) =  0.0D0
       
       ! Xi-coordinate at a tying point in a local element
       xxi_di(1, 2)  = -1.0D0
       xxi_di(2, 2)  =  0.0D0
       xxi_di(3, 2)  =  1.0D0
       xxi_di(4, 2)  =  1.0D0
       xxi_di(5, 2)  =  0.0D0
       xxi_di(6, 2)  = -1.0D0
       ! Eta-coordinate at a tying point in a local element
       eeta_di(1, 2) = -1.0D0
       eeta_di(2, 2) = -1.0D0
       eeta_di(3, 2) = -1.0D0
       eeta_di(4, 2) =  1.0D0
       eeta_di(5, 2) =  1.0D0
       eeta_di(6, 2) =  1.0D0
       
       !--------------------------------------------------------
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       !--------------------------------------------------------
        
       ! xi-coordinate at a tying point in a local element
       tpcoord(1, 1, 1)  = 0.5D0
       tpcoord(2, 1, 1)  = 0.0D0
       tpcoord(3, 1, 1)  = 0.5D0
       ! eta-coordinate at a tying point in a local element
       tpcoord(1, 2, 1) = 0.0D0
       tpcoord(2, 2, 1) = 0.5D0
       tpcoord(3, 2, 1) = 0.5D0
       
       !--------------------------------------------------------
       
      END IF
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at the center point in a local element
      ! eta-coordinate at the center point in a local element
      naturalcoord(1) = 0.0D0
      naturalcoord(2) = 0.0D0
      
      CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
      
      !--------------------------------------------------------------
      
      ! Covariant basis vector
      DO i = 1, 3
       
       g1(i) = 0.0D0
       
       DO na = 1, nn
        
        g1(i) = g1(i)+shapederiv(na, 1) &
                      *elem(i, na)      
        
       END DO
       
      END DO
      
      e_0(1) = g1(1)
      e_0(2) = g1(2)
      e_0(3) = g1(3)
      
      !--------------------------------------------------------------
      
      DO nb = 1, nn
       
       !--------------------------------------------------------
       
       naturalcoord(1) = nncoord(nb, 1)
       naturalcoord(2) = nncoord(nb, 2)
       
       CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
       
       !--------------------------------------------------------
       
       ! Covariant basis vector
       DO i = 1, 3
        
        g1(i) = 0.0D0
        g2(i) = 0.0D0
        
        DO na = 1, nn
         
         g1(i) = g1(i)+shapederiv(na, 1) &
                       *elem(i, na)      
         g2(i) = g2(i)+shapederiv(na, 2) &
                       *elem(i, na)      
         
        END DO
        
       END DO
       
       !--------------------------------------------------------
       
       det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
       det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
       det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)
       
       det_cg3_abs = DSQRT( det_cg3(1)*det_cg3(1)   &
                           +det_cg3(2)*det_cg3(2)   &
                           +det_cg3(3)*det_cg3(3) ) 
       
       v3(1, nb) = det_cg3(1)/det_cg3_abs
       v3(2, nb) = det_cg3(2)/det_cg3_abs
       v3(3, nb) = det_cg3(3)/det_cg3_abs
       
       !--------------------------------------------------------
       
       v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
       v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
       v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)
       
       v2_abs = DSQRT( v2(1, nb)*v2(1, nb)   &
                      +v2(2, nb)*v2(2, nb)   &
                      +v2(3, nb)*v2(3, nb) ) 
       
       IF( v2_abs .GT. 1.0D-15 ) THEN
        
        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs
        
        v1(1, nb) = v2(2, nb)*v3(3, nb) &
                   -v2(3, nb)*v3(2, nb) 
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
                   -v2(1, nb)*v3(3, nb) 
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
                   -v2(2, nb)*v3(1, nb) 
        
        v1_abs = DSQRT( v1(1, nb)*v1(1, nb)   &
                       +v1(2, nb)*v1(2, nb)   &
                       +v1(3, nb)*v1(3, nb) ) 
        
        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs
        
       ELSE
        
        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0
        
        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0
        
       END IF
       
       !--------------------------------------------------------
       
       v3(1, nb) = v1(2, nb)*v2(3, nb) &
                  -v1(3, nb)*v2(2, nb) 
       v3(2, nb) = v1(3, nb)*v2(1, nb) &
                  -v1(1, nb)*v2(3, nb) 
       v3(3, nb) = v1(1, nb)*v2(2, nb) &
                  -v1(2, nb)*v2(1, nb) 
       
       v3_abs = DSQRT( v3(1, nb)*v3(1, nb)   &
                      +v3(2, nb)*v3(2, nb)   &
                      +v3(3, nb)*v3(3, nb) ) 
       
       v3(1, nb) = v3(1, nb)/v3_abs
       v3(2, nb) = v3(2, nb)/v3_abs
       v3(3, nb) = v3(3, nb)/v3_abs
       
       !--------------------------------------------------------
       
       a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
       a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
       a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)
       
       !--------------------------------------------------------
       
       a_over_2_theta_cross_v3(1, nb)    &
       = theta(2, nb)*a_over_2_v3(3, nb) &
        -theta(3, nb)*a_over_2_v3(2, nb) 
       a_over_2_theta_cross_v3(2, nb)    &
       = theta(3, nb)*a_over_2_v3(1, nb) &
        -theta(1, nb)*a_over_2_v3(3, nb) 
       a_over_2_theta_cross_v3(3, nb)    &
       = theta(1, nb)*a_over_2_v3(2, nb) &
        -theta(2, nb)*a_over_2_v3(1, nb) 
       
       !--------------------------------------------------------
       
      END DO
           
!--------------------------------------------------------------------
!				Modified stress in laminated shell
!--------------------------------------------------------------------
		
!		write(*,*) 'Stress_n_total_layer',n_total_layer
		   
!--------------------------------------------------------------------
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN 
       
       zeta_ly = 0.0D0
       
      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       zeta_ly = zeta
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       zeta_ly = 0.0D0
       
      END IF
      
      !---------------------------------------------------------
      
      DO it = 1, ntying
       
       DO ip = 1, npoints_tying(it)
        
        !-------------------------------------------------
        
        naturalcoord(1) = tpcoord(ip, 1, it)
        naturalcoord(2) = tpcoord(ip, 2, it)
        
        CALL getShapeFunc(fetype, naturalcoord, shapefunc)
        
        CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
        
        !-------------------------------------------------
        
        DO na = 1, nn
         
         DO i = 1, 3
          
          u_rot(i, na)                      &
          = shapefunc(na)                   &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          
          dudxi_rot(i, na)                  &
          = shapederiv(na, 1)               &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          dudeta_rot(i, na)                 &
          = shapederiv(na, 2)               &
            *( zeta_ly*a_over_2_v3(i, na) ) 
          dudzeta_rot(i, na)                &
          = shapefunc(na)                   &
            *( a_over_2_v3(i, na) )         
          
         END DO
         
        END DO
        
        !-------------------------------------------------
        
        ! Covariant basis vector
        DO i = 1, 3
         
         g1(i) = 0.0D0
         g2(i) = 0.0D0
         g3(i) = 0.0D0
         
         DO na = 1, nn
          
          g1(i) = g1(i)+shapederiv(na, 1)  &
                        *elem(i, na)       &
                      +dudxi_rot(i, na)    
          g2(i) = g2(i)+shapederiv(na, 2)  &
                        *elem(i, na)       &
                       +dudeta_rot(i, na)  
          g3(i) = g3(i)+dudzeta_rot(i, na) 
          
         END DO
         
        END DO

        !---------------------------------------------
        
        DO i = 1, 3
         
         dudxi(i)   = 0.0D0
         dudeta(i)  = 0.0D0
         dudzeta(i) = 0.0D0
         
         DO na = 1, nn
          
          dudxi(i)                                      &
          = dudxi(i)                                    &
           +shapederiv(na, 1)                           &
            *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) ) 
          dudeta(i)                                     &
          = dudeta(i)                                   &
           +shapederiv(na, 2)                           &
            *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) ) 
          dudzeta(i)                                    &
          = dudzeta(i)                                  &
           +shapefunc(na)                               &
            *( a_over_2_theta_cross_v3(i, na) )         
          
         END DO
         
        END DO
        
        !---------------------------------------------
        
        ! Infinitesimal strain tensor
        e11_di(ip, it)                             &
        = 0.5D0                                    &
          *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
            +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
            +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) ) 
        e22_di(ip, it)                             &
        = 0.5D0                                    &
          *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
            +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
            +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) ) 
        e12_2                                   &
        = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
         +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
         +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) ) 
        e23_di_2(ip, it)                        &
        = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
         +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
         +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) ) 
        e31_di_2(ip, it)                        &
        = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
         +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
         +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) ) 
        
        !-------------------------------------------------
        
       END DO
       
      END DO
      
      !--------------------------------------------------------
      
	  	IF ( n_total_layer.ge. 2 ) THEN
			shellmatl = int(gausses(1)%pMaterial%variables(M_SHELL_MATLTYPE))
			if (shellmatl == 0)then
				sigma_layer = 0.0D0
				DO m = 1, n_layer
					sigma_layer = sigma_layer + 2 * gausses(1)%pMaterial%variables(100+3*m) / thick
				END DO
				zeta_ly = -1 + sigma_layer - gausses(1)%pMaterial%variables(100+3*n_layer) / thick * (1-zeta)
			elseif (shellmatl == 1)then
				sigma_layer = 0.0D0
				DO m = 1, n_layer
					sigma_layer = sigma_layer + 2 * gausses(1)%pMaterial%variables(100+8*m-5) / thick
				END DO
				zeta_ly = -1 + sigma_layer - gausses(1)%pMaterial%variables(100+8*n_layer-5) / thick * (1-zeta)			   
			else
				write(*,*)"ERROR : shellmatl isnot correct"; stop
			endif
		ELSE
			zeta_ly = zeta
		ENDIF

      !--------------------------------------------------------
      
      DO lx = 1, nn
       
       !--------------------------------------------------
       
       naturalcoord(1) = nncoord(lx, 1) 
       naturalcoord(2) = nncoord(lx, 2) 
       
       xi_lx  = naturalcoord(1)
       eta_lx = naturalcoord(2)
       
       CALL getShapeFunc(fetype, naturalcoord, shapefunc)
       
       CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
       
       !--------------------------------------------------
       
       DO na = 1, nn
        
        DO i = 1, 3
         
         u_rot(i, na)                      &
         = shapefunc(na)                   &
           *( zeta_ly*a_over_2_v3(i, na) ) 
         
         dudxi_rot(i, na)                  &
         = shapederiv(na, 1)               &
           *( zeta_ly*a_over_2_v3(i, na) ) 
         dudeta_rot(i, na)                 &
        = shapederiv(na, 2)                &
           *( zeta_ly*a_over_2_v3(i, na) ) 
         dudzeta_rot(i, na)                &
         = shapefunc(na)                   &
           *( a_over_2_v3(i, na) )         
         
        END DO
        
       END DO
       
       !--------------------------------------------------
       
       ! Covariant basis vector
       DO i = 1, 3
        
        g1(i) = 0.0D0
        g2(i) = 0.0D0
        g3(i) = 0.0D0
        
        DO na = 1, nn
         
         g1(i) = g1(i)+shapederiv(na, 1)  &
                       *elem(i, na)       &
                      +dudxi_rot(i, na)   
         g2(i) = g2(i)+shapederiv(na, 2)  &
                       *elem(i, na)       &
                      +dudeta_rot(i, na)  
         g3(i) = g3(i)+dudzeta_rot(i, na) 
        
        END DO
       END DO

       !--------------------------------------------------
        
       ! Jacobian
       det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
            +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
            +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) ) 
       
		if(det == 0.0d0) then
			write(*,*)"ERROR:LIB Shell in l2009 Not Jacobian"
			stop
		endif
		
       det_inv = 1.0D0/det

       !--------------------------------------------------
       
       ! Contravariant basis vector
       cg1(1) = det_inv                      &
                *( g2(2)*g3(3)-g2(3)*g3(2) ) 
       cg1(2) = det_inv                      &
                *( g2(3)*g3(1)-g2(1)*g3(3) ) 
       cg1(3) = det_inv                      &
                *( g2(1)*g3(2)-g2(2)*g3(1) ) 
       cg2(1) = det_inv                      &
                *( g3(2)*g1(3)-g3(3)*g1(2) ) 
       cg2(2) = det_inv                      &
                *( g3(3)*g1(1)-g3(1)*g1(3) ) 
       cg2(3) = det_inv                      &
                *( g3(1)*g1(2)-g3(2)*g1(1) ) 
       cg3(1) = det_inv                      &
                *( g1(2)*g2(3)-g1(3)*g2(2) ) 
       cg3(2) = det_inv                      &
                *( g1(3)*g2(1)-g1(1)*g2(3) ) 
       cg3(3) = det_inv                      &
                *( g1(1)*g2(2)-g1(2)*g2(1) ) 

       !--------------------------------------------------
       
       g3_abs = DSQRT( g3(1)*g3(1)   &
                      +g3(2)*g3(2)   &
                      +g3(3)*g3(3) ) 
       
       !--------------------------------------------------
       
       ! Orthonormal vectors
       
       e3_hat(1) = g3(1)/g3_abs
       e3_hat(2) = g3(2)/g3_abs
       e3_hat(3) = g3(3)/g3_abs
       
       e1_hat(1) = g2(2)*e3_hat(3) &
                  -g2(3)*e3_hat(2) 
       e1_hat(2) = g2(3)*e3_hat(1) &
                  -g2(1)*e3_hat(3) 
       e1_hat(3) = g2(1)*e3_hat(2) &
                  -g2(2)*e3_hat(1) 
       e1_hat_abs = DSQRT( e1_hat(1)*e1_hat(1)   &
                          +e1_hat(2)*e1_hat(2)   &
                          +e1_hat(3)*e1_hat(3) ) 
       e1_hat(1) = e1_hat(1)/e1_hat_abs
       e1_hat(2) = e1_hat(2)/e1_hat_abs
       e1_hat(3) = e1_hat(3)/e1_hat_abs
       
       e2_hat(1) = e3_hat(2)*e1_hat(3) &
                  -e3_hat(3)*e1_hat(2) 
       e2_hat(2) = e3_hat(3)*e1_hat(1) &
                  -e3_hat(1)*e1_hat(3) 
       e2_hat(3) = e3_hat(1)*e1_hat(2) &
                  -e3_hat(2)*e1_hat(1)  
       e2_hat_abs = DSQRT( e2_hat(1)*e2_hat(1)   &
                          +e2_hat(2)*e2_hat(2)   &
                          +e2_hat(3)*e2_hat(3) ) 
       e2_hat(1) = e2_hat(1)/e2_hat_abs
       e2_hat(2) = e2_hat(2)/e2_hat_abs
       e2_hat(3) = e2_hat(3)/e2_hat_abs
       
       !--------------------------------------------------
       
       DO i = 1, 3
        
        dudxi(i)   = 0.0D0
        dudeta(i)  = 0.0D0
        dudzeta(i) = 0.0D0
        
        DO na = 1, nn
         
         dudxi(i)                                      &
         = dudxi(i)                                    &
          +shapederiv(na, 1)                           &
           *( edisp(i, na)                             &
             +zeta_ly*a_over_2_theta_cross_v3(i, na) ) 
         dudeta(i)                                     &
         = dudeta(i)                                   &
          +shapederiv(na, 2)                           &
           *( edisp(i, na)                             &
             +zeta_ly*a_over_2_theta_cross_v3(i, na) ) 
         dudzeta(i)                                    &
         = dudzeta(i)                                  &
          +shapefunc(na)                               &
           *( a_over_2_theta_cross_v3(i, na) )         
         
        END DO
        
       END DO
       
       !--------------------------------------------------
       
       ! Infinitesimal strain tensor
       e11                                        &
       = 0.5D0                                    &
         *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
           +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
           +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) ) 
       e22                                        &
       = 0.5D0                                    &
         *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
           +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
           +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) ) 
       e12_2                                   &
       = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
        +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
        +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) ) 
       e23_2                                   &
       = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
        +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
        +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) ) 
       e31_2                                   &
       = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
        +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
        +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) ) 
       
       ! MITC4
       IF( etype .EQ. fe_mitc4_shell ) THEN
        
        e23_2 = 0.0D0
        e31_2 = 0.0D0
        
        ! e23_as_2
        e23_2                                   &
        = 0.5D0*( 1.0D0-xi_lx  )*e23_di_2(4, 1) &
         +0.5D0*( 1.0D0+xi_lx  )*e23_di_2(2, 1) 
        ! e31_as_2
        e31_2                                   &
        = 0.5D0*( 1.0D0-eta_lx )*e31_di_2(1, 1) &
         +0.5D0*( 1.0D0+eta_lx )*e31_di_2(3, 1) 
        
       ! MITC9
       ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
        
        xxi_lx  = xi_lx /DSQRT( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/DSQRT( 3.0D0/5.0D0 )
        
        DO ip = 1, npoints_tying(1)
         
         h(ip, 1)                                     &
         = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
           *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
              *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
             +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
              *( 1.0D0-eeta_lx*eeta_lx ) )            
         
        END DO
        
        xxi_lx  = xi_lx /DSQRT( 3.0D0/5.0D0 )
        eeta_lx = eta_lx/DSQRT( 1.0D0/3.0D0 )
        
        DO ip = 1, npoints_tying(2)
         
         h(ip, 2)                                      &
         = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
             *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
            +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
             *( 1.0D0-xxi_lx*xxi_lx ) )                &
           *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) ) 
         
        END DO
        
        xxi_lx  = xi_lx /DSQRT( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/DSQRT( 1.0D0/3.0D0 )
        
        DO ip = 1, npoints_tying(3)
         
         h(ip, 3)                                      &
         = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
           *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) ) 
         
        END DO
        
        e11   = 0.0D0
        e31_2 = 0.0D0
        
        ! e11_as, e31_as_2
        DO ip = 1, npoints_tying(1)
         
         e11   = e11  +h(ip, 1)*e11_di(ip, 1)
         e31_2 = e31_2+h(ip, 1)*e31_di_2(ip, 1)
         
        END DO
        
        e22   = 0.0D0
        e23_2 = 0.0D0
        
        ! e22_as, e23_as_2
        DO ip = 1, npoints_tying(2)
         
         e22   = e22  +h(ip, 2)*e22_di(ip, 2)
         e23_2 = e23_2+h(ip, 2)*e23_di_2(ip, 2)
         
        END DO
        
        e12_2 = 0.0D0
        
        ! e12_as_2
        DO ip = 1, npoints_tying(3)
         
         e12_2 = e12_2+h(ip, 3)*e12_di_2(ip, 3)
         
        END DO
        
       ! MITC3
       ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
        
        e23_2 = 0.0D0
        e31_2 = 0.0D0
        
        ! e23_as_2
        e23_2                                       &
         = ( 1.0D0-xi_lx )*e23_di_2(2, 1)           &
          +xi_lx *e31_di_2(1, 1)                    &
          +xi_lx *( e23_di_2(3, 1)-e31_di_2(3, 1) ) 
        ! e31_as_2
        e31_2                                       &
         = eta_lx*e23_di_2(2, 1)                    &
          +( 1.0D0-eta_lx )*e31_di_2(1, 1)          &
          -eta_lx*( e23_di_2(3, 1)-e31_di_2(3, 1) ) 
        
       END IF
       
       !--------------------------------------------------
       
       ! { E } vector
       Ev(1) = e11
       Ev(2) = e22
       Ev(3) = e12_2
       Ev(4) = e23_2
       Ev(5) = e31_2
       
       ! Infinitesimal strain tensor
       ! [ E ] matrix
       E(1, 1) = Ev(1)
       E(2, 2) = Ev(2)
       E(3, 3) = 0.0D0
       E(1, 2) = 0.5D0*Ev(3)
       E(2, 1) = 0.5D0*Ev(3)
       E(2, 3) = 0.5D0*Ev(4)
       E(3, 2) = 0.5D0*Ev(4)
       E(3, 1) = 0.5D0*Ev(5)
       E(1, 3) = 0.5D0*Ev(5)
       
       !--------------------------------------------------
    !   write(*,*) 'Stress_n_layer', n_layer
       CALL MatlMatrix_Shell                        &
            (gausses(lx), Shell, D,                 &
             e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
             alpha, n_layer)                                 
       
       !--------------------------------------------------

       Sv = MATMUL( D, Ev )
       
       ! Infinitesimal stress tensor
       ! [ S ] matrix
       S(1, 1) = Sv(1)
       S(2, 2) = Sv(2)
       S(3, 3) = 0.0D0
       S(1, 2) = Sv(3)
       S(2, 1) = Sv(3)
       S(2, 3) = Sv(4)
       S(3, 2) = Sv(4)
       S(3, 1) = Sv(5)
       S(1, 3) = Sv(5)

       !------------------------------------------------
       
       stress(lx,1)                  &
        = ( S(1, 1)*g1(1)*g1(1)   &
           +S(1, 2)*g1(1)*g2(1)   &
           +S(1, 3)*g1(1)*g3(1)   &
           +S(2, 1)*g2(1)*g1(1)   &
           +S(2, 2)*g2(1)*g2(1)   &
           +S(2, 3)*g2(1)*g3(1)   &
           +S(3, 1)*g3(1)*g1(1)   &
           +S(3, 2)*g3(1)*g2(1) ) 
        stress(lx,2)                 &
        = ( S(1, 1)*g1(2)*g1(2)   &
           +S(1, 2)*g1(2)*g2(2)   &
           +S(1, 3)*g1(2)*g3(2)   &
           +S(2, 1)*g2(2)*g1(2)   &
           +S(2, 2)*g2(2)*g2(2)   &
           +S(2, 3)*g2(2)*g3(2)   &
           +S(3, 1)*g3(2)*g1(2)   &
           +S(3, 2)*g3(2)*g2(2) ) 
        stress(lx,3)                 &
        = ( S(1, 1)*g1(3)*g1(3)   &
           +S(1, 2)*g1(3)*g2(3)   &
           +S(1, 3)*g1(3)*g3(3)   &
           +S(2, 1)*g2(3)*g1(3)   &
           +S(2, 2)*g2(3)*g2(3)   &
           +S(2, 3)*g2(3)*g3(3)   &
           +S(3, 1)*g3(3)*g1(3)   &
           +S(3, 2)*g3(3)*g2(3) ) 
        stress(lx,4)                 &
        = ( S(1, 1)*g1(1)*g1(2)   &
           +S(1, 2)*g1(1)*g2(2)   &
           +S(1, 3)*g1(1)*g3(2)   &
           +S(2, 1)*g2(1)*g1(2)   &
           +S(2, 2)*g2(1)*g2(2)   &
           +S(2, 3)*g2(1)*g3(2)   &
           +S(3, 1)*g3(1)*g1(2)   &
           +S(3, 2)*g3(1)*g2(2) ) 
        stress(lx,5)                 &
        = ( S(1, 1)*g1(2)*g1(3)   &
           +S(1, 2)*g1(2)*g2(3)   &
           +S(1, 3)*g1(2)*g3(3)   &
           +S(2, 1)*g2(2)*g1(3)   &
           +S(2, 2)*g2(2)*g2(3)   &
           +S(2, 3)*g2(2)*g3(3)   &
           +S(3, 1)*g3(2)*g1(3)   &
           +S(3, 2)*g3(2)*g2(3) ) 
        stress(lx,6)                 &
        = ( S(1, 1)*g1(3)*g1(1)   &
           +S(1, 2)*g1(3)*g2(1)   &
           +S(1, 3)*g1(3)*g3(1)   &
           +S(2, 1)*g2(3)*g1(1)   &
           +S(2, 2)*g2(3)*g2(1)   &
           +S(2, 3)*g2(3)*g3(1)   &
           +S(3, 1)*g3(3)*g1(1)   &
           +S(3, 2)*g3(3)*g2(1) ) 
        
        strain(lx,1)                   &
        = ( E(1, 1)*cg1(1)*cg1(1)   &
           +E(1, 2)*cg1(1)*cg2(1)   &
           +E(1, 3)*cg1(1)*cg3(1)   &
           +E(2, 1)*cg2(1)*cg1(1)   &
           +E(2, 2)*cg2(1)*cg2(1)   &
           +E(2, 3)*cg2(1)*cg3(1)   &
           +E(3, 1)*cg3(1)*cg1(1)   &
           +E(3, 2)*cg3(1)*cg2(1) ) 
        strain(lx,2)                   &
        = ( E(1, 1)*cg1(2)*cg1(2)   &
           +E(1, 2)*cg1(2)*cg2(2)   &
           +E(1, 3)*cg1(2)*cg3(2)   &
           +E(2, 1)*cg2(2)*cg1(2)   &
           +E(2, 2)*cg2(2)*cg2(2)   &
           +E(2, 3)*cg2(2)*cg3(2)   &
           +E(3, 1)*cg3(2)*cg1(2)   &
           +E(3, 2)*cg3(2)*cg2(2) ) 
        strain(lx,3)                   &
        = ( E(1, 1)*cg1(3)*cg1(3)   &
           +E(1, 2)*cg1(3)*cg2(3)   &
           +E(1, 3)*cg1(3)*cg3(3)   &
           +E(2, 1)*cg2(3)*cg1(3)   &
           +E(2, 2)*cg2(3)*cg2(3)   &
           +E(2, 3)*cg2(3)*cg3(3)   &
           +E(3, 1)*cg3(3)*cg1(3)   &
           +E(3, 2)*cg3(3)*cg2(3) ) 
        strain(lx,4)                   &
        = ( E(1, 1)*cg1(1)*cg1(2)   &
           +E(1, 2)*cg1(1)*cg2(2)   &
           +E(1, 3)*cg1(1)*cg3(2)   &
           +E(2, 1)*cg2(1)*cg1(2)   &
           +E(2, 2)*cg2(1)*cg2(2)   &
           +E(2, 3)*cg2(1)*cg3(2)   &
           +E(3, 1)*cg3(1)*cg1(2)   &
           +E(3, 2)*cg3(1)*cg2(2) ) 
        strain(lx,5)                   &
        = ( E(1, 1)*cg1(2)*cg1(3)   &
           +E(1, 2)*cg1(2)*cg2(3)   &
           +E(1, 3)*cg1(2)*cg3(3)   &
           +E(2, 1)*cg2(2)*cg1(3)   &
           +E(2, 2)*cg2(2)*cg2(3)   &
           +E(2, 3)*cg2(2)*cg3(3)   &
           +E(3, 1)*cg3(2)*cg1(3)   &
           +E(3, 2)*cg3(2)*cg2(3) ) 
        strain(lx,6)                   &
        = ( E(1, 1)*cg1(3)*cg1(1)   &
           +E(1, 2)*cg1(3)*cg2(1)   &
           +E(1, 3)*cg1(3)*cg3(1)   &
           +E(2, 1)*cg2(3)*cg1(1)   &
           +E(2, 2)*cg2(3)*cg2(1)   &
           +E(2, 3)*cg2(3)*cg3(1)   &
           +E(3, 1)*cg3(3)*cg1(1)   &
           +E(3, 2)*cg3(3)*cg2(1) ) 
		   
       !--------------------------------------------------
      
      END DO
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE ElementStress_Shell_MITC
!####################################################################
      
      
!####################################################################
      SUBROUTINE DL_Shell                                  &
                 (etype, nn, ndof, xx, yy, zz, rho, thick, &
                  ltype, params, vect, nsize)              
!####################################################################
      
      USE hecmw
      USE m_utilities
      USE mMechGauss
      USE gauss_integration
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype
      INTEGER(KIND = kint), INTENT(IN) :: nn
      INTEGER(KIND = kint), INTENT(IN) :: ndof
      REAL(KIND = kreal), INTENT(IN) :: xx(*), yy(*), zz(*)
      REAL(KIND = kreal), INTENT(IN) :: rho
      REAL(KIND = kreal), INTENT(IN) :: thick
      REAL(KIND = kreal), INTENT(IN) :: params(*)
      REAL(KIND = kreal), INTENT(OUT) :: vect(*)
      INTEGER(KIND = kint), INTENT(OUT) :: nsize
      
!--------------------------------------------------------------------
      
      INTEGER :: ivol, isurf
      INTEGER :: lx, ly
      INTEGER :: fetype
      INTEGER :: ny
      INTEGER :: i
      INTEGER :: na, nb
      INTEGER :: isize
      INTEGER :: jsize1, jsize2, jsize3, &
                 jsize4, jsize5, jsize6  
      INTEGER :: ltype
      
      REAL(KIND = kreal) :: elem(3, nn)
      REAL(KIND = kreal) :: val
      REAL(KIND = kreal) :: ax, ay, az
      REAL(KIND = kreal) :: rx, ry, rz
      REAL(KIND = kreal) :: xi_lx, eta_lx, zeta_ly
      REAL(KIND = kreal) :: w_w_lx, w_ly
      REAL(KIND = kreal) :: naturalcoord(2)
      REAL(KIND = kreal) :: nncoord(nn, 2)
      REAL(KIND = kreal) :: shapefunc(nn)
      REAL(KIND = kreal) :: shapederiv(nn, 2)
      REAL(KIND = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
      REAL(KIND = kreal) :: v1_abs, v2_abs, v3_abs
      REAL(KIND = kreal) :: a_over_2_v3(3, nn)
      REAL(KIND = kreal) :: u_rot(3, nn)
      REAL(KIND = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
                            dudzeta_rot(3, nn)                   
      REAL(KIND = kreal) :: g1(3), g2(3), g3(3)
      REAL(KIND = kreal) :: g1_abs, g2_abs, g3_abs
      REAL(KIND = kreal) :: g1_cross_g2(3)
      REAL(KIND = kreal) :: e_0(3)
      REAL(KIND = kreal) :: cg1(3), cg2(3), cg3(3)
      REAL(KIND = kreal) :: det
      REAL(KIND = kreal) :: det_cg3(3)
      REAL(KIND = kreal) :: det_inv
      REAL(KIND = kreal) :: det_cg3_abs
      REAL(KIND = kreal) :: w_w_w_det
      REAL(KIND = kreal) :: e3_hat(3)
      REAL(KIND = kreal) :: N(3, ndof*nn)
      REAL(KIND = kreal) :: hx, hy, hz
      REAL(KIND = kreal) :: phx, phy, phz
      REAL(KIND = kreal) :: coefx, coefy, coefz
      REAL(KIND = kreal) :: x, y, z
      
!--------------------------------------------------------------------
      
      !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
      !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
      !   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
      !   CRAV LTYPE=4  :GRAVITY FORCE
      !   CENT LTYPE=5  :CENTRIFUGAL LOAD
      !   P    LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR SHELL SURFACE
      
!--------------------------------------------------------------------
      
      nsize = ndof*nn
      
!--------------------------------------------------------------------
      
      val = params(1)
      ax  = params(2)
      ay  = params(3)
      az  = params(4)
      rx  = params(5)
      ry  = params(6)
      rz  = params(7)
      
!--------------------------------------------------------------------
      
      ! MITC4
      IF( etype .EQ. fe_mitc4_shell ) THEN
       
       fetype = fe_mitc4_shell
       
       ny = 2
       
      ! MITC9
      ELSE IF( etype .EQ. fe_mitc9_shell ) THEN
       
       fetype = fe_mitc9_shell
       
       ny = 3
       
      ! MITC3
      ELSE IF( etype .EQ. fe_mitc3_shell ) THEN
       
       fetype = fe_mitc3_shell
       
       ny = 2
       
      END IF
      
!--------------------------------------------------------------------
      
      DO na = 1, nn
       
       elem(1, na) = xx(na)
       elem(2, na) = yy(na)
       elem(3, na) = zz(na)
       
      END DO
      
!-------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      ! eta-coordinate at a node in a local element
      CALL getNodalNaturalCoord(fetype, nncoord)
      
!--------------------------------------------------------------------
      
      ! Local load vector
      DO isize = 1, ndof*nn
       
       vect(isize) = 0.0D0
       
      END DO
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at the center point in a local element
      ! eta-coordinate at the center point in a local element
      naturalcoord(1) = 0.0D0
      naturalcoord(2) = 0.0D0
      
      CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
      
      !--------------------------------------------------------------
      
      ! Covariant basis vector
      DO i = 1, 3
       
       g1(i) = 0.0D0
       
       DO na = 1, nn
        
        g1(i) = g1(i)+shapederiv(na, 1) &
                      *elem(i, na)      
        
       END DO
       
      END DO
      
      e_0(1) = g1(1)
      e_0(2) = g1(2)
      e_0(3) = g1(3)
      
      !--------------------------------------------------------------
      
      DO nb = 1, nn
       
       !--------------------------------------------------------
       
       naturalcoord(1) = nncoord(nb, 1)
       naturalcoord(2) = nncoord(nb, 2)
       
       CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
       
       !--------------------------------------------------------
       
       ! Covariant basis vector
       DO i = 1, 3
        
        g1(i) = 0.0D0
        g2(i) = 0.0D0
        
        DO na = 1, nn
         
         g1(i) = g1(i)+shapederiv(na, 1) &
                       *elem(i, na)      
         g2(i) = g2(i)+shapederiv(na, 2) &
                       *elem(i, na)      
         
        END DO
        
       END DO
       
       !--------------------------------------------------------
       
       det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
       det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
       det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)
       
       det_cg3_abs = DSQRT( det_cg3(1)*det_cg3(1)   &
                           +det_cg3(2)*det_cg3(2)   &
                           +det_cg3(3)*det_cg3(3) ) 
       
       v3(1, nb) = det_cg3(1)/det_cg3_abs
       v3(2, nb) = det_cg3(2)/det_cg3_abs
       v3(3, nb) = det_cg3(3)/det_cg3_abs
       
       !--------------------------------------------------------
       
       v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
       v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
       v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)
       
       v2_abs = DSQRT( v2(1, nb)*v2(1, nb)   &
                      +v2(2, nb)*v2(2, nb)   &
                      +v2(3, nb)*v2(3, nb) ) 
       
       IF( v2_abs .GT. 1.0D-15 ) THEN
        
        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs
        
        v1(1, nb) = v2(2, nb)*v3(3, nb) &
                   -v2(3, nb)*v3(2, nb) 
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
                   -v2(1, nb)*v3(3, nb) 
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
                   -v2(2, nb)*v3(1, nb) 
        
        v1_abs = DSQRT( v1(1, nb)*v1(1, nb)   &
                       +v1(2, nb)*v1(2, nb)   &
                       +v1(3, nb)*v1(3, nb) ) 
        
        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs
        
       ELSE
        
        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0
        
        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0
        
       END IF
       
       !--------------------------------------------------------
       
       v3(1, nb) = v1(2, nb)*v2(3, nb) &
                  -v1(3, nb)*v2(2, nb) 
       v3(2, nb) = v1(3, nb)*v2(1, nb) &
                  -v1(1, nb)*v2(3, nb) 
       v3(3, nb) = v1(1, nb)*v2(2, nb) &
                  -v1(2, nb)*v2(1, nb) 
       
       v3_abs = DSQRT( v3(1, nb)*v3(1, nb)   &
                      +v3(2, nb)*v3(2, nb)   &
                      +v3(3, nb)*v3(3, nb) ) 
       
       v3(1, nb) = v3(1, nb)/v3_abs
       v3(2, nb) = v3(2, nb)/v3_abs
       v3(3, nb) = v3(3, nb)/v3_abs
       
       !--------------------------------------------------------
       
       a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
       a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
       a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)
       
       !--------------------------------------------------------
       
      END DO
      
!--------------------------------------------------------------------
      
      ! Selction of load type
      
      ivol = 0
      isurf = 0
      
      IF( ltype .LT. 10 ) THEN
       
       ivol = 1
       
      ELSE IF( ltype .GE. 10 ) THEN
       
       isurf = 1
       
      END IF
      
!--------------------------------------------------------------------
      
      !** Surface load
      IF( isurf .EQ. 1 ) THEN
       
       !--------------------------------------------------------
       
       DO lx = 1, NumOfQuadPoints(fetype)
        
        !--------------------------------------------------
        
        CALL getQuadPoint(fetype, lx, naturalcoord)
        
        xi_lx  = naturalcoord(1)
        eta_lx = naturalcoord(2)
        
        w_w_lx = getWeight(fetype, lx)
        
        CALL getShapeFunc(fetype, naturalcoord, shapefunc)
        
        CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
        
        !--------------------------------------------------
        
        DO na = 1, nn
         
         DO i = 1, 3
          
          u_rot(i, na)                    &
          = shapefunc(na)                 &
            *( 0.0D0*a_over_2_v3(i, na) ) 
          
          dudxi_rot(i, na)                &
          = shapederiv(na, 1)             &
            *( 0.0D0*a_over_2_v3(i, na) ) 
          dudeta_rot(i, na)               &
          = shapederiv(na, 2)             &
            *( 0.0D0*a_over_2_v3(i, na) ) 
          dudzeta_rot(i, na)              &
          = shapefunc(na)                 &
            *( a_over_2_v3(i, na) )         
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        ! Covariant basis vector
        DO i = 1, 3
         
         g1(i) = 0.0D0
         g2(i) = 0.0D0
         !g3(i) = 0.0D0
         
         DO na = 1, nn
          
          g1(i) = g1(i)+shapederiv(na, 1)  &
                        *elem(i, na)       &
                       +dudxi_rot(i, na)   
          g2(i) = g2(i)+shapederiv(na, 2)  &
                        *elem(i, na)       &
                       +dudeta_rot(i, na)  
          !g3(i) = g3(i)+dudzeta_rot(i, na)
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        !g3_abs = DSQRT( g3(1)*g3(1)   &
        !               +g3(2)*g3(2)   &
        !               +g3(3)*g3(3) ) 
        
        !--------------------------------------------------
        
        !e3_hat(1) = g3(1)/g3_abs
        !e3_hat(2) = g3(2)/g3_abs
        !e3_hat(3) = g3(3)/g3_abs
        
        !--------------------------------------------------
        
        ! Jacobian
        !det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
        !     +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
        !     +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) ) 
        
        !--------------------------------------------------
        
        g1_cross_g2(1) = g1(2)*g2(3)-g1(3)*g2(2)
        g1_cross_g2(2) = g1(3)*g2(1)-g1(1)*g2(3)
        g1_cross_g2(3) = g1(1)*g2(2)-g1(2)*g2(1)
        
        !--------------------------------------------------
        
        DO nb = 1, nn
         
         jsize1 = ndof*(nb-1)+1
         jsize2 = ndof*(nb-1)+2
         jsize3 = ndof*(nb-1)+3
         jsize4 = ndof*(nb-1)+4
         jsize5 = ndof*(nb-1)+5
         jsize6 = ndof*(nb-1)+6
         
         N(1, jsize1) = shapefunc(nb)
         N(1, jsize2) = 0.0D0
         N(1, jsize3) = 0.0D0
         N(1, jsize4) = 0.0D0
         N(1, jsize5) = 0.0D0
         N(1, jsize6) = 0.0D0
         N(2, jsize1) = 0.0D0
         N(2, jsize2) = shapefunc(nb)
         N(2, jsize3) = 0.0D0
         N(2, jsize4) = 0.0D0
         N(2, jsize5) = 0.0D0
         N(2, jsize6) = 0.0D0
         N(3, jsize1) = 0.0D0
         N(3, jsize2) = 0.0D0
         N(3, jsize3) = shapefunc(nb)
         N(3, jsize4) = 0.0D0
         N(3, jsize5) = 0.0D0
         N(3, jsize6) = 0.0D0
         
        END DO
                
        DO isize = 1, ndof*nn
         
         vect(isize)                                    &
         = vect(isize)                                  &
          +w_w_lx*( N(1, isize)*g1_cross_g2(1)       &
                      +N(2, isize)*g1_cross_g2(2)       &
                      +N(3, isize)*g1_cross_g2(3) )*val 
         
        END DO
        
        !--------------------------------------------------
        
       END DO
       
       !--------------------------------------------------------
       
      END IF
      
!--------------------------------------------------------------------
      
      !** Volume load
      IF( ivol .EQ. 1 ) THEN
       
       !--------------------------------------------------------
       
       DO ly = 1, ny
        
        !--------------------------------------------------
        
        zeta_ly = xg(ny, ly)
        w_ly    = wgt(ny, ly)
        
        !--------------------------------------------------
        
        DO lx = 1, NumOfQuadPoints(fetype)
         
         !--------------------------------------------
         
         CALL getQuadPoint(fetype, lx, naturalcoord)
         
         xi_lx  = naturalcoord(1)
         eta_lx = naturalcoord(2)
         
         w_w_lx = getWeight(fetype, lx)
         
         CALL getShapeFunc(fetype, naturalcoord, shapefunc)
         
         CALL getShapeDeriv(fetype, naturalcoord, shapederiv)
         
         !--------------------------------------------
         
         DO na = 1, nn
          
          DO i = 1, 3
           
           u_rot(i, na)                      &
           = shapefunc(na)                   &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           
           dudxi_rot(i, na)                  &
           = shapederiv(na, 1)               &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           dudeta_rot(i, na)                 &
           = shapederiv(na, 2)               &
             *( zeta_ly*a_over_2_v3(i, na) ) 
           dudzeta_rot(i, na)                &
           = shapefunc(na)                   &
             *( a_over_2_v3(i, na) )         
           
          END DO
          
         END DO
         
         !--------------------------------------------
         
         ! Covariant basis vector
         DO i = 1, 3
          
          g1(i) = 0.0D0
          g2(i) = 0.0D0
          g3(i) = 0.0D0
          
          DO na = 1, nn
           
           g1(i) = g1(i)+shapederiv(na, 1)  &
                         *elem(i, na)       &
                        +dudxi_rot(i, na)   
           g2(i) = g2(i)+shapederiv(na, 2)  &
                         *elem(i, na)       &
                        +dudeta_rot(i, na)  
           g3(i) = g3(i)+dudzeta_rot(i, na) 
           
          END DO
          
         END DO
         
         !--------------------------------------------
         
         ! Jacobian
         det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
              +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
              +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) ) 
         
         !--------------------------------------------
         
         ! [ N ] matrix
         DO nb = 1, nn
          
          jsize1 = ndof*(nb-1)+1
          jsize2 = ndof*(nb-1)+2
          jsize3 = ndof*(nb-1)+3
          jsize4 = ndof*(nb-1)+4
          jsize5 = ndof*(nb-1)+5
          jsize6 = ndof*(nb-1)+6
          
          N(1, jsize1) =  shapefunc(nb)
          N(2, jsize1) =  0.0D0
          N(3, jsize1) =  0.0D0
          N(1, jsize2) =  0.0D0
          N(2, jsize2) =  shapefunc(nb)
          N(3, jsize2) =  0.0D0
          N(1, jsize3) =  0.0D0
          N(2, jsize3) =  0.0D0
          N(3, jsize3) =  shapefunc(nb)
          N(1, jsize4) =  0.0D0
          N(2, jsize4) = -u_rot(3, nb)
          N(3, jsize4) =  u_rot(2, nb)
          N(1, jsize5) =  u_rot(3, nb)
          N(2, jsize5) =  0.0D0
          N(3, jsize5) = -u_rot(1, nb)
          N(1, jsize6) = -u_rot(2, nb)
          N(2, jsize6) =  u_rot(1, nb)
          N(3, jsize6) =  0.0D0
          
         ENDDO
         
         !--------------------------------------------
         
         w_w_w_det = w_w_lx*w_ly*det
         
         !--------------------------------------------
         
         IF( ltype .EQ. 1 ) THEN
          
          DO isize = 1, ndof*nn
           
           vect(isize) = vect(isize)+w_w_w_det*N(1, isize)*val
           
          END DO
          
         ELSE IF( ltype .EQ. 2 ) THEN
          
          DO isize = 1, ndof*nn
           
           vect(isize) = vect(isize)+w_w_w_det*N(2, isize)*val
           
          END DO
          
         ELSE IF( ltype .EQ. 3 ) THEN
           
          DO isize = 1, ndof*nn
           
           vect(isize) = vect(isize)+w_w_w_det*N(3, isize)*val
           
          END DO
          
         ELSE IF( ltype .EQ. 4 ) THEN
          
          DO isize = 1, ndof*nn
           
           vect(isize) = vect(isize)+w_w_w_det*rho*ax*N(1, isize)*val
           vect(isize) = vect(isize)+w_w_w_det*rho*ay*N(2, isize)*val
           vect(isize) = vect(isize)+w_w_w_det*rho*az*N(3, isize)*val
           
          END DO
          
         ELSE IF( ltype .EQ. 5 ) THEN
          
          x = 0.0D0
          y = 0.0D0
          z = 0.0D0
          
          DO nb = 1, nn
           
           x = x+shapefunc(nb)*elem(1, nb)
           y = y+shapefunc(nb)*elem(2, nb)
           z = z+shapefunc(nb)*elem(3, nb)
           
          END DO
          
          hx = ax+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rx
          hy = ay+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*ry
          hz = az+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rz
          
          phx = x-hx
          phy = y-hy
          phz = z-hz
          
          coefx = phx*val*rho*val
          coefy = phy*val*rho*val
          coefz = phz*val*rho*val
          
          DO isize = 1, ndof*nn
           
           vect(isize)                       &
           = vect(isize)                     &
            +w_w_w_det*( N(1, isize)*coefx   &
                        +N(2, isize)*coefy   &
                        +N(3, isize)*coefz ) 
           
          END DO
          
         END IF
         
         !--------------------------------------------
         
        END DO
        
        !--------------------------------------------------
        
       END DO
       
       !--------------------------------------------------------
       
      END IF
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE DL_Shell
!####################################################################
      
!####################################################################
      END MODULE m_static_LIB_shell
