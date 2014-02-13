!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by X. YUAN (AdvanceSoft)                  !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!
!> \brief  This module provide common functions of beam elements
!
!>  \author                date              version 
!>  X.Yuan(Advancesoft)    2012/11/16        original
!>  
!======================================================================!
module m_static_LIB_beam

use hecmw
USE mMechGauss
use m_utilities
use elementinfo

implicit none

INTEGER, PARAMETER, private :: ndof=6

contains

  !< This subroutine build up coordinate frame for beam elements   
  subroutine framtr(refx, xl, le, t)
      real(kind=kreal), intent(in)  :: refx(3)      !< Reference Vector
      real(kind=kreal), intent(in)  :: xl(3,2)      !< Element coordinates
      real(kind=kreal), intent(out) :: le           !< length of the element
      real(kind=kreal), intent(out) :: t(3,3)       !< Transformation array

      real(kind=kreal) ::       dl,theta

      real(kind=kreal), parameter ::  tol = 1.d-08 

      t(1,1) = xl(1,2) - xl(1,1)
      t(1,2) = xl(2,2) - xl(2,1)
      t(1,3) = xl(3,2) - xl(3,1)
      le  = sqrt(t(1,1)*t(1,1)+t(1,2)*t(1,2)+t(1,3)*t(1,3))
      dl  = 1.0d0/le
      t(1,1) = t(1,1)*dl
      t(1,2) = t(1,2)*dl
      t(1,3) = t(1,3)*dl

      t(3,1) = refx(1)
      t(3,2) = refx(2)
      t(3,3) = refx(3)


      t(2,1) = (t(3,2)*t(1,3) - t(3,3)*t(1,2))
      t(2,2) = (t(3,3)*t(1,1) - t(3,1)*t(1,3))
      t(2,3) = (t(3,1)*t(1,2) - t(3,2)*t(1,1))
      dl  = sqrt(t(2,1)*t(2,1)+t(2,2)*t(2,2)+t(2,3)*t(2,3))
      if(dl<tol*le) then
        stop  "Bad reference for beam element!"
      else
        t(2,1) = t(2,1)/dl
        t(2,2) = t(2,2)/dl
        t(2,3) = t(2,3)/dl
        t(3,1) = t(1,2)*t(2,3) - t(1,3)*t(2,2)
        t(3,2) = t(1,3)*t(2,1) - t(1,1)*t(2,3)
        t(3,3) = t(1,1)*t(2,2) - t(1,2)*t(2,1)
      endif

  end subroutine

  !> Calculate stiff matrix of BEAM elements
   SUBROUTINE STF_Beam(etype,nn,ecoord,section,E,P,STIFF)
      INTEGER, INTENT(IN)            :: etype        !< element's type
      INTEGER, INTENT(IN)            :: nn           !< number of element's nodes
      REAL(kind=kreal), INTENT(in)   :: ecoord(3,nn) !< coordinates of elemental nodes
      REAL(kind=kreal), INTENT(in)   :: section(:)   !< section parameters
      REAL(kind=kreal), INTENT(in)   :: E,P   !< status of qudrature points
      REAL(kind=kreal), INTENT(out)  :: STIFF(nn*6,nn*6)   !< elemental stiff matrix
	  
      REAL(kind=kreal) :: le, outa(2), trans(3,3), refv(3), transt(3,3)
      REAL(kind=kreal) :: G 
      REAL(kind=kreal) :: L2, L3, A, Iy, Iz, Jx, EA, twoE, fourE, twelveE, sixE
      logical  :: ierr

      refv = section(1:3)	  
      call framtr(refv, ecoord, le, trans)
      transT= transpose(trans)
      L2 = le*le
      L3 = L2*le
	  
      G = E/(2.d0*(1.d0 + P)) 

      A = section(4);  Iy=section(5); Iz=section(6); Jx=section(7)

      EA = E*A/le
      twoE = 2.d0*E/le
      fourE = 4.d0*E/le
      twelveE = 12*E/L3
      sixE = 6*E/L2

       stiff = 0.d0
       stiff(1,1) = EA;  
       stiff(7,1) = -EA;

       stiff(2,2) = twelveE*Iz;
       stiff(6,2) = sixE*Iz;
       stiff(8,2) = -twelveE*Iz;
       stiff(12,2) = sixE*Iz;
	    
       stiff(3,3) = twelveE*Iy;
       stiff(5,3) = -sixE*Iy;
       stiff(9,3) = -twelveE*Iy;
       stiff(11,3) = -sixE*Iy;	    
	    
       stiff(4,4) = G*Jx/le;
       stiff(10,4) = -G*Jx/le;
	    
       stiff(3,5) = -sixE*Iy;
       stiff(5,5) = fourE*Iy;
       stiff(9,5) = sixE*Iy;
       stiff(11,5) = twoE*Iy;
	    
       stiff(2,6) = sixE*Iz;
       stiff(6,6) = fourE*Iz;
       stiff(8,6) = -sixE*Iz;
       stiff(12,6) = twoE*Iz;	    
	    
       stiff(1,7) = -EA;  
       stiff(7,7) = EA;

       stiff(2,8) = -twelveE*Iz;
       stiff(6,8) = -sixE*Iz;
       stiff(8,8) = twelveE*Iz;
       stiff(12,8) = -sixE*Iz;
	    
       stiff(3,9) = -twelveE*Iy;
       stiff(5,9) = sixE*Iy;
       stiff(9,9) = twelveE*Iy;
       stiff(11,9) = sixE*Iy;	    
	    
       stiff(4,10) = -G*Jx/le;
       stiff(10,10) = G*Jx/le;
	    
       stiff(3,11) = -sixE*Iy;
       stiff(5,11) = twoE*Iy;
       stiff(9,11) = sixE*Iy;
       stiff(11,11) = fourE*Iy;
	    
       stiff(2,12) = sixE*Iz;
       stiff(6,12) = twoE*Iz;
       stiff(8,12) = -sixE*Iz;
       stiff(12,12) = fourE*Iz;
	   
      stiff(1:3,:) = matmul( transT, stiff(1:3,:) )
      stiff(4:6,:) = matmul( transT, stiff(4:6,:) )
      stiff(7:9,:) = matmul( transT, stiff(7:9,:) )
      stiff(10:12,:) = matmul( transT, stiff(10:12,:) )
	  
	  stiff(:,1:3) = matmul( stiff(:,1:3), trans )
      stiff(:,4:6) = matmul( stiff(:,4:6), trans )
      stiff(:,7:9) = matmul( stiff(:,7:9), trans )
      stiff(:,10:12) = matmul( stiff(:,10:12), trans )

   END SUBROUTINE
   
   
      ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) <
      !> Calculate stiff matrix of BEAM elements
!####################################################################
      SUBROUTINE STF_Beam_641                                         &
                 (etype, nn, ecoord, gausses, section, stiff, tt, t0) 
!####################################################################
      
      USE mMechGauss
      
!--------------------------------------------------------------------
      
      INTEGER, INTENT(IN)            :: etype              !< element type
      INTEGER, INTENT(IN)            :: nn                 !< the total number of elemental nodes
      REAL(kind=kreal), INTENT(IN)   :: ecoord(3, nn)      !< coordinates of elemental nodes
      TYPE(tGaussStatus), INTENT(IN) :: gausses(:)         !< status of Gaussian qudrature points
      REAL(kind=kreal), INTENT(IN)   :: section(:)         !< section parameters
      REAL(kind=kreal), INTENT(OUT)  :: stiff(nn*3, nn*3)  !< elemental stiffness matrix
      REAL(kind=kreal), INTENT(IN), OPTIONAL :: tt(nn), t0(nn)
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: refv(3)
      REAL(KIND = kreal) :: trans(3, 3), transt(3, 3)
      REAL(kind = kreal) :: ec(3, 2)
      REAL(KIND = kreal) :: tempc
      REAL(KIND = kreal) :: ina1(1), outa1(2)
      REAL(KIND = kreal) :: ee, pp
      REAL(KIND = kreal) :: le
      REAL(KIND = kreal) :: l2, l3, g, a, iy, iz, jx
      REAL(KIND = kreal) :: ea, twoe, foure, twelvee, sixe
      
      LOGICAL :: ierr
      
!--------------------------------------------------------------------
      
      refv(1) = section(1)
      refv(2) = section(2)
      refv(3) = section(3)
      
      ec(1, 1) = ecoord(1, 1)
      ec(2, 1) = ecoord(2, 1)
      ec(3, 1) = ecoord(3, 1)
      ec(1, 2) = ecoord(1, 2)
      ec(2, 2) = ecoord(2, 2)
      ec(3, 2) = ecoord(3, 2)
      
      CALL framtr(refv, ec, le, trans)
      
      transt= TRANSPOSE( trans )
      
      l2 = le*le
      l3 = l2*le
      
!--------------------------------------------------------------------
      
      IF( PRESENT( tt ) ) THEN
       
       tempc = 0.5D0*( tt(1)+tt(2) )
       
      END IF
      
!--------------------------------------------------------------------
      
      IF( PRESENT( tt ) ) THEN
       
       ina1(1) = tempc
       
       CALL fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa1, ierr, ina1 )
       
      ELSE
       
       ierr = .TRUE.
       
      END IF
      
      !--------------------------------------------------------------
      
      IF( ierr ) THEN
       
       ee = gausses(1)%pMaterial%variables(M_YOUNGS)
       pp = gausses(1)%pMaterial%variables(M_POISSON)
       
      ELSE
       
       ee = outa1(1)
       pp = outa1(2)
       
      END IF
      
!--------------------------------------------------------------------
      
      g = ee/( 2.0D0*( 1.0D0+pp ) )
      
      a = section(4)
      
      iy = section(5)
      iz = section(6)
      jx = section(7)
      
!--------------------------------------------------------------------
      
      ea = ee*a/le
      
      twoe    = 2.0D0*ee/le
      foure   = 4.0D0*ee/le
      twelvee = 12.0D0*ee/l3
      sixe    = 6.0D0*ee/l2
      
!--------------------------------------------------------------------
      
      stiff = 0.0D0
      
      stiff(1, 1) = ea
      !stiff(7, 1) = -ea
      stiff(4, 1) = -ea
      
      stiff(2, 2)  = twelvee*iz
      !stiff(6, 2) = sixe*iz
      stiff(9, 2)  = sixe*iz
      !stiff(8, 2) = -twelvee*iz
      stiff(5, 2)  = -twelvee*iz
      stiff(12, 2) = sixe*iz
      
      stiff(3, 3)  = twelvee*iy
      !stiff(5, 3) = -sixe*iy
      stiff(8, 3)  = -sixe*iy
      !stiff(9, 3) = -twelvee*iy
      stiff(6, 3)  = -twelvee*iy
      stiff(11, 3) = -sixe*iy
      
      !stiff(4, 4) = g*jx/le
      stiff(7, 7)  = g*jx/le
      !stiff(10, 4) = -g*jx/le
      stiff(10, 7) = -g*jx/le
      
      !stiff(3, 5) = -sixe*iy
      stiff(3, 8)  = -sixe*iy
      !stiff(5, 5) = foure*iy
      stiff(8, 8)  = foure*iy
      !stiff(9, 5) = sixe*iy
      stiff(6, 8)  = sixe*iy
      !stiff(11, 5) = twoe*iy
      stiff(11, 8) = twoe*iy
       
      !stiff(2, 6) = sixe*iz
      stiff(2, 9)  = sixe*iz
      !stiff(6, 6) = foure*iz
      stiff(9, 9)  = foure*iz
      !stiff(8, 6) = -sixe*iz
      stiff(5, 9)  = -sixe*iz
      !stiff(12, 6) = twoe*iz
      stiff(12, 9) = twoe*iz
      
      !stiff(1, 7) = -ea
      stiff(1, 4) = -ea
      !stiff(7, 7) = ea
      stiff(4, 4) = ea
      
      !stiff(2, 8) = -twelvee*iz
      stiff(2, 5)  = -twelvee*iz
      !stiff(6, 8) = -sixe*iz
      stiff(9, 5)  = -sixe*iz
      !stiff(8, 8) = twelvee*iz
      stiff(5, 5)  = twelvee*iz
      !stiff(12, 8) = -sixe*iz
      stiff(12, 5) = -sixe*iz
      
      !stiff(3, 9) = -twelvee*iy
      stiff(3, 6)  = -twelvee*iy
      !stiff(5, 9) = sixe*iy
      stiff(8, 6)  = sixe*iy
      !stiff(9, 9) = twelvee*iy
      stiff(6, 6)  = twelvee*iy
      !stiff(11, 9) = sixe*iy
      stiff(11, 6) = sixe*iy
      
      !stiff(4, 10) = -g*jx/le
      stiff(7, 10)  = -g*jx/le
      stiff(10, 10) = g*jx/le
      
      stiff(3, 11) = -sixe*iy
      !stiff(5, 11) = twoe*iy
      stiff(8, 11) = twoe*iy
      !stiff(9, 11) = sixe*iy
      stiff(6, 11) = sixe*iy
      stiff(11, 11) = foure*iy
      
      stiff(2, 12)  = sixe*iz
      !stiff(6, 12) = twoe*iz
      stiff(9, 12)  = twoe*iz
      !stiff(8, 12) = -sixe*iz
      stiff(5, 12)  = -sixe*iz
      stiff(12, 12) = foure*iz
      
!--------------------------------------------------------------------
      
      stiff(  1:3, :) = MATMUL( transt, stiff(  1:3, :) )
      stiff(  4:6, :) = MATMUL( transt, stiff(  4:6, :) )
      stiff(  7:9, :) = MATMUL( transt, stiff(  7:9, :) )
      stiff(10:12, :) = MATMUL( transt, stiff(10:12, :) )
      
      stiff(:,   1:3) = MATMUL( stiff(:,   1:3), trans )
      stiff(:,   4:6) = MATMUL( stiff(:,   4:6), trans )
      stiff(:,   7:9) = MATMUL( stiff(:,   7:9), trans )
      stiff(:, 10:12) = MATMUL( stiff(:, 10:12), trans )
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE STF_Beam_641
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2014/02/06)
      
      
      ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) < 
!####################################################################
      SUBROUTINE DL_Beam_641(etype, nn, xx, yy, zz, rho, ltype, params, &
                             section, vect, nsize)                      
!####################################################################
!**
!**  SET DLOAD
!** 
!   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
!   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
!   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
!   GRAV LTYPE=4  :GRAVITY FORCE
!   CENT LTYPE=5  :CENTRIFUGAL LOAD
!   P1   LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR FACE-1
!   P2   LTYPE=20 :TRACTION IN NORMAL-DIRECTION FOR FACE-2
!   P3   LTYPE=30 :TRACTION IN NORMAL-DIRECTION FOR FACE-3
!   P4   LTYPE=40 :TRACTION IN NORMAL-DIRECTION FOR FACE-4
!   P5   LTYPE=50 :TRACTION IN NORMAL-DIRECTION FOR FACE-5
!   P6   LTYPE=60 :TRACTION IN NORMAL-DIRECTION FOR FACE-6
! I/F VARIABLES
      INTEGER(KIND = kint), INTENT(IN)  :: etype, nn
      REAL(KIND = kreal), INTENT(IN)    :: xx(:), yy(:), zz(:)
      REAL(KIND = kreal), INTENT(IN)    :: params(0:6)
      REAL(KIND = kreal), INTENT(IN)    :: section(:)
      REAL(KIND = kreal), INTENT(INOUT) :: vect(:)
      REAL(KIND = kreal) :: rho
      INTEGER(KIND = kint) :: ltype, nsize
! LOCAL VARIABLES
      INTEGER(KIND = kint) :: ndof
      PARAMETER(ndof = 3)
      REAL(KIND = kreal) :: h(nn)
      REAL(KIND = kreal) :: xj(3, 3), det, wg
      INTEGER(KIND = kint) :: ivol, isuf
      INTEGER(KIND = kint) :: nod(nn)
      INTEGER(KIND = kint) :: ig2, lx, i ,surtype, nsur
      REAL(KIND = kreal) :: vx, vy, vz, xcod, ycod, zcod
      REAL(KIND = kreal) :: ax, ay, az, rx, ry, rz, hx, hy, hz, val
      REAL(KIND = kreal) :: phx, phy, phz
      REAL(KIND = kreal) :: coefx, coefy, coefz
      REAL(KIND = kreal) :: normal(3), localcoord(3), elecoord(3, nn), deriv(nn, 3)
      REAL(KIND = kreal) :: a, aa
      
!--------------------------------------------------------------------
      
      val = params(0)
      
      !--------------------------------------------------------------
      
      ivol = 0
      isuf = 0
      
      IF( ltype .LT. 10 ) THEN
       
       ivol = 1
       
      ELSE IF( ltype .GE. 10 ) THEN
       
       isuf = 1
       
       CALL getSubFace(etype, ltype/10, surtype, nod)
       
       nsur = getNumberOfNodes(surtype)
       
      END IF
      
!--------------------------------------------------------------------
      
      nsize = nn*ndof
      
!--------------------------------------------------------------------
      
      vect(1:nsize) = 0.0D0
      
      !--------------------------------------------------------------
      
      ! Volume force
      
      IF( ivol .EQ. 1 ) THEN
       
       IF( ltype .EQ. 4 ) THEN
        
        AA = DSQRT( ( xx(2)-xx(1) )*( xx(2)-xx(1) )   &
                   +( yy(2)-yy(1) )*( yy(2)-yy(1) )   &
                   +( zz(2)-zz(1) )*( zz(2)-zz(1) ) ) 
        
        a = section(4)
        
        vx = params(1)
        vy = params(2)
        vz = params(3)
        vx = vx/DSQRT( params(1)**2+params(2)**2+params(3)**2 )
        vy = vy/DSQRT( params(1)**2+params(2)**2+params(3)**2 )
        vz = vz/DSQRT( params(1)**2+params(2)**2+params(3)**2 )
        
        DO i = 1, 2
         
         vect(3*i-2) = val*rho*a*0.5D0*AA*vx
         vect(3*i-1) = val*rho*a*0.5D0*AA*vy
         vect(3*i  ) = val*rho*a*0.5D0*AA*vz
         
        END DO
        
        DO i = 3, 4
         
         vect(3*i-2) = 0.0D0
         vect(3*i-1) = 0.0D0
         vect(3*i  ) = 0.0D0
         
        END DO
        
       END IF
       
      END IF
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE DL_Beam_641
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2014/02/06)
      
      
      ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) <
!####################################################################
      SUBROUTINE TLOAD_Beam_641                        &
                 (etype, nn, ndof, xx, yy, zz, tt, t0, &
                  gausses, section, vect)              
!####################################################################
      
      USE hecmw
      USE m_fstr
      USE m_utilities
      USE mMechGauss
      USE gauss_integration
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN) :: etype
      INTEGER(KIND = kint), INTENT(IN) :: nn
      INTEGER(KIND = kint), INTENT(IN) :: ndof
      TYPE(tGaussStatus), INTENT(IN)   :: gausses(:)
      REAL(KIND = kreal), INTENT(IN)   :: section(:)
      REAL(KIND = kreal), INTENT(IN)   :: xx(nn), yy(nn), zz(nn)
      REAL(KIND = kreal), INTENT(IN)   :: tt(nn), t0(nn)
      REAL(KIND = kreal), INTENT(OUT)  :: vect(nn*ndof)
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: tempc, temp0
      REAL(KIND = kreal) :: ecoord(3, nn)
      REAL(KIND = kreal) :: ec(3, 2)
      REAL(KIND = kreal) :: ina1(1), outa1(2)
      REAL(KIND = kreal) :: ina2(1), outa2(1)
      REAL(KIND = kreal) :: alp, alp0
      REAL(KIND = kreal) :: ee, pp
      REAL(KIND = kreal) :: a
      REAL(KIND = kreal) :: refv(3)
      REAL(KIND = kreal) :: alpha_bar
      REAL(KIND = kreal) :: g
      REAL(KIND = kreal) :: le
      REAL(KIND = kreal) :: trans(3, 3), transt(3, 3)
      
      LOGICAL :: ierr
      
!--------------------------------------------------------------------
      
      ecoord(1, 1:nn) = xx(1:nn)
      ecoord(2, 1:nn) = yy(1:nn)
      ecoord(3, 1:nn) = zz(1:nn)
      
!--------------------------------------------------------------------
      
      tempc = 0.5D0*( tt(1)+tt(2) )
      temp0 = 0.5D0*( t0(1)+t0(2) )
      
      !--------------------------------------------------------------
      
      ina1(1) = tempc
      
      CALL fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa1, ierr, ina1 )
      
      IF( ierr ) THEN
       
       ee = gausses(1)%pMaterial%variables(M_YOUNGS)
       pp = gausses(1)%pMaterial%variables(M_POISSON)
      
      ELSE
       
       ee = outa1(1)
       pp = outa1(2)
       
      END IF
      
      !--------------------------------------------------------------
      
      ina2(1) = tempc
      
      CALL fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )
       
      IF( ierr ) STOP "Fails in fetching expansion coefficient!"
       
      alp = outa2(1)
      
      !--------------------------------------------------------------
      
      ina2(1) = temp0
      
      CALL fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )
      
      IF( ierr ) STOP "Fails in fetching expansion coefficient!"
      
      alp0 = outa2(1)
      
!--------------------------------------------------------------------
      
      refv(1) = section(1)
      refv(2) = section(2)
      refv(3) = section(3)
      
      ec(1, 1) = ecoord(1, 1)
      ec(2, 1) = ecoord(2, 1)
      ec(3, 1) = ecoord(3, 1)
      ec(1, 2) = ecoord(1, 2)
      ec(2, 2) = ecoord(2, 2)
      ec(3, 2) = ecoord(3, 2)
      
      CALL framtr(refv, ec, le, trans)
      
      transt= TRANSPOSE( trans )
      
!--------------------------------------------------------------------
      
      a = section(4)
      
      g = ee/( 2.0D0*( 1.0D0+pp ))
      
!--------------------------------------------------------------------
      
      vect( 1) = -a*ee*( alp*( tempc-REF_TEMP )-alp0*( temp0-REF_TEMP ) )
      vect( 2) =  0.0D0
      vect( 3) =  0.0D0
      
      vect( 4) =  a*ee*( alp*( tempc-REF_TEMP )-alp0*( temp0-REF_TEMP ) )
      vect( 5) =  0.0D0
      vect( 6) =  0.0D0
      
      vect( 7) =  0.0D0
      vect( 8) =  0.0D0
      vect( 9) =  0.0D0
      
      vect(10) =  0.0D0
      vect(11) =  0.0D0
      vect(12) =  0.0D0
      
!--------------------------------------------------------------------
      
      vect(  1:3) = MATMUL( transt,   vect(1:3) )
      vect(  4:6) = MATMUL( transt,   vect(4:6) )
      vect(  7:9) = MATMUL( transt,   vect(7:9) )
      vect(10:12) = MATMUL( transt, vect(10:12) )
	  
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE TLOAD_Beam_641
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2013/09/13)
      
      
      ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) <
!####################################################################
      SUBROUTINE NodalStress_Beam_641                         &
                 (etype, nn, ecoord, gausses, section, edisp, &
                  ndstrain, ndstress, tt, t0, ntemp)          
!####################################################################
      
      USE m_fstr
      USE mMechGauss
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint), INTENT(IN)  :: etype
      INTEGER(KIND = kint), INTENT(IN)  :: nn
      REAL(KIND = kreal), INTENT(IN)    :: ecoord(3, nn)
      TYPE(tGaussStatus), INTENT(INOUT) :: gausses(:)
      REAL(KIND = kreal), INTENT(IN)    :: section(:)
      REAL(KIND = kreal), INTENT(IN)    :: edisp(3, nn)
      REAL(KIND = kreal), INTENT(OUT)   :: ndstrain(nn, 6)
      REAL(KIND = kreal), INTENT(OUT)   :: ndstress(nn, 6)
      REAL(kind=kreal), INTENT(IN), OPTIONAL :: tt(nn), t0(nn)
      INTEGER(KIND = kint), INTENT(IN)  :: ntemp
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint) :: i, j, k
      
      REAL(KIND = kreal) :: tempc, temp0
      REAL(KIND = kreal) :: ina1(1), outa1(2)
      REAL(KIND = kreal) :: ina2(1), outa2(1)
      REAL(KIND = kreal) :: alp, alp0
      REAL(KIND = kreal) :: ee, pp
      REAL(KIND = kreal) :: a, radius, angle(6)
      REAL(KIND = kreal) :: refv(3)
      REAL(KIND = kreal) :: le, l2, l3
      REAL(KIND = kreal) :: trans(3, 3), transT(3, 3)
      REAL(KIND = kreal) :: edisp_hat(3, nn)
      REAL(KIND = kreal) :: ec(3, 2)
      REAL(KIND = kreal) :: t(3, 3), t_hat(3, 3)
      REAL(KIND = kreal) :: t_hat_tmp(3, 3)
      REAL(KIND = kreal) :: e(3, 3), e_hat(3, 3)
      REAL(KIND = kreal) :: e_hat_tmp(3, 3)
      REAL(KIND = kreal) :: x1_hat, x2_hat, x3_hat
      REAL(KIND = kreal) :: pi
      
      LOGICAL :: ierr
      
!--------------------------------------------------------------------
      
      pi = 4.0D0*DATAN( 1.0D0 )
      
!--------------------------------------------------------------------
      
      IF( PRESENT( tt ) .AND. PRESENT( t0 ) ) THEN
       
       tempc = 0.5D0*( tt(1)+tt(2) )
       temp0 = 0.5D0*( t0(1)+t0(2) )
       
      END IF
      
!--------------------------------------------------------------------
      
      IF( ntemp .EQ. 1 ) THEN
       
       ina1(1) = tempc
       
       CALL fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa1, ierr, ina1 )
       
      ELSE
       
       ierr = .TRUE.
       
      END IF
      
      !--------------------------------------------------------------
      
      IF( ierr ) THEN
       
       ee = gausses(1)%pMaterial%variables(M_YOUNGS)
       pp = gausses(1)%pMaterial%variables(M_POISSON)
       
      ELSE
       
       ee = outa1(1)
       pp = outa1(2)
       
      END IF
      
!--------------------------------------------------------------------
      
      IF( ntemp .EQ. 1 ) THEN
       
       ina2(1) = tempc
       
       CALL fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )
        
       IF( ierr ) STOP "Fails in fetching expansion coefficient!"
       
       alp = outa2(1)
       
      END IF
      
      !--------------------------------------------------------------
      
      IF( ntemp .EQ. 1 ) THEN
       
       ina2(1) = temp0
       
       CALL fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )
       
       IF( ierr ) STOP "Fails in fetching expansion coefficient!"
       
       alp0 = outa2(1)
       
      END IF
      
!--------------------------------------------------------------------
      
      refv(1) = section(1)
      refv(2) = section(2)
      refv(3) = section(3)
      
      ec(1, 1) = ecoord(1, 1)
      ec(2, 1) = ecoord(2, 1)
      ec(3, 1) = ecoord(3, 1)
      ec(1, 2) = ecoord(1, 2)
      ec(2, 2) = ecoord(2, 2)
      ec(3, 2) = ecoord(3, 2)
      
      CALL framtr(refv, ec, le, trans)
      
      transt= TRANSPOSE( trans )
      
      l2 = le*le
      l3 = l2*le
      
!--------------------------------------------------------------------
      
      a = section(4)
      
      radius = gausses(1)%pMaterial%variables(M_BEAM_RADIUS)
      
      angle(1)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE1)
      angle(2)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE2)
      angle(3)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE3)
      angle(4)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE4)
      angle(5)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE5)
      angle(6)  = gausses(1)%pMaterial%variables(M_BEAM_ANGLE6)
      
!--------------------------------------------------------------------
      
      DO k = 1, 6
       
       !--------------------------------------------------------
       
       angle(k) = angle(k)/180.0D0*pi
       
       x2_hat = radius*DCOS( angle(k) )
       x3_hat = radius*DSIN( angle(k) )
       
       !--------------------------------------------------------
       
       DO j = 1, nn
        
        DO i = 1, 3
         
         edisp_hat(i, j) = trans(i, 1)*edisp(1, j) &
                          +trans(i, 2)*edisp(2, j) &
                          +trans(i, 3)*edisp(3, j) 
         
        END DO
        
       END DO
       
       !--------------------------------------------------------
       
       x1_hat = 0.5D0*le
       
       e_hat = 0.0D0
       e_hat(1, 1) = ( edisp_hat(1, 2)-edisp_hat(1, 1) )/le
       
       t_hat = 0.0D0
       t_hat(1, 1) = ee*( edisp_hat(1, 2)-edisp_hat(1, 1) )/le                   &
                    -ee*x2_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(2, 1)  &
                                +( -4.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(2, 2)  &
                                +( -2.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 4) ) &
                    -ee*x3_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(3, 1)  &
                                +(  4.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(3, 2)  &
                                +(  2.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 4) ) 
       
       IF( PRESENT( tt ) ) THEN
        
        t_hat(1, 1)                                             &
        = t_hat(1, 1)                                           &
         -ee*( alp*( tempc-REF_TEMP )-alp0*( temp0-REF_TEMP ) ) 
        
       END IF
       
       e_hat_tmp(1:3,:) = MATMUL( trans, e_hat(1:3,:) )
       t_hat_tmp(1:3,:) = MATMUL( trans, t_hat(1:3,:) )
       
       e(:, 1:3) = matmul( e_hat_tmp(:,1:3), transt )
       t(:, 1:3) = matmul( t_hat_tmp(:,1:3), transt )
       
       gausses(1)%strain(k) = e_hat(1, 1)
       gausses(1)%stress(k) = t_hat(1, 1)
       
       !--------------------------------------------------------
       
       ndstrain(1, k) = 0.0D0
       ndstrain(2, k) = 0.0D0
       ndstrain(3, k) = 0.0D0
       ndstrain(4, k) = 0.0D0
       
       ndstress(1, k) = 0.0D0
       ndstress(2, k) = 0.0D0
       ndstress(3, k) = 0.0D0
       ndstress(4, k) = 0.0D0
       
       !--------------------------------------------------------
       
       x1_hat = 0.0D0
       
       e_hat = 0.0D0
       e_hat(1, 1) = ( edisp_hat(1, 2)-edisp_hat(1, 1) )/le
       
       t_hat = 0.0D0
       t_hat(1, 1) = ee*( edisp_hat(1, 2)-edisp_hat(1, 1) )/le                   &
                    -ee*x2_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(2, 1)  &
                                +( -4.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(2, 2)  &
                                +( -2.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 4) ) &
                    -ee*x3_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(3, 1)  &
                                +(  4.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(3, 2)  &
                                +(  2.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 4) ) 
       
       IF( PRESENT( tt ) ) THEN
        
        t_hat(1, 1)                                             &
        = t_hat(1, 1)                                           &
         -ee*( alp*( tempc-REF_TEMP )-alp0*( temp0-REF_TEMP ) ) 
        
       END IF
       
       e_hat_tmp(1:3, :) = MATMUL( trans, e_hat(1:3, :) )
       t_hat_tmp(1:3, :) = MATMUL( trans, t_hat(1:3, :) )
       
       e(:, 1:3) = MATMUL( e_hat_tmp(:, 1:3), transt )
       t(:, 1:3) = MATMUL( t_hat_tmp(:, 1:3), transt )
       
       ndstrain(1, k) = e_hat(1, 1)
       ndstress(1, k) = t_hat(1, 1)
       
       !--------------------------------------------------------
       
       x1_hat = le
       
       e_hat = 0.0D0
       e_hat(1, 1) = ( edisp_hat(1, 2)-edisp_hat(1, 1) )/le
       
       t_hat = 0.0D0
       t_hat(1, 1) = ee*( edisp_hat(1, 2)-edisp_hat(1, 1) )/le                   &
                    -ee*x2_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(2, 1)  &
                                +( -4.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(2, 2)  &
                                +( -2.0D0/le+6.0D0*x1_hat/l2 )*edisp_hat(3, 4) ) &
                    -ee*x3_hat*( ( -6.0D0/l2+12.0D0*x1_hat/l3 )*edisp_hat(3, 1)  &
                                +(  4.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 3)   &
                                +(  6.0D0/l2-12.0D0*x1_hat/l3 )*edisp_hat(3, 2)  &
                                +(  2.0D0/le-6.0D0*x1_hat/l2 )*edisp_hat(2, 4) ) 
       
       IF( PRESENT( tt ) ) THEN
        
        t_hat(1, 1)                                             &
        = t_hat(1, 1)                                           &
         -ee*( alp*( tempc-REF_TEMP )-alp0*( temp0-REF_TEMP ) ) 
        
       END IF
       
       e_hat_tmp(1:3, :) = MATMUL( trans, e_hat(1:3, :) )
       t_hat_tmp(1:3, :) = MATMUL( trans, t_hat(1:3, :) )
       
       e(:, 1:3) = MATMUL( e_hat_tmp(:, 1:3), transt )
       t(:, 1:3) = MATMUL( t_hat_tmp(:, 1:3), transt )
       
       ndstrain(2, k) = e_hat(1, 1)
       ndstress(2, k) = t_hat(1, 1)
       
       !--------------------------------------------------------
       
      END DO
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE NodalStress_Beam_641
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2013/09/13)
      
end module
