!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.0                                   !
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
   
end module