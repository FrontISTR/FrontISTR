!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                   Written by Gaku Hashimoto, The University of Tokyo !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!                                                                      !
!> \brief  This module contains functions for interpolation in 9 node 
!!   quadrilateral element
!                                                                      !
!>  \author     Gaku Hashimoto, The University of Tokyo
!>  \date       2012/11/15
!>  \version    0.00
!======================================================================!

      MODULE shape_quad9n
!####################################################################
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      INTEGER, PARAMETER, PRIVATE :: kreal = KIND( 0.0D0 )
      
!--------------------------------------------------------------------
      
      CONTAINS
      

!####################################################################
      SUBROUTINE ShapeFunc_quad9n(lcoord, func)
!####################################################################
      
      REAL(KIND = kreal), INTENT(IN)  :: lcoord(2)
      REAL(KIND = kreal), INTENT(OUT) :: func(9)
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: xi(9), eta(9)
      REAL(KIND = kreal) :: n_xi(9), n_eta(9)
      
      INTEGER :: na
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      xi(1)  = -1.0D0
      xi(2)  =  1.0D0
      xi(3)  =  1.0D0
      xi(4)  = -1.0D0
      xi(5)  =  0.0D0
      xi(6)  =  1.0D0
      xi(7)  =  0.0D0
      xi(8)  = -1.0D0
      xi(9)  =  0.0D0
      ! eta-coordinate at a node in a local element
      eta(1) = -1.0D0
      eta(2) = -1.0D0
      eta(3) =  1.0D0
      eta(4) =  1.0D0
      eta(5) = -1.0D0
      eta(6) =  0.0D0
      eta(7) =  1.0D0
      eta(8) =  0.0D0
      eta(9) =  0.0D0
      
!--------------------------------------------------------------------
      
      DO na = 1, 9
       
       n_xi(na)  = ( 0.5D0*xi(na) *lcoord(1) )    &
                   *( 1.0D0+xi(na) *lcoord(1) )   &
                  +( 1.0D0-xi(na) *xi(na)  )      &
                   *( 1.0D0-lcoord(1)*lcoord(1) ) 
       n_eta(na) = ( 0.5D0*eta(na)*lcoord(2) )    &
                   *( 1.0D0+eta(na)*lcoord(2) )   &
                  +( 1.0D0-eta(na)*eta(na) )      &
                   *( 1.0D0-lcoord(2)*lcoord(2) ) 
       
      END DO
      
!--------------------------------------------------------------------
      
      DO na = 1, 9
       
       func(na) = n_xi(na)*n_eta(na)
       
      END DO
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE ShapeFunc_quad9n
!####################################################################

!####################################################################
      SUBROUTINE ShapeDeriv_quad9n(lcoord, func)
!####################################################################
      
      REAL(KIND = kreal), INTENT(IN)  :: lcoord(2)
      REAL(KIND = kreal), INTENT(OUT) :: func(9, 2)
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: xi(9), eta(9)
      REAL(KIND = kreal) :: n_xi(9), n_eta(9)
      REAL(KIND = kreal) :: dn_xi(9), dn_eta(9)
      
      INTEGER :: na
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      xi(1)  = -1.0D0
      xi(2)  =  1.0D0
      xi(3)  =  1.0D0
      xi(4)  = -1.0D0
      xi(5)  =  0.0D0
      xi(6)  =  1.0D0
      xi(7)  =  0.0D0
      xi(8)  = -1.0D0
      xi(9)  =  0.0D0
      ! eta-coordinate at a node in a local element
      eta(1) = -1.0D0
      eta(2) = -1.0D0
      eta(3) =  1.0D0
      eta(4) =  1.0D0
      eta(5) = -1.0D0
      eta(6) =  0.0D0
      eta(7) =  1.0D0
      eta(8) =  0.0D0
      eta(9) =  0.0D0
      
!--------------------------------------------------------------------
      
      DO na = 1, 9
       
       n_xi(na)  = ( 0.5D0*xi(na) *lcoord(1) )    &
                   *( 1.0D0+xi(na) *lcoord(1) )   &
                  +( 1.0D0-xi(na) *xi(na)  )      &
                   *( 1.0D0-lcoord(1)*lcoord(1) ) 
       n_eta(na) = ( 0.5D0*eta(na)*lcoord(2) )    &
                   *( 1.0D0+eta(na)*lcoord(2) )   &
                  +( 1.0D0-eta(na)*eta(na) )      &
                   *( 1.0D0-lcoord(2)*lcoord(2) ) 
       
       dn_xi(na)  = ( 0.5D0*xi(na)  )            &
                    *( 1.0D0+xi(na) *lcoord(1) ) &
                   +( 0.5D0*xi(na) *lcoord(1) )  &
                    *xi(na)                      &
                   +( 1.0D0-xi(na) *xi(na)  )    &
                    *( -2.0D0*lcoord(1) )        
       dn_eta(na) = ( 0.5D0*eta(na) )            &
                    *( 1.0D0+eta(na)*lcoord(2) ) &
                   +( 0.5D0*eta(na)*lcoord(2) )  &
                    *eta(na)                     &
                   +( 1.0D0-eta(na)*eta(na) )    &
                    *( -2.0D0*lcoord(2) )        
       
      END DO
      
!--------------------------------------------------------------------
      
      DO na = 1, 9
       
       func(na, 1) = dn_xi(na)*n_eta(na)
       func(na, 2) = n_xi(na) *dn_eta(na)
       
      END DO
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE ShapeDeriv_quad9n
!####################################################################

!####################################################################
      SUBROUTINE NodalNaturalCoord_quad9n(nncoord)
!####################################################################
      
      REAL(KIND = kreal), INTENT(OUT) :: nncoord(9, 2)
      
!--------------------------------------------------------------------
      
      ! xi-coordinate at a node in a local element
      nncoord(1, 1)  = -1.0D0
      nncoord(2, 1)  =  1.0D0
      nncoord(3, 1)  =  1.0D0
      nncoord(4, 1)  = -1.0D0
      nncoord(5, 1)  =  0.0D0
      nncoord(6, 1)  =  1.0D0
      nncoord(7, 1)  =  0.0D0
      nncoord(8, 1)  = -1.0D0
      nncoord(9, 1)  =  0.0D0
      ! eta-coordinate at a node in a local element
      nncoord(1, 2) = -1.0D0
      nncoord(2, 2) = -1.0D0
      nncoord(3, 2) =  1.0D0
      nncoord(4, 2) =  1.0D0
      nncoord(5, 2) = -1.0D0
      nncoord(6, 2) =  0.0D0
      nncoord(7, 2) =  1.0D0
      nncoord(8, 2) =  0.0D0
      nncoord(9, 2) =  0.0D0
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE NodalNaturalCoord_quad9n
!####################################################################

      END MODULE shape_quad9n
