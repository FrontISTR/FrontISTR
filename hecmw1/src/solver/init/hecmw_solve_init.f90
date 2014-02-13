!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!
!C*** 
!C*** hecmw_solve_init
!C***
!

module m_hecmw_solve_init

contains

      subroutine hecmw_solve_init (hecMAT) 
      use hecmw_util
      implicit REAL*8 (A-H,O-Z)
      type (hecmwST_matrix)     :: hecMAT
      
      hecMAT%Iarray= 0
      hecMAT%Rarray= 0.d0

      hecMAT%Iarray(1)= 100 
      hecMAT%Iarray(6)=  10

      hecMAT%Rarray(1)= 1.d-8
      hecMAT%Rarray(2)= 1.d0
      hecMAT%Rarray(3)= 0.d0
      hecMAT%Rarray(4)= 0.10d0
      hecMAT%Rarray(5)= 0.10d0

      end subroutine hecmw_solve_init

end module m_hecmw_solve_init


