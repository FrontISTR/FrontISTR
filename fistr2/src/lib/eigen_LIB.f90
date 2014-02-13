!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This modules just summarizes all modules used in eigen analysis
module m_eigen_LIB

use lczparm
use m_utilities
use m_solve_lineq
use m_eigen_lib_2d1mass
use m_eigen_lib_2d2mass
use m_eigen_lib_3d1mass
use m_eigen_lib_3d2mass

end module m_eigen_LIB
