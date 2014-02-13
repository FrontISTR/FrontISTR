!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> This modules just summarizes all modules used in static analysis
module m_static_LIB

use m_utilities
use m_solve_LINEQ

use m_static_LIB_1d
use m_static_LIB_2d
use m_static_LIB_3d
use m_static_LIB_C3D8
use m_static_LIB_3dIC

use m_static_LIB_beam
use m_static_LIB_shell

end module m_static_LIB
