!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Coupling Interface                                 !
!                                                                      !
!            Written by Shin'ichi Ezure (RIST)                         !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!


module hecmw_couple_define_f

  use hecmw_util

  implicit none
  private

  integer(kind=kint), parameter, public :: HECMW_COUPLE_TYPE_UNDEF  = -100
  integer(kind=kint), parameter, public :: HECMW_COUPLE_TYPE_MXN    =  101
  integer(kind=kint), parameter, public :: HECMW_COUPLE_TYPE_MAXMN  =  102
  integer(kind=kint), parameter, public :: HECMW_COUPLE_TYPE_MANUAL =  103

  integer(kind=kint), parameter, public :: HECMW_COUPLE_UNIT_UNDEF = -200
  integer(kind=kint), parameter, public :: HECMW_COUPLE_UNIT1      =  201
  integer(kind=kint), parameter, public :: HECMW_COUPLE_UNIT2      =  202

  integer(kind=kint), parameter, public :: HECMW_COUPLE_DIRECTION_UNDEF = -300
  integer(kind=kint), parameter, public :: HECMW_COUPLE_UNIT1_TO_UNIT2  =  301
  integer(kind=kint), parameter, public :: HECMW_COUPLE_UNIT2_TO_UNIT1  =  302

  integer(kind=kint), parameter, public :: HECMW_COUPLE_GROUP_UNDEF   = -400
  integer(kind=kint), parameter, public :: HECMW_COUPLE_NODE_GROUP    =  401
  integer(kind=kint), parameter, public :: HECMW_COUPLE_ELEMENT_GROUP =  402
  integer(kind=kint), parameter, public :: HECMW_COUPLE_SURFACE_GROUP =  403

  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_UNDEF        = -500
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_NODE_TO_NODE =  501
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_NODE_TO_ELEM =  502
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_NODE_TO_SURF =  503
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_ELEM_TO_NODE =  511
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_ELEM_TO_ELEM =  512
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_ELEM_TO_SURF =  513
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_SURF_TO_NODE =  521
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_SURF_TO_ELEM =  522
  integer(kind=kint), parameter, public :: HECMW_COUPLE_IP_SURF_TO_SURF =  523

  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_UNDEF        = -600
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_NODE_TO_NODE =  601
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_NODE_TO_ELEM =  602
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_NODE_TO_SURF =  603
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_ELEM_TO_NODE =  611
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_ELEM_TO_ELEM =  612
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_ELEM_TO_SURF =  613
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_SURF_TO_NODE =  621
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_SURF_TO_ELEM =  622
  integer(kind=kint), parameter, public :: HECMW_COUPLE_MAP_SURF_TO_SURF =  623

end module hecmw_couple_define_f
