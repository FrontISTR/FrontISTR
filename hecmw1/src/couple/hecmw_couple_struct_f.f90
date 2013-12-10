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


module hecmw_couple_struct_f

  use hecmw_util

  implicit none
  private

  type, public :: hecmw_couple_value
    integer(kind=kint)                        :: n
    integer(kind=kint)                        :: item_type
    integer(kind=kint)                        :: n_dof
    integer(kind=kint), dimension(:), pointer :: item
    real(kind=kreal),   dimension(:), pointer :: value
  end type hecmw_couple_value

end module hecmw_couple_struct_f
