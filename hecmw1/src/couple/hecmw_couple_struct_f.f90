!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

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
