!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : kind parameters
!!
!! Self-contained kind parameters for the SA-AMG core so that the core modules
!! can be compiled both inside hecmw and stand-alone (gfortran only) for the
!! development harness.  The values mirror hecmw_util exactly (kint=4, kreal=8).
module hecmw_precond_SAAMG_util
  implicit none
  public

  integer(kind=4), parameter :: kint  = 4
  integer(kind=4), parameter :: kreal = 8

end module hecmw_precond_SAAMG_util
