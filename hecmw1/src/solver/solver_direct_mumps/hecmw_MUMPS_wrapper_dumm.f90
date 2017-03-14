!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides wrapper for parallel sparse direct solver MUMPS
module m_hecmw_MUMPS_wrapper
  use hecmw_util
  use m_sparse_matrix

  private
  public :: hecmw_mumps_wrapper

contains

  subroutine hecmw_mumps_wrapper(spMAT, job, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(out) :: istat
    stop "MUMPS not available"
  end subroutine hecmw_mumps_wrapper

end module m_hecmw_MUMPS_wrapper
