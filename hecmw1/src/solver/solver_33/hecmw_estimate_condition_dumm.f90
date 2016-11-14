!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_estimate_condition

  public :: hecmw_estimate_cond_num_CG
  public :: hecmw_estimate_cond_num_GMRES

contains

  subroutine hecmw_estimate_condition_CG(ITER, D, E)
    use hecmw_util
    implicit none
    integer(kind=kint), intent(in) :: ITER
    real(kind=kreal), intent(in) :: D(:), E(:)
  end subroutine hecmw_estimate_condition_CG


  subroutine hecmw_estimate_condition_GMRES(I, H)
    use hecmw_util
    implicit none
    integer(kind=kint), intent(in) :: I
    real(kind=kreal), intent(in) :: H(:,:)
  end subroutine hecmw_estimate_condition_GMRES

end module hecmw_estimate_condition
