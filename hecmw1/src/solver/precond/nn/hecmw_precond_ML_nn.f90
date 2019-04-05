!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_precond_ML_nn
!C***
!C
module hecmw_precond_ML_nn
  use  hecmw_util

  private

  public:: hecmw_precond_ML_nn_setup
  public:: hecmw_precond_ML_nn_apply
  public:: hecmw_precond_ML_nn_clear

  integer(kind=kint), save :: id

  logical, save :: INITIALIZED = .false.

contains

  subroutine hecmw_precond_ML_nn_setup(hecMAT, hecMESH, sym)
    use hecmw_matrix_misc
    use hecmw_mat_id
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: sym
    integer(kind=kint) :: ierr
    integer(kind=kint), save :: n_recycle = 0
    if (INITIALIZED) then
      if (hecMAT%Iarray(98) == 1) then      ! need symbolic and numerical setup
        call hecmw_precond_ML_nn_clear()
      else if (hecMAT%Iarray(97) == 1) then ! need numerical setup only
        call hecmw_precond_ML_nn_clear()
      else                                  ! no need to setup or skip setup
        call hecmw_mat_id_clear(id)
        call hecmw_mat_id_set(hecMAT, hecMESH, id)
        return
      endif
    endif
    call hecmw_mat_id_set(hecMAT, hecMESH, id)
    call hecmw_ML_wrapper_setup(id, sym, hecMAT%NDOF, ierr)
    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0 ! symbolic setup done
    hecMAT%Iarray(97) = 0 ! numerical setup done
    n_recycle = 0
  end subroutine hecmw_precond_ML_nn_setup

  subroutine hecmw_precond_ML_nn_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: ierr
    call hecmw_ML_wrapper_apply(id, WW, ierr)
  end subroutine hecmw_precond_ML_nn_apply

  subroutine hecmw_precond_ML_nn_clear()
    use hecmw_mat_id
    implicit none
    integer(kind=kint) :: ierr
    call hecmw_ML_wrapper_clear(id, ierr)
    call hecmw_mat_id_clear(id)
    INITIALIZED = .false.
  end subroutine hecmw_precond_ML_nn_clear

end module     hecmw_precond_ML_nn
