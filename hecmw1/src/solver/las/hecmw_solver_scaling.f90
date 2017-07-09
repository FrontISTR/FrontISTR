!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_scaling
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_solver_scaling_33
  use hecmw_solver_scaling_44
  use hecmw_solver_scaling_66
  use hecmw_solver_scaling_nn
  implicit none

  private
  real(kind=kreal), private, allocatable :: SCALE(:)

  public :: hecmw_solver_scaling_fw
  public :: hecmw_solver_scaling_bk

contains

  subroutine hecmw_solver_scaling_fw(hecMESH, hecMAT, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(inout) :: COMMtime
    SELECT CASE(hecMAT%NDOF)
      CASE(3)
        call hecmw_solver_scaling_fw_33(hecMESH, hecMAT, COMMtime)
      CASE(4)
        call hecmw_solver_scaling_fw_44(hecMESH, hecMAT, COMMtime)
      CASE(6)
        call hecmw_solver_scaling_fw_66(hecMESH, hecMAT, COMMtime)
      CASE DEFAULT
        call hecmw_solver_scaling_fw_nn(hecMESH, hecMAT, COMMtime)
    END SELECT
  end subroutine hecmw_solver_scaling_fw

  subroutine hecmw_solver_scaling_bk(hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    SELECT CASE(hecMAT%NDOF)
      CASE(3)
        call hecmw_solver_scaling_bk_33(hecMAT)
      CASE(4)
        call hecmw_solver_scaling_bk_44(hecMAT)
      CASE(6)
        call hecmw_solver_scaling_bk_66(hecMAT)
      CASE DEFAULT
        call hecmw_solver_scaling_bk_nn(hecMAT)
    END SELECT
    
  end subroutine hecmw_solver_scaling_bk

end module hecmw_solver_scaling
