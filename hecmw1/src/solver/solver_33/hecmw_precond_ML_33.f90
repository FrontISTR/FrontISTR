!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/07/04                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!C
!C*** 
!C*** module hecmw_precond_ML_33
!C***
!C
module hecmw_precond_ML_33
  use  hecmw_util

  private

  public:: hecmw_precond_ML_33_setup
  public:: hecmw_precond_ML_33_apply
  public:: hecmw_precond_ML_33_clear

  integer(kind=kint) :: id

contains

  subroutine hecmw_precond_ML_33_setup(hecMAT, hecMESH, sym)
    use hecmw_matrix_misc
    use hecmw_mat_id
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: sym
    integer(kind=kint) :: ierr
    call hecmw_mat_id_set(hecMAT, hecMESH, id)
    call hecmw_ML_wrapper_setup(id, sym, ierr)
  end subroutine hecmw_precond_ML_33_setup

  subroutine hecmw_precond_ML_33_apply(WW)
    implicit none
    real(kind=kreal), intent(inout) :: WW(:)
    integer(kind=kint) :: ierr
    call hecmw_ML_wrapper_apply(id, WW, ierr)
  end subroutine hecmw_precond_ML_33_apply

  subroutine hecmw_precond_ML_33_clear()
    use hecmw_mat_id
    implicit none
    integer(kind=kint) :: ierr
    call hecmw_ML_wrapper_clear(id, ierr)
    call hecmw_mat_id_clear(id)
  end subroutine hecmw_precond_ML_33_clear

end module     hecmw_precond_ML_33
