!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize fstr_eigen and to read
!! !EIGEN into it.
module m_fstr_eigen_setup
  use m_fstr
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_ctrl_err_stop
  use fstr_ctrl_eigen
  implicit none

contains

  !> Initial setting of eigen ca;culation
  subroutine fstr_eigen_init( fstrEIG )
    implicit none
    type(fstr_eigen) :: fstrEIG

    fstrEIG%nget        = 5
    fstrEIG%maxiter     = 60
    fstrEIG%iter        = 0
    !fstrEIG%sigma       = 0.0d0
    fstrEIG%sigma       = 0.01d0
    fstrEIG%tolerance   = 1.0d-6
    fstrEIG%totalmass   = 0.0d0
  end subroutine fstr_eigen_init

  subroutine fstr_eigen_finalize( fstrEIG )
    implicit none
    type(fstr_eigen) :: fstrEIG
    if(associated(fstrEIG%eigval)) deallocate(fstrEIG%eigval)
    if(associated(fstrEIG%eigvec)) deallocate(fstrEIG%eigvec)
    if(associated(fstrEIG%resid)) deallocate(fstrEIG%resid)
    if(associated(fstrEIG%filter)) deallocate(fstrEIG%filter)
    if(associated(fstrEIG%mass)) deallocate(fstrEIG%mass)
    if(associated(fstrEIG%effmass)) deallocate(fstrEIG%effmass)
    if(associated(fstrEIG%partfactor)) deallocate(fstrEIG%partfactor)
  end subroutine fstr_eigen_finalize

  !-----------------------------------------------------------------------------!
  !> Read in !EIGEN                                                            !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_EIGEN( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_EIGEN( ctrl, P%EIGEN%nget, P%EIGEN%tolerance, P%EIGEN%maxiter, P%EIGEN%sigma)
    if( rcode /= 0) call fstr_ctrl_err_stop

  end subroutine fstr_setup_EIGEN

end module m_fstr_eigen_setup
