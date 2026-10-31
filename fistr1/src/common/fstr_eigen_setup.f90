!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine to initialize fstr_eigen.
module m_fstr_eigen_setup
  use m_fstr
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

end module m_fstr_eigen_setup
