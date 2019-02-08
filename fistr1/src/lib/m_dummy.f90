!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module defined dummy data and function
module m_dummy
  use hecmw
  use mMechGauss

  implicit none

  type tDummy
    !> DUMMY
    integer(kind=kint) :: DUMMY_egrp_tot
    integer(kind=kint), pointer :: DUMMY_egrp_GRPID  (:)  =>null()
    integer(kind=kint), pointer :: DUMMY_egrp_ID     (:)  =>null()
    integer(kind=kint), pointer :: DUMMY_egrp_amp    (:)  =>null()
    real(kind=kreal), pointer   :: DUMMY_egrp_eps    (:)  =>null()
  end type

  integer, parameter :: kDUM_UNDEFINED = -1
  integer, parameter :: kDUM_INACTIVE  = 0
  integer, parameter :: kDUM_ACTIVE    = 1

contains

  !< print dummy info
  subroutine print_dummy_info( dum )
    type(tDummy), intent(in) :: dum  !< dummy info

    integer(kind=kint) :: i, j

    write(*,'(A,I0)') 'DUMMY_egrp_tot: ', dum%DUMMY_egrp_tot
    do i=1,dum%DUMMY_egrp_tot
      write(*,*) 'DUMMY : ',i
      write(*,*) 'GRPID    : ',dum%DUMMY_egrp_GRPID
      write(*,*) 'EGRPID   : ',dum%DUMMY_egrp_ID
      write(*,*) 'AMPLITUDE: ',dum%DUMMY_egrp_amp
      write(*,*) 'EPSILON  : ',dum%DUMMY_egrp_eps
    end do
  end subroutine

end module m_dummy
