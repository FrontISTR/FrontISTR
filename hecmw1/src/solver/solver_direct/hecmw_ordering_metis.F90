!-------------------------------------------------------------------------------
! Copyright (c) 2020 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_ORDERING_METIS is a program for the Nested Dissection
!         ordering provided by Metis
!----------------------------------------------------------------------
module hecmw_ordering_metis
  use hecmw_util
  implicit none

  private
  public :: hecmw_ordering_metis_nodend

contains

  !======================================================================!
  !> @brief hecmw_ordering_metis_NodeND
  !======================================================================!
  subroutine hecmw_ordering_metis_NodeND(Neqns,Xadj,Adj0,Perm,Invp)
    implicit none
    !------
    integer(kind=kint), intent(in):: NEQns
    integer(kind=kint), intent(in):: Adj0(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(out):: Perm(:)
    integer(kind=kint), intent(out):: Invp(:)
    !------
#ifdef HECMW_WITH_METIS

#  if HECMW_METIS_VER == 5

    integer(kind=kint), allocatable:: vwght(:)
    integer(kind=kint):: options(40)
    integer(kind=kint):: ierror

    allocate(vwght(Neqns),stat=ierror)
    if ( ierror/=0 ) stop "ALLOCATION ERROR, vwght: SUB. gennd"

    vwght(:)=1

    call METIS_SetDefaultOptions(options)
    ! set fortran numbering
    options(18)=1

    call METIS_NodeND(Neqns,Xadj,Adj0,vwght,options,Perm,Invp)

    deallocate(vwght)

#  elif HECMW_METIS_VER == 4

    integer(kind=kint):: numflag
    integer(kind=kint):: options(8)

    numflag=1
    options(:)=0

    call METIS_NodeND(Neqns,Xadj,Adj0,numflag,options,Perm,Invp)

#  else
#    error unknown HECMW_METIS_VER
#  endif

#else
    stop "METIS not available"
#endif
  end subroutine hecmw_ordering_metis_NodeND

end module hecmw_ordering_metis
