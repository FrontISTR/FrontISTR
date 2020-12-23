!-------------------------------------------------------------------------------
! Copyright (c) 2020 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_ORDERING is a program for fill-reducing ordering
!         for direct solver
!----------------------------------------------------------------------
module hecmw_ordering
  use hecmw_util
  implicit none

  private
  public :: hecmw_ordering_gen

  integer(kind=kint), parameter:: ORDERING_DEFAULT = 0
  integer(kind=kint), parameter:: ORDERING_QMD     = 1
  integer(kind=kint), parameter:: ORDERING_METIS   = 2

contains

  !======================================================================!
  !> @brief hecmw_ordering_gen
  !======================================================================!
  subroutine hecmw_ordering_gen(Neqns,Nttbr,Xadj,Adj0,Perm,Invp,opt)
    use hecmw_ordering_qmd
    use hecmw_ordering_metis
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nttbr
    integer(kind=kint), intent(in):: Adj0(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(in):: opt
    integer(kind=kint), intent(out):: Perm(:)
    integer(kind=kint), intent(out):: Invp(:)
    !------
    integer(kind=kint):: ordering
    ordering = opt
    if (ordering < 0 .or. ordering > 2) then
      stop "ERROR ordering option for direct solver out of range"
    endif
    if (ordering == ORDERING_DEFAULT) then
#ifdef HECMW_WITH_METIS
      ordering = ORDERING_METIS
#else
      ordering = ORDERING_QMD
#endif
    endif
    select case (ordering)
    case(ORDERING_QMD)
      call hecmw_ordering_GENQMD(Neqns,Nttbr,Xadj,Adj0,Perm,Invp)
    case(ORDERING_METIS)
      call hecmw_ordering_METIS_NodeND(Neqns,Xadj,Adj0,Perm,Invp)
    end select
  end subroutine hecmw_ordering_gen

end module hecmw_ordering
