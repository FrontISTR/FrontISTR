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
  integer(kind=kint), parameter:: ORDERING_RCM     = 3
  integer(kind=kint), parameter:: ORDERING_NMAX    = 3

  integer(kind=kint), parameter:: ORDERING_DEBUG = 0

contains

  !======================================================================!
  !> @brief hecmw_ordering_gen
  !======================================================================!
  subroutine hecmw_ordering_gen(Neqns,Nttbr,Xadj,Adj0,Perm,Invp,opt,loglevel)
    use hecmw_ordering_qmd
    use hecmw_ordering_metis
    use hecmw_ordering_rcm
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nttbr
    integer(kind=kint), intent(in):: Adj0(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(in):: opt
    integer(kind=kint), intent(in):: loglevel
    integer(kind=kint), intent(out):: Perm(:)
    integer(kind=kint), intent(out):: Invp(:)
    !------
    integer(kind=kint):: ordering
    ordering = opt
    if (ordering < 0 .or. ordering > ORDERING_NMAX) then
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
      if (loglevel > 0) write(*,*) 'Ordering method: QMD'
      call hecmw_ordering_GENQMD(Neqns,Nttbr,Xadj,Adj0,Perm,Invp)
    case(ORDERING_METIS)
      if (loglevel > 0) write(*,*) 'Ordering method: METIS_NodeND'
      call hecmw_ordering_METIS_NodeND(Neqns,Xadj,Adj0,Perm,Invp)
    case(ORDERING_RCM)
      if (loglevel > 0) write(*,*) 'Ordering method: RCM'
      call hecmw_ordering_GENRCM(Neqns,Xadj,Adj0,Perm,Invp)
    end select
    if (ORDERING_DEBUG > 0) then
      call write_nonzero_profile(Neqns, Xadj, Adj0, perm, invp)
      call write_perm(Neqns, perm, invp)
    endif
  end subroutine hecmw_ordering_gen

  subroutine write_nonzero_profile(N, index, item, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: index(:)
    integer(kind=kint), intent(in) :: item(:)
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: F_ORG = 901
    integer(kind=kint), parameter :: F_NEW = 902
    integer(kind=kint) :: i, j, irow, jcol
    open(F_ORG, file='nzprof_org.txt', status='replace')
    do irow = 1, N
      i = irow
      do j = index(i), index(i+1)-1
        jcol = item(j)
        write(F_ORG,*) irow, jcol
      end do
    end do
    close(F_ORG)
    open(F_NEW, file='nzprof_new.txt', status='replace')
    do irow = 1, N
      i = perm(irow)
      do j = index(i), index(i+1)-1
        jcol = item(j)
        write(F_NEW,*) irow, iperm(jcol)
      end do
    end do
    close(F_NEW)
  end subroutine write_nonzero_profile

  subroutine write_perm(N, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: F_PERM = 903
    integer(kind=kint) :: i
    open(F_PERM, file='perm_iperm.txt', status='replace')
    do i = 1, N
      write(F_PERM,*) perm(i), iperm(i)
    end do
    close(F_PERM)
  end subroutine write_perm

end module hecmw_ordering
