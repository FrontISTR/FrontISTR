!-------------------------------------------------------------------------------
! Copyright (c) 2020 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_ORDERING_RCM is a program for fill-reducing ordering
!         for direct solver
!----------------------------------------------------------------------
module hecmw_ordering_rcm
  use hecmw_util
  implicit none

  private
  public :: hecmw_ordering_genrcm

contains

  subroutine hecmw_ordering_GENRCM(Neqns,Xadj,Adj0,Perm,Invp)
    use m_hecmw_matrix_ordering_CM
    implicit none
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Adj0(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(out):: Perm(:)
    integer(kind=kint), intent(out):: Invp(:)
    integer(kind=kint), allocatable:: indexL(:), indexU(:), itemU(:)
    integer(kind=kint) :: i
    allocate(indexL(0:Neqns))
    allocate(indexU(0:Neqns))
    allocate(itemU(1))
    ! 1-based numbering to 0-based numbering
    do i = 0, Neqns
      indexL(i) = Xadj(i+1)-1
    enddo
    ! empty index/item for upper triangle
    indexU(:) = 0
    itemU(:) = 0
    call hecmw_matrix_ordering_RCM(Neqns, indexL, Adj0, indexU, itemU, perm, invp)
    deallocate(indexL)
    deallocate(indexU)
    deallocate(itemU)
  end subroutine hecmw_ordering_GENRCM

end module hecmw_ordering_rcm
