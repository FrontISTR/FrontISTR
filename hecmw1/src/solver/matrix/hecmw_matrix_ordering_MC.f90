!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
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

module m_hecmw_matrix_ordering_MC
  use hecmw_util
  implicit none

  private
  public :: hecmw_matrix_ordering_MC

contains

  subroutine hecmw_matrix_ordering_MC(N, indexL, itemL, indexU, itemU, &
       perm_cur, ncolor_in, ncolor_out, COLORindex, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(in) :: perm_cur(:)
    integer(kind=kint), intent(in) :: ncolor_in
    integer(kind=kint), intent(out) :: ncolor_out
    integer(kind=kint), intent(out) :: COLORindex(0:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    integer(kind=kint), allocatable :: iwk(:)
    integer(kind=kint) :: nn_color, cntall, cnt, color
    integer(kind=kint) :: i, inode, j, jnode
    allocate(iwk(N))
    iwk = 0
    nn_color = N / ncolor_in
    cntall = 0
    COLORindex(0) = 0
    do color=1,N
      cnt = 0
      do i=1,N
        inode = perm_cur(i)
        if (iwk(inode) > 0 .or. iwk(inode) == -1) cycle
        ! if (iwk(inode) == 0)
        iwk(inode) = color
        cntall = cntall + 1
        perm(cntall) = inode
        cnt = cnt + 1
        if (cnt == nn_color) exit
        if (cntall == N) exit
        ! mark all connected and uncolored nodes
        do j = indexL(inode-1)+1, indexL(inode)
          jnode = itemL(j)
          if (iwk(jnode) == 0) iwk(jnode) = -1
        end do
        do j = indexU(inode-1)+1, indexU(inode)
          jnode = itemU(j)
          if (jnode > N) cycle
          if (iwk(jnode) == 0) iwk(jnode) = -1
        end do
      end do
      COLORindex(color) = cntall
      if (cntall == N) then
        ncolor_out = color
        exit
      end if
      ! unmark all marked nodes
      do i=1,N
        if (iwk(i) == -1) iwk(i) = 0
      end do
    end do
    deallocate(iwk)
    ! make iperm
    do i=1,N
      iperm(perm(i)) = i
    end do
  end subroutine hecmw_matrix_ordering_MC

end module m_hecmw_matrix_ordering_MC
