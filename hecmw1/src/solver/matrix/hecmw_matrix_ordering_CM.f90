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

module m_hecmw_matrix_ordering_CM
  use hecmw_util
  implicit none

  private
  public :: hecmw_matrix_ordering_CM
  public :: hecmw_matrix_ordering_RCM

contains

  subroutine hecmw_matrix_ordering_CM(N, indexL, itemL, indexU, itemU, &
       perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: NMINMAX = 5
    integer(kind=kint) :: i, nmin
    integer(kind=kint) :: mins(NMINMAX)
    integer(kind=kint), allocatable :: nlevel(:)
    integer(kind=kint), allocatable :: lv_index(:,:), lv_item(:,:)
    integer(kind=kint) :: nlevel_max, max_id
    call find_minimum_degrees(N, indexL, indexU, itemU, NMINMAX, nmin, mins)
    allocate(nlevel(nmin), lv_index(0:N,nmin), lv_item(N,nmin))
    ! perform CM ordering starting from each minimum degree node
    do i=1,nmin
      call ordering_CM_inner(N, indexL, itemL, indexU, itemU, mins(i), &
           nlevel(i), lv_index(:,i), lv_item(:,i))
    end do
    ! choose CM ordering with maximum nlevel
    nlevel_max = nlevel(1)
    max_id = 1
    do i=2,nmin
      if (nlevel(i) > nlevel_max) then
        nlevel_max = nlevel(i)
        max_id = i
      end if
    end do
    do i=1,N
      perm(i) = lv_item(i,max_id)
      iperm(perm(i)) = i
    end do
    deallocate(nlevel, lv_index, lv_item)
  end subroutine hecmw_matrix_ordering_CM

  subroutine hecmw_matrix_ordering_RCM(N, indexL, itemL, indexU, itemU, &
       perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    call hecmw_matrix_ordering_CM(N, indexL, itemL, indexU, itemU, perm, iperm)
    call reverse_ordering(N, perm, iperm)
  end subroutine hecmw_matrix_ordering_RCM

  subroutine ordering_CM_inner(N, indexL, itemL, indexU, itemU, nstart, &
       nlevel, lv_index, lv_item)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(in) :: nstart
    integer(kind=kint), intent(out) :: nlevel
    integer(kind=kint), intent(out) :: lv_index(0:)
    integer(kind=kint), intent(out) :: lv_item(:)
    integer(kind=kint), allocatable :: iwk(:)
    integer(kind=kint) :: level, cntall, cnt, j, jnode, k, knode
    allocate(iwk(N))
    iwk = 0
    lv_index(0) = 0
    ! first level
    iwk(nstart) = 1
    cntall = 1
    lv_item(1) = nstart
    lv_index(1) = 1
    ! other levels
    do level=2,N
      cnt = 0
      ! all nodes in previous level
      PRLV: do j = lv_index(level-2)+1, lv_index(level-1)
        jnode = lv_item(j)
        ! all connected nodes
        do k = indexL(jnode-1)+1, indexL(jnode)
          knode = itemL(k)
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            if (cntall == N) exit PRLV
          end if
        end do
        do k = indexU(jnode-1)+1, indexU(jnode)
          knode = itemU(k)
          if (knode > N) cycle
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            if (cntall == N) exit PRLV
          end if
        end do
      end do PRLV
      if (cnt == 0) then
        !write(*,*) 'DEBUG: choose any uncolored node..'
        do knode = 1, N
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            exit
          end if
        end do
      endif
      lv_index(level) = cntall
      if (cntall == N) then
        nlevel = level
        exit
      end if
    end do
    !write(*,*) 'DEBUG:: ordering_CM_inner: nstart, nlevel = ', nstart, nlevel
  end subroutine ordering_CM_inner

  subroutine find_minimum_degrees(N, indexL, indexU, itemU, nminmax, nmin, mins)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemU(:)
    integer(kind=kint), intent(in) :: nminmax
    integer(kind=kint), intent(out) :: nmin
    integer(kind=kint), intent(out) :: mins(nminmax)
    integer(kind=kint) :: degmin, i, deg, j
    degmin = N
    nmin = 0
    do i=1,N
      deg = indexL(i) - indexL(i-1)
      do j = indexU(i-1)+1, indexU(i)
        if (itemU(j) <= N) deg = deg + 1
      end do
      ! skip unconnected nodes
      if (deg == 0) cycle
      if (deg < degmin) then
        degmin = deg
        nmin = 1
        mins(1) = i
      else if (deg == degmin) then
        nmin = nmin + 1
        if (nmin <= nminmax) mins(nmin) = i
      end if
    end do
    !write(*,*) 'DEBUG:: find_minimum_degrees: nmin, deg = ', nmin, degmin
    if (nmin > nminmax) nmin = nminmax
  end subroutine find_minimum_degrees

  subroutine reverse_ordering(N, perm, iperm)
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(inout) :: perm(:), iperm(:)
    integer(kind=kint) :: i, N1
    N1 = N + 1
    do i=1,N
      perm(i) = N1 - perm(i)
      iperm(perm(i)) = i
    end do
  end subroutine reverse_ordering

end module m_hecmw_matrix_ordering_CM
