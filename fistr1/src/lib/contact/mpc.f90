!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provide functions of contact stiffness calculation
module m_mpc
  use elementInfo

  implicit none

  private
  public :: tMPCCond
  public :: mpc_cond_init
  public :: print_mpc_cond
  public :: get_mpc_cond

  integer, parameter :: kint = 4
  integer, parameter :: kreal = 8

  type tMPCCond
    integer                        :: nitem
    integer, allocatable           :: pid(:)
    integer, allocatable           :: dof(:)
    real(kind=kreal), allocatable  :: coeff(:)
  end type

contains

  !> Initializer
  subroutine mpc_cond_init(mpc,nitem)
    type(tMPCCond), intent(inout)  :: mpc !< mpc condition
    integer(kind=kint), intent(in) :: nitem

    mpc%nitem = nitem
    allocate(mpc%pid(nitem))
    mpc%pid = 0
    allocate(mpc%dof(nitem))
    mpc%dof = 0
    allocate(mpc%coeff(nitem))
    mpc%coeff = 0.d0
  end subroutine

  !> get mpc condition from tMPCCond
  subroutine get_mpc_cond(mpc, nitem, pid, dof, coeff)
    type(tMPCCond), intent(in)      :: mpc !< mpc condition
    integer(kind=kint), intent(out)         :: nitem
    integer(kind=kint), intent(out)         :: pid(:)
    integer(kind=kint), intent(out)         :: dof(:)
    real(kind=kreal), intent(out)           :: coeff(:)

    nitem = mpc%nitem
    pid(1:nitem) =mpc%pid(1:nitem)
    dof(1:nitem) =mpc%dof(1:nitem)
    coeff(1:nitem) =mpc%coeff(1:nitem)
  end subroutine

  !> Print out mpc condition in !EQUATION format
  subroutine print_mpc_cond(fnum, mpc, idlabel)
    integer, intent(in)             :: fnum !< file number
    type(tMPCCond), intent(in)      :: mpc !< mpc condition
    integer, optional, intent(in)   :: idlabel(:) !< id label for input file

    integer(kind=kint) :: i, id

    write(fnum,"(I0,A)") mpc%nitem,", 0.0"
    do i=1,mpc%nitem
      id = mpc%pid(i)
      if( present(idlabel) ) id = idlabel(id)
      write(fnum,"(I0,',',I0,',',f12.6)") id,mpc%dof(i),mpc%coeff(i)
    enddo

  end subroutine

  !> Print out mpc condition in !EQUATION format
  subroutine print_mpc_conditions(fnum, mpcs, idlabel)
    integer, intent(in)             :: fnum !< file number
    type(tMPCCond), intent(in)      :: mpcs(:) !< mpc condition
    integer, optional, intent(in)   :: idlabel(:) !< id label for input file

    integer(kind=kint) :: i

    write(fnum,'(A)') "!EQUATION"
    do i=1,size(mpcs)
      call print_mpc_cond(fnum, mpcs(i), idlabel)
    enddo
  end subroutine


end module m_mpc


