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
  public :: tMPCGroup
  public :: init_mpc_cond
  public :: finalize_mpc_cond
  public :: set_mpc_cond
  public :: get_mpc_cond
  public :: copy_mpc_cond
  public :: expand_mpc_cond
  public :: print_mpc_cond
  public :: print_mpc_conditions
  public :: print_full_mpc_conditions_3d

  integer, parameter :: kint = 4
  integer, parameter :: kreal = 8

  type tMPCCond
    integer                        :: nitem
    integer, allocatable           :: pid(:)
    integer, allocatable           :: dof(:)
    real(kind=kreal), allocatable  :: coeff(:)
  end type

  type tMPCGroup
    integer                     :: nitem
    integer                     :: nslaves
    type(tMPCCond), allocatable :: mpcs(:)
  end type
  
contains

  !> Initializer
  subroutine init_mpc_cond(mpc,nitem)
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

  !> Finalizer
  subroutine finalize_mpc_cond(mpc)
    type(tMPCCond), intent(inout)  :: mpc !< mpc condition

    mpc%nitem = 0
    deallocate(mpc%pid,mpc%dof,mpc%coeff)
  end subroutine

  !> set mpc condition from tMPCCond
  subroutine set_mpc_cond(mpc, nitem, pid, dof, coeff)
    type(tMPCCond), intent(inout)      :: mpc !< mpc condition
    integer(kind=kint), intent(in)         :: nitem
    integer(kind=kint), intent(in)         :: pid(:)
    integer(kind=kint), intent(in)         :: dof(:)
    real(kind=kreal), intent(in)           :: coeff(:)

    mpc%pid(1:nitem) = pid(1:nitem)
    mpc%dof(1:nitem) = dof(1:nitem)
    mpc%coeff(1:nitem) = coeff(1:nitem)
  end subroutine

  !> get mpc condition from tMPCCond
  subroutine get_mpc_cond(mpc, nitem, pid, dof, coeff)
    type(tMPCCond), intent(in)      :: mpc !< mpc condition
    integer(kind=kint), intent(out)         :: nitem
    integer(kind=kint), intent(out)         :: pid(:)
    integer(kind=kint), intent(out)         :: dof(:)
    real(kind=kreal), intent(out)           :: coeff(:)

    nitem = mpc%nitem
    pid(1:nitem) = mpc%pid(1:nitem)
    dof(1:nitem) = mpc%dof(1:nitem)
    coeff(1:nitem) = mpc%coeff(1:nitem)
  end subroutine

  !> copy mpc condition from tMPCCond
  subroutine copy_mpc_cond( mpcA, mpcB )
    type(tMPCCond), intent(in)   :: mpcA
    type(tMPCCond), intent(out)  :: mpcB

    call init_mpc_cond(mpcB,mpcA%nitem)
    call set_mpc_cond(mpcB, mpcA%nitem, mpcA%pid, mpcA%dof, mpcA%coeff)
  end subroutine

  subroutine expand_mpc_cond( mpcA, mpcB )
    type(tMPCCond), intent(inout) :: mpcA
    type(tMPCCond), intent(in)    :: mpcB

    integer(kind=kint) :: nitem_new, slave, target_id
    integer(kind=kint) :: i, j
    integer(kind=kint), allocatable :: pid_new(:), dof_new(:)
    real(kind=kreal), allocatable :: coeff_new(:)
    real(kind=kreal) :: target_coeff, ratio
    logical :: found_same_pid_dof

    !find id and coeff of target
    slave = mpcB%pid(1)
    do i=2,mpcA%nitem
      if( mpcA%pid(i) == slave ) then
        target_id = i
        target_coeff = mpcA%coeff(i)
        exit
      endif
    enddo

    nitem_new = mpcA%nitem+mpcB%nitem-2
    allocate(pid_new(nitem_new),dof_new(nitem_new),coeff_new(nitem_new))
    nitem_new = 0
    ! set mpcA
    do i=1,mpcA%nitem
      if( i == target_id ) cycle
      nitem_new = nitem_new + 1
      pid_new(nitem_new)   = mpcA%pid(i)
      dof_new(nitem_new)   = mpcA%dof(i)
      coeff_new(nitem_new) = mpcA%coeff(i)
    enddo
    ! set mpcB
    if( dabs(mpcB%coeff(1)) < 1.d-10 ) stop
    ratio = -target_coeff/mpcB%coeff(1)
    do i=2,mpcB%nitem
      ! if same pid-dof pair is already in mpcA, add coeff to it
      found_same_pid_dof = .false.
      do j=2,nitem_new
        if( pid_new(j) /= mpcB%pid(i) ) cycle
        if( dof_new(j) /= mpcB%dof(i) ) cycle
        coeff_new(j) = coeff_new(j) + ratio*mpcB%coeff(i)
        found_same_pid_dof = .true.
        exit
      enddo
      if( .not. found_same_pid_dof ) then
        nitem_new = nitem_new + 1
        pid_new(nitem_new)   = mpcB%pid(i)
        dof_new(nitem_new)   = mpcB%dof(i)
        coeff_new(nitem_new) = ratio*mpcB%coeff(i)
      endif
    enddo

    call finalize_mpc_cond(mpcA)
    call init_mpc_cond(mpcA,nitem_new)
    call set_mpc_cond(mpcA, nitem_new, pid_new, dof_new, coeff_new)
  end subroutine

  !> Print out mpc condition in !EQUATION format
  subroutine print_mpc_cond(fnum, mpc, idlabel, doflabel)
    integer, intent(in)             :: fnum !< file number
    type(tMPCCond), intent(in)      :: mpc !< mpc condition
    integer, optional, intent(in)   :: idlabel(:) !< id label for input file
    integer, optional, intent(in)   :: doflabel !< dof is set to this value if given

    integer(kind=kint) :: i, id, dof

    write(fnum,*) mpc%nitem,",",0.d0
    do i=1,mpc%nitem
      id = mpc%pid(i)
      if( present(idlabel) ) id = idlabel(id)
      dof = mpc%dof(i)
      if( present(doflabel) ) dof = doflabel
      write(fnum,*) id,",",dof,",",mpc%coeff(i)
    enddo

  end subroutine

  !> Print out mpc condition in !EQUATION format
  subroutine print_mpc_conditions(fnum, mpcs, idlabel, doflabel)
    integer, intent(in)             :: fnum !< file number
    type(tMPCCond), intent(in)      :: mpcs(:) !< mpc condition
    integer, optional, intent(in)   :: idlabel(:) !< id label for input file
    integer, optional, intent(in)   :: doflabel !< dof is set to this value if given

    integer(kind=kint) :: i

    write(fnum,'(A)') "!EQUATION"
    do i=1,size(mpcs)
      call print_mpc_cond(fnum, mpcs(i), idlabel, doflabel)
    enddo
  end subroutine

  !> Print out mpc condition in !EQUATION format
  subroutine print_full_mpc_conditions_3d( mpcs, nodeID, target_dof )
    type(tMPCCond), intent(in)   :: mpcs(:)
    integer, intent(in)          :: nodeID(:)
    integer, intent(in), optional :: target_dof

    integer(kind=kint) :: iunit, idof, idof_start, idof_end, ierror
    character(len=20) :: fname

    if( present(target_dof) ) then
      idof_start = target_dof
      idof_end = target_dof
    else
      idof_start = 1
      idof_end = 3
    endif

    do idof = idof_start, idof_end
      write(fname,'(a,i0,a)') 'tied_equation_',idof,'.dat'
      iunit = 200+idof
      open(iunit, file=fname, status='replace', iostat=ierror)
      call print_mpc_conditions(iunit, mpcs, nodeID, idof)
      close(iunit)
    enddo
  end subroutine

end module m_mpc


