!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Coupling Interface                                 !
!                                                                      !
!            Written by Shin'ichi Ezure (RIST)                         !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!


module hecmw_couple_info_f

  use hecmw_util

  implicit none
  private
  public :: hecmw_couple_get_unit_id
  public :: hecmw_couple_is_member
  public :: hecmw_couple_is_unit_member
  public :: hecmw_couple_is_unit_member_u
  public :: hecmw_couple_is_root
  public :: hecmw_couple_is_unit_root
  public :: hecmw_couple_is_unit_root_u
  public :: hecmw_intercomm_get_size
  public :: hecmw_intracomm_get_size
  public :: hecmw_intracomm_get_size_u
  public :: hecmw_intercomm_get_rank
  public :: hecmw_intracomm_get_rank
  public :: hecmw_intracomm_get_rank_u
  public :: hecmw_intercomm_get_comm
  public :: hecmw_intracomm_get_comm
  public :: hecmw_intracomm_get_comm_u
  public :: hecmw_intercomm_get_group
  public :: hecmw_intracomm_get_group
  public :: hecmw_intracomm_get_group_u

contains

!...
!...
!...
subroutine hecmw_couple_get_unit_id(boundary_id, unit_specifier, unit_id)

  character(len=HECMW_NAME_LEN), intent(in)  :: boundary_id
  integer(kind=kint),            intent(in)  :: unit_specifier
  character(len=HECMW_NAME_LEN), intent(out) :: unit_id
  integer(kind=kint)                         :: ierr

  call hecmw_couple_get_unit_id_if(boundary_id, unit_specifier, unit_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

end subroutine hecmw_couple_get_unit_id

!...
!...
!...
function hecmw_couple_is_member(boundary_id) result (is_member)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: is_member

  call hecmw_couple_is_memb_if(boundary_id, is_member)

end function hecmw_couple_is_member

!...
!...
!...
function hecmw_couple_is_unit_member(boundary_id, unit_specifier) result(is_member)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: is_member

  call hecmw_couple_is_unit_memb_if(boundary_id, unit_specifier, is_member)

end function hecmw_couple_is_unit_member

!...
!...
!...
function hecmw_couple_is_unit_member_u(unit_id) result(is_member)

  character(len=HECMW_NAME_LEN), intent(in) :: unit_id
  integer(kind=kint)                        :: is_member

  call hecmw_couple_is_unit_memb_u_if(unit_id, is_member)

end function hecmw_couple_is_unit_member_u

!...
!...
!...
function hecmw_couple_is_root(boundary_id) result(is_root)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: is_root

  call hecmw_couple_is_root_if(boundary_id, is_root)

end function hecmw_couple_is_root

!...
!...
!...
function hecmw_couple_is_unit_root(boundary_id, unit_specifier) result(is_root)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: is_root

  call hecmw_couple_is_unit_root_if(boundary_id, unit_specifier, is_root)

end function hecmw_couple_is_unit_root

!...
!...
!...
function hecmw_couple_is_unit_root_u(unit_id) result(is_root)

  character(len=HECMW_NAME_LEN), intent(in) :: unit_id
  integer(kind=kint)                        :: is_root

  call hecmw_couple_is_unit_root_u_if(unit_id, is_root)

end function hecmw_couple_is_unit_root_u

!...
!...
!...
function hecmw_intercomm_get_size(boundary_id) result(psize)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: psize

  call hecmw_intercomm_get_size_if(boundary_id, psize)

end function hecmw_intercomm_get_size

!...
!...
!...
function hecmw_intracomm_get_size(boundary_id, unit_specifier) result(psize)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: psize

  call hecmw_intracomm_get_size_if(boundary_id, unit_specifier, psize)

end function hecmw_intracomm_get_size

!...
!...
!...
function hecmw_intracomm_get_size_u(unit_id) result(psize)

  character(len=HECMW_NAME_LEN), intent(in) :: unit_id
  integer(kind=kint)                        :: psize

  call hecmw_intracomm_get_size_u_if(unit_id, psize)

end function hecmw_intracomm_get_size_u

!...
!...
!...
function hecmw_intercomm_get_rank(boundary_id) result(rank)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: rank

  call hecmw_intercomm_get_rank_if(boundary_id, rank)

end function hecmw_intercomm_get_rank

!...
!...
!...
function hecmw_intracomm_get_rank(boundary_id, unit_specifier) result(rank)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: rank

  call hecmw_intracomm_get_rank_if(boundary_id, unit_specifier, rank)

end function hecmw_intracomm_get_rank

!...
!...
!...
function hecmw_intracomm_get_rank_u(unit_id) result(rank)

  character(len=HECMW_NAME_LEN), intent(in) :: unit_id
  integer(kind=kint)                        :: rank

  call hecmw_intracomm_get_rank_u_if(unit_id, rank)

end function hecmw_intracomm_get_rank_u

!...
!...
!...
function hecmw_intercomm_get_comm(boundary_id) result(comm)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: comm

  call hecmw_intercomm_get_comm_if(boundary_id, comm)

end function hecmw_intercomm_get_comm

!...
!...
!...
function hecmw_intracomm_get_comm(boundary_id, unit_specifier) result(comm)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: comm

  call hecmw_intracomm_get_comm_if(boundary_id, unit_specifier, comm)

end function hecmw_intracomm_get_comm

!...
!...
!...
function hecmw_intracomm_get_comm_u(unit_id) result(comm)

  character(len=HECMW_NAME_LEN), intent(in) :: unit_id
  integer(kind=kint)                        :: comm

  call hecmw_intracomm_get_comm_u_if(unit_id, comm)

end function hecmw_intracomm_get_comm_u

!...
!...
!...
function hecmw_intercomm_get_group(boundary_id) result(group)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: group

  call hecmw_intercomm_get_group_if(boundary_id, group)

end function hecmw_intercomm_get_group

!...
!...
!...
function hecmw_intracomm_get_group(boundary_id, unit_specifier) result(group)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint),            intent(in) :: unit_specifier
  integer(kind=kint)                        :: group

  call hecmw_intracomm_get_group_if(boundary_id, unit_specifier, group)

end function hecmw_intracomm_get_group

!...
!...
!...
function hecmw_intracomm_get_group_u(boundary_id) result(group)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: group

  call hecmw_intracomm_get_group_u_if(boundary_id, group)

end function hecmw_intracomm_get_group_u

end module hecmw_couple_info_f

