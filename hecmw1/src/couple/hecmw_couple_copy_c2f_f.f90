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


module hecmw_couple_copy_c2f_f

  use hecmw_util
  use hecmw_couple_define_f
  use hecmw_couple_struct_f

  implicit none
  private
  public :: hecmw_couple_copy_c2f

contains

subroutine hecmw_couple_copy_c2f(couple_value, ierr)

  type(hecmw_couple_value),      intent(inout) :: couple_value
  integer(kind=kint),            intent(inout) :: ierr
  integer(kind=kint)                           :: ista, is_allocated
  character(len=HECMW_NAME_LEN)                :: sname, vname

  sname = "hecmw_couple_value"

  vname = "n"
  call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%n, ierr)
  if(ierr /= 0) return

  vname = "item_type"
  call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%item_type, ierr)
  if(ierr /= 0) return

  vname = "n_dof"
  call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%n_dof, ierr)
  if(ierr /= 0) return

  if(couple_value%n > 0) then
    vname = "item"
    call hecmw_cpl_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
    if(is_allocated == 1) then
      if(couple_value%item_type == HECMW_COUPLE_NODE_GROUP) then
        allocate(couple_value%item(couple_value%n), stat=ista)
        if(ista > 0) return
        call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%item, ierr)
        if(ierr /= 0) return
      else if(couple_value%item_type == HECMW_COUPLE_ELEMENT_GROUP) then
        allocate(couple_value%item(couple_value%n), stat=ista)
        if(ista > 0) return
        call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%item, ierr)
        if(ierr /= 0) return
      else if(couple_value%item_type == HECMW_COUPLE_SURFACE_GROUP) then
        allocate(couple_value%item(couple_value%n*2), stat=ista)
        if(ista > 0) return
        call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%item, ierr)
        if(ierr /= 0) return
      else
        return
      endif
    endif
  endif

  if(couple_value%n > 0 .AND. couple_value%n_dof > 0) then
    vname = "value"
    call hecmw_cpl_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
    if(is_allocated == 1) then
      allocate(couple_value%value(couple_value%n*couple_value%n_dof), stat=ista)
      if(ista > 0) return
      call hecmw_cpl_copy_c2f_set_if(sname, vname, couple_value%value, ierr)
      if(ierr /= 0) return
    endif
  endif

end subroutine hecmw_couple_copy_c2f

end module hecmw_couple_copy_c2f_f
