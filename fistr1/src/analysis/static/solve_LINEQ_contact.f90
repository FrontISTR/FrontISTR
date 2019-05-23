!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to solve sparse system of
!> \linear equitions in the case of contact analysis using standard
!> \Lagrange multiplier algorithm

module m_solve_LINEQ_contact

  use m_fstr
  use m_solve_LINEQ_MKL_contact
  use m_solve_LINEQ_direct_serial_lag
  use m_solve_LINEQ_MUMPS_contact
  use m_solve_LINEQ_iter_contact
  use m_fstr_mat_resid_contact

  implicit none

  private
  public :: solve_LINEQ_contact_init
  public :: solve_LINEQ_contact

contains

  !> \brief This subroutine
  subroutine solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange)
    logical :: is_sym

    if( hecMAT%Iarray(99)==1 )then
      call solve_LINEQ_iter_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
    elseif( hecMAT%Iarray(99)==2 )then
      call solve_LINEQ_serial_lag_hecmw_init(hecMAT,fstrMAT,is_sym)
    else if( hecMAT%Iarray(99)==3 )then
      call solve_LINEQ_MKL_contact_init(hecMESH,is_sym)
    elseif( hecMAT%Iarray(99)==5 ) then
      call solve_LINEQ_mumps_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
    endif
  end subroutine solve_LINEQ_contact_init


  !> \brief This subroutine
  subroutine solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,istat,rf,conMAT)

    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange)
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), optional               :: rf
    type (hecmwST_matrix),optional           :: conMAT

    real(kind=kreal)                         :: factor
    real(kind=kreal) :: resid
    real(kind=kreal) :: t1, t2
    integer(kind=kint) :: ndof

    factor = 1.0d0
    if( present(rf) )factor = rf

    t1 = hecmw_wtime()

    istat = 0
    if( hecMAT%Iarray(99)==1 )then
      if(paraContactFlag.and.present(conMAT)) then
        call solve_LINEQ_iter_contact(hecMESH,hecMAT,fstrMAT,istat,conMAT)
      else
        call solve_LINEQ_iter_contact(hecMESH,hecMAT,fstrMAT,istat)
      endif
    elseif( hecMAT%Iarray(99)==2 )then
      call solve_LINEQ_serial_lag_hecmw(hecMESH,hecMAT,fstrMAT)
    elseif( hecMAT%Iarray(99)==3 )then
      if(paraContactFlag.and.present(conMAT)) then
        call solve_LINEQ_MKL_contact(hecMESH,hecMAT,fstrMAT,istat,conMAT)
      else
        call solve_LINEQ_MKL_contact(hecMESH,hecMAT,fstrMAT,istat)
      endif
    elseif( hecMAT%Iarray(99)==5 ) then
      ! ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        call solve_LINEQ_mumps_contact(hecMESH,hecMAT,fstrMAT,istat,conMAT)
      else
        call solve_LINEQ_mumps_contact(hecMESH,hecMAT,fstrMAT,istat)
      endif
    endif

    ndof = hecMAT%NDOF
    if( ndof==3 ) then
      call hecmw_update_3_R(hecMESH,hecMAT%X,hecMESH%n_node)
    else if( ndof==2 ) then
      call hecmw_update_2_R(hecMESH,hecMAT%X,hecMESH%n_node)
    else if( ndof==4 ) then
      call hecmw_update_4_R(hecMESH,hecMAT%X,hecMESH%n_node)
    else if( ndof==6 ) then
      call hecmw_update_m_R(hecMESH,hecMAT%X,hecMESH%n_node,6)
    endif

    t2 = hecmw_wtime()
    if (hecmw_mat_get_timelog(hecMAT) .ge. 1) then
      if (myrank==0) write(*,*) ' solve time :', t2 - t1
    endif

    if(paraContactFlag.and.present(conMAT)) then
    else
      resid=fstr_get_resid_max_contact(hecMESH,hecMAT,fstrMAT)
      if (myrank==0) then
        write(*,*) ' maximum residual = ', resid
        if( hecmw_mat_get_solver_type(hecMAT) /= 1 .and. resid >= 1.0d-8) then
          write(*,*) ' ###Maximum residual exceeded 1.0d-8---Direct Solver### '
          !          stop
        endif
      endif
    endif

    hecMAT%X=factor*hecMAT%X

  end subroutine solve_LINEQ_contact

end module m_solve_LINEQ_contact
