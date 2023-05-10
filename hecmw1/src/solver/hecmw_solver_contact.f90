!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to solve sparse system of
!> \linear equitions in the case of contact analysis using standard
!> \Lagrange multiplier algorithm

module m_solve_LINEQ_contact

  use hecmw_util
  use m_solve_LINEQ_MKL_contact
  use m_solve_LINEQ_direct_serial_lag
  use m_solve_LINEQ_MUMPS_contact
  use m_solve_LINEQ_iter_contact
  use m_hecmw_mat_resid_contact
  use hecmw_matrix_misc
  use m_hecmw_comm_f

  implicit none

  private
  public :: solve_LINEQ_contact_init
  public :: solve_LINEQ_contact

contains

  !> \brief This subroutine
  subroutine solve_LINEQ_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_matrix_lagrange)           :: hecLagMAT        !< type hecmwST_matrix_lagrange)
    logical :: is_sym

    if( hecMAT%Iarray(99)==1 )then
      call solve_LINEQ_iter_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    elseif( hecMAT%Iarray(99)==2 )then
      call solve_LINEQ_serial_lag_hecmw_init(hecMAT,hecLagMAT,is_sym)
    else if( hecMAT%Iarray(99)==3 )then
      call solve_LINEQ_MKL_contact_init(hecMESH,is_sym)
    elseif( hecMAT%Iarray(99)==5 ) then
      call solve_LINEQ_mumps_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    endif
  end subroutine solve_LINEQ_contact_init


  !> \brief This subroutine
  subroutine solve_LINEQ_contact(hecMESH,hecMAT,hecLagMAT,conMAT,istat,rf,is_contact_active)

    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_matrix_lagrange)           :: hecLagMAT        !< type hecmwST_matrix_lagrange)
    type (hecmwST_matrix)                    :: conMAT
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), optional               :: rf
    logical                                  :: is_contact_active

    real(kind=kreal)                         :: factor
    real(kind=kreal) :: t1, t2
    integer(kind=kint) :: ndof

    factor = 1.0d0
    if( present(rf) )factor = rf

    t1 = hecmw_wtime()

    istat = 0
    if( hecMAT%Iarray(99)==1 )then
      call solve_LINEQ_iter_contact(hecMESH,hecMAT,hecLagMAT,istat,conMAT,is_contact_active)
    elseif( hecMAT%Iarray(99)==2 )then
      if( hecmw_comm_get_size() > 1) then
        write(*,*) 'ERROR: !SOLVER,METHOD=DIRECT not available in parallel contact analysis; please use MUMPS or DIRECTmkl instead'
        call hecmw_abort(hecmw_comm_get_comm())
      else
        call add_conMAT_to_hecMAT(hecMAT,conMAT,hecLagMat)
        call solve_LINEQ_serial_lag_hecmw(hecMESH,hecMAT,hecLagMAT)
      endif
    elseif( hecMAT%Iarray(99)==3 )then
      if( hecmw_comm_get_size() > 1) then
        call solve_LINEQ_MKL_contact(hecMESH,hecMAT,hecLagMAT,istat,conMAT)
      else
        call add_conMAT_to_hecMAT(hecMAT,conMAT,hecLagMat)
        call solve_LINEQ_MKL_contact(hecMESH,hecMAT,hecLagMAT,istat)
      endif
    elseif( hecMAT%Iarray(99)==5 ) then
      call solve_LINEQ_mumps_contact(hecMESH,hecMAT,hecLagMAT,istat,conMAT)
    endif

    ndof = hecMAT%NDOF
    call hecmw_update_R(hecMESH,hecMAT%X,hecMESH%n_node, ndof)

    t2 = hecmw_wtime()
    if (hecmw_mat_get_timelog(hecMAT) .ge. 1) then
      if ( hecmw_comm_get_rank() ==0) write(*,*) ' solve time :', t2 - t1
    endif

    hecMAT%X=factor*hecMAT%X

  end subroutine solve_LINEQ_contact


  subroutine add_conMAT_to_hecMAT(hecMAT,conMAT,hecLagMat)
    type (hecmwST_matrix),          intent(inout) :: hecMAT
    type (hecmwST_matrix),          intent(in)    :: conMAT
    type (hecmwST_matrix_lagrange), intent(in)    :: hecLagMAT


    integer(kind=kint) :: ndof,ndof2,i

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof

    do i=1,hecMAT%NP*ndof + hecLagMat%num_lagrange
      hecMAT%B(i) = hecMAT%B(i) + conMAT%B(i)
    enddo

    do i=1,hecMAT%NP*ndof2
      hecMAT%D(i) = hecMAT%D(i) + conMAT%D(i)
    enddo

    do i=1,hecMAT%NPL*ndof2
      hecMAT%AL(i) = hecMAT%AL(i) + conMAT%AL(i)
    enddo

    do i=1,hecMAT%NPU*ndof2
      hecMAT%AU(i) = hecMAT%AU(i) + conMAT%AU(i)
    enddo
  end subroutine add_conMAT_to_hecMAT


end module m_solve_LINEQ_contact
