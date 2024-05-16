!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_matrix_misc
  use hecmw_util
  use hecmw_matrix_contact
  use m_hecmw_comm_f
  implicit none

  private
  public :: hecmw_mat_clear
  public :: hecmw_mat_clear_b
  public :: hecmw_mat_init
  public :: hecmw_mat_finalize
  public :: hecmw_mat_copy_profile
  public :: hecmw_mat_copy_val

  public :: hecmw_mat_set_iter
  public :: hecmw_mat_get_iter
  public :: hecmw_mat_set_method
  public :: hecmw_mat_get_method
  public :: hecmw_mat_set_precond
  public :: hecmw_mat_get_precond
  public :: hecmw_mat_set_nset
  public :: hecmw_mat_get_nset
  public :: hecmw_mat_set_iterpremax
  public :: hecmw_mat_get_iterpremax
  public :: hecmw_mat_set_nrest
  public :: hecmw_mat_get_nrest
  public :: hecmw_mat_set_nbfgs
  public :: hecmw_mat_get_nbfgs
  public :: hecmw_mat_set_scaling
  public :: hecmw_mat_get_scaling
  public :: hecmw_mat_set_penalized
  public :: hecmw_mat_get_penalized
  public :: hecmw_mat_set_penalized_b
  public :: hecmw_mat_get_penalized_b
  public :: hecmw_mat_set_mpc_method
  public :: hecmw_mat_get_mpc_method
  public :: hecmw_mat_set_estcond
  public :: hecmw_mat_get_estcond
  public :: hecmw_mat_set_contact_elim
  public :: hecmw_mat_get_contact_elim
  public :: hecmw_mat_set_iterlog
  public :: hecmw_mat_get_iterlog
  public :: hecmw_mat_set_timelog
  public :: hecmw_mat_get_timelog
  public :: hecmw_mat_set_dump
  public :: hecmw_mat_get_dump
  public :: hecmw_mat_set_dump_exit
  public :: hecmw_mat_get_dump_exit
  public :: hecmw_mat_set_usejad
  public :: hecmw_mat_get_usejad
  public :: hecmw_mat_set_ncolor_in
  public :: hecmw_mat_get_ncolor_in
  public :: hecmw_mat_set_maxrecycle_precond
  public :: hecmw_mat_get_maxrecycle_precond
  public :: hecmw_mat_get_nrecycle_precond
  public :: hecmw_mat_reset_nrecycle_precond
  public :: hecmw_mat_incr_nrecycle_precond
  public :: hecmw_mat_set_flag_numfact
  public :: hecmw_mat_get_flag_numfact
  public :: hecmw_mat_set_flag_symbfact
  public :: hecmw_mat_get_flag_symbfact
  public :: hecmw_mat_clear_flag_symbfact
  public :: hecmw_mat_set_solver_type
  public :: hecmw_mat_get_solver_type

  public :: hecmw_mat_set_method2
  public :: hecmw_mat_get_method2
  public :: hecmw_mat_set_flag_converged
  public :: hecmw_mat_get_flag_converged
  public :: hecmw_mat_set_flag_diverged
  public :: hecmw_mat_get_flag_diverged
  public :: hecmw_mat_set_flag_mpcmatvec
  public :: hecmw_mat_get_flag_mpcmatvec

  public :: hecmw_mat_set_solver_opt
  public :: hecmw_mat_get_solver_opt

  public :: hecmw_mat_set_resid
  public :: hecmw_mat_get_resid
  public :: hecmw_mat_set_sigma_diag
  public :: hecmw_mat_get_sigma_diag
  public :: hecmw_mat_set_sigma
  public :: hecmw_mat_get_sigma
  public :: hecmw_mat_set_thresh
  public :: hecmw_mat_get_thresh
  public :: hecmw_mat_set_filter
  public :: hecmw_mat_get_filter
  public :: hecmw_mat_set_penalty
  public :: hecmw_mat_get_penalty
  public :: hecmw_mat_set_penalty_alpha
  public :: hecmw_mat_get_penalty_alpha

  public :: hecmw_mat_diag_max
  public :: hecmw_mat_recycle_precond_setting
  public :: hecmw_mat_substitute

  integer, parameter :: IDX_I_ITER               = 1
  integer, parameter :: IDX_I_METHOD             = 2
  integer, parameter :: IDX_I_PRECOND            = 3
  integer, parameter :: IDX_I_NSET               = 4
  integer, parameter :: IDX_I_ITERPREMAX         = 5
  integer, parameter :: IDX_I_NREST              = 6
  integer, parameter :: IDX_I_NBFGS              = 60
  integer, parameter :: IDX_I_SCALING            = 7
  integer, parameter :: IDX_I_PENALIZED          = 11
  integer, parameter :: IDX_I_PENALIZED_B        = 12
  integer, parameter :: IDX_I_MPC_METHOD         = 13
  integer, parameter :: IDX_I_ESTCOND            = 14
  integer, parameter :: IDX_I_CONTACT_ELIM       = 15
  integer, parameter :: IDX_I_ITERLOG            = 21
  integer, parameter :: IDX_I_TIMELOG            = 22
  integer, parameter :: IDX_I_DUMP               = 31
  integer, parameter :: IDX_I_DUMP_EXIT          = 32
  integer, parameter :: IDX_I_USEJAD             = 33
  integer, parameter :: IDX_I_NCOLOR_IN          = 34
  integer, parameter :: IDX_I_MAXRECYCLE_PRECOND = 35
  integer, parameter :: IDX_I_NRECYCLE_PRECOND   = 96
  integer, parameter :: IDX_I_FLAG_NUMFACT       = 97
  integer, parameter :: IDX_I_FLAG_SYMBFACT      = 98
  integer, parameter :: IDX_I_SOLVER_TYPE        = 99

  integer, parameter :: IDX_I_METHOD2            = 8
  integer, parameter :: IDX_I_FLAG_CONVERGED     = 81
  integer, parameter :: IDX_I_FLAG_DIVERGED      = 82
  integer, parameter :: IDX_I_FLAG_MPCMATVEC     = 83

  integer, parameter :: IDX_I_SOLVER_OPT_S       = 41
  integer, parameter :: IDX_I_SOLVER_OPT_E       = 50

  integer, parameter :: IDX_R_RESID         = 1
  integer, parameter :: IDX_R_SIGMA_DIAG    = 2
  integer, parameter :: IDX_R_SIGMA         = 3
  integer, parameter :: IDX_R_THRESH        = 4
  integer, parameter :: IDX_R_FILTER        = 5
  integer, parameter :: IDX_R_PENALTY       = 11
  integer, parameter :: IDX_R_PENALTY_ALPHA = 12

contains

  subroutine hecmw_mat_clear( hecMAT )
    type(hecmwST_matrix) :: hecMAT

    hecMAT%D = 0.0d0
    hecMAT%AL = 0.0d0
    hecMAT%AU = 0.0d0
    call hecmw_cmat_clear( hecMAT%cmat )
    call hecmw_mat_set_penalized( hecMAT, 0 )
    call hecmw_mat_set_penalty_alpha( hecMAT, 0.d0 )
  end subroutine hecmw_mat_clear

  subroutine hecmw_mat_clear_b( hecMAT )
    type(hecmwST_matrix) :: hecMAT

    hecMAT%B = 0.0d0
    call hecmw_mat_set_penalized_b( hecMAT, 0 )
  end subroutine hecmw_mat_clear_b

  subroutine hecmw_mat_init( hecMAT )
    type(hecmwST_matrix) :: hecMAT

    call hecmw_nullify_matrix( hecMAT )

    hecMAT%Iarray = 0
    hecMAT%Rarray = 0.d0

    call hecmw_mat_set_iter( hecMAT, 100 )
    call hecmw_mat_set_method( hecMAT, 1 )
    call hecmw_mat_set_precond( hecMAT, 1 )
    call hecmw_mat_set_nset( hecMAT, 0 )
    call hecmw_mat_set_iterpremax( hecMAT, 1 )
    call hecmw_mat_set_nrest( hecMAT, 10 )
    call hecmw_mat_set_nbfgs( hecMAT, 0 )
    call hecmw_mat_set_scaling( hecMAT, 0 )
    call hecmw_mat_set_iterlog( hecMAT, 0 )
    call hecmw_mat_set_timelog( hecMAT, 0 )
    call hecmw_mat_set_dump( hecMAT, 0 )
    call hecmw_mat_set_dump_exit( hecMAT, 0 )
    call hecmw_mat_set_usejad( hecMAT, 0 )
    call hecmw_mat_set_ncolor_in( hecMAT, 10 )
    call hecmw_mat_set_estcond( hecMAT, 0 )
    call hecmw_mat_set_maxrecycle_precond( hecMAT, 3 )

    call hecmw_mat_set_resid( hecMAT, 1.d-8 )
    call hecmw_mat_set_sigma_diag( hecMAT, 1.d0 )
    call hecmw_mat_set_sigma( hecMAT, 0.d0 )
    call hecmw_mat_set_thresh( hecMAT, 0.10d0 )
    call hecmw_mat_set_filter( hecMAT, 0.10d0 )

    call hecmw_mat_set_penalized( hecMAT, 0 )
    call hecmw_mat_set_penalty( hecMAT, 1.d+4 )
    call hecmw_mat_set_penalty_alpha( hecMAT, 0.d0 )
    call hecmw_mat_set_mpc_method( hecMAT, 0 )

    call hecmw_mat_reset_nrecycle_precond( hecMAT )
    call hecmw_mat_set_flag_numfact( hecMAT, 1 )
    call hecmw_mat_set_flag_symbfact( hecMAT, 1 )
    call hecmw_mat_set_solver_type( hecMAT, 1 )

    call hecmw_cmat_init( hecMAT%cmat )
  end subroutine hecmw_mat_init

  subroutine hecmw_mat_finalize( hecMAT )
    type(hecmwST_matrix) :: hecMAT
    if (associated(hecMAT%D)) deallocate(hecMAT%D)
    if (associated(hecMAT%B)) deallocate(hecMAT%B)
    if (associated(hecMAT%X)) deallocate(hecMAT%X)
    if (associated(hecMAT%AL)) deallocate(hecMAT%AL)
    if (associated(hecMAT%AU)) deallocate(hecMAT%AU)
    if (associated(hecMAT%indexL)) deallocate(hecMAT%indexL)
    if (associated(hecMAT%indexU)) deallocate(hecMAT%indexU)
    if (associated(hecMAT%itemL)) deallocate(hecMAT%itemL)
    if (associated(hecMAT%itemU)) deallocate(hecMAT%itemU)
    if (associated(hecMAT%ALU)) deallocate(hecMAT%ALU)
    call hecmw_cmat_finalize( hecMAT%cmat )
  end subroutine hecmw_mat_finalize

  subroutine hecmw_mat_copy_profile( hecMATorg, hecMAT )
    type(hecmwST_matrix), intent(in) :: hecMATorg
    type(hecmwST_matrix), intent(inout) :: hecMAT
    hecMAT%N    = hecMATorg%N
    hecMAT%NP   = hecMATorg%NP
    hecMAT%NDOF = hecMATorg%NDOF
    hecMAT%NPL  = hecMATorg%NPL
    hecMAT%NPU  = hecMATorg%NPU
    allocate(hecMAT%indexL(0:size(hecMATorg%indexL)-1))
    allocate(hecMAT%indexU(0:size(hecMATorg%indexU)-1))
    allocate(hecMAT%itemL (size(hecMATorg%itemL )))
    allocate(hecMAT%itemU (size(hecMATorg%itemU )))
    allocate(hecMAT%D (size(hecMATorg%D )))
    allocate(hecMAT%AL(size(hecMATorg%AL)))
    allocate(hecMAT%AU(size(hecMATorg%AU)))
    allocate(hecMAT%B (size(hecMATorg%B )))
    allocate(hecMAT%X (size(hecMATorg%X )))
    hecMAT%indexL = hecMATorg%indexL
    hecMAT%indexU = hecMATorg%indexU
    hecMAT%itemL  = hecMATorg%itemL
    hecMAT%itemU  = hecMATorg%itemU
    hecMAT%D  = 0.d0
    hecMAT%AL = 0.d0
    hecMAT%AU = 0.d0
    hecMAT%B  = 0.d0
    hecMAT%X  = 0.d0
  end subroutine hecmw_mat_copy_profile

  subroutine hecmw_mat_copy_val( hecMATorg, hecMAT )
    type(hecmwST_matrix), intent(in) :: hecMATorg
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: ierr
    integer(kind=kint) :: i
    ierr = 0
    if (hecMAT%N    /= hecMATorg%N) ierr = 1
    if (hecMAT%NP   /= hecMATorg%NP) ierr = 1
    if (hecMAT%NDOF /= hecMATorg%NDOF) ierr = 1
    if (hecMAT%NPL  /= hecMATorg%NPL) ierr = 1
    if (hecMAT%NPU  /= hecMATorg%NPU) ierr = 1
    if (ierr /= 0) then
      write(0,*) 'ERROR: hecmw_mat_copy_val: different profile'
      stop
    endif
    do i = 1, size(hecMAT%D)
      hecMAT%D(i)  = hecMATorg%D(i)
    enddo
    do i = 1, size(hecMAT%AL)
      hecMAT%AL(i) = hecMATorg%AL(i)
    enddo
    do i = 1, size(hecMAT%AU)
      hecMAT%AU(i) = hecMATorg%AU(i)
    enddo
  end subroutine hecmw_mat_copy_val

  subroutine hecmw_mat_set_iter( hecMAT, iter )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iter

    hecMAT%Iarray(IDX_I_ITER) = iter
  end subroutine hecmw_mat_set_iter

  function hecmw_mat_get_iter( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_iter
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iter = hecMAT%Iarray(IDX_I_ITER)
  end function hecmw_mat_get_iter

  subroutine hecmw_mat_set_method( hecMAT, method )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: method

    hecMAT%Iarray(IDX_I_METHOD) = method
  end subroutine hecmw_mat_set_method

  function hecmw_mat_get_method( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_method
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_method = hecMAT%Iarray(IDX_I_METHOD)
  end function hecmw_mat_get_method

  subroutine hecmw_mat_set_method2( hecMAT, method2 )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: method2

    hecMAT%Iarray(IDX_I_METHOD2) = method2
  end subroutine hecmw_mat_set_method2

  function hecmw_mat_get_method2( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_method2
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_method2 = hecMAT%Iarray(IDX_I_METHOD2)
  end function hecmw_mat_get_method2

  subroutine hecmw_mat_set_precond( hecMAT, precond )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: precond

    hecMAT%Iarray(IDX_I_PRECOND) = precond
  end subroutine hecmw_mat_set_precond

  function hecmw_mat_get_precond( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_precond
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_precond = hecMAT%Iarray(IDX_I_PRECOND)
  end function hecmw_mat_get_precond

  subroutine hecmw_mat_set_nset( hecMAT, nset )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nset

    hecMAT%Iarray(IDX_I_NSET) = nset
  end subroutine hecmw_mat_set_nset

  function hecmw_mat_get_nset( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_nset
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_nset = hecMAT%Iarray(IDX_I_NSET)
  end function hecmw_mat_get_nset

  subroutine hecmw_mat_set_iterpremax( hecMAT, iterpremax )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iterpremax

    if (iterpremax.lt.0) iterpremax= 0
    if (iterpremax.gt.4) iterpremax= 4

    hecMAT%Iarray(IDX_I_ITERPREMAX) = iterpremax
  end subroutine hecmw_mat_set_iterpremax

  function hecmw_mat_get_iterPREmax( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_iterPREmax
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iterPREmax = hecMAT%Iarray(IDX_I_ITERPREMAX)
  end function hecmw_mat_get_iterPREmax

  subroutine hecmw_mat_set_nrest( hecMAT, nrest )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nrest

    hecMAT%Iarray(IDX_I_NREST) = nrest
  end subroutine hecmw_mat_set_nrest

  function hecmw_mat_get_nrest( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_nrest
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_nrest = hecMAT%Iarray(IDX_I_NREST)
  end function hecmw_mat_get_nrest

  subroutine hecmw_mat_set_nbfgs( hecMAT, nbfgs )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nbfgs

    hecMAT%Iarray(IDX_I_NBFGS) = nbfgs
  end subroutine hecmw_mat_set_nbfgs

  function hecmw_mat_get_nbfgs( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_nbfgs
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_nbfgs = hecMAT%Iarray(IDX_I_NBFGS)
  end function hecmw_mat_get_nbfgs

  subroutine hecmw_mat_set_scaling( hecMAT, scaling )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: scaling

    hecMAT%Iarray(IDX_I_SCALING) = scaling
  end subroutine hecmw_mat_set_scaling

  function hecmw_mat_get_scaling( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_scaling
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_scaling = hecMAT%Iarray(IDX_I_SCALING)
  end function hecmw_mat_get_scaling

  subroutine hecmw_mat_set_penalized( hecMAT, penalized )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: penalized

    hecMAT%Iarray(IDX_I_PENALIZED) = penalized
  end subroutine hecmw_mat_set_penalized

  function hecmw_mat_get_penalized( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_penalized
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalized = hecMAT%Iarray(IDX_I_PENALIZED)
  end function hecmw_mat_get_penalized

  subroutine hecmw_mat_set_penalized_b( hecMAT, penalized_b )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: penalized_b

    hecMAT%Iarray(IDX_I_PENALIZED_B) = penalized_b
  end subroutine hecmw_mat_set_penalized_b

  function hecmw_mat_get_penalized_b( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_penalized_b
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalized_b = hecMAT%Iarray(IDX_I_PENALIZED_B)
  end function hecmw_mat_get_penalized_b

  subroutine hecmw_mat_set_mpc_method( hecMAT, mpc_method )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: mpc_method

    hecMAT%Iarray(IDX_I_MPC_METHOD) = mpc_method
  end subroutine hecmw_mat_set_mpc_method

  function hecmw_mat_get_mpc_method( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_mpc_method
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_mpc_method = hecMAT%Iarray(IDX_I_MPC_METHOD)
  end function hecmw_mat_get_mpc_method

  function hecmw_mat_get_estcond( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_estcond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_estcond = hecMAT%Iarray(IDX_I_ESTCOND)
  end function hecmw_mat_get_estcond

  subroutine hecmw_mat_set_estcond( hecMAT, estcond )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: estcond
    hecMAT%Iarray(IDX_I_ESTCOND) = estcond
  end subroutine hecmw_mat_set_estcond

  function hecmw_mat_get_contact_elim( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_contact_elim
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_contact_elim = hecMAT%Iarray(IDX_I_CONTACT_ELIM)
  end function hecmw_mat_get_contact_elim

  subroutine hecmw_mat_set_contact_elim( hecMAT, contact_elim )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: contact_elim
    hecMAT%Iarray(IDX_I_CONTACT_ELIM) = contact_elim
  end subroutine hecmw_mat_set_contact_elim

  subroutine hecmw_mat_set_iterlog( hecMAT, iterlog )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iterlog

    hecMAT%Iarray(IDX_I_ITERLOG) = iterlog
  end subroutine hecmw_mat_set_iterlog

  function hecmw_mat_get_iterlog( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_iterlog
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iterlog = hecMAT%Iarray(IDX_I_ITERLOG)
  end function hecmw_mat_get_iterlog

  subroutine hecmw_mat_set_timelog( hecMAT, timelog )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: timelog

    hecMAT%Iarray(IDX_I_TIMELOG) = timelog
  end subroutine hecmw_mat_set_timelog

  function hecmw_mat_get_timelog( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_timelog
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_timelog = hecMAT%Iarray(IDX_I_TIMELOG)
  end function hecmw_mat_get_timelog

  function hecmw_mat_get_dump( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_dump
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_dump = hecMAT%Iarray(IDX_I_DUMP)
  end function hecmw_mat_get_dump

  subroutine hecmw_mat_set_dump( hecMAT, dump_type )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: dump_type
    hecMAT%Iarray(IDX_I_DUMP) = dump_type
  end subroutine hecmw_mat_set_dump

  function hecmw_mat_get_dump_exit( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_dump_exit
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_dump_exit = hecMAT%Iarray(IDX_I_DUMP_EXIT)
  end function hecmw_mat_get_dump_exit

  subroutine hecmw_mat_set_dump_exit( hecMAT, dump_exit )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: dump_exit
    hecMAT%Iarray(IDX_I_DUMP_EXIT) = dump_exit
  end subroutine hecmw_mat_set_dump_exit

  function hecmw_mat_get_usejad( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_usejad
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_usejad = hecMAT%Iarray(IDX_I_USEJAD)
  end function hecmw_mat_get_usejad

  subroutine hecmw_mat_set_usejad( hecMAT, usejad )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: usejad
    hecMAT%Iarray(IDX_I_USEJAD) = usejad
  end subroutine hecmw_mat_set_usejad

  function hecmw_mat_get_ncolor_in( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_ncolor_in
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_ncolor_in = hecMAT%Iarray(IDX_I_NCOLOR_IN)
  end function hecmw_mat_get_ncolor_in

  subroutine hecmw_mat_set_ncolor_in( hecMAT, ncolor_in )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: ncolor_in
    hecMAT%Iarray(IDX_I_NCOLOR_IN) = ncolor_in
  end subroutine hecmw_mat_set_ncolor_in

  function hecmw_mat_get_maxrecycle_precond( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_maxrecycle_precond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_maxrecycle_precond = hecMAT%Iarray(IDX_I_MAXRECYCLE_PRECOND)
  end function hecmw_mat_get_maxrecycle_precond

  subroutine hecmw_mat_set_maxrecycle_precond( hecMAT, maxrecycle_precond )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: maxrecycle_precond
    if (maxrecycle_precond > 100) maxrecycle_precond = 100
    hecMAT%Iarray(IDX_I_MAXRECYCLE_PRECOND) = maxrecycle_precond
  end subroutine hecmw_mat_set_maxrecycle_precond

  function hecmw_mat_get_nrecycle_precond( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_nrecycle_precond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_nrecycle_precond = hecMAT%Iarray(IDX_I_NRECYCLE_PRECOND)
  end function hecmw_mat_get_nrecycle_precond

  subroutine hecmw_mat_reset_nrecycle_precond( hecMAT )
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(IDX_I_NRECYCLE_PRECOND) = 0
  end subroutine hecmw_mat_reset_nrecycle_precond

  subroutine hecmw_mat_incr_nrecycle_precond( hecMAT )
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(IDX_I_NRECYCLE_PRECOND) = hecMAT%Iarray(IDX_I_NRECYCLE_PRECOND) + 1
  end subroutine hecmw_mat_incr_nrecycle_precond

  function hecmw_mat_get_flag_numfact( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_flag_numfact
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_numfact = hecMAT%Iarray(IDX_I_FLAG_NUMFACT)
  end function hecmw_mat_get_flag_numfact

  subroutine hecmw_mat_set_flag_numfact( hecMAT, flag_numfact )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_numfact
    hecMAT%Iarray(IDX_I_FLAG_NUMFACT) = flag_numfact
  end subroutine hecmw_mat_set_flag_numfact

  function hecmw_mat_get_flag_symbfact( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_flag_symbfact
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_symbfact = hecMAT%Iarray(IDX_I_FLAG_SYMBFACT)
  end function hecmw_mat_get_flag_symbfact

  subroutine hecmw_mat_set_flag_symbfact( hecMAT, flag_symbfact )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_symbfact
    hecMAT%Iarray(IDX_I_FLAG_SYMBFACT) = flag_symbfact
  end subroutine hecmw_mat_set_flag_symbfact

  subroutine hecmw_mat_clear_flag_symbfact( hecMAT )
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(IDX_I_FLAG_SYMBFACT) = 0
  end subroutine hecmw_mat_clear_flag_symbfact

  function hecmw_mat_get_solver_type( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_solver_type
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_solver_type = hecMAT%Iarray(IDX_I_SOLVER_TYPE)
  end function hecmw_mat_get_solver_type

  subroutine hecmw_mat_set_solver_type( hecMAT, solver_type )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: solver_type
    hecMAT%Iarray(IDX_I_SOLVER_TYPE) = solver_type
  end subroutine hecmw_mat_set_solver_type

  subroutine hecmw_mat_set_flag_converged( hecMAT, flag_converged )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_converged
    hecMAT%Iarray(IDX_I_FLAG_CONVERGED) = flag_converged
  end subroutine hecmw_mat_set_flag_converged

  function hecmw_mat_get_flag_converged( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_flag_converged
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_converged = hecMAT%Iarray(IDX_I_FLAG_CONVERGED)
  end function hecmw_mat_get_flag_converged

  subroutine hecmw_mat_set_flag_diverged( hecMAT, flag_diverged )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_diverged
    hecMAT%Iarray(IDX_I_FLAG_DIVERGED) = flag_diverged
  end subroutine hecmw_mat_set_flag_diverged

  function hecmw_mat_get_flag_diverged( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_flag_diverged
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_diverged = hecMAT%Iarray(IDX_I_FLAG_DIVERGED)
  end function hecmw_mat_get_flag_diverged

  subroutine hecmw_mat_set_flag_mpcmatvec( hecMAT, flag_mpcmatvec )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_mpcmatvec
    hecMAT%Iarray(IDX_I_FLAG_MPCMATVEC) = flag_mpcmatvec
  end subroutine hecmw_mat_set_flag_mpcmatvec

  function hecmw_mat_get_flag_mpcmatvec( hecMAT )
    integer(kind=kint) :: hecmw_mat_get_flag_mpcmatvec
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_mpcmatvec = hecMAT%Iarray(IDX_I_FLAG_MPCMATVEC)
  end function hecmw_mat_get_flag_mpcmatvec

  subroutine hecmw_mat_set_solver_opt( hecMAT, solver_opt )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: solver_opt(:)
    integer(kind=kint) :: nopt
    nopt = IDX_I_SOLVER_OPT_E - IDX_I_SOLVER_OPT_S + 1
    hecMAT%Iarray(IDX_I_SOLVER_OPT_S:IDX_I_SOLVER_OPT_E) = solver_opt(1:nopt)
  end subroutine hecmw_mat_set_solver_opt

  subroutine hecmw_mat_get_solver_opt( hecMAT, solver_opt )
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: solver_opt(:)
    integer(kind=kint) :: nopt
    nopt = IDX_I_SOLVER_OPT_E - IDX_I_SOLVER_OPT_S + 1
    solver_opt(1:nopt) = hecMAT%Iarray(IDX_I_SOLVER_OPT_S:IDX_I_SOLVER_OPT_E)
  end subroutine hecmw_mat_get_solver_opt

  subroutine hecmw_mat_set_resid( hecMAT, resid )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: resid

    hecMAT%Rarray(IDX_R_RESID) = resid
  end subroutine hecmw_mat_set_resid

  function hecmw_mat_get_resid( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_resid
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_resid = hecMAT%Rarray(IDX_R_RESID)
  end function hecmw_mat_get_resid

  subroutine hecmw_mat_set_sigma_diag( hecMAT, sigma_diag )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: sigma_diag

    if( sigma_diag < 0.d0 ) then
      hecMAT%Rarray(IDX_R_SIGMA_DIAG) = -1.d0
    elseif( sigma_diag < 1.d0 ) then
      hecMAT%Rarray(IDX_R_SIGMA_DIAG) = 1.d0
    elseif( sigma_diag > 2.d0 ) then
      hecMAT%Rarray(IDX_R_SIGMA_DIAG) = 2.d0
    else
      hecMAT%Rarray(IDX_R_SIGMA_DIAG) = sigma_diag
    endif
  end subroutine hecmw_mat_set_sigma_diag

  function hecmw_mat_get_sigma_diag( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_sigma_diag
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_sigma_diag = hecMAT%Rarray(IDX_R_SIGMA_DIAG)
  end function hecmw_mat_get_sigma_diag

  subroutine hecmw_mat_set_sigma( hecMAT, sigma )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: sigma

    if (sigma < 0.d0) then
      hecMAT%Rarray(IDX_R_SIGMA) = 0.d0
    elseif (sigma > 1.d0) then
      hecMAT%Rarray(IDX_R_SIGMA) = 1.d0
    else
      hecMAT%Rarray(IDX_R_SIGMA) = sigma
    endif
  end subroutine hecmw_mat_set_sigma

  function hecmw_mat_get_sigma( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_sigma
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_sigma = hecMAT%Rarray(IDX_R_SIGMA)
  end function hecmw_mat_get_sigma

  subroutine hecmw_mat_set_thresh( hecMAT, thresh )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: thresh

    hecMAT%Rarray(IDX_R_THRESH) = thresh
  end subroutine hecmw_mat_set_thresh

  function hecmw_mat_get_thresh( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_thresh
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_thresh = hecMAT%Rarray(IDX_R_THRESH)
  end function hecmw_mat_get_thresh

  subroutine hecmw_mat_set_filter( hecMAT, filter )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: filter

    hecMAT%Rarray(IDX_R_FILTER) = filter
  end subroutine hecmw_mat_set_filter

  function hecmw_mat_get_filter( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_filter
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_filter = hecMAT%Rarray(IDX_R_FILTER)
  end function hecmw_mat_get_filter

  subroutine hecmw_mat_set_penalty( hecMAT, penalty )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: penalty

    hecMAT%Rarray(IDX_R_PENALTY) = penalty
  end subroutine hecmw_mat_set_penalty

  function hecmw_mat_get_penalty( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_penalty
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalty = hecMAT%Rarray(IDX_R_PENALTY)
  end function hecmw_mat_get_penalty

  subroutine hecmw_mat_set_penalty_alpha( hecMAT, alpha )
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: alpha

    hecMAT%Rarray(IDX_R_PENALTY_ALPHA) = alpha
  end subroutine hecmw_mat_set_penalty_alpha

  function hecmw_mat_get_penalty_alpha( hecMAT )
    real(kind=kreal) :: hecmw_mat_get_penalty_alpha
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalty_alpha = hecMAT%Rarray(IDX_R_PENALTY_ALPHA)
  end function hecmw_mat_get_penalty_alpha

  function hecmw_mat_diag_max(hecMAT, hecMESH)
    real(kind=kreal) :: hecmw_mat_diag_max
    type (hecmwST_matrix) :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: ndiag, i

    hecmw_mat_diag_max = -1.0e20
    ndiag = hecMAT%NDOF**2 * hecMAT%NP
    do i = 1, ndiag
      if( hecMAT%D(i) > hecmw_mat_diag_max ) hecmw_mat_diag_max = hecMAT%D(i)
    enddo
    call hecmw_allREDUCE_R1(hecMESH, hecmw_mat_diag_max, hecmw_max)
  end function hecmw_mat_diag_max

  subroutine hecmw_mat_recycle_precond_setting( hecMAT )
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nrecycle, maxrecycle
    if (hecMAT%Iarray(IDX_I_FLAG_SYMBFACT) >= 1) then
      hecMAT%Iarray(IDX_I_FLAG_NUMFACT)=1
      call hecmw_mat_reset_nrecycle_precond(hecMAT)
    elseif (hecMAT%Iarray(IDX_I_FLAG_NUMFACT) > 1) then
      call hecmw_mat_reset_nrecycle_precond(hecMAT)
      hecMAT%Iarray(IDX_I_FLAG_NUMFACT) = 1
    elseif (hecMAT%Iarray(IDX_I_FLAG_NUMFACT) == 1) then
      nrecycle = hecmw_mat_get_nrecycle_precond(hecMAT)
      maxrecycle = hecmw_mat_get_maxrecycle_precond(hecMAT)
      if ( nrecycle < maxrecycle ) then
        hecMAT%Iarray(IDX_I_FLAG_NUMFACT) = 0
        call hecmw_mat_incr_nrecycle_precond(hecMAT)
      else
        call hecmw_mat_reset_nrecycle_precond(hecMAT)
      endif
    endif
  end subroutine hecmw_mat_recycle_precond_setting

  subroutine hecmw_mat_substitute( dest, src )
    type (hecmwST_matrix), intent(inout) :: dest
    type (hecmwST_matrix), intent(inout) :: src
    dest%N = src%N
    dest%NP = src%NP
    dest%NPL = src%NPL
    dest%NPU = src%NPU
    dest%NDOF = src%NDOF
    dest%NPCL = src%NPCU
    if (associated(src%D)) dest%D => src%D
    if (associated(src%B)) dest%B => src%B
    if (associated(src%X)) dest%X => src%X
    if (associated(src%ALU)) dest%ALU => src%ALU
    if (associated(src%AL)) dest%AL => src%AL
    if (associated(src%AU)) dest%AU => src%AU
    if (associated(src%CAL)) dest%CAL => src%CAL
    if (associated(src%indexL)) dest%indexL => src%indexL
    if (associated(src%indexU)) dest%indexU => src%indexU
    if (associated(src%indexCL)) dest%indexCL => src%indexCL
    if (associated(src%indexCU)) dest%indexCU => src%indexCU
    if (associated(src%itemL)) dest%itemL => src%itemL
    if (associated(src%itemU)) dest%itemU => src%itemU
    if (associated(src%itemCL)) dest%itemCL => src%itemCL
    if (associated(src%itemCU)) dest%itemCU => src%itemCU
    dest%Iarray(:) = src%Iarray(:)
    dest%Rarray(:) = src%Rarray(:)
    call hecmw_cmat_substitute( dest%cmat, src%cmat )
  end subroutine hecmw_mat_substitute

end module hecmw_matrix_misc
