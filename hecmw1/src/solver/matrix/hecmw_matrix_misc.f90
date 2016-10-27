!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2014/06/11                                         !
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

module hecmw_matrix_misc
  use hecmw_util
  use hecmw_matrix_contact
  use m_hecmw_comm_f

  contains

  subroutine hecmw_mat_clear( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT

    hecMAT%D = 0.0d0
    hecMAT%AL = 0.0d0
    hecMAT%AU = 0.0d0
    call hecmw_cmat_clear( hecMAT%cmat )
    call hecmw_mat_set_penalized( hecMAT, 0 )
    call hecmw_mat_set_penalty_alpha( hecMAT, 0.d0 )
  end subroutine hecmw_mat_clear

  subroutine hecmw_mat_clear_b( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT

    hecMAT%B = 0.0d0
    call hecmw_mat_set_penalized_b( hecMAT, 0 )
  end subroutine hecmw_mat_clear_b

  subroutine hecmw_mat_init( hecMAT )
    implicit none
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
    call hecmw_mat_set_mpc_method( hecMAT, 3 )

    call hecmw_mat_reset_nrecycle_precond( hecMAT )
    call hecmw_mat_set_flag_numfact( hecMAT, 1 )
    call hecmw_mat_set_flag_symbfact( hecMAT, 1 )
    call hecmw_mat_set_solver_type( hecMAT, 1 )

    call hecmw_cmat_init( hecMAT%cmat )
  end subroutine hecmw_mat_init

  subroutine hecmw_mat_finalize( hecMAT )
    implicit none
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

  subroutine hecmw_mat_set_iter( hecMAT, iter )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iter

    hecMAT%Iarray(1) = iter
  end subroutine hecmw_mat_set_iter

  function hecmw_mat_get_iter( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_iter
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iter = hecMAT%Iarray(1)
  end function hecmw_mat_get_iter

  subroutine hecmw_mat_set_method( hecMAT, method )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: method

    hecMAT%Iarray(2) = method
  end subroutine hecmw_mat_set_method

  function hecmw_mat_get_method( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_method
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_method = hecMAT%Iarray(2)
  end function hecmw_mat_get_method

  subroutine hecmw_mat_set_precond( hecMAT, precond )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: precond

    hecMAT%Iarray(3) = precond
  end subroutine hecmw_mat_set_precond

  function hecmw_mat_get_precond( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_precond
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_precond = hecMAT%Iarray(3)
  end function hecmw_mat_get_precond

  subroutine hecmw_mat_set_nset( hecMAT, nset )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nset

    hecMAT%Iarray(4) = nset
  end subroutine hecmw_mat_set_nset

  function hecmw_mat_get_nset( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_nset
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_nset = hecMAT%Iarray(4)
  end function hecmw_mat_get_nset

  subroutine hecmw_mat_set_iterpremax( hecMAT, iterpremax )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iterpremax

        if (iterpremax.lt.0) iterpremax= 0
        if (iterpremax.gt.4) iterpremax= 4

    hecMAT%Iarray(5) = iterpremax
  end subroutine hecmw_mat_set_iterpremax

  function hecmw_mat_get_iterPREmax( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_iterPREmax
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iterPREmax = hecMAT%Iarray(5)
  end function hecmw_mat_get_iterPREmax

  subroutine hecmw_mat_set_nrest( hecMAT, nrest )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nrest

    hecMAT%Iarray(6) = nrest
  end subroutine hecmw_mat_set_nrest

  function hecmw_mat_get_nrest( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_nrest
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_nrest = hecMAT%Iarray(6)
  end function hecmw_mat_get_nrest

  subroutine hecmw_mat_set_scaling( hecMAT, scaling )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: scaling

    hecMAT%Iarray(7) = scaling
  end subroutine hecmw_mat_set_scaling

  function hecmw_mat_get_scaling( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_scaling
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_scaling = hecMAT%Iarray(7)
  end function hecmw_mat_get_scaling

  subroutine hecmw_mat_set_penalized( hecMAT, penalized )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: penalized

    hecMAT%Iarray(11) = penalized
  end subroutine hecmw_mat_set_penalized

  function hecmw_mat_get_penalized( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_penalized
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalized = hecMAT%Iarray(11)
  end function hecmw_mat_get_penalized

  subroutine hecmw_mat_set_penalized_b( hecMAT, penalized_b )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: penalized_b

    hecMAT%Iarray(12) = penalized_b
  end subroutine hecmw_mat_set_penalized_b

  function hecmw_mat_get_penalized_b( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_penalized_b
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalized_b = hecMAT%Iarray(12)
  end function hecmw_mat_get_penalized_b

  subroutine hecmw_mat_set_mpc_method( hecMAT, mpc_method )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: mpc_method

    hecMAT%Iarray(13) = mpc_method
  end subroutine hecmw_mat_set_mpc_method

  function hecmw_mat_get_mpc_method( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_mpc_method
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_mpc_method = hecMAT%Iarray(13)
  end function hecmw_mat_get_mpc_method

  function hecmw_mat_get_estcond( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_estcond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_estcond = hecMAT%Iarray(14)
  end function hecmw_mat_get_estcond

  subroutine hecmw_mat_set_estcond( hecMAT, estcond )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: estcond
    hecMAT%Iarray(14) = estcond
  end subroutine hecmw_mat_set_estcond

  subroutine hecmw_mat_set_iterlog( hecMAT, iterlog )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: iterlog

    hecMAT%Iarray(21) = iterlog
  end subroutine hecmw_mat_set_iterlog

  function hecmw_mat_get_iterlog( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_iterlog
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_iterlog = hecMAT%Iarray(21)
  end function hecmw_mat_get_iterlog

  subroutine hecmw_mat_set_timelog( hecMAT, timelog )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: timelog

    hecMAT%Iarray(22) = timelog
  end subroutine hecmw_mat_set_timelog

  function hecmw_mat_get_timelog( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_timelog
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_timelog = hecMAT%Iarray(22)
  end function hecmw_mat_get_timelog

  function hecmw_mat_get_dump( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_dump
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_dump = hecMAT%Iarray(31)
  end function hecmw_mat_get_dump

  subroutine hecmw_mat_set_dump( hecMAT, dump_type )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: dump_type
    hecMAT%Iarray(31) = dump_type
  end subroutine hecmw_mat_set_dump

  function hecmw_mat_get_dump_exit( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_dump_exit
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_dump_exit = hecMAT%Iarray(32)
  end function hecmw_mat_get_dump_exit

  subroutine hecmw_mat_set_dump_exit( hecMAT, dump_exit )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: dump_exit
    hecMAT%Iarray(32) = dump_exit
  end subroutine hecmw_mat_set_dump_exit

  function hecmw_mat_get_usejad( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_usejad
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_usejad = hecMAT%Iarray(33)
  end function hecmw_mat_get_usejad

  subroutine hecmw_mat_set_usejad( hecMAT, usejad )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: usejad
    hecMAT%Iarray(33) = usejad
  end subroutine hecmw_mat_set_usejad

  function hecmw_mat_get_ncolor_in( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_ncolor_in
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_ncolor_in = hecMAT%Iarray(34)
  end function hecmw_mat_get_ncolor_in

  subroutine hecmw_mat_set_ncolor_in( hecMAT, ncolor_in )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: ncolor_in
    hecMAT%Iarray(34) = ncolor_in
  end subroutine hecmw_mat_set_ncolor_in

  function hecmw_mat_get_maxrecycle_precond( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_maxrecycle_precond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_maxrecycle_precond = hecMAT%Iarray(35)
  end function hecmw_mat_get_maxrecycle_precond

  subroutine hecmw_mat_set_maxrecycle_precond( hecMAT, maxrecycle_precond )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: maxrecycle_precond
    if (maxrecycle_precond > 100) maxrecycle_precond = 100
    hecMAT%Iarray(35) = maxrecycle_precond
  end subroutine hecmw_mat_set_maxrecycle_precond

  function hecmw_mat_get_nrecycle_precond( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_nrecycle_precond
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_nrecycle_precond = hecMAT%Iarray(96)
  end function hecmw_mat_get_nrecycle_precond

  subroutine hecmw_mat_reset_nrecycle_precond( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(96) = 0
  end subroutine hecmw_mat_reset_nrecycle_precond

  subroutine hecmw_mat_incr_nrecycle_precond( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(96) = hecMAT%Iarray(96) + 1
  end subroutine hecmw_mat_incr_nrecycle_precond

  function hecmw_mat_get_flag_numfact( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_flag_numfact
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_numfact = hecMAT%Iarray(97)
  end function hecmw_mat_get_flag_numfact

  subroutine hecmw_mat_set_flag_numfact( hecMAT, flag_numfact )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_numfact
    hecMAT%Iarray(97) = flag_numfact
  end subroutine hecmw_mat_set_flag_numfact

  function hecmw_mat_get_flag_symbfact( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_flag_symbfact
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_flag_symbfact = hecMAT%Iarray(98)
  end function hecmw_mat_get_flag_symbfact

  subroutine hecmw_mat_set_flag_symbfact( hecMAT, flag_symbfact )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: flag_symbfact
    hecMAT%Iarray(98) = flag_symbfact
  end subroutine hecmw_mat_set_flag_symbfact

  subroutine hecmw_mat_clear_flag_symbfact( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    hecMAT%Iarray(98) = 0
  end subroutine hecmw_mat_clear_flag_symbfact

  function hecmw_mat_get_solver_type( hecMAT )
    implicit none
    integer(kind=kint) :: hecmw_mat_get_solver_type
    type(hecmwST_matrix) :: hecMAT
    hecmw_mat_get_solver_type = hecMAT%Iarray(99)
  end function hecmw_mat_get_solver_type

  subroutine hecmw_mat_set_solver_type( hecMAT, solver_type )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: solver_type
    hecMAT%Iarray(99) = solver_type
  end subroutine hecmw_mat_set_solver_type

  subroutine hecmw_mat_set_resid( hecMAT, resid )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: resid

    hecMAT%Rarray(1) = resid
  end subroutine hecmw_mat_set_resid

  function hecmw_mat_get_resid( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_resid
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_resid = hecMAT%Rarray(1)
  end function hecmw_mat_get_resid

  subroutine hecmw_mat_set_sigma_diag( hecMAT, sigma_diag )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: sigma_diag

    if( sigma_diag < 0.d0 ) then
      hecMAT%Rarray(2) = -1.d0
    elseif( sigma_diag < 1.d0 ) then
      hecMAT%Rarray(2) = 1.d0
    elseif( sigma_diag > 2.d0 ) then
      hecMAT%Rarray(2) = 2.d0
    else
      hecMAT%Rarray(2) = sigma_diag
    endif
  end subroutine hecmw_mat_set_sigma_diag

  function hecmw_mat_get_sigma_diag( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_sigma_diag
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_sigma_diag = hecMAT%Rarray(2)
  end function hecmw_mat_get_sigma_diag

  subroutine hecmw_mat_set_sigma( hecMAT, sigma )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: sigma

    if (sigma < 0.d0) then
      hecMAT%Rarray(3) = 0.d0
    elseif (sigma > 1.d0) then
      hecMAT%Rarray(3) = 1.d0
    else
      hecMAT%Rarray(3) = sigma
    endif
  end subroutine hecmw_mat_set_sigma

  function hecmw_mat_get_sigma( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_sigma
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_sigma = hecMAT%Rarray(3)
  end function hecmw_mat_get_sigma

  subroutine hecmw_mat_set_thresh( hecMAT, thresh )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: thresh

    hecMAT%Rarray(4) = thresh
  end subroutine hecmw_mat_set_thresh

  function hecmw_mat_get_thresh( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_thresh
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_thresh = hecMAT%Rarray(4)
  end function hecmw_mat_get_thresh

  subroutine hecmw_mat_set_filter( hecMAT, filter )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: filter

    hecMAT%Rarray(5) = filter
  end subroutine hecmw_mat_set_filter

  function hecmw_mat_get_filter( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_filter
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_filter = hecMAT%Rarray(5)
  end function hecmw_mat_get_filter

  subroutine hecmw_mat_set_penalty( hecMAT, penalty )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: penalty

    hecMAT%Rarray(11) = penalty
  end subroutine hecmw_mat_set_penalty

  function hecmw_mat_get_penalty( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_penalty
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalty = hecMAT%Rarray(11)
  end function hecmw_mat_get_penalty

  subroutine hecmw_mat_set_penalty_alpha( hecMAT, alpha )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    real(kind=kreal) :: alpha

    hecMAT%Rarray(12) = alpha
  end subroutine hecmw_mat_set_penalty_alpha

  function hecmw_mat_get_penalty_alpha( hecMAT )
    implicit none
    real(kind=kreal) :: hecmw_mat_get_penalty_alpha
    type(hecmwST_matrix) :: hecMAT

    hecmw_mat_get_penalty_alpha = hecMAT%Rarray(12)
  end function hecmw_mat_get_penalty_alpha

  function hecmw_mat_diag_max(hecMAT, hecMESH)
    implicit none
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
    implicit none
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: nrecycle, maxrecycle
    if (hecMAT%Iarray(98) >= 1) then
      hecMAT%Iarray(97)=1
      call hecmw_mat_reset_nrecycle_precond(hecMAT)
    elseif (hecMAT%Iarray(97) > 1) then
      call hecmw_mat_reset_nrecycle_precond(hecMAT)
      hecMAT%Iarray(97) = 1
    elseif (hecMAT%Iarray(97) == 1) then
      nrecycle = hecmw_mat_get_nrecycle_precond(hecMAT)
      maxrecycle = hecmw_mat_get_maxrecycle_precond(hecMAT)
      if ( nrecycle < maxrecycle ) then
        hecMAT%Iarray(97) = 0
        call hecmw_mat_incr_nrecycle_precond(hecMAT)
      else
        call hecmw_mat_reset_nrecycle_precond(hecMAT)
      endif
    endif
  end subroutine hecmw_mat_recycle_precond_setting

  subroutine hecmw_mat_substitute( dest, src )
    implicit none
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
