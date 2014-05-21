!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module m_hecmw_solve_error
      use hecmw_util

      integer(kind=kint), parameter :: &
           HECMW_SOLVER_ERROR_INCONS_PC = 1001, &
           HECMW_SOLVER_ERROR_ZERO_DIAG = 2001, &
           HECMW_SOLVER_ERROR_ZERO_RHS = 2002, &
           HECMW_SOLVER_ERROR_NOCONV_MAXIT = 3001, &
           HECMW_SOLVER_ERROR_DIVERGE_MAT = 3002, &
           HECMW_SOLVER_ERROR_DIVERGE_PC = 3003

contains

      subroutine hecmw_solve_error (hecMESH, IFLAG)
      use hecmw_util
      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: IFLAG

      if ( hecMESH%zero.eq.1 .and. &
           (IFLAG.eq.HECMW_SOLVER_ERROR_INCONS_PC .or. IFLAG.eq.HECMW_SOLVER_ERROR_ZERO_DIAG) ) then
        write (*,'(/a )')                                               &
     &           '###############################################'
        write (*,'( a )')                                               &
     &           '######## ERROR MESSAGE : LINEAR SOLVER ########'
        write (*,'( a/)')                                               &
     &           '###############################################'
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_INCONS_PC) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-E-1001'
          write (*,'( a/)')'    inconsistent solver/preconditioning'
        endif
!        call MPI_ABORT (hecMESH%MPI_COMM, ierr)
        call hecmw_abort( hecmw_comm_get_comm())
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_ZERO_DIAG) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-E-2001: '
          write (*,'( a/)')'    ZERO component in diagonal block'
        endif
!        call MPI_ABORT (hecMESH%MPI_COMM, ierr)
        call hecmw_abort( hecmw_comm_get_comm())
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_ZERO_RHS) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-2002: '
          write (*,'( a/)')'    ZERO RHS norm'
        endif
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_NOCONV_MAXIT) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-3001: '
          write (*,'( a/)')'    not converged within ceratin iterations'
        endif
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_DIVERGE_MAT) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-3002: '
          write (*,'( a/)')'    diverged due to indefinite or negative definite matrix'
        endif
      endif

      if (IFLAG.eq.HECMW_SOLVER_ERROR_DIVERGE_PC) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-3003: '
          write (*,'( a/)')'    diverged due to indefinite preconditioner'
        endif
      endif

      !stop " PROGRAM STOP:"

      end subroutine hecmw_solve_error

end module m_hecmw_solve_error

