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

contains

      subroutine hecmw_solve_error (hecMESH, IFLAG)
      use hecmw_util
      implicit REAL*8 (A-H,O-Z)       
      type (hecmwST_local_mesh) :: hecMESH

      if ( hecMESH%zero.eq.1 .and. &
           (IFLAG.eq.1001 .or. IFLAG.eq.2001) ) then
        write (*,'(/a )')                                               &
     &           '###############################################'
        write (*,'( a )')                                               &
     &           '######## ERROR MESSAGE : LINEAR SOLVER ########'
        write (*,'( a/)')                                               &
     &           '###############################################'
      endif

      if (IFLAG.eq.1001) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-E-1001'
          write (*,'( a/)')'    inconsistent solver/preconditioning'
        endif
!        call MPI_ABORT (hecMESH%MPI_COMM, ierr)
        call hecmw_abort( hecmw_comm_get_comm())
      endif

      if (IFLAG.eq.2001) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-E-2001: '
          write (*,'( a/)')'    ZERO component in diagonal block'
        endif
!        call MPI_ABORT (hecMESH%MPI_COMM, ierr)
        call hecmw_abort( hecmw_comm_get_comm())
      endif

      if (IFLAG.eq.2002) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-2002: '
          write (*,'( a/)')'    ZERO RHS norm'
        endif
      endif

      if (IFLAG.eq.3001) then
        if (hecMESH%zero.eq.1) then
          write (*,'(/a )')'  #### HEC-MW-SOLVER-W-3001: '
          write (*,'( a/)')'    not converged within ceratin iterations'
        endif
      endif

      !stop " PROGRAM STOP:"

      end subroutine hecmw_solve_error

end module m_hecmw_solve_error

