!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                       Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This program is a HECMW interface to a set of linear iterative and direct
!! solvers. The interface may be called from within a HECMW application, with
!! an appropriate choice of TYPE (iterative, direct), and METHOD (depending
!! on the iterative solver used).
module m_solve_LINEQ
   implicit none

   private
   public :: solve_LINEQ

   contains

   SUBROUTINE solve_LINEQ(hecMESH,hecMAT,imsg)
      USE hecmw
      USE hecmw_solver_11
      USE hecmw_solver_22
      USE hecmw_solver_33
      USE hecmw_solver_direct
      USE hecmw_solver_direct_parallel
      USE hecmw_solver_direct_MUMPS
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      INTEGER(kind=kint) imsg, i, myrank
      real(kind=kreal) :: resid
!C
      SELECT CASE(hecMAT%Iarray(99))
!C
!* Call Iterative Solver
      CASE (1)
!C
        call hecmw_mat_dump(hecMAT, hecMESH)

        SELECT CASE(hecMESH%n_dof)
        CASE(1)
!          WRITE(*,*) "Calling 1x1 Iterative Solver..."
          CALL hecmw_solve_11(hecMESH,hecMAT)
        CASE(2)
!          WRITE(*,*) "Calling 2x2 Iterative Solver..."
          CALL hecmw_solve_22(hecMESH,hecMAT)
        CASE(3)
!          WRITE(*,*) "Calling 3x3 Iterative Solver..."
          CALL hecmw_solve_33(hecMESH,hecMAT)
        CASE(4:)
          !CALL hecmw_solve_mm(hecMESH,hecMAT)
!          WRITE(*,*) "FATAL: Solve_mm not yet available..."
          call hecmw_abort( hecmw_comm_get_comm())
        END SELECT

        call hecmw_mat_dump_solution(hecMAT)
!C
!* Call Direct Solver
      CASE(2:)
!C
!* Please note the following:
!* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
!* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)

        call hecmw_mat_ass_equation( hecMESH, hecMAT )

        call hecmw_mat_dump(hecMAT, hecMESH)

        if (hecMAT%Iarray(2) .eq. 104) then
          call hecmw_solve_direct_MUMPS(hecMESH, hecMAT)
        else
          IF(hecMESH%PETOT.GT.1) THEN
            CALL hecmw_solve_direct_parallel(hecMESH,hecMAT,imsg)
          ELSE
            CALL hecmw_solve_direct(hecMESH,hecMAT,imsg)
          ENDIF
!!!       hecMAT%X = hecMAT%B -- leading stack overflow (intel9)
          do i=1,hecMAT%NP*hecMESH%n_dof
              hecMAT%X(i) = hecMAT%B(i)
          end do
        endif

        SELECT CASE(hecMESH%n_dof)
        CASE(1)
          resid=hecmw_rel_resid_L2_11(hecMESH,hecMAT)
        CASE(2)
          resid=hecmw_rel_resid_L2_22(hecMESH,hecMAT)
        CASE(3)
          resid=hecmw_rel_resid_L2_33(hecMESH,hecMAT)
        CASE(4:)
          !resid=hecmw_rel_resid_L2_mm(hecMESH,hecMAT)
          resid=0.d0 !! TEMPORARY
        END SELECT
        myrank=hecmw_comm_get_rank()
        if (myrank==0) then
          write(*,*) ' relative residual =', resid
          if( resid >= 1.0d-8) then
            write(*,*) ' ###Relative residual exceeded 1.0d-8---Direct Solver### '
!            stop
          endif
        endif

        call hecmw_mat_dump_solution(hecMAT)
!C
      END SELECT
!C
       RETURN

   end subroutine solve_LINEQ

end module m_solve_LINEQ
