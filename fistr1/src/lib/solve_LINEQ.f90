!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.0                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
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
   contains

   SUBROUTINE solve_LINEQ(hecMESH,hecMAT,imsg)
      USE hecmw
      USE hecmw_solver_11
      USE hecmw_solver_22
      USE hecmw_solver_33
      USE hecmw_solver_direct
      USE hecmw_solver_direct_parallel
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      INTEGER(kind=kint) nprocs, ierror, imsg
!C
      SELECT CASE(hecMAT%Iarray(99))
!C
!* Call Iterative Solver
      CASE (1)
!C
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
!C
!* Call Direct Solver
      CASE(2)
!C
!* Please note the following:
!* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
!* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)
        IF(hecMESH%PETOT.GT.1) THEN
          CALL hecmw_solve_direct_parallel(hecMESH,hecMAT,imsg)  
        ELSE
          CALL hecmw_solve_direct(hecMESH,hecMAT,imsg)  
        ENDIF
!!!     hecMAT%X = hecMAT%B -- leading stack overflow (intel9)
        do i=1,hecMAT%NP*hecMESH%n_dof
            hecMAT%X(i) = hecMAT%B(i)
        end do
!C
      END SELECT 
!C
       RETURN
   
   end subroutine solve_LINEQ
end module m_solve_LINEQ
