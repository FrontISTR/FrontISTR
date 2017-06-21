!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
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
      USE hecmw_solver_44
      USE hecmw_solver_66
      USE hecmw_solver_direct
      USE hecmw_solver_direct_parallel
      USE hecmw_solver_direct_MUMPS
      USE hecmw_solver_direct_clusterMKL
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
        CASE(4)
!          WRITE(*,*) "Calling 4x4 Iterative Solver..."
          CALL hecmw_solve_44(hecMESH,hecMAT)
        CASE(5)
          !CALL hecmw_solve_mm(hecMESH,hecMAT)
!          WRITE(*,*) "FATAL: Solve_mm not yet available..."
          !call hecmw_abort( hecmw_comm_get_comm() )
          call fstr_substitute_hecmw_solver(hecMESH,hecMAT,6)
        CASE(6)
!          WRITE(*,*) "Calling 6x6 Iterative Solver..."
          CALL hecmw_solve_66(hecMESH,hecMAT)
        CASE(7:)
          !CALL hecmw_solve_mm(hecMESH,hecMAT)
!          WRITE(*,*) "FATAL: Solve_mm not yet available..."
          call hecmw_abort( hecmw_comm_get_comm() )
        END SELECT
!C
!* Call Direct Solver
      CASE(2:)
!C
!* Please note the following:
!* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
!* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)

        if (hecMAT%Iarray(97) .gt. 1) hecMAT%Iarray(97)=1

        if (hecMAT%Iarray(2) .eq. 104) then
          call hecmw_solve_direct_MUMPS(hecMESH, hecMAT)
        elseif (hecMAT%Iarray(2) .eq. 105) then
          call hecmw_solve_direct_ClusterMKL(hecMESH, hecMAT)
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
          resid=hecmw_rel_resid_L2_44(hecMESH,hecMAT)
        END SELECT
        myrank=hecmw_comm_get_rank()
        if (myrank==0) then
          write(*,"(a,1pe12.5)")'### Relative residual =', resid
          if( resid >= 1.0d-8) then
            write(*,"(a)")'### Relative residual exceeded 1.0d-8---Direct Solver### '
!            stop
          endif
        endif
!C
      END SELECT
!C
      RETURN

   end subroutine solve_LINEQ
  subroutine fstr_substitute_hecmw_solver(hecMESH,hecMATorig,NDOF)
    use hecmw
    use hecmw_solver_11
    use hecmw_solver_22
    use hecmw_solver_33
    use hecmw_solver_44
    use hecmw_solver_66
    type (hecmwST_local_mesh)     :: hecMESH
    type (hecmwST_matrix)         :: hecMATorig
    type (hecmwST_matrix),pointer :: hecMAT => null()
    integer(kind=kint) NDOF
    if (NDOF == hecMATorig%NDOF) then
      call hecmw_clone_matrix(hecMATorig,hecMAT)
    else if (NDOF < hecMATorig%NDOF) then
      call hecmw_abort( hecmw_comm_get_comm() ) 
    else
      call hecmw_blockmatrix_expand(hecMATorig,hecMAT,NDOF)
      call hecmw_cmat_init(hecMAT%cmat)    
    end if
    select case(NDOF)
      case(1)
        call hecmw_solve_11(hecMESH,hecMAT)
      case(2)
        call hecmw_solve_22(hecMESH,hecMAT)
      case(3)
        call hecmw_solve_33(hecMESH,hecMAT)
      case(4)
        call hecmw_solve_44(hecMESH,hecMAT)
      case(5)
        call hecmw_abort( hecmw_comm_get_comm() )
      case(6)
        call hecmw_solve_66(hecMESH,hecMAT)
      case(7:)
        call hecmw_abort( hecmw_comm_get_comm() )
    end select
    if (NDOF /= hecMATorig%NDOF) then
      call hecmw_vector_contract(hecMATorig,hecMAT,NDOF)
    end if 
  end subroutine fstr_substitute_hecmw_solver  
end module m_solve_LINEQ
