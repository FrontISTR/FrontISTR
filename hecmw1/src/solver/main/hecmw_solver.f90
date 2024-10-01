!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver
contains

  subroutine hecmw_solve (hecMESH, hecMAT)

    use hecmw_util
    use hecmw_solver_las
    use hecmw_solver_iterative
    use hecmw_solver_direct
    use hecmw_solver_direct_parallel
    use hecmw_solver_direct_MUMPS
    use hecmw_solver_direct_MKL
    use hecmw_solver_direct_clusterMKL
    use hecmw_matrix_misc
    use m_hecmw_comm_f
    implicit none

    type (hecmwST_matrix), target :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH

    real(kind=kreal) :: resid
    integer(kind=kint) :: i, myrank, NDOF
    integer(kind=kint) :: imsg = 51
    NDOF=hecMAT%NDOF

    !C ERROR CHECK
    if(hecmw_solve_check_zerorhs(hecMESH, hecMAT))then
      hecMAT%X = 0.0d0
      return
    endif

    select case(hecMAT%Iarray(99))
        !* Call Iterative Solver
      case (1)
        call hecmw_solve_iterative(hecMESH,hecMAT)
        !* Call Direct Solver
      case(2:)
        !C
        !* Please note the following:
        !* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
        !* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)

        if (hecMAT%Iarray(97) .gt. 1) hecMAT%Iarray(97)=1

        call hecmw_mat_set_flag_converged(hecMAT, 0)
        call hecmw_mat_set_flag_diverged(hecMAT, 0)

        if (hecMAT%Iarray(2) .eq. 102) then
          if(hecMESH%PETOT.GT.1) then
            call hecmw_solve_direct_ClusterMKL(hecMESH, hecMAT)
          else
            call hecmw_solve_direct_MKL(hecMESH,hecMAT)
          endif
        elseif (hecMAT%Iarray(2) .eq. 104) then
          call hecmw_solve_direct_MUMPS(hecMESH, hecMAT)
        else
          if(hecMESH%PETOT.GT.1) then
            call hecmw_solve_direct_parallel(hecMESH,hecMAT,imsg)
          else
            call hecmw_solve_direct(hecMESH,hecMAT,imsg)
          endif
        endif

        resid=hecmw_rel_resid_L2(hecMESH,hecMAT)
        myrank=hecmw_comm_get_rank()
        if (myrank==0) then
          if (hecMAT%Iarray(21) > 0 .or. hecMAT%Iarray(22) > 0) then
            write(*,"(a,1pe12.5)")'### Relative residual =', resid
          endif
          if( resid >= 1.0d-8) then
            write(*,"(a)")'### Relative residual exceeded 1.0d-8---Direct Solver### '
            !            stop
          endif
        endif
        if (resid < hecmw_mat_get_resid(hecMAT)) then
          call hecmw_mat_set_flag_converged(hecMAT, 1)
        endif
        !C
    end select

  end subroutine hecmw_solve

  subroutine hecmw_substitute_solver(hecMESH, hecMATorig, NDOF)

    use hecmw_util
    use hecmw_solver_iterative
    use hecmw_solver_direct
    use hecmw_solver_direct_parallel
    use hecmw_solver_direct_MUMPS
    use hecmw_solver_direct_clusterMKL
    use hecmw_matrix_contact
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
    !  select case(NDOF)
    !    case(1)
    !      call hecmw_solve_11(hecMESH,hecMAT)
    !    case(2)
    !      call hecmw_solve_22(hecMESH,hecMAT)
    !    case(3)
    !      call hecmw_solve_33(hecMESH,hecMAT)
    !    case(4)
    !      call hecmw_solve_iterative(hecMESH,hecMAT)
    !    case(5)
    !      call hecmw_solve_iterative(hecMESH,hecMAT)
    !    case(6)
    !      call hecmw_solve_66(hecMESH,hecMAT)
    !    case(7:)
    !      call hecmw_solve_iterative(hecMESH,hecMAT)
    !  end select
    !call hecmw_solve_direct_MUMPS(hecMESH, hecMAT)
    call hecmw_solve_iterative(hecMESH,hecMAT)
    if (NDOF /= hecMATorig%NDOF) then
      call hecmw_vector_contract(hecMATorig,hecMAT,NDOF)
    end if
  end subroutine hecmw_substitute_solver

end module hecmw_solver
