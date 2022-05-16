!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module blopex_fortran_hold_vars
  use m_fstr
  use hecmw_util
  use iso_c_binding
  integer(c_int) :: N_hold, M_hold
  type(hecmwST_local_mesh), pointer :: hecMESH_hold
  type(hecmwST_matrix), pointer     :: hecMAT_hold
  type(fstr_eigen), pointer         :: fstrEIG_hold
end module blopex_fortran_hold_vars

module m_fstr_EIG_lobpcg
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none

contains

  subroutine fstr_solve_lobpcg(hecMESH, hecMAT, fstrSOLID, fstrEIG)
    use m_fstr
    use hecmw_util
    use hecmw_precond
    implicit none
    type(hecmwST_local_mesh), target :: hecMESH
    type(hecmwST_matrix), target     :: hecMAT
    type(fstr_solid)                 :: fstrSOLID
    type(fstr_eigen), target         :: fstrEIG
    integer(kint) :: N, n_eigs, loglevel, maxit
    real(kreal) :: tol

    if(myrank == 0)then
      write(IMSG,'(" fstr_solve_lobpcg")')
      write(*,'(" fstr_solve_lobpcg")')
    endif

    N = hecMAT%NP*hecMAT%NDOF
    n_eigs = fstrEIG%nget
    maxit = fstrEIG%maxiter
    tol = fstrEIG%tolerance
    N_hold = N
    M_hold = n_eigs
    hecMESH_hold => hecMESH
    hecMAT_hold => hecMAT
    fstrEIG_hold => fstrEIG

    call hecmw_precond_setup(hecMAT, hecMESH, 1)

    allocate(fstrEIG%eigval(n_eigs), source = 0.0d0)
    allocate(fstrEIG%eigvec(N, n_eigs), source = 0.0d0)

    loglevel = 1

    call blopex_lobpcg_solve(N, n_eigs, maxit, tol, loglevel, &
      & fstrEIG%eigval, fstrEIG%eigvec)
  end subroutine fstr_solve_lobpcg

  subroutine blopex_lobpcg_solve( &
      & N, n_eigs, maxit, tol, loglevel, eigen_val, eigen_vec)
    implicit none
    integer(c_int) :: N, NZ, i, j
    integer(c_int) :: is_B, is_T
    integer(c_int) :: n_eigs, maxit, loglevel
    real(c_double) :: tol
    real(c_double) :: eigen_val(n_eigs)
    real(c_double) :: eigen_vec(N,n_eigs)
    real(c_double) :: eigen_vec_temp(N*n_eigs)
    external :: blopex_fortran_opA
    external :: blopex_fortran_opB
    external :: blopex_fortran_opT

    if(N < n_eigs) n_eigs = N

    is_B = 1
    is_T = 0

    call blopex_lobpcg_solve_c(n_eigs, maxit, tol, N, &
      & blopex_fortran_opA, blopex_fortran_opB, blopex_fortran_opT, &
      & is_B, is_T, &
      & loglevel, eigen_val, eigen_vec_temp)

    do i = 1, n_eigs
      do j = 1, N
        eigen_vec(j,i) = eigen_vec_temp(n_eigs*(i-1) + j)
      enddo
    enddo
  end subroutine blopex_lobpcg_solve

end module m_fstr_EIG_lobpcg

subroutine blopex_fortran_opA(dum, a, b)
  use blopex_fortran_hold_vars
  use hecmw_solver_las
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    call hecmw_matvec(hecMESH_hold, hecMAT_hold, a(jS:jE), b(jS:jE))
  enddo
end subroutine blopex_fortran_opA

subroutine blopex_fortran_opB(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    do i = jS, jE
      b(i) = fstrEIG_hold%mass(i-jS+1)*a(i)
    enddo
  enddo
end subroutine blopex_fortran_opB

subroutine blopex_fortran_opT(dum, a, b)
  use blopex_fortran_hold_vars
  use iso_c_binding
  implicit none
  integer(c_int) :: dum, i, j, jS, jE
  real(c_double) :: a(N_hold*M_hold), b(N_hold*M_hold)

  do j = 1, M_hold
    jS = N_hold*(j-1) + 1
    jE = N_hold*(j)
    !call hecmw_precond_apply(hecMESH, hecMAT, WW(:,R), WW(:,Z), WW(:,WK), Tcomm)
  enddo
end subroutine blopex_fortran_opT
