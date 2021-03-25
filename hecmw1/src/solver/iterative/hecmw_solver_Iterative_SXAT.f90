!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_iterative_SXAT
!C***
!C
module hecmw_solver_iterative_SXAT

    public :: hecmw_solve_sxat_test_fortran
    public :: hecmw_solve_sxat_entry
    public :: hecmw_solve_sxat_CG_entry

contains
    subroutine hecmw_solve_sxat_CG_interface(N, & !0
            &                               NP, &  !1
            &                               NDOF, & !2
            &                               D, & !3
            &                               AL, & !4
            &                               AU, & !5
            &                               indexL, & !6
            &                               indexU, & !7
            &                               itemL, & !8
            &                               itemU, & !9
            &                               my_rank, & !10
            &                               X, & !11
            &                               B, & !12
            &                               iterlog, & !13
            &                               timelog, & !14
            &                               estcond, & !15
            &                               iter, & !16
            &                               err, & !17
            &                               time_setup, & !18
            &                               time_sol, & !19
            &                               time_comm) bind(c) !20

        use hecmw_util
        use m_hecmw_solve_error
        use m_hecmw_comm_f
        use hecmw_matrix_misc
        use hecmw_solver_misc
        use hecmw_solver_las
        use hecmw_solver_scaling
        use hecmw_precond
        use hecmw_jad_type
        use hecmw_estimate_condition
        use hecmw_solver_CG
        use hecmw_solver_CG
        use hecmw_solver_las
        use hecmw_precond
        use hecmw_solver_misc
        use Iso_C_Binding

        implicit none

        integer(kind=kint), value :: N,NP,NDOF
        real(kind=kreal), allocatable :: D(:), AL(:), AU(:)
        real(kind=kreal), allocatable :: B(:), X(:)
        real(kind=kint), allocatable :: itemL(:), itemU(:), indexL(:), indexU(:)
        integer(kind=kint), value :: iter, error
        integer(kind=kint), value :: iterlog, timelog, 
        integer(kind=kint) ::ESTCOND
        integer(kind=kint ) :: my_rank
        integer(kind=kint), value :: iter, error
        real(kind=kreal), value :: resid, time_setup, time_sol, time_comm

        type(hecmwST_matrix) :: hecMAT
        type (hecmwST_local_mesh) :: hecMESH
        hecMAT%N = N
        hecMAT%NP = NP
        hecMAT%NDOF = NDOF
        hecMAT%D => D
        hecMAT%AL => D
        hecMAT%AU => D
        hecMAT%indexL => indexL
        hecMAT%indexU => indexU
        hecMAT%itemL => itemL
        hecMAT%itemU => itemU
        hecMAT%X => X
        hecMAT%B => B
        hecMESH%my_rank = my_rank

!       print *, 'VE fortran CG now'
        call hecmw_solve_CG(hecMESH, hecMAT, iter, resid, error, Tset, Tsol, Tcomm)

    end subroutine hecmw_solve_sxat_CG_interface

!     subroutine hecmw_solve_sxat_test_fortran(val, N, array) bind(c)
!         use hecmw_util
!         use m_hecmw_solve_error
!         use m_hecmw_comm_f
!         use hecmw_matrix_misc
!         use hecmw_solver_misc
!         use hecmw_solver_las
!         use hecmw_solver_scaling
!         use hecmw_precond
!         use hecmw_jad_type
!         use hecmw_estimate_condition
!         use Iso_C_Binding
! 
!         integer(kind=kint), value :: N
!         real(kind=kreal), value :: val
!         real(kind=kreal) :: array(N)
!         type (hecmwST_local_mesh) :: hecMESH
!         Integer :: ndof
! 
!         print *, 'hello I am VE fortran'
!         print *, N, val, array(1), array(2)
! 
!         hecMESH%nn_internal = N
!         ndof = 1
! 
!         print *, array(1), array(2)
! 
!         array(1) = array(1) + 10000
! 
!         call hecmw_scale_R(hecMESH, ndof, val, array)
! 
! 
!     end subroutine hecmw_solve_sxat_test_fortran
! 
!     integer function hecmw_solve_sxat_entry( X, ndof, addr)
! 
!         use hecmw_util
!         use m_hecmw_solve_error
!         use m_hecmw_comm_f
!         use hecmw_matrix_misc
!         use hecmw_solver_misc
!         use hecmw_solver_las
!         use hecmw_solver_scaling
!         use hecmw_precond
!         use hecmw_jad_type
!         use hecmw_estimate_condition
!         use Iso_C_Binding
! 
! 
!         implicit none
!         integer, value :: x, ndof
! 
! !         type(c_ptr) :: addr
! !         real, pointer :: array(:)
!         integer(8), value :: addr
! 
!         !call c_f_pointer(addr, array, [100])
! 
!         print *, X, ndof, addr
! 
!         hecmw_solve_sxat_entry = 0
!     end function hecmw_solve_sxat_entry
end module hecmw_solver_iterative_SXAT
