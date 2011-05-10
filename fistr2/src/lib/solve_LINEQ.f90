!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.0                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Xi YUAN (Advancesoft)                          !
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

   use iso_c_binding
   use m_fstr

   interface
     integer(C_INT) function mw_solve_( iter_max, tolerance, method, pre_condition ) bind(c)
	    use iso_c_binding
        implicit none
        type(C_PTR), VALUE :: iter_max
        type(C_PTR), VALUE :: tolerance
        type(C_PTR), VALUE :: method
        type(C_PTR), VALUE :: pre_condition
	 END FUNCTION mw_solve_
   end interface
	 
   contains

   SUBROUTINE solve_LINEQ(imsg, ierr)
     integer, intent(in) :: imsg
	 integer(C_INT) :: ierr
     integer(C_INT), target :: iter, med, pre
     real(C_DOUBLE), target :: tol
     tol = svRarray(1)
     ierr = mw_solve_(c_loc(iter), c_loc(tol), c_loc(med), c_loc(pre))
   end subroutine solve_LINEQ
   
end module m_solve_LINEQ
