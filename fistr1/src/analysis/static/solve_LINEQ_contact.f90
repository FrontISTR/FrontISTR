!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : m_solve_LINEQ_contact                             !
!                                                                      !
!            Written by Z. Sun(ASTOM)                                  !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief This module provides functions to solve sparse system of 
!> \linear equitions in the case of contact analysis using standard
!> \Lagrange multiplier algorithm   
!>
!>  \author     Z. Sun(ASTOM)
!>  \date       2010/11   
!>  \version    0.00
!!
!======================================================================!

module m_solve_LINEQ_contact                                 

   use m_fstr                                                  
   use m_solve_LINEQ_mkl                             
   use m_solve_LINEQ_direct_serial_lag                    
   use m_solve_LINEQ_MUMPS_contact
   use m_fstr_mat_resid_contact
   
  implicit none   
  
  private
  public :: solve_LINEQ_contact_init
  public :: solve_LINEQ_contact

  contains
    
!> \brief This subroutine    
    subroutine solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
      type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
      type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
      type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange)
      logical :: is_sym

      if( hecMAT%Iarray(99)==3 )then
        call solve_LINEQ_mkl_init(hecMAT,fstrMAT,is_sym)
      elseif( hecMAT%Iarray(99)==4 )then
        call solve_LINEQ_serial_lag_hecmw_init(hecMAT,fstrMAT,is_sym)
      elseif( hecMAT%Iarray(99)==5 ) then
        call solve_LINEQ_mumps_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
      endif
    end subroutine solve_LINEQ_contact_init


!> \brief This subroutine
    subroutine solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,rf)
  
      type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
      type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
      type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange)
      real(kind=kreal), optional              :: rf

      real(kind=kreal)                         :: factor
      real(kind=kreal) :: resid

      factor = 1.0d0
      if( present(rf) )factor = rf
         
      if( hecMAT%Iarray(99)==3 )then                   
        call solve_LINEQ_mkl(hecMAT,fstrMAT)
      elseif( hecMAT%Iarray(99)==4 )then
        call solve_LINEQ_serial_lag_hecmw(hecMAT,fstrMAT)
      elseif( hecMAT%Iarray(99)==5 ) then
        call solve_LINEQ_mumps_contact(hecMESH,hecMAT,fstrMAT)
      endif
     
      resid=fstr_get_resid_max_contact(hecMESH,hecMAT,fstrMAT)
      if (myrank==0) then
        write(*,*) ' maximum residual = ', resid
        if( resid >= 1.0d-8) then
          write(*,*) ' ###Maximum residual exceeded 1.0d-8---Direct Solver### '
!          stop
        endif
      endif
      
      hecMAT%X=factor*hecMAT%X

    end subroutine solve_LINEQ_contact  
    
    
end module m_solve_LINEQ_contact
