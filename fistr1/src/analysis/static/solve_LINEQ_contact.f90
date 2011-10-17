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
   use m_set_arrays_directsolver_contact              
   use m_solve_LINEQ_mkl                             
   use m_solve_LINEQ_direct_serial_lag                    
   
  implicit none   
  
  contains
    
!> \brief This subroutine    
    subroutine solve_LINEQ_contact(hecMAT,fstrMAT)          
  
      type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
      type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange)
   
      integer(kind=kint)                       :: ntdf   
      real(kind=kreal), allocatable            :: y(:)           !< right-hand side vector
      real(kind=kreal)                          :: residual_Max   !< maximum residual
     
      ntdf = hecMAT%NP*hecMAT%NDOF + fstrMAT%num_lagrange 
         
      if( hecMAT%Iarray(99)==3 )then                   
        call solve_LINEQ_mkl(hecMAT,fstrMAT)
      elseif( hecMAT%Iarray(99)==4 )then
        call solve_LINEQ_serial_lag_hecmw(hecMAT,fstrMAT)
      endif
     
      allocate(y(size(hecMAT%B)))
      y = 0.0d0  
      call getApproximateB(ntdf,hecMAT%X,y)                                     
      residual_Max=MAXVAL(DABS(y-hecMAT%B))  
      write(*,*)' maximum residual = ',residual_Max                   
      if( dabs(residual_Max) >= 1.0d-8) then                          
        write(*,*) ' ###Maximum residual exceeded 1.0d-8---Direct Solver### '
!        stop 
      endif                                                  
      deallocate(y)
      
      deallocate(values)
   
    end subroutine solve_LINEQ_contact  
    
    
!> \brief This subroutine gets the residual vector  
    subroutine getApproximateB(ntdf,x,y)          
   
     integer(kind=kint)     :: ntdf                       !< total degree of freedom
     integer(kind=kint)     :: i, j, k                    
     real(kind=kreal)        :: x(ntdf)                    !< solution vector
     real(kind=kreal)        :: y(ntdf)                    !< residual vector
     
       y = 0.0d0   
       do i = 1, ntdf
         do j = pointers(i), pointers(i+1)-1
           k = indices(j)
           y(i) = y(i) + values(j)*x(k)
           if( symmetricMatrixStruc .and. k/=i )&         
           y(k)=y(k)+values(j)*x(i)   
         enddo           
       enddo 
                 
     end subroutine getApproximateB   
     
   
   
end module m_solve_LINEQ_contact
