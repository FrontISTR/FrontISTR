!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : m_dynamic_init_variables                          !
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

!C================================================================C
!> \brief This module provides functions to initialize variables 
!> when initial velocity or acceleration boundary conditions are given.
!> attention : just for rigid motion in the initial state.   
!C================================================================C

module m_dynamic_init_variables                                  

use m_fstr
use lczparm
use m_dynamic_mat_ass_load

contains

   subroutine dynamic_init_varibles( hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC )
   
    implicit none

    type ( hecmwST_local_mesh  ) :: hecMESH
    type ( hecmwST_matrix      ) :: hecMAT
    type ( lczparam            ) :: myEIG
    type ( fstr_solid          ) :: fstrSOLID
    type ( fstr_dynamic        ) :: fstrDYNAMIC
    
    integer(kind=kint) :: j    
     
    call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
    
    if( fstrSOLID%VELOCITY_type == kbcInitial ) then   
      do j = 1, hecMESH%n_node*hecMESH%n_dof  
        fstrDYNAMIC%ACC(j,1)=(hecMAT%B(j)-fstrDYNAMIC%ray_m*myEIG%mass(j)*fstrDYNAMIC%VEL(j,1))/&         
                              myEIG%mass(j) 
      enddo
    elseif( fstrSOLID%ACCELERATION_type == kbcInitial ) then
      do j = 1, hecMESH%n_node*hecMESH%n_dof
        fstrDYNAMIC%VEL(j,1)=(hecMAT%B(j)-myEIG%mass(j)*fstrDYNAMIC%ACC(j,1))/&         
                             (fstrDYNAMIC%ray_m*myEIG%mass(j))              
      enddo    
    endif  

   end subroutine dynamic_init_varibles

end module m_dynamic_init_variables
