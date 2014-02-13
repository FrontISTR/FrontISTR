!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provide a function to prepare output of static analysis
module m_static_make_result

use hecmw_result
implicit none

contains
!>  MAKE RESULT for static analysis(WITHOUT ELEMENTAL RESULTS)
      subroutine FSTR_MAKE_RESULT(fstrSOLID,fstrRESULT)
      use m_fstr
      type (fstr_solid)         :: fstrSOLID
      type (hecmwST_result_data):: fstrRESULT   

      integer :: i,k,ndof	  

      call hecmw_nullify_result_data( fstrRESULT )
      ndof = assDOF(1)                            
      if( NDOF==3  .or. NDOF==2 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(16*total_node))
        fstrRESULT%nn_dof(1)=3
        fstrRESULT%nn_dof(2)=6
        fstrRESULT%nn_dof(3)=7
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'
        fstrRESULT%node_label(3)='STRESS'
  
        do i= 1, total_node
          do k=1,NDOF
            fstrRESULT%node_val_item(16*(i-1)+k)     &
            =fstrSOLID%unode(NDOF*i-NDOF+k)
          enddo
          do k=1,6
            fstrRESULT%node_val_item(16*(i-1)+k+3)   &                  
            =fstrSOLID%STRAIN(6*i-6+k)
          enddo
          do k=1,7
            fstrRESULT%node_val_item(16*(i-1)+k+9)   &
            =fstrSOLID%STRESS(7*i-7+k)
          enddo
        enddo
      
      elseif( NDOF == 6 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(28*total_node))
        fstrRESULT%nn_dof(1)=6
        fstrRESULT%nn_dof(2)=10
        fstrRESULT%nn_dof(3)=12
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'

        do i= 1, total_node
          do k=1,6
            fstrRESULT%node_val_item(28*(i-1)+k)       &          
            =fstrSOLID%unode(6*(i-1)+k)
          enddo
          do k=1,10
            fstrRESULT%node_val_item(28*(i-1)+k+6)     &                
            =fstrSOLID%STRAIN(10*(i-1)+k)
          enddo
          do k=1,12
            fstrRESULT%node_val_item(28*(i-1)+k+16)    &                
            =fstrSOLID%STRESS(12*(i-1)+k)
          enddo
        enddo
      endif

      end subroutine FSTR_MAKE_RESULT


end module m_static_make_result
