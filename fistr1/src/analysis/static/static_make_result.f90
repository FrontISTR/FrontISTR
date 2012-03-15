!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
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
contains
!C***
!>  MAKE RESULT for static analysis(WITHOUT ELEMENTAL RESULTS)
!C***
      subroutine FSTR_MAKE_RESULT(hecMESH,fstrSOLID,fstrRESULT)
      use m_fstr
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (hecmwST_result_data):: fstrRESULT      
!C***
      call hecmw_nullify_result_data( fstrRESULT )
      NDOF = hecMESH%n_dof                             
      if( NDOF.eq.3  .or. NDOF.eq.2 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(16*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=3
        fstrRESULT%nn_dof(2)=6
        fstrRESULT%nn_dof(3)=7
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'
        fstrRESULT%node_label(3)='STRESS'
!C*** Set Displacement/Strain/Stress @node
        do i= 1, hecMESH%n_node
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
      
      elseif( NDOF .eq. 6 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(28*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=6
        fstrRESULT%nn_dof(2)=10
        fstrRESULT%nn_dof(3)=12
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'
        fstrRESULT%node_label(3)='STRESS'
!C*** Set Displacement/Strain/Stress @node
        do i= 1, hecMESH%n_node
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
	  
	  subroutine FSTR_MAKE_EIGENRESULT(hecMESH,unode,fstrRESULT)
      use m_fstr
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal)          :: unode(:)
      type (hecmwST_result_data):: fstrRESULT      
!C***
      call hecmw_nullify_result_data( fstrRESULT )
      NDOF = hecMESH%n_dof                             
      if( NDOF==3  .or. NDOF==2 ) then
        fstrRESULT%nn_component=1
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(1) )
        allocate( fstrRESULT%node_label(1) )
        allocate( fstrRESULT%node_val_item(3*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=3
        fstrRESULT%node_label(1)='DISPLACEMENT'
        do i= 1, hecMESH%n_node
          do k=1,NDOF
            fstrRESULT%node_val_item(3*(i-1)+k)     &
            =unode(NDOF*i-NDOF+k)
          enddo
        enddo
      
      elseif( NDOF == 6 ) then
        fstrRESULT%nn_component=1
        fstrRESULT%ne_component=0
        allocate( fstrRESULT%nn_dof(1) )
        allocate( fstrRESULT%node_label(1) )
        allocate( fstrRESULT%node_val_item(6*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=6
        fstrRESULT%node_label(1)='DISPLACEMENT'
        do i= 1, hecMESH%n_node
          do k=1,6
            fstrRESULT%node_val_item(6*(i-1)+k)       &          
            =unode(6*(i-1)+k)
          enddo
        enddo
      endif

      end subroutine FSTR_MAKE_EIGENRESULT

!C***
!>  FSTR MAKE RESULT for FUTURE VERTION (WITH ELEMENTAL RESULTS)
!C***
      subroutine FSTR_MAKE_RESULT_FUTURE(hecMESH,fstrSOLID,fstrRESULT)
      use m_fstr
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (hecmwST_result_data):: fstrRESULT
!C***
      NDOF = hecMESH%n_dof
      if( NDOF.eq.3  .or. NDOF.eq.2 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=2
!C @node
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(16*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=3
        fstrRESULT%nn_dof(2)=6
        fstrRESULT%nn_dof(3)=7
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'
        fstrRESULT%node_label(3)='STRESS'
!C @element
        allocate( fstrRESULT%ne_dof(2) )
        allocate( fstrRESULT%elem_label(2) )
        allocate( fstrRESULT%elem_val_item(13*hecMESH%n_elem))
        fstrRESULT%ne_dof(1)=6
        fstrRESULT%ne_dof(2)=7
        fstrRESULT%elem_label(1)='ESTRAIN'
        fstrRESULT%elem_label(2)='ESTRESS'
        fstrRESULT%elem_val_item=0.0
!C*** Set Displacement/Strain/Stress 
!C @node
        do i= 1, hecMESH%n_node
          do k=1,NDOF
            fstrRESULT%node_val_item(16*(i-1)+k)    &
            =fstrSOLID%unode(NDOF*i-NDOF+k)
          enddo
          do k=1,6
            fstrRESULT%node_val_item(16*(i-1)+k+3)  &                   
            =fstrSOLID%STRAIN(6*i-6+k)
          enddo
          do k=1,7
            fstrRESULT%node_val_item(16*(i-1)+k+9)  &        
            =fstrSOLID%STRESS(7*i-7+k)
          enddo
        enddo
!C @element
        do i= 1, hecMESH%n_elem
!!!          ielem   = hecMESH%elem_ID( i*2-1 )
          ielem   = i
          ID_area = hecMESH%elem_ID( i*2 )
          if( ID_area /= hecMESH%my_rank ) cycle
          do k=1,6
            fstrRESULT%elem_val_item(13*(i-1)+k)     &                
            =fstrSOLID%ESTRAIN(6*ielem-6+k)
          enddo
          do k=1,7
            fstrRESULT%elem_val_item(13*(i-1)+k+6)   &                 
            =fstrSOLID%ESTRESS(7*ielem-7+k)
          enddo
        enddo
      
      elseif( NDOF .eq. 6 ) then
        fstrRESULT%nn_component=3
        fstrRESULT%ne_component=2
!C @node
        allocate( fstrRESULT%nn_dof(3) )
        allocate( fstrRESULT%node_label(3) )
        allocate( fstrRESULT%node_val_item(28*hecMESH%n_node))
        fstrRESULT%nn_dof(1)=6
        fstrRESULT%nn_dof(2)=10
        fstrRESULT%nn_dof(3)=12
        fstrRESULT%node_label(1)='DISPLACEMENT'
        fstrRESULT%node_label(2)='STRAIN'
        fstrRESULT%node_label(3)='STRESS'
!C @elem
        allocate( fstrRESULT%ne_dof(2) )
        allocate( fstrRESULT%elem_label(2) )
        allocate( fstrRESULT%elem_val_item(22*hecMESH%n_elem))
        fstrRESULT%ne_dof(1)=10
        fstrRESULT%ne_dof(2)=12
        fstrRESULT%elem_label(1)='ESTRAIN'
        fstrRESULT%elem_label(2)='ESTRESS'
        fstrRESULT%elem_val_item=0.0
!C*** Set Displacement/Strain/Stress 
!C @node
        do i= 1, hecMESH%n_node
          do k=1,6
            fstrRESULT%node_val_item(28*(i-1)+k)     &
            =fstrSOLID%unode(6*(i-1)+k)
          enddo
          do k=1,10
            fstrRESULT%node_val_item(28*(i-1)+k+6)   &                  
            =fstrSOLID%STRAIN(10*(i-1)+k)
          enddo
          do k=1,12
            fstrRESULT%node_val_item(28*(i-1)+k+16)  &            
            =fstrSOLID%STRESS(12*(i-1)+k)
          enddo
        enddo
!C @element
        do i= 1, hecMESH%n_elem
!!!          ielem   = hecMESH%elem_ID( i*2-1 )
          ielem   = i
          ID_area = hecMESH%elem_ID( i*2 )
          if( ID_area == hecMESH%my_rank ) then
            do k=1,10
              fstrRESULT%elem_val_item(22*(i-1)+k)      &
              =fstrSOLID%ESTRAIN(10*(ielem-1)+k)
            enddo
            do k=1,12
              fstrRESULT%elem_val_item(22*(i-1)+k+10)   &
              =fstrSOLID%ESTRESS(12*(ielem-1)+k)
            enddo
          else
            do k=1,10
              fstrRESULT%elem_val_item(22*(i-1)+k) = 0.0
            enddo
            do k=1,12
              fstrRESULT%elem_val_item(22*(i-1)+k+10) = 0.0
            enddo
          endif
        enddo
      endif

      end subroutine FSTR_MAKE_RESULT_FUTURE
end module m_static_make_result
