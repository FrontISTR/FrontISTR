!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : fstr_matrix_con_contact                           !
!                                                                      !
!      Written by Z. Sun(ASTOM)                                        !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief This module provides functions of reconstructing 
!         stiffness matrix structure for the contact analysis 
!         employing standard Lagrange multiplier algorithm
!
!>  \author     Z. Sun(ASTOM)
!>  \date       2010/11
!>  \version    0.00
!!
!======================================================================!
module fstr_matrix_con_contact     

      use m_fstr
      use elementInfo 
      use mContact                                                                       

      implicit none

!> Structure for defining stiffness matrix structure
      type nodeRelated
        integer (kind=kint)             :: num_node=0, num_lagrange=0       !< total number of related nodes and Lagrange multipliers          
        integer (kind=kint), pointer   :: id_node(:) => NULL()         !< list of related nodes
        integer (kind=kint), pointer   :: id_lagrange(:) => NULL()     !< list of related Lagrange multipliers           
      end type

!> Structure for Lagrange multiplier-related part of stiffness matrix
!> (Lagrange multiplier-related matrix) 
      type fstrST_matrix_contact_lagrange
        integer (kind=kint)            :: num_lagrange=0                               !< total number of Lagrange multipliers
        integer (kind=kint)            :: numL_lagrange=0, numU_lagrange=0               !< node-based number of non-zero items in lower triangular half of matrix
                                                                                      !< node-based number of non-zero items in upper triangular half of matrix
        integer (kind=kint), pointer  :: indexL_lagrange(:) => NULL(), &     
                                           indexU_lagrange(:) => NULL()               !< node-based index of first non-zero item of each row 
        integer (kind=kint), pointer  :: itemL_lagrange(:) => NULL(), &      
                                           itemU_lagrange(:) => NULL()                !< node-based column number of non-zero items
        real (kind=kreal),    pointer  :: AL_lagrange(:) => NULL(), &         
                                           AU_lagrange(:) => NULL()                   !< values of non-zero items  
        real (kind=kreal),    pointer  :: Lagrange(:) => NULL()                      !< values of Lagrange multipliers                    
      end type fstrST_matrix_contact_lagrange
   
      integer (kind=kint),          save    :: NPL_org, NPU_org                      !< original number of non-zero items    
      type (nodeRelated), pointer, save    :: list_nodeRelated_org(:) => null()     !< original structure of matrix   
                  
      type (nodeRelated), pointer          :: list_nodeRelated(:) => null()          !< current structure of matrix
  
      logical                               :: permission = .false.               
  
      private                               :: insert_lagrange, insert_node, find_locationINarray, reallocate_memory
   
 contains 
   

!> \brief This subroutine saves original matrix structure constructed originally by hecMW_matrix    
    subroutine fstr_save_originalMatrixStructure(hecMAT)
       
      type (hecmwST_matrix)           :: hecMAT                                !< type hecmwST_matrix
   
      integer (kind=kint)             :: numL, numU, num_nodeRelated           !< original number of nodes related to each node
      integer (kind=kint)             :: i, j     
      integer (kind=kint)             :: ierr
  
       NPL_org = hecMAT%NPL;  NPU_org = hecMAT%NPU
   
      allocate(list_nodeRelated_org(hecMAT%NP), stat=ierr)   
       if( ierr /= 0) stop " Allocation error, list_nodeRelated_org "

      do i = 1, hecMAT%NP
       
        numL = hecMAT%indexL(i) - hecMAT%indexL(i-1) 
        numU = hecMAT%indexU(i) - hecMAT%indexU(i-1)

        num_nodeRelated = numL + numU + 1      
                             
        allocate(list_nodeRelated_org(i)%id_node(num_nodeRelated), stat=ierr)
          if( ierr /= 0) stop " Allocation error, list_nodeRelated_org%id_node "      
            
        list_nodeRelated_org(i)%num_node = num_nodeRelated
                  
        do j = 1, numL
          list_nodeRelated_org(i)%id_node(j) = hecMAT%itemL(hecMAT%indexL(i-1)+j) 
        enddo
        list_nodeRelated_org(i)%id_node(numL+1) = i
        do j = 1, numU         
          list_nodeRelated_org(i)%id_node(numL+1+j) = hecMAT%itemU(hecMAT%indexU(i-1)+j)   
        enddo
           
      enddo    
                 
    end subroutine fstr_save_originalMatrixStructure
 
!> \brief this subroutine reconstructs node-based (stiffness) matrix structure 
!> \corresponding to contact state       
    subroutine fstr_mat_con_contact(cstep,hecMAT,fstrSOLID,fstrMAT,infoCTChange,conMAT)    
  
      integer (kind=kint)                     :: cstep                                  !< current loading step     
      type (hecmwST_matrix)                    :: hecMAT                                !< type hecmwST_matrix
      type (fstr_solid)                        :: fstrSOLID                             !< type fstr_solid
      type (fstrST_matrix_contact_lagrange)    :: fstrMAT                               !< type fstrST_matrix_contact_lagrange 
      type (fstr_info_contactChange)           :: infoCTChange                          !< type fstr_contactChange           

      integer (kind=kint)                     :: num_lagrange                           !< number of Lagrange multipliers
      integer (kind=kint)                     :: countNon0LU_node, countNon0LU_lagrange !< counter of node-based number of non-zero items
      integer (kind=kint)                     :: numNon0_node, numNon0_lagrange         !< node-based number of displacement-related non-zero items in half of the matrix
                                                                                         !< node-based number of Lagrange multiplier-related non-zero items in half of the matrix  
      type (hecmwST_matrix),optional          :: conMAT 

      num_lagrange = infoCTChange%contactNode_current                      
      fstrMAT%num_lagrange = num_lagrange  
      
! Get original list of related nodes                                 
      call getOriginalListOFrelatedNodes(hecMAT%NP,num_lagrange)     
       
! Construct new list of related nodes and Lagrange multipliers
      countNon0LU_node = NPL_org + NPU_org
      countNon0LU_lagrange = 0 
      if( fstr_is_contact_active() ) &                            
      call getNewListOFrelatednodesANDLagrangeMultipliers(cstep,hecMAT%NP,fstrSOLID,countNon0LU_node,countNon0LU_lagrange)      

! Construct new matrix structure(hecMAT&fstrMAT)
      numNon0_node = countNon0LU_node/2                                     
      numNon0_lagrange = countNon0LU_lagrange/2
!     ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        call constructNewMatrixStructure(hecMAT,fstrMAT,numNon0_node,numNon0_lagrange,conMAT)
      else
        call constructNewMatrixStructure(hecMAT,fstrMAT,numNon0_node,numNon0_lagrange)
      endif
         
! Copy Lagrange multipliers
      if( fstr_is_contact_active() ) &                                 
      call fstr_copy_lagrange_contact(fstrSOLID,fstrMAT)

    end subroutine fstr_mat_con_contact

!> Get original list of related nodes
    subroutine getOriginalListOFrelatedNodes(np,num_lagrange)
    
      integer (kind=kint)            :: np, num_lagrange                    !< total number of nodes
      integer (kind=kint)            :: num_nodeRelated_org                 !< original number of nodes related to each node
      integer (kind=kint)            :: i, ierr   
            
      allocate(list_nodeRelated(np+num_lagrange), stat=ierr)  
      if( ierr /= 0) stop " Allocation error, list_nodeRelated "
   
      do i = 1, np  !hecMAT%NP    
        num_nodeRelated_org = list_nodeRelated_org(i)%num_node
        allocate(list_nodeRelated(i)%id_node(num_nodeRelated_org), stat=ierr) 
        if( ierr /= 0) stop " Allocation error, list_nodeRelated%id_node "
        list_nodeRelated(i)%num_node = num_nodeRelated_org                                                
        list_nodeRelated(i)%id_node(1:num_nodeRelated_org) = list_nodeRelated_org(i)%id_node(1:num_nodeRelated_org)
      enddo
    
      if( fstr_is_contact_active() ) then                
        do i = np+1, np+num_lagrange  !hecMAT%NP+1, hecMAT%NP+num_lagrange
          allocate(list_nodeRelated(i)%id_lagrange(5), stat=ierr)                       
          if( ierr /= 0) stop " Allocation error, list_nodeRelated%id_lagrange "
          list_nodeRelated(i)%num_lagrange = 0                                  
          list_nodeRelated(i)%id_lagrange = 0   
        enddo  
      endif                                      
    
    end subroutine getOriginalListOFrelatedNodes
    

!> Construct new list of related nodes and Lagrange multipliers. Here, a procedure similar to HEC_MW is used.      
    subroutine getNewListOFrelatednodesANDLagrangeMultipliers(cstep,np,fstrSOLID,countNon0LU_node,countNon0LU_lagrange)   
   
      type (fstr_solid)              :: fstrSOLID                                           !< type fstr_solid
      
      integer (kind=kint)            :: cstep                                               !< current loading step 
      integer (kind=kint)            :: np                                                  !< total number of nodes 
      integer (kind=kint)            :: countNon0LU_node, countNon0LU_lagrange              !< counters of node-based number of non-zero items
      integer (kind=kint)            :: grpid                                               !< contact pairs group ID     
      integer (kind=kint)            :: count_lagrange                                      !< counter of Lagrange multiplier
      integer (kind=kint)            :: ctsurf, etype, nnode, ndLocal(l_max_surface_node+1) !< contants of type tContact       
      integer (kind=kint)            :: i, j, k, l, num, num_nodeRelated_org, ierr 
      real (kind=kreal)               :: fcoeff                                              !< friction coefficient
    
      count_lagrange = 0       
      do i = 1, size(fstrSOLID%contacts)                        
       
        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle                  
                                                           
        fcoeff = fstrSOLID%contacts(i)%fcoeff          
          
        do j = 1, size(fstrSOLID%contacts(i)%slave)
        
          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle  
          ctsurf = fstrSOLID%contacts(i)%states(j)%surface         
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          if( etype/=fe_tri3n .and. etype/=fe_quad4n ) &                                                     
          stop " ##Error: This element type is not supported in contact analysis !!! "   
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)
                       
          count_lagrange = count_lagrange + 1

          do k = 1, nnode+1
          
            if( .not. associated(list_nodeRelated(ndLocal(k))%id_lagrange) )then
              num = 10                                                             
 !             if( k == 1 ) num = 1                              
              allocate(list_nodeRelated(ndLocal(k))%id_lagrange(num),stat=ierr)
              if( ierr /= 0) stop " Allocation error, list_nodeRelated%id_lagrange "
              list_nodeRelated(ndLocal(k))%num_lagrange = 0
              list_nodeRelated(ndLocal(k))%id_lagrange = 0 
            endif
             
            if( fcoeff /= 0.0d0 ) then  
              num_nodeRelated_org = list_nodeRelated_org(ndLocal(k))%num_node                   
              if( list_nodeRelated(ndLocal(k))%num_node == num_nodeRelated_org )then   
                num = 10                    
                if(k==1) num = 4                                                                        
                call reallocate_memory(num,list_nodeRelated(ndLocal(k)))   
              endif
            endif

            call insert_lagrange(k,count_lagrange,list_nodeRelated(ndLocal(k)),countNon0LU_lagrange)  
 
            do l = k, nnode+1                                                           
              if( fcoeff /= 0.0d0 ) then                         
                if( k /= l) then 
                  num_nodeRelated_org = list_nodeRelated_org(ndLocal(k))%num_node     
                  call insert_node(ndLocal(l),list_nodeRelated(ndLocal(k)),countNon0LU_node)  
                  num_nodeRelated_org = list_nodeRelated_org(ndLocal(l))%num_node        
                  if( list_nodeRelated(ndLocal(l))%num_node == num_nodeRelated_org )then    
                    num = 10                                                 
                    call reallocate_memory(num,list_nodeRelated(ndLocal(l)))  
                  endif
                  call insert_node(ndLocal(k),list_nodeRelated(ndLocal(l)),countNon0LU_node)    
                endif
              endif                       

              if(k == 1) &             
              call insert_lagrange(0,ndLocal(l),list_nodeRelated(np+count_lagrange),countNon0LU_lagrange)         
               
            enddo                 
            
          enddo                     
          
        enddo     
        
      enddo
    
    end subroutine getNewListOFrelatednodesANDLagrangeMultipliers

     
!> Construct new stiffness matrix structure    
    subroutine constructNewMatrixStructure(hecMAT,fstrMAT,numNon0_node,numNon0_lagrange,conMAT)
    
      type (hecmwST_matrix)                    :: hecMAT                                  !< type hecmwST_matrix
      type (fstrST_matrix_contact_lagrange)    :: fstrMAT                                 !< type fstrST_matrix_contact_lagrange
      
      integer (kind=kint)                      :: numNon0_node, numNon0_lagrange          !< node-based number of non-zero items in half of the matrix
      integer (kind=kint)                      :: countNon0L_node, countNon0U_node, countNon0U_lagrange, countNon0L_lagrange  !< counters of node-based number ofnon-zero items
      integer (kind=kint)                      :: i, j, ierr  
      integer (kind=kint)                      :: numI_node, numI_lagrange 
      integer (kind=kint)                      :: ndof, nn
      type (hecmwST_matrix),optional           :: conMAT

!     ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        conMAT%N  = hecMAT%N
        conMAT%NP = hecMAT%NP
        conMAT%ndof = hecMAT%ndof
        if(associated(conMAT%indexL).and.associated(conMAT%indexU))deallocate(conMAT%indexL,conMAT%indexU) 
        allocate(conMAT%indexL(0:conMAT%NP), conMAT%indexU(0:conMAT%NP), stat=ierr) 
        if ( ierr /= 0) stop " Allocation error, conMAT%indexL-conMAT%indexU "
        conMAT%indexL = 0 ; conMAT%indexU = 0
        if(associated(conMAT%itemL).and.associated(conMAT%itemU))deallocate(conMAT%itemL,conMAT%itemU) 
        allocate(conMAT%itemL(numNon0_node), conMAT%itemU(numNon0_node), stat=ierr)
        if ( ierr /= 0) stop " Allocation error, conMAT%itemL-conMAT%itemU "
        conMAT%itemL = 0 ; conMAT%itemU = 0
!
        conMAT%NPL = numNon0_node
        conMAT%NPU = numNon0_node
      endif
    
      if(associated(hecMAT%indexL).and.associated(hecMAT%indexU))deallocate(hecMAT%indexL,hecMAT%indexU) 
      allocate(hecMAT%indexL(0:hecMAT%NP), hecMAT%indexU(0:hecMAT%NP), stat=ierr) 
      if ( ierr /= 0) stop " Allocation error, hecMAT%indexL-hecMAT%indexU "
      hecMAT%indexL = 0 ; hecMAT%indexU = 0
      if(associated(hecMAT%itemL).and.associated(hecMAT%itemU))deallocate(hecMAT%itemL,hecMAT%itemU) 
      allocate(hecMAT%itemL(numNon0_node), hecMAT%itemU(numNon0_node), stat=ierr)
      if ( ierr /= 0) stop " Allocation error, hecMAT%itemL-hecMAT%itemU "
      hecMAT%itemL = 0 ; hecMAT%itemU = 0  
          
      if(associated(fstrMAT%indexL_lagrange).and.associated(fstrMAT%indexU_lagrange)) &
                                                deallocate(fstrMAT%indexL_lagrange,fstrMAT%indexU_lagrange)
      if(associated(fstrMAT%itemL_lagrange).and.associated(fstrMAT%itemU_lagrange)) &
                                                deallocate(fstrMAT%itemL_lagrange,fstrMAT%itemU_lagrange)      
      if( fstr_is_contact_active() ) then                                                                       
        allocate(fstrMAT%indexL_lagrange(0:fstrMAT%num_lagrange), fstrMAT%indexU_lagrange(0:hecMAT%NP), stat=ierr) 
        if ( ierr /= 0) stop " Allocation error, fstrMAT%indexL_lagrange-fstrMAT%indexU_lagrange "
        fstrMAT%indexL_lagrange = 0 ; fstrMAT%indexU_lagrange = 0  
        allocate(fstrMAT%itemL_lagrange(numNon0_lagrange), fstrMAT%itemU_lagrange(numNon0_lagrange), stat=ierr) 
        if ( ierr /= 0) stop " Allocation error, fstrMAT%itemL_lagrange-fstrMAT%itemU_lagrange "
        fstrMAT%itemL_lagrange = 0 ; fstrMAT%itemU_lagrange = 0   
      endif                                                              
              
      hecMAT%NPL = numNon0_node                                
      hecMAT%NPU = numNon0_node                                 
      
      fstrMAT%numL_lagrange = numNon0_lagrange
      fstrMAT%numU_lagrange = numNon0_lagrange

      countNon0L_node = 0 
      countNon0U_node = 0
      countNon0U_lagrange = 0     
      do i = 1, hecMAT%NP                                          
          
        list_nodeRelated(i)%num_node = count(list_nodeRelated(i)%id_node /= 0)        
        numI_node = list_nodeRelated(i)%num_node
        if( fstr_is_contact_active() ) &                                     
        numI_lagrange = list_nodeRelated(i)%num_lagrange
        
        do j = 1, numI_node 
          if( list_nodeRelated(i)%id_node(j) < i ) then 
            countNon0L_node = countNon0L_node + 1
            hecMAT%itemL(countNon0L_node) = list_nodeRelated(i)%id_node(j)      
          elseif( list_nodeRelated(i)%id_node(j) > i ) then
            countNon0U_node = countNon0U_node + 1
            hecMAT%itemU(countNon0U_node) = list_nodeRelated(i)%id_node(j)          
          endif   
        enddo
        hecMAT%indexL(i) = countNon0L_node
        hecMAT%indexU(i) = countNon0U_node         
      
        if( fstr_is_contact_active() ) then                                 
          do j = 1, numI_lagrange
            countNon0U_lagrange = countNon0U_lagrange + 1 
            fstrMAT%itemU_lagrange(countNon0U_lagrange) = list_nodeRelated(i)%id_lagrange(j)
          enddo          
          fstrMAT%indexU_lagrange(i) = countNon0U_lagrange  
        endif                                                                                                                          
       
        deallocate(list_nodeRelated(i)%id_node)  
        if(associated(list_nodeRelated(i)%id_lagrange)) deallocate(list_nodeRelated(i)%id_lagrange)   
      
      end do
      
!     ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        conMAT%itemL(:)   = hecMAT%itemL(:) 
        conMAT%indexL(:)  = hecMAT%indexL(:)
        conMAT%itemU(:)   = hecMAT%itemU(:) 
        conMAT%indexU(:)  = hecMAT%indexU(:)
      endif  
  
      if( fstr_is_contact_active() ) then     
        countNon0L_lagrange = 0                    
        do i = 1, fstrMAT%num_lagrange            
          numI_lagrange = list_nodeRelated(hecMAT%NP+i)%num_lagrange     
          do j = 1, numI_lagrange   
            countNon0L_lagrange = countNon0L_lagrange + 1
            fstrMAT%itemL_lagrange(countNon0L_lagrange) = list_nodeRelated(hecMAT%NP+i)%id_lagrange(j) 
          enddo  
          fstrMAT%indexL_lagrange(i) = countNon0L_lagrange  
          deallocate(list_nodeRelated(hecMAT%NP+i)%id_lagrange)    
        enddo 
      endif                                   

      deallocate(list_nodeRelated)
    
      ndof = hecMAT%NDOF
      nn = ndof*ndof
      if(associated(hecMAT%AL)) deallocate(hecMAT%AL) 
      allocate(hecMAT%AL(nn*hecMAT%NPL), stat=ierr) 
      if ( ierr /= 0 ) stop " Allocation error, hecMAT%AL "
      hecMAT%AL = 0.0D0  

      if(associated(hecMAT%AU)) deallocate(hecMAT%AU) 
      allocate(hecMAT%AU(nn*hecMAT%NPU), stat=ierr) 
      if ( ierr /= 0 ) stop " Allocation error, hecMAT%AU "
      hecMAT%AU = 0.0D0   
      
      if(associated(fstrMAT%AL_lagrange)) deallocate(fstrMAT%AL_lagrange)
      if(associated(fstrMAT%AU_lagrange)) deallocate(fstrMAT%AU_lagrange)
      if(associated(fstrMAT%Lagrange)) deallocate(fstrMAT%Lagrange)
                
      if( fstr_is_contact_active() ) then                 
        allocate(fstrMAT%AL_lagrange(ndof*fstrMAT%numL_lagrange), stat=ierr)
        if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AL_lagrange "
        fstrMAT%AL_lagrange = 0.0D0      
        allocate(fstrMAT%AU_lagrange(ndof*fstrMAT%numU_lagrange), stat=ierr)
        if ( ierr /= 0 ) stop " Allocation error, fstrMAT%AU_lagrange "
        fstrMAT%AU_lagrange = 0.0D0      
        allocate(fstrMAT%Lagrange(fstrMAT%num_lagrange))      
        fstrMAT%Lagrange = 0.0D0 
      endif                                                                                                                     

      if(associated(hecMAT%B)) deallocate(hecMAT%B)  
      allocate(hecMAT%B(hecMAT%NP*ndof+fstrMAT%num_lagrange))      
      hecMAT%B = 0.0D0                                                              

      if(associated(hecMAT%X)) deallocate(hecMAT%X) 
      allocate(hecMAT%X(hecMAT%NP*ndof+fstrMAT%num_lagrange))      
      hecMAT%X = 0.0D0

      if(associated(hecMAT%D)) deallocate(hecMAT%D) 
      allocate(hecMAT%D(hecMAT%NP*ndof**2+fstrMAT%num_lagrange))      
      hecMAT%D = 0.0D0
!
!     ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        if(associated(conMAT%AL)) deallocate(conMAT%AL) 
        allocate(conMAT%AL(nn*conMAT%NPL), stat=ierr) 
        if ( ierr /= 0 ) stop " Allocation error, conMAT%AL "
        conMAT%AL = 0.0D0  
    
        if(associated(conMAT%AU)) deallocate(conMAT%AU) 
        allocate(conMAT%AU(nn*conMAT%NPU), stat=ierr) 
        if ( ierr /= 0 ) stop " Allocation error, conMAT%AU "
        conMAT%AU = 0.0D0   
    
        if(associated(conMAT%B)) deallocate(conMAT%B)  
        allocate(conMAT%B(conMAT%NP*ndof+fstrMAT%num_lagrange))      
        conMAT%B = 0.0D0                                                              
    
        if(associated(conMAT%X)) deallocate(conMAT%X) 
        allocate(conMAT%X(conMAT%NP*ndof+fstrMAT%num_lagrange))      
        conMAT%X = 0.0D0
        
        if(associated(conMAT%D)) deallocate(conMAT%D) 
        allocate(conMAT%D(conMAT%NP*ndof**2+fstrMAT%num_lagrange))      
        conMAT%D = 0.0D0
      endif   
    
    end subroutine ConstructNewMatrixStructure
    

!> Insert a Lagrange multiplier in list of related Lagrange multipliers
    subroutine insert_lagrange(i,id_lagrange,list_node,countNon0_lagrange)      
   
      type (nodeRelated)             :: list_node                    !< type nodeRelated
      integer (kind=kint)            :: i, id_lagrange               !< local number of node in current contact pair
                                                                      !< Lagrange multiplier ID
      integer (kind=kint)            :: countNon0_lagrange           !< counter of node-based number of non-zero items 
                                                                      !< in Lagrange multiplier-related matrix
      integer (kind=kint)            :: ierr, num_lagrange, location
      integer (kind=kint)            :: id_lagrange_save(1000)  
      
      character(len=1)               :: answer                                               

      ierr = 0  
      
      num_lagrange = count(list_node%id_lagrange /= 0 )                   

!      if( i == 1 .and. num_lagrange /= 0) return                                       
      if( i == 1 .and. num_lagrange /= 0 .and. .not. permission) then                  
     1  write(*,*) '##Error: node is both slave and master node simultaneously !'
        write(*,*) '         Please check contact surface definition !'
        write(*,'(''          Do you want to continue(y/n)):'',$)')
        read(*,'(A1)',err=1) answer 
        if(answer == 'Y' .OR. answer == 'y')then  
          permission = .true.          
        else
          stop
        endif
      endif                                                                            
        
      if (num_lagrange == 0)then
        list_node%num_lagrange = 1
        list_node%id_lagrange(1) = id_lagrange 
        countNon0_lagrange = countNon0_lagrange + 1  
      else 
        id_lagrange_save(1:num_lagrange) = list_node%id_lagrange(1:num_lagrange)
        location = find_locationINarray(id_lagrange,num_lagrange,list_node%id_lagrange)
        if(location /= 0)then
          num_lagrange = num_lagrange + 1
          if( num_lagrange > size(list_node%id_lagrange)) then
            deallocate(list_node%id_lagrange)                                         
            allocate(list_node%id_lagrange(num_lagrange),stat=ierr)                    
            if( ierr /= 0 ) stop " Allocation error, list_nodeRelated%id_lagrange "     
          endif    
          list_node%num_lagrange = num_lagrange 
          list_node%id_lagrange(location) = id_lagrange
          if(location /= 1) list_node%id_lagrange(1:location-1) = id_lagrange_save(1:location-1)       
          if(location /= num_lagrange) list_node%id_lagrange(location+1:num_lagrange) = id_lagrange_save(location:num_lagrange-1)
          countNon0_lagrange = countNon0_lagrange + 1
        endif
      endif  

    end subroutine insert_lagrange
   
!> Insert a node in list of related nodes
    subroutine insert_node(id_node,list_node,countNon0_node)    
   
      type (nodeRelated)             :: list_node                        !< type nodeRelated
      integer (kind=kint)            :: id_node                          !< local number of node in current contact pair 
                                                                          !< global number of node 
      integer (kind=kint)            :: countNon0_node                   !< counter of node-based number of non-zero items in displacement-related matrix
      integer (kind=kint)            :: ierr, num_node, location
      integer (kind=kint)            :: id_node_save(1000)                                  
     
      ierr = 0
     
      num_node = list_node%num_node 
      
      id_node_save(1:num_node) = list_node%id_node(1:num_node)  
      location = find_locationINarray(id_node,num_node,list_node%id_node)
      if(location /= 0)then
        num_node = num_node + 1
        if( num_node > size(list_node%id_node)) then    
           deallocate(list_node%id_node)                                             
           allocate(list_node%id_node(num_node),stat=ierr)                            
           if( ierr /= 0) stop " Allocation error, list_nodeRelated%id_node "                 
        endif  
        list_node%num_node = num_node        
        list_node%id_node(location) = id_node
        if(location /= 1) list_node%id_node(1:location-1) = id_node_save(1:location-1)       
        if(location /= num_node) list_node%id_node(location+1:num_node) = id_node_save(location:num_node-1)
        countNon0_node = countNon0_node + 1
      endif      
      
    end subroutine insert_node
   

!> Find location of an item in an array by bisection method   
    integer function find_locationINarray(item,n,a)
 
      integer (kind=kint)            :: item, n           !< item to be found; length of array
      integer (kind=kint), pointer  :: a(:)              !< array               
      integer (kind=kint)            :: l, r, m
     
      find_locationINarray = 0
       
      l = 1 ; r = n ; m = (l+r)/2 
      if( item == a(l) .or. item == a(r) )then
        return
      elseif( item < a(l) )then
        find_locationINarray = 1
        return
      elseif( item > a(r) )then
        find_locationINarray = n + 1
        return
      endif 
        
      do while ( l <= r)
        if( item > a(m) ) then
          l = m + 1
          m = (l + r)/2
        elseif( item < a(m) ) then
          r = m - 1
          m = (l + r)/2
        elseif( item == a(m) )then
          return
        endif    
      enddo  
      
      find_locationINarray = m + 1     
      
    end function find_locationINarray
    

!> Reallocate memory for list_relatedNodes    
    subroutine reallocate_memory(num,list_node)                     
  
      type(nodeRelated)       :: list_node            !< type nodeRelated
      integer (kind=kint)     :: num                  !< length to be added 
      integer (kind=kint)     :: num_node_org         !< original number of related nodes 
                                                       !< before reallocation
      integer (kind=kint)     :: id_save(1000)
      integer (kind=kint)     :: ierr    
     
      num_node_org = size(list_node%id_node)
      id_save(1:num_node_org) = list_node%id_node(1:num_node_org) 
      deallocate(list_node%id_node)
      allocate(list_node%id_node(num_node_org+num),stat=ierr)  
      if( ierr /= 0) stop " reAllocation error, list_nodeRelated%id_node "  
      list_node%id_node = 0     
      list_node%id_node(1:num_node_org) = id_save(1:num_node_org)
    
    end subroutine reallocate_memory  
    
    
!> Copy Lagrange multipliers    
    subroutine fstr_copy_lagrange_contact(fstrSOLID,fstrMAT)                 
    
      type(fstr_solid)                        :: fstrSOLID                !< type fstr_solid
      type(fstrST_matrix_contact_lagrange)    :: fstrMAT                  !< fstrST_matrix_contact_lagrange
      integer (kind=kint)                    :: id_lagrange, i, j     
      
      id_lagrange = 0
        
      do i = 1, size(fstrSOLID%contacts)               
        do j = 1, size(fstrSOLID%contacts(i)%slave)
          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle  
          id_lagrange = id_lagrange + 1
          fstrMAT%Lagrange(id_lagrange)=fstrSOLID%contacts(i)%states(j)%multiplier(1)
        enddo
      enddo    
    
    end subroutine fstr_copy_lagrange_contact 
 

!> \brief this function judges whether sitiffness matrix is symmetric or not
    logical function fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)

      type(fstr_solid )        :: fstrSOLID
      type(hecmwST_local_mesh) :: hecMESH
      integer (kind=kint)      :: is_in_contact

      is_in_contact = 0
      if( any(fstrSOLID%contacts(:)%fcoeff /= 0.0d0) )  &
           is_in_contact = 1
      call hecmw_allreduce_I1(hecMESH, is_in_contact, HECMW_MAX)
      fstr_is_matrixStruct_symmetric = (is_in_contact == 0)

    end function fstr_is_matrixStruct_symmetric


end module fstr_matrix_con_contact
