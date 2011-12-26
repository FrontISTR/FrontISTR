!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato(Advancesoft), X. YUAN(AdavanceSoft)    !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to read in and write out
!>         restart fiels
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/12/26
!>  \version    0.00
!>  \author     Z. Sun(ASTOM)
!>  \date       2011/11   
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Restart
   use m_utilities
   use m_fstr
   implicit none
   
   integer, parameter :: fnum = 60
   character(len=HECMW_FILENAME_LEN):: fname 

   contains
   
!> Read in restart file
!----------------------------------------------------------------------*
      subroutine fstr_read_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM,contactNode)
!----------------------------------------------------------------------*
      integer, intent(out)                  :: cstep       !< current step
      integer, intent(out)                  :: substep     !< current sub step
      integer, intent(out)                  :: step_count  !< current through step
      integer, intent(out), optional       :: contactNode !< total number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
      type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
      type(fstr_param), intent(in)          :: fstrPARAM     

      integer :: i,j,ierr,restrt_step(3), nif(2)
      integer :: nslave,k
      real(kind=kreal) :: fdum(1)
      character(len=1)::dummy                   

      fname=trim(restartfilNAME)
      call append_int2name( hecMESH%my_rank,fname,fstrSOLID%restart_nin )  
      open(fnum,status='OLD',file=fname,form='UNFORMATTED',err=1)
	  
      read(fnum) restrt_step
      read(fnum) fstrSOLID%unode
      read(fnum) fstrSOLID%QFORCE

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            read(fnum) nif
            read(fnum) fstrSOLID%elements(i)%gausses(j)%strain
            read(fnum) fstrSOLID%elements(i)%gausses(j)%stress
            if( nif(1)>0 )        &
              read(fnum) fstrSOLID%elements(i)%gausses(j)%istatus
            if( nif(2)>0 )        &
              read(fnum) fstrSOLID%elements(i)%gausses(j)%fstatus
        enddo
      enddo
	  
      if(present(contactNode))then   !!associated(fstrSOLID%contacts))then       
        read(fnum)contactNode  
        do i=1,size(fstrSOLID%contacts)
          nslave=size(fstrSOLID%contacts(i)%slave)
          read(fnum)fstrSOLID%contacts(i)%states%surface
          read(fnum)fstrSOLID%contacts(i)%states%state 
          read(fnum)((fstrSOLID%contacts(i)%states(j)%lpos(k),k=1,2),j=1,nslave)
          read(fnum)((fstrSOLID%contacts(i)%states(j)%direction(k),k=1,3),j=1,nslave)
          read(fnum)((fstrSOLID%contacts(i)%states(j)%multiplier(k),k=1,3),j=1,nslave) 
          read(fnum)((fstrSOLID%contacts(i)%states(j)%tangentForce_trial(k),k=1,3),j=1,nslave)                   
          read(fnum)((fstrSOLID%contacts(i)%states(j)%tangentForce_final(k),k=1,3),j=1,nslave)   
        enddo   
      endif  
      cstep = restrt_step(1) 
      substep = restrt_step(2) + 1
      step_count = restrt_step(3)   
      if(fstrPARAM%solution_type==kstNLSTATIC)then      
        if(restrt_step(2)==fstrSOLID%step_ctrl(cstep)%num_substep)then   
          cstep = cstep + 1                                                          
          substep = 1
        endif 
      endif
      read(fnum,end=6)dummy                                                         
      return                                           
    6 close( fnum)                                     
 
      return       
    1 stop "Cannot open restart file"
	  
      end subroutine fstr_read_restart

!> write out restart file
!----------------------------------------------------------------------*
      subroutine fstr_write_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM,contactNode)  
!----------------------------------------------------------------------*
      integer, intent(in)                   :: cstep      !< current step
      integer, intent(in)                   :: substep    !< current sub step
      integer, intent(in)                   :: step_count !< current through step
      integer, intent(in),optional         :: contactNode   !< number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
      type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
      type(fstr_param), intent(in)          :: fstrPARAM           

      integer :: i,j,restrt_step(3)
      integer :: nif(2)
      integer :: nslave,k
      real(kind=kreal) :: fdum(1)

      fname=trim(restartfilNAME)                       
      
      if(fstrPARAM%restart_out_Type==restart_outLast)then                            
        call append_int2name( hecMESH%my_rank,fname )  
      elseif(fstrPARAM%restart_out_Type==restart_outAll)then
        call append_int2name( hecMESH%my_rank,fname,step_count )  
      endif                                                                                               
      open(fnum,status='unknown',file=fname,form='UNFORMATTED')        

      restrt_step(1)=cstep
      restrt_step(2)=substep
	  restrt_step(3)=step_count
      write(fnum) restrt_step
      write(fnum) fstrSOLID%unode
      write(fnum) fstrSOLID%QFORCE

      do i= 1, hecMESH%n_elem
        do j= 1, size(fstrSOLID%elements(i)%gausses)
            nif = 0
            if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) ) nif(1)=size(fstrSOLID%elements(i)%gausses(j)%istatus)
            if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) ) nif(2)=size(fstrSOLID%elements(i)%gausses(j)%fstatus)
            write(fnum) nif
            write(fnum) fstrSOLID%elements(i)%gausses(j)%strain
            write(fnum) fstrSOLID%elements(i)%gausses(j)%stress
            if( associated(fstrSOLID%elements(i)%gausses(j)%istatus) )        &
              write(fnum) fstrSOLID%elements(i)%gausses(j)%istatus
            if( associated(fstrSOLID%elements(i)%gausses(j)%fstatus) )        &
               write(fnum) fstrSOLID%elements(i)%gausses(j)%fstatus
        enddo
      enddo
      
      if(present(contactNode))then   
        write(fnum)contactNode
        do i=1,size(fstrSOLID%contacts)
          nslave=size(fstrSOLID%contacts(i)%slave)
          write(fnum)fstrSOLID%contacts(i)%states%surface
          write(fnum)fstrSOLID%contacts(i)%states%state 
          write(fnum)((fstrSOLID%contacts(i)%states(j)%lpos(k),k=1,2),j=1,nslave)
          write(fnum)((fstrSOLID%contacts(i)%states(j)%direction(k),k=1,3),j=1,nslave)
          write(fnum)((fstrSOLID%contacts(i)%states(j)%multiplier(k),k=1,3),j=1,nslave)
          write(fnum)((fstrSOLID%contacts(i)%states(j)%tangentForce_trial(k),k=1,3),j=1,nslave)      
          write(fnum)((fstrSOLID%contacts(i)%states(j)%tangentForce_final(k),k=1,3),j=1,nslave)     
        enddo   
      endif   
	  
      close(fnum)

      end subroutine fstr_write_restart
	  
!> Read in restart file for dynamic analysis
!----------------------------------------------------------------------*
      subroutine fstr_read_restart_dyna(cstep,substep,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,contactNode)
      integer, intent(out)                  :: cstep       !< current step
      integer, intent(out)                  :: substep     !< current sub step
      integer                               :: step_count  !< current through step
      integer, intent(out), optional      :: contactNode !< number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH     !< hecmw mesh
      type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid
      type ( fstr_dynamic), intent(inout)   :: fstrDYNAMIC
      type(fstr_param), intent(in)          :: fstrPARAM  
      integer :: i                                     
	  
	  if(present(contactNode))then
        call fstr_read_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM,contactNode)
      else
        call fstr_read_restart(cstep,substep,step_count,hecMESH,fstrSOLID,fstrPARAM)  
      endif
      backspace(fnum) 
      read(fnum) fstrDYNAMIC%t_curr                                                   
      read(fnum) fstrDYNAMIC%idx_eqa
      if( fstrDYNAMIC%idx_eqa == 1 ) then
        read(fnum) fstrDYNAMIC%DISP(:,1)
        read(fnum) fstrDYNAMIC%VEL(:,1)
        read(fnum) fstrDYNAMIC%ACC(:,1)
      else	
        read(fnum) fstrDYNAMIC%DISP(:,1)
        read(fnum) fstrDYNAMIC%DISP(:,3)
      endif
      read(fnum)fstrDYNAMIC%strainEnergy                 
      do i= 1, hecMESH%n_elem                                                             
        read(fnum)fstrSOLID%elements(i)%equiForces                    
      enddo  
      cstep=cstep+1                                                
      close(fnum)
      end subroutine
	  
!> write out restart file for dynamic analysis
!----------------------------------------------------------------------*
      subroutine fstr_write_restart_dyna(cstep,substep,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,contactNode) 
!----------------------------------------------------------------------*
      integer, intent(in)                   :: cstep      !< current step
      integer, intent(in)                   :: substep    !< current sub step
      integer                               :: step_count !< current through step
      integer, intent(in), optional       :: contactNode !< number of contact nodes
      type (hecmwST_local_mesh), intent(in) :: hecMESH    !< hecmw mesh
      type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid
      type ( fstr_dynamic), intent(in)      :: fstrDYNAMIC
      type(fstr_param), intent(in)          :: fstrPARAM            
      integer :: i                                                   
	  
      if(present(contactNode))then
        call fstr_write_restart(cstep,1,cstep,hecMESH,fstrSOLID,fstrPARAM,contactNode) 
      else
        call fstr_write_restart(cstep,1,cstep,hecMESH,fstrSOLID,fstrPARAM)   
      endif  
      OPEN(fnum,status='OLD',file=fname,form='UNFORMATTED',position='APPEND')
      write(fnum) fstrDYNAMIC%t_curr                   
      write(fnum) fstrDYNAMIC%idx_eqa
      if( fstrDYNAMIC%idx_eqa == 1 ) then
        write(fnum) fstrDYNAMIC%DISP(:,2)
        write(fnum) fstrDYNAMIC%VEL(:,2)
        write(fnum) fstrDYNAMIC%ACC(:,2)
      else
        write(fnum) fstrDYNAMIC%DISP(:,1)
        write(fnum) fstrDYNAMIC%DISP(:,3)
      endif
      write(fnum)fstrDYNAMIC%strainEnergy                        
      do i= 1, hecMESH%n_elem                                                                 
        write(fnum)fstrSOLID%elements(i)%equiForces                    
      enddo                                                                                                 
      close(fnum)
      end subroutine

end module m_fstr_restart
