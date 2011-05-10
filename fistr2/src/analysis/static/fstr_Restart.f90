!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.0                                   !
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
!>  \date       2010/01/14
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Restart
   use m_utilities
   implicit none
!
   contains
   
!> Read in restart file
!----------------------------------------------------------------------*
      subroutine fstr_read_restart(cstep,substep,fstrSOLID)
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(out)                  :: cstep       !< current step
      integer, intent(out)                  :: substep     !< current sub step
      type (fstr_solid),intent(inout)       :: fstrSOLID   !< fstr_solid

      integer :: i,j,fnum,ierr,restrt_step(3), nif(2)
      real(kind=kreal) :: fdum(1)
      character(len=HECMW_FILENAME_LEN):: fname

      fname=trim(restartfilNAME)
      call append_int2name( myrank, fname )
      fnum =60
      OPEN(fnum,status='OLD',file=fname,form='UNFORMATTED',err=1)
	  
      read(fnum) restrt_step
      read(fnum) fstrSOLID%unode
      read(fnum) fstrSOLID%QFORCE

      do i= 1, total_elem
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
	  
      do i=1,restrt_step(3)
        if( .not. associated(fstrSOLID%contacts(i)%states)) cycle
        do j=1,size(fstrSOLID%contacts(i)%states)
          read(fnum) nif
          fstrSOLID%contacts(i)%states(j)%state=nif(1)
          fstrSOLID%contacts(i)%states(j)%surface=nif(2)
          read(fnum) fdum
          fstrSOLID%contacts(i)%states(j)%distance=fdum(1)
          read(fnum) fstrSOLID%contacts(i)%states(j)%lpos
          read(fnum) fstrSOLID%contacts(i)%states(j)%gpos
          read(fnum) fstrSOLID%contacts(i)%states(j)%direction
          read(fnum) fstrSOLID%contacts(i)%states(j)%multiplier
        enddo
      enddo

      close( fnum)
      cstep = restrt_step(1)
      substep = restrt_step(2) + 1
      return
	  
    1 stop "Cannot open restart file"
	  
      end subroutine fstr_read_restart

!> write out restart file
!----------------------------------------------------------------------*
      subroutine fstr_write_restart(cstep,substep,fstrSOLID)
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: cstep      !< current step
      integer, intent(in)                   :: substep    !< current sub step
      type (fstr_solid), intent(in)         :: fstrSOLID  !< fstr_solid

      integer :: i,j,fnum,restrt_step(3)
      integer :: nif(2)
      real(kind=kreal) :: fdum(1)
      character(len=HECMW_FILENAME_LEN):: fname

      fname=trim(restartfilNAME)
      call append_int2name( myrank, fname )
      fnum =60
      OPEN(fnum,status='unknown',file=fname,form='UNFORMATTED')

      restrt_step(1)=cstep
      restrt_step(2)=substep
	  restrt_step(3)=0
      if( associated(fstrSOLID%contacts) ) restrt_step(3)=size(fstrSOLID%contacts)
      write(fnum) restrt_step
      write(fnum) fstrSOLID%unode
      write(fnum) fstrSOLID%QFORCE

      do i= 1, total_elem
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
	  
      do i=1,restrt_step(3)
        if( .not. associated(fstrSOLID%contacts(i)%states)) cycle
        do j=1,size(fstrSOLID%contacts(i)%states)
          nif(1)= fstrSOLID%contacts(i)%states(j)%state
          nif(2)= fstrSOLID%contacts(i)%states(j)%surface
          fdum(1)= fstrSOLID%contacts(i)%states(j)%distance
          write(fnum) nif
          write(fnum) fdum
          write(fnum) fstrSOLID%contacts(i)%states(j)%lpos
          write(fnum) fstrSOLID%contacts(i)%states(j)%gpos
          write(fnum) fstrSOLID%contacts(i)%states(j)%direction
          write(fnum) fstrSOLID%contacts(i)%states(j)%multiplier
        enddo
      enddo
	    
      close(fnum)

      end subroutine fstr_write_restart

end module m_fstr_restart
