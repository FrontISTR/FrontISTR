!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Subroutine Name : pardiso                                       !
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
!> \brief This module provides functions to stop contact calculations  
!> \with Lagrage multipliers. 
!>
!>  \author     Z. Sun(ASTOM)
!>  \date       2010/11   
!>  \version    0.00
!!
!======================================================================!
 
        
   subroutine pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices, &                                 
                       idum, nrhs, iparm, msglvl, ddum1, ddum2, ierr)

   use m_fstr

!< ------------------------------------------------ input parameters for MKL solver(see Intel(R) MKL Reference Manual)
   integer (kind=8)              :: pt(64)                 
   integer (kind=kint)           :: maxfct
   integer (kind=kint)           :: mnum
   integer (kind=kint)           :: mtype
   integer (kind=kint)           :: phase
   integer (kind=kint)           :: ntdf
   integer (kind=kint)           :: pointers(:)        !< ia     
   integer (kind=kint)           :: indices(:)         !< ja    
   real     (kind=kreal)          :: values(:)          !< a    
   integer (kind=kint)           :: idum
   integer (kind=kint)           :: nrhs
   integer (kind=kint)           :: iparm(64)
   integer (kind=kint)           :: msglvl   
   real(kind=kreal)               :: ddum1, ddum2
   integer (kind=kint)           :: ierr
   
   write(*,*) "paradiso: ERROR was detected. Please install MKL library and try again."
   stop
  
   end subroutine
