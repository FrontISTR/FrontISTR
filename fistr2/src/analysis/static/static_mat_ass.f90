!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Xi YUAN (Advancesoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!C
!C***
!> Matrix assemble for FSTR solver
!C***
!C
module m_static_mat_ass
   contains

   SUBROUTINE fstr_mat_ass (myEIG, fstrSOLID)

      USE m_fstr
      USE lczparm
      use m_static_lib
      use m_static_mat_ass_main
      use m_fstr_ass_load
      use m_fstr_AddBC

      implicit none
      type (fstr_solid)  :: fstrSOLID
      type(lczparam)     :: myEIG


      IF(myrank == 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' ****   STAGE Stiffness Matrix assembly  **'
      ENDIF


      call fstr_mat_ass_main (fstrSOLID)
      write(IDBG,*) 'fstr_mat_ass_main: OK'

      IF(myEIG%eqset==0) THEN
        call fstr_ass_load (1, fstrSOLID,1.d0)
        write(IDBG,*) 'fstr_mat_ass_load: OK'

      ELSE IF(myEIG%eqset==1) THEN
        IF(myrank == 0) THEN
          WRITE(IMSG,*) '*-------------------------------------------*'
          WRITE(IMSG,*) 'NOTE: Loads ignored for eigenvalue analysis.' 
          WRITE(IMSG,*) '*-------------------------------------------*'
        ENDIF

      ENDIF

      call fstr_AddBC(1,fstrSOLID,1)
      write(IDBG,*) 'fstr_mat_ass_bc: OK'

      end subroutine FSTR_MAT_ASS
	  
end module m_static_mat_ass
