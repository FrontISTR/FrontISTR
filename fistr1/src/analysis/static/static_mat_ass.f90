!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!C
!C***
!> Matrix assemble for FSTR solver
!C***
!C
module m_static_mat_ass
   contains

   SUBROUTINE fstr_mat_ass (hecMESH, hecMAT, myEIG, fstrSOLID)

      USE m_fstr
      USE lczparm
      use m_static_lib
      use m_static_mat_ass_main
      use m_fstr_ass_load
      use m_fstr_AddBC
      use fstr_matrix_con_contact                           

      implicit none
      integer(kind=kint) :: IFLAG, numnp, ndof, i
      real(kind=kreal) :: bsize
      data IFLAG/0/

      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (fstr_param       )  :: fstrPARAM          
      type (fstrST_matrix_contact_lagrange)  :: fstrMAT                                                  
!* LCZ
      type(lczparam) :: myEIG

      INTEGER(kind=kint) :: n1, n2

      iexit = 0
      IF(myrank == 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' ****   STAGE Stiffness Matrix assembly  **'
      ENDIF

      hecMAT%NDOF = hecMESH%n_dof

      if (IFLAG==0) then
!        CALL memget(masbr,n2*hecMAT%NPL,8)
!        CALL memget(masbr,n2*hecMAT%NPU,8)
!        CALL memget(masbr,(n1*2+n2)*hecMAT%NP,8)
!        CALL memget(masbr,n2*hecMAT%N,8)

        IFLAG= 1
      endif

      fstrSOLID%factor(1)=0.d0; fstrSOLID%factor(2)=1.d0
      call fstr_mat_ass_main (hecMESH, hecMAT, fstrSOLID)
      call fstr_AddSPRING(1, 1, hecMESH, hecMAT, fstrSOLID, fstrPARAM)
      write(IDBG,*) 'fstr_mat_ass_main: OK'

      IF(myEIG%eqset==0) THEN
        call fstr_ass_load (1, hecMESH, hecMAT, fstrSOLID, fstrPARAM)
        write(IDBG,*) 'fstr_mat_ass_load: OK'

      ELSE IF(myEIG%eqset==1) THEN
        IF(myrank == 0) THEN
          WRITE(IMSG,*) '*-------------------------------------------*'
          WRITE(IMSG,*) 'NOTE: Loads ignored for eigenvalue analysis.' 
          WRITE(IMSG,*) '*-------------------------------------------*'
        ENDIF

      ENDIF

      call fstr_AddBC(1, 1, hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrMAT,1)
      write(IDBG,*) 'fstr_mat_ass_bc: OK'
!C
!C  RHS LOAD VECTOR CHECK
!C
      IF(myEIG%eqset==0) THEN
        numnp=hecMAT%NP
        ndof =hecMESH%n_dof
        bsize=0.0
        do i=1,numnp*ndof
          bsize=bsize+hecMAT%B(i)**2
        enddo
!C
!C Gather RHS vector
!C
        call hecmw_allREDUCE_R1( hecMESH,bsize,hecmw_sum )

        if( hecMESH%my_rank == 0 ) then
          write(IMSG,*) 'Total RHS size=',bsize
          write(IMSG,*) 'Total number of equations=',hecMESH%mpc%n_mpc
        endif

        if( bsize < 1.0e-31 .and. hecMESH%mpc%n_mpc==0 ) then
          iexit = 1
          WRITE(IMSG,*) '###Load Vector Error!'
          stop
        endif

      ENDIF

      return
      end subroutine FSTR_MAT_ASS
end module m_static_mat_ass
