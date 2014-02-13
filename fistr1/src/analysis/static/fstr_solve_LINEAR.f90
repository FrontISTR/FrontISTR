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

!> SOLVE STATIC LINEAR SOLID MECHANICS
module m_fstr_solve_LINEAR
   contains

   subroutine FSTR_SOLVE_LINEAR ( hecMESH,hecMAT,myEIG,fstrSOLID,fstrPARAM )

      use m_fstr
      use m_static_mat_ass
      use hecmw_solver_11
      use hecmw_solver_22
      use hecmw_solver_33
      use hecmw_solver_direct
      use m_solve_lineq
      use m_fstr_Update
      use m_static_output
      use m_hecmw2fstr_mesh_conv

      implicit REAL(kind=kreal) (A-H,O-Z)
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM

      integer(kind=kint) :: i

      if( fstrSOLID%TEMP_ngrp_tot>0 .and. hecMESH%hecmw_flag_initcon==1 ) then
        do j=1, hecMESH%n_node
          fstrSOLID%last_temp(j) = hecMESH%node_init_val_item(j)
          fstrSOLID%temperature(j) = hecMESH%node_init_val_item(j)
        end do
      endif
!C
!C-- MATRIX ASSEMBLING
!C
      call fstr_mat_ass(hecMESH, hecMAT, myEIG, fstrSOLID)

      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *****  STAGE Solve static problem    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- LINEAR SOLVER
!C
      hecMAT%Iarray(98) = 1   !Assmebly complete
      hecMAT%Iarray(97) = 1   !Need numerical factorization
      CALL solve_LINEQ(hecMESH,hecMAT,imsg)

      IF(myrank .EQ. 0) THEN
        IF(hecMAT%Iarray(99) .EQ. 1) THEN
          write(IMSG,*) '*----------------------------*'
          write(IMSG,*) '## No.of ITER:',hecMAT%Iarray(1)
          write(IMSG,*) '*----------------------------*'
          call flush(IMSG)
        ENDIF
      ENDIF

!C
!C-- UPDATE DISPLACEMENT, STRAIN, STRESS
!C
      do i = 1, hecMESH%n_node*hecMESH%n_dof
        fstrSOLID%unode(i) = hecMAT%X(i)
      enddo
      if( hecMESH%n_dof==3 ) then
        call hecmw_update_3_R ( hecMESH, fstrSOLID%unode, hecMAT%NP )
        call fstr_Update3D ( hecMESH, fstrSOLID )
      else if( hecMESH%n_dof==2 ) then
        call hecmw_update_2_R ( hecMESH, fstrSOLID%unode, hecMAT%NP )
        call fstr_Update2D ( hecMESH, fstrSOLID )
      else if( hecMESH%n_dof==6) THEN
        call hecmw_update_m_R ( hecMESH, fstrSOLID%unode, hecMAT%NP, hecMESH%n_dof )
        call fstr_Update6D ( hecMESH, fstrSOLID )
      endif 

      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- POST PROCESSING
!C
      call fstr_static_Output( 1, 1, hecMESH, fstrSOLID, fstrPR%solution_type )

      end subroutine FSTR_SOLVE_LINEAR
end module m_fstr_solve_LINEAR
