!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.0                                   !
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
      use m_static_make_result
      use m_static_output
      use m_static_post
      use m_hecmw2fstr_mesh_conv

      implicit REAL(kind=kreal) (A-H,O-Z)
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM

      integer(kind=kint)    :: nstep, istep
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
!C-- OUTPUT FOR NODES
      if( hecMESH%n_dof==3 ) then
        call hecmw_update_3_R ( hecMESH, hecMAT%X, hecMAT%NP )
        call fstr_output_3d(hecMESH,hecMAT,fstrSOLID,fstrPARAM)
      else if( hecMESH%n_dof==2 ) then
        call hecmw_update_2_R (hecMESH, hecMAT%X, hecMAT%NP)
        call fstr_output_2d(hecMESH,hecMAT,fstrSOLID,fstrPARAM)
      else if( hecMESH%n_dof==6) THEN
        call hecmw_update_m_R (hecMESH,hecMAT%X,hecMAT%NP,hecMESH%n_dof)
        call fstr_output_6d(hecMESH,hecMAT,fstrSOLID,fstrPARAM)
      endif 
	  
!C-- OUTPUT FOR ELEMENTS
      call fstr_output_elem(hecMESH,fstrSOLID)

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
      nstep = 1
      istep = 1
      call fstr_post(hecMESH,hecMAT,fstrSOLID,istep)
!C
!C-- POST PROCESSING VIA MEMORY
!C
      if( IVISUAL==1 ) then
        call fstr_make_result(hecMESH,fstrSOLID,fstrRESULT)
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        istep =0
        idummy=0
        call hecmw_visualize(hecMESH,fstrRESULT,istep,istep,idummy)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif

      end subroutine FSTR_SOLVE_LINEAR
end module m_fstr_solve_LINEAR
