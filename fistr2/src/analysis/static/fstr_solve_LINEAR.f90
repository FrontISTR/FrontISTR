!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.0                                   !
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

!> SOLVE STATIC LINEAR SOLID MECHANICS
module m_fstr_solve_LINEAR
   implicit none
   contains

   subroutine FSTR_SOLVE_LINEAR( myEIG,fstrSOLID )
      use m_fstr
      use m_fstr_NodalStress
      use m_static_make_result
      use m_static_mat_ass
      use m_static_output
      use m_static_post
      use m_fstr_Update
      use hecmw_result

      include "HEC_MW3_For.h"

      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID

      type ( hecmwST_result_data ) :: fstrRESULT
      integer :: level, partID
      character(len=HECMW_FILENAME_LEN) :: ctrlfile

      integer :: iAss, iPart, iElem, iNode, iGrp, ndID, iErr, ik, ndof
      integer :: snode, enode
	  
      ndof = assDOF(1)
!C
!C-- MATRIX ASSEMBLING
!C
      call fstr_mat_ass(myEIG, fstrSOLID)

      IF(myrank == 0) THEN
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
      iErr=mw_solve(svIarray(1), svRarray(1), svIarray(2), svIarray(3))
      if( iErr==0 ) then
        write(*,*) ' Fatal error: Fails in solver! '
        write(IMSG,*) ' Fatal error: Fails in solver! '
        call hecmw_abort(hecmw_comm_get_comm())
      endif
!C
!C-- OUTPUT 
!C
!      call mw_send_recv()
	  
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            snode = part_nodes(iAss+1,iPart+1)+1
            enode = part_nodes(iAss+1,iPart+2)
            call mw_get_solution_vector(fstrSOLID%unode((snode-1)*ndof+1:enode*ndof), iPart)
         enddo
      enddo
      call fstr_UpdateNewton( fstrSOLID, 1,0.d0,1.d0,1)

      call fstr_NodalStress3D( fstrSOLID, fstrSOLID%STRAIN, fstrSOLID%STRESS )

      call solid_output(fstrSOLID)

      IF(myrank == 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- POST PROCESSING
!C
      call fstr_post(fstrSOLID,1)

!C
!C-- POST PROCESSING VIA MEMORY
!C
      if( IVISUAL==1 ) then
        level = 0
        partID = -1
        ctrlfile = 'hecmw_ctrl.dat'
        call fstr_make_result(fstrSOLID,fstrRESULT)
        call mw3_visualize_init_ex(ctrlfile,level,partID)
        call hecmw_result_copy_f2c(fstrRESULT,iErr)
        if( iErr/=0 ) then
          write(*,*) ' Fatal error: Fails in copying data for visualizer! '
          write(IMSG,*) ' Fatal error: Fails in copying data for visualizer! '
          call hecmw_abort(hecmw_comm_get_comm())
        endif
        call mw3_visualize(1,1,0)
        call mw3_visualize_finalize
        call hecmw_result_free(fstrRESULT)
      endif

      end subroutine FSTR_SOLVE_LINEAR
end module m_fstr_solve_LINEAR
