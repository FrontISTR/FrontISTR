!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines controlling dynamic calculation

module fstr_solver_dynamic

use m_fstr
use fstr_dynamic_implicit
use fstr_dynamic_explicit

contains

!C================================================================C
!> Master subroutine for dynamic analysis
!C================================================================C
  subroutine fstr_solve_DYNAMIC(hecMESH,hecMAT,fstrSOLID,myEIG   &
                          ,fstrDYNAMIC,fstrPARAM,fstrCPL )
      use m_fstr_setup
      implicit none
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      type ( fstr_couple         ) :: fstrCPL         !for COUPLE


      if(dabs(fstrDYNAMIC%t_delta) .lt. 1.0e-20) then
              if( hecMESH%my_rank .eq. 0 ) then
                 write(imsg,*) 'stop due to fstrDYNAMIC%t_delta = 0'
              end if
              call hecmw_abort( hecmw_comm_get_comm())
      end if
      call fstr_dynamic_alloc( fstrDYNAMIC )

      if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                 call fstr_solve_dynamic_implicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL )
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if

      else if(fstrDYNAMIC%idx_eqa == 11) then  ! explicit dynamic analysis

             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                call fstr_solve_dynamic_explicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL )

             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if
      end if
		 
      call fstr_dynamic_finalize( fstrDYNAMIC )
      call hecMAT_finalize( imsg, hecMAT )

  end subroutine fstr_solve_DYNAMIC


end module fstr_solver_dynamic
