!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains functions to print out calculation settings

module fstr_debug_dump
  use m_fstr

contains

  !> This subroutine prints out global control parameters
  subroutine dump_fstr_global
    implicit none

    write(*,*) 'global parameters dump ***********'
    write(*,*)
    write(*,*) 'IECHO   ',IECHO
    write(*,*) 'IRESULT ',IRESULT
    write(*,*) 'IVISUAL ',IVISUAL
    write(*,*)
    write(*,*) 'for heat ...'
    write(*,*) 'INEUTRAL ', INEUTRAL
    write(*,*) 'IRRES    ', IRRES
    write(*,*) 'IWRES    ', IWRES
    write(*,*) 'NRRES    ', NRRES
    write(*,*) 'NPRINT   ', NPRINT
    write(*,*)
    write(*,*) 'REF_TEMP ', REF_TEMP
    write(*,*)
    write(*,*) 'ANALYSIS CONTROL for NLGEOM and HEAT'
    write(*,*) 'DT     ',DT
    write(*,*) 'ETIME  ',ETIME
    write(*,*) 'ITMAX  ',ITMAX
    write(*,*) 'EPS    ',EPS
    write(*,*)
  end subroutine dump_fstr_global

  !> This subroutine prints out solution control parameters
  subroutine dump_fstr_param( p )
    implicit none
    type( fstr_param ) :: p

    write(*,*) 'fstrPARAM dump ********************'
    write(*,*)
    write(*,*) 'solution_type ',p%solution_type
    write(*,*) 'solver_method ',p%solver_method
    write(*,*)
    write(*,*) '!!STATIC !HEAT'
    write(*,*) p%analysis_n
    if( associated( P%dtime))  write(*,*) 'dtime ', p%dtime
    if( associated( P%etime))  write(*,*) 'etime', p%etime
    if( associated( P%dtmin))  write(*,*) 'dtmin ', p%dtmin
    if( associated( P%delmax))   write(*,*) 'delmax ', p%delmax
    if( associated( P%itmax))  write(*,*) 'itmax ', p%itmax
    if( associated( P%eps))    write(*,*) 'eps ', p%eps
    write(*,*) 'ref_temp ', p%ref_temp
    write(*,*)
    write(*,*) 'output control'
    write(*,*) 'fg_echo ', p%fg_echo
    write(*,*) 'fg_result ', p%fg_result
    write(*,*) 'fg_visual ', p%fg_visual
    write(*,*)
    write(*,*) 'for heat ...'
    write(*,*) 'fg_neutral ', p%fg_neutral
    write(*,*) 'fg_irres ', p%fg_irres
    write(*,*) 'fg_iwres ', p%fg_iwres
    write(*,*) 'nrres ', p%nrres
    write(*,*) 'nprint ', p%nprint
    write(*,*)
    write(*,*) 'for couple ...'
    write(*,*) 'fg_couple ',p%fg_couple
    write(*,*)
    write(*,*) 'ndex table for global node ID sorting'
    write(*,*) 'n_node ', p%n_node
    if( associated( P%global_local_ID)) write(*,*) 'global_local_ID ', p%global_local_ID
  end subroutine dump_fstr_param

  !> This subroutine prints out data for static analysis
  subroutine dump_fstr_solid( s )
    implicit none
    type( fstr_solid ) :: s

    write(*,*) 'fstrSOLID dump ********************'
    write(*,*)
    write(*,*) 'file_type ', s%file_type
    write(*,*)
    write(*,*) '!!BOUNDARY'
    write(*,*) 'BOUNDARY_ngrp_tot ', s%BOUNDARY_ngrp_tot
    if( s%BOUNDARY_ngrp_tot /= 0 ) then
      write(*,*) 'BOUNDARY_ngrp_ID ', s%BOUNDARY_ngrp_ID
      write(*,*) 'BOUNDARY_ngrp_type ',s%BOUNDARY_ngrp_type
      write(*,*) 'BOUNDARY_ngrp_val ',s%BOUNDARY_ngrp_val
      write(*,*) 'BOUNDARY_ngrp_amp ',s%BOUNDARY_ngrp_amp
    end if
    write(*,*)
    write(*,*) '!!VELOCITY'
    write(*,*) 'VELOCITY_ngrp_tot ', s%VELOCITY_ngrp_tot
    if( s%VELOCITY_ngrp_tot /= 0 ) then
      write(*,*) 'VELOCITY_ngrp_ID ', s%VELOCITY_ngrp_ID
      write(*,*) 'VELOCITY_ngrp_type ',s%VELOCITY_ngrp_type
      write(*,*) 'VELOCITY_ngrp_val ',s%VELOCITY_ngrp_val
      write(*,*) 'VELOCITY_ngrp_amp ',s%VELOCITY_ngrp_amp
    end if
    write(*,*)
    write(*,*) '!!ACCELERATION'
    write(*,*) 'ACCELERATION_ngrp_tot ', s%ACCELERATION_ngrp_tot
    if( s%ACCELERATION_ngrp_tot /= 0 ) then
      write(*,*) 'ACCELERATION_ngrp_ID ', s%ACCELERATION_ngrp_ID
      write(*,*) 'ACCELERATION_ngrp_type ',s%ACCELERATION_ngrp_type
      write(*,*) 'ACCELERATION_ngrp_val ',s%ACCELERATION_ngrp_val
      write(*,*) 'ACCELERATION_ngrp_amp ',s%ACCELERATION_ngrp_amp
    end if
    write(*,*)
    write(*,*) '!!CLOAD'
    write(*,*) 'CLOAD_ngrp_tot ', s%CLOAD_ngrp_tot
    if( s%CLOAD_ngrp_tot /= 0 ) then
      write(*,*) 'CLOAD_ngrp_ID ', s%CLOAD_ngrp_ID
      write(*,*) 'CLOAD_ngrp_DOF ', s%CLOAD_ngrp_DOF
      write(*,*) 'CLOAD_ngrp_val ',s%CLOAD_ngrp_val
      write(*,*) 'CLOAD_ngrp_amp ',s%CLOAD_ngrp_amp
    end if
    write(*,*)
    write(*,*) '!!DLOAD'
    write(*,*) 'DLOAD_ngrp_tot ', s%DLOAD_ngrp_tot
    if( s%DLOAD_ngrp_tot/= 0 ) then
      write(*,*) 'DLOAD_ngrp_ID ',s%DLOAD_ngrp_ID
      write(*,*) 'DLOAD_ngrp_LID ',s%DLOAD_ngrp_LID
      write(*,*) 'DLOAD_ngrp_params ', s%DLOAD_ngrp_params
      write(*,*) 'DLOAD_ngrp_amp ',s%DLOAD_ngrp_amp
    end if
    write(*,*)
    write(*,*) '!!TEMPERATURE'
    write(*,*) 'TEMP_ngrp_tot ',s%TEMP_ngrp_tot
    if( s%TEMP_ngrp_tot/= 0 ) then
      write(*,*) 'TEMP_ngrp_ID ',s%TEMP_ngrp_ID
      write(*,*) 'TEMP_ngrp_val ', s%TEMP_ngrp_val
    end if
    write(*,*)
    write(*,*) '!!STATIC'
    write(*,*) 'restart_nout ',s%restart_nout
    write(*,*)
    write(*,*) '!!COUPLE'
    write(*,*) 'COUPLE_ngrp_tot ',s%COUPLE_ngrp_tot
    if( s%COUPLE_ngrp_tot>0 ) then
      write(*,*) 'COUPLE_ngrp_ID ', s%COUPLE_ngrp_ID
    endif
    write(*,*)
  end subroutine dump_fstr_solid

  !> This subroutine prints out data for heat conductive analysis
  subroutine dump_fstr_heat( h )
    implicit none
    type( fstr_heat ) :: h

    write(*,*) 'fstrHEAT dump ********************'
    write(*,*)
    write(*,*) 'TIME CONTROL'
    write(*,*) 'STEPtot ', h%STEPtot
    if( h%STEPtot /= 0 ) then
      write(*,*) 'STEP_DLTIME ', h%STEP_DLTIME
      write(*,*) 'STEP_EETIME ', h%STEP_EETIME
      write(*,*) 'STEP_DELMIN ', h%STEP_DELMIN
      write(*,*) 'STEP_DELMAX ', h%STEP_DELMAX
    end if
    write(*,*)
    write(*,*) 'MATERIAL'
    write(*,*) 'ATERIALtot ', h%MATERIALtot
    if( h%MATERIALtot /= 0 ) then
      write(*,*) 'RHO ', h%RHO
      write(*,*) 'RHOtemp ', h%RHOtemp
      write(*,*) 'CP ',h%CP
      write(*,*) 'CPtemp ', h%CPtemp
      write(*,*) 'COND ', h%COND
      write(*,*) 'CONDtemp ',h%CONDtemp
      write(*,*)
      write(*,*) 'RHOtab ', h%RHOtab
      write(*,*) 'CPtab ', h%CPtab
      write(*,*) 'CONDtab ',h%CONDtab
      write(*,*)
      write(*,*) 'RHOfuncA ', h%RHOfuncA
      write(*,*) 'RHOfuncB ', h%RHOfuncB
      write(*,*) 'CPfuncA ',h%CPfuncA
      write(*,*) 'CPfuncB ',h%CPfuncB
      write(*,*) 'CONDfuncA ',h%CONDfuncA
      write(*,*) 'CONDfuncB ',h%CONDfuncB
    end if
    write(*,*)
    write(*,*) 'AMPLITUDE'
    write(*,*) 'AMPLITUDEtot ',h%AMPLITUDEtot
    if( h%AMPLITUDEtot /=0 ) then
      write(*,*) 'AMPL ',h%AMPL
      write(*,*) 'AMPLtime ',h%AMPLtime
      write(*,*) 'AMPLtab ', h%AMPLtab
      write(*,*) 'AMPLfuncA ', h%AMPLfuncA
      write(*,*) 'AMPLfuncB ', h%AMPLfuncB
    end if
    write(*,*)
    write(*,*) 'VALUE'
    if( associated(h%TEMP)) then
      write(*,*) 'TEMP ',h%TEMP
      write(*,*) 'TEMP0 ',h%TEMP0
      write(*,*) 'TEMPC ',h%TEMPC
    end if
    write(*,*)
    write(*,*) 'BOUNDARY CONDTIONS -------'
    write(*,*)
    write(*,*) '!FIXTEMP '
    write(*,*) 'T_FIX_tot ',h%T_FIX_tot
    if( h%T_FIX_tot /= 0 ) then
      write(*,*) 'T_FIX_node ',h%T_FIX_node
      write(*,*) 'T_FIX_ampl ',h%T_FIX_ampl
      write(*,*) 'T_FIX_val ',h%T_FIX_val
    end if
    write(*,*)
    write(*,*) '!CFLUX'
    write(*,*) 'Q_NOD_tot ',h%Q_NOD_tot
    if( h%Q_NOD_tot /= 0 ) then
      write(*,*) 'Q_NOD_node ', h%Q_NOD_node
      write(*,*) 'Q_NOD_ampl ',h%Q_NOD_ampl
      write(*,*) 'Q_NOD_val ', h%Q_NOD_val
    end if
    write(*,*)
    write(*,*) '!DFLUX (not used)'
    write(*,*) 'Q_VOL_tot ',h%Q_VOL_tot
    if( h%Q_VOL_tot /= 0 ) then
      write(*,*) 'Q_VOL_elem ',h%Q_VOL_elem
      write(*,*) 'Q_VOL_ampl ',h%Q_VOL_ampl
      write(*,*) 'Q_VOL_val ',h%Q_VOL_val
    end if
    write(*,*)
    write(*,*) '!DFLUX, !SFLUX'
    write(*,*) 'Q_SUF_tot ', h%Q_SUF_tot
    if( h%Q_SUF_tot /= 0 ) then
      write(*,*) 'Q_SUF_elem ', h%Q_SUF_elem
      write(*,*) 'Q_SUF_ampl ',h%Q_SUF_ampl
      write(*,*) 'Q_SUF_surf ',h%Q_SUF_surf
      write(*,*) 'Q_SUF_val ',h%Q_SUF_val
    end if
    write(*,*)
    write(*,*) '!RADIATE, !SRADIATE'
    write(*,*) 'R_SUF_tot ',h%R_SUF_tot
    if( h%R_SUF_tot /= 0 ) then
      write(*,*) 'R_SUF_elem ',h%R_SUF_elem
      write(*,*) 'R_SUF_ampl ', h%R_SUF_ampl
      write(*,*) 'R_SUF_surf ',h%R_SUF_surf
      write(*,*) 'R_SUF_val ',h%R_SUF_val
    end if
    write(*,*)
    write(*,*) '!FILM, SFILM'
    write(*,*) 'H_SUF_tot ',h%H_SUF_tot
    if( h%H_SUF_tot /= 0 ) then
      write(*,*) 'H_SUF_elem ',h%H_SUF_elem
      write(*,*) 'H_SUF_ampl ',h%H_SUF_ampl
      write(*,*) 'H_SUF_surf ',h%H_SUF_surf
      write(*,*) 'H_SUF_val ',h%H_SUF_val
    end if
    write(*,*)
  end subroutine dump_fstr_heat

  !> This subroutine prints out parameters for eigen analysis
  subroutine dump_fstr_eigen( e )
    implicit none
    type( fstr_eigen ) :: e

    write(*,*) 'lczparam dump ********************'
    write(*,*)
    write(*,*) 'nget    ', e%nget
    write(*,*)
  end subroutine dump_fstr_eigen

  !> This subroutine prints out data for dynamic analysis
  subroutine dump_fstr_dynamic( d )
    implicit none
    type( fstr_dynamic ) :: d

    write(*,*) 'fstrDYNAMIC dump ********************'
    write(*,*)
    write(*,*) 'idx_eqa  ', d%idx_eqa
    write(*,*) 'idx_resp ', d%idx_resp
    write(*,*) 'n_step   ', d%n_step
    write(*,*) 't_start  ', d%t_start
    write(*,*) 't_end    ', d%t_end
    write(*,*) 't_delta  ', d%t_delta
    write(*,*) 'ganma    ', d%ganma
    write(*,*) 'beta     ', d%beta
    write(*,*) 'idx_mas  ', d%idx_mas
    write(*,*) 'idx_dmp  ', d%idx_dmp
    write(*,*) 'ray_m    ', d%ray_m
    write(*,*) 'ray_k    ', d%ray_k
    write(*,*) 'restart_nout',d%restart_nout
    write(*,*) 'nout     ', d%nout
    write(*,*) 'ngrp_monit  ', d%ngrp_monit
    write(*,*) 'nout_monit  ', d%nout_monit
    write(*,*) 'iout_list   ', d%iout_list
    write(*,*)
  end subroutine dump_fstr_dynamic

  !> This subroutine prints out coupleing analysis
  subroutine dump_fstr_couple( c )
    implicit none
    type( fstr_couple ) :: c
    integer( kind=kint) :: i,j
    write(*,*) 'fstrCLP dump ********************'
    write(*,*)
    write(*,*) 'dof ', c%dof
    write(*,*) 'ndof ', c%ndof
    write(*,*) 'coupled_node_n ', c%coupled_node_n
    if( c%coupled_node_n >0 ) then
      write(*,*) 'coupled_node'
      write(*,*) c%coupled_node
    endif
    write(*,*) 'trac'
    do i=1, c%coupled_node_n
      j = c%dof * i;
      write(*,*)  c%trac(j-2),' ',c%trac(j-1),' ',c%trac(j)
    end do
    write(*,*) 'velo'
    do i=1, c%coupled_node_n
      j = c%dof * i;
      write(*,*)  c%velo(j-2),' ',c%velo(j-1),' ',c%velo(j)
    end do
  end subroutine dump_fstr_couple


end module fstr_debug_dump


