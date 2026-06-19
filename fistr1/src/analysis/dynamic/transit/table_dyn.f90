!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Amplitude evaluation for loading conditions in dynamic analysis
module m_table_dyn

  use m_fstr
  use hecmw

  implicit none

contains

  !C================================================================C
  !C-- subroutine fstr_get_amplitude_dyn
  !C================================================================C
  !> \brief Evaluate the amplitude-scaled target value a(t) for dynamic analysis.
  !!
  !! On entry `value` holds the nominal card value; on exit it holds the full
  !! target value a(t). When no per-card amplitude is attached the value is
  !! returned unchanged (equivalent to a factor of 1.0). The elapsed-time basis
  !! of the dynamic solver is preserved (the explicit scheme evaluates one step
  !! behind, at t_curr - t_delta); the actual interpolation, end-point clamping
  !! and RELATIVE/ABSOLUTE handling are delegated to the canonical evaluator
  !! hecmw_get_amplitude_value, shared with static analysis.
  subroutine fstr_get_amplitude_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, t_curr, value, flag_u)
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC
    integer(kind=kint)       :: ig0
    real(kind=kreal)         :: t_curr
    real(kind=kreal)         :: value    !< in: nominal value, out: a(t)
    integer(kind=kint)       :: flag_u

    integer(kind=kint) :: jj_n_amp
    real(kind=kreal)   :: t_eval

    jj_n_amp = 0
    if( flag_u == 1 ) then
      jj_n_amp = fstrSOLID%BOUNDARY_ngrp_amp(ig0)
    else if( flag_u == 2 ) then
      jj_n_amp = fstrSOLID%VELOCITY_ngrp_amp(ig0)
    else if( flag_u == 3 ) then
      jj_n_amp = fstrSOLID%ACCELERATION_ngrp_amp(ig0)
    else if( flag_u == 0 ) then
      jj_n_amp = fstrSOLID%CLOAD_ngrp_amp(ig0)
    else if( flag_u == 10 ) then
      jj_n_amp = fstrSOLID%DLOAD_ngrp_amp(ig0)
    end if

    if( jj_n_amp <= 0 ) return  ! no amplitude: value unchanged (factor 1.0)

    t_eval = t_curr
    if( fstrDYNAMIC%idx_eqa == 11 ) t_eval = t_curr - fstrDYNAMIC%t_delta

    call hecmw_get_amplitude_value(hecMESH%amp, jj_n_amp, t_eval, value)

  end subroutine fstr_get_amplitude_dyn

end module m_table_dyn
