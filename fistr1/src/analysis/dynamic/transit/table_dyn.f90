!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Table of lading step in dynamic analysis
module m_table_dyn

  use m_fstr
  use hecmw

  implicit none

contains

  !C================================================================C
  !C-- subroutine table_dyn
  !C================================================================C
  subroutine table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC

    integer(kind=kint) :: i, ig0
    integer(kind=kint) :: jj_n_amp, jj1, jj2
    integer(kind=kint) :: s1, s2, flag_u
    real(kind=kreal)   :: t_1, t_2, t_t, f_1, f_2, f_t

    jj_n_amp = 0
    s1 = 0; s2 = 0
    t_1 = 0.0d0; t_2 = 0.0d0; t_t = 0.0d0; f_1 = 0.0d0; f_2 = 0.0d0; f_t = 0.0d0

    if( flag_u .eq. 1 ) then
      jj_n_amp = fstrSOLID%BOUNDARY_ngrp_amp(ig0)
    else if( flag_u .eq. 2 ) then
      jj_n_amp = fstrSOLID%VELOCITY_ngrp_amp(ig0)
    else if( flag_u .eq. 3 ) then
      jj_n_amp = fstrSOLID%ACCELERATION_ngrp_amp(ig0)
    else if( flag_u .eq. 0 ) then
      jj_n_amp = fstrSOLID%CLOAD_ngrp_amp(ig0)
    else if( flag_u .eq. 10 ) then
      jj_n_amp = fstrSOLID%DLOAD_ngrp_amp(ig0)
    end if

    if( jj_n_amp == 0 ) then
      f_t = 1.d0
    else

      jj1 = hecMESH%amp%amp_index(jj_n_amp - 1)
      jj2 = hecMESH%amp%amp_index(jj_n_amp)

      jj1 = jj1 + 2
      if( fstrDYNAMIC%idx_eqa == 1 ) then
        t_t = fstrDYNAMIC%t_curr

      else if( fstrDYNAMIC%idx_eqa == 11 ) then
        select case (flag_u)
          case (0)
            t_t = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
          case (10)
            t_t = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
          case (1)
            t_t = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
          case (2)
            t_t = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
          case (3)
            t_t = fstrDYNAMIC%t_curr - fstrDYNAMIC%t_delta
        end select
      end if

      if( fstrDYNAMIC%i_step == 0 ) then
        t_t = fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step
      end if

      !     if(jj2 .eq. 0) then
      !        f_t = 1.0

      if(t_t .gt. hecMESH%amp%amp_table(jj2)) then
        f_t = hecMESH%amp%amp_val(jj2)
      else if(t_t .le. hecMESH%amp%amp_table(jj2)) then
        do i = jj1, jj2
          if(t_t .le. hecMESH%amp%amp_table(i)) then
            s2 = i
            s1 = i - 1
            exit
          end if
        end do

        t_2 = hecMESH%amp%amp_table(s2)
        t_1 = hecMESH%amp%amp_table(s1)
        f_2 = hecMESH%amp%amp_val(s2)
        f_1 = hecMESH%amp%amp_val(s1)
        if( t_2-t_1 .lt. 1.0e-20) then
          if( hecMESH%my_rank.eq.0) then
            write(imsg,*) 'stop due to t_2-t_1 <= 0'
          end if
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        f_t = ((t_2*f_1 - t_1*f_2) + (f_2 - f_1)*t_t) / (t_2 - t_1)
      end if

    end if

  end subroutine table_dyn

end module m_table_dyn
