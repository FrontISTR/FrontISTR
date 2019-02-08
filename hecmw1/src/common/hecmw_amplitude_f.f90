!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!> \brief Handle amplitude

module hecmw_amplitude
  use hecmw_util
  implicit none

contains

  subroutine hecmw_get_amplitude_value(amp, amp_id, time, value)
    type(hecmwST_amplitude), intent(in) :: amp
    integer(kind=kint), intent(in)      :: amp_id
    real(kind=kreal), intent(in)        :: time
    real(kind=kreal), intent(inout)     :: value

    integer(kind=kint) :: i, jj1, jj2, s1, s2
    real(kind=kreal)   :: t_1, t_2
    real(kind=kreal)   :: f_t, f_1, f_2

    if( amp_id <= 0 .or. amp_id > amp%n_amp ) return

    jj1 = amp%amp_index(amp_id-1)+1
    jj2 = amp%amp_index(amp_id)

    if(time >= amp%amp_table(jj2)) then
      f_t = amp%amp_val(jj2)
    else if(time <= amp%amp_table(jj1)) then
      f_t = amp%amp_val(jj1)
    else
      do i = jj1+1, jj2
        if(time <= amp%amp_table(i)) then
          s2 = i
          s1 = s2 - 1
          exit
        endif
      end do

      t_2 = amp%amp_table(s2)
      t_1 = amp%amp_table(s1)
      f_2 = amp%amp_val(s2)
      f_1 = amp%amp_val(s1)
      if( t_2-t_1 <= 1.0d-20) then
        f_t = f_2
      else
        f_t = ((t_2*f_1 - t_1*f_2) + (f_2 - f_1)*time) / (t_2 - t_1)
      endif
    endif

    if( amp%amp_type_value(amp_id) == HECMW_AMP_TYPEVAL_RELATIVE ) then
      value = f_t*value
    else
      value = f_t
    end if

  end subroutine hecmw_get_amplitude_value

end module hecmw_amplitude
