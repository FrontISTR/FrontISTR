module hecmw_api_common
  use iso_c_binding
  implicit none

contains

  ! hecmw_lib_fc.c の HECMW_strcpy_f2c と同等だがメモリ確保していることを前提
  ! c_str は c_len の大きさで呼び出し側で確保済み
  subroutine f_c_str_copy(f_str, c_str, c_len)
    implicit none
    character(len=*), intent(in)  :: f_str
    character(kind=c_char), intent(inout) :: c_str(*)
    integer(c_int), value, intent(in) :: c_len
    integer :: f_len, i

    f_len = len_trim(f_str)
    if (f_len>=c_len) then
      f_len = c_len-1
    end if

    do i = 1, f_len
       c_str(i) = f_str(i:i)
    end do
    c_str(f_len+1) = c_null_char
  end subroutine

  ! hecmw_lib_fc.c の HECMW_strcpy_c2f と同等
  ! f_str は f_len の大きさで呼び出し側で確保済み
  subroutine c_f_str_copy(c_str, f_str, f_len)
    implicit none
    character(kind=c_char), intent(in) :: c_str(*)
    character(len=*), intent(inout) :: f_str
    integer(c_int), value, intent(in) :: f_len

    integer :: c_len, i

    c_len = 0
    do while (c_str(c_len + 1) /= c_null_char)
      c_len = c_len + 1
    end do

    if (c_len > f_len) then
      c_len = f_len
    end if

    do i = 1, c_len
       f_str(i:i) = c_str(i)
    end do

    if (c_len < f_len) then
       f_str(c_len + 1 : f_len) = ' '
    end if
  end subroutine c_f_str_copy

end module
