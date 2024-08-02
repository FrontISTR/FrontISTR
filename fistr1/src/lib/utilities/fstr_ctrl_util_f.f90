module fstr_ctrl_util_f
  use hecmw
  implicit none

  ! external functions implemented in fstr_ctr_util.c

  interface

    function fstr_ctrl_rewind( ctrl ) &
        bind(c,name='fstr_ctrl_rewind')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_rewind
      integer(c_int) :: ctrl
    end function fstr_ctrl_rewind

    subroutine fstr_ctrl_get_err_msg_c( buf, len ) &
        bind(c,name='fstr_ctrl_get_err_msg')
      use iso_c_binding
      type(c_ptr),value :: buf
      integer(c_int) :: len
    end subroutine fstr_ctrl_get_err_msg_c

    function fstr_ctrl_open_c( filename ) &
        bind(c,name='fstr_ctrl_open')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_open_c
      type(c_ptr),value :: filename
    end function fstr_ctrl_open_c

    function fstr_ctrl_get_rec_number( ctrl ) &
        bind(c,name='fstr_ctrl_get_rec_number')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_rec_number
      integer(c_int) :: ctrl
    end function fstr_ctrl_get_rec_number

    function fstr_ctrl_get_line_c( ctrl, rec_no, buff, buff_size ) &
        bind(c,name='fstr_ctrl_get_line')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_line_c
      integer(c_int) :: ctrl
      integer(c_int) :: rec_no
      type(c_ptr),value :: buff
      integer(c_int) :: buff_size
    end function fstr_ctrl_get_line_c

    function fstr_ctrl_seek_header_c( ctrl, header_name ) &
        bind(c,name='fstr_ctrl_seek_header')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_seek_header_c
      integer(c_int) :: ctrl
      type(c_ptr),value :: header_name
    end function fstr_ctrl_seek_header_c

    function fstr_ctrl_seek_next_header( ctrl ) &
        bind(c,name='fstr_ctrl_seek_next_header')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_seek_next_header
      integer(c_int) :: ctrl
    end function fstr_ctrl_seek_next_header

    function fstr_ctrl_get_c_h_name_c( ctrl, header_name, buf_size ) &
        bind(c,name='fstr_ctrl_get_c_h_name')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_c_h_name_c
      integer(c_int) :: ctrl
      type(c_ptr),value :: header_name
      integer(c_int) :: buf_size
    end function fstr_ctrl_get_c_h_name_c

    !!! not used
    function fstr_ctrl_get_c_h_line_no( ctrl ) &
        bind(c,name='fstr_ctrl_get_c_h_line_no')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_c_h_line_no
      integer(c_int) :: ctrl
    end function fstr_ctrl_get_c_h_line_no

    function fstr_ctrl_get_c_h_pos( ctrl ) &
        bind(c,name='fstr_ctrl_get_c_h_pos')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_c_h_pos
      integer(c_int) :: ctrl
    end function fstr_ctrl_get_c_h_pos

    !!! not used
    function fstr_ctrl_get_param_c( ctrl, param_name, value_list, type, val ) &
        bind(c,name='fstr_ctrl_get_param')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_param_c
      integer(c_int) :: ctrl
      type(c_ptr),value :: param_name, value_list
      type(c_ptr),value :: type
      type(c_ptr),value :: val
    end function fstr_ctrl_get_param_c

    function fstr_ctrl_get_param_ex_c( ctrl, param_name, value_list, necessity, type, val ) &
        bind(c,name='fstr_ctrl_get_param_ex')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_param_ex_c
      integer(c_int) :: ctrl
      type(c_ptr),value :: param_name, value_list
      integer(c_int) :: necessity
      type(c_ptr),value :: type
      type(c_ptr),value :: val
    end function fstr_ctrl_get_param_ex_c

    function fstr_ctrl_get_data_line_n( ctrl ) &
        bind(c,name='fstr_ctrl_get_data_line_n')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_line_n
      integer(c_int) :: ctrl
    end function fstr_ctrl_get_data_line_n

    !!! not used
    function fstr_ctrl_get_data_n_in_line_c( ctrl, line_no, delim ) &
        bind(c,name='fstr_ctrl_get_data_n_in_line')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_n_in_line_c
      integer(c_int) :: ctrl, line_no
      type(c_ptr),value :: delim
    end function fstr_ctrl_get_data_n_in_line_c

    !!! not used
    function fstr_ctrl_get_data_error_pos() &
        bind(c,name='fstr_ctrl_get_data_error_pos')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_error_pos
    end function fstr_ctrl_get_data_error_pos

    function fstr_ctrl_get_data_error_line() &
        bind(c,name='fstr_ctrl_get_data_error_line')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_error_line
    end function fstr_ctrl_get_data_error_line

    !!! not used
    function fstr_ctrl_get_data_c( ctrl, line_no, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) &
        bind(c,name='fstr_ctrl_get_data')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_c
      integer(c_int) :: ctrl
      integer(c_int) :: line_no
      type(c_ptr),value :: format
      type(c_ptr),value :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    end function fstr_ctrl_get_data_c

    function fstr_ctrl_get_data_ex_c( ctrl, line_no, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) &
        bind(c,name='fstr_ctrl_get_data_ex')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_ex_c
      integer(c_int) :: ctrl
      integer(c_int) :: line_no
      type(c_ptr),value :: format
      type(c_ptr),value :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    end function fstr_ctrl_get_data_ex_c

    function fstr_ctrl_get_data_array_ex_c( ctrl, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) &
        bind(c,name='fstr_ctrl_get_data_array_ex')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_get_data_array_ex_c
      integer(c_int) :: ctrl
      type(c_ptr),value :: format
      type(c_ptr),value :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    end function fstr_ctrl_get_data_array_ex_c

    function fstr_ctrl_close( ctrl ) &
        bind(c,name='fstr_ctrl_close')
      use iso_c_binding
      integer(c_int) :: fstr_ctrl_close
      integer(c_int) :: ctrl
    end function fstr_ctrl_close

    subroutine fstr_ctrl_dump( ctrl ) &
        bind(c,name='fstr_ctrl_dump')
      use iso_c_binding
      integer(c_int) :: ctrl
    end subroutine fstr_ctrl_dump

  end interface

contains

  subroutine fstr_ctrl_get_err_msg( buf, len )
    use iso_c_binding
    character(len=*),target :: buf
    integer(c_int) :: len
    call fstr_ctrl_get_err_msg_c( c_loc(buf), len )
  end subroutine fstr_ctrl_get_err_msg

  function fstr_ctrl_open( filename )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_open
    character(len=*),target :: filename
    fstr_ctrl_open = fstr_ctrl_open_c( c_loc(filename) )
  end function fstr_ctrl_open

  function fstr_ctrl_get_line( ctrl, rec_no, buff, buff_size )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_line
    integer(c_int) :: ctrl
    integer(c_int) :: rec_no
    character(len=*),target :: buff
    integer(c_int) :: buff_size
    fstr_ctrl_get_line = fstr_ctrl_get_line_c( ctrl, rec_no, c_loc(buff), buff_size )
  end function fstr_ctrl_get_line

  function fstr_ctrl_seek_header( ctrl, header_name )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_seek_header
    integer(c_int) :: ctrl
    character(len=*),target :: header_name
    fstr_ctrl_seek_header = fstr_ctrl_seek_header_c( ctrl, c_loc(header_name) )
  end function fstr_ctrl_seek_header

  function fstr_ctrl_get_c_h_name( ctrl, header_name, buf_size )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_c_h_name
    integer(c_int) :: ctrl
    character(len=*),target :: header_name
    integer(c_int) :: buf_size
    fstr_ctrl_get_c_h_name = fstr_ctrl_get_c_h_name_c( ctrl, c_loc(header_name), buf_size )
  end function fstr_ctrl_get_c_h_name

  function fstr_ctrl_get_param( ctrl, param_name, value_list, type, val )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_param
    integer(c_int) :: ctrl
    character(len=*),target :: param_name, value_list
    character(c_char),target :: type
    type(*),dimension(..),target :: val
    fstr_ctrl_get_param = fstr_ctrl_get_param_c( ctrl, c_loc(param_name), c_loc(value_list), c_loc(type), c_loc(val) )
  end function fstr_ctrl_get_param

  function fstr_ctrl_get_param_ex( ctrl, param_name, value_list, necessity, type, val )
    use iso_c_binding
    use hecmw
    integer(c_int) :: fstr_ctrl_get_param_ex
    integer(c_int) :: ctrl
    character(len=*),target :: param_name, value_list
    integer(c_int) :: necessity
    character(c_char),target :: type
    type(*),dimension(..),target :: val
    fstr_ctrl_get_param_ex = fstr_ctrl_get_param_ex_c( ctrl,c_loc(param_name),c_loc(value_list),necessity,c_loc(type),c_loc(val) )
  end function fstr_ctrl_get_param_ex

  function fstr_ctrl_get_data_n_in_line( ctrl, line_no, delim )
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_data_n_in_line
    integer(c_int) :: ctrl, line_no
    character(len=*),target :: delim
    fstr_ctrl_get_data_n_in_line = fstr_ctrl_get_data_n_in_line_c( ctrl, line_no, c_loc(delim) )
  end function fstr_ctrl_get_data_n_in_line

  function fstr_ctrl_get_data( ctrl, line_no, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_data
    integer(c_int) :: ctrl
    integer(c_int) :: line_no
    character(len=*),target :: format
    type(*),dimension(..),target,optional :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    type(c_ptr) :: pv(10)
    pv(:) = C_NULL_PTR
    if( present(v1) ) pv(1) = c_loc(v1)
    if( present(v2) ) pv(2) = c_loc(v2)
    if( present(v3) ) pv(3) = c_loc(v3)
    if( present(v4) ) pv(4) = c_loc(v4)
    if( present(v5) ) pv(5) = c_loc(v5)
    if( present(v6) ) pv(6) = c_loc(v6)
    if( present(v7) ) pv(7) = c_loc(v7)
    if( present(v8) ) pv(8) = c_loc(v8)
    if( present(v9) ) pv(9) = c_loc(v9)
    if( present(v10) ) pv(10) = c_loc(v10)
    fstr_ctrl_get_data = fstr_ctrl_get_data_c( ctrl, line_no, c_loc(format), &
        pv(1), pv(2), pv(3), pv(4), pv(5), pv(6), pv(7), pv(8), pv(9), pv(10) )
  end function fstr_ctrl_get_data

  function fstr_ctrl_get_data_ex( ctrl, line_no, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_data_ex
    integer(c_int) :: ctrl
    integer(c_int) :: line_no
    character(len=*),target :: format
    type(*),dimension(..),target,optional :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    type(c_ptr) :: pv(10)
    pv(:) = C_NULL_PTR
    if( present(v1) ) pv(1) = c_loc(v1)
    if( present(v2) ) pv(2) = c_loc(v2)
    if( present(v3) ) pv(3) = c_loc(v3)
    if( present(v4) ) pv(4) = c_loc(v4)
    if( present(v5) ) pv(5) = c_loc(v5)
    if( present(v6) ) pv(6) = c_loc(v6)
    if( present(v7) ) pv(7) = c_loc(v7)
    if( present(v8) ) pv(8) = c_loc(v8)
    if( present(v9) ) pv(9) = c_loc(v9)
    if( present(v10) ) pv(10) = c_loc(v10)
    fstr_ctrl_get_data_ex = fstr_ctrl_get_data_ex_c( ctrl, line_no, c_loc(format), &
        pv(1), pv(2), pv(3), pv(4), pv(5), pv(6), pv(7), pv(8), pv(9), pv(10) )
  end function fstr_ctrl_get_data_ex

  function fstr_ctrl_get_data_array_ex( ctrl, format, v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
    use iso_c_binding
    integer(c_int) :: fstr_ctrl_get_data_array_ex
    integer(c_int) :: ctrl
    character(len=*),target :: format
    type(*),dimension(..),target,optional :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10
    type(c_ptr) :: pv(10)
    pv(:) = C_NULL_PTR
    if( present(v1) ) pv(1) = c_loc(v1)
    if( present(v2) ) pv(2) = c_loc(v2)
    if( present(v3) ) pv(3) = c_loc(v3)
    if( present(v4) ) pv(4) = c_loc(v4)
    if( present(v5) ) pv(5) = c_loc(v5)
    if( present(v6) ) pv(6) = c_loc(v6)
    if( present(v7) ) pv(7) = c_loc(v7)
    if( present(v8) ) pv(8) = c_loc(v8)
    if( present(v9) ) pv(9) = c_loc(v9)
    if( present(v10) ) pv(10) = c_loc(v10)
    fstr_ctrl_get_data_array_ex = fstr_ctrl_get_data_array_ex_c( ctrl, c_loc(format), &
        pv(1), pv(2), pv(3), pv(4), pv(5), pv(6), pv(7), pv(8), pv(9), pv(10) )
  end function fstr_ctrl_get_data_array_ex

end module fstr_ctrl_util_f
