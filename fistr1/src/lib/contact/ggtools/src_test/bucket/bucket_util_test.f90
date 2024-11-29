!> バケットセルテストモジュール
module mod_ggtools_bucket_util_test
  use mod_monolis_utils
  use mod_ggtools
  implicit none

contains

  subroutine ggtools_bucket_util_test()
    implicit none
    integer(kint) :: index(3), nx(3), id
    real(kdouble) :: xmin(3), dx(3), pos(3)

    !ggtools_bucket_get_integer_coordinate_from_real_coordinate
    call monolis_std_global_log_string("ggtools_bucket_get_integer_coordinate_from_real_coordinate")

    xmin(1) = 1.0d0
    xmin(2) = 2.0d0
    xmin(3) = 3.0d0

    dx(1) = 0.1d0
    dx(2) = 0.2d0
    dx(3) = 0.3d0

    pos(1) = 1.0d0
    pos(2) = 2.0d0
    pos(3) = 3.0d0

    index = ggtools_bucket_get_integer_coordinate_from_real_coordinate(xmin, dx, pos)

    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 1-1", index(1), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 1-2", index(2), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 1-3", index(3), 1)

    pos(1) = 2.0d0
    pos(2) = 4.0d0
    pos(3) = 6.0d0

    index = ggtools_bucket_get_integer_coordinate_from_real_coordinate(xmin, dx, pos)

    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 2-1", index(1), 11)
    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 2-2", index(2), 11)
    call monolis_test_check_eq_I1("integer_coordinate_from_real_coordinate 2-3", index(3), 11)

    !ggtools_bucket_get_integer_coordinate_from_1d_array_index
    call monolis_std_global_log_string("ggtools_bucket_get_integer_coordinate_from_1d_array_index")

    nx(1) = 2
    nx(2) = 3
    nx(3) = 4

    index = ggtools_bucket_get_integer_coordinate_from_1d_array_index(1, nx)

    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 1-1", index(1), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 1-2", index(2), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 1-3", index(3), 1)

    index = ggtools_bucket_get_integer_coordinate_from_1d_array_index(6, nx)

    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 2-1", index(1), 2)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 2-2", index(2), 3)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 2-3", index(3), 1)

    index = ggtools_bucket_get_integer_coordinate_from_1d_array_index(7, nx)

    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 3-1", index(1), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 3-2", index(2), 1)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 3-3", index(3), 2)

    index = ggtools_bucket_get_integer_coordinate_from_1d_array_index(24, nx)

    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 4-1", index(1), 2)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 4-2", index(2), 3)
    call monolis_test_check_eq_I1("integer_coordinate_from_1d_array_index 4-3", index(3), 4)

    !ggtools_bucket_get_1d_array_index_from_integer_coodinate
    call monolis_std_global_log_string("ggtools_bucket_get_1d_array_index_from_integer_coodinate")

    nx(1) = 2
    nx(2) = 3
    nx(3) = 4

    index(1) = 1
    index(2) = 1
    index(3) = 1

    id =  ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, index)

    call monolis_test_check_eq_I1("ggtools_bucket_get_1d_array_index_from_integer_coodinate 1", id, 1)

    index(1) = 2
    index(2) = 1
    index(3) = 1

    id =  ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, index)

    call monolis_test_check_eq_I1("ggtools_bucket_get_1d_array_index_from_integer_coodinate 2", id, 2)

    index(1) = 1
    index(2) = 2
    index(3) = 1

    id =  ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, index)

    call monolis_test_check_eq_I1("ggtools_bucket_get_1d_array_index_from_integer_coodinate 3", id, 3)

    index(1) = 1
    index(2) = 3
    index(3) = 1

    id =  ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, index)

    call monolis_test_check_eq_I1("ggtools_bucket_get_1d_array_index_from_integer_coodinate 4", id, 5)

    index(1) = 1
    index(2) = 1
    index(3) = 2

    id =  ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, index)

    call monolis_test_check_eq_I1("ggtools_bucket_get_1d_array_index_from_integer_coodinate 5", id, 7)
  end subroutine ggtools_bucket_util_test

end module mod_ggtools_bucket_util_test
