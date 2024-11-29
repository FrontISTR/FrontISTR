!> ローカルバケットセルモジュール
module mod_ggtools_bucket_util
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup dev
  !> 実数座標を入れるとバケットセル整数座標を取得
  function ggtools_bucket_get_integer_coordinate_from_real_coordinate(xmin, dx, pos)
    implicit none
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バケットセルひとつのサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    !> 入力座標
    real(kdouble) :: pos(3)
    !> バケットセル整数座標
    integer(kint) :: ggtools_bucket_get_integer_coordinate_from_real_coordinate(3)
    integer(kint) :: id(3)

    ! if pos(i) = xmin(i), set id = 1(not 0)
    ! else ceiling((pos(1) - xmin(1))/dx(1))
    id(1) = max(ceiling((pos(1) - xmin(1))/dx(1)),1)
    id(2) = max(ceiling((pos(2) - xmin(2))/dx(2)),1)
    id(3) = max(ceiling((pos(3) - xmin(3))/dx(3)),1)

    !if(id(1) < 1)
    !if(id(2) < 1)
    !if(id(3) < 1)
    !if(id(1) > nx(1))
    !if(id(2) > nx(2))
    !if(id(3) > nx(3))

    ggtools_bucket_get_integer_coordinate_from_real_coordinate = id
  end function ggtools_bucket_get_integer_coordinate_from_real_coordinate

  !> @ingroup dev
  !> バケットセルの 1 次元配列通し番号を入れるとバケットセル整数座標を取得
  function ggtools_bucket_get_integer_coordinate_from_1d_array_index(array_index, nx)
    implicit none
    !> バケットセルの 1 次元配列通し番号
    integer(kint), intent(in) :: array_index
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    !> バケットセル整数座標
    integer(kint) :: ggtools_bucket_get_integer_coordinate_from_1d_array_index(3)
    integer(kint) :: id(3), in

    !> 割り切れない場合は切り捨て（言語仕様）
    in = array_index - 1
    id(3) = in/(nx(1)*nx(2))
    in = mod(in, nx(1)*nx(2))
    id(2) = in/nx(1)
    id(1) = mod(in, nx(1))

    ggtools_bucket_get_integer_coordinate_from_1d_array_index = id + 1
  end function ggtools_bucket_get_integer_coordinate_from_1d_array_index

  !> @ingroup dev
  !> バケットセル整数座標からバケットセルの 1 次元配列通し番号を取得
  function ggtools_bucket_get_1d_array_index_from_integer_coodinate(nx, id)
    implicit none
    !> バケットセル分割数
    integer(kint) :: nx(3)
    !> バケットセル整数座標
    integer(kint) :: id(3)
    !> バケットセル id
    integer(kint) :: ggtools_bucket_get_1d_array_index_from_integer_coodinate
    integer(kint) :: index

    !if(id(1) < 1)
    !if(id(2) < 1)
    !if(id(3) < 1)
    !if(id(1) > nx(1))
    !if(id(2) > nx(2))
    !if(id(3) > nx(3))

    index = id(1) + (id(2)-1)*nx(1) + (id(3)-1)*nx(1)*nx(2)
    ggtools_bucket_get_1d_array_index_from_integer_coodinate = index
  end function ggtools_bucket_get_1d_array_index_from_integer_coodinate

end module mod_ggtools_bucket_util
