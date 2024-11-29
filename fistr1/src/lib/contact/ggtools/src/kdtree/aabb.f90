!> ローカルバウンディングボックスモジュール
module mod_ggtools_lbb
  use mod_monolis_utils
  use mod_ggtools_def_local_bucket

  implicit none

contains

  !> @ingroup aabb
  !> 座標の集合からバウンディングボックスを作成
  subroutine ggtools_get_lbb(coord, xmin, xmax)
    implicit none
    !> [in] 座標の集合（[1,:] に x 座標、[2,:] に y 座標、[3,:] に z 座標を格納）
    real(kdouble), intent(in) :: coord(:,:)
    !> [out] バウンディングボックス
    real(kdouble), intent(out) :: xmin(3)
    real(kdouble), intent(out) :: xmax(3)

    xmin(1) = minval(coord(1,:))
    xmin(2) = minval(coord(2,:))
    xmin(3) = minval(coord(3,:))

    xmax(1) = maxval(coord(1,:))
    xmax(2) = maxval(coord(2,:))
    xmax(3) = maxval(coord(3,:))
  end subroutine ggtools_get_lbb

  !> @ingroup aabb
  !> 座標がバウンディングボックスに含まれているか判定
  subroutine monolis_check_coordinates_inside_in_lbb(coord, xmin, xmax, is_inside)
    implicit none
    !> [in] 座標
    real(kdouble), intent(in) :: coord(3)
    !> [in] バウンディングボックス（x_min, x_max, y_min, y_max, z_min, z_max の順に格納）
    real(kdouble), intent(in) :: BB(6)
    !> [out] 内包判定フラグ（内側に含まれていれば .true. を返す）
    logical, intent(out) :: is_inside

    is_inside = .false.

    if(xmin(1) - ths <= coord(1) .and. coord(1) <= xmax(1) + ths .and. &
       xmin(2) - ths <= coord(2) .and. coord(2) <= xmax(2) + ths .and. &
       xmin(3) - ths <= coord(3) .and. coord(3) <= xmax(3) + ths)then
      is_inside = .true.
    endif
  end subroutine monolis_check_coordinates_inside_in_lbb

end module mod_ggtools_lbb
