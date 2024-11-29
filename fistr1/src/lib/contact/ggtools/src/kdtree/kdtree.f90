!> k-d ツリーモジュール
module mod_monolis_utils_kdtree
  use mod_monolis_utils_define_prm
  use mod_monolis_utils_error
  use mod_monolis_utils_aabb
  use mod_monolis_utils_alloc
  use mod_monolis_utils_std_sort_R
  use mod_monolis_utils_std_sort_I

  implicit none

  private

  public :: monolis_kdtree_structure
  public :: monolis_kdtree_init_by_BB
  public :: monolis_kdtree_get_BB_including_coordinates
  public :: monolis_kdtree_finalize

  !> k-d ツリー構造体
  type monolis_kdtree_structure_main
    !> 評価・分割する軸番号（1:x 軸、2:y 軸、3:z 軸）
    integer(kint) :: dim
    !> 格納する id 番号
    integer(kint) :: id
    !> root のバウンディングボックス（x_min, x_max, y_min, y_max, z_min, z_max の順に格納）
    real(kdouble) :: BB(6)
    !> root 以下全てのバウンディングボックス（x_min, x_max, y_min, y_max, z_min, z_max の順に格納）
    real(kdouble) :: BB_all(6)
    !> BB の中間座標
    real(kdouble) :: mid_coord(3)
    !> 中間点より座標が小さい集合の k-d ツリー構造体
    type(monolis_kdtree_structure_main), pointer :: left => null()
    !> 中間点より座標が大きい集合の k-d ツリー構造体
    type(monolis_kdtree_structure_main), pointer :: right => null()
  end type monolis_kdtree_structure_main

  !> k-d ツリー構造体
  type monolis_kdtree_structure
    !> k-d ツリー構造体（メイン）
    type(monolis_kdtree_structure_main), pointer :: kdtree => null()
  end type monolis_kdtree_structure

contains

  !> @ingroup kdtree
  !> バウンディングボックス情報を入力して k-d ツリーを構築
  subroutine monolis_kdtree_init_by_BB(monolis_kdtree, n_BB, BB_id, BB)
    implicit none
    !> [in,out] k-d ツリー構造体
    type(monolis_kdtree_structure), intent(inout) :: monolis_kdtree
    !> [in,out] バウンディングボックスの入力数
    integer(kint), intent(inout) :: n_BB
    !> [in,out] バウンディングボックスの id（サイズ [n_BB]）
    integer(kint), intent(inout) :: BB_id(:)
    !> [in,out] バウンディングボックスの入力座標（サイズ [6, n_BB]、x_min, x_max, y_min, y_max, z_min, z_max の順に格納）
    real(kdouble), intent(inout) :: BB(:,:)

    if(associated(monolis_kdtree%kdtree))then
      call monolis_std_error_string("monolis_kdtree_init_by_BB")
      call monolis_std_error_string("monolis_kdtree is already allocated or initialized")
      call monolis_std_error_stop()
    endif

    call monolis_kdtree_init_by_BB_main(monolis_kdtree%kdtree, n_BB, BB_id, BB, 1)
  end subroutine monolis_kdtree_init_by_BB

  !> @ingroup kdtree
  !> 座標を入力して座標を含むバウンディングボックスを取得
  subroutine monolis_kdtree_get_BB_including_coordinates(monolis_kdtree, pos, n_BB, BB_id)
    implicit none
    !> [in,out] k-d ツリー構造体
    type(monolis_kdtree_structure), intent(inout) :: monolis_kdtree
    !> [in] 入力座標
    real(kdouble), intent(in) :: pos(3)
    !> [out] 座標を含むバウンディングボックスの数
    integer(kint), intent(out) :: n_BB
    !> [out] 座標を含むバウンディングボックスの id（サイズ [n_BB]）
    integer(kint), allocatable, intent(out) :: BB_id(:)

    n_BB = 0
    call monolis_dealloc_I_1d(BB_id)

    call monolis_kdtree_get_BB_including_coordinates_main(monolis_kdtree%kdtree, pos, n_BB, BB_id)
  end subroutine monolis_kdtree_get_BB_including_coordinates

  !> @ingroup kdtree
  !> k-d ツリー構造体の終了処理
  subroutine monolis_kdtree_finalize(monolis_kdtree)
    implicit none
    !> [in,out] k-d ツリー構造体
    type(monolis_kdtree_structure), intent(inout) :: monolis_kdtree

    call monolis_kdtree_finalize_main(monolis_kdtree%kdtree)
  end subroutine monolis_kdtree_finalize

  !> @ingroup dev_kdtree
  !> バウンディングボックス情報を入力して k-d ツリーを構築（メイン関数）
  recursive subroutine monolis_kdtree_init_by_BB_main(kdtree, n_BB, BB_id, BB, depth)
    implicit none
    !> [out] k-d ツリー構造体
    type(monolis_kdtree_structure_main), pointer, intent(out) :: kdtree
    !> [in] バウンディングボックスの入力数
    integer(kint), intent(in) :: n_BB
    !> [in,out] バウンディングボックスの id（サイズ [n_BB]）
    integer(kint), intent(inout) :: BB_id(:)
    !> [in,out] バウンディングボックスの入力座標（サイズ [6, n_BB]、x_min, x_max, y_min, y_max, z_min, z_max の順に格納）
    real(kdouble), intent(inout) :: BB(:,:)
    !> [in] k-d ツリーの深さ
    integer(kint), intent(in) :: depth
    integer(kint) :: i, mid_id, iSl, iEl, iSr, iEr
    integer(kint), allocatable :: perm(:)
    real(kdouble), allocatable :: mid_coord(:,:)

    if(n_BB <= 0) return

    allocate(kdtree)

    !> 評価軸の決定
    kdtree%dim = 3 - mod(depth, 3)

    !> 保持している全てのデータ全体の BB を取得
    call monolis_get_aabb_from_BB(BB, kdtree%BB_all)

    !> 保持しているデータごとに中間座標を取得しソート
    call monolis_alloc_R_2d(mid_coord, 3, n_BB)

    do i = 1, n_BB
      mid_coord(1,i) = 0.5d0*(BB(1,i) + BB(2,i))
      mid_coord(2,i) = 0.5d0*(BB(3,i) + BB(4,i))
      mid_coord(3,i) = 0.5d0*(BB(5,i) + BB(6,i))
    enddo

    call monolis_alloc_I_1d(perm, n_BB)
    call monolis_get_sequence_array_I(perm, n_BB, 1, 1)
    call monolis_qsort_R_1d_I_1d(mid_coord(kdtree%dim,:), perm, 1, n_BB)

    !> BB id、BB の置換
    call monolis_perm_array_I(BB_id, perm, n_BB)
    call monolis_perm_array_R(BB(1,:), perm, n_BB)
    call monolis_perm_array_R(BB(2,:), perm, n_BB)
    call monolis_perm_array_R(BB(3,:), perm, n_BB)
    call monolis_perm_array_R(BB(4,:), perm, n_BB)
    call monolis_perm_array_R(BB(5,:), perm, n_BB)
    call monolis_perm_array_R(BB(6,:), perm, n_BB)

    !> 中点座標を計算し、左右で再帰的に k-d ツリーを構築
    mid_id = max(1, n_BB / 2)
    kdtree%id = BB_id(mid_id)
    kdtree%BB = BB(:,mid_id)

    iSl = 1
    iEl = mid_id - 1
    iSr = mid_id + 1
    iEr = n_BB

    call monolis_kdtree_init_by_BB_main(kdtree%left,  iEl - iSl + 1, BB_id(iSl:iEl), BB(:,iSl:iEl), depth + 1)
    call monolis_kdtree_init_by_BB_main(kdtree%right, iEr - iSr + 1, BB_id(iSr:iEr), BB(:,iSr:iEr), depth + 1)
  end subroutine monolis_kdtree_init_by_BB_main

  !> @ingroup dev_kdtree
  !> 座標を入力して座標を含むバウンディングボックスを取得（メイン関数）
  recursive subroutine monolis_kdtree_get_BB_including_coordinates_main(kdtree, pos, n_BB, BB_id)
    implicit none
    !> [in,out] k-d ツリー構造体
    type(monolis_kdtree_structure_main), pointer, intent(inout) :: kdtree
    !> [in] 入力座標
    real(kdouble), intent(in) :: pos(3)
    !> [in,out] 座標を含むバウンディングボックスの数
    integer(kint), intent(inout) :: n_BB
    !> [in,out] 座標を含むバウンディングボックスの id（サイズ [n_BB]）
    integer(kint), intent(inout), allocatable :: BB_id(:)
    integer(kint) :: iadd(1)
    real(kdouble) :: ths = 1.0d-8
    logical :: is_inside

    if(.not. associated(kdtree)) return

    !> k-d ツリーの root のバウンディングボックスに内包されていれば、戻り値に追加
    call monolis_check_inside_in_aabb(pos, kdtree%BB, ths, is_inside)

    if(is_inside)then
      n_BB = n_BB + 1
      iadd(1) = kdtree%id
      call monolis_append_I_1d(BB_id, 1, iadd)
    endif

    !> 入力座標が child の全てのバウンディングボックスに内包されていれば検索を継続
    if(associated(kdtree%left))then
      if(kdtree%left%BB_all(2*kdtree%dim - 1) - ths <= pos(kdtree%dim) .and. &
         pos(kdtree%dim) <= kdtree%left%BB_all(2*kdtree%dim) + ths)then
        call monolis_kdtree_get_BB_including_coordinates_main(kdtree%left , pos, n_BB, BB_id)
      endif
    endif

    if(associated(kdtree%right))then
      if(kdtree%right%BB_all(2*kdtree%dim - 1) - ths <= pos(kdtree%dim) .and. &
         pos(kdtree%dim) <= kdtree%right%BB_all(2*kdtree%dim) + ths)then
        call monolis_kdtree_get_BB_including_coordinates_main(kdtree%right, pos, n_BB, BB_id)
      endif
    endif
  end subroutine monolis_kdtree_get_BB_including_coordinates_main

  !> @ingroup dev_kdtree
  !> k-d ツリー構造体の終了処理（メイン関数）
  recursive subroutine monolis_kdtree_finalize_main(kdtree)
    implicit none
    !> [in,out] k-d ツリー構造体
    type(monolis_kdtree_structure_main), pointer, intent(inout) :: kdtree

    if(associated(kdtree%left))then
      call monolis_kdtree_finalize_main(kdtree%left)
      deallocate(kdtree%left)
      kdtree%left => null()
    endif

    if(associated(kdtree%right))then
      call monolis_kdtree_finalize_main(kdtree%right)
      deallocate(kdtree%right)
      kdtree%right => null()
    endif
  end subroutine monolis_kdtree_finalize_main
end module mod_monolis_utils_kdtree
