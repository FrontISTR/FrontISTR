module mod_ggtools_element
  use mod_ggtools_def_element
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup element
  !> エレメントidのリストからエレメント構造体配列のindex配列を取得する関数（順番は保存されない）
subroutine ggtools_element_grep_index_by_id_array(ggtools_element_all,nid,id_array,lid_array)
    implicit none
    !> エレメント構造体
    type(type_ggtools_element), intent(in) :: ggtools_element_all(:)
    !> 取得するエレメントの個数
    integer(kint) :: nid
    !> 取得するエレメントidの 1 次元整数配列
    integer(kint) :: id_array(:)
    !> エレメント構造体配列のindex
    integer(kint) :: lid_array(:)

    integer(kint) :: i, n_element_all, idx, element_id
    integer(kint), allocatable :: index(:), id_array_all(:), id_array_tmp(:)

    n_element_all = size(ggtools_element_all)

    call monolis_alloc_I_1d(index, n_element_all)
    call monolis_alloc_I_1d(id_array_all, n_element_all)
    do i=1,n_element_all
      index(i) = i
      id_array_all(i) = ggtools_element_get_element_id(ggtools_element_all(i))
    enddo

    ! id_array_allをソート
    call monolis_qsort_I_2d(id_array_all, index, 1, n_element_all)

    ! 入力されたid_arrayをソート
    call monolis_alloc_I_1d(id_array_tmp, nid)
    id_array_tmp(:) = id_array(:)
    call monolis_qsort_I_1d(id_array_tmp, 1, nid)

    idx = 1
    do i=1,n_element_all
      element_id = id_array_all(i)
      if ( element_id == id_array_tmp(idx)) then
        lid_array(idx) = index(i)
        idx = idx + 1
        if( idx > nid ) exit
      endif
    enddo
    ! indexの中に見つからないものがあった場合は、その個数を減らしたnidとlid_arrayを返す
    if( idx <= nid ) lid_array(idx:nid) = 0
    nid = idx - 1
  end subroutine

  !> @ingroup element
  !> エレメント構造体をエレメントidのリストで取得する関数
  !> 呼び出しの度にエレメント配列ggtools_element_allのindexをidでソートするので
  !> 巨大なggtools_element_allに対して少ないnidで頻繁に呼び出すことは推奨しない
  subroutine ggtools_element_grep_element_by_id_array(ggtools_element_all,nid,id_array,ggtools_element_out)
    implicit none
    !> エレメント構造体
    type(type_ggtools_element), intent(in) :: ggtools_element_all(:)
    !> 取得するエレメントの個数
    integer(kint) :: nid
    !> 取得するエレメントidの 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
    !> 取得したエレメント構造体配列
    type(type_ggtools_element), allocatable :: ggtools_element_out(:)

    integer(kint) :: i
    integer(kint), allocatable :: lid_array(:)

    call monolis_alloc_I_1d(lid_array, nid)   
    call ggtools_element_grep_index_by_id_array(ggtools_element_all,nid,id_array,lid_array)

    allocate(ggtools_element_out(nid))
    do i=1,nid
      call ggtools_element_copy(ggtools_element_all(lid_array(i)), ggtools_element_out(i))
    enddo
  end subroutine ggtools_element_grep_element_by_id_array

  subroutine ggtools_neighbor_local_bucket_get_elements(n_neighbor, neighbor_pids, &
             & id_array_index, id_array_item, nid, id_array, ggtools_element_out, comm)
    !> 探索半径内の隣接PEの個数
    integer(kint) :: n_neighbor
    !> 探索半径内の隣接PEの番号
    integer(kint), allocatable :: neighbor_pids(:)
    !> 探索半径内の隣接PEに含まれる要素番号配列の添字番号（差分が領域ごとの個数）
    integer(kint), allocatable :: id_array_index(:)
    !> 探索半径内の隣接PEに含まれる要素番号配列
    integer(kint), allocatable :: id_array_item(:)
    !> 取得するエレメントの個数
    integer(kint) :: nid
    !> 取得するエレメントidの 1 次元整数配列
    integer(kint) :: id_array(:)
    !> 取得したエレメント構造体配列
    type(type_ggtools_element), allocatable :: ggtools_element_out(:)
    !> MPI コミュニケータ
    integer(kint), intent(in) :: comm

    integer(kint) :: i, n_element_all, idx, element_id
    integer(kint), allocatable :: index(:), id_array_all(:), id_array_tmp(:)

    call monolis_alloc_I_1d(index, n_element_all)
    call monolis_alloc_I_1d(id_array_all, n_element_all)
    do i=1,n_element_all
      index(i) = i
    enddo

    ! id_array_allをソート
    call monolis_qsort_I_2d(id_array_all, index, 1, n_element_all)

    ! 入力されたid_arrayをソート
    call monolis_alloc_I_1d(id_array_tmp, nid)
    id_array_tmp(:) = id_array(:)
    call monolis_qsort_I_1d(id_array_tmp, 1, nid)

    idx = 1
    do i=1,n_element_all
      element_id = id_array_all(i)
      if ( element_id == id_array_tmp(idx)) then
        idx = idx + 1
        if( idx > nid ) exit
      endif
    enddo
    ! indexの中に見つからないものがあった場合は、その個数を減らしたnidとlid_arrayを返す
    nid = idx - 1
  
  end subroutine

end module
