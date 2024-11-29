!> ローカルバケットセルモジュール
module mod_ggtools_local_bucket
  use mod_monolis_utils
  use mod_ggtools_def_local_bucket
  use mod_ggtools_bucket_util
  use mod_monolis_utils
  implicit none

contains

  !> @ingroup bucket
  !> 座標で指定したLBBを番号の範囲に変換する
  subroutine ggtools_local_bucket_xminmax_to_iminmax(ggtools_bucket, xmin_local, xmax_local, imin, imax)
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> ローカルバウンダリボックス LBB の最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンダリボックス LBB の最大座標
    real(kdouble) :: xmax_local(3)
    !> ローカルバウンダリボックス LBB の最小座標のインデックス（x,y,z形式）
    integer(kint) :: imin(3)
    !> ローカルバウンダリボックス LBB の最大座標のインデックス（x,y,z形式）
    integer(kint) :: imax(3)

    real(kdouble) :: x_point(3)

    x_point(1) = xmin_local(1) - ths
    x_point(2) = xmin_local(2) - ths
    x_point(3) = xmin_local(3) - ths
    imin = ggtools_bucket_get_integer_coordinate_from_real_coordinate(ggtools_bucket%xmin, ggtools_bucket%dx, x_point)
    ! thsの影響でx_pointがバケット範囲外に出た場合は端のバケットに引き戻す
    imin(1) = min(imin(1),1)
    imin(2) = min(imin(2),1)
    imin(3) = min(imin(3),1)

    x_point(1) = xmax_local(1) + ths
    x_point(2) = xmax_local(2) + ths
    x_point(3) = xmax_local(3) + ths
    imax = ggtools_bucket_get_integer_coordinate_from_real_coordinate(ggtools_bucket%xmin, ggtools_bucket%dx, x_point)
    ! thsの影響でx_pointがバケット範囲外に出た場合は端のバケットに引き戻す
    imax(1) = min(imax(1),ggtools_bucket%nx(1))
    imax(2) = min(imax(2),ggtools_bucket%nx(2))
    imax(3) = min(imax(3),ggtools_bucket%nx(3))

  end subroutine

  !> @ingroup bucket
  !> LBB を指定してバケットセルに id を登録する関数
  subroutine ggtools_local_bucket_set_id_by_lbb(ggtools_bucket, ggtools_bucket_cell, xmin_local, xmax_local, id)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> ローカルバウンダリボックス LBB の最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンダリボックス LBB の最大座標
    real(kdouble) :: xmax_local(3)
    !> バケットセルに登録する id
    integer(kint) :: id, x, y, z
    integer(kint) :: int_id(3), in, imin(3), imax(3)

    call ggtools_local_bucket_xminmax_to_iminmax(ggtools_bucket, xmin_local, xmax_local, imin, imax)

    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      int_id(1) = x
      int_id(2) = y
      int_id(3) = z
      in = ggtools_bucket_get_1d_array_index_from_integer_coodinate(ggtools_bucket%nx, int_id)
      call ggtools_local_bucket_set_id_main(ggtools_bucket_cell(in), id)
    enddo
    enddo
    enddo
  end subroutine ggtools_local_bucket_set_id_by_lbb

  !> @ingroup bucket
  !> 座標を指定して 1 つのバケットセルに id を登録する関数
  subroutine ggtools_local_bucket_set_id_by_point(ggtools_bucket, ggtools_bucket_cell, x_point, id)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バケットセルに登録する id
    integer(kint) :: id
    integer(kint) :: in
    integer(kint) :: int_id(3)

    int_id = ggtools_bucket_get_integer_coordinate_from_real_coordinate(ggtools_bucket%xmin, ggtools_bucket%dx, x_point)
    in = ggtools_bucket_get_1d_array_index_from_integer_coodinate(ggtools_bucket%nx, int_id)
    call ggtools_local_bucket_set_id_main(ggtools_bucket_cell(in), id)
  end subroutine ggtools_local_bucket_set_id_by_point

  !> @ingroup dev
  !> バケットへの情報登録（メイン関数）
  subroutine ggtools_local_bucket_set_id_main(ggtools_bucket_cell, data)
    implicit none
    !> バケット検索構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell
    !> 登録する要素領域 id
    integer(kint) :: data
    integer(kint) :: add(1)

    add = data
    call monolis_append_I_1d(ggtools_bucket_cell%id_array, 1, add)
    ggtools_bucket_cell%nid = ggtools_bucket_cell%nid + 1
  end subroutine ggtools_local_bucket_set_id_main

  !> @ingroup bucket
  !> バケットに登録した情報の取得
  subroutine ggtools_local_bucket_get_id_main(ggtools_bucket_cell, nid, id_array)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)

    nid = ggtools_bucket_cell%nid
    call monolis_realloc_I_1d(id_array, nid)
    id_array(1:nid) = ggtools_bucket_cell%id_array(1:nid)
  end subroutine

  !> @ingroup bucket
  !> 入力した座標を内包するバケットセルから、nid 個の整数配列 id_array を取得する関数
  subroutine ggtools_local_bucket_get_id_by_point(ggtools_bucket, ggtools_bucket_cell, x_point, nid, id_array)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
    integer(kint) :: int_id(3), index

    int_id = ggtools_bucket_get_integer_coordinate_from_real_coordinate(ggtools_bucket%xmin, ggtools_bucket%dx, x_point)
    index = ggtools_bucket_get_1d_array_index_from_integer_coodinate(ggtools_bucket%nx, int_id)
    call ggtools_local_bucket_get_id_main(ggtools_bucket_cell(index), nid, id_array)
  end subroutine ggtools_local_bucket_get_id_by_point

  !> @ingroup bucket
  !> 入力した LBB を内包するバケットセルから、nid 個の整数配列 id_array を取得する関数
  subroutine ggtools_local_bucket_get_id_by_lbb(ggtools_bucket, ggtools_bucket_cell, xmin_local, xmax_local, nid, id_array)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> ローカルバウンダリボックス LBB の最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンダリボックス LBB の最大座標
    real(kdouble) :: xmax_local(3)
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
    integer(kint) :: x, y, z, int_id(3), imin(3), imax(3), nid_temp, nid_all, index
    integer(kint), allocatable :: id_array_temp(:), id_array_all(:)

    call ggtools_local_bucket_xminmax_to_iminmax(ggtools_bucket, xmin_local, xmax_local, imin, imax)

    nid_all = 0

    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      int_id(1) = x
      int_id(2) = y
      int_id(3) = z
      index = ggtools_bucket_get_1d_array_index_from_integer_coodinate(ggtools_bucket%nx, int_id)
      call ggtools_local_bucket_get_id_main(ggtools_bucket_cell(index), nid_temp, id_array_temp)
      call monolis_append_I_1d(id_array_all, nid_temp, id_array_temp)
      nid_all = nid_all + nid_temp
    enddo
    enddo
    enddo

    if(nid_all == 0)then
      nid = 0
      return
    endif

    call monolis_qsort_I_1d(id_array_all, 1, nid_all)
    call monolis_get_uniq_array_I(id_array_all, nid_all, nid)

    call monolis_alloc_I_1d(id_array, nid)
    id_array(:) = id_array_all(1:nid)
  end subroutine ggtools_local_bucket_get_id_by_lbb

  !> @ingroup bucket
  !> 入力した座標を内包するバケットセルから、nid 個の整数配列 id_array を取得する関数
  subroutine ggtools_local_bucket_get_id_by_1d_array_index(ggtools_bucket_cell, index, nid, id_array)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> バケットセルの 1 次元配列通し番号
    integer(kint) :: index
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)

    if(ggtools_bucket_cell(index)%nid == 0)then
      nid = 0
      return
    endif

    nid = ggtools_bucket_cell(index)%nid
    call monolis_alloc_I_1d(id_array, nid)
    id_array(:) = ggtools_bucket_cell(index)%id_array(:)
  end subroutine ggtools_local_bucket_get_id_by_1d_array_index

  !> @ingroup bucket
  !> 検索座標・検索半径を入力し検索半径内に含まれるローカルバケットセルの個数とローカルバケットセル番号配列を取得する関数
  subroutine ggtools_local_bucket_get_id_by_radius(ggtools_bucket, x_point, radius, &
      & n_bucketcell, bucketcell_id_array)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バゲットセルの探索半径
    real(kdouble) :: radius
    !> 探索半径内のバケットセルの個数
    integer(kint) :: n_bucketcell
    !> 探索半径内のバケットセルの番号
    integer(kint), allocatable :: bucketcell_id_array(:)

    integer(kint) :: idof, x, y, z
    integer(kint) :: int_id(3), index, imin(3), imax(3)
    real(kdouble)  :: xmin_local(3), xmax_local(3)
    real(kdouble)  :: xmin(3), xmax(3),dx(3), center(3)
    real(kdouble)  :: bkradius, distance

    !> 検索球のLBB（x_point ± radiusの範囲）でローカルバケットセルを絞り込む
    xmin_local(1:3) = x_point(1:3) - radius
    xmax_local(1:3) = x_point(1:3) + radius
    call ggtools_local_bucket_xminmax_to_iminmax(ggtools_bucket, xmin_local, xmax_local, imin, imax)    

    !> 絞り込んだLBBのうち、x_pointからのユークリッド距離がradiusに届かないものを除外する
    !!> バケットの座標計算用
    call ggtools_local_bucket_get_size_of_one_bucket_cell(ggtools_bucket, dx)
    call ggtools_local_bucket_get_bucket_size(ggtools_bucket, xmin, xmax)
    bkradius = 0.5d0*dsqrt(dot_product(dx,dx))
    !!> 対象バケットで距離判定。球体と直方体の厳密な距離判定ではなく、高速に枝刈りする目的で
    !!> 自明に遠いバケット＝x_pointから中心への距離がradius+バケット外接球の半径に満たないものを除外する
    n_bucketcell=0
    call monolis_realloc_I_1d(bucketcell_id_array, 100)
    do z = imin(3), imax(3)
    do y = imin(2), imax(2)
    do x = imin(1), imax(1)
      int_id(1) = x
      int_id(2) = y
      int_id(3) = z
      do idof = 1, 3
        center(idof) = xmin(idof) + (dble(int_id(idof))-0.5d0)*dx(idof)
      enddo
      distance = dsqrt(dot_product(center-x_point,center-x_point))
      if( distance > radius + bkradius ) cycle

      index = ggtools_bucket_get_1d_array_index_from_integer_coodinate(ggtools_bucket%nx, int_id)
      n_bucketcell = n_bucketcell + 1
      if( n_bucketcell > size(bucketcell_id_array) ) then
        call monolis_realloc_I_1d(bucketcell_id_array, 2*size(bucketcell_id_array))
      endif
      bucketcell_id_array(n_bucketcell) = index
    enddo
    enddo
    enddo

  end subroutine ggtools_local_bucket_get_id_by_radius

  subroutine ggtools_local_bucket_get_id_by_bucketids(ggtools_bucket_cell, n_bucketcell, bucketcell_id_array, nid, id_array)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> 探索半径内のバケットセルの個数
    integer(kint) :: n_bucketcell
    !> 探索半径内のバケットセルの番号
    integer(kint), allocatable :: bucketcell_id_array(:)
    !> バケットセルに登録された個数
    integer(kint) :: nid
    !> バケットセルに登録された nid 個の 1 次元整数配列（引数が与えられたときだけ返す）
    integer(kint), allocatable, optional :: id_array(:)

    integer(kint) :: i
    integer(kint) :: nid_temp, nid_all, index
    integer(kint), allocatable :: id_array_temp(:), id_array_all(:)

    nid_all = 0
    do i = 1, n_bucketcell
      index = bucketcell_id_array(i)
      call ggtools_local_bucket_get_id_main(ggtools_bucket_cell(index), nid_temp, id_array_temp)
      call monolis_append_I_1d(id_array_all, nid_temp, id_array_temp)
      nid_all = nid_all + nid_temp
    enddo

    if(nid_all == 0)then
      nid = 0
      return
    endif

    call monolis_qsort_I_1d(id_array_all, 1, nid_all)
    call monolis_get_uniq_array_I(id_array_all, nid_all, nid)

    if( .not. present(id_array) ) return
    call monolis_realloc_I_1d(id_array, nid)
    id_array(:) = id_array_all(1:nid)
  end subroutine ggtools_local_bucket_get_id_by_bucketids

  subroutine ggtools_neighbor_local_bucket_get_id_by_radius(gbucket, lbucket, x_point, radius, n_neighbor, &
    & neighbor_pids, nbucketcell_neighbor_index, nbucketcell_neighbor_item, comm)
    !> グローバルバケット構造体
    type(type_ggtools_global_bucket) :: gbucket
    !> ローカルバケット構造体
    type(type_ggtools_local_bucket) :: lbucket
    !> バゲットセルを指定する座標
    real(kdouble) :: x_point(3)
    !> バゲットセルの探索半径
    real(kdouble) :: radius
    !> 探索半径内の隣接PEの個数
    integer(kint) :: n_neighbor
    !> 探索半径内の隣接PEの番号
    integer(kint), allocatable :: neighbor_pids(:)
    !> 探索半径内の隣接PEに含まれるバケットセルの番号の添字番号（差分が領域ごとの個数）
    integer(kint), allocatable :: nbucketcell_neighbor_index(:)
    !> 探索半径内の隣接PEに含まれるバケットセルの番号
    integer(kint), allocatable :: nbucketcell_neighbor_item(:)
    !> MPI コミュニケータ
    integer(kint), intent(in) :: comm

    integer(kint) :: pid, my_rank, n_local
    integer(kint) :: int_id(3), index, imin(3), imax(3)
    real(kdouble)  :: xmin_local(3), xmax_local(3)
    !> ローカルバケット構造体（他領域探索用）
    type(type_ggtools_local_bucket) :: lbucket_neighbor
    integer(kint) :: n_bucketcell
    integer(kint), allocatable :: bucketcell_id_array(:)

    n_local = monolis_mpi_get_global_comm_size()
    my_rank = monolis_mpi_get_global_my_rank()

    call monolis_realloc_I_1d(neighbor_pids, n_local)
    call monolis_realloc_I_1d(nbucketcell_neighbor_index, n_local+1)

    n_neighbor = 0
    nbucketcell_neighbor_index(1) = 0
    do pid=0,n_local-1
      if(pid == my_rank) cycle
      call ggtools_global_bucket_get_local_bb(gbucket, pid, xmin_local, xmax_local)
      call ggtools_local_bucket_init_by_global_bucket(gbucket, lbucket_neighbor)
      call ggtools_local_bucket_get_id_by_radius(lbucket_neighbor, x_point, radius, &
          & n_bucketcell, bucketcell_id_array)
      if( n_bucketcell == 0 ) cycle
      n_neighbor = n_neighbor + 1
      nbucketcell_neighbor_index(n_neighbor+1) = nbucketcell_neighbor_index(n_neighbor) + n_bucketcell
      call monolis_append_I_1d(nbucketcell_neighbor_item, n_bucketcell, bucketcell_id_array)
    enddo
  end subroutine

  subroutine ggtools_neighbor_local_bucket_get_id_by_bucketids(gbucket, lbucket, n_neighbor, &
    & neighbor_pids, nbucketcell_neighbor_index, nbucketcell_neighbor_item, &
    & id_array_index, comm, id_array_item)
    !> グローバルバケット構造体
    type(type_ggtools_global_bucket) :: gbucket
    !> ローカルバケット構造体
    type(type_ggtools_local_bucket) :: lbucket
    !> 探索半径内の隣接PEの個数
    integer(kint) :: n_neighbor
    !> 探索半径内の隣接PEの番号
    integer(kint), allocatable :: neighbor_pids(:)
    !> 探索半径内の隣接PEに含まれるバケットセルの番号の添字番号（差分が領域ごとの個数）
    integer(kint), allocatable :: nbucketcell_neighbor_index(:)
    !> 探索半径内の隣接PEに含まれるバケットセルの番号
    integer(kint), allocatable :: nbucketcell_neighbor_item(:)
    !> 探索半径内の隣接PEに含まれる要素番号配列の添字番号（差分が領域ごとの個数）
    integer(kint), allocatable :: id_array_index(:)
    !> MPI コミュニケータ
    integer(kint), intent(in) :: comm
    !> 探索半径内の隣接PEに含まれる要素番号配列
    integer(kint), allocatable, optional :: id_array_item(:)

    integer(kint) :: i, pid, my_rank, n_local
    integer(kint) :: int_id(3), index, imin(3), imax(3)
    real(kdouble)  :: xmin_local(3), xmax_local(3)
    !> ローカルバケット構造体（他領域探索用）
    type(type_ggtools_local_bucket) :: lbucket_neighbor
    integer(kint) :: n_bucketcell
    integer(kint), allocatable :: bucketcell_id_array(:)

    n_local = monolis_mpi_get_global_comm_size()
    my_rank = monolis_mpi_get_global_my_rank()

    n_neighbor = 0
    nbucketcell_neighbor_index(1) = 0
    do i=1,n_neighbor
      if(pid == my_rank) cycle
      call ggtools_global_bucket_get_local_bb(gbucket, pid, xmin_local, xmax_local)
      call ggtools_local_bucket_init_by_global_bucket(gbucket, lbucket_neighbor)
      n_neighbor = n_neighbor + 1
      nbucketcell_neighbor_index(n_neighbor+1) = nbucketcell_neighbor_index(n_neighbor) + n_bucketcell
      call monolis_append_I_1d(nbucketcell_neighbor_item, n_bucketcell, bucketcell_id_array)
    enddo
  end subroutine


end module mod_ggtools_local_bucket
