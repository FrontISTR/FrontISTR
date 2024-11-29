!> グローバルバケットセルテストモジュール
module mod_ggtools_global_bucket_test
    use mod_monolis_utils
    use mod_ggtools
    implicit none
  
  contains
  
    subroutine ggtools_global_bucket_test()
      implicit none
  
      !# バケット検索用テスト
      call ggtools_global_bucket_test_search_serial()
      call ggtools_global_bucket_test_search_2mpi()
      !call ggtools_bucket_initialize_test()
      !call ggtools_bucket_search_test()
  
      !# 一時的なテスト
      !call monolis_utils_aabb_test()
      !call monolis_utils_kdtree_test()
    end subroutine ggtools_global_bucket_test

    subroutine ggtools_global_bucket_test_search_serial()
      implicit none

      !# 節点要素
      type(type_ggtools_element), allocatable :: slaves(:)
      type(type_ggtools_element), allocatable :: masters(:)
      type(type_ggtools_element), allocatable :: masters_out(:)

      !# バケット
      type(type_ggtools_global_bucket) :: gbucket
      type(type_ggtools_local_bucket)  :: lbucket
      type(type_ggtools_local_bucket_cells), allocatable :: lbucket_cells(:)
      real(kdouble)  :: x_min(3), x_max(3)
      real(kdouble)  :: dx(3)
      integer(kint)  :: nx(3), icel, id, n_local
      real(kdouble)  :: xmin_local(3), xmax_local(3)
      integer(kint)  :: n_bucketcell, nid
      integer(kint), allocatable  :: bucketcell_id_array(:), id_array(:)
      real(kdouble)  :: x_point(3), radius

      n_local = monolis_mpi_get_global_comm_size()
      if( n_local > 1 ) return

      !# 0. 節点要素の定義
      !# 1. 各分割領域で、ローカルバケットのx_min、x_max を取得
      call ggtools_global_bucket_test_search_setelements(slaves, masters, x_min, x_max)

      !# 2. ローカルバケットのx_min、x_max を⽤いてグローバルバケット構造体を初期化
      dx(1) = 5.2d0/5.d0  ! 5分割
      dx(2) = 1.2d0/2.d0  ! 2分割
      dx(3) = 0.2d0       ! 1分割
      call ggtools_global_bucket_init_by_local_bucket(gbucket, x_min, x_max, dx, MPI_COMM_WORLD)

      !# 3. ローカルバケットのx_min、x_max、グローバルバケット構造体を⽤いて
      !#    ローカルバケット構造体を初期化
      call ggtools_local_bucket_init_by_global_bucket(gbucket, lbucket)

      !# 4. ローカルバケット構造体を⽤いて、ローカルバケットセル構造体を初期化
      call ggtools_local_bucket_get_number_of_bucket_divisions(lbucket, nx)
      call ggtools_local_bucket_cells_init(lbucket_cells, nx)

      !# 5-A. ローカルバケットの登録先指定に利⽤する座標と登録する整数id を⼊⼒し、
      !#      ローカルバケットセルに情報を登録
      !# slave要素（節点）
      do icel=1,size(slaves)
        call ggtools_element_get_element_lbb(slaves(icel),xmin_local,xmax_local)
        id = ggtools_element_get_element_id(slaves(icel))
        call ggtools_local_bucket_set_id_by_point(lbucket, lbucket_cells, xmin_local, id)
      enddo
      !# 5-B. ローカルバケットの登録先指定に利用する
      !#      エレメントバウンディングボックスの座標x_min_EBB、x_max_EBB を取得。
      !#      登録する整数id とともに⼊⼒し、ローカルバケットセルに情報を登録
      !# master要素
      do icel=1,size(masters)
        call ggtools_element_get_element_lbb(masters(icel),xmin_local,xmax_local)
        id = ggtools_element_get_element_id(masters(icel))
        call ggtools_local_bucket_set_id_by_lbb(lbucket, lbucket_cells, xmin_local, xmax_local, id)
      enddo
      
      !自領域のデータ取得関数（基本関数）
      !前提：隣接するローカルバケットの個数・ランク番号、検索半径内の登録id 個数・id 配列、
      !エレメントおよびエレメントに付随する整数・実数情報の取得が最⼩構成
      !2-A. 検索座標・検索半径を入力し、
      !     検索半径内に含まれるローカルバケットセルの個数とローカルバケットセル番号配列を取得
      x_point(1:3) = [1.d0, 0.5d0, 0.d0]
      radius = 1.d0
      call ggtools_local_bucket_get_id_by_radius(lbucket, x_point, radius, n_bucketcell, bucketcell_id_array)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius nid", n_bucketcell, 9)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius bucketcell_id_array(1)", bucketcell_id_array(1), 1)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius bucketcell_id_array(9)", bucketcell_id_array(9), 13)
      !2-B. ローカルバケットセルの個数とローカルバケットセル番号配列を入力し、
      !     バケットセル内に含まれる重複のない登録 id の個数を取得
      call ggtools_local_bucket_get_id_by_bucketids(lbucket_cells, n_bucketcell, bucketcell_id_array, nid)
      !2-C. ローカルバケットセルの個数とローカルバケットセル番号配列を入力し、
      !     バケットセル内に含まれる重複のない登録 id の個数と登録id 配列を取得
      call ggtools_local_bucket_get_id_by_bucketids(lbucket_cells, n_bucketcell, bucketcell_id_array, nid, id_array)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius nid", nid, 8)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius bucketcell_id_array(1)", id_array(1), 1)
      call monolis_test_check_eq_I1("ggtools_local_bucket_get_id_by_radius bucketcell_id_array(6)", id_array(8), 10)
      !2-D. データが格納されたエレメント構造体配列、登録 id の個数と登録id 配列を入力し、
      !     指定id 配列で指定されたエレメント構造体配列を取得
      nid = 2
      id_array = [8, 6]
      call ggtools_element_grep_element_by_id_array(masters,nid,id_array,masters_out)
      call monolis_test_check_eq_I1("ggtools_element_grep_element_by_id_array eid1", masters_out(1)%element_id, 6)
      call monolis_test_check_eq_I1("ggtools_element_grep_element_by_id_array eid2", masters_out(2)%element_id, 8)

    end subroutine ggtools_global_bucket_test_search_serial

    subroutine ggtools_global_bucket_test_search_setelements(slaves,masters,xmin,xmax)
      type(type_ggtools_element), allocatable :: slaves(:)
      type(type_ggtools_element), allocatable :: masters(:)
      real(kdouble)  :: xmin(3)
      real(kdouble)  :: xmax(3)

      real(kdouble)  :: slave_nodes(3,5)
      real(kdouble)  :: master_nodes(3,12)
      real(kdouble), allocatable  :: coordinate(:,:)
      integer(kint)  :: n_slaves, n_masters, conn(4,5)
      integer(kint)  :: slave_list(5), master_list(5)
      integer(kint)  :: slave_ids(5), master_ids(5)
      integer(kint)  :: i, icel, inod, my_rank, n_local, idof

      n_local = monolis_mpi_get_global_comm_size()
      my_rank = monolis_mpi_get_global_my_rank()

      if( n_local == 1 ) then
        n_slaves = 5
        slave_list = [1,2,3,4,5]
        n_masters = 5
        master_list = [1,2,3,4,5]
      else if( n_local == 2 ) then
        if( my_rank == 0 ) then
          n_slaves = 2
          slave_list = [1,2,0,0,0]
          n_masters = 2
          master_list = [1,2,0,0,0]
        else if( my_rank == 1 ) then
          n_slaves = 3
          slave_list = [3,4,5,0,0]
          n_masters = 3
          master_list = [3,4,5,0,0]
        endif
      else
        call monolis_std_error_stop()
      endif

      !# 節点要素の定義
      slave_ids(1:5) = [1,2,3,4,5] 
      slave_nodes(1:3,1) = [0.5d0, 0.d0, 0.d0]
      slave_nodes(1:3,2) = [1.5d0, 0.d0, 0.d0]
      slave_nodes(1:3,3) = [2.5d0, 0.d0, 0.d0]
      slave_nodes(1:3,4) = [3.5d0, 0.d0, 0.d0]
      slave_nodes(1:3,5) = [4.5d0, 0.d0, 0.d0]

      master_ids(1:5) = [6,7,8,9,10] 
      master_nodes(1:3, 1) = [0.0d0, 0.0d0, 0.d0]
      master_nodes(1:3, 2) = [0.0d0, 1.0d0, 0.d0]
      master_nodes(1:3, 3) = [1.0d0, 0.0d0, 0.d0]
      master_nodes(1:3, 4) = [1.0d0, 1.0d0, 0.d0]
      master_nodes(1:3, 5) = [2.0d0, 0.0d0, 0.d0]
      master_nodes(1:3, 6) = [2.0d0, 1.0d0, 0.d0]
      master_nodes(1:3, 7) = [3.0d0, 0.0d0, 0.d0]
      master_nodes(1:3, 8) = [3.0d0, 1.0d0, 0.d0]
      master_nodes(1:3, 9) = [4.0d0, 0.0d0, 0.d0]
      master_nodes(1:3,10) = [4.0d0, 1.0d0, 0.d0]
      master_nodes(1:3,11) = [5.0d0, 0.0d0, 0.d0]
      master_nodes(1:3,12) = [5.0d0, 1.0d0, 0.d0]

      conn(1:4,1) = [1, 2, 4, 3]
      conn(1:4,2) = [3, 4, 6, 5]
      conn(1:4,3) = [5, 6, 8, 7]
      conn(1:4,4) = [7, 8,10, 9]
      conn(1:4,5) = [9,10,12,11]

      call monolis_alloc_R_2d(coordinate,3,4)

      xmin(:) = 9999.d0
      xmax(:) = -9999.d0

      !slave要素（点要素）の定義
      allocate(slaves(n_slaves))
      do i=1,n_slaves
        icel = slave_list(i)
        coordinate(1:3,1) = slave_nodes(1:3,icel)
        call ggtools_element_init(slaves(i), 1, coordinate, my_rank, slave_ids(icel))
        do idof=1,3
          if( coordinate(idof,1) < xmin(idof) ) xmin(idof) = coordinate(idof,1)
          if( coordinate(idof,1) > xmax(idof) ) xmax(idof) = coordinate(idof,1)
        enddo
      enddo

      !master要素（四角形面要素）の定義
      allocate(masters(n_masters))
      do i=1,n_masters
        icel = master_list(i)
        do inod=1,4
          coordinate(1:3,inod) = master_nodes(1:3, conn(inod,icel))
          do idof=1,3
            if( coordinate(idof,inod) < xmin(idof) ) xmin(idof) = coordinate(idof,inod)
            if( coordinate(idof,inod) > xmax(idof) ) xmax(idof) = coordinate(idof,inod)
          enddo
        enddo
        call ggtools_element_init(masters(i), 4, coordinate, my_rank, master_ids(icel))
      enddo

      ! xmin, xmaxを拡大
      xmin(:) = xmin(:) - 0.1d0
      xmax(:) = xmax(:) + 0.1d0
    end subroutine

    subroutine ggtools_global_bucket_test_search_2mpi()
      implicit none

      !# 節点要素
      type(type_ggtools_element), allocatable :: slaves(:)
      type(type_ggtools_element), allocatable :: masters(:)
      type(type_ggtools_element), allocatable :: masters_out(:)

      !# バケット
      type(type_ggtools_global_bucket) :: gbucket
      type(type_ggtools_local_bucket)  :: lbucket
      type(type_ggtools_local_bucket_cells), allocatable :: lbucket_cells(:)
      real(kdouble)  :: x_min(3), x_max(3)
      real(kdouble)  :: dx(3)
      integer(kint)  :: nx(3), icel, id, n_local, my_rank
      real(kdouble)  :: xmin_local(3), xmax_local(3)
      integer(kint)  :: n_bucketcell, nid
      integer(kint), allocatable  :: bucketcell_id_array(:), id_array(:)
      integer(kint)  :: n_neighbor
      integer(kint), allocatable  :: neighbor_pids(:), nbucketcell_neighbor_index(:), nbucketcell_neighbor_item(:)
      real(kdouble)  :: x_point(3), radius
      integer(kint), allocatable  :: id_array_index(:), id_array_item(:)

      n_local = monolis_mpi_get_global_comm_size()
      if( n_local /= 2 ) return
      my_rank = monolis_mpi_get_global_my_rank()
      
      !# 0. 節点要素の定義
      !# 1. 各分割領域で、ローカルバケットのx_min、x_max を取得
      call ggtools_global_bucket_test_search_setelements(slaves, masters, x_min, x_max)

      !# 2. ローカルバケットのx_min、x_max を⽤いてグローバルバケット構造体を初期化
      dx(1) = 5.2d0/5.d0  ! 5分割
      dx(2) = 1.2d0/2.d0  ! 2分割
      dx(3) = 0.2d0       ! 1分割
      call ggtools_global_bucket_init_by_local_bucket(gbucket, x_min, x_max, dx, MPI_COMM_WORLD)

      !# 3. ローカルバケットのx_min、x_max、グローバルバケット構造体を⽤いて
      !#    ローカルバケット構造体を初期化
      call ggtools_local_bucket_init_by_global_bucket(gbucket, lbucket)

      !# 4. ローカルバケット構造体を⽤いて、ローカルバケットセル構造体を初期化
      call ggtools_local_bucket_get_number_of_bucket_divisions(lbucket, nx)
      call ggtools_local_bucket_cells_init(lbucket_cells, nx)

      !# 5-A. ローカルバケットの登録先指定に利⽤する座標と登録する整数id を⼊⼒し、
      !#      ローカルバケットセルに情報を登録
      !# slave要素（節点）
      do icel=1,size(slaves)
        call ggtools_element_get_element_lbb(slaves(icel),xmin_local,xmax_local)
        id = ggtools_element_get_element_id(slaves(icel))
        call ggtools_local_bucket_set_id_by_point(lbucket, lbucket_cells, xmin_local, id)
      enddo
      !# 5-B. ローカルバケットの登録先指定に利用する
      !#      エレメントバウンディングボックスの座標x_min_EBB、x_max_EBB を取得。
      !#      登録する整数id とともに入力し、ローカルバケットセルに情報を登録
      !# master要素
      do icel=1,size(masters)
        call ggtools_element_get_element_lbb(masters(icel),xmin_local,xmax_local)
        id = ggtools_element_get_element_id(masters(icel))
        call ggtools_local_bucket_set_id_by_lbb(lbucket, lbucket_cells, xmin_local, xmax_local, id)
      enddo

      !自領域のデータ取得関数（基本関数）
      !前提：隣接するローカルバケットの個数・ランク番号、検索半径内の登録id 個数・id 配列、
      !エレメントおよびエレメントに付随する整数・実数情報の取得が最⼩構成
      !2-A. 検索座標・検索半径を入力し、
      !     検索半径内に含まれるローカルバケットセルの個数とローカルバケットセル番号配列を取得
      x_point(1:3) = [1.d0, 0.5d0, 0.d0]
      radius = 1.d0
      call ggtools_local_bucket_get_id_by_radius(lbucket, x_point, radius, n_bucketcell, bucketcell_id_array)
      !2-B. ローカルバケットセルの個数とローカルバケットセル番号配列を入力し、
      !     バケットセル内に含まれる重複のない登録 id の個数を取得
      call ggtools_local_bucket_get_id_by_bucketids(lbucket_cells, n_bucketcell, bucketcell_id_array, nid)
      !2-C. ローカルバケットセルの個数とローカルバケットセル番号配列を入力し、
      !     バケットセル内に含まれる重複のない登録 id の個数と登録id 配列を取得
      call ggtools_local_bucket_get_id_by_bucketids(lbucket_cells, n_bucketcell, bucketcell_id_array, nid, id_array)
      !2-D. データが格納されたエレメント構造体配列、登録 id の個数と登録id 配列を入力し、
      !     指定id 配列で指定されたエレメント構造体配列を取得
      nid = 2
      id_array = [8, 6]
      call ggtools_element_grep_element_by_id_array(masters,nid,id_array,masters_out)

      !2-E. 検索座標・検索半径を入力し、
      !     検索半径内に含まれる隣接ローカルバケットの個数とランク番号配列を取得
      x_point(1:3) = [1.d0, 0.5d0, 0.d0]
      radius = 1.d0
      call ggtools_neighbor_local_bucket_get_id_by_radius(gbucket, lbucket, x_point, radius, n_neighbor, &
                       & neighbor_pids, nbucketcell_neighbor_index, nbucketcell_neighbor_item, MPI_COMM_WORLD)
      !2-F. 検索座標・検索半径を入力し、隣接ローカルバケット数の大きさの配列に、
      !     それぞれの隣接ローカルバケットの検索半径内に含まれる登録 id の個数が出力
      call ggtools_neighbor_local_bucket_get_id_by_bucketids(gbucket, lbucket, n_neighbor, &
                       & neighbor_pids, nbucketcell_neighbor_index, nbucketcell_neighbor_item, &
                       & id_array_index, MPI_COMM_WORLD)
      !2-G. 検索座標・検索半径を入力し、2-F で得られた登録 id の個数の合計値の大きさの配列に、
      !     検索半径内に含まれる登録id 配列を出力
      call ggtools_neighbor_local_bucket_get_id_by_bucketids(gbucket, lbucket, n_neighbor, &
                       & neighbor_pids, nbucketcell_neighbor_index, nbucketcell_neighbor_item, &
                       & id_array_index, MPI_COMM_WORLD, id_array_item)
      !2-H. 送受信したいデータ配列（整数型・実数型・エレメント構造体配列型）・取得データの個数・
      !     取得したいデータに対応する登録id 配列を入力し、隣接領域で保存されているデータ配列を取得
      call ggtools_neighbor_local_bucket_get_elements(n_neighbor, neighbor_pids, &
                       & id_array_index, id_array_item, nid, id_array, masters_out, MPI_COMM_WORLD)

    end subroutine

  end module mod_ggtools_global_bucket_test
  