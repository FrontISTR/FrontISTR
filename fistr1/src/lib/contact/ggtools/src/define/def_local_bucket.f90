!> ローカルバケットセルモジュール
module mod_ggtools_def_local_bucket
  use mod_monolis_utils
  use mod_ggtools_def_global_bucket
  implicit none

  !> @ingroup bucket
  !> 第 i 番目のバケットセルの基本情報構造体
  type type_ggtools_local_bucket_cells
    !> バケットセルに登録された個数
    integer(kint) :: nid = 0
    !> バケットセルに登録された nid 個の 1 次元整数配列
    integer(kint), allocatable :: id_array(:)
  end type type_ggtools_local_bucket_cells

  !> @ingroup bucket
  !> バケットの基本情報構造体
  type type_ggtools_local_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    !> バケットセルひとつのサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
  end type type_ggtools_local_bucket

contains

  !> @ingroup bucket
  !> バケットの初期化処理（バケットセルサイズを入力）
  !> @details ローカルバケットのx_min、x_max、グローバルバケット構造体を⽤いて
  !> @details ローカルバケット構造体を初期化する
  subroutine ggtools_local_bucket_init_by_global_bucket(ggtools_bucket_global, ggtools_bucket_local)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_bucket_global
    !> バケット構造体
    type(type_ggtools_local_bucket)  :: ggtools_bucket_local

    integer(kint) :: my_rank, nx(3)
    real(kdouble) :: dx(3)

    my_rank = monolis_mpi_get_global_my_rank()

    ggtools_bucket_local%xmin(1:3) = ggtools_bucket_global%x_min_local(1:3,my_rank+1)
    ggtools_bucket_local%xmax(1:3) = ggtools_bucket_global%x_max_local(1:3,my_rank+1)

    dx(1:3) = ggtools_bucket_global%dx(1:3)
    ggtools_bucket_local%dx = dx(1:3)

    nx(1) = ceiling( (ggtools_bucket_local%xmax(1) - ggtools_bucket_local%xmin(1))/dx(1) -ths )
    nx(2) = ceiling( (ggtools_bucket_local%xmax(2) - ggtools_bucket_local%xmin(2))/dx(2) -ths )
    nx(3) = ceiling( (ggtools_bucket_local%xmax(3) - ggtools_bucket_local%xmin(3))/dx(3) -ths )
    ggtools_bucket_local%nx(1:3) = nx(1:3)

    if(nx(1) < 1 .or. nx(2) < 1 .or. nx(3) < 1)then
      call monolis_std_error_string("the number of local bucket cells is less than 1")
      call monolis_std_error_stop()
    endif

  end subroutine ggtools_local_bucket_init_by_global_bucket

  !> @ingroup bucket
  !> バケットの終了処理
  subroutine ggtools_local_bucket_finalize(ggtools_bucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket

    ggtools_bucket%xmin = 0.0d0
    ggtools_bucket%xmax = 0.0d0
    ggtools_bucket%dx = 0.0d0
    ggtools_bucket%nx = 0
  end subroutine ggtools_local_bucket_finalize

  !> @ingroup bucket
  !> バケットのプリント
  subroutine ggtools_local_bucket_print(ggtools_bucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket

    write(*,*) "print local bucket"
    write(*,*) "xmin"
    write(*,*) ggtools_bucket%xmin
    write(*,*) "xmax"
    write(*,*) ggtools_bucket%xmax
    write(*,*) "nx"
    write(*,*) ggtools_bucket%nx
    write(*,*) "dx"
    write(*,*) ggtools_bucket%dx

  end subroutine ggtools_local_bucket_print

  !> @ingroup bucket
  !> バケットセルの初期化処理
  subroutine ggtools_local_bucket_cells_init(ggtools_bucket_cell, nx)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells), allocatable :: ggtools_bucket_cell(:)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    integer(kint) :: n_total

    n_total = nx(1)*nx(2)*nx(3)
    allocate(ggtools_bucket_cell(n_total))
  end subroutine ggtools_local_bucket_cells_init

  !> @ingroup bucket
  !> バケットセルの終了処理
  subroutine ggtools_local_bucket_cells_finalize(ggtools_bucket_cell, nx)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    integer(kint) :: i, n_total

    n_total = nx(1)*nx(2)*nx(3)

    do i = 1, n_total
      if(ggtools_bucket_cell(i)%nid == 0) cycle
      deallocate(ggtools_bucket_cell(i)%id_array)
    enddo
  end subroutine ggtools_local_bucket_cells_finalize

  !> @ingroup bucket
  !> バケットサイズを取得
  subroutine ggtools_local_bucket_get_bucket_size(ggtools_bucket, xmin, xmax)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: xmin(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: xmax(3)
    xmin = ggtools_bucket%xmin
    xmax = ggtools_bucket%xmax
  end subroutine ggtools_local_bucket_get_bucket_size

  !> @ingroup bucket
  !> バケットセル分割数を取得
  subroutine ggtools_local_bucket_get_number_of_bucket_divisions(ggtools_bucket, nx)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    nx = ggtools_bucket%nx
  end subroutine ggtools_local_bucket_get_number_of_bucket_divisions

  !> @ingroup bucket
  !> バケットセルサイズを取得
  subroutine ggtools_local_bucket_get_size_of_one_bucket_cell(ggtools_bucket, dx)
    implicit none
    !> バケット構造体
    type(type_ggtools_local_bucket) :: ggtools_bucket
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    dx = ggtools_bucket%dx
  end subroutine ggtools_local_bucket_get_size_of_one_bucket_cell

  !> @ingroup bucket
  !> バケットセルのプリント
  subroutine ggtools_local_bucket_cell_print(ggtools_bucket_cell)
    implicit none
    !> バケットセル構造体
    type(type_ggtools_local_bucket_cells) :: ggtools_bucket_cell(:)

    integer(kint) :: i, n_bucket_cell

    write(*,*) "print local bucket cell"
    n_bucket_cell = size(ggtools_bucket_cell(:))
    do i=1, n_bucket_cell
      write(*,*) "bucket cell index, nid:", i, ggtools_bucket_cell(i)%nid
      write(*,*) "id_array:"
      write(*,'(10(I9,","))') ggtools_bucket_cell(i)%id_array
    enddo
  end subroutine ggtools_local_bucket_cell_print

end module mod_ggtools_def_local_bucket
