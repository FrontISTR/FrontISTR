!> バケットセルモジュール
module mod_ggtools_def_global_bucket
  use mod_monolis_utils
  implicit none

  real(kdouble) :: ths = 1.0d-6

  !> @ingroup bucket
  !> グローバルバケットの基本情報構造体
  type type_ggtools_global_bucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: x_min(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: x_max(3)
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    !> バケットセルひとつのサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    !> バウンダリボックスの最小座標：x_min_local[3, n_local]
    real(kdouble), allocatable :: x_min_local(:,:)
    !> バウンダリボックスの最大座標：x_max_local[3, n_local]
    real(kdouble), allocatable :: x_max_local(:,:)
  end type type_ggtools_global_bucket

contains

  !> @ingroup bucket
  !> グローバルバケットの初期化処理（バケットセルサイズを入力）
  !> @details バウンダリボックスの最小座標 x_min を起点に、バウンディングボックスを覆う最小のバケットセルを確保する。
  !> @details バケットセル確保後に、バウンダリボックスの最大座標 x_max が更新される。
  subroutine ggtools_global_bucket_init(ggtools_gbucket, x_min, x_max, dx, x_min_locals, x_max_locals )
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: x_min(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: x_max(3)
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    !> コミュニケータのproc数
    integer(kint) :: n_local
    real(kdouble) :: x_min_locals(:,:)
    real(kdouble) :: x_max_locals(:,:)

    !# 2-1. MPI 関数を⽤いてn_local を取得（n_local = 並列プロセス数）
    n_local = monolis_mpi_get_global_comm_size() 
       
    ggtools_gbucket%x_min = x_min
    ggtools_gbucket%x_max = x_max

    ggtools_gbucket%dx = dx

    ggtools_gbucket%nx(1) = 1
    ggtools_gbucket%nx(2) = 2
    ggtools_gbucket%nx(3) = 3

    call monolis_alloc_R_2d(ggtools_gbucket%x_min_local, 3, n_local)
    call monolis_alloc_R_2d(ggtools_gbucket%x_max_local, 3, n_local)

    ggtools_gbucket%x_min_local(:,:) = x_min_locals(:,:)
    ggtools_gbucket%x_max_local(:,:) = x_max_locals(:,:)
  end subroutine ggtools_global_bucket_init

  !> @ingroup bucket
  !> バケットの終了処理
  subroutine ggtools_global_bucket_finalize(ggtools_gbucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket

    ggtools_gbucket%x_min = 0.0d0
    ggtools_gbucket%x_max = 0.0d0
    ggtools_gbucket%dx = 0.0d0
    ggtools_gbucket%nx = 0
    deallocate(ggtools_gbucket%x_min_local)
    deallocate(ggtools_gbucket%x_max_local)
  end subroutine ggtools_global_bucket_finalize

  !> @ingroup bucket
  !> バケットのプリント
  subroutine ggtools_global_bucket_print(ggtools_gbucket)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket

    write(*,*) "ggtools_gbucket%x_min"
    write(*,*) ggtools_gbucket%x_min
    write(*,*) "ggtools_gbucket%x_max"
    write(*,*) ggtools_gbucket%x_max
    write(*,*) "ggtools_gbucket%dx"
    write(*,*) ggtools_gbucket%dx
    write(*,*) "ggtools_gbucket%nx"
    write(*,*) ggtools_gbucket%nx
    write(*,*) "ggtools_gbucket%x_min_local"
    write(*,*) ggtools_gbucket%x_min_local
    write(*,*) "ggtools_gbucket%x_max_local"
    write(*,*) ggtools_gbucket%x_max_local

  end subroutine ggtools_global_bucket_print

  !> @ingroup bucket
  !> バケットサイズを取得
  subroutine ggtools_global_bucket_get_bucket_size(ggtools_gbucket, x_min, x_max)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket
    !> バウンダリボックスの最小座標
    real(kdouble) :: x_min(3)
    !> バウンダリボックスの最大座標
    real(kdouble) :: x_max(3)
    x_min = ggtools_gbucket%x_min
    x_max = ggtools_gbucket%x_max
  end subroutine ggtools_global_bucket_get_bucket_size

  !> @ingroup bucket
  !> バケットセル分割数を取得
  subroutine ggtools_global_bucket_get_number_of_bucket_divisions(ggtools_gbucket, nx)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket
    !> バケットセル分割数（nx, ny, nz）
    integer(kint) :: nx(3)
    nx = ggtools_gbucket%nx
  end subroutine ggtools_global_bucket_get_number_of_bucket_divisions

  !> @ingroup bucket
  !> バケットセルサイズを取得
  subroutine ggtools_global_bucket_get_size_of_one_bucket_cell(ggtools_gbucket, dx)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket
    !> バケットセルサイズ（dx, dy, dz）
    real(kdouble) :: dx(3)
    dx = ggtools_gbucket%dx
  end subroutine ggtools_global_bucket_get_size_of_one_bucket_cell

  !> @ingroup bucket
  !> ローカルバウンディングボックスを取得
  subroutine ggtools_global_bucket_get_local_bb(ggtools_gbucket, pid, xmin_local, xmax_local)
    implicit none
    !> バケット構造体
    type(type_ggtools_global_bucket) :: ggtools_gbucket
    !> バウンディングボックスのpid
    integer(kint) :: pid
    !> ローカルバウンディングボックスの最小座標
    real(kdouble) :: xmin_local(3)
    !> ローカルバウンディングボックスの最大座標
    real(kdouble) :: xmax_local(3)
    xmin_local(1:3) = ggtools_gbucket%x_min_local(1:3,pid+1)
    xmax_local(1:3) = ggtools_gbucket%x_max_local(1:3,pid+1)
  end subroutine ggtools_global_bucket_get_local_bb


end module mod_ggtools_def_global_bucket
