!> ローカルバケットセルモジュール
module mod_ggtools_global_bucket
  use mod_monolis_utils
  use mod_ggtools_def_global_bucket
  implicit none

contains

  !> @ingroup bucket
  !> ローカルバケットのx_min、x_max を⽤いてグローバルバケット構造体を初期化する関数
  subroutine ggtools_global_bucket_init_by_local_bucket(gbucket, x_min_local, x_max_local, dx_in, comm)
    implicit none
    !> グローバルバケット構造体
    type(type_ggtools_global_bucket) :: gbucket
    !> ローカルバケットのバウンダリボックスの最小座標、実行後はグローバルバケットの格⼦点座標に調整される
    real(kdouble) :: x_min_local(3)
    !> ローカルバケットのバウンダリボックスの最大座標、実行後はグローバルバケットの格⼦点座標に調整される
    real(kdouble) :: x_max_local(3)
    !> グローバルバケットの目標セルサイズ
    real(kdouble) :: dx_in(3)
    !> MPI コミュニケータ
    integer(kint), intent(in) :: comm

    integer(kint) :: n_local, my_rank, iproc, idof
    real(kdouble), allocatable :: x_min_locals(:,:), x_max_locals(:,:)
    real(kdouble) :: x_min_global(3), x_max_global(3), dx(3)
    integer(kint) :: nx(3), idx

    !# 2-1. MPI 関数を⽤いてn_local を取得（n_local = 並列プロセス数）
    n_local = monolis_mpi_get_global_comm_size()

    !# 2-2. MPI 関数を⽤いてx_min_local[3, n_local]、x_max_local[3, n_local] を取得
    call monolis_alloc_R_2d(x_min_locals,3,n_local)
    call monolis_alloc_R_2d(x_max_locals,3,n_local)
    my_rank = monolis_mpi_get_global_my_rank()
    x_min_locals(1:3,my_rank+1) = x_min_local(1:3)
    x_max_locals(1:3,my_rank+1) = x_max_local(1:3)
    call monolis_allreduce_R(n_local, x_min_locals(1,1:n_local), monolis_mpi_sum, comm)
    call monolis_allreduce_R(n_local, x_min_locals(2,1:n_local), monolis_mpi_sum, comm)
    call monolis_allreduce_R(n_local, x_min_locals(3,1:n_local), monolis_mpi_sum, comm)
    call monolis_allreduce_R(n_local, x_max_locals(1,1:n_local), monolis_mpi_sum, comm)
    call monolis_allreduce_R(n_local, x_max_locals(2,1:n_local), monolis_mpi_sum, comm)
    call monolis_allreduce_R(n_local, x_max_locals(3,1:n_local), monolis_mpi_sum, comm)

    !# 2-3. 手順2-2 の最⼤値・最⼩値を計算しglobal%x_min 配列、global%x_max 配列、global%nx 配列、
    !#      global%dx 配列の値を計算
    !## 最大最小値の計算
    x_min_global(1:3) = x_min_locals(1:3,1)
    x_max_global(1:3) = x_max_locals(1:3,1)
    do iproc=2,n_local
      do idof=1,3
        if( x_min_locals(idof,iproc) < x_min_global(idof) ) x_min_global(idof) = x_min_locals(idof,iproc)
        if( x_max_locals(idof,iproc) > x_max_global(idof) ) x_max_global(idof) = x_max_locals(idof,iproc)
      enddo
    enddo
    !## nxとdxの辻褄が合うように入力されたdxを調整する
    do idof=1,3
      nx(idof) = max(ceiling((x_max_global(idof) - x_min_global(idof)) / dx_in(idof)),1)
      dx(idof) = (x_max_global(idof) - x_min_global(idof)) / dble(nx(idof))
    enddo

    !# 2-4. 手順2-2 で得られたglobal%x_min_local[3, n_local]、global%x_max_local[3, n_local] を
    !#      グローバルバケットの格⼦点座標に調整
    do iproc=1,n_local
      do idof=1,3
        idx = floor((x_min_locals(idof,iproc)-x_min_global(idof))/dx(idof))
        x_min_locals(idof,iproc) = x_min_global(idof) + dble(idx)*dx(idof)
        idx = ceiling((x_max_locals(idof,iproc)-x_min_global(idof))/dx(idof))
        x_max_locals(idof,iproc) = x_min_global(idof) + dble(idx)*dx(idof)
      enddo
    enddo

    !グローバルバケットの初期化
    call ggtools_global_bucket_init(gbucket, x_min_global, x_max_global, dx, x_min_locals, x_max_locals )

  end subroutine


end module mod_ggtools_global_bucket
  