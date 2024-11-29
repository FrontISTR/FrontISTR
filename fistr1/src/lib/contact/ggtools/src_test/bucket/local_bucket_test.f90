!> バケットセルテストモジュール
module mod_ggtools_local_bucket_test
  use mod_monolis_utils
  use mod_ggtools
  implicit none

contains

  subroutine ggtools_local_bucket_test()
    implicit none

    !# バケット検索用テスト
    !call ggtools_bucket_initialize_test()
    !call ggtools_bucket_search_test()

    !# 一時的なテスト
    !call monolis_utils_aabb_test()
    !call monolis_utils_kdtree_test()
  end subroutine ggtools_local_bucket_test

end module mod_ggtools_local_bucket_test
