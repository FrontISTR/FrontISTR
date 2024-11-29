!> エレメントモジュール
module mod_ggtools_def_element
    use mod_monolis_utils
    implicit none
  
    !> @ingroup element
    !> エレメント構造体
    type type_ggtools_element
      !> エレメントを構成する計算点数
      integer(kint) :: n_point
      !> 計算点の座標
      real(kdouble), allocatable :: coordinate(:,:)
      !> 領域番号（MPI ランク番号）
      integer(kint) :: domain_id
      !> 領域番号domain_id におけるローカルエレメント番号
      integer(kint) :: element_id   
    end type type_ggtools_element
  
    contains

    !> @ingroup element
    !> エレメント構造体の初期化処理
    subroutine ggtools_element_init(ggtools_element, n_point, coordinate, domain_id, element_id)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element) :: ggtools_element
      !> エレメントを構成する計算点数
      integer(kint) :: n_point
      !> 計算点の座標
      real(kdouble), allocatable :: coordinate(:,:)
      !> 領域番号（MPI ランク番号）
      integer(kint) :: domain_id
      !> 領域番号domain_id におけるローカルエレメント番号
      integer(kint) :: element_id

      ggtools_element%n_point = n_point
      call monolis_alloc_R_2d(ggtools_element%coordinate, 3, n_point)
      ggtools_element%coordinate(1:3,1:n_point) = coordinate(1:3,1:n_point)
      ggtools_element%domain_id = domain_id
      ggtools_element%element_id = element_id
    end subroutine ggtools_element_init

    !> @ingroup element
    !> エレメント構造体のコピー
    subroutine ggtools_element_copy(ggtools_element1, ggtools_element2)
      implicit none
      !> エレメント構造体（コピー元）
      type(type_ggtools_element) :: ggtools_element1
      !> エレメント構造体（コピー先）
      type(type_ggtools_element) :: ggtools_element2

      call ggtools_element_init(ggtools_element2, ggtools_element1%n_point, &
        & ggtools_element1%coordinate, ggtools_element1%domain_id, ggtools_element1%element_id)
    end subroutine ggtools_element_copy
    
    !> @ingroup element
    !> エレメント構造体の初期化処理
    subroutine ggtools_element_finalize(ggtools_element)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element) :: ggtools_element

      ggtools_element%n_point = 0
      call monolis_dealloc_R_2d(ggtools_element%coordinate)
      ggtools_element%domain_id = 0
      ggtools_element%element_id = 0
    end subroutine ggtools_element_finalize

    !> @ingroup element
    !> エレメント構造体のプリント
    subroutine ggtools_element_print(ggtools_element)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element) :: ggtools_element

      integer(kint) :: i

      write(*,*) "ggtools_element%n_point",ggtools_element%n_point
      write(*,*) "ggtools_element%coordinate"
      do i=1,ggtools_element%n_point
        write(*,*) ggtools_element%coordinate(1:3,i)
      enddo
      write(*,*) "ggtools_element%domain_id",ggtools_element%domain_id
      write(*,*) "ggtools_element%element_id",ggtools_element%element_id
    end subroutine ggtools_element_print

    !> @ingroup element
    !> エレメント構造体の構成節点数取得
    function ggtools_element_get_num_points(ggtools_element)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element) :: ggtools_element
      !> エレメントを構成する計算点数
      integer(kint) :: ggtools_element_get_num_points

      ggtools_element_get_num_points = ggtools_element%n_point
    end function ggtools_element_get_num_points

    !> @ingroup element
    !> エレメント構造体の領域番号取得
    function ggtools_element_get_domain_id(ggtools_element) result(domain_id)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element), intent(in) :: ggtools_element
      !> 領域番号
      integer :: domain_id
  
      domain_id = ggtools_element%domain_id
    end function ggtools_element_get_domain_id
  
    !> @ingroup element
    !> エレメント構造体のローカルエレメント番号取得
    function ggtools_element_get_element_id(ggtools_element) result(element_id)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element), intent(in) :: ggtools_element
      !> ローカルエレメント番号
      integer :: element_id
  
      element_id = ggtools_element%element_id
    end function ggtools_element_get_element_id

    !> @ingroup element
    !> エレメント構造体のLBB取得
    subroutine ggtools_element_get_element_lbb(ggtools_element,xmin_local,xmax_local)
      implicit none
      !> エレメント構造体
      type(type_ggtools_element), intent(in) :: ggtools_element
      !> エレメントLBB
      real(kdouble)  :: xmin_local(3), xmax_local(3)
      real(kdouble)  :: tmpcoord(3)

      integer(kint) :: i, idof, n_point

      n_point = ggtools_element_get_num_points(ggtools_element)

      xmin_local = ggtools_element%coordinate(1:3,1)
      xmax_local = xmin_local

      do i=2,n_point
        tmpcoord(1:3) = ggtools_element%coordinate(1:3,i)
        do idof = 1,3
          if( tmpcoord(idof) < xmin_local(idof) ) xmin_local(idof) = tmpcoord(idof)
          if( tmpcoord(idof) > xmax_local(idof) ) xmax_local(idof) = tmpcoord(idof)
        enddo
      enddo

    end subroutine ggtools_element_get_element_lbb

end module mod_ggtools_def_element