      module JAD_TYPE
      use hecmw_util
      implicit none
!C---------------------- AU&AL
      real(kind=kreal), allocatable      :: AJAD(:)
      integer(kind=kint), allocatable    :: JAJAD(:)
      integer(kind=kint), allocatable    :: JADORD(:)
      integer(kind=kint), allocatable    :: IAJAD(:)
      integer(kind=kint) :: MJAD
      real(kind=kreal), allocatable  :: WP1(:), WP2(:), WP3(:)

      logical:: JAD_ON_F=.FALSE.
      end module JAD_TYPE
