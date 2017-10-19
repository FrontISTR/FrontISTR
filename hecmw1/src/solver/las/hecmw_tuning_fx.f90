!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_tuning_fx
  use hecmw_util
  implicit none

  private

  public :: hecmw_tuning_fx_calc_sector_cache

  !!
  !! Please set TotalSectorCacheSize to
  !!  (on K-computer) : 12
  !!  (on FX10)       : 24
  !!
  integer, parameter :: TotalSectorCacheSize = 12

contains

  subroutine hecmw_tuning_fx_calc_sector_cache( N, NDOF, &
      sectorCacheSize0, sectorCacheSize1 )
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(out) :: sectorCacheSize0, sectorCacheSize1
    ! calculate sector cache size
    sectorCacheSize1 = int((dble(N) * NDOF * kreal / (4096 * 128)) + 0.999)
    if (sectorCacheSize1 > TotalSectorCacheSize / 2 ) &
      sectorCacheSize1 = TotalSectorCacheSize / 2
    sectorCacheSize0 = TotalSectorCacheSize - sectorCacheSize1
    ! write(*,*) 'Vector size =', N * NDOF * kreal, '[byte]  ', &
      !            'sectorCache0 =', sectorCacheSize0, '[way]  ', &
      !            'sectorCache1 =', sectorCacheSize1, '[way]'
  end subroutine hecmw_tuning_fx_calc_sector_cache

end module hecmw_tuning_fx
