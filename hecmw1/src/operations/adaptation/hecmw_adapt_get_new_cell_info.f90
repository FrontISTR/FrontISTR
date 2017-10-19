!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

!C
!C***
!C*** hecmw_adapt_GET_NEW_CELL_INFO
!C***
!C
!C    get new CELL info.
!C
subroutine hecmw_adapt_GET_NEW_CELL_INFO (hecMESH)

  use hecmw_util
  implicit real*8 (A-H,O-Z)

  integer(kind=kint), dimension(:), allocatable :: IW1
  type (hecmwST_local_mesh) :: hecMESH

  allocate (IW1(hecMESH%n_elem))
  do neib= 1, hecMESH%n_neighbor_pe
    icou= 0
    neib0= hecMESH%neighbor_pe(neib)
    do icel= 1, hecMESH%n_elem
      ic1= hecMESH%elem_ID(2*icel)
      if (ic1.eq.neib0 .and. hecMESH%adapt_type(icel).ne.0) then
        is= hecMESH%adapt_children_index(icel-1) + 1
        iE= hecMESH%adapt_children_index(icel)
        icS= hecMESH%adapt_children_local(is)
        if (hecMESH%when_i_was_refined_elem(icS).eq.                &
            &          hecMESH%n_adapt) then
          do k= is, iE
            if (hecMESH%adapt_children_item(2*k-1).ne.0) then
              icou = icou + 1
              IW1(icou)= hecMESH%adapt_children_item(2*k-1)
            endif
          enddo
        endif
      endif
    enddo

    icou= 0
    do icel= 1, hecMESH%n_elem
      ic1= hecMESH%elem_ID(2*icel)
      ic2= hecMESH%when_i_was_refined_elem(icel)
      if (ic1.eq.neib0 .and. ic2.eq.hecMESH%n_adapt) then
        icou= icou + 1
        hecMESH%elem_ID(2*icel-1)= IW1(icou)
      endif
    enddo
  enddo

  deallocate (IW1)

end subroutine hecmw_adapt_GET_NEW_CELL_INFO


