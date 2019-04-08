!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

!C
!C***
!C*** hecmw_adapt_REPRO_COMM_TABLE
!C***
!C
!C    reproduce COMMUNICATION table
!C
subroutine hecmw_adapt_REPRO_COMM_TABLE (hecMESH)

  use  hecmw_util
  use  hecmw_adapt_STACK_SR
  use  hecmw_adapt_ITEM_SR
  use  hecmw_adapt_INT_SR

  integer(kind=kint ), pointer     :: wSI(:), wSE(:), wiI(:), wiE(:)
  integer(kind=kint ), allocatable :: IW1(:), IW2(:)

  type (hecmwST_local_mesh) :: hecMESH

  !C
  !C-- init.
  allocate (IW1(hecMESH%n_neighbor_pe))
  IW1 = 0

  allocate (wSE(0:hecMESH%n_neighbor_pe), wSI(0:hecMESH%n_neighbor_pe))
  wSE = 0
  wSI = 0

  !C
  !C-- search IMPORT items
  do in= 1, hecMESH%n_adapt_node_cur
    is= hecMESH%adapt_NEWtoOLD_node(in)
    ih= hecMESH%node_ID(2*is)
    if (ih.ne.hecMESH%my_rank) then
      ihR = hecMESH%rev_neighbor_pe(ih)
      wSI(ihR)= wSI(ihR) + 1
    endif
  enddo

  !C
  !C-- exchange INFO. on IMPORT item #
  call hecmw_adapt_STACK_SEND_RECV                                  &
    &   ( hecMESH%n_neighbor_pe, hecMESH%neighbor_pe, wSI, wSE,        &
    &     hecMESH%MPI_COMM, hecMESH%my_rank)

  !C
  !C-- IMPORT/EXPORT item #
  do neib= 1, hecMESH%n_neighbor_pe
    wSI(neib)= wSI(neib-1) + wSI(neib)
    wSE(neib)= wSE(neib-1) + wSE(neib)
  enddo

  !C
  !C-- send as IMPORT/recv. as EXPORT

  allocate (wiI(wSI(hecMESH%n_neighbor_pe)))
  allocate (wiE(wSE(hecMESH%n_neighbor_pe)))
  wiI = 0
  wiE = 0

  do in= 1, hecMESH%n_adapt_node_cur
    is= hecMESH%adapt_NEWtoOLD_node(in)
    ih= hecMESH%node_ID(2*is)
    if (ih.ne.hecMESH%my_rank) then
      ihR             = hecMESH%rev_neighbor_pe(ih)
      IW1(ihR  )          = IW1(ihR) + 1
      wiI(wSI(ihR-1)+IW1(ihR))= hecMESH%node_ID(2*is-1)
    endif
  enddo
  deallocate (IW1)

  N= hecMESH%n_adapt_node_cur
  len= max(wSI(hecMESH%n_neighbor_pe),wSE(hecMESH%n_neighbor_pe),N)
  allocate (IW1(len), IW2(len))
  IW1 = 0
  IW2 = 0

  call hecmw_adapt_ITEM_SEND_RECV                                   &
    &    (len, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
    &     wSI, wiI, wSE, wiE, IW1, IW2,                                &
    &     hecMESH%MPI_COMM, hecMESH%my_rank, 1)

  !C
  !C-- reconstruct IMPORT table

  deallocate (IW1)
  allocate (IW1(hecMESH%n_neighbor_pe))
  IW1 = 0

  do in= 1, hecMESH%n_adapt_node_cur
    is= hecMESH%adapt_NEWtoOLD_node(in)
    ih= hecMESH%node_ID(2*is)
    if (ih.ne.hecMESH%my_rank) then
      ihR             = hecMESH%rev_neighbor_pe(ih)
      IW1(ihR  )          = IW1(ihR) + 1
      wiI(wSI(ihR-1)+IW1(ihR))= in
    endif
  enddo

  !C
  !C-- new ARRAY
  allocate (hecMESH%adapt_import_new_index(0:hecMESH%n_neighbor_pe))
  allocate (hecMESH%adapt_export_new_index(0:hecMESH%n_neighbor_pe))

  do neib= 0, hecMESH%n_neighbor_pe
    hecMESH%adapt_import_new_index(neib)= wSI(neib)
    hecMESH%adapt_export_new_index(neib)= wSE(neib)
  enddo

  MAXimport= wSI(hecMESH%n_neighbor_pe)
  MAXexport= wSE(hecMESH%n_neighbor_pe)

  allocate (hecMESH%adapt_import_new_item(MAXimport))
  allocate (hecMESH%adapt_export_new_item(MAXexport))

  do k= 1, MAXimport
    hecMESH%adapt_import_new_item(k)= wiI(k)
  enddo
  do k= 1, MAXexport
    hecMESH%adapt_export_new_item(k)= wiE(k)
  enddo

  deallocate (IW1,IW2)
  deallocate (wSE,wSI,wiE,wiI)
  !C===

end subroutine hecmw_adapt_REPRO_COMM_TABLE
