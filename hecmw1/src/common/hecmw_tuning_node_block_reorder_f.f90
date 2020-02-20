!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module hecmw_tuning_node_block_reorder
  use hecmw_util
  implicit none


  private ! default

  public :: hecmw_tuning_reorder_init
  public :: hecmw_tuning_reorder_do


  !C
  !C***
  !C*** TYPE definition
  !C***
  !C

  type tVector
    integer(kind=kint) :: n
    integer(kind=kint) :: maxn
    logical            :: sorted
    integer(kind=kint), allocatable :: vec(:)
  end type


contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !C
  !C***
  !C*** SUBROUTINES
  !C***
  !C
  !!!!!!! --- add connectivity_reorder start --- !!!!!!!
  subroutine tVecInit(tvec,init_maxn)
    type(tVector), intent(inout)   :: tvec
    integer(kind=kint),intent(in)  :: init_maxn

    tvec%n = 0
    tvec%maxn = init_maxn
    allocate(tvec%vec(init_maxn))
    tvec%vec(1:init_maxn) = 0
    tvec%sorted = .true.
  end subroutine

  subroutine tVecPrint(tvec)
    type(tVector), intent(inout)   :: tvec

    write(*,*) 'n,maxn',tvec%n,tvec%maxn
    write(*,*) tvec%vec(1:tvec%n)
  end subroutine

  subroutine tVecAdd(tvec,item)
    type(tVector), intent(inout)   :: tvec
    integer(kind=kint),intent(in)  :: item

    integer(kind=kint) :: maxn_old
    integer(kind=kint), allocatable :: tmp(:)

    if( tvec%n == tvec%maxn ) then
      maxn_old = tvec%maxn
      allocate(tmp(maxn_old))
      tmp(1:maxn_old) = tvec%vec(1:maxn_old)
      deallocate(tvec%vec)

      if( maxn_old < 1000 ) then
        tvec%maxn = min(maxn_old*2,1000)
      else
        tvec%maxn = maxn_old+1000
      end if

      allocate(tvec%vec(tvec%maxn))
      tvec%vec(1:maxn_old) = tmp(1:maxn_old)
      tvec%vec(maxn_old+1:tvec%maxn) = 0
      deallocate(tmp)
    end if

    tvec%n = tvec%n + 1
    tvec%vec(tvec%n) = item

  end subroutine

  recursive subroutine sort_int_array(array, istart, iend)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint) :: left, right, center
    integer(kind=kint) :: pivot, tmp
    if (istart >= iend) return
    center = (istart + iend) / 2
    pivot = array(center)
    left = istart
    right = iend
    do
      do while (array(left) < pivot)
        left = left + 1
      end do
      do while (pivot < array(right))
        right = right - 1
      end do
      if (left >= right) exit
      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp
      left = left + 1
      right = right - 1
    end do
    if (istart < left-1) call sort_int_array(array, istart, left-1)
    if (right+1 < iend) call sort_int_array(array, right+1, iend)
  end subroutine sort_int_array

  recursive subroutine sort_int_array2(array, istart, iend)
    implicit none
    integer(kind=kint), intent(inout) :: array(:,:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint) :: left, right, center
    integer(kind=kint) :: pivot, tmp(2)
    if (istart >= iend) return
    center = (istart + iend) / 2
    pivot = array(2,center)
    left = istart
    right = iend
    do
      do while (array(2,left) < pivot)
        left = left + 1
      end do
      do while (pivot < array(2,right))
        right = right - 1
      end do
      if (left >= right) exit
      tmp(1:2) = array(1:2,left)
      array(1:2,left) = array(1:2,right)
      array(1:2,right) = tmp
      left = left + 1
      right = right - 1
    end do
    if (istart < left-1) call sort_int_array2(array, istart, left-1)
    if (right+1 < iend) call sort_int_array2(array, right+1, iend)
  end subroutine sort_int_array2

  subroutine tVecCompress(tvec)
    type(tVector), intent(inout)   :: tvec

    integer(kind=kint) :: i, idx

    call sort_int_array(tvec%vec, 1, tvec%n)
    idx = 1
    do i=1,tvec%n-1
      if( tvec%vec(i+1) == tvec%vec(i) ) cycle
      idx = idx + 1
      tvec%vec(idx) = tvec%vec(i+1)
    end do
    tvec%n = idx
  end subroutine

  subroutine tVecFinilize(tvec)
    type(tVector), intent(inout)   :: tvec
    deallocate(tvec%vec)
  end subroutine
  !!!!!!! --- add connectivity_reorder end --- !!!!!!!

  !C
  !C    INIT variables.
  !C
  subroutine hecmw_tuning_reorder_init(hecMESH, ierr)

    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint) :: ierr

    integer :: i

    allocate(hecMESH%tuning_block_reorder_old2new(hecMESH%n_node))
    allocate(hecMESH%tuning_block_reorder_new2old(hecMESH%n_node))

    do i=1, hecMESH%n_node
      hecMESH%tuning_block_reorder_new2old(i)=i
      hecMESH%tuning_block_reorder_old2new(i)=i
    end do

    hecMESH%tuning_block_reorder_on   = .FALSE.

    ierr=0
  end subroutine hecmw_tuning_reorder_init


  subroutine hecmw_tuning_reorder_do(hecMESH,    &
                                     block_numx, &
                                     block_numy, &
                                     block_numz, &
                                     block_inout_ratio)

    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint), intent(in), optional :: block_numx
    integer(kind=kint), intent(in), optional :: block_numy
    integer(kind=kint), intent(in), optional :: block_numz
    real(kind=kreal),   intent(in), optional :: block_inout_ratio

    integer(kind=kint) :: numx, numy, numz
    real(kind=kreal)   :: ratio

    hecMESH%tuning_block_reorder_on = .TRUE.

    ! default block reordering parameter
    numx  = 10
    numy  = 10
    numz  = 10
    ratio = 0.1

    if (present(block_numx)) then
      if (block_numx <= 1) then
        numx = 1
      else
        numx = block_numx
      end if
    end if

    if (present(block_numy)) then
      if (block_numy <= 1) then
        numy = 1
      else
        numy = block_numy
      end if
    end if

    if (present(block_numz)) then
      if (block_numz <= 1) then
        numz = 1
      else
        numz = block_numz
      end if
    end if

    if (present(block_inout_ratio)) then
      if (block_inout_ratio <= 0.001) then
        ratio = 0.001
      else if (block_inout_ratio >= 0.999) then
        ratio = 0.999
      else
        ratio = block_inout_ratio
      end if
    end if

    call make_reorder_table(hecMESH, numx, numy, numz, ratio)

    call reorder_local_node_ID(hecMESH)

  end subroutine hecmw_tuning_reorder_do


  subroutine reorder_local_node_ID(hecMESH)

    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer :: i, ndof, new_id, nnode

    real(kind=kreal), pointer   :: new_node(:)
    real(kind=kreal), pointer   :: old_node(:)

    integer(kind=kint), pointer :: new_global_node_ID(:)
    integer(kind=kint), pointer :: old_global_node_ID(:)

    integer(kind=kint), pointer :: new_elem_node_item(:)

    integer(kind=kint), pointer :: new_node_group_grp_item(:)

    integer(kind=kint), pointer :: old2new(:)

    nnode = hecMESH%n_node
    old2new => hecMESH%tuning_block_reorder_old2new

    ! reorder node coordinate
    allocate(new_node(nnode*3))
    do i=1, nnode
      new_id = old2new(i)
      new_node(new_id*3  )=hecMESH%node(i*3  )
      new_node(new_id*3-1)=hecMESH%node(i*3-1)
      new_node(new_id*3-2)=hecMESH%node(i*3-2)
    end do
    deallocate(hecMESH%node)
    hecMESH%node => new_node

    ! reorder global node ID
    allocate(new_global_node_ID(nnode))
    do i=1, nnode
      new_id = old2new(i)
      new_global_node_ID(new_id)=hecMESH%global_node_ID(i)
    end do
    deallocate(hecMESH%global_node_ID)
    hecMESH%global_node_ID => new_global_node_ID

    ! reorder local node ID on each element
    allocate(new_elem_node_item(hecMESH%elem_node_index(hecMESH%n_elem)))
    do i=1, hecMESH%elem_node_index(hecMESH%n_elem)
      new_id = old2new(hecMESH%elem_node_item(i))
      new_elem_node_item(i) = new_id
    end do
    deallocate(hecMESH%elem_node_item)
    hecMESH%elem_node_item => new_elem_node_item

    ! reorder local node ID on node group
    allocate(new_node_group_grp_item(size(hecMESH%node_group%grp_item)))
    do i=1, size(hecMESH%node_group%grp_item)
      new_id = old2new(hecMESH%node_group%grp_item(i))
      new_node_group_grp_item(i)=new_id
    end do
    deallocate(hecMESH%node_group%grp_item)
    hecMESH%node_group%grp_item => new_node_group_grp_item

    ! reorder communication table
    if (size(hecMESH%export_item) > 1) then
      do i=1, size(hecMESH%export_item)
        new_id = old2new(hecMESH%export_item(i))
        hecMESH%export_item(i) = new_id
      end do
    end if

  end subroutine reorder_local_node_ID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_reorder_table(hecMESH, numx, numy, numz, ratio)

    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint),        intent(in)    :: numx, numy, numz
    real(kind=kreal),          intent(in)    :: ratio

    integer(kind=kint) :: num
    integer(kind=kint), allocatable :: num_reorder(:)
    real(kind=kreal), allocatable   :: pos(:,:)

    integer :: i

    num = hecMESH%nn_internal ! only internal nodes are reorderd.

    allocate(pos(3,num))
    do i=1, num
      pos(1,i) = hecMESH%node(i*3-2)
      pos(2,i) = hecMESH%node(i*3-1)
      pos(3,i) = hecMESH%node(i*3  )
    end do

    allocate(num_reorder(num))

    !call block_reorder(num,numx,numy,numz,ratio,num_reorder,pos)
    !call connectivity_reorder(hecMESH,num_reorder)
    call metis_reorder(hecMESH,1200,num_reorder)

    ! Because reorder table for external nodes are already initialized,
    ! only internal nodes are updated.
    do i=1,num
      hecMESH%tuning_block_reorder_old2new(i) = num_reorder(i)
      hecMESH%tuning_block_reorder_new2old(num_reorder(i)) = i
    end do

    deallocate(num_reorder)
    deallocate(pos)

  end subroutine make_reorder_table

  subroutine metis_reorder(hecMESH,num_per_block,num_reorder)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in)        :: num_per_block
    integer(kind=kint), intent(inout)     :: num_reorder(:)

    integer              :: i, numpart
    type(tVector), allocatable :: n2nlist(:)

    !for metis
    integer, allocatable :: npart(:)

    !for numbering
    integer              :: region
    integer, allocatable :: id_inblock(:)
    integer, allocatable :: nregion(:), offset(:)


    !create node to node list
    call make_n2nlist(hecMESH,n2nlist)

    ! partition by metis
    numpart = hecMESH%n_node/num_per_block+1
    allocate(npart(hecMESH%n_node))
    call part_by_metis(hecMESH,n2nlist,numpart,npart)

    ! reorder inter-blocks
    call reorder_interblocks(n2nlist,numpart,npart)

    ! call reorder inblocks
    allocate(id_inblock(hecMESH%n_node))
    call reorder_inblocks(n2nlist,hecMESH%node,numpart,npart,id_inblock)

    ! numbering
    allocate(nregion(0:numpart),offset(0:numpart))
    nregion(:) = 0
    offset(:) = 0
    do i=1,hecMESH%nn_internal
      nregion(npart(i)) = nregion(npart(i)) + 1
    end do
    do i=1,numpart-1
      offset(i) = nregion(i-1) + offset(i-1)
    end do

    nregion(:) = 0
    do i=1,hecMESH%nn_internal
      region = npart(i)
      num_reorder(i) = offset(region)+id_inblock(i)
      nregion(region) = nregion(region) + 1
    end do


    deallocate(nregion,offset)

    deallocate(npart,id_inblock)

  end subroutine

  subroutine reorder_interblocks(n2nlist,numpart,npart)
    type(tVector), allocatable, intent(in)         :: n2nlist(:)
    integer(kind=kint), intent(in)                 :: numpart
    integer(kind=kint), allocatable, intent(inout) :: npart(:)

    integer(kind=kint) :: i,j,k,iter,nnode,con(numpart,numpart)
    integer(kind=kint) :: pt1,pt2,threshold
    type(tVector)      :: dcon(numpart),weight(numpart)
    logical :: used(numpart)
    integer(kind=kint) :: mini, minw, maxi, maxw, counter, next
    integer(kind=kint) :: neworder(numpart), reverse_order(numpart)

    nnode = size(n2nlist)

    con(:,:) = 0
    do i=1,nnode
      pt1 = npart(i)+1
      if( pt1 < 1 ) cycle
      do k=1,n2nlist(i)%n
        j = n2nlist(i)%vec(k)
        pt2 = npart(j)+1
        if( pt2 < 1 ) cycle
        con(pt1,pt2) = con(pt1,pt2) + 1
      end do
    end do

    threshold = 0
    do i=1,numpart
      do j=i+1,numpart
        threshold = threshold + 2*con(i,j)
      end do
    end do
    threshold = threshold/(numpart*numpart*10)
    !write(*,*) threshold !tmp

    do i=1,numpart
      call tVecInit(dcon(i),numpart)
      call tVecInit(weight(i),numpart)
      do j=1,numpart

        if( i == j ) cycle
        if( con(i,j) < threshold ) cycle
        call tVecAdd(dcon(i),j)
        call tVecAdd(weight(i),con(i,j))
      end do
      !write(*,*) "print",i
      !call tVecPrint(dcon(i))
      !call tVecPrint(weight(i))
    end do

    used(:) = .false.
    next = 0
    do iter=1,numpart
      if( next == 0 ) then
        mini = 0
        minw = numpart+1
        do i=1,numpart
          if( used(i) ) cycle
          counter = 0
          do j=1,dcon(i)%n
            if( used(dcon(i)%vec(j)) ) cycle
            counter = counter + 1
          end do
          if( counter < minw ) then
            mini = i
            minw = dcon(i)%n
          end if
        end do
      else
        mini = next
      end if
      neworder(iter) = mini
      used(mini) = .true.

      ! set next
      maxi = 0
      maxw = 0
      do j=1,dcon(mini)%n
        if( used(dcon(mini)%vec(j)) ) cycle
        if( weight(mini)%vec(j) > maxw ) then
          maxw = weight(mini)%vec(j)
          maxi = dcon(mini)%vec(j)
        end if
      end do
      next = maxi
    end do

    do i=1,numpart
      reverse_order(neworder(i)) = i
    end do

    do i=1,nnode
      pt1 = npart(i)+1
      if( pt1 < 1 ) cycle
      pt2 = reverse_order(pt1)-1
      npart(i) = pt2
    end do

    do i=1,numpart
      call tVecFinilize(dcon(i))
      call tVecFinilize(weight(i))
    end do

  end subroutine

  subroutine reorder_inblocks(n2nlist,nodes,numpart,npart,id_inblock)
    type(tVector), allocatable, intent(in)         :: n2nlist(:)
    real(kind=kreal), pointer, intent(in  )        :: nodes(:) !< coordinate for inter ordering
    integer(kind=kint), intent(in)                 :: numpart
    integer(kind=kint), allocatable, intent(in)    :: npart(:)
    integer(kind=kint), allocatable, intent(inout) :: id_inblock(:)


    integer(kind=kint) :: i, j, k, nid1, nid2, pt1, pt2, nnode, lid
    type(tVector) :: inblock_list

    integer(kind=kint), allocatable :: id_local(:,:)

    nnode = size(n2nlist)
    allocate(id_local(2,nnode))
    call tVecInit(inblock_list,2*nnode/numpart)

    do i=0,numpart-1
      !make node list
      inblock_list%n = 0 !clear inblock_list
      do j=1,nnode
        if( npart(j) /= i ) cycle
        call tVecAdd(inblock_list,j)
        id_local(1,inblock_list%n) = inblock_list%n
      end do

      ! label inter-block node
      do j=1,inblock_list%n
        nid1 = inblock_list%vec(j)
        pt1 = npart(nid1)
        lid = 0
        do k=1,n2nlist(nid1)%n
          nid2 = n2nlist(nid1)%vec(k)
          pt2 = npart(nid2)
          lid = lid + pt2
          if( pt2 == pt1-1 ) then
            lid = lid - numpart
          else if( pt2 == pt1+1 ) then
            lid = lid + numpart
          end if
        end do
        id_local(2,j) = 10000*lid/n2nlist(nid1)%n
      end do

      call sort_int_array2(id_local, 1, inblock_list%n)

      ! reorder
      do j=1,inblock_list%n
        nid1 = inblock_list%vec(id_local(1,j))
        id_inblock(nid1) = j
      end do

    end do
  end subroutine

  subroutine reorder_inblocks0(n2nlist,nodes,numpart,npart,id_inblock)
    type(tVector), allocatable, intent(in)         :: n2nlist(:)
    real(kind=kreal), pointer, intent(in  )        :: nodes(:) !< coordinate for inter ordering
    integer(kind=kint), intent(in)                 :: numpart
    integer(kind=kint), allocatable, intent(in)    :: npart(:)
    integer(kind=kint), allocatable, intent(inout) :: id_inblock(:)


    integer(kind=kint) :: i, j, k, nid1, nid2, pt1, pt2, nnode, lid, currnnode, currnnode0
    type(tVector) :: inblock_list, lower_list_local, upper_list_local
    type(tVector) :: inter_lower_local, inter_upper_local, stock
    type(tVector) :: inter_lower_tmp, inter_upper_tmp

    integer(kind=kint), allocatable :: internal_id_local(:), inorder(:)
    logical, allocatable :: used(:)

    nnode = size(n2nlist)

    call tVecInit(inblock_list,2*nnode/numpart)
    call tVecInit(lower_list_local,inblock_list%maxn)
    call tVecInit(upper_list_local,inblock_list%maxn)
    call tVecInit(inter_lower_local,inblock_list%maxn)
    call tVecInit(inter_upper_local,inblock_list%maxn)
    call tVecInit(inter_lower_tmp,inblock_list%maxn)
    call tVecInit(inter_upper_tmp,inblock_list%maxn)
    call tVecInit(stock,inblock_list%maxn)

    id_inblock(:) = 0
    allocate(internal_id_local(nnode),used(inblock_list%maxn),inorder(inblock_list%maxn))
    internal_id_local(:) = 0

    do i=0,numpart-1
      !make node list
      inblock_list%n = 0 !clear inblock_list
      do j=1,nnode
        if( npart(j) /= i ) cycle
        call tVecAdd(inblock_list,j)
        internal_id_local(j) = inblock_list%n
      end do

      ! label inter-block node
      lower_list_local%n = 0
      upper_list_local%n = 0
      used(:) = .false.
      do j=1,inblock_list%n
        nid1 = inblock_list%vec(j)
        pt1 = npart(nid1)
        lid = 0
        do k=1,n2nlist(nid1)%n
          nid2 = n2nlist(nid1)%vec(k)
          pt2 = npart(nid2)
          if( abs(pt1-pt2) /= 1 ) cycle
          lid = pt1 - pt2
          exit
        end do
        if( lid > 0 ) then
          call tVecAdd(lower_list_local,j)
          used(j) = .true.
        else if( lid < 0 ) then
          call tVecAdd(upper_list_local,j)
          used(j) = .true.
        end if
      end do
      if( lower_list_local%n + upper_list_local%n == 0 ) then
        call tVecAdd(lower_list_local,1)
        used(1) = .true.
      end if


      currnnode0 = lower_list_local%n+upper_list_local%n
      currnnode = currnnode0
      inter_lower_local%n = 0 !clear
      inter_upper_local%n = 0 !clear
      inter_lower_tmp%n = lower_list_local%n
      inter_lower_tmp%vec(1:inter_lower_tmp%n) = lower_list_local%vec(1:lower_list_local%n)
      inter_upper_tmp%n = upper_list_local%n
      do j=1,upper_list_local%n
         inter_upper_tmp%vec(j) = upper_list_local%vec(upper_list_local%n-j+1)
      end do
      do while( currnnode < inblock_list%n )
        !expand lower
        stock%n = 0
        if( inter_lower_tmp%n > 0 ) then
          do j=1,inter_lower_tmp%n
            !write(*,*) j,inter_lower_tmp%n,inter_lower_tmp%vec(j)
            nid1 = inblock_list%vec(inter_lower_tmp%vec(j))
            pt1 = npart(nid1)
            do k=1,n2nlist(nid1)%n
              nid2 = n2nlist(nid1)%vec(k)
              lid = internal_id_local(nid2)
              pt2 = npart(nid2)
              if( pt2 /= i ) cycle
              if( used(lid) ) cycle
              call tVecAdd(stock,lid)
              used(lid) = .true.
            end do
          end do
          inter_lower_tmp%n = 0
          do k=1,stock%n
            call tVecAdd(inter_lower_tmp,stock%vec(k))
            call tVecAdd(inter_lower_local,stock%vec(k))
          end do
        end if

        !expand upper
        stock%n = 0
        if( inter_upper_tmp%n > 0 ) then
          do j=inter_upper_tmp%n,1,-1
            nid1 = inblock_list%vec(inter_upper_tmp%vec(j))
            pt1 = npart(nid1)
            do k=n2nlist(nid1)%n,1,-1
              nid2 = n2nlist(nid1)%vec(k)
              lid = internal_id_local(nid2)
              pt2 = npart(nid2)
              if( pt2 /= i ) cycle
              if( used(lid) ) cycle
              call tVecAdd(stock,lid)
              used(lid) = .true.
            end do
          end do
          inter_upper_tmp%n = 0
          do k=stock%n,1,-1
            call tVecAdd(inter_upper_tmp,stock%vec(k))
            call tVecAdd(inter_upper_local,stock%vec(k))
          end do
        end if

        currnnode = currnnode0 + inter_upper_tmp%n+inter_lower_tmp%n
        if( currnnode0 == currnnode ) then
          do j=1,inblock_list%n
            if( used(j) ) cycle
            call tVecAdd(inter_lower_tmp,j)
            call tVecAdd(inter_lower_local,j)
            used(j) = .true.
            currnnode = currnnode + 1
            exit
          end do
        end if
        currnnode0 = currnnode
      end do

      inorder(:) = 0
      do j=1,lower_list_local%n
        inorder(j) = lower_list_local%vec(j)
      end do
      do j=1,inter_lower_local%n
        inorder(lower_list_local%n+j) = inter_lower_local%vec(j)
      end do

      do j=1,upper_list_local%n
        inorder(inblock_list%n-j+1) = upper_list_local%vec(j)
      end do
      do j=1,inter_upper_local%n
        inorder(inblock_list%n-upper_list_local%n-j+1) = inter_upper_local%vec(j)
      end do

      !check inorder
      !do k=1,inblock_list%n
      !  if( .not. used(k) ) write(*,*) "notused:",k
      !end do
      !used(:) = .false.
      !do k=1,inblock_list%n
      !  if( inorder(k) > inblock_list%n ) write(*,*) "inorder(k) > inblock_list%n:",inorder(k),inblock_list%n
      !  if( inorder(k) < 1 ) write(*,*) "inorder(k) < 1:",inorder(k)
      !  if( used(inorder(k)) ) write(*,*) "duplicate!:",inorder(k)
      !  used(inorder(k)) = .true.
      !end do

      ! reorder
      do k=1,inblock_list%n
        nid1 = inblock_list%vec(inorder(k))
        id_inblock(nid1) = k
        internal_id_local(nid1) = 0
      end do

    end do

    deallocate(internal_id_local)
    call tVecFinilize(inblock_list)
    call tVecFinilize(lower_list_local)
    call tVecFinilize(upper_list_local)
    call tVecFinilize(inter_lower_local)
    call tVecFinilize(inter_upper_local)

  end subroutine


  subroutine reorder_inblocks2(n2nlist,nodes,numpart,npart,id_inblock) !tmp
    type(tVector), allocatable, intent(in)         :: n2nlist(:)
    real(kind=kreal), pointer, intent(in  )        :: nodes(:) !< coordinate for inter ordering
    integer(kind=kint), intent(in)                 :: numpart
    integer(kind=kint), allocatable, intent(in)    :: npart(:)
    integer(kind=kint), allocatable, intent(inout) :: id_inblock(:)


    type(tVector) :: inblock_list, lower_list_local, upper_list_local

    integer(kind=kint) :: i, j, k, nid1, nid2, pt1, pt2, nnode, lid
    integer(kind=kint), allocatable :: flag(:,:)

    !for inblock ordering
    integer(kind=kint)              :: insize
    integer(kind=kint), allocatable :: inorder(:), gid(:)
    real(kind=kreal), allocatable   :: pos(:,:)

    nnode = size(n2nlist)
    allocate(pos(3,nnode),inorder(2*nnode/numpart),gid(2*nnode/numpart))

    call tVecInit(inblock_list,2*nnode/numpart)
    call tVecInit(lower_list_local,inblock_list%n)
    call tVecInit(upper_list_local,inblock_list%n)

    id_inblock(:) = 0

    do i=0,numpart-1

      !make node list
      inblock_list%n = 0 !clear inblock_list
      do j=1,nnode
        if( npart(j) /= i ) cycle
        call tVecAdd(inblock_list,j)
      end do

      ! label inter-block node

      allocate(flag(2,inblock_list%n))
      do j=1,inblock_list%n
        flag(1,j) = j
        flag(2,j) = 0
      end do
      do j=1,inblock_list%n
        nid1 = inblock_list%vec(j)
        pt1 = npart(nid1)
        do k=1,n2nlist(nid1)%n
          nid2 = n2nlist(nid1)%vec(k)
          pt2 = npart(nid2)
          if( pt1 == pt2 ) cycle
          flag(2,j) = flag(2,j)+100*numpart/(pt2-pt1)
        end do
      end do
      call sort_int_array2(flag, 1, inblock_list%n)

      ! two-sided CM
      call two_sided_CM(flag,n2nlist,inblock_list,i,npart)

      ! reorder
      do k=1,inblock_list%n
        nid1 = inblock_list%vec(flag(1,k))
        id_inblock(nid1) = k
      end do

      deallocate(flag)
    end do

    deallocate(pos,inorder,gid)
    call tVecFinilize(inblock_list)
  end subroutine

  subroutine two_sided_CM(flag,n2nlist,inblock_list,pid,npart)
    integer(kind=kint), allocatable, intent(in)    :: flag(:,:)
    type(tVector), allocatable, intent(in)         :: n2nlist(:)
    type(tVector), intent(in)                      :: inblock_list
    integer(kind=kint), intent(in)                 :: pid
    integer(kind=kint), allocatable, intent(in)    :: npart(:)


  end subroutine

  subroutine part_by_metis(hecMESH,n2nlist,nparts,part)
    use iso_c_binding
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type(tVector), allocatable, intent(in) :: n2nlist(:)
    integer(kind=kint), intent(in)        :: nparts
    integer(kind=kint), allocatable, intent(inout)     :: part(:)

    integer              :: i, j, nid
    logical              :: is_external

    !for metis
    integer              :: nvtxs, ncon
    integer, allocatable :: xadj(:), adjncy(:)
    integer, allocatable :: vwgt(:)
    !type(c_ptr)          :: vwgt = c_null_ptr
    type(c_ptr)          :: vsize = c_null_ptr
    integer, allocatable :: adjwgt(:)
    type(c_ptr)          :: adjwgt_dum = c_null_ptr
    real, allocatable    :: tpwgts(:)
    real                 :: ubvec
    integer, allocatable :: mopts(:)
    integer              :: objval

    nvtxs = hecMESH%nn_internal
    ncon = 0
    do i=1,nvtxs
      do j=1,n2nlist(i)%n
        nid = n2nlist(i)%vec(j)
        if( nid > hecMESH%nn_internal ) cycle
        ncon = ncon + 1
      end do
      !ncon = ncon + n2nlist(i)%n
    end do
    allocate(xadj(nvtxs+1),adjncy(ncon),adjwgt(ncon),tpwgts(nparts),mopts(41))
    allocate(vwgt(nvtxs))
    vwgt = 1
    xadj(1) = 0
    ncon = 0
    do i=1,nvtxs
      is_external = ( i > hecMESH%nn_internal )
      do j=1,n2nlist(i)%n
        nid = n2nlist(i)%vec(j)
        if( nid > hecMESH%nn_internal ) cycle
        ncon = ncon + 1
        adjncy(ncon) = nid-1
        if( nid > nvtxs ) write(*,*) "nid > nvtxs!: ",nid,nvtxs
        if( nid < 1 ) write(*,*) "nid < 1!: ",nid
        if( is_external .or. ( nid> hecMESH%nn_internal ) ) then
          adjwgt(ncon) = 4
        else
          adjwgt(ncon) = 1
        end if
      end do
      xadj(i+1) = ncon
    end do
    ubvec = 1.001
    tpwgts(:) = 1.0 / real(nparts)
    mopts(:)  = -1
    ncon = 1

    do i=1,nvtxs
      part(i) = 0
    end do
    if( nparts > 1 ) then
      call METIS_PartGraphRecursive(nvtxs,     &  ! (in), int
                                    ncon,      &  ! (in), int
                                    xadj,      &  ! (in), int(:)
                                    adjncy,    &  ! (in), int(:)
                                    vwgt,      &  ! (in), int(:)
                                    vsize,     &  ! (in), int(:)
                                    adjwgt,    &  ! (in), int(:)
                                    nparts,    &  ! (in), int(:)
                                    tpwgts,    &  ! (in), real(:)
                                    ubvec,     &  ! (in), real(:)
                                    mopts,     &  ! (in), int(:)
                                    objval,    &  ! (out) int(:)
                                    part)         ! (out) int(:)

      !call METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt_dum,nparts,tpwgts,ubvec,mopts,objval,part)
      !call METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize,adjwgt,tpwgts,nparts,ubvec,mopts,objval,part)
    end if

    do i=hecMESH%nn_internal+1,hecMESH%n_node
      part(i) = -1
    end do

  end subroutine

  subroutine part_by_metis_graph(hecMESH,n2nlist,numpart,npart)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type(tVector), allocatable, intent(in) :: n2nlist(:)
    integer(kind=kint), intent(in)        :: numpart
    integer(kind=kint), allocatable, intent(inout)     :: npart(:)

    !for metis
    integer              :: nels, nnds, nitem, i, n
    integer, allocatable :: eptr(:), nodes(:), epart(:)
    integer, allocatable :: vwgt, vsize, mopts
    real(8), allocatable :: tpwgts

    nnds = hecMESH%n_node
    nels = hecMESH%n_elem
    nitem = hecMESH%elem_node_index(nels)

    allocate(nodes(nitem),epart(nels))

    if( numpart > 1 ) then
      do i=1,nitem
        nodes(i) = hecMESH%elem_node_item(i)-1
      end do
      call METIS_PartMeshNodal(nels,nnds,hecMESH%elem_node_index,nodes,vwgt,vsize,numpart,tpwgts,mopts,n,epart,npart)
    else
      do i=1,nnds
        npart(i) = 0
      end do
      do i=1,nels
        epart(i) = 0
      end do
    end if

    do i=hecMESH%nn_internal+1,nnds
      npart(i) = -1
    end do

    deallocate(nodes,epart)

  end subroutine

  !!!!!!! --- add connectivity_reorder start --- !!!!!!!
  subroutine connectivity_reorder(hecMESH,num_reorder)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(inout)     :: num_reorder(:)

    integer(kind=kint) :: i
    type(tVector), allocatable :: n2nlist(:)
    integer(kind=kint) :: n_groups
    type(tVector), allocatable :: groups(:)

    !create node to node list
    call make_n2nlist(hecMESH,n2nlist)

    !create group info
    call make_grouping(hecMESH,n2nlist,n_groups,groups)

    !create new node number
    call make_connectivity_reorder(hecMESH,n2nlist,n_groups,groups,num_reorder)

    !finalize n2nlist
    do i=1,hecMESH%n_node
      !call tVecPrint(n2nlist(i))
      call tVecFinilize(n2nlist(i))
    end do
    deallocate(n2nlist)
    !finalize grouplist
    do i=1,n_groups
      call tVecFinilize(groups(i))
    end do
    deallocate(groups)

  end subroutine

  subroutine make_connectivity_reorder(hecMESH,n2nlist,n_groups,groups,num_reorder)
    type (hecmwST_local_mesh), intent(in)     :: hecMESH
    type(tVector), allocatable, intent(in)    :: n2nlist(:)
    integer(kind=kint), intent(in)            :: n_groups
    type(tVector), allocatable, intent(in)    :: groups(:)
    integer(kind=kint), intent(inout)         :: num_reorder(:)

    integer(kind=kint) :: i, j, k, nid, nid2, counter
    integer(kind=kint), allocatable :: newnumber(:)
    integer(kind=kint), allocatable :: sortarray(:,:)

    allocate(newnumber(hecMESH%n_node))
    newnumber(1:hecMESH%nn_internal) = 0
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      newnumber(i) = i
    end do

    allocate(sortarray(2,hecMESH%nn_internal))
    counter = 0
    do i=2,n_groups
      sortarray(1:2,1:groups(i)%n) = 0
      do j=1,groups(i)%n
        nid = groups(i)%vec(j)
        sortarray(1,j) = nid
        do k=1,n2nlist(nid)%n
          nid2 = n2nlist(nid)%vec(k)
          sortarray(2,j) = max(newnumber(nid2),sortarray(2,j))
        end do
      end do

      call sort_int_array2(sortarray, 1, groups(i)%n)
      do j=1,groups(i)%n
        nid = sortarray(1,groups(i)%n-j+1)
        if( newnumber(nid) > 0 ) cycle
        newnumber(nid) = hecMESH%nn_internal-counter
        num_reorder(nid) = newnumber(nid)
        counter = counter + 1
      end do

    end do

    deallocate(newnumber,sortarray)
  end subroutine

  subroutine make_grouping(hecMESH,n2nlist,n_groups,groups)
    type (hecmwST_local_mesh), intent(in)     :: hecMESH
    type(tVector), allocatable, intent(in)    :: n2nlist(:)
    integer(kind=kint), intent(out)           :: n_groups
    type(tVector), allocatable, intent(inout) :: groups(:)

    integer(kind=kint), allocatable :: grpids(:)
    integer(kind=kint) :: i, j, k, nid0, nid1, nnode_rest
    integer(kind=kint) :: counter

    allocate(grpids(hecMESH%n_node))
    allocate(groups(hecMESH%nn_internal+1))

    !initial group is import nodes
    grpids(1:hecMESH%nn_internal) = -1
    grpids(hecMESH%nn_internal+1:hecMESH%n_node) = 1
    n_groups = 1
    call tVecInit(groups(n_groups),hecMESH%n_node-hecMESH%nn_internal)
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      call tVecAdd(groups(n_groups),i)
    end do
    n_groups = n_groups + 1
    nnode_rest = hecMESH%nn_internal

    do k=2,hecMESH%nn_internal
      call tVecInit(groups(n_groups),1000)
      do i=1,groups(n_groups-1)%n
        nid0 = groups(n_groups-1)%vec(i)

        ! connectivity loop
        do j=1,n2nlist(nid0)%n
          nid1 = n2nlist(nid0)%vec(j)
          if( grpids(nid1) > 0 ) cycle

          grpids(nid1) = n_groups
          call tVecAdd(groups(n_groups),nid1)
        end do
      end do

      if( groups(n_groups)%n == 0 ) then
        do i=1,hecMESH%nn_internal
          if( grpids(i) > 0 ) cycle
          call tVecAdd(groups(n_groups),i)
        end do
      end if

      nnode_rest = nnode_rest - groups(n_groups)%n
      if( nnode_rest == 0 ) exit
      n_groups = n_groups + 1
    end do

  end subroutine

  subroutine make_n2nlist(hecMESH,n2nlist)
    use hecmw_etype
    type (hecmwST_local_mesh), intent(in)     :: hecMESH
    type(tVector), allocatable, intent(inout) :: n2nlist(:)

    integer(kind=kint) :: itype, iS, iE, ic_type, nn, icel, iiS
    integer(kind=kint) :: i, j, nodlocal(20)

    allocate(n2nlist(hecMESH%n_node))
    do i=1,hecMESH%n_node
      call tVecInit(n2nlist(i),100)
    end do

    ! element loop
    do itype= 1, hecMESH%n_elem_type
      iS= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      nn = hecmw_get_max_node(ic_type)

      do icel= is, iE
        iiS= hecMESH%elem_node_index(icel-1)
        nodLOCAL(1:nn)= hecMESH%elem_node_item (iiS+1:iiS+nn)
        do i=1,nn
          do j=1,nn
            call tVecAdd(n2nlist(nodLOCAL(i)),nodLOCAL(j))
          end do
        end do
      end do
    end do

    do i=1,hecMESH%n_node
      call tVecCompress(n2nlist(i))
    end do
  end subroutine

  !!!!!!! --- add connectivity_reorder end --- !!!!!!!

  subroutine block_reorder(num,numx,numy,numz,ratio,num_reorder,pos)
    implicit none

    integer(kind=kint), intent(IN)  :: num,numx,numy,numz
    real(kind=kreal),   intent(IN)  :: ratio
    integer(kind=kint), intent(OUT) :: num_reorder(:)
    real(kind=kreal),   intent(IN)  :: pos(:,:) ! pos(3,num)

    integer(kind=kint), dimension(:,:,:),  allocatable :: icount_in
    integer(kind=kint), dimension(:,:,:),  allocatable :: icount_out
    integer(kind=kint), dimension(:,:,:,:),allocatable :: num_new_in
    integer(kind=kint), dimension(:,:,:,:),allocatable :: num_new_out

    real(kind=kreal) :: xunit, yunit, zunit
    real(kind=kreal) :: xl_ratio, yl_ratio, zl_ratio
    real(kind=kreal) :: posx,     posy,     posz
    real(kind=kreal) :: postx,    posty,    postz
    real(kind=kreal) :: posx_max, posy_max, posz_max
    real(kind=kreal) :: posx_min, posy_min, posz_min
    real(kind=kreal) :: avec, sigma

    integer(kind=kint)  :: numx_pos, numy_pos, numz_pos
    integer(kind=kint)  :: max_num_box
    integer(kind=kint)  :: iout,iin

    integer :: i
    integer :: ii, jj, kk
    integer :: ic, icc, ic1, ic2
    integer :: nblock, minc, maxc, ll, iold


    !! param
    !     parameter(xl_ratio=0.1,yl_ratio=0.1,zl_ratio=0.1)

    xl_ratio = ratio
    yl_ratio = ratio
    zl_ratio = ratio

    !! allocate section
    allocate(icount_in(numx,numy,numz))
    allocate(icount_out(numx,numy,numz))
    !! search max,min position for x-y-z direction
    posx_max=pos(1,1)
    posy_max=pos(2,1)
    posz_max=pos(3,1)
    posx_min=pos(1,1)
    posy_min=pos(2,1)
    posz_min=pos(3,1)
    do i = 1,num
      posx=pos(1,i)
      posy=pos(2,i)
      posz=pos(3,i)
      if(posx < posx_min ) then
        posx_min = posx
      endif
      if(posy < posy_min ) then
        posy_min = posy
      endif
      if(posz < posz_min ) then
        posz_min = posz
      endif
      if(posx_max < posx ) then
        posx_max = posx
      endif
      if(posy_max < posy ) then
        posy_max = posy
      endif
      if(posz_max < posz ) then
        posz_max = posz
      endif
    end do
    !! calculation x-y-z unit
    xunit=(posx_max - posx_min)/numx
    yunit=(posy_max - posy_min)/numy
    zunit=(posz_max - posz_min)/numz
    !    9 format(a,4(1x,a10))
    !      write(*, 9) "       ", "min", "max", "range", "pitch"
    !   10 format(a,4(1x,f10.5))
    !      write(*,10) "posx = ", posx_min,posx_max,posx_max-posx_min,xunit
    !      write(*,10) "posy = ", posy_min,posy_max,posy_max-posy_min,yunit
    !      write(*,10) "posz = ", posz_min,posz_max,posz_max-posz_min,zunit


    !! count number in box points
    icount_in = 0
    do i = 1,num
      posx=pos(1,i) - posx_min
      posy=pos(2,i) - posy_min
      posz=pos(3,i) - posz_min
      if(posx-1.0e-10 <= 0.0) then
        numx_pos=1
      else
        numx_pos=int((posx-1.0e-10)/xunit)+1
      end if
      if(numx_pos < 0 ) then
        write(6,*) 'minus x_number o!!urs i=',i
      endif
      if(posy-1.0e-10 <= 0.0) then
        numy_pos=1
      else
        numy_pos=int((posy-1.0e-10)/yunit)+1
      end if
      if(numy_pos < 0 ) then
        write(6,*) 'minus y_number o!!urs i=',i
      endif
      if(posz-1.0e-10 <= 0.0) then
        numz_pos=1
      else
        numz_pos=int((posz-1.0e-10)/zunit)+1
      end if
      if(numz_pos < 0 ) then
        write(6,*) 'minus z_number o!!urs i=',i
      endif
      icount_in(numx_pos,numy_pos,numz_pos) = &
          icount_in(numx_pos,numy_pos,numz_pos)+1
    end do


    !! check
    !      do kk=1,numz
    !        do jj=1,numy
    !          do ii=1,numx
    !   11 format(4(a,i0))
    !        if( icount_in(ii,jj,kk) > 0 ) then
    !        write(6,11) 'icount(',ii,',',jj,',',kk,')=',icount_in(ii,jj,kk)
    !         endif
    !          end do
    !        end do
    !      end do

    ! search max number of box points
    max_num_box = 0
    do kk=1,numz
      do jj=1,numy
        do ii=1,numx
          ic=icount_in(ii,jj,kk)
          if(max_num_box < ic ) then
            max_num_box = ic
          end if
        end do
      end do
    end do


    !! allocate section
    allocate(num_new_in(max_num_box,numx,numy,numz))
    allocate(num_new_out(max_num_box,numx,numy,numz))
    !! store new number
    icount_in = 0
    icount_out = 0
    num_new_in = 0
    num_new_out = 0
    iin = 0
    iout = 0
    do i = 1,num
      posx=pos(1,i) - posx_min
      posy=pos(2,i) - posy_min
      posz=pos(3,i) - posz_min
      if(posx-1.0e-10 <= 0.0) then
        numx_pos=1
        postx=0.0
      else
        numx_pos=int((posx-1.0e-10)/xunit)+1
        postx=posx-(numx_pos-1)*xunit
      end if
      if(posy-1.0e-10 <= 0.0) then
        numy_pos=1
        posty=0.0
      else
        numy_pos=int((posy-1.0e-10)/yunit)+1
        posty=posy-(numy_pos-1)*yunit
      end if
      if(posz-1.0e-10 <= 0.0) then
        numz_pos=1
        postz=0.0
      else
        numz_pos=int((posz-1.0e-10)/zunit)+1
        postz=posz-(numz_pos-1)*zunit
      end if
      if((xunit*xl_ratio <= postx)        .and. &
         (postx <= xunit*(1.0-xl_ratio))  .and. &
         (yunit*yl_ratio <= posty)        .and. &
         (posty <= yunit*(1.0-yl_ratio))  .and. &
         (zunit*zl_ratio <= postz)        .and. &
         (postz <= zunit*(1.0-zl_ratio))) then
        icount_in(numx_pos,numy_pos,numz_pos) = &
            icount_in(numx_pos,numy_pos,numz_pos)+1
        num_new_in(icount_in(numx_pos,numy_pos,numz_pos) &
            ,numx_pos,numy_pos,numz_pos)=i
        iin= iin+1
      else
        icount_out(numx_pos,numy_pos,numz_pos) = &
            icount_out(numx_pos,numy_pos,numz_pos)+1
        num_new_out(icount_out(numx_pos,numy_pos,numz_pos) &
            ,numx_pos,numy_pos,numz_pos)=i
        iout= iout+1
      end if
    end do


    ! make renumber list
    !      num_reorder = 0
    icc=0
    nblock=0
    minc=num
    maxc=0
    do kk=1,numz
      do jj=1,numy
        do ii=1,numx
          ic1=icount_in(ii,jj,kk)
          do ll=1,ic1
            iold=num_new_in(ll,ii,jj,kk)
            icc = icc + 1
            num_reorder(iold)=icc
            !              write(6,*) icc,iold !DEBUG
          end do
          ic2=icount_out(ii,jj,kk)
          do ll=1,ic2
            iold=num_new_out(ll,ii,jj,kk)
            icc = icc + 1
            num_reorder(iold)=icc
            !              write(6,*) icc,iold !DEBUG
          end do
          ic=ic1+ic2
          if(ic>0) then
            nblock = nblock + 1
            minc = min(ic,minc)
            maxc = max(ic,maxc)
          endif
        end do
      end do
    end do


    !! count
    !      avec = dble(num)/dble(nblock)
    !      if(.true.)then
    !        sigma = 0.d0
    !        do kk=1,numz
    !          do jj=1,numy
    !            do ii=1,numx
    !              ic1=icount_in(ii,jj,kk)
    !              do ll=1,ic1
    !                iold=num_new_in(ll,ii,jj,kk)
    !              end do
    !              ic2=icount_out(ii,jj,kk)
    !              do ll=1,ic2
    !                iold=num_new_out(ll,ii,jj,kk)
    !              end do
    !              ic=ic1+ic2
    !              if(ic>0) then
    !                sigma = sigma + (ic-avec)**2
    !              endif
    !            end do
    !          end do
    !        end do
    !        sigma = sqrt(sigma/nblock)
    !        write(*,*) 'number of block  = ', nblock
    !        write(*,*) '  ave.grid/block = ', avec
    !        write(*,*) '  max            = ', maxc
    !        write(*,*) '  min            = ', minc
    !        write(*,*) '  sigma          = ', sigma
    !      endif

    ! check
    !      write(6,*) 'i!!=',icc
    !      call qsort(num_reorder,num,4,compare4)
    !      do i = 1,num
    !        num_chk=num_reorder(i)
    !        if(i/=num_chk) then
    !          write(6,*) i,num_chk,'+++ error co!!ured +++'
    !        end if
    !      end do

    !      write(6,*) '+++ execution OK +++'
    !      write(6,*) '## in - out ##',iin,iout

    deallocate(icount_in)
    deallocate(icount_out)
    deallocate(num_new_in)
    deallocate(num_new_out)
    return
  end subroutine block_reorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! For dubug.
  ! Make reorder table as inverse order.
  subroutine make_reorder_table_inverse_order(hecMESH)

    type (hecmwST_local_mesh), intent(inout) :: hecMESH

    integer :: i

    do i=1, hecMESH%nn_internal
      hecMESH%tuning_block_reorder_new2old(i)=hecMESH%nn_internal - i + 1
      hecMESH%tuning_block_reorder_old2new(hecMESH%tuning_block_reorder_new2old(i))=i
    end do
    do i=hecMESH%nn_internal + 1, hecMESH%n_node
      hecMESH%tuning_block_reorder_new2old(i)=i
      hecMESH%tuning_block_reorder_old2new(hecMESH%tuning_block_reorder_new2old(i))=i
    end do
  end subroutine make_reorder_table_inverse_order

end module hecmw_tuning_node_block_reorder
