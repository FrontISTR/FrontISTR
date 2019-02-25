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


contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !C
  !C***
  !C*** SUBROUTINES
  !C***
  !C


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
    allocate(new_elem_node_item(size(hecMESH%elem_node_item)))
    do i=1, size(hecMESH%elem_node_item)
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

    call block_reorder(num,numx,numy,numz,ratio,num_reorder,pos)

    ! Because reorder table for external nodes are already initialized, 
    ! only internal nodes are updated. 
    do i=1,num
      hecMESH%tuning_block_reorder_old2new(i) = num_reorder(i)
      hecMESH%tuning_block_reorder_new2old(num_reorder(i)) = i
    end do

    deallocate(num_reorder)
    deallocate(pos)

  end subroutine make_reorder_table



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
