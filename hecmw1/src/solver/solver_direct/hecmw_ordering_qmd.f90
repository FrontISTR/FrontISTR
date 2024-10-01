!-------------------------------------------------------------------------------
! Copyright (c) 2020 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------
!> @brief HECMW_ORDERING_QMD is a program for the minimum degree
!         ordering using quotient graphs
!----------------------------------------------------------------------
module hecmw_ordering_qmd
  use hecmw_util
  implicit none

  private
  public :: hecmw_ordering_genqmd

contains

  !======================================================================!
  !> @brief hecmw_ordering_GENQMD
  !======================================================================!
  subroutine hecmw_ordering_GENQMD(Neqns,Nttbr,Xadj,Adj0,Perm,Invp)
    implicit none
    !------
    integer(kind=kint), intent(in):: Neqns
    integer(kind=kint), intent(in):: Nttbr
    integer(kind=kint), intent(in):: Adj0(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(out):: Perm(:)
    integer(kind=kint), intent(out):: Invp(:)
    !------
    integer(kind=kint), allocatable:: Deg(:)
    integer(kind=kint), allocatable:: Marker(:)
    integer(kind=kint), allocatable:: Rchset(:)
    integer(kind=kint), allocatable:: Nbrhd(:)
    integer(kind=kint), allocatable:: Qsize(:)
    integer(kind=kint), allocatable:: Qlink(:)
    integer(kind=kint), allocatable:: Adjncy(:)
    integer(kind=kint):: inode
    integer(kind=kint):: ip
    integer(kind=kint):: irch
    integer(kind=kint):: j
    integer(kind=kint):: mindeg
    integer(kind=kint):: ndeg
    integer(kind=kint):: nhdsze
    integer(kind=kint):: node
    integer(kind=kint):: np
    integer(kind=kint):: num
    integer(kind=kint):: nump1
    integer(kind=kint):: nxnode
    integer(kind=kint):: rchsze
    integer(kind=kint):: search
    integer(kind=kint):: thresh
    integer(kind=kint):: ierror
    logical:: found

    allocate (DEG(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, deg: SUB. genqmd"
    allocate (MARker(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, marker: SUB. genqmd"
    allocate (RCHset(NEQns+2),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, rchset: SUB. genqmd"
    allocate (NBRhd(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, nbrhd: SUB. genqmd"
    allocate (QSIze(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, qsize: SUB. genqmd"
    allocate (QLInk(NEQns+1),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, qlink: SUB. genqmd"
    allocate (ADJncy(2*NTTbr),stat=IERror)
    if ( IERror/=0 ) stop "ALLOCATION ERROR, adjncy: SUB. genqmd"

    mindeg = Neqns
    Adjncy(1:Xadj(Neqns+1) - 1) = Adj0(1:Xadj(Neqns+1) - 1)
    do node = 1, Neqns
      Perm(node) = node
      Invp(node) = node
      Marker(node) = 0
      Qsize(node) = 1
      Qlink(node) = 0
      ndeg = Xadj(node+1) - Xadj(node)
      Deg(node) = ndeg
      if ( ndeg<mindeg ) mindeg = ndeg
    enddo

    num = 0
    loop1: do
      search = 1
      thresh = mindeg
      mindeg = Neqns
      loop2: do
        nump1 = num + 1
        if ( nump1>search ) search = nump1
        found = .false.
        do j = search, Neqns
          node = Perm(j)
          if ( Marker(node)>=0 ) then
            ndeg = Deg(node)
            if ( ndeg<=thresh ) then
              found = .true.
              exit
            endif
            if ( ndeg<mindeg ) mindeg = ndeg
          endif
        enddo
        if (.not. found) cycle loop1

        search = j
        Marker(node) = 1
        call QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
        nxnode = node
        do
          num = num + 1
          np = Invp(nxnode)
          ip = Perm(num)
          Perm(np) = ip
          Invp(ip) = np
          Perm(num) = nxnode
          Invp(nxnode) = num
          Deg(nxnode) = -1
          nxnode = Qlink(nxnode)
          if ( nxnode<=0 ) then
            if ( rchsze>0 ) then
              call QMDUPD(Xadj,Adjncy,rchsze,Rchset,Deg,Qsize,Qlink,Marker,Rchset(rchsze+1:),Nbrhd(nhdsze+1:))
              Marker(node) = 0
              do irch = 1, rchsze
                inode = Rchset(irch)
                if ( Marker(inode)>=0 ) then
                  Marker(inode) = 0
                  ndeg = Deg(inode)
                  if ( ndeg<mindeg ) mindeg = ndeg
                  if ( ndeg<=thresh ) then
                    mindeg = thresh
                    thresh = ndeg
                    search = Invp(inode)
                  endif
                endif
              enddo
              if ( nhdsze>0 ) call QMDOT(node,Xadj,Adjncy,Marker,rchsze,Rchset,Nbrhd)
            endif
            if ( num>=Neqns ) exit
            cycle loop2
          endif
        enddo
        exit
      enddo loop2
      exit
    enddo loop1

    deallocate (DEG)
    deallocate (MARker)
    deallocate (RCHset)
    deallocate (NBRhd)
    deallocate (QSIze)
    deallocate (QLInk)
    deallocate (ADJncy)
  end subroutine HECMW_ORDERING_GENQMD

  !======================================================================!
  !> @brief QMDRCH
  !======================================================================!
  subroutine QMDRCH(Root,Xadj,Adjncy,Deg,Marker,Rchsze,Rchset,Nhdsze,Nbrhd)
    implicit none
    !------
    integer(kind=kint), intent(in):: Root
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: Deg(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(out):: Nhdsze
    integer(kind=kint), intent(out):: Rchsze
    integer(kind=kint), intent(inout):: Marker(:)
    integer(kind=kint), intent(out):: Rchset(:)
    integer(kind=kint), intent(out):: Nbrhd(:)
    !------
    integer(kind=kint):: i
    integer(kind=kint):: istrt
    integer(kind=kint):: istop
    integer(kind=kint):: j
    integer(kind=kint):: jstrt
    integer(kind=kint):: jstop
    integer(kind=kint):: nabor
    integer(kind=kint):: node

    Nhdsze = 0
    Rchsze = 0
    istrt = Xadj(Root)
    istop = Xadj(Root+1) - 1
    if ( istop<istrt ) return
    do i = istrt, istop
      nabor = Adjncy(i)
      if ( nabor==0 ) return
      if ( Marker(nabor)==0 ) then
        if ( Deg(nabor)<0 ) then
          Marker(nabor) = -1
          Nhdsze = Nhdsze + 1
          Nbrhd(Nhdsze) = nabor
          loop1: do
            jstrt = Xadj(nabor)
            jstop = Xadj(nabor+1) - 1
            do j = jstrt, jstop
              node = Adjncy(j)
              nabor = -node
              if ( node<0 ) cycle loop1
              if ( node==0 ) exit
              if ( Marker(node)==0 ) then
                Rchsze = Rchsze + 1
                Rchset(Rchsze) = node
                Marker(node) = 1
              endif
            enddo
            exit
          enddo loop1
        else
          Rchsze = Rchsze + 1
          Rchset(Rchsze) = nabor
          Marker(nabor) = 1
        endif
      endif
    enddo
  end subroutine QMDRCH

  !======================================================================!
  !> @brief QMDUPD
  !======================================================================!
  subroutine QMDUPD(Xadj,Adjncy,Nlist,List,Deg,Qsize,Qlink,Marker,Rchset,Nbrhd)
    implicit none
    !------
    integer(kind=kint), intent(in):: Nlist
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: List(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(inout):: Deg(:)
    integer(kind=kint), intent(inout):: Marker(:)
    integer(kind=kint), intent(out):: Rchset(:)
    integer(kind=kint), intent(out):: Nbrhd(:)
    integer(kind=kint), intent(inout):: Qsize(:)
    integer(kind=kint), intent(inout):: Qlink(:)
    !------
    integer(kind=kint):: deg0
    integer(kind=kint):: deg1
    integer(kind=kint):: il
    integer(kind=kint):: inhd
    integer(kind=kint):: inode
    integer(kind=kint):: irch
    integer(kind=kint):: j
    integer(kind=kint):: jstrt
    integer(kind=kint):: jstop
    integer(kind=kint):: mark
    integer(kind=kint):: nabor
    integer(kind=kint):: nhdsze
    integer(kind=kint):: node
    integer(kind=kint):: rchsze

    if ( Nlist<=0 ) return
    deg0 = 0
    nhdsze = 0
    do il = 1, Nlist
      node = List(il)
      deg0 = deg0 + Qsize(node)
      jstrt = Xadj(node)
      jstop = Xadj(node+1) - 1
      do j = jstrt, jstop
        nabor = Adjncy(j)
        if ( Marker(nabor)==0 .and. Deg(nabor)<0 ) then
          Marker(nabor) = -1
          nhdsze = nhdsze + 1
          Nbrhd(nhdsze) = nabor
        endif
      enddo
    enddo

    if ( nhdsze>0 ) call QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,deg0,nhdsze,Nbrhd,Rchset,Nbrhd(nhdsze+1:))
    do il = 1, Nlist
      node = List(il)
      mark = Marker(node)
      if ( mark<=1 .and. mark>=0 ) then
        call QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
        deg1 = deg0
        if ( rchsze>0 ) then
          do irch = 1, rchsze
            inode = Rchset(irch)
            deg1 = deg1 + Qsize(inode)
            Marker(inode) = 0
          enddo
        endif
        Deg(node) = deg1 - 1
        if ( nhdsze>0 ) then
          do inhd = 1, nhdsze
            inode = Nbrhd(inhd)
            Marker(inode) = 0
          enddo
        endif
      endif
    enddo
  end subroutine QMDUPD

  !======================================================================!
  !> @brief QMDMRG
  !======================================================================!
  subroutine QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,Deg0,Nhdsze,Nbrhd,Rchset,Ovrlp)
    implicit none
    !------
    integer(kind=kint), intent(in):: Deg0
    integer(kind=kint), intent(in):: Nhdsze
    integer(kind=kint), intent(in):: Adjncy(:)
    integer(kind=kint), intent(in):: Nbrhd(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(inout):: Deg(:)
    integer(kind=kint), intent(inout):: Qsize(:)
    integer(kind=kint), intent(inout):: Qlink(:)
    integer(kind=kint), intent(inout):: Marker(:)
    integer(kind=kint), intent(out):: Rchset(:)
    integer(kind=kint), intent(out):: Ovrlp(:)
    !------
    integer(kind=kint):: deg1
    integer(kind=kint):: head
    integer(kind=kint):: inhd
    integer(kind=kint):: iov
    integer(kind=kint):: irch
    integer(kind=kint):: j
    integer(kind=kint):: jstrt
    integer(kind=kint):: jstop
    integer(kind=kint):: link
    integer(kind=kint):: lnode
    integer(kind=kint):: mark
    integer(kind=kint):: mrgsze
    integer(kind=kint):: nabor
    integer(kind=kint):: node
    integer(kind=kint):: novrlp
    integer(kind=kint):: rchsze
    integer(kind=kint):: root

    if ( Nhdsze<=0 ) return
    do inhd = 1, Nhdsze
      root = Nbrhd(inhd)
      Marker(root) = 0
    enddo
    do inhd = 1, Nhdsze
      root = Nbrhd(inhd)
      Marker(root) = -1
      rchsze = 0
      novrlp = 0
      deg1 = 0
      loop1: do
        jstrt = Xadj(root)
        jstop = Xadj(root+1) - 1
        do j = jstrt, jstop
          nabor = Adjncy(j)
          root = -nabor
          if ( nabor<0 ) cycle loop1
          if ( nabor==0 ) exit
          mark = Marker(nabor)

          if ( mark>=0 ) then
            if ( mark<=0 ) then
              rchsze = rchsze + 1
              Rchset(rchsze) = nabor
              deg1 = deg1 + Qsize(nabor)
              Marker(nabor) = 1
            elseif ( mark<=1 ) then
              novrlp = novrlp + 1
              Ovrlp(novrlp) = nabor
              Marker(nabor) = 2
            endif
          endif
        enddo
        exit
      enddo loop1
      head = 0
      mrgsze = 0
      loop2: do iov = 1, novrlp
        node = Ovrlp(iov)
        jstrt = Xadj(node)
        jstop = Xadj(node+1) - 1
        do j = jstrt, jstop
          nabor = Adjncy(j)
          if ( Marker(nabor)==0 ) then
            Marker(node) = 1
            cycle loop2
          endif
        enddo
        mrgsze = mrgsze + Qsize(node)
        Marker(node) = -1
        lnode = node
        do
          link = Qlink(lnode)
          if ( link<=0 ) then
            Qlink(lnode) = head
            head = node
            exit
          else
            lnode = link
          endif
        enddo
      enddo loop2
      if ( head>0 ) then
        Qsize(head) = mrgsze
        Deg(head) = Deg0 + deg1 - 1
        Marker(head) = 2
      endif
      root = Nbrhd(inhd)
      Marker(root) = 0
      if ( rchsze>0 ) then
        do irch = 1, rchsze
          node = Rchset(irch)
          Marker(node) = 0
        enddo
      endif
    enddo
  end subroutine QMDMRG

  !======================================================================!
  !> @brief QMDOT
  !======================================================================!
  subroutine QMDOT(Root,Xadj,Adjncy,Marker,Rchsze,Rchset,Nbrhd)
    implicit none
    !------
    integer(kind=kint), intent(in):: Rchsze
    integer(kind=kint), intent(in):: Root
    integer(kind=kint), intent(in):: Marker(:)
    integer(kind=kint), intent(in):: Rchset(:)
    integer(kind=kint), intent(in):: Nbrhd(:)
    integer(kind=kint), intent(in):: Xadj(:)
    integer(kind=kint), intent(inout):: Adjncy(:)
    !------
    integer(kind=kint):: inhd
    integer(kind=kint):: irch
    integer(kind=kint):: j
    integer(kind=kint):: jstrt
    integer(kind=kint):: jstop
    integer(kind=kint):: link
    integer(kind=kint):: nabor
    integer(kind=kint):: node

    irch = 0
    inhd = 0
    node = Root
    loop1: do
      jstrt = Xadj(node)
      jstop = Xadj(node+1) - 2
      if ( jstop>=jstrt ) then
        do j = jstrt, jstop
          irch = irch + 1
          Adjncy(j) = Rchset(irch)
          if ( irch>=Rchsze ) exit loop1
        enddo
      endif
      link = Adjncy(jstop+1)
      node = -link
      if ( link>=0 ) then
        inhd = inhd + 1
        node = Nbrhd(inhd)
        Adjncy(jstop+1) = -node
      endif
    enddo loop1

    Adjncy(j+1) = 0
    do irch = 1, Rchsze
      node = Rchset(irch)
      if ( Marker(node)>=0 ) then
        jstrt = Xadj(node)
        jstop = Xadj(node+1) - 1
        do j = jstrt, jstop
          nabor = Adjncy(j)
          if ( Marker(nabor)<0 ) then
            Adjncy(j) = Root
            exit
          endif
        enddo
      endif
    enddo
  end subroutine QMDOT

end module hecmw_ordering_qmd
