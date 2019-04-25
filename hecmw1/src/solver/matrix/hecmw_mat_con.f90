!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_matrix_con
  use hecmw_util
  implicit none

  private

  public :: hecmw_mat_con

  integer(kind=kint) :: NU, NL
  integer(kind=kint), pointer :: INL(:), INU(:)
  integer(kind=kint), pointer :: IAL(:,:), IAU(:,:)

contains
  !C***
  !C*** MAT_CON for solver
  !C***
  !C
  subroutine hecmw_mat_con ( hecMESH, hecMAT )

    use hecmw_util
    use hecmw_matrix_contact

    implicit none
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH

    call hecmw_mat_con0  (hecMESH, hecMAT)
    call hecmw_mat_con1  (         hecMAT)
    call hecmw_cmat_init (hecMAT%cmat)

  end subroutine hecmw_mat_con
  !C
  !C***
  !C*** MAT_CON0 for solver
  !C***
  !C
  subroutine hecmw_mat_con0 ( hecMESH, hecMAT )

    use hecmw_util
    use hecmw_etype

    implicit none
    integer(kind=kint) ierr,itype,is,iE,ic_type,nn,icel,j,k,inod
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) nid(20)

    integer(kind=kint), dimension(2048) :: NCOL1, NCOL2
    !C
    !C +-------+
    !C | INIT. |
    !C +-------+
    !C===
    hecMAT%NP= hecMESH%n_node
    hecMAT%N = hecMESH%nn_internal

    NU= 10
    NL= 10

    allocate (INL(hecMAT%NP), IAL(hecMAT%NP,NL))
    allocate (INU(hecMAT%NP), IAU(hecMAT%NP,NU))

    INL= 0
    IAL= 0
    INU= 0
    IAU= 0
    !C===
    !C
    !C +----------------------------------------+
    !C | CONNECTIVITY according to ELEMENT TYPE |
    !C +----------------------------------------+
    !C===
    do
      ierr = 0
      do itype= 1, hecMESH%n_elem_type
        is= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if( hecmw_is_etype_patch(ic_type) ) cycle
        !C Set number of nodes
        nn = hecmw_get_max_node(ic_type)
        !C element loop
        do icel= is, iE
          is= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nid(j)= hecMESH%elem_node_item (is+j)
          enddo
          do j=1,nn
            do k=1,nn
              if( k .ne. j ) then
                call hecmw_FIND_NODE( hecMAT,nid(j),nid(k), ierr  )
                if( ierr.ne.0 ) then
                  call hecmw_mat_con0_clear (ierr)
                  exit
                endif
              endif
            enddo
            if( ierr.ne.0 ) exit
          enddo
          if( ierr.ne.0 ) exit
        enddo
        if( ierr.ne.0 ) exit
      enddo
      if( ierr.eq.0 ) exit
    enddo
    !C===
    !C
    !C +---------+
    !C | SORTING |
    !C +---------+
    !C===
    do inod= 1, hecMAT%NP
      NN= INL(inod)
      do k= 1, NN
        NCOL1(k)= IAL(inod,k)
      enddo
      call hecmw_mSORT (NCOL1, NCOL2, NN)
      do k= NN, 1, -1
        IAL(inod,NN-k+1)= NCOL1(NCOL2(k))
      enddo
      NN= INU(inod)
      do k= 1, NN
        NCOL1(k)= IAU(inod,k)
      enddo
      call hecmw_mSORT (NCOL1, NCOL2, NN)
      do k= NN, 1, -1
        IAU(inod,NN-k+1)= NCOL1(NCOL2(k))
      enddo
    enddo
    !C===
  contains
    !C
    !C*** MAT_CON0_CLEAR
    !C
    subroutine hecmw_mat_con0_clear (IERR)

      implicit none
      integer(kind=kint) IERR

      deallocate (INL, IAL, INU, IAU)

      if (IERR.eq.1) NL= NL + 5
      if (IERR.eq.2) NU= NU + 5
      allocate (INL(hecMAT%NP),IAL(hecMAT%NP,NL))
      allocate (INU(hecMAT%NP),IAU(hecMAT%NP,NU))

      INL= 0
      IAL= 0
      INU= 0
      IAU= 0

    end subroutine hecmw_mat_con0_clear
  end subroutine hecmw_mat_con0
  !C
  !C***
  !C*** FIND_TS_NODE
  !C***
  !C
  subroutine hecmw_FIND_NODE ( hecMAT, ip1,ip2, IERR )

    use hecmw_util

    implicit none
    integer(kind=kint) ip1,ip2,IERR
    integer(kind=kint) kk,icou
    type (hecmwST_matrix)     :: hecMAT

    if (ip1.gt.ip2) then
      do kk= 1, INL(ip1)
        if (ip2.eq.IAL(ip1,kk)) return
      enddo
      icou= INL(ip1) + 1
      if (icou.gt.NL) then
        IERR= 1
        return
      endif
      IAL(ip1,icou)= ip2
      INL(ip1     )= icou
      return
    endif

    if (ip2.gt.ip1) then
      do kk= 1, INU(ip1)
        if (ip2.eq.IAU(ip1,kk)) return
      enddo
      icou= INU(ip1) + 1
      if (icou.gt.NU) then
        IERR= 2
        return
      endif
      IAU(ip1,icou)= ip2
      INU(ip1     )= icou
      return
    endif

  end subroutine hecmw_FIND_NODE
  !C
  !C***
  !C*** fstr_mSORT
  !C***
  !C
  subroutine hecmw_mSORT (STEM,INUM,NN)
    use hecmw_util

    implicit none
    integer(kind=kint) NN
    integer(kind=kint) STEM(NN), INUM(NN)
    integer(kind=kint) ii,jj,ITEM
    do ii = 1,NN
      INUM(ii)= ii
    enddo
    do ii= 1,NN-1
      !CDIR NOVECTOR
      do jj= 1,NN-ii
        if (STEM(INUM(jj)) .lt. STEM(INUM(jj+1))) then
          ITEM      = INUM(jj+1)
          INUM(jj+1)= INUM(jj)
          INUM(jj)  = ITEM
        endif
      enddo
    enddo
    return
  end subroutine hecmw_mSORT
  !C
  !C***
  !C*** MAT_CON1 for solver
  !C***
  !C
  subroutine hecmw_mat_con1 (hecMAT)

    use hecmw_util

    implicit none
    integer(kind=kint) i,k,kk
    type (hecmwST_matrix    ) :: hecMAT

    allocate (hecMAT%indexL(0:hecMAT%NP), hecMAT%indexU(0:hecMAT%NP))

    hecMAT%indexL = 0
    hecMAT%indexU = 0
    do i = 1, hecMAT%NP
      hecMAT%indexL(i) = hecMAT%indexL(i-1) + INL(i)
      hecMAT%indexU(i) = hecMAT%indexU(i-1) + INU(i)
    enddo

    hecMAT%NPL = hecMAT%indexL(hecMAT%NP)
    hecMAT%NPU = hecMAT%indexU(hecMAT%NP)

    allocate (hecMAT%itemL(hecMAT%NPL), hecMAT%itemU(hecMAT%NPU))

    do i = 1, hecMAT%NP
      do k = 1, INL(i)
        kk  = k + hecMAT%indexL(i-1)
        hecMAT%itemL(kk) =     IAL(i,k)
      enddo
      do k= 1, INU(i)
        kk  = k + hecMAT%indexU(i-1)
        hecMAT%itemU(kk) =     IAU(i,k)
      enddo
    enddo

    deallocate (INL, INU, IAL, IAU)

  end subroutine hecmw_mat_con1
end module hecmw_matrix_con
