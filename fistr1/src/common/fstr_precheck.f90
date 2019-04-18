!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides function to check input data of IFSTR solver
module m_fstr_precheck
contains

  subroutine fstr_get_thickness(hecMESH,mid,thick)
    use hecmw
    use m_fstr
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: mid, ihead
    real(kind=kreal)   :: thick

    ihead = hecMESH%section%sect_R_index(mid-1)
    thick = hecMESH%section%sect_R_item(ihead+1)
    !  if(thick.LE.0.0) then
    !      write(*,*) "Zero thickness <= 0 is illegal"
    !      call hecmw_abort( hecmw_comm_get_comm())
    !  endif
  end subroutine fstr_get_thickness
  !C
  !C***
  !C*** Pre Check for FSTR solver
  !C***
  !C
  subroutine fstr_precheck ( hecMESH, hecMAT, soltype )
    use m_fstr
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix )    :: hecMAT
    integer(kind=kint)        :: soltype

    if(myrank .EQ. 0) then
      write(IMSG,*)
      write(IMSG,*) ' ****   STAGE PreCheck  **'
    endif

    call fstr_precheck_elem ( hecMESH, hecMAT )
    write(IDBG,*) 'fstr_precheck_elem: OK'

    !C
    !C Output sparsity pattern
    !C
    if( soltype == kstNZPROF )then
      call hecmw_nonzero_profile(hecMESH, hecMAT)
    endif

  end subroutine fstr_precheck
  !C
  !C
  subroutine fstr_precheck_elem ( hecMESH, hecMAT )
    use m_fstr
    use m_precheck_LIB_2d
    use m_precheck_LIB_3d
    use m_precheck_LIB_shell

    implicit none

    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID

    !** Local variables
    integer(kind=kint) :: nelem, mid, j, isect, nline, tline, icel, iiS
    integer(kind=kint) :: ndof2, nelem_wo_mpc
    integer(kind=kint) :: ie, ia, jelem, ic_type, nn, jS, jE, itype
    integer(kind=kint) :: nodLOCAL(20),NTOTsum(1)
    integer(kind=8)    :: ntdof2, nonzero
    real(kind=kreal)   :: al, almin, almax, AA, thick, vol, avvol
    real(kind=kreal)   :: tvol, tvmax, tvmin, tlmax, tlmin, asp, aspmax
    real(kind=kreal)   :: xx(20),yy(20),zz(20)
    real(kind=kreal)   :: TOTsum(1),TOTmax(3),TOTmin(2)
    !C
    !C INIT
    !C
    nelem  = 0
    tvol   = 0.0
    tvmax  = 0.0
    tvmin  = 1.0e+20
    tlmax  = 0.0
    tlmin  = 1.0e+20
    aspmax = 0.0

    !C
    !C Mesh Summary
    !C
    write(ILOG,"(a)") '###  Mesh Summary  ###'
    write(ILOG,"(a,i12)") '  Num of node:',hecMESH%n_node
    write(*   ,"(a,i12)") '  Num of node:',hecMESH%n_node
    write(ILOG,"(a,i12)") '  Num of DOF :',hecMESH%n_node*hecMESH%n_dof
    write(*   ,"(a,i12)") '  Num of DOF :',hecMESH%n_node*hecMESH%n_dof
    ndof2  = hecMESH%n_dof**2
    ntdof2 = (hecMESH%n_node*hecMESH%n_dof)**2
    write(ILOG,"(a,i12)") '  Num of elem:',hecMESH%n_elem
    write(*   ,"(a,i12)") '  Num of elem:',hecMESH%n_elem
    nelem_wo_mpc = 0
    do itype = 1, hecMESH%n_elem_type
      jS = hecMESH%elem_type_index(itype-1)
      jE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(.not. (hecmw_is_etype_link(ic_type) .or. hecmw_is_etype_patch(ic_type))) &
           nelem_wo_mpc = nelem_wo_mpc + jE-jS
    enddo
    write(ILOG,"(a,i12)") '   ** w/o MPC/Patch:',nelem_wo_mpc
    write(*   ,"(a,i12)") '   ** w/o MPC/Patch:',nelem_wo_mpc
    do itype = 1, hecMESH%n_elem_type
      jS = hecMESH%elem_type_index(itype-1)
      jE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      write(ILOG,"(a,i4,a,i12)") '  Num of ',ic_type,':',jE-jS
      write(*   ,"(a,i4,a,i12)") '  Num of ',ic_type,':',jE-jS
    enddo
    nonzero = ndof2*(hecMAT%NP + hecMAT%NPU + hecMAT%NPL)
    write(ILOG,"(a,i12)") '  Num of NZ  :',nonzero
    write(*   ,"(a,i12)") '  Num of NZ  :',nonzero
    write(ILOG,"(a,1pe12.5,a)") '  Sparsity   :',100.0d0*dble(nonzero)/dble(ntdof2),"[%]"
    write(*   ,"(a,1pe12.5,a)") '  Sparsity   :',100.0d0*dble(nonzero)/dble(ntdof2),"[%]"

    !C
    !C 3D
    !C
    if( hecMESH%n_dof .eq. 3 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        !          je = hecMESH%elem_ID(ie*2-1)
        jelem = hecMESH%global_elem_ID(ie)
        !C
        ic_type = hecMESH%elem_type(ie)
        !C
        if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_solid(ic_type) &
            & .or. HECMW_is_etype_beam(ic_type) .or. HECMW_is_etype_shell(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        nn = hecmw_get_max_node(ic_type)
        !C
        jS = hecMESH%elem_node_index(ie-1)
        jE = hecMESH%elem_node_index(ie)

        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
          zz(j) = hecMESH%node(3*nodLOCAL(j)  )
        enddo
        !C
        if    ( ic_type.eq.111 ) then
          isect = hecMESH%section_ID(ie)
          mid = hecMESH%section%sect_mat_ID_item(isect)
          call fstr_get_thickness( hecMESH,mid,AA )
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
          nline = 1
          tline = al
          vol = AA*al
          almax = al
          almin = al
        elseif( ic_type.eq.341 ) then
          call PRE_341 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.351 ) then
          call PRE_351 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.361 ) then
          call PRE_361 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.342 ) then
          call PRE_342 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.352 ) then
          call PRE_352 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.362 ) then
          call PRE_362 ( xx,yy,zz,vol,almax,almin )
        elseif( ic_type.eq.641 ) then
          vol = 1.0d-12
        elseif( ic_type.eq.761 ) then
          vol = 1.0d-12
        elseif( ic_type.eq.781 ) then
          vol = 1.0d-12
        endif
        !C
        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo
      !C
      !C 2D
      !C
    elseif( hecMESH%n_dof .eq. 2 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        !          je = hecMESH%elem_ID(ie*2-1)
        jelem = hecMESH%global_elem_ID(ie)
        !C
        ic_type = hecMESH%elem_type(ie)
        !C
        if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_surface(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        nn = hecmw_get_max_node(ic_type)
        !C
        jS = hecMESH%elem_node_index(ie-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
        enddo
        !C
        isect = hecMESH%section_ID(ie)
        mid = hecMESH%section%sect_mat_ID_item(isect)
        call fstr_get_thickness( hecMESH,mid,AA )
        !C
        if    ( ic_type.eq.111 ) then
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2 )
          vol = AA*al
          if( al.gt.tlmax ) tlmax = al
          if( al.lt.tlmin ) tlmin = al
          aspmax = 1.0
        elseif( ic_type.eq.231 ) then
          call PRE_231 ( xx,yy,AA,vol,almax,almin )
        elseif( ic_type.eq.241 ) then
          call PRE_241 ( xx,yy,AA,vol,almax,almin )
        elseif( ic_type.eq.232 ) then
          call PRE_232 ( xx,yy,AA,vol,almax,almin )
        elseif( ic_type.eq.242 ) then
          call PRE_242 ( xx,yy,AA,vol,almax,almin )
        else
          vol = 0.0
        endif
        !C
        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo
      !C
      !C SHELL
      !C
    elseif( hecMESH%n_dof .eq. 6 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        !          je = hecMESH%elem_ID(ie*2-1)
        jelem = hecMESH%global_elem_ID(ie)
        !C
        ic_type = hecMESH%elem_type(ie)
        !C
        if (.not. (hecmw_is_etype_beam(ic_type) .or. hecmw_is_etype_shell(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        nn = hecmw_get_max_node(ic_type)
        !C
        jS = hecMESH%elem_node_index(ie-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
          zz(j) = hecMESH%node(3*nodLOCAL(j)  )
        enddo
        !C
        isect = hecMESH%section_ID(ie)
        mid = hecMESH%section%sect_mat_ID_item(isect)
        call fstr_get_thickness( hecMESH,mid,AA )
        !C
        if    ( ic_type.eq.111 ) then
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
          nline = nline + 1
          tline = tline + al
          vol = AA*al
          if( al.gt.tlmax ) tlmax = al
          if( al.lt.tlmin ) tlmin = al
          aspmax = 1.0
        elseif( ic_type.eq.731 ) then
          call PRE_731 ( xx,yy,zz,AA,vol,almax,almin )
        elseif( ic_type.eq.741 ) then
          call PRE_741 ( xx,yy,zz,AA,vol,almax,almin )
        endif
        !C
        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo
    endif
    !C
    avvol = tvol / nelem
    write(ILOG,*) '###  Sub Summary  ###'
    write(ILOG,*) ' Total Volumes in this region        = ',tvol
    write(ILOG,*) ' Average Volume of elements          = ',avvol
    write(ILOG,*) ' Maximum Volume of elements          = ',tvmax
    write(ILOG,*) ' Minimum Volume of elements          = ',tvmin
    write(ILOG,*) ' Maximum length of element edges     = ',tlmax
    write(ILOG,*) ' Minimum length of element edges     = ',tlmin

    write(ILOG,*) ' Maximum aspect ratio in this region = ',aspmax
    TOTsum(1) = tvol
    call hecmw_allREDUCE_R(hecMESH,TOTsum,1,hecmw_sum)
    NTOTsum(1) = nelem
    call hecmw_allREDUCE_I(hecMESH,NTOTsum,1,hecmw_sum)
    TOTmax(1) = tvmax
    TOTmax(2) = tlmax
    TOTmax(3) = aspmax
    call hecmw_allREDUCE_R(hecMESH,TOTmax,3,hecmw_max)
    TOTmin(1) = tvmin
    TOTmin(2) = tlmin
    call hecmw_allREDUCE_R(hecMESH,TOTmin,2,hecmw_min)
    if( hecMESH%my_rank .eq. 0 ) then
      avvol = TOTsum(1) / NTOTsum(1)
      write(ILOG,*) '###  Global Summary  ###'
      write(ILOG,*) ' TOTAL VOLUME = ',TOTsum(1)
      write(*,*)    ' TOTAL VOLUME = ',TOTsum(1)
      write(ILOG,*) ' AVERAGE VOLUME OF ELEMENTS = ',avvol
      write(*,*)    ' AVERAGE VOLUME OF ELEMENTS = ',avvol
      write(ILOG,*) ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
      write(*,*)    ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
      write(ILOG,*) ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
      write(*,*)    ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
      write(ILOG,*) ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
      write(*,*)    ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
      write(ILOG,*) ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
      write(*,*)    ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
      write(ILOG,*) ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
      write(*,*)    ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
    endif
    !C
  end subroutine fstr_precheck_elem

  subroutine hecmw_nonzero_profile(hecMESH, hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix)     :: hecMAT

    integer(kind=kint) :: i, j, in, jS, jE, ftype, n, ndof, nnz, fio
    real(kind=kreal) :: rnum, dens, cond
    character :: fileid*3

    fio = 70 + hecMESH%my_rank
    write(fileid,"(i3.3)")hecMESH%my_rank

    !ftype = 2: eps
    !ftype = 4: png
    ftype = 4

    n = hecMAT%N
    ndof = 3*n
    nnz = 9*n + 9*2*hecMAT%indexL(n)
    dens = 100*dble(nnz)/dble(9*n*n)
    rnum = (7.21d0+0.01*dlog10(dble(hecMAT%N)))*10.0d0/dble(hecMAT%N)
    cond = 1.0d0
    !rnum = (7.25d0)*10.0d0/dble(hecMAT%N)

    open(fio,file='nonzero.dat.'//fileid, status='replace')
    !write(fio,"(a,f12.5,i0)")"##magic number 10 : 7.2, ",rnum,hecMAT%N
    do i= 1, n
      jS= hecMAT%indexL(i-1) + 1
      jE= hecMAT%indexL(i  )
      write(fio,"(i0,a,i0)")i,"  ",i
      do j= jS, jE
        in = hecMAT%itemL(j)
        write(fio,"(i0,a,i0)")i, "  ",in
        write(fio,"(i0,a,i0)")in,"  ",i
      enddo
    enddo
    close(fio)

    open(fio,file='nonzero.plt.'//fileid, status='replace')
    if(ftype == 4)then
      write(fio,"(a)")'set terminal png size 1500,1500'
    else
      write(fio,"(a)")'set terminal postscript eps enhanced color solid "TimesNewRomanPSMT" 20'
    endif

    write(fio,"(a)")'unset key'
    write(fio,"(a)")'unset xtics'
    write(fio,"(a)")'unset ytics'
    write(fio,"(a)")'set size ratio 1.0'
    write(fio,"(a)")'set border lw 1.0'
    write(fio,"(a,i0,a)")'set xrange[0.5:',n,'.5]'
    write(fio,"(a,i0,a)")'set yrange[0.5:',n,'.5] reverse '

    if(ftype == 4)then
      write(fio,"(a)")'set out "image.'//fileid//'.png"'
    else
      write(fio,"(a)")'set out "image.'//fileid//'.eps"'
      write(fio,"(a)"     )'set label 1 "Name" at graph 1.1,0.9'
      write(fio,"(a)")'set label 2 "N" at graph 1.1,0.85'
      write(fio,"(a)")'set label 3 "Non-Zero Elem." at graph 1.1,0.8'
      write(fio,"(a)")'set label 4 "Density [%]" at graph 1.1,0.75'
      write(fio,"(a)")'set label 9 "Condition Num." at graph 1.1,0.7'
      write(fio,"(a)"     )'set label 5 ":  matrix" at graph 1.4,0.9'
      write(fio,"(a,i0,a)")'set label 6 ":  ',ndof,'" at graph 1.4,0.85'
      write(fio,"(a,i0,a)")'set label 7 ":  ',nnz,'" at graph 1.4,0.8'
      write(fio,"(a,1pe9.2,a)")'set label 8 ": ',dens,'" at graph 1.4,0.75'
      write(fio,"(a,1pe9.2,a)")'set label 10 ": ',cond,'" at graph 1.4,0.7'
    endif

    write(fio,"(a,f12.5,a)")'plot "nonzero.dat.'//fileid//'" pointtype 5 pointsize ',rnum,' linecolor rgb "#F96566"'
    close(fio)

    write(*,*)''
    write(*,*)' ### Command recommendation'
    write(*,*)' gnuplot -persist "nonzero.plt"'

    !call system('gnuplot -persist "nonzero.plt"')
    !open(fio,file='nonzero.dat',status='old')
    !close(fio,status="delete")
  end subroutine hecmw_nonzero_profile
end module m_fstr_precheck
