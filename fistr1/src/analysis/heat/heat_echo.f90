!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> ECHO for HEAT solver
module m_heat_echo
contains

  subroutine heat_echo ( p, hecMESH, fstrHEAT )

    use m_fstr
    use m_hecmw2fstr_mesh_conv

    implicit none

    type(fstr_param)         :: p
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_heat)          :: fstrHEAT

    !** Local variables
    integer(kind=kint) :: i, j, itype, is, iE, ic_type, nn, icel, isect, mid, ic, im, jm, jS, jE, kc, km
    real(kind=kreal)   :: val, temp, aa, bb, time, x, y, z
    integer(kind=kint) :: ig1, iS0, iE0, ik, in, inod, nid
    integer(kind=kint) :: nids(20)


    !C +-------------------------------+
    !C | GLOBAL PARAMETERS             |
    !C +-------------------------------+

    write(ILOG,*) 'global parameters  ***********'
    write(ILOG,*)
    write(ILOG,*) 'IECHO   ',IECHO
    write(ILOG,*) 'IRESULT ',IRESULT
    write(ILOG,*) 'IVISUAL ',IVISUAL
    write(ILOG,*)
    write(ILOG,*) 'for heat ...'
    write(ILOG,*) 'INEUTRAL ', INEUTRAL
    write(ILOG,*) 'IRRES    ', IRRES
    write(ILOG,*) 'IWRES    ', IWRES
    write(ILOG,*) 'NRRES    ', NRRES
    write(ILOG,*) 'NPRINT   ', NPRINT
    write(ILOG,*)
    write(ILOG,*) 'REF_TEMP ', REF_TEMP
    write(ILOG,*)
    write(ILOG,*) 'ANALYSIS CONTROL for HEAT'
    write(ILOG,*) 'DT     ',DT
    write(ILOG,*) 'ETIME  ',ETIME
    write(ILOG,*) 'ITMAX  ',ITMAX
    write(ILOG,*) 'EPS    ',EPS

    !C +-------------------------------+
    !C | fstrPARAM                     |
    !C +-------------------------------+

    write(ILOG,*) 'fstrPARAM  ********************'
    write(ILOG,*)
    write(ILOG,*) 'solution_type ',p%solution_type
    write(ILOG,*) 'solver_method ',p%solver_method
    write(ILOG,*)
    write(ILOG,*) '!!STATIC !HEAT'
    write(ILOG,*) p%analysis_n
    if( associated( P%dtime))  write(ILOG,*) 'dtime  ', p%dtime
    if( associated( P%etime))  write(ILOG,*) 'etime  ', p%etime
    if( associated( P%dtmin))  write(ILOG,*) 'dtmin  ', p%dtmin
    if( associated( P%delmax)) write(ILOG,*) 'delmax ', p%delmax
    if( associated( P%itmax))  write(ILOG,*) 'itmax  ', p%itmax
    if( associated( P%eps))    write(ILOG,*) 'eps    ', p%eps
    write(ILOG,*) 'ref_temp ', p%ref_temp
    write(ILOG,*)
    write(ILOG,*) 'output control'
    write(ILOG,*) 'fg_echo   ', p%fg_echo
    write(ILOG,*) 'fg_result ', p%fg_result
    write(ILOG,*) 'fg_visual ', p%fg_visual
    write(ILOG,*)
    write(ILOG,*) 'for heat ...'
    write(ILOG,*) 'fg_neutral ', p%fg_neutral
    write(ILOG,*) 'fg_irres   ', p%fg_irres
    write(ILOG,*) 'fg_iwres   ', p%fg_iwres
    write(ILOG,*) 'nrres   ',    p%nrres
    write(ILOG,*) 'nprint  ',    p%nprint
    write(ILOG,*)
    write(ILOG,*) 'index table for global node ID sorting'
    write(ILOG,*) 'n_node ', p%n_node
    if( associated( P%global_local_ID)) write(ILOG,*) 'global_local_ID ', p%global_local_ID

    !C +-------------------------------+
    !C | NODE                          |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Nodes'
    write(ILOG,*) '### Number of nodes',hecMESH%n_node
    write(ILOG,*) 'ID X Y Z'
    do i=1,hecMESH%n_node
      nid = hecMESH%global_node_ID(i)
      x = hecMESH%node(3*i-2)
      y = hecMESH%node(3*i-1)
      z = hecMESH%node(3*i)
      write(ILOG,*) nid,x,y,z
    enddo


    !C +-------------------------------+
    !C | ELEMENT                       |
    !C +-------------------------------+

    call fstr2hecmw_mesh_conv( hecMESH )
    write(ILOG,*)
    write(ILOG,*) '### Elements'

    do itype= 1, hecMESH%n_elem_type
      is= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)

      !C** Set number of nodes
      nn = hecmw_get_max_node(ic_type)
      !C element loop
      do icel= is, iE
        !C** node ID
        is= hecMESH%elem_node_index(icel-1)
        do j=1,nn
          if( hecMESH%n_refine > 0 ) then
            nids(j)= hecMESH%elem_node_item (is+j)
          else
            nids(j)= hecMESH%global_node_ID( hecMESH%elem_node_item (is+j))
          endif
        enddo
        !C** section  ID
        isect= hecMESH%section_ID(icel)
        !C** material ID
        mid= hecMESH%section%sect_mat_ID_item(isect)
        write(ILOG,*) '### Element ID=',hecMESH%global_elem_ID(icel)
        write(ILOG,*) ic_type,isect,mid
        write(ILOG,*) (nids(j),j=1,nn)
      enddo
    enddo
    call hecmw2fstr_mesh_conv(hecMESH)
    !C +-------------------------------+
    !C | Material                      |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Material'

    ic = 0
    do im = 1, hecMESH%material%n_mat
      write(ILOG,*)
      write(ILOG,*) '  Material No. =', im
      do jm = 1, 3
        ic = ic  + 1
        jS= hecMESH%material%mat_TABLE_index(ic-1) + 1
        jE= hecMESH%material%mat_TABLE_index(ic  )
        nn= jE - jS + 1
        if (jm.eq.1) write(ILOG,*) ' Density       Temperature   functionA     functionB'
        if (jm.eq.2) write(ILOG,*) ' Spcific heat  Temperature   functionA     functionB'
        if (jm.eq.3) write(ILOG,*) ' Conductivity  Temperature   functionA     functionB'
        kc = 0
        do km = jS, jE
          kc = kc + 1
          val  = hecMESH%material%mat_VAL (km)
          temp = hecMESH%material%mat_TEMP(km)
          if (jm.eq.1) aa = fstrHEAT%RHOfuncA(im,kc)
          if (jm.eq.1) bb = fstrHEAT%RHOfuncB(im,kc)
          if (jm.eq.2) aa = fstrHEAT%CPfuncA(im,kc)
          if (jm.eq.2) bb = fstrHEAT%CPfuncB(im,kc)
          if (jm.eq.3) aa = fstrHEAT%CONDfuncA(im,kc)
          if (jm.eq.3) bb = fstrHEAT%CONDfuncB(im,kc)
          write(ILOG,'(1p4e13.4)') val,temp,aa,bb
        enddo
        if (jm.eq.1) aa = fstrHEAT%RHOfuncA(im,kc+1)
        if (jm.eq.1) bb = fstrHEAT%RHOfuncB(im,kc+1)
        if (jm.eq.2) aa = fstrHEAT%CPfuncA(im,kc+1)
        if (jm.eq.2) bb = fstrHEAT%CPfuncB(im,kc+1)
        if (jm.eq.3) aa = fstrHEAT%CONDfuncA(im,kc+1)
        if (jm.eq.3) bb = fstrHEAT%CONDfuncB(im,kc+1)
        write(ILOG,'(26x,1p2e13.4)') aa,bb
      enddo
    enddo

    !C +-------------------------------+
    !C | NODE GROUP                    |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Ngroup'

    do ig1= 1, hecMESH%node_group%n_grp
      write(ILOG,*)
      write(ILOG,'(a80)') hecMESH%node_group%grp_name(ig1)
      iS0= hecMESH%node_group%grp_index(ig1-1) + 1
      iE0= hecMESH%node_group%grp_index(ig1  )
      do ik= iS0, iE0
        in   = hecMESH%node_group%grp_item(ik)
        write(ILOG,*) in
      enddo
    enddo

    !C +-------------------------------+
    !C | ELEMEN GROUP                  |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Egroup'

    do ig1= 1, hecMESH%elem_group%n_grp
      write(ILOG,*)
      write(ILOG,'(a80)') hecMESH%elem_group%grp_name(ig1)
      iS0= hecMESH%elem_group%grp_index(ig1-1) + 1
      iE0= hecMESH%elem_group%grp_index(ig1  )
      do ik= iS0, iE0
        in   = hecMESH%elem_group%grp_item(ik)
        write(ILOG,*) in
      enddo
    enddo

    !C +-------------------------------+
    !C | BOUNDARY                      |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Boundary'

    write(ILOG,*) '  T_FIX_tot :', fstrHEAT%T_FIX_tot
    write(ILOG,*) '      No./ NID/  amp/ TEMP-BOUNDARY '
    do i = 1, fstrHEAT%T_FIX_tot
      inod = hecMESH%global_node_ID( fstrHEAT%T_FIX_node(i) )
      write(ILOG,'(2i10,i5,1PE15.7)') i, inod, fstrHEAT%T_FIX_ampl(i)        &
        , fstrHEAT%T_FIX_val (i)
    enddo

    write(ILOG,*) '  Q_NOD_tot :', fstrHEAT%Q_NOD_tot
    write(ILOG,*) '      No./ NID/ amp/ Q-POINT '
    do i = 1, fstrHEAT%Q_NOD_tot
      in = hecMESH%global_node_ID( fstrHEAT%Q_NOD_node(i) )
      write(ILOG,'(2i10,i5,1PE15.7)') i, inod, fstrHEAT%Q_NOD_ampl(i)        &
        , fstrHEAT%Q_NOD_val (i)
    enddo

    write(ILOG,*) '  Q_VOL_tot :', fstrHEAT%Q_VOL_tot
    write(ILOG,*) '      No./ EID/ SID/ amp/ Q-VOL '
    do i = 1, fstrHEAT%Q_VOL_tot
      ie = hecMESH%global_elem_ID( fstrHEAT%Q_VOL_elem(i) )
      write(ILOG,'(2i10,i5,1PE15.7)') i, ie, fstrHEAT%Q_VOL_ampl(i)          &
        , fstrHEAT%Q_VOL_val (i)
    enddo

    write(ILOG,*) '  Q_SUF_tot :', fstrHEAT%Q_SUF_tot
    write(ILOG,*) '      No./ EID/ SID/ amp/ Q-SURF '
    do i = 1, fstrHEAT%Q_SUF_tot
      ie = hecMESH%global_elem_ID( fstrHEAT%Q_SUF_elem(i) )
      write(ILOG,'(2i10,2i5,1PE15.7)') i, ie, fstrHEAT%Q_SUF_surf(i)         &
        , fstrHEAT%Q_SUF_ampl(i), fstrHEAT%Q_SUF_val (i)

    enddo

    write(ILOG,*) '  H_SUF_tot :', fstrHEAT%H_SUF_tot
    write(ILOG,*) '      No./ EID/ SID/ H_amp/ T_amp/ HH/ Sink/ '
    do i = 1, fstrHEAT%H_SUF_tot
      ie = hecMESH%global_elem_ID( fstrHEAT%H_SUF_elem(i) )
      write(ILOG,'(2i10,3i5,1P2E15.7)') i, ie, fstrHEAT%H_SUF_surf(i)        &
        , fstrHEAT%H_SUF_ampl(i,1), fstrHEAT%H_SUF_ampl(i,2)       &
        , fstrHEAT%H_SUF_val (i,1), fstrHEAT%H_SUF_val (i,2)
    enddo

    write(ILOG,*) '  R_SUF_tot :', fstrHEAT%R_SUF_tot
    write(ILOG,*) '      No./ EID/ SID/ R_amp/ T_amp/ RR/ Sink/ '
    do i = 1, fstrHEAT%R_SUF_tot
      ie = hecMESH%global_elem_ID( fstrHEAT%R_SUF_elem(i) )
      write(ILOG,'(2i10,3i5,1P2E15.7)') i, ie, fstrHEAT%R_SUF_surf(i)        &
        , fstrHEAT%R_SUF_ampl(i,1), fstrHEAT%R_SUF_ampl(i,2)       &
        , fstrHEAT%R_SUF_val (i,1), fstrHEAT%R_SUF_val (i,2)
    enddo

    !C +-------------------------------+
    !C | Amplitude                     |
    !C +-------------------------------+

    write(ILOG,*)
    write(ILOG,*) '### Amplitude'

    write(ILOG,*) ' AMPLITUDEtot  :', fstrHEAT%AMPLITUDEtot
    do i = 1, fstrHEAT%AMPLITUDEtot
      nn = fstrHEAT%AMPLtab(i)
      write(ILOG,'(2i5,a,a10)') i, nn,' : name=', hecMESH%amp%amp_name(i)

      do j = 1, nn
        time = fstrHEAT%AMPLtime(i,j)
        val  = fstrHEAT%AMPL    (i,j)
        aa   = fstrHEAT%AMPLfuncA(i,j)
        bb   = fstrHEAT%AMPLfuncB(i,j)
        write(ILOG,'(i5,1p4e12.4)') j,time,val,aa,bb
      enddo
      aa   = fstrHEAT%AMPLfuncA(i,nn+1)
      bb   = fstrHEAT%AMPLfuncB(i,nn+1)
      write(ILOG,'(i5,1p4e12.4)') nn+1,time,val,aa,bb

    enddo

  end subroutine heat_echo
end module m_heat_echo
