!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides functions to initialize heat analysis
module m_heat_init
contains

  subroutine heat_init(hecMESH, fstrHEAT)
    use m_fstr
    implicit none
    type(fstr_heat) :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j

    allocate(fstrHEAT%TEMP0(hecMESH%n_node))
    allocate(fstrHEAT%TEMPC(hecMESH%n_node))
    allocate(fstrHEAT%TEMP (hecMESH%n_node))
    fstrHEAT%TEMP0 = 0.0d0
    fstrHEAT%TEMPC = 0.0d0
    fstrHEAT%TEMP  = 0.0d0

    if(hecMESH%hecmw_flag_initcon == 1)then
      do i = 1, hecMESH%n_node
        j = hecMESH%node_init_val_index(i)
        fstrHEAT%TEMP0(i) = hecMESH%node_init_val_item(j)
        fstrHEAT%TEMPC(i) = fstrHEAT%TEMP0(i)
        fstrHEAT%TEMP (i) = fstrHEAT%TEMP0(i)
      enddo
      write(ILOG,*) ' Initial condition of temperatures: OK'
    endif
    if( associated(g_InitialCnd) ) then
        do j=1,size(g_InitialCnd)
          if( g_InitialCnd(j)%cond_name=="temperature" ) then
            do i= 1, hecMESH%n_node
              fstrHEAT%TEMP0(i)= g_InitialCnd(j)%realval(i)
              fstrHEAT%TEMPC(i)= fstrHEAT%TEMP0(i)
              fstrHEAT%TEMP (i)= fstrHEAT%TEMP0(i)
            enddo
            exit
          endif
        enddo
        write(ILOG,*) ' Initial condition of temperatures: OK'
    endif
  end subroutine heat_init

  subroutine heat_init_log(hecMESH)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh) :: hecMESH

    if(hecMESH%my_rank == 0)then
      write(IMSG,*) '============================='
      write(IMSG,*) '  H E A T   T R A N S F E R  '
      write(IMSG,*) '============================='
      write(ISTA,*)
      write(ISTA,*)'  ISTEP    INCR    ITER     RESIDUAL     IITER   '
      write(ISTA,*)'-------------------------------------------------'
    endif
  end subroutine heat_init_log

  subroutine heat_finalize(fstrHEAT)
    use m_fstr
    implicit none
    type(fstr_heat) :: fstrHEAT
    deallocate(fstrHEAT%TEMP0)
    deallocate(fstrHEAT%TEMPC)
    deallocate(fstrHEAT%TEMP )
  end subroutine heat_finalize

  !C***
  !C*** INIT_AMPLITUDE
  !C***
  subroutine heat_init_amplitude (hecMESH, fstrHEAT)

    use m_fstr

    implicit none
    integer(kind=kint) :: namax, i, nn, is, iE, icou, j, k
    real(kind=kreal)   :: x1, y1, x2, y2
    type(fstr_heat)    :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH
    !C
    !C===
    namax = 0
    do i = 1, hecMESH%amp%n_amp
      nn = hecMESH%amp%amp_index(i) - hecMESH%amp%amp_index(i-1)
      namax = max(nn,namax)
    enddo

    fstrHEAT%AMPLITUDEtot= hecMESH%amp%n_amp

    allocate (fstrHEAT%AMPLtab (fstrHEAT%AMPLITUDEtot) )
    allocate (fstrHEAT%AMPL    (fstrHEAT%AMPLITUDEtot,namax),                &
      fstrHEAT%AMPLtime(fstrHEAT%AMPLITUDEtot,namax) )

    fstrHEAT%AMPLtab  = 0
    fstrHEAT%AMPL     = 0.d0
    fstrHEAT%AMPLtime = 0.d0

    do i = 1, fstrHEAT%AMPLITUDEtot
      is = hecMESH%amp%amp_index(i-1) + 1
      iE = hecMESH%amp%amp_index(i)
      nn = iE - is + 1

      fstrHEAT%AMPLtab(i) = nn

      icou = 0
      do j = is, iE
        icou = icou + 1
        fstrHEAT%AMPL    (i,icou) = hecMESH%amp%amp_val  (j)
        fstrHEAT%AMPLtime(i,icou) = hecMESH%amp%amp_table(j)
      enddo
    enddo
    !C===

    !C
    !C +-----------+
    !C | AMP-TABLE |
    !C +-----------+
    !C===
    allocate ( fstrHEAT%AMPLfuncA( fstrHEAT%AMPLITUDEtot,namax+1 ) )
    allocate ( fstrHEAT%AMPLfuncB( fstrHEAT%AMPLITUDEtot,namax+1 ) )

    fstrHEAT%AMPLfuncA  = 0.d0
    fstrHEAT%AMPLfuncB  = 0.d0
    !C
    !C--
    do i = 1, fstrHEAT%AMPLITUDEtot

      fstrHEAT%AMPLfuncA(i,1) = 0.d0
      fstrHEAT%AMPLfuncB(i,1) = fstrHEAT%AMPL(i,1)

      nn = fstrHEAT%AMPLtab(i)

      do k = 2, nn

        x1 = fstrHEAT%AMPLtime(i,k-1)
        y1 = fstrHEAT%AMPL    (i,k-1)
        x2 = fstrHEAT%AMPLtime(i,k)
        y2 = fstrHEAT%AMPL    (i,k)

        fstrHEAT%AMPLfuncA(i,k) =  (y2-y1)/(x2-x1)
        fstrHEAT%AMPLfuncB(i,k) = -(y2-y1)/(x2-x1)*x1 + y1

      enddo

      fstrHEAT%AMPLfuncA(i,nn+1) = 0.d0
      fstrHEAT%AMPLfuncB(i,nn+1) = fstrHEAT%AMPL(i,nn)

    enddo

    !C===
  end subroutine heat_init_amplitude
  !C***
  !C*** INIT_MATERIAL
  !C***
  subroutine heat_init_material (hecMESH, fstrHEAT)

    use m_fstr

    implicit none
    integer(kind=kint) :: m1max, m2max, m3max, icou, im, jm, nn, ic, jS, jE, kc, km, k
    real(kind=kreal)   :: aa, bb
    type(fstr_heat)    :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C +----------+
    !C | MATERIAL |
    !C +----------+
    !C===
    fstrHEAT%MATERIALtot= hecMESH%material%n_mat

    m1max = 0
    m2max = 0
    m3max = 0
    icou  = 0
    do im = 1, hecMESH%material%n_mat
      do jm = 1, 3
        icou = icou + 1
        nn = hecMESH%material%mat_TABLE_index(icou) - hecMESH%material%mat_TABLE_index(icou-1)
        if( jm.eq.1 ) m1max = max(nn,m1max)
        if( jm.eq.2 ) m2max = max(nn,m2max)
        if( jm.eq.3 ) m3max = max(nn,m3max)
      enddo
    enddo

    allocate (fstrHEAT%RHOtab  (fstrHEAT%MATERIALtot),               &
      fstrHEAT%CPtab   (fstrHEAT%MATERIALtot),                       &
      fstrHEAT%CONDtab (fstrHEAT%MATERIALtot))
    allocate (fstrHEAT%RHO     (fstrHEAT%MATERIALtot,m1max),         &
      fstrHEAT%RHOtemp (fstrHEAT%MATERIALtot,m1max))
    allocate (fstrHEAT%CP      (fstrHEAT%MATERIALtot,m2max),         &
      fstrHEAT%CPtemp  (fstrHEAT%MATERIALtot,m2max))
    allocate (fstrHEAT%COND    (fstrHEAT%MATERIALtot,m3max),         &
      fstrHEAT%CONDtemp(fstrHEAT%MATERIALtot,m3max))

    fstrHEAT%RHO  = 0.d0
    fstrHEAT%CP   = 0.d0
    fstrHEAT%COND = 0.d0

    fstrHEAT%RHOtemp  = 0.d0
    fstrHEAT%CPtemp   = 0.d0
    fstrHEAT%CONDtemp = 0.d0
    fstrHEAT%RHOtab   = 0
    fstrHEAT%CPtab    = 0
    fstrHEAT%CONDtab  = 0

    ic = 0
    do im = 1, fstrHEAT%MATERIALtot
      do jm = 1, 3
        ic = ic + 1
        jS = hecMESH%material%mat_TABLE_index(ic-1) + 1
        jE = hecMESH%material%mat_TABLE_index(ic  )
        nn = jE - jS + 1
        if( jm.eq.1 ) fstrHEAT%RHOtab (im) = nn
        if( jm.eq.2 ) fstrHEAT%CPtab  (im) = nn
        if( jm.eq.3 ) fstrHEAT%CONDtab(im) = nn

        kc = 0
        do km = jS, jE
          kc = kc + 1
          if( jm.eq.1 ) then
            fstrHEAT%RHO     (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%RHOtemp (im,kc) = hecMESH%material%mat_TEMP(km)
          endif
          if( jm.eq.2 ) then
            fstrHEAT%CP      (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%CPtemp  (im,kc) = hecMESH%material%mat_TEMP(km)
          endif
          if( jm.eq.3 ) then
            fstrHEAT%COND    (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%CONDtemp(im,kc) = hecMESH%material%mat_TEMP(km)
          endif
        enddo
      enddo
    enddo
    !C===

    !C
    !C +-----------+
    !C | MAT-TABLE |
    !C +-----------+
    !C===
    allocate (fstrHEAT%RHOfuncA (fstrHEAT%MATERIALtot, m1max+1)              &
      ,fstrHEAT%RHOfuncB (fstrHEAT%MATERIALtot, m1max+1))
    allocate (fstrHEAT%CPfuncA  (fstrHEAT%MATERIALtot, m2max+1)              &
      ,fstrHEAT%CPfuncB  (fstrHEAT%MATERIALtot, m2max+1))
    allocate (fstrHEAT%CONDfuncA(fstrHEAT%MATERIALtot, m3max+1)              &
      ,fstrHEAT%CONDfuncB(fstrHEAT%MATERIALtot, m3max+1))

    fstrHEAT%RHOfuncA  = 0.d0
    fstrHEAT%RHOfuncB  = 0.d0
    fstrHEAT%CPfuncA   = 0.d0
    fstrHEAT%CPfuncB   = 0.d0
    fstrHEAT%CONDfuncA = 0.d0
    fstrHEAT%CONDfuncB = 0.d0
    !C
    !C--RHO
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%RHOfuncB(im,1) = fstrHEAT%RHO(im,1)
      do k = 2, fstrHEAT%RHOtab(im)
        bb= fstrHEAT%RHO    (im,k) - fstrHEAT%RHO    (im,k-1)
        aa= fstrHEAT%RHOtemp(im,k) - fstrHEAT%RHOtemp(im,k-1)
        fstrHEAT%RHOfuncA(im,k) =   bb/aa
        fstrHEAT%RHOfuncB(im,k) = -(bb/aa)*fstrHEAT%RHOtemp(im,k-1) + fstrHEAT%RHO(im,k-1)
      enddo
      fstrHEAT%RHOfuncB(im,fstrHEAT%RHOtab(im)+1) = fstrHEAT%RHO(im,fstrHEAT%RHOtab(im))
    enddo
    !C
    !C-- CP
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%CPfuncB(im,1) = fstrHEAT%CP(im,1)
      do k = 2, fstrHEAT%CPtab(im)
        bb= fstrHEAT%CP    (im,k) - fstrHEAT%CP    (im,k-1)
        aa= fstrHEAT%CPtemp(im,k) - fstrHEAT%CPtemp(im,k-1)
        fstrHEAT%CPfuncA(im,k) =   bb/aa
        fstrHEAT%CPfuncB(im,k) = -(bb/aa)*fstrHEAT%CPtemp(im,k-1) + fstrHEAT%CP(im,k-1)
      enddo
      fstrHEAT%CPfuncB(im,fstrHEAT%CPtab(im)+1) =  fstrHEAT%CP(im,fstrHEAT%CPtab(im))
    enddo
    !C
    !C-- COND.
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%CONDfuncB(im,1)= fstrHEAT%COND(im,1)
      do k = 2, fstrHEAT%CONDtab(im)
        bb = fstrHEAT%COND    (im,k) - fstrHEAT%COND    (im,k-1)
        aa = fstrHEAT%CONDtemp(im,k) - fstrHEAT%CONDtemp(im,k-1)
        fstrHEAT%CONDfuncA(im,k) =   bb/aa
        fstrHEAT%CONDfuncB(im,k) = -(bb/aa)*fstrHEAT%CONDtemp(im,k-1) + fstrHEAT%COND(im,k-1)
      enddo
      fstrHEAT%CONDfuncB(im,fstrHEAT%CONDtab(im)+1) = fstrHEAT%COND(im,fstrHEAT%CONDtab(im))
    enddo
    !C===
  end subroutine heat_init_material

end module m_heat_init
