!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine for setting heat radiate
!! boundary conditions
module m_heat_mat_ass_bc_RADIATE
contains
  !C
  !C***
  !C*** MAT_ASS_RADIATE
  !C***
  !C
  subroutine heat_mat_ass_bc_RADIATE( hecMESH, hecMAT, fstrHEAT, CTIME, DTIME, beta )

    use m_fstr
    use m_heat_get_amplitude
    use m_heat_LIB_RADIATE

    implicit none
    integer(kind=kint) :: k, icel, isuf, iam1, iam2, ic_type, isect, nn, is, j, mm, m, ic, ip
    integer(kind=kint) :: inod, jp, jnod, isU, ieU, ik, isL, ieL
    real(kind=kreal)   :: CTIME, DTIME, TZERO, QQ, RR, SINK, thick, beta
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH

    real(kind=kreal)   :: xx(20), yy(20), zz(20), tt(20)
    real(kind=kreal)   :: term1(64), term2(8)
    integer(kind=kint) :: nodLocal(20), nsuf(8)

    TZERO = hecMESH%zero_temp
    !C
    do k = 1, fstrHEAT%R_SUF_tot
      icel    = fstrHEAT%R_SUF_elem(k)
      isuf    = fstrHEAT%R_SUF_surf(k)
      iam1    = fstrHEAT%R_SUF_ampl(k,1)
      iam2    = fstrHEAT%R_SUF_ampl(k,2)
      call heat_get_amplitude ( fstrHEAT, iam1, CTIME, QQ )
      RR      = fstrHEAT%R_SUF_val (k,1) * QQ
      call heat_get_amplitude ( fstrHEAT, iam2, CTIME, QQ )
      SINK    = fstrHEAT%R_SUF_val (k,2) * QQ
      ic_type = hecMESH%elem_type(icel)
      isect   = hecMESH%section_ID(icel)
      !C**
      nn = hecmw_get_max_node(ic_type)
      !C**
      is = hecMESH%elem_node_index(icel-1)
      do j = 1, nn
        nodLOCAL(j) = hecMESH%elem_node_item(is+j)
        xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
        yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
        zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
        tt(j) = fstrHEAT%TEMP( nodLOCAL(j) )
      enddo
      !C**
      if    ( ic_type.eq.231 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm=2
        call heat_RADIATE_231(nn,xx,yy,zz,thick,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.232 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm=3
        call heat_RADIATE_232(nn,xx,yy,zz,thick,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.241 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm=2
        call heat_RADIATE_241(nn,xx,yy,zz,thick,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.242 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm=3
        call heat_RADIATE_242(nn,xx,yy,zz,thick,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.341 ) then
        mm=3
        call heat_RADIATE_341(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.342 ) then
        mm=6
        call heat_RADIATE_342(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.351 ) then
        mm=4
        if( isuf.eq.1 .or. isuf.eq.2 ) mm=3
        call heat_RADIATE_351(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.352 ) then
        mm=8
        if( isuf.eq.1 .or. isuf.eq.2 ) mm=6
        call heat_RADIATE_352(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.361 ) then
        mm=4
        call heat_RADIATE_361(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.362 ) then
        mm=8
        call heat_RADIATE_362(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,mm,term1,term2,nsuf)
      elseif( ic_type.eq.731 ) then
        call heat_RADIATE_731(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,term1,term2)
        mm=3
        do m = 1, mm
          nsuf(m) = m
        enddo
      elseif( ic_type.eq.741 ) then
        call heat_RADIATE_741(nn,xx,yy,zz,tt,isuf,RR,SINK,TZERO,term1,term2)
        mm=4
        do m = 1, mm
          nsuf(m) = m
        enddo
      endif
      !C
      ic = 0
      do ip = 1, mm
        inod = nodLOCAL(nsuf(ip))
        do jp = 1, mm
          jnod = nodLOCAL(nsuf(jp))
          ic = ic + 1
          if( jnod.gt.inod ) then
            isU = hecMAT%indexU(inod-1) + 1
            ieU = hecMAT%indexU(inod)
            do ik = isU, ieU
              if( hecMAT%itemU(ik).eq.jnod ) hecMAT%AU(ik) = hecMAT%AU(ik) - term1(ic)
            enddo
          elseif( jnod.lt.inod ) then
            isL = hecMAT%indexL(inod-1) + 1
            ieL = hecMAT%indexL(inod)
            do ik = isL, ieL
              if( hecMAT%itemL(ik).eq.jnod ) hecMAT%AL(ik) = hecMAT%AL(ik) - term1(ic)
            enddo
          else
            hecMAT%D(inod) = hecMAT%D(inod) - term1(ic)
            hecMAT%B(jnod) = hecMAT%B(jnod) - term2(jp)
          endif
        enddo
      enddo
      !C
    enddo

  end subroutine heat_mat_ass_bc_RADIATE
end module m_heat_mat_ass_bc_RADIATE
