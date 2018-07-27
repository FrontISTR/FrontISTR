!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine to assemble heat capacity matrix
module m_heat_mat_ass_capacity
contains
  !C***
  !C*** MAT_ASS_CAPACITY
  !C***
  subroutine heat_mat_ass_capacity ( hecMESH,hecMAT,fstrHEAT,DTIME )

    use m_fstr
    use m_heat_LIB_CAPACITY

    implicit none
    integer(kind=kint) :: itype, is, iE, ic_type, icel, isect, IMAT, in0, nn, i, nodLOCAL, ip, inod
    real(kind=kreal)   :: DTIME, XX, YY, ZZ, TT, ASECT, S0, THICK

    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH

    dimension nodLOCAL(20), XX(20), YY(20), ZZ(20), TT(20), S0(20)


    !C +-------------------------------+
    !C | ELEMENT-by-ELEMENT ASSEMBLING |
    !C | according to ELEMENT TYPE     |
    !C +-------------------------------+

    do itype= 1, hecMESH%n_elem_type
      is= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if (hecmw_is_etype_link(ic_type)) cycle

      do icel = is, iE
        isect = hecMESH%section_ID(icel)
        IMAT = hecMESH%section%sect_mat_ID_item(isect)
        in0 = hecMESH%elem_node_index(icel-1)

        nn = hecmw_get_max_node(ic_type)
        do i = 1, nn
          nodLOCAL(i) = hecMESH%elem_node_item(in0+i)
          XX(i) = hecMESH%node ( 3*nodLOCAL(i)-2 )
          YY(i) = hecMESH%node ( 3*nodLOCAL(i)-1 )
          ZZ(i) = hecMESH%node ( 3*nodLOCAL(i)   )
          TT(i) = fstrHEAT%TEMP0(   nodLOCAL(i)   )
        enddo
        do i = 1, nn
          S0(i) = 0.0
        enddo

        if( ic_type.eq.111 ) then
          is = hecMesh%section%sect_R_index(isect)
          ASECT = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_111 ( nn,XX,YY,ZZ,TT,IMAT,ASECT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.231 ) then
          is = hecMesh%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_231 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.232 ) then
          is = hecMesh%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_232 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.241 ) then
          is = hecMesh%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_241 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.242 ) then
          is = hecMesh%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_242 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.341 ) then
          call heat_CAPACITY_341 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.342 ) then
          call heat_CAPACITY_342 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.351 ) then
          call heat_CAPACITY_351 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.352 ) then
          call heat_CAPACITY_352 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.361 ) then
          call heat_CAPACITY_361 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.362 ) then
          call heat_CAPACITY_362 ( nn,XX,YY,ZZ,TT,IMAT,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.731 ) then
          is = hecMesh%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_731 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        elseif( ic_type.eq.741 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_CAPACITY_741 ( nn,XX,YY,ZZ,TT,IMAT,THICK,S0 &
            ,fstrHEAT%CPtab(IMAT)     ,fstrHEAT%CPtemp(IMAT,:) &
            ,fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:) &
            ,fstrHEAT%RHOtab(IMAT)    ,fstrHEAT%RHOtemp(IMAT,:) &
            ,fstrHEAT%RHOfuncA(IMAT,:),fstrHEAT%RHOfuncB(IMAT,:) )
        endif
        !C
        do ip = 1, nn
          inod = nodLOCAL(ip)
          hecMAT%D(inod) = hecMAT%D(inod) + S0(ip) / DTIME
          hecMAT%B(inod) = hecMAT%B(inod) + S0(ip) * TT(ip) / DTIME
        enddo
        !C
      enddo
    enddo

  end subroutine heat_mat_ass_capacity
end module m_heat_mat_ass_capacity
