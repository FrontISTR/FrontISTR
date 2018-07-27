!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine to assemble heat conductance
!! attribution
module m_heat_mat_ass_conductivity
contains
  !C***
  !C*** MAT_ASS_CONDUCTIVITY
  !C***
  subroutine heat_mat_ass_conductivity( hecMESH,hecMAT,fstrHEAT,BETA )

    use m_fstr
    use m_heat_lib

    implicit none
    integer(kind=kint) :: itype, is, iE, ic_type, icel, isect, IMAT, ntab, itab
    integer(kind=kint) :: in0, nn, i, nodLOCAL, jsect, ic, ip, inod, jp, jnod, isU, ieU, ik, isL, ieL
    real(kind=kreal)   :: BETA, TZERO, ALFA, temp, funcA, funcB, XX, YY, ZZ, TT, T0, SS
    real(kind=kreal)   :: asect, thick, GTH, GHH, GR1, GR2
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH

    !!    real (kind=kreal),pointer :: temp(:),funcA(:),funcB(:)
    dimension nodLOCAL(20),XX(20),YY(20),ZZ(20),TT(20),T0(20)                &
      ,SS(400),temp(1000),funcA(1000),funcB(1000)
    !C +-------+
    !C | INIT. |
    !C +-------+

    TZERO = hecMESH%zero_temp
    ALFA = 1.0 - BETA

    call hecmw_mat_clear(hecMAT)
    call hecmw_mat_clear_b(hecMAT)

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

        if( hecMESH%section%sect_type(isect) .ne. 4 ) then
          ntab = fstrHEAT%CONDtab(IMAT)
          do itab = 1, ntab
            temp(itab)  = fstrHEAT%CONDtemp (IMAT,itab)
            funcA(itab) = fstrHEAT%CONDfuncA(IMAT,itab)
            funcB(itab) = fstrHEAT%CONDfuncB(IMAT,itab)
          enddo
          funcA(ntab+1) = fstrHEAT%CONDfuncA(IMAT,ntab+1)
          funcB(ntab+1) = fstrHEAT%CONDfuncB(IMAT,ntab+1)
        endif

        in0 = hecMESH%elem_node_index(icel-1)

        nn = hecmw_get_max_node(ic_type)
        do i = 1, nn
          nodLOCAL(i) = hecMESH%elem_node_item(in0+i)
          XX(i) = hecMESH%node  ( 3*nodLOCAL(i)-2 )
          YY(i) = hecMESH%node  ( 3*nodLOCAL(i)-1 )
          ZZ(i) = hecMESH%node  ( 3*nodLOCAL(i)   )
          TT(i) = fstrHEAT%TEMP (   nodLOCAL(i)   )
          T0(i) = fstrHEAT%TEMP0(   nodLOCAL(i)   )
        enddo
        do i = 1, nn*nn
          SS(i) = 0.0
        enddo

        if( ic_type.eq.111 ) then
          is = hecMESH%section%sect_R_index(isect)
          ASECT = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_111 ( nn,XX,YY,ZZ,TT,IMAT,ASECT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.231 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_231 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.232 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_232 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.241 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_241 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.242 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_242 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.341 ) then
          call heat_THERMAL_341 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.342 ) then
          call heat_THERMAL_342 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.351 ) then
          call heat_THERMAL_351 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.352 ) then
          call heat_THERMAL_352 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.361 ) then
          call heat_THERMAL_361 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.362 ) then
          call heat_THERMAL_362 ( nn,XX,YY,ZZ,TT,IMAT,SS,ntab,temp,funcA,funcB )

        elseif (ic_type.eq.541) then
          jsect = hecMESH%section%sect_R_index(isect-1)+1
          GTH = hecMESH%section%sect_R_item(jsect)
          GHH = hecMESH%section%sect_R_item(jsect+1)
          GR1 = hecMESH%section%sect_R_item(jsect+2)
          GR2 = hecMESH%section%sect_R_item(jsect+3)
          call heat_THERMAL_541 ( nn,XX,YY,ZZ,TT,TZERO,GTH,GHH,GR1,GR2,SS )

        elseif( ic_type.eq.731 ) then
          nn = 4
          nodLOCAL(nn) = hecMESH%elem_node_item(in0+nn-1)
          XX(nn) = XX(nn-1)
          YY(nn) = YY(nn-1)
          ZZ(nn) = ZZ(nn-1)
          TT(nn) = TT(nn-1)
          T0(nn) = T0(nn-1)
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          do i = (nn-1)*(nn-1)+1, nn*nn
            SS(i) = 0.0
          enddo
          call heat_THERMAL_731 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )

        elseif( ic_type.eq.741 ) then
          is = hecMESH%section%sect_R_index(isect)
          THICK = hecMESH%section%sect_R_item(is)
          call heat_THERMAL_741 ( nn,XX,YY,ZZ,TT,IMAT,THICK,SS,ntab,temp,funcA,funcB )
        endif

        ic = 0
        do ip = 1, nn
          inod = nodLOCAL(ip)
          do jp = 1, nn
            ic = ic + 1
            jnod = nodLOCAL(jp)

            if( jnod.gt.inod ) then
              isU = hecMAT%indexU(inod-1) + 1
              ieU = hecMAT%indexU(inod)
              do ik = isU, ieU
                if( hecMAT%itemU(ik).eq.jnod ) hecMAT%AU(ik) = hecMAT%AU(ik) + SS(ic)*BETA
              enddo
            elseif( jnod.lt.inod ) then
              isL = hecMAT%indexL(inod-1) + 1
              ieL = hecMAT%indexL(inod)
              do ik = isL, ieL
                if( hecMAT%itemL(ik).eq.jnod ) hecMAT%AL(ik) = hecMAT%AL(ik) + SS(ic)*BETA
              enddo
            else
              hecMAT%D(inod) = hecMAT%D(inod) + SS(ic)*BETA
            endif

            hecMAT%B(inod) = hecMAT%B(inod) - SS(ic)*T0(jp)*ALFA
          enddo
        enddo
        !C
      enddo
    enddo

  end subroutine heat_mat_ass_conductivity
end module m_heat_mat_ass_conductivity
