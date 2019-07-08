!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_heat_mat_ass_conductivity
contains

  subroutine heat_mat_ass_conductivity( hecMESH,hecMAT,fstrHEAT,beta )
    use hecmw
    use m_fstr
    use m_heat_lib
    use m_static_LIB
    implicit none
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: itype, is, iE, ic_type, icel, isect, IMAT, ntab, itab, NDOF
    integer(kind=kint) :: in0, nn, i, in, j, nodLOCAL(20), jsect, ic, ip, inod, jp, jnod, isU, ieU, ik, isL, ieL
    real(kind=kreal)   :: beta, TZERO, ALFA, temp(1000), funcA(1000), funcB(1000), TT(20), T0(20), SS(400)
    real(kind=kreal)   :: asect, thick, GTH, GHH, GR1, GR2
    real(kind=kreal) :: lumped(20), stiff(20, 20), ecoord(3,20)
    real(kind=kreal), allocatable :: S(:)

    NDOF = hecMESH%n_dof
    TZERO = hecMESH%zero_temp
    ALFA = 1.0 - beta

    call hecmw_mat_clear(hecMAT)
    call hecmw_mat_clear_b(hecMAT)

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if (hecmw_is_etype_link(ic_type)) cycle
      if (hecmw_is_etype_patch(ic_type)) cycle

      !$omp parallel default(none), &
        !$omp&  private(icel,isect,IMAT,ntab,itab,temp,funcA,funcB,in0,nn,i,j,nodLOCAL,SS,TT,GTH,GHH,GR1,GR2,ecoord,&
        !$omp&          stiff,in,ASECT,thick,jsect), &
        !$omp&  shared(iS,iE,hecMESH,ic_type,hecMAT,fstrHEAT,TZERO)
      !$omp do
      do icel = iS, iE
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
          TT(i) = fstrHEAT%TEMP (   nodLOCAL(i)   )
          !T0(i) = fstrHEAT%TEMP0(   nodLOCAL(i)   )
          do j = 1, 3
            ecoord(j,i) = hecMESH%node(3*(nodLOCAL(i)-1)+j)
          enddo
        enddo
        do i = 1, nn*nn
          SS(i) = 0.0
        enddo

        if(ic_type == 111) then
          in = hecMESH%section%sect_R_index(isect)
          ASECT = hecMESH%section%sect_R_item(in)
          call heat_conductivity_C1(ic_type, nn, ecoord(1:2,1:nn), TT, IMAT, ASECT, stiff, &
            fstrHEAT%CONDtab(IMAT), fstrHEAT%CONDtemp(IMAT,:), fstrHEAT%CONDfuncA(IMAT,:) ,fstrHEAT%CONDfuncB(IMAT,:))

        elseif(ic_type == 231 .or. ic_type == 232 .or. ic_type == 241 .or. ic_type == 242)then
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          call heat_conductivity_C2(ic_type, nn, ecoord(1:2,1:nn), TT, IMAT, thick, stiff, &
            fstrHEAT%CONDtab(IMAT), fstrHEAT%CONDtemp(IMAT,:), fstrHEAT%CONDfuncA(IMAT,:) ,fstrHEAT%CONDfuncB(IMAT,:))

        elseif(ic_type == 341 .or. ic_type == 342 .or. ic_type == 351 .or. ic_type == 352 .or. &
             & ic_type == 361 .or. ic_type == 362)then
          call heat_conductivity_C3(ic_type, nn, ecoord(1:3,1:nn), TT, IMAT, stiff, &
            fstrHEAT%CONDtab(IMAT), fstrHEAT%CONDtemp(IMAT,:), fstrHEAT%CONDfuncA(IMAT,:) ,fstrHEAT%CONDfuncB(IMAT,:))

        elseif (ic_type == 541) then
          jsect = hecMESH%section%sect_R_index(isect-1)+1
          GTH = hecMESH%section%sect_R_item(jsect)
          GHH = hecMESH%section%sect_R_item(jsect+1)
          GR1 = hecMESH%section%sect_R_item(jsect+2)
          GR2 = hecMESH%section%sect_R_item(jsect+3)
          call heat_conductivity_541(nn, ecoord, TT, TZERO, GTH, GHH, GR1, GR2, SS, stiff)

        elseif(ic_type == 731)then
          nn = 4
          nodLOCAL(nn) = hecMESH%elem_node_item(in0+nn-1)
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          SS = 0.0d0
          call heat_conductivity_shell_731(ic_type, nn, ecoord(1:3,1:nn), TT, IMAT, thick, SS, stiff, &
            fstrHEAT%CONDtab(IMAT), fstrHEAT%CONDtemp(IMAT,:), fstrHEAT%CONDfuncA(IMAT,:) ,fstrHEAT%CONDfuncB(IMAT,:))

        elseif(ic_type == 741)then
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          call heat_conductivity_shell_741(ic_type, nn, ecoord(1:3,1:nn), TT, IMAT, thick, SS, stiff, &
            fstrHEAT%CONDtab(IMAT), fstrHEAT%CONDtemp(IMAT,:), fstrHEAT%CONDfuncA(IMAT,:) ,fstrHEAT%CONDfuncB(IMAT,:))

        else
          write(*,*)"** error setMASS"
        endif

        if(ic_type == 541 .or. ic_type == 731 .or. ic_type == 741)then
          stiff = 0.0d0
          in = 1
          do i = 1, nn
            do j = 1, nn
              stiff(j,i) = SS(in)
              in = in + 1
            enddo
          enddo
        endif

        call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiff)

      enddo
      !$omp end do
      !$omp end parallel
    enddo

    allocate(S(hecMAT%NP))
    S = 0.0d0

    call hecmw_matvec(hecMESH, hecMAT, fstrHEAT%TEMP0, S)

    hecMAT%D  = beta*hecMAT%D
    hecMAT%AU = beta*hecMAT%AU
    hecMAT%AL = beta*hecMAT%AL
    hecMAT%B  = hecMAT%B - ALFA*S

    deallocate(S)

  end subroutine heat_mat_ass_conductivity
end module m_heat_mat_ass_conductivity
