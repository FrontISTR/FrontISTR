!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine to assemble heat capacity matrix
module m_heat_mat_ass_capacity
contains

  subroutine heat_mat_ass_capacity(hecMESH, hecMAT, fstrHEAT, delta_time)
    use m_fstr
    use m_heat_LIB_CAPACITY
    implicit none
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, in, j, nodLOCAL(20), ip, inod
    integer(kind=kint) :: itype, iS, iE, ic_type, icel, isect, IMAT, in0, nn, NDOF
    real(kind=kreal) :: temp(20), lumped(120), mass(20*6, 20*6), ecoord(3,20)
    real(kind=kreal) :: delta_time, surf, THICK, ALFA, BETA

    NDOF = hecMESH%n_dof
    beta = fstrHEAT%beta

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)

      if (hecmw_is_etype_link(ic_type)) cycle
      if (hecmw_is_etype_patch(ic_type)) cycle
      if(ic_type == 3414) cycle

      !$omp parallel default(none), &
        !$omp&  private(icel,isect,IMAT,nn,temp,in0,i,j,nodLOCAL,ecoord,in,thick,surf,inod,lumped,mass), &
        !$omp&  shared(iS,iE,hecMESH,ic_type,hecMAT,fstrHEAT,ndof,delta_time)
      !$omp do
      do icel = iS, iE
        isect = hecMESH%section_ID(icel)
        IMAT = hecMESH%section%sect_mat_ID_item(isect)
        in0 = hecMESH%elem_node_index(icel-1)
        nn = hecmw_get_max_node(ic_type)
        do i = 1, nn
          nodLOCAL(i) = hecMESH%elem_node_item(in0+i)
          do j = 1, NDOF
            ecoord(j,i) = hecMESH%node(3*(nodLOCAL(i)-1)+j)
          enddo
          temp(i) = fstrHEAT%TEMP0(nodLOCAL(i))
        enddo

        lumped = 0.0d0
        if(ic_type == 231 .or. ic_type == 232 .or. ic_type == 241 .or. ic_type == 242)then
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          call heat_capacity_C2(ic_type, nn, ecoord(1:2,1:nn), temp, IMAT, thick, lumped, mass, &
            fstrHEAT%CPtab(IMAT), fstrHEAT%CPtemp(IMAT,:), fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:), &
            fstrHEAT%RHOtab(IMAT), fstrHEAT%RHOtemp(IMAT,:), fstrHEAT%RHOfuncA(IMAT,:), fstrHEAT%RHOfuncB(IMAT,:))

        elseif(ic_type == 341 .or. ic_type == 342 .or. ic_type == 351 .or. ic_type == 352 .or. &
             & ic_type == 361 .or. ic_type == 362)then
          call heat_capacity_C3(ic_type, nn, ecoord(1:3,1:nn), temp, IMAT, lumped, mass, &
            fstrHEAT%CPtab(IMAT), fstrHEAT%CPtemp(IMAT,:), fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:), &
            fstrHEAT%RHOtab(IMAT), fstrHEAT%RHOtemp(IMAT,:), fstrHEAT%RHOfuncA(IMAT,:), fstrHEAT%RHOfuncB(IMAT,:))

        elseif(ic_type == 731)then
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          call heat_capacity_shell_731(ic_type, nn, ecoord(1:3,1:nn), temp, IMAT, thick, lumped, mass, &
            fstrHEAT%CPtab(IMAT), fstrHEAT%CPtemp(IMAT,:), fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:), &
            fstrHEAT%RHOtab(IMAT), fstrHEAT%RHOtemp(IMAT,:), fstrHEAT%RHOfuncA(IMAT,:), fstrHEAT%RHOfuncB(IMAT,:))

        elseif(ic_type == 741)then
          in = hecMesh%section%sect_R_index(isect)
          thick = hecMESH%section%sect_R_item(in)
          call heat_capacity_shell_741(ic_type, nn, ecoord(1:3,1:nn), temp, IMAT, thick, lumped, mass, &
            fstrHEAT%CPtab(IMAT), fstrHEAT%CPtemp(IMAT,:), fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:), &
            fstrHEAT%RHOtab(IMAT), fstrHEAT%RHOtemp(IMAT,:), fstrHEAT%RHOfuncA(IMAT,:), fstrHEAT%RHOfuncB(IMAT,:))

        elseif(ic_type == 111 .or. ic_type == 301 .or. ic_type == 611 .or. ic_type == 641)then
          in = hecMesh%section%sect_R_index(isect)
          surf = hecMESH%section%sect_R_item(in)
          call heat_capacity_C1(ic_type, nn, ecoord(1:3,1:nn), temp, IMAT, surf, lumped, mass, &
            fstrHEAT%CPtab(IMAT), fstrHEAT%CPtemp(IMAT,:), fstrHEAT%CPfuncA(IMAT,:) ,fstrHEAT%CPfuncB(IMAT,:), &
            fstrHEAT%RHOtab(IMAT), fstrHEAT%RHOtemp(IMAT,:), fstrHEAT%RHOfuncA(IMAT,:), fstrHEAT%RHOfuncB(IMAT,:))

        else
          write(*,*)"** error setMASS"
        endif

        do ip = 1, nn
          inod = nodLOCAL(ip)
          !$omp atomic
          hecMAT%D(inod) = hecMAT%D(inod) + lumped(ip) / delta_time
          !$omp atomic
          hecMAT%B(inod) = hecMAT%B(inod) + lumped(ip)*temp(ip) / delta_time
        enddo
      enddo
      !$omp end do
      !$omp end parallel
    enddo
  end subroutine heat_mat_ass_capacity
end module m_heat_mat_ass_capacity
