!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Set up lumped mass matrix
module m_fstr_EIG_setMASS

  use m_eigen_lib
  use m_static_get_prop
  use m_dynamic_mass
  use m_fstr

contains

  subroutine setMASS(fstrSOLID, hecMESH, hecMAT, fstrEIG)
    use hecmw_util
    use m_fstr
    use mMaterial
    use m_static_LIB_shell
    implicit none
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_eigen)         :: fstrEIG
    integer(kind=kint) :: i, iS, iE, ii, nn, j, jn, jS, k
    integer(kind=kint) :: N, NP, NDOF
    integer(kind=kint) :: icel, ic_type, itype, isect, ihead, sec_opt, cid
    integer(kind=kint) :: nodLOCAL(20)
    real(kind=kreal) :: val, surf, chkmass
    real(kind=kreal) :: rho, thick, length
    real(kind=kreal) :: lumped(120), mass(20*6, 20*6), ecoord(3,20)

    if(myrank == 0)then
      write(IMSG,*)"* mass matrix generation"
    endif

    N  = hecMAT%N
    NP = hecMAT%NP
    NDOF = hecMESH%n_dof

    allocate(fstrEIG%mass(NP*NDOF))
    fstrEIG%mass = 0.0d0
    val    = 0.0d0

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      nn = hecmw_get_max_node(ic_type)

      if(hecmw_is_etype_link(ic_type)) cycle
      if(hecmw_is_etype_patch(ic_type)) cycle
      if(ic_type == 3414) cycle

      do icel = iS, iE
        jS = hecMESH%elem_node_index(icel-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
          do i = 1, 3
            ecoord(i,j) = hecMESH%node(3*(nodLOCAL(j)-1)+i)
          enddo
        enddo

        isect = hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        cid   = hecMESH%section%sect_mat_ID_item(isect)
        sec_opt = hecMESH%section%sect_opt(isect)
        rho   = fstrSOLID%materials(cid)%variables(M_DENSITY)
        thick = fstrSOLID%materials(cid)%variables(M_THICK)

        lumped = 0.0d0
        if(ic_type == 231 .or. ic_type == 232 .or. ic_type == 241 .or. ic_type == 242)then
          call mass_C2(ic_type, nn, ecoord(1:2,1:nn), fstrSOLID%elements(icel)%gausses, sec_opt, thick, mass, lumped)

        elseif(ic_type == 341 .or. ic_type == 342 .or. ic_type == 351 .or. ic_type == 352 .or. &
             & ic_type == 361 .or. ic_type == 362 )then
          call mass_C3(ic_type, nn, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses, mass, lumped)

        elseif(ic_type==731 .or. ic_type==741 .or. ic_type==743) then
          rho = fstrSOLID%materials(cid)%variables(M_DENSITY)
          thick = fstrSOLID%materials(cid)%variables(M_THICK)
          call mass_shell(ic_type, nn, ecoord(1:3,1:nn), rho, thick, fstrSOLID%elements(icel)%gausses, mass, lumped)

        elseif(ic_type == 761)then
          surf = get_face3(ecoord(1:3,1:nn))
          rho = fstrSOLID%materials(cid)%variables(M_DENSITY)
          thick = fstrSOLID%materials(cid)%variables(M_THICK)
          val = surf*thick*rho/3.0d0

        elseif(ic_type == 781)then
          surf = get_face4(ecoord(1:3,1:nn))
          rho = fstrSOLID%materials(cid)%variables(M_DENSITY)
          thick = fstrSOLID%materials(cid)%variables(M_THICK)
          val = surf*thick*rho/4.0d0

        elseif(ic_type == 611 .or. ic_type == 641)then
          surf = hecMESH%section%sect_R_item(ihead+4)
          length = get_length(ecoord(1:3,1:nn))
          rho = fstrSOLID%materials(cid)%variables(M_DENSITY)
          val = 0.5d0*surf*length*rho

        elseif(ic_type == 301)then
          surf = hecMESH%section%sect_R_item(ihead+1)
          length = get_length(ecoord(1:3,1:nn))
          rho = fstrSOLID%materials(cid)%variables(M_DENSITY)
          val = 0.5d0*surf*length*rho

        elseif( ic_type/100 == 5 )then !skip interface element
        else
          write(*,*)"** error setMASS"
        endif

        do j = 1,nn
          jn = nodLOCAL(j)
          js = NDOF*(jn-1)

          if(ic_type == 611)then
            fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + val
            fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + val
            fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + val
            fstrEIG%mass(js+4) = fstrEIG%mass(js+4) + 0.0d0
            fstrEIG%mass(js+5) = fstrEIG%mass(js+5) + 0.0d0
            fstrEIG%mass(js+6) = fstrEIG%mass(js+6) + 0.0d0

          elseif(ic_type == 641)then
            if(j == 1 .or. j == 2)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + val
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + val
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + val
            elseif(j == 3 .or. j == 4)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + 0.0d0
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + 0.0d0
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + 0.0d0
            end if

          elseif(ic_type == 761)then
            if(j == 1 .or. j == 2 .or. j == 3)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + val
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + val
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + val
            elseif(j == 4 .or. j == 5 .or. j == 6)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + 0.0d0
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + 0.0d0
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + 0.0d0
            endif

          else if( ic_type == 781 ) then
            if(j == 1 .or. j == 2 .or. j == 3 .or. j == 4)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + val
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + val
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + val
            elseif(j == 5 .or. j == 6 .or. j == 7 .or. j == 8)then
              fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + 0.0d0
              fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + 0.0d0
              fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + 0.0d0
            endif

          elseif(ic_type == 301) then
            fstrEIG%mass(js+1) = fstrEIG%mass(js+1) + val
            fstrEIG%mass(js+2) = fstrEIG%mass(js+2) + val
            fstrEIG%mass(js+3) = fstrEIG%mass(js+3) + val

          else
            do k = 1, NDOF
              fstrEIG%mass(js+k) = fstrEIG%mass(js+k) + lumped(NDOF*(j-1)+k)
            enddo
          endif
        enddo
      enddo
    enddo

    call hecmw_update_R(hecMESH, fstrEIG%mass, NP, NDOF)

    chkmass = 0.0d0
    do i = 1, N
      ii = NDOF*(i-1) + 1
      chkmass = chkmass + fstrEIG%mass(ii)
    end do
    call hecmw_allreduce_R1(hecMESH, chkmass, hecmw_sum)
    fstrEIG%totalmass = chkmass

    if(myrank == 0)then
      write(IMSG,"(a,1pe12.5)")"** Total mass: ", chkmass
    endif
  end subroutine setMASS

end module m_fstr_EIG_setMASS

