!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_matrix_ass
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_matrix_contact
  implicit none

  private

  public :: hecmw_mat_ass_elem
  public :: hecmw_mat_add_node
  public :: hecmw_array_search_i
  public :: hecmw_mat_ass_equation
  public :: hecmw_mat_ass_equation_rhs
  public :: hecmw_mat_add_dof
  public :: hecmw_mat_ass_bc
  public :: hecmw_mat_ass_bc_contactlag
  public :: hecmw_mat_ass_contact
  public :: hecmw_mat_ass_contactlag
  public :: stf_get_block

contains

  !C
  !C***
  !C*** MAT_ASS_ELEM
  !C***
  !C
  subroutine hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiffness)
    type (hecmwST_matrix)     :: hecMAT
    integer(kind=kint) :: nn
    integer(kind=kint) :: nodLOCAL(:)
    real(kind=kreal) :: stiffness(:, :)
    !** Local variables
    integer(kind=kint) :: ndof, inod_e, jnod_e, inod, jnod
    real(kind=kreal) :: a(6,6)

    ndof = hecMAT%NDOF

    do inod_e = 1, nn
      inod = nodLOCAL(inod_e)
      do jnod_e = 1, nn
        jnod = nodLOCAL(jnod_e)
        !***** Add components
        call stf_get_block(stiffness, ndof, inod_e, jnod_e, a)
        call hecmw_mat_add_node(hecMAT, inod, jnod, a)
      enddo
    enddo

  end subroutine hecmw_mat_ass_elem

  subroutine stf_get_block(stiffness, ndof, inod, jnod, a)
    real(kind=kreal) :: stiffness(:, :), a(:, :)
    integer(kind=kint) :: ndof, inod, jnod
    !** Local variables
    integer(kind=kint) :: row_offset, col_offset, i, j

    row_offset = ndof*(inod-1)
    do i = 1, ndof

      col_offset = ndof*(jnod-1)
      do j = 1, ndof

        a(i, j) = stiffness(i + row_offset, j + col_offset)
      enddo
    enddo
  end subroutine stf_get_block


  subroutine hecmw_mat_add_node(hecMAT, inod, jnod, a)
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: inod, jnod
    real(kind=kreal) :: a(:, :)
    !** Local variables
    integer(kind=kint) :: NDOF, is, iE, k, idx_base, idx, idof, jdof

    NDOF = hecMAT%NDOF

    if (inod < jnod) then
      is = hecMAT%indexU(inod-1)+1
      iE = hecMAT%indexU(inod)
      k = hecmw_array_search_i(hecMAT%itemU, is, iE, jnod)

      if (k < is .or. iE < k) then
        write(*,*) '###ERROR### : cannot find connectivity (1)'
        write(*,*) ' myrank = ', hecmw_comm_get_rank(), ', inod = ', inod, ', jnod = ', jnod
        call hecmw_abort(hecmw_comm_get_comm())
      endif

      idx_base = NDOF**2 * (k-1)
      do idof = 1, NDOF
        do jdof = 1, NDOF
          idx = idx_base + jdof
          !$omp atomic
          hecMAT%AU(idx) = hecMAT%AU(idx) + a(idof, jdof)
        enddo
        idx_base = idx_base + NDOF
      enddo

    else if (inod > jnod) then
      is = hecMAT%indexL(inod-1)+1
      iE = hecMAT%indexL(inod)
      k = hecmw_array_search_i(hecMAT%itemL, is, iE, jnod)

      if (k < is .or. iE < k) then
        write(*,*) '###ERROR### : cannot find connectivity (2)'
        write(*,*) ' myrank = ', hecmw_comm_get_rank(), ', inod = ', inod, ', jnod = ', jnod
        call hecmw_abort(hecmw_comm_get_comm())
      endif

      idx_base = NDOF**2 * (k-1)
      do idof = 1, NDOF
        do jdof = 1, NDOF
          idx = idx_base + jdof
          !$omp atomic
          hecMAT%AL(idx) = hecMAT%AL(idx) + a(idof, jdof)
        enddo
        idx_base = idx_base + NDOF
      enddo

    else
      idx_base = NDOF**2 * (inod - 1)
      do idof = 1, NDOF
        do jdof = 1, NDOF
          idx = idx_base + jdof
          !$omp atomic
          hecMAT%D(idx) = hecMAT%D(idx) + a(idof, jdof)
        enddo
        idx_base = idx_base + NDOF
      enddo
    endif
  end subroutine hecmw_mat_add_node


  function hecmw_array_search_i(array, is, iE, ival)
    integer(kind=kint) :: hecmw_array_search_i
    integer(kind=kint) :: array(:)
    integer(kind=kint) :: is, iE, ival
    !** Local variables
    integer(kind=kint) :: left, right, center, cval

    left = is
    right = iE
    do
      if (left > right) then
        center = -1
        exit
      endif

      center = (left + right) / 2
      cval = array(center)

      if (ival < cval) then
        right = center - 1
      else if (cval < ival) then
        left = center + 1
      else
        exit
      endif
    enddo

    hecmw_array_search_i = center

  end function hecmw_array_search_i


  !C
  !C***
  !C*** MAT_ASS_EQUATION
  !C***
  !C
  subroutine hecmw_mat_ass_equation ( hecMESH, hecMAT )
    type (hecmwST_matrix), target :: hecMAT
    type (hecmwST_local_mesh)     :: hecMESH
    !** Local variables
    real(kind=kreal), pointer :: penalty
    real(kind=kreal) :: ALPHA, a1_2inv, ai, aj, factor
    integer(kind=kint) :: impc, is, iE, i, j, inod, idof, jnod, jdof
    logical :: is_internal_i

    if( hecmw_mat_get_penalized(hecMAT) == 1 ) return

    ! write(*,*) "INFO: imposing MPC by penalty"

    penalty => hecMAT%Rarray(11)

    if (penalty < 0.0) stop "ERROR: negative penalty"
    if (penalty < 1.0) write(*,*) "WARNING: penalty ", penalty, " smaller than 1"

    ALPHA= hecmw_mat_diag_max(hecMAT, hecMESH) * penalty
    call hecmw_mat_set_penalty_alpha(hecMAT, ALPHA)

    OUTER: do impc = 1, hecMESH%mpc%n_mpc
      is = hecMESH%mpc%mpc_index(impc-1) + 1
      iE = hecMESH%mpc%mpc_index(impc)

      do i = is, iE
        if (hecMESH%mpc%mpc_dof(i) > hecMAT%NDOF) cycle OUTER
      enddo

      a1_2inv = 1.0 / hecMESH%mpc%mpc_val(is)**2


      do i = is, iE
        inod = hecMESH%mpc%mpc_item(i)

        is_internal_i = (hecMESH%node_ID(2*inod) == hecmw_comm_get_rank())
        if (.not. is_internal_i) cycle

        idof = hecMESH%mpc%mpc_dof(i)
        ai = hecMESH%mpc%mpc_val(i)
        factor = ai * a1_2inv

        do j = is, iE
          jnod = hecMESH%mpc%mpc_item(j)
          jdof = hecMESH%mpc%mpc_dof(j)
          aj = hecMESH%mpc%mpc_val(j)

          call hecmw_mat_add_dof(hecMAT, inod, idof, jnod, jdof, aj*factor*ALPHA)
        enddo

      enddo
    enddo OUTER

    call hecmw_mat_set_penalized(hecMAT, 1)

  end subroutine hecmw_mat_ass_equation


  subroutine hecmw_mat_ass_equation_rhs ( hecMESH, hecMAT )
    type (hecmwST_matrix), target :: hecMAT
    type (hecmwST_local_mesh)     :: hecMESH
    !** Local variables
    real(kind=kreal) :: ALPHA, a1_2inv, ai, factor, ci
    integer(kind=kint) :: ndof, impc, iS, iE, i, inod, idof

    if( hecmw_mat_get_penalized_b(hecMAT) == 1) return

    ALPHA = hecmw_mat_get_penalty_alpha(hecMAT)
    if (ALPHA <= 0.0) stop "ERROR: penalty applied on vector before matrix"

    ndof = hecMAT%NDOF

    OUTER: do impc = 1, hecMESH%mpc%n_mpc
      iS = hecMESH%mpc%mpc_index(impc-1) + 1
      iE = hecMESH%mpc%mpc_index(impc)

      do i = is, iE
        if (hecMESH%mpc%mpc_dof(i) > ndof) cycle OUTER
      enddo

      a1_2inv = 1.0 / hecMESH%mpc%mpc_val(iS)**2

      do i = iS, iE
        inod = hecMESH%mpc%mpc_item(i)

        idof = hecMESH%mpc%mpc_dof(i)
        ai = hecMESH%mpc%mpc_val(i)
        factor = ai * a1_2inv

        ci = hecMESH%mpc%mpc_const(impc)
        !$omp atomic
        hecMAT%B(ndof*(inod-1)+idof) = hecMAT%B(ndof*(inod-1)+idof) + ci*factor*ALPHA
      enddo
    enddo OUTER

    call hecmw_mat_set_penalized_b(hecMAT, 1)

  end subroutine hecmw_mat_ass_equation_rhs


  subroutine hecmw_mat_add_dof(hecMAT, inod, idof, jnod, jdof, val)
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: inod, idof, jnod, jdof
    real(kind=kreal) :: val
    !** Local variables
    integer(kind=kint) :: NDOF, is, iE, k, idx

    NDOF = hecMAT%NDOF
    if (inod < jnod) then
      is = hecMAT%indexU(inod-1)+1
      iE = hecMAT%indexU(inod)
      k = hecmw_array_search_i(hecMAT%itemU, is, iE, jnod)

      if (k < is .or. iE < k) then
        write(*,*) '###ERROR### : cannot find connectivity (3)'
        write(*,*) '  myrank = ', hecmw_comm_get_rank(), ', inod = ', inod, ', jnod = ', jnod
        call hecmw_abort(hecmw_comm_get_comm())
        return
      endif

      idx = NDOF**2 * (k-1) + NDOF * (idof-1) + jdof
      !$omp atomic
      hecMAT%AU(idx) = hecMAT%AU(idx) + val

    else if (inod > jnod) then
      is = hecMAT%indexL(inod-1)+1
      iE = hecMAT%indexL(inod)
      k = hecmw_array_search_i(hecMAT%itemL, is, iE, jnod)

      if (k < is .or. iE < k) then
        write(*,*) '###ERROR### : cannot find connectivity (4)'
        write(*,*) ' myrank = ', hecmw_comm_get_rank(), ', inod = ', inod, ', jnod = ', jnod
        call hecmw_abort(hecmw_comm_get_comm())
        return
      endif

      idx = NDOF**2 * (k-1) + NDOF * (idof-1) + jdof
      !$omp atomic
      hecMAT%AL(idx) = hecMAT%AL(idx) + val

    else
      idx = NDOF**2 * (inod - 1) + NDOF * (idof - 1) + jdof
      !$omp atomic
      hecMAT%D(idx) = hecMAT%D(idx) + val
    endif

  end subroutine hecmw_mat_add_dof

  !C
  !C***
  !C*** MAT_ASS_BC
  !C***
  !C
  subroutine hecmw_mat_ass_bc(hecMAT, inode, idof, RHS, conMAT)
    type (hecmwST_matrix)     :: hecMAT
    integer(kind=kint) :: inode, idof
    real(kind=kreal) :: RHS, val
    type (hecmwST_matrix),optional     :: conMAT
    integer(kind=kint) :: NDOF, in, i, ii, iii, ndof2, k, is, iE, iiS, iiE, ik, idx

    NDOF = hecMAT%NDOF
    if( NDOF < idof ) return

    !C-- DIAGONAL block

    hecMAT%B(NDOF*inode-(NDOF-idof)) = RHS
    if(present(conMAT)) conMAT%B(NDOF*inode-(NDOF-idof)) = 0.0D0
    ndof2 = NDOF*NDOF
    ii  = ndof2 - idof

    do i = NDOF-1,0,-1
      if( i .NE. NDOF-idof ) then
        idx = NDOF*inode-i
        val = hecMAT%D(ndof2*inode-ii)*RHS
        !$omp atomic
        hecMAT%B(idx) = hecMAT%B(idx) - val
        if(present(conMAT)) then
          val = conMAT%D(ndof2*inode-ii)*RHS
          !$omp atomic
          conMAT%B(idx) = conMAT%B(idx) - val
        endif
      endif
      ii = ii - NDOF
    end do

    !*Set diagonal row to zero
    ii  = ndof2-1 - (idof-1)*NDOF

    do i = 0, NDOF - 1
      hecMAT%D(ndof2*inode-ii+i)=0.d0
      if(present(conMAT)) conMAT%D(ndof2*inode-ii+i)=0.d0
    end do

    !*Set diagonal column to zero
    ii = ndof2 - idof
    do i = 1, NDOF
      if( i.NE.idof ) then
        hecMAT%D(ndof2*inode-ii) = 0.d0
        if(present(conMAT)) conMAT%D(ndof2*inode-ii) = 0.d0
      else
        hecMAT%D(ndof2*inode-ii) = 1.d0
        if(present(conMAT)) conMAT%D(ndof2*inode-ii) = 0.d0
      endif
      ii = ii - NDOF
    end do

    !C-- OFF-DIAGONAL blocks

    ii  = ndof2-1 - (idof-1)*NDOF
    is = hecMAT%indexL(inode-1) + 1
    iE = hecMAT%indexL(inode  )

    do k= is, iE

      !*row (left)
      do i = 0, NDOF - 1
        hecMAT%AL(ndof2*k-ii+i) = 0.d0
        if(present(conMAT)) conMAT%AL(ndof2*k-ii+i) = 0.d0
      end do

      !*column (upper)
      in = hecMAT%itemL(k)
      iiS = hecMAT%indexU(in-1) + 1
      iiE = hecMAT%indexU(in  )
      do ik = iiS, iiE
        if (hecMAT%itemU(ik) .eq. inode) then
          iii = ndof2 - idof
          do i = NDOF-1,0,-1
            idx = NDOF*in-i
            val = hecMAT%AU(ndof2*ik-iii)*RHS
            !$omp atomic
            hecMAT%B(idx) = hecMAT%B(idx) - val
            hecMAT%AU(ndof2*ik-iii)= 0.d0
            if(present(conMAT)) then
              val = conMAT%AU(ndof2*ik-iii)*RHS
              !$omp atomic
              conMAT%B(idx) = conMAT%B(idx) - val
              conMAT%AU(ndof2*ik-iii)= 0.d0
            endif
            iii = iii - NDOF
          end do
          exit
        endif
      enddo

    enddo

    ii = ndof2-1 - (idof-1)*NDOF
    is = hecMAT%indexU(inode-1) + 1
    iE = hecMAT%indexU(inode  )

    do k= is, iE

      !*row (right)
      do i = 0,NDOF-1
        hecMAT%AU(ndof2*k-ii+i) = 0.d0
        if(present(conMAT)) conMAT%AU(ndof2*k-ii+i) = 0.d0
      end do

      !*column (lower)
      in = hecMAT%itemU(k)
      iiS = hecMAT%indexL(in-1) + 1
      iiE = hecMAT%indexL(in  )
      do ik= iiS, iiE
        if (hecMAT%itemL(ik) .eq. inode) then
          iii  = ndof2 - idof

          do i = NDOF-1, 0, -1
            idx = NDOF*in-i
            val = hecMAT%AL(ndof2*ik-iii)*RHS
            !$omp atomic
            hecMAT%B(idx) = hecMAT%B(idx) - val
            hecMAT%AL(ndof2*ik-iii) = 0.d0
            if(present(conMAT)) then
              val = conMAT%AL(ndof2*ik-iii)*RHS
              !$omp atomic
              conMAT%B(idx) = conMAT%B(idx) - val
              conMAT%AL(ndof2*ik-iii) = 0.d0
            endif
            iii = iii - NDOF
          end do
          exit
        endif
      enddo

    enddo
    !*End off - diagonal blocks

    call hecmw_cmat_ass_bc(hecMAT, inode, idof, RHS)

  end subroutine hecmw_mat_ass_bc


  !> Modify Lagrange multiplier-related part of stiffness matrix and right-hand side vector
  !> for dealing with prescribed displacement boundary condition
  subroutine hecmw_mat_ass_bc_contactlag(hecMAT,hecLagMAT,inode,idof,RHS)

    type(hecmwST_matrix)                 :: hecMAT !< hecmwST_matrix
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< hecmwST_matrix_lagrange
    integer(kind=kint) :: inode, idof !< number of node; degree of freedom
    integer(kind=kint) :: isU, ieU, isL, ieL, i, l, k
    real(kind=kreal)   :: RHS !< value of prescribed displacement

    isU = hecLagMAT%indexU_lagrange(inode-1)+1
    ieU = hecLagMAT%indexU_lagrange(inode)
    do i = isU, ieU
      hecLagMAT%AU_lagrange((i-1)*3+idof) = 0.0d0
      l = hecLagMAT%itemU_lagrange(i)
      isL = hecLagMAT%indexL_lagrange(l-1)+1
      ieL = hecLagMAT%indexL_lagrange(l)
      k = hecmw_array_search_i(hecLagMAT%itemL_lagrange,isL,ieL,inode)
      if(k < isL .or. k > ieL) cycle
      hecMAT%B(hecMAT%NP*hecMAT%NDOF+l) = hecMAT%B(hecMAT%NP*hecMAT%NDOF+l) - hecLagMAT%AL_lagrange((k-1)*3+idof)*RHS
      hecLagMAT%AL_lagrange((k-1)*3+idof) = 0.0d0
    enddo

  end subroutine hecmw_mat_ass_bc_contactlag

  !C
  !C***
  !C*** MAT_ASS_CONTACT
  !C***
  !C
  subroutine hecmw_mat_ass_contact(hecMAT, nn, nodLOCAL, stiffness)
    type (hecmwST_matrix)     :: hecMAT
    integer(kind=kint) :: nn
    integer(kind=kint) :: nodLOCAL(:)
    real(kind=kreal) :: stiffness(:, :)
    !** Local variables
    integer(kind=kint) :: ndof, inod_e, jnod_e, inod, jnod
    real(kind=kreal) :: a(3,3)

    ndof = hecMAT%NDOF
    if( ndof .ne. 3 ) then
      write(*,*) '###ERROR### : ndof=',ndof,'; contact matrix supports only ndof==3'
      call hecmw_abort(hecmw_comm_get_comm())
      return
    endif

    do inod_e = 1, nn
      inod = nodLOCAL(inod_e)
      do jnod_e = 1, nn
        jnod = nodLOCAL(jnod_e)
        !***** Add components
        call stf_get_block(stiffness, ndof, inod_e, jnod_e, a)
        call hecmw_cmat_add(hecMAT%cmat, inod, jnod, a)
      enddo
    enddo
    call hecmw_cmat_pack(hecMAT%cmat)

  end subroutine hecmw_mat_ass_contact

  !> \brief This subroutine assembles contact stiffness matrix of a contact pair into global stiffness matrix
  subroutine hecmw_mat_ass_contactlag(nnode,ndLocal,id_lagrange,fcoeff,stiffness,hecMAT,hecLagMAT)

    type(hecmwST_matrix)                 :: hecMAT !< type hecmwST_matrix
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint) :: nnode, ndLocal(nnode + 1), id_lagrange !< total number of nodes of master segment
!< global number of nodes of contact pair
!< number of Lagrange multiplier
    integer(kind=kint) :: i, j, inod, jnod, l
    integer(kind=kint) :: isL, ieL, idxL_base, kL, idxL, isU, ieU, idxU_base, kU, idxU
    real(kind=kreal)   :: fcoeff !< friction coefficient
    real(kind=kreal)   :: stiffness(:,:) !< contact stiffness matrix
    real(kind=kreal)   :: a(3, 3)

    i = nnode + 1 + 1
    inod = id_lagrange
    isL = hecLagMAT%indexL_lagrange(inod-1)+1
    ieL = hecLagMAT%indexL_lagrange(inod)

    do j = 1, nnode + 1
      jnod = ndLocal(j)
      isU = hecLagMAT%indexU_lagrange(jnod-1)+1
      ieU = hecLagMAT%indexU_lagrange(jnod)

      kL = hecmw_array_search_i(hecLagMAT%itemL_lagrange,isL,ieL,jnod)
      if( kL<isL .or. kL>ieL ) then
        write(*,*) '###ERROR### : cannot find connectivity (Lagrange1)'
        stop
      endif
      kU = hecmw_array_search_i(hecLagMAT%itemU_lagrange,isU,ieU,inod)
      if( kU<isU .or. kU>ieU ) then
        write(*,*) '###ERROR### : cannot find connectivity (Lagrange2)'
        stop
      endif

      idxL_base = (kL-1)*3
      idxU_base = (kU-1)*3

      do l = 1, 3
        idxL = idxL_base + l
        hecLagMAT%AL_lagrange(idxL) = hecLagMAT%AL_lagrange(idxL) + stiffness((i-1)*3+1,(j-1)*3+l)
        idxU = idxU_base + l
        hecLagMAT%AU_lagrange(idxU) = hecLagMAT%AU_lagrange(idxU) + stiffness((j-1)*3+l,(i-1)*3+1)
      enddo
    enddo


    if(fcoeff /= 0.0d0)then

      do i = 1, nnode + 1
        inod = ndLocal(i)
        do j = 1, nnode + 1
          jnod = ndLocal(j)
          call stf_get_block(stiffness(1:(nnode+1)*3,1:(nnode+1)*3), 3, i, j, a)
          call hecmw_mat_add_node(hecMAT, inod, jnod, a)
        enddo
      enddo

    endif

  end subroutine hecmw_mat_ass_contactlag

end module hecmw_matrix_ass
