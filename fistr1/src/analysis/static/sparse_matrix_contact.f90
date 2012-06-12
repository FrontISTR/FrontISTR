!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by K. Goto (VINAS)                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provides conversion routines between HEC and FISTR data
!! structures and DOF based sparse matrix data structures (CSR/COO).
module m_sparse_matrix_contact
  use m_fstr
  use fstr_matrix_con_contact
  use m_sparse_matrix
  use m_sparse_matrix_hec
  implicit none

  private
  public :: sparse_matrix_contact_init_prof
  public :: sparse_matrix_contact_set_vals
  public :: sparse_matrix_contact_set_rhs
  public :: sparse_matrix_contact_get_rhs

contains

!!$#define NEED_DIAG

  subroutine sparse_matrix_contact_init_prof(spMAT, hecMAT, fstrMAT, hecMESH)
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint) :: ndof, ndof2, N_loc, NL, NU, NZ
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    N_loc=hecMAT%N*ndof+fstrMAT%num_lagrange
    if (sparse_matrix_is_sym(spMAT)) then
       NU=hecMAT%indexU(hecMAT%N)
       NZ=hecMAT%N*(ndof2+ndof)/2+NU*ndof2 &
            +fstrMAT%numU_lagrange*ndof
!!$#ifdef NEED_DIAG
!!$       NZ=NZ+fstrMAT%num_lagrange
!!$#endif

    else
       NL=hecMAT%indexL(hecMAT%N)
       NU=hecMAT%indexU(hecMAT%N)
       NZ=(hecMAT%N+NU+NL)*ndof2 &
            +(fstrMAT%numL_lagrange+fstrMAT%numU_lagrange)*ndof
    endif
    call sparse_matrix_init(spMAT, N_loc, NZ)
    call sparse_matrix_hec_set_conv_ext(spMAT, hecMESH, ndof)
    spMAT%timelog = hecMAT%Iarray(22)
    call sparse_matrix_contact_set_prof(spMAT, hecMAT, fstrMAT)
  end subroutine sparse_matrix_contact_init_prof

  subroutine sparse_matrix_contact_set_prof(spMAT, hecMAT, fstrMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    type(hecmwST_matrix), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange)
    integer(kind=kint) :: ndof, ndof2
    integer(kind=kint) :: m, i, idof, i0, ii, ls, le, l, j, j0, jdof, jdofs
    !integer(kind=kint) :: offset_l, offset_d, offset_u
    ! CONVERT TO CSR or COO STYLE
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    m=1
    do i=1,hecMAT%N
       do idof=1,ndof
          i0=spMAT%OFFSET+ndof*(i-1)
          ii=i0+idof
          if (spMAT%type==SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii-spMAT%OFFSET)=m
          ! Lower
          if (.not. sparse_matrix_is_sym(spMAT)) then
             ls=hecMAT%indexL(i-1)+1
             le=hecMAT%indexL(i)
             do l=ls,le
                j=hecMAT%itemL(l)
                !if (j <= hecMAT%N) then
                j0=spMAT%OFFSET+ndof*(j-1)
                !else
                !   j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
                !endif
                !offset_l=ndof2*(l-1)+ndof*(idof-1)
                do jdof=1,ndof
                   if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
                   spMAT%JCN(m)=j0+jdof
                   !spMAT%A(m)=hecMAT%AL(offset_l+jdof)
                   m=m+1
                enddo
             enddo
          endif
          ! Diag
          !offset_d=ndof2*(i-1)+ndof*(idof-1)
          if (sparse_matrix_is_sym(spMAT)) then; jdofs=idof; else; jdofs=1; endif
          do jdof=jdofs,ndof
             if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
             spMAT%JCN(m)=i0+jdof
             !spMAT%A(m)=hecMAT%D(offset_d+jdof)
             m=m+1
          enddo
          ! Upper
          ls=hecMAT%indexU(i-1)+1
          le=hecMAT%indexU(i)
          do l=ls,le
             j=hecMAT%itemU(l)
             if (j <= hecMAT%N) then
                j0=spMAT%OFFSET+ndof*(j-1)
             else
                j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
                if (sparse_matrix_is_sym(spMAT) .and. j0 < i0) cycle
             endif
             !offset_u=ndof2*(l-1)+ndof*(idof-1)
             do jdof=1,ndof
                if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
                spMAT%JCN(m)=j0+jdof
                !spMAT%A(m)=hecMAT%AU(offset_u+jdof)
                m=m+1
             enddo
          enddo
          ! Upper Lagrange
          if (fstrMAT%num_lagrange > 0) then
            j0=spMAT%OFFSET+ndof*hecMAT%N
            ls=fstrMAT%indexU_lagrange(i-1)+1
            le=fstrMAT%indexU_lagrange(i)
            do l=ls,le
              j=fstrMAT%itemU_lagrange(l)
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
              spMAT%JCN(m)=j0+j
              m=m+1
            enddo
          endif
       enddo
    enddo
    ! Lower Lagrange
    if (fstrMAT%num_lagrange > 0) then
      i0=spMAT%OFFSET+ndof*hecMAT%N
      if (.not. sparse_matrix_is_sym(spMAT)) then
        do i=1,fstrMAT%num_lagrange
          ii=i0+i
          if (spMAT%type==SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii-spMAT%OFFSET)=m
          ls=fstrMAT%indexL_lagrange(i-1)+1
          le=fstrMAT%indexL_lagrange(i)
          do l=ls,le
            j=fstrMAT%itemL_lagrange(l)
            j0=spMAT%OFFSET+ndof*(j-1)
            do jdof=1,ndof
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
              spMAT%JCN(m)=j0+jdof
              m=m+1
            enddo
          enddo
        enddo
!!$#ifdef NEED_DIAG
!!$      else
!!$        do i=1,fstrMAT%num_lagrange
!!$          ii=i0+i
!!$          if (spMAT%type==SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii-spMAT%OFFSET)=m
!!$          if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
!!$          spMAT%JCN(m)=ii
!!$          m=m+1
!!$        enddo
!!$#endif
      endif
    endif
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii+1-spMAT%OFFSET)=m
    if (sparse_matrix_is_sym(spMAT) .and. m-1 < spMAT%NZ) spMAT%NZ=m-1
    if (m-1 /= spMAT%NZ) then
       write(*,*) 'ERROR: sparse_matrix_contact_set_prof on rank ',myrank
       write(*,*) 'm-1 = ',m-1,', NZ=',spMAT%NZ,',num_lagrange=',fstrMAT%num_lagrange
       call hecmw_abort(hecmw_comm_get_comm())
    endif
  end subroutine sparse_matrix_contact_set_prof

  subroutine sparse_matrix_contact_set_vals(spMAT, hecMAT, fstrMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    type(hecmwST_matrix), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange)
    integer(kind=kint) :: ndof, ndof2
    integer(kind=kint) :: m, i, idof, i0, ii, ls, le, l, j, j0, jdof, jdofs
    integer(kind=kint) :: offset_l, offset_d, offset_u
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    m=1
    do i=1,hecMAT%N
       do idof=1,ndof
          i0=spMAT%OFFSET+ndof*(i-1)
          ii=i0+idof
          if (spMAT%type==SPARSE_MATRIX_TYPE_CSR .and. spMAT%IRN(ii-spMAT%OFFSET)/=m) &
               stop "ERROR: sparse_matrix_contact_set_a"
          ! Lower
          if (.not. sparse_matrix_is_sym(spMAT)) then
             ls=hecMAT%indexL(i-1)+1
             le=hecMAT%indexL(i)
             do l=ls,le
                j=hecMAT%itemL(l)
                !if (j <= hecMAT%N) then
                j0=spMAT%OFFSET+ndof*(j-1)
                !else
                !   j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
                !endif
                offset_l=ndof2*(l-1)+ndof*(idof-1)
                do jdof=1,ndof
                   if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
                        stop "ERROR: sparse_matrix_contact_set_a"
                   if (spMAT%JCN(m)/=j0+jdof) stop "ERROR: sparse_matrix_contact_set_a"
                   spMAT%A(m)=hecMAT%AL(offset_l+jdof)
                   m=m+1
                enddo
             enddo
          endif
          ! Diag
          offset_d=ndof2*(i-1)+ndof*(idof-1)
          if (sparse_matrix_is_sym(spMAT)) then; jdofs=idof; else; jdofs=1; endif
          do jdof=jdofs,ndof
             if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
                  stop "ERROR: sparse_matrix_contact_set_a"
             if (spMAT%JCN(m)/=i0+jdof) stop "ERROR: sparse_matrix_contact_set_a"
             spMAT%A(m)=hecMAT%D(offset_d+jdof)
             m=m+1
          enddo
          ! Upper
          ls=hecMAT%indexU(i-1)+1
          le=hecMAT%indexU(i)
          do l=ls,le
             j=hecMAT%itemU(l)
             if (j <= hecMAT%N) then
                j0=spMAT%OFFSET+ndof*(j-1)
             else
                j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
                if (sparse_matrix_is_sym(spMAT) .and. j0 < i0) cycle
             endif
             offset_u=ndof2*(l-1)+ndof*(idof-1)
             do jdof=1,ndof
                if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
                     stop "ERROR: sparse_matrix_contact_set_a"
                if (spMAT%JCN(m)/=j0+jdof) stop "ERROR: sparse_matrix_contact_set_a"
                spMAT%A(m)=hecMAT%AU(offset_u+jdof)
                m=m+1
             enddo
          enddo
          ! Upper Lagrange
          if (fstrMAT%num_lagrange > 0) then
            j0=spMAT%OFFSET+ndof*hecMAT%N
            ls=fstrMAT%indexU_lagrange(i-1)+1
            le=fstrMAT%indexU_lagrange(i)
            do l=ls,le
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
                   stop "ERROR: sparse_matrix_contact_set_a"
              if (spMAT%JCN(m)/=j0+fstrMAT%itemU_lagrange(l)) &
                   stop "ERROR: sparse_matrix_contact_set_a"
              spMAT%A(m)=fstrMAT%AU_lagrange((l-1)*ndof+idof)
              m=m+1
            enddo
          endif
       enddo
    enddo
    ! Lower Lagrange
    if (fstrMAT%num_lagrange > 0) then
      i0=spMAT%OFFSET+ndof*hecMAT%N
      if (.not. sparse_matrix_is_sym(spMAT)) then
        do i=1,fstrMAT%num_lagrange
          ii=i0+i
          if (spMAT%type == SPARSE_MATRIX_TYPE_CSR .and. spMAT%IRN(ii-spMAT%OFFSET)/=m) &
               stop "ERROR: sparse_matrix_contact_set_a"
          ls=fstrMAT%indexL_lagrange(i-1)+1
          le=fstrMAT%indexL_lagrange(i)
          do l=ls,le
            j=fstrMAT%itemL_lagrange(l)
            j0=spMAT%OFFSET+ndof*(j-1)
            do jdof=1,ndof
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
                   stop "ERROR: sparse_matrix_contact_set_a"
              if (spMAT%JCN(m)/=j0+jdof) &
                   stop "ERROR: sparse_matrix_contact_set_a"
              spMAT%A(m)=fstrMAT%AL_lagrange((l-1)*ndof+jdof)
              m=m+1
            enddo
          enddo
        enddo
!!$#ifdef NEED_DIAG
!!$      else
!!$        do i=1,fstrMAT%num_lagrange
!!$          ii=i0+i
!!$          if (spMAT%type==SPARSE_MATRIX_TYPE_CSR .and. spMAT%IRN(ii-spMAT%OFFSET)/=m) &
!!$               stop "ERROR: sparse_matrix_contact_set_a"
!!$          if (spMAT%type==SPARSE_MATRIX_TYPE_COO .and. spMAT%IRN(m)/=ii) &
!!$               stop "ERROR: sparse_matrix_contact_set_a"
!!$          if (spMAT%JCN(m)/=ii) &
!!$               stop "ERROR: sparse_matrix_contact_set_a"
!!$          spMAT%A(m)=0.0d0
!!$          m=m+1
!!$        enddo
!!$#endif
      endif
    endif
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR .and. spMAT%IRN(ii+1-spMAT%OFFSET)/=m) &
         stop "ERROR: sparse_matrix_contact_set_a"
    if (m-1 /= spMAT%NZ) stop "ERROR: sparse_matrix_contact_set_a"
  end subroutine sparse_matrix_contact_set_vals

  subroutine sparse_matrix_contact_set_rhs(spMAT, hecMAT, fstrMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: nndof,npndof,ierr
    allocate(spMAT%rhs(spMAT%N_loc), stat=ierr)
    if (ierr /= 0) then
      write(*,*) " Allocation error, spMAT%rhs"
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    nndof=hecMAT%N*hecMAT%ndof
    npndof=hecMAT%NP*hecMAT%ndof
    spMAT%rhs(1:nndof) = hecMAT%b(1:nndof)
    ! skip external DOF in hecMAT%b
    if (fstrMAT%num_lagrange > 0) &
         spMAT%rhs(nndof+1:spMAT%N_loc) = hecMAT%b(npndof+1:npndof+fstrMAT%num_lagrange)
  end subroutine sparse_matrix_contact_set_rhs

  subroutine sparse_matrix_contact_get_rhs(spMAT, hecMAT, fstrMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: nndof,npndof
    nndof=hecMAT%N*hecMAT%ndof
    npndof=hecMAT%NP*hecMAT%ndof
    hecMAT%x(1:nndof) = spMAT%rhs(1:nndof)
    ! TODO: call update to get external DOF from neighboring domains???
    if (fstrMAT%num_lagrange > 0) &
         hecMAT%x(npndof+1:npndof+fstrMAT%num_lagrange) = spMAT%rhs(nndof+1:spMAT%N_loc)
    deallocate(spMAT%rhs)
  end subroutine sparse_matrix_contact_get_rhs

end module m_sparse_matrix_contact
