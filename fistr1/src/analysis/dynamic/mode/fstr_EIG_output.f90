!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_fstr_EIG_output

contains

  subroutine fstr_eigen_output(hecMESH, hecMAT, fstrEIG)
    use m_fstr
    use m_fstr_EIG_lanczos_util
    use hecmw_solver_las
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG

    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: i, j, k, in, jn, nget
    real(kind=kreal)   :: chk, gm
    real(kind=kreal), allocatable :: s(:), t(:), u(:), r(:)
    real(kind=kreal), pointer     :: mass(:), eigval(:), eigvec(:,:)

    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    nget   =  fstrEIG%nget
    mass   => fstrEIG%mass
    eigval => fstrEIG%eigval
    eigvec => fstrEIG%eigvec

    allocate(fstrEIG%effmass(NDOF*nget))
    allocate(fstrEIG%partfactor(NDOF*nget))
    allocate(s(NPNDOF))
    allocate(t(NPNDOF))
    allocate(u(NPNDOF))
    allocate(r(NDOF))
    fstrEIG%effmass    = 0.0d0
    fstrEIG%partfactor = 0.0d0

    do i = 1, nget
      r  = 0.0d0
      gm = 0.0d0
      do j = 1, N
        do k = 1, NDOF
          in = NDOF*(j-1) + k
          r(k) = r(k) + mass(in)*eigvec(in,i)
          gm   = gm   + mass(in)*eigvec(in,i)*eigvec(in,i)
        enddo
      enddo

      call hecmw_allreduce_R(hecMESH, r, NDOF, hecmw_sum)
      call hecmw_allreduce_R1(hecMESH, gm, hecmw_sum)

      gm = 1.0d0/gm
      do j = 1, NDOF
        in = NDOF*(i-1) + j
        fstrEIG%partfactor(in) = gm*r(j)
        fstrEIG%effmass(in)    = gm*r(j)*r(j)
      enddo
    enddo

    call EGLIST(hecMESH, hecMAT, fstrEIG)

    if(myrank == 0)then
      write(IMSG,*) ''
      write(IMSG,*) '*----------------------------------------------*'
      write(IMSG,*) '*--E I G E N V A L U E  C O N V E R G E N C E--*'
      write(IMSG,*) 'Absolute residual =  |(||Kx - lambda*Mx||)|'
      write(IMSG,*) '   Iter.#   Eigenvalue    Abs. Residual  '
      write(ILOG,*) '   Iter.#   Eigenvalue    Abs. Residual  '
      write(IMSG,*) '   *-----*  *---------*  *--------------*'
    endif

    do j = 1, fstrEIG%nget
      do i = 1, NNDOF
        u(i) = eigvec(i,j)
      enddo
      call hecmw_matvec(hecMESH, hecMAT, u, t)

      s = 0.0d0
      do i = 1, NNDOF
        s(i) = mass(i)*eigvec(i,j)
      enddo

      chk = 0.0d0
      do i = 1,NNDOF
        chk = chk + (t(i) - (eigval(j)-fstrEIG%sigma)*s(i))**2
      enddo
      call hecmw_allreduce_R1(hecMESH, chk, hecmw_sum)
      chk = dsqrt(chk)

      if(myrank == 0)then
        write(IMSG,'(2x,i5,2x,1p5e15.6)') j, eigval(j), chk
        write(ILOG,'(i5,1p5e12.4)') j, eigval(j), chk
      endif
    enddo

    if(myrank == 0)then
      write(IMSG,*)'*    ---END Eigenvalue listing---     *'
    endif

    deallocate(s)
    deallocate(t)
    deallocate(u)
    deallocate(r)
  end subroutine fstr_eigen_output

  subroutine fstr_eigen_make_result(hecMESH, hecMAT, fstrEIG, fstrRESULT)
    use m_fstr
    use m_hecmw2fstr_mesh_conv
    use hecmw_util
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(fstr_eigen)          :: fstrEIG
    type(hecmwST_result_data) :: fstrRESULT

    integer(kind=kint) :: i, istep, nget, NP, NDOF, NPNDOF, totalmpc, MPC_METHOD
    real(kind=kreal)   :: t1
    real(kind=kreal), pointer :: eigvec(:,:)
    real(kind=kreal), allocatable :: X(:), egval(:)
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN)    :: comment
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    nget   = fstrEIG%nget
    NDOF   = hecMAT%NDOF
    NPNDOF = hecMAT%NP*hecMAT%NDOF
    !totalmpc = hecMESH%mpc%n_mpc
    !call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    eigvec => fstrEIG%eigvec

    allocate(X(NPNDOF))
    X = 0.0d0
    allocate(egval(1))

    do istep = 1, nget
      egval(1) = fstrEIG%eigval(istep)
      do i=1,NPNDOF
        X(i) = eigvec(i,istep)
      enddo

      !if (totalmpc > 0) then
      !  MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
      !  if (MPC_METHOD < 1 .or. 3 < MPC_METHOD) MPC_METHOD = 3
      !  if (MPC_METHOD == 3) then  ! elimination
      !    call hecmw_tback_x_33(hecMESH, X, t1)
      !  else
      !    if (hecMESH%my_rank.eq.0) write(0,*) "### ERROR: MPC_METHOD must set to 3"
      !    stop
      !  endif
      !endif

      call hecmw_update_R(hecMESH, X, hecMAT%NP, NDOF)

      if( IRESULT.eq.1 ) then
        header = "*fstrresult"
        comment = "eigen_result"
        call hecmw_result_init(hecMESH,istep,header,comment)
        label = "EIGENVALUE"
        call hecmw_result_add(HECMW_RESULT_DTYPE_GLOBAL,1,label,egval)
        label = "DISPLACEMENT"
        call hecmw_result_add(HECMW_RESULT_DTYPE_NODE,NDOF,label,X)
        nameID = "fstrRES"
        call hecmw_result_write_by_name(nameID)
        call hecmw_result_finalize
      endif

      if( IVISUAL.eq.1 ) then
        call hecmw_nullify_result_data(fstrRESULT)
        fstrRESULT%ng_component = 1
        fstrRESULT%nn_component = 1
        fstrRESULT%ne_component = 0
        allocate(fstrRESULT%ng_dof(1))
        allocate(fstrRESULT%global_label(1))
        allocate(fstrRESULT%global_val_item(1))
        fstrRESULT%ng_dof(1) = 1
        fstrRESULT%global_label(1) = 'EIGENVALUE'
        fstrRESULT%global_val_item(1) = egval(1)
        allocate(fstrRESULT%nn_dof(1))
        allocate(fstrRESULT%node_label(1))
        allocate(fstrRESULT%node_val_item(NDOF*hecMAT%NP))
        fstrRESULT%nn_dof(1) = NDOF
        fstrRESULT%node_label(1) = 'DISPLACEMENT'
        fstrRESULT%node_val_item = X
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize( hecMESH, fstrRESULT, istep )
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif
    enddo

    deallocate(X)

  end subroutine fstr_eigen_make_result

  !> Output eigenvalues and vectors
  subroutine EGLIST(hecMESH, hecMAT, fstrEIG)
    use m_fstr
    use hecmw_util
    use m_fstr_EIG_lanczos_util
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG

    integer(kind=kint) :: NDOF
    integer(kind=kint) :: i, j, in, iter ,nget
    real(kind=kreal)   :: pi, angle, freq, pf(3), em(3)
    real(kind=kreal), pointer :: eigval(:)

    NDOF = hecMAT%NDOF
    nget = fstrEIG%nget
    iter = fstrEIG%iter
    PI   = 4.0d0 * datan(1.0d0)
    eigval => fstrEIG%eigval

    if(myrank == 0)then
      write(ILOG,*)""
      write(ILOG,"(a)")"********************************"
      write(ILOG,"(a)")"*RESULT OF EIGEN VALUE ANALYSIS*"
      write(ILOG,"(a)")"********************************"
      write(ILOG,"(a)")""
      write(ILOG,"(a,i8)")"NUMBER OF ITERATIONS = ",iter
      write(ILOG,"(a,1pe12.4)")"TOTAL MASS = ",fstrEIG%totalmass
      write(ILOG,"(a)")""
      write(ILOG,"(3a)")"                   ANGLE       FREQUENCY   ",&
        "PARTICIPATION FACTOR                EFFECTIVE MASS"
      write(ILOG,"(3a)")"  NO.  EIGENVALUE  FREQUENCY   (HZ)        ",&
        "X           Y           Z           X           Y           Z"
      write(ILOG,"(3a)")"  ---  ----------  ----------  ----------  ",&
        "----------  ----------  ----------  ----------  ----------  ----------"
      write(*,*)""
      write(*,"(a)")"#----------------------------------#"
      write(*,"(a)")"#  RESULT OF EIGEN VALUE ANALYSIS  #"
      write(*,"(a)")"#----------------------------------#"
      write(*,"(a)")""
      write(*,"(a,i8)")"### NUMBER OF ITERATIONS = ",iter
      write(*,"(a,1pe12.4)")"### TOTAL MASS = ",fstrEIG%totalmass
      write(*,"(a)")""
      write(*,"(3a)")"       PERIOD     FREQUENCY  ",&
        "PARTICIPATION FACTOR             EFFECTIVE MASS"
      write(*,"(3a)")"  NO.  [Sec]      [HZ]       ",&
        "X          Y          Z          X          Y          Z"
      write(*,"(3a)")"  ---  ---------  ---------  ",&
        "---------  ---------  ---------  ---------  ---------  ---------"

      in = 0
      do i = 1, nget
        in = in + 1
        if(eigval(i) < 0.0d0) eigval(i) = 0.0d0
        angle = dsqrt(eigval(i))
        freq = angle*0.5d0/PI

        pf = 0.0d0
        em = 0.0d0
        do j = 1, min(NDOF, 3)
          pf(j) = fstrEIG%partfactor(NDOF*(i-1) + j)
          em(j) = fstrEIG%effmass(NDOF*(i-1) + j)
        enddo

        write(ILOG,'(I5,1P9E12.4)') in, eigval(i), angle, freq, pf(1), pf(2), pf(3), em(1), em(2), em(3)
        write(*   ,'(I5,1P8E11.3)') in, 1.0d0/freq, freq, pf(1), pf(2), pf(3), em(1), em(2), em(3)
      enddo
      write(ILOG,*)""
      write(*,*)""
    endif

  end subroutine EGLIST
end module m_fstr_EIG_output


