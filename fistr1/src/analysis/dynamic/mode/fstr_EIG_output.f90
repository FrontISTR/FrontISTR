!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
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
    integer(kind=kint) :: i, j, jn, nget
    integer(kind=kint) :: in1, in2, in3
    real(kind=kreal)   :: chk, gm, r1, r2, r3
    real(kind=kreal), allocatable :: s(:), t(:), u(:)
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

    allocate(fstrEIG%effmass(3*NGET))
    allocate(fstrEIG%partfactor(3*NGET))
    allocate(s(NPNDOF))
    allocate(t(NPNDOF))
    allocate(u(NPNDOF))
    fstrEIG%effmass    = 0.0d0
    fstrEIG%partfactor = 0.0d0

    do i=1, NGET
      r1 = 0.0d0
      r2 = 0.0d0
      r3 = 0.0d0
      gm = 0.0d0
      if(NDOF == 3)then
        do j=1, N
          in1 = 3*j-2
          in2 = 3*j-1
          in3 = 3*j
          r1 = r1 + mass(in1)*eigvec(in1,i)
          r2 = r2 + mass(in2)*eigvec(in2,i)
          r3 = r3 + mass(in3)*eigvec(in3,i)
          gm = gm + mass(in1)*eigvec(in1,i)*eigvec(in1,i) &
            & + mass(in2)*eigvec(in2,i)*eigvec(in2,i) &
            & + mass(in3)*eigvec(in3,i)*eigvec(in3,i)
        enddo
      elseif(NDOF == 2)then
        do j=1, N
          in1 = 2*j-1
          in2 = 2*j
          r1 = r1 + mass(in1)*eigvec(in1,i)
          r2 = r2 + mass(in2)*eigvec(in2,i)
          gm = gm + mass(in1)*eigvec(in1,i)*eigvec(in1,i) &
            & + mass(in2)*eigvec(in2,i)*eigvec(in2,i)
        enddo
      endif
      call hecmw_allreduce_R1(hecMESH,r1,hecmw_sum)
      call hecmw_allreduce_R1(hecMESH,r2,hecmw_sum)
      call hecmw_allreduce_R1(hecMESH,r3,hecmw_sum)
      call hecmw_allreduce_R1(hecMESH,gm,hecmw_sum)

      fstrEIG%partfactor(3*i-2) = r1/gm
      fstrEIG%partfactor(3*i-1) = r2/gm
      fstrEIG%partfactor(3*i  ) = r3/gm
      fstrEIG%effmass(3*i-2) = r1*r1/gm
      fstrEIG%effmass(3*i-1) = r2*r2/gm
      fstrEIG%effmass(3*i  ) = r3*r3/gm
    enddo

    call EGLIST(hecMESH, hecMAT, fstrEIG)

    if(myrank == 0)then
      write(IMSG,*) ''
      write(IMSG,*) '*----------------------------------------------*'
      write(IMSG,*) '*--E I G E N V A L U E  C O N V E R G E N C E--*'
      write(IMSG,*) 'Absolute residual =  |(||Kx - lambda*Mx||)|'
      write(IMSG,*) '   Iter.#   Eigenvalue    Abs. Residual  '
      write(IMSG,*) '   *-----*  *---------*  *--------------*'
    endif

    do j = 1,fstrEIG%nget
      jn = j

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
        chk = chk + (t(i) - (eigval(jn)-fstrEIG%sigma)*s(i))**2
      enddo
      call hecmw_allreduce_R1(hecMESH, chk, hecmw_sum)
      chk = dsqrt(chk)

      if(myrank == 0)then
        write(IMSG,'(2X,I5,2X,5E15.6)') jn, eigval(jn), chk
      endif
    enddo

    if(myrank == 0)then
      write(IMSG,*)'*    ---END Eigenvalue listing---     *'
    endif

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
    real(kind=kreal), allocatable :: X(:)
    character(len=HECMW_HEADER_LEN) :: header
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

    do istep = 1, nget
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

      call hecmw_update_m_R(hecMESH, X, hecMAT%NP, NDOF)

      if( IRESULT.eq.1 ) then
        header = "*fstrresult"
        call hecmw_result_init(hecMESH,nget,istep,header)
        label = "DISPLACEMENT"
        call hecmw_result_add(1,NDOF,label,X)
        nameID = "fstrRES"
        call hecmw_result_write_by_name(nameID)
        call hecmw_result_finalize
      endif

      if( IVISUAL.eq.1 ) then
        call hecmw_nullify_result_data(fstrRESULT)
        fstrRESULT%nn_component = 1
        fstrRESULT%ne_component = 0
        allocate(fstrRESULT%nn_dof(1))
        allocate(fstrRESULT%node_label(1))
        allocate(fstrRESULT%node_val_item(NDOF*hecMAT%NP))
        fstrRESULT%nn_dof(1) = NDOF
        fstrRESULT%node_label(1) = 'DISPLACEMENT'
        fstrRESULT%node_val_item = X
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize(hecMESH,fstrRESULT,istep,nget,1)
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

    integer(kind=kint) :: i, j, in, iter ,nget
    real(kind=kreal)   :: pi, EEE, WWW, FFF, PFX, PFY, PFZ, EMX, EMY, EMZ
    real(kind=kreal), pointer :: eigval(:)

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
      do i=1, nget
        in = in + 1
        EEE = eigval(i)
        if(EEE < 0.0d0) EEE = 0.0d0
        WWW = dsqrt(EEE)
        FFF = WWW*0.5/PI
        PFX = fstrEIG%partfactor(3*i-2)
        PFY = fstrEIG%partfactor(3*i-1)
        PFZ = fstrEIG%partfactor(3*i  )
        EMX = fstrEIG%effmass(3*i-2)
        EMY = fstrEIG%effmass(3*i-1)
        EMZ = fstrEIG%effmass(3*i  )
        write(ILOG,'(I5,1P9E12.4)') in, EEE,  WWW, FFF, PFX, PFY, PFZ, EMX, EMY, EMZ
        write(*   ,'(I5,1P8E11.3)') in, 1.0d0/FFF, FFF, PFX, PFY, PFZ, EMX, EMY, EMZ
      enddo
      write(ILOG,*)""
      write(*,*)""
    endif

  end subroutine EGLIST
end module m_fstr_EIG_output


