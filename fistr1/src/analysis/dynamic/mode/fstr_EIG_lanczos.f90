!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos
  contains

!> SOLVE EIGENVALUE PROBLEM
  subroutine fstr_solve_lanczos(hecMESH, hecMAT, fstrEIG, maxItr)
    use m_fstr
    use m_fstr_StiffMatrix
    use m_fstr_AddBC
    use fstr_matrix_con_contact
    use hecmw_util
    use lczeigen
    use m_eigen_lib
    use m_fstr_EIG_getamat
    use m_fstr_EIG_lanczos_util
    use m_fstr_EIG_matmult
    use m_fstr_EIG_mgs1
    use m_fstr_EIG_output
    use m_fstr_EIG_setMASS
    use m_fstr_EIG_tridiag
    use m_static_lib
    use m_static_make_result
    use m_hecmw2fstr_mesh_conv

    implicit none

    type (hecmwST_local_mesh ) :: hecMESH
    type (hecmwST_matrix     ) :: hecMAT
    type (hecmwST_result_data) :: fstrRESULT
    type (fstr_param         ) :: fstrPARAM
    type (fstrST_matrix_contact_lagrange)  :: fstrMAT
    type (fstr_eigen)          :: fstrEIG

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    type lczvec
      real(kind=kreal), pointer :: q(:) => null()
    end type lczvec
    TYPE(lczvec), pointer :: lvecq(:)       !< Array of Q vectors

    integer(kind=kint) :: i, j, k , ii, iii, ik, in, in1, in2, in3, nstep, istep, maxItr
    integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0, jiter, iiter, kiter
    integer(kind=kint) :: kk, jjiter, ppc
    integer(kind=kint) IOUT,IREOR,eITMAX,itype,iS,iE,ic_type,icel,jS,nn
    real(kind=kreal)   :: t1, t2, aalf, tmp, tmp2, gm, gm2, r1, r2, r3, r4, r5, r6

    integer(kind=kint), allocatable :: isnode33(:)
    real(kind=kreal), allocatable   :: gmass(:)

    numnp  = hecMAT%NP
    numn   = hecMAT%N
    NDOF   = hecMESH%n_dof
    ntotal = numnp*NDOF

    allocate(lvecq(0:lvecq_size))
    allocate( mass( ntotal )            )
    allocate( EM( ntotal )              )
    allocate( ewk( ntotal, neig)        )
    allocate( eval( neig )              )
    allocate( modal( neig )             )
    allocate( new( neig )               )
    allocate( work( neig*(3*neig + 5) ) )
    allocate( LVECP(NTOTAL)             )
    allocate( LVECPP(NTOTAL)            )
    allocate( lvecq(0)%q(ntotal)        )
    allocate( lvecq(1)%q(ntotal)        )
    allocate( LWRK(NTOTAL)              )
    allocate( LLWRK(NTOTAL)             )
    allocate( LLLWRK(NTOTAL)            )
    allocate( ALF(NEIG+2)               )
    allocate( BTA(NEIG+2)               )
    ALLOCATE(my_ntotal(0:nprocs), STAT=ierror)

    mass       = 0.0
    ewk        = 0.0
    eval       = 0.0
    new        = 0.0
    work       = 0.0
    lwrk       = 0.0
    LVECP      = 0.0
    LVECPP     = 0.0
    lvecq(0)%q = 0.0
    lvecq(1)%q = 0.0
    LWRK       = 0.0
    LLWRK      = 0.0
    ALF        = 0.0
    BTA        = 0.0

    do i=1,ntotal
      mass(i) = fstrEIG%mass(i)
    enddo

    my_ntotal(myrank) = Gntotal
      DO I = 0,nprocs-1
        call hecmw_bcast_I1(hecMESH,my_ntotal(I),I)
      enddo

    call SETIVL(mass,EM,EFILT,ewk,LVECP,lvecq(0)%q,lvecq(1)%q,BTA,ntotal,&
   &                         neig,my_ntotal,hecMESH,hecMAT,NDOF,Gtotal)

    CONV = .FALSE.

    do i=1,ntotal
      hecMAT%B(i) = EM(i)
    enddo
    hecMAT%X = 0.
    hecMAT%Iarray(98) = 1   !Assmebly complete
    hecMAT%Iarray(97) = 1   !Need numerical factorization
    hecMAT%Rarray(2)  = 1.0
    hecMAT%Rarray(1)  = fstrEIG%iluetol

    IF(myrank .EQ. 0) THEN
      WRITE(IMSG,*)
      WRITE(IMSG,*) ' *****   STAGE Begin Lanczos loop     **'
    ENDIF

    DO ITER=1,maxItr
      call DUPL(EM,LVECP,ntotal)
      call hecmw_mat_clear_b(hecMAT)

      do i=1,ntotal
        hecMAT%B(i) = EM(i)
      enddo
      hecMAT%X = 0.0d0

      call solve_LINEQ( hecMESH,hecMAT,IMSG )

      do i=1,ntotal
        EM(i) = hecMAT%X(i)
      enddo

      i = 0
      do ii = 1,ntotal
        IF(efilt(ii) .EQ. 0) i = i + 1
      enddo

      do ii = 1,ntotal
        EM(ii) = EM(ii)*EFILT(ii)
      enddo

      ALLOCATE( lvecq(iter+1)%q(ntotal), STAT=ierror )

      call SCSHFT(EM,lvecq(iter-1)%q,BTA(ITER),ntotal)

      call VECPRO1(AALF,LVECP,EM,Gntotal)
      ALF(ITER)=AALF

      call hecmw_allreduce_R1(hecMESH,ALF(ITER),hecmw_sum)

      call SCSHFT(EM,lvecq(iter)%q,ALF(ITER),ntotal)

      LWRK = 0.
      call MATPRO(LWRK,mass,EM,ntotal,1)
      call VECPRO1(prechk,LWRK,LWRK,Gntotal)
        call hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)

      prechk = sqrt(prechk)
      IF(prechk.NE.0.0D0) LWRK = LWRK/prechk
      IREOR = ITER*(1.0 - fstrEIG%lczrod)

      DO KK = IREOR,ITER
        prechk = 0.0D0
        call VECPRO1(prechk,lvecq(kk)%q,LWRK,Gntotal)
          call hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)

        prechk1 = 0.0D0
        call VECPRO1(prechk1,lvecq(kk)%q,lvecq(kk)%q,Gntotal)
          call hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)

        prechk1 = sqrt(prechk1)
        if(prechk1.ne.0.0D0) prechk = prechk/prechk1
        IF(abs(prechk).GT.fstrEIG%lczrot) THEN
          call MGS1(lvecq(kk)%q,EM,mass,ntotal,my_ntotal,myrank,&
   &                hecMESH,Gntotal)
        ENDIF
      enddo

      call MATPRO(LVECPP,mass,EM,ntotal,1)

      call VECPRO1(AALF,LVECPP,EM,Gntotal)
      BTA(ITER+1) = AALF
        call hecmw_allreduce_R1(hecMESH,BTA(ITER+1),hecmw_sum)

      BTA(ITER+1) = SQRT(BTA(ITER+1))

      call hecmw_barrier(hecMESH)
      call UPLCZ(LVECP,lvecq(iter+1)%q,LVECPP,EM,BTA(ITER+1),ntotal)

      LTRIAL = ITER
      ALLOCATE( LLDIAG(LTRIAL)        )
      ALLOCATE( LNDIAG(LTRIAL)        ) !Unordered
      ALLOCATE( LSUB(LTRIAL)          )
      ALLOCATE( LZMAT(LTRIAL,LTRIAL)  )
      ALLOCATE( LNZMAT(LTRIAL,LTRIAL) ) !Unordered

      DO JITER = 1,LTRIAL
        DO IITER = 1,LTRIAL
          LZMAT(JITER,IITER)  = 0.0D0
          LNZMAT(JITER,IITER) = 0.0D0 !Unordered
        enddo
      enddo

      DO IITER = 1,LTRIAL
        LLDIAG(IITER) = ALF(IITER)
        LZMAT(IITER,IITER) = 1.0D0
      enddo

      LSUB(1) = 0.0
      IITER   = 0
      DO IITER = 2,LTRIAL
        LSUB(IITER)  = BTA(IITER)
      enddo

      call TRIDIAG(LTRIAL,LTRIAL,LLDIAG,LNDIAG,&
   &                      LSUB,LZMAT,LNZMAT,IERROR) !Unordered

      DO IITER = 1, LTRIAL
        IF( LLDIAG(IITER).NE.0.0D0 ) THEN
          EVAL(IITER) = 1.0D0/LLDIAG(IITER) + fstrEIG%lczsgm
        ENDIF
      enddo
      call EVSORT(EVAL,NEW,LTRIAL)

      CCHK  = 0.0
      KITER = NGET+2            !Extra 2 values as a safety feature
      IF(LTRIAL .LT. KITER) KITER = LTRIAL

      DO JJITER = 1,KITER
        JITER = NEW(JJITER)
        call VECPRO1(prechk1,LZMAT(1,JITER),LZMAT(1,JITER),LTRIAL)

        IF(prechk1 .GT. 0.) THEN
          prechk1 = SQRT(prechk1)
        ELSE
          prechk1 = 1.
        ENDIF

        CCHK1 = (BTA(JITER))*(LZMAT(LTRIAL,JITER)/prechk)
        IF(cchk .LT. ABS(cchk1)) THEN
          cchk = ABS(cchk1)
          iiter = jiter
          ppc = prechk1
        ENDIF
      enddo

      IF(CCHK.LE.CTOL.AND.ITER.GE.NGET) THEN
        WRITE(IDBG,*) '*=====Desired convergence was obtained =====*'
        CONV = .TRUE.
        GO TO 25
      ENDIF

      IF(ITER.LT.maxItr) THEN
        if( allocated(LLDIAG) ) DEALLOCATE(LLDIAG)
        if( allocated(LNDIAG) ) DEALLOCATE(LNDIAG)
        if( allocated(LSUB) )   DEALLOCATE(LSUB)
        if( allocated(LZMAT) )  DEALLOCATE(LZMAT)
        if( allocated(LNZMAT) ) DEALLOCATE(LNZMAT)
      ENDIF

    enddo

25  CONTINUE
    CERR = CCHK

    DO IITER = 1, LTRIAL
      IF(LLDIAG(IITER).NE.0.0D0) THEN
        EVAL(IITER) = 1.0D0/LLDIAG(IITER) + fstrEIG%lczsgm
      ENDIF
    enddo
    call evsort(eval,new,ltrial)

    ewk = 0.0
    k = nget
    IF(k .GT. ltrial) k = ltrial
    DO kk = 1,k
      kiter = NEW(kk)
      DO jiter = 1,ltrial
        DO iiter =1,ntotal
          ewk(iiter,kk) = ewk(iiter,kk) &
   &                    + lvecq(jiter)%q(iiter)*LZMAT(jiter,kiter)
        enddo
      enddo
    enddo

    DO iiter=0,ltrial
      if( associated(lvecq(iiter)%q) ) DEALLOCATE(lvecq(iiter)%q)
    enddo

    modal = 0
    do i = 1,LTRIAL
      if(EVAL(i).NE.0) modal(i) = 1
    enddo

    lczmult = .FALSE.

    t2 = hecmw_Wtime() !DEBUG elap
    write(idbg,'(a,f10.2)') 'Lanczos loop (sec) :', T2 - T1 ! elap

    IF(myrank .EQ. 0) THEN
      WRITE(IMSG,*)
      WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
    ENDIF

    DO JITER=1,NGET
      !prechk1 = 0.0
!C        call VECPRO1(prechk1,ewk(1:,JITER:),ewk(1:,JITER:),Gntotal)
      !do iii = 1, Gntotal
      !    prechk1 = prechk1 + mass(iii)*ewk(iii,JITER)*ewk(iii,JITER)
      !enddo
      prechk1 = maxval(ewk(:,JITER))
        call hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)

      !prechk1 = sqrt(prechk1)
      if(prechk1.NE.0.0D0)then
        do i = 1, Gntotal
          ewk(i,JITER) = ewk(i,JITER)/prechk1
        enddo
      endif
    enddo

!C***** compute effective mass and participation factor
    allocate(fstrEIG%effmass(3*NGET))
    allocate(fstrEIG%partfactor(3*NGET))
    fstrEIG%effmass    = 0.0d0
    fstrEIG%partfactor = 0.0d0

    if(NDOF == 3)then
      DO i=1,NGET
        r1 = 0.0d0
        r2 = 0.0d0
        r3 = 0.0d0
        gm = 0.0d0
        do j = 1, numn
          in1 = 3*j-2
          in2 = 3*j-1
          in3 = 3*j
          r1 = r1 + mass(in1)*ewk(in1,i)
          r2 = r2 + mass(in2)*ewk(in2,i)
          r3 = r3 + mass(in3)*ewk(in3,i)
          gm = gm + mass(in1)*ewk(in1,i)*ewk(in1,i) &
          & + mass(in2)*ewk(in2,i)*ewk(in2,i) &
          & + mass(in3)*ewk(in3,i)*ewk(in3,i)
        enddo
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

    elseif(NDOF == 2)then
      DO i=1,NGET
        r1 = 0.0d0
        r2 = 0.0d0
        gm = 0.0d0
        do j = 1, numn
          in1 = 2*j-1
          in2 = 2*j
          r1 = r1 + mass(in1)*ewk(in1,i)
          r2 = r2 + mass(in2)*ewk(in2,i)
          gm = gm + r1*ewk(in1,i) + r2*ewk(in2,i)
        enddo
        call hecmw_allreduce_R1(hecMESH,r1,hecmw_sum)
        call hecmw_allreduce_R1(hecMESH,r2,hecmw_sum)
        call hecmw_allreduce_R1(hecMESH,gm,hecmw_sum)
        fstrEIG%partfactor(3*i-2) = r1/gm
        fstrEIG%partfactor(3*i-1) = r2/gm
        fstrEIG%effmass(3*i-2) = r1*r1/gm
        fstrEIG%effmass(3*i-1) = r2*r2/gm
      enddo
    endif

  end subroutine fstr_solve_lanczos

end module m_fstr_EIG_lanczos
