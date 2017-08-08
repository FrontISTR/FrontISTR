!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos
  contains

!> SOLVE EIGENVALUE PROBLEM
  subroutine fstr_solve_lanczos(hecMESH, hecMAT, fstrSOLID, fstrEIG)
    use m_fstr
    use hecmw_util
    use lczeigen
    use m_eigen_lib
    use m_fstr_EIG_getamat
    use m_fstr_EIG_lanczos_util
    use m_fstr_EIG_matmult
    use m_fstr_EIG_mgs1
    use m_fstr_EIG_tridiag

    implicit none

    type (hecmwST_local_mesh ) :: hecMESH
    type (hecmwST_matrix     ) :: hecMAT
    type (fstr_solid         ) :: fstrSOLID
    type (fstr_eigen         ) :: fstrEIG

    type fstr_eigen_vec
      real(kind=kreal), pointer :: q(:) => null()
    end type fstr_eigen_vec

    TYPE(fstr_eigen_vec), pointer :: lvecq(:) !< Array of Q vectors

    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: iter, ierr

    integer(kind=kint) :: i, j, k , ii, iii, ik, in, in1, in2, in3, nstep, istep, maxItr
    integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0, jiter, iiter, kiter, it
    integer(kind=kint) :: kk, jjiter, ppc, LTRIAL, nget, NEIG
    integer(kind=kint) :: IOUT,IREOR,eITMAX,itype,iS,iE,ic_type,icel,jS,nn
    real(kind=kreal)   :: CTOL, prechk
      REAL(KIND=KREAL) :: PRECHK1,PRECHK2,cchk0,cchk1,cchk,CERR
    real(kind=kreal)   :: t1, t2, aalf, tmp, tmp2, gm, gm2, r1, r2, r3, r4, r5, r6


      REAL(KIND=KREAL), ALLOCATABLE ::  LLDIAG(:), LNDIAG(:), LSUB(:)
      REAL(KIND=KREAL), ALLOCATABLE ::  LZMAT(:,:), LNZMAT(:,:)


    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    allocate(fstrEIG%filter(NPNDOF))
    fstrEIG%filter = 1.0d0

    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0  = hecMESH%node_group%grp_index(ig-1) + 1
      iE0  = hecMESH%node_group%grp_index(ig  )
      it   = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      itS0 = (it - mod(it,10))/10
      itE0 = mod(it,10)

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        do i = itS0,itE0
          fstrEIG%filter((in-1)*NDOF+i) = 0.0d0
        enddo
      enddo
    enddo

    call hecmw_update_m_R(hecMESH, fstrEIG%filter, NP, NDOF)


    in = 0
    do i = 1, NNDOF
      IF(fstrEIG%filter(i) == 1.0d0) in = in + 1
    enddo
    call hecmw_allreduce_I1(hecMESH, in, hecmw_sum)

    eITMAX  = fstrEIG%lczmax
    i = NPNDOF - in
    IF(eITMAX.GT.i) THEN
      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*) '*-------------------------------------------*'
        WRITE(IMSG,*) '  WARNING: LCZMAX exceeds system matrix size. '
        WRITE(IMSG,*) '  Resetting LCZMAX to system matrix size.'
        WRITE(IMSG,*) '*-------------------------------------------*'
      endif
      eITMAX = i
    endif

    CTOL = fstrEIG%lcztol
    NGET = fstrEIG%nget
    neig = eITMAX + NGET
    !maxItr = i
    maxItr = neig - 1

    allocate(lvecq(0:lvecq_size))
    allocate( EM( NPNDOF )              )
    allocate( ewk( NPNDOF, neig)        )
    allocate( eval( neig )              )
    allocate( modal( neig )             )
    allocate( new( neig )               )
    allocate( work( neig*(3*neig + 5) ) )
    allocate( LVECP(NPNDOF)             )
    allocate( LVECPP(NPNDOF)            )
    allocate( lvecq(0)%q(NPNDOF)        )
    allocate( lvecq(1)%q(NPNDOF)        )
    allocate( LWRK(NPNDOF)              )
    allocate( LLWRK(NPNDOF)             )
    allocate( LLLWRK(NPNDOF)            )
    allocate( ALF(NEIG+2)               )
    allocate( BTA(NEIG+2)               )
    ALLOCATE(my_ntotal(0:nprocs) )

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

    my_ntotal(myrank) = NNDOF
      DO I = 0,nprocs-1
        call hecmw_bcast_I1(hecMESH,my_ntotal(I),I)
      enddo

    call SETIVL(fstrEIG%mass,EM,fstrEIG%filter,ewk,LVECP,lvecq(0)%q,lvecq(1)%q,BTA,NPNDOF,&
   &                         neig,my_ntotal,hecMESH,hecMAT,NDOF,NPNDOF)

    do i=1,NPNDOF
      hecMAT%B(i) = EM(i)
    enddo
    hecMAT%X = 0.
    hecMAT%Iarray(98) = 1   !Assmebly complete
    hecMAT%Iarray(97) = 1   !Need numerical factorization
    hecMAT%Rarray(2)  = 1.0

    IF(myrank .EQ. 0) THEN
      WRITE(IMSG,*)
      WRITE(IMSG,*) ' *****   STAGE Begin Lanczos loop     **'
    ENDIF

    DO ITER=1,maxItr
      call DUPL(EM,LVECP,NPNDOF)
      call hecmw_mat_clear_b(hecMAT)

      do i=1,NPNDOF
        hecMAT%B(i) = EM(i)
      enddo
      hecMAT%X = 0.0d0

      call solve_LINEQ( hecMESH,hecMAT,IMSG )

      do i=1,NPNDOF
        EM(i) = hecMAT%X(i)
      enddo

      i = 0
      do ii = 1,NPNDOF
        IF(fstrEIG%filter(ii) .EQ. 0) i = i + 1
      enddo

      do ii = 1,NPNDOF
        EM(ii) = EM(ii)*fstrEIG%filter(ii)
      enddo

      ALLOCATE( lvecq(iter+1)%q(NPNDOF) )

      call SCSHFT(EM,lvecq(iter-1)%q,BTA(ITER),NPNDOF)

      call VECPRO1(AALF,LVECP,EM,NNDOF)
      ALF(ITER)=AALF

      call hecmw_allreduce_R1(hecMESH,ALF(ITER),hecmw_sum)

      call SCSHFT(EM,lvecq(iter)%q,ALF(ITER),NPNDOF)

      LWRK = 0.
      call MATPRO(LWRK,fstrEIG%mass,EM,NPNDOF,1)
      call VECPRO1(prechk,LWRK,LWRK,NNDOF)
        call hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)

      prechk = sqrt(prechk)
      IF(prechk.NE.0.0D0) LWRK = LWRK/prechk
      !IREOR = ITER*(1.0 - fstrEIG%lczrod)
      IREOR = 0

      DO KK = IREOR,ITER
        prechk = 0.0D0
        call VECPRO1(prechk,lvecq(kk)%q,LWRK,NNDOF)
          call hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)

        prechk1 = 0.0D0
        call VECPRO1(prechk1,lvecq(kk)%q,lvecq(kk)%q,NNDOF)
          call hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)

        prechk1 = sqrt(prechk1)
        if(prechk1.ne.0.0D0) prechk = prechk/prechk1
          call MGS1(lvecq(kk)%q,EM,fstrEIG%mass,NPNDOF,my_ntotal,myrank,&
   &                hecMESH,NNDOF)
      enddo

      call MATPRO(LVECPP,fstrEIG%mass,EM,NPNDOF,1)

      call VECPRO1(AALF,LVECPP,EM,NNDOF)
      BTA(ITER+1) = AALF
        call hecmw_allreduce_R1(hecMESH,BTA(ITER+1),hecmw_sum)

      BTA(ITER+1) = SQRT(BTA(ITER+1))

      call hecmw_barrier(hecMESH)
      call UPLCZ(LVECP,lvecq(iter+1)%q,LVECPP,EM,BTA(ITER+1),NPNDOF)

      LTRIAL = ITER
      fstrEIG%iter = ITER

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
   &                      LSUB,LZMAT,LNZMAT,ierr) !Unordered

      DO IITER = 1, LTRIAL
        IF( LLDIAG(IITER).NE.0.0D0 ) THEN
          EVAL(IITER) = 1.0D0/LLDIAG(IITER) + fstrEIG%lczsgm
        ENDIF
      enddo
      call EVSORT(EVAL,NEW,LTRIAL)

      cchk  = 0.0
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

        cchk1 = (BTA(JITER))*(LZMAT(LTRIAL,JITER)/prechk)
        IF(cchk .LT. ABS(cchk1)) THEN
          cchk = ABS(cchk1)
          iiter = jiter
          ppc = prechk1
        ENDIF
      enddo

      IF(cchk.LE.CTOL.AND.ITER.GE.NGET) THEN
        WRITE(IDBG,*) '*=====Desired convergence was obtained =====*'
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
        DO iiter =1,NPNDOF
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

    t2 = hecmw_Wtime() !DEBUG elap
    write(idbg,'(a,f10.2)') 'Lanczos loop (sec) :', T2 - T1 ! elap

    IF(myrank .EQ. 0) THEN
      WRITE(IMSG,*)
      WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
    ENDIF

    DO JITER=1,NGET
      prechk1 = maxval(ewk(:,JITER))
      call hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)

      if(prechk1.NE.0.0D0)then
        do i = 1, NNDOF
          ewk(i,JITER) = ewk(i,JITER)/prechk1
        enddo
      endif
    enddo

  end subroutine fstr_solve_lanczos

end module m_fstr_EIG_lanczos
