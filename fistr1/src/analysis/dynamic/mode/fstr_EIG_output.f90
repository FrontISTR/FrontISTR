!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Display results in .out files
module m_fstr_EIG_output
contains

!C======================================================================
!C----------------------------------------------------------------------
      subroutine fstr_eigen_output(hecMESH,hecMAT,fstrEIG)
!C----------------------------------------------------------------------
!C*
!C* SHOW RESULTS
!C*
      use m_fstr
      use lczeigen
      use m_fstr_EIG_lanczos_util
      use m_fstr_EIG_matmult
      implicit none
      integer(kind=kint) :: I,IREOR,JJITER,JITER,IITER
      real(kind=kreal) :: prechk0(1)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      type (fstr_solid       )  :: fstrSOLID
      type (hecmwST_result_data):: fstrRESULT
      type (fstr_eigen) :: fstrEIG

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

      CALL EGLIST(hecMESH,hecMAT,fstrEIG)
      IF(myrank.EQ.0) THEN
        WRITE(IMSG,*) ''
        WRITE(IMSG,*) '*----------------------------------------------*'
        WRITE(IMSG,*) '*--E I G E N V A L U E  C O N V E R G E N C E--*'
        WRITE(IMSG,*) '*----------------------------------------------*'
        WRITE(IMSG,*) 'Absolute residual =  |(||Kx - lambda*Mx||)|'
      ENDIF
!C
!Normalize eigenvectors
      IF(myrank==0) THEN
        WRITE(IMSG,*) '*------------- Residual Check ----------------*'
        WRITE(IMSG,*) '   Iter.#   Eigenvalue    Abs. Residual  '
        WRITE(IMSG,*) '   *-----*  *---------*  *--------------*'
      ENDIF
!C
      LLLWRK = 0.0D0
      DO JJITER = 1,fstrEIG%nget
        JITER = NEW(JJITER)

        if( modal(JITER).eq.1 ) then
          LWRK = 0.0D0
!!!          CALL VECPRO( prechk1,ewk(1,JJITER),ewk(1,JJITER),Gntotal,1 )
          prechk1=0.0
          do i = 1, Gntotal
            prechk1=prechk1+ewk(i,JJITER)**2
          enddo
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1( hecMESH,prechk1,hecmw_sum )
          end if
          prechk1 = sqrt(prechk1)
          if( prechk1.NE.0.0D0 ) ewk(:,JJITER) = ewk(:,JJITER)/prechk1
!C
          DO IITER = 1, ntotal
            LWRK(IITER) = ewk(IITER,JJITER)
          ENDDO
!C
          IF(NDOF.EQ.3) THEN
            CALL MATMULT3( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ELSE IF(NDOF.EQ.2) THEN
            CALL MATMULT2( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ELSE IF(NDOF.EQ.6) THEN
            CALL MATMULT6( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ENDIF
!C
          LLWRK = 0.0D0
!!!          CALL MATPRO(LLWRK,mass,ewk(1,JJITER),ntotal,1)
          do i = 1, ntotal
            LLWRK(i) = mass(i)*ewk(i,JJITER)
          enddo
!C
          CCHK1 = 0.0D0
          DO IITER = 1,Gntotal
            CCHK1 = CCHK1 + ( LWRK(IITER) - (eval(JITER) &
     &                      - fstrEIG%lczsgm)*LLWRK(IITER) )**2
          END DO
!C
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1(hecMESH,CCHK1,hecmw_sum)
          end if
          CCHK1 = SQRT(CCHK1)
          IF(myrank.EQ.0) THEN
             WRITE(IMSG,'(2X,I5,2X,5E15.6)') JITER,eval(JITER),CCHK1
          ENDIF
        endif
      ENDDO
!C
      IF(myrank==0) THEN
        WRITE(IMSG,*) ''
        WRITE(IMSG,*)'*    ---END Eigenvalue listing---     *'
        WRITE(IMSG,*) ''
      ENDIF

      end subroutine fstr_eigen_output

      subroutine fstr_eigen_make_result(hecMESH, hecMAT, fstrEIG, fstrRESULT, ewk)
        use m_fstr
        use m_hecmw2fstr_mesh_conv
        use hecmw_util
        implicit none
        type (hecmwST_local_mesh ) :: hecMESH
        type (hecmwST_matrix     ) :: hecMAT
        type (fstr_eigen         ) :: fstrEIG
        type (hecmwST_result_data) :: fstrRESULT

        integer(kind=kint) :: i, istep, nget, NP, NDOF, NPNDOF, totalmpc, MPC_METHOD
        real(kind=kreal)   :: t1, ewk(:,:)
        real(kind=kreal), allocatable :: X(:)
        character(len=HECMW_HEADER_LEN) :: header
        character(len=HECMW_NAME_LEN)   :: label
        character(len=HECMW_NAME_LEN)   :: nameID

        nget   = fstrEIG%nget
        NDOF   = hecMAT%NDOF
        NPNDOF = hecMAT%NP*hecMAT%NDOF
        !totalmpc = hecMESH%mpc%n_mpc
        !call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

        allocate(X(NPNDOF))
        X = 0.0d0

        do istep = 1, nget
          do i=1,NPNDOF
            !X(i) = fstrEIG%evec(i,istep)
            X(i) = ewk(i,istep)
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

end module m_fstr_EIG_output


