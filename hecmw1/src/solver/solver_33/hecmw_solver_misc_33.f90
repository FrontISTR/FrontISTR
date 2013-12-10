module hecmw_solver_misc_33

      contains

!C
!C***
!C*** hecmw_matvec_33
!C***
!C
      subroutine hecmw_matvec_33 (hecMESH, hecMAT, X, Y, COMMtime)
      use hecmw_util
      use m_hecmw_comm_f
      use hecmw_matrix_contact

      implicit none
      real(kind=kreal) :: X(:), Y(:)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix)     :: hecMAT
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: START_TIME, END_TIME
      integer(kind=kint) :: i, j, jS, jE, in
      real(kind=kreal) :: YV1, YV2, YV3, X1, X2, X3

      if (present(COMMtime)) then
      START_TIME= HECMW_WTIME()
      call hecmw_update_3_R (hecMESH, X, hecMAT%NP)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      else
      call hecmw_update_3_R (hecMESH, X, hecMAT%NP)
      endif

      do i= 1, hecMAT%N
        X1= X(3*i-2)
        X2= X(3*i-1)
        X3= X(3*i  )
        YV1= hecMAT%D(9*i-8)*X1 + hecMAT%D(9*i-7)*X2                    &
     &                          + hecMAT%D(9*i-6)*X3
        YV2= hecMAT%D(9*i-5)*X1 + hecMAT%D(9*i-4)*X2                    &
     &                          + hecMAT%D(9*i-3)*X3
        YV3= hecMAT%D(9*i-2)*X1 + hecMAT%D(9*i-1)*X2                    &
     &                          + hecMAT%D(9*i  )*X3

        jS= hecMAT%indexL(i-1) + 1
        jE= hecMAT%indexL(i  )
        do j= jS, jE
          in  = hecMAT%itemL(j)
          X1= X(3*in-2)
          X2= X(3*in-1)
          X3= X(3*in  )
          YV1= YV1 + hecMAT%AL(9*j-8)*X1 + hecMAT%AL(9*j-7)*X2          &
     &                                   + hecMAT%AL(9*j-6)*X3
          YV2= YV2 + hecMAT%AL(9*j-5)*X1 + hecMAT%AL(9*j-4)*X2          &
     &                                   + hecMAT%AL(9*j-3)*X3
          YV3= YV3 + hecMAT%AL(9*j-2)*X1 + hecMAT%AL(9*j-1)*X2          &
     &                                   + hecMAT%AL(9*j  )*X3
        enddo
        jS= hecMAT%indexU(i-1) + 1
        jE= hecMAT%indexU(i  )
        do j= jS, jE
          in  = hecMAT%itemU(j)
          X1= X(3*in-2)
          X2= X(3*in-1)
          X3= X(3*in  )
          YV1= YV1 + hecMAT%AU(9*j-8)*X1 + hecMAT%AU(9*j-7)*X2          &
     &                                   + hecMAT%AU(9*j-6)*X3
          YV2= YV2 + hecMAT%AU(9*j-5)*X1 + hecMAT%AU(9*j-4)*X2          &
     &                                   + hecMAT%AU(9*j-3)*X3
          YV3= YV3 + hecMAT%AU(9*j-2)*X1 + hecMAT%AU(9*j-1)*X2          &
     &                                   + hecMAT%AU(9*j  )*X3
        enddo
        Y(3*i-2)= YV1
        Y(3*i-1)= YV2
        Y(3*i  )= YV3
      enddo

      call hecmw_cmat_multvec_add( hecMAT%cmat, X, Y, hecMAT%NP * hecMAT%NDOF )

      end subroutine hecmw_matvec_33

!C
!C***
!C*** hecmw_matresid_33
!C***
!C
      subroutine hecmw_matresid_33 (hecMESH, hecMAT, X, B, R, COMMtime)
      use hecmw_util

      implicit none
      real(kind=kreal) :: X(:), B(:), R(:)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i

      if (present(COMMtime)) then
      call hecmw_matvec_33 (hecMESH, hecMAT, X, R, COMMtime)
      else
      call hecmw_matvec_33 (hecMESH, hecMAT, X, R)
      endif
      do i = 1, hecMAT%N * 3
        R(i) = B(i) - R(i)
      enddo

      end subroutine hecmw_matresid_33

!C
!C***
!C*** hecmw_rel_resid_L2_33
!C***
!C
      function hecmw_rel_resid_L2_33 (hecMESH, hecMAT, COMMtime)
      use hecmw_util
      use hecmw_solver_misc

      implicit none
      real(kind=kreal) :: hecmw_rel_resid_L2_33
      type ( hecmwST_local_mesh ), intent(in) :: hecMESH
      type ( hecmwST_matrix     ), intent(in) :: hecMAT
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal), allocatable :: r(:)
      real(kind=kreal) :: bnorm2, rnorm2

      allocate(r(hecMAT%NDOF*hecMAT%NP))

      if (present(COMMtime)) then
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, COMMtime)
      call hecmw_matresid_33(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, COMMtime)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, COMMtime)
      else
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2)
      call hecmw_matresid_33(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2)
      endif
      hecmw_rel_resid_L2_33 = sqrt(rnorm2 / bnorm2)

      deallocate(r)
      end function hecmw_rel_resid_L2_33

!C
!C***
!C*** hecmw_Tvec_33
!C***
!C
      subroutine hecmw_Tvec_33 (hecMESH, X, Y, COMMtime)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      real(kind=kreal) :: X(:), Y(:)
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: START_TIME, END_TIME
      integer(kind=kint) :: i, j, jj, k, kk

      if (present(COMMtime)) then
      START_TIME= HECMW_WTIME()
      call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      else
      call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
      endif

      do i= 1, hecMESH%nn_internal * hecMESH%n_dof
        Y(i)= X(i)
      enddo

      do i= 1, hecMESH%mpc%n_mpc
        k = hecMESH%mpc%mpc_index(i-1) + 1
        kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
        Y(kk) = 0.d0
        do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
          jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
          Y(kk) = Y(kk) - hecMESH%mpc%mpc_val(j) * X(jj)
        enddo
      enddo

      end subroutine hecmw_Tvec_33

!C
!C***
!C*** hecmw_Ttvec_33
!C***
!C
      subroutine hecmw_Ttvec_33 (hecMESH, X, Y, COMMtime)
      use hecmw_util
      use m_hecmw_comm_f

      implicit none
      real(kind=kreal) :: X(:), Y(:)
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal) :: START_TIME, END_TIME
      integer(kind=kint) :: i, j, jj, k, kk

      if (present(COMMtime)) then
      START_TIME= HECMW_WTIME()
      call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      else
      call hecmw_update_3_R (hecMESH, X, hecMESH%n_node)
      endif

      do i= 1, hecMESH%nn_internal * hecMESH%n_dof
        Y(i)= X(i)
      enddo

      do i= 1, hecMESH%mpc%n_mpc
        k = hecMESH%mpc%mpc_index(i-1) + 1
        kk = 3 * (hecMESH%mpc%mpc_item(k) - 1) + hecMESH%mpc%mpc_dof(k)
        Y(kk) = 0.d0
        do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
          jj = 3 * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
          Y(jj) = Y(jj) - hecMESH%mpc%mpc_val(j) * X(kk)
        enddo
      enddo

      end subroutine hecmw_Ttvec_33

!C
!C***
!C*** hecmw_TtmatTvec_33
!C***
!C
      subroutine hecmw_TtmatTvec_33 (hecMESH, hecMAT, X, Y, W, COMMtime)
      use hecmw_util

      implicit none
      real(kind=kreal) :: X(:), Y(:), W(:)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix)     :: hecMAT
      real(kind=kreal), optional :: COMMtime

      if (present(COMMtime)) then
      call hecmw_Tvec_33(hecMESH, X, Y, COMMtime)
      call hecmw_matvec_33 (hecMESH, hecMAT, Y, W, COMMtime)
      call hecmw_Ttvec_33(hecMESH, W, Y, COMMtime)
      else
      call hecmw_Tvec_33(hecMESH, X, Y)
      call hecmw_matvec_33 (hecMESH, hecMAT, Y, W)
      call hecmw_Ttvec_33(hecMESH, W, Y)
      endif

      end subroutine hecmw_TtmatTvec_33

!C
!C***
!C*** hecmw_precond_33
!C***
!C
      subroutine hecmw_precond_33(hecMESH, hecMAT, R, Z, ZP, COMMtime)
      use hecmw_util
      use hecmw_matrix_misc

      implicit none
      real(kind=kreal) :: R(:), Z(:), ZP(:)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix)     :: hecMAT
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint ) :: PRECOND, iterPREmax
      integer(kind=kint ) :: i, j, k, iterPRE, isL, ieL, isU, ieU
      real   (kind=kreal) ::X1,X2,X3
      real   (kind=kreal) ::SW1,SW2,SW3

      PRECOND = hecmw_mat_get_precond( hecMAT )
      iterPREmax = hecmw_mat_get_iterpremax( hecMAT )

      if (iterPREmax.le.0) then
        do i= 1, hecMAT%N * 3
          Z(i)= R(i)
        enddo
        return
      endif

!C===
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

!C
!C== Block SSOR
      if (PRECOND.le.2) then

      do i= 1, hecMAT%NP * 3
        ZP(i)= R(i)
      enddo
      do i= 1, hecMAT%NP * 3
        Z(i)= 0.d0
      enddo

      do iterPRE= 1, iterPREmax

!C-- FORWARD

        do i= 1, hecMAT%N
          SW1= ZP(3*i-2)
          SW2= ZP(3*i-1)
          SW3= ZP(3*i  )
          isL= hecMAT%indexL(i-1)+1
          ieL= hecMAT%indexL(i)
          do j= isL, ieL
              k= hecMAT%itemL(j)
             X1= ZP(3*k-2)
             X2= ZP(3*k-1)
             X3= ZP(3*k  )
            SW1= SW1 - hecMAT%AL(9*j-8)*X1 - hecMAT%AL(9*j-7)*X2 - hecMAT%AL(9*j-6)*X3
            SW2= SW2 - hecMAT%AL(9*j-5)*X1 - hecMAT%AL(9*j-4)*X2 - hecMAT%AL(9*j-3)*X3
            SW3= SW3 - hecMAT%AL(9*j-2)*X1 - hecMAT%AL(9*j-1)*X2 - hecMAT%AL(9*j  )*X3
          enddo

          if (hecMAT%cmat%n_val.ne.0) then
            isL= hecMAT%indexCL(i-1)+1
            ieL= hecMAT%indexCL(i)
            do j= isL, ieL
                k= hecMAT%itemCL(j)
               X1= ZP(3*k-2)
               X2= ZP(3*k-1)
               X3= ZP(3*k  )
              SW1= SW1 - hecMAT%CAL(9*j-8)*X1 - hecMAT%CAL(9*j-7)*X2 - hecMAT%CAL(9*j-6)*X3
              SW2= SW2 - hecMAT%CAL(9*j-5)*X1 - hecMAT%CAL(9*j-4)*X2 - hecMAT%CAL(9*j-3)*X3
              SW3= SW3 - hecMAT%CAL(9*j-2)*X1 - hecMAT%CAL(9*j-1)*X2 - hecMAT%CAL(9*j  )*X3
            enddo
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - hecMAT%ALU(9*i-5)*X1
          X3= X3 - hecMAT%ALU(9*i-2)*X1 - hecMAT%ALU(9*i-1)*X2
          X3= hecMAT%ALU(9*i  )*  X3
          X2= hecMAT%ALU(9*i-4)*( X2 - hecMAT%ALU(9*i-3)*X3 )
          X1= hecMAT%ALU(9*i-8)*( X1 - hecMAT%ALU(9*i-6)*X3 - hecMAT%ALU(9*i-7)*X2)
          ZP(3*i-2)= X1
          ZP(3*i-1)= X2
          ZP(3*i  )= X3
        enddo

!C-- BACKWARD

        do i= hecMAT%N, 1, -1
          isU= hecMAT%indexU(i-1) + 1
          ieU= hecMAT%indexU(i)
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= ieU, isU, -1
              k= hecMAT%itemU(j)
             X1= ZP(3*k-2)
             X2= ZP(3*k-1)
             X3= ZP(3*k  )
            SW1= SW1 + hecMAT%AU(9*j-8)*X1 + hecMAT%AU(9*j-7)*X2 + hecMAT%AU(9*j-6)*X3
            SW2= SW2 + hecMAT%AU(9*j-5)*X1 + hecMAT%AU(9*j-4)*X2 + hecMAT%AU(9*j-3)*X3
            SW3= SW3 + hecMAT%AU(9*j-2)*X1 + hecMAT%AU(9*j-1)*X2 + hecMAT%AU(9*j  )*X3
          enddo

          if (hecMAT%cmat%n_val.gt.0) then
            isU= hecMAT%indexCU(i-1) + 1
            ieU= hecMAT%indexCU(i)
            do j= ieU, isU, -1
                k= hecMAT%itemCU(j)
               X1= ZP(3*k-2)
               X2= ZP(3*k-1)
               X3= ZP(3*k  )
              SW1= SW1 + hecMAT%CAU(9*j-8)*X1 + hecMAT%CAU(9*j-7)*X2 + hecMAT%CAU(9*j-6)*X3
              SW2= SW2 + hecMAT%CAU(9*j-5)*X1 + hecMAT%CAU(9*j-4)*X2 + hecMAT%CAU(9*j-3)*X3
              SW3= SW3 + hecMAT%CAU(9*j-2)*X1 + hecMAT%CAU(9*j-1)*X2 + hecMAT%CAU(9*j  )*X3
            enddo
          endif

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - hecMAT%ALU(9*i-5)*X1
          X3= X3 - hecMAT%ALU(9*i-2)*X1 - hecMAT%ALU(9*i-1)*X2
          X3= hecMAT%ALU(9*i  )*  X3
          X2= hecMAT%ALU(9*i-4)*( X2 - hecMAT%ALU(9*i-3)*X3 )
          X1= hecMAT%ALU(9*i-8)*( X1 - hecMAT%ALU(9*i-6)*X3 - hecMAT%ALU(9*i-7)*X2)
          ZP(3*i-2)=  ZP(3*i-2) - X1
          ZP(3*i-1)=  ZP(3*i-1) - X2
          ZP(3*i  )=  ZP(3*i  ) - X3
        enddo
!C
!C-- additive Schwartz
!C
        do i= 1, hecMAT%N * 3
          Z(i)= Z(i) + ZP(i)
        enddo

        if (iterPRE.eq.iterPREmax) exit

!C--    {ZP} = {R} - [A] {Z}

        if (present(COMMtime)) then
        call hecmw_matresid_33 (hecMESH, hecMAT, Z, R, ZP, COMMtime)
        else
        call hecmw_matresid_33 (hecMESH, hecMAT, Z, R, ZP)
        endif

      enddo
      endif

!C
!C== Block SCALING
      if (PRECOND.eq.3) then
      do i= 1, hecMAT%N
        X1= R(3*i-2)
        X2= R(3*i-1)
        X3= R(3*i  )
        X2= X2 - hecMAT%ALU(9*i-5)*X1
        X3= X3 - hecMAT%ALU(9*i-2)*X1 - hecMAT%ALU(9*i-1)*X2
        X3= hecMAT%ALU(9*i  )*  X3
        X2= hecMAT%ALU(9*i-4)*( X2 - hecMAT%ALU(9*i-3)*X3 )
        X1= hecMAT%ALU(9*i-8)*( X1 - hecMAT%ALU(9*i-6)*X3 - hecMAT%ALU(9*i-7)*X2)
        Z(3*i-2)= X1
        Z(3*i-1)= X2
        Z(3*i  )= X3
      enddo
      endif

      end subroutine hecmw_precond_33

!C
!C***
!C*** hecmw_mpc_scale
!C***
!C
      subroutine hecmw_mpc_scale(hecMESH)
      use hecmw_util

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: i, j, k
      real(kind=kreal) :: WVAL

      do i = 1, hecMESH%mpc%n_mpc
        k = hecMESH%mpc%mpc_index(i-1)+1
        WVAL = 1.d0 / hecMESH%mpc%mpc_val(k)
        hecMESH%mpc%mpc_val(k) = 1.d0
        do j = hecMESH%mpc%mpc_index(i-1)+2, hecMESH%mpc%mpc_index(i)
          hecMESH%mpc%mpc_val(j) = hecMESH%mpc%mpc_val(j) * WVAL
        enddo
        hecMESH%mpc%mpc_const(i) = hecMESH%mpc%mpc_const(i) * WVAL
      enddo

      end subroutine hecmw_mpc_scale

!C
!C***
!C*** hecmw_trans_b_33
!C***
!C
      subroutine hecmw_trans_b_33(hecMESH, hecMAT, B, BT, W, COMMtime)
      use hecmw_util

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix)     :: hecMAT
      real(kind=kreal) :: B(:), W(:)
      real(kind=kreal), target :: BT(:)
      real(kind=kreal), optional :: COMMtime

      real(kind=kreal), pointer :: XG(:)
      integer(kind=kint) :: i, k, kk

!C===
!C +---------------------------+
!C | {bt}= [T']({b} - [A]{xg}) |
!C +---------------------------+
!C===
      XG => BT
      XG = 0.d0

!C-- Generate {xg} from mpc_const
      do i = 1, hecMESH%mpc%n_mpc
        k = hecMESH%mpc%mpc_index(i-1) + 1
        kk = 3 * hecMESH%mpc%mpc_item(k) + hecMESH%mpc%mpc_dof(k) - 3
        XG(kk) = hecMESH%mpc%mpc_const(i)
      enddo

!C-- {w} = {b} - [A]{xg}
      if (present(COMMtime)) then
      call hecmw_matresid_33 (hecMESH, hecMAT, XG, B, W, COMMtime)
      else
      call hecmw_matresid_33 (hecMESH, hecMAT, XG, B, W)
      endif

!C-- {bt} = [T'] {w}
      if (present(COMMtime)) then
      call hecmw_Ttvec_33(hecMESH, W, BT, COMMtime)
      else
      call hecmw_Ttvec_33(hecMESH, W, BT)
      endif

      end subroutine hecmw_trans_b_33

!C
!C***
!C*** hecmw_tback_x_33
!C***
!C
      subroutine hecmw_tback_x_33(hecMESH, X, W, COMMtime)
      use hecmw_util

      implicit none
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal) :: X(:), W(:)
      real(kind=kreal), optional :: COMMtime

      integer(kind=kint) :: i, k, kk

!C-- {tx} = [T]{x}
      if (present(COMMtime)) then
      call hecmw_Tvec_33(hecMESH, X, W, COMMtime)
      else
      call hecmw_Tvec_33(hecMESH, X, W)
      endif

!C-- {x} = {tx} + {xg}

      do i= 1, hecMESH%nn_internal * 3
        X(i)= W(i)
      enddo

      do i = 1, hecMESH%mpc%n_mpc
        k = hecMESH%mpc%mpc_index(i-1) + 1
        kk = 3 * hecMESH%mpc%mpc_item(k) + hecMESH%mpc%mpc_dof(k) - 3
        X(kk) = X(kk) + hecMESH%mpc%mpc_const(i)
      enddo

      end subroutine hecmw_tback_x_33

end module hecmw_solver_misc_33
