!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : m_addContactStiffness                             !
!                                                                      !
!            Written by Z. Sun(ASTOM)                                  !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief This module provides functions:
!!  1) obtain contact stiffness matrix of each contact pair and assemble
!!     it into global stiffness matrix.
!!  2) obtain contact nodal force vector of each contact pair and assemble
!!      it into right-hand side vector to update non-equilibrated nodal force vector.
!!  3) Modify Lagrange multiplier-related part of stiffness matrix and right-hand side
!!     vector for dealing with prescribed displacement boundary condition.
!
!>  \author     Z. Sun(ASTOM)
!>  \date       2010/11
!>  \version    0.00
!!
!======================================================================!
module m_addContactStiffness

use m_fstr
use elementInfo
use m_contact_lib
use fstr_matrix_con_contact
use hecmw_matrix_ass
use m_fstr_Residual

implicit none

public :: fstr_AddContactStiffness
public :: fstr_Update_NDForce_contact
public :: update_NDForce_contact
public :: fstr_ass_load_contact
public :: fstr_mat_ass_bc_contact

private :: getContactStiffness
private :: fstr_mat_ass_contact
private :: getContactNodalForce
private :: getTrialFricForceANDcheckFricState

contains

!> \brief This subroutine obtains contact stiffness matrix of each contact pair
!!  and assembles it into global stiffness matrix.
   subroutine fstr_AddContactStiffness(cstep,iter,hecMAT,fstrMAT,fstrSOLID)

      integer (kind=kint)                    :: cstep                                 !< current loding step
      integer (kind=kint)                    :: iter
      type(hecmwST_matrix)                    :: hecMAT                               !< type hecmwST_matrix
      type(fstr_solid)                        :: fstrSOLID                            !< type fstr_solid
      type(fstrST_matrix_contact_lagrange)    :: fstrMAT                              !< type fstrST_matrix_contact_lagrange
      integer (kind=kint)                    :: ctsurf, etype, nnode, ndLocal(9)      !< contants of type tContact
      integer (kind=kint)                    :: i, j, id_lagrange, grpid
      real(kind=kreal)                        :: lagrange
      real(kind=kreal)                        :: stiffness(9*3+1,9*3+1)

      fstrMAT%AL_lagrange = 0.0d0
      fstrMAT%AU_lagrange = 0.0d0

      id_lagrange = 0

      do i = 1, size(fstrSOLID%contacts)

        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

        do j = 1, size(fstrSOLID%contacts(i)%slave)

          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

          id_lagrange = id_lagrange + 1
          lagrange = fstrMAT%Lagrange(id_lagrange)

          ctsurf = fstrSOLID%contacts(i)%states(j)%surface
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)

! Obtain contact stiffness matrix of contact pair
          call getContactStiffness(iter,etype,nnode,fstrSOLID%contacts(i)%states(j),  &
                                   fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,stiffness)

! Assemble contact stiffness matrix of contact pair into global stiffness matrix
          call fstr_mat_ass_contact(nnode,ndLocal,id_lagrange,fstrSOLID%contacts(i)%fcoeff,stiffness,hecMAT,fstrMAT)

        enddo

      enddo


   end subroutine fstr_AddContactStiffness


!> \brief This subroutine obtains contact stiffness matrix of contact pair
   subroutine getContactStiffness(iter,etype,nnode,ctState,tPenalty,fcoeff,lagrange,stiffness)

      type(tContactState)   :: ctState                                        !< type tContactState
      integer (kind=kint)   :: iter
      integer (kind=kint)   :: etype, nnode                                   !< type of master segment; number of nodes of master segment
      integer (kind=kint)   :: i, j, k, l
      real(kind=kreal)      :: normal(3), shapefunc(nnode)                     !< normal vector at target point; shape functions
      real(kind=kreal)      :: nTm((nnode+1)*3)                                !< vector
      real(kind=kreal)      :: fcoeff, tPenalty                                !< friction coefficient; tangential penalty
      real(kind=kreal)      :: lagrange                                        !< value of Lagrange multiplier
      real(kind=kreal)      :: tf_trial(3), length_tft
      real(kind=kreal)      :: tangent(3), tTm((nnode+1)*3)
      real(kind=kreal)      :: stiffness(9*3+1,9*3+1)                          !< contact stiffness matrix

      stiffness = 0.0d0

      call getShapeFunc( etype, ctState%lpos(:), shapefunc )

      normal(1:3) = ctState%direction(1:3)

      nTm(1:3) = normal(1:3)
      do i = 1, nnode
        nTm(i*3+1:i*3+3) = -shapefunc(i)*normal(1:3)
      enddo

      i = (nnode+1)*3 + 1
      do j = 1, (nnode+1)*3
        stiffness(i,j) = nTm(j);  stiffness(j,i) = nTm(j)
      enddo


     if( fcoeff /= 0.0d0 ) then
       if( lagrange>0.0d0 .or. iter==1 ) then

          do i = 1, nnode+1
            do j = 1, i
              do k = 1, 3
                do l = 1, k
                  stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - tPenalty*nTm((i-1)*3+k)*nTm((j-1)*3+l)
                  if( k==l ) then
                    if(i==1 .and. j==1)then
                      stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) + tPenalty
                    elseif(i>1 .and. j==1)then
                      stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - tPenalty*shapefunc(i-1)
                    elseif(i>1 .and. j>1)then
                      stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) + tPenalty*shapefunc(i-1)*shapefunc(j-1)
                    endif
                  endif
                  if(i==j .and. k==l) cycle
                  stiffness((j-1)*3+l,(i-1)*3+k) = stiffness((i-1)*3+k,(j-1)*3+l)
                enddo
              enddo
            enddo
          enddo

          if( ctstate%state == contactSlip ) then

            tf_trial(1:3) = ctstate%tangentForce_trial(1:3)
            length_tft = dsqrt(dot_product(tf_trial,tf_trial))
            tangent(1:3) = tf_trial(1:3)/length_tft

            tTm(1:3) = -tangent(1:3)
            do i = 1, nnode
              tTm(i*3+1:i*3+3) = shapefunc(i)*tangent(1:3)
            enddo

            do i = 1, nnode+1
              do j = 1, nnode+1
                do k = 1, 3
                  do l = 1, 3
                    stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l)       &
                                                   + tPenalty*(-tTm((i-1)*3+k)*tTm((j-1)*3+l)  &
                                                               +tTm((i-1)*3+k)*nTm((j-1)*3+l)*dot_product(tangent,normal))
                  enddo
                enddo
              enddo
            enddo
            stiffness(1:(nnode+1)*3,1:(nnode+1)*3) = (fcoeff*lagrange/length_tft)*stiffness(1:(nnode+1)*3,1:(nnode+1)*3)

!
!        do i = 1, (nnode + 1)*3
!          do j = 1, i
!            do k = 1, 3
!              do l = 1, k
!                stiffness((i-1)*3+k,(j-1)*3+l) = stiffness((i-1)*3+k,(j-1)*3+l) - mut*tTm((i-1)*3+k)*tTm((j-1)*3+l))
!                if( k/=l )stiffness((j-1)*3+l,(i-1)*3+k) = stiffness((i-1)*3+k,(j-1)*3+l)
!              enddo !- l
!            enddo !- k
!          enddo !- j
!        enddo !- i
!        stiffness(1:(nnode+1)*3,1:(nnode+1)*3) = (fcoeff*dabs(lagrange)/length_tft)*stiffness(1:(nnode+1)*3,1:(nnode+1)*3)

!         j = (nnode+1)*3 + 1
!         do i = 1, (nnode+1)*3
!           stiffness(i,j) = stiffness(i,j) - fcoeff*tTm(i)
!         enddo
         stiffness(1:(nnode+1)*3,(nnode+1)*3+1) = stiffness(1:(nnode+1)*3,(nnode+1)*3+1) - fcoeff*tTm(1:(nnode+1)*3)

         endif
       endif
     endif

   end subroutine getContactStiffness


!> \brief This subroutine assembles contact stiffness matrix of a contact pair into global stiffness matrix
   subroutine fstr_mat_ass_contact(nnode,ndLocal,id_lagrange,fcoeff,stiffness,hecMAT,fstrMAT)

     type(hecmwST_matrix)                    :: hecMAT                                     !< type hecmwST_matrix
     type(fstrST_matrix_contact_lagrange)    :: fstrMAT                                    !< type fstrST_matrix_contact_lagrange
     integer (kind=kint)                     :: nnode, ndLocal(nnode+1),id_lagrange        !< total number of nodes of master segment
                                                                                            !< global number of nodes of contact pair
                                                                                            !< number of Lagrange multiplier
     integer (kind=kint)                     :: i, j, inod, jnod, l
     integer (kind=kint)                     :: isL, ieL, idxL_base, kL, idxL, isU, ieU, idxU_base, kU, idxU
     real(kind=kreal)                         :: fcoeff                                    !< friction coefficient
     real(kind=kreal)                         :: stiffness(9*3+1,9*3+1)                    !< contact stiffness matrix
     real(kind=kreal)                         :: a(3,3)

     i = nnode + 1 + 1
     inod = id_lagrange
     isL = fstrMAT%indexL_lagrange(inod-1)+1
     ieL = fstrMAT%indexL_lagrange(inod)

     do j = 1, nnode + 1
       jnod = ndLocal(j)
       isU = fstrMAT%indexU_lagrange(jnod-1)+1
       ieU = fstrMAT%indexU_lagrange(jnod)

       kL = hecmw_array_search_i(fstrMAT%itemL_lagrange,isL,ieL,jnod)
       if( kL<isL .or. kL>ieL ) then
         write(*,*) '###ERROR### : cannot find connectivity (Lagrange1)'
         stop
       endif
       kU = hecmw_array_search_i(fstrMAT%itemU_lagrange,isU,ieU,inod)
       if( kU<isU .or. kU>ieU ) then
         write(*,*) '###ERROR### : cannot find connectivity (Lagrange2)'
         stop
       endif

       idxL_base = (kL-1)*3
       idxU_base = (kU-1)*3

       do l = 1, 3
         idxL = idxL_base + l
         fstrMAT%AL_lagrange(idxL) = fstrMAT%AL_lagrange(idxL) + stiffness((i-1)*3+1,(j-1)*3+l)
         idxU = idxU_base + l
         fstrMAT%AU_lagrange(idxU) = fstrMAT%AU_lagrange(idxU) + stiffness((j-1)*3+l,(i-1)*3+1)
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

   end subroutine fstr_mat_ass_contact


!> \brief This subroutine obtains contact nodal force vector of each contact pair
!! and assembles it into right-hand side vector to update non-equilibrated nodal force vector.
   subroutine fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID,conMAT)

      type(hecmwST_local_mesh)                :: hecMESH                                  !< type hecmwST_local_mesh
      type(hecmwST_matrix)                    :: hecMAT                                   !< type hecmwST_matrix
      type(fstr_solid)                        :: fstrSOLID                                !< type fstr_solid
      type(fstrST_matrix_contact_lagrange)    :: fstrMAT                                  !< type fstrST_matrix_contact_lagrange
      type(hecmwST_matrix),optional           :: conMAT                                   !< type hecmwST_matrix for contact part only
      integer (kind=kint)                    :: ctsurf, etype, nnode, ndLocal(9)          !< contants of type tContact
      integer (kind=kint)                    :: i, j, k, id_lagrange
      real(kind=kreal)                        :: ndCoord(9*3), ndDu(9*3)                  !< nodal coordinates ; nodal displacement increment
      real(kind=kreal)                        :: lagrange                                 !< value of Lagrange multiplier
      real(kind=kreal)                        :: ctForce(9*3+1)                           !< nodal contact force vector

      integer(kind=kint)                     :: cstep                                    !< current calculation step
      integer(kind=kint)                     :: ig0, grpid, ig, ityp, is0, ie0, ik, in, idof1, idof2, idof, ndof
      real(kind=kreal)                        :: rhs                                      !< value of prescribed displacement

      id_lagrange = 0

      do i = 1, size(fstrSOLID%contacts)

        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

        do j = 1, size(fstrSOLID%contacts(i)%slave)

          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

          id_lagrange = id_lagrange + 1
          lagrange = fstrMAT%Lagrange(id_lagrange)

          fstrSOLID%contacts(i)%states(j)%multiplier(1) = fstrMAT%Lagrange(id_lagrange)

          ctsurf = fstrSOLID%contacts(i)%states(j)%surface
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)
          do k = 1, nnode+1
            ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)          &
                                         + fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)       &
                                         + fstrSOLID%dunode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
            ndDu((k-1)*3+1:(k-1)*3+3) = fstrSOLID%dunode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          enddo
! Obtain contact nodal force vector of contact pair
          call getContactNodalForce(etype,nnode,ndCoord,ndDu,fstrSOLID%contacts(i)%states(j),    &
                                    fstrSOLID%contacts(i)%tPenalty,fstrSOLID%contacts(i)%fcoeff,lagrange,ctForce)
! Update non-eqilibrited force vector
          if(present(conMAT)) then
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,ctForce,conMAT)
          else
            call update_NDForce_contact(nnode,ndLocal,id_lagrange,ctForce,hecMAT)
          endif

        enddo

      enddo

!    Consider SPC condition
      call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, hecMAT%B)
      if(present(conMAT)) call fstr_Update_NDForce_SPC(cstep, hecMESH, fstrSOLID, conMAT%B)

   end subroutine fstr_Update_NDForce_contact

!> \brief This subroutine obtains contact nodal force vector of contact pair
   subroutine getContactNodalForce(etype,nnode,ndCoord,ndDu,ctState,tPenalty,fcoeff,lagrange,ctForce)

      type(tContactState)   :: ctState                                    !< type tContactState
      integer (kind=kint)   :: etype, nnode                               !< type of master segment; number of nodes of master segment
      integer (kind=kint)   :: i, j
      real(kind=kreal)       :: fcoeff, tPenalty                           !< friction coefficient; tangential penalty
      real(kind=kreal)       :: lagrange                                   !< value of Lagrange multiplier
      real(kind=kreal)       :: ndCoord((nnode+1)*3), ndDu((nnode+1)*3)    !< nodal coordinates; nodal displacement increment
      real(kind=kreal)       :: normal(3), shapefunc(nnode)                !< normal vector at target point; shape functions
      real(kind=kreal)       :: nTm((nnode+1)*3)                           !< vector
      real(kind=kreal)       :: tf_trial(3), length_tft, tangent(3), tf_final(3)
      real(kind=kreal)       :: ctForce((nnode+1)*3+1)                     !< contact force vector

      ctForce = 0.0d0

      call getShapeFunc( etype, ctState%lpos(:), shapefunc )

      normal(1:3) = ctState%direction(1:3)

      nTm(1:3) = -normal(1:3)
      do i = 1, nnode
        nTm(i*3+1:i*3+3) = shapefunc(i)*normal(1:3)
      enddo

      do j = 1, (nnode+1)*3
        ctForce(j) = lagrange*nTm(j)
      enddo
      j = (nnode+1)*3 + 1
      ctForce(j) = dot_product(nTm,ndCoord)

      if(fcoeff /= 0.0d0 .and. lagrange > 0.0d0)then

        call getTrialFricForceANDcheckFricState(nnode,tPenalty,fcoeff,lagrange,normal,shapefunc,nTm,ndDu,ctstate)

        if( ctstate%state == contactStick ) then
          tf_final(1:3) = ctstate%tangentForce_trial(1:3)
        elseif( ctstate%state == contactSlip ) then
          tf_trial(1:3) = ctstate%tangentForce_trial(1:3)
          length_tft = dsqrt(dot_product(tf_trial,tf_trial))
          tangent(1:3) = tf_trial(1:3)/length_tft
          tf_final(1:3) = fcoeff*dabs(lagrange)*tangent(1:3)
        endif

        ctForce(1:3) = ctForce(1:3) - tf_final(1:3)
        do j = 1, nnode
          ctForce(j*3+1:j*3+3) = ctForce(j*3+1:j*3+3) + shapefunc(j)*tf_final(1:3)
        enddo

        ctstate%tangentForce_final(1:3) = tf_final(1:3)

      endif

   end subroutine getContactNodalForce


!> \brief This subroutine calculates trial friction force and checks friction state
   subroutine getTrialFricForceANDcheckFricState(nnode,tPenalty,fcoeff,lagrange,normal,shapefunc,nTm,ndDu,ctstate)

     type(tContactState)    :: ctState                                        !< type tContactState
     integer (kind=kint)   :: nnode                                           !< number of nodes of master segment
     integer (kind=kint)   :: i, j
     real(kind=kreal)       :: fcoeff, tPenalty                                !< friction coefficient; tangential penalty
     real(kind=kreal)       :: lagrange                                        !< value of Lagrange multiplier
     real(kind=kreal)       :: ndDu((nnode+1)*3)                               !< nodal displacement increment
     real(kind=kreal)       :: normal(3), shapefunc(nnode)                     !< normal vector at target point; shape functions
     real(kind=kreal)       :: nTm((nnode+1)*3)                                !< vector
     real(kind=kreal)       :: dotP
     real(kind=kreal)       :: relativeDisp(3)                                 !< relative displacement
     real(kind=kreal)       :: tf_yield

     relativeDisp = 0.0d0

     dotP = dot_product(nTm,ndDu)
     do i = 1, 3
       relativeDisp(i) = - ndDu(i)
       do j = 1, nnode
         relativeDisp(i) = relativeDisp(i) + shapefunc(j)*ndDu(j*3+i)
       enddo
       relativeDisp(i) = relativeDisp(i) - dotP*normal(i)
       ctstate%tangentForce_trial(i) = ctstate%tangentForce(i) -tPenalty*relativeDisp(i)
     enddo

     tf_yield = fcoeff*dabs(lagrange)
     if(ctstate%state == contactSlip) tf_yield =0.99d0*tf_yield
     if( dsqrt(dot_product(ctstate%tangentForce_trial,ctstate%tangentForce_trial)) <= tf_yield ) then
       ctstate%state = contactStick
     else
       ctstate%state = contactSlip
     endif

   end subroutine getTrialFricForceANDcheckFricState


!> \brief This subroutine assembles contact nodal force vector into right-hand side vector
!! to update non-equilibrated nodal force vector.
   subroutine update_NDForce_contact(nnode,ndLocal,id_lagrange,ctForce,hecMAT)

     type(hecmwST_matrix)                    :: hecMAT                          !< type hecmwST_matrix
     type(fstrST_matrix_contact_lagrange)    :: fstrMAT                         !< type fstrST_matrix_contact_lagrange
     integer (kind=kint)                     :: nnode, ndLocal(nnode+1)         !< number of nodes of master segment
                                                                                 !< global number of nodes of contact pair
     integer (kind=kint)                     :: id_lagrange                     !< number of Lagrange multiplier
     integer (kind=kint)                     :: np, ndof                        !< total number of nodes; degree of freedom
     integer (kind=kint)                     :: i, j, inod
     real(kind=kreal)                         :: ctForce(:)                      !< nodal contact force vector

     np = hecMAT%NP ; ndof = hecMAT%NDOF

     do i = 1, nnode + 1
       inod = ndLocal(i)
       hecMAT%B((inod-1)*3+1:(inod-1)*3+3) = hecMAT%B((inod-1)*3+1:(inod-1)*3+3) + ctForce((i-1)*3+1:(i-1)*3+3)
     enddo

     hecMAT%B(np*ndof+id_lagrange) = ctForce((nnode+1)*3+1)

   end subroutine update_NDForce_contact


!> \brief This subroutine adds initial contact force vector to the right-hand side vector
!> \at the beginning of each substep calculation
   subroutine fstr_ass_load_contact(cstep, hecMESH, hecMAT, fstrSOLID, fstrMAT)

      type(hecmwST_local_mesh)                :: hecMESH                                  !< type hecmwST_local_mesh
      type(hecmwST_matrix)                    :: hecMAT                                   !< type hecmwST_matrix
      type(fstr_solid)                        :: fstrSOLID                                !< type fstr_solid
      type(fstrST_matrix_contact_lagrange)    :: fstrMAT                                  !< type fstrST_matrix_contact_lagrange
      integer (kind=kint)                    :: cstep                                     !< current step
      integer (kind=kint)                    :: np, ndof                                  !< total number of nodes; degree of freedom
      integer (kind=kint)                    :: i, j, k, l, id_lagrange, lnod, grpid
      integer (kind=kint)                    :: ctsurf, etype, nnode, ndLocal(9)          !< contants of type tContact
      real(kind=kreal)                        :: ndCoord(9*3), lagrange                    !< nodal coordinates; value of Lagrange mutiplier
      real(kind=kreal)                        :: normal(3), shapefunc(9)                   !< normal vector; shape functions
      real(kind=kreal)                        :: nTm(10*3)                                 !< vector
      real(kind=kreal)                        :: tf_final(3)                               !< final friciton force vector
      real(kind=kreal)                        :: ctForce(9*3+1)                            !< initial nodal contact force vector

      np = hecMAT%NP ; ndof = hecMAT%NDOF

      id_lagrange = 0

      do i = 1, size(fstrSOLID%contacts)

        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

        do j = 1, size(fstrSOLID%contacts(i)%slave)

          if( fstrSOLID%contacts(i)%states(j)%state == CONTACTFREE ) cycle

          id_lagrange = id_lagrange + 1
          lagrange = fstrSOLID%contacts(i)%states(j)%multiplier(1)

          ctsurf = fstrSOLID%contacts(i)%states(j)%surface
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          ndLocal(2:nnode+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(1:nnode)
          do k = 1, nnode+1
            ndCoord((k-1)*3+1:(k-1)*3+3) = hecMESH%node((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)          &
                                         + fstrSOLID%unode((ndLocal(k)-1)*3+1:(ndLocal(k)-1)*3+3)
          enddo

          ctForce = 0.0d0

          call getShapeFunc( etype, fstrSOLID%contacts(i)%states(j)%lpos(:), shapefunc )
          normal(1:3) = fstrSOLID%contacts(i)%states(j)%direction(1:3)
          nTm(1:3) = -normal(1:3)
          do k = 1, nnode
            nTm(k*3+1:k*3+3) = shapefunc(k)*normal(1:3)
          enddo
          do l = 1, (nnode+1)*3
            ctForce(l) = lagrange*nTm(l)
          enddo
          l = (nnode+1)*3 + 1
          ctForce(l) = dot_product(nTm(1:(nnode+1)*3),ndCoord(1:(nnode+1)*3))

          if( fstrSOLID%contacts(i)%fcoeff/=0.0d0 .and. lagrange>0.0d0 )then
            tf_final(1:3) = fstrSOLID%contacts(i)%states(j)%tangentForce_final(1:3)
            ctForce(1:3) = ctForce(1:3) - tf_final(1:3)
            do l = 1, nnode
              ctForce(l*3+1:l*3+3) = ctForce(l*3+1:l*3+3) + shapefunc(l)*tf_final(1:3)
            enddo
          endif

          do l = 1, nnode + 1
            lnod = ndLocal(l)
            hecMAT%B((lnod-1)*3+1:(lnod-1)*3+3) = hecMAT%B((lnod-1)*3+1:(lnod-1)*3+3) + ctForce((l-1)*3+1:(l-1)*3+3)
          enddo
          hecMAT%B(np*ndof+id_lagrange) = ctForce((nnode+1)*3+1)

        enddo

      enddo

   end subroutine fstr_ass_load_contact


!> Modify Lagrange multiplier-related part of stiffness matrix and right-hand side vector
!> for dealing with prescribed displacement boundary condition
    subroutine fstr_mat_ass_bc_contact(hecMAT,fstrMAT,inode,idof,RHS)

      type(hecmwST_matrix)                   :: hecMAT                            !< hecmwST_matrix
      type(fstrST_matrix_contact_lagrange)   :: fstrMAT                           !< fstrST_matrix_contact_lagrange
      integer(kind=kint)    :: inode, idof                                        !< number of node; degree of freedom
      integer(kind=kint)    :: isU, ieU, isL, ieL, i, l, k
      real(kind=kreal)       :: RHS                                                !< value of prescribed displacement

      isU = fstrMAT%indexU_lagrange(inode-1)+1
      ieU = fstrMAT%indexU_lagrange(inode)
      do i = isU, ieU
        fstrMAT%AU_lagrange((i-1)*3+idof) = 0.0d0
        l = fstrMAT%itemU_lagrange(i)
        isL = fstrMAT%indexL_lagrange(l-1)+1
        ieL = fstrMAT%indexL_lagrange(l)
        k = hecmw_array_search_i(fstrMAT%itemL_lagrange,isL,ieL,inode)
        if(k < isL .or. k > ieL) cycle
        hecMAT%B(hecMAT%NP*hecMAT%NDOF+l) = hecMAT%B(hecMAT%NP*hecMAT%NDOF+l) - fstrMAT%AL_lagrange((k-1)*3+idof)*RHS
        fstrMAT%AL_lagrange((k-1)*3+idof) = 0.0d0
      enddo

    end subroutine fstr_mat_ass_bc_contact


end module m_addContactStiffness



