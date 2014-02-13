!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by X. Yuan (Advancesoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!C
!C***
!> CONSTRUCT the GLOBAL STIFF MATRIX 
!C***
!C
module m_static_mat_ass_main
   implicit none

   contains

   subroutine FSTR_MAT_ASS_MAIN (fstrSOLID)
      use m_fstr
      type (fstr_solid)         :: fstrSOLID
	  
	  include "HEC_MW3_For.h"

      real(kind=kreal) :: x,y,z,xx(20), yy(20), zz(20), stiffness(20*6, 20*6)
      integer(kind=kint) :: nids(0:20)
      integer(kind=kint) :: ndof, itype, cid, ic_type, nn, j, icel
      integer :: iAss, iPart, iElem, iNode, iGrp, iErr

      ndof = assDOF(1)
	  
      icel = 0	  
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         call mw_select_algebra(0)
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            call mw_matrix_clear( iPart )
            do iElem = 0, mw_get_num_of_element()-1
              icel = icel+1
              call mw_select_element( iElem )
              call  mw_get_element_vert_node_id( nids )
              nn = mw_get_num_of_element_vert()
              ic_type = mw_get_element_type()
              ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
              do iNode = 1, nn
                cid = nids(iNode-1)
                call mw_get_node_coord( cid, x,y,z )
                xx(iNode)=x
                yy(iNode)=y
                zz(iNode)=z
              enddo
              call fstr_local_stf_create(ndof, ic_type, xx, yy, zz, fstrSOLID%elements(icel)%gausses, &
                                     fstrSOLID%elements(icel)%iset, stiffness)
              iErr = mw_matrix_add_elem( iPart, iElem, stiffness(1:nn*ndof,1:nn*ndof) )
            enddo
         enddo
      enddo

   end subroutine FSTR_MAT_ASS_MAIN


   !> Calculate stiff matrix of current element
   subroutine FSTR_LOCAL_STF_CREATE(ndof, ic_type, xx, yy, zz, gausses, iset, stiffness)
      use m_fstr
      use m_static_lib
      use mMechGauss

      integer(kind=kint) :: ndof, ic_type, iset
      real(kind=kreal) :: xx(:), yy(:), zz(:), stiffness(:, :)
      type( tGaussStatus ), intent(in) :: gausses(:)
      real(kind=kreal) :: ee, pp, thick,ecoord(3,20)
      type( tMaterial ), pointer :: material

!** Local variables
      real(kind=kreal) :: local_stf(1830)
      integer(kind=kint) :: nn


      nn = getNumberOfNodes(ic_type)
      ecoord(1,1:nn) = xx(1:nn)
      ecoord(2,1:nn) = yy(1:nn)
      ecoord(3,1:nn) = zz(1:nn)
      material => gausses(1)%pMaterial
      ee = material%variables(M_YOUNGS)
      pp = material%variables(M_POISSON)
      thick = material%variables(M_THICK)
      if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
      if ( ic_type==241 .or. ic_type==242 .or.    &
           ic_type==231 .or. ic_type==232 ) then
        call STF_C2( ic_type,nn,ecoord(1:2,1:nn),gausses(:),thick,stiffness(1:nn*ndof,1:nn*ndof),iset)

      else if (ic_type==361) then
        call STF_C3D8IC( ic_type,nn,ecoord(:,1:nn),gausses(:),stiffness(1:nn*ndof,1:nn*ndof))
      else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.     &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
        call STF_C3( ic_type,nn,ecoord(:,1:nn),gausses(:),stiffness(1:nn*ndof,1:nn*ndof),1.d0)
!
      else if ( ic_type==731) then
        call STF_S3(xx,yy,zz,ee,pp,thick,local_stf)
        call fstr_local_stf_restore(local_stf, nn*ndof, stiffness)
      else if ( ic_type==741) then
        call STF_S4(xx,yy,zz,ee,pp,thick,local_stf)
        call fstr_local_stf_restore(local_stf, nn*ndof, stiffness)
      else
        write(*,*) '###ERROR### : Element type not supported for linear static analysis'
        write(*,*) ' ic_type = ', ic_type
        call hecmw_abort(hecmw_comm_get_comm())
      endif

   end subroutine FSTR_LOCAL_STF_CREATE

   !> restore stiff matrix in 1-dimensional form into 2-dimwnsion array
   subroutine FSTR_LOCAL_STF_RESTORE(local_stf, mat_size, stiffness)
      use m_fstr

      implicit none

      real(kind=kreal) :: local_stf(:), stiffness(:, :)
      integer(kind=kint) :: mat_size

!** Local variables
      integer(kind=kint) :: i, j, num

      do i = 1, mat_size
        do j = 1, i-1
          num = (i-1)*i/2 + j
          stiffness(i, j) = local_stf(num)
          stiffness(j, i) = local_stf(num)
        enddo
        num = i*(i+1)/2
        stiffness(i, i) = local_stf(num)
      enddo
   end subroutine FSTR_LOCAL_STF_RESTORE


end module m_static_mat_ass_main
