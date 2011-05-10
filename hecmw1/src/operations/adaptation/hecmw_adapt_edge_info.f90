!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

!C
!C***
!C*** hecmw_adapt_EDGE_INFO
!C***      

      subroutine hecmw_adapt_EDGE_INFO (hecMESH, nod1, nod2, iedge, NFLAG)
      use  hecmw_util

      integer(kind=4), save :: INITflag, nbuckets
      integer(kind=4), dimension(:), allocatable, save :: ieaddrs
      
      data INITflag/0/
      type      (hecmwST_local_mesh) :: hecMESH

!C
!C-- init.
      if (INITflag.eq.0) then
        INITflag= 1
        nbuckets= 2*max(hecMESH%n_elem,hecMESH%n_node)
        allocate (ieaddrs(-nbuckets:+nbuckets))
        ieaddrs=   0
      endif      

!C
!C    NFALG= 0 : CREATE NEW EDGEs
!C    NFLAG= 1 : REFER  the EDGE INFORMATION
!C    NFLAG= 2 : DEALLOCATE ieaddrs
!C

      if (NFLAG.eq.2) then
        INITflag = 0
        deallocate( ieaddrs )
        return
      endif

      iedge= 0
       
      nn1 = mod(nod1, nbuckets) * mod(nod2, nbuckets)
      iarg= mod( nn1, nbuckets)

      if (NFLAG.eq.0) then
        if (ieaddrs (iarg).gt.hecMESH%n_adapt_edge) then
          ieaddrs (iarg)= 0
        endif
      endif

   50 continue


!C
!C-- NEW EDGE

      if (ieaddrs (iarg).eq.0) then
        hecMESH%n_adapt_edge= hecMESH%n_adapt_edge + 1
        iedge= hecMESH%n_adapt_edge
        hecMESH%adapt_edge_node (2*iedge-1)= nod1
        hecMESH%adapt_edge_node (2*iedge  )= nod2

!      if (iarg.gt.nbuckets) write (*,*) nod1,nod2,iarg
        ieaddrs (iarg)= hecMESH%n_adapt_edge
        return
      else

!      if (iarg.gt.nbuckets) write (*,*) nod1,nod2,iarg
        iedge= ieaddrs (iarg)
        in1= hecMESH%adapt_edge_node (2*iedge-1)
        in2= hecMESH%adapt_edge_node (2*iedge  )

!C
!C-- EXISTING EDGE
        if (in1.eq.nod1 .and. in2.eq.nod2  .or.                         &
     &      in1.eq.nod2 .and. in2.eq.nod1) return

        incr= 1
        ioldadd= iarg
  100   continue
        inewadd= mod (ioldadd + incr**3, nbuckets)

        if (inewadd .eq. ioldadd) then
          icount= icount+ 1
          ioldadd= ioldadd + 1
          inewadd= ioldadd
        endif

        if (NFLAG .eq. 0) then
        if (ieaddrs (inewadd).gt.hecMESH%n_adapt_edge) then
            ieaddrs (inewadd)= 0
            goto 50
        endif
        endif

        if (ieaddrs (inewadd) .ne. 0) then
          iedge= ieaddrs (inewadd)
          in1= hecMESH%adapt_edge_node (2*iedge-1)
          in2= hecMESH%adapt_edge_node (2*iedge  )
!C
!C-- EXISTING EDGE
          if (in1.eq.nod1 .and. in2.eq.nod2  .or.                       &
     &        in1.eq.nod2 .and. in2.eq.nod1) return
          incr= incr + 1
          go to 100

         else
!C
!C-- NEW EDGE
          hecMESH%n_adapt_edge= hecMESH%n_adapt_edge + 1
          iedge= hecMESH%n_adapt_edge
          hecMESH%adapt_edge_node (2*iedge-1)= nod1
          hecMESH%adapt_edge_node (2*iedge  )= nod2

!      if (inewadd.gt.nbuckets) write (*,*) nod1,nod2,inewadd
          ieaddrs (inewadd)= iedge
          return
        endif
      endif

      return
      end


