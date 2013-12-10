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

      subroutine hecmw_adapt_active (hecMESH)
!C
!C***
!C*** hecmw_adapt_active
!C***
!C
!C    FIND active NODEs, EDGEs, ELEMs
!C
      use  hecmw_util
      type (hecmwST_local_mesh) :: hecMESH
!C
!C-- arrays
      allocate (hecMESH%adapt_edge_home(hecMESH%n_adapt_edge))
      allocate (hecMESH%adapt_act_edge (hecMESH%n_adapt_edge))

      hecMESH%adapt_edge_home= hecMESH%my_rank
      hecMESH%adapt_act_edge = 0

!C
!C-- EDGE home
      icou= 0
      do ie= 1, hecMESH%n_adapt_edge
        in1= hecMESH%adapt_edge_node(2*ie-1)
        in2= hecMESH%adapt_edge_node(2*ie  )
        ip1= hecMESH%node_ID(2*in1)
        ip2= hecMESH%node_ID(2*in2)
        hecMESH%adapt_edge_home(ie)= min(ip1,ip2)

        if ((ip1.eq.hecMESH%my_rank.and.ip2.eq.hecMESH%my_rank).or.     &
     &      (ip1.eq.hecMESH%my_rank.and.ip1.lt.ip2            ).or.     &     
     &      (ip2.eq.hecMESH%my_rank.and.ip2.lt.ip1    )) then
                 icou = icou + 1
          hecMESH%adapt_act_edge(icou)= ie
        endif
      enddo
      hecMESH%n_adapt_act_edge= icou

      return
      end






