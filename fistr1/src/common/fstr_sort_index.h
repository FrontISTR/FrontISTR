/*=====================================================================*
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : I/O and Utility                                   !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
 *=====================================================================*/



/*
	Index Sorting for FSTR
	2004.10.18 by N.Imai
	-------------------------
	[Fortran]
	integer(kind=4) :: index_data(2,:), n
	call fstr_sort_index( index_data, n )
*/

#ifndef fstr_sort_indexH
#define fstr_sort_indexH

void c_fstr_sort_index( int* index_data, int n);

#endif
