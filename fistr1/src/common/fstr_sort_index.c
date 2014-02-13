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



#include <stdlib.h>


static int my_comp( const void* a, const void* b)
{
	return (*(int*)a - *(int*)b);
}


void c_fstr_sort_index( int* index_data, int n)
{
	qsort( index_data, n, sizeof(int)*2, my_comp);
}



/*----------- Fortran Interface ---------------*/

void fstr_sort_index( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

void fstr_sort_index_( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

void fstr_sort_index__( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

void FSTR_SORT_INDEX( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

void FSTR_SORT_INDEX_( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

void FSTR_SORT_INDEX__( int* index_data, int* n)
{
	c_fstr_sort_index( index_data, *n);
}

