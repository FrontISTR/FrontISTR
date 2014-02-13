/*=====================================================================*
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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

              
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define Table232_Size 6
#define Table342_Size 10
#define Table352_Size 15

/* for HECMW to FSTR */
int Table232[] = { 1,2,3,  6,4,5 };
int Table342[] = { 1,2,3,4, 7,5,6,  8,9,10 };
int Table352[] = { 1,2,3,4,5,6,  9,7,8,  12,10,11,  13,14,15 };
/*int Table232[] = { 1,2,3,  4,5,6 };
int Table342[] = { 1,2,3,4, 5,6,7,  8,9,10 };
int Table352[] = { 1,2,3,4,5,6,  7,8,9,  10,11,12,  13,14,15 };*/


#define MaxItemNumber 15


/*======================================================================*/
/*                                                                      */
/* HECMW to FSTR Mesh Data Converter                                    */
/* Convering Conectivity of Element Type 232, 342 and 352               */
/*                                                                      */
/*======================================================================*/


void hecmw2fstr_connect_conv( int n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	int i, j;
	int start_no;
	int type;
	int table_n;
	int* table;
	int* item;
	int buff[MaxItemNumber];

	for(i=0; i<n_elem; i++ ){
		type = elem_type[i];
		if( type != 232 && type != 342 && type != 352 )
			continue;

		switch(type){
		case 232:
			table_n = Table232_Size;
			table = Table232;
			break;
		case 342:
			table_n = Table342_Size;
			table = Table342;
			break;
		case 352:
			table_n = Table352_Size;
			table = Table352;
			break;
		}
		start_no = elem_node_index[i];
		item = &elem_node_item[start_no];
		memcpy( buff, item, sizeof(int)*table_n );

		for(j=0; j<table_n; j++){
			item[j] = buff[table[j]-1];
		}
	}
}


/*======================================================================*/
/*                                                                      */
/* FSTR to HECMW Mesh Data Converter                                    */
/* Convering Conectivity of Element Type 232, 342 and 352               */
/*                                                                      */
/*======================================================================*/


void fstr2hecmw_connect_conv( int n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	int i, j;
	int start_no;
	int type;
	int table_n;
	int* table;
	int* item;
	int buff[MaxItemNumber];

	for(i=0; i<n_elem; i++ ){
		type = elem_type[i];
		if( type != 232 && type != 342 && type != 352 )
			continue;

		switch(type){
		case 232:
			table_n = Table232_Size;
			table = Table232;
			break;
		case 342:
			table_n = Table342_Size;
			table = Table342;
			break;
		case 352:
			table_n = Table352_Size;
			table = Table352;
			break;
		}
		start_no = elem_node_index[i];
		item = &elem_node_item[start_no];
		memcpy( buff, item, sizeof(int)*table_n );

		for(j=0; j<table_n; j++){
			item[table[j]-1] = buff[j];
		}
	}
}


void fstr2hecmw_elem_conv( int type, int node[] )
{
	int j;
	int table_n;
	int* table;
	int buff[MaxItemNumber];

	if( type != 232 && type != 342 && type != 352 ) return;

	switch(type){
	case 232:
		table_n = Table232_Size;
		table = Table232;
		break;
	case 342:
		table_n = Table342_Size;
		table = Table342;
		break;
	case 352:
		table_n = Table352_Size;
		table = Table352;
		break;
	}
	memcpy( buff, node, sizeof(int)*table_n );

	for(j=0; j<table_n; j++){
		node[table[j]-1] = buff[j];
	}
}



/********************************************************************************************************/
/* Fortran Interface */

void hecmw2fstr_connect_conv_( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	hecmw2fstr_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void hecmw2fstr_connect_conv__( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	hecmw2fstr_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void HECMW2FSTR_CONNECT_CONV( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	hecmw2fstr_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void HECMW2FSTR_CONNECT_CONV_( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	hecmw2fstr_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void HECMW2FSTR_CONNECT_CONV__( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	hecmw2fstr_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

/*---------------------------------------------------------------------------------------------------------*/


void fstr2hecmw_connect_conv_( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	fstr2hecmw_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void fstr2hecmw_connect_conv__( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	fstr2hecmw_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void FSTR2HECMW_CONNECT_CONV( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	fstr2hecmw_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void FSTR2HECMW_CONNECT_CONV_( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	fstr2hecmw_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}

void FSTR2HECMW_CONNECT_CONV__( int* n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] )
{
	fstr2hecmw_connect_conv( *n_elem, elem_type, elem_node_index, elem_node_item );
}
