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



/*======================================================================*/
/*                                                                      */
/* HECMW to FSTR Mesh Data Converter                                    */
/* Convering Conectivity of Element Type 232, 342 and 352               */
/*                                                                      */
/*======================================================================*/

void c_hecmw2fstr_connect_conv( int n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] );

/*======================================================================*/
/*                                                                      */
/* FSTR to HECMW Mesh Data Converter                                    */
/* Convering Conectivity of Element Type 232, 342 and 352               */
/*                                                                      */
/*======================================================================*/

void c_fstr2hecmw_connect_conv( int n_elem, int elem_type[], int elem_node_index[], int elem_node_item[] );
void c_fstr2hecmw_elem_conv( int elem_type, int node[] );
