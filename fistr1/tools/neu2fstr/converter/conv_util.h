/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : FrontSTR Input File Converter                     *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/


#ifndef conv_utilH
#define conv_utilH

#include <stdio.h>

namespace n2h_util {

inline void create_ngrp_name( int id, char* name ) 		{ sprintf( name, "NGRP%d",    id ); }
inline void create_boundary_ngrp_name( int id, char* name )	{ sprintf( name, "BOUNDARY%d",id ); }
inline void create_cload_ngrp_name( int id, char* name )	{ sprintf( name, "CLOAD%d",   id ); }
inline void create_dload_ngrp_name( int id, char* name )	{ sprintf( name, "DLOAD%d",   id ); }
inline void create_egrp_name( int id, char* name )		{ sprintf( name, "EGRP%d",    id ); }
inline void create_sgrp_name( int id, char* name )		{ sprintf( name, "SGRP%d",    id ); }
inline void create_mat_name ( int id, char* name )		{ sprintf( name, "MAT%d",     id ); }
inline void create_egrp_name_for_sec( int id, char* name )	{ sprintf( name, "SECT%d",    id ); }


inline int hec_face_no( int hec_etype, int neu_face, int& fg_front)
{
	fg_front = 1; 
	switch( hec_etype ){
	case 231: case 232: case 241: case 242:
		return (neu_face - 2);
	case 731: case 732: case 741: case 742:
		if( neu_face == 2 ) fg_front = false;
		return 0;
	default:
		return neu_face;
	}
}


}




#endif
