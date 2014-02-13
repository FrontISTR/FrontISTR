/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef HECMW_CTRLLEX_INCLUDED
#define HECMW_CTRLLEX_INCLUDED

#include <stdio.h>

enum {
	HECMW_CTRLLEX_NL = 1000,
	HECMW_CTRLLEX_INT,
	HECMW_CTRLLEX_DOUBLE,
	HECMW_CTRLLEX_NAME,
	HECMW_CTRLLEX_FILENAME,

	HECMW_CTRLLEX_H_CONTROL = 2000,
	HECMW_CTRLLEX_H_MESH,
	HECMW_CTRLLEX_H_MESH_GROUP,
	HECMW_CTRLLEX_H_RESULT,
	HECMW_CTRLLEX_H_RESTART,
	HECMW_CTRLLEX_H_SUBDIR,

	HECMW_CTRLLEX_K_ABAQUS = 3000,
	HECMW_CTRLLEX_K_DIR,
	HECMW_CTRLLEX_K_FEMAP,
	HECMW_CTRLLEX_K_GEOFEM,
	HECMW_CTRLLEX_K_HECMW_DIST,
	HECMW_CTRLLEX_K_HECMW_ENTIRE,
	HECMW_CTRLLEX_K_IN,
	HECMW_CTRLLEX_K_INOUT,
	HECMW_CTRLLEX_K_IO,
	HECMW_CTRLLEX_K_LIMIT,
	HECMW_CTRLLEX_K_NAME,
	HECMW_CTRLLEX_K_NASTRAN,
	HECMW_CTRLLEX_K_ON,
	HECMW_CTRLLEX_K_OFF,
	HECMW_CTRLLEX_K_OUT,
	HECMW_CTRLLEX_K_REFINE,
	HECMW_CTRLLEX_K_TYPE,
};

extern double HECMW_ctrllex_get_number(void);

extern char *HECMW_ctrllex_get_text(void);

extern int HECMW_ctrllex_get_lineno(void);

extern int HECMW_ctrllex_next_token(void);

extern int HECMW_ctrllex_next_token_skip(int skip_token);

extern int HECMW_ctrllex_set_input(FILE *fp);

extern int HECMW_ctrllex_skip_line(void);

extern int HECMW_ctrllex_unput_token(void);

#endif
