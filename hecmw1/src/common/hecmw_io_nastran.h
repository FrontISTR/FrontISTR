/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

/* debug vertion */
/*
  2004.06.07
*/

#ifndef HECMW_IO_NASTRAN_INCLUDED
#define HECMW_IO_NASTRAN_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "hecmw_util.h"

extern int HECMW_read_nastran_mesh(const char *filename);

extern struct hecmwST_local_mesh *HECMW_get_nastran_mesh(const char *filename);

#endif

