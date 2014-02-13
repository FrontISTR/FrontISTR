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





#ifndef HECMW_IO_ABAQUS_INCLUDED
#define HECMW_IO_ABAQUS_INCLUDED

#include <stdio.h>
#include "hecmw_struct.h"


extern int HECMW_read_abaqus_mesh(const char *filename);


extern struct hecmwST_local_mesh *HECMW_get_abaqus_mesh(const char *filename);

#endif
