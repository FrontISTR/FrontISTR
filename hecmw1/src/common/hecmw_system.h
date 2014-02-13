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



#ifndef HECMW_SYSTEM_INCLUDED
#define HECMW_SYSTEM_INCLUDED

#include "hecmw_geometric.h"


struct hecmw_system_param {
	double xa;
	double ya;
	double za;
	double xb;
	double yb;
	double zb;
	double xc;
	double yc;
	double zc;
};


extern int HECMW_system(struct hecmw_system_param *param,
						struct hecmw_coord *coord, struct hecmw_coord *result);

#endif
