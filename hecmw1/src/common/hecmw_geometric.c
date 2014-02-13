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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_geometric.h"


#define HECMW_PI 3.1415926535897932

/* Degree to Radian Transformation */
double
HECMW_degree_to_radian(double deg)
{
	return deg * HECMW_PI / 180;
}


/* Radian to Degree Transformation */
double
HECMW_radian_to_degree(double rad)
{
	return rad * 180 / HECMW_PI;
}


int
HECMW_cylindrical_to_cartesian(
		const struct hecmw_coord *coord, struct hecmw_coord *result)
{
	double r;
	double rtheta;
	double z;

	if(result == NULL) return -1;

	r      = coord->x;
	rtheta = coord->y;	/* radian */
	z      = coord->z;

	result->x = r * cos(rtheta);
	result->y = r * sin(rtheta);
	result->z = z;

	return 0;
}

