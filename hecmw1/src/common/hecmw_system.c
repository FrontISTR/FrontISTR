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
#include "hecmw_system.h"

int
HECMW_system(struct hecmw_system_param *param,
			struct hecmw_coord *coord, struct hecmw_coord *result)
{
	if(param == NULL) {
		/* do nothing */
		*result = *coord;
		return 0;
	}
	if(coord == NULL) return -1;
	if(result == NULL) return -1;

	*result = *coord;

	return 0;
}
