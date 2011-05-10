/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#include <stdlib.h>

#include "hecmw_config.h"
#include "hecmw_lib_fc.h"

#include "hecmw_couple_finalize.h"

/*================================================================================================*/

extern void
hecmw_couple_finalize_if(char *boundary_id, int *err, int len)
{
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(boundary_id, len, cname, sizeof(cname)) == NULL) return;
	if(HECMW_couple_finalize(cname)) return;

	*err = 0;
}



extern void
hecmw_couple_finalize_if_(char *boundary_id, int *err, int len)
{
	hecmw_couple_finalize_if(boundary_id, err, len);
}



extern void
hecmw_couple_finalize_if__(char *boundary_id, int *err, int len)
{
	hecmw_couple_finalize_if(boundary_id, err, len);
}



extern void
HECMW_COUPLE_FINALIZE_IF(char *boundary_id, int *err, int len)
{
	hecmw_couple_finalize_if(boundary_id, err, len);
}
