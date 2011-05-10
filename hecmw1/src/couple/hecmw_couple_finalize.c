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

#include "hecmw_msgno.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_init.h"
#include "hecmw_couple_info.h"
#include "hecmw_couple_finalize.h"



extern int
HECMW_couple_finalize(char *boundary_id)
{
	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple_finalize(): 'boundary_id' is NULL");
		return -1;
	}

	HECMW_couple_free_init(boundary_id);
	HECMW_couple_free_couple_info();

	return 0;
}
