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




#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hecmw_struct.h"
#include "hecmw_lib_fc.h"

#include "hecmw_couple_copy_c2f.h"
#include "hecmw_couple_startup.h"


static struct hecmw_couple_value *couple_value;


/*================================================================================================*/

extern void
hecmw_couple_startup_init_if(char *boundary_id, int *err, int len)
{
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(boundary_id, len, cname, sizeof(cname)) == NULL) return;
	if((couple_value = HECMW_couple_startup(cname)) == NULL) return;
	if(HECMW_couple_copy_c2f_init(couple_value)) return;

	*err = 0;
}



extern void
hecmw_couple_startup_init_if_(char *boundary_id, int *err, int len)
{
	hecmw_couple_startup_init_if(boundary_id, err, len);
}



extern void
hecmw_couple_startup_init_if__(char *boundary_id, int *err, int len)
{
	hecmw_couple_startup_init_if(boundary_id, err, len);
}



extern void
HECMW_COUPLE_STARTUP_INIT_IF(char *boundary_id, int *err, int len)
{
	hecmw_couple_startup_init_if(boundary_id, err, len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_startup_final_if(int *err)
{
	*err = 1;

	if(HECMW_couple_copy_c2f_finalize()) return;

	HECMW_couple_free_couple_value(couple_value);
	couple_value = NULL;

	*err = 0;
}



extern void
hecmw_couple_startup_final_if_(int *err)
{
	hecmw_couple_startup_final_if(err);
}



extern void
hecmw_couple_startup_final_if__(int *err)
{
	hecmw_couple_startup_final_if(err);
}



extern void
HECMW_COUPLE_STARTUP_FINAL_IF(int *err)
{
	hecmw_couple_startup_final_if(err);
}
