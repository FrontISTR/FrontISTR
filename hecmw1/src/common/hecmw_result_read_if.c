/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
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
#include <errno.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_result.h"
#include "hecmw_result_copy_c2f.h"


static struct hecmwST_result_data *result;



void
hecmw_result_read_initbyname_if(char *name_ID, int *tstep, int *n_node, int *n_elem, int *err, int len)
{
	int nnode, nelem;
	char *name = NULL;
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(name_ID) {
		if(HECMW_strcpy_f2c_r(name_ID, len, cname, sizeof(cname)) == NULL) {
			return;
		}
		name = cname;
	}

	result = HECMW_result_read_by_name(name, *tstep);
	if(result == NULL) {
		return;
	}

	*n_node = nnode = HECMW_result_get_nnode();
	*n_elem = nelem = HECMW_result_get_nelem();
	if(HECMW_result_copy_c2f_init(result, nnode, nelem)) {
		return;
	}

	*err = 0;
}



void
hecmw_result_read_initbyname_if_(char *name_ID, int *tstep, int *n_node, int *n_elem, int *err, int len)
{
	hecmw_result_read_initbyname_if(name_ID, tstep, n_node, n_elem, err, len);
}



void
hecmw_result_read_initbyname_if__(char *name_ID, int *tstep, int *n_node, int *n_elem, int *err, int len)
{
	hecmw_result_read_initbyname_if(name_ID, tstep, n_node, n_elem, err, len);
}



void
HECMW_RESULT_READ_INITBYNAME_IF(char *name_ID, int *tstep, int *n_node, int *n_elem, int *err, int len)
{
	hecmw_result_read_initbyname_if(name_ID, tstep, n_node, n_elem, err, len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_result_read_init_if(int *tstep, int *n_node, int *n_elem, int *err)
{
	hecmw_result_read_initbyname_if(NULL, tstep, n_node, n_elem, err, 0);
}



void
hecmw_result_read_init_if_(int *tstep, int *n_node, int *n_elem, int *err)
{
	hecmw_result_read_init_if(tstep, n_node, n_elem, err);
}



void
hecmw_result_read_init_if__(int *tstep, int *n_node, int *n_elem, int *err)
{
	hecmw_result_read_init_if(tstep, n_node, n_elem, err);
}



void
HECMW_RESULT_READ_INIT_IF(int *tstep, int *n_node, int *n_elem, int *err)
{
	hecmw_result_read_init_if(tstep, n_node, n_elem, err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_result_read_finalize_if(int *err)
{
	*err = 1;

	if(HECMW_result_copy_c2f_finalize()) {
		return;
	}

	HECMW_result_free(result);
	result = NULL;

	*err = 0;
}



void
hecmw_result_read_finalize_if_(int *err)
{
	hecmw_result_read_finalize_if(err);
}



void
hecmw_result_read_finalize_if__(int *err)
{
	hecmw_result_read_finalize_if(err);
}



void
HECMW_RESULT_READ_FINALIZE_IF(int *err)
{
	hecmw_result_read_finalize_if(err);
}

