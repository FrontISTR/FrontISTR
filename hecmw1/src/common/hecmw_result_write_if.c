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
#include "hecmw_result_copy_f2c.h"


static struct hecmwST_result_data *result;



void
hecmw_result_write_st_byname_if(char *name_ID, int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	char *name = NULL;
	char nameid[HECMW_NAME_LEN+1];
	char head[HECMW_HEADER_LEN+1];

	*err = 1;

	if(name_ID) {
		if(HECMW_strcpy_f2c_r(name_ID, name_len, nameid, sizeof(nameid)) == NULL) {
			return;
		}
		name = nameid;
	}
	if(HECMW_strcpy_f2c_r(header, header_len, head, sizeof(head)) == NULL) {
		return;
	}

	if(HECMW_result_write_ST_by_name(name, result, *n_node, *n_elem, *tstep, head)) {
		return;
	}

	*err = 0;
}



void
hecmw_result_write_st_byname_if_(char *name_ID, int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_byname_if(name_ID, n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
hecmw_result_write_st_byname_if__(char *name_ID, int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_byname_if(name_ID, n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
HECMW_RESULT_WRITE_ST_BYNAME_IF(char *name_ID, int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_byname_if(name_ID, n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
hecmw_result_write_st_if(int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_byname_if(NULL, n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
hecmw_result_write_st_if_(int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_if(n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
hecmw_result_write_st_if__(int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_if(n_node, n_elem, tstep, header, err, name_len, header_len);
}



void
HECMW_RESULT_WRITE_ST_IF(int *n_node, int *n_elem, int *tstep, char *header, int *err, int name_len, int header_len)
{
	hecmw_result_write_st_if(n_node, n_elem, tstep, header, err, name_len, header_len);
}


/*----------------------------------------------------------------------------*/


void
hecmw_result_write_st_init_if(int *nnode, int *nelem,  int *err)
{
	result = HECMW_calloc(1, sizeof(*result));
	if(result == NULL) {
		HECMW_set_error(errno, "");
		*err = 1;
		return;
	}
	if(HECMW_result_copy_f2c_init(result, *nnode, *nelem)) {
		*err = 1;
		return;
	}
	*err = 0;
}



void
hecmw_result_write_st_init_if_(int *nnode, int *nelem, int *err)
{
	hecmw_result_write_st_init_if(nnode, nelem, err);
}



void
hecmw_result_write_st_init_if__(int *nnode, int *nelem, int *err)
{
	hecmw_result_write_st_init_if(nnode, nelem, err);
}



void
HECMW_RESULT_WRITE_ST_INIT_IF(int *nnode, int *nelem, int *err)
{
	hecmw_result_write_st_init_if(nnode, nelem, err);
}


/*----------------------------------------------------------------------------*/


void
hecmw_result_w_st_finalize_if(int *err)
{
	*err = 1;

	if(HECMW_result_copy_f2c_finalize()) {
		return;
	}

	HECMW_result_free(result);
	result = NULL;

	*err = 0;
}



void
hecmw_result_w_st_finalize_if_(int *err)
{
	hecmw_result_w_st_finalize_if(err);
}



void
hecmw_result_w_st_finalize_if__(int *err)
{
	hecmw_result_w_st_finalize_if(err);
}



void
HECMW_RESULT_W_ST_FINALIZE_IF(int *err)
{
	hecmw_result_w_st_finalize_if(err);
}

