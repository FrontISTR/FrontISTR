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
#include <errno.h>

#include "hecmw_malloc.h"
#include "hecmw_lib_fc.h"
#include "hecmw_msgno.h"
#include "hecmw_struct.h"
#include "hecmw_couple_define.h"
#include "hecmw_couple_startup.h"

static struct hecmw_couple_value *couple_value;


/*------------------------------------------------------------------------------------------------*/
/*
 * SetFunc
 */

static int
set_n(void *src)
{
	couple_value->n = *((int *)src);
	return 0;
}


static int
set_item_type(void *src)
{
	couple_value->item_type = *((int *)src);
	return 0;
}


static int
set_n_dof(void *src)
{
	couple_value->n_dof = *((int *)src);
	return 0;
}


static int
set_item(void *src)
{
	int size;

	if(couple_value->n <= 0) return 0;

	if(couple_value->item_type == HECMW_COUPLE_NODE_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n;
	} else if(couple_value->item_type == HECMW_COUPLE_ELEMENT_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n;
	} else if(couple_value->item_type == HECMW_COUPLE_SURFACE_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n * 2;
	} else {
		return 0;
	}
	couple_value->item = HECMW_malloc(size);
	if(couple_value->item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(couple_value->item, src, size);
	return 0;
}


static int
set_value(void *src)
{
	int size;

	if(couple_value->n <= 0 || couple_value->n_dof <= 0) return 0;

	size = sizeof(*couple_value->value) * couple_value->n * couple_value->n_dof;
	couple_value->value = HECMW_malloc(size);
	if(couple_value->value == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	memcpy(couple_value->value, src, size);
	return 0;
}


/*------------------------------------------------------------------------------------------------*/
/*
 * SetFunc table
 */

typedef int (*SetFunc)(void *);

static struct func_table {
	char *struct_name;
	char *var_name;
	SetFunc set_func;
} functions[] = {
/*  { Struct_name, Variable name, memcpy function } */
	{"hecmw_couple_value", "n",         set_n},
	{"hecmw_couple_value", "item_type", set_item_type},
	{"hecmw_couple_value", "n_dof",     set_n_dof},
	{"hecmw_couple_value", "item",      set_item},
	{"hecmw_couple_value", "value",     set_value},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);



static SetFunc
get_set_func(char *struct_name, char *var_name)
{
	int i;

	for(i=0; i<NFUNC; i++) {
		if(strcmp(functions[i].struct_name, struct_name) == 0 &&
				strcmp(functions[i].var_name, var_name) == 0) {
			return functions[i].set_func;
		}
	}
	return NULL;
}


/*------------------------------------------------------------------------------------------------*/

extern int
HECMW_couple_copy_f2c_init(struct hecmw_couple_value *_couple_value)
{
	couple_value = _couple_value;
	return 0;
}



extern int
HECMW_couple_copy_f2c_finalize(void)
{
	couple_value = NULL;
	return 0;
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_cpl_copy_f2c_set_if(char *struct_name, char *var_name,
		void *src, int *err, int slen, int vlen)
{
	SetFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(couple_value == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_f2c_set_if(): 'couple_value' has not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_f2c_set_if(): 'struct_name' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_f2c_set_if(): 'var_name' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, slen, sname, sizeof(sname)) == NULL) return;
	if(HECMW_strcpy_f2c_r(var_name, vlen, vname, sizeof(vname)) == NULL) return;

	if((func = get_set_func(sname, vname)) == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_f2c_set_if(): SetFunc not found");
		return;
	}

	if((*func)(src)) return;

	*err = 0;
}



extern void
hecmw_cpl_copy_f2c_set_if_(char *struct_name, char *var_name,
		void *src, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}



extern void
hecmw_cpl_copy_f2c_set_if__(char *struct_name, char *var_name,
		void *src, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}



extern void
HECMW_CPL_COPY_F2C_SET_IF(char *struct_name, char *var_name,
		void *src, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_f2c_set_if(struct_name, var_name, src, err, slen, vlen);
}
