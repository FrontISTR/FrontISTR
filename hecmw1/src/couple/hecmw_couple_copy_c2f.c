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
#include "hecmw_couple_define.h"
#include "hecmw_couple_startup.h"

static struct hecmw_couple_value *couple_value;


/*================================================================================================*/
static int
set_n(void *dst)
{
	void *src;
	int size;

	src = &couple_value->n;
	size = sizeof(couple_value->n);
	memcpy(dst, src, size);

	return 0;
}


static int
set_item_type(void *dst)
{
	void *src;
	int size;

	src = &couple_value->item_type;
	size = sizeof(couple_value->item_type);
	memcpy(dst, src, size);

	return 0;
}


static int
set_n_dof(void *dst)
{
	void *src;
	int size;

	src = &couple_value->n_dof;
	size = sizeof(couple_value->n_dof);
	memcpy(dst, src, size);

	return 0;
}


static int
set_item(void *dst)
{
	void *src;
	int size;

	if(couple_value->n <= 0) return 0;

	src = couple_value->item;
	if(couple_value->item_type == HECMW_COUPLE_NODE_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n;
	} else if(couple_value->item_type == HECMW_COUPLE_ELEMENT_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n;
	} else if(couple_value->item_type == HECMW_COUPLE_SURFACE_GROUP) {
		size = sizeof(*couple_value->item) * couple_value->n * 2;
	} else {
		return 0;
	}
	memcpy(dst, src, size);

	return 0;
}


static int
set_value(void *dst)
{
	void *src;
	int size;

	if(couple_value->n <= 0 && couple_value->n_dof <= 0) return 0;

	src = couple_value->value;
	size = sizeof(*couple_value->value) * couple_value->n * couple_value->n_dof;
	memcpy(dst, src, size);

	return 0;
}


/*------------------------------------------------------------------------------------------------*/
static int
is_alloc_item(void)
{
	return couple_value->item ? 1 : 0;
}


static int
is_alloc_value(void)
{
	return couple_value->value ? 1 : 0;
}


/*------------------------------------------------------------------------------------------------*/
typedef int (*SetFunc)(void *);
typedef int (*IsAllocatedFunc)(void);

static struct func_table {
	char *struct_name;
	char *var_name;
	SetFunc set_func;
	IsAllocatedFunc is_allocated_func;
} functions[] = {
/* { Struct name, Variable name, memcpy function, check allocation function } */
	{"hecmw_couple_value", "n",         set_n,         NULL          },
	{"hecmw_couple_value", "item_type", set_item_type, NULL          },
	{"hecmw_couple_value", "n_dof",     set_n_dof,     NULL          },
	{"hecmw_couple_value", "item",      set_item,      is_alloc_item },
	{"hecmw_couple_value", "value",     set_value,     is_alloc_value},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);



static IsAllocatedFunc
get_is_allocated_func(char *struct_name, char *var_name)
{
	int i;

	for(i=0; i<NFUNC; i++) {
		if(strcmp(functions[i].struct_name, struct_name) == 0 &&
				strcmp(functions[i].var_name, var_name) == 0) {
			return functions[i].is_allocated_func;
		}
	}
	return NULL;
}



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
HECMW_couple_copy_c2f_init(struct hecmw_couple_value *_couple_value)
{
	couple_value = _couple_value;
	return 0;
}



extern int
HECMW_couple_copy_c2f_finalize(void)
{
	couple_value = NULL;
	return 0;
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_cpl_copy_c2f_isalloc_if(char *struct_name, char *var_name,
		int *is_allocated, int *err, int slen, int vlen)
{
	IsAllocatedFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(couple_value == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_c2f_isalloc_if(): 'couple_value' is not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_isalloc_if(): 'struct_name' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_isalloc_if(): 'var_name' is NULL");
		return;
	}
	if(is_allocated == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_isalloc_if(): 'is_allocated' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, slen, sname, sizeof(sname)) == NULL) return;
	if(HECMW_strcpy_f2c_r(var_name, vlen, vname, sizeof(vname)) == NULL) return;

	if((func = get_is_allocated_func(sname, vname)) == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_c2f_isalloc_if(): IsAllocatedFunc not found");
		return;
	}

	if((*func)()) {
		*is_allocated = 1;
	} else {
		*is_allocated = 0;
	}

	*err = 0;
}



extern void
hecmw_cpl_copy_c2f_isalloc_if_(char *struct_name, char *var_name,
		int *is_allocated, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_isalloc_if(struct_name, var_name, is_allocated, err, slen, vlen);
}



extern void
hecmw_cpl_copy_c2f_isalloc_if__(char *struct_name, char *var_name,
		int *is_allocated, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_isalloc_if(struct_name, var_name, is_allocated, err, slen, vlen);
}



extern void
HECMW_CPL_COPY_C2F_ISALLOC_IF(char *struct_name, char *var_name,
		int *is_allocated, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_isalloc_if(struct_name, var_name, is_allocated, err, slen, vlen);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_cpl_copy_c2f_set_if(char *struct_name, char *var_name,
		void *dst, int *err, int slen, int vlen)
{
	SetFunc func;
	char sname[HECMW_NAME_LEN+1];
	char vname[HECMW_NAME_LEN+1];

	*err = 1;

	if(couple_value == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_c2f_set_if(): 'couple_value' has not initialized yet");
		return;
	}
	if(struct_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_set_if(): 'struct_name' is NULL");
		return;
	}
	if(var_name == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_set_if(): 'var_name' is NULL");
		return;
	}
	if(dst == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"hecmw_cpl_copy_c2f_set_if(): 'dst' is NULL");
		return;
	}

	if(HECMW_strcpy_f2c_r(struct_name, slen, sname, sizeof(sname)) == NULL) return;
	if(HECMW_strcpy_f2c_r(var_name, vlen, vname, sizeof(vname)) == NULL) return;

	if((func = get_set_func(sname, vname)) == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_NULL_PTR,
				"hecmw_cpl_copy_c2f_set_if(): SetFunc not found");
		return;
	}

	if((*func)(dst)) return;

	*err = 0;
}



extern void
hecmw_cpl_copy_c2f_set_if_(char *struct_name, char *var_name,
		void *dst, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_set_if(struct_name, var_name, dst, err, slen, vlen);
}



extern void
hecmw_cpl_copy_c2f_set_if__(char *struct_name, char *var_name,
		void *dst, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_set_if(struct_name, var_name, dst, err, slen, vlen);
}



extern void
HECMW_CPL_COPY_C2F_SET_IF(char *struct_name, char *var_name,
		void *dst, int *err, int slen, int vlen)
{
	hecmw_cpl_copy_c2f_set_if(struct_name, var_name, dst, err, slen, vlen);
}
