/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hecmw_struct.h"
#include "hecmw_lib_fc.h"

#include "hecmw_couple_copy_c2f.h"
#include "hecmw_couple_copy_f2c.h"
#include "hecmw_couple_startup.h"

static struct hecmw_couple_value *couple_value;


/*================================================================================================*/

extern void
hecmw_couple_exec_init_if(int *err)
{
	*err = 1;

	if((couple_value = HECMW_couple_alloc_couple_value()) == NULL) return;
	if(HECMW_couple_copy_f2c_init(couple_value)) return;

	*err = 0;
}



extern void
hecmw_couple_exec_init_if_(int *err)
{
	hecmw_couple_exec_init_if(err);
}



extern void
hecmw_couple_exec_init_if__(int *err)
{
	hecmw_couple_exec_init_if(err);
}



extern void
HECMW_COUPLE_EXEC_INIT_IF(int *err)
{
	hecmw_couple_exec_init_if(err);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_if(char *boundary_id, int *err, int len)
{
	char cname[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(boundary_id, len, cname, sizeof(cname)) == NULL) return;
	if(HECMW_couple(cname, couple_value) != HECMW_SUCCESS) return;
	if(HECMW_couple_copy_c2f_init(couple_value)) return;

	*err = 0;
}



extern void
hecmw_couple_if_(char *boundary_id, int *err, int len)
{
	hecmw_couple_if(boundary_id, err, len);
}



extern void
hecmw_couple_if__(char *boundary_id, int *err, int len)
{
	hecmw_couple_if(boundary_id, err, len);
}



extern void
HECMW_COUPLE_IF(char *boundary_id, int *err, int len)
{
	hecmw_couple_if(boundary_id, err, len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_exec_finalize_if(int *err)
{
	*err = 1;

	if(HECMW_couple_copy_f2c_finalize()) return;
	HECMW_couple_free_couple_value(couple_value);
	couple_value = NULL;

	*err = 0;
}



extern void
hecmw_couple_exec_finalize_if_(int *err)
{
	hecmw_couple_exec_finalize_if(err);
}



extern void
hecmw_couple_exec_finalize_if__(int *err)
{
	hecmw_couple_exec_finalize_if(err);
}



extern void
HECMW_COUPLE_EXEC_FINALIZE_IF(int *err)
{
	hecmw_couple_exec_finalize_if(err);
}
