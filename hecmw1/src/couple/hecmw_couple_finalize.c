/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

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
