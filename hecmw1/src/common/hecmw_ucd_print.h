/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_UCD_PRINT
#define INC_HECMW_UCD_PRINT

#include "hecmw_struct.h"
#include "hecmw_result.h"

extern int
HECMW_ucd_print( const struct hecmwST_local_mesh *mesh,
                 const struct hecmwST_result_data *result,
                 const char *ofname );
extern int
HECMW_ucd_legacy_print( const struct hecmwST_local_mesh *mesh,
                        const struct hecmwST_result_data *result,
                        const char *ofname );

#endif

