/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_IO_DIST
#define INC_HECMW_IO_DIST

#include "hecmw_struct.h"

extern struct hecmwST_local_mesh *HECMW_get_dist_mesh( char *fname );

extern int HECMW_put_dist_mesh( const struct hecmwST_local_mesh *mesh, char *fname );

#endif
