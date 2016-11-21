/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_DIST_INIT
#define INC_HECMW_DIST_INIT

#include "hecmw_struct.h"


extern struct hecmwST_local_mesh *HECMW_dist_alloc( void );


extern int HECMW_dist_init( struct hecmwST_local_mesh *mesh );

#endif
