/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_IO_ABAQUS_INCLUDED
#define HECMW_IO_ABAQUS_INCLUDED

#include <stdio.h>
#include "hecmw_struct.h"


extern int HECMW_read_abaqus_mesh(const char *filename);


extern struct hecmwST_local_mesh *HECMW_get_abaqus_mesh(const char *filename);

#endif
