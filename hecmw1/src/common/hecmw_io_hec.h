/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_IO_HEC_INCLUDED
#define HECMW_IO_HEC_INCLUDED

#include <stdio.h>
#include "hecmw_struct.h"

extern int HECMW_read_entire_mesh(const char *filename);

extern struct hecmwST_local_mesh *HECMW_get_entire_mesh(const char *filename);

#endif
