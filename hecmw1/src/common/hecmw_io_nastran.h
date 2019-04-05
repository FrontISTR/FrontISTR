/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_IO_NASTRAN_INCLUDED
#define HECMW_IO_NASTRAN_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "hecmw_util.h"

extern int HECMW_read_nastran_mesh(const char *filename);

extern struct hecmwST_local_mesh *HECMW_get_nastran_mesh(const char *filename);

#endif
