/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_RESTART_INCLUDED
#define HECMW_RESTART_INCLUDED

#include <stdio.h>

extern int HECMW_restart_open_by_name(char *name_ID);

extern int HECMW_restart_open(void);

extern int HECMW_restart_close(void);

extern void *HECMW_restart_read(void *addr);

extern int HECMW_restart_add(void *data, size_t size);

extern int HECMW_restart_add_int(int *data, int n_data);

extern int HECMW_restart_add_double(double *data, int n_data);

extern int HECMW_restart_write_by_name(char *name_ID);

extern int HECMW_restart_write(void);

#endif
