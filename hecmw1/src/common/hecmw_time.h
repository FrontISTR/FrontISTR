/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_TIME_INCLUDED
#define HECMW_TIME_INCLUDED

#ifndef HECMW_SERIAL
#include "mpi.h"
#else
#include <time.h>
#include <sys/timeb.h>
#endif

double HECMW_Wtime(void);
double HECMW_Wtick(void);

#endif
