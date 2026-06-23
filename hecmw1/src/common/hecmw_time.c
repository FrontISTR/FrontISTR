/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_time.h"

double HECMW_Wtime(void) {
#ifndef HECMW_SERIAL
  return MPI_Wtime();
#else
  #if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
      return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
    }
  #endif
  /* Fallback for systems without CLOCK_MONOTONIC */
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
}

double HECMW_Wtick(void) {
#ifndef HECMW_SERIAL
  return MPI_Wtick();
#else
  #if defined(CLOCK_MONOTONIC)
    struct timespec res;
    if (clock_getres(CLOCK_MONOTONIC, &res) == 0) {
      return (double)res.tv_sec + (double)res.tv_nsec * 1e-9;
    }
  #endif
  return 1e-6;  /* gettimeofday microsecond resolution */
#endif
}

/* interface for fortran */

double hecmw_wtime_fi(void) { return HECMW_Wtime(); }
double hecmw_wtime_fi_(void) { return HECMW_Wtime(); }
double hecmw_wtime_fi__(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI_(void) { return HECMW_Wtime(); }
double HECMW_WTIME_FI__(void) { return HECMW_Wtime(); }

double hecmw_wtick_fi(void) { return HECMW_Wtick(); }
double hecmw_wtick_fi_(void) { return HECMW_Wtick(); }
double hecmw_wtick_fi__(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI_(void) { return HECMW_Wtick(); }
double HECMW_WTICK_FI__(void) { return HECMW_Wtick(); }
