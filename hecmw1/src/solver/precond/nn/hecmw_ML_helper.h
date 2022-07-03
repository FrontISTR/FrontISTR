/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_ML_HELPER_INCLUDED
#define HECMW_ML_HELPER_INCLUDED

/*
 * prototype of helper functions in hecmw_ML_helper.f90
 */
extern void hecmw_ml_get_nlocal_(int *id, int *nlocal, int *nlocal_allcolumns, int *ierr);
extern void hecmw_ml_get_coord_(int *id, double x[], double y[], double z[], int *ierr);
extern void hecmw_ml_get_rbm_(int *id, double rbm[], int *ierr);
extern void hecmw_ml_get_loglevel_(int *id, int *level);
extern void hecmw_ml_get_opt_(int *id, int opt[], int *ierr);
extern void hecmw_ml_set_opt_(int *id, int opt[], int *ierr);

#endif
