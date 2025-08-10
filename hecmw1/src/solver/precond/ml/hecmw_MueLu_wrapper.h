/*****************************************************************************
 * Copyright (c) 2021 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_MUELU_WRAPPER_INCLUDED
#define HECMW_MUELU_WRAPPER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/*
 * MueLu wrapper functions
 */
void hecmw_MueLu_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr);
void hecmw_MueLu_wrapper_apply(int *id, double rhs[], int *ierr);
void hecmw_MueLu_wrapper_clear(int *id, int *ierr);

/* Fortran interface */
void hecmw_muelu_wrapper_setup_(int *id, int *sym, int *ndof, int *ierr);
void hecmw_muelu_wrapper_setup__(int *id, int *sym, int *ndof, int *ierr);
void HECMW_MUELU_WRAPPER_SETUP(int *id, int *sym, int *ndof, int *ierr);

void hecmw_muelu_wrapper_apply_(int *id, double rhs[], int *ierr);
void hecmw_muelu_wrapper_apply__(int *id, double rhs[], int *ierr);
void HECMW_MUELU_WRAPPER_APPLY(int *id, double rhs[], int *ierr);

void hecmw_muelu_wrapper_clear_(int *id, int *ierr);
void hecmw_muelu_wrapper_clear__(int *id, int *ierr);
void HECMW_MUELU_WRAPPER_CLEAR(int *id, int *ierr);

#ifdef __cplusplus
}
#endif

#endif /* HECMW_MUELU_WRAPPER_INCLUDED */
