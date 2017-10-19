/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_vis_subimage_composite_vr.h"

void composite_subimage_vr(int pesize, int *pe_id, int pixn, double *n_subimage,
                           double *n_subopa, double *subimage) {
  int i, j, pe_no;
  double acc_opa, r, g, b;

  for (j = 0; j < pixn; j++) {
    pe_no               = pe_id[0];
    subimage[j * 3 + 0] = n_subimage[pe_no * pixn * 3 + j * 3 + 0];
    subimage[j * 3 + 1] = n_subimage[pe_no * pixn * 3 + j * 3 + 1];
    subimage[j * 3 + 2] = n_subimage[pe_no * pixn * 3 + j * 3 + 2];
    acc_opa             = n_subopa[pe_no * pixn + j];
    i                   = 0;
    while ((acc_opa < 0.99) && (i < pesize - 1)) {
      i++;
      pe_no = pe_id[i];
      r     = n_subimage[pe_no * pixn * 3 + j * 3];
      g     = n_subimage[pe_no * pixn * 3 + j * 3 + 1];
      b     = n_subimage[pe_no * pixn * 3 + j * 3 + 2];
      subimage[j * 3] += r * (1.0 - acc_opa);
      subimage[j * 3 + 1] += g * (1.0 - acc_opa);
      subimage[j * 3 + 2] += b * (1.0 - acc_opa);
      acc_opa += n_subopa[pe_no * pixn + j] * (1.0 - acc_opa);
    }
  }
  return;
}
