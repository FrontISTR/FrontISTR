/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_geometric.h"
#include "hecmw_system.h"

int HECMW_system(struct hecmw_system_param *param, struct hecmw_coord *coord,
                 struct hecmw_coord *result) {
  if (param == NULL) {
    /* do nothing */
    *result = *coord;
    return 0;
  }
  if (coord == NULL) return -1;
  if (result == NULL) return -1;

  *result = *coord;

  return 0;
}
