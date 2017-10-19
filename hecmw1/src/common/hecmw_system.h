/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_SYSTEM_INCLUDED
#define HECMW_SYSTEM_INCLUDED

#include "hecmw_geometric.h"

struct hecmw_system_param {
  double xa;
  double ya;
  double za;
  double xb;
  double yb;
  double zb;
  double xc;
  double yc;
  double zc;
};

extern int HECMW_system(struct hecmw_system_param *param,
                        struct hecmw_coord *coord, struct hecmw_coord *result);

#endif
