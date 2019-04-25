/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_GEOMETRIC_INCLUDED
#define HECMW_GEOMETRIC_INCLUDED

struct hecmw_coord {
  double x;
  double y;
  double z;
};

extern double HECMW_degree_to_radian(double deg);

extern double HECMW_radian_to_degree(double rad);

extern int HECMW_cylindrical_to_cartesian(const struct hecmw_coord *coord,
                                          struct hecmw_coord *result);

#endif
