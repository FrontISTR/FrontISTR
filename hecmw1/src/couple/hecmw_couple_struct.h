/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_STRUCT
#define INC_HECMW_COUPLE_STRUCT

#include "hecmw_config.h"
#include "hecmw_struct.h"

struct hecmw_couple_comm {
  int psize;
  int rank;
  int *ranks;
  HECMW_Comm comm;
  HECMW_Group group;
  int root;
  int is_root;
  int is_member;
};

#endif /* INC_HECMW_COUPLE_STRUCT */
