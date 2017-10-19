/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdlib.h>

#include "hecmw_config.h"
#include "hecmw_lib_fc.h"

#include "hecmw_couple_finalize.h"

/*================================================================================================*/

extern void hecmw_couple_finalize_if(char *boundary_id, int *err, int len) {
  char cname[HECMW_NAME_LEN + 1];

  *err = 1;

  if (HECMW_strcpy_f2c_r(boundary_id, len, cname, sizeof(cname)) == NULL)
    return;
  if (HECMW_couple_finalize(cname)) return;

  *err = 0;
}

extern void hecmw_couple_finalize_if_(char *boundary_id, int *err, int len) {
  hecmw_couple_finalize_if(boundary_id, err, len);
}

extern void hecmw_couple_finalize_if__(char *boundary_id, int *err, int len) {
  hecmw_couple_finalize_if(boundary_id, err, len);
}

extern void HECMW_COUPLE_FINALIZE_IF(char *boundary_id, int *err, int len) {
  hecmw_couple_finalize_if(boundary_id, err, len);
}
