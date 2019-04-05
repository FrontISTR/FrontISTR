/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_repart.h"
void HECMW_dlb_memory_exit(char *var) {
  fprintf(stderr,
          "#### HEC-MW-VIS-E0001:There is no enough memory allocated for "
          "variable %s\n",
          var);
  HECMW_Finalize();
  exit(0);
}

void HECMW_dlb_print_exit(char *var) {
  fprintf(stderr, "%s\n", var);
  HECMW_Finalize();
  exit(0);
}
