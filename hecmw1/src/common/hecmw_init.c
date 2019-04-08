/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/**
 * @brief I/O and Utility 
 */

#include <stdio.h>
#include "hecmw_util.h"
#include "hecmw_init.h"
/* #include "hecmw_couple_info.h"  2007/12/27 S.Ito   */

int HECMW_init_ex(int *argc, char ***argv, const char *ctrlfile) {
  if (HECMW_comm_init(argc, argv)) return -1;
  HECMW_log(HECMW_LOG_DEBUG, "Initilalizing...");
  if (ctrlfile == NULL) ctrlfile = HECMW_CTRL_FILE;
  if (HECMW_ctrl_init_ex(ctrlfile)) return -1;
  /*     if(HECMW_couple_comm_init() != HECMW_SUCCESS) return -1;  2007/12/27
   * S.Ito */
  return 0;
}

int HECMW_init(int *argc, char ***argv) {
  return HECMW_init_ex(argc, argv, HECMW_CTRL_FILE);
}
