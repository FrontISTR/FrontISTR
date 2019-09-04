/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include "hecmw_finalize.h"
#include "hecmw_util.h"

int HECMW_finalize(void) {
  HECMW_log(HECMW_LOG_DEBUG, "Finalizing...");

  HECMW_ctrl_finalize();

#ifndef HECMW_SERIAL
  MPI_Finalize();

  HECMW_log(HECMW_LOG_DEBUG, "MPI finalized");
#endif

  return 0;
}
