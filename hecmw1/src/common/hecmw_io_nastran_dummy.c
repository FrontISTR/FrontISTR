/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/* JP-0 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_io_mesh.h"
#include "hecmw_io_struct.h"
#include "hecmw_struct.h"
#include "hecmw_config.h"
#include "hecmw_system.h"
#include "hecmw_dist.h"
#include "hecmw_dist_print.h"
#include "hecmw_common.h"
#include "hecmw_path.h"
#include "hecmw_conn_conv.h"
#include "hecmw_io_nastran.h"

/* DUMMY */
int HECMW_read_nastran_mesh(const char *filename) {
  fprintf(stdout,
          "##FATAL : HEC-MW IO ERROR : Nastran data is not supported\n");
  fflush(stdout);
  HECMW_abort(HECMW_comm_get_comm());
  return -1;
}

struct hecmwST_local_mesh *HECMW_get_nastran_mesh(const char *filename) {
  fprintf(stdout,
          "##FATAL : HEC-MW IO ERROR : Nastran data is not supported\n");
  fflush(stdout);
  HECMW_abort(HECMW_comm_get_comm());
  return NULL;
}
