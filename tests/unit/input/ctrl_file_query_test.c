/*****************************************************************************
 * Copyright (c) 2026 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_control.h"

static int write_file(const char *name, const char *text) {
  FILE *fp = fopen(name, "w");
  if (fp == NULL) return 1;
  fputs(text, fp);
  return fclose(fp) != 0;
}

int main(int argc, char **argv) {
  int failed = 0;
  if (argc != 2) return EXIT_FAILURE;
  if (HECMW_ctrl_is_control_file(NULL) ||
      HECMW_ctrl_is_control_file(argv[1])) return EXIT_FAILURE;
  if (write_file("ctrl_query.inp", "** input\n") ||
      write_file("ctrl_other.cnt", "!END\n") ||
      write_file("ctrl_query.dat", "!CONTROL, NAME=analysis_input\n./ctrl_query.inp\n") ||
      write_file("ctrl_other.dat", "!CONTROL, NAME=other_input\nctrl_other.cnt\n"))
    return EXIT_FAILURE;
  if (HECMW_ctrl_init_ex("ctrl_query.dat")) return EXIT_FAILURE;
  failed |= !HECMW_ctrl_is_control_file("ctrl_query.inp");
  failed |= !HECMW_ctrl_is_control_file("./ctrl_query.inp");
  failed |= !HECMW_ctrl_is_control_file(argv[1]);
  failed |= HECMW_ctrl_is_control_file("ctrl_other.cnt");
  failed |= HECMW_ctrl_is_control_file("missing_input.inp");
  failed |= HECMW_ctrl_is_control_file(NULL);
  HECMW_ctrl_finalize();
  failed |= HECMW_ctrl_is_control_file(argv[1]);
  if (HECMW_ctrl_init_ex("ctrl_other.dat")) return EXIT_FAILURE;
  failed |= HECMW_ctrl_is_control_file(argv[1]);
  failed |= !HECMW_ctrl_is_control_file("ctrl_other.cnt");
  HECMW_ctrl_finalize();
  remove("ctrl_query.inp");
  remove("ctrl_other.cnt");
  remove("ctrl_query.dat");
  remove("ctrl_other.dat");
  if (failed) fputs("Control file identity query failed\n", stderr);
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
