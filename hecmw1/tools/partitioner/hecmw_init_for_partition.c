/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "hecmw_msgno.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"

#include "hecmw_part_define.h"
#include "hecmw_part_get_control.h"
#include "hecmw_init_for_partition.h"

#define DEFAULT_CONTROL_FILE_NAME "hecmw_part_ctrl.dat"

static void print_usage(void) {
  fprintf(stderr, "Usage: hecmw_part1 [-f filename] [-v] [-d number] [-m KMETIS|PMETIS|RCB ] [-e number] \n");
  fprintf(stderr, "         [ -t NODE-BASED|ELEMENT-BASED ] [ -u filename ] [ -c DEFAULT|AGGREGATE|DISTRIBUTE|SIMPLE]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -f  specify control file name\n");
  fprintf(stderr, "  -v  print verbose messages\n");
  fprintf(stderr, "  -h  print usage\n");
  fprintf(stderr, "*** Following option is set parameter instead of hecmw_part_ctrl.dat ***\n");  
  fprintf(stderr, "  -d  number of sub-domains \n");
  fprintf(stderr, "  -t  partitioning type            (unimplemented, default: NODE-BASED) \n");
  fprintf(stderr, "  -m  partitioning method          (unimplemented, default: KMETIS) \n");
  fprintf(stderr, "  -e  depth of overlapping zone    (unimplemented, default: 1) \n");
  fprintf(stderr, "  -u  UCD file.                    (unimplemented, default: #no output#)  \n");
  fprintf(stderr, "  -c  partitioning contact         (unimplemented, default: DEFAULT)  \n");

}

extern int HECMW_init_for_partition(int argc, char **argv) {
  int counter;
  char control_file_name[HECMW_FILENAME_LEN + 1];

  strcpy(control_file_name, DEFAULT_CONTROL_FILE_NAME);

  if (argc > 1) {
    counter = 1;
    while (counter < argc) {
      if (!strcmp(argv[counter], "-f")) {
        counter++;
        if (counter >= argc) {
          print_usage();
          goto error;
        }
        if (strlen(argv[counter]) > HECMW_FILENAME_LEN) {
          HECMW_set_error(HECMW_PART_E_TOO_LONG_FNAME, "%s",
                          "control file for partitioner");
          goto error;
        }
        strcpy(control_file_name, argv[counter]);
        counter++;

      }else if (!strcmp(argv[counter], "-d")) {
        counter++;
        if (counter >= argc) {
          print_usage();
          goto error;
        }
        if (HECMW_part_set_subdomains(atoi(argv[counter]))) goto error;
        counter++;
      } else if (!strcmp(argv[counter], "-v")) {
        counter++;
        HECMW_setloglv(HECMW_LOG_DEBUG);

      } else {
        print_usage();
        goto error;
      }
    }
  }

  if (HECMW_part_set_ctrl_file_name(control_file_name)) goto error;

  return 0;

error:
  return -1;
}
