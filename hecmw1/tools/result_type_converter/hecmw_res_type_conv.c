/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_result_io_txt.h"
#include "hecmw_util.h"
#include "hecmw_io.h"

int strid                         = 1;
int endid                         = 1;
int intid                         = 1;
char out_file[HECMW_NAME_LEN + 1] = "";

void help(void) {
  printf(" HECMW Result File Type Converter\n");
  printf("usage)  rconv [options]\n");
  printf("[option]\n");
  printf(" -h             : help\n");
  printf(" -o [file]      : output file name without rank and step number\n");
  printf(" -s [step]      : start step number (default:%d)\n", strid);
  printf(" -e [step]      : end step number (default:%d)\n", endid);
  printf(" -i [step]      : interval step number (default:%d)\n", intid);
}

int set_params(int argc, char **argv) {
  int i;

  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      help();
      return -1;
    } else if (strcmp(argv[i], "-o") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        return -1;
      }
      i++;
      strcpy(out_file, argv[i]);
    } else if (strcmp(argv[i], "-s") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        return -1;
      }
      i++;
      if (sscanf(argv[i], "%d", &strid) != 1) {
        fprintf(
            stderr,
            "Error : parameter %s cannot be converted to start step number\n",
            argv[i]);
        return -1;
      }
    } else if (strcmp(argv[i], "-e") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        return -1;
      }
      i++;
      if (sscanf(argv[i], "%d", &endid) != 1) {
        fprintf(stderr,
                "Error : parameter %s cannot be converted to end step number\n",
                argv[i]);
        return -1;
      }
    } else if (strcmp(argv[i], "-i") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        return -1;
      }
      i++;
      if (sscanf(argv[i], "%d", &intid) != 1) {
        fprintf(stderr,
                "Error : parameter %s cannot be converted to interval step "
                "number\n",
                argv[i]);
        return -1;
      }
    } else {
      fprintf(stderr, "Error : invalid parameter %s\n", argv[i]);
      help();
      return -1;
    }
  }

  return 0;
}

int main(int argc, char **argv) {
  struct hecmwST_result_data *data;
  char *fileheader, resultfile[HECMW_FILENAME_LEN + 1];
  char header[HECMW_HEADER_LEN + 1];
  char comment[HECMW_MSG_LEN + 1];
  char dirname[HECMW_HEADER_LEN + 1];
  char buff[HECMW_HEADER_LEN + 1];
  char *ptoken, *ntoken;
  int n_node, n_elem, rcode, fg_text;
  int i, mynode;

  if (HECMW_init(&argc, &argv)) {
    HECMW_abort(HECMW_comm_get_comm());
  }

  if (set_params(argc, argv)) {
    HECMW_abort(HECMW_comm_get_comm());
  }

  mynode = HECMW_comm_get_rank();

  for (i = strid; i <= endid; i++) {
    if ((i % intid) != 0 && i != endid) continue;

    fileheader =
        HECMW_ctrl_get_result_fileheader("fstrRES", i, &fg_text);
    sprintf(resultfile, "%s.%d.%d", fileheader, mynode, i);
    fprintf(stdout, "Input file : %s\n", resultfile);
    data = HECMW_result_read_by_fname(resultfile);
    if (!data) {
      HECMW_abort(HECMW_comm_get_comm());
    }

    if (out_file[0]) {
      strcpy(buff, resultfile);
      strcpy(dirname, "");
      ptoken = strtok(buff, "/");
      ntoken = strtok(NULL, "/");
      while (ntoken) {
        strcat(dirname, ptoken);
        strcat(dirname, "/");
        ptoken = ntoken;
        ntoken = strtok(NULL, "/");
      }
      sprintf(resultfile, "%s%s.%d.%d", dirname, out_file, mynode, i);
    }
    fprintf(stdout, "Output file : %s\n", resultfile);

    n_node = HECMW_result_get_nnode();
    n_elem = HECMW_result_get_nelem();
    HECMW_result_get_header(header);
    HECMW_result_get_comment(comment);
    rcode = HECMW_result_io_txt_write_ST_by_fname(resultfile, data, n_node, n_elem,
                                               header, comment);
    if (rcode) {
      HECMW_abort(HECMW_comm_get_comm());
    }

    HECMW_result_free(data);
    HECMW_result_free_nodeID();
    HECMW_result_free_elemID();
  }

  HECMW_finalize();

  return 0;
}
