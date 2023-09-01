/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/*
 * Merge FSTR result output files into a single file.
 * Control file and mesh file used in calculations are required.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "hecmw_util.h"
#include "hecmw_result_io_bin.h"
#include "hecmw_result_io_txt.h"

#include "fstr_rmerge_util.h"

#define BINARY_DEFAULT 0
#define NRANK_DEFAULT 0
#define STRID_DEFAULT 0
#define ENDID_DEFAULT -1
#define INTID_DEFAULT 1

void error_stop(void) {
  char* msg;
  HECMW_get_error(&msg);
  fputs(msg, stderr);
  exit(-1);
}

void help(void) {
  printf(" HECMW Result File Merger\n");
  printf("usage)  rmerge [options] out_fileheader\n");
  printf("[option]\n");
  printf(" -h             : help\n");
  printf(" -o [type]      : output file type (binary/text, default:%s)\n",
         (BINARY_DEFAULT == 1) ? "binary" : "text");
  printf(" -n [rank]      : number of ranks (default:%d)\n", NRANK_DEFAULT);
  printf(" -s [step]      : start step number (default:%d)\n", STRID_DEFAULT);
  printf(" -e [step]      : end step number (default:%d)\n", ENDID_DEFAULT);
  printf(" -i [step]      : interval step number (default:%d)\n", INTID_DEFAULT);
}

void set_fname(int argc, char** argv, char* out_fheader, int* binary,
               int *nrank, int *strid, int *endid, int *intid) {
  int i;
  char* fheader = NULL;

  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h")) {
      help();
      exit(-1);
    } else if (!strcmp(argv[i], "-o")) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        exit(-1);
      }
      i++;
      if (!strcmp(argv[i], "binary")) {
        *binary = 1;
      } else if (!strcmp(argv[i], "text")) {
        *binary = 0;
      } else {
        fprintf(stderr, "Error : text or binary is required after -o\n");
        exit(-1);
      }
    } else if (strcmp(argv[i], "-n") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        exit(-1);
      }
      i++;
      if (sscanf(argv[i], "%d", nrank) != 1) {
        fprintf(stderr,
                "Error : parameter %s cannot be converted to number of ranks\n",
                argv[i]);
        exit(-1);
      }
    } else if (strcmp(argv[i], "-s") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        exit(-1);
      }
      i++;
      if (sscanf(argv[i], "%d", strid) != 1) {
        fprintf(
            stderr,
            "Error : parameter %s cannot be converted to start step number\n",
            argv[i]);
        exit(-1);
      }
    } else if (strcmp(argv[i], "-e") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        exit(-1);
      }
      i++;
      if (sscanf(argv[i], "%d", endid) != 1) {
        fprintf(stderr,
                "Error : parameter %s cannot be converted to end step number\n",
                argv[i]);
        exit(-1);
      }
    } else if (strcmp(argv[i], "-i") == 0) {
      if (argc == i + 1) {
        fprintf(stderr, "Error : parameter required after %s\n", argv[i]);
        exit(-1);
      }
      i++;
      if (sscanf(argv[i], "%d", intid) != 1) {
        fprintf(stderr,
                "Error : parameter %s cannot be converted to interval step "
                "number\n",
                argv[i]);
        exit(-1);
      }
    } else {
      fheader = argv[i];
    }
  }

  if (!fheader) {
    help();
    exit(-1);
  } else {
    strcpy(out_fheader, fheader);
  }
}

int main(int argc, char** argv) {
  int binary = BINARY_DEFAULT;
  int nrank = NRANK_DEFAULT;
  int strid = STRID_DEFAULT;
  int endid = ENDID_DEFAULT;
  int intid = INTID_DEFAULT;
  int area_n, step_n, refine, fg_text;
  int step, rcode;
  char out_fheader[HECMW_HEADER_LEN + 1];
  char out_fname[HECMW_FILENAME_LEN + 1];
  struct hecmwST_local_mesh** mesh;
  struct hecmwST_local_mesh* glmesh;
  struct hecmwST_result_data* data;
  fstr_res_info** res;
  fstr_glt* glt;
  char header[HECMW_HEADER_LEN + 1];
  char comment[HECMW_MSG_LEN + 1];
  char* fileheader;
  char dirname[HECMW_HEADER_LEN + 1];
  char buff[HECMW_HEADER_LEN + 1];
  char *ptoken, *ntoken;

  fstr_set_log_fp(stderr);

  if (HECMW_init(&argc, &argv)) error_stop();

  set_fname(argc, argv, out_fheader, &binary, &nrank, &strid, &endid, &intid);
  fstr_out_log("out file name header is %s\n", out_fheader);

  mesh = fstr_get_all_local_mesh("fstrMSH", nrank, &area_n, &refine);
  if (!mesh) error_stop();

  fstr_out_log("table creating .. \n");
  glt = fstr_create_glt(mesh, area_n);
  if (!glt) {
    fprintf(stderr, "ERROR : Cannot create global_local table.\n");
    fstr_free_mesh(mesh, area_n);
    exit(-1);
  }

  glmesh = fstr_create_glmesh(glt);
  if (!glmesh) {
    fprintf(stderr, "ERROR : Cannot create global table.\n");
    fstr_free_mesh(mesh, area_n);
    fstr_free_glt(glt);
    exit(-1);
  }

  fstr_free_mesh(mesh, area_n);

  step_n = (endid > -1) ? endid : fstr_get_step_n("fstrRES", nrank);

  for (step = strid; step <= step_n; step++) {
    if ((step % intid) != 0 && step != step_n) continue;

    fstr_out_log("step:%d .. reading .. ", step);
    res = fstr_get_all_result("fstrRES", step, area_n, refine, nrank);
    if (!res) {
      fstr_free_result(res, area_n);
      continue;
    }
    fstr_out_log("end\n");

    fstr_out_log("step:%d .. combining .. ", step);
    data = fstr_all_result(glt, res, refine);
    if (!data) {
      fprintf(stderr, "ERROR : Cannot combine result structure.\n");
      fstr_free_glt(glt);
      fstr_free_glmesh(glmesh);
      fstr_free_result(res, area_n);
      exit(-1);
    }
    fstr_out_log("end\n");

    if (nrank == 0) {
      if ((fileheader = HECMW_ctrl_get_result_fileheader("fstrRES", step,
                                                         &fg_text)) == NULL)
        return 0;
    } else {
      if ((fileheader = HECMW_ctrl_get_result_fileheader_sub(
               "fstrRES", step, nrank, 0, &fg_text)) == NULL)
        return 0;
    }
    strcpy(buff, fileheader);
    strcpy(dirname, "");
    ptoken = strtok(buff, "/");
    ntoken = strtok(NULL, "/");
    while (ntoken) {
      strcat(dirname, ptoken);
      strcat(dirname, "/");
      ptoken = ntoken;
      ntoken = strtok(NULL, "/");
    }
    sprintf(out_fname, "%s%s.%d", dirname, out_fheader, step);
    fstr_out_log("output to %s .. ", out_fname);
    HECMW_result_get_header(header);
    HECMW_result_get_comment(comment);
    HECMW_result_init(glmesh, step, header, comment);
    if (binary) {
      rcode = HECMW_result_io_bin_write_ST_by_fname(out_fname, data, glt->node_n,
                                                 glt->elem_n, header, comment);
    } else {
      rcode = HECMW_result_io_txt_write_ST_by_fname(out_fname, data, glt->node_n,
                                                 glt->elem_n, header, comment);
    }
    if (rcode) {
      fprintf(stderr, "ERROR : Cannot open/write file %s\n", out_fname);
      fstr_free_glt(glt);
      fstr_free_glmesh(glmesh);
      fstr_free_result(res, area_n);
      HECMW_result_free(data);
      exit(-1);
    }
    fstr_out_log("end\n");

    fstr_free_result(res, area_n);
    HECMW_result_free(data);
  }

  fstr_free_glt(glt);
  fstr_free_glmesh(glmesh);

  HECMW_finalize();

  exit(0);
}
