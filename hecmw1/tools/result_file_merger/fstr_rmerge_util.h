/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/**
 * @brief Utility for reading and processing results computed in parallel
 */

#include <stdio.h>
#include "hecmw_struct.h"
#include "hecmw_result.h"

/**
 * @struct fstr_res_info
 * @brief Analysis result file info
 */
typedef struct {
  int nnode_gid;
  int nelem_gid;
  int* node_gid;
  int* elem_gid;
  struct hecmwST_result_data* result;
} fstr_res_info;

/**
 * @struct fstr_gl_rec
 * @brief Global ID, local ID and belonging area record
 */
typedef struct {
  int global;
  int local;
  int area;
} fstr_gl_rec;

/**
 * @struct fstr_glt
 * @brief Table for global ID, local ID and belonging area records
 */
typedef struct {
  fstr_gl_rec* nrec;
  fstr_gl_rec* erec;
  int node_n;
  int elem_n;
} fstr_glt;

/**
 * @brief Set file pointer for log output
 */
void fstr_set_log_fp(FILE *log_fp);

/**
 * @brief Log output
 */
void fstr_out_log(const char* fmt, ...);

/**
 * @brief Read all distributed meshes
 */
struct hecmwST_local_mesh** fstr_get_all_local_mesh(char* name_ID, int nrank,
                                                    int* area_n, int* refine);

/**
 * @brief Delete mesh
 */
void fstr_free_mesh(struct hecmwST_local_mesh** mesh, int area_n);

/**
 * @brief Check the number of steps (check for the existence of files)
 */
int fstr_get_step_n(char* name_ID, int nrank);

/**
 * @brief Read all area data of step
 */
fstr_res_info** fstr_get_all_result(char* name_ID, int step, int area_n,
                                    int refine, int nrank);

/**
 * @brief Combine data in all areas of the step
 */
struct hecmwST_result_data* fstr_all_result(fstr_glt* glt, fstr_res_info** res,
                                            int refine);

/**
 * @biref Delete fstr_res_info
 */
void fstr_free_result(fstr_res_info** res, int area_n);

/**
 * @brief Create table for global ID, local ID and belonging area records fstr_glt
 */
fstr_glt* fstr_create_glt(struct hecmwST_local_mesh** mesh, int area_n);

/**
 * @brief Delete fstr_glt
 */
void fstr_free_glt(fstr_glt* glt);

/**
 * @brief Create a single region mesh
 */
struct hecmwST_local_mesh* fstr_create_glmesh(fstr_glt* glt);

/**
 * @brief Delete a single region mesh
 */
void fstr_free_glmesh(struct hecmwST_local_mesh* mp);
