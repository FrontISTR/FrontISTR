/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/**
 * @brief 分散で計算された結果を読込み処理するためのユーティリティ
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "hecmw_util.h"
#include "hecmw_io_mesh.h"
#include "hecmw_io_struct.h"
#include "hecmw_struct.h"
#include "hecmw_config.h"
#include "hecmw_dist.h"
#include "hecmw_dist_free.h"
#include "hecmw_common.h"

#include "hecmw_control.h"
#include "hecmw_result.h"
#include "hecmw_io_dist.h"
#include "hecmw_io_get_mesh.h"

/**
 * @struct fstr_res_info
 * @brief 解析結果ファイル情報
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
 * @brief グローバル・ローカルID対応
 */
typedef struct {
  int global;
  int local;
  int area;
} fstr_gl_rec;

/**
 * @struct fstr_glt
 * @brief グローバル・ローカル対応表
 */
typedef struct {
  fstr_gl_rec* nrec;
  fstr_gl_rec* erec;
  int node_n;
  int elem_n;
} fstr_glt;

/**
 * @brief 全分散メッシュの読込み
 */
struct hecmwST_local_mesh** fstr_get_all_local_mesh(char* name_ID, int* area_n,
                                                    int* refine);

/**
 * @brief メッシュの削除
 */
void fstr_free_mesh(struct hecmwST_local_mesh** mesh, int area_n);

/**
 * @brief ステップ数を調べる(ファイルの存在を調べる)
 */
int fstr_get_step_n(char* name_ID);

/**
 * @brief ステップの全領域データの読み込み
 */
fstr_res_info** fstr_get_all_result(char* name_ID, int step, int area_n,
                                    int refine);

/**
 * @brief ステップの全領域のデータの結合
 */
struct hecmwST_result_data* fstr_all_result(fstr_glt* glt, fstr_res_info** res,
                                            int refine);

/**
 * @biref fstr_res_info の削除
 */
void fstr_free_result(fstr_res_info** res, int area_n);

/**
 * @brief テーブル fstr_glt の作成
 */
fstr_glt* fstr_create_glt(struct hecmwST_local_mesh** mesh, int area_n);

/**
 * @brief fstr_glt の削除
 */
void fstr_free_glt(fstr_glt* glt);

/**
 * @brief 単一領域メッシュの作成
 */
struct hecmwST_local_mesh* fstr_create_glmesh(fstr_glt* glt);

/**
 * @brief 単一領域メッシュの削除
 */
void fstr_free_glmesh(struct hecmwST_local_mesh* mp);
