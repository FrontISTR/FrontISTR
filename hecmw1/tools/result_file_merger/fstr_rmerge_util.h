/*
	fstr_res_util Ver.1.3
	2010.09.27 by K.Suemitsu (AdvanceSoft)
	2004.10.26 by N.Imai (RIST)
	--------------------------------------------------------
	分散で計算された結果を読込み処理するためのユーティリティ
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


/* 解析結果ファイル情報 */
typedef struct {
	int nnode_gid;
	int nelem_gid;
	int* node_gid;
	int* elem_gid;
	struct hecmwST_result_data* result;
} fstr_res_info;

/* グローバル⇔ローカルID 対応 */
typedef struct {
	int global;
	int local;
	int area;
} fstr_gl_rec;

/* グローバル⇔ローカル対応表 */
/* 注）グルーバルIDが重複する可能性あり */
typedef struct {
	fstr_gl_rec* nrec;
	fstr_gl_rec* erec;
	int node_n;
	int elem_n;
} fstr_gl_t;


/* 全分散メッシュの読込み( area_n : 分散領域の数) */
struct hecmwST_local_mesh** fstr_get_all_local_mesh( char* name_ID, int* area_n );

/* メッシュの削除 */
void fstr_free_mesh( struct hecmwST_local_mesh** mesh, int area_n );

/* ステップ数を調べる（ファイルの存在を調べる） */
int fstr_get_step_n( char* name_ID );

/* ステップ step の全領域のデータを得る */
fstr_res_info** fstr_get_all_result( char* name_ID, int step, int area_n );

/* ステップ step の全領域のデータを結合する */
struct hecmwST_result_data* fstr_all_result( fstr_gl_t* glt,
	fstr_res_info** res, int* n_node, int* n_elem );

/* result の削除 */
void fstr_free_result( fstr_res_info** res, int area_n );

/* テーブル fstr_gl_t の作成 */
fstr_gl_t* fstr_create_glt( struct hecmwST_local_mesh** mesh, int area_n );

/* リファイン時要素テーブルの作成 */
fstr_gl_t* fstr_refine_glt( fstr_gl_t* glt, fstr_res_info** res, int area_n );

/* fstr_gl_t の削除 */
void fstr_free_gl_t( fstr_gl_t* glt );

/* グローバルIDメッシュの作成 */
struct hecmwST_local_mesh* fstr_create_glmesh( fstr_gl_t* glt );

/* グローバルIDメッシュの削除 */
void fstr_free_glmesh( struct hecmwST_local_mesh* mp );

