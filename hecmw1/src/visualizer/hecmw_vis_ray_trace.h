/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Visualization                                     *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef HECMW_VIS_RAY_TRACE_H_INCLUDED
#define HECMW_VIS_RAY_TRACE_H_INCLUDED

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include    <memory.h>
#include <ctype.h>
#include "hecmw_vis_resampling.h"
#include "hecmw_vis_bmp.h"
#include "hecmw_vis_comm_util.h"
*/
#define MASTER_PE 	0
/*
#define CONVERSE_ORDER

#include    <sys/types.h>
#include    <sys/timeb.h>
#include    <time.h>
#include    "glos.h"
#include    <GL/gl.h>
#include    <GL/glu.h>
#include    <GL/glaux.h>
 */


#define EPSILON	 	0.00000001
#define PI  3.1415926
/*#define TABLE_SIZE 	100
#define	VERTEX_PACK	50
#define POLYGON_PACK	100
#define VERTEX_KIND	27

#define HEX_N_NODE	8
#define HEX_N_FACE	6
#define HEX_NODE_INDEX	255
#define HEX_FACE_INDEX	63
 */
#define SQR(x) (x)*(x)
#define HASH_TABLE_SIZE	10000
/*#define  ResSize 1
#define  PixelSize 1
 */
#define  DIGN_PE  2

#define BAR_WIDTH 10
#define NUM_CONTROL_PVR 51
#define MAX_N_NODE      20




typedef struct _vr_parameter_struct {
	int max_level;
	int xr;
	int yr;
	int projection_style;
	int num_of_lights;
	double *light_point;

	double view_point_d[3];
	double screen_point[3];

	double up[3];

	double k_ads[3];
	int surface_on;
	double surface_opacity;

	int color_mapping_style;

	int interval_mapping_num;


	double *interval_point; /* 2:mincolor, maxcolor 3: interval_mapping_num*2 (value, mark_value) */


	int transfer_function_style;
	/*  1: constant    input: value
    2: first-order derivatives  input: none
	3: feature points  input:  num_of_featurepoints, point[num]
	4: feature intervals  input: num_of_intervals  point[num*2]
	5: distance inverse
	6: distance proportional
	7: look-up table      input: name of the look-up table file
	 */
	double opa_value;
	int num_of_features;
	double *fea_point;
	char name_lookup[128];

	int  rotate_style;
	int color_mapping_bar_on;
	int scale_marking_on;
	int num_of_frames;
	char name_voxelfile[128];
	double background_color[3];
	double font_color[3];
	int   color_system_type;
	double font_size;
	int   color_bar_style;
	int   fixed_range_on;
	double range_value[2];
	int    num_of_scale;
	int    mark_0_on;
	int    histogram_on;


	int    remove_0_display_on;
	int    specified_level[3];
	char      color_comp_name[100];
	int       color_comp;
	char      color_subcomp_name[5];
	int       color_subcomp;
	int    nv_xyz[3];
	double display_range[6];
	int    time_mark_on;
	int    fixed_scale_mark;
} Parameter_vr;



typedef struct _pvr_link_struct {

	int                      num_of_pvr;

	Parameter_vr             *vr;
	struct _pvr_link_struct  *next_pvr;
	int                      stat_para[NUM_CONTROL_PVR];
	int                      visual_start_step;
	int                      visual_end_step;
	int                      visual_interval_step;
} PVR_link;


typedef struct surface_info_struct {
	int num;
	double *surf_data;
} Surface_info;

typedef struct _data_structured_vr_struct {
	double dxyz[3];
	int leveltot;
	int varnumtot;
	char **varname;
	int nxyz[3];
	double xyz0[3];
	int voxtotadd;
	int voxtotall;
	/*	int *rlevel;
	int *parent;
	 */
	int r_nxyz[3];
	double r_dxyz[3];
	int *empty_flag;
	double *var;
	double *grad_var;
	Surface_info *surface;
} VR_data;

typedef struct _tree_pointer_struct {
	int     cell_id[8];
	int     surf_id;
	int     level;
	double  bound_box[6];
	int     local_child_no;
	int     local_face_in;
	int     local_face_out;
	struct _tree_pointer_struct *child;
	struct _tree_pointer_struct *parent;
} Tree_pointer;

typedef Tree_pointer *Tree_pointer_ptr;
/*
typedef struct _ray_volume_struct {
	Elem_no elem_id;
	int     face_id;
    double  p[3];
	struct _ray_volume_struct *next_elem;
	} Ray_volume;

typedef struct _head_ray_volume_struct {
	int elem_num;
	Ray_volume *next_elem;
} Head_ray_volume;
 */

int find_first_inter(double point_o[3], double view_point_d[3],  int r_level[3], double orig_xyz[3], double dxyz[3],  double r_dxyz[3], double ray_direction[3], double first_p[3], int ijk[3]);
void ray_trace(int remove_0_display_on, int color_mapping_style, double *interval_point,int transfer_function_style,
		double opa_value, int num_of_features, double *fea_point,
		double view_point_d[3],  int interval_mapping_num, int color_system_type, int num_of_lights,
		double *light_point, double k_ads[3], double orig_xyz[3], double dxyz[3],double r_dxyz[3], int r_level[3], int *empty_flag, double *var, double *grad_var, double first_p[3], int first_ijk[3], double ray_direction[3], double mincolor,
		double maxcolor, double accum_rgba[4], double grad_minmax[2], double feap_minmax[2],
		double feai_minmax[2], double dis_minmax[2], double *opa_table, double tav_length, int time_step, int test_i, int test_j);

#endif /* HECMW_VIS_RAY_TRACE_H_INCLUDED */
