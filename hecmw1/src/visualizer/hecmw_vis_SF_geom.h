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


#ifndef HECMW_VIS_SF_GEOM_H_INCLUDED
#define HECMW_VIS_SF_GEOM_H_INCLUDED

/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <ctype.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_io.h"
*/
#include "hecmw_vis_psf_rendering.h"
/*
#include "hecmw_vis_comm_util.h"
*/

#define MASTER_PE 	0

#define EPSILON	 	0.00000001


#define TABLE_SIZE 	100


#define	VERTEX_PACK	50

#define POLYGON_PACK	100

#define VERTEX_KIND	27

#define HEX_N_NODE	8

#define HEX_N_FACE	6

#define HEX_NODE_INDEX	255	/* 2^8 */


#define HEX_FACE_INDEX	63	/* 2^6 */

#define PRISM_N_NODE	6

#define PRISM_N_FACE	5

#define PRISM_NODE_INDEX	63	/* 2^6 */

#define PRISM_FACE_INDEX	31	/* 2^5 */

#define TETRA_N_NODE	4

#define TETRA_N_FACE	4

#define TETRA_NODE_INDEX	15	/* 2^4 */

#define TETRA_FACE_INDEX	15	/* 2^4 */



#define MAX_N_NODE      20

#define HASH_TABLE_SIZE	10000


#define NUM_CONTROL_PSF     73
/*
#define MAX_LINE_LEN   256
#define buffer_size  300
#define UCD_NUM_CELL_TYPES	8
#define UCD_LABEL_LEN		1024
 */



typedef struct _psf_link_struct {

	int                      num_of_psf;

	struct surface_module    *sf;

	Parameter_rendering      *sr;
	struct _psf_link_struct  *next_psf;
	int                      stat_para[NUM_CONTROL_PSF];
	int                      visual_start_step;
	int                      visual_end_step;

	int                      visual_interval_step;

} PSF_link;


struct surface_module {
	int		surface_style;
	char      group_name[128];
	int       defined_style;
	char      data_comp_name[128];
	int       data_comp;
	char      data_subcomp_name[128];
	int       data_subcomp;
	double    iso_value;
	int       method;
	double    point[3];
	double    radius;
	double    length[3];
	double	coef[10];
	int		display_method;
	char      color_comp_name[128];
	int       color_comp;
	char      color_subcomp_name[128];
	int       color_subcomp;
	int       isoline_number;
	double    specified_color;
	int       normalize_flag;
	int       range_output;
	char      range_filename[128];

	int       deform_display_on;

	int       disp_comp;

	char      disp_comp_name[128];

	double    disp_scale;
	double    real_disp_scale;

	int       initial_style;

	int       deform_style;

	double    initial_line_color[3];

	double    deform_line_color[3];

	int       output_type;
};


typedef struct _surface_module_struct {
	int		surface_style;
	char      *group_name;
	int       defined_style;
	int       data_comp;
	int       data_subcomp;
	double    iso_value;
	int		color_comp;
	int       color_subcomp;
	int       cross_type;
	double    cont_equ[10];
	int       display_way;
	double    rgbrange[3];
	int       isonumber;
	double    specified_color;

	int       output_type;
	int opacity_choice;
	double    opacity;

	int       deform_display_on;

	int       disp_comp;

	char      disp_comp_name[128];

	double    disp_scale;
	double    real_disp_scale;
	int       initial_style;

	int       deform_style;

	double    initial_line_color[3];

	double    deform_line_color[3];
} Surface;


typedef struct _result_struct {
	int    n_vertex;
	int    n_patch;
	double *vertex;
	int    *patch;
	double *color;
	double *disp;
} Result;


typedef struct _fgeom_struct {
	double 	x;
	double 	y;
	double	z;
} Fgeom;


typedef struct _triangle_struct {
	int vertex[3];
} Triangle;

typedef struct _isoline_struct {
	Fgeom point[2];
	struct _isoline_struct *nextline;
} Isoline;

typedef struct _isohead_struct {
	int linenum;
	struct _isoline_struct *nextline;
} Isohead;

typedef struct _point_struct {
	int			ident;
	double			field;
	double			cdata;
	double          disp[3];
	Fgeom           	geom;
	int	                locator;
	int			bdflag;
	int			info;
	struct _point_struct    *nextpoint;
} Point;

typedef struct _polygon_struct {

	int		   	type;
	int     *plist;
	struct _polygon_struct  *nextpolygon;
	int			elem_id[2];
	int			bdflag;
} Polygon;
/* type: 0 then polygon is owned by alpha isosurface
         1 then polygon is owned by beta isosurface
         2 then polygon is owned by cross section   	*/
/* *plist = {a,b,c,d,e,....} : a is number of vertex,
		       		b,c,d,c,... is array of vertex(left turn) */
/* flag : 境界の要素かどうか */

typedef struct _polygon_obj_struct {
	int		verts_num;
	double		*field;
	double		*verts;
	unsigned long	*colors;
	int		*plist;
} Polygon_obj;


typedef struct _CS_polygon_struct {
	Polygon_obj	*alpha_obj;
	Polygon_obj	*cross_obj;
	double		area;
	double		volume;
	double		integral;
} CS_polygon_obj;


typedef struct _cube_polygons_struct {
	int		verts_num;
	int		verts[VERTEX_KIND];	/* 8 + 12 + 7 */
	int		**isosurf;		/* terminator is -1 */
} Cube_polygons;


typedef struct _rotation_info_struct {
	unsigned char	face[6];
	unsigned char	grid_point[8];
	unsigned char	edge[12];
	unsigned char	inside[7];
} Rotation_info;

typedef struct _cell_struct {
	double	axis[3*8];
	double	s_data[8];
	double	c_data[8];
	double     disp[3*8];
	double     v_data[3*8];
	int		elem_id[2];
} Cell;



typedef struct _tetra_struct {
	double	axis[3*4];
	double	s_data[4];
	double	c_data[4];
	double    disp[3*4];
	double    v_data[3*4];
	int		elem_id[2];
	int       local_vid[4];
} Tetra;

typedef struct _prism_struct {
	double	axis[3*6];
	double	s_data[6];
	double	c_data[6];
	double    disp[3*6];
	double    v_data[3*6];
	int		elem_id[2];
	int       local_vid[6];
} Prism;

/*
typedef struct _overlap_struct {
  int		index;
  int		elem_id[2];
  int		verts_num;
  int		*verts;
  struct _overlap_struct	*next_elem;
} Overlap;
 */



typedef struct _vertex_struct {
	int		index;
	double	x;
	double	y;
	double	z;
	double	color;
	double    disp[3];
	struct _vertex_struct	*next_verts;
} Vertex;
typedef struct _overlap_struct {
	int		index;
	int		elem_id[2];
	int		verts_num;
	int		*verts;
	struct _overlap_struct	*next_elem;
} Overlap;

/*
typedef struct _elem_no {
  int	peID;
  int	localID;
} Elem_no;
 */

typedef struct _hash_table_struct {
	int elemID;
	int faceID;
	struct _hash_table_struct *next_elem;
}  Hash_table;


typedef struct _boundary_patch_struct {
	int   type;
	/* 3== tri  4==quad */
	double vertex[4*3];
	double color[4];
	double disp[4*3];
	struct _boundary_patch_struct *next_patch;
} Boundary_patch;



typedef struct _point_tetra_struct {
	int			ident;
	double			cdata;
	double          disp[3];
	double          geom[3];

	struct _point_tetra_struct    *nextpoint;
} Tetra_point;


typedef struct _patch_tetra_struct {
	int  patch[3];
	struct _patch_tetra_struct *next_patch;
} Patch_tetra;

typedef struct _head_patch_tetra_struct {
	int  num_patch;
	Patch_tetra *patch_link;
} Head_patch_tetra;

typedef struct _hash_vertex_struct {
	int ident;
	double geom[3];
	struct _hash_vertex_struct *next_vertex;
} Hash_vertex;


typedef struct _connect_inf {
	int *index_connect;
	int *connect;
} Connect_inf;

#endif /* HECMW_VIS_SF_GEOM_H_INCLUDED */
