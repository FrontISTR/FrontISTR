#ifndef HECMW_VIS_GENERATE_HISTOGRAM_SF_H_INCLUDED
#define HECMW_VIS_GENERATE_HISTOGRAM_SF_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"

void generate_histogram_graph_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue, Result *result, int mynode,
		int pesize, HECMW_Comm VIS_COMM, int color_system_type);
void generate_interval_point_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue,  Result *result,
		int mynode, int pesize, HECMW_Comm VIS_COMM, double *interval_point);
void output_histogram_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue, Result *result,int mynode, int pesize, HECMW_Comm VIS_COMM);
void find_minmax_sf(struct hecmwST_local_mesh *mesh, int mynode, double range[6]);
void find_patch_minmax_sf(Result *result, struct surface_module *sf, double range[6]);
void find_new_patch_minmax_sf(Result *result, struct surface_module *sf, double range[6]);
void get_deform_scale(struct surface_module *sf, int ii, double range_x, double range_y, double range_z,
		struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int pesize, HECMW_Comm VIS_COMM);

#endif /* HECMW_VIS_GENERATE_HISTOGRAM_SF_H_INCLUDED */
