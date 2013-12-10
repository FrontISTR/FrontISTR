#ifndef HECMW_VIS_CASE_TABLE_H_INCLUDED
#define HECMW_VIS_CASE_TABLE_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_vis_SF_geom.h"

int p2h(int i);
int get_data(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int elemID,
		Cell *cell, int tn_component);
double get_value_equ(double cont_equ[10],int cross_type,double x,double y,double z);
int make_tile(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, Cube_polygons *cube_polys, int disamb_flag);
int get_config_info(int index);

#endif /* HECMW_VIS_CASE_TABLE_H_INCLUDED */




