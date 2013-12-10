#ifndef HECMW_VIS_INTERSECTION_FIND_H_INCLUDED
#define HECMW_VIS_INTERSECTION_FIND_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_vis_SF_geom.h"


/*----------------------------------------------------------------------
#     Subroutines in this file on isosurface generation by Marching Cubes is based
	  on the revision of Dr. Yuriko Takeshima's codes when she was working part time in RIST
#---------------------------------------------------------------------- */


int get_hex_rotation_info(int rot_type, Rotation_info *rot_info);
int point_convert(int point_index, Rotation_info *rot_info);
int separation_test(struct hecmwST_local_mesh *mesh, Cell *cell, int amb,
		double fvalue, int CS_type, int disamb_flag);

void case0opp_tiler_hex(Cube_polygons *cube_polys);
void case1_tiler_hex(int index, Cube_polygons *cube_polys);
void case1opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case2_tiler_hex(int index, Cube_polygons *cube_polys);
void case2opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case3_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case3opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case4_tiler_hex(int index, Cube_polygons *cube_polys);
void case4opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case5_tiler_hex(int index, Cube_polygons *cube_polys);
void case5opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case6_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case6opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case7_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case7opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case8_tiler_hex(int index, Cube_polygons *cube_polys);
void case8opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case9_tiler_hex(int index, Cube_polygons *cube_polys);
void case9opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case10_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case10opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case11_tiler_hex(int index, Cube_polygons *cube_polys);
void case11opp_tiler_hex(int index, Cube_polygons *cube_polys);
void case12_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case12opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case13_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case13opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag);
void case14_tiler_hex(int index, Cube_polygons *cube_polys);
void case14opp_tiler_hex(int index, Cube_polygons *cube_polys);
int choice_disambiguation(int disamb_flag, double fvalue, int CS_type,
		double voxel0, double voxel1, double voxel2,
		double voxel3);

#endif /* HECMW_VIS_INTERSECTION_FIND_H_INCLUDED */
