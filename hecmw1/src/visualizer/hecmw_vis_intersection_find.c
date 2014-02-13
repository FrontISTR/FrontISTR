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

#include "hecmw_vis_intersection_find.h"

#include "hecmw_vis_calc_attr.h"


/*  provide information for rotation of cube  */
int get_hex_rotation_info(int rot_type, Rotation_info *rot_info)
{
	static Rotation_info rotation[24] = {
			{
					{0,1,2,3,4,5}, {0,1,2,3,4,5,6,7},  /* rotation No.0 */
					{0,1,2,3,4,5,6,7,8,9,10,11}, {0,1,2,3,4,5,6}
			},{
					{5,2,4,0,1,3}, {0,4,5,1,3,7,6,2},  /* rotation No.1 */
					{8,4,9,0,10,6,11,2,3,7,1,5}, {0,2,3,1,6,4,5}
			},{
					{3,4,1,5,2,0}, {0,3,7,4,1,2,6,5},  /* rotation No.2 */
					{3,10,7,8,1,11,5,9,0,2,4,6}, {0,3,1,2,5,6,4}
			},{
					{1,2,3,0,4,5}, {1,5,6,2,0,4,7,3},  /* rotation No.3 */
					{9,5,11,1,8,7,10,3,0,4,2,6}, {1,3,0,2,5,4,6}
			},{
					{5,3,4,1,2,0}, {1,0,4,5,2,3,7,6},  /* rotation No.4 */
					{0,8,4,9,2,10,6,11,1,3,5,7}, {1,0,2,3,6,5,4}
			},{
					{0,4,2,5,3,1}, {1,2,3,0,5,6,7,4},  /* rotation No.5 */
					{1,2,3,0,5,6,7,4,9,11,8,10}, {1,2,3,0,4,6,5}
			},{
					{1,5,3,4,2,0}, {2,1,5,6,3,0,4,7},  /* rotation No.6 */
					{1,9,5,11,3,8,7,10,2,0,6,4}, {2,1,3,0,5,6,4}
			},{
					{0,3,2,1,5,4}, {2,3,0,1,6,7,4,5},  /* rotation No.7 */
					{2,3,0,1,6,7,4,5,11,10,9,8}, {2,3,0,1,4,5,6}
			},{
					{4,2,5,0,3,1}, {2,6,7,3,1,5,4,0},  /* rotation No.8 */
					{11,6,10,2,9,4,8,0,1,5,3,7}, {2,0,1,3,6,4,5}
			},{
					{4,1,5,3,2,0}, {3,2,6,7,0,1,5,4},  /* rotation No.9 */
					{2,11,6,10,0,9,4,8,3,1,7,5}, {3,2,0,1,6,5,4}
			},{
					{0,5,2,4,1,3}, {3,0,1,2,7,4,5,6},  /* rotation No.10 */
					{3,0,1,2,7,4,5,6,10,8,11,9}, {3,0,1,2,4,6,5}
			},{
					{3,2,1,0,5,4}, {3,7,4,0,2,6,5,1},  /* rotation No.11 */
					{10,7,8,3,11,5,9,1,2,6,0,4}, {3,1,2,0,5,4,6}
			},{
					{5,1,4,3,0,2}, {4,5,1,0,7,6,2,3},  /* rotation No.12 */
					{4,9,0,8,6,11,2,10,7,5,3,1}, {2,3,1,0,6,5,4}
			},{
					{2,4,0,5,1,3}, {4,7,6,5,0,3,2,1},  /* rotation No.13 */
					{7,6,5,4,3,2,1,0,8,10,9,11}, {2,1,0,3,4,6,5}
			},{
					{3,0,1,2,4,5}, {4,0,3,7,5,1,2,6},  /* rotation No.14 */
					{8,3,10,7,9,1,11,5,4,0,6,2}, {2,0,3,1,5,4,6}
			},{
					{2,3,0,1,4,5}, {5,4,7,6,1,0,3,2},  /* rotation No.15 */
					{4,7,6,5,0,3,2,1,9,8,11,10}, {3,2,1,0,4,5,6}
			},{
					{5,0,4,2,3,1}, {5,1,0,4,6,2,3,7},  /* rotation No.16 */
					{9,0,8,4,11,2,10,6,5,1,7,3}, {3,1,0,2,6,4,5}
			},{
					{1,4,3,5,0,2}, {5,6,2,1,4,7,3,0},  /* rotation No.17 */
					{5,11,1,9,7,10,3,8,4,6,0,2}, {3,0,2,1,5,6,4}
			},{
					{1,0,3,2,5,4}, {6,2,1,5,7,3,0,4},  /* rotation No.18 */
					{11,1,9,5,10,3,8,7,6,2,4,0}, {0,2,1,3,5,4,6}
			},{
					{4,3,5,1,0,2}, {6,7,3,2,5,4,0,1},  /* rotation No.19 */
					{6,10,2,11,4,8,0,9,5,7,1,3}, {0,1,3,2,6,5,4}
			},{
					{2,5,0,4,3,1}, {6,5,4,7,2,1,0,3},  /* rotation No.20 */
					{5,4,7,6,1,0,3,2,11,9,10,8}, {0,3,2,1,4,6,5}
			},{
					{4,0,5,2,1,3}, {7,3,2,6,4,0,1,5},  /* rotation No.21 */
					{10,2,11,6,8,0,9,4,7,3,5,1}, {1,3,2,0,6,4,5}
			},{
					{3,5,1,4,0,2}, {7,4,0,3,6,5,1,2},  /* rotation No.22 */
					{7,8,3,10,5,9,1,11,6,4,2,0}, {1,2,0,3,5,6,4}
			},{
					{2,1,0,3,5,4}, {7,6,5,4,3,2,1,0},  /* rotation No.23 */
					{6,5,4,7,2,1,0,3,10,11,8,9}, {1,0,3,2,4,5,6}
			}
	};

	*rot_info = rotation[rot_type];

	return 1;
}

/*  decide the point location after the rotation  */
int point_convert(int point_index, Rotation_info *rot_info)
{
	int 		converted_index;

	if (point_index < 100) {
		converted_index = (int) rot_info->edge[point_index];
	} else if (point_index < 200) {
		converted_index = 100 + (int) rot_info->inside[point_index-100];
	} else {
		converted_index = 200 + (int) rot_info->grid_point[point_index-200];
	}

	return converted_index;
}

/*  decide the separation pattern on ambiguous face  */
int separation_test(struct hecmwST_local_mesh *mesh, Cell *cell, int amb,
		double fvalue, int CS_type, int disamb_flag)
{
	int 		sep;

	switch(amb) {
	case  0:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[0], cell->s_data[1],
				cell->s_data[2], cell->s_data[3]);
		break;
	case  1:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[1], cell->s_data[2],
				cell->s_data[3], cell->s_data[0]);
		break;
	case  2:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[1], cell->s_data[5],
				cell->s_data[6], cell->s_data[2]);
		break;
	case  3:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[5], cell->s_data[6],
				cell->s_data[2], cell->s_data[1]);
		break;
	case  4:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[5], cell->s_data[4],
				cell->s_data[7], cell->s_data[6]);
		break;
	case  5:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[4], cell->s_data[7],
				cell->s_data[6], cell->s_data[5]);
		break;
	case  6:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[4], cell->s_data[0],
				cell->s_data[3], cell->s_data[7]);
		break;
	case  7:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[0], cell->s_data[3],
				cell->s_data[7], cell->s_data[4]);
		break;
	case  8:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[3], cell->s_data[2],
				cell->s_data[6], cell->s_data[7]);
		break;
	case  9:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[2], cell->s_data[6],
				cell->s_data[7], cell->s_data[3]);
		break;
	case 10:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[4], cell->s_data[5],
				cell->s_data[1], cell->s_data[0]);
		break;
	case 11:
		sep = choice_disambiguation(disamb_flag, fvalue, CS_type,
				cell->s_data[5], cell->s_data[1],
				cell->s_data[0], cell->s_data[4]);
		break;
	}

	return sep;
}


/*  make the boundary polygons of (alpha || beta) interval volume in cube  */
void case0opp_tiler_hex(Cube_polygons *cube_polys)
{
	int 		i;
	static int 	verts[8] = {200, 201, 202, 203, 204, 205, 206, 207};
	static int 	isosurf[1] = {-1};
	/*  static int 	bounds[6][6] = {{4, 0,1,2,3, -1},{4, 1,5,6,2, -1},
				{4, 5,4,7,6, -1},{4, 4,0,3,7, -1},
				{4, 3,2,6,7, -1},{4, 4,5,1,0, -1}};
	 */

	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++)
		cube_polys->verts[i] = verts[i];

	*cube_polys->isosurf = isosurf;
	/*  for (i = 0; i < 6; i++) {
	 *cube_polys->bounds[i] = bounds[i];
  }
	 */
}

void case1_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info	rot_info;
	static int 	verts[4] = {0, 3, 8, 200};
	static int 	isosurf[4] = {0,2,1, -1};
	/*  static int 	bounds[6][5] = {{3, 3,0,1, -1},{-1, -1,-1,-1,-1},
				{-1, -1,-1,-1,-1},{3, 3,1,2, -1},
				{-1, -1,-1,-1,-1},{3, 3,2,0, -1}};
	 */
	switch(index) {
	case   1: get_hex_rotation_info( 0, &rot_info);  break;
	case   2: get_hex_rotation_info( 3, &rot_info);  break;
	case   4: get_hex_rotation_info( 6, &rot_info);  break;
	case   8: get_hex_rotation_info( 9, &rot_info);  break;
	case  16: get_hex_rotation_info(12, &rot_info);  break;
	case  32: get_hex_rotation_info(15, &rot_info);  break;
	case  64: get_hex_rotation_info(18, &rot_info);  break;
	case 128: get_hex_rotation_info(21, &rot_info);  break;
	}

	cube_polys->verts_num = 4;
	for (i = 0; i < 4; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	/*  for (i = 0; i < 6; i++) {
	 *cube_polys->bounds[rot_info.face[i]] = bounds[i];
  }
	 */
}

void case1opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 3, 8, 201, 202, 203, 204, 205, 206, 207};
	static int 	isosurf[4] = {0,1,2, -1};
	/*  static int 	bounds[6][7] = {{5, 0,3,4,5,1, -1},{4, 3,7,8,4, -1,-1},
				{4, 7,6,9,8, -1,-1},{5, 1,5,9,6,2, -1},
				{4, 5,4,8,9, -1,-1},{5, 2,6,7,3,0, -1}};
	 */
	switch(index) {
	case 254: get_hex_rotation_info( 0, &rot_info);  break;
	case 253: get_hex_rotation_info( 3, &rot_info);  break;
	case 251: get_hex_rotation_info( 6, &rot_info);  break;
	case 247: get_hex_rotation_info( 9, &rot_info);  break;
	case 239: get_hex_rotation_info(12, &rot_info);  break;
	case 223: get_hex_rotation_info(15, &rot_info);  break;
	case 191: get_hex_rotation_info(18, &rot_info);  break;
	case 127: get_hex_rotation_info(21, &rot_info);  break;
	}

	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	/*  for (i = 0; i < 6; i++) {
	 *cube_polys->bounds[rot_info.face[i]] = bounds[i];
  }
	 */
}

void case2_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[6] = {1, 3, 8, 9, 200, 201};
	static int 	isosurf[7] = {1,0,2,  3,2,0, -1};
	/*  static int 	bounds[6][6] = {{4, 4,5,0,1, -1},{3, 5,3,0, -1,-1},
				{-1, -1,-1,-1,-1,-1},{3, 4,1,2, -1,-1},
				{-1, -1,-1,-1,-1,-1},{4, 4,2,3,5, -1}};
	 */
	switch(index) {
	case   3: get_hex_rotation_info( 0, &rot_info);  break;
	case   6: get_hex_rotation_info( 5, &rot_info);  break;
	case   9: get_hex_rotation_info( 2, &rot_info);  break;
	case  12: get_hex_rotation_info( 7, &rot_info);  break;
	case  17: get_hex_rotation_info( 1, &rot_info);  break;
	case  34: get_hex_rotation_info( 3, &rot_info);  break;
	case  48: get_hex_rotation_info(12, &rot_info);  break;
	case  68: get_hex_rotation_info( 8, &rot_info);  break;
	case  96: get_hex_rotation_info(17, &rot_info);  break;
	case 136: get_hex_rotation_info(11, &rot_info);  break;
	case 144: get_hex_rotation_info(13, &rot_info);  break;
	case 192: get_hex_rotation_info(19, &rot_info);  break;
	}

	cube_polys->verts_num = 6;
	for (i = 0; i < 6; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	/*  for (i = 0; i < 6; i++) {
	 *cube_polys->bounds[rot_info.face[i]] = bounds[i];
  }
	 */
}

void case2opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {1, 3, 8, 9, 202, 203, 204, 205, 206, 207};
	static int 	isosurf[7] = {0,1,2,  2,3,0, -1};
	/*  static int 	bounds[6][7] = {{4, 5,1,0,4, -1,-1},{5, 4,0,3,7,8, -1},
				{4, 7,6,9,8, -1,-1},{5, 1,5,9,6,2, -1},
				{4, 5,4,8,9, -1,-1},{4, 2,6,7,3, -1,-1}};
	 */
	switch(index) {
	case 252: get_hex_rotation_info( 0, &rot_info);  break;
	case 249: get_hex_rotation_info( 5, &rot_info);  break;
	case 246: get_hex_rotation_info( 2, &rot_info);  break;
	case 243: get_hex_rotation_info( 7, &rot_info);  break;
	case 238: get_hex_rotation_info( 1, &rot_info);  break;
	case 221: get_hex_rotation_info( 3, &rot_info);  break;
	case 207: get_hex_rotation_info(12, &rot_info);  break;
	case 187: get_hex_rotation_info( 8, &rot_info);  break;
	case 159: get_hex_rotation_info(17, &rot_info);  break;
	case 119: get_hex_rotation_info(11, &rot_info);  break;
	case 111: get_hex_rotation_info(13, &rot_info);  break;
	case  63: get_hex_rotation_info(19, &rot_info);  break;
	}

	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	/*  for (i = 0; i < 6; i++) {
	 *cube_polys->bounds[rot_info.face[i]] = bounds[i];
  }
	 */
}

void case3_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[8] = {0, 1, 2, 3, 8, 11, 200, 202};
	static int 	isosurf_A[7] = {0,4,3,  2,5,1, -1};
	static int 	isosurf_B[13] = {3,2,4,  2,5,4,  4,5,1,  4,1,0, -1};
	/*  static int 	bounds_0[2][9] = {{3, 6,0,3,  3, 7,2,1, -1},
				  {6, 6,0,1,7,2,3, -1,-1}};
  static int 	bounds_1[5] = {3, 7,1,5, -1};
  static int 	bounds_2[1] = {-1};
  static int 	bounds_3[5] = {3, 6,3,4, -1};
  static int 	bounds_4[5] = {3, 7,5,2, -1};
  static int 	bounds_5[5] = {3, 6,4,0, -1};
	 */
	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;	/* isosurf_A */
		switch(index) {
		case   5: get_hex_rotation_info( 0, &rot_info); break;
		case  10: get_hex_rotation_info(10, &rot_info); break;
		case  18: get_hex_rotation_info( 4, &rot_info); break;
		case  24: get_hex_rotation_info(11, &rot_info); break;
		case  33: get_hex_rotation_info( 1, &rot_info); break;
		case  36: get_hex_rotation_info( 6, &rot_info); break;
		case  66: get_hex_rotation_info(18, &rot_info); break;
		case  72: get_hex_rotation_info( 9, &rot_info); break;
		case  80: get_hex_rotation_info(13, &rot_info); break;
		case 129: get_hex_rotation_info(22, &rot_info); break;
		case 132: get_hex_rotation_info(21, &rot_info); break;
		case 160: get_hex_rotation_info(15, &rot_info); break;
		}
	} else {
		switch(index) {
		case   5: amb_index = separation_test(mesh, cell, 0,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		case  10: amb_index = separation_test(mesh, cell, 1,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(10, &rot_info);
		break;
		case  18: amb_index = separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 4, &rot_info);
		break;
		case  24: amb_index = separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(11, &rot_info);
		break;
		case  33: amb_index = separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 1, &rot_info);
		break;
		case  36: amb_index = separation_test(mesh, cell, 3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 6, &rot_info);
		break;
		case  66: amb_index = separation_test(mesh, cell, 2,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(18, &rot_info);
		break;
		case  72: amb_index = separation_test(mesh, cell, 8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 9, &rot_info);
		break;
		case  80: amb_index = separation_test(mesh, cell, 5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(13, &rot_info);
		break;
		case 129: amb_index = separation_test(mesh, cell, 7,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(22, &rot_info);
		break;
		case 132: amb_index = separation_test(mesh, cell, 9,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(21, &rot_info);
		break;
		case 160: amb_index = separation_test(mesh, cell, 4,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(15, &rot_info);
		break;
		}
	}

	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	if (!amb_index) {
		*cube_polys->isosurf = isosurf_A;
		/*    *cube_polys->bounds[rot_info.face[0]] = bounds_0[0];*/
	} else {
		*cube_polys->isosurf = isosurf_B;
		/*    *cube_polys->bounds[rot_info.face[0]] = bounds_0[1];*/
	}

	/*  *cube_polys->bounds[rot_info.face[1]] = bounds_1;
	 *cube_polys->bounds[rot_info.face[2]] = bounds_2;
	 *cube_polys->bounds[rot_info.face[3]] = bounds_3;
	 *cube_polys->bounds[rot_info.face[4]] = bounds_4;
	 *cube_polys->bounds[rot_info.face[5]] = bounds_5;
	 */
}

void case3opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[12] = {0, 1, 2, 3, 8, 11, 201, 203, 204, 205, 206, 207};
	static int 	isosurf_A[7] = {0,3,4,  2,1,5, -1};
	static int 	isosurf_B[13] = {3,4,2,  2,4,5,  4,1,5,  4,0,1, -1};
	/*  static int 	bounds_0[2][9] = {{6, 7,3,0,6,1,2, -1,-1},
				  {3, 6,1,0,  3, 7,3,2, -1}};
  static int 	bounds_1[7] = {5, 6,9,10,5,1, -1};
  static int 	bounds_2[6] = {4, 9,8,11,10, -1};
  static int 	bounds_3[7] = {5, 7,11,8,4,3, -1};
  static int 	bounds_4[7] = {5, 10,11,7,2,5, -1};
  static int 	bounds_5[7] = {5, 8,9,6,0,4, -1};
	 */
	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;	/* isosurf_A */
		switch(index) {
		case 250: get_hex_rotation_info( 0, &rot_info); break;
		case 245: get_hex_rotation_info(10, &rot_info); break;
		case 237: get_hex_rotation_info( 4, &rot_info); break;
		case 231: get_hex_rotation_info(11, &rot_info); break;
		case 222: get_hex_rotation_info( 1, &rot_info); break;
		case 219: get_hex_rotation_info( 6, &rot_info); break;
		case 189: get_hex_rotation_info(18, &rot_info); break;
		case 183: get_hex_rotation_info( 9, &rot_info); break;
		case 175: get_hex_rotation_info(13, &rot_info); break;
		case 126: get_hex_rotation_info(22, &rot_info); break;
		case 123: get_hex_rotation_info(21, &rot_info); break;
		case  95: get_hex_rotation_info(15, &rot_info); break;
		}
	} else {
		switch(index) {
		case 250: amb_index = separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		case 245: amb_index = separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(10, &rot_info);
		break;
		case 237: amb_index = separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 4, &rot_info);
		break;
		case 231: amb_index = separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(11, &rot_info);
		break;
		case 222: amb_index = separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 1, &rot_info);
		break;
		case 219: amb_index = separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 6, &rot_info);
		break;
		case 189: amb_index = separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(18, &rot_info);
		break;
		case 183: amb_index = separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 9, &rot_info);
		break;
		case 175: amb_index = separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(13, &rot_info);
		break;
		case 126: amb_index = separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(22, &rot_info);
		break;
		case 123: amb_index = separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(21, &rot_info);
		break;
		case  95: amb_index = separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(15, &rot_info);
		break;
		}
	}

	cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	if (!amb_index) {
		*cube_polys->isosurf = isosurf_A;
		/*    *cube_polys->bounds[rot_info.face[0]] = bounds_0[0];*/
	} else {
		*cube_polys->isosurf = isosurf_B;
		/*    *cube_polys->bounds[rot_info.face[0]] = bounds_0[1];*/
	}

	/*  *cube_polys->bounds[rot_info.face[1]] = bounds_1;
	 *cube_polys->bounds[rot_info.face[2]] = bounds_2;
	 *cube_polys->bounds[rot_info.face[3]] = bounds_3;
	 *cube_polys->bounds[rot_info.face[4]] = bounds_4;
	 *cube_polys->bounds[rot_info.face[5]] = bounds_5;
	 */
}

void case4_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[8] = {0, 3, 5, 6, 8, 11, 200, 206};
	static int 	isosurf[7] = {0,4,1,  5,3,2, -1};
	/*  static int 	bounds[6][5] = {{3, 6,0,1, -1},{3, 7,5,2, -1},
				{3, 7,2,3, -1},{3, 6,1,4, -1},
				{3, 7,3,5, -1},{3, 6,4,0, -1}};
	 */
	switch(index) {
	case  20: get_hex_rotation_info(6, &rot_info);  break;
	case  40: get_hex_rotation_info(9, &rot_info);  break;
	case  65: get_hex_rotation_info(0, &rot_info);  break;
	case 130: get_hex_rotation_info(4, &rot_info);  break;
	}

	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
		/*    *cube_polys->bounds[rot_info.face[i]] = bounds[i];*/
	}
}

void case4opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[12] = {0, 3, 5, 6, 8, 11, 201, 202, 203, 204, 205, 207};
	static int 	isosurf[7] = {0,1,4,  5,2,3, -1};
	/*  static int 	bounds[6][7] = {{5, 6,7,8,1,0, -1},{5, 7,6,10,2,5, -1},
				{5, 10,9,11,3,2, -1},{5, 8,11,9,4,1, -1},
				{5, 11,8,7,5,3, -1},{5, 9,10,6,0,4, -1}};
	 */
	switch(index) {
	case 235: get_hex_rotation_info( 6, &rot_info);  break;
	case 215: get_hex_rotation_info( 9, &rot_info);  break;
	case 190: get_hex_rotation_info( 0, &rot_info);  break;
	case 125: get_hex_rotation_info( 4, &rot_info);  break;
	}

	cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
		/*    *cube_polys->bounds[rot_info.face[i]] = bounds[i];*/
	}
}

void case5_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[8] = {0, 1, 5, 7, 8, 201, 204, 205};
	static int 	isosurf[10] = {0,1,4,  1,3,4,  1,2,3, -1};
	/*  static int 	bounds[6][7] = {{3, 5,1,0, -1,-1,-1},{4, 5,7,2,1, -1,-1},
				{4, 7,6,3,2, -1,-1},{3, 6,4,3, -1,-1,-1},
				{-1, -1,-1,-1,-1,-1,-1},{5, 6,7,5,0,4, -1}};
	 */
	if (index < 80) {
		switch(index) {
		case  7: get_hex_rotation_info( 9, &rot_info);  break;
		case 11: get_hex_rotation_info( 6, &rot_info);  break;
		case 13: get_hex_rotation_info( 4, &rot_info);  break;
		case 14: get_hex_rotation_info( 2, &rot_info);  break;
		case 19: get_hex_rotation_info(15, &rot_info);  break;
		case 25: get_hex_rotation_info(21, &rot_info);  break;
		case 35: get_hex_rotation_info(14, &rot_info);  break;
		case 38: get_hex_rotation_info(20, &rot_info);  break;
		case 49: get_hex_rotation_info( 3, &rot_info);  break;
		case 50: get_hex_rotation_info( 0, &rot_info);  break;
		case 70: get_hex_rotation_info(16, &rot_info);  break;
		case 76: get_hex_rotation_info(23, &rot_info);  break;
		}
	} else {
		switch(index) {
		case  98: get_hex_rotation_info( 8, &rot_info);  break;
		case 100: get_hex_rotation_info( 5, &rot_info);  break;
		case 112: get_hex_rotation_info(22, &rot_info);  break;
		case 137: get_hex_rotation_info(13, &rot_info);  break;
		case 140: get_hex_rotation_info(18, &rot_info);  break;
		case 145: get_hex_rotation_info(10, &rot_info);  break;
		case 152: get_hex_rotation_info( 1, &rot_info);  break;
		case 176: get_hex_rotation_info(19, &rot_info);  break;
		case 196: get_hex_rotation_info(11, &rot_info);  break;
		case 200: get_hex_rotation_info( 7, &rot_info);  break;
		case 208: get_hex_rotation_info(17, &rot_info);  break;
		case 224: get_hex_rotation_info(12, &rot_info);  break;
		}
	}

	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
		/*    *cube_polys->bounds[rot_info.face[i]] = bounds[i];*/
	}
}

void case5opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 1, 5, 7, 8, 200, 202, 203, 206, 207};
	static int 	isosurf[10] = {0,4,1,  1,4,3,  1,3,2, -1};
	/*  static int 	bounds[6][7] = {{5, 6,7,5,0,1, -1},{4, 8,6,1,2, -1,-1},
				{4, 9,8,2,3, -1,-1},{5, 5,7,9,3,4, -1},
				{4, 6,8,9,7, -1,-1},{3, 5,4,0, -1,-1,-1}};
	 */
	if (index > 170) {
		switch(index) {
		case 248: get_hex_rotation_info( 9, &rot_info);  break;
		case 244: get_hex_rotation_info( 6, &rot_info);  break;
		case 242: get_hex_rotation_info( 4, &rot_info);  break;
		case 241: get_hex_rotation_info( 2, &rot_info);  break;
		case 236: get_hex_rotation_info(15, &rot_info);  break;
		case 230: get_hex_rotation_info(21, &rot_info);  break;
		case 220: get_hex_rotation_info(14, &rot_info);  break;
		case 217: get_hex_rotation_info(20, &rot_info);  break;
		case 206: get_hex_rotation_info( 3, &rot_info);  break;
		case 205: get_hex_rotation_info( 0, &rot_info);  break;
		case 185: get_hex_rotation_info(16, &rot_info);  break;
		case 179: get_hex_rotation_info(23, &rot_info);  break;
		}
	} else {
		switch(index) {
		case 157: get_hex_rotation_info( 8, &rot_info);  break;
		case 155: get_hex_rotation_info( 5, &rot_info);  break;
		case 143: get_hex_rotation_info(22, &rot_info);  break;
		case 118: get_hex_rotation_info(13, &rot_info);  break;
		case 115: get_hex_rotation_info(18, &rot_info);  break;
		case 110: get_hex_rotation_info(10, &rot_info);  break;
		case 103: get_hex_rotation_info( 1, &rot_info);  break;
		case  79: get_hex_rotation_info(19, &rot_info);  break;
		case  59: get_hex_rotation_info(11, &rot_info);  break;
		case  55: get_hex_rotation_info( 7, &rot_info);  break;
		case  47: get_hex_rotation_info(17, &rot_info);  break;
		case  31: get_hex_rotation_info(12, &rot_info);  break;
		}
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
		/*    *cube_polys->bounds[rot_info.face[i]] = bounds[i];*/
	}
}

void case6_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[10] = {1, 3, 5, 6, 8, 9, 11, 200, 201, 206};
	static int 	isosurf_A[16] = {0,6,1,  1,6,3,  1,3,4,  4,3,2,  4,2,5, -1};
	static int 	isosurf_B[10] = {1,0,4,  4,0,5,  6,3,2, -1};
	/*  static int 	bounds_0[6]  = {4, 7,8,0,1, -1};
  static int 	bounds_1[2][9] = {{6, 8,5,2,9,6,0, -1,-1},
				  {3, 8,5,0,  3, 9,6,2, -1}};
  static int 	bounds_2[5] = {3, 9,2,3, -1};
  static int 	bounds_3[5] = {3, 7,1,4, -1};
  static int 	bounds_4[5] = {3, 9,3,6, -1};
  static int 	bounds_5[6] = {4, 8,7,4,5, -1};
	 */
	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;	/* isosurf_B */
		if (index < 82) {
			switch(index) {
			case 21: get_hex_rotation_info(14, &rot_info); break;
			case 22: get_hex_rotation_info( 6, &rot_info); break;
			case 28: get_hex_rotation_info( 7, &rot_info); break;
			case 41: get_hex_rotation_info(10, &rot_info); break;
			case 42: get_hex_rotation_info(16, &rot_info); break;
			case 44: get_hex_rotation_info( 9, &rot_info); break;
			case 52: get_hex_rotation_info(12, &rot_info); break;
			case 56: get_hex_rotation_info(15, &rot_info); break;
			case 67: get_hex_rotation_info( 0, &rot_info); break;
			case 69: get_hex_rotation_info(18, &rot_info); break;
			case 73: get_hex_rotation_info( 2, &rot_info); break;
			case 81: get_hex_rotation_info( 1, &rot_info); break;
			}
		} else {
			switch(index) {
			case  84: get_hex_rotation_info( 8, &rot_info); break;
			case  97: get_hex_rotation_info(20, &rot_info); break;
			case 104: get_hex_rotation_info(17, &rot_info); break;
			case 131: get_hex_rotation_info( 4, &rot_info); break;
			case 134: get_hex_rotation_info( 5, &rot_info); break;
			case 138: get_hex_rotation_info(21, &rot_info); break;
			case 146: get_hex_rotation_info(22, &rot_info); break;
			case 148: get_hex_rotation_info(13, &rot_info); break;
			case 162: get_hex_rotation_info( 3, &rot_info); break;
			case 168: get_hex_rotation_info(11, &rot_info); break;
			case 193: get_hex_rotation_info(19, &rot_info); break;
			case 194: get_hex_rotation_info(23, &rot_info); break;
			}
		}
	} else {
		if (index < 82) {
			switch(index) {
			case 21: amb_index = separation_test(mesh, cell,  0,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(14, &rot_info);
			break;
			case 22: amb_index = separation_test(mesh, cell, 10,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 6, &rot_info);
			break;
			case 28: amb_index = separation_test(mesh, cell,  6,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 7, &rot_info);
			break;
			case 41: amb_index = separation_test(mesh, cell, 11,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(10, &rot_info);
			break;
			case 42: amb_index = separation_test(mesh, cell,  1,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(16, &rot_info);
			break;
			case 44: amb_index = separation_test(mesh, cell,  3,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 9, &rot_info);
			break;
			case 52: amb_index = separation_test(mesh, cell,  3,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(12, &rot_info);
			break;
			case 56: amb_index = separation_test(mesh, cell,  6,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(15, &rot_info);
			break;
			case 67: amb_index = separation_test(mesh, cell,  2,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 0, &rot_info);
			break;
			case 69: amb_index = separation_test(mesh, cell,  0,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(18, &rot_info);
			break;
			case 73: amb_index = separation_test(mesh, cell,  8,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 2, &rot_info);
			break;
			case 81: amb_index = separation_test(mesh, cell,  5,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 1, &rot_info);
			break;
			}
		} else {
			switch(index) {
			case  84: amb_index = separation_test(mesh, cell,  5,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 8, &rot_info);
			break;
			case  97: amb_index = separation_test(mesh, cell, 11,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(20, &rot_info);
			break;
			case 104: amb_index = separation_test(mesh, cell,  8,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(17, &rot_info);
			break;
			case 131: amb_index = separation_test(mesh, cell,  7,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 4, &rot_info);
			break;
			case 134: amb_index = separation_test(mesh, cell,  9,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 5, &rot_info);
			break;
			case 138: amb_index = separation_test(mesh, cell,  1,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(21, &rot_info);
			break;
			case 146: amb_index = separation_test(mesh, cell, 10,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(22, &rot_info);
			break;
			case 148: amb_index = separation_test(mesh, cell,  9,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(13, &rot_info);
			break;
			case 162: amb_index = separation_test(mesh, cell,  4,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 3, &rot_info);
			break;
			case 168: amb_index = separation_test(mesh, cell,  4,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(11, &rot_info);
			break;
			case 193: amb_index = separation_test(mesh, cell,  7,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(19, &rot_info);
			break;
			case 194: amb_index = separation_test(mesh, cell,  2,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(23, &rot_info);
			break;
			}
		}
	}

	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	/*  *cube_polys->bounds[rot_info.face[0]] = bounds_0;*/

	if (amb_index) {
		*cube_polys->isosurf = isosurf_A;
		/*    *cube_polys->bounds[rot_info.face[1]] = bounds_1[0];*/
	} else {
		*cube_polys->isosurf = isosurf_B;
		/*    *cube_polys->bounds[rot_info.face[1]] = bounds_1[1];*/
	}

	/*  *cube_polys->bounds[rot_info.face[2]] = bounds_2;
	 *cube_polys->bounds[rot_info.face[3]] = bounds_3;
	 *cube_polys->bounds[rot_info.face[4]] = bounds_4;
	 *cube_polys->bounds[rot_info.face[5]] = bounds_5;
	 */
}

void case6opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[12] = {1, 3, 5, 6, 8, 9, 11, 202, 203, 204, 205, 207};
	static int 	isosurf_A[16] = {0,1,6,  1,3,6,  1,4,3,  4,2,3,  4,5,2, -1};
	static int 	isosurf_B[10] = {1,4,0,  4,5,0,  6,2,3, -1};
	/*  static int 	bounds_0[6]  = {4, 7,8,1,0, -1};
  static int 	bounds_1[2][9] = {{3, 7,0,6,  3, 10,2,5, -1},
				  {6, 7,0,5,10,2,6, -1,-1}};
  static int 	bounds_2[7] = {5, 10,9,11,3,2, -1};
  static int 	bounds_3[7] = {5, 8,11,9,4,1, -1};
  static int 	bounds_4[7] = {5, 11,8,7,6,3, -1};
  static int 	bounds_5[6] = {4, 9,10,5,4, -1};
	 */
	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;	/* isosurf_B */
		if (index > 173) {
			switch(index) {
			case 234: get_hex_rotation_info(14, &rot_info); break;
			case 233: get_hex_rotation_info( 6, &rot_info); break;
			case 227: get_hex_rotation_info( 7, &rot_info); break;
			case 214: get_hex_rotation_info(10, &rot_info); break;
			case 213: get_hex_rotation_info(16, &rot_info); break;
			case 211: get_hex_rotation_info( 9, &rot_info); break;
			case 203: get_hex_rotation_info(12, &rot_info); break;
			case 199: get_hex_rotation_info(15, &rot_info); break;
			case 188: get_hex_rotation_info( 0, &rot_info); break;
			case 186: get_hex_rotation_info(18, &rot_info); break;
			case 182: get_hex_rotation_info( 2, &rot_info); break;
			case 174: get_hex_rotation_info( 1, &rot_info); break;
			}
		} else {
			switch(index) {
			case 171: get_hex_rotation_info( 8, &rot_info); break;
			case 158: get_hex_rotation_info(20, &rot_info); break;
			case 151: get_hex_rotation_info(17, &rot_info); break;
			case 124: get_hex_rotation_info( 4, &rot_info); break;
			case 121: get_hex_rotation_info( 5, &rot_info); break;
			case 117: get_hex_rotation_info(21, &rot_info); break;
			case 109: get_hex_rotation_info(22, &rot_info); break;
			case 107: get_hex_rotation_info(13, &rot_info); break;
			case  93: get_hex_rotation_info( 3, &rot_info); break;
			case  87: get_hex_rotation_info(11, &rot_info); break;
			case  62: get_hex_rotation_info(19, &rot_info); break;
			case  61: get_hex_rotation_info(23, &rot_info); break;
			}
		}
	} else {
		if (index > 173) {
			switch(index) {
			case 234: amb_index = separation_test(mesh, cell,  0,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(14, &rot_info);
			break;
			case 233: amb_index = separation_test(mesh, cell, 10,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 6, &rot_info);
			break;
			case 227: amb_index = separation_test(mesh, cell,  6,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 7, &rot_info);
			break;
			case 214: amb_index = separation_test(mesh, cell, 11,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(10, &rot_info);
			break;
			case 213: amb_index = separation_test(mesh, cell,  1,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(16, &rot_info);
			break;
			case 211: amb_index = separation_test(mesh, cell,  3,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 9, &rot_info);
			break;
			case 203: amb_index = separation_test(mesh, cell,  3,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(12, &rot_info);
			break;
			case 199: amb_index = separation_test(mesh, cell,  6,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(15, &rot_info);
			break;
			case 188: amb_index = separation_test(mesh, cell,  2,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 0, &rot_info);
			break;
			case 186: amb_index = separation_test(mesh, cell,  0,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(18, &rot_info);
			break;
			case 182: amb_index = separation_test(mesh, cell,  8,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 2, &rot_info);
			break;
			case 174: amb_index = separation_test(mesh, cell,  5,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 1, &rot_info);
			break;
			}
		} else {
			switch(index) {
			case 171: amb_index = separation_test(mesh, cell,  5,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 8, &rot_info);
			break;
			case 158: amb_index = separation_test(mesh, cell, 11,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(20, &rot_info);
			break;
			case 151: amb_index = separation_test(mesh, cell,  8,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(17, &rot_info);
			break;
			case 124: amb_index = separation_test(mesh, cell,  7,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 4, &rot_info);
			break;
			case 121: amb_index = separation_test(mesh, cell,  9,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 5, &rot_info);
			break;
			case 117: amb_index = separation_test(mesh, cell,  1,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(21, &rot_info);
			break;
			case 109: amb_index = separation_test(mesh, cell, 10,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(22, &rot_info);
			break;
			case 107: amb_index = separation_test(mesh, cell,  9,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(13, &rot_info);
			break;
			case  93: amb_index = separation_test(mesh, cell,  4,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info( 3, &rot_info);
			break;
			case  87: amb_index = separation_test(mesh, cell,  4,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(11, &rot_info);
			break;
			case  62: amb_index = separation_test(mesh, cell,  7,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(19, &rot_info);
			break;
			case  61: amb_index = separation_test(mesh, cell,  2,
					fvalue, CS_type, disamb_flag);
			get_hex_rotation_info(23, &rot_info);
			break;
			}
		}
	}

	cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	/*  *cube_polys->bounds[rot_info.face[0]] = bounds_0;*/

	if (amb_index) {
		*cube_polys->isosurf = isosurf_A;
		/*    *cube_polys->bounds[rot_info.face[1]] = bounds_1[0];*/
	} else {
		*cube_polys->isosurf = isosurf_B;
		/*    *cube_polys->bounds[rot_info.face[1]] = bounds_1[1];*/
	}

	/*  *cube_polys->bounds[rot_info.face[2]] = bounds_2;
	 *cube_polys->bounds[rot_info.face[3]] = bounds_3;
	 *cube_polys->bounds[rot_info.face[4]] = bounds_4;
	 *cube_polys->bounds[rot_info.face[5]] = bounds_5;
	 */
}

void case7_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[13] = {0,1,2,3,5,6,9,10,11,201,203,206,101};
	static int 	isosurf_A[10] = {0,1,6,  3,7,2,  8,5,4, -1};
	static int 	isosurf_B[16] = {8,5,4,  3,7,0,  0,7,6,  6,7,2,  6,2,1, -1};
	static int 	isosurf_C[28] = {6,0,12,  0,3,12,  3,7,12,  7,2,12,  2,1,12,
			1,8,12,  8,5,12,  5,4,12,  4,6,12, -1};
	static int 	isosurf_D[16] = {8,2,1,  0,3,7,  0,7,5,  0,5,6,  6,5,4, -1};



	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;
		switch(index) {
		case  26: get_hex_rotation_info(16, &rot_info); break;
		case  37: get_hex_rotation_info(10, &rot_info); break;
		case  74: get_hex_rotation_info( 0, &rot_info); break;
		case  82: get_hex_rotation_info( 6, &rot_info); break;
		case  88: get_hex_rotation_info( 2, &rot_info); break;
		case 133: get_hex_rotation_info(14, &rot_info); break;
		case 161: get_hex_rotation_info( 4, &rot_info); break;
		case 164: get_hex_rotation_info( 9, &rot_info); break;
		}
	} else {
		amb_index = 0;
		switch(index) {
		case  26: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(16, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 7, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(22, &rot_info);
			break;
		}
		break;
		case  37: amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<=1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(10, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(12, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(18, &rot_info);
			break;
		}
		break;
		case  74: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 0, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(17, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(21, &rot_info);
			break;
		}
		break;
		case  82: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		amb_index <<=1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 6, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 1, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(23, &rot_info);
			break;
		}
		break;
		case  88: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 2, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 8, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(15, &rot_info);
			break;
		}
		break;
		case 133: amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(14, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 5, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(19, &rot_info);
			break;
		}
		break;
		case 161: amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 4, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(11, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(20, &rot_info);
			break;
		}
		break;
		case 164: amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 9, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 3, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(13, &rot_info);
			break;
		}
		break;
		}
	}

	switch(amb_index) {
	case 0:
		cube_polys->verts_num = 12;
		for (i = 0; i < 12; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_A;

		break;
	case 1: case 2: case 4:
		cube_polys->verts_num = 12;
		for (i = 0; i < 12; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_B;

		break;
	case 3: case 5: case 6:
		cube_polys->verts_num = 13;
		for (i = 0; i < 13; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_C;

		break;
	case 7:
		cube_polys->verts_num = 12;
		for (i = 0; i < 12; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_D;

		break;
	}
}

void case7opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[15] = {0,1,2,3,5,6,9,10,11,200,202,204,205,207,101};
	static int 	isosurf_A[10] = {0,6,1,  3,2,7,  8,4,5, -1};
	static int 	isosurf_B[16] = {8,4,5,  3,0,7,  0,6,7,  6,2,7,  6,1,2, -1};
	static int 	isosurf_C[28] = {0,6,14,  3,0,14,  7,3,14,  2,7,14,  1,2,14,
			8,1,14,  5,8,14,  4,5,14,  6,4,14, -1};
	static int 	isosurf_D[16] = {8,1,2,  0,7,3,  0,5,7,  0,6,5,  6,4,5, -1};


	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;
		switch(index) {
		case 229: get_hex_rotation_info(16, &rot_info); break;
		case 218: get_hex_rotation_info(10, &rot_info); break;
		case 181: get_hex_rotation_info( 0, &rot_info); break;
		case 173: get_hex_rotation_info( 6, &rot_info); break;
		case 167: get_hex_rotation_info( 2, &rot_info); break;
		case 122: get_hex_rotation_info(14, &rot_info); break;
		case  94: get_hex_rotation_info( 4, &rot_info); break;
		case  91: get_hex_rotation_info( 9, &rot_info); break;
		}
	} else {
		amb_index = 0;
		switch(index) {
		case 229: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(16, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 7, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(22, &rot_info);
			break;
		}
		break;
		case 218: amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(10, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(12, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(18, &rot_info);
			break;
		}
		break;
		case 181: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 0, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(17, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(21, &rot_info);
			break;
		}
		break;
		case 173: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 6, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 1, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(23, &rot_info);
			break;
		}
		break;
		case 167: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 2, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 8, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(15, &rot_info);
			break;
		}
		break;
		case 122: amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info(14, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 5, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(19, &rot_info);
			break;
		}
		break;
		case  94: amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 4, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info(11, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(20, &rot_info);
			break;
		}
		break;
		case  91: amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);
		switch(amb_index) {
		case 0: case 1: case 3: case 7:
			get_hex_rotation_info( 9, &rot_info);
			break;
		case 2: case 6:
			get_hex_rotation_info( 3, &rot_info);
			break;
		case 4: case 5:
			get_hex_rotation_info(13, &rot_info);
			break;
		}
		break;
		}
	}

	switch(amb_index) {
	case 0:
		cube_polys->verts_num = 14;
		for (i = 0; i < 14; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_A;

		break;
	case 1: case 2: case 4:
		cube_polys->verts_num = 14;
		for (i = 0; i < 14; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_B;

		break;
	case 3: case 5: case 6:
		cube_polys->verts_num = 15;
		for (i = 0; i < 15; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_C;

		break;
	case 7:
		cube_polys->verts_num = 14;
		for (i = 0; i < 14; i++) {
			cube_polys->verts[i] = point_convert(verts[i], &rot_info);
		}
		*cube_polys->isosurf = isosurf_D;

		break;
	}
}


void case8_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[8] = {1, 3, 5, 7, 200, 201, 204, 205};
	static int 	isosurf[7] = {0,2,1,  1,2,3, -1};

	switch(index) {
	case  15: get_hex_rotation_info( 2, &rot_info);  break;
	case  51: get_hex_rotation_info( 0, &rot_info);  break;
	case 102: get_hex_rotation_info( 5, &rot_info);  break;
	}
	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {

	}
}


void case8opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[8] = {1, 3, 5, 7, 202, 203, 206, 207};
	static int 	isosurf[7] = {0,1,2,  1,3,2, -1};


	switch(index) {
	case 240: get_hex_rotation_info( 2, &rot_info);  break;
	case 204: get_hex_rotation_info( 0, &rot_info);  break;
	case 153: get_hex_rotation_info( 5, &rot_info);  break;
	}
	cube_polys->verts_num = 8;
	for (i = 0; i < 8; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
	}
}

void case9_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 3, 5, 6, 9, 10, 200, 204, 205, 207};
	static int 	isosurf[13] = {1,3,5,  0,3,1,  0,2,3,  0,4,2, -1};

	switch(index) {
	case 141: get_hex_rotation_info( 1, &rot_info);  break;
	case 177: get_hex_rotation_info( 0, &rot_info);  break;
	case 216: get_hex_rotation_info(10,&rot_info);  break;
	case 228: get_hex_rotation_info( 7, &rot_info);  break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
	}
}


void case9opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 3, 5, 6, 9, 10, 201, 202, 203, 206};
	static int 	isosurf[13] = {1,5,3,  0,1,3,  0,3,2,  0,2,4, -1};


	switch(index) {
	case 114: get_hex_rotation_info( 1, &rot_info);  break;
	case  78: get_hex_rotation_info( 0, &rot_info);  break;
	case  39: get_hex_rotation_info(10,&rot_info);  break;
	case  27: get_hex_rotation_info( 7, &rot_info);  break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {
	}
}

void case10_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[13] = {0, 2, 4, 6, 8, 9, 10, 11, 200, 203, 205, 206, 106};
	static int 	isosurf_A[13] = {4,1,0,  1,4,6,  5,7,3,  3,2,5, -1};
	static int 	isosurf_B[25] = {1,0,12,  7,1,12,  5,7,12,  2,5,12,
			3,2,12,  6,3,12,  4,6,12,  0,4,12, -1};
	static int 	isosurf_C[13] = {6,3,2,  6,2,4,  7,1,0,  7,0,5, -1};
	static int 	isosurf_D[25] = {0,5,12,  5,7,12,  7,3,12,  3,2,12,
			2,4,12,  4,6,12,  6,1,12,  1,0,12, -1};



	if (disamb_flag == 0) {  /* disambiguation : none */
		amb_index = 0;
	switch(index) {
	case  60: get_hex_rotation_info(10, &rot_info); break;
	case  85: get_hex_rotation_info( 2, &rot_info); break;
	case 105: get_hex_rotation_info( 0, &rot_info); break;
	}
	}else {
		amb_index = 0;
		switch(index) {
		case  60: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(10, &rot_info);
		break;
		case  85: amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 2, &rot_info);
		break;
		case 105: amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		}
	}

	switch(amb_index) {
	case 0: cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_A;

	break;
	case 1: cube_polys->verts_num = 13;
	for (i = 0; i < 13; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_B;

	break;
	case 3: cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_C;

	break;
	case 2: cube_polys->verts_num = 13;
	for (i = 0; i < 13; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_D;

	break;
	}
}

void case10opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[13] = {0, 2, 4, 6, 8, 9, 10, 11, 201, 202, 204, 207, 106};
	static int 	isosurf_A[13] = {4,0,1,  1,6,4,  5,3,7,  3,5,2, -1};
	static int 	isosurf_B[25] = {0,1,12,  1,7,12,  7,5,12,  5,2,12,
			2,3,12,  3,6,12,  6,4,12,  4,0,12, -1};
	static int 	isosurf_C[13] = {6,2,3,  6,4,2,  7,0,1,  7,5,0, -1};
	static int 	isosurf_D[25] = {5,0,12,  7,5,12,  3,7,12,  2,3,12,
			4,2,12,  6,4,12,  1,6,12,  0,1,12, -1};



	if (disamb_flag == 0) {  /* disambiguation : none */
		amb_index = 0;
		switch(index) {
		case 195: get_hex_rotation_info(10, &rot_info); break;
		case 170: get_hex_rotation_info( 2, &rot_info); break;
		case 150: get_hex_rotation_info( 0, &rot_info); break;
		}
	} else {
		amb_index = 0;
		switch(index) {
		case 195: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(10, &rot_info);
		break;
		case 170: amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 2, &rot_info);
		break;
		case 150: amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		}
	}

	switch(amb_index) {
	case 0: cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_A;

	break;
	case 1: cube_polys->verts_num = 13;
	for (i = 0; i < 13; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_B;

	break;
	case 3: cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_C;

	break;
	case 2: cube_polys->verts_num = 13;
	for (i = 0; i < 13; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_D;

	break;
	}
}

void case11_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 3, 6, 7, 9, 11, 200, 204, 205, 206};
	static int 	isosurf[13] = {0,3,1,  4,3,0,  4,2,3,  4,5,2, -1};


	switch(index) {
	case  29: get_hex_rotation_info(13, &rot_info); break;
	case  43: get_hex_rotation_info( 9, &rot_info); break;
	case  54: get_hex_rotation_info( 8, &rot_info); break;
	case  71: get_hex_rotation_info(20, &rot_info); break;
	case 108: get_hex_rotation_info(11, &rot_info); break;
	case 113: get_hex_rotation_info( 0, &rot_info); break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {

	}
}


void case11opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 3, 6, 7, 9, 11, 201, 202, 203, 207};
	static int 	isosurf[13] = {0,1,3,  4,0,3,  4,3,2,  4,2,5, -1};


	switch(index) {
	case 226: get_hex_rotation_info(13, &rot_info); break;
	case 212: get_hex_rotation_info( 9, &rot_info); break;
	case 201: get_hex_rotation_info( 8, &rot_info); break;
	case 184: get_hex_rotation_info(20, &rot_info); break;
	case 147: get_hex_rotation_info(11, &rot_info); break;
	case 142: get_hex_rotation_info( 0, &rot_info); break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {

	}
}

void case12_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[12] = {0, 1, 2, 3, 5, 7, 8, 10, 201, 203, 204, 205};
	static int 	isosurf_A[13] = {7,2,3,  0,1,6,  1,5,6,  1,4,5, -1};
	static int 	isosurf_B[19] = {4,2,1,  4,7,2,  4,0,7,  0,3,7,
			4,6,0,  4,5,6, -1};
	static int 	isosurf_C[13] = {0,3,6,  2,5,7,  1,5,2,  4,5,1, -1};
	static int 	isosurf_D[19] = {2,3,6,  4,5,7,  4,7,2,  4,2,6,
			4,6,0,  4,0,1, -1};


	if (disamb_flag == 0) {  /* disambiguation : none */
		amb_index = 0;	/* isosurf_A */
		switch(index) {
		case  30: get_hex_rotation_info( 2, &rot_info); break;
		case  45: get_hex_rotation_info( 4, &rot_info); break;
		case  53: get_hex_rotation_info( 3, &rot_info); break;
		case  58: get_hex_rotation_info( 0, &rot_info); break;
		case  75: get_hex_rotation_info( 6, &rot_info); break;
		case  83: get_hex_rotation_info(15, &rot_info); break;
		case  86: get_hex_rotation_info(16, &rot_info); break;
		case  89: get_hex_rotation_info(21, &rot_info); break;
		case  92: get_hex_rotation_info(23, &rot_info); break;
		case 101: get_hex_rotation_info( 5, &rot_info); break;
		case 106: get_hex_rotation_info( 8, &rot_info); break;
		case 120: get_hex_rotation_info(22, &rot_info); break;
		}
	} else {
		amb_index = 0;
		switch(index) {
		case  30: amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 2, &rot_info);
		break;
		case  45: amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 4, &rot_info);
		break;
		case  53: amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 3, &rot_info);
		break;
		case  58: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		case  75: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 6, &rot_info);
		break;
		case  83: amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(15, &rot_info);
		break;
		case  86: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(16, &rot_info);
		break;
		case  89: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(21, &rot_info);
		break;
		case  92: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(23, &rot_info);
		break;
		case 101: amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 5, &rot_info);
		break;
		case 106: amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 8, &rot_info);
		break;
		case 120: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(22, &rot_info);
		break;
		}
	}

	cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}

	switch(amb_index) {
	case 0: *cube_polys->isosurf = isosurf_A;

	break;
	case 1:     *cube_polys->isosurf = isosurf_B;

	break;
	case 3: *cube_polys->isosurf = isosurf_C;

	break;
	case 2:
		*cube_polys->isosurf = isosurf_D;

		break;
	}
}

void case12opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[12] = {0, 1, 2, 3, 5, 7, 8, 10, 200, 202, 206, 207};
	static int 	isosurf_A[13] = {7,3,2,  0,6,1,  1,6,5,  1,5,4, -1};
	static int 	isosurf_B[19] = {4,1,2,  4,2,7,  4,7,0,  0,7,3,
			4,0,6,  4,6,5, -1};
	static int 	isosurf_C[13] = {0,6,3,  2,7,5,  1,2,5,  4,1,5, -1};
	static int 	isosurf_D[19] = {2,6,3,  4,7,5,  4,2,7,  4,6,2,
			4,0,6,  4,1,0, -1};



	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 0;
		switch(index) {
		case 225: get_hex_rotation_info( 2, &rot_info); break;
		case 210: get_hex_rotation_info( 4, &rot_info); break;
		case 202: get_hex_rotation_info( 3, &rot_info); break;
		case 197: get_hex_rotation_info( 0, &rot_info); break;
		case 180: get_hex_rotation_info( 6, &rot_info); break;
		case 172: get_hex_rotation_info(15, &rot_info); break;
		case 169: get_hex_rotation_info(16, &rot_info); break;
		case 166: get_hex_rotation_info(21, &rot_info); break;
		case 163: get_hex_rotation_info(23, &rot_info); break;
		case 154: get_hex_rotation_info( 5, &rot_info); break;
		case 149: get_hex_rotation_info( 8, &rot_info); break;
		case 135: get_hex_rotation_info(22, &rot_info); break;
		}
	} else {
		amb_index = 0;
		switch(index) {
		case 225: amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 2, &rot_info);
		break;
		case 210: amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 4, &rot_info);
		break;
		case 202: amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 3, &rot_info);
		break;
		case 197: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 0, &rot_info);
		break;
		case 180: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 6, &rot_info);
		break;
		case 172: amb_index |= separation_test(mesh, cell,  2,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(15, &rot_info);
		break;
		case 169: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell, 10,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(16, &rot_info);
		break;
		case 166: amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(21, &rot_info);
		break;
		case 163: amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  5,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(23, &rot_info);
		break;
		case 154: amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 5, &rot_info);
		break;
		case 149: amb_index |= separation_test(mesh, cell,  1,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info( 8, &rot_info);
		break;
		case 135: amb_index |= separation_test(mesh, cell,  8,
				fvalue, CS_type, disamb_flag);
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  6,
				fvalue, CS_type, disamb_flag);
		get_hex_rotation_info(22, &rot_info);
		break;
		}
	}

	cube_polys->verts_num = 12;
	for (i = 0; i < 12; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}

	switch(amb_index) {
	case 0: *cube_polys->isosurf = isosurf_A;

	break;
	case 1:
		*cube_polys->isosurf = isosurf_B;

		break;
	case 3: *cube_polys->isosurf = isosurf_C;

	break;
	case 2:
		*cube_polys->isosurf = isosurf_D;

		break;
	}
}


void case13_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
			200, 202, 205, 207};
	static int 	isosurf_A[13] = {2,10,3,  0,9,1,  8,7,4,  5,6,11, -1};
	static int 	isosurf_B[19] = {8,7,4,  5,6,11,  9,1,10,  1,2,10,
			9,10,3,  9,3,0, -1};
	static int 	isosurf_C[25] = {0,9,3,  9,10,3,  10,9,1,  1,2,10,
			5,4,11,  11,4,8,  11,8,6,  6,8,7, -1};
	static int 	isosurf_D[31] = {0,9,16,  9,5,16,  5,6,16,  6,11,16,
			11,1,16,  1,2,16,  2,10,16,  10,3,16,
			3,0,16,  4,8,7, -1};
	static int 	isosurf_E[37] = {0,9,16,  9,5,16,  5,6,16,  6,11,16,
			11,1,16,  1,2,16,  2,10,16,  10,7,16,
			7,4,16,  4,8,16,  8,3,16,  3,0,16, -1};
	static int 	isosurf_F[37] = {0,8,16,  8,7,16,  7,4,16,  4,9,16,
			9,5,16,  5,6,16,  6,11,16,  11,1,16,
			1,2,16,  2,10,16,  10,3,16,  3,0,16, -1};
	static int 	isosurf_G[25] = {0,8,1,  1,8,11,  11,8,3,  3,2,11,
			9,7,4,  9,10,7,  9,5,10,  5,6,10, -1};
	static int 	isosurf_H[31] = {0,8,16,  8,3,16,  3,2,16,  2,11,16,
			11,5,16,  5,4,16,  4,9,16,  9,1,16,
			1,0,16,  6,10,7, -1};
	static int 	isosurf_I[19] = {8,3,2,  8,2,11,  8,11,0, 0,11,1,
			9,5,4,  10,7,6, -1};
	static int 	isosurf_J[13] = {0,8,3,  4,9,5,  2,11,1,  10,7,6, -1};
	static int 	isosurf_K[37] = {0,8,16,  8,7,16,  7,4,16,  4,9,16,
			9,1,16,  1,2,16,  2,11,16,  11,5,16,
			5,6,16,  6,10,16,  10,3,16,  3,0,16, -1};
	static int 	isosurf_L[19] = {2,11,1,  0,10,3,  0,6,10,  9,6,0,
			9,5,6,  8,7,4, -1};
	/*
  static int 	sep_bd[6][9] = {{3, 12,0,3,  3, 13,2,1, -1},
				{3, 14,5,9,  3, 13,1,11, -1},
				{3, 14,4,5,  3, 15,6,7, -1},
				{3, 12,3,8,  3, 15,7,10, -1},
				{3, 13,11,2,  3, 15,10,6, -1},
				{3, 14,9,4,  3, 12,8,0, -1}};
  static int 	not_sep_bd[6][8] = {{6, 12,0,1,13,2,3, -1},
				    {6, 14,5,11,13,1,9, -1},
				    {6, 14,4,7,15,6,5, -1},
				    {6, 12,3,10,15,7,8, -1},
				    {6, 13,11,6,15,10,2, -1},
				    {6, 14,9,0,12,8,4, -1}};
	 */
	static int 	cube_patch[64] = { 9, 8, 8, 7, 8, 6, 7, 4, 8, 7,
			6,10, 7, 4,10, 2, 8, 7, 7, 5,
			7,10,11, 3, 7,11, 4, 3, 5, 3,
			3, 1, 8, 7, 7,11, 7,10, 5, 3,
			7, 5, 4, 3,11, 3, 3, 1, 6, 4,
			10, 3, 4, 2, 3, 1,10, 3, 2, 1,
			3, 1, 1, 0};
	static int 	cube_rotation[64] = { 0, 0, 6, 0,15, 0,23,22, 2, 7,
			6,15,15, 6, 0, 1, 8,21,17,15,
			8, 1, 7,22, 2,23, 1, 1, 0, 6,
			16, 1,16,16, 6,15, 1,21, 7, 2,
			22,17, 8, 8, 0,17,21, 8,21,15,
			22,15, 0, 6, 7,22,17,23, 0,15,
			0, 6, 0, 0};


	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 9;	/* isosurf_J */
	} else {
		amb_index = 0;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);  /* bottom */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);  /* top    */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);  /* left   */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);  /* back   */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);  /* right  */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);  /* front  */
	}
	get_hex_rotation_info(cube_rotation[amb_index], &rot_info);

	switch(cube_patch[amb_index]) {
	case  0: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_A;

	break;
	case  1: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_B;

	break;
	case  2:
		cube_polys->verts_num = 16;
		for (i = 0; i < 16; i++) {
			cube_polys->verts[i] = point_convert(verts[i],&rot_info);
		}
		*cube_polys->isosurf = isosurf_C;

		break;
	case  3: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_D;

	break;
	case  4: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_E;

	break;
	case  5: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_F;

	break;
	case  6:
		cube_polys->verts_num = 16;
		for (i = 0; i < 16; i++) {
			cube_polys->verts[i] = point_convert(verts[i],&rot_info);
		}
		*cube_polys->isosurf = isosurf_G;

		break;
	case  7: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(102,&rot_info);
	*cube_polys->isosurf = isosurf_H;

	break;
	case  8: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_I;

	break;
	case  9: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_J;

	break;
	case 10: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info); }
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_K;

	break;
	case 11: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_L;

	break;
	}
}


void case13opp_tiler_hex(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, int index, Cube_polygons *cube_polys,
		int disamb_flag)
{
	int 		i;
	int 		amb_index;
	Rotation_info rot_info;
	static int 	verts[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
			201, 203, 204, 206};
	static int 	isosurf_A[13] = {2,3,10,  0,1,9,  8,4,7,  5,11,6, -1};
	static int 	isosurf_B[19] = {8,4,7,  5,11,6,  9,10,1,  1,10,2,
			9,3,10,  9,0,3, -1};
	static int 	isosurf_C[25] = {0,3,9,  9,3,10,  10,1,9,  1,10,2,
			5,11,4,  11,8,4,  11,6,8,  6,7,8, -1};
	static int 	isosurf_D[31] = {9,0,16,  5,9,16,  6,5,16,  11,6,16,
			1,11,16,  2,1,16,  10,2,16,  3,10,16,
			0,3,16,  4,7,8, -1};
	static int 	isosurf_E[37] = {9,0,16,  5,9,16,  6,5,16,  11,6,16,
			1,11,16,  2,1,16,  10,2,16,  7,10,16,
			4,7,16,  8,4,16,  3,8,16,  0,3,16, -1};
	static int 	isosurf_F[37] = {8,0,16,  7,8,16,  4,7,16,  9,4,16,
			5,9,16,  6,5,16,  11,6,16,  1,11,16,
			2,1,16,  10,2,16,  3,10,16,  0,3,16, -1};
	static int 	isosurf_G[25] = {0,1,8,  1,11,8,  11,3,8,  3,11,2,
			9,4,7,  9,7,10,  9,10,5,  5,10,6, -1};
	static int 	isosurf_H[31] = {8,0,16,  3,8,16,  2,3,16,  11,2,16,
			5,11,16,  4,5,16,  9,4,16,  1,9,16,
			0,1,16,  6,7,10, -1};
	static int 	isosurf_I[19] = {8,2,3,  8,11,2,  8,0,11,  0,1,11,
			9,4,5,  10,6,7, -1};
	static int 	isosurf_J[13] = {0,3,8,  4,5,9,  2,1,11,  10,6,7, -1};
	static int 	isosurf_K[37] = {8,0,16,  7,8,16,  4,7,16,  9,4,16,
			1,9,16,  2,1,16,  11,2,16,  5,11,16,
			6,5,16,  10,6,16,  3,10,16,  0,3,16, -1};
	static int 	isosurf_L[19] = {2,1,11,  0,3,10,  0,10,6,  9,0,6,
			9,6,5,  8,4,7, -1};
	/*
  static int 	sep_bd[6][9] = {{3, 12,1,0,  3, 13,3,2, -1},
				{3, 12,9,1,  3, 15,11,5, -1},
				{3, 14,7,4,  3, 15,5,6, -1},
				{3, 14,8,7,  3, 13,10,3, -1},
				{3, 13,2,10,  3, 15,6,11, -1},
				{3, 14,4,8,  3, 12,0,9, -1}};
  static int 	not_sep_bd[6][8] = {{6, 12,1,2,13,3,0, -1},
				    {6, 12,9,5,15,11,1, -1},
				    {6, 14,7,6,15,5,4, -1},
				    {6, 14,8,3,13,10,7, -1},
				    {6, 13,2,11,15,6,10, -1},
				    {6, 14,4,9,12,0,8, -1}};
	 */
	static int 	cube_patch[64] = { 9, 8, 8, 7, 8, 6, 7, 4, 8, 7,
			6,10, 7, 4,10, 2, 8, 7, 7, 5,
			7,10,11, 3, 7,11, 4, 3, 5, 3,
			3, 1, 8, 7, 7,11, 7,10, 5, 3,
			7, 5, 4, 3,11, 3, 3, 1, 6, 4,
			10, 3, 4, 2, 3, 1,10, 3, 2, 1,
			3, 1, 1, 0};
	static int 	cube_rotation[64] = { 0, 0, 6, 0,15, 0,23,22, 2, 7,
			6,15,15, 6, 0, 1, 8,21,17,15,
			8, 1, 7,22, 2,23, 1, 1, 0, 6,
			16, 1,16,16, 6,15, 1,21, 7, 2,
			22,17, 8, 8, 0,17,21, 8,21,15,
			22,15, 0, 6, 7,22,17,23, 0,15,
			0, 6, 0, 0};


	if (disamb_flag == 0) { /* disambiguation : none */
		amb_index = 9;
	} else {
		amb_index = 0;
		amb_index |= separation_test(mesh, cell, 11,
				fvalue, CS_type, disamb_flag);  /* bottom */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  9,
				fvalue, CS_type, disamb_flag);  /* top    */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  7,
				fvalue, CS_type, disamb_flag);  /* left   */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  4,
				fvalue, CS_type, disamb_flag);  /* back   */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  3,
				fvalue, CS_type, disamb_flag);  /* right  */
		amb_index <<= 1;
		amb_index |= separation_test(mesh, cell,  0,
				fvalue, CS_type, disamb_flag);  /* front  */
	}

	get_hex_rotation_info(cube_rotation[amb_index], &rot_info);

	switch(cube_patch[amb_index]) {
	case  0: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_A;

	break;
	case  1: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_B;

	break;
	case  2:
		cube_polys->verts_num = 16;
		for (i = 0; i < 16; i++) {
			cube_polys->verts[i] = point_convert(verts[i],&rot_info);
		}
		*cube_polys->isosurf = isosurf_C;

		break;
	case  3: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_D;

	break;
	case  4: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_E;

	break;
	case  5: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_F;

	break;
	case  6:
		cube_polys->verts_num = 16;
		for (i = 0; i < 16; i++) {
			cube_polys->verts[i] = point_convert(verts[i],&rot_info);
		}
		*cube_polys->isosurf = isosurf_G;

		break;
	case  7: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(102,&rot_info);
	*cube_polys->isosurf = isosurf_H;

	break;
	case  8: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_I;

	break;
	case  9: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_J;

	break;
	case 10: cube_polys->verts_num = 17;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	cube_polys->verts[16] = point_convert(101,&rot_info);
	*cube_polys->isosurf = isosurf_K;


	break;
	case 11: cube_polys->verts_num = 16;
	for (i = 0; i < 16; i++) {
		cube_polys->verts[i] = point_convert(verts[i],&rot_info);
	}
	*cube_polys->isosurf = isosurf_L;

	break;
	}
}


void case14_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 1, 5, 6, 8, 10, 201, 204, 205, 207};
	static int 	isosurf[13] = {1,2,0,  0,2,4,  4,2,3,  4,3,5, -1};


	switch(index) {
	case 139: get_hex_rotation_info(13, &rot_info); break;
	case 156: get_hex_rotation_info(18, &rot_info); break;
	case 178: get_hex_rotation_info( 0, &rot_info); break;
	case 198: get_hex_rotation_info(11, &rot_info); break;
	case 209: get_hex_rotation_info(10, &rot_info); break;
	case 232: get_hex_rotation_info( 7, &rot_info); break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {

	}
}


void case14opp_tiler_hex(int index, Cube_polygons *cube_polys)
{
	int 		i;
	Rotation_info rot_info;
	static int 	verts[10] = {0, 1, 5, 6, 8, 10, 200, 202, 203, 206};
	static int 	isosurf[13] = {1,0,2,  0,4,2,  4,3,2,  4,5,3, -1};


	switch(index) {
	case 116: get_hex_rotation_info(13, &rot_info); break;
	case  99: get_hex_rotation_info(18, &rot_info); break;
	case  77: get_hex_rotation_info( 0, &rot_info); break;
	case  57: get_hex_rotation_info(11, &rot_info); break;
	case  46: get_hex_rotation_info(10, &rot_info); break;
	case  23: get_hex_rotation_info( 7, &rot_info); break;
	}
	cube_polys->verts_num = 10;
	for (i = 0; i < 10; i++) {
		cube_polys->verts[i] = point_convert(verts[i], &rot_info);
	}

	*cube_polys->isosurf = isosurf;
	for (i = 0; i < 6; i++) {

	}
}

int choice_disambiguation(int disamb_flag, double fvalue, int CS_type,
		double voxel0, double voxel1, double voxel2,
		double voxel3)
{
	int 		sep;
	double		field;

	switch(disamb_flag) {
	case 1: /* asymptotic decider */
	field = calc_cross_field(voxel0, voxel1, voxel2, voxel3);
	break;
	case 2: /* facial average */
		field = facial_average(voxel0, voxel1, voxel2, voxel3);
		break;
	}

	if (!CS_type) {
		/*  when alpha volume  */
		if (voxel0 >= fvalue) { sep = (field < fvalue) ? 0 : 1; }
		else { sep = (field >= fvalue) ? 0 : 1; }
	} else {
		/*  when beta volume  */
		if (voxel0 <= fvalue) { sep = (field > fvalue) ? 0 : 1; }
		else { sep = (field <= fvalue) ? 0 : 1; }
	}

	return sep;
}
