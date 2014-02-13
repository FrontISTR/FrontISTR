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

#include "hecmw_vis_patch_const.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_calc_attr.h"
#include "hecmw_malloc.h"

int merge_vol_iso(int iso_type, Cell *cell,
		double falpha, Cube_polygons *alpha_cube, double fbeta,
		Cube_polygons *beta_cube, int bdflag, int *sum_verts,
		Point **CS_verts_tail, Point **CS_verts_refer,
		Point **CS_verts_head, Polygon **CS_polys_tail)
{
	int 		i;
	int 		*alpha_verts, *beta_verts;
	int 		*alpha_iso, *beta_iso;
	int 		alpha_vident[VERTEX_KIND], beta_vident[VERTEX_KIND];

	alpha_verts = alpha_cube->verts;
	alpha_iso = *alpha_cube->isosurf;
	for (i = 0; i < VERTEX_KIND; i++) { alpha_vident[i] = 0; }
	beta_verts = beta_cube->verts;
	beta_iso = *beta_cube->isosurf;
	for (i = 0; i < VERTEX_KIND; i++) { beta_vident[i] = 0; }

	if (!iso_type) {
		/*  make isosurface , iso_type = 0  */
		if (alpha_iso[0] != -1) {
			for (i = 0; i < alpha_cube->verts_num; i++) {
				if (alpha_verts[i] < 200) {
					if (!(alpha_vident[i]
					                   = get_vert_ident(alpha_verts[i], cell, falpha, sum_verts,
					                		   CS_verts_tail, CS_verts_refer,
					                		   CS_verts_head, bdflag))) {
						fprintf(stderr, "Error: Cannot get vertex index.\n");

						return 0;
					}
				}
			}
			i = 0;
			do {
				(*CS_polys_tail)->type = 0;
				(*CS_polys_tail)->plist = (int *) HECMW_calloc(4, sizeof(int));
				(*CS_polys_tail)->plist[0] = 3;
				(*CS_polys_tail)->plist[1] = alpha_vident[alpha_iso[i++]];
				(*CS_polys_tail)->plist[2] = alpha_vident[alpha_iso[i++]];
				(*CS_polys_tail)->plist[3] = alpha_vident[alpha_iso[i++]];
				(*CS_polys_tail)->bdflag = bdflag;
				(*CS_polys_tail)->elem_id[0] = cell->elem_id[0];
				(*CS_polys_tail)->elem_id[1] = cell->elem_id[1];

				if ((*CS_polys_tail)->nextpolygon == NULL) {
					if (((*CS_polys_tail)->nextpolygon =
						alloc_polygons(POLYGON_PACK)) == NULL) {
						fprintf(stderr, "Cannot allocate memory.\n");
						return 0;
					}
				}
				*CS_polys_tail = (*CS_polys_tail)->nextpolygon;
			} while (alpha_iso[i] != -1);
		}
#ifdef BCUBE
		if (beta_iso[0] != -1) {
			for (i = 0; i < beta_cube->verts_num; i++) {
				if (beta_verts[i] < 100) {
					beta_vident[i] = alpha_vident[i];
				}
			}
			i = 0;
			do {
				(*CS_polys_tail)->type = 1;
				(*CS_polys_tail)->plist = (int *) HECMW_calloc(4, sizeof(int));
				(*CS_polys_tail)->plist[0] = 3;
				(*CS_polys_tail)->plist[1] = beta_vident[beta_iso[i++]];
				(*CS_polys_tail)->plist[2] = beta_vident[beta_iso[i++]];
				(*CS_polys_tail)->plist[3] = beta_vident[beta_iso[i++]];
				(*CS_polys_tail)->bdflag = bdflag;
				(*CS_polys_tail)->elem_id[0] = cell->elem_id[0];
				(*CS_polys_tail)->elem_id[1] = cell->elem_id[1];
				if ((*CS_polys_tail)->nextpolygon == NULL) {
					if (((*CS_polys_tail)->nextpolygon =
						alloc_polygons(POLYGON_PACK)) == NULL) {
						return 0;
					}
				}
				*CS_polys_tail = (*CS_polys_tail)->nextpolygon;
			} while (beta_iso[i] != -1);
		}
#endif
		return 1;
	}
	return 1;    /* *** 2007/12/26 S. Ito *** */
}


int get_edge_index(int bound_index, int vert_index)
{
	int edge_index;

	switch(bound_index) {
	case 0:
		switch(vert_index) {
		case 200:
			edge_index = 0;
			break;
		case 201:
			edge_index = 1;
			break;
		case 202:
			edge_index = 2;
			break;
		case 203:
			edge_index = 3;
			break;
		}
		break;
	case 1:
		switch(vert_index) {
		case 201:
			edge_index = 9;
			break;
		case 205:
			edge_index = 5;
			break;
		case 206:
			edge_index = 11;
			break;
		case 202:
			edge_index = 1;
			break;
		}
		break;
	case 2:
		switch(vert_index) {
		case 205:
			edge_index = 4;
			break;
		case 204:
			edge_index = 7;
			break;
		case 207:
			edge_index = 6;
			break;
		case 206:
			edge_index = 5;
			break;
		}
		break;
	case 3:
		switch(vert_index) {
		case 204:
			edge_index = 8;
			break;
		case 200:
			edge_index = 3;
			break;
		case 203:
			edge_index = 10;
			break;
		case 207:
			edge_index = 7;
			break;
		}
		break;
	case 4:
		switch(vert_index) {
		case 203:
			edge_index = 2;
			break;
		case 202:
			edge_index = 11;
			break;
		case 206:
			edge_index = 6;
			break;
		case 207:
			edge_index = 10;
			break;
		}
		break;
	case 5:
		switch(vert_index) {
		case 204:
			edge_index = 4;
			break;
		case 205:
			edge_index = 9;
			break;
		case 201:
			edge_index = 0;
			break;
		case 200:
			edge_index = 8;
			break;
		}
		break;
	}

	return edge_index;
}

int get_vert_ident(int pindex, Cell *cell, double fvalue, int *sum_verts,
		Point **CS_verts_tail, Point **CS_verts_refer,
		Point **CS_verts_head, int bdflag)
{
	int 		table_no;
	int 		pident;
	double	cdata;
	Fgeom		pgeom;

	get_point_geom(pindex, cell, fvalue, &pgeom, &cdata, 1);
	table_no = 0;

	CS_verts_refer[table_no] = CS_verts_head[table_no];

	if (pindex < 100) {
		pident = search_verts_table(&pgeom, CS_verts_refer, bdflag, table_no);
		if (pident) { return pident; }
	} else if (pindex >= 200) {
		fvalue = cell->s_data[pindex-200];
		cdata = cell->c_data[pindex-200];
		pident = search_verts_table(&pgeom, CS_verts_refer, bdflag, table_no);
		if (pident) { return pident; }
	}
	pident = ++(*sum_verts);

	add_verts_table(CS_verts_tail, table_no, pident, fvalue, cdata,
			&pgeom, bdflag);

	return pident;
}


int search_verts_table(Fgeom *pgeom, Point **CS_verts_refer, int bdflag,
		int table_no)
{
	Point 	*CS_verts;

	CS_verts = CS_verts_refer[table_no];

	while (CS_verts->ident) {
		if ((fabs(CS_verts->geom.x - pgeom->x) < EPSILON) &&
				(fabs(CS_verts->geom.y - pgeom->y) < EPSILON) &&
				(fabs(CS_verts->geom.z - pgeom->z) < EPSILON)) {
			if ((CS_verts->info == 1) || (CS_verts->info == 0)) {
				if (bdflag < 0) {
					CS_verts->info = -1;
				} else if (bdflag >= 1024) {
					CS_verts->info = 2;
				} else if (bdflag != HEX_FACE_INDEX) {
					CS_verts->info = 1;
				} else {
					CS_verts->info = 0;
				}
			} else if (CS_verts->info == 2) {
				if (bdflag < 0) CS_verts->info = -1;
			}
			return CS_verts->ident; }
		CS_verts = CS_verts->nextpoint;
	}
	return 0;
}

int add_verts_table(Point **CS_verts_tail, int table_no, int pident,
		double pfield, double cdata, Fgeom *pgeom, int bdflag)
{
	CS_verts_tail[table_no]->ident = pident;
	CS_verts_tail[table_no]->field = pfield;
	CS_verts_tail[table_no]->cdata = cdata;
	CS_verts_tail[table_no]->geom = *pgeom;

	if (bdflag < 0) {
		CS_verts_tail[table_no]->info = -1;
	} else if (bdflag >= 1024) {
		CS_verts_tail[table_no]->info = 2;
	} else if (bdflag != HEX_FACE_INDEX) {
		CS_verts_tail[table_no]->info = 1;
	} else {
		CS_verts_tail[table_no]->info = 0;
	}

	if (CS_verts_tail[table_no]->nextpoint == NULL) {
		if ((CS_verts_tail[table_no]->nextpoint =
			alloc_verts(VERTEX_PACK)) == NULL) {
			fprintf(stderr, "Cannot allocate memory.\n");
			return 0; }
	}
	CS_verts_tail[table_no] = CS_verts_tail[table_no]->nextpoint;

	return pident;
}

int check_vertex(Cell *cell, double fvalue, Cube_polygons *cube_polys)
{
	int 	i;
	int 	flag;
	int 	verts_num;
	int 	*verts;

	flag = 0;
	verts_num = cube_polys->verts_num;
	verts = cube_polys->verts;

	if (cell->s_data[0] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 0)||(verts[i] == 3)
					||(verts[i] == 8)||(verts[i] == 100)) {
				verts[i] = 200;
				flag++;
			}
		}
	}
	if (cell->s_data[1] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 0)||(verts[i] == 1)
					||(verts[i] == 9)||(verts[i] == 101)) {
				verts[i] = 201;
				flag++;
			}
		}
	}
	if (cell->s_data[2] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 1)||(verts[i] == 2)
					||(verts[i] == 11)||(verts[i] == 102)) {
				verts[i] = 202;
				flag++;
			}
		}
	}
	if (cell->s_data[3] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 2)||(verts[i] == 3)
					||(verts[i] == 10)||(verts[i] == 103)) {
				verts[i] = 203;
				flag++;
			}
		}
	}
	if (cell->s_data[4] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 4)||(verts[i] == 7)
					||(verts[i] == 8)||(verts[i] == 102)) {
				verts[i] = 204;
				flag++;
			}
		}
	}
	if (cell->s_data[5] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 4)||(verts[i] == 5)
					||(verts[i] == 9)||(verts[i] == 103)) {
				verts[i] = 205;
				flag++;
			}
		}
	}
	if (cell->s_data[6] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 5)||(verts[i] == 6)
					||(verts[i] == 11)||(verts[i] == 100)) {
				verts[i] = 206;
				flag++;
			}
		}
	}
	if (cell->s_data[7] == fvalue) {
		for (i = 0; i < verts_num; i++) {
			if ((verts[i] == 6)||(verts[i] == 7)
					||(verts[i] == 10)||(verts[i] == 101)) {
				verts[i] = 207;
				flag++;
			}
		}
	}

	return flag;
}

