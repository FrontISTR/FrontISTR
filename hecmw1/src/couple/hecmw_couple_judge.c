/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#include "hecmw_struct.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_judge.h"


#define FRAC_1_2 (0.5)


#define FRAC_1_3 (0.33333333333333333)


#define FRAC_1_4 (0.25)

#define EPS_PRODUCT (1.0E-06)
#define EPS_DISTANCE (1.0E-02)

/*================================================================================================*/

static void
calc_gravity_tri(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z,
		double p3_x, double p3_y, double p3_z, double *g_x, double *g_y, double *g_z)
{
	*g_x = (p1_x + p2_x + p3_x) * FRAC_1_3;
	*g_y = (p1_y + p2_y + p3_y) * FRAC_1_3;
	*g_z = (p1_z + p2_z + p3_z) * FRAC_1_3;
}



static void
calc_gravity_quad(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z,
		double p3_x, double p3_y, double p3_z, double p4_x, double p4_y, double p4_z,
		double *g_x, double *g_y, double *g_z)
{
	*g_x = (p1_x + p2_x + p3_x + p4_x) * FRAC_1_4;
	*g_y = (p1_y + p2_y + p3_y + p4_y) * FRAC_1_4;
	*g_z = (p1_z + p2_z + p3_z + p4_z) * FRAC_1_4;
}



static void
calc_normal_vector_tri(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z,
		double p3_x, double p3_y, double p3_z, double *vn_x, double *vn_y, double *vn_z)
{
	double v1_x, v1_y, v1_z;
	double v2_x, v2_y, v2_z;

	v1_x = p2_x - p1_x;
	v1_y = p2_y - p1_y;
	v1_z = p2_z - p1_z;

	v2_x = p3_x - p1_x;
	v2_y = p3_y - p1_y;
	v2_z = p3_z - p1_z;

	*vn_x = v1_y * v2_z - v1_z * v2_y;
	*vn_y = v1_z * v2_x - v1_x * v2_z;
	*vn_z = v1_x * v2_y - v1_y * v2_x;
}



static void
calc_normal_vector_quad(double p1_x, double p1_y, double p1_z,
		double p2_x, double p2_y, double p2_z, double p3_x, double p3_y, double p3_z,
		double p4_x, double p4_y, double p4_z, double *vn_x, double *vn_y, double *vn_z)
{
	double vn1_x, vn1_y, vn1_z;
	double vn2_x, vn2_y, vn2_z;

	calc_normal_vector_tri(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p4_x, p4_y, p4_z,
			&vn1_x, &vn1_y, &vn1_z);
	calc_normal_vector_tri(p3_x, p3_y, p3_z, p4_x, p4_y, p4_z, p2_x, p2_y, p2_z,
			&vn2_x, &vn2_y, &vn2_z);

	*vn_x = (vn1_x + vn2_x) * FRAC_1_2;
	*vn_y = (vn1_y + vn2_y) * FRAC_1_2;
	*vn_z = (vn1_z + vn2_z) * FRAC_1_2;
}



static double
calc_dot_product(double v1_x, double v1_y, double v1_z, double v2_x, double v2_y, double v2_z,
		double *dot_product)
{
	*dot_product = v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;

	return *dot_product;
}



static int
calc_dot_product_tri(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z,
		double p3_x, double p3_y, double p3_z, double p_x, double p_y, double p_z,
		double *dot_product, double *distance)
{
	double vn_x, vn_y, vn_z;
	double g_x, g_y, g_z;
	double vp_x, vp_y, vp_z;

	calc_normal_vector_tri(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z,
			&vn_x, &vn_y, &vn_z);
	calc_gravity_tri(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, &g_x, &g_y, &g_z);

	vp_x = p_x - g_x;
	vp_y = p_y - g_y;
	vp_z = p_z - g_z;

	calc_dot_product(vn_x, vn_y, vn_z, vp_x, vp_y, vp_z, dot_product);

	*distance = fabs(*dot_product / sqrt(vn_x * vn_x + vn_y * vn_y + vn_z * vn_z));

	return *dot_product;
}



static int
calc_dot_product_quad(double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z,
		double p3_x, double p3_y, double p3_z, double p4_x, double p4_y, double p4_z,
		double p_x, double p_y, double p_z, double *dot_product, double *distance)
{
	double vn_x, vn_y, vn_z;
	double g_x, g_y, g_z;
	double vp_x, vp_y, vp_z;

	calc_normal_vector_quad(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z,
			&vn_x, &vn_y, &vn_z);

	calc_gravity_quad(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z, p4_x, p4_y, p4_z,
			&g_x, &g_y, &g_z);

	vp_x = p_x - g_x;
	vp_y = p_y - g_y;
	vp_z = p_z - g_z;

	calc_dot_product(vn_x, vn_y, vn_z, vp_x, vp_y, vp_z, dot_product);

	*distance = fabs(*dot_product / sqrt(vn_x * vn_x + vn_y * vn_y + vn_z * vn_z));

	return *dot_product;
}



extern int
HECMW_couple_judge_tet1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance)
{
	double coord_x[4], coord_y[4], coord_z[4], _dot_product[4], _distance[4];
	int node_index, node_id, n_positive;
	int i;

	node_index = local_mesh->elem_node_index[elem-1];
	for(i=0; i<4; i++) {
		node_id    = local_mesh->elem_node_item[node_index+i];
		coord_x[i] = local_mesh->node[3*(node_id-1)  ];
		coord_y[i] = local_mesh->node[3*(node_id-1)+1];
		coord_z[i] = local_mesh->node[3*(node_id-1)+2];
	}

	calc_dot_product_tri(coord_x[1], coord_y[1], coord_z[1], coord_x[2], coord_y[2], coord_z[2],
			coord_x[3], coord_y[3], coord_z[3], coord_px, coord_py, coord_pz,
			&_dot_product[0], &_distance[0]);
	calc_dot_product_tri(coord_x[0], coord_y[0], coord_z[0], coord_x[3], coord_y[3], coord_z[3],
			coord_x[2], coord_y[2], coord_z[2], coord_px, coord_py, coord_pz,
			&_dot_product[1], &_distance[1]);
	calc_dot_product_tri(coord_x[0], coord_y[0], coord_z[0], coord_x[1], coord_y[1], coord_z[1],
			coord_x[3], coord_y[3], coord_z[3], coord_px, coord_py, coord_pz,
			&_dot_product[2], &_distance[2]);
	calc_dot_product_tri(coord_x[0], coord_y[0], coord_z[0], coord_x[2], coord_y[2], coord_z[2],
			coord_x[1], coord_y[1], coord_z[1], coord_px, coord_py, coord_pz,
			&_dot_product[3], &_distance[3]);

	*dot_product = _dot_product[surf_id-1];
	*distance    = _distance[surf_id-1];

	n_positive = 0;
	for(i=0; i<4; i++) {
		if(_dot_product[i] > EPS_PRODUCT) {
			if(_distance[i] > EPS_DISTANCE) {
				return -1;
			}
			n_positive++;
		}
	}

	return n_positive;
}


extern int
HECMW_couple_judge_hex1(const struct hecmwST_local_mesh *local_mesh, int elem, int surf_id,
		double coord_px, double coord_py, double coord_pz, double *dot_product, double *distance)
{
	double coord_x[8], coord_y[8], coord_z[8], _dot_product[8], _distance[8];
	int node_index, node_id, n_positive, i;

	node_index = local_mesh->elem_node_index[elem-1];
	for(i=0; i<8; i++) {
		node_id    = local_mesh->elem_node_item[node_index+i];
		coord_x[i] = local_mesh->node[3*(node_id-1)  ];
		coord_y[i] = local_mesh->node[3*(node_id-1)+1];
		coord_z[i] = local_mesh->node[3*(node_id-1)+2];
	}

	calc_dot_product_quad(coord_x[3], coord_y[3], coord_z[3], coord_x[0], coord_y[0], coord_z[0],
			coord_x[4], coord_y[4], coord_z[4], coord_x[7], coord_y[7], coord_z[7],
			coord_px, coord_py, coord_pz, &_dot_product[0], &_distance[0]);
	calc_dot_product_quad(coord_x[1], coord_y[1], coord_z[1], coord_x[2], coord_y[2], coord_z[2],
			coord_x[6], coord_y[6], coord_z[6], coord_x[5], coord_y[5], coord_z[5],
			coord_px, coord_py, coord_pz, &_dot_product[1], &_distance[1]);
	calc_dot_product_quad(coord_x[0], coord_y[0], coord_z[0], coord_x[1], coord_y[1], coord_z[1],
			coord_x[5], coord_y[5], coord_z[5], coord_x[4], coord_y[4], coord_z[4],
			coord_px, coord_py, coord_pz, &_dot_product[2], &_distance[2]);
	calc_dot_product_quad(coord_x[2], coord_y[2], coord_z[2], coord_x[3], coord_y[3], coord_z[3],
			coord_x[7], coord_y[7], coord_z[7], coord_x[6], coord_y[6], coord_z[6],
			coord_px, coord_py, coord_pz, &_dot_product[3], &_distance[3]);
	calc_dot_product_quad(coord_x[3], coord_y[3], coord_z[3], coord_x[2], coord_y[2], coord_z[2],
			coord_x[1], coord_y[1], coord_z[1], coord_x[0], coord_y[0], coord_z[0],
			coord_px, coord_py, coord_pz, &_dot_product[4], &_distance[4]);
	calc_dot_product_quad(coord_x[4], coord_y[4], coord_z[4], coord_x[5], coord_y[5], coord_z[5],
			coord_x[6], coord_y[6], coord_z[6], coord_x[7], coord_y[7], coord_z[7],
			coord_px, coord_py, coord_pz, &_dot_product[5], &_distance[5]);

	*dot_product = _dot_product[surf_id-1];
	*distance    = _distance[surf_id-1];

	n_positive = 0;
	for(i=0; i<6; i++) {
		if(_dot_product[i] > EPS_PRODUCT) {
			if(_distance[i] > EPS_DISTANCE) {
				return -1;
			}
			n_positive++;
		}
	}

	return n_positive;
}
