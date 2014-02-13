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

#include "hecmw_vis_calc_attr.h"

#include <math.h>

/*----------------------------------------------------------------------
#     Subroutines in this file on isosurface generation by Marching Cubes is based
	  on the revision of Dr. Yuriko Takeshima's codes when she was working part time in RIST
#---------------------------------------------------------------------- */


/*  calculate the geometry at the point in cube  */
int get_point_geom(int point_index, Cell *cell, double fvalue,
		Fgeom *point_geom, double *cdata, int disamb_flag)
{
	if (point_index < 100) {
		get_edgepoint(point_index, cell, fvalue, point_geom, cdata);
	} else if (point_index < 200) {
		get_insidepoint((point_index - 100), cell, fvalue, point_geom, cdata,
				disamb_flag);
	} else {
		get_gridpoint((point_index - 200), cell, point_geom, cdata);
	}

	return 1;
}

/*  return the geometry at the vertex of cube  */
void get_gridpoint(int voxel_index, Cell *cell, Fgeom *vert_geom,
		double *cdata)
{
	vert_geom->x = cell->axis[voxel_index*3];
	vert_geom->y = cell->axis[voxel_index*3+1];
	vert_geom->z = cell->axis[voxel_index*3+2];
	*cdata = cell->c_data[voxel_index];

}

void get_edgepoint(int edge_index, Cell *cell, double fvalue,
		Fgeom *vert_geom, double *cdata)
{
	double		flip_ratio, ratio;

	switch(edge_index) {
	case  0:
		ratio = linear_interpolate(cell->s_data[0], cell->s_data[1], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[0] + ratio * cell->axis[3];
		vert_geom->y = flip_ratio * cell->axis[1] + ratio * cell->axis[4];
		vert_geom->z = flip_ratio * cell->axis[2] + ratio * cell->axis[5];
		*cdata = flip_ratio * cell->c_data[0] + ratio * cell->c_data[1];
		break;
	case  1:
		ratio = linear_interpolate(cell->s_data[1], cell->s_data[2], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[3] + ratio * cell->axis[6];
		vert_geom->y = flip_ratio * cell->axis[4] + ratio * cell->axis[7];
		vert_geom->z = flip_ratio * cell->axis[5] + ratio * cell->axis[8];
		*cdata = flip_ratio * cell->c_data[1] + ratio * cell->c_data[2];
		break;
	case  2:
		ratio = linear_interpolate(cell->s_data[3], cell->s_data[2], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[9]  + ratio * cell->axis[6];
		vert_geom->y = flip_ratio * cell->axis[10] + ratio * cell->axis[7];
		vert_geom->z = flip_ratio * cell->axis[11] + ratio * cell->axis[8];
		*cdata = flip_ratio * cell->c_data[3] + ratio * cell->c_data[2];
		break;
	case  3:
		ratio = linear_interpolate(cell->s_data[0], cell->s_data[3], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[0] + ratio * cell->axis[9];
		vert_geom->y = flip_ratio * cell->axis[1] + ratio * cell->axis[10];
		vert_geom->z = flip_ratio * cell->axis[2] + ratio * cell->axis[11];
		*cdata = flip_ratio * cell->c_data[0] + ratio * cell->c_data[3];
		break;
	case  4:
		ratio = linear_interpolate(cell->s_data[4], cell->s_data[5], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[12] + ratio * cell->axis[15];
		vert_geom->y = flip_ratio * cell->axis[13] + ratio * cell->axis[16];
		vert_geom->z = flip_ratio * cell->axis[14] + ratio * cell->axis[17];
		*cdata = flip_ratio * cell->c_data[4] + ratio * cell->c_data[5];
		break;
	case  5:
		ratio = linear_interpolate(cell->s_data[5], cell->s_data[6], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[15] + ratio * cell->axis[18];
		vert_geom->y = flip_ratio * cell->axis[16] + ratio * cell->axis[19];
		vert_geom->z = flip_ratio * cell->axis[17] + ratio * cell->axis[20];
		*cdata = flip_ratio * cell->c_data[5] + ratio * cell->c_data[6];
		break;
	case  6:
		ratio = linear_interpolate(cell->s_data[7], cell->s_data[6], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[21] + ratio * cell->axis[18];
		vert_geom->y = flip_ratio * cell->axis[22] + ratio * cell->axis[19];
		vert_geom->z = flip_ratio * cell->axis[23] + ratio * cell->axis[20];
		*cdata = flip_ratio * cell->c_data[7] + ratio * cell->c_data[6];
		break;
	case  7:
		ratio = linear_interpolate(cell->s_data[4], cell->s_data[7], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[12] + ratio * cell->axis[21];
		vert_geom->y = flip_ratio * cell->axis[13] + ratio * cell->axis[22];
		vert_geom->z = flip_ratio * cell->axis[14] + ratio * cell->axis[23];
		*cdata = flip_ratio * cell->c_data[4] + ratio * cell->c_data[7];
		break;
	case  8:
		ratio = linear_interpolate(cell->s_data[0], cell->s_data[4], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[0] + ratio * cell->axis[12];
		vert_geom->y = flip_ratio * cell->axis[1] + ratio * cell->axis[13];
		vert_geom->z = flip_ratio * cell->axis[2] + ratio * cell->axis[14];
		*cdata = flip_ratio * cell->c_data[0] + ratio * cell->c_data[4];
		break;
	case  9:
		ratio = linear_interpolate(cell->s_data[1], cell->s_data[5], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[3] + ratio * cell->axis[15];
		vert_geom->y = flip_ratio * cell->axis[4] + ratio * cell->axis[16];
		vert_geom->z = flip_ratio * cell->axis[5] + ratio * cell->axis[17];
		*cdata = flip_ratio * cell->c_data[1] + ratio * cell->c_data[5];
		break;
	case 10:
		ratio = linear_interpolate(cell->s_data[3], cell->s_data[7], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[9]  + ratio * cell->axis[21];
		vert_geom->y = flip_ratio * cell->axis[10] + ratio * cell->axis[22];
		vert_geom->z = flip_ratio * cell->axis[11] + ratio * cell->axis[23];
		*cdata = flip_ratio * cell->c_data[3] + ratio * cell->c_data[7];
		break;
	case 11:
		ratio = linear_interpolate(cell->s_data[2], cell->s_data[6], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[6] + ratio * cell->axis[18];
		vert_geom->y = flip_ratio * cell->axis[7] + ratio * cell->axis[19];
		vert_geom->z = flip_ratio * cell->axis[8] + ratio * cell->axis[20];
		*cdata = flip_ratio * cell->c_data[2] + ratio * cell->c_data[6];
		break;
	}
}

/* return the geometry at the interpolate point in cube  */
void get_insidepoint(int inside_index, Cell *cell, double fvalue,
		Fgeom *vert_geom, double *cdata, int disamb_flag)
{
	double 	ratio, flip_ratio;

	switch(inside_index) {
	case 0:
		ratio = linear_interpolate(cell->s_data[0], cell->s_data[6], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[0] + ratio * cell->axis[18];
		vert_geom->y = flip_ratio * cell->axis[1] + ratio * cell->axis[19];
		vert_geom->z = flip_ratio * cell->axis[2] + ratio * cell->axis[20];
		*cdata = flip_ratio * cell->c_data[0] + ratio * cell->c_data[6];
		break;
	case 1:
		ratio = linear_interpolate(cell->s_data[1], cell->s_data[7], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[3] + ratio * cell->axis[21];
		vert_geom->y = flip_ratio * cell->axis[4] + ratio * cell->axis[22];
		vert_geom->z = flip_ratio * cell->axis[5] + ratio * cell->axis[23];
		*cdata = flip_ratio * cell->c_data[1] + ratio * cell->c_data[7];
		break;
	case 2:
		ratio = linear_interpolate(cell->s_data[2], cell->s_data[4], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[6] + ratio * cell->axis[12];
		vert_geom->y = flip_ratio * cell->axis[7] + ratio * cell->axis[13];
		vert_geom->z = flip_ratio * cell->axis[8] + ratio * cell->axis[14];
		*cdata = flip_ratio * cell->c_data[2] + ratio * cell->c_data[4];
		break;
	case 3:
		ratio = linear_interpolate(cell->s_data[3], cell->s_data[5], fvalue);
		flip_ratio = 1 - ratio;
		vert_geom->x = flip_ratio * cell->axis[9]  + ratio * cell->axis[15];
		vert_geom->y = flip_ratio * cell->axis[10] + ratio * cell->axis[16];
		vert_geom->z = flip_ratio * cell->axis[11] + ratio * cell->axis[17];
		*cdata = flip_ratio * cell->c_data[3] + ratio * cell->c_data[5];
		break;
	case 4: case 5: case 6:
		vert_geom->x = (cell->axis[0] + cell->axis[3] + cell->axis[6] +
				cell->axis[9] + cell->axis[12] + cell->axis[15] +
				cell->axis[18] + cell->axis[21])/8.0;
		vert_geom->y = (cell->axis[1] + cell->axis[4] + cell->axis[7] +
				cell->axis[10] + cell->axis[13] + cell->axis[16] +
				cell->axis[19] + cell->axis[22])/8.0;
		vert_geom->z = (cell->axis[2] + cell->axis[5] + cell->axis[8] +
				cell->axis[11] + cell->axis[14] + cell->axis[17] +
				cell->axis[20] + cell->axis[23])/8.0;
		*cdata = (cell->c_data[0] + cell->c_data[1] + cell->c_data[2] +
				cell->c_data[3] + cell->c_data[4] + cell->c_data[5] +
				cell->c_data[6] + cell->c_data[7])/8.0;
		break;
	}
}

/*  calculate linear interpolate value  */
double linear_interpolate(double left, double right, double fvalue)
{
	double ratio;

	if(fabs(left-right)<EPSILON) ratio=0.0;
	else ratio=((fvalue - left)/(right - left));
	return ratio;
}


/*  calculate the field value at the intersection point of the asymptotes  */
double calc_cross_field(double f00, double f10, double f11, double f01)
{
	double value;
	if(fabs(f00+f11-f01-f10)<EPSILON)
		value=(f00+f11+f01+f10)/4.0;
	else
		value=((f00*f11 - f10*f01)/(f00 + f11 - f01 - f10));
	return value;
}

double facial_average(double f00, double f10, double f11, double f01)
{
	return ((f00+f11+ f10+f01)/4);
}

