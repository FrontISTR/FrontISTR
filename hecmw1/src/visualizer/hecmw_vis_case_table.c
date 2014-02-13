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

#include "hecmw_vis_case_table.h"

#include <math.h>
#include "hecmw_vis_intersection_find.h"


int p2h(int i) {
	if(i==0) return 0;
	else if(i==1) return 1;
	else if(i==2) return 2;
	else if(i==3) return 2;
	else if(i==4) return 3;
	else if(i==5) return 4;
	else if(i==6) return 5;
	else if(i==7) return 5;
	return 0;
}

int get_data(Surface *sff,struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int elemID, Cell *cell, int tn_component)
{
	int		i, j;
	int		nodeID;
	int		c_base, s_base;
	double    c[13], s[13];
	int       tmp_int;
	double     tmp_d;

	c_base = 0;

	if(sff->display_way!=4) {
		for (i = 0; i < sff->color_comp; i++) {
			c_base += data->nn_dof[i];
		}
	}
	if(sff->surface_style==2) {
		s_base=0;
		for(i=0;i<sff->data_comp;i++)
			s_base+=data->nn_dof[i];
	}

	/*  set field data of voxel in cube  */
	if((mesh->elem_type[elemID]==361) || (mesh->elem_type[elemID]==362)){
		for (i = 0; i < HEX_N_NODE; i++) {
			nodeID = mesh->elem_node_item[mesh->elem_node_index[elemID]+i];
			if(sff->display_way!=4) {
				if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
					tmp_d=0.0;
					for(j=0;j<data->nn_dof[sff->color_comp];j++) {
						c[j]=data->node_val_item[(nodeID-1)*tn_component+c_base+j];
						tmp_d+=c[j]*c[j];
					}
					cell->c_data[i]=sqrt(tmp_d);
				}
				else if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp!=0))
					cell->c_data[i]=data->node_val_item[(nodeID-1)*tn_component+c_base+(sff->color_subcomp-1)];
				else if(data->nn_dof[sff->color_comp]==1)
					cell->c_data[i] = data->node_val_item[(nodeID-1)*tn_component+c_base];

			}
			else if(sff->display_way==4)
				cell->c_data[i]=sff->specified_color;

			cell->axis[i*3] = mesh->node[(nodeID-1)*3];
			cell->axis[i*3+1] = mesh->node[(nodeID-1)*3+1];
			cell->axis[i*3+2] = mesh->node[(nodeID-1)*3+2];

			if(sff->surface_style==2) {

				if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp==0)) {
					tmp_d=0.0;
					for(j=0;j<data->nn_dof[sff->data_comp];j++) {
						s[j]=data->node_val_item[(nodeID-1)*tn_component+s_base+j];
						tmp_d+=s[j]*s[j];
					}
					cell->s_data[i]=sqrt(tmp_d);
				}
				else if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp!=0))
					cell->s_data[i] = data->node_val_item[(nodeID-1)*tn_component+s_base+(sff->data_subcomp-1)];
				else if(data->nn_dof[sff->data_comp]==1)
					cell->s_data[i]=data->node_val_item[(nodeID-1)*tn_component+s_base];
			}
			else if(sff->surface_style==3)

				cell->s_data[i]=get_value_equ(sff->cont_equ,sff->cross_type,cell->axis[i*3],
						cell->axis[i*3+1], cell->axis[i*3+2]);
		}
	}
	else if((mesh->elem_type[elemID]==351) || (mesh->elem_type[elemID]==352)) {
		for (i = 0; i < 8; i++) {
			tmp_int=p2h(i);
			nodeID = mesh->elem_node_item[mesh->elem_node_index[elemID]+tmp_int];
			if(sff->display_way!=4) {
				if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
					tmp_d=0.0;
					for(j=0;j<data->nn_dof[sff->color_comp];j++) {
						c[j]=data->node_val_item[(nodeID-1)*tn_component+c_base+j];
						tmp_d+=c[j]*c[j];
					}
					cell->c_data[i]=sqrt(tmp_d);
				}
				else if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp!=0))
					cell->c_data[i]=data->node_val_item[(nodeID-1)*tn_component+c_base+(sff->color_subcomp-1)];
				else if(data->nn_dof[sff->color_comp]==1)
					cell->c_data[i] = data->node_val_item[(nodeID-1)*tn_component+c_base];

			}
			else if(sff->display_way==4)
				cell->c_data[i]=sff->specified_color;

			cell->axis[i*3] = mesh->node[(nodeID-1)*3];
			cell->axis[i*3+1] = mesh->node[(nodeID-1)*3+1];
			cell->axis[i*3+2] = mesh->node[(nodeID-1)*3+2];

			if(sff->surface_style==2) {

				if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp==0)) {
					tmp_d=0.0;
					for(j=0;j<data->nn_dof[sff->data_comp];j++) {
						s[j]=data->node_val_item[(nodeID-1)*tn_component+s_base+j];
						tmp_d+=s[j]*s[j];
					}
					cell->s_data[i]=sqrt(tmp_d);
				}
				else if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp!=0))
					cell->s_data[i] = data->node_val_item[(nodeID-1)*tn_component+s_base+(sff->data_subcomp-1)];
				else if(data->nn_dof[sff->data_comp]==1)
					cell->s_data[i]=data->node_val_item[(nodeID-1)*tn_component+s_base];
			}
			else if(sff->surface_style==3)

				cell->s_data[i]=get_value_equ(sff->cont_equ,sff->cross_type,cell->axis[i*3],
						cell->axis[i*3+1], cell->axis[i*3+2]);
		}




	}

	cell->elem_id[0] = mesh->elem_ID[elemID*2];
	cell->elem_id[1] = mesh->elem_ID[elemID*2+1];


	return 1;
}

double get_value_equ(double cont_equ[10],int cross_type,double x,double y,double z)

{
	double value;
	if(cross_type==1)
		value=cont_equ[0]*x+cont_equ[1]*y+cont_equ[2]*z+cont_equ[3];
	else if(cross_type==2)
		value=cont_equ[0]*x*x+cont_equ[1]*y*y+cont_equ[2]*z*z+cont_equ[3]*x*y+cont_equ[4]*y*z
		+cont_equ[5]*x*z+cont_equ[6]*x+cont_equ[7]*y+cont_equ[8]*z+cont_equ[9];
	else {
		fprintf(stderr,"The value of cross type is wrong\n");
		return(-1);
	}
	return(value);
}

/*  make the boundary patches of interval volumes in each cube  */
int make_tile(struct hecmwST_local_mesh *mesh, Cell *cell, double fvalue,
		int CS_type, Cube_polygons *cube_polys, int disamb_flag)
{
	int 		index, config;

	/*  decide the configuration of volume in cube  */
	index = 0;
	if (!CS_type) {
		/*  when alpha volume  */
		index |= (cell->s_data[0] >= fvalue) ? (int)   1 : (int) 0;
		index |= (cell->s_data[1] >= fvalue) ? (int)   2 : (int) 0;
		index |= (cell->s_data[2] >= fvalue) ? (int)   4 : (int) 0;
		index |= (cell->s_data[3] >= fvalue) ? (int)   8 : (int) 0;
		index |= (cell->s_data[4] >= fvalue) ? (int)  16 : (int) 0;
		index |= (cell->s_data[5] >= fvalue) ? (int)  32 : (int) 0;
		index |= (cell->s_data[6] >= fvalue) ? (int)  64 : (int) 0;
		index |= (cell->s_data[7] >= fvalue) ? (int) 128 : (int) 0;
	} else {
		/*  when beta volume  */
		index |= (cell->s_data[0] < fvalue) ? (int)   1 : (int) 0;
		index |= (cell->s_data[1] < fvalue) ? (int)   2 : (int) 0;
		index |= (cell->s_data[2] < fvalue) ? (int)   4 : (int) 0;
		index |= (cell->s_data[3] < fvalue) ? (int)   8 : (int) 0;
		index |= (cell->s_data[4] < fvalue) ? (int)  16 : (int) 0;
		index |= (cell->s_data[5] < fvalue) ? (int)  32 : (int) 0;
		index |= (cell->s_data[6] < fvalue) ? (int)  64 : (int) 0;
		index |= (cell->s_data[7] < fvalue) ? (int) 128 : (int) 0;
	}

	/*  make tile data of cube  */
	if (!index) return 0;
	if ((config = get_config_info(index)) > 0) {
		switch(config) {
		case   1: case1_tiler_hex(index, cube_polys);
		break;
		case   2: case2_tiler_hex(index, cube_polys);
		break;
		case   3: case3_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case   4: case4_tiler_hex(index, cube_polys);
		break;
		case   5: case5_tiler_hex(index, cube_polys);
		break;
		case   6: case6_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case   7: case7_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case   8: case8_tiler_hex(index, cube_polys);
		break;
		case   9: case9_tiler_hex(index, cube_polys);
		break;
		case  10: case10_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  11: case11_tiler_hex(index, cube_polys);
		break;
		case  12: case12_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  13: case13_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  14: case14_tiler_hex(index, cube_polys);
		break;
		}
	} else {
		switch(config) {
		case   0: case0opp_tiler_hex(cube_polys);
		break;
		case  -1: case1opp_tiler_hex(index, cube_polys);
		break;
		case  -2: case2opp_tiler_hex(index, cube_polys);
		break;
		case  -3: case3opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  -4: case4opp_tiler_hex(index, cube_polys);
		break;
		case  -5: case5opp_tiler_hex(index, cube_polys);
		break;
		case  -6: case6opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  -7: case7opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case  -8: case8opp_tiler_hex(index, cube_polys);
		break;
		case  -9: case9opp_tiler_hex(index, cube_polys);
		break;
		case -10: case10opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case -11: case11opp_tiler_hex(index,cube_polys);
		break;
		case -12: case12opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case -13: case13opp_tiler_hex(mesh, cell, fvalue, CS_type, index,
				cube_polys, disamb_flag);
		break;
		case -14: case14opp_tiler_hex(index, cube_polys);
		break;
		}
	}

	return index;
}

/*  decide the cube configuration  */
int get_config_info(int index)
{
	static int	config_list[256] = {
			0,  1,  1,  2,  1,  3,  2,  5,  1,  2,
			3,  5,  2,  5,  5,  8,  1,  2,  3,  5,
			4,  6,  6,-14,  3,  5,  7, -9,  6, 11,
			12, -5,  1,  3,  2,  5,  3,  7,  5, -9,
			4,  6,  6, 11,  6, 12,-14, -5,  2,  5,
			5,  8,  6, 12, 11, -5,  6,-14, 12, -5,
			10, -6, -6, -2,  1,  4,  3,  6,  2,  6,
			5, 11,  3,  6,  7, 12,  5,-14, -9, -5,
			3,  6,  7, 12,  6, 10, 12, -6,  7, 12,
			-13, -7, 12, -6, -7, -3,  2,  6,  5,-14,
			5, 12,  8, -5,  6, 10, 12, -6, 11, -6,
			-5, -2,  5, 11, -9, -5,-14, -6, -5, -2,
			12, -6, -7, -3, -6, -4, -3, -1,  1,  3,
			4,  6,  3,  7,  6,-12,  2,  5,  6, 14,
			5,  9,-11, -5,  2,  5,  6,-11,  6,-12,
			-10, -6,  5, -8,-12, -5, 14, -5, -6, -2,
			3,  7,  6,-12,  7, 13,-12, -7,  6,-12,
			-10, -6,-12, -7, -6, -3,  5,  9, 14, -5,
			-12, -7, -6, -3,-11, -5, -6, -2, -6, -3,
			-4, -1,  2,  6,  6,-10,  5,-12, 14, -6,
			5,-11,-12, -6, -8, -5, -5, -2,  5, 14,
			-12, -6,-11, -6, -6, -4,  9, -5, -7, -3,
			-5, -2, -3, -1,  5,-12,-11, -6,  9, -7,
			-5, -3, 14, -6, -6, -4, -5, -3, -2, -1,
			-8, -5, -5, -2, -5, -3, -2, -1, -5, -2,
			-3, -1, -2, -1, -1,  0};

	return config_list[index];
}




