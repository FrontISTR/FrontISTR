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

#include "hecmw_vis_ucd_trans.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_voxel_gen.h"
#include "hecmw_malloc.h"


/* This program is used to transform GeoFEM UCD files into PVR data format */


void transform_ucd_pvr(struct hecmwST_result_data *data, double *node1,  struct hecmwST_local_mesh *mesh,
		Parameter_vr *vr,int my_rank, int pe_size, HECMW_Comm VIS_COMM, double *voxel_dxyz, double *voxel_orig_xyz,
		int *level, int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, int voxel_on, int display_range_on,
		int init_flag, int num_of_pvr)

{
	int          i,j;
	int          c_base;
	double        range[6], c_range[2];
	int           name_len, find;
	double      x, y, z;
	int          nodeid;


	int           tn_component;
	double        tmp;

	if((init_flag==1) || (num_of_pvr>1)) {
		if(my_rank==0) {
			fprintf(stderr, "start ucd-pvr transformation\n");
			/*      fprintf(stderr, "------in ucd_trans  color_comp_name  = %s\n", vr->color_comp_name);
			 */
		}


		if(vr->color_comp==-1) {
			find=0;
			if(strncmp(vr->color_comp_name, "NULL", 4)==0) {
				vr->color_comp=0;
				find=1;
			}
			else {

				for(j=0;j<data->nn_component;j++) {
					name_len=strlen(data->node_label[j]);
					if(strncmp(vr->color_comp_name, data->node_label[j], name_len)==0) {
						vr->color_comp=j;
						find=1;
						break;
					}
				}
			}
			if(find==0) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1053:the name for color component is not correct:%s\n", vr->color_comp_name);
				HECMW_vis_print_exit("please check it again");
			}
		}
		if(data->nn_dof[vr->color_comp]>1) {
			if(vr->color_subcomp==-1) {
				find=0;
				if(strncmp(vr->color_subcomp_name, "norm", 4)==0) {
					vr->color_subcomp=0;
					find=1;
				}
				else if(strncmp(vr->color_subcomp_name, "x", 1)==0) {
					vr->color_subcomp=1;
					find=1;
				}
				else if(strncmp(vr->color_subcomp_name, "y", 1)==0) {
					vr->color_subcomp=2;
					find=1;
				}
				else if(strncmp(vr->color_subcomp_name, "z", 1)==0) {
					vr->color_subcomp=3;
					find=1;
				}
				else
					vr->color_subcomp=1;
			}
		}
	}
	/*  fprintf(stderr, "*****ok here color_comp=%d  color_subcomp=%d\n", vr->color_comp, vr->color_subcomp);
	 */
	range[0]=range[2]=range[4]=1.0E17;
	range[1]=range[3]=range[5]=-1.0E17;
	c_range[0]=1.0E17;
	c_range[1]=-1.0E17;
	/*  for(i=0;i<nstep;i++) {
	  tstep=step_start+step_increment*i;
	 sprintf(infile, "%s.%d.%d.inp", resultfile, tstep, my_rank);
     if ((infp = fopen(infile, "r")) == NULL) {
        fprintf(stderr, "There is no such an input data file %s:\n", infile);
        exit (0);
	 }

     v = (struct visual_buf *)HECMW_malloc(sizeof(struct visual_buf));
     color_comp=read_ucd(infp, v, range, c_range, str_color, color_subcomp);
     fclose(infp);
  fprintf(outfp, "%d\n", v->mesh->n_node);
  fprintf(outfp, "%d\n", tstep);
  fprintf(outfp, "%lf\n", 0.0);
	 */

	if(vr->color_comp>=data->nn_component)
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1054: color_comp is wrong: >nn_component");
	if(vr->color_subcomp>data->nn_dof[vr->color_comp])
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1055: color_subcomp is wrong: >dof");
	if(data->nn_dof[vr->color_comp]==1) {
		c_base = 0;
		for(i=0;i<vr->color_comp;i++)
			c_base+=data->nn_dof[i];
		vr->color_subcomp=1;
	}
	else if(data->nn_dof[vr->color_comp]>1) {
		c_base=0;
		for(i=0;i<vr->color_comp;i++)
			c_base+=data->nn_dof[i];

	}
	tn_component=0;
	for(i=0;i<data->nn_component;i++)
		tn_component+=data->nn_dof[i];

	if((data->nn_dof[vr->color_comp]>1) && (vr->color_subcomp==0)) {
		for(i=0;i<mesh->n_node;i++) {
			node1[i]=0.0;
			for(j=0;j<data->nn_dof[vr->color_comp];j++) {
				tmp=data->node_val_item[c_base+i*tn_component+j];
				node1[i]+=tmp*tmp;
			}
			node1[i]=sqrt(node1[i]);
			if(node1[i]<c_range[0]) c_range[0]=node1[i];
			if(node1[i]>c_range[1]) c_range[1]=node1[i];
		}
	}

	else if(data->nn_dof[vr->color_comp]>1) {
		for(i=0;i<mesh->n_node;i++) {
			node1[i]= data->node_val_item[c_base+i*tn_component+(vr->color_subcomp-1)];
			if(node1[i]<c_range[0]) c_range[0]=node1[i];
			if(node1[i]>c_range[1]) c_range[1]=node1[i];
		}
	}
	else if(data->nn_dof[vr->color_comp]==1) {
		for(i=0;i<mesh->n_node;i++) {
			node1[i]= data->node_val_item[c_base + i*tn_component];
			if(node1[i]<c_range[0]) c_range[0]=node1[i];
			if(node1[i]>c_range[1]) c_range[1]=node1[i];
		}
	}
	if(my_rank==0)
		fprintf(stderr, " colorminmax=%lf %lf\n", c_range[0], c_range[1]);
	/*  free_v_node(v);
	 */
	 HECMW_Barrier(VIS_COMM);
	 /*---------find minmax value of mesh  ----------------*/
	 for(i=0;i<mesh->n_elem;i++) {
		 if(mesh->elem_type[i]<400) {
			 for(j=mesh->elem_node_index[i];j<mesh->elem_node_index[i+1];j++) {
				 nodeid=mesh->elem_node_item[j];
				 x=mesh->node[(nodeid-1)*3];
				 y=mesh->node[(nodeid-1)*3+1];
				 z=mesh->node[2+(nodeid-1)*3];
				 if(x>range[1]) range[1]=x;
				 if(x<range[0]) range[0]=x;
				 if(y>range[3]) range[3]=y;
				 if(y<range[2]) range[2]=y;
				 if(z>range[5]) range[5]=z;
				 if(z<range[4]) range[4]=z;
			 }
		 }
	 }

	 voxel_gen(range, c_range, vr->nv_xyz, voxel_dxyz, voxel_orig_xyz, level,
			 voxel_n_neighbor_pe, voxel_neighbor_pe, VIS_COMM, voxel_on,display_range_on, vr->display_range);


	 return;
}


