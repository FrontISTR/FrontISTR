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

#include "hecmw_vis_tetra_intersect.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_case_table.h"
#include "hecmw_malloc.h"

int get_tetra_data(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int elemID, Tetra *tetra, int tn_component)
{
	int		i, j;
	int		nodeID;
	int		c_base, s_base;
	double   tmp;
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
	for (i = 0; i < 4; i++) {
		nodeID = mesh->elem_node_item[mesh->elem_node_index[elemID]+i];
		tetra->local_vid[i]=nodeID-1;
		if(sff->display_way!=4) {
			if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
				tetra->c_data[i]=0.0;
				for(j=0;j<data->nn_dof[sff->color_comp];j++) {
					tmp=data->node_val_item[(nodeID-1)*tn_component+c_base+j];
					tetra->c_data[i]+=tmp*tmp;
				}
				tetra->c_data[i]=sqrt(tetra->c_data[i]);
			}
			else if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp!=0))
				tetra->c_data[i]=data->node_val_item[(nodeID-1)*tn_component+c_base+(sff->color_subcomp-1)];
			else if(data->nn_dof[sff->color_comp]==1)
				tetra->c_data[i] = data->node_val_item[(nodeID-1)*tn_component+c_base];

		}
		else if(sff->display_way==4)
			tetra->c_data[i]=sff->specified_color;

		tetra->axis[i*3] = mesh->node[(nodeID-1)*3];
		tetra->axis[i*3+1] = mesh->node[(nodeID-1)*3+1];
		tetra->axis[i*3+2] = mesh->node[(nodeID-1)*3+2];

		if(sff->surface_style==2) {

			if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp==0)) {
				tetra->s_data[i]=0.0;
				for(j=0;j<data->nn_dof[sff->data_comp];j++) {
					tmp=data->node_val_item[(nodeID-1)*tn_component+s_base+j];
					tetra->s_data[i]+=tmp*tmp;
				}
				tetra->s_data[i]=sqrt(tetra->s_data[i]);
			}
			else if((data->nn_dof[sff->data_comp]>1) && (sff->data_subcomp!=0))
				tetra->s_data[i] = data->node_val_item[(nodeID-1)*tn_component+s_base+(sff->data_subcomp-1)];
			else if(data->nn_dof[sff->data_comp]==1)
				tetra->s_data[i]=data->node_val_item[(nodeID-1)*tn_component+s_base];
		}
		else if(sff->surface_style==3)

			tetra->s_data[i]=get_value_equ(sff->cont_equ,sff->cross_type,tetra->axis[i*3],
					tetra->axis[i*3+1], tetra->axis[i*3+2]);
	}

	tetra->elem_id[0] = mesh->elem_ID[elemID*2];
	tetra->elem_id[1] = mesh->elem_ID[elemID*2+1];


	return 1;
}


static double intersect_line(int v0, int v1, double isovalue, Tetra *tetra, double point[3]) {
	int j;
	double rate, color;

	if(fabs(tetra->s_data[v1]-tetra->s_data[v0])<EPSILON)
		HECMW_vis_print_exit("There is something wrong in data precison");
	else {
		rate=(isovalue-tetra->s_data[v0])/(tetra->s_data[v1]-tetra->s_data[v0]);
		for(j=0;j<3;j++)
			point[j]=rate*(tetra->axis[v1*3+j]-tetra->axis[v0*3+j])+tetra->axis[v0*3+j];
		color=rate*(tetra->c_data[v1]-tetra->c_data[v0])+tetra->c_data[v0];
	}
	return (color);
}





void find_intersection_tetra(Tetra *tetra, double isovalue, Tetra_point *tetra_point, Head_patch_tetra *head_patch_tetra,
		Hash_vertex *vertex_hash_table)
{
	int  i,j,k;
	int  num_gt_0, v0_id;
	int  tmp_int;
	double point[3], color;
	Hash_vertex *h1, *h2;
	Tetra_point  *p1, *p2;
	Patch_tetra *t1, *t2;
	int  flag_existing, index_patch, patch[4];
	int  v1, v2, v[4], count_minus, count_plus;
	double  fp[4][3], n_f[3], n_norm, f_cen_point[3], sign_f, n_c, c_c[3];
	num_gt_0=0;
	for(i=0;i<4;i++) {
		if(tetra->s_data[i]>isovalue)
			num_gt_0++;
	}

	if((num_gt_0!=0) && (num_gt_0!=4)) { /* There are patches inside the tetra */

		if((num_gt_0==1) || (num_gt_0==3)) { /* There is one patch inside the tetra */
			t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
			if(t1==NULL)
				HECMW_vis_memory_exit("t1");
			t2=head_patch_tetra->patch_link;
			head_patch_tetra->num_patch++;
			head_patch_tetra->patch_link=t1;
			t1->next_patch=t2;

			if(num_gt_0==1) {
				v0_id=-1;
				for(i=0;i<4;i++) {
					if(tetra->s_data[i]>isovalue)
						v0_id=i;
				}
			}
			else if(num_gt_0==3) {
				v0_id=-1;
				for(i=0;i<4;i++) {
					if(tetra->s_data[i]<=isovalue)
						v0_id=i;
				}
			}
			/* find intersection */
			index_patch=-1;
			for(i=0;i<4;i++) {
				if(i!=v0_id) {
					color=intersect_line(v0_id, i, isovalue,  tetra, point);
					tmp_int=tetra->local_vid[v0_id]+tetra->local_vid[i];
					/*					if(vertex_hash_table[tmp_int].ident==0) {
						vertex_hash_table[tmp_int].ident++;
						h1=(Hash_vertex *)HECMW_malloc(sizeof(Hash_vertex));
						if(h1==NULL)
							HECMW_vis_memory_exit("h1");
					    vertex_hash_table[tmp_int].next_vertex=h1;
						h1->next_vertex=NULL;
						h1->ident=tetra_point->ident;
						for(j=0;j<3;j++)
							h1->geom[j]=point[j];
						p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
						if(p1==NULL)
							HECMW_vis_memory_exit("p1");
						tetra_point->ident++;
						p2=tetra_point->nextpoint;
						tetra_point->nextpoint=p1;
						p1->cdata=color;
						p1->ident=tetra_point->ident-1;
						for(j=0;j<3;j++)
							p1->geom[j]=point[j];
						p1->nextpoint=p2;
						index_patch++;
						t1->patch[index_patch]=p1->ident;
					}
					else if(vertex_hash_table[tmp_int].ident>0) { */
					flag_existing=0;
					h1=vertex_hash_table[tmp_int].next_vertex;
					for(j=0;j<vertex_hash_table[tmp_int].ident;j++) {
						if((fabs(h1->geom[0]-point[0])<EPSILON) && (fabs(h1->geom[1]-point[1])<EPSILON) &&
								(fabs(h1->geom[2]-point[2])<EPSILON)) {
							flag_existing=1;
							index_patch++;
							patch[index_patch]=h1->ident;
							for(k=0;k<3;k++)
								fp[index_patch][k]=point[k];
							break;
						}
						h1=h1->next_vertex;
					}
					if(flag_existing==0) { /*adding new vertex */
						vertex_hash_table[tmp_int].ident++;
						h1=(Hash_vertex *)HECMW_malloc(sizeof(Hash_vertex));
						if(h1==NULL)
							HECMW_vis_memory_exit("h1");
						h2=vertex_hash_table[tmp_int].next_vertex;
						vertex_hash_table[tmp_int].next_vertex=h1;
						h1->next_vertex=h2;
						h1->ident=tetra_point->ident;
						for(j=0;j<3;j++)
							h1->geom[j]=point[j];
						p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
						if(p1==NULL)
							HECMW_vis_memory_exit("p1");
						tetra_point->ident++;
						p2=tetra_point->nextpoint;
						tetra_point->nextpoint=p1;
						p1->cdata=color;
						p1->ident=tetra_point->ident-1;
						for(j=0;j<3;j++)
							p1->geom[j]=point[j];
						p1->nextpoint=p2;
						index_patch++;
						patch[index_patch]=p1->ident;
						for(k=0;k<3;k++)
							fp[index_patch][k]=point[k];
					}
					/*					}
					 */				}

			}/* end for 0, 4 */
			/*  judge whether the rotation direction pointed out to isosurface */

			n_f[0]=(fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
			n_f[1]=-(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
			n_f[2]=(fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
			for(j=0;j<3;j++)
				fp[3][j]=tetra->axis[v0_id*3+j];
			n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
			if(fabs(n_norm)>EPSILON) {
				for(j=0;j<3;j++)
					n_f[j]/=n_norm;
			}
			/*selce the direction point to inside the element */
			for(j=0;j<3;j++)
				f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j])/3.0;
			for(j=0;j<3;j++)
				c_c[j]=fp[3][j]-f_cen_point[j];
			n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
			if(fabs(n_c)>EPSILON) {
				for(j=0;j<3;j++)
					c_c[j]/=n_c;
			}
			sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
			if(((tetra->s_data[v0_id]>isovalue) && (sign_f<-EPSILON)) || ((tetra->s_data[v0_id]<=isovalue) && (sign_f>EPSILON))) {
				tmp_int=patch[1];
				patch[1]=patch[2];
				patch[2]=tmp_int;
			}
			for(j=0;j<3;j++)
				t1->patch[j]=patch[j];



		} /* end if (num_gt_0==1 or 3) */

		else if (num_gt_0==2) { /* Two patches inside the tetra */
			index_patch=-1;
			count_minus=0; count_plus=0;
			for(i=0;i<4;i++) {
				if(tetra->s_data[i]<=isovalue) {
					v[count_minus]=i;
					count_minus++;
				}
				else {
					v[2+count_plus]=i;
					count_plus++;
				}
			}
#ifdef old
			/*   judge the tetrahedra whether in the right rotation side */
			for(i=0;i<4;i++)
				for(j=0;j<3;j++)
					fp[i][j]=tetra->axis[v[i]*3+j];
			n_f[0]=(fp[2][1]-fp[1][1])*(fp[1][2]-fp[0][2])-(fp[1][1]-fp[0][1])*(fp[2][2]-fp[1][2]);
			n_f[1]=-(fp[2][0]-fp[1][0])*(fp[1][2]-fp[0][2])+(fp[1][0]-fp[0][0])*(fp[2][2]-fp[1][2]);
			n_f[2]=(fp[2][0]-fp[1][0])*(fp[1][1]-fp[0][1])-(fp[1][0]-fp[0][0])*(fp[2][1]-fp[1][1]);
			n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
			if(fabs(n_norm)>EPSILON) {
				for(j=0;j<3;j++)
					n_f[j]/=n_norm;
			}
			/*selce the direction point to inside the element */
			for(j=0;j<3;j++)
				f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j])/3.0;
			for(j=0;j<3;j++)
				c_c[j]=fp[3][j]-f_cen_point[j];
			n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
			if(fabs(n_c)>EPSILON) {
				for(j=0;j<3;j++)
					c_c[j]/=n_c;
			}
			sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
			if(sign_f<-EPSILON) {
				tmp_int=v[2];
				v[2]=v[3];
				v[3]=tmp_int;
			}
#endif

			for(i=0;i<4;i++) {
				if(i==0) {
					v1=v[0];  v2=v[2];
				}
				else if(i==1) {
					v1=v[1];  v2=v[2];
				}
				else if (i==2) {
					v1=v[1]; v2=v[3];
				}
				else if(i==3) {
					v1=v[0];  v2=v[3];
				}


				color=intersect_line(v1, v2, isovalue,  tetra, point);
				tmp_int=tetra->local_vid[v1]+tetra->local_vid[v2];
				/*					if(vertex_hash_table[tmp_int].ident==0) {
						vertex_hash_table[tmp_int].ident++;
						h1=(Hash_vertex *)HECMW_malloc(sizeof(Hash_vertex));
						if(h1==NULL)
							HECMW_vis_memory_exit("h1");
					    vertex_hash_table[tmp_int].next_vertex=h1;
						h1->next_vertex=NULL;
						h1->ident=tetra_point->ident;
						for(j=0;j<3;j++)
							h1->geom[j]=point[j];
						p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
						if(p1==NULL)
							HECMW_vis_memory_exit("p1");
						tetra_point->ident++;
						p2=tetra_point->nextpoint;
						tetra_point->nextpoint=p1;
						p1->cdata=color;
						p1->ident=tetra_point->ident-1;
						for(j=0;j<3;j++)
							p1->geom[j]=point[j];
						p1->nextpoint=p2;
						index_patch++;
						patch[index_patch]=p1->ident;
					}
					else if(vertex_hash_table[tmp_int].ident>0) {
				 */
				flag_existing=0;
				h1=vertex_hash_table[tmp_int].next_vertex;
				for(j=0;j<vertex_hash_table[tmp_int].ident;j++) {
					if((fabs(h1->geom[0]-point[0])<EPSILON) && (fabs(h1->geom[1]-point[1])<EPSILON) &&
							(fabs(h1->geom[2]-point[2])<EPSILON)) {
						flag_existing=1;
						index_patch++;
						patch[index_patch]=h1->ident;
						for(k=0;k<3;k++)
							fp[index_patch][k]=point[k];
						break;
					}
					h1=h1->next_vertex;
				}
				if(flag_existing==0) { /*adding new vertex */
					vertex_hash_table[tmp_int].ident++;
					h1=(Hash_vertex *)HECMW_malloc(sizeof(Hash_vertex));
					if(h1==NULL)
						HECMW_vis_memory_exit("h1");
					h2=vertex_hash_table[tmp_int].next_vertex;
					vertex_hash_table[tmp_int].next_vertex=h1;
					h1->next_vertex=h2;
					h1->ident=tetra_point->ident;
					for(j=0;j<3;j++)
						h1->geom[j]=point[j];
					p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
					if(p1==NULL)
						HECMW_vis_memory_exit("p1");
					tetra_point->ident++;
					p2=tetra_point->nextpoint;
					tetra_point->nextpoint=p1;
					p1->cdata=color;
					p1->ident=tetra_point->ident-1;
					for(j=0;j<3;j++)
						p1->geom[j]=point[j];
					p1->nextpoint=p2;
					index_patch++;
					patch[index_patch]=p1->ident;
					for(k=0;k<3;k++)
						fp[index_patch][k]=point[k];
					/*						}
					}
					 */
				}
			} /* end for i*/

			/*  judge whether the rotation direction pointed out to isosurface */

			n_f[0]=(fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
			n_f[1]=-(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
			n_f[2]=(fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
			for(j=0;j<3;j++)
				fp[3][j]=tetra->axis[0*3+j];
			n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
			if(fabs(n_norm)>EPSILON) {
				for(j=0;j<3;j++)
					n_f[j]/=n_norm;
			}
			/*selce the direction point to inside the element */
			for(j=0;j<3;j++)
				f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j])/3.0;
			for(j=0;j<3;j++)
				c_c[j]=fp[3][j]-f_cen_point[j];
			n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
			if(fabs(n_c)>EPSILON) {
				for(j=0;j<3;j++)
					c_c[j]/=n_c;
			}
			sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
			if(((tetra->s_data[0]>isovalue) && (sign_f<-EPSILON)) || ((tetra->s_data[0]<=isovalue) && (sign_f>EPSILON))) {
				tmp_int=patch[1];
				patch[1]=patch[3];
				patch[3]=tmp_int;
			}


			t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
			if(t1==NULL)
				HECMW_vis_memory_exit("t1");
			t2=head_patch_tetra->patch_link;
			head_patch_tetra->num_patch++;
			head_patch_tetra->patch_link=t1;
			t1->next_patch=t2;
			for(j=0;j<3;j++)
				t1->patch[j]=patch[j];
			t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
			if(t1==NULL)
				HECMW_vis_memory_exit("t1");
			t2=head_patch_tetra->patch_link;
			head_patch_tetra->num_patch++;
			head_patch_tetra->patch_link=t1;
			t1->next_patch=t2;
			t1->patch[0]=patch[0];
			t1->patch[1]=patch[2];
			t1->patch[2]=patch[3];


		} /* end if (num_gt_0==2) */
	}
	return;
}

