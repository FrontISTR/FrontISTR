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

/*----------------------------------------------------------------------
#     Subroutines in this file on isosurface generation for hexahedra by Marching Cubes is based
	  on the revision of Dr. Yuriko Takeshima's codes when she was working part time in RIST
#---------------------------------------------------------------------- */
#include "hecmw_vis_surface_compute.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_case_table.h"
#include "hecmw_vis_tetra_intersect.h"
#include "hecmw_vis_patch_const.h"
#include "hecmw_malloc.h"


int HECMW_vis_surface_compute(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *bdflag, int *sum_v, int *sum_t,int *tvertex, int *tpatch,
		double *minc, double *maxc, Result *result, int sf_i, int mynode,  HECMW_Comm VIS_COMM )

{
	/* in_volume */
	int		i, j,k, mm;
	Point **CS_verts_head;
	Polygon **CS_polys_head;
	int 		alpha_index, beta_index;
	int 		sum_verts;
	int 		sum_table;
	Point 	**CS_verts_tail, **CS_verts_refer;
	Polygon 	**CS_polys_tail;
	Cube_polygons alpha_cube, beta_cube;

	int		n_elem, tmp_int;
	Cell		*cell;
	Tetra     *tetra;


	int		disamb_flag;

	Polygon	*CS_polys, *CS_polys_tmp;
	Point		*CS_verts, *CS_verts_tmp;
	int		sum_polys, poly_num;
	int    aplist_size, bplist_size, cplist_size;
	int       num_nh_verts, num_nh_patch; /* vertex on patches generated in non-hexahedra */
	int       flag_hexa, flag_tetra;
	Hash_vertex *vertex_hash_table;
	Tetra_point  *tetra_point;
	Head_patch_tetra  *head_patch_tetra;
	Tetra_point  *p1, *p2;
	Patch_tetra *t1, *t2;
	Hash_vertex *h1, *h2;
	int    minus_patch;
	double  tmp_data[9];
	int   tn_component, c_base;
	double tmp;

	sum_verts=0;
	num_nh_verts=0;
	num_nh_patch=0;
	sum_polys=0;

	n_elem = mesh->n_elem;
	tn_component=0;
	for(i=0;i<data->nn_component;i++)
		tn_component+=data->nn_dof[i];

	/*  prism = (Prism *)HECMW_malloc(sizeof(Prism));
	 */
	if(mesh->elem_type[0]>300) {
		flag_hexa=0;
		flag_tetra=0;

		disamb_flag = 1;
		/*   fprintf(stderr, "sff inf: color=%d color_sub=%d con=%lf %lf %lf% lf\n", sff->color_comp, sff->color_subcomp, sff->cont_equ[6],
	   sff->cont_equ[7], sff->cont_equ[8],sff->cont_equ[9]);

   sprintf(test_file, "test_ucd.%d.%d.inp", time_step, mynode);
   write_mesh_display(test_file, mesh, data);
   HECMW_Barrier(VIS_COMM);
		 */
		/*
    for(i = 0; i < mesh->ne_internal; i++) {
		tmp_int=mesh->elem_internal_list[i]-1;
		 */
		for(i=0;i<mesh->n_elem;i++)
			if(mesh->elem_ID[i*2+1]==mynode) {
				tmp_int=i;
				if((mesh->elem_type[tmp_int]==361) || (mesh->elem_type[tmp_int]==362)
						|| (mesh->elem_type[tmp_int]==351) 	|| (mesh->elem_type[tmp_int]==352)) {
					if(flag_hexa==0) {
						flag_hexa=1;

						CS_polys_head = (Polygon **)HECMW_malloc(sizeof(Polygon *));
						CS_verts_head = (Point **)HECMW_calloc(TABLE_SIZE, sizeof(Point *));
						sum_table = TABLE_SIZE;
						CS_verts_tail = (Point **)HECMW_calloc(sum_table, sizeof(Point *));
						CS_verts_refer = (Point **)HECMW_calloc(sum_table, sizeof(Point *));
						for (j= 0; j < sum_table; j++) {
							if ((CS_verts_refer[j] = CS_verts_tail[j] = CS_verts_head[j]
							                                                          = alloc_verts(VERTEX_PACK)) == NULL)
								HECMW_vis_memory_exit("verts for hexahedra");
						}

						CS_polys_tail = (Polygon **)HECMW_malloc(sizeof(Polygon *));
						if ((*CS_polys_tail = *CS_polys_head =
							alloc_polygons(POLYGON_PACK)) == NULL)
							HECMW_vis_memory_exit("CS_polys_tail");

						alpha_cube.isosurf = (int **)HECMW_malloc(sizeof(int *));
						beta_cube.isosurf = (int **)HECMW_malloc(sizeof(int *));
						cell = (Cell *)HECMW_malloc(sizeof(Cell));
					}
					get_data(sff,mesh, data, tmp_int, cell, tn_component);
					if(sff->surface_style==2) {
						alpha_index=make_tile(mesh, cell, sff->iso_value, 0, &alpha_cube, disamb_flag);
						beta_index=make_tile(mesh, cell, sff->iso_value, 1, &beta_cube, disamb_flag);
						if (alpha_index && beta_index) {
							if ((alpha_index + beta_index) == HEX_NODE_INDEX) {
								if (!merge_vol_iso(0, cell, sff->iso_value, &alpha_cube, 0,
										&beta_cube, bdflag[i], &sum_verts,
										CS_verts_tail, CS_verts_refer,
										CS_verts_head, CS_polys_tail)) {
									return 0;
								}
							}
						}

					}
					else if(sff->surface_style==3) {
						alpha_index = make_tile(mesh, cell, 0, 0, &alpha_cube, disamb_flag);
						beta_index = make_tile(mesh, cell, 0, 1, &beta_cube, disamb_flag);
						if (alpha_index && beta_index) {
							if ((alpha_index + beta_index) == HEX_NODE_INDEX) {
								if (!merge_vol_iso(0, cell, 0, &alpha_cube, 0,
										&beta_cube, bdflag[i], &sum_verts,
										CS_verts_tail, CS_verts_refer,
										CS_verts_head, CS_polys_tail)) {
									return 0;
								}
							}
						}
					}
				}
				else if((mesh->elem_type[tmp_int]==341) || (mesh->elem_type[tmp_int]==342)) {
					/* tetrahedra */
					if(flag_tetra==0) {
						flag_tetra=1;
						/* initialize  */
						tetra = (Tetra *)HECMW_malloc(sizeof(Tetra));
						vertex_hash_table=(Hash_vertex *)HECMW_calloc(mesh->n_node*2, sizeof(Hash_vertex));
						tetra_point=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
						head_patch_tetra=(Head_patch_tetra *)HECMW_malloc(sizeof(Head_patch_tetra));
						if((vertex_hash_table==NULL) || (tetra==NULL) || (tetra_point==NULL) || (head_patch_tetra==NULL))
							HECMW_vis_memory_exit("initialize tetra");
						for(j=0;j<mesh->n_node*2;j++) {
							vertex_hash_table[j].ident=0;
							vertex_hash_table[j].next_vertex=NULL;
						}
						tetra_point->ident=0;
						head_patch_tetra->num_patch=0;
					}
					get_tetra_data(sff, mesh, data, tmp_int, tetra, tn_component);
					if(sff->surface_style==3)
						sff->iso_value=0.0;
					find_intersection_tetra(tetra, sff->iso_value, tetra_point, head_patch_tetra, vertex_hash_table);
				}
			}
		if(flag_tetra==1) {
			num_nh_verts=tetra_point->ident;
			num_nh_patch=head_patch_tetra->num_patch;
			HECMW_free(tetra);
			for(j=0;j<mesh->n_node*2;j++) {
				h1=vertex_hash_table[j].next_vertex;
				for(i=0;i<vertex_hash_table[j].ident;i++) {
					h2=h1->next_vertex;
					HECMW_free(h1);
					h1=h2;
				}
			}
			HECMW_free(vertex_hash_table);
		}



		if(flag_hexa>0) {
			mfree(alpha_cube.isosurf);  mfree(beta_cube.isosurf);
			mfree(CS_verts_tail);   mfree(CS_verts_refer);
			mfree(CS_polys_tail);


			*sum_v = sum_verts;
			*sum_t = sum_table;
			HECMW_free(cell);

			CS_polys = *CS_polys_head;

			sum_polys=0;

			aplist_size =  bplist_size =  cplist_size = 1;
			while (CS_polys->plist != NULL) {
				switch(CS_polys->type) {
				case 0: aplist_size += CS_polys->plist[0] + 1;
				sum_polys++;
				break;
				case 1: bplist_size += CS_polys->plist[0] + 1;
				sum_polys++;
				break;
				case 2: cplist_size += (CS_polys->plist[0] - 2) * 4;
				sum_polys += CS_polys->plist[0] - 2;
				break;
				}
				CS_polys = CS_polys->nextpolygon;
			}
		}
		/*sum_polys=sum_polys*2;*/
		if((sum_verts+num_nh_verts)>0) {
			result[sf_i].n_vertex=sum_verts+num_nh_verts;
			result[sf_i].n_patch=sum_polys+num_nh_patch;
			result[sf_i].vertex=(double *)HECMW_calloc(result[sf_i].n_vertex*3, sizeof(double));
			result[sf_i].patch=(int *)HECMW_calloc(result[sf_i].n_patch*3, sizeof(int));
			result[sf_i].color=(double *)HECMW_calloc(result[sf_i].n_vertex, sizeof(double));
			if((result[sf_i].vertex==NULL) || (result[sf_i].patch==NULL) || (result[sf_i].color==NULL))
				HECMW_vis_memory_exit("result");
		}
		mm=0;
		if(sum_verts>0) {
			/*  make main vertex table of CS  */
			for (i = 0; i < sum_table; i++) {
				CS_verts_tmp = CS_verts = CS_verts_head[i];
				j = 0;
				while (CS_verts->ident != 0) {
					/*      verts_geom[CS_verts->ident] = CS_verts->geom;
      verts_field[CS_verts->ident] = CS_verts->field;
      verts_color[CS_verts->ident] = CS_verts->cdata;

      fprintf(vfile, " %d %lf %lf %lf\n",  *tvertex+CS_verts->ident,
	      verts_geom[CS_verts->ident].x,
	      verts_geom[CS_verts->ident].y,
	      verts_geom[CS_verts->ident].z);
					 */
					result[sf_i].vertex[mm*3]=CS_verts->geom.x;
					result[sf_i].vertex[mm*3+1]=CS_verts->geom.y;
					result[sf_i].vertex[mm*3+2]=CS_verts->geom.z;
					result[sf_i].color[mm]=CS_verts->cdata;
					mm++;
					CS_verts = CS_verts->nextpoint;
					if (!((++j)%VERTEX_PACK)) {
						mfree(CS_verts_tmp);
						CS_verts_tmp = CS_verts;
					}
				}
				mfree(CS_verts_tmp);
			}
			mfree(CS_verts_head);

			/*  make polygon table and decide vertex ID
      of each object (alpha,beta,cross)  */
			/*  trilist=(Triangle *)HECMW_calloc(sum_polys+1,sizeof(Triangle));
			 */
			CS_polys_tmp = CS_polys = *CS_polys_head;
			i = 0;
			poly_num = 1;
			minus_patch=0;

			while (CS_polys->plist != NULL) {
				switch(CS_polys->type) {
				case 0:
					/*	      fprintf(pfile, "%d ", poly_num+(*tpatch));
					 */
					k=-1;
					for ( j =1; j <= CS_polys->plist[0];j++) {
						k++;
						/*		fprintf(pfile, "%d ", (*tvertex)+CS_polys->plist[j]);

		trilist[poly_num -1].vertex[k]=CS_polys->plist[j];
						 */
						if((sff->output_type==1) || (sff->output_type==2))
							result[sf_i].patch[(poly_num-1)*3+k]=*tvertex+CS_polys->plist[j];
						else
							result[sf_i].patch[(poly_num-1)*3+k]=CS_polys->plist[j];

						tmp_data[k*3]=result[sf_i].vertex[(CS_polys->plist[j]-1)*3];
						tmp_data[k*3+1]=result[sf_i].vertex[(CS_polys->plist[j]-1)*3+1];
						tmp_data[k*3+2]=result[sf_i].vertex[(CS_polys->plist[j]-1)*3+2];
					}
					if(((fabs(tmp_data[0]-tmp_data[3])<EPSILON) && (fabs(tmp_data[1]-tmp_data[4])<EPSILON) &&
							(fabs(tmp_data[2]-tmp_data[5])<EPSILON)) ||
							((fabs(tmp_data[0]-tmp_data[6])<EPSILON) && (fabs(tmp_data[1]-tmp_data[7])<EPSILON) &&
									(fabs(tmp_data[2]-tmp_data[8])<EPSILON)) ||
									((fabs(tmp_data[6]-tmp_data[3])<EPSILON) && (fabs(tmp_data[7]-tmp_data[4])<EPSILON) &&
											(fabs(tmp_data[8]-tmp_data[5])<EPSILON)))
						minus_patch++;
					else
						poly_num++;

					break;

				}

				CS_polys = CS_polys->nextpolygon;
				if (!((++i)%POLYGON_PACK)) {
					mfree(CS_polys_tmp->plist);

					mfree(CS_polys_tmp);
					CS_polys_tmp = CS_polys;
				}
			}
			/*#ifdef DEBUG
 fprintf(stderr, "the previous patch num is %d now is %d\n", result[sf_i].n_patch, result[sf_i].n_patch-minus_patch);
#endif
			 */
			result[sf_i].n_patch-=minus_patch;

			mfree(CS_polys_tmp);
		}
		if(num_nh_verts>0) {
			p1=tetra_point->nextpoint;
			for(i=0;i<tetra_point->ident;i++) {
				for(j=0;j<3;j++)
					result[sf_i].vertex[(sum_verts+p1->ident)*3+j]=p1->geom[j];
				result[sf_i].color[sum_verts+p1->ident]=p1->cdata;
				p2=p1;
				p1=p1->nextpoint;
				HECMW_free(p2);
			}
			HECMW_free(tetra_point);
		}
		if(num_nh_patch>0) {
			t1=head_patch_tetra->patch_link;
			for(i=0;i<head_patch_tetra->num_patch;i++) {
				for(j=0;j<3;j++) {
					if(sff->output_type==3)
						result[sf_i].patch[(sum_polys-minus_patch+i)*3+j]=t1->patch[j]+1+sum_verts;
					else
						result[sf_i].patch[(sum_polys-minus_patch+i)*3+j]=*tvertex+t1->patch[j]+1+sum_verts;
				}
				t2=t1;
				t1=t1->next_patch;
				HECMW_free(t2);
			}
			HECMW_free(head_patch_tetra);
		}


		/*
 fprintf(stderr, "On surface %d PE %d: n_vertex is %d  n_patch is %d\n", sf_i, mynode, result[sf_i].n_vertex, result[sf_i].n_patch);
		 */
		if(result[sf_i].n_vertex>0) {
			*minc=*maxc=result[sf_i].color[0];
			for(i=1;i<=result[sf_i].n_vertex;i++) {
				if(result[sf_i].color[i-1]<(*minc)) (*minc)=result[sf_i].color[i-1];
				if(result[sf_i].color[i-1]>(*maxc)) (*maxc)=result[sf_i].color[i-1];
			}
			/*
#ifdef DEBUG
  fprintf(stderr,"On surface %d PE %d: minimum color=%lf maximum color=%lf\n", sf_i, mynode, *minc,*maxc);
#endif
			 */
		}
	} /* endof if elem_type>300 */
	else if((mesh->elem_type[0]>200) && (mesh->elem_type[0]<300)) {
		result[sf_i].n_patch=0;
		for(i=0;i<n_elem;i++) {
			if(mesh->elem_type[i]==231)
				result[sf_i].n_patch++;
			else if(mesh->elem_type[i]==241)
				result[sf_i].n_patch+=2;
		}
		result[sf_i].n_vertex=mesh->n_node;
		result[sf_i].vertex=(double *)HECMW_calloc(mesh->n_node*3, sizeof(double));
		result[sf_i].color=(double *)HECMW_calloc(mesh->n_node, sizeof(double));
		result[sf_i].patch=(int *)HECMW_calloc(result[sf_i].n_patch*3, sizeof(int));
		if((result[sf_i].vertex==NULL) || (result[sf_i].color==NULL) || (result[sf_i].patch==NULL))
			HECMW_vis_memory_exit("result: vertex, color and patch");
		for(i=0;i<mesh->n_node;i++) {
			for(j=0;j<3;j++)
				result[sf_i].vertex[i*3+j]=mesh->node[i*3+j];
		}
		if(data->nn_dof[sff->color_comp]==1) {
			c_base = 0;
			for(i=0;i<sff->color_comp;i++)
				c_base+=data->nn_dof[i];
		}
		else if(data->nn_dof[sff->color_comp]>1) {
			c_base=0;
			for(i=0;i<sff->color_comp;i++)
				c_base+=data->nn_dof[i];

		}

		if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
			for(i=0;i<mesh->n_node;i++) {
				result[sf_i].color[i]=0.0;
				for(j=0;j<data->nn_dof[sff->color_comp];j++) {
					tmp=data->node_val_item[c_base+i*tn_component+j];
					result[sf_i].color[i]+=tmp*tmp;
				}
				result[sf_i].color[i]=sqrt(result[sf_i].color[i]);
			}
		}

		else if(data->nn_dof[sff->color_comp]>1) {
			for(i=0;i<mesh->n_node;i++) {
				result[sf_i].color[i]= data->node_val_item[c_base+i*tn_component+(sff->color_subcomp-1)];
			}
		}
		else if(data->nn_dof[sff->color_comp]==1) {
			for(i=0;i<mesh->n_node;i++) {
				result[sf_i].color[i]= data->node_val_item[c_base + i*tn_component];
			}
		}
		poly_num=0;
		for(i=0;i<n_elem;i++) {
			if(mesh->elem_type[i]==231) {
				for(j=0;j<3;j++)
					result->patch[poly_num*3+j]=mesh->elem_node_item[mesh->elem_node_index[i]+j]+*tvertex;
				poly_num++;
			}
			else if(mesh->elem_type[i]==241) {
				for(j=0;j<3;j++)
					result->patch[poly_num*3+j]=mesh->elem_node_item[mesh->elem_node_index[i]+j]+*tvertex;
				poly_num++;
				result->patch[poly_num*3+0]=mesh->elem_node_item[mesh->elem_node_index[i]+0]+*tvertex;
				result->patch[poly_num*3+1]=mesh->elem_node_item[mesh->elem_node_index[i]+2]+*tvertex;
				result->patch[poly_num*3+2]=mesh->elem_node_item[mesh->elem_node_index[i]+3]+*tvertex;
				poly_num++;
			}
		}
		/*
 fprintf(stderr, "n_vertex is %d  n_patch is %d\n", result[sf_i].n_vertex, result[sf_i].n_patch);
		 */
		if(result[sf_i].n_vertex>0) {
			*minc=*maxc=result[sf_i].color[0];
			for(i=1;i<=result[sf_i].n_vertex;i++) {
				if(result[sf_i].color[i-1]<(*minc)) (*minc)=result[sf_i].color[i-1];
				if(result[sf_i].color[i-1]>(*maxc)) (*maxc)=result[sf_i].color[i-1];
			}
			/*
#ifdef DEBUG
  fprintf(stderr,"On surface %d PE %d: minimum color=%lf maximum color=%lf\n", sf_i, mynode, *minc,*maxc);
#endif
			 */
		}
	}








	(*tvertex)+=result[sf_i].n_vertex;
	(*tpatch)+=result[sf_i].n_patch;
	/*  if(sum_verts>0) {
	 *minv=mincolor;
	 *maxv=maxcolor;
  }

 mfree(verts_info);
  mfree(verts_geom);
  mfree(verts_field);
  mfree(verts_color);
  mfree(trilist);
	 */

	return (1);
}



int HECMW_vis_chk_bounds(struct hecmwST_local_mesh *mesh, int *bdflag)
{
	int		i;
	int		n_elem;


	n_elem = mesh->n_elem;
	for (i = 0; i < n_elem; i++) {


		*(bdflag + i) =64;
	}
	return 1;

}

#ifdef old
int chk_node_data(struct visual_buf *v, int s_comp, int c_comp)
{
	int		*global_node_id;
	double	*shape_data;
	double	*color_data;
	int		imp_node;
	int		s_base, c_base;
	int		i, j, k;
	int 		ne;
	int		n_export_node, export_base;
	int 		n_import_node, import_base;

	HECMW_Status	stat;

	s_base = c_base = 0;

	for(i = 0; i < s_comp; i++) {
		s_base += v->node.n_free[i]*v->node.n;
	}
	for(i = 0; i < c_comp; i++) {
		c_base += v->node.n_free[i]*v->node.n;
	}

	export_base = 0;
	import_base = 0;
	for (i = 0; i < v->mesh->n_neighbor_pe; i++) {
		ne = v->mesh->neighbor_pe[i];
		n_export_node = v->mesh->export_index[i+1] - v->mesh->export_index[i];
		global_node_id = (int *)HECMW_calloc(n_export_node, sizeof(int));
		shape_data = (double *)HECMW_calloc(n_export_node, sizeof(double));
		color_data = (double *)HECMW_calloc(n_export_node, sizeof(double));
		for (j = 0; j < n_export_node; j++) {
			global_node_id[j]
			               = v->mesh->global_node_id[v->mesh->export_node[export_base + j] - 1];
			shape_data[j] = (double)*(v->node.data + s_base
					+ v->mesh->export_node[export_base + j] - 1);
			color_data[j] = (double)*(v->node.data + c_base
					+ v->mesh->export_node[export_base + j] - 1);
		}
		HECMW_Send(&n_export_node, 1, HECMW_INT, ne, 0, geofem_app_comm);
		HECMW_Send(global_node_id, n_export_node, HECMW_INT, ne, 0, geofem_app_comm);
		HECMW_Send(shape_data, n_export_node, HECMW_DOUBLE, ne, 0, geofem_app_comm);
		HECMW_Send(color_data, n_export_node, HECMW_DOUBLE, ne, 0, geofem_app_comm);

		HECMW_free(global_node_id);
		HECMW_free(shape_data);
		HECMW_free(color_data);
		export_base = v->mesh->export_index[i+1];

		n_import_node = v->mesh->import_index[i+1] - v->mesh->import_index[i];
		HECMW_Recv(&imp_node, 1, HECMW_INT, ne, HECMW_ANY_TAG, geofem_app_comm, &stat);

		if (imp_node != n_import_node) {
			GEOFEM_abort(711, "chk_node_data error\n");
		}
		global_node_id = (int *)HECMW_calloc(n_import_node, sizeof(int));
		shape_data = (double *)HECMW_calloc(n_import_node, sizeof(double));
		color_data = (double *)HECMW_calloc(n_import_node, sizeof(double));
		HECMW_Recv(global_node_id, n_import_node, HECMW_INT, ne, HECMW_ANY_TAG,
				geofem_app_comm, &stat);
		HECMW_Recv(shape_data, n_import_node, HECMW_DOUBLE, ne, HECMW_ANY_TAG,
				geofem_app_comm, &stat);
		HECMW_Recv(color_data, n_import_node, HECMW_DOUBLE, ne, HECMW_ANY_TAG,
				geofem_app_comm, &stat);
		for (k = import_base; k < import_base + n_import_node; k++) {
			for (j = 0; j < n_import_node; j++) {
				if (v->mesh->global_node_id[v->mesh->import_node[k] - 1]
				                            == global_node_id[j]) {
					v->node.data[s_base + v->mesh->import_node[k] - 1] = shape_data[j];
					v->node.data[c_base + v->mesh->import_node[k] - 1] = color_data[j];
					break;
				}
			}
		}
		HECMW_free(global_node_id);
		HECMW_free(shape_data);
		HECMW_free(color_data);
		import_base = v->mesh->import_index[i+1];
	}

	return 1;
}
#endif

/*
int chk_bounds(struct visual_buf *vol, int *bdflag)
{
  int	i, j, k, l;
  int	bd1, bd2;
  int	n_elem;

  n_elem = vol->mesh->n_elem;
  for (i = 0; i < n_elem; i++) {
    for (j = 0; j < i; j++) {
      if (*(bdflag + j) < HEX_FACE_INDEX) {
	bd1 = bd2 = 0;
	for (k = 0; k < HEX_N_NODE; k++) {
	  for (l = 0; l < HEX_N_NODE; l++) {
	    if ((*(vol->mesh->elem + i + n_elem * k))
		== (*(vol->mesh->elem + j + n_elem * l))) {
	      bd1 += 1 << k;
	      bd2 += 1 << l;
	    }
	  }
	}
	switch (bd1) {
	  case  15: *(bdflag + i) +=  1; break;
	  case 102: *(bdflag + i) +=  2; break;
	  case 240: *(bdflag + i) +=  4; break;
	  case 153: *(bdflag + i) +=  8; break;
	  case 204: *(bdflag + i) += 16; break;
	  case  51: *(bdflag + i) += 32; break;
	  default: break;
	}
	switch (bd2) {
	  case  15: *(bdflag + j) +=  1; break;
	  case 102: *(bdflag + j) +=  2; break;
	  case 240: *(bdflag + j) +=  4; break;
	  case 153: *(bdflag + j) +=  8; break;
	  case 204: *(bdflag + j) += 16; break;
	  case  51: *(bdflag + j) += 32; break;
	  default: break;
	}
      }

      if (*(bdflag + i) == HEX_FACE_INDEX) break;
    }
  }

  return 1;

}

 */
void find_isoline(int isonum, int sum_polys,double mincolor, double maxcolor, double *vcoord,
		int *plist, double *vcolor, Isohead *isohead)
{
	int i,j,k;
	double c[3];
	double isocolor,deltac;
	Fgeom g[3];

	if(isonum<=0) {
		fprintf(stderr, "isonumber input wrong\n");
		exit(0);
	}
	deltac=(maxcolor-mincolor)/(isonum+1);
	for(k=0;k<isonum;k++) {
		isocolor=mincolor+(k+1)*deltac;
		for(i=0;i<sum_polys;i++) {
			for(j=0;j<3;j++) {
				c[j]=vcolor[plist[i*3+j]-1];
				g[j].x=vcoord[(plist[i*3+j]-1)*3];
				g[j].y=vcoord[(plist[i*3+j]-1)*3+1];
				g[j].z=vcoord[(plist[i*3+j]-1)*3+2];
			}
			line_find(isocolor, c, g, k, isohead);
		}
	}
}

int find_line_segment(double f[3][3], double c[3], double isocolor, double iso_p[6])
{
	int j,m,j1;
	int pnum, flag;

	flag=0;
	pnum=-1;
	for(j=0;j<3;j++) {
		j1=j+1;
		if(j1==3) j1=0;
		if((fabs(c[j]-c[j1])<0.0000001) && (fabs(c[j]-isocolor)<0.0000001)) {
			/*this edge is isoline*/
			flag=1;

			for(m=0;m<3;m++)
				iso_p[m]=f[j][m];
			for(m=0;m<3;m++)
				iso_p[m+3]=f[j1][m];
			return 1;
		}
		else if(((c[j]>=isocolor) && (c[j1]<isocolor)) || ((c[j]<isocolor) && (c[j1]>=isocolor))){
			pnum++;
			for(m=0;m<3;m++)
				iso_p[pnum*3+m]=f[j][m]+(isocolor-c[j])/(c[j1]-c[j])*(f[j1][m]-f[j][m]);
		}
	}
	if(pnum==1)
		flag=1;


	return flag;
}



void line_find(double isocolor, double c[3], Fgeom g[3], int k, Isohead *isohead)
{
	Isoline *p1, *p2;
	int j,m,j1;
	int pnum;
	Fgeom isopoint[2];

	pnum=-1;
	for(j=0;j<2;j++){
		isopoint[j].x=isopoint[j].y=isopoint[j].z=0.0;
	}
	for(j=0;j<3;j++) {
		j1=j+1;
		if(j1==3) j1=0;
		if((fabs(c[j]-c[j1])<0.0000001) && (fabs(c[j]-isocolor)<0.0000001)) {
			/*this edge is isoline*/
			isohead[k].linenum++;
			p1=isohead[k].nextline;
			p2=(Isoline *)HECMW_malloc(sizeof(Isoline));

			p2->point[0].x=g[j].x;
			p2->point[0].y=g[j].y;
			p2->point[0].z=g[j].z;
			p2->point[1].x=g[j1].x;
			p2->point[1].y=g[j1].y;
			p2->point[1].z=g[j1].z;

			p2->nextline=p1;
			isohead[k].nextline=p2;
			return;
		}

		else if(((c[j]>=isocolor) && (c[j1]<isocolor)) || ((c[j]<isocolor) && (c[j1]>=isocolor))){
			pnum++;
			isopoint[pnum].x=g[j].x+(isocolor-c[j])/(c[j1]-c[j])*(g[j1].x-g[j].x);
			isopoint[pnum].y=g[j].y+(isocolor-c[j])/(c[j1]-c[j])*(g[j1].y-g[j].y);
			isopoint[pnum].z=g[j].z+(isocolor-c[j])/(c[j1]-c[j])*(g[j1].z-g[j].z);
		}
	}
	if(pnum==1) {
		isohead[k].linenum++;
		p1=isohead[k].nextline;
		p2=(Isoline *)HECMW_malloc(sizeof(Isoline));
		for(m=0;m<2;m++) {
			p2->point[m].x=isopoint[m].x;
			p2->point[m].y=isopoint[m].y;
			p2->point[m].z=isopoint[m].z;
		}
		p2->nextline=p1;
		isohead[k].nextline=p2;
	}

	return;
}

void isoline_free(int isonum,Isohead *isohead)
{
	Isoline *p1, *p2;
	int i;
	for(i=0;i<isonum;i++) {
		p1=isohead[i].nextline;
		while(p1!=NULL) {
			p2=p1;
			p1=p1->nextline;
			HECMW_free(p2);
		}
	}
	HECMW_free(isohead);
	return;
}









