/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Dynamic Load Balancing                            *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_repart.h"
/*
extern void write_mesh_display(char *outfile, struct local_mesh *mesh, struct node_elem_data *node); 
extern void write2_mesh_display(char *outfile, struct local_mesh *mesh, struct node_elem_data *node); 
*/

extern struct hecmwST_local_mesh *mesh;
extern struct hecmwST_result_data *data;
extern struct hecmwST_local_mesh *new_mesh;
extern struct hecmwST_result_data *new_data;

extern void mesh_migration(int mynode, int pesize, Result_part *result, int *vtxdist);
extern void mesh_migration_adapt(int mynode, int pesize, Result_part *result, int *vtxdist);
extern void write_dist_mesh_display(char *outfile, struct hecmwST_local_mesh *new_mesh, struct hecmwST_result_data *new_data); 
extern void write_one_mesh_display(char *outfile, struct hecmwST_local_mesh *new_mesh, struct hecmwST_result_data *new_data,
							HECMW_Comm VIS_COMM, int mynode, int pesize); 
extern void HECMW_put_result_from_structure(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, 
											char *resultfile_dist);


/*
extern void mesh_migration(int pesize, int mynode, struct local_mesh *mesh, struct node_elem_data *node, struct local_mesh *new_mesh, 
					struct node_elem_data *new_node,  Result_part *result, 
					int *vtxdist, int *new_vtxdist, int *global_new2old, HECMW_Comm repart_comm);
extern void grp_migration(int pesize, int mynode, int t_node, struct local_mesh *new_mesh, struct local_mesh *mesh, struct grp_data *grp, 
				   struct grp_data *new_grp, int *vtxdist,int *new_vtxdist, 
				   int *global_new2old, HECMW_Comm repart_comm);



void redistribute_mesh(struct local_mesh *mesh, struct grp_data *grp, GraphType *graph, 
				Result_part *result, struct local_mesh *new_mesh, struct grp_data *new_grp, struct node_elem_data *new_node, 
				 HECMW_Comm VIS_COMM)
				 */



void redistribute_mesh(GraphType *graph, Result_part *result, int mynode, int pesize)

{
	int i,j,k,m;
	Adj_find *p1, *p2;
  HECMW_Status	stat;
  int *tmp_index, tmp_nvtxs, tmp_sum_node, tmp_sum_adj, tmp_sum_elem, tmp_pe, tmp_lid, tmp_sum;
   FILE *fp_test, *fp_orig, *fp;
   char test_file[128], orig_file[128];
    int *send_node_num, *send_node, *count_node, *send_adj_num, *count_adj, *send_adj, *send_elem_num, *send_elem, *count_elem;
	int *recv_node_num, *recv_elem_num, *recv_elem, *recv_node;
	int *flag_hit;
	int n_elem_type;
	int *new_part;
	int tmp_int, *tmp_send, *tmp_recv;
	int max_global_elem_num, new_n_elem;
	int *elem_keep_flag, *elem_hit_flag, *node_global_hit;
	Import_link_struct *import_link, *p3, *p4;
	int *recv_import_num, *recv_import, *send_import_num, *send_import, *send_recv_num, import_num;
	double *tmp_send_d, *tmp_recv_d;
	int *global_new2old, *new_vtxdist, *vtxdist;
	int  tn_component, ne_internal, tmp_count;
	char *resultfile, buf2[HECMW_FILENAME_LEN], resultfile_dist[HECMW_FILENAME_LEN];

	if(mynode==0) 
		fprintf(stderr, "Start migration for generating new mesh among PEs \n");
/*    update_ne_internal(mynode, mesh);
	fp_orig=fopen(orig_filename, "r");
	for(i=0;i<mesh->n_elem;i++) 
		fprintf(fp_orig, "n_elem=%d  ne_internal=%d\n", i+1, mesh->ne_internal_list[i]);
	fclose(fp_orig);
	*/
/*	 

	write2_mesh_display(orig_filename, mesh, node);


	sprintf(orig_filename, "orig_mesh.%d.inp", mynode);

	write_mesh_display(orig_filename, mesh, node);
*/
	vtxdist=(int *)calloc(pesize+1, sizeof(int));
		for(i=0;i<pesize+1;i++)
			vtxdist[i]=graph->vtxdist[i];
/*	new_vtxdist=(int *)calloc(pesize+1, sizeof(int));

	global_new2old=(int *)calloc(result->t_node, sizeof(int));
	if((new_vtxdist==NULL) || (global_new2old==NULL))
		HECMW_dlb_memory_exit("global_new2old");
		*/
/*	if(mesh->hecmw_flag_adapt==0)
       mesh_migration(mynode,pesize, result, vtxdist);
	else if(mesh->hecmw_flag_adapt==1)
	*/
    mesh_migration_adapt(mynode, pesize, result, vtxdist);
    if(mynode==0)
       fprintf(stderr, "Start output balanced mesh and result to files\n");
	HECMW_put_mesh(new_mesh, "mesh-dlb-out");

	resultfile=HECMW_ctrl_get_result_fileheader("result-out", buf2, HECMW_FILENAME_LEN);
	sprintf(resultfile_dist, "%s.%d", resultfile, mynode);
	HECMW_put_result_from_structure(new_mesh, new_data, resultfile_dist);

	   /*   free(node->data);
   free(node->n_free);
	
	grp_migration(pesize, mynode, result->t_node, new_mesh, mesh, grp, new_grp,  graph->vtxdist, new_vtxdist, global_new2old,
		repart_comm);

   free(global_new2old);
   free(new_vtxdist);
   */
/*   sprintf(orig_file, "orig_mesh.%d.inp", mynode);
   write_dist_mesh_display(orig_file, mesh, data); 

   sprintf(test_file, "bala_mesh.%d.inp", mynode);
   write_dist_mesh_display(test_file, new_mesh, new_data); 
*/
/*  sprintf(test_file, "test_result.%d", mynode);

  if ((fp = fopen(test_file, "w")) == NULL) {
    fprintf(stderr, "output file error: %s\n", test_file);
    exit (1);
  }
  */
/*  tn_component=0;
  for(i=0;i<new_data->nn_component;i++)
	  tn_component+=new_data->nn_dof[i];
  ne_internal=0;
  for(i=0;i<new_mesh->ne_internal;i++) {
	  if(new_mesh->adapt_type[new_mesh->elem_internal_list[i]-1]==0)
		  ne_internal++;
  }

  fprintf(fp, "%d %d %d %d %d\n", new_mesh->n_node, ne_internal, tn_component, 0, 0);
  for(i=0;i<new_mesh->n_node;i++)
	  fprintf(fp, "%d %lf %lf %lf\n", i+1, new_mesh->node[i*3], new_mesh->node[i*3+1], new_mesh->node[i*3+2]);
  tmp_count=0;
  for(i=0;i<new_mesh->ne_internal;i++) {
	  if(new_mesh->adapt_type[new_mesh->elem_internal_list[i]-1]==0) {
			tmp_count++;
			fprintf(fp, "%d 1  ", tmp_count);
			tmp_int=new_mesh->elem_internal_list[i]-1;
			if((new_mesh->elem_node_index[tmp_int+1]-new_mesh->elem_node_index[tmp_int]) == 4)
				fprintf(fp, " tet ");
			else
				fprintf(fp, " prism ");
			for(j=new_mesh->elem_node_index[tmp_int];j<new_mesh->elem_node_index[tmp_int+1];j++)
				fprintf(fp, "%d ", new_mesh->elem_node_item[j]);
			fprintf(fp,"\n");		
		}
  }
  fprintf(fp, "%d   ", new_data->nn_component);
  for(i=0;i<new_data->nn_component;i++)
    fprintf(fp, "%d ", new_data->nn_dof[i]);
  fprintf(fp, "\n");
  for(i=0;i<new_data->nn_component;i++)
	  fprintf(fp, "%s,\n", new_data->node_label[i]);
  for(i=0;i<new_mesh->n_node;i++) {
	  fprintf(fp, "%d   ", i+1);
	  for(j=0;j<tn_component;j++)
		  fprintf(fp, "%lf  ", new_data->node_val_item[i*tn_component+j]);
	  fprintf(fp, "\n");
   }
 */
/*  fprintf(fp, "%d\n", new_mesh->n_neighbor_pe);
  for(i=0;i<new_mesh->n_neighbor_pe;i++)
	  fprintf(fp, "%d  ", new_mesh->neighbor_pe[i]);
  fprintf(fp, "\n");
  for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
	  fprintf(fp, "%d ", new_mesh->import_index[i]);
  fprintf(fp, "\n");
  for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
	  fprintf(fp, "%d ", new_mesh->export_index[i]);
  fprintf(fp, "\n");
  fprintf(fp, "=====import nodes ====\n");
  for(i=0;i<new_mesh->import_index[new_mesh->n_neighbor_pe];i++)
	  fprintf(fp, "%d\n", new_mesh->import_item[i]);
  fprintf(fp, "=====export nodes ====\n");
  for(i=0;i<new_mesh->export_index[new_mesh->n_neighbor_pe];i++)
	  fprintf(fp, "%d\n", new_mesh->export_item[i]);
  
  fclose(fp);
  */



	return;
	}
	
	














		   







		





