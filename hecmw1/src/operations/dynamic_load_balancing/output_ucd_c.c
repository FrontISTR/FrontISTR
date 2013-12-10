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




void write_dist_mesh_display(char *outfile, struct hecmwST_local_mesh *new_mesh, struct hecmwST_result_data *new_data) 
{
	int i,j,k;
	FILE *fp;
	int  id_elem;
	int flag_hit;
	int tmp_count, n_elem_type;
	int  tn_component, ne_internal, tmp_int;

  if ((fp = fopen(outfile, "w")) == NULL) {
    fprintf(stderr, "output file error: %s\n", outfile);
    exit (1);
  }

  tn_component=0;
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
  fclose(fp);
  return;
}

/*
void write_one_mesh_display(char *outfile, struct hecmwST_local_mesh *new_mesh, struct hecmwST_result_data *new_data,
							MPI_Comm VIS_COMM, int mynode, int pesize) 
{
	int i,j,k;
	FILE *fp;
	int  id_elem;
	int flag_hit;
	int tmp_count, n_elem_type;
	int  tn_component, ne_internal, tmp_int;
	int n_node, n_elem;
	int  nvtxs, *new_vtxdist, tmp_sum, tmp_nvtxs, *tmp_data;

	if(mynode==0) {
  if ((fp = fopen(outfile, "w")) == NULL) {
    fprintf(stderr, "output file error: %s\n", outfile);
    exit (1);
  }
    new_vtxdist=(int *)calloc(pesize+1, sizeof(int));
	nvtxs=new_mesh->nn_internal;
	if(mynode==0) {
		new_vtxdist[0]=0;
		new_vtxdist[1]=nvtxs;
		tmp_sum=nvtxs;
		for(i=1;i<pesize;i++) {
			MPI_Recv(&tmp_nvtxs, 1, MPI_INT, i, MPI_ANY_TAG, VIS_COMM, &stat);
			tmp_sum+=tmp_nvtxs;
			new_vtxdist[i+1]=tmp_sum;
		}
		for(i=1;i<pesize;i++)
	       MPI_Send(new_vtxdist,pesize+1,MPI_INT, i, 0, VIS_COMM);
	}
	else {
		MPI_Send(&nvtxs, 1, MPI_INT, 0, 0, VIS_COMM);
        MPI_Recv(new_vtxdist, pesize+1, MPI_INT, 0, MPI_ANY_TAG, VIS_COMM, &stat);
	}

  tn_component=0;
  for(i=0;i<new_data->nn_component;i++)
	  tn_component+=new_data->nn_dof[i];
	}
  ne_internal=0;
  for(i=0;i<new_mesh->ne_internal;i++) {
	  if(new_mesh->adapt_type[new_mesh->elem_internal_list[i]-1]==0)
		  ne_internal++;
  }
	if(pesize>1)
	  MPI_Allreduce(&new_mesh->nn_internal, &n_node, 1, MPI_INT, MPI_SUM, VIS_COMM);
	else
		n_node=mesh->nn_internal;
	if(pesize>1)
	  MPI_Allreduce(&ne_internal, &n_elem, 1, MPI_INT, MPI_SUM, VIS_COMM);
	else
		n_elem=ne_internal;

	if(mynode==0) {
  fprintf(fp, "%d %d %d %d %d\n", n_node, n_elem, tn_component, 0, 0);
  for(i=0;i<new_mesh->nn_internal;i++)
	  fprintf(fp, "%d %lf %lf %lf\n", i+1, new_mesh->node[i*3], new_mesh->node[i*3+1], new_mesh->node[i*3+2]);
  for(i=1;i<pesize;i++) {
	  tmp_int=new_vtxdist[i]-new_vtxdist[i-1];
	  tmp_data=(double *)calloc(tmp_int*3, sizeof(double));
	  if(tmp_data==NULL)
		  memory_exit("tmp_data");
	   MPI_Recv(tmp_data, tmp_int*3, MPI_DOUBLE, i, MPI_ANY_TAG, VIS_COMM, &stat);
	   for(j=0;j<tmp_int;j++)
		   fprintf(fp, "%d %lf %lf %lf\n", new_vtxdist[i-1]+j+1, tmp_data[j*3], tmp_data[j*3+1], tmp_data[j*3+2]);
	   free(tmp_data);
  }
	}
	else if(mynode!=0) {
		tmp_data=(double *)calloc(new_mesh->nn_internal*3, sizeof(double));
		if(tmp_data==NULL)
			memory_exit("tmp_data");
		for(j=0;j<new_mesh->nn_internal;j++) {
			for(k=0;k<3;k++)
				tmp_data[j*3+k]=new_mesh->node[j*3+k]);
		}
		MPI_Send(tmp_data, new_mesh->nn_internal*3, MPI_DOUBLE, 0, 0, VIS_COMM);
		free(tmp_data);
	}







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
  fclose(fp);
  return;
}
*/
