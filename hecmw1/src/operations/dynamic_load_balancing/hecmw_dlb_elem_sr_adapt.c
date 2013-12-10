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

extern struct hecmwST_local_mesh *mesh;
extern struct hecmwST_result_data *data;
extern struct hecmwST_local_mesh *new_mesh;
extern struct hecmwST_result_data *new_data;
extern int  stack_part_send_recv(int neibpetot, int *neibpe, int *stack_import,  int *stack_export, 
								 HECMW_Comm repart_comm, int my_rank); 
extern int  stack_whole_send_recv(int pesize, int *stack_import,  int *stack_export, HECMW_Comm repart_comm, int my_rank);                                    
              
extern int int_part_send_recv(int n, int neibpetot, int *neibpe,int *stack_import, int *nod_import,int *stack_export, int *nod_export, 
		     int *x,  HECMW_Comm repart_comm, int my_rank);
extern int double_part_send_recv(int n, int neibpetot, int *neibpe, int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export, double *x, HECMW_Comm repart_comm, int my_rank);
extern void int_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export,int *x, int *y, HECMW_Comm repart_comm, int my_rank);
extern void double_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, int *nod_import,
		     int *stack_export, int *nod_export, double *x, double *y, HECMW_Comm repart_comm, int my_rank);

extern void int2_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, 
		     int *x, int *y, 
		     HECMW_Comm repart_comm, int my_rank);
extern void int3_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, 
		     int *x, int *y, 
		     HECMW_Comm repart_comm, int my_rank);
extern void double2_whole_send_recv(int n1, int n2, int pesize,  int *stack_import, 
		     int *stack_export, double *x, double *y, 
		     HECMW_Comm repart_comm, int my_rank);
extern void whole_copy_array(int *recv_elem_num, int *global_recv_elem_num, int mynode, int pesize, HECMW_Comm repart_comm);



void mesh_migration_adapt(int mynode, int pesize, Result_part *result, int *vtxdist)

{
    int *send_elem_num,  *send_elem, *recv_elem_num, *recv_elem, *recv_parent_num, *count_elem, *send_parent, *send_parent_num, *count_parent;
	int *flag_hit, *flag_elem_hit;
	int  ncon, id_elem;
	int i,j,k, m;
	int  *tmp_int_send, *tmp_int_recv, *tmp2_int_send, *tmp2_int_recv;
	int *new_part;
	int *tmp_stack_send, *tmp_stack_recv, *tmp_neibpe, tmp_int, tmp_sum;
	int *send_ptr_num, *send_ptr_parent_num, *recv_ptr_num, *recv_ptr_parent_num;
	int *global_index, *global_index_hit;
	int *send_node_num, *recv_node_num, *send_node, *recv_node;
	int  *tmp_node, *count_node, *tmp_send, *tmp_recv;
	int *import_index, *export_index, *count_index, *export_node, *import_node;
	int import_n_neighbor_pe, export_n_neighbor_pe, *import_neighbor_pe, *export_neighbor_pe;
/*	int *import2_index, *export2_index, *export2_node, *import2_node;
*/
	int *tmp_int_nodeid, *tmp2_int_nodeid;
	double *tmp_node_d, *recv_node_d, *tmp2_node_d, *recv2_node_d, *tmp_send_d, *tmp_recv_d, *tmp2_send_d, *tmp2_recv_d;
	int  *tmp_node_i, *recv_node_i, *tmp2_node_i, *recv2_node_i;
	 HECMW_Status stat;
	 FILE *test_fp, *test_fp2;
	 char test_file[128];
	 int min_pe, tmp_pe, local_nid;
	int  nvtxs, tmp_nvtxs, tmp_peid, tmp_lid;
	int *new_elem, *global_recv_elem_num, *global_recv_parent_num;
	int *send_adapt_parent, *recv_adapt_parent, *send_tmp, *recv_tmp, *send_tmp_num, *recv_tmp_num, *send_adapt_ptype, *count_num;
	int t_elem, l_child, new_l_child;
	int *inter_elem, *send_adapt_child, *recv_adapt_child, *send_adapt_child_num, *recv_adapt_child_num, *send_index_child, *recv_index_child;
	int *send2_adapt_child_num, *recv2_adapt_child_num;
	int  *send_inter, *recv_inter, *send_inter_num, *recv_inter_num, *global_recv_inter_num;
	int  *send_parent_inter, *recv_parent_inter, *send_parent_inter_num, *recv_parent_inter_num, *global_recv_parent_inter_num;
	int init_flag, new_elemid;
	int *tmp_elem_grp, num_grp_item, *tmp_surf_grp, *tmp_surf_id, *tmp2_surf_id;
    Tmp_grp_inf *tmp_grp;
	int  *new_nelem_dist;
	int  tn_component;
	int   *new_vtxdist;
    int  *new2old, *count_elem_index, *new_tmp, *new_tmp2, *old2new;
    double  t1, t2;

  if(mynode==0) {
      t1=HECMW_Wtime();
      }

#ifdef TEST
    sprintf(test_file, "test3.%d", mynode);
	test_fp=fopen(test_file, "w");
	if(test_fp==NULL) {
		fprintf(stderr, "Cannot open test_file\n");
		exit(0);
	}
#endif
/*	fprintf(test_fp, "%d\n", mesh->ne_internal);
	for(i=0;i<mesh->ne_internal;i++)
		fprintf(test_fp, "%d  %d\n", mesh->elem_internal_list[i], mesh->adapt_type[mesh->elem_internal_list[i]]);
	fprintf(test_fp, "%d  %d\n", mesh->n_node, mesh->nn_internal);
	
	fprintf(test_fp, "%d \n", mesh->n_neighbor_pe);
	for(i=0;i<mesh->n_neighbor_pe;i++)
		fprintf(test_fp, "%d  ", mesh->neighbor_pe[i]);
	fprintf(test_fp, "\n");
	for(i=0;i<mesh->n_neighbor_pe+1;i++)
		fprintf(test_fp, "%d ", mesh->import_index[i]);
	fprintf(test_fp, "\n");
	for(i=0;i<mesh->n_neighbor_pe+1;i++)
		fprintf(test_fp, "%d ", mesh->export_index[i]);
	fprintf(test_fp, "\n");
	for(i=0;i<mesh->import_index[mesh->n_neighbor_pe];i++)
		fprintf(test_fp, "%d\n", mesh->import_item[i]);
	fprintf(test_fp, "export node*****\n");
	for(i=0;i<mesh->export_index[mesh->n_neighbor_pe];i++)
		fprintf(test_fp, "%d\n", mesh->export_item[i]);
*/
/*
        fprintf(test_fp, "in PE %d n_elem=%d ne_internal=%d\n", mynode, mesh->n_elem, mesh->ne_internal);
        for(i=0;i<mesh->n_elem;i++) 
          fprintf(test_fp, "i=%d %d  %d  %d  %d %d %d\n", i, mesh->elem_ID[i*2], mesh->elem_ID[i*2+1], mesh->adapt_parent[i*2], mesh->adapt_parent[i*2+1], mesh->adapt_type[i], mesh->elem_node_index[i]);
	fclose(test_fp);
*/
#ifdef test
    sprintf(test_file, "test5.%d", mynode);
	test_fp2=fopen(test_file, "w");
	if(test_fp==NULL) 
		HECMW_dlb_print_exit("Cannot open test_file\n");
#endif	
	
/* 1: first update new partition to include import nodes */

	  new_part=(int *)calloc(mesh->n_node, sizeof(int));
	  if(new_part==NULL) 
	   	HECMW_dlb_memory_exit("new_part");
	  for(i=0;i<mesh->nn_internal;i++)
		new_part[i]=result->part[i];
      int_part_send_recv(mesh->n_node, mesh->n_neighbor_pe, mesh->neighbor_pe, 
		   mesh->import_index, mesh->import_item,
		   mesh->export_index, mesh->export_item,
		   new_part, mesh->HECMW_COMM, mynode);
	  free(result->part);
	  result->part=new_part;
      new_elem=(int *)calloc(mesh->n_elem, sizeof(int));
	  if(new_elem==NULL)
	    	HECMW_dlb_memory_exit("new_elem");
	  inter_elem=(int *)calloc(mesh->n_elem, sizeof(int));
/*	  tmp_int=0;
*/
      for(i=0;i<mesh->n_elem;i++) {
	     inter_elem[i]=-1;
	     if(mesh->elem_ID[i*2+1]==mynode) {
		    inter_elem[mesh->elem_ID[i*2]-1]=i;
/*		    tmp_int++;
*/
		 }
	  }
/*	  fprintf(test_fp, "ne_internal====%d\n", mesh->ne_internal);
	for(i=0;i<mesh->n_elem;i++)
		fprintf(test_fp, "%d   %d  %d   %d\n", i, inter_elem[i], mesh->elem_ID[i*2], mesh->elem_ID[i*2+1]);  
	fclose(test_fp);
*/
	   for(i=0;i<mesh->n_elem;i++)
		 new_elem[i]=-1;
	   if(pesize>1)
	     HECMW_Allreduce(&mesh->n_elem, &t_elem, 1, HECMW_INT, HECMW_SUM, mesh->HECMW_COMM);
	   else
		 t_elem=mesh->n_elem;


/* send_parent_num=(int *)calloc(pesize+1, sizeof(int)); */
	send_elem_num=(int *)calloc(pesize+1, sizeof(int));
	send_inter_num=(int *)calloc(pesize+1, sizeof(int));
	recv_elem_num=(int *)calloc(pesize+1, sizeof(int));
	recv_inter_num=(int *)calloc(pesize+1, sizeof(int));
	if((send_elem_num==NULL) || (send_inter_num==NULL) || (recv_elem_num==NULL) || (recv_inter_num==NULL))
		HECMW_dlb_memory_exit("send_elem_num ");
	  flag_hit=(int *)calloc(pesize, sizeof(int));
	  flag_elem_hit=(int *)calloc(mesh->n_elem, sizeof(int));
	  if((flag_hit==NULL) || (flag_elem_hit==NULL))
		 HECMW_dlb_memory_exit("flag_hit");

	for(i=0;i<pesize+1;i++) {
		send_elem_num[i]=0;
		send_inter_num[i]=0;
		recv_elem_num[i]=0;
		recv_inter_num[i]=0;
	}

	   for(i=0;i<mesh->n_elem;i++) {
		   for(j=0;j<pesize;j++)
		     	flag_hit[j]=0;
		   init_flag=pesize;
		if((mesh->elem_ID[i*2+1]==mynode) && (mesh->adapt_type[i]==0)) {
		   ncon=mesh->elem_node_index[i+1]-mesh->elem_node_index[i];
		      for(j=mesh->elem_node_index[i];j<mesh->elem_node_index[i+1];j++) {
			       if(flag_hit[result->part[mesh->elem_node_item[j]-1]]==0) {
				       flag_hit[result->part[mesh->elem_node_item[j]-1]]=1;
				       send_elem_num[result->part[mesh->elem_node_item[j]-1]+1]++;
					   if(result->part[mesh->elem_node_item[j]-1]<init_flag)
						   init_flag=result->part[mesh->elem_node_item[j]-1];
				   }                      
			  }
			  send_inter_num[init_flag+1]++;
		}
	   }
/*		for(i=1;i<pesize+1;i++)
			fprintf(stderr, "in PE %d the number send to PE %d is %d\n", mynode, i, send_elem_num[i]);
			*/




	for(i=1;i<pesize+1;i++) {
		send_elem_num[i]=send_elem_num[i-1]+send_elem_num[i];
		send_inter_num[i]=send_inter_num[i-1]+send_inter_num[i];
	}
/*
	for(i=0;i<pesize+1;i++) 
		fprintf(stderr, "PE %d: send %d is %d\n", mynode, i, send_elem_num[i]);
  */  
    HECMW_Barrier(mesh->HECMW_COMM);
	
	stack_whole_send_recv(pesize, send_elem_num, recv_elem_num, mesh->HECMW_COMM, mynode);
	stack_whole_send_recv(pesize, send_inter_num, recv_inter_num, mesh->HECMW_COMM, mynode);
/*
	for(i=0;i<pesize+1;i++) 
		fprintf(stderr, "PE %d: recv %d is %d\n", mynode, i, recv_elem_num[i]);
*/

	global_recv_elem_num=(int *)calloc(pesize*(pesize+1), sizeof(int));
	global_recv_inter_num=(int *)calloc(pesize*(pesize+1), sizeof(int));
	if((global_recv_elem_num==NULL) || (global_recv_inter_num==NULL))
		HECMW_dlb_memory_exit("global_recv_elem_num");
	whole_copy_array(recv_elem_num, global_recv_elem_num, mynode, pesize, mesh->HECMW_COMM);
	whole_copy_array(recv_inter_num, global_recv_inter_num, mynode, pesize, mesh->HECMW_COMM);
/*		for(i=1;i<pesize+1;i++)
			fprintf(stderr, "AFTER--in PE0 the number send to PE %d is %d\n", i, send_elem_num[i]);
			*/
	if(send_elem_num[pesize]>0) {
	send_elem=(int *)calloc(send_elem_num[pesize], sizeof(int));
	count_elem=(int *)calloc(pesize, sizeof(int));
    count_num=(int *)calloc(pesize, sizeof(int));
	send_inter=(int *)calloc(send_elem_num[pesize], sizeof(int));

	if((send_elem==NULL) || (count_elem==NULL)) 
		HECMW_dlb_memory_exit("send_elem and count_elem");
	for(i=0;i<pesize;i++) {
		count_elem[i]=0;
		count_num[i]=0;
	}
	}
	for(i=0;i<mesh->n_elem;i++) {
		for(j=0;j<pesize;j++)
			flag_hit[j]=0;
		init_flag=pesize;
		if((mesh->elem_ID[i*2+1]==mynode) && (mesh->adapt_type[i]==0)) {
		   ncon=mesh->elem_node_index[i+1]-mesh->elem_node_index[i];
		      for(j=mesh->elem_node_index[i];j<mesh->elem_node_index[i+1];j++) {
			       if(flag_hit[result->part[mesh->elem_node_item[j]-1]]==0) {
				       flag_hit[result->part[mesh->elem_node_item[j]-1]]=1;
				  send_elem[send_elem_num[result->part[mesh->elem_node_item[j]-1]]
					  +count_elem[result->part[mesh->elem_node_item[j]-1]]]=i;
				  count_elem[result->part[mesh->elem_node_item[j]-1]]++;
					   if(result->part[mesh->elem_node_item[j]-1]<init_flag)
						   init_flag=result->part[mesh->elem_node_item[j]-1];
				   }
			  }
		  new_elemid=init_flag*t_elem+global_recv_inter_num[init_flag*(pesize+1)+mynode]
			  +count_num[init_flag];
			  new_elem[i]=new_elemid;
/*			  new_elemid=init_flag*t_elem+global_recv_elem_num[init_flag*(pesize+1)+mynode]+count_elem[init_flag]-1;
			  new_elem[i]=new_elemid;
			  */
			  count_num[init_flag]++;

		}
	}
    for(i=0;i<send_elem_num[pesize];i++)
		send_inter[i]=new_elem[send_elem[i]];


/*  -------start finding parent information ----------------*/
/* find all the parents which need be sent */
	send_parent_num=(int *)calloc(pesize+1, sizeof(int));
	send_parent_inter_num=(int *)calloc(pesize+1, sizeof(int));
	recv_parent_num=(int *)calloc(pesize+1, sizeof(int));
	recv_parent_inter_num=(int *)calloc(pesize+1, sizeof(int));

	for(i=0;i<pesize+1;i++) {
       send_parent_num[i]=0;
	   send_parent_inter_num[i]=0;
	   }
		for(i=0;i<mesh->n_elem;i++)
			flag_elem_hit[i]=0;
	for(j=0;j<pesize;j++) {
		send_parent_num[j+1]=0;
		for(k=send_elem_num[j];k<send_elem_num[j+1];k++) {
			id_elem=send_elem[k];
			while((mesh->adapt_parent[id_elem*2]>0) && (mesh->adapt_parent[id_elem*2+1]==mynode) && 
				(flag_elem_hit[mesh->adapt_parent[id_elem*2]-1]==0)) {
				flag_elem_hit[mesh->adapt_parent[id_elem*2]-1]=1;
				send_parent_num[j+1]++;
				id_elem=inter_elem[mesh->adapt_parent[id_elem*2]-1];
			}
		}
	}

	for(i=1;i<pesize+1;i++) {
		send_parent_num[i]=send_parent_num[i-1]+send_parent_num[i];
	}

	recv_parent_num=(int *)calloc(pesize+1, sizeof(int));
	if(recv_parent_num==NULL)
		HECMW_dlb_memory_exit("recv_parent_num");
    stack_whole_send_recv(pesize, send_parent_num, recv_parent_num,  mesh->HECMW_COMM, mynode);  
	global_recv_parent_num=(int *)calloc(pesize*(pesize+1), sizeof(int));
	if(global_recv_parent_num==NULL)
		HECMW_dlb_memory_exit("global_recv_parent_num");
	whole_copy_array(recv_parent_num, global_recv_parent_num, mynode, pesize, mesh->HECMW_COMM);

	send_parent=(int *)calloc(send_parent_num[pesize], sizeof(int));
	send_parent_inter=(int *)calloc(send_parent_num[pesize], sizeof(int));
	count_parent=(int *)calloc(pesize, sizeof(int));

	if((send_parent==NULL) || (count_parent==NULL)) 
		HECMW_dlb_memory_exit("send_parent and count_parent");
	/* to keep each parent just occur once in all PEs */
	for(i=0;i<mesh->n_elem;i++)
		flag_elem_hit[i]=0;
	for(j=0;j<pesize;j++) {
		count_parent[j]=0;
		for(k=send_elem_num[j];k<send_elem_num[j+1];k++) {
			id_elem=send_elem[k];
			while((mesh->adapt_parent[id_elem*2]>0) && (mesh->adapt_parent[id_elem*2+1]==mynode)
				&& (flag_elem_hit[mesh->adapt_parent[id_elem*2]-1]==0)) {
				flag_elem_hit[mesh->adapt_parent[id_elem*2]-1]=1;
				send_parent[send_parent_num[j]+count_parent[j]]=inter_elem[mesh->adapt_parent[id_elem*2]-1];
		/*		send_parent_inter[send_parent_num[j]+count_parent[j]]=j*t_elem+global_recv_inter_num[j*(pesize+1)+pesize]+global_recv_parent_num[j*(pesize+1)+
					mynode]+count_parent[j];
					*/
					
				count_parent[j]++;

				new_elem[inter_elem[mesh->adapt_parent[id_elem*2]-1]]=j*t_elem+global_recv_inter_num[j*(pesize+1)+pesize]+
					global_recv_parent_num[j*(pesize+1)+mynode]+count_parent[j]-1;
			
				id_elem=inter_elem[mesh->adapt_parent[id_elem*2]-1];
			}
		}
	}
	for(i=0;i<send_parent_num[pesize];i++) {
		send_parent_inter[i]=new_elem[send_parent[i]];
		if(send_parent_inter[i]<0) 
			HECMW_dlb_print_exit("There is something wrong with send_parent_inter");
	}
	/* ------------------- start sending parent information ----------------*/

    send_tmp_num=(int *)calloc(pesize+1, sizeof(int));
	recv_tmp_num=(int *)calloc(pesize+1, sizeof(int));
	if((send_tmp_num==NULL) || (recv_tmp_num==NULL))
		HECMW_dlb_memory_exit("send_tmp_num");
	if(send_elem_num[pesize]>0) {
	send_adapt_parent=(int *)calloc(send_elem_num[pesize], sizeof(int));
	send_adapt_ptype=(int *)calloc(send_elem_num[pesize], sizeof(int));
	if((send_adapt_parent==NULL) || (send_adapt_ptype==NULL))
		HECMW_dlb_memory_exit("send_adapt_parent");
	}

	for(i=0;i<pesize+1;i++) 
		send_tmp_num[i]=0;

	for(i=0;i<send_elem_num[pesize];i++) {
		if(mesh->adapt_parent[send_elem[i]*2+1]<0)
            send_adapt_parent[i]=-1;
		else if(mesh->adapt_parent[send_elem[i]*2+1]!=mynode) 
			send_tmp_num[mesh->adapt_parent[send_elem[i]*2+1]+1]++;
		else if(new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]]!=-1)
			send_adapt_parent[i]=new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]];
		else if(new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]]==-1){
			fprintf(stderr, "There is something wrong with parent information\n");
			fprintf(stderr, "i=%d, send_elem[i]=%d parent is %d PE=%d\n", i, send_elem[i], mesh->adapt_parent[send_elem[i]*2]-1,
               mesh->adapt_parent[send_elem[i]*2]);
			HECMW_dlb_print_exit("Error in parent finding");
		}
	}
	for(i=1;i<pesize+1;i++) {
		send_tmp_num[i]=send_tmp_num[i-1]+send_tmp_num[i];
	}

	stack_whole_send_recv(pesize, send_tmp_num, recv_tmp_num, mesh->HECMW_COMM, mynode);

    for(i=0;i<pesize;i++)
		count_num[i]=0;
	if(send_tmp_num[pesize]>0) {
	send_tmp=(int *)calloc(send_tmp_num[pesize], sizeof(int));
	if(send_tmp==NULL) 
		HECMW_dlb_memory_exit("send_tmp");
	}
	if(recv_tmp_num[pesize]>0) {
	recv_tmp=(int *)calloc(recv_tmp_num[pesize], sizeof(int));
	if (recv_tmp==NULL)
		HECMW_dlb_memory_exit("send_tmp, recv_tmp");
	}
	for(i=0;i<send_elem_num[pesize];i++) {
		if(mesh->adapt_parent[send_elem[i]*2+1]<0)
            send_adapt_parent[i]=-1;
		else if(mesh->adapt_parent[send_elem[i]*2+1]!=mynode) {
            send_tmp[send_tmp_num[mesh->adapt_parent[send_elem[i]*2+1]]+
				count_num[mesh->adapt_parent[send_elem[i]*2+1]]]
				=inter_elem[mesh->adapt_parent[send_elem[i]*2]-1];
			count_num[mesh->adapt_parent[send_elem[i]*2+1]]++;
		}
	}

	if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(send_tmp_num[pesize], recv_tmp_num[pesize], pesize, recv_tmp_num, 
		send_tmp_num, send_tmp, recv_tmp, mesh->HECMW_COMM, mynode);
	HECMW_Barrier(mesh->HECMW_COMM);
    for(i=0;i<recv_tmp_num[pesize];i++) {
		if(new_elem[recv_tmp[i]]==-1) 
			HECMW_dlb_print_exit("There is something wrong in parent inf");
		else
			recv_tmp[i]=new_elem[recv_tmp[i]];
	}
    if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(recv_tmp_num[pesize], send_tmp_num[pesize], pesize, send_tmp_num, recv_tmp_num, 
		 recv_tmp, send_tmp, mesh->HECMW_COMM, mynode);

	for(i=0;i<pesize;i++)
		count_num[i]=0;
	for(i=0;i<send_elem_num[pesize];i++) {
		if(mesh->adapt_parent[send_elem[i]*2+1]<0) {
            send_adapt_parent[i]=-1;
			send_adapt_ptype[i]=-1;
		}
		else if(mesh->adapt_parent[send_elem[i]*2+1]!=mynode) {
			send_adapt_parent[i]=send_tmp[send_tmp_num[mesh->adapt_parent[send_elem[i]*2+1]]
				+count_num[mesh->adapt_parent[send_elem[i]*2+1]]];
			send_adapt_ptype[i]=mesh->adapt_parent_type[send_elem[i]];

			count_num[mesh->adapt_parent[send_elem[i]*2+1]]++;
		}
		else if(new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]]!=-1) {
			send_adapt_parent[i]=new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]];
			send_adapt_ptype[i]=mesh->adapt_parent_type[send_elem[i]];
		}
		else if(new_elem[inter_elem[mesh->adapt_parent[send_elem[i]*2]-1]]==-1)
			HECMW_dlb_print_exit("There is something wrong with parent information");

	}
    if(recv_elem_num[pesize]>0) {
	tmp_int_recv=(int *)calloc(recv_elem_num[pesize], sizeof(int));
	if(tmp_int_recv==NULL)
		HECMW_dlb_memory_exit("tmp_int_recv");
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, send_adapt_parent, tmp_int_recv, mesh->HECMW_COMM, mynode);
	}
/* new mesh whole inf. copy  */
/*        for(i=0;i<mesh->n_elem;i++) {
          if(new_elem[i]==-1)
             fprintf(stderr, "Wrong in new elem: %d %d\n", mynode, i);
             HECMW_dlb_print_exit("check");
        }
*/
	new_mesh->hecmw_n_file=mesh->hecmw_n_file;
    new_mesh->files=(char **)calloc(mesh->hecmw_n_file, sizeof(char *));
	for(i=9;i<mesh->hecmw_n_file;i++)
    new_mesh->files[i]=(char *)calloc(128, sizeof(char));

	sprintf(new_mesh->header, "%s", mesh->header);
	sprintf(new_mesh->gridfile, "%s", mesh->gridfile);
	new_mesh->hecmw_flag_adapt=mesh->hecmw_flag_adapt;
	new_mesh->hecmw_flag_initcon=mesh->hecmw_flag_initcon;
	new_mesh->hecmw_flag_parttype=mesh->hecmw_flag_parttype;
	new_mesh->hecmw_flag_partdepth=mesh->hecmw_flag_partdepth;
	new_mesh->hecmw_flag_version=mesh->hecmw_flag_version;
	new_mesh->zero_temp=mesh->zero_temp;
	new_mesh->n_dof=mesh->n_dof;
	new_mesh->n_dof_grp=mesh->n_dof_grp;

    new_mesh->n_elem=recv_elem_num[pesize]+recv_parent_num[pesize];
	new_mesh->n_elem_type=mesh->n_elem_type;
	new_mesh->n_elem_mat_ID=mesh->n_elem_mat_ID;
	if(new_mesh->n_elem<=0) 
		HECMW_dlb_print_exit("Error: New mesh: n_elem==0");
    new_mesh->adapt_parent=(int *)calloc(new_mesh->n_elem*2, sizeof(int));
	if(new_mesh->adapt_parent==NULL)
		HECMW_dlb_memory_exit("new_mesh: adaption_parent");
    new_mesh->adapt_parent_type=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_mesh->adapt_parent_type==NULL)
		HECMW_dlb_memory_exit("new_mesh: adapt_parent_type");
	for(i=0;i<recv_elem_num[pesize];i++) {
		if(tmp_int_recv[i]==-1) {
		new_mesh->adapt_parent[i*2]=0;
		new_mesh->adapt_parent[i*2+1]=-1;
		}
		else {
		new_mesh->adapt_parent[i*2]=(tmp_int_recv[i] % t_elem) +1;
		new_mesh->adapt_parent[i*2+1]=tmp_int_recv[i] / t_elem;
		}
	}
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, send_adapt_ptype, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->adapt_parent_type[i]=tmp_int_recv[i];
	if(send_tmp_num[pesize]>0)
	free(send_tmp);
	if(recv_tmp_num[pesize]>0)
	free(recv_tmp);
	if(send_elem_num[pesize]>0) {
	free(send_adapt_parent);
	free(send_adapt_ptype);
	}
	if(send_parent_num[pesize]>0) {
	send_adapt_parent=(int *)calloc(send_parent_num[pesize], sizeof(int));
	send_adapt_ptype=(int *)calloc(send_parent_num[pesize], sizeof(int));
	if((send_adapt_parent==NULL) || (send_adapt_ptype==NULL))
		HECMW_dlb_memory_exit("send_adapt_parent");
	}

	for(i=0;i<pesize+1;i++) 
		send_tmp_num[i]=0;

	for(i=0;i<send_parent_num[pesize];i++) {
		if(mesh->adapt_parent[send_parent[i]*2+1]<0)
            send_adapt_parent[i]=-1;
		else if(mesh->adapt_parent[send_parent[i]*2+1]!=mynode) 
			send_tmp_num[mesh->adapt_parent[send_parent[i]*2+1]+1]++;
		else if(new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]]!=-1)
			send_adapt_parent[i]=new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]];
		else if(new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]]==-1)
			HECMW_dlb_print_exit("There is something wrong with parent information");
	}
	for(i=1;i<pesize+1;i++) {
		send_tmp_num[i]=send_tmp_num[i-1]+send_tmp_num[i];
	}
	
	stack_whole_send_recv(pesize, send_tmp_num, recv_tmp_num, mesh->HECMW_COMM, mynode);

    for(i=0;i<pesize;i++)
		count_num[i]=0;
	if(send_tmp_num[pesize]>0) {
	send_tmp=(int *)calloc(send_tmp_num[pesize], sizeof(int));
	if(send_tmp==NULL) 
		HECMW_dlb_memory_exit("send_tmp, recv_tmp");
	}
	if(recv_tmp_num[pesize]>0) {
	recv_tmp=(int *)calloc(recv_tmp_num[pesize], sizeof(int));
	if(recv_tmp==NULL)
		HECMW_dlb_memory_exit("send_tmp, recv_tmp");
	}
	for(i=0;i<send_parent_num[pesize];i++) {
		if(mesh->adapt_parent[send_parent[i]*2+1]<0)
            send_adapt_parent[i]=-1;
		else if(mesh->adapt_parent[send_parent[i]*2+1]!=mynode) {
            send_tmp[send_tmp_num[mesh->adapt_parent[send_parent[i]*2+1]]+
				count_num[mesh->adapt_parent[send_parent[i]*2+1]]]
				=inter_elem[mesh->adapt_parent[send_parent[i]*2]-1];
			count_num[mesh->adapt_parent[send_parent[i]*2+1]]++;
		}
	}
	if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(send_tmp_num[pesize], recv_tmp_num[pesize], pesize, recv_tmp_num, 
		send_tmp_num, send_tmp, recv_tmp, mesh->HECMW_COMM, mynode);
    for(i=0;i<recv_tmp_num[pesize];i++) {
		if(new_elem[recv_tmp[i]]==-1) 
			HECMW_dlb_print_exit("There is something wrong in parent inf");
		else
			recv_tmp[i]=new_elem[recv_tmp[i]];
	}
    if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(recv_tmp_num[pesize], send_tmp_num[pesize], pesize, send_tmp_num, recv_tmp_num, 
		 recv_tmp, send_tmp, mesh->HECMW_COMM, mynode);

	for(i=0;i<pesize;i++)
		count_num[i]=0;
	for(i=0;i<send_parent_num[pesize];i++) {
		if(mesh->adapt_parent[send_parent[i]*2+1]<0) {
            send_adapt_parent[i]=-1;
			send_adapt_ptype[i]=-1;
		}
		else if(mesh->adapt_parent[send_parent[i]*2+1]!=mynode) {
			send_adapt_parent[i]=send_tmp[send_tmp_num[mesh->adapt_parent[send_parent[i]*2+1]]
				+count_num[mesh->adapt_parent[send_parent[i]*2+1]]];
			send_adapt_ptype[i]=mesh->adapt_parent_type[send_parent[i]];

			count_num[mesh->adapt_parent[send_parent[i]*2+1]]++;
		}
		else if(new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]]!=-1) {
			send_adapt_parent[i]=new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]];
			send_adapt_ptype[i]=mesh->adapt_parent_type[send_parent[i]];
		}
		else if(new_elem[inter_elem[mesh->adapt_parent[send_parent[i]*2]-1]]==-1)
			HECMW_dlb_print_exit("There is something wrong with parent information");

	}

	tmp2_int_recv=(int *)calloc(recv_parent_num[pesize], sizeof(int));
	if(tmp2_int_recv==NULL)
		HECMW_dlb_memory_exit("tmp2_int_recv");
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, send_adapt_parent, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	
	for(i=0;i<recv_parent_num[pesize];i++) {
		if(tmp2_int_recv[i]==-1) {
		new_mesh->adapt_parent[(i+recv_elem_num[pesize])*2]=0;
		new_mesh->adapt_parent[(i+recv_elem_num[pesize])*2+1]=-1;
		}
		else {

		new_mesh->adapt_parent[(i+recv_elem_num[pesize])*2]=(tmp2_int_recv[i] % t_elem) +1;
		new_mesh->adapt_parent[(i+recv_elem_num[pesize])*2+1]=tmp2_int_recv[i] / t_elem;
		}
	}
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, send_adapt_ptype, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->adapt_parent_type[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(send_tmp_num[pesize]>0)
	free(send_tmp);
	if(recv_tmp_num[pesize]>0)
	free(recv_tmp);
	if(send_parent_num[pesize]>0) {
	free(send_adapt_parent);
	free(send_adapt_ptype);
	}
    HECMW_Barrier(mesh->HECMW_COMM);




/* ------------- start migration child information ------------------  */
   send_adapt_child_num=(int *)calloc(pesize+1, sizeof(int));
   for(j=0;j<pesize+1;j++) 
	   send_adapt_child_num[j]=0;
   for(j=0;j<pesize;j++) {
	   for(i=send_elem_num[j];i<send_elem_num[j+1];i++) 
		   send_adapt_child_num[j+1]+=mesh->adapt_children_index[send_elem[i]+1]-mesh->adapt_children_index[send_elem[i]];
   }
	for(i=1;i<pesize+1;i++) {
		send_adapt_child_num[i]=send_adapt_child_num[i-1]+send_adapt_child_num[i];
	}
	   
	send_adapt_child=(int *)calloc(send_adapt_child_num[pesize], sizeof(int));
	send_index_child=(int *)calloc(send_elem_num[pesize], sizeof(int));
	if((send_adapt_child==NULL) || (send_index_child==NULL)) 
		HECMW_dlb_memory_exit("send_adapt_child, send_index_child");

	for(i=0;i<pesize+1;i++) 
		send_tmp_num[i]=0;
    tmp_int=-1;
	for(i=0;i<send_elem_num[pesize];i++) {
		for(j=mesh->adapt_children_index[send_elem[i]];j<mesh->adapt_children_index[send_elem[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0)
               send_adapt_child[tmp_int]=-1;
		    else if(mesh->adapt_children_item[j*2+1]!=mynode) 
		    	send_tmp_num[mesh->adapt_children_item[j*2+1]+1]++;
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]!=-1)
		    	send_adapt_child[tmp_int]=new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]];
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]==-1) {
			    fprintf(stderr, "There is something wrong with child information\n");
		    	fprintf(stderr, "i=%d, send_elem[i]=%d child is %d PE=%d\n", i, send_elem[i], mesh->adapt_children_item[j*2]-1,
               mesh->adapt_children_item[j*2+1]);
		    	HECMW_dlb_print_exit("Error: in finding children inf");
			}
		}
	}
	for(i=1;i<pesize+1;i++) {
		send_tmp_num[i]=send_tmp_num[i-1]+send_tmp_num[i];
	}

	stack_whole_send_recv(pesize, send_tmp_num, recv_tmp_num, mesh->HECMW_COMM, mynode);
    for(i=0;i<pesize;i++)
		count_num[i]=0;
	if(send_tmp_num[pesize]>0) {
	send_tmp=(int *)calloc(send_tmp_num[pesize], sizeof(int));
	if(send_tmp==NULL) 
		HECMW_dlb_memory_exit("send_tmp");
	}
	if(recv_tmp_num[pesize]>0) {
	recv_tmp=(int *)calloc(recv_tmp_num[pesize], sizeof(int));
	if (recv_tmp==NULL)
		HECMW_dlb_memory_exit("recv_tmp");
	}
	tmp_int=-1;
	for(i=0;i<send_elem_num[pesize];i++) {
		for(j=mesh->adapt_children_index[send_elem[i]];j<mesh->adapt_children_index[send_elem[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0)
               send_adapt_child[tmp_int]=-1;
			else if(mesh->adapt_children_item[j*2+1]!=mynode) {
            send_tmp[send_tmp_num[mesh->adapt_children_item[j*2+1]]+count_num[mesh->adapt_children_item[j*2+1]]]
				=inter_elem[mesh->adapt_children_item[j*2]-1];
			count_num[mesh->adapt_children_item[j*2+1]]++;
			}
		}
	}
	if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(send_tmp_num[pesize], recv_tmp_num[pesize], pesize, recv_tmp_num, 
		send_tmp_num, send_tmp, recv_tmp, mesh->HECMW_COMM, mynode);
    for(i=0;i<recv_tmp_num[pesize];i++) {
		if(new_elem[recv_tmp[i]]==-1) 
			HECMW_dlb_print_exit("There is something wrong in children inf");
		else
			recv_tmp[i]=new_elem[recv_tmp[i]];
	}
    if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(recv_tmp_num[pesize], send_tmp_num[pesize], pesize, send_tmp_num, recv_tmp_num, 
		 recv_tmp, send_tmp, mesh->HECMW_COMM, mynode);

	for(i=0;i<pesize;i++)
		count_num[i]=0;
	tmp_int=-1;
	for(i=0;i<send_elem_num[pesize];i++) {
		send_index_child[i]=mesh->adapt_children_index[send_elem[i]+1]-mesh->adapt_children_index[send_elem[i]];
		for(j=mesh->adapt_children_index[send_elem[i]];j<mesh->adapt_children_index[send_elem[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0) {
            send_adapt_child[tmp_int]=-1;
			}
		    else if(mesh->adapt_children_item[j*2+1]!=mynode) {
				send_adapt_child[tmp_int]= send_tmp[send_tmp_num[mesh->adapt_children_item[j*2+1]]
					+count_num[mesh->adapt_children_item[j*2+1]]];
			count_num[mesh->adapt_children_item[j*2+1]]++;
			}
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]!=-1) 
		    	send_adapt_child[tmp_int]=new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]];
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]==-1) {
			fprintf(stderr, "There is something wrong with child information\n");
			fprintf(stderr, "i=%d, send_elem[i]=%d child is %d PE=%d\n", i, send_elem[i], mesh->adapt_children_item[j*2]-1,
               mesh->adapt_children_item[j*2+1]);
			HECMW_dlb_print_exit("There is something wrong in children inf");
			}
		}
	}

	recv_adapt_child_num=(int *)calloc(pesize+1, sizeof(int));
	if(recv_adapt_child_num==NULL)
		HECMW_dlb_memory_exit("recv_adapt_child_num");

	stack_whole_send_recv(pesize, send_adapt_child_num, recv_adapt_child_num, mesh->HECMW_COMM, mynode);
	recv_adapt_child=(int *)calloc(recv_adapt_child_num[pesize], sizeof(int));
	if(recv_adapt_child==NULL)
		HECMW_dlb_memory_exit("recv_adapt_child");
	int2_whole_send_recv(send_adapt_child_num[pesize], recv_adapt_child_num[pesize], pesize, recv_adapt_child_num, 
		send_adapt_child_num, send_adapt_child, recv_adapt_child, mesh->HECMW_COMM, mynode);

	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, send_index_child, tmp_int_recv, mesh->HECMW_COMM, mynode);

/*--------first find child number in send_parent --------------- */

   send2_adapt_child_num=(int *)calloc(pesize+1, sizeof(int));
   for(j=0;j<pesize+1;j++) 
	   send2_adapt_child_num[j]=0;
   for(j=0;j<pesize;j++) {
	   for(i=send_parent_num[j];i<send_parent_num[j+1];i++) 
		   send2_adapt_child_num[j+1]+=mesh->adapt_children_index[send_parent[i]+1]-mesh->adapt_children_index[send_parent[i]];
   }
	for(i=1;i<pesize+1;i++) {
		send2_adapt_child_num[i]=send2_adapt_child_num[i-1]+send2_adapt_child_num[i];
	}
	recv2_adapt_child_num=(int *)calloc(pesize+1, sizeof(int));
	if(recv2_adapt_child_num==NULL)
		HECMW_dlb_memory_exit("recv2_adapt_child_num");

	stack_whole_send_recv(pesize, send2_adapt_child_num, recv2_adapt_child_num, mesh->HECMW_COMM, mynode);
	new_mesh->adapt_children_item=(int *)calloc(2*(recv_adapt_child_num[pesize]+recv2_adapt_child_num[pesize]), sizeof(int));
	if(new_mesh->adapt_children_item==NULL)
		HECMW_dlb_memory_exit("new_mesh: adapt_children_item");
    new_mesh->adapt_children_index=(int *)calloc(new_mesh->n_elem+1, sizeof(int));
	if(new_mesh->adapt_children_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: adapt_children_index");
/*
	new_l_child=recv_adapt_child_num[pesize]+recv2_adapt_child_num[pesize];
	*/
	
	for(i=0;i<recv_adapt_child_num[pesize];i++) {
		if(recv_adapt_child[i]==-1) {
		new_mesh->adapt_children_item[i*2]=0;
		new_mesh->adapt_children_item[i*2+1]=-1;
		}
		else {
		new_mesh->adapt_children_item[i*2]=(recv_adapt_child[i] % t_elem) +1;
		new_mesh->adapt_children_item[i*2+1]=recv_adapt_child[i] / t_elem;
		}
	}
	new_mesh->adapt_children_index[0]=0;
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->adapt_children_index[i+1]=new_mesh->adapt_children_index[i]+tmp_int_recv[i];
	if(send_tmp_num[pesize]>0)
	free(send_tmp);
	if(recv_tmp_num[pesize]>0)
	free(recv_tmp);
	if(send_elem_num[pesize]>0) {
	free(send_adapt_child);
	free(recv_adapt_child);
	free(send_index_child);
	}

	send_adapt_child=(int *)calloc(send2_adapt_child_num[pesize], sizeof(int));
	send_index_child=(int *)calloc(send_parent_num[pesize], sizeof(int));
	if((send_adapt_child==NULL) || (send_index_child==NULL)) 
		HECMW_dlb_memory_exit("send_adapt_child, send_index_child");

	for(i=0;i<pesize+1;i++) 
		send_tmp_num[i]=0;
    tmp_int=-1;
	for(i=0;i<send_parent_num[pesize];i++) {
		for(j=mesh->adapt_children_index[send_parent[i]];j<mesh->adapt_children_index[send_parent[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0)
               send_adapt_child[tmp_int]=-1;
		    else if(mesh->adapt_children_item[j*2+1]!=mynode) 
		    	send_tmp_num[mesh->adapt_children_item[j*2+1]+1]++;
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]!=-1)
		    	send_adapt_child[tmp_int]=new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]];
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]==-1) {
			    fprintf(stderr, "There is something wrong with child information\n");
		    	fprintf(stderr, "in PE %d i=%d, send_parent[i]=%d child is %d PE=%d\n", mynode, i, send_parent[i], 
					mesh->adapt_children_item[j*2]-1,mesh->adapt_children_item[j*2+1]);
		    	HECMW_dlb_print_exit("Error in finding children inf");
			}
		}
	}
	for(i=1;i<pesize+1;i++) {
		send_tmp_num[i]=send_tmp_num[i-1]+send_tmp_num[i];
	}

	stack_whole_send_recv(pesize, send_tmp_num, recv_tmp_num, mesh->HECMW_COMM, mynode);
    for(i=0;i<pesize;i++)
		count_num[i]=0;
	if(send_tmp_num[pesize]>0) {
	send_tmp=(int *)calloc(send_tmp_num[pesize], sizeof(int));
	if(send_tmp==NULL) 
		HECMW_dlb_memory_exit("send_tmp");
	}
	if(recv_tmp_num[pesize]>0) {
	recv_tmp=(int *)calloc(recv_tmp_num[pesize], sizeof(int));
	if (recv_tmp==NULL)
		HECMW_dlb_memory_exit("recv_tmp");
	}
    tmp_int=-1;
	for(i=0;i<send_parent_num[pesize];i++) {
		for(j=mesh->adapt_children_index[send_parent[i]];j<mesh->adapt_children_index[send_parent[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0)
               send_adapt_child[tmp_int]=-1;
			else if(mesh->adapt_children_item[j*2+1]!=mynode) {
            send_tmp[send_tmp_num[mesh->adapt_children_item[j*2+1]]+count_num[mesh->adapt_children_item[j*2+1]]]
				=inter_elem[mesh->adapt_children_item[j*2]-1];
			count_num[mesh->adapt_children_item[j*2+1]]++;
			}
		}
	}
	if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(send_tmp_num[pesize], recv_tmp_num[pesize], pesize, recv_tmp_num, 
		send_tmp_num, send_tmp, recv_tmp, mesh->HECMW_COMM, mynode);
    for(i=0;i<recv_tmp_num[pesize];i++) {
		if(new_elem[recv_tmp[i]]==-1) 
			HECMW_dlb_print_exit("There is something wrong in parents' children inf");
		else
			recv_tmp[i]=new_elem[recv_tmp[i]];
	}
    if(send_tmp_num[pesize]>0)
	int2_whole_send_recv(recv_tmp_num[pesize], send_tmp_num[pesize], pesize, send_tmp_num, recv_tmp_num, 
		 recv_tmp, send_tmp, mesh->HECMW_COMM, mynode);

	for(i=0;i<pesize;i++)
		count_num[i]=0;
	tmp_int=-1;
	for(i=0;i<send_parent_num[pesize];i++) {
		send_index_child[i]=mesh->adapt_children_index[send_parent[i]+1]-mesh->adapt_children_index[send_parent[i]];
		for(j=mesh->adapt_children_index[send_parent[i]];j<mesh->adapt_children_index[send_parent[i]+1];j++) {
			tmp_int++;
			if(mesh->adapt_children_item[j*2+1]<0) {
            send_adapt_child[tmp_int]=-1;
			}
		    else if(mesh->adapt_children_item[j*2+1]!=mynode) {
				send_adapt_child[tmp_int]= send_tmp[send_tmp_num[mesh->adapt_children_item[j*2+1]]
					+count_num[mesh->adapt_children_item[j*2+1]]];
			count_num[mesh->adapt_children_item[j*2+1]]++;
			}
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]!=-1) 
		    	send_adapt_child[tmp_int]=new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]];
			else if(new_elem[inter_elem[mesh->adapt_children_item[j*2]-1]]==-1) {
			fprintf(stderr, "There is something wrong with child information\n");
			fprintf(stderr, "i=%d, send_elem[i]=%d child is %d PE=%d\n", i, send_parent[i], mesh->adapt_children_item[j*2]-1,
               mesh->adapt_children_item[j*2+1]);
			HECMW_dlb_print_exit("ERROR in finding parent elements' children inf");
			}
		}
	}

	recv_adapt_child=(int *)calloc(recv2_adapt_child_num[pesize], sizeof(int));
	if(recv_adapt_child==NULL)
		HECMW_dlb_memory_exit("recv_adapt_child");
	int2_whole_send_recv(send2_adapt_child_num[pesize], recv2_adapt_child_num[pesize], pesize, recv2_adapt_child_num, 
		send2_adapt_child_num, send_adapt_child, recv_adapt_child, mesh->HECMW_COMM, mynode);

	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, send_index_child, tmp2_int_recv, mesh->HECMW_COMM, mynode);

	
	for(i=0;i<recv2_adapt_child_num[pesize];i++) {
		if(recv_adapt_child[i]==-1) {
		new_mesh->adapt_children_item[(i+recv_adapt_child_num[pesize])*2]=0;
		new_mesh->adapt_children_item[(i+recv_adapt_child_num[pesize])*2+1]=-1;
		}
		else {
		new_mesh->adapt_children_item[(i+recv_adapt_child_num[pesize])*2]=(recv_adapt_child[i] % t_elem) +1;
		new_mesh->adapt_children_item[(i+recv_adapt_child_num[pesize])*2+1]=recv_adapt_child[i] / t_elem;
		}
	}
	new_mesh->adapt_children_index[recv_elem_num[pesize]+1]=new_mesh->adapt_children_index[recv_elem_num[pesize]]
		+tmp2_int_recv[0];
	for(i=1;i<recv_parent_num[pesize];i++)
		new_mesh->adapt_children_index[i+1+recv_elem_num[pesize]]=new_mesh->adapt_children_index[i+recv_elem_num[pesize]]
		+tmp2_int_recv[i];
	if(send_tmp_num[pesize]>0)
	free(send_tmp);
	if(recv_tmp_num[pesize]>0)
	free(recv_tmp);
	if(send_parent_num[pesize]>0) {
	free(send_adapt_child);
	free(recv_adapt_child);
	free(send_index_child);
	}



/*---------------send elem_id --------------------------------*/

	if(new_mesh->n_elem>0) {
	   new_mesh->elem_ID=(int *)calloc(2*new_mesh->n_elem, sizeof(int));
	   if(new_mesh->elem_ID==NULL)
		  HECMW_dlb_memory_exit("new_mesh: elem_id");
	}
	    
	tmp_int_send=(int *)calloc(send_elem_num[pesize]+1, sizeof(int));

	if(tmp_int_send==NULL)
		HECMW_dlb_memory_exit("tmp_int_send ");
	if(send_elem_num[pesize]>0) {
    	for(i=0;i<send_elem_num[pesize];i++)
		    tmp_int_send[i]=send_inter[i];
	}
	tmp_int_recv=(int *)calloc(recv_elem_num[pesize]+1, sizeof(int));
	if(tmp_int_recv==NULL)
		HECMW_dlb_memory_exit("tmp_int_recv");

	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	if(recv_elem_num[pesize]>0) {
    	for(i=0;i<recv_elem_num[pesize];i++) {
     		new_mesh->elem_ID[i*2]=(tmp_int_recv[i] % t_elem) +1;
    		new_mesh->elem_ID[i*2+1]=tmp_int_recv[i] / t_elem;
		}
	}
	tmp2_int_send=(int *)calloc(send_parent_num[pesize], sizeof(int));
	if(tmp2_int_send==NULL)
		HECMW_dlb_memory_exit("tmp2_int_send");

	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=send_parent_inter[i];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++) {
		new_mesh->elem_ID[(i+recv_elem_num[pesize])*2]=(tmp2_int_recv[i] % t_elem) +1;
		new_mesh->elem_ID[(i+recv_elem_num[pesize])*2+1]=tmp2_int_recv[i] / t_elem;
	}



    if(mynode==0)
    fprintf(stderr, "Finish sending elem_id\n");
	if(mesh->n_elem>0) {
    	free(new_elem);
    	free(send_inter);
	    free(send_parent_inter);
	}
/*  ------  generate elem_internal_list and global_elem_id -------------- */
	new_mesh->ne_internal=recv_inter_num[pesize]+recv_parent_num[pesize];
	new_mesh->elem_internal_list=(int *)calloc(new_mesh->ne_internal, sizeof(int));
	if(new_mesh->elem_internal_list==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_internal_list");
	new_nelem_dist=(int *)calloc(pesize+1, sizeof(int));
	if(new_nelem_dist==NULL)
		HECMW_dlb_memory_exit("new_nelem_dist");
	if(mynode==0) {
		new_nelem_dist[0]=0;
		new_nelem_dist[1]=new_mesh->ne_internal;
		tmp_sum=new_mesh->ne_internal;
		for(i=1;i<pesize;i++) {
			HECMW_Recv(&tmp_int, 1, HECMW_INT, i, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
			tmp_sum+=tmp_int;
			new_nelem_dist[i+1]=tmp_sum;
		}
		for(i=1;i<pesize;i++)
	       HECMW_Send(new_nelem_dist,pesize+1,HECMW_INT, i, 0, mesh->HECMW_COMM);
	}
	else {
		HECMW_Send(&new_mesh->ne_internal, 1, HECMW_INT, 0, 0, mesh->HECMW_COMM);
        HECMW_Recv(new_nelem_dist, pesize+1, HECMW_INT, 0, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
	}
	new_mesh->global_elem_ID=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_mesh->global_elem_ID==NULL)
		HECMW_dlb_memory_exit("new_mesh: global_elem_ID");

	tmp_int=0;
	for(i=0;i<new_mesh->n_elem;i++) {
		new_mesh->global_elem_ID[i]=new_nelem_dist[new_mesh->elem_ID[i*2+1]]+new_mesh->elem_ID[i*2];
		if(new_mesh->elem_ID[i*2+1]==mynode) {
			new_mesh->elem_internal_list[tmp_int]=i+1;
			tmp_int++;
		}
	}
		



/*---------------send elem_type --------------------------------*/
	if(new_mesh->n_elem>0) {
	    new_mesh->elem_type=(int *)calloc(new_mesh->n_elem, sizeof(int));
	    if(new_mesh->elem_type==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_type");
	}
    if(send_elem_num[pesize]>0) {
	for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->elem_type[send_elem[i]];
	}
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->elem_type[i]=tmp_int_recv[i];
	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->elem_type[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->elem_type[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(mynode==0)
      fprintf(stderr, "Finish sending elem_type\n");

	if(new_mesh->n_elem>0) {
	    new_mesh->section_ID=(int *)calloc(new_mesh->n_elem, sizeof(int));
	    if(new_mesh->section_ID==NULL)
		HECMW_dlb_memory_exit("new_mesh: section_ID");
	}
    if(send_elem_num[pesize]>0) {
	for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->section_ID[send_elem[i]];
	}
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->section_ID[i]=tmp_int_recv[i];
	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->section_ID[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->section_ID[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(mynode==0)
      fprintf(stderr, "Finish sending section_ID\n");

/* ----------  send material inf. ---------------- */
    if(new_mesh->n_elem>0) {
		if(mesh->elem_mat_int_val!=NULL) {
	    new_mesh->elem_mat_int_val=(double *)calloc(new_mesh->n_elem, sizeof(double));
		if(new_mesh->elem_mat_int_val==NULL)
			HECMW_dlb_memory_exit("new_mesh: elem_mat_int_val");
	tmp_send_d=(double *)calloc(send_elem_num[pesize]+1, sizeof(double));
	tmp_recv_d=(double *)calloc(recv_elem_num[pesize]+1, sizeof(double));
	if(tmp_send_d==NULL)
		HECMW_dlb_memory_exit("tmp_send_d");
	if(send_elem_num[pesize]>0) {
	for(i=0;i<send_elem_num[pesize];i++)
		tmp_send_d[i]=mesh->elem_mat_int_val[send_elem[i]];
	}
	double2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_send_d, tmp_recv_d, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->elem_mat_int_val[i]=tmp_recv_d[i];
	free(tmp_send_d);
	free(tmp_recv_d);
	tmp2_send_d=(double *)calloc(send_parent_num[pesize]+1, sizeof(double));
	tmp2_recv_d=(double *)calloc(recv_parent_num[pesize]+1, sizeof(double));
	if(tmp2_send_d==NULL)
		HECMW_dlb_memory_exit("tmp2_send_d");
	
	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_send_d[i]=mesh->elem_mat_int_val[send_parent[i]];
	double2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_send_d, tmp2_recv_d, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->elem_mat_int_val[i+recv_elem_num[pesize]]=tmp2_recv_d[i];
	free(tmp2_send_d);
	free(tmp2_recv_d);

	}



	}
	if(new_mesh->n_elem>0) {
		new_mesh->elem_mat_ID_index=(int *)calloc(new_mesh->n_elem+1, sizeof(int));
		new_mesh->elem_mat_ID_item=(int *)calloc(new_mesh->n_elem, sizeof(int));
		if((new_mesh->elem_mat_ID_index==NULL) || (new_mesh->elem_mat_ID_item==NULL))
			HECMW_dlb_memory_exit("new_mesh: elem_mat_ID_index, elem_mat_ID_item");
			new_mesh->elem_mat_ID_index[0]=0;
		for(i=0;i<new_mesh->n_elem;i++) {
			new_mesh->elem_mat_ID_item[i]=1;
			new_mesh->elem_mat_ID_index[i+1]=i+1;
		}
	}

		
			


	if(mynode==0)
		fprintf(stderr, "Finish sending elem_material inf\n");

/*---------------send adaptation_level --------------------------------*/
	new_mesh->coarse_grid_level=mesh->coarse_grid_level;
	new_mesh->n_adapt=mesh->n_adapt;
	new_mesh->adapt_level=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_mesh->adapt_level==NULL)
		HECMW_dlb_memory_exit("new_mesh: adapt_level");

	for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->adapt_level[send_elem[i]];
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->adapt_level[i]=tmp_int_recv[i];

	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->adapt_level[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->adapt_level[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(mynode==0)
       fprintf(stderr, "Finish sending adaptation_level\n");
	/*---------------send adaptation_type --------------------------------*/
	new_mesh->adapt_type=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_mesh->adapt_type==NULL)
		HECMW_dlb_memory_exit("new_mesh: adapt_type");

	for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->adapt_type[send_elem[i]];
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->adapt_type[i]=tmp_int_recv[i];

	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->adapt_type[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->adapt_type[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(mynode==0)
      fprintf(stderr, "Finish sending adaptation_type\n");


/*  send wheniwasrefined elem -------------  */
	new_mesh->when_i_was_refined_elem=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_mesh->when_i_was_refined_elem==NULL)
		HECMW_dlb_memory_exit("new_mesh: when_i_was_refined_elem");

	for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->when_i_was_refined_elem[send_elem[i]];
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_elem_num[pesize];i++)
		new_mesh->when_i_was_refined_elem[i]=tmp_int_recv[i];

	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->when_i_was_refined_elem[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->when_i_was_refined_elem[i+recv_elem_num[pesize]]=tmp2_int_recv[i];
	if(mynode==0)
      fprintf(stderr, "Finish sending when_i_was_refined_elem\n");

	

	/*  ------- send elem_index --------  */
    if(new_mesh->n_elem>0) {
	    new_mesh->elem_node_index=(int *)calloc(new_mesh->n_elem+1, sizeof(int));
		if(new_mesh->elem_node_index==NULL)
			HECMW_dlb_memory_exit("new_mesh: elem_node_index");
	}
    for(i=0;i<send_elem_num[pesize];i++)
		tmp_int_send[i]=mesh->elem_node_index[send_elem[i]+1]-mesh->elem_node_index[send_elem[i]];
	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode); 
	for(i=0;i<send_parent_num[pesize];i++)
		tmp2_int_send[i]=mesh->elem_node_index[send_parent[i]+1]-mesh->elem_node_index[send_parent[i]];
	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);

	if(new_mesh->n_elem>0) {
	    new_mesh->elem_node_index[0]=0;
    	for(i=0;i<recv_elem_num[pesize];i++) 
		new_mesh->elem_node_index[i+1]=new_mesh->elem_node_index[i]+tmp_int_recv[i];
	for(i=0;i<recv_parent_num[pesize];i++)
		new_mesh->elem_node_index[recv_elem_num[pesize]+1+i]=tmp2_int_recv[i]+new_mesh->elem_node_index[recv_elem_num[pesize]+i];
	}

	if(mynode==0)
		fprintf(stderr, "Finish sending elem_node_index\n");
/*  -------------- send elem_node_item ------------------------------------- */
    send_ptr_num=(int *)calloc(pesize+1, sizeof(int));
	send_ptr_parent_num=(int *)calloc(pesize+1, sizeof(int));
    recv_ptr_num=(int *)calloc(pesize+1, sizeof(int));
	recv_ptr_parent_num=(int *)calloc(pesize+1, sizeof(int));
	if((send_ptr_num==NULL) || (send_ptr_parent_num==NULL) || (recv_ptr_num==NULL) || (recv_ptr_parent_num==NULL)) 
		HECMW_dlb_memory_exit("send_recv_ptr_num, send_recv_parent_ptr_num");
	for(i=0;i<pesize+1;i++) {
		send_ptr_num[i]=0;
	    send_ptr_parent_num[i]=0;
	}
	for(i=1;i<pesize+1;i++) {
		for(j=send_elem_num[i-1]; j<send_elem_num[i];j++)
			send_ptr_num[i]+=tmp_int_send[j];
	}
	for(i=1;i<pesize+1;i++)
		send_ptr_num[i]=send_ptr_num[i-1]+send_ptr_num[i];
	for(i=1;i<pesize+1;i++) {
		for(j=send_parent_num[i-1]; j<send_parent_num[i];j++)
			send_ptr_parent_num[i]+=tmp2_int_send[j];
	}
	for(i=1;i<pesize+1;i++)
		send_ptr_parent_num[i]=send_ptr_parent_num[i-1]+send_ptr_parent_num[i];
	
	stack_whole_send_recv(pesize, send_ptr_num, recv_ptr_num, mesh->HECMW_COMM, mynode);
    stack_whole_send_recv(pesize, send_ptr_parent_num, recv_ptr_parent_num, mesh->HECMW_COMM, mynode);


	   new_mesh->node_group=(struct hecmwST_node_grp *)calloc(1,sizeof(struct hecmwST_node_grp));
       new_mesh->elem_group=(struct hecmwST_elem_grp *)calloc(1,sizeof(struct hecmwST_elem_grp));
       new_mesh->surf_group=(struct hecmwST_surf_grp *)calloc(1,sizeof(struct hecmwST_surf_grp));

		new_mesh->node_group->n_grp=mesh->node_group->n_grp;
		new_mesh->elem_group->n_grp=mesh->elem_group->n_grp;
		new_mesh->surf_group->n_grp=mesh->surf_group->n_grp;
		if(mesh->node_group->n_grp>0) {
			new_mesh->node_group->grp_name=(char **)calloc(mesh->node_group->n_grp, sizeof(char*));
			for(m=0;m<mesh->node_group->n_grp;m++) {
				new_mesh->node_group->grp_name[m]=(char *)calloc(128, sizeof(char));
				if(new_mesh->node_group->grp_name[m]==NULL)
					HECMW_dlb_memory_exit("new_mesh: grp_name");
				strcpy(new_mesh->node_group->grp_name[m], mesh->node_group->grp_name[m]);
			}
		}
		if(mesh->elem_group->n_grp>0) {
			new_mesh->elem_group->grp_name=(char **)calloc(mesh->elem_group->n_grp, sizeof(char*));
			for(m=0;m<mesh->elem_group->n_grp;m++) {
				new_mesh->elem_group->grp_name[m]=(char *)calloc(128, sizeof(char));
				if(new_mesh->elem_group->grp_name[m]==NULL)
					HECMW_dlb_memory_exit("new_mesh: grp_name");
				strcpy(new_mesh->elem_group->grp_name[m], mesh->elem_group->grp_name[m]);
			}
		}
		if(mesh->surf_group->n_grp>0) {
			new_mesh->surf_group->grp_name=(char **)calloc(mesh->surf_group->n_grp, sizeof(char*));
			for(m=0;m<mesh->surf_group->n_grp;m++) {
				new_mesh->surf_group->grp_name[m]=(char *)calloc(128, sizeof(char));
				if(new_mesh->surf_group->grp_name[m]==NULL)
					HECMW_dlb_memory_exit("new_mesh: grp_name");
				strcpy(new_mesh->surf_group->grp_name[m], mesh->surf_group->grp_name[m]);
			}
		}



	if(new_mesh->elem_group->n_grp>0) {	
		tmp_grp=(Tmp_grp_inf *)calloc(new_mesh->elem_group->n_grp, sizeof(Tmp_grp_inf));
		if(tmp_grp==NULL)
			HECMW_dlb_memory_exit("tmp_grp");
/*			
		tmp_int_send=(int *)calloc(send_elem_num[pesize]+1, sizeof(int));
		if(tmp_int_send==NULL)
		    HECMW_dlb_memory_exit("tmp_int_send ");
		tmp_int_recv=(int *)calloc(recv_elem_num[pesize]+1, sizeof(int));
		if(tmp_int_recv==NULL)
		  HECMW_dlb_memory_exit("tmp_int_recv");
		tmp2_int_send=(int *)calloc(send_parent_num[pesize]+1, sizeof(int));
		if(tmp2_int_send==NULL)
		    HECMW_dlb_memory_exit("tmp2_int_send ");
		tmp2_int_recv=(int *)calloc(recv_parent_num[pesize]+1, sizeof(int));
		if(tmp2_int_recv==NULL)
		  HECMW_dlb_memory_exit("tmp2_int_recv");
		  */

		for(m=0;m<new_mesh->elem_group->n_grp;m++)
			tmp_grp[m].num_of_item=0;
    	for(m=0;m<new_mesh->elem_group->n_grp;m++) {
				if(m==0) {
        		tmp_elem_grp=(int *)calloc(mesh->n_elem, sizeof(int));
        		if(tmp_elem_grp==NULL)
        			HECMW_dlb_memory_exit("tmp_elem_grp");
				}
    			for(i=0;i<mesh->n_elem;i++)
	    			tmp_elem_grp[i]=0;
	    		if((mesh->elem_group->grp_index[m+1]-mesh->elem_group->grp_index[m])>0) {
     				for(i=mesh->elem_group->grp_index[m]; i<mesh->elem_group->grp_index[m+1];i++) 
    					tmp_elem_grp[mesh->elem_group->grp_item[i]-1]=1;
				}
            if(send_elem_num[pesize]>0) {
	           for(i=0;i<send_elem_num[pesize];i++)
		           tmp_int_send[i]=tmp_elem_grp[send_elem[i]];
        	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		          send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
			}
            if(send_parent_num[pesize]>0) {
	           for(i=0;i<send_parent_num[pesize];i++)
		           tmp2_int_send[i]=tmp_elem_grp[send_parent[i]];
        	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		          send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
			}
			num_grp_item=0;
			for(i=0;i<recv_elem_num[pesize];i++) {
				if(tmp_int_recv[i]==1)
					num_grp_item++;
			}
			for(i=0;i<recv_parent_num[pesize];i++) {
				if(tmp2_int_recv[i]==1)
					num_grp_item++;
			}
			tmp_grp[m].num_of_item=num_grp_item;
			if(num_grp_item>0) {
    			tmp_grp[m].item=(int *)calloc(num_grp_item, sizeof(int));
				if(tmp_grp[m].item==NULL)
					HECMW_dlb_memory_exit("tmp_grp:item");
				tmp_int=0;
				for(i=0;i<recv_elem_num[pesize];i++) {
					if(tmp_int_recv[i]==1) {
						tmp_grp[m].item[tmp_int]=i+1;
						tmp_int++;
					}
				}
				for(i=0;i<recv_parent_num[pesize];i++) {
					if(tmp2_int_recv[i]==1) {
						tmp_grp[m].item[tmp_int]=recv_elem_num[pesize]+i+1;
						tmp_int++;
					}
				}
			}
		}

    		free(tmp_elem_grp);
		num_grp_item=0;
		for(m=0;m<new_mesh->elem_group->n_grp;m++)
			num_grp_item+=tmp_grp[m].num_of_item;
		new_mesh->elem_group->grp_index=(int *)calloc(new_mesh->elem_group->n_grp+1, sizeof(int));
		if(new_mesh->elem_group->grp_index==NULL)
			HECMW_dlb_memory_exit("new_mesh: elem_grp: grp_index");
		new_mesh->elem_group->grp_index[0]=0;
		for(m=0;m<new_mesh->elem_group->n_grp;m++)
			new_mesh->elem_group->grp_index[m+1]=new_mesh->elem_group->grp_index[m]+tmp_grp[m].num_of_item;
		if(num_grp_item>0) {
			new_mesh->elem_group->grp_item=(int *)calloc(num_grp_item, sizeof(int));
			if(new_mesh->elem_group->grp_item==NULL)
				HECMW_dlb_memory_exit("new_mesh: elem_grp: grp_item");
			tmp_int=0;
			for(m=0;m<new_mesh->elem_group->n_grp;m++) {
				for(i=0;i<tmp_grp[m].num_of_item;i++) {
					new_mesh->elem_group->grp_item[tmp_int]=tmp_grp[m].item[i];
					tmp_int++;
				}
			}
		}
	    for(m=0;m<new_mesh->elem_group->n_grp;m++) {
			if(tmp_grp[m].num_of_item>0)
				free(tmp_grp[m].item);
		}
		free(tmp_grp);
	if(mynode==0) 
		fprintf(stderr, "Finish generating new elem_grp inf.\n");
	}


	if(new_mesh->surf_group->n_grp>0) {	
		tmp_grp=(Tmp_grp_inf *)calloc(new_mesh->surf_group->n_grp, sizeof(Tmp_grp_inf));
		if(tmp_grp==NULL)
			HECMW_dlb_memory_exit("tmp_grp");

		tmp_surf_id=(int *)calloc(recv_elem_num[pesize]+1, sizeof(int));
		if(tmp_surf_id==NULL)
		  HECMW_dlb_memory_exit("tmp_surf_id");
		tmp2_surf_id=(int *)calloc(recv_parent_num[pesize]+1, sizeof(int));
		if(tmp2_surf_id==NULL)
		  HECMW_dlb_memory_exit("tmp2_surf_id");

		for(m=0;m<new_mesh->surf_group->n_grp;m++)
			tmp_grp[m].num_of_item=0;
    	for(m=0;m<new_mesh->surf_group->n_grp;m++) {
				if(m==0) {
        		tmp_elem_grp=(int *)calloc(mesh->n_elem, sizeof(int));
        		if(tmp_elem_grp==NULL)
        			HECMW_dlb_memory_exit("tmp_elem_grp");
				tmp_surf_grp=(int *)calloc(mesh->n_elem, sizeof(int));
        		if(tmp_surf_grp==NULL)
        			HECMW_dlb_memory_exit("tmp_surf_grp");
				}
    			for(i=0;i<mesh->n_elem;i++)
	    			tmp_elem_grp[i]=0;
    			for(i=0;i<mesh->n_elem;i++)
	    			tmp_surf_grp[i]=-1;
	    		if((mesh->surf_group->grp_index[m+1]-mesh->surf_group->grp_index[m])>0) {
					for(i=mesh->surf_group->grp_index[m]; i<mesh->surf_group->grp_index[m+1];i++) { 
    					tmp_elem_grp[mesh->surf_group->grp_item[i*2]-1]=1;
    					tmp_surf_grp[mesh->surf_group->grp_item[i*2]-1]=mesh->surf_group->grp_item[i*2+1];
					}

				}
            if(send_elem_num[pesize]>0) {
	           for(i=0;i<send_elem_num[pesize];i++)
		           tmp_int_send[i]=tmp_elem_grp[send_elem[i]];
        	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		          send_elem_num, tmp_int_send, tmp_int_recv, mesh->HECMW_COMM, mynode);
	           for(i=0;i<send_elem_num[pesize];i++)
		           tmp_int_send[i]=tmp_surf_grp[send_elem[i]];
        	int2_whole_send_recv(send_elem_num[pesize], recv_elem_num[pesize], pesize, recv_elem_num, 
		          send_elem_num, tmp_int_send, tmp_surf_id, mesh->HECMW_COMM, mynode);
			}
            if(send_parent_num[pesize]>0) {
	           for(i=0;i<send_parent_num[pesize];i++)
		           tmp2_int_send[i]=tmp_elem_grp[send_parent[i]];
        	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		          send_parent_num, tmp2_int_send, tmp2_int_recv, mesh->HECMW_COMM, mynode);
	           for(i=0;i<send_parent_num[pesize];i++)
		           tmp2_int_send[i]=tmp_surf_grp[send_parent[i]];
        	int2_whole_send_recv(send_parent_num[pesize], recv_parent_num[pesize], pesize, recv_parent_num, 
		          send_parent_num, tmp2_int_send, tmp2_surf_id, mesh->HECMW_COMM, mynode);
			}
			num_grp_item=0;
			for(i=0;i<recv_elem_num[pesize];i++) {
				if(tmp_int_recv[i]==1)
					num_grp_item++;
			}
			for(i=0;i<recv_parent_num[pesize];i++) {
				if(tmp2_int_recv[i]==1)
					num_grp_item++;
			}
			tmp_grp[m].num_of_item=num_grp_item;
			if(num_grp_item>0) {
    			tmp_grp[m].item=(int *)calloc(num_grp_item*2, sizeof(int));
				if(tmp_grp[m].item==NULL)
					HECMW_dlb_memory_exit("tmp_grp:item");
				tmp_int=0;
				for(i=0;i<recv_elem_num[pesize];i++) {
					if(tmp_int_recv[i]==1) {
						tmp_grp[m].item[tmp_int*2]=i+1;
						tmp_grp[m].item[tmp_int*2+1]=tmp_surf_id[i];
						tmp_int++;
					}
				}
				for(i=0;i<recv_parent_num[pesize];i++) {
					if(tmp2_int_recv[i]==1) {
						tmp_grp[m].item[tmp_int*2]=recv_elem_num[pesize]+i+1;
						tmp_grp[m].item[tmp_int*2+1]=tmp2_surf_id[i];

						tmp_int++;
					}
				}
			}
		}

    		free(tmp_elem_grp);
			free(tmp_surf_grp);
			free(tmp2_surf_id);
			free(tmp_surf_id);
		num_grp_item=0;
		for(m=0;m<new_mesh->surf_group->n_grp;m++)
			num_grp_item+=tmp_grp[m].num_of_item;
		new_mesh->surf_group->grp_index=(int *)calloc(new_mesh->surf_group->n_grp+1, sizeof(int));
		if(new_mesh->surf_group->grp_index==NULL)
			HECMW_dlb_memory_exit("new_mesh: surf_grp: grp_index");
		new_mesh->surf_group->grp_index[0]=0;
		for(m=0;m<new_mesh->surf_group->n_grp;m++)
			new_mesh->surf_group->grp_index[m+1]=new_mesh->surf_group->grp_index[m]+tmp_grp[m].num_of_item;
		if(num_grp_item>0) {
			new_mesh->surf_group->grp_item=(int *)calloc(num_grp_item*2, sizeof(int));
			if(new_mesh->surf_group->grp_item==NULL)
				HECMW_dlb_memory_exit("new_mesh: surf_grp: grp_item");
			tmp_int=0;
			for(m=0;m<new_mesh->surf_group->n_grp;m++) {
				for(i=0;i<tmp_grp[m].num_of_item;i++) {
					new_mesh->surf_group->grp_item[tmp_int*2]=tmp_grp[m].item[i*2];
					new_mesh->surf_group->grp_item[tmp_int*2+1]=tmp_grp[m].item[i*2+1];
					tmp_int++;
				}
			}
		}
	    for(m=0;m<new_mesh->surf_group->n_grp;m++) {
			if(tmp_grp[m].num_of_item>0)
				free(tmp_grp[m].item);
		}
		free(tmp_grp);
	if(mynode==0) 
		fprintf(stderr, "Finish generating new surf_grp inf.\n");
	}




    if(tmp_int_send!=NULL)
	free(tmp_int_send);  
    if(tmp2_int_send!=NULL)
	free(tmp2_int_send);

	/*	recv_index_int=0;
	for(i=0;i<recv_elem_num[pesize];i++)
		recv_index_int+=tmp_int_recv[i];
	recv2_index_int=0;
	for(i=0;i<recv_parent_num[pesize];i++)
		recv2_index_int+=tmp2_int_recv[i];
*/		
    if(tmp_int_recv!=NULL)
	free(tmp_int_recv);
    if(tmp2_int_recv!=NULL)
	free(tmp2_int_recv);
	tmp_int_send=(int *)calloc(send_ptr_num[pesize], sizeof(int));
	tmp_int_nodeid=(int *)calloc(recv_ptr_num[pesize], sizeof(int));

	if((tmp_int_send==NULL) || (tmp_int_nodeid==NULL))
		HECMW_dlb_memory_exit("tmp_int_send, tmp_int_nodeid for ptr_elem sending");
/*	for(i=0;i<pesize+1;i++)
		fprintf(stderr, "vtxdist=== %d  ", vtxdist[i]);
	fprintf(stderr, "\n");
	*/
	tmp_int=0;
	for(i=0;i<send_elem_num[pesize];i++) {
		for(j=mesh->elem_node_index[send_elem[i]];j<mesh->elem_node_index[send_elem[i]+1];j++) {
			tmp_int_send[tmp_int]=mesh->node_ID[(mesh->elem_node_item[j]-1)*2]-1+
				vtxdist[mesh->node_ID[(mesh->elem_node_item[j]-1)*2+1]];
			tmp_int++;
		}
	}

	int2_whole_send_recv(send_ptr_num[pesize], recv_ptr_num[pesize], pesize, recv_ptr_num, 
		send_ptr_num, tmp_int_send, tmp_int_nodeid, mesh->HECMW_COMM, mynode);
	tmp2_int_send=(int *)calloc(send_ptr_parent_num[pesize], sizeof(int));
	tmp2_int_nodeid=(int *)calloc(recv_ptr_parent_num[pesize], sizeof(int));

	if((tmp2_int_send==NULL) || (tmp2_int_nodeid==NULL))
		HECMW_dlb_memory_exit("tmp_int_send, tmp_int_nodeid for ptr_elem sending");
	tmp_int=0;
	for(i=0;i<send_parent_num[pesize];i++) {
		for(j=mesh->elem_node_index[send_parent[i]];j<mesh->elem_node_index[send_parent[i]+1];j++) {
			tmp2_int_send[tmp_int]=mesh->node_ID[(mesh->elem_node_item[j]-1)*2]-1+
				vtxdist[mesh->node_ID[(mesh->elem_node_item[j]-1)*2+1]];
			tmp_int++;
		}
	}

	int2_whole_send_recv(send_ptr_parent_num[pesize], recv_ptr_parent_num[pesize], pesize, recv_ptr_parent_num, 
		send_ptr_parent_num, tmp2_int_send, tmp2_int_nodeid, mesh->HECMW_COMM, mynode);
	free(tmp_int_send);
	free(tmp2_int_send);

   HECMW_Barrier(mesh->HECMW_COMM);


/* --------- start to find import_nodes according to tmp_int_nodeid, tmp_int_peid, tmp_int_repart -------*/
	/* first building global_index according to node redistribution */
	send_node_num=(int *)calloc(pesize+1, sizeof(int));
	count_node=(int *)calloc(pesize, sizeof(int));
    send_node=(int *)calloc(mesh->nn_internal, sizeof(int));
	recv_node_num=(int *)calloc(pesize+1, sizeof(int));
	if((send_node_num==NULL) || (count_node==NULL) || (send_node==NULL) || (recv_node_num==NULL))
		HECMW_dlb_memory_exit("send_node_num, count_node, send_node, and recv_node_num");
	for(i=0;i<pesize+1;i++) {
		send_node_num[i]=0;
	}
    for(i=0;i<mesh->nn_internal;i++) { 
		send_node_num[result->part[i]+1]++;
	}
	for(i=1;i<pesize+1;i++) 
		send_node_num[i]=send_node_num[i-1]+send_node_num[i];
	for(i=0;i<pesize;i++) 
		count_node[i]=0;
    for(i=0;i<mesh->nn_internal;i++) {
		send_node[send_node_num[result->part[i]]+count_node[result->part[i]]]=i;
		count_node[result->part[i]]++;
	}
	tmp_node=(int *)calloc(send_node_num[pesize]+1, sizeof(int));
	for(i=0;i<send_node_num[pesize];i++)
		tmp_node[i]=vtxdist[mynode]+send_node[i];

	stack_whole_send_recv(pesize, send_node_num, recv_node_num, mesh->HECMW_COMM, mynode);
	recv_node=(int *)calloc(recv_node_num[pesize]+1, sizeof(int));
	if(recv_node==NULL)
		HECMW_dlb_memory_exit("recv_elem");
	int2_whole_send_recv(send_node_num[pesize], recv_node_num[pesize], pesize, recv_node_num, 
		send_node_num, tmp_node, recv_node, mesh->HECMW_COMM, mynode);
    global_index=(int *)calloc(result->t_node, sizeof(int));
    global_index_hit=(int *)calloc(result->t_node, sizeof(int));
    if((global_index==NULL) || (global_index_hit==NULL))
    	HECMW_dlb_memory_exit("global_index");
    for(i=0;i<result->t_node;i++) {
     	global_index[i]=-1;
    	global_index_hit[i]=-1;
	}
    for(i=0;i<recv_node_num[pesize];i++) {
		if(recv_node[i]>=result->t_node)
		   HECMW_dlb_print_exit("!!!! wrong in recv_node");
    	global_index[recv_node[i]]=mynode*result->t_node+i;
    	global_index_hit[recv_node[i]]=i;
	}
	free(tmp_node);


/*
    if(mynode>=0) {
	*/
    	if(mynode==0) {
     		tmp_recv=(int *)calloc(result->t_node, sizeof(int));
    		if(tmp_recv==NULL) 
     			HECMW_dlb_memory_exit("tmp_recv");
    		for(i=1;i<pesize;i++) {
    			HECMW_Recv(tmp_recv, result->t_node, HECMW_INT, i, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
    			for(j=0;j<result->t_node;j++) {
    				if(tmp_recv[j]>=0)
	    				global_index[j]=tmp_recv[j];
				}
			}
            free(tmp_recv);
    		for(i=1;i<pesize;i++)
    	       HECMW_Send(global_index,result->t_node,HECMW_INT, i, 0, mesh->HECMW_COMM);
		}
    	else {
    		HECMW_Send(global_index, result->t_node, HECMW_INT, 0, 0, mesh->HECMW_COMM);
            HECMW_Recv(global_index, result->t_node, HECMW_INT, 0, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
		}


/*   for(i=0;i<result->t_node;i++)
	   fprintf(test_fp, "%d\n", global_index[i]);
   fclose(test_fp);
   */
      import_index=(int *)calloc(pesize+1, sizeof(int));
      if(import_index==NULL)
    	   HECMW_dlb_memory_exit("import_index");
      for(i=0;i<pesize+1;i++)
    	   import_index[i]=0;
      for(i=0;i<recv_ptr_num[pesize];i++) {
    	   tmp_int=tmp_int_nodeid[i];
	       if(tmp_int>=result->t_node) {
		       fprintf(stderr, "There is somethign wrong with data: i=%d tmp_int=%d\n", i, tmp_int);
		       HECMW_dlb_print_exit("Please check again");
		   }
    	   if(global_index_hit[tmp_int]==-1) {
    		   global_index_hit[tmp_int]=-2;
    		   m=global_index[tmp_int]/result->t_node;
    		   import_index[m+1]++;
		   }
	  }
      for(i=0;i<recv_ptr_parent_num[pesize];i++) {
    	   tmp_int=tmp2_int_nodeid[i];
	       if(tmp_int>=result->t_node) {
		       fprintf(stderr, "There is somethign wrong with data: i=%d tmp_int=%d\n", i, tmp_int);
		       HECMW_dlb_print_exit("Please check again");
		   }
    	   if(global_index_hit[tmp_int]==-1) {
    		   global_index_hit[tmp_int]=-2;
    		   m=global_index[tmp_int]/result->t_node;
    		   import_index[m+1]++;
		   }
	  }
      for(i=1;i<pesize+1;i++)
	      import_index[i]=import_index[i-1]+import_index[i];
      count_index=(int *)calloc(pesize, sizeof(int));
      for(i=0;i<pesize;i++)
	      count_index[i]=0;
	  for(i=0;i<result->t_node;i++) {
	    	global_index_hit[i]=-1;
	  }
	  for(i=0;i<recv_node_num[pesize];i++) {
		global_index_hit[recv_node[i]]=i;
	  }
    	free(recv_node);
/*  ********origianl mesh free
	free(mesh->node_id);


*/
        new_mesh->n_node=recv_node_num[pesize]+import_index[pesize];
		new_mesh->node_dof_index=(int *)calloc(new_mesh->n_dof_grp+1, sizeof(int));
		if(new_mesh->node_dof_index==NULL)
			HECMW_dlb_memory_exit("new_mesh: node_dof_index");
		new_mesh->node_dof_index[0]=0;
		new_mesh->node_dof_index[1]=new_mesh->n_node;
		new_mesh->node_dof_item=(int *)calloc(new_mesh->n_dof_grp, sizeof(int));
		if(new_mesh->node_dof_item==NULL)
			HECMW_dlb_memory_exit("new_mesh: node_dof_item");
		for(i=0;i<new_mesh->n_dof_grp;i++)
			new_mesh->node_dof_item[i]=mesh->node_dof_item[i];


    	new_mesh->nn_internal=recv_node_num[pesize];
		new_mesh->n_dof=mesh->n_dof;
		new_mesh->n_dof_grp=mesh->n_dof_grp;
    	new_mesh->node_ID=(int *)calloc(2*new_mesh->n_node, sizeof(int));
    	if(new_mesh->node_ID==NULL)
    		HECMW_dlb_memory_exit("new_mesh: node_ID");
     	for(i=0;i<recv_node_num[pesize];i++) {
    		new_mesh->node_ID[i*2+1]=mynode;
    		new_mesh->node_ID[i*2]=i+1;
		}
		if((recv_ptr_num[pesize]+recv_ptr_parent_num[pesize])>0) {
        	new_mesh->elem_node_item=(int *)calloc(recv_ptr_num[pesize]+recv_ptr_parent_num[pesize], sizeof(int));

        	if(new_mesh->elem_node_item==NULL)
        		HECMW_dlb_memory_exit("new_mesh: elem_node_item");
            for(i=0;i<recv_ptr_num[pesize];i++) {
        	   tmp_int=tmp_int_nodeid[i];
        	   if(global_index_hit[tmp_int]!=-1) {
        		   new_mesh->elem_node_item[i]=global_index_hit[tmp_int];
			   }
               else {
        		   m=global_index[tmp_int]/result->t_node;
        		   global_index_hit[tmp_int]=new_mesh->nn_internal+import_index[m]+count_index[m];
                   new_mesh->elem_node_item[i]=global_index_hit[tmp_int];
         		   new_mesh->node_ID[global_index_hit[tmp_int]*2]=(global_index[tmp_int] % result->t_node)+1;
        		   new_mesh->node_ID[global_index_hit[tmp_int]*2+1]=global_index[tmp_int] / result->t_node;
         		   count_index[m]++;

			   }
			}
            for(i=0;i<recv_ptr_parent_num[pesize];i++) {
        	   tmp_int=tmp2_int_nodeid[i];
        	   if(global_index_hit[tmp_int]!=-1) {
        		   new_mesh->elem_node_item[i+recv_ptr_num[pesize]]=global_index_hit[tmp_int];
			   }
               else {
        		   m=global_index[tmp_int]/result->t_node;
        		   global_index_hit[tmp_int]=new_mesh->nn_internal+import_index[m]+count_index[m];
                   new_mesh->elem_node_item[i+recv_ptr_num[pesize]]=global_index_hit[tmp_int];
         		   new_mesh->node_ID[global_index_hit[tmp_int]*2]=(global_index[tmp_int] % result->t_node)+1;
        		   new_mesh->node_ID[global_index_hit[tmp_int]*2+1]=global_index[tmp_int] / result->t_node;
         		   count_index[m]++;

			   }
			}

		}
	    for(i=0;i<recv_ptr_num[pesize]+recv_ptr_parent_num[pesize];i++)
           	   new_mesh->elem_node_item[i]+=1;


   	  new_mesh->elem_internal_list=(int *)calloc(new_mesh->ne_internal, sizeof(int));
	  if(new_mesh->elem_internal_list==NULL) 
		  HECMW_dlb_memory_exit("new_mesh: elem_internal_list");
    	for(i=0;i<new_mesh->ne_internal;i++) 
	    	new_mesh->elem_internal_list[i]=0;
    	new_mesh->ne_internal=0;
/*	for(i=0;i<recv_elem_num[pesize];i++) {
			min_pe=mynode;
			for(j=new_smesh->index_elem[i];j<new_smesh->index_elem[i+1];j++) {
				local_nid=new_smesh->ptr_elem[j]-1;
				tmp_pe=new_smesh->node_id[new_smesh->n_node+local_nid];
				if(tmp_pe<min_pe)
					min_pe=tmp_pe;
			}
			if(min_pe==mynode) {
				new_smesh->ne_internal++;
				new_smesh->ne_internal_list[i]=1;
			}
		}
		*/
		
	for(i=0;i<new_mesh->n_elem;i++) {
		if(new_mesh->elem_ID[i*2+1]==mynode) {
				new_mesh->elem_internal_list[new_mesh->ne_internal]=i+1;
            	new_mesh->ne_internal++;
		}
	}
/*
     fprintf(stderr, "In pe %d ne_internal=%d  n_elem=%d t_elem=%d\n", mynode, new_mesh->ne_internal, new_mesh->n_elem,
		 t_elem);
*/
	 tmp_int=0;
     for(i=1;i<pesize+1;i++) {
	    if((import_index[i]-import_index[i-1])>0)
		   tmp_int++;
     }
     import_n_neighbor_pe=tmp_int;
     import_neighbor_pe=(int *)calloc(tmp_int, sizeof(int));

      tmp_int=-1;
      for(i=1;i<pesize+1;i++) {
    	   if((import_index[i]-import_index[i-1])>0) {
     		   tmp_int++;
    	       import_neighbor_pe[tmp_int]=i-1;
		   }
	  }
      export_index=(int *)calloc(pesize+1, sizeof(int));
      stack_whole_send_recv(pesize, import_index, export_index, mesh->HECMW_COMM, mynode);
      export_node=(int *)calloc(export_index[pesize]+1, sizeof(int));
    	tmp_send=(int *)calloc(import_index[pesize]+1, sizeof(int));
	if(tmp_send==NULL)
		HECMW_dlb_memory_exit("tmp_send");
	for(i=0;i<import_index[pesize];i++)
		tmp_send[i]=new_mesh->node_ID[(new_mesh->nn_internal+i)*2];
/*
	for(i=0;i<vis_pesize+1;i++)
		fprintf(test_fp, "%d  ", import_index[i]);
	fprintf(test_fp, "\n");
    for(i=0;i<import_index[vis_pesize];i++)
		fprintf(test_fp, "%d  %d\n", tmp_send[i], new_smesh->node_ID[(new_smesh->nn_internal+i)*2+1]);
	fclose(test_fp);
*/	

	int2_whole_send_recv(import_index[pesize], export_index[pesize], pesize, export_index, 
		import_index, tmp_send, export_node, mesh->HECMW_COMM, mynode);

/*	for(i=0;i<pesize+1;i++)
		fprintf(test_fp2, "%d  ", export_index[i]);
	fprintf(test_fp2, "\n");
    for(i=0;i<export_index[pesize];i++)
		fprintf(test_fp2, "%d \n", export_node[i]);
		*/


	 tmp_int=0;
     for(i=1;i<pesize+1;i++) {
	   if((export_index[i]-export_index[i-1])>0) 
		   tmp_int++;
	 }
     export_n_neighbor_pe=tmp_int;
     export_neighbor_pe=(int *)calloc(tmp_int, sizeof(int));
     tmp_int=-1;
     for(i=1;i<pesize+1;i++) {
	   if((export_index[i]-export_index[i-1])>0) {
		   tmp_int++;
	       export_neighbor_pe[tmp_int]=i-1;
	   }
	 }

     if(export_n_neighbor_pe>import_n_neighbor_pe) {
  	     new_mesh->n_neighbor_pe=export_n_neighbor_pe;
		 new_mesh->neighbor_pe=(int *)calloc(new_mesh->n_neighbor_pe, sizeof(int));
		 for(i=0;i<new_mesh->n_neighbor_pe;i++)
    	 new_mesh->neighbor_pe[i]=export_neighbor_pe[i];
	 }
     else {
    	   new_mesh->n_neighbor_pe=import_n_neighbor_pe;
	       new_mesh->neighbor_pe=(int *)calloc(new_mesh->n_neighbor_pe, sizeof(int));
		   for(i=0;i<new_mesh->n_neighbor_pe;i++)
    	     new_mesh->neighbor_pe[i]=import_neighbor_pe[i];
	 }
    /* free(import_neighbor_pe);
	 free(export_neighbor_pe);
*/
	 new_mesh->export_index=(int *)calloc(new_mesh->n_neighbor_pe+1, sizeof(int));
     new_mesh->import_index=(int *)calloc(new_mesh->n_neighbor_pe+1, sizeof(int));
	 new_mesh->export_index[0]=0;
     new_mesh->import_index[0]=0;
     for(i=0;i<new_mesh->n_neighbor_pe;i++) {
	       new_mesh->export_index[i+1]=export_index[new_mesh->neighbor_pe[i]+1];
	       new_mesh->import_index[i+1]=import_index[new_mesh->neighbor_pe[i]+1];
	 }
	 if(import_index[pesize]>0) {
     new_mesh->import_item=(int *)calloc(import_index[pesize], sizeof(int));
      if(new_mesh->import_item==NULL)
	     HECMW_dlb_memory_exit("new_mesh: import_item");
	 }
     for(i=0;i<import_index[pesize];i++)
	     new_mesh->import_item[i]=new_mesh->nn_internal+i+1;
	 if(export_index[pesize]>0) {
     new_mesh->export_item=(int *)calloc(export_index[pesize], sizeof(int));
	 if(new_mesh->export_item==NULL)
		 HECMW_dlb_memory_exit("new_mesh: export_item");
	 }
     for(i=0;i<export_index[pesize];i++)
	    new_mesh->export_item[i]=export_node[i];
     free(global_index_hit);

     HECMW_Barrier(mesh->HECMW_COMM);



/*
	} 
	 *//* end of if(mynode in VIS_COMM) */

/* ----------------sending node information --------------------- */

	 if(new_mesh->n_node>0) {
   new_mesh->node=(double *)calloc(3*new_mesh->n_node, sizeof(double));
   if(new_mesh->node==NULL)
	   HECMW_dlb_memory_exit("new_mesh: node");
	 }

    

	tmp_node_d=(double *)calloc(send_node_num[pesize]+1, sizeof(double));
	recv_node_d=(double *)calloc(recv_node_num[pesize]+1, sizeof(double));
	tmp2_node_d=(double *)calloc(export_index[pesize]+1, sizeof(double));
	recv2_node_d=(double *)calloc(import_index[pesize]+1, sizeof(double));
	if((tmp_node_d==NULL) || (recv_node_d==NULL) || (tmp2_node_d==NULL) || (recv2_node_d==NULL))
		HECMW_dlb_memory_exit("new_mesh: recv_node");
	for(j=0;j<3;j++) {
	    for(i=0;i<recv_node_num[pesize];i++)
		   recv_node_d[i]=0.0;
	    for(i=0;i<send_node_num[pesize];i++)
		   tmp_node_d[i]=mesh->node[send_node[i]*3+j];
	    double2_whole_send_recv(send_node_num[pesize], recv_node_num[pesize], pesize, recv_node_num, 
		  send_node_num, tmp_node_d, recv_node_d, mesh->HECMW_COMM, mynode);
     	for(i=0;i<export_index[pesize];i++)
     	   tmp2_node_d[i]=recv_node_d[export_node[i]-1];
	
        double2_whole_send_recv(export_index[pesize], import_index[pesize], pesize, 
		      import_index, export_index, tmp2_node_d, recv2_node_d, mesh->HECMW_COMM, mynode);
        for(i=0;i<new_mesh->nn_internal;i++)
     	    new_mesh->node[i*3+j]=recv_node_d[i];
     	for(i=0;i<import_index[pesize];i++)
	        new_mesh->node[(i+new_mesh->nn_internal)*3+j]=recv2_node_d[i];
		}
		

/*	
	fprintf(stderr, "n_node=%d nn_internal=%d recv_node_num=%d import_index_num=%d\n", new_mesh->n_node, new_mesh->nn_internal,
		recv_node_num[pesize], new_mesh->import_index[new_mesh->n_neighbor_pe]);
*/
/*
	global_comm_table->send_node_num=send_node_num;
	global_comm_table->recv_node_num=recv_node_num;
	
	global_comm_table->send_node=send_node;
	if(mynode>=mesh_pesize) {
	global_comm_table->import_index=import_index;
	global_comm_table->export_index=export_index;
	global_comm_table->export_node=export_node;
	}
	*/
/*    free(node->data);
*/
	free(tmp_int_nodeid);
	free(tmp2_int_nodeid);

/*  --------------------start sending result data ----------------------   */
	new_data->nn_component=data->nn_component;
	new_data->ne_component=data->ne_component;
	if(new_data->nn_component>0) {
	new_data->nn_dof=(int *)calloc(new_data->nn_component, sizeof(int));
	new_data->node_label=(char **)calloc(new_data->nn_component, sizeof(char *));
	for(i=0;i<new_data->nn_component;i++)
		new_data->node_label[i]=(char *)calloc(128, sizeof(char));
	if(new_data->nn_dof==NULL)
		HECMW_dlb_memory_exit("new_data: nn_dof");
		for(i=0;i<new_data->nn_component;i++) {
			new_data->nn_dof[i]=data->nn_dof[i];
			strcpy(new_data->node_label[i], data->node_label[i]);
		}
	}
	if(new_data->ne_component>0) {
	new_data->ne_dof=(int *)calloc(new_data->ne_component, sizeof(int));
	new_data->elem_label=(char **)calloc(new_data->ne_component, sizeof(char *));
	for(i=0;i<new_data->ne_component;i++)
		new_data->elem_label[i]=(char *)calloc(128, sizeof(char));
	if(new_data->ne_dof==NULL)
		HECMW_dlb_memory_exit("new_data: ne_dof");
		for(i=0;i<new_data->ne_component;i++) {
			new_data->ne_dof[i]=data->ne_dof[i];
			strcpy(new_data->elem_label[i], data->elem_label[i]);
		}
	}
	if(new_data->nn_component>0) {
		tn_component=0;
		for(i=0;i<new_data->nn_component;i++)
			tn_component+=new_data->nn_dof[i];
	new_data->node_val_item=(double *)calloc(tn_component*new_mesh->n_node, sizeof(double));
	if(new_data->node_val_item==NULL)
		HECMW_dlb_memory_exit("new_data: node_val_item");
	for(j=0;j<tn_component;j++) {
	for(i=0;i<send_node_num[pesize];i++)
		tmp_node_d[i]=data->node_val_item[send_node[i]*tn_component+j];

	double2_whole_send_recv(send_node_num[pesize], recv_node_num[pesize], pesize, recv_node_num, 
		send_node_num, tmp_node_d, recv_node_d, mesh->HECMW_COMM, mynode);
	for(i=0;i<export_index[pesize];i++)
		tmp2_node_d[i]=recv_node_d[export_node[i]-1];
	
    double2_whole_send_recv(export_index[pesize], import_index[pesize], pesize, 
		   import_index, export_index, tmp2_node_d, 
		   recv2_node_d, mesh->HECMW_COMM, mynode);
	for(i=0;i<new_mesh->nn_internal;i++)
		new_data->node_val_item[i*tn_component+j]=recv_node_d[i];
	for(i=0;i<import_index[pesize];i++)
		new_data->node_val_item[(i+new_mesh->nn_internal)*tn_component+j]=recv2_node_d[i];
		   

	}
	}


	free(tmp_node_d);
	free(recv_node_d);
	free(tmp2_node_d);
	free(recv2_node_d);
	new_mesh->when_i_was_refined_node=(int *)calloc(new_mesh->n_node, sizeof(int));
	if(new_mesh->when_i_was_refined_node==NULL)
		HECMW_dlb_memory_exit("new_mesh: when_i_was_refined_node");
	tmp_node_i=(int *)calloc(send_node_num[pesize]+1, sizeof(int));
	recv_node_i=(int *)calloc(recv_node_num[pesize]+1, sizeof(int));
	tmp2_node_i=(int *)calloc(export_index[pesize]+1, sizeof(int));
	recv2_node_i=(int *)calloc(import_index[pesize]+1, sizeof(int));
	if((tmp_node_i==NULL) || (recv_node_i==NULL))
		HECMW_dlb_memory_exit("new_mesh: when_i_was_refined_node");
	for(i=0;i<recv_node_num[pesize];i++)
		recv_node_i[i]=-1;
	for(i=0;i<send_node_num[pesize];i++)
		tmp_node_i[i]=mesh->when_i_was_refined_node[send_node[i]];

	int2_whole_send_recv(send_node_num[pesize], recv_node_num[pesize], pesize, recv_node_num, 
		send_node_num, tmp_node_i, recv_node_i,  mesh->HECMW_COMM, mynode);
	for(i=0;i<export_index[pesize];i++)
		tmp2_node_i[i]=recv_node_i[export_node[i]-1];
	
    int2_whole_send_recv(export_index[pesize], import_index[pesize], pesize, 
		   import_index, export_index, tmp2_node_i, 
		   recv2_node_i,  mesh->HECMW_COMM, mynode);
	for(i=0;i<new_mesh->nn_internal;i++)
		new_mesh->when_i_was_refined_node[i]=recv_node_i[i];
	for(i=0;i<import_index[pesize];i++)
		new_mesh->when_i_was_refined_node[i+new_mesh->nn_internal]=recv2_node_i[i];

	if(new_mesh->node_group->n_grp>0) {	
		tmp_grp=(Tmp_grp_inf *)calloc(new_mesh->node_group->n_grp, sizeof(Tmp_grp_inf));
		if(tmp_grp==NULL)
			HECMW_dlb_memory_exit("tmp_grp");

		for(m=0;m<new_mesh->node_group->n_grp;m++)
			tmp_grp[m].num_of_item=0;
    	for(m=0;m<new_mesh->node_group->n_grp;m++) {
				if(m==0) {
        		tmp_elem_grp=(int *)calloc(mesh->nn_internal, sizeof(int));
        		if(tmp_elem_grp==NULL)
        			HECMW_dlb_memory_exit("tmp_elem_grp");
				}
    			for(i=0;i<mesh->nn_internal;i++)
	    			tmp_elem_grp[i]=0;
	    		if((mesh->node_group->grp_index[m+1]-mesh->node_group->grp_index[m])>0) {
					for(i=mesh->node_group->grp_index[m]; i<mesh->node_group->grp_index[m+1];i++) {
						if(mesh->node_group->grp_item[i]<=mesh->nn_internal)
    					tmp_elem_grp[mesh->node_group->grp_item[i]-1]=1;
					}
				}
	for(i=0;i<recv_node_num[pesize];i++)
		recv_node_i[i]=-1;
	for(i=0;i<send_node_num[pesize];i++)
		tmp_node_i[i]=tmp_elem_grp[send_node[i]];

	int2_whole_send_recv(send_node_num[pesize], recv_node_num[pesize], pesize, recv_node_num, 
		send_node_num, tmp_node_i, recv_node_i,  mesh->HECMW_COMM, mynode);
	for(i=0;i<export_index[pesize];i++)
		tmp2_node_i[i]=recv_node_i[export_node[i]-1];
	
    int2_whole_send_recv(export_index[pesize], import_index[pesize], pesize, 
		   import_index, export_index, tmp2_node_i, 
		   recv2_node_i,  mesh->HECMW_COMM, mynode);
	num_grp_item=0;
	for(i=0;i<new_mesh->nn_internal;i++) {
		if(recv_node_i[i]==1)
			num_grp_item++;
	}
	for(i=0;i<import_index[pesize];i++) {
		if(recv2_node_i[i]==1)
			num_grp_item++;
	}
			tmp_grp[m].num_of_item=num_grp_item;
			if(num_grp_item>0) {
    			tmp_grp[m].item=(int *)calloc(num_grp_item, sizeof(int));
				if(tmp_grp[m].item==NULL)
					HECMW_dlb_memory_exit("tmp_grp:item");
				tmp_int=0;
				for(i=0;i<new_mesh->nn_internal;i++) {
					if(recv_node_i[i]==1) {
						tmp_grp[m].item[tmp_int]=i+1;
						tmp_int++;
					}
				}
				for(i=0;i<import_index[pesize];i++) {
		          if(recv2_node_i[i]==1) {
						tmp_grp[m].item[tmp_int]=new_mesh->nn_internal+i+1;
						tmp_int++;
					}
				}
			}
		}

    		free(tmp_elem_grp);
		num_grp_item=0;
		for(m=0;m<new_mesh->node_group->n_grp;m++)
			num_grp_item+=tmp_grp[m].num_of_item;
		new_mesh->node_group->grp_index=(int *)calloc(new_mesh->node_group->n_grp+1, sizeof(int));
		if(new_mesh->node_group->grp_index==NULL)
			HECMW_dlb_memory_exit("new_mesh: node_grp: grp_index");
		new_mesh->node_group->grp_index[0]=0;
		for(m=0;m<new_mesh->node_group->n_grp;m++)
			new_mesh->node_group->grp_index[m+1]=new_mesh->node_group->grp_index[m]+tmp_grp[m].num_of_item;
		if(num_grp_item>0) {
			new_mesh->node_group->grp_item=(int *)calloc(num_grp_item, sizeof(int));
			if(new_mesh->node_group->grp_item==NULL)
				HECMW_dlb_memory_exit("new_mesh: node_grp: grp_item");
			tmp_int=0;
			for(m=0;m<new_mesh->node_group->n_grp;m++) {
				for(i=0;i<tmp_grp[m].num_of_item;i++) {
					new_mesh->node_group->grp_item[tmp_int]=tmp_grp[m].item[i];
					tmp_int++;
				}
			}
		}
	    for(m=0;m<new_mesh->node_group->n_grp;m++) {
			if(tmp_grp[m].num_of_item>0)
				free(tmp_grp[m].item);
		}
		free(tmp_grp);
	if(mynode==0) 
		fprintf(stderr, "Finish generating new node_grp inf.\n");
	}




	free(tmp_node_i);
	free(tmp2_node_i);
	free(recv_node_i);
	free(recv2_node_i);

	nvtxs=new_mesh->nn_internal;
    new_vtxdist=(int *)calloc(pesize+1, sizeof(int));
	
	if(mynode==0) {
		new_vtxdist[0]=0;
		new_vtxdist[1]=nvtxs;
		tmp_sum=nvtxs;
		for(i=1;i<pesize;i++) {
			HECMW_Recv(&tmp_nvtxs, 1, HECMW_INT, i, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
			tmp_sum+=tmp_nvtxs;
			new_vtxdist[i+1]=tmp_sum;
		}
		for(i=1;i<pesize;i++)
	       HECMW_Send(new_vtxdist,pesize+1,HECMW_INT, i, 0, mesh->HECMW_COMM);
	}
	else {
		HECMW_Send(&nvtxs, 1, HECMW_INT, 0, 0, mesh->HECMW_COMM);
        HECMW_Recv(new_vtxdist, pesize+1, HECMW_INT, 0, HECMW_ANY_TAG, mesh->HECMW_COMM, &stat);
	}
	   new_mesh->global_node_ID=(int *)calloc(new_mesh->n_node, sizeof(int));
	   if(new_mesh->global_node_ID==NULL)
		   HECMW_dlb_print_exit("new_mesh: global_node_ID");
	   for(i=0;i<new_mesh->n_node;i++)
		   new_mesh->global_node_ID[i]=new_vtxdist[new_mesh->node_ID[i*2+1]]+new_mesh->node_ID[i*2];
	   free(new_vtxdist);
/*
   for(i=0;i<result->t_node;i++)
	   global_new2old[i]=-1;

	for(i=0; i<result->t_node;i++) {
		tmp_peid=global_index[i] / result->t_node;
		tmp_lid=global_index[i] % result->t_node;
		global_new2old[new_vtxdist[tmp_peid]+tmp_lid]=i;
	}
	*/
	free(global_index);
	/* set new data structure */
/*	if(mynode>=mesh_pesize) {
		new_data->node_val_item=(double *)calloc(new_smesh->n_node, sizeof(double));
		if(new_data->node_val_item==NULL)
			HECMW_dlb_memory_exit("new_data: node_val_item");
	}
  */ 
	
    if(new_mesh->n_elem_type>1) {
	new2old=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new2old==NULL)
		HECMW_dlb_memory_exit("new2old");
    new_mesh->elem_type_index=(int *)calloc(new_mesh->n_elem_type+1, sizeof(int));
	if(new_mesh->elem_type_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_type_index");
	for(i=0;i<new_mesh->n_elem_type+1;i++)
		new_mesh->elem_type_index[i]=0;
	new_mesh->elem_type_item=(int *)calloc(new_mesh->n_elem_type, sizeof(int));
	if(new_mesh->elem_type_item==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_type_item");
	for(i=0;i<new_mesh->n_elem_type;i++)
		new_mesh->elem_type_item[i]=mesh->elem_type_item[i];

	for(i=0;i<new_mesh->n_elem;i++) {
		for(j=0;j<new_mesh->n_elem_type;j++) {
			if(new_mesh->elem_type[i]==new_mesh->elem_type_item[j]) {
				new_mesh->elem_type_index[j+1]++;
				break;
			}
		}
	}
	for(j=1;j<new_mesh->n_elem_type+1;j++)
		new_mesh->elem_type_index[j]+=new_mesh->elem_type_index[j-1];
/*
	fprintf(stderr, "new_mesh: elem_type_index= %d %d %d\n", new_mesh->elem_type_index[0], new_mesh->elem_type_index[1], 
		new_mesh->elem_type_index[2]);
*/
	count_elem_index=(int *)calloc(new_mesh->n_elem_type, sizeof(int));
	if(count_elem_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: count_elem_index");
	for(i=0;i<new_mesh->n_elem_type;i++)
		count_elem_index[i]=0;
	for(i=0;i<new_mesh->n_elem;i++) {
		for(j=0;j<new_mesh->n_elem_type;j++) {
			if(new_mesh->elem_type[i]==new_mesh->elem_type_item[j]) {
				new2old[i]=new_mesh->elem_type_index[j]+count_elem_index[j];
				count_elem_index[j]++;
				break;
			}
		}
	}
	free(count_elem_index);
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) 
		new_tmp[i]=new_mesh->elem_type[new2old[i]];
	free(new_mesh->elem_type);
	new_mesh->elem_type=new_tmp;
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) 
		new_tmp[i]=new_mesh->section_ID[new2old[i]];
	free(new_mesh->section_ID);
	new_mesh->section_ID=new_tmp;
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) 
		new_tmp[i]=new_mesh->elem_mat_ID_item[new2old[i]];
	free(new_mesh->elem_mat_ID_item);
	new_mesh->elem_mat_ID_item=new_tmp;

/*
    fprintf(test_fp, "n_node=%d n_elem=%d nn_internal=%d ne_internal=%d\n", new_mesh->n_node, new_mesh->n_elem,
		new_mesh->nn_internal, new_mesh->ne_internal);
    for(i=0;i<new_mesh->n_elem;i++)
		fprintf(test_fp, "%d %d\n", i, new_mesh->elem_node_index[i+1]-new_mesh->elem_node_index[i]);
	fclose(test_fp);
*/
	
	new_tmp2=(int *)calloc(new_mesh->n_elem+1, sizeof(int));
	for(i=0;i<new_mesh->n_elem+1;i++)
		new_tmp2[i]=new_mesh->elem_node_index[i];
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) 
		new_tmp[i]=new_mesh->elem_node_index[new2old[i]+1]-new_mesh->elem_node_index[new2old[i]];
    for(i=1;i<new_mesh->n_elem+1;i++)
    	new_mesh->elem_node_index[i]=new_mesh->elem_node_index[i-1]+new_tmp[i-1];
	free(new_tmp);
	new_tmp=(int *)calloc(new_mesh->elem_node_index[new_mesh->n_elem], sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++){
		for(j=new_tmp2[new2old[i]];j<new_tmp2[new2old[i]+1];j++) {
			new_tmp[new_mesh->elem_node_index[i]+j-new_tmp2[new2old[i]]]=new_mesh->elem_node_item[j];
		}
	}
	free(new_mesh->elem_node_item);
	new_mesh->elem_node_item=new_tmp;
	free(new_tmp2);
	new_tmp=(int *)calloc(new_mesh->n_elem*2, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) {
		new_tmp[i*2]=new_mesh->elem_ID[new2old[i]*2];
		new_tmp[i*2+1]=new_mesh->elem_ID[new2old[i]*2+1];
	}
	free(new_mesh->elem_ID);
	new_mesh->elem_ID=new_tmp;
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++)
		new_tmp[i]=new_mesh->global_elem_ID[new2old[i]];
	free(new_mesh->global_elem_ID);
	new_mesh->global_elem_ID=new_tmp;
#ifdef TEST
	for(i=0;i<new_mesh->n_elem;i++)
		fprintf(test_fp, "i= %d  elem_ID=%d %d\n", i, new_mesh->elem_ID[i*2], new_mesh->elem_ID[i*2+1]);
    fclose(test_fp);
#endif
	old2new=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(old2new==NULL)
		HECMW_dlb_memory_exit("old2new");
	for(i=0;i<new_mesh->n_elem;i++) 
		old2new[new2old[i]]=i;
	new_mesh->ne_internal=0;
	for(i=0;i<new_mesh->n_elem;i++) {
		if(new_mesh->elem_ID[i*2+1]==mynode) {
				new_mesh->elem_internal_list[new_mesh->ne_internal]=i+1;
            	new_mesh->ne_internal++;
		}
	}

	new_mesh->my_rank=mesh->my_rank;
	new_mesh->zero=mesh->zero;
	new_mesh->PETOT=mesh->PETOT;
	new_mesh->PEsmpTOT=mesh->PEsmpTOT;
	new_mesh->errnof=mesh->errnof;
	HECMW_Comm_dup(mesh->HECMW_COMM, &new_mesh->HECMW_COMM);
	new_mesh->n_subdomain=mesh->n_subdomain;

	new_mesh->shared_index=(int *)calloc(new_mesh->n_neighbor_pe+1, sizeof(int));
	if(new_mesh->shared_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: shared_index");
	for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
		new_mesh->shared_index[i]=0;
   
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++)
		new_tmp[i]=new_mesh->when_i_was_refined_elem[new2old[i]];
	free(new_mesh->when_i_was_refined_elem);
	new_mesh->when_i_was_refined_elem=new_tmp;

	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++)
		new_tmp[i]=new_mesh->adapt_parent_type[new2old[i]];
	free(new_mesh->adapt_parent_type);
	new_mesh->adapt_parent_type=new_tmp;

	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++)
		new_tmp[i]=new_mesh->adapt_type[new2old[i]];
	free(new_mesh->adapt_type);
	new_mesh->adapt_type=new_tmp;

	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++)
		new_tmp[i]=new_mesh->adapt_level[new2old[i]];
	free(new_mesh->adapt_level);
	new_mesh->adapt_level=new_tmp;

	new_tmp=(int *)calloc(new_mesh->n_elem*2, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) {
		new_tmp[i*2]=new_mesh->adapt_parent[new2old[i]*2];
		new_tmp[i*2+1]=new_mesh->adapt_parent[new2old[i]*2+1];
	}
	free(new_mesh->adapt_parent);
	new_mesh->adapt_parent=new_tmp;

	new_tmp2=(int *)calloc(new_mesh->n_elem+1, sizeof(int));
	for(i=0;i<new_mesh->n_elem+1;i++)
		new_tmp2[i]=new_mesh->adapt_children_index[i];
	new_tmp=(int *)calloc(new_mesh->n_elem, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++) 
		new_tmp[i]=new_mesh->adapt_children_index[new2old[i]+1]-new_mesh->adapt_children_index[new2old[i]];
    for(i=1;i<new_mesh->n_elem+1;i++)
    	new_mesh->adapt_children_index[i]=new_mesh->adapt_children_index[i-1]+new_tmp[i-1];
	free(new_tmp);
	new_tmp=(int *)calloc(new_mesh->adapt_children_index[new_mesh->n_elem]*2, sizeof(int));
	if(new_tmp==NULL)
		HECMW_dlb_memory_exit("new_tmp");
	for(i=0;i<new_mesh->n_elem;i++){
		for(j=new_tmp2[new2old[i]];j<new_tmp2[new2old[i]+1];j++) {
			new_tmp[(new_mesh->adapt_children_index[i]+j-new_tmp2[new2old[i]])*2]=new_mesh->adapt_children_item[j*2];
			new_tmp[(new_mesh->adapt_children_index[i]+j-new_tmp2[new2old[i]])*2+1]=new_mesh->adapt_children_item[j*2+1];

		}
	}
	free(new_mesh->adapt_children_item);
	new_mesh->adapt_children_item=new_tmp;
	free(new_tmp2);
    if(mesh->section!=NULL) {
    new_mesh->section=(struct hecmwST_section *)calloc(1, sizeof(struct hecmwST_section));
	if(new_mesh->section==NULL)
		HECMW_dlb_memory_exit("new_mesh: section");
	new_mesh->section->n_sect=mesh->section->n_sect;
	new_mesh->section->sect_type=(int *)calloc(new_mesh->section->n_sect, sizeof(int));
	new_mesh->section->sect_opt=(int *)calloc(new_mesh->section->n_sect, sizeof(int));
	if((new_mesh->section->sect_type==NULL) || (new_mesh->section->sect_opt==NULL))
		HECMW_dlb_memory_exit("new_mesh: section");
	new_mesh->section->sect_mat_ID_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));
	new_mesh->section->sect_mat_ID_item=(int *)calloc(mesh->section->sect_mat_ID_index[mesh->section->n_sect], sizeof(int));
	for(i=0;i<mesh->section->n_sect;i++)
		new_mesh->section->sect_type[i]=mesh->section->sect_type[i];
	for(i=0;i<mesh->section->n_sect;i++)
		new_mesh->section->sect_opt[i]=mesh->section->sect_opt[i];
	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_mat_ID_index[i]=mesh->section->sect_mat_ID_index[i];
	for(i=0;i<mesh->section->sect_mat_ID_index[mesh->section->n_sect];i++)
		new_mesh->section->sect_mat_ID_item[i]=mesh->section->sect_mat_ID_item[i];
	new_mesh->section->sect_I_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));
	new_mesh->section->sect_R_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));

	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_I_index[i]=mesh->section->sect_I_index[i];
	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_R_index[i]=mesh->section->sect_R_index[i];
	new_mesh->section->sect_I_item=(int *)calloc(new_mesh->section->sect_I_index[new_mesh->section->n_sect], sizeof(int));
	new_mesh->section->sect_R_item=(double *)calloc(new_mesh->section->sect_R_index[new_mesh->section->n_sect], sizeof(double));
	for(i=0;i<new_mesh->section->sect_I_index[new_mesh->section->n_sect];i++)
		new_mesh->section->sect_I_item[i]=mesh->section->sect_I_item[i];
	for(i=0;i<new_mesh->section->sect_R_index[new_mesh->section->n_sect];i++)
		new_mesh->section->sect_R_item[i]=mesh->section->sect_R_item[i];
	}
    if(mesh->material!=NULL) {
    new_mesh->material=(struct hecmwST_material *)calloc(1, sizeof(struct hecmwST_material));
	if(new_mesh->material==NULL)
		HECMW_dlb_memory_exit("new_mesh: material");
	new_mesh->material->n_mat=mesh->material->n_mat;
	new_mesh->material->n_mat_item=mesh->material->n_mat_item;
	new_mesh->material->n_mat_subitem=mesh->material->n_mat_subitem;
	new_mesh->material->n_mat_table=mesh->material->n_mat_table;
	new_mesh->material->mat_name=(char **)calloc(new_mesh->material->n_mat, sizeof(char *));
	if(new_mesh->material->mat_name==NULL)
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->n_mat;i++) {
		new_mesh->material->mat_name[i]=(char *)calloc(128, sizeof(char));
		sprintf(new_mesh->material->mat_name[i], "%s", mesh->material->mat_name[i]);
	}
	new_mesh->material->mat_item_index=(int *)calloc(new_mesh->material->n_mat+1, sizeof(int));
	new_mesh->material->mat_subitem_index=(int *)calloc(new_mesh->material->n_mat_item+1, sizeof(int));
	new_mesh->material->mat_table_index=(int *)calloc(new_mesh->material->n_mat_subitem+1, sizeof(int));
	if((new_mesh->material->mat_item_index==NULL) || (new_mesh->material->mat_subitem_index==NULL) || 
		(new_mesh->material->mat_table_index==NULL))
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->n_mat+1;i++)
		new_mesh->material->mat_item_index[i]=mesh->material->mat_item_index[i];
	for(i=0;i<new_mesh->material->n_mat_item+1;i++)
		new_mesh->material->mat_subitem_index[i]=mesh->material->mat_subitem_index[i];
	for(i=0;i<new_mesh->material->n_mat_subitem+1;i++)
		new_mesh->material->mat_table_index[i]=mesh->material->mat_table_index[i];
	new_mesh->material->mat_val=(double *)calloc(new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem], 
		sizeof(double));
	new_mesh->material->mat_temp=(double *)calloc(new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem], 
		sizeof(double));
	if((new_mesh->material->mat_val==NULL) || (new_mesh->material->mat_temp==NULL))
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem];i++) {
		new_mesh->material->mat_val[i]=mesh->material->mat_val[i];
		new_mesh->material->mat_temp[i]=mesh->material->mat_temp[i];
	}
	}
	if(mesh->mpc!=NULL) {
		new_mesh->mpc=(struct hecmwST_mpc *)calloc(1, sizeof(struct hecmwST_mpc));
		if(new_mesh->mpc==NULL)
			HECMW_dlb_memory_exit("new_mesh: mpc");
		new_mesh->mpc->n_mpc=mesh->mpc->n_mpc;
		if(new_mesh->mpc->n_mpc!=0) {
			new_mesh->mpc->mpc_index=(int *)calloc(new_mesh->mpc->n_mpc+1, sizeof(int));
			new_mesh->mpc->mpc_item=(int *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(int));
			new_mesh->mpc->mpc_dof=(int *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(int));
			new_mesh->mpc->mpc_val=(double *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(double));
			for(i=0;i<mesh->mpc->n_mpc+1;i++)
				new_mesh->mpc->mpc_index[i]=mesh->mpc->mpc_index[i];
			for(i=0;i<mesh->mpc->n_mpc;i++) {
				new_mesh->mpc->mpc_item[i]=mesh->mpc->mpc_item[i];
				new_mesh->mpc->mpc_dof[i]=mesh->mpc->mpc_dof[i];
				new_mesh->mpc->mpc_val[i]=mesh->mpc->mpc_val[i];
			}
		}
	}
	if(mesh->amp!=NULL) {
		new_mesh->amp=(struct hecmwST_amplitude *)calloc(1, sizeof(struct hecmwST_amplitude));
		if(new_mesh->amp==NULL)
			HECMW_dlb_memory_exit("new_mesh: amp");
		new_mesh->amp->n_amp=mesh->amp->n_amp;
		if(new_mesh->amp->n_amp!=0) {
			new_mesh->amp->amp_type_definition=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_type_time=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_type_value=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_index=(int *)calloc(mesh->amp->n_amp+1, sizeof(int));
			new_mesh->amp->amp_val=(double *)calloc(mesh->amp->amp_index[mesh->amp->n_amp], sizeof(double));
			new_mesh->amp->amp_table=(double *)calloc(mesh->amp->amp_index[mesh->amp->n_amp], sizeof(double));
			if((new_mesh->amp->amp_type_definition==NULL) || (new_mesh->amp->amp_type_time==NULL) || 
				(new_mesh->amp->amp_type_value==NULL) || (new_mesh->amp->amp_index==NULL) || 
				(new_mesh->amp->amp_val==NULL) || (new_mesh->amp->amp_table==NULL))
				HECMW_dlb_memory_exit("new_mesh: amp");
			for(i=0;i<mesh->amp->n_amp;i++) {
				new_mesh->amp->amp_type_definition[i]=mesh->amp->amp_type_definition[i];
				new_mesh->amp->amp_type_time[i]=mesh->amp->amp_type_time[i];
				new_mesh->amp->amp_type_value[i]=mesh->amp->amp_type_value[i];
			}
			for(i=0;i<mesh->amp->n_amp+1;i++) 
				new_mesh->amp->amp_index[i]=mesh->amp->amp_index[i];
			for(i=0;i<mesh->amp->amp_index[mesh->amp->n_amp];i++) {
				new_mesh->amp->amp_val[i]=mesh->amp->amp_val[i];
				new_mesh->amp->amp_table[i]=mesh->amp->amp_table[i];
			}
		}
	}

	if(new_mesh->elem_group->n_grp>0) {
		if(new_mesh->elem_group->grp_index[new_mesh->elem_group->n_grp]>0) {
			new_tmp=(int *)calloc(new_mesh->elem_group->grp_index[new_mesh->elem_group->n_grp], sizeof(int));
			for(i=0;i<new_mesh->elem_group->grp_index[new_mesh->elem_group->n_grp];i++)
				new_tmp[i]=old2new[new_mesh->elem_group->grp_item[i]-1]+1;
			free(new_mesh->elem_group->grp_item);
			new_mesh->elem_group->grp_item=new_tmp;
		}
	}

	if(new_mesh->surf_group->n_grp>0) {
		if(new_mesh->surf_group->grp_index[new_mesh->surf_group->n_grp]>0) {
			new_tmp=(int *)calloc(new_mesh->surf_group->grp_index[new_mesh->surf_group->n_grp]*2, sizeof(int));
			for(i=0;i<new_mesh->surf_group->grp_index[new_mesh->surf_group->n_grp];i++) {
				new_tmp[i*2]=old2new[new_mesh->surf_group->grp_item[i*2]-1]+1;
				new_tmp[i*2+1]=new_mesh->surf_group->grp_item[i*2+1];
			}
			free(new_mesh->surf_group->grp_item);
			new_mesh->surf_group->grp_item=new_tmp;
		}
	}
	}
else {

    new_mesh->elem_type_index=(int *)calloc(2, sizeof(int));
	if(new_mesh->elem_type_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_type_index");
	for(i=0;i<new_mesh->n_elem_type+1;i++)
		new_mesh->elem_type_index[i]=0;
	new_mesh->elem_type_item=(int *)calloc(new_mesh->n_elem_type, sizeof(int));
	if(new_mesh->elem_type_item==NULL)
		HECMW_dlb_memory_exit("new_mesh: elem_type_item");
	for(i=0;i<new_mesh->n_elem_type;i++)
		new_mesh->elem_type_item[i]=mesh->elem_type_item[i];
	new_mesh->elem_type_index[1]=new_mesh->n_elem;

/*
	fprintf(stderr, "new_mesh: elem_type_index= %d %d %d\n", new_mesh->elem_type_index[0], new_mesh->elem_type_index[1], 
		new_mesh->elem_type_index[2]);
*/

/*
	for(i=0;i<new_mesh->n_elem;i++)
		fprintf(test_fp, "i= %d  elem_ID=%d %d\n", i, new_mesh->elem_ID[i*2], new_mesh->elem_ID[i*2+1]);
    fclose(test_fp);
*/
	new_mesh->my_rank=mesh->my_rank;
	new_mesh->zero=mesh->zero;
	new_mesh->PETOT=mesh->PETOT;
	new_mesh->PEsmpTOT=mesh->PEsmpTOT;
	new_mesh->errnof=mesh->errnof;
	HECMW_Comm_dup(mesh->HECMW_COMM, &new_mesh->HECMW_COMM);
	new_mesh->n_subdomain=mesh->n_subdomain;

	new_mesh->shared_index=(int *)calloc(new_mesh->n_neighbor_pe+1, sizeof(int));
	if(new_mesh->shared_index==NULL)
		HECMW_dlb_memory_exit("new_mesh: shared_index");
	for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
		new_mesh->shared_index[i]=0;
   
    if(mesh->section!=NULL) {
    new_mesh->section=(struct hecmwST_section *)calloc(1, sizeof(struct hecmwST_section));
	if(new_mesh->section==NULL)
		HECMW_dlb_memory_exit("new_mesh: section");
	new_mesh->section->n_sect=mesh->section->n_sect;
	new_mesh->section->sect_type=(int *)calloc(new_mesh->section->n_sect, sizeof(int));
	new_mesh->section->sect_opt=(int *)calloc(new_mesh->section->n_sect, sizeof(int));
	if((new_mesh->section->sect_type==NULL) || (new_mesh->section->sect_opt==NULL))
		HECMW_dlb_memory_exit("new_mesh: section");
	new_mesh->section->sect_mat_ID_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));
	new_mesh->section->sect_mat_ID_item=(int *)calloc(mesh->section->sect_mat_ID_index[mesh->section->n_sect], sizeof(int));
	for(i=0;i<mesh->section->n_sect;i++)
		new_mesh->section->sect_type[i]=mesh->section->sect_type[i];
	for(i=0;i<mesh->section->n_sect;i++)
		new_mesh->section->sect_opt[i]=mesh->section->sect_opt[i];
	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_mat_ID_index[i]=mesh->section->sect_mat_ID_index[i];
	for(i=0;i<mesh->section->sect_mat_ID_index[mesh->section->n_sect];i++)
		new_mesh->section->sect_mat_ID_item[i]=mesh->section->sect_mat_ID_item[i];
	new_mesh->section->sect_I_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));
	new_mesh->section->sect_R_index=(int *)calloc(new_mesh->section->n_sect+1, sizeof(int));

	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_I_index[i]=mesh->section->sect_I_index[i];
	for(i=0;i<mesh->section->n_sect+1;i++)
		new_mesh->section->sect_R_index[i]=mesh->section->sect_R_index[i];
	}
    if(mesh->material!=NULL) {
    new_mesh->material=(struct hecmwST_material *)calloc(1, sizeof(struct hecmwST_material));
	if(new_mesh->material==NULL)
		HECMW_dlb_memory_exit("new_mesh: material");
	new_mesh->material->n_mat=mesh->material->n_mat;
	new_mesh->material->n_mat_item=mesh->material->n_mat_item;
	new_mesh->material->n_mat_subitem=mesh->material->n_mat_subitem;
	new_mesh->material->n_mat_table=mesh->material->n_mat_table;
	new_mesh->material->mat_name=(char **)calloc(new_mesh->material->n_mat, sizeof(char *));
	if(new_mesh->material->mat_name==NULL)
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->n_mat;i++) {
		new_mesh->material->mat_name[i]=(char *)calloc(128, sizeof(char));
		sprintf(new_mesh->material->mat_name[i], "%s", mesh->material->mat_name[i]);
	}
	new_mesh->material->mat_item_index=(int *)calloc(new_mesh->material->n_mat+1, sizeof(int));
	new_mesh->material->mat_subitem_index=(int *)calloc(new_mesh->material->n_mat_item+1, sizeof(int));
	new_mesh->material->mat_table_index=(int *)calloc(new_mesh->material->n_mat_subitem+1, sizeof(int));
	if((new_mesh->material->mat_item_index==NULL) || (new_mesh->material->mat_subitem_index==NULL) || 
		(new_mesh->material->mat_table_index==NULL))
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->n_mat+1;i++)
		new_mesh->material->mat_item_index[i]=mesh->material->mat_item_index[i];
	for(i=0;i<new_mesh->material->n_mat_item+1;i++)
		new_mesh->material->mat_subitem_index[i]=mesh->material->mat_subitem_index[i];
	for(i=0;i<new_mesh->material->n_mat_subitem+1;i++)
		new_mesh->material->mat_table_index[i]=mesh->material->mat_table_index[i];
	new_mesh->material->mat_val=(double *)calloc(new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem], 
		sizeof(double));
	new_mesh->material->mat_temp=(double *)calloc(new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem], 
		sizeof(double));
	if((new_mesh->material->mat_val==NULL) || (new_mesh->material->mat_temp==NULL))
		HECMW_dlb_memory_exit("new_mesh: material");
	for(i=0;i<new_mesh->material->mat_table_index[new_mesh->material->n_mat_subitem];i++) {
		new_mesh->material->mat_val[i]=mesh->material->mat_val[i];
		new_mesh->material->mat_temp[i]=mesh->material->mat_temp[i];
	}
	}
	if(mesh->mpc!=NULL) {
		new_mesh->mpc=(struct hecmwST_mpc *)calloc(1, sizeof(struct hecmwST_mpc));
		if(new_mesh->mpc==NULL)
			HECMW_dlb_memory_exit("new_mesh: mpc");
		new_mesh->mpc->n_mpc=mesh->mpc->n_mpc;
		if(new_mesh->mpc->n_mpc!=0) {
			new_mesh->mpc->mpc_index=(int *)calloc(new_mesh->mpc->n_mpc+1, sizeof(int));
			new_mesh->mpc->mpc_item=(int *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(int));
			new_mesh->mpc->mpc_dof=(int *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(int));
			new_mesh->mpc->mpc_val=(double *)calloc(mesh->mpc->mpc_index[new_mesh->mpc->n_mpc], sizeof(double));
			for(i=0;i<mesh->mpc->n_mpc+1;i++)
				new_mesh->mpc->mpc_index[i]=mesh->mpc->mpc_index[i];
			for(i=0;i<mesh->mpc->n_mpc;i++) {
				new_mesh->mpc->mpc_item[i]=mesh->mpc->mpc_item[i];
				new_mesh->mpc->mpc_dof[i]=mesh->mpc->mpc_dof[i];
				new_mesh->mpc->mpc_val[i]=mesh->mpc->mpc_val[i];
			}
		}
	}
	if(mesh->amp!=NULL) {
		new_mesh->amp=(struct hecmwST_amplitude *)calloc(1, sizeof(struct hecmwST_amplitude));
		if(new_mesh->amp==NULL)
			HECMW_dlb_memory_exit("new_mesh: amp");
		new_mesh->amp->n_amp=mesh->amp->n_amp;
		if(new_mesh->amp->n_amp!=0) {
			new_mesh->amp->amp_type_definition=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_type_time=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_type_value=(int *)calloc(mesh->amp->n_amp, sizeof(int));
			new_mesh->amp->amp_index=(int *)calloc(mesh->amp->n_amp+1, sizeof(int));
			new_mesh->amp->amp_val=(double *)calloc(mesh->amp->amp_index[mesh->amp->n_amp], sizeof(double));
			new_mesh->amp->amp_table=(double *)calloc(mesh->amp->amp_index[mesh->amp->n_amp], sizeof(double));
			if((new_mesh->amp->amp_type_definition==NULL) || (new_mesh->amp->amp_type_time==NULL) || 
				(new_mesh->amp->amp_type_value==NULL) || (new_mesh->amp->amp_index==NULL) || 
				(new_mesh->amp->amp_val==NULL) || (new_mesh->amp->amp_table==NULL))
				HECMW_dlb_memory_exit("new_mesh: amp");
			for(i=0;i<mesh->amp->n_amp;i++) {
				new_mesh->amp->amp_type_definition[i]=mesh->amp->amp_type_definition[i];
				new_mesh->amp->amp_type_time[i]=mesh->amp->amp_type_time[i];
				new_mesh->amp->amp_type_value[i]=mesh->amp->amp_type_value[i];
			}
			for(i=0;i<mesh->amp->n_amp+1;i++) 
				new_mesh->amp->amp_index[i]=mesh->amp->amp_index[i];
			for(i=0;i<mesh->amp->amp_index[mesh->amp->n_amp];i++) {
				new_mesh->amp->amp_val[i]=mesh->amp->amp_val[i];
				new_mesh->amp->amp_table[i]=mesh->amp->amp_table[i];
			}
		}
	}
	}
  new_mesh->node_internal_list=(int *)calloc(new_mesh->nn_internal, sizeof(int));
  if(new_mesh->node_internal_list==NULL)
    HECMW_dlb_memory_exit("node_internal_list");
  for(i=0;i<new_mesh->nn_internal;i++)
    new_mesh->node_internal_list[i]=i;
  if(mynode==0) {
      t2=HECMW_Wtime();
      fprintf(stderr, "Finish migration now. The time cost for migration and generating new mesh is %lf\n", t2-t1);
      }
	return;
	}





			
		




