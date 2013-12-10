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

extern void mesh2graph(struct hecmwST_local_mesh *mesh, GraphType *graph, Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], 
					   Result_part *result, HECMW_Comm repart_comm);
extern void redistribute_mesh(GraphType *graph, 
				Result_part *result, int mynode, int pesize);
/*extern void write_mesh_display(char *outfile, struct local_mesh *mesh, struct node_elem_data *node); 
extern void write3_mesh_display(char *outfile, struct local_mesh *mesh, struct node_elem_data *node); 
*/
extern void hecmw_dlb_read_control(char *contfile, Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], int pesize);
extern void hecmw_dlb_set_default_control(Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], int pesize);


void hecmw_dynamic_load_balancing_()
{
	int mynode, pesize;
	FILE  *test_fp;
	char  test_filename[128], out_filename[128];
	int  i, j;
    int  num_leaf_elem;

   GraphType *graph;
  int maxnvtxs = -1, maxnedges = -1;
  int readew = -1, readvw = -1, dummy, edge;
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt;
  idxtype *your_xadj, *your_adjncy, *your_vwgt, *your_adjwgt, graphinfo[4];
  int fmt, ncon, nobj;

  struct Global_inf *global_index;
  Result_part *result;
  double t1, t2, t3;

   int   *num_move;
/*   int  l_child, tmp_int, *inter_elem;
*/

   int cglevel, adaplevel, tmp, neibpetot, *neibpe, allnodtotcur, intnodetotcur, *idnode, *whennode, *xyz;
   int num_char;
   char *contfile, buf[HECMW_FILENAME_LEN];
   Control_para *ctl_para;
   int stat_para[NUM_CONTROL_PARAS];
   int flag_control;
  HECMW_Comm_rank(mesh->HECMW_COMM, &mynode);
  HECMW_Comm_size(mesh->HECMW_COMM, &pesize);
  if(mynode==0) {
      t1=HECMW_Wtime();
      fprintf(stderr, "Start graph repartition now ...\n");
      }
  
    ctl_para=(Control_para *)malloc(sizeof(Control_para));
	if(ctl_para==NULL)
		  HECMW_dlb_memory_exit("ctl_para");
    flag_control=HECMW_ctrl_is_exists_control("dlb-ctrl");
    if(flag_control==0)
	   hecmw_dlb_set_default_control(ctl_para,stat_para, pesize);
    else {
    contfile=HECMW_ctrl_get_control_file("dlb-ctrl", buf, HECMW_FILENAME_LEN); 
   if(contfile!=NULL)
       hecmw_dlb_read_control(contfile, ctl_para, stat_para, pesize);
    }  
  
   

  graph=(GraphType *)malloc(sizeof(GraphType));
  if(graph==NULL) 
     HECMW_dlb_memory_exit("graph");
/*  global_index=(struct Global_inf *)malloc(sizeof(struct Global_inf));
*/
  result=(Result_part *)malloc(sizeof(Result_part));
  if(result==NULL) 
	  HECMW_dlb_memory_exit("result");
  

  mesh2graph(mesh, graph,  ctl_para, stat_para,result, mesh->HECMW_COMM);
  if(mynode==0)
     fprintf(stderr, "the edgecut number is %d\n", result->edgecut);
/*
  num_move=(int *)calloc(pesize, sizeof(int));
  for(i=0;i<pesize;i++) {
	  num_move[i]=0;
  }
  for(i=0;i<mesh->nn_internal;i++) {
	  num_move[result->part[i]]++;
  }
  for(i=0;i<pesize;i++)
      fprintf(stderr, "in PE %d  the number of nodes belong to PE %d is %d\n", mynode, i, num_move[i]);
  */
  if(mynode==0) {
      t2=HECMW_Wtime();
      fprintf(stderr, "Finish repartition now. \n");
	  fprintf(stderr, "The time for repartition is %lf\n", t2-t1);
	  fprintf(stderr, "Start migration...\n");
  }
 


  HECMW_Barrier(mesh->HECMW_COMM);
  new_mesh=(struct hecmwST_local_mesh *)calloc(1, sizeof(struct hecmwST_local_mesh));
  if(new_mesh==NULL)
	  HECMW_dlb_memory_exit("new_mesh");
  new_data=(struct hecmwST_result_data *)calloc(1, sizeof(struct hecmwST_result_data));
  if(new_data==NULL)
	  HECMW_dlb_memory_exit("new_data");
/*  global_comm_table=(Comm_table *)malloc(sizeof(Comm_table));
  if(global_comm_table==NULL)
	  HECMW_dlb_memory_exit("global_comm_table");
	  */

  redistribute_mesh(graph, result, mynode, pesize);
/*
  if(mynode==0) {
      t3=HECMW_Wtime();
      fprintf(stderr, "Finish migration\n");
	  fprintf(stderr, "The time for migration is %lf\n", t3-t2);
	  t3=t2;
      }
*/
  return;
}


