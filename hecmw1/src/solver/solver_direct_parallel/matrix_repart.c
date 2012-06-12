/*=====================================================================*
 *                                                                     *
 *  Software Name : HPC-MW Library for PC-cluster                      *
 *        Version : 1.00                                               *
 *                                                                     *
 *    Last Update : 2010/06/16                                         *
 *    Modified by : KITAYAMA Takeshi (Univ. of Tokyo)                  *
 *                                                                     *
 *       Category : Dynamic Load Balancing                             *
 *         Author : Li Chen (RIST)                                     *
 *                                                                     *
 *                                                                     *
 *                The University of Tokyo, FSIS Project                *
 *                                                                     *
 *       "High-Performance Computing Middleware (HPC-MW)" Group.       *
 *                                                                     *
 *=====================================================================*/

#include "separator.h"
Separator_result *separator;

extern void matrix2graph(int num_of_row, int num_of_col, int num_of_nzero, int *irow, int *jcol, GraphType *graph);


void bi_part_directive(int *neqns, int *nttbr, int *irow, int *jcol, int *num_graph1, int *num_graph2, int *num_separator)
{
	FILE  *test_fp;
	char  test_filename[128], out_filename[128];
	int  i, j;

   GraphType *graph;
  int maxnvtxs = -1, maxnedges = -1;
  int readew = -1, readvw = -1, dummy, edge;
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt;
  idxtype *your_xadj, *your_adjncy, *your_vwgt, *your_adjwgt, graphinfo[4];
  int fmt, ncon, nobj;
  int num_of_row, num_of_col, num_of_nzero;
  double t1, t2, t3;

   int   *num_move;

   int cglevel, adaplevel, tmp, neibpetot, *neibpe, allnodtotcur, intnodetotcur, *idnode, *whennode, *xyz;
   int num_char;
   Control_para *ctl_para;
   int stat_para[NUM_CONTROL_PARAS];
   int flag_control;
   int options[8];
   int *perm, *iperm;
   int num_flag;
  
  
	num_of_row=*neqns;
	num_of_col=*neqns;
	num_of_nzero=*nttbr;
   
  
  graph=(GraphType *)malloc(sizeof(GraphType));
  if(graph==NULL) 
     separator_memory_exit("graph");
  
  fprintf(stderr, "Start transforming matrix to graph\n");
  matrix2graph(num_of_row, num_of_col, num_of_nzero, irow, jcol, graph);
  fprintf(stderr, "Graph Information ---------------------------------------------------\n");
  fprintf(stderr, "#Vertices: %d, #Edges: %d\n\n",  graph->nvtxs, graph->nedges/2);
  perm=(int *)calloc(num_of_col, sizeof(int));
  iperm=(int *)calloc(num_of_col, sizeof(int));
  if((perm==NULL) || (iperm==NULL))
	  separator_memory_exit("matrix_repart: perm, iperm");

  /* followings are options. see METIS manual for METIS_NODEND() */
	options[0]=1; /* specify parameters */
	options[1]=3; /* default */
	options[2]=1; /* default */
	options[3]=2; /* default */
	options[4]=0; /* default */
	options[5]=0; /* do not try to compress the matrix */
	options[6]=0; /* default */
	options[7]=1; /* default */

  num_flag=0;
  METIS_NodeND(&graph->nvtxs, graph->xadj, graph->adjncy, &num_flag, options, perm, iperm);
  /* copy to separator */
  
	*num_graph1=separator->num_of_lgraph;
	*num_graph2=separator->num_of_rgraph;
	*num_separator=separator->num_of_separator;


  return;
}


void bi_part_directive_(int *neqns, int *nttbr, int *irow, int *jcol, int *num_graph1, int *num_graph2, int *num_separator)
{
  bi_part_directive(neqns, nttbr, irow, jcol, num_graph1, num_graph2, num_separator);
}

void bi_part_directive__(int *neqns, int *nttbr, int *irow, int *jcol, int *num_graph1, int *num_graph2, int *num_separator)
{
  bi_part_directive(neqns, nttbr, irow, jcol, num_graph1, num_graph2, num_separator);
}

void BI_PART_DIRECTIVE(int *neqns, int *nttbr, int *irow, int *jcol, int *num_graph1, int *num_graph2, int *num_separator)
{
  bi_part_directive(neqns, nttbr, irow, jcol, num_graph1, num_graph2, num_separator);
}
