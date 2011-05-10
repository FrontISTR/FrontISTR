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

void find_new_v8(int new_v8[24])
{

    new_v8[0]=1;  new_v8[1]=3;  new_v8[2]=4;
	new_v8[1*3]=0; new_v8[1*3+1]=2; new_v8[1*3+2]=5;
	new_v8[2*3]=1; new_v8[2*3+1]=3; new_v8[2*3+2]=6;
	new_v8[3*3]=0; new_v8[3*3+1]=2; new_v8[3*3+2]=7;
	new_v8[4*3]=0; new_v8[4*3+1]=5; new_v8[4*3+2]=7;
	new_v8[5*3]=1; new_v8[5*3+1]=4; new_v8[5*3+2]=6;
	new_v8[6*3]=2; new_v8[6*3+1]=5; new_v8[6*3+2]=7;
	new_v8[7*3]=3; new_v8[7*3+1]=4; new_v8[7*3+2]=6;
	return;
}

void find_new_v6(int new_v6[18])
{

    new_v6[0]=1;  new_v6[1]=2;  new_v6[2]=3;
	new_v6[1*3]=0; new_v6[1*3+1]=2; new_v6[1*3+2]=4;
	new_v6[2*3]=0; new_v6[2*3+1]=1; new_v6[2*3+2]=5;
	new_v6[3*3]=0; new_v6[3*3+1]=4; new_v6[3*3+2]=5;
	new_v6[4*3]=1; new_v6[4*3+1]=3; new_v6[4*3+2]=5;
	new_v6[5*3]=2; new_v6[5*3+1]=3; new_v6[5*3+2]=4;
	return;
}

void find_new_v4(int new_v4[12])
{
	new_v4[0]=1;   new_v4[1]=2;   new_v4[2]=3;
	new_v4[1*3]=0; new_v4[1*3+1]=2;  new_v4[1*3+2]=3;
	new_v4[2*3]=0; new_v4[2*3+1]=1;  new_v4[2*3+2]=3;
	new_v4[3*3]=0; new_v4[3*3+1]=1;  new_v4[3*3+2]=2;
	return;
}
		


void add8_adj_link(Adj_find *adj_link, int id_elem,  struct hecmwST_local_mesh *mesh, int new_v8[24])
{
    Adj_find *p1, *p2;
	int  i,j, k, m;
	int  new_v, id_node, flag_hit;
	int  ncon;

    ncon=mesh->elem_node_index[id_elem+1]-mesh->elem_node_index[id_elem];
	if(ncon!=8) 
		HECMW_dlb_print_exit("The type of element is conflict with its index: data error\n");
    for(j=mesh->elem_node_index[id_elem];j<mesh->elem_node_index[id_elem+1];j++) {
		id_node=mesh->elem_node_item[j]-1;
		if((id_node>mesh->n_node) || (id_node<0))
			HECMW_dlb_print_exit("There is something wrong in index_elem\n");
		if(id_node<mesh->nn_internal) {
			for(k=0;k<3;k++) {
				new_v=mesh->elem_node_item[mesh->elem_node_index[id_elem]+new_v8[(j-mesh->elem_node_index[id_elem])*3+k]]-1;
				flag_hit=0;
				p1=adj_link[id_node].next_vertex;
				for(m=0;m<adj_link[id_node].vertex_num;m++) {
					if(new_v==p1->vertex_num) {
						flag_hit=1;
						break;
					}
				    p1=p1->next_vertex;
				}
				if(flag_hit==0) { /* adding the vertex new_v */
					p1=(Adj_find *)malloc(sizeof(Adj_find));
					if(p1==NULL) {
						fprintf(stderr, "There is no enough memory for p1 in adj_link\n");
						exit(0);
					}
					p2=adj_link[id_node].next_vertex;
					adj_link[id_node].vertex_num++;
					adj_link[id_node].next_vertex=p1;
					p1->next_vertex=p2;
					p1->vertex_num=new_v;
				}
			}
		}
	}
	return;
}

void add6_adj_link(Adj_find *adj_link, int id_elem,  struct hecmwST_local_mesh *mesh, int new_v6[18])
{
    Adj_find *p1, *p2;
	int  i,j, k, m;
	int  new_v, id_node, flag_hit;
	int  ncon;

    ncon=mesh->elem_node_index[id_elem+1]-mesh->elem_node_index[id_elem];
	if(ncon!=6) 
		HECMW_dlb_print_exit("The type of element is conflict with its index: data error\n");
    for(j=mesh->elem_node_index[id_elem];j<mesh->elem_node_index[id_elem+1];j++) {
		id_node=mesh->elem_node_item[j]-1;
		if((id_node>=mesh->n_node) || (id_node<0)) 
			HECMW_dlb_print_exit("There is something wrong in index_elem\n");
		if(id_node<mesh->nn_internal) {
			for(k=0;k<3;k++) {
				new_v=mesh->elem_node_item[mesh->elem_node_index[id_elem]+new_v6[(j-mesh->elem_node_index[id_elem])*3+k]]-1;
				flag_hit=0;
				p1=adj_link[id_node].next_vertex;
				for(m=0;m<adj_link[id_node].vertex_num;m++) {
					if(new_v==p1->vertex_num) {
						flag_hit=1;
						break;
					}
				    p1=p1->next_vertex;
				}
				if(flag_hit==0) { /* adding the vertex new_v */
					p1=(Adj_find *)malloc(sizeof(Adj_find));
					if(p1==NULL) {
						fprintf(stderr, "There is no enough memory for p1 in adj_link\n");
						exit(0);
					}
					p2=adj_link[id_node].next_vertex;
					adj_link[id_node].vertex_num++;
					adj_link[id_node].next_vertex=p1;
					p1->next_vertex=p2;
					p1->vertex_num=new_v;
				}
			}
		}
	}
	return;
}

void add4_adj_link(Adj_find *adj_link, int id_elem,  struct hecmwST_local_mesh *mesh, int new_v4[12])
{
    Adj_find *p1, *p2;
	int  i,j, k, m;
	int  new_v, id_node, flag_hit;
	int  ncon;

    ncon=mesh->elem_node_index[id_elem+1]-mesh->elem_node_index[id_elem];
	if(ncon!=4) 
		HECMW_dlb_print_exit("The type of element is conflict with its index: data error\n");
    for(j=mesh->elem_node_index[id_elem];j<mesh->elem_node_index[id_elem+1];j++) {
		id_node=mesh->elem_node_item[j]-1;
		if((id_node>=mesh->n_node) || (id_node<0)) 
			HECMW_dlb_print_exit("There is something wrong in index_elem\n");
		if(id_node<mesh->nn_internal) {
			for(k=0;k<3;k++) {
				new_v=mesh->elem_node_item[mesh->elem_node_index[id_elem]+new_v4[(j-mesh->elem_node_index[id_elem])*3+k]]-1;
				flag_hit=0;
				p1=adj_link[id_node].next_vertex;
				for(m=0;m<adj_link[id_node].vertex_num;m++) {
					if(new_v==p1->vertex_num) {
						flag_hit=1;
						break;
					}
				    p1=p1->next_vertex;
				}
				if(flag_hit==0) { /* adding the vertex new_v */
					p1=(Adj_find *)malloc(sizeof(Adj_find));
					if(p1==NULL) 
						HECMW_dlb_memory_exit("Adj_find: p1");
					p2=adj_link[id_node].next_vertex;
					adj_link[id_node].vertex_num++;
					adj_link[id_node].next_vertex=p1;
					p1->next_vertex=p2;
					p1->vertex_num=new_v;
				}
			}
		}
	}
	return;
}

void adj_link_free(Adj_find *adj_link, int num)
{
	Adj_find *p1, *p2;
	int i,j;
	for(i=0;i<num;i++) {
		p1=adj_link[i].next_vertex;
		for(j=0;j<adj_link[i].vertex_num;j++) {
			p2=p1;
			p1=p1->next_vertex;
			free(p2);
		}
	}
	free(adj_link);
	return;
}
		


void mesh2graph(struct hecmwST_local_mesh *mesh, GraphType *graph, Control_para *ctl_para, int stat_para[NUM_CONTROL_PARAS], 
				Result_part *result, HECMW_Comm repart_comm)
{
	int i,j,k,m;
	int mynode, pesize;
	Adj_find *adj_link;
	Adj_find *p1, *p2;
	int new_v8[3*8], new_v4[3*4], new_v6[3*6];
  int nvtxs  = 0, nedges = 0, global_num_node=0;
  
  int readew = -1, readvw = -1, dummy, edge;
  int *vtxdist, *xadj, *adjncy, *vwgts, *adjwgts;
  int wgtflag, numflag, ncon, nparts, edgecut, options[4];
  float *tpwgts, *ubvec, itr;
  HECMW_Status	stat;
  int *tmp_index, tmp_nvtxs, tmp_sum, tmp_pe, tmp_lid;
   FILE *fp_test, *fp_wgts;
   char test_file[128];
   float *xyz;
   int ndim;
   int tmp_int;
   int flag_count;

    HECMW_Comm_rank(repart_comm, &mynode);
    HECMW_Comm_size(repart_comm, &pesize);
/*
   mynode=mesh->my_rank;
   pesize=mesh->PEtotMESH;
   */
/*
   fprintf(stderr, "In PE %d  pesize=%d\n", mynode, pesize);
    fprintf(stderr, "wgtflag=%d \n", ctl_para->wgtflag);
    fprintf(stderr, "ncon=%d\n", ctl_para->num_repartition);
*/
	if(mynode==0) 
		fprintf(stderr, "Start transform original mesh data into graph structure of ParMetis\n");

	if(pesize>1)
	  HECMW_Allreduce(&mesh->nn_internal, &global_num_node, 1, HECMW_INT, HECMW_SUM, repart_comm);
	else
		global_num_node=mesh->nn_internal;
/*
	fprintf(stderr, "the global node number is %d\n", global_num_node);
*/  
  result->t_node=global_num_node;
	nvtxs=mesh->nn_internal;	
	graph->vtxdist=(int *)calloc((pesize+1),sizeof(int));
	if(graph->vtxdist==NULL) 
		HECMW_dlb_memory_exit("graph->vtxdist");
	if(mynode==0) {
		graph->vtxdist[0]=0;
		graph->vtxdist[1]=nvtxs;
		tmp_sum=nvtxs;
		for(i=1;i<pesize;i++) {
			HECMW_Recv(&tmp_nvtxs, 1, HECMW_INT, i, HECMW_ANY_TAG, repart_comm, &stat);
			tmp_sum+=tmp_nvtxs;
			graph->vtxdist[i+1]=tmp_sum;
		}
		for(i=1;i<pesize;i++)
	       HECMW_Send(graph->vtxdist,pesize+1,HECMW_INT, i, 0, repart_comm);
	}
	else {
		HECMW_Send(&nvtxs, 1, HECMW_INT, 0, 0, repart_comm);
        HECMW_Recv(graph->vtxdist, pesize+1, HECMW_INT, 0, HECMW_ANY_TAG, repart_comm, &stat);
	}
/*
	if(mynode==0) {
		for(i=0;i<pesize+1;i++)
			fprintf(stderr, "vtxdist=%d ", graph->vtxdist[i]);
		fprintf(stderr, "\n");
	}
*/	

	adj_link=(Adj_find *)calloc(mesh->nn_internal, sizeof(Adj_find));
    if(adj_link==NULL) 
		HECMW_dlb_memory_exit("adj_link");
    for(i=0;i<mesh->nn_internal;i++)
		adj_link[i].vertex_num=0;
	find_new_v8(new_v8);
	find_new_v6(new_v6);
	find_new_v4(new_v4);
	for(i=0;i<mesh->n_elem;i++) {
		flag_count=1;
		if(mesh->hecmw_flag_adapt==1) {
			if(mesh->adapt_type[i]!=0)
				flag_count=0;
		}
		if(flag_count==1) {

            if(mesh->elem_type[i]==341)
			     add4_adj_link(adj_link, i, mesh, new_v4);
		    else if(mesh->elem_type[i]==351)
				add6_adj_link(adj_link, i, mesh, new_v6);
			else if(mesh->elem_type[i]==361)
		        add8_adj_link(adj_link, i, mesh, new_v8);
		}

	}
	for(i=0;i<mesh->nn_internal;i++) 
		nedges+=adj_link[i].vertex_num;
	graph->xadj=(int *)calloc(nvtxs+1, sizeof(int));
	graph->adjncy=(int *)calloc(nedges, sizeof(int));
	if((graph->xadj==NULL) || (graph->adjncy==NULL)) 
		HECMW_dlb_memory_exit("graph: xadj and adjncy");
	m=0;
	graph->xadj[0]=0;
	for(i=0;i<mesh->nn_internal;i++) {
		graph->xadj[i+1]=m+adj_link[i].vertex_num;
		p1=adj_link[i].next_vertex;
		for(j=0;j<adj_link[i].vertex_num;j++) {
			if(p1->vertex_num<mesh->nn_internal)
				graph->adjncy[m+j]=graph->vtxdist[mynode]+p1->vertex_num;
			else {
				tmp_pe=mesh->node_ID[p1->vertex_num*2+1];
				tmp_lid=mesh->node_ID[p1->vertex_num*2]-1;
				graph->adjncy[m+j]=graph->vtxdist[tmp_pe]+tmp_lid;
			}
			p1=p1->next_vertex;
		}
		m+=adj_link[i].vertex_num;
	}
    adj_link_free(adj_link, mesh->nn_internal);
/*
	sprintf(test_file, "test.%d", mynode);
	fp_test=fopen(test_file, "w");
	tmp_int=0;
    for(i=0;i<mesh->n_internal;i++) {
		fprintf(fp_test, "%d  %d\n", i+1, graph->xadj[i+1]-graph->xadj[i]);
		if((graph->xadj[i+1]-graph->xadj[i])==27) {
			tmp_int++;
			fprintf(stderr, "the num of 27 edges is in PE %d nodeid is %d\n", mynode, i+1);
		}
		for(j=graph->xadj[i];j<graph->xadj[i+1];j++)
			fprintf(fp_test, "%d  ", graph->adjncy[j]+1);
		fprintf(fp_test, "\n");
	}
		fprintf(stderr, "the num of 27 edges is %d in PE %d\n", tmp_int, mynode);
	fclose(fp_test);
	*/
	graph->gnvtxs=global_num_node;
	graph->ncon=0; graph->nobj=0;
	wgtflag=ctl_para->wgtflag;
	numflag=0;
	ncon=ctl_para->num_criteria;
	nparts=ctl_para->num_repartition;
/*	graph->vwgt=(int *)calloc(nvtxs, sizeof(int));
	if(graph->vwgt==NULL) {
		fprintf(stderr, "There is no enough memory for vwgt\n");
        HECMW_Finalize();
		exit(0);
	}
	for(i=0;i<nvtxs;i++)
		graph->vwgt[i]=1;
*/	
/*	tpwgts=(float *)calloc(pesize, sizeof(float));
	if(tpwgts==NULL) {
		fprintf(stderr, "There is no enough memory for tpwgts\n");
        HECMW_Finalize();
		exit(0);
	}
    for(i=0;i<pesize;i++)
		tpwgts[i]=ctl_para;
		
	ubvec=(float *)calloc(ncon, sizeof(float));
	if(ubvec==NULL) 
		HECMW_dlb_memory_exit("ubvec");
	for(i=0;i<ncon;i++)
		ubvec[i]=1.05;
	*/
	options[0]=1;
	options[1]=3;
	options[2]=1;
        options[3]=1;
	edgecut=0;
	result->part=(int *)calloc(nvtxs, sizeof(int));
	if(result->part==NULL) 
		HECMW_dlb_memory_exit("result: part");

/*	xyz=(float *)calloc(mesh->n_internal*3, sizeof(float));
	for(i=0;i<mesh->n_internal;i++) {
		xyz[i*3]=(float)mesh->node[i];
		xyz[i*3+1]=(float)mesh->node[i+mesh->n_node];
		xyz[i*3+2]=(float)mesh->node[i+2*mesh->n_node];
	}
	ndim=3;
    ParMETIS_V3_PartGeomKway(graph->vtxdist, graph->xadj, graph->adjncy, NULL, NULL, &wgtflag,
        &numflag, &ndim, xyz, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, result->part, &repart_comm);
	result->edgecut=edgecut;
*/
	if(strncmp(ctl_para->adaptive_repartition, "off", 3)==0) {
	if(ctl_para->wgtflag==0)
       ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, NULL, NULL, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		options, &edgecut, result->part, &repart_comm);
	else if((stat_para[6]!=0) && (stat_para[7]==0)){
		vwgts=(int *)calloc(mesh->nn_internal, sizeof(int));
		if(vwgts==NULL)
			HECMW_dlb_memory_exit("vwgts");
		fp_wgts=fopen(ctl_para->vwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: vwgt_filename\n");
			exit(0);
		}
		for(i=0;i<mesh->nn_internal;i++)
			fscanf(fp_wgts, "%d", &vwgts[i]);
		fclose(fp_wgts);
       ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, vwgts, NULL, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		options, &edgecut, result->part, &repart_comm);
		free(vwgts);
	   }
	else if ((stat_para[6]==0) && (stat_para[7]!=0)){
		adjwgts=(int *)calloc(graph->xadj[mesh->nn_internal], sizeof(int));
		if(adjwgts==NULL)
			HECMW_dlb_memory_exit("adjwgts");
		fp_wgts=fopen(ctl_para->adjwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: adjwgt_filename\n");
			exit(0);
		}
		for(i=0;i<graph->xadj[mesh->nn_internal];i++)
			fscanf(fp_wgts, "%d", &adjwgts[i]);
		fclose(fp_wgts);
       ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, NULL, adjwgts, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		options, &edgecut, result->part, &repart_comm);
		free(adjwgts);
	   }
	else if ((stat_para[6]!=0) && (stat_para[7]!=0)){
		vwgts=(int *)calloc(mesh->nn_internal, sizeof(int));
		if(vwgts==NULL)
			HECMW_dlb_memory_exit("vwgts");
		fp_wgts=fopen(ctl_para->vwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: vwgt_filename\n");
			exit(0);
		}
		for(i=0;i<mesh->nn_internal;i++)
			fscanf(fp_wgts, "%d", &vwgts[i]);
		fclose(fp_wgts);
		adjwgts=(int *)calloc(graph->xadj[mesh->nn_internal], sizeof(int));
		if(adjwgts==NULL)
			HECMW_dlb_memory_exit("adjwgts");
		fp_wgts=fopen(ctl_para->adjwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: adjwgt_filename\n");
			exit(0);
		}
		for(i=0;i<graph->xadj[mesh->nn_internal];i++)
			fscanf(fp_wgts, "%d", &adjwgts[i]);
		fclose(fp_wgts);
       ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, vwgts, adjwgts, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		options, &edgecut, result->part, &repart_comm);
		free(adjwgts);
		free(vwgts);
	   }

	result->edgecut=edgecut;
	}
	else if(strncmp(ctl_para->adaptive_repartition, "on", 2)==0) {
	if(ctl_para->wgtflag==0)
    ParMETIS_V3_AdaptiveRepart(graph->vtxdist, graph->xadj, graph->adjncy, NULL, NULL, NULL, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		&(ctl_para->itr_rate),  options, &edgecut, result->part, &repart_comm);
	

	else if((stat_para[6]!=0) && (stat_para[7]==0)){
		vwgts=(int *)calloc(mesh->nn_internal, sizeof(int));
		if(vwgts==NULL)
			HECMW_dlb_memory_exit("vwgts");
		fp_wgts=fopen(ctl_para->vwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: vwgt_filename\n");
			exit(0);
		}
		for(i=0;i<mesh->nn_internal;i++)
			fscanf(fp_wgts, "%d", &vwgts[i]);
		fclose(fp_wgts);
    ParMETIS_V3_AdaptiveRepart(graph->vtxdist, graph->xadj, graph->adjncy, vwgts, NULL, NULL, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		&(ctl_para->itr_rate),  options, &edgecut, result->part, &repart_comm);
	free(vwgts);
	   }
	else if ((stat_para[6]==0) && (stat_para[7]!=0)){
		adjwgts=(int *)calloc(graph->xadj[mesh->nn_internal], sizeof(int));
		if(adjwgts==NULL)
			HECMW_dlb_memory_exit("adjwgts");
		fp_wgts=fopen(ctl_para->adjwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: adjwgt_filename\n");
			exit(0);
		}
		for(i=0;i<graph->xadj[mesh->nn_internal];i++)
			fscanf(fp_wgts, "%d", &adjwgts[i]);
		fclose(fp_wgts);
    ParMETIS_V3_AdaptiveRepart(graph->vtxdist, graph->xadj, graph->adjncy, NULL, NULL, adjwgts, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		&(ctl_para->itr_rate),  options, &edgecut, result->part, &repart_comm);
		free(adjwgts);
	   }
	else if ((stat_para[6]!=0) && (stat_para[7]!=0)){
		vwgts=(int *)calloc(mesh->nn_internal, sizeof(int));
		if(vwgts==NULL)
			HECMW_dlb_memory_exit("vwgts");
		fp_wgts=fopen(ctl_para->vwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: vwgt_filename\n");
			exit(0);
		}
		for(i=0;i<mesh->nn_internal;i++)
			fscanf(fp_wgts, "%d", &vwgts[i]);
		fclose(fp_wgts);
		adjwgts=(int *)calloc(graph->xadj[mesh->nn_internal], sizeof(int));
		if(adjwgts==NULL)
			HECMW_dlb_memory_exit("adjwgts");
		fp_wgts=fopen(ctl_para->adjwgt_filename, "r");
		if(fp_wgts==NULL) {
			fprintf(stderr, "control file wrong: adjwgt_filename\n");
			exit(0);
		}
		for(i=0;i<graph->xadj[mesh->nn_internal];i++)
			fscanf(fp_wgts, "%d", &adjwgts[i]);
		fclose(fp_wgts);
    ParMETIS_V3_AdaptiveRepart(graph->vtxdist, graph->xadj, graph->adjncy, vwgts, NULL, adjwgts, &(ctl_para->wgtflag),
        &numflag, &(ctl_para->num_criteria), &(ctl_para->num_repartition), ctl_para->machine_wgt, ctl_para->balance_rate,
		&(ctl_para->itr_rate),  options, &edgecut, result->part, &repart_comm);
		free(adjwgts);
		free(vwgts);
	   }


	result->edgecut=edgecut;
	}
	/* write test file */
/*    for(i=0;i<mesh->n_internal;i++) {
		fprintf(fp_test, "%d \n", result->part[i]);
	}
	
	fclose(fp_test);
*/	
	free(graph->xadj);
	free(graph->adjncy);

	return;
	}
	
	














		   







		





