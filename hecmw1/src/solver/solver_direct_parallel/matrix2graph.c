#include "matrix2graph.h"

#include <stdio.h>
#include <stdlib.h>
#include "separator.h"

void matrix2graph(int num_of_row, int num_of_col, int num_of_nzero, int *irow, int *jcol, graph_type *graph)
{
	int i;
	int *count, tmp_nodeid;

	graph->xadj=(int *)calloc(num_of_col+1, sizeof(int));
	if(graph->xadj==NULL)
		separator_memory_exit("matrix2graph: graph");
	graph->xadj[0]=0;
	for(i=0;i<num_of_nzero;i++) {
			if(irow[i]!=jcol[i]) {
		graph->xadj[irow[i]]++;
		graph->xadj[jcol[i]]++;
			}
	}
	for(i=0;i<num_of_col;i++) 
		graph->xadj[i+1]+=graph->xadj[i];
	graph->adjncy=(int *)calloc(graph->xadj[num_of_col],sizeof(int));
	if(graph->adjncy==NULL)
		separator_memory_exit("matrix2graph: graph");
	count=(int *)calloc(num_of_row, sizeof(int));
	if(count==NULL)
		separator_memory_exit("tmp: count");
	for(i=0;i<num_of_row;i++)
		count[i]=0;
	for(i=0;i<num_of_nzero;i++) {
			if(irow[i]!=jcol[i]) {
				tmp_nodeid=irow[i]-1;
				graph->adjncy[graph->xadj[tmp_nodeid]+count[tmp_nodeid]]=jcol[i]-1;
				count[tmp_nodeid]++;
				tmp_nodeid=jcol[i]-1;
				graph->adjncy[graph->xadj[tmp_nodeid]+count[tmp_nodeid]]=irow[i]-1;
				count[tmp_nodeid]++;
			}
		}
	free(count);
	graph->nvtxs=num_of_col;
	graph->nedges=graph->xadj[num_of_col];
	graph->ncon=0; 

	return;
}



