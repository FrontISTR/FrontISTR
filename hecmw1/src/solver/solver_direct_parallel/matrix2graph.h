#ifndef MATRIX2GRAPH_H
#define MATRIX2GRAPH_H

typedef struct {
	int nvtxs;
	int nedges;
	int ncon;
	int *xadj;
	int *adjncy;
} graph_type;

extern void matrix2graph(int num_of_row, int num_of_col, int num_of_nzero, int *irow, int *jcol, graph_type *graph);

#endif /* MATRIX2GRAPH_H */
