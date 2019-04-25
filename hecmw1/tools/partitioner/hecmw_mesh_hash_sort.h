/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_MESH_HASH_SORT
#define INC_HECMW_MESH_HASH_SORT

/* edge */
extern int HECMW_mesh_hsort_edge_init(int n_node, int n_elem);
extern int HECMW_mesh_hsort_edge_realloc(void);
extern void HECMW_mesh_hsort_edge_final(void);
extern long long int HECMW_mesh_hsort_edge_get_n(void);
extern int *HECMW_mesh_hsort_edge_get_v(void);
extern long long int HECMW_mesh_hsort_edge(int node1, int node2);

/* triangular surface */
extern int HECMW_mesh_hsort_tsuf_init(int n_node, int n_elem);
extern int HECMW_mesh_hsort_tsuf_realloc(void);
extern void HECMW_mesh_hsort_tsuf_final(void);
extern int HECMW_mesh_hsort_tsuf_get_n(void);
extern int *HECMW_mesh_hsort_tsuf_get_v(void);
extern int HECMW_mesh_hsort_tsuf(int node1, int node2, int node3);

/* quadrilateral surface */
extern int HECMW_mesh_hsort_qsuf_init(int n_node, int n_elem);
extern int HECMW_mesh_hsort_qsuf_realloc(void);
extern void HECMW_mesh_hsort_qsuf_final(void);
extern int HECMW_mesh_hsort_qsuf_get_n(void);
extern int *HECMW_mesh_hsort_qsuf_get_v(void);
extern int HECMW_mesh_hsort_qsuf(int node1, int node2, int node3, int node4);

#endif /* INC_HECMW_MESH_HASH_SORT */
