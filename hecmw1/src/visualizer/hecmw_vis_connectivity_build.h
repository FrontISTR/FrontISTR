#ifndef HECMW_VIS_CONNECTIVITY_BUILD_H_INCLUDED
#define HECMW_VIS_CONNECTIVITY_BUILD_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_vis_SF_geom.h"


void find_index_connectivity(struct hecmwST_local_mesh *mesh, int *index_connect);
void find_index_a_connect(struct hecmwST_local_mesh *mesh, int num_export, int pe_no, int *export_element,
		int *index_a_connect);
void generate_face(int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24]);
void add_to_hash(int elemID, int faceID, int hashID, Hash_table *h_table);
void build_hash_table(struct hecmwST_local_mesh *mesh, int *index_connect, Hash_table  *h_table,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24]);
int is_equal_array(int n[4],int nn[4], int num);
int is_connect(int elemID1, int faceID1, int elemID2, int faceID2, struct hecmwST_local_mesh *mesh,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],int index_face_hexa[7],
		int face_hexa[24]);
int find_to_hash(int elemID1,int faceID1, int hashID,  Hash_table *h_table, struct hecmwST_local_mesh *mesh, int tmp_connect[2],
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],int index_face_hexa[7],
		int face_hexa[24]);
void h_free(Hash_table *h_table,int maxadd);
void free_b_patch(Boundary_patch *b_patch);
void  build_connectivity(struct hecmwST_local_mesh *mesh, Hash_table *h_table, int *index_connect, int *connect,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24]);
void add_one_patch(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *node_hit, Tetra_point  *b_point, Head_patch_tetra  *head_b_patch, int node[4], int c_base, int d_base,
		int tn_component);
void add_two_patch(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *node_hit, Tetra_point  *b_point, Head_patch_tetra  *head_b_patch, int node[4], int c_base, int d_base,
		int tn_component);
void HECMW_vis_find_boundary_surface(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *tvertex, int *tpatch,
		double *minc, double *maxc, Result *result, int sf_i, HECMW_Comm VIS_COMM, int init_flag,
		Connect_inf *global_connect);

#endif /* HECMW_VIS_CONNECTIVITY_BUILD_H_INCLUDED */
