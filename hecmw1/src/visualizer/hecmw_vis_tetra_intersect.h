#ifndef HECMW_VIS_TETRA_INTERSECT_H_INCLUDED
#define HECMW_VIS_TETRA_INTERSECT_H_INCLUDED

#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_vis_SF_geom.h"

extern int get_tetra_data(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int elemID, Tetra *tetra, int tn_component);
extern void find_intersection_tetra(Tetra *tetra, double isovalue, Tetra_point *tetra_point, Head_patch_tetra *head_patch_tetra,
		Hash_vertex *vertex_hash_table);

#endif /* HECMW_VIS_TETRA_INTERSECT_H_INCLUDED */







