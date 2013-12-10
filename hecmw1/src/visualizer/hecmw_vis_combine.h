#ifndef HECMW_VIS_COMBINE_H_INCLUDED
#define HECMW_VIS_COMBINE_H_INCLUDED

#include <stdio.h>
#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_vis_SF_geom.h"

void HECMW_vis_combine(struct surface_module *sf, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int tvertex, int tpatch, int *color_list, double *minvalue, double *maxvalue,
		Result *result, char *outfile,  HECMW_Comm VIS_COMM);
void put_neutral_head(FILE *outfp);
void put_neutral_601(FILE *outfp, struct hecmwST_local_mesh *mesh);
void put_neutral_402(FILE *outfp, struct hecmwST_local_mesh *mesh);
void put_neutral_middle(FILE *outfp);
void put_neutral_409(FILE *outfp);

#endif /* HECMW_VIS_COMBINE_H_INCLUDED */

