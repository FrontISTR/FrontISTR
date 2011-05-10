/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#ifndef INC_HECMW_COUPLE_INTER_IFTABLE
#define INC_HECMW_COUPLE_INTER_IFTABLE

#include "hecmw_struct.h"
#include "hecmw_couple_comm.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"
#include "hecmw_couple_background_cell.h"
#include "hecmw_couple_mapped_point.h"


struct hecmw_couple_inter_iftable {
	int n_neighbor_pe_import;	
	int *neighbor_pe_import;    
	int *import_index;			
	int *import_item;			
	int n_neighbor_pe_export;	
	int *neighbor_pe_export;	
	int *export_index;			
	int *export_item;			
};

extern void
HECMW_couple_free_inter_iftable(struct hecmw_couple_inter_iftable *p);
extern struct hecmw_couple_inter_iftable *
HECMW_couple_alloc_inter_iftable(void);
extern void
HECMW_couple_print_inter_iftable(const struct hecmw_couple_inter_iftable *p, FILE *fp);

extern struct hecmw_couple_inter_iftable *
HECMW_couple_set_map_data(
		const struct hecmwST_local_mesh *mesh_src,
		const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_bounding_box *bbox_dst,
		const struct hecmw_couple_background_cell *bgcell_src,
		const struct hecmw_couple_mapped_point *mapped_point);

#endif	/* INC_HECMW_COUPLE_INTER_IFTABLE */
