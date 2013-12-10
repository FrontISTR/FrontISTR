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




#ifndef INC_HECMW_COUPLE_INTRA_IFTABLE
#define INC_HECMW_COUPLE_INTRA_IFTABLE

#include "hecmw_struct.h"
#include "hecmw_couple_comm.h"
#include "hecmw_couple_boundary_info.h"


struct hecmw_couple_intra_iftable {
	int n_neighbor_pe;	
	int *neighbor_pe;	
	int *import_index;	
	int *import_item;	
	int *export_index;	
	int *export_item;	
};


extern void
HECMW_couple_free_intra_iftable(struct hecmw_couple_intra_iftable *intra_tbl);

extern struct hecmw_couple_intra_iftable *
HECMW_couple_alloc_intra_iftable(void);

extern struct hecmw_couple_intra_iftable *
HECMW_couple_make_intra_iftable(const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_boundary *boundary,
		const struct hecmw_couple_comm *intracomm);

#endif	/* INC_HECMW_COUPLE_INTRA_IFTABLE */
