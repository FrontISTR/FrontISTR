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




#ifndef INC_HECMW_COUPLE_MAPPED_NODE
#define INC_HECMW_COUPLE_MAPPED_NODE

#include "hecmw_struct.h"
#include "hecmw_couple_boundary_info.h"


struct hecmw_couple_mapped_point {
	int n;			
	int type;		
	int *item;		
	int *id;		
	double *coord;	
};

extern struct hecmw_couple_mapped_point *
HECMW_couple_alloc_mapped_point(void);
extern void
HECMW_couple_free_mapped_point(struct hecmw_couple_mapped_point *mapped_point);
extern struct hecmw_couple_mapped_point *
HECMW_couple_set_mapped_point(const char *boundary_id, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_boundary *boundary_dst);

#endif	/* INC_HECMW_COUPLE_MAPPED_NODE */
