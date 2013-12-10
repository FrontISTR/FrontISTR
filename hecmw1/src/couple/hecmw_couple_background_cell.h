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




#ifndef INC_HECMW_COUPLE_BACKGROUND_CELL
#define INC_HECMW_COUPLE_BACKGROUND_CELL

#include "hecmw_struct.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"


struct hecmw_couple_background_cell {
	int n;			
	double coef;	
	int nx;			
	int ny;			
	int nz;			
	double dx;		
	double dy;		
	double dz;		
};

extern void
HECMW_couple_free_background_cell(struct hecmw_couple_background_cell *bgcell);

extern struct hecmw_couple_background_cell *
HECMW_couple_set_background_cell(const char *boundary_id,
		const struct hecmwST_local_mesh *mesh,
		const struct hecmw_couple_bounding_box *bbox,
		const struct hecmw_couple_boundary *boundary);

#endif	/* INC_HECMW_COUPLE_BACKGROUND_CELL */
