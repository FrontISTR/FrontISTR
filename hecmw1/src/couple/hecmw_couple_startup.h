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




#ifndef INC_HECMW_COUPLE_STARTUP
#define INC_HECMW_COUPLE_STARTUP


struct hecmw_couple_value {
	int n;			
	int item_type;	
	int n_dof;		
	int *item;		
	double *value;	
};


extern void
HECMW_couple_free_couple_value(struct hecmw_couple_value *couple_value);
extern struct hecmw_couple_value *
HECMW_couple_alloc_couple_value(void);
extern void
HECMW_couple_print_couple_value(const struct hecmw_couple_value *couple_value, FILE *fp);

extern struct hecmw_couple_value *
HECMW_couple_startup(const char *boundary_id);
extern void
HECMW_couple_cleanup(struct hecmw_couple_value *couple_value);

#endif	/* INC_HECMW_COUPLE_STARTUP */
