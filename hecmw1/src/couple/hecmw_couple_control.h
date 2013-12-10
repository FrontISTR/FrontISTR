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




#ifndef INC_HECMW_COUPLE_CONTROL
#define INC_HECMW_COUPLE_CONTROL

#include <stdio.h>


struct hecmw_couple_ctrl_unit_ids {
	int n_unit;		
	char **ids;		
};


struct hecmw_couple_ctrl_couple_ids {
	int n_couple;	
	char **ids;		
};


struct hecmw_couple_ctrl_boundary_ids {
	int n_boundary;	
	char **ids;		
};


struct hecmw_couple_ctrl_proc {
	int n_proc;
	int is_specified_ranks;
	int *ranks;
};


struct hecmw_couple_group {
	int n_grp;
	int geom_type;
	int data_type;
	char **grp_name;
};

extern void
HECMW_couple_free_unit_ids(struct hecmw_couple_ctrl_unit_ids *unit_ids);
extern void
HECMW_couple_free_couple_ids(struct hecmw_couple_ctrl_couple_ids *couple_ids);
extern void
HECMW_couple_free_boundary_ids(struct hecmw_couple_ctrl_boundary_ids *boundary_ids);
extern void
HECMW_couple_ctrl_free_proc(struct hecmw_couple_ctrl_proc *proc_info);
extern void
HECMW_couple_ctrl_free_group(struct hecmw_couple_group *grp_info);

extern void
HECMW_couple_ctrl_free_couplemesh(void);
extern void
HECMW_couple_ctrl_free_couple(void);
extern void
HECMW_couple_ctrl_free(void);

extern void
HECMW_couple_ctrl_print_unit(FILE *fp);
extern void
HECMW_couple_ctrl_print_couple(FILE *fp);
extern void
HECMW_couple_ctrl_print_boundary(FILE *fp);

extern int
HECMW_couple_ctrl_unit(void);
extern int
HECMW_couple_ctrl_couple(void);
extern int
HECMW_couple_ctrl_boundary(void);

extern int
HECMW_couple_ctrl_get_n_unit(void);
extern int
HECMW_couple_ctrl_get_n_couple(void);
extern int
HECMW_couple_ctrl_get_n_boundary(void);
extern struct hecmw_couple_ctrl_unit_ids *
HECMW_couple_get_unit_ids(void);
extern struct hecmw_couple_ctrl_couple_ids *
HECMW_couple_get_couple_ids(void);
extern struct hecmw_couple_ctrl_boundary_ids *
HECMW_couple_get_boundary_ids(void);
extern char *
HECMW_couple_ctrl_get_unit_id(const char *couple_id, int unit_specifier, char *buf, int bufsize);
extern char *
HECMW_couple_ctrl_get_couple_id(const char *boundary_id, char *buf, int bufsize);
extern struct hecmw_couple_ctrl_proc *
HECMW_couple_ctrl_get_proc(const char *unit_id);
extern int
HECMW_couple_ctrl_get_type(const char *couple_id, int *couple_type);
extern int
HECMW_couple_ctrl_get_direction(const char *boundary_id, int *direction);
extern int
HECMW_couple_ctrl_get_tolerance(const char *boundary_id, double *tolerance);
extern int
HECMW_couple_ctrl_get_bbcoef(const char *boundary_id, double *bbcoef);
extern int
HECMW_couple_ctrl_get_bgcoef(const char *boundary_id, double *bgcoef);
extern struct hecmw_couple_group *
HECMW_couple_ctrl_get_group(const char *boundary_id, int unit_specifier);

#endif	/* INC_HECMW_COUPLE_CONTROL */
