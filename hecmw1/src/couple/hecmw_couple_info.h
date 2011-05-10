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




#ifndef INC_HECMW_COUPLE_INFO
#define INC_HECMW_COUPLE_INFO

#include "hecmw_config.h"
#include "hecmw_couple_struct.h"

extern void
HECMW_couple_free_comm(struct hecmw_couple_comm *comm);
extern void
HECMW_couple_free_couple_info(void);
extern int
HECMW_couple_comm_init(void);

extern char *
HECMW_couple_get_unit_id(const char *boundary_id, int unit_specifier, char *buf, int bufsize);
extern int
HECMW_couple_is_member(const char *boundary_id);
extern int
HECMW_couple_is_unit_member(const char *boundary_id, int unit_specifier);
extern int
HECMW_couple_is_unit_member_u(const char *unit_id);
extern int
HECMW_couple_is_root(const char *boundary_id);
extern int
HECMW_couple_is_unit_root(const char *boundary_id, int unit_specifier);
extern int
HECMW_couple_is_unit_root_u(const char *unit_id);
extern int
HECMW_intercomm_get_size(const char *boundary_id);
extern int
HECMW_intracomm_get_size(const char *boundary_id, int unit_specifier);
extern int
HECMW_intracomm_get_size_u(const char *unit_id);
extern int
HECMW_intercomm_get_rank(const char *boundary_id);
extern int
HECMW_intracomm_get_rank(const char *boundary_id, int unit_specifier);
extern int
HECMW_intracomm_get_rank_u(const char *unit_id);
extern HECMW_Comm
HECMW_intercomm_get_comm(const char *boundary_id);
extern HECMW_Comm
HECMW_intracomm_get_comm(const char *boundary_id, int unit_specifier);
extern HECMW_Comm
HECMW_intracomm_get_comm_u(const char *unit_id);
extern HECMW_Group
HECMW_intercomm_get_group(const char *boundary_id);
extern HECMW_Group
HECMW_intracomm_get_group(const char *boundary_id, int unit_specifier);
extern HECMW_Group
HECMW_intracomm_get_group_u(const char *unit_id);
extern struct hecmw_couple_comm *
HECMW_couple_get_intracomm(const char *boundary_id, int unit_specifier);
extern struct hecmw_couple_comm *
HECMW_couple_get_intracomm_u(const char *unit_id);
extern struct hecmw_couple_comm *
HECMW_couple_get_intercomm(const char *boundary_id);

#endif	/* INC_HECMW_COUPLE_INFO */
