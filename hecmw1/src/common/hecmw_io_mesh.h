/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef HECMW_IO_MESH_INCLUDED
#define HECMW_IO_MESH_INCLUDED

#include <stdio.h>
#include "hecmw_util.h"
#include "hecmw_io_struct.h"
#include "hecmw_system.h"


extern int HECMW_io_get_version(void);


extern int HECMW_io_init(void);


extern int HECMW_io_finalize(void);


extern void HECMW_io_print_all(FILE *fp);


extern int HECMW_io_free_all(void);


extern int HECMW_io_set_gridfile(char *gridfile);


extern struct hecmw_io_amplitude *
HECMW_io_add_amp(const char *name, int definition, int time, int value, double val, double t);


extern struct hecmw_io_initial *HECMW_io_get_initial(int node);


extern struct hecmw_io_initial *
HECMW_io_add_initial(int type, int node, const char *ngrp , double val);


extern struct hecmw_io_element *HECMW_io_get_elem(int id);


extern int HECMW_io_get_n_elem(void);


extern int HECMW_io_get_elem_max_id(void);


extern struct hecmw_io_element *
HECMW_io_add_elem(int id, int type, int *node, int nmatitem, double *matitem);


struct hecmw_io_id_array *HECMW_io_get_elem_in_egrp(const char *name);


extern struct hecmw_io_egrp *HECMW_io_get_egrp(const char *name);


extern int HECMW_io_add_egrp(const char *name, int nelem, int *elem);


extern struct hecmw_io_node *HECMW_io_get_node(int id);


extern int HECMW_io_get_n_node(void);


extern struct hecmw_io_node *HECMW_io_add_node(int id, double x, double y, double z);


extern int HECMW_io_get_nnode_in_ngrp(const char *name);


extern int HECMW_io_remove_node(int id);


extern struct hecmw_io_ngrp *HECMW_io_get_ngrp(const char *name);


extern struct hecmw_io_id_array *HECMW_io_get_node_in_ngrp(const char *name);


extern int HECMW_io_add_ngrp(const char *name, int nnode, int *node);


extern int HECMW_io_add_sgrp(const char *name, int n, int *elem, int *surf);


extern struct hecmw_io_mpc *HECMW_io_add_mpc(int neq, const struct hecmw_io_mpcitem *mpcitem, double cnst);


extern struct hecmw_io_section *HECMW_io_add_sect(struct hecmw_io_section *sect);


extern struct hecmw_io_material *HECMW_io_get_mat(const char *name);


extern struct hecmw_io_material *HECMW_io_add_mat(struct hecmw_io_material *mat);


extern void HECMW_io_set_header(struct hecmw_io_header *header);


extern struct hecmw_system_param *HECMW_io_get_system(void);


extern void HECMW_io_set_system(struct hecmw_system_param *system);


extern void HECMW_io_set_zero(struct hecmw_io_zero *zero);


extern struct hecmw_io_contact *HECMW_io_add_contact(const char *name, int type, const char *slave_grp, const char *master_grp);

extern int HECMW_io_post_process(void);


extern int HECMW_io_pre_process(void);


extern int HECMW_io_check_mpc_dof(int dof);


extern int HECMW_io_is_reserved_name(const char *name);


extern struct hecmwST_local_mesh *HECMW_io_make_local_mesh(void);

#endif
