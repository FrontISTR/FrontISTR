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



#ifndef HECMW_DIST_INCLUDED
#define HECMW_DIST_INCLUDED

#include "hecmw_struct.h"


extern int HECMW_dist_get_mat_id(const struct hecmwST_material *mat, const char *name);


extern int HECMW_dist_get_ngrp_id(const struct hecmwST_node_grp *ngrp, const char *name);


extern int HECMW_dist_get_egrp_id(const struct hecmwST_elem_grp *egrp, const char *name);


extern int HECMW_dist_get_sgrp_id(const struct hecmwST_surf_grp *sgrp, const char *name);


extern int HECMW_dist_gid2lid_node(const struct hecmwST_local_mesh *mesh, int gid);


extern int HECMW_dist_gid2lid_elem(const struct hecmwST_local_mesh *mesh, int gid);

#endif
