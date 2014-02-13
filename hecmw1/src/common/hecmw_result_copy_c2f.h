/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef HECMW_RESULT_COPY_C2F_INCLUDED
#define HECMW_RESULT_COPY_C2F_INCLUDED

#include "hecmw_result.h"

extern int HECMW_result_copy_c2f_init(struct hecmwST_result_data *result_data, int n_node, int n_elem);
extern int HECMW_result_copy_c2f_finalize(void);

#endif
