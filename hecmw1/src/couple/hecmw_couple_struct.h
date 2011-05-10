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




#ifndef INC_HECMW_COUPLE_STRUCT
#define INC_HECMW_COUPLE_STRUCT

#include "hecmw_config.h"
#include "hecmw_struct.h"

struct hecmw_couple_comm {
  int psize;
  int rank;
  int *ranks;
  HECMW_Comm comm;
  HECMW_Group group;
  int root;
  int is_root;
  int is_member;
};

#endif	/* INC_HECMW_COUPLE_STRUCT */
