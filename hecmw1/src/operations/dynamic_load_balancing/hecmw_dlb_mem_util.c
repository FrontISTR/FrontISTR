/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Dynamic Load Balancing                            *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_repart.h"
void HECMW_dlb_memory_exit(char *var) 
{
	fprintf(stderr, "#### HEC-MW-VIS-E0001:There is no enough memory allocated for variable %s\n", var);
    HECMW_Finalize();
    exit(0);
	}

void HECMW_dlb_print_exit(char *var) 
{
	fprintf(stderr, "%s\n", var);
    HECMW_Finalize();
    exit(0);
	}


