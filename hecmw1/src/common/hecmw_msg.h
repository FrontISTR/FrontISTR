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



#ifndef HECMW_MSG_INCLUDED
#define HECMW_MSG_INCLUDED

#include  "hecmw_msgno.h"


struct hecmw_msgent {
	int msgno;			
	char *msgno_str;	
	char *msg;			
};


extern struct hecmw_msgent hecmw_msg_table[];


extern char *HECMW_strmsg(int msgno);


extern int HECMW_is_syserr(int msgno);

#endif
