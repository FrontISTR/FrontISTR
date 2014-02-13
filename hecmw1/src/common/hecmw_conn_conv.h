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



#ifndef HECMW_CONN_CONV_INCLUDED
#define HECMW_CONN_CONV_INCLUDED

#define HECMW_CONNTYPE_HECMW 1 
#define HECMW_CONNTYPE_ABAQUS 2 
#define HECMW_CONNTYPE_NASTRAN 3 


extern int HECMW_convert_connectivity(int from, int hecmw_etype, int *conn);

#endif
