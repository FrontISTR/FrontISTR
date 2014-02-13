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



#ifndef HECMW_LIB_FC_INCLUDED
#define HECMW_LIB_FC_INCLUDED


extern char *HECMW_strcpy_f2c(const char *fstr, int flen);


extern char *HECMW_strcpy_f2c_r(const char *fstr, int flen, char *buf, int bufsize);


extern int HECMW_strcpy_c2f(const char *cstr, char *fstr, int flen);

#endif
