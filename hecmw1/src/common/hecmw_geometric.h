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



#ifndef HECMW_GEOMETRIC_INCLUDED
#define HECMW_GEOMETRIC_INCLUDED


struct hecmw_coord {
	double x;	
	double y;	
	double z;	
};


extern double HECMW_degree_to_radian(double deg);


extern double HECMW_radian_to_degree(double rad);


extern int HECMW_cylindrical_to_cartesian(
		const struct hecmw_coord *coord, struct hecmw_coord *result);

#endif
