/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/* matrix for converting ver.1.0 */


#ifndef cconv_matH
#define cconv_matH

#include <math.h>



enum coord_type {
	coord_t_cartesian,
	coord_t_cylinder,
	coord_t_sphere
};

class cconv_mat {
public:
	double e[4][4];
	coord_type type;

	cconv_mat( coord_type t=coord_t_cartesian );
	cconv_mat( const cconv_mat& m, coord_type t=coord_t_cartesian );
	cconv_mat& operator=( const cconv_mat& m );
	cconv_mat& operator*=( const cconv_mat& m );

	void zero();
	void unit();
	void transfer( double x, double y, double z );
	void rotate( char axis, double angle );

	void convert( double x, double y, double z,  double& X, double& Y, double& Z  );

	double& operator()(int i, int j) { return e[i][j];}
	double operator()(int i, int j) const { return e[i][j];}

	void cartesian_convert( double x, double y, double z,  double& X, double& Y, double& Z  );
	void cylinder2cartesian( double r, double a, double h, double& x, double& y, double& z );
	void sphere2cartesian( double r, double a, double b, double& x, double& y, double& z );
};


cconv_mat operator * ( cconv_mat& a, cconv_mat& b );


#endif
