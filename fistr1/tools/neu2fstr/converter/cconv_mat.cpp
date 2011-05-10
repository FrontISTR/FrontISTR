/*=====================================================================*
 *                                                                     *
 *   Software Name : neu2fstr                                          *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : FrontSTR Input File Converter                     *
 *                                                                     *
 *            Written by Noboru Imai (Univ. of Tokyo)                  *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

/* matrix for converting ver.1.0 */


#include <assert.h>
#include "cconv_mat.h"


cconv_mat::cconv_mat( coord_type t )
 : type(t)
{
	unit();
}



cconv_mat::cconv_mat( const cconv_mat& m, coord_type t )
 : type(t)
{
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			e[i][j] = m(i,j);
		}
	}
}

cconv_mat& cconv_mat::operator=( const cconv_mat& m )
{
	type = m.type;
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			e[i][j] = ((cconv_mat&)m)(i,j);
		}
	}
	return *this;
}

cconv_mat& cconv_mat::operator*=( const cconv_mat& m )
{
	type = m.type;
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			e[i][j] *= ((cconv_mat&)m)(i,j);
		}
	}
	return *this;
}


void cconv_mat::zero()
{
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			e[i][j] = 0.0;
		}
	}
}


void cconv_mat::unit()
{
	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			if(i==j)
				e[i][j] = 1.0;
			else
				e[i][j] = 0.0;
		}
	}
}


void cconv_mat::transfer( double x, double y, double z )
{
	unit();
	e[0][3] = x;
	e[1][3] = y;
	e[2][3] = z;
}


void cconv_mat::rotate( char axis, double angle )
{
	double s = sin(angle);
	double c = cos(angle);

	e[3][0] = e[3][1] = e[3][2] = 0.0;
	e[3][3] = 1.0;

	switch(axis){
	case 'x':case 'X':
			e[0][0] = 1;	e[0][1] = 0;	e[0][2] = 0;
			e[1][0] = 0;	e[1][1] = c;	e[1][2] =-s;
			e[2][0] = 0;	e[2][1] = s;	e[2][2] = c;
			break;
	case 'y':case 'Y':
			e[0][0] = c;	e[0][1] = 0;	e[0][2] = s;
			e[1][0] = 0;	e[1][1] = 1;	e[1][2] = 0;
			e[2][0] =-s;	e[2][1] = 0;	e[2][2] = c;
			break;
	case 'z':case 'Z':
			e[0][0] = c;	e[0][1] =-s;	e[0][2] = 0;
			e[1][0] = s;	e[1][1] = c;	e[1][2] = 0;
			e[2][0] = 0;	e[2][1] = 0;	e[2][2] = 1;
			break;
	}
}


void cconv_mat::convert( double x, double y, double z,  double& X, double& Y, double& Z  )
{
	double xx, yy, zz;

	switch( type ){
	case coord_t_cartesian:
		cartesian_convert( x, y, z, X,Y,Z );
		break;
	case coord_t_cylinder:
		cylinder2cartesian( x, y, z, xx,yy,zz );
		cartesian_convert( xx, yy, zz, X,Y,Z );
		break;
	case coord_t_sphere:
		sphere2cartesian( x, y, z, xx,yy,zz );
		cartesian_convert( xx, yy, zz, X,Y,Z );
		break;
	default:
		assert(0);
	}
}


void cconv_mat::cartesian_convert( double x, double y, double z,  double& X, double& Y, double& Z  )
{
	X = e[0][0] * x + e[0][1] * y + e[0][2] * z + e[0][3];
	Y = e[1][0] * x + e[1][1] * y + e[1][2] * z + e[1][3];
	Z = e[2][0] * x + e[2][1] * y + e[2][2] * z + e[2][3];
}


void cconv_mat::cylinder2cartesian( double r, double a, double h, double& x, double& y, double& z )
{
	x = r * cos(a);
	y = r * sin(a);
	z = h;
}


void cconv_mat::sphere2cartesian( double r, double a, double b, double& x, double& y, double& z )
{
	double R = r * cos(b);
	x = R * cos(a);
	y = R * sin(a);
	z = r * sin(b);
}

cconv_mat operator * ( cconv_mat& a, cconv_mat& b )
{
	cconv_mat m;

	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			double r = 0.0;
			for(int k=0; k<4; k++) {
				r += a(i,k)* b(k,j);
			}
			m(i,j) = r;
		}
	}
	return m;
}

