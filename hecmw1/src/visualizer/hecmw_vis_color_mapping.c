/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Visualization                                     *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_vis_color_mapping.h"

#include <math.h>

#define EPSILON	 	0.00000001


void value_to_rgb(double value, double color[3], double mincolor, double maxcolor, int color_mapping_style,
		double *interval_point, int interval_mapping_num, int color_system_type) {
	double r, g, b;
	int i;
	if(color_mapping_style==1) {
		if(fabs(maxcolor-mincolor)>EPSILON)
			value=(value-mincolor)/(maxcolor-mincolor);
	}

	if(color_mapping_style==2) {
		mincolor=interval_point[0];
		maxcolor=interval_point[1];
		if(fabs(maxcolor-mincolor)>EPSILON)
			value=(value-mincolor)/(maxcolor-mincolor);
	}
	if((color_mapping_style==3) || (color_mapping_style==4)) {
		if(value<interval_point[0])
			value=0.0;
		else if(value>interval_point[interval_mapping_num*2])
			value=1.0;
		else {
			for(i=1;i<interval_mapping_num+1;i++) {
				if((value<=interval_point[i*2]) && (value>interval_point[(i-1)*2])) {
					value=(value-interval_point[(i-1)*2])/(interval_point[i*2]-interval_point[(i-1)*2])*
					(interval_point[i*2+1]-interval_point[(i-1)*2+1])+interval_point[(i-1)*2+1];
				}
			}
		}
	}

	if(color_system_type==1) {

		if(value<0.0) value=0.0;
		if(value>1.0) value=1.0;
		if(value<=0.25) {
			r=0.0;
			g=value*4.0;
			b=1.0;
		}
		else if((value>0.25) && (value<=0.5)) {
			r=0.0;
			g=1.0;
			b=(0.5-value)*4.0;
		}
		else if((value>0.5) && (value<=0.75)) {
			r=(value-0.5)*4.0;
			g=1.0;
			b=0.0;
		}
		else if(value>0.75) {
			r=1.0;
			g=(1.0-value)*4.0;
			b=0.0;
		}

	}
	else if(color_system_type==2) {
		if(value<0.0)
			value=0.0;
		if(value>1.0)
			value=1.0;
		if(value<=0.2) {
			g=0.0; b=1.0; r=(0.2-value)*5.0;
		}
		else if((value>0.2) && (value<=0.4)) {
			r=0.0; b=1.0; g=(value-0.2)*5.0;
		}
		else if((value>0.4) && (value<=0.6)) {
			r=0.0; g=1.0; b=1.0-(value-0.4)*5.0;
		}
		else if((value>0.6) && (value<=0.8)) {
			r=(value-0.6)*5.0; g=1.0; b=0.0;
		}
		else if(value>0.0) {
			r=1.0; g=1.0-(value-0.8)*5.0;  b=0.0;
		}
	}
	else if(color_system_type==3) {
		r=g=b=value;
	}

	color[0]=r;
	color[1]=g;  color[2]=b;
	return;
}

