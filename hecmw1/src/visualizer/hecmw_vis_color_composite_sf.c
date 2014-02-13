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

#include "hecmw_vis_color_composite_sf.h"

#include <math.h>
#include "hecmw_vis_SF_geom.h"

int find_surface_point(double fp[3][3], double point_o[3], double view_point_d[3],double point[3], double n_f[4],
		double v_f[9], double c_value[3], double *value, double normal[9], int normal_flag, int smooth_shading)
{
	int i,j,n1, n2, flag_inside;
	double s[3], sina, cosa, ss, t, dis[3], dist, n_norm, v_norm, nv_sign;
	double  tmp_p[7];

	/*find the equation of the patch first */
	flag_inside=0;
	if(normal_flag==1) {
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
	}
	for(i=0;i<3;i++)
		v_f[i]=view_point_d[i]-point_o[i];
	v_norm=sqrt(v_f[0]*v_f[0]+v_f[1]*v_f[1]+v_f[2]*v_f[2]);
	if(fabs(v_norm)>EPSILON) {
		for(i=0;i<3;i++)
			v_f[i]/=v_norm;
	}
	nv_sign=n_f[0]*v_f[0]+n_f[1]*v_f[1]+n_f[2]*v_f[2];
	if(nv_sign<0.0) {
		for(i=0;i<3;i++)
			n_f[i]=-n_f[i];
	}
	if(normal_flag==1)
		n_f[3]=-n_f[0]*fp[0][0]-n_f[1]*fp[0][1]-n_f[2]*fp[0][2];
	if(smooth_shading==0)
		for(i=0;i<3;i++)
			normal[i]=n_f[i];
	/*find intersection point*/
	if(fabs(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+
			n_f[2]*(point_o[2]-view_point_d[2]))>EPSILON) {
		t=(-n_f[3]-n_f[0]*view_point_d[0]-n_f[1]*view_point_d[1]-n_f[2]*view_point_d[2])
		/(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+n_f[2]*(point_o[2]-view_point_d[2]));
		if(t>-EPSILON) {
			for(j=0;j<3;j++)
				point[j]=view_point_d[j]+t*(point_o[j]-view_point_d[j]);
			/*judge whether it is inside the range of patch by equal area method*/
			for(j=0;j<3;j++) {
				n1=j;
				n2=j+1;
				if(n2==3) n2=0;
				if((fabs(point[0]-fp[n1][0])<EPSILON) && (fabs(point[1]-fp[n1][1])<EPSILON) && (fabs(point[2]-fp[n1][2])<EPSILON)) {
					flag_inside=1;
					*value=c_value[n1];
					if(smooth_shading==1) {
						tmp_p[0]=normal[n1*3+0];
						tmp_p[1]=normal[n1*3+1];
						tmp_p[2]=normal[n1*3+2];
						normal[0]=tmp_p[0];
						normal[1]=tmp_p[1];
						normal[2]=tmp_p[2];
					}
					break;
				}
				cosa=((point[0]-fp[n1][0])*(fp[n2][0]-fp[n1][0])+(point[1]-fp[n1][1])*(fp[n2][1]-fp[n1][1])+
						(point[2]-fp[n1][2])*(fp[n2][2]-fp[n1][2]))/(sqrt(SQR(point[0]-fp[n1][0])+SQR(point[1]-fp[n1][1])
								+SQR(point[2]-fp[n1][2]))*sqrt(SQR(fp[n2][0]-fp[n1][0])+SQR(fp[n2][1]-fp[n1][1])+SQR(fp[n2][2]-fp[n1][2])));
				sina=sqrt(1-cosa*cosa);
				s[j]=sqrt(SQR(fp[n2][0]-fp[n1][0])+SQR(fp[n2][1]-fp[n1][1])+SQR(fp[n2][2]-fp[n1][2]))*sqrt(SQR(point[0]-fp[n1][0])+SQR(point[1]-fp[n1][1])
						+SQR(point[2]-fp[n1][2]))*sina/2.0;
			}
			if(flag_inside==0) {
				cosa=((fp[2][0]-fp[0][0])*(fp[1][0]-fp[0][0])+(fp[2][1]-fp[0][1])*(fp[1][1]-fp[0][1])+
						(fp[2][2]-fp[0][2])*(fp[1][2]-fp[0][2]))/(sqrt(SQR(fp[2][0]-fp[0][0])+SQR(fp[2][1]-fp[0][1])
								+SQR(fp[2][2]-fp[0][2]))*sqrt(SQR(fp[1][0]-fp[0][0])+SQR(fp[1][1]-fp[0][1])+SQR(fp[1][2]-fp[0][2])));
				sina=sqrt(1-cosa*cosa);
				ss=sqrt(SQR(fp[1][0]-fp[0][0])+SQR(fp[1][1]-fp[0][1])+SQR(fp[1][2]-fp[0][2]))*sqrt(SQR(fp[2][0]-fp[0][0])+SQR(fp[2][1]-fp[0][1])
						+SQR(fp[2][2]-fp[0][2]))*sina/2.0;
				/*		 if(fabs(ss-s[0]-s[1]-s[2])<EPSILON*ss/100.0)
				 */
				if(fabs(ss-s[0]-s[1]-s[2])<EPSILON) {
					flag_inside=1;
					for(j=0;j<3;j++)
						dis[j]=sqrt(SQR(point[0]-fp[j][0])+SQR(point[1]-fp[j][1])+SQR(point[2]-fp[j][2]));
					dist=1.0/dis[0]+1.0/dis[1]+1.0/dis[2];
					*value=c_value[0]/(dis[0]*dist)+c_value[1]/(dis[1]*dist)+c_value[2]/(dis[2]*dist);
					if(smooth_shading==1) {
						for(j=0;j<3;j++) {
							tmp_p[j]=normal[0*3+j]/(dis[0]*dist)+normal[1*3+j]/(dis[1]*dist)+normal[2*3+j]/(dis[2]*dist);
						}
						for(j=0;j<3;j++)
							normal[j]=tmp_p[j];

					}

				}
			}
		}
	}

	return(flag_inside);
}

int find_point_depth(double fp[3][3], double point_o[3], double view_point_d[3], double n_f[4], double point[3],
		int normal_flag)
{
	int i,j, flag_inside;
	double v_f[3], t, n_norm, v_norm, nv_sign;

	/*find the equation of the patch first */
	flag_inside=0;
	if(normal_flag==1) {
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
	}
	for(i=0;i<3;i++)
		v_f[i]=view_point_d[i]-point_o[i];
	v_norm=sqrt(v_f[0]*v_f[0]+v_f[1]*v_f[1]+v_f[2]*v_f[2]);
	if(fabs(v_norm)>EPSILON) {
		for(i=0;i<3;i++)
			v_f[i]/=v_norm;
	}
	nv_sign=n_f[0]*v_f[0]+n_f[1]*v_f[1]+n_f[2]*v_f[2];
	if(nv_sign<0.0) {
		for(i=0;i<3;i++)
			n_f[i]=-n_f[i];
	}
	if(normal_flag==1)
		n_f[3]=-n_f[0]*fp[0][0]-n_f[1]*fp[0][1]-n_f[2]*fp[0][2];

	/*find intersection point*/
	if(fabs(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+
			n_f[2]*(point_o[2]-view_point_d[2]))>EPSILON) {
		t=(-n_f[3]-n_f[0]*view_point_d[0]-n_f[1]*view_point_d[1]-n_f[2]*view_point_d[2])
		/(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+n_f[2]*(point_o[2]-view_point_d[2]));
		if(t>-EPSILON) {
			for(j=0;j<3;j++)
				point[j]=view_point_d[j]+t*(point_o[j]-view_point_d[j]);
			flag_inside=1;

		}
	}

	return(flag_inside);
}


int find2_point_depth(double fp[3][3], double point_o[3], double view_point_d[3], double point[3])
{
	int i,j,flag_inside;
	double v_f[3], t, n_norm, v_norm, nv_sign;
	double  n_f[4];

	/*find the equation of the patch first */
	flag_inside=0;
	n_f[0]=  (fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
	n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
	n_f[2]=  (fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
	n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
	if(fabs(n_norm)>EPSILON) {
		for(j=0;j<3;j++)
			n_f[j]/=n_norm;
	}
	for(i=0;i<3;i++)
		v_f[i]=view_point_d[i]-point_o[i];
	v_norm=sqrt(v_f[0]*v_f[0]+v_f[1]*v_f[1]+v_f[2]*v_f[2]);
	if(fabs(v_norm)>EPSILON) {
		for(i=0;i<3;i++)
			v_f[i]/=v_norm;
	}
	nv_sign=n_f[0]*v_f[0]+n_f[1]*v_f[1]+n_f[2]*v_f[2];
	if(nv_sign<0.0) {
		for(i=0;i<3;i++)
			n_f[i]=-n_f[i];
	}
	n_f[3]=-n_f[0]*fp[0][0]-n_f[1]*fp[0][1]-n_f[2]*fp[0][2];
	/*find intersection point*/
	if(fabs(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+
			n_f[2]*(point_o[2]-view_point_d[2]))>EPSILON) {
		t=(-n_f[3]-n_f[0]*view_point_d[0]-n_f[1]*view_point_d[1]-n_f[2]*view_point_d[2])
		/(n_f[0]*(point_o[0]-view_point_d[0])+n_f[1]*(point_o[1]-view_point_d[1])+n_f[2]*(point_o[2]-view_point_d[2]));
		if(t>-EPSILON) {
			for(j=0;j<3;j++)
				point[j]=view_point_d[j]+t*(point_o[j]-view_point_d[j]);
			flag_inside=1;

		}
	}

	return(flag_inside);
}


void compute_color_sf(double p[3], double value, double n[3], int color_mapping_style, double *interval_point,
		double view_point_d[3],  int interval_mapping_num, int color_system_type, int num_of_lights,
		double *light_point, double k_ads[3], double accum_rgba[4],
		double mincolor, double maxcolor, int display_method)
{
	int i, j;
	double cosalpha, costheta;
	double color[3];
	double lp[3], vp[3], lp_norm, vp_norm, norm, hp[3], hp_norm;
	double inprodLN, inprodVN, inprodHN;
	double r, g, b;




	/*--------------map value to rgb -------------------   */
	if(display_method!=4) {
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
						break;
					}
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
	/* --------------end of mapping value to rgb -----------------------  */

	r=0.0;   g=0.0;    b=0.0;
	for(j=0;j<num_of_lights;j++) {

		for(i=0;i<3;i++) {
			lp[i]=light_point[j*3+i]-p[i];
			vp[i]=view_point_d[i]-p[i];
			hp[i]=(lp[i]+vp[i])/2.0;
		}
		lp_norm=sqrt(SQR(lp[0])+SQR(lp[1])+SQR(lp[2]));
		vp_norm=sqrt(SQR(vp[0])+SQR(vp[1])+SQR(vp[2]));
		hp_norm=sqrt(SQR(hp[0])+SQR(hp[1])+SQR(hp[2]));
		if(fabs(lp_norm)>EPSILON) {
			for(i=0;i<3;i++)
				lp[i]/=lp_norm;
		}
		if(fabs(vp_norm)>EPSILON) {
			for(i=0;i<3;i++)
				vp[i]/=vp_norm;
		}
		if(fabs(hp_norm)>EPSILON) {
			for(i=0;i<3;i++)
				hp[i]/=hp_norm;
		}
		norm=sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
		if(fabs(norm)>EPSILON) {
			for(i=0;i<3;i++)
				n[i]/=norm;
		}
		inprodLN=n[0]*lp[0]+n[1]*lp[1]+n[2]*lp[2];
		inprodVN=n[0]*vp[0]+n[1]*vp[1]+n[2]*vp[2];

		inprodHN=n[0]*hp[0]+n[1]*hp[1]+n[2]*hp[2];
		/*	a_current=opacity_decision(in_voxel, vr, vd, grad_minmax, feap_minmax, feai_minmax,
		dis_minmax, opa_table, mincolor, maxcolor, time_step);
	cosalpha=sqrt(1.0-pow(inprodLN,2));
		 */
		cosalpha=inprodLN;
		costheta=inprodLN*inprodVN-sqrt(1.0-inprodLN*inprodLN)*sqrt(1.0-inprodVN*inprodVN);


		cosalpha=fabs(cosalpha);


		if(cosalpha>0) {

			r+=color[0]*(k_ads[0]+k_ads[1]*cosalpha+k_ads[2]*costheta*costheta*costheta*costheta*costheta*costheta);
			g+=color[1]*(k_ads[0]+k_ads[1]*cosalpha+k_ads[2]*costheta*costheta*costheta*costheta*costheta*costheta);
			b+=color[2]*(k_ads[0]+k_ads[1]*cosalpha+k_ads[2]*costheta*costheta*costheta*costheta*costheta*costheta);
		}

		else {
			r+=k_ads[0]*color[0];
			g+=k_ads[0]*color[1];
			b+=k_ads[0]*color[2];

		}

	}
	/*	r*=a_current; g*=a_current;  b*=a_current;
	if(accum_rgba[3]<0.99) {
	accum_rgba[0]+=r*(1.0-accum_rgba[3]);
	accum_rgba[1]+=g*(1.0-accum_rgba[3]);
	accum_rgba[2]+=b*(1.0-accum_rgba[3]);
	accum_rgba[3]+=a_current*(1.0-accum_rgba[3])*coff_i;
        }
	 */
	accum_rgba[0]=r;
	accum_rgba[1]=g;
	accum_rgba[2]=b;
	accum_rgba[3]=1.0;
	return;
}





