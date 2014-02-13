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

#include "hecmw_vis_color_composite_vr.h"

#include <math.h>
#include "hecmw_vis_ray_trace.h"
#include "hecmw_malloc.h"


#if 0
void find_sample_point(double in_point[3], double out_point[3], int num_sample, double *sample_point)
{
	int i,j;
	double t;

	t=1.0/(double)(num_sample+1);
	for(i=0;i<num_sample;i++) {
		for(j=0;j<3;j++)
			sample_point[i*3+j]=in_point[j]+t*(i+1)*(out_point[j]-in_point[j]);
	}
	return;
}

double value_compute(double point[3], int current_ijk[3], double orig_xyz[3], double r_dxyz[3], int r_level[3], double *var)
{
	int i, j,k, index, m;
	double value[8], vv[8*3], v;
	double dis[8], d;
	int  return_flag;
	return_flag=0;
	vv[0*3]=vv[4*3]=vv[7*3]=vv[3*3]=orig_xyz[0]+current_ijk[0]*r_dxyz[0];
	vv[1*3]=vv[5*3]=vv[6*3]=vv[2*3]=orig_xyz[0]+(current_ijk[0]+1)*r_dxyz[0];
	vv[0*3+2]=vv[1*3+2]=vv[2*3+2]=vv[3*3+2]=orig_xyz[2]+current_ijk[2]*r_dxyz[2];
	vv[4*3+2]=vv[5*3+2]=vv[6*3+2]=vv[7*3+2]=orig_xyz[2]+(current_ijk[2]+1)*r_dxyz[2];
	vv[0*3+1]=vv[1*3+1]=vv[5*3+1]=vv[4*3+1]=orig_xyz[1]+current_ijk[1]*r_dxyz[1];
	vv[3*3+1]=vv[2*3+1]=vv[6*3+1]=vv[7*3+1]=orig_xyz[1]+(current_ijk[1]+1)*r_dxyz[1];


	for(m=0;m<8;m++) {
		i=m % 2;
		j=(m/2) % 2;
		k=(m/2) /2;
		index=(current_ijk[2]+k)*(r_level[0]+1)*(r_level[1]+1)+(current_ijk[1]+j)*(r_level[0]+1)+
		current_ijk[0]+i;
		value[k*4+j*2+i]=var[index];
	}

	for(i=0;i<8;i++) {
		dis[i]=sqrt(SQR(point[0]-vv[i*3])+SQR(point[1]-vv[i*3+1])+SQR(point[2]-vv[i*3+2]));
		if(fabs(dis[i])<EPSILON) {
			v=value[i];
			return_flag=1;
		}
	}
	if(return_flag==0) {
		d=0.0;
		for(i=0;i<8;i++)
			d+=1.0/dis[i];
		v=0.0;
		for(i=0;i<8;i++) {
			v+=value[i]/(dis[i]*d);

		}
	}
	return v;
}


/*
void find_grad_minmax(int r_nxyz[3], grad_var, double grad_minmax[2])
{
	int i;
	double grad;
	grad_minmax[0]=grad_minmax[1]=sqrt(SQR(grad_var[0])+SQR(grad_var[1])
	+SQR(grad_var[2]));
	for(i=0;i<(r_nxyz[0]+1)*(r_nxyz[1]+1)*(r_nxyz[2]+1);i++) {
		grad=sqrt(SQR(grad_var[i*3])+SQR(grad_var[i*3+1])
			+SQR(grad_var[i*3+2]));
		if(grad<grad_minmax[0]) grad_minmax[0]=grad;
		if(grad>grad_minmax[1]) grad_minmax[1]=grad;
	}
	return;
}
 */





double opacity_decision(int current_ijk[3],  int transfer_function_style, double opa_value, int num_of_features, double *fea_point,
		double view_point_d[3], int color_mapping_style, double *interval_point, int interval_mapping_num,
		double orig_xyz[3], double r_dxyz[3],
		double value, double n[3], double grad_minmax[2],
		double feap_minmax[2], double feai_minmax[2],
		double dis_minmax[2], double *opa_table, double mincolor, double maxcolor, int time_step)
{
	double opacity, grad;
	double t, mint;
	int i, /* cell_id, */ level, min_type;
	double cp[3], del_l, dist;
	/*	cell_id=in_voxel->cell_id;
	 */
	if(transfer_function_style==1)
		opacity=opa_value;
	else if(transfer_function_style==2) {
		grad=sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
		opacity=(grad-grad_minmax[0])/(grad_minmax[1]-grad_minmax[0])/200.0+0.0002;
	}
	else if(transfer_function_style==3) {
		/*		value=vd->var[cell_id*vd->varnumtot+vr->color_comp];
		 */
		mint=1.0E17;
		for(i=0;i<num_of_features;i++) {
			t=fabs(value- fea_point[i*3]);
			if(t<mint) {
				mint=t;
				min_type=i;
			}
		}
		/*		if(mint<0.00001) {
                t=feap_minmax[1]-mint+feap_minmax[0];
		opacity=(t-feap_minmax[0])/(feap_minmax[1]-feap_minmax[0])*3.05+0.15;
        }
		 */
		if(mint<fea_point[min_type*3+1]) {
			opacity=fea_point[min_type*3+2]*(fea_point[min_type*3+1]-mint)/fea_point[min_type*3+1]+opa_value;
		}
		else opacity=opa_value;
	}
	else if(transfer_function_style==4) {
		/*		value=vd->var[cell_id*vd->varnumtot+vr->color_comp];
		 */
		/*		mint=1.0E17;
		for(i=0;i<vr->num_of_features;i++) {
			if((value>=vr->fea_point[i*3]) && (value<=vr->fea_point[i*3+1]))
				t=0.0;
			else {
				t1=fabs(value-vr->fea_point[i*2]);
				t2=fabs(value-vr->fea_point[i*2+1]);
				if(t1<t2) t=t1;
				else t=t2;
			}
			if(t<mint) mint=t;
		}
		t=feai_minmax[1]-mint+feai_minmax[0];
		opacity=(t-feai_minmax[0])/(feai_minmax[1]-feai_minmax[0])*0.4+0.1;
		 */
		opacity=opa_value;
		for(i=0;i<num_of_features;i++) {
			if((value>=fea_point[i*3]) && (value<=fea_point[i*3+1])) {
				opacity=fea_point[i*3+2];
				break;
			}
		}

	}
	else if(transfer_function_style==5) {
		/*		for(i=0;i<3;i++)
			cp[i]=(in_voxel->bound_box[i*2+1]+in_voxel->bound_box[i*2])/2.0;
		 */
		for(i=0;i<3;i++)
			cp[i]=(current_ijk[i]+0.5)*r_dxyz[i]+orig_xyz[i];
		dist=sqrt(SQR(cp[0]-view_point_d[0])+SQR(cp[1]-view_point_d[1])
				+SQR(cp[2]-view_point_d[2]));
		dist=dis_minmax[1]-dist+dis_minmax[0];
		opacity=(dist-dis_minmax[0])/(dis_minmax[1]-dis_minmax[0])/200.0+0.0002;
	}
	else if(transfer_function_style==6) {
		/*		for(i=0;i<3;i++)
			cp[i]=(in_voxel->bound_box[i*2+1]+in_voxel->bound_box[i*2])/2.0;
		 */
		for(i=0;i<3;i++)

			cp[i]=(current_ijk[i]+0.5)*r_dxyz[i]+orig_xyz[i];

		dist=sqrt(SQR(cp[0]-view_point_d[0])+SQR(cp[1]-view_point_d[1])
				+SQR(cp[2]-view_point_d[2]));
		opacity=(dist-dis_minmax[0])/(dis_minmax[1]-dis_minmax[0])/200.0+0.0002;
	}
	else if(transfer_function_style==7) {


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
			for(i=1;i<interval_mapping_num+1;i++) {
				if((value<=interval_point[i*2]) && (value>interval_point[(i-1)*2])) {
					value=(value-interval_point[(i-1)*2])/(interval_point[i*2]-interval_point[(i-1)*2])*
					(interval_point[i*2+1]-interval_point[(i-1)*2+1])+interval_point[(i-1)*2+1];
					break;
				}
			}
		}

		if(value>1.0) value=1.0;
		if(value<0.0) value=0.0;
		time_step+=0;
		/*      ss=(sqrt((double)time_step)-6.0)/60.0;
   	opacity=((value-0.6)/6.0*((double)time_step-40.0)/7.0)+0.23-ss;

        ss=((double)time_step)/((time_step+280));
		 */
		opacity=value/200.0+0.0002;;
		if(opacity<0) opacity=0.0;
	}
	else if(transfer_function_style==8) {
		del_l=(maxcolor-mincolor)/255.0;
		level=(int)(value-mincolor)/del_l;
		if(level<0) level=0;
		if(level>255) level=255;
		opacity=opa_table[level];
	}

	return(opacity);
}
#endif

#ifdef surface_vr
int find_surface_point(int index,VR_data *vd,double in[3], double out[3], double *surf_p)
{
	int i,j,k, m, n1, n2, flag, flag_inside;
	double fp[3][3], n_f[4], v_f[3], d, point[3], s[3], sina, cosa, ss, t, color, dis[3], dist, n_norm, v_norm, nv_sign;
	int   num_surp, num_surp1;
	double *dd, dd_norm;
	int    *o_flag;
	double *surf_p1;
	double  tmp_dd, tmp_p[7];

	num_surp=0;
	surf_p1=(double *)HECMW_calloc(7*vd->surface[index].num, sizeof(double));
	if(surf_p1==NULL) {
		fprintf(stderr, "There is no enough memory for surf_p\n"),
		exit(0);
	}

	for(m=0;m<vd->surface[index].num;m++) {
		/*find the equation of the patch first */

		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				fp[i][j]=vd->surface[index].surf_data[m*12+i*3+j];
		n_f[0]=(fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=(fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
		for(i=0;i<3;i++)
			v_f[i]=in[i]-out[i];
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
		if(fabs(n_f[0]*(out[0]-in[0])+n_f[1]*(out[1]-in[1])+n_f[2]*(out[2]-in[2]))>EPSILON) {
			t=(-n_f[3]-n_f[0]*in[0]-n_f[1]*in[1]-n_f[2]*in[2])/(n_f[0]*(out[0]-in[0])+n_f[1]*(out[1]-in[1])+n_f[2]*(out[2]-in[2]));
			if((t>-EPSILON) && (t<1.0)) {
				for(j=0;j<3;j++)
					point[j]=in[j]+t*(out[j]-in[j]);
				/*judge whether it is inside the range of patch by equal area method*/
				flag_inside=0;
				for(j=0;j<3;j++) {
					n1=j;
					n2=j+1;
					if(n2==3) n2=0;
					if((fabs(point[0]-fp[n1][0])<EPSILON) && (fabs(point[1]-fp[n1][1])<EPSILON) && (fabs(point[2]-fp[n1][2])<EPSILON)) {
						flag_inside=1;
						break;
					}
					else
						cosa=((point[0]-fp[n1][0])*(fp[n2][0]-fp[n1][0])+(point[1]-fp[n1][1])*(fp[n2][1]-fp[n1][1])+
								(point[2]-fp[n1][2])*(fp[n2][2]-fp[n1][2]))/(sqrt(SQR(point[0]-fp[n1][0])+SQR(point[1]-fp[n1][1])
										+SQR(point[2]-fp[n1][2]))*sqrt(SQR(fp[n2][0]-fp[n1][0])+SQR(fp[n2][1]-fp[n1][1])+SQR(fp[n2][2]-fp[n1][2])));
					sina=sqrt(1-cosa*cosa);
					s[j]=sqrt(SQR(fp[n2][0]-fp[n1][0])+SQR(fp[n2][1]-fp[n1][1])+SQR(fp[n2][2]-fp[n1][2]))*sqrt(SQR(point[0]-fp[n1][0])+SQR(point[1]-fp[n1][1])
							+SQR(point[2]-fp[n1][2]))*sina/2.0;
				}
				cosa=((fp[2][0]-fp[0][0])*(fp[1][0]-fp[0][0])+(fp[2][1]-fp[0][1])*(fp[1][1]-fp[0][1])+
						(fp[2][2]-fp[0][2])*(fp[1][2]-fp[0][2]))/(sqrt(SQR(fp[2][0]-fp[0][0])+SQR(fp[2][1]-fp[0][1])
								+SQR(fp[2][2]-fp[0][2]))*sqrt(SQR(fp[1][0]-fp[0][0])+SQR(fp[1][1]-fp[0][1])+SQR(fp[1][2]-fp[0][2])));
				sina=sqrt(1-cosa*cosa);
				ss=sqrt(SQR(fp[1][0]-fp[0][0])+SQR(fp[1][1]-fp[0][1])+SQR(fp[1][2]-fp[0][2]))*sqrt(SQR(fp[2][0]-fp[0][0])+SQR(fp[2][1]-fp[0][1])
						+SQR(fp[2][2]-fp[0][2]))*sina/2.0;
				if(flag_inside==0) {
					/*		 if(fabs(ss-s[0]-s[1]-s[2])<EPSILON*ss/100.0)
					 */
					if(fabs(ss-s[0]-s[1]-s[2])<EPSILON*ss/100.0)
						flag_inside=1;
				}
				/*this intersection point on the patch*/
				if(flag_inside==1) {
					flag=1;
					for(j=0;j<3;j++) {
						dis[j]=sqrt(SQR(point[0]-fp[j][0])+SQR(point[1]-fp[j][1])+SQR(point[2]-fp[j][2]));
						if(dis[j]<EPSILON) {
							flag=0;
							color=vd->surface[index].surf_data[m*12+9+j];
							break;
						}
					}
					if(flag==1) {
						dist=1.0/(1.0/dis[0]+1.0/dis[1]+1.0/dis[2]);
						color=vd->surface[index].surf_data[m*12+9]*dist/dis[0]+vd->surface[index].surf_data[m*12+10]*dist/dis[1]
						                                                                                                      +vd->surface[index].surf_data[m*12+11]*dist/dis[2];
					}



					for(j=0;j<3;j++) {
						surf_p1[num_surp*7+j]=point[j];
						surf_p1[num_surp*7+3+j]=n_f[j];
						surf_p1[num_surp*7+6]=color;

					}
					num_surp++;
				}
			}
		} /* end of if fabs >epsilon */
	}  /*end of m*/
	/*delete overlapped point*/
	num_surp1=0;
	if(num_surp>0) {
		dd=(double *)HECMW_calloc(num_surp, sizeof(double));
		o_flag=(int *)HECMW_calloc(num_surp, sizeof(int));
		if((dd==NULL) || (o_flag==NULL)) {
			fprintf(stderr, "There is no enough memory for dd and overlap_flag\n");
			exit(0);
		}
		for(i=0;i<num_surp;i++)
			dd[i]=sqrt(SQR(surf_p1[i*7+0]-in[0])+SQR(surf_p1[i*7+1]-in[1])+SQR(surf_p1[i*7+2]-in[2]));
		dd_norm=sqrt(SQR(out[0]-in[0])+SQR(out[1]-in[1])+SQR(out[2]-in[2]));
		for(i=0;i<num_surp-1;i++)
			for(j=i+1;j<num_surp;j++) {
				if(dd[i]>dd[j]) {
					tmp_dd=dd[i];
					dd[i]=dd[j];
					dd[j]=tmp_dd;
					for(k=0;k<7;k++)
						tmp_p[k]=surf_p1[i*7+k];
					for(k=0;k<7;k++)
						surf_p1[i*7+k]=surf_p1[j*7+k];
					for(k=0;k<7;k++)
						surf_p1[j*7+k]=tmp_p[k];
				}
			}
		for(i=0;i<num_surp;i++)
			o_flag[i]=1;
		for(i=1;i<num_surp;i++) {
			if(fabs(dd[i]-dd[i-1])<EPSILON*dd_norm*100.0)
				o_flag[i]=0;
		}
		/*        for(i=1;i<num_surp;i++)
           o_flag[i]=0;
		 */
		num_surp1=0;
		for(i=0;i<num_surp;i++) {
			if(o_flag[i]==1) {
				for(j=0;j<7;j++)
					surf_p[num_surp1*7+j]=surf_p1[i*7+j];
				num_surp1++;
			}
		}
		HECMW_free(dd);
		HECMW_free(o_flag);
		HECMW_free(surf_p1);
	}
	return(num_surp1);
}

#endif




void compute_color_vr(int current_ijk[3], int color_mapping_style, double *interval_point,int transfer_function_style,
		double opa_value, int num_of_features, double *fea_point,
		double view_point_d[3],  int interval_mapping_num, int color_system_type, int num_of_lights,
		double *light_point, double k_ads[3], int r_level[3], double orig_xyz[3],
		double r_dxyz[3], double *var, double *grad_var, double accum_rgba[4],
		double mincolor, double maxcolor, double grad_minmax[2], double feap_minmax[2],
		double feai_minmax[2], double dis_minmax[2], double *opa_table,
		double in_point[3], double out_point[3], double tav_length,
		int time_step, int print_flag)
{
	int i, j, k;
	double a_current;
	double color[3];
	double value;
	double length, coff_i;
	int   index;
	int   m;
	double value2[8], vv[8*3], v;
	double dis[8], d;
	int  return_flag;
	double t;
	double r, g, b;



	length=sqrt(SQR(out_point[0]-in_point[0])+SQR(out_point[1]-in_point[1])+SQR(out_point[2]-
			in_point[2]));
	coff_i=length/tav_length;
	index=current_ijk[2]*r_level[0]*r_level[1]+current_ijk[1]*r_level[0]+current_ijk[0];





	/*---------------------start computing value of in_point --------------------*/
	return_flag=0;
	vv[0*3]=vv[4*3]=vv[7*3]=vv[3*3]=orig_xyz[0]+current_ijk[0]*r_dxyz[0];
	vv[1*3]=vv[5*3]=vv[6*3]=vv[2*3]=orig_xyz[0]+(current_ijk[0]+1)*r_dxyz[0];
	vv[0*3+2]=vv[1*3+2]=vv[2*3+2]=vv[3*3+2]=orig_xyz[2]+current_ijk[2]*r_dxyz[2];
	vv[4*3+2]=vv[5*3+2]=vv[6*3+2]=vv[7*3+2]=orig_xyz[2]+(current_ijk[2]+1)*r_dxyz[2];
	vv[0*3+1]=vv[1*3+1]=vv[5*3+1]=vv[4*3+1]=orig_xyz[1]+current_ijk[1]*r_dxyz[1];
	vv[3*3+1]=vv[2*3+1]=vv[6*3+1]=vv[7*3+1]=orig_xyz[1]+(current_ijk[1]+1)*r_dxyz[1];


	for(m=0;m<8;m++) {
		i=m % 2;
		j=(m/2) % 2;
		k=(m/2) /2;
		index=(current_ijk[2]+k)*(r_level[0]+1)*(r_level[1]+1)+(current_ijk[1]+j)*(r_level[0]+1)+
		current_ijk[0]+i;
		value2[k*4+j*2+i]=var[index];
	}

	for(i=0;i<8;i++){
		dis[i]=sqrt(SQR(in_point[0]-vv[i*3])+SQR(in_point[1]-vv[i*3+1])+SQR(in_point[2]-vv[i*3+2]));
	}
	/* if(return_flag==0) {
	 */
	 d=0.0;
	 for(i=0;i<8;i++)
		 d+=1.0/(dis[i]+EPSILON);
	 v=0.0;
	 for(i=0;i<8;i++) {
		 v+=value2[i]/((dis[i]+EPSILON)*d);

	 }
	 /* }*/
	 /*---------------------end computing value of in_point --------------------*/

	 value=v;


#ifdef transfer_change

	 /*	if(transfer_function_style==1)
		opacity=opa_value;
	if(transfer_function_style==2) {
		grad=sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
             opacity=(grad-grad_minmax[0])/(grad_minmax[1]-grad_minmax[0])/200.0+0.0002;
        }

	  */
#endif
	 a_current=opa_value;
	 if(transfer_function_style==3) {
		 a_current=opa_value;

		 for(i=0;i<num_of_features;i++) {
			 t=fabs(value- fea_point[i*3]);
			 if(t<fea_point[i*3+1])
				 a_current+=fea_point[i*3+2]*(fea_point[i*3+1]-t)/fea_point[i*3+1];
		 }

	 }

	 if(transfer_function_style==4) {
		 a_current=opa_value;
		 for(i=0;i<num_of_features;i++) {
			 if((value>=fea_point[i*3]) && (value<=fea_point[i*3+1])) {
				 a_current=fea_point[i*3+2];
			 }
		 }

	 }

	 /*--------------map value to rgb -------------------   */
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

	 if(color_mapping_style==3) {

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
	 /*	for(j=0;j<num_of_lights;j++) {
	  */
#ifdef slow
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
#endif

	 /*		r=color[0]*(ka+kd*cosalpha+ks*pow(costheta,6));
		g=color[1]*(ka+kd*cosalpha+ks*pow(costheta,6));
		b=color[2]*(ka+kd*cosalpha+ks*pow(costheta,6));
	  */

	 r+=color[0]*k_ads[0]*coff_i;
	 g+=color[1]*k_ads[0]*coff_i;
	 b+=color[2]*k_ads[0]*coff_i;

	 /*
	}
	  */
	 r*=a_current; g*=a_current;  b*=a_current;
	 if(accum_rgba[3]<0.99) {
		 accum_rgba[0]+=r*(1.0-accum_rgba[3]);
		 accum_rgba[1]+=g*(1.0-accum_rgba[3]);
		 accum_rgba[2]+=b*(1.0-accum_rgba[3]);
		 accum_rgba[3]+=a_current*(1.0-accum_rgba[3])*coff_i;
	 }


	 return;
}





