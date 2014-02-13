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

#include "hecmw_vis_define_parameters.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"

#define EPSILON	 	0.00000001
#define PI  3.1415926


void transform_range_vertex(double range[6], double vertex[24])
{
	vertex[0*3]=vertex[4*3]=vertex[7*3]=vertex[3*3]=range[0];
	vertex[1*3]=vertex[5*3]=vertex[6*3]=vertex[2*3]=range[1];
	vertex[0*3+1]=vertex[1*3+1]=vertex[5*3+1]=vertex[4*3+1]=range[2];
	vertex[3*3+1]=vertex[2*3+1]=vertex[6*3+1]=vertex[7*3+1]=range[3];
	vertex[0*3+2]=vertex[1*3+2]=vertex[2*3+2]=vertex[3*3+2]=range[4];
	vertex[4*3+2]=vertex[5*3+2]=vertex[6*3+2]=vertex[7*3+2]=range[5];
	return;
}

static void normalize_f(double vector[3])
{
	int i;
	double norm_v;
	norm_v=sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
	if(fabs(norm_v)>EPSILON) {
		for(i=0;i<3;i++)
			vector[i]/=norm_v;
	}
	return;
}

void get_frame_transform_matrix(double view_point_d[3], double screen_point[3], double up[3], double coff_matrix[3][3])
{
	double U[3], V[3], N[3];
	int i;

	for(i=0;i<3;i++) {
		N[i]=-(view_point_d[i]-screen_point[i]);
	}
	normalize_f(N);
	/* find the direction of axis U */
	U[0]=up[1]*N[2]-N[1]*up[2];
	U[1]=-up[0]*N[2]+N[0]*up[2];
	U[2]=up[0]*N[1]-N[0]*up[1];
	normalize_f(U);
	/*find the direction of axix V */
	V[0]=N[1]*U[2]-U[1]*N[2];
	V[1]=-N[0]*U[2]+U[0]*N[2];
	V[2]=N[0]*U[1]-U[0]*N[1];
	normalize_f(V);
	for(i=0;i<3;i++) {
		coff_matrix[i][0]=U[i];
		coff_matrix[i][1]=V[i];
		coff_matrix[i][2]=N[i];
	}
	return;
}

void find_inverse_matrix(double coff_matrix[3][3], double inv_matrix[3][3])
{
	int i, j;
	double a[3][3], norm_a, aa[3][3];

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			a[i][j]=coff_matrix[i][j];
	norm_a=a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0]
	                                                                                                    -a[0][1]*a[1][0]*a[2][2]-a[0][0]*a[1][2]*a[2][1];
	if(fabs(norm_a)<1.0E-7)
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2001: There is something wrong with transform matrix, invers =0");

	aa[0][0]=  a[1][1]*a[2][2]-a[1][2]*a[2][1];
	aa[0][1]=-(a[0][1]*a[2][2]-a[2][1]*a[0][2]);
	aa[0][2]=  a[0][1]*a[1][2]-a[1][1]*a[0][2];
	aa[1][0]=-(a[1][0]*a[2][2]-a[2][0]*a[1][2]);
	aa[1][1]=  a[0][0]*a[2][2]-a[2][0]*a[0][2];
	aa[1][2]=-(a[0][0]*a[1][2]-a[1][0]*a[0][2]);
	aa[2][0]=  a[1][0]*a[2][1]-a[2][0]*a[1][1];
	aa[2][1]=-(a[0][0]*a[2][1]-a[2][0]*a[0][1]);
	aa[2][2]=  a[0][0]*a[1][1]-a[1][0]*a[0][1];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			inv_matrix[i][j]=aa[i][j]/norm_a;
	return;
}


void transform_frame(double screen_point[3], double vertex[24], double coff_matrix[3][3],
		double n_vertex[24])
{
	int i, j;
	double xx, yy, zz;

	for(i=0;i<24;i++)
		n_vertex[i]=vertex[i];

	for(i=0;i<8;i++)
		for(j=0;j<3;j++)
			n_vertex[i*3+j]-=screen_point[j];

	for(i=0;i<8;i++) {
		xx=n_vertex[i*3];
		yy=n_vertex[i*3+1];
		zz=n_vertex[i*3+2];
		n_vertex[i*3]=xx*coff_matrix[0][0]+yy*coff_matrix[1][0]+zz*coff_matrix[2][0];
		n_vertex[i*3+1]=xx*coff_matrix[0][1]+yy*coff_matrix[1][1]+zz*coff_matrix[2][1];
		n_vertex[i*3+2]=xx*coff_matrix[0][2]+yy*coff_matrix[1][2]+zz*coff_matrix[2][2];
	}
	return;
}

void transform_frame3(double screen_point[3], double f[3][3], double coff_matrix[3][3],
		double n_f[3][3])
{
	int i, j;
	double xx, yy, zz;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			n_f[i][j]=f[i][j];

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			n_f[i][j]-=screen_point[j];

	for(i=0;i<3;i++) {
		xx=n_f[i][0];
		yy=n_f[i][1];
		zz=n_f[i][2];
		n_f[i][0]=xx*coff_matrix[0][0]+yy*coff_matrix[1][0]+zz*coff_matrix[2][0];
		n_f[i][1]=xx*coff_matrix[0][1]+yy*coff_matrix[1][1]+zz*coff_matrix[2][1];
		n_f[i][2]=xx*coff_matrix[0][2]+yy*coff_matrix[1][2]+zz*coff_matrix[2][2];
	}
	return;
}


void transform2_frame(double coff_matrix[3][3], double view_point[3])
{
	double xx, yy, zz;

	xx=view_point[0];
	yy=view_point[1];
	zz=view_point[2];
	view_point[0]=xx*coff_matrix[0][0]+yy*coff_matrix[1][0]+zz*coff_matrix[2][0];
	view_point[1]=xx*coff_matrix[0][1]+yy*coff_matrix[1][1]+zz*coff_matrix[2][1];
	view_point[2]=xx*coff_matrix[0][2]+yy*coff_matrix[1][2]+zz*coff_matrix[2][2];


	return;
}

void tranverse_transform(double screen_point[3], double point_s[3], double inv_matrix[3][3], double point_o[3])
{
	int i;
	double xx, yy,zz;

	xx=point_s[0];  yy=point_s[1];  zz=point_s[2];
	point_o[0]=xx*inv_matrix[0][0]+yy*inv_matrix[1][0]+zz*inv_matrix[2][0];
	point_o[1]=xx*inv_matrix[0][1]+yy*inv_matrix[1][1]+zz*inv_matrix[2][1];
	point_o[2]=xx*inv_matrix[0][2]+yy*inv_matrix[1][2]+zz*inv_matrix[2][2];


	/* transform (x0, y0, z0)*/
	for(i=0;i<3;i++)
		point_o[i]+=screen_point[i];
	return;
}
/*
void output_frame(Parameter_rendering *vr, double view_point[3], double n_vertex[24], double scr_area[4],
				  double xd, double yd)
{

	int i, pp[8][2];
	int start_x, start_y;
    double p[8][2];
	FILE *f1;

    f1=fopen("frame.dat", "w");
	if(f1==NULL) {
		fprintf(stderr, "cannot open the frame output file\n");
		exit(0);
	}
	for(i=0;i<8;i++) {
		if(fabs(n_vertex[i*3+2]-view_point[2])<EPSILON) {
			fprintf(stderr, " The viewpoint position is not correct\n");
			exit(0);
		}
		p[i][0]=view_point[0]-view_point[2]/(n_vertex[i*3+2]-view_point[2])*
			(n_vertex[i*3]-view_point[0]);
		p[i][1]=view_point[1]-view_point[2]/(n_vertex[i*3+2]-view_point[2])*
			(n_vertex[i*3+1]-view_point[1]);
	}
	if((vr->color_mapping_bar_on==0) && (vr->scale_marking_on==0)) {
		start_x=10; start_y=10;
	}
	else if((vr->color_mapping_bar_on==1) && (vr->scale_marking_on==0)) {
		start_x=30;  start_y=10;
	}
	else if((vr->color_mapping_bar_on==1) && (vr->scale_marking_on==1)) {
		start_x=45; start_y=10;
	}
	for(i=0;i<8;i++) {
		pp[i][0]=(int)((p[i][0]-scr_area[0])/xd)+start_x;
		pp[i][1]=(int)((p[i][1]-scr_area[2])/yd)+start_y;
		fprintf(f1, "%d %d\n", pp[i][0], pp[i][1]);
	}
	fclose(f1);
	return;
}

 */
void transform_frame4(double screen_point[3], double iso_p[6], double coff_matrix[3][3], double n_iso[6])
{
	int i, j;
	double xx, yy, zz;

	for(i=0;i<6;i++)
		n_iso[i]=iso_p[i];

	for(i=0;i<2;i++)
		for(j=0;j<3;j++)
			n_iso[i*3+j]-=screen_point[j];

	for(i=0;i<2;i++) {
		xx=n_iso[i*3+0];
		yy=n_iso[i*3+1];
		zz=n_iso[i*3+2];
		n_iso[i*3+0]=xx*coff_matrix[0][0]+yy*coff_matrix[1][0]+zz*coff_matrix[2][0];
		n_iso[i*3+1]=xx*coff_matrix[0][1]+yy*coff_matrix[1][1]+zz*coff_matrix[2][1];
		n_iso[i*3+2]=xx*coff_matrix[0][2]+yy*coff_matrix[1][2]+zz*coff_matrix[2][2];
	}
	return;
}


void find_projection_range3(double view_point[3],double n_iso[6], double pixel_d[2][2], double iso_p[6])
{

	int i;
	double tmp_d;


	for(i=0;i<2;i++) {
		if(fabs(n_iso[i*3+2]-view_point[2])<EPSILON)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2002: The viewpoint position is not correct");
		pixel_d[i][0]=view_point[0]-view_point[2]/(n_iso[i*3+2]-view_point[2])*
		(n_iso[i*3+0]-view_point[0]);
		pixel_d[i][1]=view_point[1]-view_point[2]/(n_iso[i*3+2]-view_point[2])*
		(n_iso[i*3+1]-view_point[1]);
	}
	if(pixel_d[0][0]>pixel_d[1][0]) {
		for(i=0;i<2;i++) {
			tmp_d=pixel_d[1][i];
			pixel_d[1][i]=pixel_d[0][i];
			pixel_d[0][i]=tmp_d;
		}
		for(i=0;i<3;i++) {
			tmp_d=n_iso[3+i];
			n_iso[3+i]=n_iso[i];
			n_iso[i]=tmp_d;
		}
		for(i=0;i<3;i++) {
			tmp_d=iso_p[3+i];
			iso_p[3+i]=iso_p[i];
			iso_p[i]=tmp_d;
		}
	}


	return;
}


void find_projection_range2(double view_point[3],  double n_f[3][3],
		double scr_area[4])
{
	int i;
	double p[3][2];


	for(i=0;i<3;i++) {
		if(fabs(n_f[i][2]-view_point[2])<EPSILON)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2002: The viewpoint position is not correct");
		p[i][0]=view_point[0]-view_point[2]/(n_f[i][2]-view_point[2])*
		(n_f[i][0]-view_point[0]);
		p[i][1]=view_point[1]-view_point[2]/(n_f[i][2]-view_point[2])*
		(n_f[i][1]-view_point[1]);
	}
	scr_area[0]=scr_area[1]=p[0][0];
	scr_area[2]=scr_area[3]=p[0][1];
	for(i=0;i<3;i++) {
		if(p[i][0]<scr_area[0]) scr_area[0]=p[i][0];
		if(p[i][0]>scr_area[1]) scr_area[1]=p[i][0];
		if(p[i][1]<scr_area[2]) scr_area[2]=p[i][1];
		if(p[i][1]>scr_area[3]) scr_area[3]=p[i][1];
	}
	return;
}

void find_projection_range(double view_point[3],  double n_vertex[24],
		double scr_area[4])
{
	int i;
	double p[8][2];


	for(i=0;i<8;i++) {
		if(fabs(n_vertex[i*3+2]-view_point[2])<EPSILON)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2002: The viewpoint position is not correct");
		p[i][0]=view_point[0]-view_point[2]/(n_vertex[i*3+2]-view_point[2])*
		(n_vertex[i*3]-view_point[0]);
		p[i][1]=view_point[1]-view_point[2]/(n_vertex[i*3+2]-view_point[2])*
		(n_vertex[i*3+1]-view_point[1]);
	}
	scr_area[0]=scr_area[1]=p[0][0];
	scr_area[2]=scr_area[3]=p[0][1];
	for(i=0;i<8;i++) {
		if(p[i][0]<scr_area[0]) scr_area[0]=p[i][0];
		if(p[i][0]>scr_area[1]) scr_area[1]=p[i][0];
		if(p[i][1]<scr_area[2]) scr_area[2]=p[i][1];
		if(p[i][1]>scr_area[3]) scr_area[3]=p[i][1];
	}
	return;
}


void view_parameter_define(int ii, int num_of_frames, int rotate_style, double view_point_d[3], double screen_point[3], double up[3],
		int num_of_lights, double *light_point , double trange[6])
{
	int i,j;
	double center[3], t[3], angle, tminx, tmaxx, tminy, tmaxy, tminz, tmaxz;

	for(i=0;i<3;i++)
		center[i]=(trange[i*2]+trange[i*2+1])/2.0;
	angle=2.0*PI/(double)num_of_frames;
	if(rotate_style!=0) {
		if(rotate_style==1) { /*rotate 30 along x axis */
			/* rotate viewpoint */
			up[0]=1.0;
			up[1]=0.0;
			up[2]=0.0;
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0];
			view_point_d[1]=t[1]*cos(angle)+t[2]*sin(angle);
			view_point_d[2]=-t[1]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			/*			if((fabs(screen_point[1]-center[1])>EPSILON) || (fabs(screen_point[2]-center[2])>EPSILON)) {
			for(i=0;i<3;i++)
			  t[i]=screen_point[i]-center[i];
			screen_point[0]=t[0];
			screen_point[1]=t[1]*cos(angle)+t[2]*sin(angle);
			screen_point[2]=-t[1]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				screen_point[i]+=center[i];
			}
			 */
			/* rotate light_position */
			for(j=0;j<num_of_lights;j++) {
				for(i=0;i<3;i++)
					t[i]=light_point[j*3+i]-center[i];
				light_point[j*3+0]=t[0];
				light_point[j*3+1]=t[1]*cos(angle)+t[2]*sin(angle);
				light_point[j*3+2]=-t[1]*sin(angle)+t[2]*cos(angle);
				for(i=0;i<3;i++)
					light_point[j*3+i]+=center[i];
			}
			/* rotate up_direction */
			for(i=0;i<3;i++)
				t[i]=up[i];
			up[0]=t[0];
			up[1]=t[1]*cos(angle)+t[2]*sin(angle);
			up[2]=-t[1]*sin(angle)+t[2]*cos(angle);

		}
		if(rotate_style==2) { /*rotate 30 along y axis */
			/* rotate viewpoint */
			up[0]=0.0;
			up[1]=1.0;
			up[2]=0.0;
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0]*cos(angle)+t[2]*sin(angle);
			view_point_d[1]=t[1];
			view_point_d[2]=-t[0]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			/*			if((fabs(screen_point[0]-center[0])>EPSILON) || (fabs(screen_point[2]-center[2])>EPSILON)) {
			for(i=0;i<3;i++)
			  t[i]=screen_point[i]-center[i];
			screen_point[0]=t[0]*cos(angle)+t[2]*sin(angle);
			screen_point[1]=t[1];
			screen_point[2]=-t[0]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				screen_point[i]+=center[i];
			}
			 */
			/* rotate light_position */
			for(j=0;j<num_of_lights;j++) {
				for(i=0;i<3;i++)
					t[i]=light_point[j*3+i]-center[i];
				light_point[j*3+0]=t[0]*cos(angle)+t[2]*sin(angle);
				light_point[j*3+1]=t[1];
				light_point[j*3+2]=-t[0]*sin(angle)+t[2]*cos(angle);
				for(i=0;i<3;i++)
					light_point[j*3+i]+=center[i];
			}
			/* rotate up direction */
			for(i=0;i<3;i++)
				t[i]=up[i];
			up[0]=t[0]*cos(angle)+t[2]*sin(angle);
			up[1]=t[1];
			up[2]=-t[0]*sin(angle)+t[2]*cos(angle);
		}
		if(rotate_style==3) { /*rotate 30 along z axis */
			/* rotate viewpoint */
			/*			up[0]=0.0;
			up[1]=0.0;
			up[2]=1.0;
			 */
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0]*cos(angle)+t[1]*sin(angle);
			view_point_d[1]=-t[0]*sin(angle)+t[1]*cos(angle);
			view_point_d[2]=t[2];
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			/*			if((fabs(screen_point[0]-center[0])>EPSILON) || (fabs(screen_point[1]-center[1])>EPSILON)) {
			for(i=0;i<3;i++)
			  t[i]=screen_point[i]-center[i];
		screen_point[0]=t[0]*cos(angle)+t[1]*sin(angle);
			screen_point[1]=-t[0]*sin(angle)+t[1]*cos(angle);
			screen_point[2]=t[2];
			for(i=0;i<3;i++)
				screen_point[i]+=center[i];
			}
			 */
			/* rotate light_position */
			for(j=0;j<num_of_lights;j++) {
				for(i=0;i<3;i++)
					t[i]=light_point[j*3+i]-center[i];
				light_point[j*3+0]=t[0]*cos(angle)+t[1]*sin(angle);
				light_point[j*3+1]=-t[0]*sin(angle)+t[1]*cos(angle);
				light_point[j*3+2]=t[2];
				for(i=0;i<3;i++)
					light_point[j*3+i]+=center[i];
			}
			/* rotate up direction */
			for(i=0;i<3;i++)
				t[i]=up[i];
			up[0]=t[0]*cos(angle)+t[1]*sin(angle);
			up[1]=-t[0]*sin(angle)+t[1]*cos(angle);
			up[2]=t[2];

		}
		if(rotate_style==4) {
			if(ii>0) {
				tminx=trange[0]; tmaxx=trange[1];
				tminy=trange[2]; tmaxy=trange[3];
				tminz=trange[4]; tmaxz=trange[5];
				for(i=0;i<3;i++)
					screen_point[i]=center[i];
				num_of_lights=1;
				if(ii==1) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					light_point[0]=(tminx+tmaxx)/2.0;
					light_point[1]=tmaxy+0.1*(tmaxy-tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
				else if(ii==2) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=(tmaxz+tminz)/2.0;
					light_point[0]=(tminx+tmaxx)/2.0;
					light_point[1]=tmaxy+0.1*(tmaxy-tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*0.5;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
				else if(ii==3) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=tminz-1.5*(tmaxz-tminz);
					light_point[0]=(tminx+tmaxx)/2.0;
					light_point[1]=tmaxy+0.1*(tmaxy-tminy);
					light_point[2]=tminz-(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
				else if(ii==4) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					light_point[0]=(tminx+tmaxx)/2.0;
					light_point[1]=0.7*(tmaxy+tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=-1.0;
					up[2]=0.0;
				}
				else if(ii==5) {
					view_point_d[0]=tmaxx+1.5*(tmaxx-tminx);
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					light_point[0]=tmaxx+0.5*(tmaxx-tminx);
					light_point[1]=0.5*(tmaxy+tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
				else if(ii==6) {
					view_point_d[0]=tminx-1.5*(tmaxx-tminx);
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					light_point[0]=tminx-0.5*(tmaxx-tminx);
					light_point[1]=0.5*(tmaxy+tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
				else if(ii==7) {
					view_point_d[0]=(tmaxx+tminx)/2.0;
					view_point_d[1]=tminy-1.5*(tmaxy-tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					light_point[0]=(tmaxx+tminx)/2.0;
					light_point[1]=tminy-0.5*(tmaxy-tminy);
					light_point[2]=tmaxz+(tmaxz-tminz)*2.0;
					up[0]=0.0;
					up[1]=0.0;
					up[2]=1.0;
				}
			}

		}
	}
	return;
}


void view1_parameter_define(int ii, int num_of_frames, int rotate_style, double view_point_d[3], double screen_point[3], int num_of_lights, double *light_point, double up[3], double trange[6])

{
	int i;
	double center[3], t[3], angle, tminx, tmaxx, tminy, tmaxy, tminz, tmaxz;

	for(i=0;i<3;i++)
		center[i]=(trange[i*2]+trange[i*2+1])/2.0;
	angle=2.0*PI/(double)num_of_frames;
	if(rotate_style!=0) {
		if(rotate_style==1) { /*rotate 30 along x axis */
			/* rotate viewpoint */
			up[0]=1.0;
			up[1]=0.0;
			up[2]=0.0;
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0];
			view_point_d[1]=t[1]*cos(angle)+t[2]*sin(angle);
			view_point_d[2]=-t[1]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			/*			if((fabs(screen_point[1]-center[1])>EPSILON) || (fabs(screen_point[2]-center[2])>EPSILON)) {
			for(i=0;i<3;i++)
			  t[i]=screen_point[i]-center[i];
			screen_point[0]=t[0];
			screen_point[1]=t[1]*cos(angle)+t[2]*sin(angle);
			screen_point[2]=-t[1]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				screen_point[i]+=center[i];
			}
			 */
			/* rotate light_position */
			/*			for(j=0;j<num_of_lights;j++) {
			for(i=0;i<3;i++)
			  t[i]=light_point[j*3+i]-center[i];
			light_point[j*3+0]=t[0];
			light_point[j*3+1]=t[1]*cos(angle)+t[2]*sin(angle);
			light_point[j*3+2]=-t[1]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				light_point[j*3+i]+=center[i];
			}
			 */
		}
		if(rotate_style==2) { /*rotate 30 along y axis */
			/* rotate viewpoint */
			up[0]=0.0;
			up[1]=1.0;
			up[2]=0.0;
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0]*cos(angle)+t[2]*sin(angle);
			view_point_d[1]=t[1];
			view_point_d[2]=-t[0]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			if((fabs(screen_point[0]-center[0])>EPSILON) || (fabs(screen_point[2]-center[2])>EPSILON)) {
				for(i=0;i<3;i++)
					t[i]=screen_point[i]-center[i];
				screen_point[0]=t[0]*cos(angle)+t[2]*sin(angle);
				screen_point[1]=t[1];
				screen_point[2]=-t[0]*sin(angle)+t[2]*cos(angle);
				for(i=0;i<3;i++)
					screen_point[i]+=center[i];
			}
			/* rotate light_position */
			/*			for(j=0;j<num_of_lights;j++) {
			for(i=0;i<3;i++)
			  t[i]=light_point[j*3+i]-center[i];
			light_point[j*3+0]=t[0]*cos(angle)+t[2]*sin(angle);
			light_point[j*3+1]=t[1];
			light_point[j*3+2]=-t[0]*sin(angle)+t[2]*cos(angle);
			for(i=0;i<3;i++)
				light_point[j*3+i]+=center[i];
			}
			 */
		}
		if(rotate_style==3) { /*rotate 30 along z axis */
			/*			up[0]=0.0;
			up[1]=0.0;
			up[2]=1.0;
			 */
			/* rotate viewpoint */
			for(i=0;i<3;i++)
				t[i]=view_point_d[i]-center[i];
			view_point_d[0]=t[0]*cos(angle)+t[1]*sin(angle);
			view_point_d[1]=-t[0]*sin(angle)+t[1]*cos(angle);
			view_point_d[2]=t[2];
			for(i=0;i<3;i++)
				view_point_d[i]+=center[i];
			/* rotate screen_point if it is not in center */
			if((fabs(screen_point[0]-center[0])>EPSILON) || (fabs(screen_point[1]-center[1])>EPSILON)) {
				for(i=0;i<3;i++)
					t[i]=screen_point[i]-center[i];
				screen_point[0]=t[0]*cos(angle)+t[1]*sin(angle);
				screen_point[1]=-t[0]*sin(angle)+t[1]*cos(angle);
				screen_point[2]=t[2];
				for(i=0;i<3;i++)
					screen_point[i]+=center[i];
			}
			/* rotate light_position */
			/*			for(j=0;j<num_of_lights;j++) {
			for(i=0;i<3;i++)
			  t[i]=light_point[j*3+i]-center[i];
			light_point[j*3+0]=t[0]*cos(angle)+t[1]*sin(angle);
			light_point[j*3+1]=-t[0]*sin(angle)+t[1]*cos(angle);
			light_point[j*3+2]=t[2];
			for(i=0;i<3;i++)
				light_point[j*3+i]+=center[i];
			}
			 */
		}
		if(rotate_style==4) {
			if(ii>0) {
				tminx=trange[0]; tmaxx=trange[1];
				tminy=trange[2]; tmaxy=trange[3];
				tminz=trange[4]; tmaxz=trange[5];
				if(ii==1) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
				else if(ii==2) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=(tmaxz+tminz)/2.0;
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
				else if(ii==3) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
					view_point_d[2]=tminz-1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
				else if(ii==4) {
					view_point_d[0]=(tminx+tmaxx)/2.0;
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=-1.0;  up[2]=0.0;
				}
				else if(ii==5) {
					view_point_d[0]=tmaxx+1.5*(tmaxx-tminx);
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
				else if(ii==6) {
					view_point_d[0]=tminx-1.5*(tmaxx-tminx);
					view_point_d[1]=0.5*(tmaxy+tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
				else if(ii==7) {
					view_point_d[0]=(tmaxx+tminx)/2.0;
					view_point_d[1]=tminy-1.5*(tmaxy-tminy);
					view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
					up[0]=0.0; up[1]=0.0;  up[2]=1.0;
				}
			}

		}
	}
	return;
}
