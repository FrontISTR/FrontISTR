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

#include "hecmw_vis_generate_histogram_vr.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_vis_bmp.h"
#include "hecmw_vis_font_texture.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_malloc.h"


void find_color_minmax_vr(double *var, int *empty_flag, int nx, int ny, int nz, double *mincolor, double *maxcolor)
{
	int i;
	double value;

	for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
		if(empty_flag[i]==1) {
			value=var[i];
			if(value<*mincolor) *mincolor=value;
			if(value>*maxcolor) *maxcolor=value;
		}
	}
	return;
}
void generate_histogram_graph_vr(double tmincolor, double tmaxcolor, double *var, int *empty_flag, int nx, int ny, int nz,int mynode,
		int pesize, HECMW_Comm VIS_COMM, int color_system_type)
{
	int i, j, k, m;
	double delta, value, color[3];
	int count[500], t_count[500], max_number, max_length, start_x, end_x, start_y, end_y;
	FILE *fp;
	double *graph;
	BITMAPFILEHEADER header;       /* File header */

	unsigned char r, g, b;
	int  ri, gi, bi;
	BITMAPINFOHEADER info;

	int start_xs, start_ys;
	char buf[128];
	int output7[7][7];

	delta=(tmaxcolor-tmincolor)/500.0;
	for(i=0;i<500;i++) {
		count[i]=0;
		t_count[i]=0;
	}
	for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
		if(empty_flag[i]==1) {
			j=(int)((var[i]-tmincolor)/delta);
			if(j<0) j=0;
			if(j>499) j=499;
			count[j]++;
		}
	}
	if(pesize>1)
		HECMW_Allreduce(count, t_count, 500, HECMW_INT, HECMW_SUM, VIS_COMM);
	else {
		for(i=0;i<500;i++)
			t_count[i]=count[i];
	}
	if(mynode==0) {
		fp=fopen("histogram.bmp", "wb");
		if(fp==NULL)
			HECMW_vis_print_exit("Cannot generate the histogram output file");
		graph=(double *)HECMW_calloc(400*530*3, sizeof(double));
		if(graph==NULL)
			HECMW_vis_memory_exit("graph");
		for(i=0;i<400*530*3;i++)
			graph[i]=0.0;
		max_number=0;
		for(i=0;i<500;i++) {
			if(t_count[i]>max_number)
				max_number=t_count[i];
		}
		if(max_number==0)
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2003:Cannot generate histogram graph, the number of voxels is 0");
		max_length=(int)(400-30-5-45*1.5);
		start_x=(int)(5+45*1.5+15);
		start_y=15;
		for(j=0;j<500;j++) {
			end_x=(int)((double)t_count[j]*(double)max_length/(double)max_number+start_x)+2;
			value=j/500.0;
			value2_to_rgb(value, color, color_system_type);
			for(i=start_x;i<end_x;i++) {
				graph[((j+15)*400+i)*3]=color[0];
				graph[((j+15)*400+i)*3+1]=color[1];
				graph[((j+15)*400+i)*3+2]=color[2];
			}
		}
		/*start mark scales */
		start_y=15;
		end_y=515;
		for(k=0;k<11; k++) {
			value=tmincolor+(tmaxcolor-tmincolor)/10.0*k;
			start_ys=start_y+(int)((double)500.0/10*k)-(int)7/2;
			start_xs=15;
			sprintf(buf, "%3.2E", value);
			for(m=0;m<9;m++) {
				font7_generate(buf[8-m], output7);
				for(j=0;j<7;j++)
					for(i=0;i<7;i++) {
						graph[((start_ys+j)*400+start_xs-i)*3]=(double)output7[6-j][i];
						graph[((start_ys+j)*400+start_xs-i)*3+1]=(double)output7[6-j][i];
						graph[((start_ys+j)*400+start_xs-i)*3+2]=(double)output7[6-j][i];
					}
				start_xs+=7;
				if((value>=0) && (m==0))
					start_xs-=7;
			}
			if((k!=0) && (k!=10)) {
				start_ys+=(int)7/2;

				for(i=start_x;i<start_x+5;i++)
					for(j=0;j<3;j++)
						graph[(start_ys*400+i)*3+j]=1.0;
			}

		}
		header.bfSize=54+3*400*530;
#ifdef CONVERSE_ORDER
		header.bfSize=change_unsigned_int_order(54+3*400*530);
#endif
		header.bfReserved1=0;
#ifdef CONVERSE_ORDER
		header.bfReserved1=change_short_int_order(0);
#endif

		header.bfReserved2=0;
#ifdef CONVERSE_ORDER
		header.bfReserved2=change_short_int_order(0);
#endif

		header.bfOffBits=54;
#ifdef CONVERSE_ORDER
		header.bfOffBits=change_unsigned_int_order(54);
#endif


		info.biBitCount=24;
#ifdef CONVERSE_ORDER
		info.biBitCount=change_short_int_order(24);
#endif

		info.biSize=40;
#ifdef CONVERSE_ORDER
		info.biSize=change_unsigned_int_order(40);
#endif

		info.biWidth=400;
#ifdef CONVERSE_ORDER
		info.biWidth=change_int_order(400);
#endif

		info.biHeight=530;
#ifdef CONVERSE_ORDER
		info.biHeight=change_int_order(530);
#endif

		info.biSizeImage=3*400*530;
#ifdef CONVERSE_ORDER
		info.biSizeImage=change_unsigned_int_order(3*400*530);
#endif

		info.biClrImportant=0;
#ifdef CONVERSE_ORDER
		info.biClrImportant=change_unsigned_int_order(0);
#endif

		info.biClrUsed=0;
#ifdef CONVERSE_ORDER
		info.biClrUsed=change_unsigned_int_order(0);
#endif

		info.biCompression=0;
#ifdef CONVERSE_ORDER
		info.biCompression=change_unsigned_int_order(0);
#endif

		info.biPlanes=1;
#ifdef CONVERSE_ORDER
		info.biPlanes=change_short_int_order(1);
#endif

		info.biXPelsPerMeter=3780;
#ifdef CONVERSE_ORDER
		info.biXPelsPerMeter=change_int_order(3780);
#endif

		info.biYPelsPerMeter=3780;
#ifdef CONVERSE_ORDER
		info.biYPelsPerMeter=change_int_order(3780);
#endif

		putc('B', fp);
		putc('M', fp);
		fwrite(&(header.bfSize), sizeof(unsigned int), 1,fp);
		fwrite(&header.bfReserved1, sizeof(unsigned short int), 1, fp);
		fwrite(&header.bfReserved2, sizeof(unsigned short int), 1, fp);
		fwrite(&header.bfOffBits, sizeof(unsigned int), 1, fp);
		fwrite(&info, 40, 1, fp);
		for(j=0;j<530;j++)
			for(i=400-1;i>=0;i--) {
				ri=(int)(graph[(j*400+i)*3]*256);
				gi=(int)(graph[(j*400+i)*3+1]*256);
				bi=(int)(graph[(j*400+i)*3+2]*256);
				if(ri<0) ri=0;
				if(ri>255) ri=255;
				if(gi<0) gi=0;
				if(gi>255) gi=255;
				if(bi<0) bi=0;
				if(bi>255) bi=255;
				r=ri;  g=gi; b=bi;
				putc(b, fp);
				putc(g, fp);
				putc(r, fp);
			}


		fclose(fp);
		HECMW_free(graph);
	}
	return;
}



void generate_interval_point_vr(double tmincolor, double tmaxcolor, double *var, int *empty_flag,
		int nx, int ny, int nz,int mynode, int pesize, HECMW_Comm VIS_COMM, double *interval_point)
{
	int i, j;
	double delta;
	int count[500], t_count[500], tmp_count[500], sum_count, interv_count, sum, current_j;
	delta=(tmaxcolor-tmincolor)/500.0;
	for(i=0;i<500;i++) {
		count[i]=0;
		t_count[i]=0;
	}
	for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
		if(empty_flag[i]==1) {
			j=(int)((var[i]-tmincolor)/delta);
			if(j<0) j=0;
			if(j>499) j=499;
			count[j]++;
		}
	}
	if(pesize>1)
		HECMW_Allreduce(count, t_count, 500, HECMW_INT, HECMW_SUM, VIS_COMM);
	else {
		for(i=0;i<500;i++)
			t_count[i]=count[i];
	}
	sum_count=0;
	for(i=0;i<500;i++) {
		sum_count+=t_count[i];
		tmp_count[i]=t_count[i];
	}
	interv_count=(int)sum_count/10;
	current_j=0;
	interval_point[0]=tmincolor;
	interval_point[1]=0.0;
	for(i=1;i<10;i++) {
		interval_point[i*2+1]=(double)i/10.0;
		sum=0;
		j=current_j;
		while((j<500) && (sum<interv_count)) {
			sum+=t_count[j];
			j++;
		}
		j--;
		interval_point[i*2]=((double)j/500.0+1.0/500.0*(1.0-(double)(sum-interv_count)/(double)tmp_count[j]))
		*(tmaxcolor-tmincolor)+tmincolor;
		t_count[j]=sum-interv_count;
		current_j=j;
	}
	interval_point[20]=tmaxcolor;
	interval_point[21]=1.0;
	fprintf(stderr, "The automatic color mapping set is :\n");
	for(i=0;i<11;i++)
		fprintf(stderr, "%lf    %lf   \n", interval_point[i*2], interval_point[i*2+1]);
	return;
}






void output_histogram_vr(double tmincolor, double tmaxcolor, double *var, int *empty_flag, int nx, int ny, int nz,int mynode, int pesize, HECMW_Comm VIS_COMM)
{
	int i, j;
	double delta;
	int count[100], t_count[100];
	FILE *fp;

	delta=(tmaxcolor-tmincolor)/100.0;
	for(i=0;i<100;i++) {
		count[i]=0;
		t_count[i]=0;
	}
	for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
		if(empty_flag[i]==1) {
			j=(int)((var[i]-tmincolor)/delta);
			if(j<0) j=0;
			if(j>99) j=99;
			count[j]++;
		}
	}
	if(pesize>1)
		HECMW_Allreduce(count, t_count, 100, HECMW_INT, HECMW_SUM, VIS_COMM);
	else {
		for(i=0;i<100;i++)
			t_count[i]=count[i];
	}
	if(mynode==0) {
		fp=fopen("histogram.file", "w");
		if(fp==NULL)
			HECMW_vis_print_exit("Cannot generate the histogram output file");
		for(i=0;i<100;i++)
			fprintf(fp, "%d   %d   -----(%lf --- %lf)\n", i, t_count[i], tmincolor+i*delta,  tmincolor+(i+1)*delta);
		fclose(fp);
	}
	return;
}

void find_minmax_vr(double *voxel_dxyz, double *voxel_orig_xyz, int mynode, double range[6])
{
	int  i;
	for(i=0;i<6;i++)
		range[i]=voxel_orig_xyz[mynode*3+i/2]+(i % 2)*voxel_dxyz[mynode*3+i/2];
	return;
}


void find_dis_minmax(double view_point_d[3], double vertex[24], double dis_minmax[2])
{
	int i;
	double dis;
	dis_minmax[0]=dis_minmax[1]=sqrt(SQR(vertex[0])+SQR(vertex[1])+SQR(vertex[2]));
	for(i=0;i<8;i++) {
		dis=sqrt(SQR(vertex[i*3]-view_point_d[0])+SQR(vertex[i*3+1]-view_point_d[1])
				+SQR(vertex[i*3+2]-view_point_d[2]));
		if(dis<dis_minmax[0]) dis_minmax[0]=dis;
		if(dis>dis_minmax[1]) dis_minmax[1]=dis;
	}
	return;
}

void find_feap_minmax(int num_of_features, double *fea_point, double mincolor, double maxcolor, double feap_minmax[2])
{
	int i,j;
	double t, mint, color;

	for(i=0;i<=255;i++) {
		color=mincolor+(maxcolor-mincolor)*i/255.0;
		mint=1.0E17;
		for(j=0;j<num_of_features;j++) {
			t=fabs(color-fea_point[j]);
			if(t<mint) mint=t;
		}
		if(mint>feap_minmax[1]) feap_minmax[1]=mint;
	}
	feap_minmax[0]=0.0;
	return;
}

void find_feai_minmax(int num_of_features, double *fea_point, double mincolor, double maxcolor, double feai_minmax[2])
{
	int i, j;
	double t,t1, t2, mint, color;

	for(i=0;i<=255;i++) {
		color=mincolor+(maxcolor-mincolor)*i/255.0;
		mint=1.0E17;
		for(j=0;j<num_of_features;j++) {
			if((t>=fea_point[j*2]) && (t<=fea_point[j*2+1]))
				t=0.0;
			else {
				t1=fabs(color-fea_point[j*2]);
				t2=fabs(color-fea_point[j*2+1]);
				if(t1<t2) t=t1;
				else t=t2;
			}
			if(t<mint) mint=t;
		}
		if(mint>feai_minmax[1]) feai_minmax[1]=mint;
	}
	feai_minmax[0]=0.0;
	return;
}
