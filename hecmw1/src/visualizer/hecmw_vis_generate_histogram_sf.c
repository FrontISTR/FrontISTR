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

#include "hecmw_vis_generate_histogram_sf.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_vis_bmp.h"
#include "hecmw_vis_font_texture.h"
#include "hecmw_vis_mem_util.h"
#include "hecmw_malloc.h"

void generate_histogram_graph_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue, Result *result, int mynode,
		int pesize, HECMW_Comm VIS_COMM, int color_system_type)
{
	int i, j, k, m, ii;
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
	int  color_id;


	for(i=0;i<data->nn_component;i++) {
		if(color_list[i]==1) {
			color_id=i;
			delta=(mavalue[i]-mivalue[i])/500.0;
		}
	}

	for(i=0;i<500;i++) {
		count[i]=0;
		t_count[i]=0;
	}

	for(ii=1;ii<=sf[0].surface_style;ii++) {
		if(sf[ii].display_method!=4) {
			for(i=0;i<result[ii-1].n_vertex;i++) {
				j=(int)((result[ii-1].color[i]-mivalue[color_id])/delta);
				if(j<0) j=0;
				if(j>499) j=499;
				count[j]++;
			}
		}
	}
	if(pesize>1)
		HECMW_Allreduce(count, t_count, 500, HECMW_INT, HECMW_SUM, VIS_COMM);
	else {
		for(i=0;i<500;i++)
			t_count[i]=count[i];
	}
	/*   fprintf(stderr, "count[2]=%d  t_count[2]=%d\n", count[2], t_count[2]);
	 */
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
			value=mivalue[color_id]+(mavalue[color_id]-mivalue[color_id])/10.0*k;
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



void generate_interval_point_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue,  Result *result,
		int mynode, int pesize, HECMW_Comm VIS_COMM, double *interval_point)
{
	int i, j, ii;
	double delta;
	int count[500], t_count[500], tmp_count[500], sum_count, interv_count, sum, current_j;
	int   color_id;

	for(i=0;i<data->nn_component;i++) {
		if(color_list[i]==1) {
			color_id=i;
			delta=(mavalue[i]-mivalue[i])/500.0;
		}
	}

	for(i=0;i<500;i++) {
		count[i]=0;
		t_count[i]=0;
	}

	for(ii=1;ii<=sf[0].surface_style;ii++) {
		if(sf[ii].display_method!=4) {
			for(i=0;i<result[ii-1].n_vertex;i++) {
				j=(int)((result[ii-1].color[i]-mivalue[color_id])/delta);
				if(j<0) j=0;
				if(j>499) j=499;
				count[j]++;
			}
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
	interval_point[0]=mivalue[color_id];
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
		*(mavalue[color_id]-mivalue[color_id])+mivalue[color_id];
		t_count[j]=sum-interv_count;
		current_j=j;
	}
	interval_point[20]=mavalue[color_id];
	interval_point[21]=1.0;
	fprintf(stderr, "The automatic color mapping set is :\n");
	for(i=0;i<11;i++)
		fprintf(stderr, "%lf    %lf   \n", interval_point[i*2], interval_point[i*2+1]);
	return;
}






void output_histogram_sf(struct surface_module *sf, int *color_list, struct hecmwST_result_data *data,
		double *mivalue, double *mavalue, Result *result,int mynode, int pesize, HECMW_Comm VIS_COMM)
{
	int i, j, ii;
	double delta;
	int count[100], t_count[100];
	FILE *fp;
	int   color_id;
	for(i=0;i<data->nn_component;i++) {
		if(color_list[i]==1) {
			color_id=i;
			delta=(mavalue[i]-mivalue[i])/100.0;
		}
	}

	for(i=0;i<100;i++) {
		count[i]=0;
		t_count[i]=0;
	}

	for(ii=1;ii<=sf[0].surface_style;ii++) {
		if(sf[ii].display_method!=4) {
			for(i=0;i<result[ii-1].n_vertex;i++) {
				j=(int)((result[ii-1].color[i]-mivalue[color_id])/delta);
				if(j<0) j=0;
				if(j>99) j=99;
				count[j]++;
			}
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
			fprintf(fp, "%d   %d   -----(%lf --- %lf)\n", i, t_count[i], mivalue[color_id]+i*delta, mivalue[color_id]+(i+1)*delta);
		fclose(fp);
	}
	return;
}


void find_minmax_sf(struct hecmwST_local_mesh *mesh, int mynode, double range[6])
{
	int  i;
	for(i=0;i<3;i++) {
		range[i*2]=1.0E17;
		range[i*2+1]=-1.0E17;
	}
	for(i=0;i<mesh->n_node;i++) {
		if(mesh->node[i*3]<range[0])
			range[0]=mesh->node[i*3];
		if(mesh->node[i*3]>range[1])
			range[1]=mesh->node[i*3];
		if(mesh->node[i*3+1]<range[2])
			range[2]=mesh->node[i*3+1];
		if(mesh->node[i*3+1]>range[3])
			range[3]=mesh->node[i*3+1];
		if(mesh->node[i*3+2]<range[4])
			range[4]=mesh->node[i*3+2];
		if(mesh->node[i*3+2]>range[5])
			range[5]=mesh->node[i*3+2];
	}


	return;
}

void find_patch_minmax_sf(Result *result, struct surface_module *sf, double range[6])
{
	int i, ii, j;

	for(i=0;i<3;i++) {
		range[i*2]=1.0E17;
		range[i*2+1]=-1.0E17;
	}

	for(ii=0;ii<sf[0].surface_style;ii++) {
		for(i=0;i<result[ii].n_vertex;i++) {
			for(j=0;j<3;j++) {
				if(result[ii].vertex[i*3+j]<range[j*2])
					range[j*2]=result[ii].vertex[i*3+j];
				if(result[ii].vertex[i*3+j]>range[j*2+1])
					range[j*2+1]=result[ii].vertex[i*3+j];
			}
		}


	}
	return;
}

void find_new_patch_minmax_sf(Result *result, struct surface_module *sf, double range[6])
{
	int i, ii, j;
	double new_v;


	for(ii=0;ii<sf[0].surface_style;ii++) {
		if(sf[ii+1].deform_display_on==1) {
			for(i=0;i<result[ii].n_vertex;i++) {
				for(j=0;j<3;j++) {
					new_v=result[ii].vertex[i*3+j]+result[ii].disp[i*3+j]*sf[ii+1].disp_scale;
					if(new_v<range[j*2])
						range[j*2]=new_v;
					if(new_v>range[j*2+1])
						range[j*2+1]=new_v;
				}
			}
		}


	}
	return;
}

void get_deform_scale(struct surface_module *sf, int ii, double range_x, double range_y, double range_z,
		struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int pesize, HECMW_Comm VIS_COMM)
{
	double max_disp, t_max_disp, tmp[3], disp, max_range, s_scale;
	int  i,j;
	int  tn_component, d_base;
	int mynode;

	tn_component=0;
	for(i=0;i<data->nn_component;i++)
		tn_component+=data->nn_dof[i];
	d_base=0;
	for(i=0;i<sf[ii].disp_comp;i++)
		d_base+=data->nn_dof[i];

	max_disp=0.0;

	for(i=0;i<mesh->n_node;i++) {
		for(j=0;j<3;j++)
			tmp[j]=data->node_val_item[i*tn_component+d_base+j];
		disp=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
		if(disp>max_disp)
			max_disp=disp;
	}
	HECMW_Comm_rank(VIS_COMM, &mynode);
	if(pesize>1)
		HECMW_Allreduce(&max_disp, &t_max_disp, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	else
		t_max_disp=max_disp;


	max_range=sqrt(range_x*range_x+range_y*range_y+range_z*range_z);
	if(fabs(t_max_disp)<EPSILON*0.0001)
		HECMW_vis_print_exit("No deformation in this dataset");
	s_scale=max_range*0.1/t_max_disp;
	if(sf[ii].disp_scale<0.0)
		sf[ii].disp_scale=s_scale;
	else
		sf[ii].disp_scale*=s_scale;
	if(mynode==0)
		fprintf(stderr, "The real deformation scale = %lf\n", sf[ii].disp_scale);
	return;
}
