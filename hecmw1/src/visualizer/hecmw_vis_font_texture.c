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

#include "hecmw_vis_font_texture.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_color_mapping.h"

#define BAR_WIDTH 10
#define EPSILON	 	0.00000001


static void font5_generate(char input, int output[5][5])
{
	int output0[5][5]={
			{0,1,1,1,0},
			{0,1,0,1,0},
			{0,1,0,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0}
	}, output1[5][5]={
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0}
	}, output2[5][5]={
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0}
	}, output3[5][5]={
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0}
	}, output4[5][5]={
			{0,1,0,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,0,0,1,0}
	}, output5[5][5]={
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0}
	}, output6[5][5]={
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0}
	}, output7[5][5]={
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,0,0,1,0},
			{0,0,0,1,0},
			{0,0,0,1,0}
	}, output8[5][5]={
			{0,1,1,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0}
	}, output9[5][5]={
			{0,1,1,1,0},
			{0,1,0,1,0},
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0}
	}, outpute[5][5]={
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0}
	}, outputd[5][5]={
			{0,0,0,0,0},
			{0,0,0,0,0},
			{0,0,0,0,0},
			{0,1,1,0,0},
			{0,1,1,0,0}
	}, outputa[5][5]={
			{0,0,0,0,0},
			{0,0,1,0,0},
			{0,1,1,1,0},
			{0,0,1,0,0},
			{0,0,0,0,0}
	}, outputs[5][5]={
			{0,0,0,0,0},
			{0,0,0,0,0},
			{0,1,1,1,0},
			{0,0,0,0,0},
			{0,0,0,0,0}
	}, outputt[5][5]={
			{1,1,1,1,1},
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0},
			{0,0,1,0,0}
	}, outputq[5][5]={
			{0,0,0,0,0},
			{0,1,1,1,0},
			{0,0,0,0,0},
			{0,1,1,1,0},
			{0,0,0,0,0}
	}, outputf[5][5]={
			{0,1,1,1,0},
			{0,1,0,0,0},
			{0,1,1,1,0},
			{0,0,0,1,0},
			{0,1,1,1,0}
	};
	int i,j;

	switch (input)
	{
	case '0':
	for (i = 0; i < 5; i++)
		for (j = 0; j < 5; j++)
			output[i][j] = output0[i][j];
	break;
	case '1':
	for (i = 0; i < 5; i++)
		for (j = 0; j < 5; j++)
			output[i][j] = output1[i][j];
	break;
	case '2':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output2[i][j];
		break;
	case '3':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output3[i][j];
		break;
	case '4':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output4[i][j];
		break;
	case '5':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output5[i][j];
		break;
	case '6':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output6[i][j];
		break;
	case '7':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output7[i][j];
		break;
	case '8':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output8[i][j];
		break;
	case '9':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = output9[i][j];
		break;
	case 'E':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outpute[i][j];
		break;
	case '.':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputd[i][j];
		break;
	case '+':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputa[i][j];
		break;
	case '-':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputs[i][j];
		break;
	case 'T':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputt[i][j];
		break;
	case '=':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputq[i][j];
		break;
	case 's':
		for (i = 0; i < 5; i++)
			for (j = 0; j < 5; j++)
				output[i][j] = outputf[i][j];
		break;
	default:
		for(i=0;i<5;i++)
			for(j=0;j<5;j++)
				output[i][j]=0;
		break;
	}
	return;
}


void font7_generate(char input, int output[7][7])
{
	int output0[7][7]={
			{0,0,1,1,1,0,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,0,1,1,1,0,0}
	}, output1[7][7]={
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0}
	}, output2[7][7]={
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{0,1,1,1,1,1,0}
	}, output3[7][7]={
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,1,1,1,1,1,0}
	}, output4[7][7]={
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0}
	}, output5[7][7]={
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,1,1,1,1,1,0}
	}, output6[7][7]={
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,1,1,1,1,0}
	}, output7[7][7]={
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0}
	}, output8[7][7]={
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,1,1,1,1,0}
	}, output9[7][7]={
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,1,0},
			{0,1,0,0,0,1,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,1,0},
			{0,0,0,0,0,1,0},
			{0,1,1,1,1,1,0}
	}, outpute[7][7]={
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,1,0,0,0,0,0},
			{0,1,0,0,0,0,0},
			{0,1,1,1,1,1,0}
	}, outputd[7][7]={
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,1,1,0,0,0},
			{0,0,1,1,0,0,0},
			{0,0,0,0,0,0,0}
	}, outputa[7][7]={
			{0,0,0,0,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,1,1,1,1,1,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,0,0,0,0}
	}, outputs[7][7]={
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0}
	}, outputt[7][7]={
			{1,1,1,1,1,1,1},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,1,0,0,0}
	}, outputq[7][7]={
			{0,0,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,0,0,0,0,0,0},
			{0,1,1,1,1,1,0},
			{0,0,0,0,0,0,0}
	}, outputf[7][7]={
			{0,0,1,1,1,0,0},
			{0,1,0,0,0,1,0},
			{0,0,1,0,0,0,0},
			{0,0,0,1,0,0,0},
			{0,0,0,0,1,0,0},
			{0,1,0,0,0,1,0},
			{0,0,1,1,1,0,0}
	};
	int i,j;

	switch (input)
	{
	case '0':
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
			output[i][j] = output0[i][j];
	break;
	case '1':
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
			output[i][j] = output1[i][j];
	break;
	case '2':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output2[i][j];
		break;
	case '3':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output3[i][j];
		break;
	case '4':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output4[i][j];
		break;
	case '5':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output5[i][j];
		break;
	case '6':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output6[i][j];
		break;
	case '7':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output7[i][j];
		break;
	case '8':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output8[i][j];
		break;
	case '9':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = output9[i][j];
		break;
	case 'E':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outpute[i][j];
		break;
	case '.':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputd[i][j];
		break;
	case '+':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputa[i][j];
		break;
	case '-':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputs[i][j];
		break;
	case 'T':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputt[i][j];
		break;
	case '=':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputq[i][j];
		break;
	case 's':
		for (i = 0; i < 7; i++)
			for (j = 0; j < 7; j++)
				output[i][j] = outputf[i][j];
		break;
	default:
		for(i=0;i<7;i++)
			for(j=0;j<7;j++)
				output[i][j]=0;
		break;
	}
	return;
}

void mark_time_label(double font_size,
		int xr, int yr,
		double font_color[3], double background_color[3],
		double start_time, double time_interval,
		int timestep, int max_len_step, double *image)
{
	int i,j,m;
	int start_ys, start_xs, scale, len_step;
	double vv;
	char buf[128], buf1[128], buf2[128], buf3[128];
	int output7[7][7], tmp1, tmp2;

	font_size=2.5;
	scale=(int)font_size;
	start_xs=30;
	start_ys=yr-25;
	vv=start_time+(double)timestep*time_interval;
	tmp1=(int)vv;
	tmp2=(int)((vv-(double)tmp1)*100.0);
	sprintf(buf1, "%d", tmp1);
	len_step=strlen(buf1);

	if(len_step<max_len_step)
	{
		for(i=0;i<max_len_step-len_step;i++)
			buf2[i]=' ';

		for(i=max_len_step-len_step; i<max_len_step;i++)
			buf2[i]=buf1[i-(max_len_step-len_step)];

		buf2[max_len_step]='\0';
	}
	else
	{
		sprintf(buf2, "%s", buf1);
	}

	sprintf(buf3, "%d", tmp2);

	if(strlen(buf3)<2)
	{
		buf3[2]='\0';
		buf3[1]=buf3[0];
		buf3[0]='0';
	}

	sprintf(buf, "%s%s.%s%s", "T=",buf2, buf3, "s");

	for(m=0;m<max_len_step+6;m++)
	{
		font7_generate(buf[max_len_step+5-m], output7);
		for(j=0;j<7*scale;j++)
			for(i=0;i<7*scale;i++)
			{
				if(output7[6-(int)(j/scale)][(int)(i/scale)]==1)
				{
					image[((start_ys+j)*xr+start_xs-i)*3]=font_color[0];
					image[((start_ys+j)*xr+start_xs-i)*3+1]=font_color[1];
					image[((start_ys+j)*xr+start_xs-i)*3+2]=font_color[2];
				}
				else
				{
					image[((start_ys+j)*xr+start_xs-i)*3]=background_color[0];
					image[((start_ys+j)*xr+start_xs-i)*3+1]=background_color[1];
					image[((start_ys+j)*xr+start_xs-i)*3+2]=background_color[2];
				}
			}
		start_xs+=7*scale;
	}
	return;
}

void value2_to_rgb(double value, double color[3], int color_system_type)
{
	double r, g, b;

	if(color_system_type==1)
	{
		if(value<0.0) value=0.0;
		if(value>1.0) value=1.0;

		if(value<=0.25)
		{
			r=0.0;
			g=value*4.0;
			b=1.0;
		}
		else if((value>0.25) && (value<=0.5))
		{
			r=0.0;
			g=1.0;
			b=(0.5-value)*4.0;
		}
		else if((value>0.5) && (value<=0.75))
		{
			r=(value-0.5)*4.0;
			g=1.0;
			b=0.0;
		}
		else if(value>0.75)
		{
			r=1.0;
			g=(1.0-value)*4.0;
			b=0.0;
		}

	}
	else if(color_system_type==2)
	{
		if(value<0.0) value=0.0;
		if(value>1.0) value=1.0;

		if(value<=0.2)
		{
			r=(0.2-value)*5.0;
			g=0.0;
			b=1.0;
		}
		else if((value>0.2) && (value<=0.4))
		{
			r=0.0;
			g=(value-0.2)*5.0;
			b=1.0;
		}
		else if((value>0.4) && (value<=0.6))
		{
			r=0.0;
			g=1.0;
			b=1.0-(value-0.4)*5.0;
		}
		else if((value>0.6) && (value<=0.8))
		{
			r=(value-0.6)*5.0;
			g=1.0;
			b=0.0;
		}
		else if(value>0.0)
		{
			r=1.0;
			g=1.0-(value-0.8)*5.0;
			b=0.0;
		}
	}
	else if(color_system_type==3)
	{
		r=g=b=value;
	}

	color[0]=r;
	color[1]=g;
	color[2]=b;

	return;
}


static double rgb_to_value(double value_rgb,
		double mincolor, double maxcolor,
		int color_mapping_style,
		double *interval_point,
		int interval_mapping_num)
{
	double value;
	int i;

	if(color_mapping_style==1)
	{
		if(fabs(maxcolor-mincolor)>EPSILON)
			value=value_rgb*(maxcolor-mincolor)+mincolor;
	}
	else if(color_mapping_style==2)
	{
		mincolor=interval_point[0];
		maxcolor=interval_point[1];

		if(fabs(maxcolor-mincolor)>EPSILON)
			value=value_rgb*(maxcolor-mincolor)+mincolor;
	}
	else if((color_mapping_style==3) || (color_mapping_style==4))
	{
		if(value_rgb<=interval_point[1])
		{
			value=mincolor;
		}
		else if(value_rgb>=interval_point[interval_mapping_num*2+1])
		{
			value=maxcolor;
		}
		else
		{
			for(i=1;i<interval_mapping_num+1;i++)
			{
				if((value_rgb<=interval_point[i*2+1]) &&
						(value_rgb>interval_point[(i-1)*2+1]))
				{
					value=(value_rgb         -interval_point[(i-1)*2+1])
					/(interval_point[i*2+1]-interval_point[(i-1)*2+1])
					*(interval_point[i*2  ]-interval_point[(i-1)*2  ])
					+interval_point[(i-1)*2];
					break;
				}
			}
		}
	}
	return value;
}

void generate_color_bar(int scale_marking_on, double font_size, int color_bar_style,
		int mark_0_on, int color_mapping_bar_on,
		int xr, int yr, double font_color[3], int color_system_type,
		int color_mapping_style, double *interval_point,
		int interval_mapping_num, int num_of_scale,
		double tmincolor, double tmaxcolor,
		double org_mincolor, double org_maxcolor, double *image)
{
	int i,j,k, m;
	int start_x, start_y, end_x, end_y, start_ys, start_xs, scale;
	double value, color[3], delta_y;
	double vv;
	char buf[128];
	int output5[5][5], output7[7][7], type;

	if(scale_marking_on==1)
	{
		if((font_size-(int)font_size)<0.5-EPSILON)
			type=1;
		else
			type=2;
	}

	if(color_bar_style==1)
	{
		if(scale_marking_on==0)
		{
			start_x=10;
		}
		else if(scale_marking_on==1)
		{
			if(type==1)
			{
				start_x=25;
				if(mark_0_on==1)
					start_x=40;

				if(font_size>=2.0)
				{
					start_x=25+(int)(font_size-1)*20;
					if(mark_0_on==1)
						start_x=40+(int)(font_size-1)*20;
				}
			}
			else if(type==2)
			{
				start_x=35;
				if(mark_0_on==1)
					start_x=55;

				if(font_size>=2.0)
				{
					start_x=35+(int)(font_size-1)*30;
					if(mark_0_on==1)
						start_x=55+(int)(font_size-1)*30;
				}
			}
		}
		start_y=(yr-20)/10+10;
		end_y=(yr-20)/10*5+start_y;
		end_x=start_x+BAR_WIDTH;
		delta_y=1.0/(double)(end_y-start_y);

		if(color_mapping_bar_on==1)
		{
			for(j=start_y;j<end_y;j++)
			{
				value=delta_y*(double)(j-start_y)*(tmaxcolor-tmincolor)+tmincolor;
				value_to_rgb(value, color, tmincolor, tmaxcolor, color_mapping_style,
						interval_point, interval_mapping_num, color_system_type);

				for(i=start_x;i<end_x;i++)
					for(k=0;k<3;k++)
						image[(j*xr+i)*3+k]=color[k];
			}
		}

		if(scale_marking_on==1)
		{
			/* transform tmincolor and tmaxcolor into scientific count */
			scale=(int)font_size;
			for(k=0;k<2; k++)
			{
				if(k==0)
				{
					vv=org_mincolor;
					if(type==1)
						start_ys=start_y-4-(int)font_size*5;
					else if(type==2)
						start_ys=start_y-4-(int)font_size*7;
				}
				else if(k==1)
				{
					vv=org_maxcolor; start_ys=end_y+4;
				}

				start_xs=10;
				sprintf(buf, "%10.2E", vv);

				for(m=0;m<10;m++)
				{
					if(type==1)
					{
						font5_generate(buf[9-m], output5);
						for(j=0;j<5*scale;j++)
							for(i=0;i<5*scale;i++)
							{
								image[((start_ys+j)*xr+start_xs-i)*3]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[0];
								image[((start_ys+j)*xr+start_xs-i)*3+1]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[1];
								image[((start_ys+j)*xr+start_xs-i)*3+2]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[2];
							}
						start_xs+=5*scale;
					}
					else if(type==2)
					{
						font7_generate(buf[9-m], output7);
						for(j=0;j<7*scale;j++)
							for(i=0;i<7*scale;i++)
							{
								image[((start_ys+j)*xr+start_xs-i)*3]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[0];
								image[((start_ys+j)*xr+start_xs-i)*3+1]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[1];
								image[((start_ys+j)*xr+start_xs-i)*3+2]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[2];
							}
						start_xs+=7*scale;
					}
				}
			}

			if(mark_0_on==1)
			{
				if((tmincolor<0) && (tmaxcolor>0))
				{
					start_ys=(int)((0-tmincolor)/(tmaxcolor-tmincolor)*(end_y-start_y)
							+start_y);
					start_xs=10;
					buf[0]='0';
					buf[1]='.';
					buf[2]='0';

					for(m=0;m<3;m++)
					{
						if(type==1)
						{
							font5_generate(buf[2-m], output5);
							for(j=0;j<5*scale;j++)
								for(i=0;i<5*scale;i++) {
									image[((start_ys+j)*xr+start_xs-i)*3]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[0];
									image[((start_ys+j)*xr+start_xs-i)*3+1]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[1];
									image[((start_ys+j)*xr+start_xs-i)*3+2]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[2];
								}
							start_xs+=5*scale;
						}
						else if(type==2)
						{
							font7_generate(buf[2-m], output7);
							for(j=0;j<7*scale;j++)
								for(i=0;i<7*scale;i++)
								{
									image[((start_ys+j)*xr+start_xs-i)*3]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[0];
									image[((start_ys+j)*xr+start_xs-i)*3+1]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[1];
									image[((start_ys+j)*xr+start_xs-i)*3+2]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[2];
								}
							start_xs+=7*scale;
						}
					}

					for(i=start_xs;i<start_x;i++)
						for(j=0;j<3;j++)
							image[(start_ys*xr+i)*3+j]=font_color[j];
				}
			}
		}
	}
	else if(color_bar_style==2)
	{
		if(scale_marking_on==0)
		{
			start_x=10;
		}
		else if(scale_marking_on==1)
		{
			scale=(int)font_size;
			if(type==1)
				start_x=45*scale+15;
			else if(type==2)
				start_x=63*scale+15;
		}

		start_y=(yr-20)/10+10;
		end_y=(yr-20)/10*5+start_y;
		end_x=start_x+BAR_WIDTH;
		delta_y=1.0/(double)(end_y-start_y);

		if(color_mapping_bar_on==1)
		{
			for(j=start_y;j<end_y;j++)
			{
				value=((double)(j-start_y))/(double)(end_y-start_y);
				value2_to_rgb(value, color, color_system_type);

				for(i=start_x;i<end_x;i++)
					for(k=0;k<3;k++)
						image[(j*xr+i)*3+k]=color[k];
			}
		}

		if(scale_marking_on==1)
		{
			/* transform tmincolor and tmaxcolor into scientific count */
			scale=(int)font_size;
			for(k=0;k<num_of_scale; k++)
			{
				if((k!=0) && (k!=num_of_scale-1))
					vv=rgb_to_value(1.0/(num_of_scale-1)*k, tmincolor, tmaxcolor,
							color_mapping_style, interval_point,
							interval_mapping_num);
				else if(k==0)
					vv=org_mincolor;
				else if(k==num_of_scale-1)
					vv=org_maxcolor;

				if(type==1)
					start_ys=start_y+(int)((double)(end_y-start_y)/(num_of_scale-1)*k)
					-(int)5*scale/2;
				else if(type==2)
					start_ys=start_y+(int)((double)(end_y-start_y)/(num_of_scale-1)*k)
					-(int)7*scale/2;

				start_xs=10;
				sprintf(buf, "%10.2E", vv);
				for(m=0;m<10;m++)
				{
					if(type==1)
					{
						font5_generate(buf[9-m], output5);
						for(j=0;j<5*scale;j++)
							for(i=0;i<5*scale;i++)
							{
								image[((start_ys+j)*xr+start_xs-i)*3]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[0];
								image[((start_ys+j)*xr+start_xs-i)*3+1]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[1];
								image[((start_ys+j)*xr+start_xs-i)*3+2]=
									(double)output5[4-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[2];
							}
						start_xs+=5*scale;
						if((vv>=0) && (m==9))
							start_xs-=5*scale;
					}
					else if(type==2)
					{
						font7_generate(buf[9-m], output7);
						for(j=0;j<7*scale;j++)
							for(i=0;i<7*scale;i++)
							{
								image[((start_ys+j)*xr+start_xs-i)*3]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[0];
								image[((start_ys+j)*xr+start_xs-i)*3+1]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[1];
								image[((start_ys+j)*xr+start_xs-i)*3+2]=
									(double)output7[6-(int)(j/scale)][(int)(i/scale)]
									                                  *font_color[2];
							}
						start_xs+=7*scale;
						if((vv>=0) && (m==9))
							start_xs-=7*scale;
					}
				}

				if((k!=0) && (k!=num_of_scale-1))
				{
					if(type==1)
						start_ys+=(int)5*scale/2;
					else if(type==2)
						start_ys+=(int)7*scale/2;

					for(i=start_x;i<end_x;i++)
						for(j=0;j<3;j++)
							image[(start_ys*xr+i)*3+j]=font_color[j];
				}
			}

			if(mark_0_on==1)
			{
				if((tmincolor<0) && (tmaxcolor>0))
				{
					start_ys=(int)((0-tmincolor)/(tmaxcolor-tmincolor)*(end_y-start_y)
							+start_y);
					start_xs=10;
					sprintf(buf, "%10.2E", 0.0);

					for(m=0;m<10;m++)
					{
						if(type==1)
						{
							font5_generate(buf[9-m], output5);
							for(j=0;j<5*scale;j++)
								for(i=0;i<5*scale;i++)
								{
									image[((start_ys+j)*xr+start_xs-i)*3]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[0];
									image[((start_ys+j)*xr+start_xs-i)*3+1]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[1];
									image[((start_ys+j)*xr+start_xs-i)*3+2]=
										(double)output5[4-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[2];
								}
							start_xs+=5*scale;
						}
						else if(type==2)
						{
							font7_generate(buf[9-m], output7);
							for(j=0;j<7*scale;j++)
								for(i=0;i<7*scale;i++)
								{
									image[((start_ys+j)*xr+start_xs-i)*3]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[0];
									image[((start_ys+j)*xr+start_xs-i)*3+1]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[1];
									image[((start_ys+j)*xr+start_xs-i)*3+2]=
										(double)output7[6-(int)(j/scale)][(int)(i/scale)]
										                                  *font_color[2];
								}
							start_xs+=7*scale;
						}
					}

					for(i=start_xs+1;i<start_xs+5;i++)
						for(j=0;j<3;j++)
							image[(start_ys*xr+i)*3+j]=font_color[j];
				}
			}
		}
	}

	return;
}


static unsigned short int transform_hex(char a)
{
	unsigned short int aa;

	switch (a)
	{
	case 'a':
		aa=10;
		break;
	case 'b':
		aa=11;
		break;
	case 'c':
		aa=12;
		break;
	case 'd':
		aa=13;
		break;
	case 'e':
		aa=14;
		break;
	case 'f':
		aa=15;
		break;
	case '0':
		aa=0;
		break;
	case '1':
		aa=1;
		break;
	case '2':
		aa=2;
		break;
	case '3':
		aa=3;
		break;
	case '4':
		aa=4;
		break;
	case '5':
		aa=5;
		break;
	case '6':
		aa=6;
		break;
	case '7':
		aa=7;
		break;
	case '8':
		aa=8;
		break;
	case '9':
		aa=9;
		break;
	}
	return (aa);
}



unsigned short int change_short_int_order(unsigned short int n)
{
	char c_buf[10], newc_buf[10], changed_buf[10];
	unsigned short int  i,m, digit;

	sprintf(c_buf, "%x", n);
	m=0;

	while(c_buf[m]!='\0')
		m++;

	if(m>4)
		HECMW_vis_print_exit("there is something wrong for the unsigned short int");

	for(i=0;i<4-m;i++)
		newc_buf[i]='0';

	for(i=0;i<m;i++)
		newc_buf[4-m+i]=c_buf[i];

	newc_buf[4]='\0';
	changed_buf[0]=newc_buf[2];
	changed_buf[1]=newc_buf[3];
	changed_buf[2]=newc_buf[0];
	changed_buf[3]=newc_buf[1];
	changed_buf[4]='\0';
	digit=0;

	for(i=0;i<4;i++)
	{
		digit*=16;
		digit+=transform_hex(changed_buf[i]);
	}
	return (digit);
}

unsigned int change_unsigned_int_order(unsigned int n)
{
	char c_buf[10], newc_buf[10], changed_buf[10], changed2_buf[10];
	unsigned int i, m, digit;

	sprintf(c_buf, "%x", n);
	m=0;

	for(i=0;i<8;i++)
	{
		if(c_buf[i]=='\0')
			break;
		else
			m++;
	}

	if(m>8)
		HECMW_vis_print_exit("there is something wrong for the unsigned short int");

	for(i=0;i<8-m;i++)
		newc_buf[i]='0';

	for(i=0;i<m;i++)
		newc_buf[8-m+i]=c_buf[i];

	newc_buf[8]='\0';

	for(i=0;i<4;i++)
		changed_buf[i]=newc_buf[i+4];

	for(i=0;i<4;i++)
		changed_buf[i+4]=newc_buf[i];

	changed2_buf[0]=changed_buf[2];
	changed2_buf[1]=changed_buf[3];
	changed2_buf[2]=changed_buf[0];
	changed2_buf[3]=changed_buf[1];
	changed2_buf[4]=changed_buf[6];
	changed2_buf[5]=changed_buf[7];
	changed2_buf[6]=changed_buf[4];
	changed2_buf[7]=changed_buf[5];
	changed2_buf[8]='\0';
	digit=0;

	for(i=0;i<8;i++)
	{
		digit*=16;
		digit+=transform_hex(changed2_buf[i]);
	}

	return (digit);
}



int change_int_order(int n)
{
	char c_buf[10], newc_buf[10], changed_buf[10], changed2_buf[10];
	int i, m, digit;

	sprintf(c_buf, "%x", n);
	m=0;

	for(i=0;i<8;i++)
	{
		if(c_buf[i]=='\0')
			break;
		else
			m++;
	}

	if(m>8)
		HECMW_vis_print_exit("there is something wrong for the unsigned short int");

	for(i=0;i<8-m;i++)
		newc_buf[i]='0';

	for(i=0;i<m;i++)
		newc_buf[8-m+i]=c_buf[i];

	newc_buf[8]='\0';

	for(i=0;i<4;i++)
		changed_buf[i]=newc_buf[i+4];

	for(i=0;i<4;i++)
		changed_buf[i+4]=newc_buf[i];

	changed2_buf[0]=changed_buf[2];
	changed2_buf[1]=changed_buf[3];
	changed2_buf[2]=changed_buf[0];
	changed2_buf[3]=changed_buf[1];
	changed2_buf[4]=changed_buf[6];
	changed2_buf[5]=changed_buf[7];
	changed2_buf[6]=changed_buf[4];
	changed2_buf[7]=changed_buf[5];
	changed2_buf[8]='\0';
	digit=0;

	for(i=0;i<8;i++)
	{
		digit*=16;
		digit+=transform_hex(changed2_buf[i]);
	}
	return (digit);
}
