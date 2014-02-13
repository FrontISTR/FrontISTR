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

#include "hecmw_vis_read_control.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_malloc.h"


int is_blank_line(char *buf){
	int i;
	int flag;

	i=0;
	flag=1;
	while(buf[i]!='\n') {
		if(buf[i]!=' ') {
			flag=0;
			break;
		}
		i++;
	}
	return(flag);
}

int is_comment_line(char *buf) {
	int  flag;
	flag=0;
	if(buf[0]=='#')
		flag=1;
	else if((buf[0]=='!') && (buf[1]=='!'))
		flag=1;
	return(flag);
}


void get_string_item(char *para, char *buf, int *start_location, char para2[128])
{
	int  i,j;

	i=*start_location;
	while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
		i++;
	if(buf[i]=='\n') {
		fprintf(stderr, "No string value for %s\n", para);
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	}
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para2[j]=buf[i];
		i++;
		j++;
	}
	para2[j]='\0';
	*start_location=i;
	return;
}

int get_int_item(char *para, char *buf, int *start_location)
{
	int value;
	int  i,j;
	char para2[128];

	i=*start_location;
	while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
		i++;
	if(buf[i]=='\n') {
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0003:The control parameter format error:No integer value for %s");
	}
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para2[j]=buf[i];
		i++;
		j++;
	}
	para2[j]='\0';
	if((isdigit(para2[0])==0) && (para2[0]!='+') && (para2[0]!='-')) {
		fprintf(stderr, "ERROR: HEC-MW-VIS-E0004: %s should be integer \n", para);
		HECMW_vis_print_exit("Please re-input and run again!");
	}
	value=atoi(para2);
	*start_location=i;
	return(value);
}


double get_double_item(char *para, char *buf, int *start_location)
{
	double value;
	int  i,j;
	char para2[128];

	i=*start_location;
	while((buf[i]==',') || (buf[i]==' ') || (buf[i]=='='))
		i++;
	if(buf[i]=='\n') {
		fprintf(stderr, "No integer value for %s\n", para);
		HECMW_vis_print_exit("The control parameter format error:!");
	}
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para2[j]=buf[i];
		i++;
		j++;
	}
	para2[j]='\0';
	if((isdigit(para2[0])==0) && (para2[0]!='+') && (para2[0]!='-')) {
		fprintf(stderr, "ERROR: HEC-MW-VIS-E0005:%s should be a real \n", para);
		HECMW_vis_print_exit("Please re-input and run again!");
	}
	value=atof(para2);
	*start_location=i;
	return(value);
}

int get_keyword_item(char *buf, char *para) {
	int i,j;
	i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!') {
		fprintf(stderr, "Please check the line %s\n", buf);
		HECMW_vis_print_exit("The control parameter format error:!");
	}
	i++;
	j=0;
	while((buf[i]!=' ') && (buf[i]!='=') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	return (i);
}

static int identify_surface(char *buf) {
	int i,j, ii, len_tmp;
	int flag;
	char para[128], para1[128];
	i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	i=i+1;  j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	flag=1;
	len_tmp=strlen(para);
	sprintf(para1, "%s", para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para1[ii]);
	}
	if(strlen(para)>8)
		flag=0;
	else if((strncmp(para, "SURFACE", 7)!=0) && (strncmp(para, "surface", 7)!=0))
		flag=0;
	return (flag);
}

#if 0
int identify_rendering(char *buf) {
	int i,j, ii, len_tmp;
	int flag;
	char para[128];
	i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	i=i+1;  j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	flag=1;
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}

	if(strlen(para)>8)
		flag=0;
	else if((strncmp(para, "SURFACE", 7)!=0) && (strncmp(para, "surface", 7)!=0))
		flag=0;
	return (flag);
}
#endif

static int get_keyword_visual(char *buf) {
	int i,j, ii, len_tmp;
	int flag;
	char para[128];

	i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	i=i+1;  j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	flag=0;
	para[j]='\0';
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}

	if((strncmp(para, "VISUAL", 6)!=0) && (strncmp(para, "visual", 6)!=0)) {
		flag=0;
		return (flag);
	}
	while((buf[i]!='=') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0006:The control parameter format error: visual ID");
	i++;
	while((buf[i]==' ') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0006:The control parameter format error: visual ID");
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}

	if((strncmp(para, "PSR", 3)==0) || (strncmp(para, "psr", 3)==0)) {
		flag=1;
		return (flag);
	}

	if((strncmp(para, "PVR", 3)==0) || (strncmp(para, "pvr", 3)==0)) {
		flag=2;
		return (flag);
	}

	HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0007:The control parameter format error: method only can be PSR or PVR");
	return (flag);
}


static int set_keyword_visual(char *buf, PSF_link *psf, PVR_link *pvr) {
	int i,j,ii,len_tmp;
	int flag;
	char para[128];
	PSF_link  *p1, *p2;
	PVR_link  *t1, *t2;
	i=0;
	while(buf[i]==' ')
		i++;
	if(buf[i]!='!')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	i=i+1;  j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	flag=0;
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}

	if((strncmp(para, "VISUAL", 6)!=0) && (strncmp(para, "visual", 6)!=0)) {
		flag=0;
		return (flag);
	}
	while((buf[i]!='=') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0006:The control parameter format error: visual ID");
	i++;
	while((buf[i]==' ') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0006:The control parameter format error: visual ID");
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}

	if((strncmp(para, "PSR", 3)==0) || (strncmp(para, "psr", 3)==0)) {
		flag=1;
		p1=(PSF_link *)HECMW_calloc(1,sizeof(PSF_link));
		if(p1==NULL)
			HECMW_vis_memory_exit("PSF_link: p1");
		p2=psf->next_psf;
		psf->next_psf=p1;
		p1->next_psf=p2;
		psf->num_of_psf++;
		p1->visual_start_step=-1;
		p1->visual_end_step=-1;
		p1->visual_interval_step=1;


		while((buf[i]==',') || (buf[i]==' '))
			i++;
		while(buf[i]!='\n') {
			j=0;
			while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
				para[j]=buf[i];
				i++;
				j++;
			}
			para[j]='\0';
			len_tmp=strlen(para);
			for(ii=0;ii<len_tmp;ii++) {
				para[ii]=tolower(para[ii]);
			}

			if(strncmp(para, "visual_start_step", 12)==0)
				p1->visual_start_step=get_int_item(para, buf, &i);
			else if(strncmp(para, "visual_end_step", 10)==0)
				p1->visual_end_step=get_int_item(para, buf, &i);
			else if(strncmp(para, "visual_interval_step", 14)==0)
				p1->visual_interval_step=get_int_item(para, buf, &i);
			while((buf[i]==',') || (buf[i]==' '))
				i++;
		}

		if((p1->visual_end_step>=1) && (p1->visual_start_step==-1))
			p1->visual_start_step=1;
		return (flag);
	}
	if((strncmp(para, "PVR", 3)==0) || (strncmp(para, "pvr", 3)==0)) {
		flag=2;
		t1=(PVR_link *)HECMW_malloc(sizeof(PVR_link));
		if(t1==NULL)
			HECMW_vis_memory_exit("PVR_link: t1");
		t2=pvr->next_pvr;
		pvr->next_pvr=t1;
		t1->next_pvr=t2;
		pvr->num_of_pvr++;
		t1->visual_start_step=-1;
		t1->visual_end_step=-1;
		t1->visual_interval_step=1;


		while((buf[i]==',') || (buf[i]==' '))
			i++;
		while(buf[i]!='\n') {
			j=0;
			while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='=') && (buf[i]!='\n')) {
				para[j]=buf[i];
				i++;
				j++;
			}
			para[j]='\0';
			len_tmp=strlen(para);
			for(ii=0;ii<len_tmp;ii++) {
				para[ii]=tolower(para[ii]);
			}

			if(strncmp(para, "visual_start_step", 12)==0)
				t1->visual_start_step=get_int_item(para, buf, &i);
			else if(strncmp(para, "visual_end_step", 10)==0)
				t1->visual_end_step=get_int_item(para, buf, &i);
			else if(strncmp(para, "visual_interval_step", 14)==0)
				t1->visual_interval_step=get_int_item(para, buf, &i);
			while((buf[i]==',') || (buf[i]==' '))
				i++;
		}
		if((t1->visual_end_step>=1) && (t1->visual_start_step==-1))
			t1->visual_start_step=1;
		return (flag);
	}

	HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0007:The control parameter format error: method only can be PSR or PVR");
	return (flag);
}


static int get_keyword_surface_num(char *buf) {
	int i,j,ii, len_tmp;
	char para[128];
	int surface_num;

	if(buf[0]!='!')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0002: The control parameter format error: should start from !");
	i=1;  j=0;
	while((buf[i]!=' ') && (buf[i]!='=') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	len_tmp=strlen(para);
	for(ii=0;ii<len_tmp;ii++) {
		para[ii]=tolower(para[ii]);
	}
	if((strncmp(para, "surface_num", 11)!=0) && (strncmp(para, "SURFACE_NUM", 6)!=0)) {
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0008:The control parameter format error: surface_num should be defined");
	}
	while((buf[i]!='=') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0008:The control parameter format error: surface_num should be defined");
	i++;
	while((buf[i]==' ') && (buf[i]!='\n'))
		i++;
	if(buf[i]=='\n')
		HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0008:The control parameter format error: surface_num should be defined");
	j=0;
	while((buf[i]!=' ') && (buf[i]!=',') && (buf[i]!='\n')) {
		para[j]=buf[i];
		i++;
		j++;
	}
	para[j]='\0';
	surface_num=atoi(para);

	return (surface_num);
}


void HECMW_vis_read_control(FILE *fp, int pesize, int mynode, PSF_link *psf, PVR_link *pvr)
{
	int i, j, k,ii, len_tmp;
	char		buf[MAX_LINE_LEN];
	char**      parameters;
	int      *len_para;

	char     para[128], para1[128];
	int      hit;
	int      surface_num;
	int     *stat_para;
	int     location, visual_method;
	PSF_link  *p1;
	struct surface_module	*sf;
	Parameter_rendering      *sr;
	PVR_link  *t1;
	Parameter_vr      *vr;
	int    cont_flag;

	cont_flag=1;
	while (cont_flag) {
		if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
			if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
		}
		else
			cont_flag=0;
	}
	/*  fseek(fp, offset, SEEK_SET);
	 */

	hit=0;
	visual_method=0;
	if(cont_flag==1)
		visual_method=get_keyword_visual(buf);
	while((visual_method>0) && (cont_flag==1)){
		set_keyword_visual(buf, psf, pvr);
		while (cont_flag) {
			if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
				if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
			}
			else
				cont_flag=0;
		}
		if((visual_method==1) && (cont_flag==1)) {
			visual_method=0;
			surface_num=get_keyword_surface_num(buf);
			sf=(struct surface_module *)HECMW_calloc(surface_num+1, sizeof(struct surface_module));
			if(sf==NULL)
				HECMW_vis_memory_exit("surface parameters: sf");
			sf[0].surface_style=surface_num;
			sr=(Parameter_rendering *)HECMW_calloc(1, sizeof(Parameter_rendering));
			if(sr==NULL)
				HECMW_vis_memory_exit("rendering parameters: sr");
			sr->projection_style=1;
			parameters=(char **)HECMW_calloc(NUM_CONTROL_PSF, sizeof(char *));
			len_para=(int *)HECMW_calloc(NUM_CONTROL_PSF, sizeof(int));
			stat_para=(int *)HECMW_calloc(NUM_CONTROL_PSF, sizeof(int));
			for(i=0;i<NUM_CONTROL_PSF;i++) {
				parameters[i]=(char *)HECMW_calloc(128, sizeof(char));
				if(parameters[i]==NULL)
					HECMW_vis_memory_exit("tempory parameters: Parameter");
			}
			strcpy(parameters[0], "surface_style");
			len_para[0]=12;
			strcpy(parameters[1], "group_name");
			len_para[1]=10;
			strcpy(parameters[2], "defined_style");
			len_para[2]=12;
			strcpy(parameters[3], "data_comp_name");
			len_para[3]=12;
			strcpy(parameters[4], "data_comp");
			len_para[4]=9;
			strcpy(parameters[5], "data_subcomp_name");
			len_para[5]=14;
			strcpy(parameters[6], "data_subcomp");
			len_para[6]=12;
			strcpy(parameters[7], "iso_value");
			len_para[7]=9;
			strcpy(parameters[8], "method");
			len_para[8]=6;
			strcpy(parameters[9], "point");
			len_para[9]=5;
			strcpy(parameters[10], "radius");
			len_para[10]=6;
			strcpy(parameters[11], "length");
			len_para[11]=6;
			strcpy(parameters[12], "coef");
			len_para[12]=4;
			strcpy(parameters[13], "display_method");
			len_para[13]=14;
			strcpy(parameters[14], "color_comp_name");
			len_para[14]=14;
			strcpy(parameters[15], "color_comp");
			len_para[15]=10;
			strcpy(parameters[16], "color_subcomp_name");
			len_para[16]=18;
			strcpy(parameters[17], "color_subcomp");
			len_para[17]=13;
			strcpy(parameters[18], "isoline_number");
			len_para[18]=14;
			strcpy(parameters[19], "specified_color");
			len_para[19]=14;
			strcpy(parameters[20], "output_type");
			len_para[20]=11;
			strcpy(parameters[21], "range_filename");
			len_para[21]=10;
			strcpy(parameters[22], "normalize_on");
			len_para[22]=12;

			strcpy(parameters[23], "x_resolution");
			strcpy(parameters[24], "y_resolution");
			len_para[23]=len_para[24]=12;
			strcpy(parameters[25], "num_of_lights");
			len_para[25]=13;
			strcpy(parameters[26], "position_of_lights");
			len_para[26]=18;
			strcpy(parameters[27], "viewpoint");
			len_para[27]=9;
			strcpy(parameters[28], "look_at_point");
			len_para[28]=13;
			strcpy(parameters[29], "up_direction");
			len_para[29]=12;
			strcpy(parameters[30], "ambient_coef");
			len_para[30]=12;
			strcpy(parameters[31], "diffuse_coef");
			len_para[31]=12;
			strcpy(parameters[32], "specular_coef");
			len_para[32]=13;
			strcpy(parameters[33], "color_mapping_style");
			len_para[33]=19;
			strcpy(parameters[34], "interval_mapping_num");
			len_para[34]=20;
			strcpy(parameters[35], "interval_mapping");
			len_para[35]=16;
			strcpy(parameters[36], "rotate_style");
			len_para[36]=12;
			strcpy(parameters[37], "rotate_num_of_frames");
			len_para[37]=15;
			strcpy(parameters[38], "color_mapping_bar_on");
			len_para[38]=17;
			strcpy(parameters[39], "scale_marking_on");
			len_para[39]=15;
			strcpy(parameters[40], "background_color");
			len_para[40]=12;
			strcpy(parameters[41], "font_color");
			len_para[41]=9;
			strcpy(parameters[42], "color_system_type");
			len_para[42]=12;
			strcpy(parameters[43], "font_size");
			len_para[43]=8;
			strcpy(parameters[44], "color_bar_style");
			len_para[44]=13;
			strcpy(parameters[45], "fixed_range_on");
			len_para[45]=11;
			strcpy(parameters[46], "range_value");
			len_para[46]=11;
			strcpy(parameters[47], "num_of_scale");
			len_para[47]=11;
			strcpy(parameters[48], "mark_0_on");
			len_para[48]=9;
			strcpy(parameters[49], "opacity_mapping_style");
			len_para[49]=21;
			strcpy(parameters[50], "opa_value");
			len_para[50]=8;
			strcpy(parameters[51], "num_of_features");
			len_para[51]=14;
			strcpy(parameters[52], "fea_point");
			len_para[52]=8;
			strcpy(parameters[53], "lookup_filename");
			len_para[53]=15;
			strcpy(parameters[54], "histogram_on");
			len_para[54]=11;
			strcpy(parameters[55], "boundary_line_on");
			len_para[55]=11;
			strcpy(parameters[56], "isoline_color");
			len_para[56]=12;
			strcpy(parameters[57], "time_mark_on");
			len_para[57]=10;
			strcpy(parameters[58], "fixed_scale_mark");
			len_para[58]=15;
			strcpy(parameters[59], "deform_display_on");
			len_para[59]=13;
			strcpy(parameters[60], "deform_scale");
			len_para[60]=11;
			strcpy(parameters[61], "initial_style");
			len_para[61]=12;
			strcpy(parameters[62], "deform_style");
			len_para[62]=11;
			strcpy(parameters[63], "deform_comp_name");
			len_para[63]=15;
			strcpy(parameters[64], "deform_comp");
			len_para[64]=11;
			strcpy(parameters[65], "initial_line_color");
			len_para[65]=15;
			strcpy(parameters[66], "deform_line_color");
			len_para[66]=15;
			strcpy(parameters[67], "deform_num_of_frames");
			len_para[67]=15;
			strcpy(parameters[68], "smooth_shading_on");
			len_para[68]=15;
			strcpy(parameters[69], "real_deform_scale");
			len_para[69]=16;
			strcpy(parameters[70], "fixed_mesh_range");
			len_para[70]=14;
			strcpy(parameters[71], "start_time");
			len_para[71]=9;
			strcpy(parameters[72], "time_interval");
			len_para[72]=9;

			for(i=0;i<NUM_CONTROL_PSF;i++) {
				stat_para[i]=0;
			}


			while (cont_flag==1) {
				if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
					if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
				}
				else
					cont_flag=0;
			}
			for(k=1;k<surface_num+1;k++) {
				sprintf(sf[k].data_comp_name, "%s", "NULL");
				sprintf(sf[k].color_comp_name, "%s", "NULL");
				sprintf(sf[k].data_subcomp_name, "%s", "NULL");
				sprintf(sf[k].color_subcomp_name, "%s", "NULL");
				sprintf(sf[k].disp_comp_name, "%s", "NULL");
				sf[k].data_comp=-1;
				sf[k].data_subcomp=-1;
				sf[k].color_comp=-1;
				sf[k].color_subcomp=-1;
				sf[k].isoline_number=0;
				sf[k].range_output=0;
				sf[k].normalize_flag=1;
				sf[k].disp_comp=-1;
				for(i=0;i<23;i++) {
					stat_para[i]=0;
				}
				for(i=59;i<67;i++)
					stat_para[i]=0;
				if(cont_flag==1)
					if((get_keyword_visual(buf)==0) && (identify_surface(buf)==1)) {
						while (cont_flag==1) {
							if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
								if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
							}
							else
								cont_flag=0;
						}
						if(cont_flag==1)
							while((get_keyword_visual(buf)==0) && (identify_surface(buf)==0)) {
								hit=-1;
								location=get_keyword_item(buf, para);
								len_tmp=strlen(para);
								for(ii=0;ii<len_tmp;ii++) {
									para[ii]=tolower(para[ii]);
								}
								for(i=0;i<NUM_CONTROL_PSF;i++) {
									if((strncmp(para, parameters[i], len_para[i]))==0) {
										hit=i;
										stat_para[i]=1;
										break;
									}
								}
								/*	 fprintf(stderr, "para=%s hit=%d\n", para, hit);
								 */
								if((hit>=0) && (hit<NUM_CONTROL_PSF)) {
									switch (hit) {
									case 0:
										sf[k].surface_style=get_int_item(para, buf, &location);
										break;
									case 1:
										get_string_item(para, buf, &location, sf[k].group_name);
										break;
									case 2:
										sf[k].defined_style=get_int_item(para, buf, &location);
										break;
									case 3:
										get_string_item(para, buf, &location,sf[k].data_comp_name);
										break;
									case 4:
										sf[k].data_comp=get_int_item(para, buf, &location);
										break;
									case 5:
										get_string_item(para, buf, &location, sf[k].data_subcomp_name);
										break;

									case 6:
										sf[k].data_subcomp=get_int_item(para, buf, &location);
										break;
									case 7:
										sf[k].iso_value=get_double_item(para, buf, &location);
										break;
									case 8:

										sf[k].method=get_int_item(para, buf, &location);
										break;
									case 9:
										for(j=0;j<3;j++)  {
											sf[k].point[j]=get_double_item(para, buf, &location);
										}
										break;

									case 10:

										sf[k].radius=get_double_item(para, buf, &location);
										break;
									case 11:
										for(j=0;j<3;j++)  {
											sf[k].length[j]=get_double_item(para, buf, &location);
										}
										break;
									case 12:
										for(j=0;j<10;j++) {

											sf[k].coef[j]=get_double_item(para, buf, &location);
										}
										break;
									case 13:
										sf[k].display_method=get_int_item(para, buf, &location);
										break;
									case 14:
										get_string_item(para, buf, &location, sf[k].color_comp_name);
										break;
									case 15:
										sf[k].color_comp=get_int_item(para, buf, &location);
										break;
									case 16:
										get_string_item(para, buf, &location, sf[k].color_subcomp_name);
										break;
									case 17:
										sf[k].color_subcomp=get_int_item(para, buf, &location);
										break;
									case 18:
										sf[k].isoline_number=get_int_item(para, buf, &location);
										/*                  fprintf(stderr, "isoline_number is %d\n", sf[k].isoline_number);
										 */
										break;
									case 19:
										sf[k].specified_color=get_double_item(para, buf, &location);

										break;
									case 20:
										get_string_item(para, buf, &location, para1);
										len_tmp=strlen(para1);
										for(ii=0;ii<len_tmp;ii++) {
											para1[ii]=toupper(para1[ii]);
										}
										if(strncmp(para1, "AVS", 3)==0)
											sf[k].output_type=1;
										else if(strncmp(para1, "GPPVIEW", 7)==0)
											sf[k].output_type=2;
										else if(strncmp(para1, "BMP", 3)==0)
											sf[k].output_type=3;
										else if(strncmp(para1, "VIS_FEMAP_NEUTRAL", 17)==0)
											sf[k].output_type=4;
										else if(strncmp(para1, "FSTR_FEMAP_NEUTRAL", 18)==0)
											sf[k].output_type=5;
										else if(strncmp(para1, "COMPLETE_AVS", 12)==0)
											sf[k].output_type=6;
										else if(strncmp(para1, "BIN_COMPLETE_AVS", 16)==0)
											sf[k].output_type=7;
										else if(strncmp(para1, "COMPLETE_REORDER_AVS", 20)==0)
											sf[k].output_type=8;
										else if(strncmp(para1, "COMPLETE_MICROAVS", 17)==0)
											sf[k].output_type=9;
										else if(strncmp(para1, "SEPARATE_COMPLETE_AVS", 21)==0)
											sf[k].output_type=10;
										else {
											fprintf(stderr, "ERROR: HEC-MW-VIS-E1001: output_type is not correct\n");
											HECMW_vis_print_exit("AVS or BMP or VIS_FEMAP_NEUTRAL or FSTR_FEMAP_NEUTRAL");
										}
										break;
									case 21:
										sf[k].range_output=1;
										get_string_item(para, buf, &location,  sf[k].range_filename);
										break;
									case 22:
										sf[k].normalize_flag=get_int_item(para, buf, &location);
										break;
									case 23:
										sr->xr=get_int_item(para, buf, &location);
										break;
									case 24:
										sr->yr=get_int_item(para, buf, &location);
										break;
									case 25:
										sr->num_of_lights=get_int_item(para, buf, &location);
										break;
									case 26:
										sr->light_point=(double *)HECMW_calloc(sr->num_of_lights*3, sizeof(double));
										if(sr->light_point==NULL)
											HECMW_vis_memory_exit("sr: light_point");

										for(i=0;i<sr->num_of_lights;i++) {
											sr->light_point[i*3]=get_double_item(para, buf, &location);
											sr->light_point[i*3+1]=get_double_item(para, buf, &location);
											sr->light_point[i*3+2]=get_double_item(para, buf, &location);
										}
										break;
									case 27:
										for(i=0;i<3;i++) {
											sr->view_point_d[i]=get_double_item(para, buf, &location);
										}

										break;
									case 28:
										for(i=0;i<3;i++) {
											sr->screen_point[i]=get_double_item(para, buf, &location);
										}
										break;
									case 29:
										for(i=0;i<3;i++) {
											sr->up[i]=get_double_item(para, buf, &location);
										}
										break;
									case 30:
										sr->k_ads[0]=get_double_item(para, buf, &location);
										break;
									case 31:
										sr->k_ads[1]=get_double_item(para, buf, &location);
										break;
									case 32:

										sr->k_ads[2]=get_double_item(para, buf, &location);
										break;
									case 33:
										sr->color_mapping_style=get_int_item(para, buf, &location);
										if(sr->color_mapping_style==1) {
											sr->interval_mapping_num=1;
										}
										else if(sr->color_mapping_style==2) {
											sr->interval_mapping_num=1;
										}
										else if(sr->color_mapping_style==4) {
											sr->interval_mapping_num=10;
										}
										break;
									case 34:

										sr->interval_mapping_num=get_int_item(para, buf, &location);
										break;
									case 35:
										if(sr->color_mapping_style==2) {
											sr->interval_point=(double *)HECMW_calloc(2, sizeof(double));

											if(sr->interval_point==NULL)
												HECMW_vis_memory_exit("interval_point");

											sr->interval_point[0]=get_double_item(para, buf, &location);
											sr->interval_point[1]=get_double_item(para, buf, &location);
										}
										else if(sr->color_mapping_style==3) {
											sr->interval_point=(double *)HECMW_calloc(2*(sr->interval_mapping_num+1), sizeof(double));
											if(sr->interval_point==NULL)
												HECMW_vis_memory_exit("interval_point");

											for(i=0;i<2*((sr->interval_mapping_num)+1);i++) {
												sr->interval_point[i]=get_double_item(para, buf, &location);

											}
										}
										break;
									case 36:
										sr->rotate_style=get_int_item(para, buf, &location);
										break;
									case 37:
										sr->rotate_num_of_frames=get_int_item(para, buf, &location);
										break;
									case 38:
										sr->color_mapping_bar_on=get_int_item(para, buf, &location);
										break;
									case 39:
										sr->scale_marking_on=get_int_item(para, buf, &location);
										break;
									case 40:
										for(i=0;i<3;i++) {
											sr->background_color[i]=get_double_item(para, buf, &location);
										}
										break;
									case 41:
										for(i=0;i<3;i++) {
											sr->font_color[i]=get_double_item(para, buf, &location);
										}
										break;
									case 42:
										sr->color_system_type=get_int_item(para, buf, &location);
										break;
									case 43:
										sr->font_size=get_double_item(para, buf, &location);
										break;
									case 44:
										sr->color_bar_style=get_int_item(para, buf, &location);
										break;
									case 45:
										sr->fixed_range_on=get_int_item(para, buf, &location);
										break;
									case 46:
										for(i=0;i<2;i++) {
											sr->range_value[i]=get_double_item(para, buf, &location);
										}

										break;
									case 47:
										sr->num_of_scale=get_int_item(para, buf, &location);
										break;
									case 48:
										sr->mark_0_on=get_int_item(para, buf, &location);
										break;
									case 49:
										sr->transfer_function_style=get_int_item(para, buf, &location);

										if(sr->transfer_function_style!=1)
											sr->opa_value=0.0;
										if((sr->transfer_function_style!=3) && (sr->transfer_function_style!=4)) {
											sr->num_of_features=0;
										}
										break;
									case 50:
										sr->opa_value=get_double_item(para, buf, &location);

										break;
									case 51:
										sr->num_of_features=get_int_item(para, buf, &location);
										break;
									case 52:
										if(sr->transfer_function_style==3)
											sr->fea_point=(double *)HECMW_calloc(sr->num_of_features*3, sizeof(double));

										else if(sr->transfer_function_style==4)
											sr->fea_point=(double *)HECMW_calloc(sr->num_of_features*3, sizeof(double));
										if(sr->fea_point==NULL)
											HECMW_vis_memory_exit("sr: fea_point");

										if(sr->transfer_function_style==3) {
											for(i=0;i<sr->num_of_features*3;i++) {
												sr->fea_point[i]=get_double_item(para, buf, &location);

											}
										}
										if(sr->transfer_function_style==4) {
											for(i=0;i<sr->num_of_features*3;i++) {
												sr->fea_point[i]=get_double_item(para, buf, &location);

											}
										}
										break;
									case 53:
										get_string_item(para, buf, &location, sr->name_lookup);
										break;
									case 54:
										sr->histogram_on=get_int_item(para, buf, &location);
										break;
									case 55:
										sr->boundary_line_on=get_int_item(para, buf, &location);
										break;
									case 56:
										for(i=0;i<3;i++) {
											sr->isoline_color[i]=get_double_item(para, buf, &location);
										}
										break;
									case 57:
										sr->time_mark_on=get_int_item(para, buf, &location);
										break;
									case 58:
										sr->fixed_scale_mark=get_int_item(para, buf, &location);
										break;
									case 59:
										sf[k].deform_display_on=get_int_item(para, buf, &location);
										break;
									case 60:
										sf[k].disp_scale=get_double_item(para, buf, &location);
										break;
									case 61:
										sf[k].initial_style=get_int_item(para, buf, &location);
										break;
									case 62:
										sf[k].deform_style=get_int_item(para, buf, &location);
										break;
									case 63:
										get_string_item(para, buf, &location,sf[k].disp_comp_name);
										break;
									case 64:
										sf[k].disp_comp=get_int_item(para, buf, &location);
										break;
									case 65:
										for(i=0;i<3;i++) {
											sf[k].initial_line_color[i]=get_double_item(para, buf, &location);
										}
										break;
									case 66:
										for(i=0;i<3;i++) {
											sf[k].deform_line_color[i]=get_double_item(para, buf, &location);
										}
										break;
									case 67:
										sr->deform_num_of_frames=get_int_item(para, buf, &location);
										break;
									case 68:
										sr->smooth_shading=get_int_item(para, buf, &location);
										break;
									case 69:
										sf[k].real_disp_scale=get_double_item(para, buf, &location);
										break;
									case 70:
										for(i=0;i<6;i++) {
											sr->fixed_mesh_range[i]=get_double_item(para, buf, &location);
										}
										break;
									case 71:
										sr->start_time=get_double_item(para, buf, &location);
										break;
									case 72:
										sr->time_interval=get_double_item(para, buf, &location);
										break;
									}
								}
								while (cont_flag) {
									if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
										if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
									}
									else
										cont_flag=0;
								}
								if(cont_flag==0)
									break;
							}
						/* check the parameters for the surface k */
						if(stat_para[0]==0)
							sf[k].surface_style=1;
						if((sf[k].surface_style<1) || (sf[k].surface_style>3))
							HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1002: surface_style only can be 1, 2 or 3");

						if(sf[k].surface_style==1) { /* boundary and group surface */
							/*		  if(stat_para[1]==0) {
                      fprintf(stderr, "Please specify the group name for the boundary surace generateion; if to find the whole boundary, please input boundary\n");
                      exit(0);
                  }
							 */
							if(stat_para[2]==0)  { /* default value */
								sf[k].defined_style=2;
							}
						}
						if(sf[k].surface_style==2) { /* iso-surfaces */
							if((stat_para[3]==0) && (stat_para[4]==0))
								sf[k].data_comp=0;
							if((stat_para[5]==0) && (stat_para[6]==0))
								sf[k].data_subcomp=1;
							if(stat_para[7]==0)
								HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1002: iso_value should be defined for isosurfaces");
						}
						if(sf[k].surface_style==3) { /* arbitrary surfaces defined by equation */
							if(stat_para[8]==0)
								sf[k].method=5;
							if(sf[k].method==1) {
								if(stat_para[9]==0) {
									sf[k].point[0]=0.0, sf[k].point[1]=0.0; sf[k].point[2]=0.0;
									if(mynode==0)
										fprintf(stderr, "The default center point (0.0, 0.0, 0.0) is used\n");

								}
								if(stat_para[10]==0) {
									sf[k].radius=1.0;
									if(mynode==0)
										fprintf(stderr, "The default radius 1.0 is used\n");

								}
							}
							else if((sf[k].method==2) || (sf[k].method==3) || (sf[k].method==4)) {
								if(stat_para[9]==0) {
									sf[k].point[0]=0.0; sf[k].point[1]=0.0; sf[k].point[2]=0.0;
									if(mynode==0)
										fprintf(stderr, "The default center point (0.0, 0.0, 0.0) is used\n");

								}
								if(stat_para[11]==0) {
									sf[k].length[0]=1.0; sf[k].length[1]=1.0; sf[k].length[2]=1.0;
									if(mynode==0)
										fprintf(stderr, "The default length (1.0, 1.0, 1.0) is used\n");
								}
							}
							else if (sf[k].method==5) {
								if(stat_para[12]==0) {
									if(mynode==0)
										fprintf(stderr, "ERROR: HEC-MW-VIS-E1003:The 10 coefficients of the equation should be defined \n");
									HECMW_vis_print_exit("Please re-input and run again");
								}
							}
						}

						if(stat_para[13]==0)
							sf[k].display_method=1;
						if((sf[k].display_method<=0) || (sf[k].display_method>5)) {
							if(mynode==0)
								fprintf(stderr, "ERROR: HEC-MW-VIS-E1004:display_method is not in the reasonable range\n");
							HECMW_vis_print_exit("Please re-input and run again");
						}
						if(sf[k].display_method!=4) {
							if((stat_para[14]==0) && (stat_para[15]==0))
								sf[k].color_comp=0;
							if((stat_para[16]==0) && (stat_para[17]==0))
								sf[k].color_subcomp=1;
						}
						if((sf[k].display_method==2) || (sf[k].display_method==3) || (sf[k].display_method==5)) {
							if(stat_para[18]==0) {
								if(mynode==0)
									fprintf(stderr, "No value for the number of isolines. Now using the default value 10");
								sf[k].isoline_number=10;
							}
						}
						if(sf[k].display_method==4) {
							if(stat_para[19]==0) {
								if(mynode==0) {
									fprintf(stderr, "The number of specified_color has not been defined\n");
									fprintf(stderr, "Now using the default value 0.5\n");
								}
								sf[k].specified_color=0.5;
							}
						}
						if(stat_para[20]==0)
							sf[k].output_type=1;
						if((sf[k].output_type<1) || (sf[k].output_type>10)) {
							if(mynode==0)
								fprintf(stderr, "the output_type only can be 1 -- 10\n");
							HECMW_vis_print_exit( "pls input and run again");
						}
						if(stat_para[22]==0)
							sf[k].normalize_flag = 0;
						if(stat_para[59]==0)
							sf[k].deform_display_on=0;
						if((sf[k].deform_display_on<0) || (sf[k].deform_display_on>1)) {
							fprintf(stderr, "ERROR: HEC-MW-VIS-E1037: deform_display_on should be 0 or 1");
							HECMW_vis_print_exit("Please re-input again");
						}
						if(sf[k].deform_display_on==0)
							sr->deform_num_of_frames=1;
						if((sf[k].deform_display_on==1) && (stat_para[67]==0))
							sr->deform_num_of_frames=8;
						/*		  if((sf[k].deform_display_on==1) && (sr->deform_num_of_frames==1)) {
			  fprintf(stderr, "For deformation display, the deform_num_of_frames should be greater than 1\n");
			  fprintf(stderr, "The default value 8 will be used\n");
			  sr->deform_num_of_frames = 8;
		  }
						 */
						if(stat_para[60]==0)
							sf[k].disp_scale=-1.0;
						/*		  if(sf[k].disp_scale<0.0) {
			  fprintf(stderr, "ERROR: HEC-MW-VIS-E1037: deform_scale should be greater than 0.0");
			  HECMW_vis_print_exit("Please re-input again");
		  }
						 */
						if(stat_para[61]==0)
							sf[k].initial_style=4;
						if(stat_para[62]==0)
							sf[k].deform_style=1;
						if((sf[k].initial_style<0) || (sf[k].initial_style>4)){
							fprintf(stderr, "ERROR: HEC-MW-VIS-E1037: initial_display_style should be in (0, 4)");
							HECMW_vis_print_exit("Please re-input again");
						}
						if((sf[k].deform_style<0) || (sf[k].deform_style>4)){
							fprintf(stderr, "ERROR: HEC-MW-VIS-E1037: deform_display_style should be in (0, 4)");
							HECMW_vis_print_exit("Please re-input again");
						}
						if((stat_para[63]==0) && (stat_para[64]==0))
							sprintf(sf[k].disp_comp_name, "%s", "DISPLACEMENT");
						if(stat_para[65]==0) {
							sf[k].initial_line_color[0]=0.2;
							sf[k].initial_line_color[1]=0.2;
							sf[k].initial_line_color[2]=1.0;
						}
						if(stat_para[66]==0) {
							sf[k].deform_line_color[0]=1.0;
							sf[k].deform_line_color[1]=1.0;
							sf[k].deform_line_color[2]=0.0;
						}
						for(i=0;i<3;i++) {
							if(sf[k].initial_line_color[i]<0.0)
								sf[k].initial_line_color[i]=0.0;
							if(sf[k].initial_line_color[i]>1.0)
								sf[k].initial_line_color[i]=1.0;
							if(sf[k].deform_line_color[i]<0.0)
								sf[k].deform_line_color[i]=0.0;
							if(sf[k].deform_line_color[i]>1.0)
								sf[k].deform_line_color[i]=1.0;
						}


						if(cont_flag==0)
							break;
						else if(get_keyword_visual(buf)>0) {
							visual_method=get_keyword_visual(buf);
							break;
						}
					}
			}  /* end of loop k */
			/*     fprintf(stderr, "the group name of three is 1: %s 2: %s 3:%s\n", sf[1].group_name, sf[2].group_name, sf[3].group_name);
			 */

			/*  start reading rendering parameters if output_type==BMP */

			/*  check and set default */
			if(sf[1].output_type==3) {
				if(stat_para[23]==0) {
					/* use default value 256 */
					if(mynode==0)
						fprintf(stderr, "No value for xr. Now use the default value 512\n");
					sr->xr=512;
				}
				/* check xr whether can be divided by 8 */
				sr->xr=(int)(sr->xr/8)*8;

				if(stat_para[24]==0) {
					if(mynode==0)
						fprintf(stderr, "No value for yr. Now use the default value 512\n");
					sr->yr=512;
				}
				if(sr->xr<=20) {
					HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1005: The x_resolution should be greater than 20");
				}
				if(sr->yr<=20) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1006: The y_resolution should be greater than 20\n");
					HECMW_vis_print_exit("Please re-input a correct one");
				}
				if(stat_para[25]==0) {
					if(mynode==0)
						fprintf(stderr, "No value for num_of_lights. Now use the default value 1\n");
					sr->num_of_lights=1;
					stat_para[25]=1;
				}
				if((stat_para[25]==1) && (sr->num_of_lights<=0)) {
					HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1007: The number of light sources should be greater than 0");
				}
				if(stat_para[29]==0) {
					stat_para[29]=1;
					if(mynode==0)
						fprintf(stderr, "No value for up_direction. The default value (0.0, 0.0, 1.0) is used\n");
					sr->up[0]=0.0;
					sr->up[1]=0.0;
					sr->up[2]=1.0;
				}
				if(stat_para[30]==0) {
					stat_para[30]=1;
					sr->k_ads[0]=0.5;
				}
				if((stat_para[30]==1) && (sr->k_ads[0]<0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1008: The ambient_coef is not correct. Should be >=0.0\n");
					HECMW_vis_print_exit("Please re-input the ambient_coef in your control file");
				}
				if(stat_para[31]==0) {
					stat_para[31]=1;
					sr->k_ads[1]=0.5;
				}
				if((stat_para[31]==1) && (sr->k_ads[1]<0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1009: The diffuse_coef is not correct. Should be >=0.0\n");
					HECMW_vis_print_exit("Please re-input the diffuse_coef in your control file");
				}
				if(stat_para[32]==0) {
					stat_para[32]=1;
					sr->k_ads[2]=0.6;
				}
				if((stat_para[32]==1) && (sr->k_ads[2]<0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1010: The specular_coef is not correct. Should be >=0.0\n");
					HECMW_vis_print_exit("Please re-input the specular_coef in your control file");
				}
				if(stat_para[33]==0) {
					stat_para[33]=1;
					sr->color_mapping_style=1;
				}
				if((sr->color_mapping_style<1) || (sr->color_mapping_style>4)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1011: color_mapping_style should be between 1 and 4\n");
					HECMW_vis_print_exit("Please re-input it and run again\n");
				}
				if((sr->color_mapping_style==3) && (stat_para[34]==0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1012:For color_mapping_style 3, interval_mapping_num is required\n");
					HECMW_vis_print_exit("Please re-input the value of interval_mapping_num");
				}
				if((stat_para[34]==1) && (sr->interval_mapping_num<=0)) {
					HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1013: For color_mapping_style 3, the interval_mapping_num should be greater than 0");
				}
				if(((sr->color_mapping_style==2) || (sr->color_mapping_style==3)) && (stat_para[35]==0)) {
					fprintf(stderr,"ERROR: HEC-MW-VIS-E1014: For color_mapping_style =2 or 3, the interval_point should be defined\n");
					HECMW_vis_print_exit("Please re-input the interval_point");
				}
				if(stat_para[36]==0) {
					sr->rotate_style=0;
					sr->rotate_num_of_frames=1;
				}
				if((sr->rotate_style>=1) && (sr->rotate_style<=4)) {
					if(stat_para[37]==0)
						sr->rotate_num_of_frames=8;
					if(sr->rotate_num_of_frames<=0) {
						fprintf(stderr, "ERROR: HEC-MW-VIS-E1015: The parameter rotate_num_of_frames cannot be less than 1.\n");
						HECMW_vis_print_exit("Please re-input and run again");
					}
					if(sr->rotate_style==4)
						sr->rotate_num_of_frames=8;
				}

				if(stat_para[38]==0)
					sr->color_mapping_bar_on=0;
				if((sr->color_mapping_bar_on<0) || (sr->color_mapping_bar_on>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1016: color_mapping_bar_on parameter only can be defined as 0 or 1\n");
					HECMW_vis_print_exit("Please re-input it and run again");
				}
				if(stat_para[39]==0)
					sr->scale_marking_on=0;
				if((sr->scale_marking_on<0) || (sr->scale_marking_on>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1017: scale_marking_on parameter only can be defined as 0 or 1\n");
					HECMW_vis_print_exit("Please re-input it and run again");
				}

				if((sr->color_mapping_bar_on==1) && (sr->xr<40)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1018: x resolution should be larger than 40 for adding color_mapping_bar case\n");
					HECMW_vis_print_exit( "please reinput x resolution");
				}
				if((sr->scale_marking_on==1) && (sr->xr<65)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1019: x resolution should be larger than 65 for adding color_mapping and scale marking case\n");
					HECMW_vis_print_exit( "please re-input x resolution again");
				}
				if(stat_para[40]==0)
					sr->background_color[0]=sr->background_color[1]=sr->background_color[2]=0.0;
				if(stat_para[41]==0)
					sr->font_color[0]=sr->font_color[1]=sr->font_color[2]=1.0;
				if(stat_para[42]==0)
					sr->color_system_type=1;
				if((sr->color_system_type<=0) || (sr->color_system_type>3)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1020:color_system_type should be between 1 and 3\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if((sr->background_color[0]<-EPSILON) || (sr->background_color[0]>1.0+EPSILON) ||
						(sr->background_color[1]<-EPSILON) || (sr->background_color[1]>1.0+EPSILON) ||
						(sr->background_color[2]<-EPSILON) || (sr->background_color[2]>1.0+EPSILON)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1021:The background color should be in the interval of  (0.0, 1.0)\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if((sr->font_color[0]<-EPSILON) || (sr->font_color[0]>1.0+EPSILON) ||
						(sr->font_color[1]<-EPSILON) || (sr->font_color[1]>1.0+EPSILON) ||
						(sr->font_color[2]<-EPSILON) || (sr->font_color[2]>1.0+EPSILON)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1022: The font color should be in the interval of  (0.0, 1.0)\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if(stat_para[43]==0)
					sr->font_size=1.0;
				if(sr->font_size<1.0-EPSILON) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1023: font_size paramters cannot be between 1.0 and 4.0\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if(sr->font_size>4.0)
					sr->font_size=4.0;
				if(stat_para[44]==0)
					sr->color_bar_style=2;
				if((sr->color_bar_style<1) || (sr->color_bar_style>2)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1024: color_bar_style only can be 1 or 2\n");
					HECMW_vis_print_exit( "Please input and run again");
				}
				if(stat_para[45]==0)
					sr->fixed_range_on=0;
				if((sr->fixed_range_on<0) || (sr->fixed_range_on>1)) {
					fprintf(stderr, "fixed_range_on only can be 0 or 1\n");
					HECMW_vis_print_exit("Please input and run again");
				}
				if(stat_para[47]==0)
					sr->num_of_scale=3;
				if(sr->num_of_scale<=0) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1025: num_of_scale only can be greater than 0\n");
					HECMW_vis_print_exit("Please input and run again");
				}
				if((sr->fixed_range_on==1) && (stat_para[46]==0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1026: range_value is required for fixed_range_on style\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if(stat_para[48]==0)
					sr->mark_0_on=0;
				if((sr->mark_0_on<0) || (sr->mark_0_on>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1027: mark_0_on only can be 0 or 1\n");
					HECMW_vis_print_exit("Please input and run again");
				}

				if(stat_para[49]==0) {
					stat_para[49]=1;
					sr->transfer_function_style=1;
					stat_para[50]=1;
					sr->opa_value=1.0;
				}
				if((sr->transfer_function_style<1) || (sr->transfer_function_style>9)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1028: transfer_function_style should be between 1 and 8\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if((stat_para[50]==1) && (sr->opa_value<0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1029: opacity_value cannot be less than 0.0\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if(((sr->transfer_function_style==3) || (sr->transfer_function_style==4)) && (stat_para[51]==0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1030:When transfer_function_style =3 or 4, num_of_features should be defined\n");
					HECMW_vis_print_exit("Please re-input the num_of_features");
				}
				if((stat_para[51]==1) && (sr->num_of_features<=0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1031: When transfer_function_style =3 or 4, num_of_features should be greater than 0\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if((stat_para[51]==1) && (stat_para[52]==0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1032:For transfer_function_style =3 or 4, feature_points should be defined\n");
					HECMW_vis_print_exit("Please re-input your control file");
				}
				if((sr->transfer_function_style==8) && (stat_para[53]==0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1033:For transfer_function_style=8, lookup_filename should be specified\n");
					HECMW_vis_print_exit("Please re-input the filename");
				}
				if(stat_para[54]==0)
					sr->histogram_on=0;
				if((sr->histogram_on<0) || (sr->histogram_on>2)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1034: histogram_on value should be 0,1, or 2\n");
					HECMW_vis_print_exit("Please re-input again");
				}
				if(stat_para[55]==0)
					sr->boundary_line_on=0;
				if((sr->boundary_line_on<0) || (sr->boundary_line_on>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1035: histogram_on value should be 0 or 1");
					HECMW_vis_print_exit("Please re-input again");
				}
				if(stat_para[56]==0) {
					sr->isoline_color[0]=0.0;
					sr->isoline_color[1]=0.0;
					sr->isoline_color[2]=0.0;
				}
				if(stat_para[57]==0)
					sr->time_mark_on=0;
				if((sr->time_mark_on<0) || (sr->time_mark_on>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1036: time_mark_on value should be 0 or 1");
					HECMW_vis_print_exit("Please re-input again");
				}
				if(stat_para[58]==0)
					sr->fixed_scale_mark=0;
				if((sr->fixed_scale_mark<0) || (sr->fixed_scale_mark>1)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1037: fixed_scale_mark value should be 0 or 1");
					HECMW_vis_print_exit("Please re-input again");
				}
				if(stat_para[68]==0)
					sr->smooth_shading=0;

				/*		  fprintf(stderr, "resolution is %d  %d\n", sr->xr, sr->yr);
		  fprintf(stderr, "color_mapping style is %d\n", sr->color_mapping_style);
		  fprintf(stderr, "view_point==%lf %lf %lf\n", sr->view_point_d[0], sr->view_point_d[1], sr->view_point_d[2]);
				 */
			}
			/* copy the parameters for PSF into psf link */
			p1=psf->next_psf;
			p1->num_of_psf=sf[0].surface_style;
			p1->sr=sr;
			p1->sf=sf;
			/*		  fprintf(stderr, "******** color_comp_name =%s\n", sf[1].color_comp_name);
		  fprintf(stderr, "******** surface_style=%d\n", sf[1].surface_style);
		  fprintf(stderr, "*********current PSF number= %d\n", psf->num_of_psf);
			 */

			for(i=0;i<NUM_CONTROL_PSF;i++)
				p1->stat_para[i]=stat_para[i];
			for(i=0;i<NUM_CONTROL_PSF;i++)
				HECMW_free(parameters[i]);
			HECMW_free(parameters);
			HECMW_free(stat_para);
			HECMW_free(len_para);

		} /*end if visual_method=1 */
		else if((visual_method==2) && (cont_flag==1)) {
			visual_method=0;
			vr=(Parameter_vr *)HECMW_malloc(sizeof(Parameter_vr));
			if(vr==NULL)
				HECMW_vis_memory_exit("PVR parameters: vr");
			vr->projection_style=1;


			for(i=0;i<NUM_CONTROL_PVR;i++) {
				parameters=(char **)HECMW_calloc(NUM_CONTROL_PVR, sizeof(char *));
				len_para=(int *)HECMW_calloc(NUM_CONTROL_PVR, sizeof(int));
				stat_para=(int *)HECMW_calloc(NUM_CONTROL_PVR, sizeof(int));
				parameters[i]=(char *)HECMW_calloc(128, sizeof(char));
				if(parameters[i]==NULL)
					HECMW_vis_memory_exit("tempory variable parameters");
			}
			vr->opa_value=0.0;
			vr->color_comp=-1;
			vr->color_subcomp=-1;
			strcpy(parameters[0], "maximum_refinement");
			len_para[0]=18;
			strcpy(parameters[1], "x_resolution");
			strcpy(parameters[2], "y_resolution");
			len_para[1]=len_para[2]=12;
			strcpy(parameters[3], "num_of_lights");
			len_para[3]=13;
			strcpy(parameters[4], "position_of_lights");
			len_para[4]=18;
			strcpy(parameters[5], "viewpoint");
			len_para[5]=9;
			strcpy(parameters[6], "look_at_point");
			len_para[6]=13;
			strcpy(parameters[7], "up_direction");
			len_para[7]=12;
			strcpy(parameters[8], "ambient_coef");
			len_para[8]=12;
			strcpy(parameters[9], "diffuse_coef");
			len_para[9]=12;
			strcpy(parameters[10], "specular_coef");
			len_para[10]=13;
			strcpy(parameters[11], "surface_on");
			len_para[11]=10;
			strcpy(parameters[12], "surface_opacity");
			len_para[12]=15;
			strcpy(parameters[13], "color_mapping_style");
			len_para[13]=19;
			strcpy(parameters[14], "interval_mapping_num");
			len_para[14]=20;
			strcpy(parameters[15], "interval_mapping");
			len_para[15]=16;
			strcpy(parameters[16], "opacity_mapping_style");
			len_para[16]=21;
			strcpy(parameters[17], "opacity_value");
			len_para[17]=13;
			strcpy(parameters[18], "num_of_features");
			len_para[18]=14;
			strcpy(parameters[19], "feature_points");
			len_para[19]=14;
			strcpy(parameters[20], "lookup_filename");
			len_para[20]=15;
			strcpy(parameters[21], "rotate_style");
			len_para[21]=12;
			strcpy(parameters[22], "voxel_filename");
			len_para[22]=10;
			strcpy(parameters[23], "color_mapping_bar_on");
			len_para[23]=17;
			strcpy(parameters[24], "scale_marking_on");
			len_para[24]=15;
			strcpy(parameters[25], "color_comp_name");
			len_para[25]=14;
			strcpy(parameters[26], "color_subcomp_name");
			len_para[26]=17;
			strcpy(parameters[27], "n_voxel_x");
			len_para[27]=9;
			strcpy(parameters[28], "n_voxel_y");
			len_para[28]=9;
			strcpy(parameters[29], "n_voxel_z");
			len_para[29]=9;
			strcpy(parameters[30], "surface_filename");
			len_para[30]=12;
			strcpy(parameters[31], "num_of_frames");
			len_para[31]=12;
			strcpy(parameters[32], "background_color");
			len_para[32]=12;
			strcpy(parameters[33], "font_color");
			len_para[33]=9;
			strcpy(parameters[34], "color_system_type");
			len_para[34]=12;
			strcpy(parameters[35], "font_size");
			len_para[35]=8;
			strcpy(parameters[36], "color_bar_style");
			len_para[36]=13;
			strcpy(parameters[37], "fixed_range_on");
			len_para[37]=11;
			strcpy(parameters[38], "range_value");
			len_para[38]=11;
			strcpy(parameters[39], "num_of_scale");
			len_para[39]=11;
			strcpy(parameters[40], "mark_0_on");
			len_para[40]=9;
			strcpy(parameters[41], "remove_0_display_on");
			len_para[41]=17;
			strcpy(parameters[42], "x_specified_level");
			len_para[42]=15;
			strcpy(parameters[43], "y_specified_level");
			len_para[43]=15;
			strcpy(parameters[44], "z_specified_level");
			len_para[44]=15;
			strcpy(parameters[45], "histogram_on");
			len_para[45]=10;
			strcpy(parameters[46], "display_range");
			len_para[46]=13;
			strcpy(parameters[47], "time_mark_on");
			len_para[47]=10;
			strcpy(parameters[48], "fixed_scale_mark");
			len_para[48]=15;
			strcpy(parameters[49], "color_comp");
			len_para[49]=9;
			strcpy(parameters[50], "color_subcomp");
			len_para[50]=12;


			for(i=0;i<NUM_CONTROL_PVR;i++) {
				stat_para[i]=0;
			}
			vr->projection_style=1;
			vr->color_comp=-1;
			vr->color_subcomp=-1;
			sprintf(vr->color_comp_name, "%s", "NULL");
			sprintf(vr->color_subcomp_name, "%s", "NULL");

			if(cont_flag==1)
				while(get_keyword_visual(buf)==0) {
					hit=-1;
					location=get_keyword_item(buf, para);
					len_tmp=strlen(para);
					for(ii=0;ii<len_tmp;ii++) {
						para[ii]=tolower(para[ii]);
					}
					for(i=0;i<NUM_CONTROL_PVR;i++) {
						if((strncmp(para, parameters[i], len_para[i]))==0) {
							hit=i;
							stat_para[i]=1;
							break;
						}
					}
					/*	 fprintf(stderr, "para=%s hit=%d\n", para, hit);
					 */

					if((hit>=0) && (hit<NUM_CONTROL_PVR)) {
						switch (hit) {
						case 0:
							vr->max_level=get_int_item(para, buf, &location);
							break;
						case 1:
							vr->xr=get_int_item(para, buf, &location);
							break;
						case 2:
							vr->yr=get_int_item(para, buf, &location);
							break;
						case 3:
							vr->num_of_lights=get_int_item(para, buf, &location);
							break;
						case 4:
							vr->light_point=(double *)HECMW_calloc(vr->num_of_lights*3, sizeof(double));
							if(vr->light_point==NULL)
								HECMW_vis_memory_exit("Parameter vr: light_point");

							for(i=0;i<vr->num_of_lights;i++) {
								vr->light_point[i*3]=get_double_item(para, buf, &location);
								vr->light_point[i*3+1]=get_double_item(para, buf, &location);
								vr->light_point[i*3+2]=get_double_item(para, buf, &location);
							}
							break;
						case 5:
							for(i=0;i<3;i++) {
								vr->view_point_d[i]=get_double_item(para, buf, &location);
							}

							break;
						case 6:
							for(i=0;i<3;i++) {
								vr->screen_point[i]=get_double_item(para, buf, &location);
							}
							break;
						case 7:
							for(i=0;i<3;i++) {

								vr->up[i]=get_double_item(para, buf, &location);
							}
							break;
						case 8:

							vr->k_ads[0]=get_double_item(para, buf, &location);
							break;
						case 9:
							vr->k_ads[1]=get_double_item(para, buf, &location);
							break;
						case 10:

							vr->k_ads[2]=get_double_item(para, buf, &location);
							break;
						case 11:

							break;
						case 12:
							break;
						case 13:
							vr->color_mapping_style=get_int_item(para, buf, &location);
							if(vr->color_mapping_style==1) {
								vr->interval_mapping_num=1;
							}
							else if(vr->color_mapping_style==2) {
								vr->interval_mapping_num=1;
							}
							else if(vr->color_mapping_style==4) {
								vr->interval_mapping_num=10;
							}
							break;
						case 14:

							vr->interval_mapping_num=get_int_item(para, buf, &location);
							break;
						case 15:
							if(vr->color_mapping_style==2) {
								vr->interval_point=(double *)HECMW_calloc(2, sizeof(double));

								if(vr->interval_point==NULL)
									HECMW_vis_memory_exit("interval point");

								vr->interval_point[0]=get_double_item(para, buf, &location);
								vr->interval_point[1]=get_double_item(para, buf, &location);
							}
							else if(vr->color_mapping_style==3) {
								vr->interval_point=(double *)HECMW_calloc(2*(vr->interval_mapping_num+1), sizeof(double));
								if(vr->interval_point==NULL)
									HECMW_vis_memory_exit("interval point");

								for(i=0;i<2*((vr->interval_mapping_num)+1);i++) {
									vr->interval_point[i]=get_double_item(para, buf, &location);

								}
							}
							break;
						case 16:
							vr->transfer_function_style=get_int_item(para, buf, &location);
							if(vr->transfer_function_style!=1)
								vr->opa_value=0.0;
							if((vr->transfer_function_style!=3) && (vr->transfer_function_style!=4)) {
								vr->num_of_features=0;
							}
							break;
						case 17:
							vr->opa_value=get_double_item(para, buf, &location);

							break;
						case 18:
							vr->num_of_features=get_int_item(para, buf, &location);
							break;
						case 19:
							if(vr->transfer_function_style==3)
								vr->fea_point=(double *)HECMW_calloc(vr->num_of_features*3, sizeof(double));

							else if(vr->transfer_function_style==4)
								vr->fea_point=(double *)HECMW_calloc(vr->num_of_features*3, sizeof(double));
							if(vr->fea_point==NULL)
								HECMW_vis_memory_exit("vr: fea_point");

							if(vr->transfer_function_style==3) {
								for(i=0;i<vr->num_of_features*3;i++) {
									vr->fea_point[i]=get_double_item(para, buf, &location);

								}
							}
							if(vr->transfer_function_style==4) {
								for(i=0;i<vr->num_of_features*3;i++) {
									vr->fea_point[i]=get_double_item(para, buf, &location);

								}
							}
							break;
						case 20:
							get_string_item(para, buf, &location, vr->name_lookup);
							break;
						case 21:
							vr->rotate_style=get_int_item(para, buf, &location);
							break;
						case 22:
							get_string_item(para, buf, &location, vr->name_voxelfile);
							break;
						case 23:
							vr->color_mapping_bar_on=get_int_item(para, buf, &location);
							break;
						case 24:
							vr->scale_marking_on=get_int_item(para, buf, &location);
							break;
						case 25:
							get_string_item(para, buf, &location, vr->color_comp_name);
							break;
						case 26:
							get_string_item(para, buf, &location, vr->color_subcomp_name);
							break;
						case 27:
							vr->nv_xyz[0]=get_int_item(para, buf, &location);
							break;
						case 28:
							vr->nv_xyz[1]=get_int_item(para, buf, &location);
							break;
						case 29:
							vr->nv_xyz[2]=get_int_item(para, buf, &location);
							break;
						case 30:
							break;
						case 31:
							vr->num_of_frames=get_int_item(para, buf, &location);
							break;
						case 32:
							for(i=0;i<3;i++) {
								vr->background_color[i]=get_double_item(para, buf, &location);
							}
							break;
						case 33:
							for(i=0;i<3;i++) {
								vr->font_color[i]=get_double_item(para, buf, &location);
							}
							break;
						case 34:
							vr->color_system_type=get_int_item(para, buf, &location);
							break;
						case 35:
							vr->font_size=get_double_item(para, buf, &location);
							break;
						case 36:
							vr->color_bar_style=get_int_item(para, buf, &location);
							break;
						case 37:
							vr->fixed_range_on=get_int_item(para, buf, &location);
							break;
						case 38:
							for(i=0;i<2;i++) {
								vr->range_value[i]=get_double_item(para, buf, &location);
							}

							break;
						case 39:
							vr->num_of_scale=get_int_item(para, buf, &location);
							break;
						case 40:
							vr->mark_0_on=get_int_item(para, buf, &location);
							break;
						case 41:
							vr->remove_0_display_on=get_int_item(para, buf, &location);
							break;
						case 42:
							vr->specified_level[0]=get_int_item(para, buf, &location);
							break;
						case 43:
							vr->specified_level[1]=get_int_item(para, buf, &location);
							break;
						case 44:
							vr->specified_level[2]=get_int_item(para, buf, &location);
							break;
						case 45:
							vr->histogram_on=get_int_item(para, buf, &location);
							break;
						case 46:
							for(i=0;i<6;i++) {
								vr->display_range[i]=get_double_item(para, buf, &location);
							}
							break;
						case 47:
							vr->time_mark_on=get_int_item(para, buf, &location);
							break;
						case 48:
							vr->fixed_scale_mark=get_int_item(para, buf, &location);
							break;
						case 49:
							vr->color_comp=get_int_item(para, buf, &location);
							break;
						case 50:
							vr->color_subcomp=get_int_item(para, buf, &location);
							break;
						}
					}
					while (cont_flag) {
						if(fgets(buf, MAX_LINE_LEN, fp) != NULL) {
							if ((is_blank_line(buf)==0) && (is_comment_line(buf)==0)) break;
						}
						else
							cont_flag=0;
					}
					if(cont_flag==0)
						break;
					else if(get_keyword_visual(buf)>0) {
						visual_method=get_keyword_visual(buf);
						break;
					}
				}
			/* check the parameters */
			if(stat_para[0]==0) {
				/*  set default value */
				vr->max_level=100;
			}
			if((vr->max_level<=0) && (stat_para[0]==1)){
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1038: maximum_refinement parameter should be greater than 0\n");
				HECMW_vis_print_exit("Please re-input the value");
			}
			if(stat_para[1]==0) {
				/* use default value 256 */
				vr->xr=256;
			}
			/* check xr whether can be divided by 8 */
			vr->xr=(int)(vr->xr/8)*8;

			if(stat_para[2]==0) {
				vr->yr=256;
			}
			if(vr->xr<=20) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1005: The x_resolution should be greater than 20\n");
				HECMW_vis_print_exit("Please re-input a correct one");
			}
			if(vr->yr<=20) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1006: The y_resolution should be greater than 20\n");
				HECMW_vis_print_exit("Please re-input a correct one");
			}
			if(stat_para[3]==0) {
				vr->num_of_lights=1;
				stat_para[3]=1;
			}
			if((stat_para[3]==1) && (vr->num_of_lights<=0)) {
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1007: The number of light sources should be greater than 0");
			}
			if(stat_para[7]==0) {
				stat_para[7]=1;
				vr->up[0]=0.0;
				vr->up[1]=1.0;
				vr->up[2]=0.0;
			}
			if(stat_para[8]==0) {
				stat_para[8]=1;
				vr->k_ads[0]=0.5;
			}
			if((stat_para[8]==1) && (vr->k_ads[0]<0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1008: The ambient_coef is not correct. Should be >=0.0\n");
				HECMW_vis_print_exit("Please re-input the ambient_coef in your control file");
			}
			if(stat_para[9]==0) {
				stat_para[9]=1;
				vr->k_ads[1]=0.8;
			}
			if((stat_para[9]==1) && (vr->k_ads[1]<0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1009: The diffuse_coef is not correct. Should be >=0.0\n");
				HECMW_vis_print_exit("Please re-input the diffuse_coef in your control file");
			}
			if(stat_para[10]==0) {
				stat_para[10]=1;
				vr->k_ads[2]=0.6;
			}
			if((stat_para[10]==1) && (vr->k_ads[2]<0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1010: The specular_coef is not correct. Should be >=0.0\n");
				HECMW_vis_print_exit("Please re-input the specular_coef in your control file");
			}
			/*		  if(stat_para[11]==0) {
			  stat_para[11]=1;
			  vr->surface_on=0;
			  pvr->surface_on=0;
		  }

        	  if((stat_para[11]==1) && ((vr->surface_on<0) || (vr->surface_on>1))) {
			  fprintf(stderr, "surface_on parameter only can be defined as 0 or 1\n");
			  fprintf(stderr, "Please re-input it and run again\n");
			  exit(0);
		  }
		  if((vr->surface_on==1) && (stat_para[12]==0))
			  vr->surface_opacity=0.2;
		  if((stat_para[12]==1) && (vr->surface_opacity<0.0)) {
			  fprintf(stderr, "surface_opacity parameter cannot be less than 0.0\n");
			  fprintf(stderr, "Please re-input it and run again\n");
			  exit(0);
		  }
			 */
			if(stat_para[13]==0) {
				stat_para[13]=1;
				vr->color_mapping_style=1;
			}
			if((vr->color_mapping_style<1) || (vr->color_mapping_style>4)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1011:color_mapping_style should be between 1 and 4\n");
				HECMW_vis_print_exit("Please re-input it and run again");
			}
			if((vr->color_mapping_style==3) && (stat_para[14]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1012:For color_mapping_style 3, the parameter interval_mapping_num is required\n");
				HECMW_vis_print_exit("Please re-input the value of interval_mapping_num");
			}
			if((stat_para[14]==1) && (vr->interval_mapping_num<=0)) {
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1013:For color_mapping_style 3, the interval_mapping_num should be greater than 0");
			}
			if(((vr->color_mapping_style==2) || (vr->color_mapping_style==3)) && (stat_para[15]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1014:For color_mapping_style =2 or 3, the interval_point should be defined\n");
				HECMW_vis_print_exit("Please re-input the interval_point");
			}
			if(stat_para[16]==0) {
				stat_para[16]=1;
				vr->transfer_function_style=1;
				stat_para[17]=1;
				vr->opa_value=0.05;
			}
			if((vr->transfer_function_style<1) || (vr->transfer_function_style>8)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1028:transfer_function_style should be between 1 and 8\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if((stat_para[17]==1) && (vr->opa_value<0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1029:opacity_value cannot be less than 0.0\n");
				fprintf(stderr, "Please re-input and run again\n");
				exit(0);
			}
			if(((vr->transfer_function_style==3) || (vr->transfer_function_style==4)) && (stat_para[18]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1030:When transfer_function_style =3 or 4, num_of_features should be defined\n");
				HECMW_vis_print_exit("Please re-input the num_of_features");
			}
			if((stat_para[18]==1) && (vr->num_of_features<=0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1031:When transfer_function_style =3 or 4, num_of_features should be greater than 0\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if((stat_para[18]==1) && (stat_para[19]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1032:For transfer_function_style =3 or 4, feature_points should be defined\n");
				HECMW_vis_print_exit("Please re-input your control file");
			}
			if((vr->transfer_function_style==8) && (stat_para[20]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1033:For transfer_function_style=8, lookup_filename should be specified\n");
				HECMW_vis_print_exit("Please re-input the filename");
			}
			if(stat_para[21]==0)
				vr->rotate_style=0;
			if((vr->rotate_style<0) || (vr->rotate_style>4)) {
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1039:The rotate_style parameter is wrong. Please input one within (0, 4)");
			}
			if(stat_para[23]==0)
				vr->color_mapping_bar_on=0;
			if((vr->color_mapping_bar_on<0) || (vr->color_mapping_bar_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1016:color_mapping_bar_on parameter only can be defined as 0 or 1\n");
				HECMW_vis_print_exit("Please re-input it and run again");
			}
			if(stat_para[24]==0)
				vr->scale_marking_on=0;
			if((vr->scale_marking_on<0) || (vr->scale_marking_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1017:scale_marking_on parameter only can be defined as 0 or 1\n");
				HECMW_vis_print_exit("Please re-input it and run again");
			}

			if((vr->color_mapping_bar_on==1) && (vr->xr<40)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1018:x resolution should be larger than 40 for adding color_mapping_bar case\n");
				HECMW_vis_print_exit("please reinput x resolution");
			}
			if((vr->scale_marking_on==1) && (vr->xr<65)) {
				fprintf(stderr, " ERROR: HEC-MW-VIS-E1019: x resolution should be larger than 65 for adding color_mapping and scale marking case\n");
				HECMW_vis_print_exit("please re-input x resolution again");
			}
			if(stat_para[25]==0) {
				sprintf(vr->color_comp_name, "%s","NULL");
			}
			if(stat_para[26]==0) {
				sprintf(vr->color_subcomp_name,"%s", "norm");
			}
			if(stat_para[22]==1) {
				vr->nv_xyz[0]=pesize;
				vr->nv_xyz[1]=1;
				vr->nv_xyz[2]=1;
			}
			else if(stat_para[22]==0) {
				if((stat_para[27]==0) && (stat_para[28]==0) && (stat_para[29]==0)) {
					vr->nv_xyz[0]=pesize;

					vr->nv_xyz[1]=1;
					vr->nv_xyz[2]=1;
				}
				else if((stat_para[27]==0) || (stat_para[28]==0) || (stat_para[29]==0)) {
					HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1040: please define all the n_voxel_x, n_voxel_y, n_voxel_z parameters");
				}
				else if(vr->nv_xyz[0]*vr->nv_xyz[1]*vr->nv_xyz[2]!=pesize) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1041: n_voxel_x*n_voxel_y*n_voxel_z should be equal to the number of pes\n");
					HECMW_vis_print_exit("Please re-input again");
				}
				else if((vr->nv_xyz[0]<=0) || (vr->nv_xyz[1]<=0) || (vr->nv_xyz[2]<=0)) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1042: n_voxel_x,n_voxel_y, and n_voxel_z cannot be less or equal to 0\n");
					HECMW_vis_print_exit("Please re-input again");
				}

			}
			/*		  if((vr->surface_on==1) && (stat_para[30]==0)) {
			  fprintf(stderr, "The surface_filename should be defined in the surface on case\n");
			  fprintf(stderr, "Please input again\n");
			  exit(0);
		  }
			 */
			if((vr->rotate_style>=1) && (vr->rotate_style<=4)) {
				if(stat_para[31]==0)
					vr->num_of_frames=8;
				if(vr->num_of_frames<=0) {
					fprintf(stderr, "ERROR: HEC-MW-VIS-E1015: The parameter num_of_frames cannot be less than 1.\n");
					HECMW_vis_print_exit("Please re-input and run again");
				}
				if(vr->rotate_style==4)
					vr->num_of_frames=8;
			}
			if(stat_para[32]==0)
				vr->background_color[0]=vr->background_color[1]=vr->background_color[2]=0.0;
			if(stat_para[33]==0)
				vr->font_color[0]=vr->font_color[1]=vr->font_color[2]=1.0;
			if(stat_para[34]==0)
				vr->color_system_type=1;
			if((vr->color_system_type<=0) || (vr->color_system_type>3)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1020: color_system_type should be between 1 and 3\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if((vr->background_color[0]<-EPSILON) || (vr->background_color[0]>1.0+EPSILON) ||
					(vr->background_color[1]<-EPSILON) || (vr->background_color[1]>1.0+EPSILON) ||
					(vr->background_color[2]<-EPSILON) || (vr->background_color[2]>1.0+EPSILON)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1021: The background color should be in the interval of  (0.0, 1.0)\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if((vr->font_color[0]<-EPSILON) || (vr->font_color[0]>1.0+EPSILON) ||
					(vr->font_color[1]<-EPSILON) || (vr->font_color[1]>1.0+EPSILON) ||
					(vr->font_color[2]<-EPSILON) || (vr->font_color[2]>1.0+EPSILON)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1022: The font color should be in the interval of  (0.0, 1.0)\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if(stat_para[35]==0)
				vr->font_size=1.0;
			if(vr->font_size<1.0-EPSILON) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1023: font_size paramters cannot be between 1.0 and 4.0\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if(vr->font_size>4.0)
				vr->font_size=4.0;
			if(stat_para[36]==0)
				vr->color_bar_style=2;
			if((vr->color_bar_style<1) || (vr->color_bar_style>2)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1024: color_bar_style only can be 1 or 2\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[37]==0)
				vr->fixed_range_on=0;
			if((vr->fixed_range_on<0) || (vr->fixed_range_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1043: fixed_range_on only can be 0 or 1\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[39]==0)
				vr->num_of_scale=3;
			if(vr->num_of_scale<=0) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1025: num_of_scale only can be greater than 0\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if((stat_para[37]==1) && (stat_para[38]==0)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1026:Please input range_value for fixed_range_on style\n");
				HECMW_vis_print_exit("Please re-input and run again");
			}
			if(stat_para[40]==0)
				vr->mark_0_on=0;
			if((vr->mark_0_on<0) || (vr->mark_0_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1044:mark_0_on only can be 0 or 1\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[41]==0)
				vr->remove_0_display_on=0;
			if((vr->remove_0_display_on<0) || (vr->remove_0_display_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1045:remove_0_display_on only can be 0 or 1\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[42]==0)
				vr->specified_level[1]=0;
			if(vr->specified_level[0]<0) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1046: x_specified_level should be greater than 0\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[43]==0)
				vr->specified_level[1]=0;
			if(vr->specified_level[1]<0) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1047: y_specified_level should be greater than 0\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[44]==0)
				vr->specified_level[2]=0;
			if(vr->specified_level[2]<0) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1048: z_specified_level should be greater than 0\n");
				HECMW_vis_print_exit("Please input and run again");
			}
			if(stat_para[45]==0)
				vr->histogram_on=0;
			if(stat_para[47]==0)
				vr->time_mark_on=0;
			if((vr->time_mark_on<0) || (vr->time_mark_on>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1050: time_mark_on value should be 0 or 1");
				HECMW_vis_print_exit("Please re-input again");
			}
			if(stat_para[48]==0)
				vr->fixed_scale_mark=0;
			if((vr->fixed_scale_mark<0) || (vr->fixed_scale_mark>1)) {
				fprintf(stderr, "ERROR: HEC-MW-VIS-E1049: fixed_scale_mark value should be 0 or 1");
				HECMW_vis_print_exit("Please re-input again");
			}
			t1=pvr->next_pvr;
			t1->vr=vr;

			for(i=0;i<NUM_CONTROL_PVR;i++)
				t1->stat_para[i]=stat_para[i];
			for(i=0;i<NUM_CONTROL_PVR;i++)
				HECMW_free(parameters[i]);
			HECMW_free(parameters);
			HECMW_free(stat_para);
			HECMW_free(len_para);
		} /*end if visual_method=2 */




	} /*end if visual_method>0 */
	/* print some control information */
	return;
}




#if 0
void set_default_vr( Parameter_vr *vr,int stat_para[NUM_CONTROL_PVR],  int pesize)
{
	int i;
	stat_para[0]=stat_para[1]=stat_para[2]=0;
	vr->max_level=64;
	vr->xr=256;
	vr->yr=256;
	stat_para[3]=1;
	vr->num_of_lights=1;
	stat_para[4]=stat_para[5]=stat_para[6]=0;
	stat_para[7]=1;
	vr->up[0]=0.0;
	vr->up[1]=1.0;
	vr->up[2]=0.0;
	stat_para[8]=stat_para[9]=stat_para[10]=1;
	vr->k_ads[0]=0.5;
	vr->k_ads[1]=0.8;
	vr->k_ads[2]=0.6;
	stat_para[11]=1;
	/*	vr->surface_on=0;
	 */	stat_para[12]=0;
	 stat_para[13]=1;
	 vr->color_mapping_style=1;
	 stat_para[14]=stat_para[15]=0;
	 stat_para[16]=stat_para[17]=1;
	 stat_para[46]=0;
	 vr->transfer_function_style=1;
	 vr->opa_value=0.1;
	 for(i=18;i<NUM_CONTROL_PVR;i++)
		 stat_para[i]=0;
	 vr->color_mapping_bar_on=0;
	 vr->scale_marking_on=0;
	 vr->rotate_style=0;
	 vr->color_system_type=1;
	 sprintf(vr->color_comp_name,"%s", "NULL");
	 sprintf(vr->color_subcomp_name,"%s", "norm");
	 vr->nv_xyz[0]=pesize;
	 vr->nv_xyz[1]=1;
	 vr->nv_xyz[2]=1;
	 vr->font_size=1.0;
	 vr->color_bar_style=2;
	 vr->fixed_range_on=0;
	 vr->num_of_scale=3;
	 vr->mark_0_on=0;
	 vr->remove_0_display_on=0;
	 for(i=0;i<3;i++)
		 vr->specified_level[i]=0;
	 vr->histogram_on=0;
	 return;
}


void read_lookup_table(Parameter_vr *vr, double *opa_table)
{
	int i;
	FILE *opa_fp;
	if((opa_fp=fopen(vr->name_lookup,"r"))== NULL) {
		fprintf(stderr, "There is not such a opacity file:\n");
		exit (0);
	}
	for(i=0;i<256;i++)
		fscanf(opa_fp, "%lf", &(opa_table[i]));
	return;
}
#endif
