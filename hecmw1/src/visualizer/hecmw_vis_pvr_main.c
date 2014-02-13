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

#include "hecmw_vis_pvr_main.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_bmp.h"
#include "hecmw_vis_read_voxel.h"
#include "hecmw_vis_ucd_trans.h"
#include "hecmw_vis_mesh_extent.h"
#include "hecmw_vis_new_refine.h"
#include "hecmw_vis_generate_histogram_vr.h"
#include "hecmw_vis_define_parameters.h"
#include "hecmw_vis_subimage_composite_vr.h"
#include "hecmw_vis_font_texture.h"
#include "hecmw_vis_resampling.h"
#include "hecmw_malloc.h"


double *voxel_dxyz, *voxel_orig_xyz;
int *level,*voxel_n_neighbor_pe, **voxel_neighbor_pe;
double *extent;


void HECMW_vis_pvr_rendering(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *timestep, int *init_flag,
		int num_of_pvr, Parameter_vr *vr,
		int stat_para[NUM_CONTROL_PVR], char *outfile,  HECMW_Comm VIS_COMM)
{

	int			mynode, pesize;
	char 		fname[HECMW_FILENAME_LEN];

	FILE			*FP;

	In_surface    *surface;
	int           time_step;
	double        t1, t2, t3;
	int         n_voxel;
	/*
  int max_level, xr, yr, projection_style, num_of_lights, color_mapping_style, interval_mapping_num,
	  num_of_features, rotate_style, transfer_function_style, color_mapping_bar_on, scale_marking_on, num_of_frames,
	  color_system_type,color_bar_style, fixed_range_on, num_of_scale, mark_0_on, specified_level[3],
	  nv_xyz[3], remove_0_display_on, histogram_on;
char name_lookup[128], name_voxelfile[128],color_comp_name[128], color_subcomp_name[5];

double *light_point, view_point_d[3], screen_point[3], up[3], k_ads[3], *interval_point, opa_value, *fea_point,
        background_color[3], font_color[3], font_size, range_value[2], display_range[6];
	 */
	double view_n[3], nview_n;
	double vertex[3*8];
	double point_s[3], point_o[3], ray_direction[3];
	int i, j, k, m;
	double range[6], xd, yd;
	double coff_matrix[3][3], n_vertex[3*8], inv_matrix[3][3];
	double mincolor, maxcolor, view_point[3], scr_area[4], x, y, nn;
	double tmincolor, tmaxcolor, trange[6], svertex[3*8], sn_vertex[3*8], sscr_area[4];
	int intersection;
	double accum_rgba[4];
	double grad_minmax[2], dis_minmax[2], feap_minmax[2], feai_minmax[2];
	double tfeap_minmax[2];
	double *opa_table;
	extern void HECMW_vis_memory_exit(char *var);
	extern void HECMW_vis_print_exit(char *var);
	double *accum_opa;
	double *image;
	double center[3], dis_c_v, *ndis_c_v, tmpd;
	int *pe_id, tmpi, pix_num, pix0_num, pixn;
	double *n_subimage, *n_subopa, *subimage, *subopa;
	double *subimage1, *subopa1;
	double first_p[3];
	HECMW_Status	stat;
	int starti, startj, endi, endj;
	int start_x, start_y;
	double minx, miny, minz, maxx, maxy, maxz, tminx, tminy, tminz, tmaxx, tmaxy, tmaxz;
	double av_length, tav_length;
	int first_ijk[3];
	int  num_img, ii, base_x, base_y;
	double o_v_point[3],  m_scr[4];
	int  trace_flag;
	BITMAPFILEHEADER header;       /* File header */

	unsigned char r, g, b;
	int  ri, gi, bi;
	BITMAPINFOHEADER info;
	int scale, scale_type;

	char timestep_str[128], rotate_str[128];

	/* --------new variables for vectorization --------*/

	double *node1;
	/*
  double	*voxel_dxyz;
  double	*voxel_orig_xyz;
  int		*level;
  int		*voxel_n_neighbor_pe;
  int		**voxel_neighbor_pe;
	 */
	int       x_specified_l,y_specified_l,z_specified_l;

	int   nx, ny, nz, *empty_flag;
	double *var, *grad_var;
	int   r_level[3];
	double r_dxyz[3], dxyz[3], orig_xyz[3];
	double org_mincolor, org_maxcolor;



	HECMW_Comm_size(VIS_COMM, &pesize);
	HECMW_Comm_rank(VIS_COMM, &mynode);
	n_voxel=pesize;
	/*  if(*init_flag==1) {
      contfile1[0] = '\0';
	  j=0;
      while (contfile[j] == ' ')
	    j++;
      k = 0;
      while (contfile[j] != ' ') {
    	contfile1[k] = contfile[j];
    	k++;
    	j++;
      }
      contfile1[k] = '\0';
      if((contfp=fopen(contfile1,"r"))== NULL) {
         fprintf(stderr, "There is not such a control file:%s\n", contfile1);
          exit (0);
	  }
      vr=(Parameter_vr *)HECMW_malloc(sizeof(Parameter_vr));
	  if(vr==NULL)
		  HECMW_vis_memory_exit("vr");


  if(strncmp(contfile1, "NULL", 4)==0)
	  set_default_vr(vr, stat_para, pesize);
  else {
   if((contfp=fopen(contfile1,"r"))== NULL)
	   HECMW_vis_print_exit("There is not such a control file:\n");
    read_control_file(contfp, vr, stat_para, pesize);
    fprintf(stderr, "max_level=%d xr=%d yr=%d num_of_light=%d\n", vr->max_level, vr->xr, vr->yr, vr->num_of_lights);
   fclose(contfp);
  }
	 */

	if((*init_flag==1) || (num_of_pvr>1)) {
		voxel_dxyz=(double *)HECMW_calloc(pesize*3, sizeof(double));
		voxel_orig_xyz=(double *)HECMW_calloc(pesize*3, sizeof(double));
		level=(int *) HECMW_calloc(pesize*3, sizeof(int));
		voxel_n_neighbor_pe=(int *)HECMW_calloc(pesize, sizeof(int));
		voxel_neighbor_pe=(int **)HECMW_calloc(pesize, sizeof(int *));
		if((voxel_dxyz==NULL) || (voxel_orig_xyz==NULL) || (level==NULL) || (voxel_n_neighbor_pe==NULL) || (voxel_neighbor_pe==NULL))
			HECMW_vis_memory_exit("voxel information");

		for(j=0;j<pesize;j++) {
			voxel_neighbor_pe[j]=(int *)HECMW_calloc(pesize, sizeof(int));
			if(voxel_neighbor_pe[j]==NULL)
				HECMW_vis_memory_exit("voxel_neighbor_pe[j]");
			for(k=0;k<pesize;k++)
				voxel_neighbor_pe[j][k]=0;
		}
		n_voxel=pesize;

		if(stat_para[22]==1) {
			/*	  fprintf(stderr, "voxel file is %s\n", vr->name_voxelfile);
			 */
			read_voxel_file(vr->name_voxelfile, pesize, voxel_dxyz, voxel_orig_xyz, level, voxel_n_neighbor_pe, voxel_neighbor_pe);
		}
	} /* end of *init_flag=1 */


	node1 = (double *)HECMW_calloc(mesh->n_node, sizeof(double));

	transform_ucd_pvr(data, node1, mesh, vr,mynode, pesize,
			VIS_COMM,voxel_dxyz, voxel_orig_xyz, level, voxel_n_neighbor_pe, voxel_neighbor_pe, stat_para[22], stat_para[46],
			*init_flag, num_of_pvr);
	surface=NULL;
	/*  if(pvr->surface_on==1) {
  surface=(In_surface *)HECMW_malloc(sizeof(In_surface));
  if(surface==NULL) {
	  fprintf(stderr, "There is no enough memory for surface\n");
	  exit(0);
  }
  }
	 */
	if(*init_flag==1) {
		extent = (double *)HECMW_calloc(mesh->n_elem*6, sizeof(double));
		if(extent==NULL)
			HECMW_vis_memory_exit("extent");
		calc_extent(mesh, extent);

		calc_voxel_level(n_voxel, mesh, voxel_dxyz, voxel_orig_xyz, extent, level, VIS_COMM);
		if(mynode==0)
			for(i=0;i<n_voxel;i++)
				fprintf(stderr, "voxel %d level: %d %d %d \n", i, level[i*3], level[i*3+1], level[i*3+2]);
		if(vr->specified_level[0]>0) {
			x_specified_l=vr->specified_level[0];
			for(i=0;i<n_voxel;i++)
				level[i*3]=x_specified_l;
		}
		else if(vr->max_level>0) {
			for(i=0;i<n_voxel;i++) {
				if(level[i*3]>vr->max_level)
					level[i*3]=vr->max_level;
			}
		}
		if(vr->specified_level[1]>0) {
			y_specified_l=vr->specified_level[1];
			for(i=0;i<n_voxel;i++)
				level[i*3+1]=y_specified_l;
		}
		else if(vr->max_level>0) {
			for(i=0;i<n_voxel;i++) {
				if(level[i*3+1]>vr->max_level)
					level[i*3+1]=vr->max_level;
			}
		}
		if(vr->specified_level[2]>0) {
			z_specified_l=vr->specified_level[2];
			for(i=0;i<n_voxel;i++)
				level[i*3+2]=z_specified_l;
		}
		else if(vr->max_level>0) {
			for(i=0;i<n_voxel;i++) {
				if(level[i*3+2]>vr->max_level)
					level[i*3+2]=vr->max_level;
			}
		}
	}

#ifdef slow
	phi_x = (double *)HECMW_calloc(n_node, sizeof(double));
	phi_y = (double *)HECMW_calloc(n_node, sizeof(double));
	phi_z = (double *)HECMW_calloc(n_node, sizeof(double));
	for(i=0;i<n_node;i++) {
		phi_x[i]=0.0;
		phi_y[i]=0.0;
		phi_z[i]=0.0;
	}

#endif
	HECMW_Barrier(VIS_COMM);
	/*      t2=HECMW_Wtime();
    if (mynode == 0){
      fprintf(stderr, "finish reading files\n");
      fprintf(stderr, "the time for reading file is %lf the total time up to now is %lf\n", t2-t1, t2-t1);
    }
	 */
	t1=HECMW_Wtime();
	if (mynode == 0){
		fprintf(stderr, "Start parallel volume rendering module\n");
	}
	t2=t1;

	if(mynode==0)
		fprintf(stderr, "calc_gradient\n");
#ifdef slow
	calc_gradient(n_node, n_internal, n_elem, node, elem, node1,
			n_neighbor_pe, neighbor_pe,
			import_index, import_node,
			export_index, export_node,  phi_x,  phi_y, phi_z,
			mynode, VIS_COMM);
	t3=HECMW_Wtime();
	if(mynode==0) {
		fprintf(stderr, "The time for computing gradient is %lf the total time up to now is %lf\n", t3-t2, t3-t1);
		t2=t3;
	}
	HECMW_Barrier(VIS_COMM);
	fprintf(stderr, "refinement\n");
#endif
	nx=level[mynode*3];
	ny=level[mynode*3+1];
	nz=level[mynode*3+2];
	empty_flag=(int *)HECMW_calloc((nx+1)*(ny+1)*(nz+1), sizeof(int));
	var=(double *)HECMW_calloc((nx+1)*(ny+1)*(nz+1), sizeof(double));
#ifdef slow
	grad_var=(double *)HECMW_calloc((nx+1)*(ny+1)*(nz+1)*3, sizeof(double));
	if((grad_var==NULL) || (var==NULL) || (empty_flag==NULL)) {
		fprintf(stderr, "There is no enough memory: empty_flag, var, gradient\n");
		exit(0);
	}

#endif
	refinement(mesh, node1, n_voxel, voxel_dxyz, voxel_orig_xyz, level,
			voxel_n_neighbor_pe, voxel_neighbor_pe, extent,  mynode, VIS_COMM, empty_flag, var);
#ifdef slow
	HECMW_free(phi_x);
	HECMW_free(phi_y);
	HECMW_free(phi_z);
#endif
	HECMW_Barrier(VIS_COMM);

	t3=HECMW_Wtime();
	if(mynode==0) {
		fprintf(stderr, " Finish refinement now, refinement level: \n");
		fprintf(stderr, " The time for refinement is %lf the total time up to now is %lf\n", t3-t2, t3-t1);
		t2=t3;
	}
	HECMW_free(node1);
	mincolor=tmincolor=1.0E17;
	maxcolor=tmaxcolor=-1.0E17;
	find_color_minmax_vr(var, empty_flag, nx, ny, nz, &mincolor, &maxcolor);
	HECMW_Barrier(VIS_COMM);
	if(pesize>1) {
		HECMW_Allreduce(&mincolor, &tmincolor, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&maxcolor, &tmaxcolor, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	}
	else {
		tmincolor=mincolor; tmaxcolor=maxcolor;
	}
	org_mincolor=tmincolor; org_maxcolor=tmaxcolor;
	if(vr->fixed_range_on==1) {
		tmincolor=vr->range_value[0];
		tmaxcolor=vr->range_value[1];
	}

	/*   fprintf(stderr, "mincolor= %lf   maxcolor= %lf\n", tmincolor, tmaxcolor);
	 */
	if(vr->histogram_on==1)
		output_histogram_vr(tmincolor, tmaxcolor, var, empty_flag, nx, ny, nz, mynode, pesize, VIS_COMM);
	if(vr->histogram_on==2)
		generate_histogram_graph_vr(tmincolor, tmaxcolor, var, empty_flag, nx, ny, nz, mynode, pesize, VIS_COMM, vr->color_system_type);
	if(vr->color_mapping_style==4)
		generate_interval_point_vr(tmincolor, tmaxcolor, var, empty_flag, nx, ny, nz, mynode, pesize, VIS_COMM, vr->interval_point);

	find_minmax_vr(voxel_dxyz, voxel_orig_xyz, mynode, range);
	for(i=0;i<6;i++)
		trange[i]=0.0;
	minx=range[0]; maxx=range[1]; miny=range[2]; maxy=range[3]; minz=range[4]; maxz=range[5];
	HECMW_Barrier(VIS_COMM);
	if(pesize>1) {
		HECMW_Allreduce(&minx, &tminx, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);

		HECMW_Allreduce(&maxx, &tmaxx, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
		HECMW_Allreduce(&miny, &tminy, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&maxy, &tmaxy, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
		HECMW_Allreduce(&minz, &tminz, 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&maxz, &tmaxz, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	}
	else {
		tminx=minx;
		tmaxx=maxx;
		tminy=miny;
		tmaxy=maxy;
		tminz=minz;
		tmaxz=maxz;
	}
	if(mynode==0)
		fprintf(stderr, "the range of the field is %lf %lf %lf %lf %lf %lf\n", tminx, tmaxx, tminy, tmaxy, tminz, tmaxz);
	trange[0]=tminx; trange[1]=tmaxx; trange[2]=tminy; trange[3]=tmaxy;
	trange[4]=tminz; trange[5]=tmaxz;
	if((*init_flag==1) || (num_of_pvr>1)) {
		if(stat_para[4]==0) {
			vr->light_point=(double *)HECMW_calloc(3, sizeof(double));
			if(vr->light_point==NULL)
				HECMW_vis_memory_exit("light_point");

			vr->light_point[0]=(tminx+tmaxx)/2.0;
			vr->light_point[1]=tmaxy+0.1*(tmaxy-tminy);
			vr->light_point[2]=tmaxz+(tmaxz-tminz)*2.0;

		}
		if(stat_para[5]==0) {
			vr->view_point_d[0]=(tminx+tmaxx)/2.0;
			vr->view_point_d[1]=tmaxy+1.5*(tmaxy-tminy);
			vr->view_point_d[2]=tmaxz+1.5*(tmaxz-tminz);
		}
		if(stat_para[6]==0) {
			vr->screen_point[0]=(tminx+tmaxx)/2.0;
			vr->screen_point[1]=(tminy+tmaxy)/2.0;
			vr->screen_point[2]=(tminz+tmaxz)/2.0;
		}

	}
	HECMW_Barrier(VIS_COMM);

	/*----------------------------------------
   Define the sampling and projection parameters
	 */
	/* First, define the equation of screen */
	if(vr->rotate_style==0)
		num_img=1;
	else if((vr->rotate_style>=1) && (vr->rotate_style<=4)) {
		num_img=vr->num_of_frames;
		for(i=0;i<3;i++)
			vr->screen_point[i]=(trange[i*2]+trange[i*2+1])/2.0;
		/*find the range of projection */
		for(i=0;i<3;i++)
			o_v_point[i]=vr->view_point_d[i];
		m_scr[0]=m_scr[2]=1.0E+17;
		m_scr[1]=m_scr[3]=-1.0E+17;

		for(ii=0;ii<num_img;ii++) {
			if(ii!=0)
				view1_parameter_define(ii, vr->num_of_frames, vr->rotate_style, vr->view_point_d, vr->screen_point, vr->num_of_lights,
						vr->light_point, vr->up, trange);
			for(i=0;i<3;i++)
				view_n[i]=vr->screen_point[i]-vr->view_point_d[i];
			nview_n=sqrt(SQR(view_n[0])+SQR(view_n[1])+SQR(view_n[2]));
			if(fabs(nview_n)>EPSILON)
				for(i=0;i<3;i++)
					view_n[i]/=nview_n;
			transform_range_vertex(trange, vertex);

			get_frame_transform_matrix(vr->view_point_d, vr->screen_point, vr->up, coff_matrix);
			transform_frame(vr->screen_point, vertex, coff_matrix,  n_vertex);
			for(i=0;i<3;i++)
				view_point[i]=vr->view_point_d[i]-vr->screen_point[i];
			transform2_frame(coff_matrix, view_point);
			find_projection_range(view_point, n_vertex, scr_area);
			if(scr_area[0]<m_scr[0])
				m_scr[0]=scr_area[0];
			if(scr_area[1]>m_scr[1])
				m_scr[1]=scr_area[1];
			if(scr_area[2]<m_scr[2])
				m_scr[2]=scr_area[2];
			if(scr_area[3]>m_scr[3])
				m_scr[3]=scr_area[3];
		}

		for(i=0;i<3;i++)
			vr->view_point_d[i]=o_v_point[i];
	}


	for(ii=0; ii<num_img;ii++) {
		if(ii!=0)
			view_parameter_define(ii, vr->num_of_frames, vr->rotate_style, vr->view_point_d, vr->screen_point, vr->up, vr->num_of_lights,
					vr->light_point,  trange);
		for(i=0;i<3;i++)
			view_n[i]=vr->screen_point[i]-vr->view_point_d[i];
		nview_n=sqrt(SQR(view_n[0])+SQR(view_n[1])+SQR(view_n[2]));
		if(fabs(nview_n)>EPSILON)
			for(i=0;i<3;i++)
				view_n[i]/=nview_n;
		/*   for(i=0;i<3;i++)
	   p_screen[i]=view_n[i];
   p_screen[3]=-view_n[0]*vr->screen_point[0]-view_n[1]*vr->screen_point[1]-
	   view_n[2]*vr->screen_point[2];
		 */

		/* Second, find the projection of the dataset in each PE */
		transform_range_vertex(trange, vertex);
		if(mynode==0)
			fprintf(stderr, "viewpoint=%lf %lf %lf screen_point=%lf %lf %lfup=%lf %lf %lf\n", vr->view_point_d[0], vr->view_point_d[1],
					vr->view_point_d[2], vr->screen_point[0], vr->screen_point[1], vr->screen_point[2], vr->up[0], vr->up[1], vr->up[2]);
		get_frame_transform_matrix(vr->view_point_d, vr->screen_point, vr->up, coff_matrix);
		transform_frame(vr->screen_point, vertex, coff_matrix,  n_vertex);
		/*  if(mynode==0){
   for(i=0;i<3;i++)
    for(j=0;j<3;j++)
  }
		 */
		for(i=0;i<3;i++)
			view_point[i]=vr->view_point_d[i]-vr->screen_point[i];
		transform2_frame(coff_matrix, view_point);
		find_projection_range(view_point, n_vertex, scr_area);
		/* for multi-PEs, should find minx, maxx, miny, maxy in scr_area for all PE
		 */
		if(ii==0) {
			if((vr->color_mapping_bar_on==0) && (vr->scale_marking_on==0)) {
				start_x=10; start_y=10;
				/*  if((vr->xr<=20) || (vr->yr<=20)) {
       fprintf(stderr, "The resolution cannot be less than 20\n");
       exit(0);
   }
				 */
				if(num_img==1) {
					xd=(scr_area[1]-scr_area[0])/(vr->xr-20);
					yd=(scr_area[3]-scr_area[2])/(vr->yr-20);
				}
				if(num_img>1) {
					xd=(m_scr[1]-m_scr[0])/(vr->xr-20);
					yd=(m_scr[3]-m_scr[2])/(vr->yr-20);
				}
				if(xd>=yd) {
					yd=xd;
					if(num_img==1)
						vr->yr=(int)((int)(20+(scr_area[3]-scr_area[2])/yd)/10)*10+10;
					else if(num_img>1)
						vr->yr=(int)((int)(20+(m_scr[3]-m_scr[2])/yd)/10)*10+10;
				}
				if(xd<yd) {
					xd=yd;
					if(num_img==1)
						vr->xr=(int)((int)(20+(scr_area[1]-scr_area[0])/xd)/8)*8+8;
					else if(num_img>1)
						vr->xr=(int)((int)(20+(m_scr[1]-m_scr[0])/xd)/8)*8+8;
				}
			}
			else if((vr->color_mapping_bar_on==1) && (vr->scale_marking_on==0)) {
				start_x=30; start_y=10;
				if(num_img==1) {
					xd=(scr_area[1]-scr_area[0])/(vr->xr-40);
					yd=(scr_area[3]-scr_area[2])/(vr->yr-20);
				}
				else if(num_img>1) {
					xd=(m_scr[1]-m_scr[0])/(vr->xr-40);
					yd=(m_scr[3]-m_scr[2])/(vr->yr-20);
				}

				if(xd>=yd) {
					yd=xd;
					if(num_img==1)
						vr->yr=(int)((int)(20+(scr_area[3]-scr_area[2])/yd)/10)*10+10;
					else if(num_img>1)
						vr->yr=(int)((int)(20+(m_scr[3]-m_scr[2])/yd)/10)*10+10;
				}
				else if(xd<yd) {
					xd=yd;
					if(num_img==1)
						vr->xr=(int)((int)(40+(scr_area[1]-scr_area[0])/xd)/8)*8+8;
					else if(num_img>1)
						vr->xr=(int)((int)(40+(m_scr[1]-m_scr[0])/xd)/8)*8+8;

				}
			}
			else if((vr->color_mapping_bar_on==1) && (vr->scale_marking_on==1)) {
				scale=(int)vr->font_size;
				if((vr->font_size-scale)<0.5-EPSILON)
					scale_type=1;
				else
					scale_type=2;
				if(scale_type==1) {
					start_x=10+45*scale;
					if(vr->color_bar_style==2)
						start_x=10+45*scale+15;
				}
				else if(scale_type==2) {
					start_x=10+63*scale;
					if(vr->color_bar_style==2)
						start_x=10+63*scale+15;
				}

				start_y=10;
				if(vr->xr<start_x+10) {
					fprintf(stderr, "The image x_resolution cannot write such size charaters\n");
					fprintf(stderr, "Please reduce the font size or enlarge x_resolution and run again\n");
					exit(0);
				}
				if(num_img==1) {
					xd=(scr_area[1]-scr_area[0])/(vr->xr-(start_x+10));
					yd=(scr_area[3]-scr_area[2])/(vr->yr-20);
				}
				else if(num_img>1) {
					xd=(m_scr[1]-m_scr[0])/(vr->xr-(start_x+10));
					yd=(m_scr[3]-m_scr[2])/(vr->yr-20);
				}

				if(xd>=yd) {
					yd=xd;
					if(num_img==1)
						vr->yr=(int)((int)(20+(scr_area[3]-scr_area[2])/yd)/10)*10+10;
					else if(num_img>1)
						vr->yr=(int)((int)(20+(m_scr[3]-m_scr[2])/yd)/10)*10+10;
				}
				else if(xd<yd) {
					xd=yd;
					if(num_img==1)
						vr->xr=(int)((int)(start_x+10+(scr_area[1]-scr_area[0])/xd)/8)*8+8;
					else if(num_img>1)
						vr->xr=(int)((int)(start_x+10+(m_scr[1]-m_scr[0])/xd)/8)*8+8;

				}
			}


		}
		if(ii>0) {
			start_x-=base_x;
			start_y-=base_y;
		}
		if(num_img>1) {
			base_x=(int)((scr_area[0]-m_scr[0])/xd);
			base_y=(int)((scr_area[2]-m_scr[2])/yd);

			start_x+=base_x;
			start_y+=base_y;
		}
		/* find the subimage range for this PE */
		transform_range_vertex(range, svertex);
		/*   output_frame(vr, view_point, n_vertex, scr_area, xd, yd);
		 */
		transform_frame(vr->screen_point, svertex, coff_matrix, sn_vertex);
		find_projection_range(view_point, sn_vertex, sscr_area);
		/*   starti=(int)((sscr_area[0]-scr_area[0]+EPSILON)/xd);
   endi=(int)((sscr_area[1]-scr_area[0]-EPSILON)/xd)+1;
   startj=(int)((sscr_area[2]-scr_area[2]+EPSILON)/yd);
   endj=(int)((sscr_area[3]-scr_area[2]-EPSILON)/yd)+1;
   if(starti<0) starti=0;
   if(startj<0) startj=0;
   if(endi>vr->xr) endi=vr->xr;
   if(endj>vr->yr) endj=vr->yr;
		 */
		starti=(int)((sscr_area[0]-scr_area[0])/xd)+start_x-1;
		endi=(int)((sscr_area[1]-scr_area[0])/xd)+start_x+1;
		startj=(int)((sscr_area[2]-scr_area[2])/yd)+start_y-1;
		endj=(int)((sscr_area[3]-scr_area[2])/yd)+start_y+1;

		HECMW_Barrier(VIS_COMM);
		grad_minmax[0]=1.0E17;
		grad_minmax[1]=-1.0E17;
		if(vr->transfer_function_style==2) {
			/*	   find_grad_minmax(vr, grad_var, nx, ny, nz, grad_minmax);
	   HECMW_Barrier(VIS_COMM);
   HECMW_Allreduce(&grad_minmax[0], &tgrad_minmax[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
   HECMW_Allreduce(&grad_minmax[1], &tgrad_minmax[1], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
   grad_minmax[0]=tgrad_minmax[0]; grad_minmax[1]=tgrad_minmax[1];
			 */   }
		feap_minmax[0]=1.0E17;
		feap_minmax[1]=-1.0E17;
		if(vr->transfer_function_style==3) {
			find_feap_minmax(vr->num_of_features, vr->fea_point, tmincolor, tmaxcolor, feap_minmax);
			HECMW_Barrier(VIS_COMM);
			HECMW_Allreduce(&feap_minmax[0], &tfeap_minmax[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
			HECMW_Allreduce(&feap_minmax[1], &tfeap_minmax[1], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
			feap_minmax[0]=tfeap_minmax[0];  feap_minmax[1]=tfeap_minmax[1];
		}
		feai_minmax[0]=1.0E17;
		feai_minmax[1]=-1.0E17;
		if(vr->transfer_function_style==4) {
			/*	   find_feai_minmax(num_of_features, fea_point, tmincolor, tmaxcolor, feai_minmax);
	   HECMW_Barrier(VIS_COMM);
       HECMW_Allreduce(&feai_minmax[0], &tfeai_minmax[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
       HECMW_Allreduce(&feai_minmax[1], &tfeai_minmax[1], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	   feai_minmax[0]=tfeai_minmax[0];  feai_minmax[1]=tfeai_minmax[1];
			 */
		}
		dis_minmax[0]=1.0E17;
		dis_minmax[1]=-1.0E17;
		if((vr->transfer_function_style==5) || (vr->transfer_function_style==6))
			find_dis_minmax(vr->view_point_d, vertex, dis_minmax);
		/*       HECMW_Allreduce(&dis_minmax[0], &dis_minmax[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
	   if(time_step==0) {

   if(transfer_function_style==7) {
	   opa_table=(double *)HECMW_calloc(256, sizeof(double));
	   if(opa_table==NULL) {
		   fprintf(stderr, " There is no enough memory: opa_table\n");
		   exit(0);
	   }
	   if(mynode==0)
		   read_lookup_table(name_lookup, opa_table);
       HECMW_Bcast(opa_table, 256, HECMW_DOUBLE, 0, VIS_COMM);
   }

} */
		/*if(time_step==0) {
      num_refine=1;
   for(i=0;i<vd->leveltot;i++)
	   num_refine*=2;
		 */
		av_length=sqrt(SQR(voxel_dxyz[mynode*3]/nx)+SQR(voxel_dxyz[mynode*3+1]/ny)
				+SQR(voxel_dxyz[mynode*3+2]/nz));

		HECMW_Barrier(VIS_COMM);
		if(pesize>1)
			HECMW_Allreduce(&av_length, &tav_length, 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
		else
			tav_length=av_length;

		/* Build hierachical tree*/
		image=(double *)HECMW_calloc(vr->xr*vr->yr*3, sizeof(double));
		if(image==NULL)
			HECMW_vis_memory_exit("image");
		accum_opa=(double *)HECMW_calloc(vr->xr*vr->yr, sizeof(double));
		if(accum_opa==NULL)
			HECMW_vis_memory_exit("accum_opa");
		/*}
		 */
		/*---------new change for vec ---------*/
		for(j=0;j<vr->yr*vr->xr*3;j++)
			image[j]=0.0;
		for(j=0;j<vr->yr*vr->xr;j++)
			accum_opa[j]=0.0;

		HECMW_Barrier(VIS_COMM);
		t3=HECMW_Wtime();
		if(mynode==0) {
			fprintf(stderr, "The time for preprocessing of pvr is %lf, the total time up to now is %lf\n", t3-t2, t3-t1);
			t2=t3;
		}
		find_inverse_matrix(coff_matrix, inv_matrix);
		for(i=0;i<3;i++) {
			r_level[i]=level[mynode*3+i];
			dxyz[i]=voxel_dxyz[mynode*3+i];
			r_dxyz[i]=dxyz[i]/r_level[i];
			orig_xyz[i]=voxel_orig_xyz[mynode*3+i];
		}
		/*   fprintf(stderr, "projection range is  %d %d %d %d", starti, endi, startj, endj);
		 */
		mincolor=tmincolor;
		maxcolor=tmaxcolor;
		if(vr->color_mapping_style==2) {
			mincolor=vr->interval_point[0];
			maxcolor=vr->interval_point[1];
		}
		if(vr->color_mapping_style==3) {
			mincolor=vr->interval_point[0];
			maxcolor=vr->interval_point[2*vr->interval_mapping_num];
		}

		for(j=startj;j<endj;j++)
			for(i=starti;i<endi;i++) {
				/*   for(j=0;j<yr;j++)
	   for(i=0;i<vr;i++) {

				 */
				/* for Multi-PEs, first must judge whether it is within the PE projection range

				 */
				x=xd*(i-start_x+0.5)+scr_area[0];
				y=yd*(j-start_y+0.5)+scr_area[2];
				/*transform to the original frame */
				point_s[0]=x; point_s[1]=y; point_s[2]=0.0;
				tranverse_transform(vr->screen_point, point_s, inv_matrix, point_o);
				/*find the first intersected voxel*/
				for(m=0;m<3;m++)
					ray_direction[m]=point_o[m]-vr->view_point_d[m];
				/* new---here----to remove the overlap pixel computation */
				trace_flag=1;
				if(fabs(ray_direction[0])<EPSILON) {
					if(fabs(orig_xyz[0]+dxyz[0]-vr->view_point_d[0])<EPSILON) {
						trace_flag=0;
					}
				}
				else if(fabs(ray_direction[1])<EPSILON) {
					if(fabs(orig_xyz[1]+dxyz[1]-vr->view_point_d[1])<EPSILON) {
						trace_flag=0;
					}
				}
				else if(fabs(ray_direction[2])<EPSILON) {
					if(fabs(orig_xyz[2]+dxyz[2]-vr->view_point_d[2])<EPSILON) {
						trace_flag=0;
					}
				}
				if(trace_flag==1) {
					nn=sqrt(SQR(ray_direction[0])+SQR(ray_direction[1])+SQR(ray_direction[2]));
					if(fabs(nn)>EPSILON) {
						for(m=0;m<3;m++)
							ray_direction[m]/=nn;
					}
					/*
	   intersection=find_first_inter(point_o, vr, vd, root_tree, &voxel_p, ray_direction, first_p);
	   if(intersection==1) {
		   int connect[8][6];
		   ray_trace(vr, vd, root_tree, voxel_p, first_p, ray_direction, tmincolor, tmaxcolor, connect,
		   local_face, accum_rgba, grad_minmax, feap_minmax, feai_minmax, dis_minmax,
		   opa_table, tav_length,time_step);
					 */
					intersection=find_first_inter(point_o, vr->view_point_d, r_level, orig_xyz, dxyz, r_dxyz, ray_direction,
							first_p, first_ijk);
					if(intersection==1) {
						ray_trace(vr->remove_0_display_on, vr->color_mapping_style, vr->interval_point, vr->transfer_function_style,
								vr->opa_value, vr->num_of_features, vr->fea_point,
								vr->view_point_d,  vr->interval_mapping_num, vr->color_system_type, vr->num_of_lights,
								vr->light_point, vr->k_ads, orig_xyz, dxyz, r_dxyz, r_level, empty_flag, var, grad_var, first_p,
								first_ijk, ray_direction, mincolor, maxcolor, accum_rgba, grad_minmax, feap_minmax,
								feai_minmax, dis_minmax, opa_table, tav_length, time_step, i,j);

						for(m=0;m<3;m++)
							image[(j*vr->xr+i)*3+m]=accum_rgba[m];
						accum_opa[j*vr->xr+i]=accum_rgba[3];
					}
				}

			}

		HECMW_Barrier(VIS_COMM);
		if(mynode==0) {
			fprintf(stderr, "Finish the computing on each PE\n");
			t3=HECMW_Wtime();
			fprintf(stderr, "The  volume rendering on each PE is %lf the total time up to now is %lf\n", t3-t2, t3-t1);
			t2=t3;
		}


		if(pesize>1) {
			pix_num=(int)(vr->xr*vr->yr/pesize);
			pix0_num=vr->xr*vr->yr-pix_num*(pesize-1);
			if(mynode==0) pixn=pix0_num;
			if(mynode!=0) pixn=pix_num;
			n_subimage=(double *)HECMW_calloc(pesize*pixn*3, sizeof(double));
			n_subopa=(double *)HECMW_calloc(pesize*pixn, sizeof(double));
			if((n_subimage==NULL) || (n_subopa==NULL))
				HECMW_vis_memory_exit("n_subimage, n_subopa");


			for(j=0;j<pixn*3;j++) {
				if(mynode==0)
					n_subimage[j]=image[j];
				else if(mynode!=0)
					n_subimage[mynode*pixn*3+j]=image[pix0_num*3+(mynode-1)*pix_num*3+j];
			}
			if(mynode==0)
				for(j=0;j<pixn;j++)
					n_subopa[j]=accum_opa[j];
			else if(mynode!=0)
				for(j=0;j<pixn;j++)
					n_subopa[mynode*pixn+j]=accum_opa[pix0_num+(mynode-1)*pix_num+j];

			subimage=(double *)HECMW_calloc(pixn*3, sizeof(double));
			subopa=(double *)HECMW_calloc(pixn, sizeof(double));
			if((subimage==NULL) || (subopa==NULL))
				HECMW_vis_memory_exit("subimage, subopa");
			HECMW_Barrier(VIS_COMM);

			for(i=0;i<pesize;i++) {
				if(mynode==i) {
					/*			fprintf(stderr, "in PE=%d \n", mynode);
					 */
					for(k=0;k<pesize;k++) {
						if(k!=i) {
							if(k==0) {
								subimage1=(double *)HECMW_calloc(pix0_num*3, sizeof(double));
								subopa1=(double *)HECMW_calloc(pix0_num, sizeof(double));
							}
							else if(k!=0) {
								subimage1=(double *)HECMW_calloc(pix_num*3, sizeof(double));
								subopa1=(double *)HECMW_calloc(pix_num, sizeof(double));
							}

							if((subimage1==NULL) || (subopa1==NULL))
								HECMW_vis_memory_exit("subimage1, subopa1");
							if(k==0) {
								for(j=0;j<pix0_num*3;j++)
									subimage1[j]=image[j];
								for(j=0;j<pix0_num;j++)
									subopa1[j]=accum_opa[j];
							}
							else if(k!=0) {
								for(j=0;j<pix_num*3;j++)
									subimage1[j]=image[pix0_num*3+(k-1)*pix_num*3+j];
								for(j=0;j<pix_num;j++)
									subopa1[j]=accum_opa[pix0_num+(k-1)*pix_num+j];
							}
							if(k==0) {
								HECMW_Send(subimage1, pix0_num*3, HECMW_DOUBLE, k, 0, VIS_COMM);
								HECMW_Send(subopa1, pix0_num, HECMW_DOUBLE, k, 0, VIS_COMM);
							}
							else if(k!=0) {
								HECMW_Send(subimage1, pix_num*3, HECMW_DOUBLE, k, 0, VIS_COMM);
								HECMW_Send(subopa1, pix_num, HECMW_DOUBLE, k, 0, VIS_COMM);
							}

							HECMW_free(subimage1);
							HECMW_free(subopa1);
						} /*if k!=i */
					} /*loop k*/


				}
				if(mynode!=i) {
					HECMW_Recv(subimage, pixn*3, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);
					HECMW_Recv(subopa, pixn, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);
					for(j=0;j<pixn*3;j++)
						n_subimage[i*pixn*3+j]=subimage[j];
					for(j=0;j<pixn;j++)
						n_subopa[i*pixn+j]=subopa[j];
				}
				HECMW_Barrier(VIS_COMM);

			}

			/*	if(time_step==0) {
			 */
			for(i=0;i<3;i++)
				center[i]=orig_xyz[i]+dxyz[i]/2.0;
			dis_c_v=sqrt(SQR(center[0]-vr->view_point_d[0])+SQR(center[1]-vr->view_point_d[1])
					+SQR(center[2]-vr->view_point_d[2]));
			if(mynode!=0) {
				HECMW_Send(&dis_c_v, 1, HECMW_DOUBLE, MASTER_PE, 0, VIS_COMM);
				pe_id=(int *)HECMW_calloc(pesize, sizeof(int));
				if(pe_id==NULL)
					HECMW_vis_memory_exit("pe_id");
				HECMW_Recv(pe_id, pesize, HECMW_INT, MASTER_PE, HECMW_ANY_TAG, VIS_COMM, &stat);
			}
			if(mynode==0) {
				ndis_c_v=(double *)HECMW_calloc(pesize, sizeof(double));
				if(ndis_c_v==NULL)
					HECMW_vis_memory_exit("ndis_c_v");
				ndis_c_v[0]=dis_c_v;
				for(i=1;i<pesize;i++) {
					HECMW_Recv(&dis_c_v, 1, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);
					ndis_c_v[i]=dis_c_v;
				}
				pe_id=(int *)HECMW_calloc(pesize, sizeof(int));
				if(pe_id==NULL)
					HECMW_vis_memory_exit("pe_id");
				for(i=0;i<pesize;i++)
					pe_id[i]=i;
				for(i=0;i<pesize-1;i++)
					for(j=i+1;j<pesize;j++) {
						if(ndis_c_v[i]>ndis_c_v[j]) {
							tmpd=ndis_c_v[i];
							ndis_c_v[i]=ndis_c_v[j];
							ndis_c_v[j]=tmpd;
							tmpi=pe_id[i];
							pe_id[i]=pe_id[j];
							pe_id[j]=tmpi;
						}
					}
				for(i=1;i<pesize;i++)
					HECMW_Send(pe_id,pesize,HECMW_INT, i, 0, VIS_COMM);
			}
			/*	}
			 */

			composite_subimage_vr(pesize, pe_id, pixn, n_subimage, n_subopa, subimage);
			HECMW_free(n_subimage);
			HECMW_free(n_subopa);
			HECMW_free(subopa);
			/*send subimage to master PE */
			if(mynode!=0)
				HECMW_Send(subimage, pixn*3, HECMW_DOUBLE, 0, 0, VIS_COMM);
			if(mynode==0) {
				for(j=0;j<pix0_num*3;j++)
					image[j]=subimage[j];
				for(i=1;i<pesize;i++) {
					HECMW_Recv(subimage, pix_num*3, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);
					for(j=0;j<pix_num*3;j++)
						image[pix0_num*3+(i-1)*pix_num*3+j]=subimage[j];
				}
			}
			HECMW_free(subimage);
		}
		if(mynode==0) {
			fprintf(stderr, " Finish combining \n");
			t3=HECMW_Wtime();
			fprintf(stderr, "The combining time is %lf the total time up to now is %lf\n", t3-t2, t3-t1);
			t2=t3;
			if(vr->fixed_scale_mark==1)
				generate_color_bar(vr->scale_marking_on, vr->font_size, vr->color_bar_style, vr->mark_0_on, vr->color_mapping_bar_on,
						vr->xr, vr->yr, vr->font_color, vr->color_system_type, vr->color_mapping_style, vr->interval_point,
						vr->interval_mapping_num, vr->num_of_scale, tmincolor, tmaxcolor, tmincolor, tmaxcolor,image);
			else if(vr->fixed_scale_mark==0)
				generate_color_bar(vr->scale_marking_on, vr->font_size, vr->color_bar_style, vr->mark_0_on, vr->color_mapping_bar_on,
						vr->xr, vr->yr, vr->font_color, vr->color_system_type, vr->color_mapping_style, vr->interval_point,
						vr->interval_mapping_num, vr->num_of_scale, tmincolor, tmaxcolor,org_mincolor, org_maxcolor, image);

			/*      outfile1[0] = '\0';
	  j=0;
      while (outfile[j] == ' ')
	    j++;
      k = 0;
      while (outfile[j] != ' ') {
    	outfile1[k] = outfile[j];
    	k++;
    	j++;
      }

      outfile1[k] = '\0';

	  sprintf(fname, "%s.%d.%d.bmp", outfile1, *timestep, ii);
			 */
			if(*timestep>=1000)
				sprintf(timestep_str, "%d", *timestep);
			else if(*timestep>=100)
				sprintf(timestep_str, "0%d", *timestep);
			else if(*timestep>=10)
				sprintf(timestep_str, "00%d", *timestep);
			else
				sprintf(timestep_str, "000%d", *timestep);
			if(ii>=10)
				sprintf(rotate_str, "%d", ii);
			else
				sprintf(rotate_str, "0%d", ii);
			sprintf(fname, "%s.%s.%s.bmp", outfile, timestep_str, rotate_str);

			FP = fopen(fname, "wb");

			if(!FP)
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E0009: Cannot open output file");


			header.bfSize=54+3*vr->xr*vr->yr;
#ifdef CONVERSE_ORDER
			header.bfSize=change_unsigned_int_order(54+3*vr->xr*vr->yr);
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

			info.biWidth=vr->xr;
#ifdef CONVERSE_ORDER
			info.biWidth=change_int_order(vr->xr);
#endif

			info.biHeight=vr->yr;
#ifdef CONVERSE_ORDER
			info.biHeight=change_int_order(vr->yr);
#endif

			info.biSizeImage=3*vr->xr*vr->yr;
#ifdef CONVERSE_ORDER
			info.biSizeImage=change_unsigned_int_order(3*vr->xr*vr->yr);
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

			putc('B', FP);
			putc('M', FP);
			fwrite(&(header.bfSize), sizeof(unsigned int), 1,FP);
			fwrite(&header.bfReserved1, sizeof(unsigned short int), 1, FP);
			fwrite(&header.bfReserved2, sizeof(unsigned short int), 1, FP);
			fwrite(&header.bfOffBits, sizeof(unsigned int), 1, FP);
			fwrite(&info, 40, 1, FP);
			for(j=0;j<vr->yr;j++)
				for(i=vr->xr-1;i>=0;i--) {
					if((image[(j*vr->xr+i)*3]<0.1) && (image[(j*vr->xr+i)*3+1]<0.1) && (image[(j*vr->xr+i)*3+2]<0.1)) {
						ri=(int)(vr->background_color[0]*256);
						gi=(int)(vr->background_color[1]*256);
						bi=(int)(vr->background_color[2]*256);
					}
					else {
						ri=(int)(image[(j*vr->xr+i)*3]*256);
						gi=(int)(image[(j*vr->xr+i)*3+1]*256);
						bi=(int)(image[(j*vr->xr+i)*3+2]*256);
					}
					if(ri<0) ri=0;
					if(ri>255) ri=255;
					if(gi<0) gi=0;
					if(gi>255) gi=255;
					if(bi<0) bi=0;
					if(bi>255) bi=255;
					r=ri;  g=gi; b=bi;
					putc(b, FP);
					putc(g, FP);
					putc(r, FP);
				}


			fclose(FP);

		}
	}
	HECMW_free(var);
#ifdef slow
	HECMW_free(grad_var);
#endif
	HECMW_free(empty_flag);
	HECMW_free(image);
	HECMW_free(accum_opa);
	if(num_of_pvr>1) {
		HECMW_free(voxel_dxyz);
		HECMW_free(voxel_orig_xyz);
		HECMW_free(level);
		HECMW_free(voxel_n_neighbor_pe);
		HECMW_free(voxel_neighbor_pe);
		HECMW_free(extent);
	}
	if(mynode==0) {
		fprintf(stderr, " Finish the module\n");
		t3=HECMW_Wtime();
		fprintf(stderr, " the time for output is %lf the total time of the module is %lf\n", t3-t2, t3-t1);
	}
	return;
}

