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

#include "hecmw_vis_voxel_gen.h"

#include <stdio.h>
#include <stdlib.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_ray_trace.h"
#include "hecmw_malloc.h"


void voxel_gen(double range[6], double c_range[2], int nv[3], double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, HECMW_Comm VIS_COMM, int vox_on, int display_range_on, double display_range[6])
{

	int			my_rank;
	int			pe_size;

	FILE			*fp;

	int			i, j,k, i1;

	double         trange[6], tc_range[2];
	double        *s_range, tmp_range[6];
	int           n_voxel;
	double        dx, dy, dz;
	int           flag, num_neipe, *neipe;
	double        minx, miny, minz, maxx, maxy, maxz;
	HECMW_Status	stat;



	HECMW_Comm_rank(VIS_COMM, &my_rank);
	HECMW_Comm_size(VIS_COMM, &pe_size);


	s_range=(double *)HECMW_calloc(6*pe_size, sizeof(double));
	if(s_range==NULL)
		HECMW_vis_memory_exit("s_range");
	if(display_range_on==0) {
		HECMW_Allreduce(&range[0], &trange[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&range[1], &trange[1], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
		HECMW_Allreduce(&range[2], &trange[2], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&range[3], &trange[3], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
		HECMW_Allreduce(&range[4], &trange[4], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
		HECMW_Allreduce(&range[5], &trange[5], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	}
	else if(display_range_on==1) {
		for(i=0;i<6;i++)
			trange[i]=display_range[i];
	}
	if(my_rank!=0) {
		HECMW_Send(range, 6, HECMW_DOUBLE, 0, 0, VIS_COMM);
		/*       HECMW_Recv(s_range, 6*pe_size, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
		 */
	}
	if(my_rank==0) {
		for(i=0;i<6;i++)
			s_range[i]=range[i];
		for(j=1;j<pe_size;j++) {
			HECMW_Recv(tmp_range, 6, HECMW_DOUBLE, j, HECMW_ANY_TAG, VIS_COMM, &stat);
			for(i=0;i<6;i++)
				s_range[j*6+i]=tmp_range[i];
		}
		/*	   for(j=1;j<pe_size;j++)
		   HECMW_Send(s_range, 6*pe_size,HECMW_DOUBLE, j, 0, VIS_COMM);
		 */
	}

	HECMW_Allreduce(&c_range[0], &tc_range[0], 1, HECMW_DOUBLE, HECMW_MIN, VIS_COMM);
	HECMW_Allreduce(&c_range[1], &tc_range[1], 1, HECMW_DOUBLE, HECMW_MAX, VIS_COMM);
	if(my_rank==0) {
		if ((fp = fopen("extent.file", "w")) == NULL)
			HECMW_vis_print_exit("output file name error: extent.file");
		fprintf(fp, "The range of the whole data field is\n");
		fprintf(fp, "Minimum x= %lf    Maximum x=%lf\n", trange[0], trange[1]);
		fprintf(fp, "Minimum y= %lf    Maximum y=%lf\n", trange[2], trange[3]);
		fprintf(fp, "Minimum z= %lf    Maximum z=%lf\n", trange[4], trange[5]);
		fprintf(fp, "The range of color component is\n");
		fprintf(fp, "Minimum color=%lf     maximum color=%lf\n", tc_range[0], tc_range[1]);
		fprintf(stderr, "Minimum color=%lf     maximum color=%lf\n", tc_range[0], tc_range[1]);
		fprintf(fp, "The subrange of each PE\n");
		for(j=0;j<pe_size;j++) {
			fprintf(fp, "PE %d:\n", j);
			fprintf(fp, "Minimum x= %lf    Maximum x=%lf\n", s_range[j*6+0], s_range[j*6+1]);
			fprintf(fp, "Minimum y= %lf    Maximum y=%lf\n", s_range[j*6+2], s_range[j*6+3]);
			fprintf(fp, "Minimum z= %lf    Maximum z=%lf\n", s_range[j*6+4], s_range[j*6+5]);
		}

		fclose(fp);
		if(vox_on==0) {
			n_voxel=nv[0]*nv[1]*nv[2];
			/*   fprintf(stderr, "the number of voxel is %d\n", n_voxel);
			 */
			if ((fp = fopen("voxel.file", "w")) == NULL)
				HECMW_vis_print_exit("output voxel file error: voxel.file");
			dx=(trange[1]-trange[0])/nv[0];
			dy=(trange[3]-trange[2])/nv[1];
			dz=(trange[5]-trange[4])/nv[2];
			neipe=(int *)HECMW_calloc(pe_size, sizeof(int));
			if(neipe==NULL)
				HECMW_vis_memory_exit("in voxel_gen: neipe");
			if(n_voxel==0)
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1042: n_voxel_x,n_voxel_y, and n_voxel_z cannot be less or equal to 0");
			/*  vox->info = (Voxel_info *)HECMW_calloc(n_voxel, sizeof(Voxel_info));
  if(vox->info==NULL) {
	  fprintf(stderr, "There is no enough memory for vox\n");
	  exit(0);
  }
  vox->n_voxel=n_voxel;
			 */

			for(k=0;k<nv[2];k++)
				for(j=0;j<nv[1];j++)
					for(i=0;i<nv[0];i++) {
						fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", trange[0]+i*dx, trange[2]+j*dy, trange[4]+k*dz, dx,dy,dz);
						voxel_orig_xyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3]=trange[0]+i*dx;
						voxel_orig_xyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3+1]=trange[2]+j*dy;
						voxel_orig_xyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3+2]=trange[4]+k*dz;
						voxel_dxyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3]=dx;
						voxel_dxyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3+1]=dy;
						voxel_dxyz[(k*nv[1]*nv[0]+j*nv[0]+i)*3+2]=dz;


						/*			   fprintf(stderr, "i j k=%d %d %d\n", i,j,k);
			   start_i=i-1; end_i=i+1;
			   start_j=j-1; end_j=j+1;
			   start_k=k-1; end_k=k+1;
			   if(start_i<0) start_i=0;
			   if(end_i>nv[0]-1) end_i=nv[0]-1;
			   if(start_j<0) start_j=0;
			   if(end_j>nv[1]-1) end_j=nv[1]-1;
			   if(start_k<0) start_k=0;
			   if(end_k>nv[2]-1) end_k=nv[2]-1;
			   fprintf(stderr, "%d %d %d %d %d %d\n", start_i, end_i, start_j, end_j, start_k, end_k);
			   fprintf(fp, "%d\n", (end_i-start_i+1)*(end_j-start_j+1)*(end_k-start_k+1));
			   for(k1=start_k;k1<=end_k;k1++)
				   for(j1=start_j;j1<=end_j;j1++)
					   for(i1=start_i;i1<=end_i;i1++)
						   fprintf(fp, "%d\n", k1*nv[0]*nv[1]+j1*nv[0]+i1);

						 */
						minx=trange[0]+i*dx;  maxx=minx+dx;
						miny=trange[2]+j*dy;  maxy=miny+dy;
						minz=trange[4]+k*dz;  maxz=minz+dz;
						num_neipe=0;
						for(i1=0;i1<pe_size;i1++) {
							flag=1;
							neipe[i1]=-1;
							if((s_range[i1*6]>maxx+EPSILON*dx) || (s_range[i1*6+1]<minx-EPSILON*dx) || (s_range[i1*6+2]>maxy+EPSILON*dy)
									|| (s_range[i1*6+3]<miny-EPSILON*dy) || (s_range[i1*6+4]>maxz+EPSILON*dz)
									|| (s_range[i1*6+5]<minz-EPSILON*dz))
								flag=0;
							if(flag==1) {
								neipe[num_neipe]=i1;
								num_neipe++;
							}
						}
						fprintf(fp, "%d\n", num_neipe);
						voxel_n_neighbor_pe[k*nv[1]*nv[0]+j*nv[0]+i]=num_neipe;
						if(num_neipe>0) {

							for(i1=0;i1<num_neipe;i1++) {
								fprintf(fp, "%d ", neipe[i1]);

								voxel_neighbor_pe[k*nv[1]*nv[0]+j*nv[0]+i][i1]=neipe[i1];

								if((i1 % 8) == 7)
									fprintf(fp, "\n");
							}
							fprintf(fp, "\n");
						}

					}
			fclose(fp);
			/* send vox information to all other PEs */
			for(i=1;i<pe_size;i++) {
				for(j=0;j<n_voxel;j++) {
					HECMW_Send(&(voxel_orig_xyz[j*3]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_orig_xyz[j*3+1]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_orig_xyz[j*3+2]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_dxyz[j*3]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_dxyz[j*3+1]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_dxyz[j*3+2]),1,HECMW_DOUBLE, i, 0, VIS_COMM);
					HECMW_Send(&(voxel_n_neighbor_pe[j]),1,HECMW_INT, i, 0, VIS_COMM);
					if(voxel_n_neighbor_pe[j]>0)
						HECMW_Send(voxel_neighbor_pe[j],n_voxel,HECMW_INT, i, 0, VIS_COMM);
				}
			}
		}


	}
	if(my_rank!=0) {
		if(vox_on==0) {
			n_voxel=nv[0]*nv[1]*nv[2];
			if(n_voxel==0)
				HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E1042: n_voxel_x,n_voxel_y, and n_voxel_z cannot be less or equal to 0");
			for(j=0;j<n_voxel;j++) {
				HECMW_Recv(&(voxel_orig_xyz[j*3]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_orig_xyz[j*3+1]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_orig_xyz[j*3+2]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_dxyz[j*3]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_dxyz[j*3+1]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_dxyz[j*3+2]), 1, HECMW_DOUBLE, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				HECMW_Recv(&(voxel_n_neighbor_pe[j]), 1, HECMW_INT, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				if(voxel_n_neighbor_pe[j]>0) {
					HECMW_Recv(voxel_neighbor_pe[j], n_voxel, HECMW_INT, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
				}
			}
		}
	}



	return;
}
