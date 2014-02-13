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

#include "hecmw_vis_new_refine.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_vis_resampling.h"
#include "hecmw_vis_ray_trace.h"
#include "hecmw_malloc.h"


/* static void free_headpatch(Head_surfacep_info *head_patch, int nx, int ny, int nz); */
static void free_cubehead(Cube_head *cubehead, int voxn);

void transform_face_node(int face, int node[4])
{

	switch (face) {
	case 0:
		node[0]=0;
		node[1]=4;
		node[2]=7;
		node[3]=3;
		break;
	case 1:
		node[0]=1; node[1]=5;
		node[2]=6; node[3]=2;
		break;
	case 2:
		node[0]=0; node[1]=1;
		node[2]=5; node[3]=4;
		break;
	case 3:
		node[0]=3; node[1]=2;
		node[2]=6; node[3]=7;
		break;
	case 4:
		node[0]=0; node[1]=1;
		node[2]=2; node[3]=3;
		break;
	case 5:
		node[0]=4; node[1]=5;
		node[2]=6; node[3]=7;
		break;
	}
	return;
}

static void transform_face_node_351(int face, int node[4])
{

	switch (face) {
	case 0:
		node[0]=0;
		node[1]=1;
		node[2]=4;
		node[3]=3;
		break;
	case 1:
		node[0]=1; node[1]=2;
		node[2]=5; node[3]=4;
		break;
	case 2:
		node[0]=2; node[1]=0;
		node[2]=3; node[3]=5;
		break;
	case 3:
		node[0]=0; node[1]=2;
		node[2]=1; node[3]=1;
		break;
	case 4:
		node[0]=3; node[1]=4;
		node[2]=5; node[3]=5;
		break;
	}
	return;
}


static void value2_compute(struct hecmwST_local_mesh *mesh, double *node1, double x, double y, double z,
		int elem_ID, double *field)
{
	int		i, j, k;
	double	dist[MAX_N_NODE], coef[MAX_N_NODE];
	double	coord_x, coord_y, coord_z;
	double	v_data[MAX_N_NODE];
	int		nodeID;
	double	sum_coef;
	int           return_flag;

	*field = 0.0;
	return_flag=0;

	i=elem_ID;
	for (j = 0; j < mesh->elem_node_index[i+1]-mesh->elem_node_index[i]; j++) {
		dist[j] = 0.0;
		nodeID = mesh->elem_node_item[mesh->elem_node_index[i]+j]-1;
		coord_x = mesh->node[nodeID*3];
		coord_y = mesh->node[nodeID*3+1];
		coord_z = mesh->node[nodeID*3+2];
		v_data[j] = node1[nodeID];

		dist[j] = (x - coord_x)*(x - coord_x)
		+ (y - coord_y)*(y - coord_y)
		+ (z - coord_z)*(z - coord_z);
		if (fabs(dist[j]) <EPSILON) {
			*field = v_data[j];
			return_flag=1;
		}
	}
	if(return_flag==0) {
		sum_coef = 0.0;
		for (j = 0; j <mesh->elem_node_index[i+1]-mesh->elem_node_index[i]; j++) {
			coef[j] = 1.0;
			for (k = 0; k < mesh->elem_node_index[i+1]-mesh->elem_node_index[i]; k++) {
				if (j != k) coef[j] = coef[j]*dist[k];
			}
			sum_coef += coef[j];
		}

		for (j = 0; j < mesh->elem_node_index[i+1]-mesh->elem_node_index[i]; j++) {
			*field = *field + v_data[j] * coef[j] / sum_coef;
		}
	}
	return;
}


static int judge_inner_voxel_361(double vv[MAX_N_NODE*3], double point[3])
{
	double fp[4][3], n_f[3], cen_point[3], sign_f, sign_p, f_cen_point[3], c_c[3], p_c[3];
	double n_norm, n_c, n_p;
	int i,j, m, inside,node[4];

	/*find the center point of the element */
	for(j=0;j<3;j++) {
		cen_point[j]=0.0;
		for(i=0;i<8;i++)
			cen_point[j]+=vv[i*3+j];
		cen_point[j]/=8.0;
	}
	inside=1;
	/*for each face of the element */
	for(m=0;m<6;m++) {
		transform_face_node(m,node);
		for(j=0;j<3;j++)
			for(i=0;i<4;i++){
				fp[i][j]=vv[node[i]*3+j];
			}
		/*find the normal of the face */
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[3][2]-fp[0][2])-(fp[3][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[3][2]-fp[0][2])+(fp[3][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[3][1]-fp[0][1])-(fp[3][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
		/*selce the direction point to inside the element */
		for(j=0;j<3;j++)
			f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j]+fp[3][j])/4.0;
		for(j=0;j<3;j++)
			c_c[j]=cen_point[j]-f_cen_point[j];
		n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
		if(fabs(n_c)>EPSILON) {
			for(j=0;j<3;j++)
				c_c[j]/=n_c;
		}
		sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
		if(sign_f<-EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]=-n_f[j];
		}
		for(j=0;j<3;j++)
			p_c[j]=point[j]-f_cen_point[j];
		n_p=sqrt(p_c[0]*p_c[0]+p_c[1]*p_c[1]+p_c[2]*p_c[2]);
		if(fabs(n_p)>EPSILON) {
			for(j=0;j<3;j++)
				p_c[j]/=n_p;
		}
		sign_p=p_c[0]*n_f[0]+p_c[1]*n_f[1]+p_c[2]*n_f[2];
		if(sign_p<-EPSILON) {
			inside=0;
			break;
		}
	}
	return(inside);
}

static int judge_inner_voxel_351(double vv[MAX_N_NODE*3], double point[3])
{
	double fp[4][3], n_f[3], cen_point[3], sign_f, sign_p, f_cen_point[3], c_c[3], p_c[3];
	double n_norm, n_c, n_p;
	int i,j, m, inside,node[4];

	/*find the center point of the element */
	for(j=0;j<3;j++) {
		cen_point[j]=0.0;
		for(i=0;i<6;i++)
			cen_point[j]+=vv[i*3+j];
		cen_point[j]/=6.0;
	}
	inside=1;
	/*for each face of the element */
	for(m=0;m<3;m++) {
		transform_face_node_351(m,node);
		for(j=0;j<3;j++)
			for(i=0;i<4;i++){
				fp[i][j]=vv[node[i]*3+j];
			}
		/*find the normal of the face */
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[3][2]-fp[0][2])-(fp[3][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[3][2]-fp[0][2])+(fp[3][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[3][1]-fp[0][1])-(fp[3][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
		/*selce the direction point to inside the element */
		for(j=0;j<3;j++)
			f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j]+fp[3][j])/4.0;
		for(j=0;j<3;j++)
			c_c[j]=cen_point[j]-f_cen_point[j];
		n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
		if(fabs(n_c)>EPSILON) {
			for(j=0;j<3;j++)
				c_c[j]/=n_c;
		}
		sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
		if(sign_f<-EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]=-n_f[j];
		}
		for(j=0;j<3;j++)
			p_c[j]=point[j]-f_cen_point[j];
		n_p=sqrt(p_c[0]*p_c[0]+p_c[1]*p_c[1]+p_c[2]*p_c[2]);
		if(fabs(n_p)>EPSILON) {
			for(j=0;j<3;j++)
				p_c[j]/=n_p;
		}
		sign_p=p_c[0]*n_f[0]+p_c[1]*n_f[1]+p_c[2]*n_f[2];
		if(sign_p<-EPSILON) {
			inside=0;
			break;
		}
	}
	for(m=3;m<5;m++) {
		for(j=0;j<3;j++)
			for(i=0;i<3;i++){
				fp[i][j]=vv[node[i]*3+j];
			}
		/*find the normal of the face */
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
		/*selce the direction point to inside the element */
		for(j=0;j<3;j++)
			f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j])/3.0;
		for(j=0;j<3;j++)
			c_c[j]=cen_point[j]-f_cen_point[j];
		n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
		if(fabs(n_c)>EPSILON) {
			for(j=0;j<3;j++)
				c_c[j]/=n_c;
		}
		sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
		if(sign_f<-EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]=-n_f[j];
		}
		for(j=0;j<3;j++)
			p_c[j]=point[j]-f_cen_point[j];
		n_p=sqrt(p_c[0]*p_c[0]+p_c[1]*p_c[1]+p_c[2]*p_c[2]);
		if(fabs(n_p)>EPSILON) {
			for(j=0;j<3;j++)
				p_c[j]/=n_p;
		}
		sign_p=p_c[0]*n_f[0]+p_c[1]*n_f[1]+p_c[2]*n_f[2];
		if(sign_p<-EPSILON) {
			inside=0;
			break;
		}
	}

	return(inside);
}

static int judge_inner_voxel_341(double vv[MAX_N_NODE*3], double point[3])
{
	double fp[4][3], n_f[3], cen_point[3], sign_f, sign_p, f_cen_point[3], c_c[3], p_c[3];
	double n_norm, n_c, n_p;
	int i,j, m, inside,node[4];

	/*find the center point of the element */
	for(j=0;j<3;j++) {
		cen_point[j]=0.0;
		for(i=0;i<4;i++)
			cen_point[j]+=vv[i*3+j];
		cen_point[j]/=4.0;
	}
	inside=1;
	/*for each face of the element */
	for(m=0;m<4;m++) {
		if(m==0) {
			node[0]=0;  node[1]=2;  node[2]=1;
		}
		else if(m==1) {
			node[0]=0;  node[1]=1;  node[2]=3;
		}
		else if(m==2) {
			node[0]=1;  node[1]=2;  node[2]=3;
		}
		else if(m==3) {
			node[0]=2;  node[1]=0;  node[2]=3;
		}

		for(j=0;j<3;j++)
			for(i=0;i<3;i++){
				fp[i][j]=vv[node[i]*3+j];
			}
		/*find the normal of the face */
		n_f[0]=  (fp[1][1]-fp[0][1])*(fp[2][2]-fp[0][2])-(fp[2][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		n_f[1]= -(fp[1][0]-fp[0][0])*(fp[2][2]-fp[0][2])+(fp[2][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		n_f[2]=  (fp[1][0]-fp[0][0])*(fp[2][1]-fp[0][1])-(fp[2][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		n_norm=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
		if(fabs(n_norm)>EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]/=n_norm;
		}
		/*selce the direction point to inside the element */
		for(j=0;j<3;j++)
			f_cen_point[j]=(fp[0][j]+fp[1][j]+fp[2][j])/3.0;
		for(j=0;j<3;j++)
			c_c[j]=cen_point[j]-f_cen_point[j];
		n_c=sqrt(c_c[0]*c_c[0]+c_c[1]*c_c[1]+c_c[2]*c_c[2]);
		if(fabs(n_c)>EPSILON) {
			for(j=0;j<3;j++)
				c_c[j]/=n_c;
		}
		sign_f=n_f[0]*c_c[0]+n_f[1]*c_c[1]+n_f[2]*c_c[2];
		if(sign_f<-EPSILON) {
			for(j=0;j<3;j++)
				n_f[j]=-n_f[j];
		}
		for(j=0;j<3;j++)
			p_c[j]=point[j]-f_cen_point[j];
		n_p=sqrt(p_c[0]*p_c[0]+p_c[1]*p_c[1]+p_c[2]*p_c[2]);
		if(fabs(n_p)>EPSILON) {
			for(j=0;j<3;j++)
				p_c[j]/=n_p;
		}
		sign_p=p_c[0]*n_f[0]+p_c[1]*n_f[1]+p_c[2]*n_f[2];
		if(sign_p<-EPSILON) {
			inside=0;
			break;
		}
	}
	return(inside);
}

/*
void find2_minmax(In_surface *surface, int i, double minmax[6])
{
	int j;


	for(j=0;j<3;j++)
		minmax[j*2]=minmax[j*2+1]=surface->vertex[(surface->patch[i*3]-1)*3+j];
	for(j=0;j<3;j++) {
		if(surface->vertex[(surface->patch[i*3+1]-1)*3+j]<minmax[j*2])
			minmax[j*2]=surface->vertex[(surface->patch[i*3+1]-1)*3+j];
		if(surface->vertex[(surface->patch[i*3+1]-1)*3+j]>minmax[j*2+1])
	        minmax[j*2+1]=surface->vertex[(surface->patch[i*3+1]-1)*3+j];
		if(surface->vertex[(surface->patch[i*3+2]-1)*3+j]<minmax[j*2])
			minmax[j*2]=surface->vertex[(surface->patch[i*3+2]-1)*3+j];
		if(surface->vertex[(surface->patch[i*3+2]-1)*3+j]>minmax[j*2+1])
	        minmax[j*2+1]=surface->vertex[(surface->patch[i*3+2]-1)*3+j];
	}

	return;
}

 */



void refinement(struct hecmwST_local_mesh *mesh, double *node1,
		int n_voxel, double *voxel_dxyz, double *voxel_orig_xyz, int *level,
		int *voxel_n_neighbor_pe, int **voxel_neighbor_pe, double *extent,
		int my_rank, HECMW_Comm VIS_COMM,  int *empty_flag, double *var)

{
	int		i, j;

	int		pe_size;
	Cube_head *cubehead;
	Cube_point *p1;
	int m1, m2, n1,n2, n3, m;
	int nx, ny, nz, flag;
	int i1, j1,k1, i2, j2,k2;
	double dx, dy, dz, minx, miny,minz, maxx,maxy,maxz, x,y,z;
	double field;
	int *code;
	double *value;
	int vox_id, point_num;
	int pe_id;
	HECMW_Status	stat;
	double vv[MAX_N_NODE*3], in_point[3];
	int inside, nodeID;

	HECMW_Comm_size(VIS_COMM, &pe_size);
	cubehead=(Cube_head *)HECMW_calloc(n_voxel, sizeof(Cube_head));
	if(cubehead==NULL)
		HECMW_vis_memory_exit("In refinement: cubehead");
	for(i=0;i<n_voxel;i++) {
		cubehead[i].point_num=0;
		cubehead[i].cube_link=NULL;
	}


	for(m1=0;m1<n_voxel;m1++)
		for(m2=0;m2<voxel_n_neighbor_pe[m1];m2++) {
			if(voxel_neighbor_pe[m1][m2]==my_rank) {
				nx=level[m1*3];
				ny=level[m1*3+1];
				nz=level[m1*3+2];
				dx=voxel_dxyz[m1*3]/(double)nx;
				dy=voxel_dxyz[m1*3+1]/(double)ny;
				dz=voxel_dxyz[m1*3+2]/(double)nz;
				/*                fprintf(stderr, "in voxel %d: nx, ny, nz=%d %d %d dx, dy, dz=%lf %lf %lf\n", m1,nx, ny, nz, dx, dy, dz);
				 */
				minx=voxel_orig_xyz[m1*3];
				maxx=voxel_orig_xyz[m1*3]+voxel_dxyz[m1*3];
				miny=voxel_orig_xyz[m1*3+1];
				maxy=voxel_orig_xyz[m1*3+1]+voxel_dxyz[m1*3+1];
				minz=voxel_orig_xyz[m1*3+2];
				maxz=voxel_orig_xyz[m1*3+2]+voxel_dxyz[m1*3+2];
				for (i = 0; i < mesh->n_elem; i++) {
					if(mesh->elem_type[i]<400) {
						flag=1;
						if((extent[i*6]>maxx) || (extent[i*6+1]<minx) || (extent[i*6+2]>maxy)
								|| (extent[i*6+3]<miny) || (extent[i*6+4]>maxz)
								|| (extent[i*6+5]<minz))
							flag=0;
						if(flag==1) {
							i1=(int)((extent[i*6]-EPSILON-minx)/dx);
							j1=(int)((extent[i*6+2]-EPSILON-miny)/dy);
							k1=(int)((extent[i*6+4]-EPSILON-minz)/dz);
							i2=(int)((extent[i*6+1]+EPSILON-minx)/dx);
							j2=(int)((extent[i*6+3]+EPSILON-miny)/dy);
							k2=(int)((extent[i*6+5]+EPSILON-minz)/dz);
							if(i1<0) i1=0;
							if(i1>nx) i1=i2+1;
							if(j1<0) j1=0;
							if(j1>ny) j1=j2+1;
							if(k1<0) k1=0;
							if(k1>nz) k1=k2+1;
							if(i2<0) i2=i1-1;
							if(i2>nx) i2=nx;
							if(j2<0) j2=j1-1;
							if(j2>ny) j2=ny;
							if(k2<0) k2=k1-1;
							if(k2>nz) k2=nz;
							for(n3=k1;n3<=k2;n3++)
								for(n2=j1;n2<=j2;n2++)
									for(n1=i1;n1<=i2;n1++) {
										x=minx+dx*n1;
										y=miny+dy*n2;
										z=minz+dz*n3;
										/*								if((x>=extent[i].x_min-EPSILON) && (x<=extent[i].x_max+EPSILON) &&
									(y>=extent[i].y_min-EPSILON) && (y<=extent[i].y_max+EPSILON) &&
									(z>=extent[i].z_min-EPSILON) && (z<=extent[i].z_max+EPSILON)) {
										 */
										in_point[0]=x;
										in_point[1]=y;
										in_point[2]=z;
										for(m=mesh->elem_node_index[i];m<mesh->elem_node_index[i+1];m++) {
											nodeID = mesh->elem_node_item[m] - 1;
											vv[(m-mesh->elem_node_index[i])*3]=mesh->node[nodeID*3];
											vv[(m-mesh->elem_node_index[i])*3+1]=mesh->node[nodeID*3+1];
											vv[(m-mesh->elem_node_index[i])*3+2]=mesh->node[nodeID*3+2];
										}
										if((mesh->elem_type[i]==361) || (mesh->elem_type[i]==362))
											inside=judge_inner_voxel_361(vv, in_point);
										else if((mesh->elem_type[i]==341) || (mesh->elem_type[i]==342))
											inside=judge_inner_voxel_341(vv, in_point);
										else if((mesh->elem_type[i]==351) || (mesh->elem_type[i]==352))
											inside=judge_inner_voxel_351(vv, in_point);

										/*------------need add judge inner for tetrahedras later !!!!!!!!!!!!!!!!!!!!
-----------------*/
										if(inside==1) {

											cubehead[m1].point_num++;
											p1=(Cube_point *)HECMW_malloc(sizeof(Cube_point));
											if(p1==NULL)
												HECMW_vis_memory_exit("new_refine: p1");
											p1->code[0]=n1;
											p1->code[1]=n2;
											p1->code[2]=n3;

											value2_compute(mesh, node1, x, y, z, i, &field);
											p1->field=field;
											p1->next_point=cubehead[m1].cube_link;
											cubehead[m1].cube_link=p1;
										} /*end if(xyz) */
									} /*end of n1, n2, n3 loop */
						} /*end of flag=1 */
					} /*end of element type */
				} /* end of element loop */
			} /* end of if */
		} /* end of m1, m2 loop */
	/*		  pe_vox=(int *)HECMW_calloc(pe_size, sizeof(int));
		  if(pe_vox==NULL) {
			  fprintf(stderr, "There is no enough memory for pe_vox\n");
			  exit(0);
		  }
		  pe_vox[0]=7; pe_vox[1]=6; pe_vox[2]=5; pe_vox[3]=4; pe_vox[4]=3;
		  pe_vox[5]=2; pe_vox[6]=1; pe_vox[7]=0;
	 */
	HECMW_Barrier(VIS_COMM);

	/* start send cube points for each voxel in this PE */

	/*		    if (nflag == 0) {
    if ((sta1 = (HECMW_Status *)HECMW_calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status)))
	== NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
    if ((sta2 = (HECMW_Status *)HECMW_calloc(HECMW_STATUS_SIZE, sizeof(HECMW_Status)))
	== NULL) {
    fprintf(stderr, "Not enough memory\n");
    exit(1);
    }
    if ((req1 = (HECMW_Request *)HECMW_calloc(pe_size-1, sizeof(HECMW_Request)))
	== NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
    if ((req2 = (HECMW_Request *)HECMW_calloc(pe_size-1, sizeof(HECMW_Request)))
	== NULL) {
      fprintf(stderr, "Not enough memory\n");
      exit(1);
    }
    nflag = 1;
  }
	 */
	for(m1=0;m1<n_voxel;m1++) {
		pe_id=m1 % pe_size;
		/*                  fprintf(stderr, "m1=%d pe_size=%d\n", m1, pe_size);
		 */
		/*		  pe_id=pe_vox[m1];
		 */
		if(pe_id!=my_rank) {
			if(cubehead[m1].point_num==0) {
				HECMW_Send(&m1,1,HECMW_INT, pe_id, 0, VIS_COMM);
				HECMW_Send(&(cubehead[m1].point_num),1,HECMW_INT, pe_id, 0, VIS_COMM);

				/*			HECMW_Isend(&m1,1,HECMW_INT, pe_id, 0, VIS_COMM,&req1[pe_id]);
            HECMW_Isend(&(cubehead[m1].point_num),1,HECMW_INT, pe_id, 0, VIS_COMM,&req1[pe_id]);
				 */
			}
			else if(cubehead[m1].point_num>0) {
				code=(int *)HECMW_calloc(3*cubehead[m1].point_num, sizeof(int));
				value=(double *)HECMW_calloc(cubehead[m1].point_num, sizeof(double));
				if((code==NULL) || (value==NULL))
					HECMW_vis_memory_exit("new_refine: code, value");
				p1=cubehead[m1].cube_link;
				for(i=0;i<cubehead[m1].point_num;i++) {
					for(j=0;j<3;j++)
						code[i*3+j]=p1->code[j];
					value[i]=p1->field;
					p1=p1->next_point;
				}
				HECMW_Send(&m1, 1, HECMW_INT, pe_id, 0, VIS_COMM);
				HECMW_Send(&(cubehead[m1].point_num),1,HECMW_INT, pe_id, 0, VIS_COMM);
				HECMW_Send(code,3*cubehead[m1].point_num,HECMW_INT, pe_id, 0, VIS_COMM);
				HECMW_Send(value, cubehead[m1].point_num,HECMW_DOUBLE, pe_id, 0, VIS_COMM);

				/*			HECMW_Isend(&m1, 1, HECMW_INT, pe_id, 0, VIS_COMM,&req1[pe_id]);
            HECMW_Isend(&(cubehead[m1].point_num),1,HECMW_INT, pe_id, 0, VIS_COMM,&req1[pe_id]);
            HECMW_Isend(code,3*cubehead[m1].point_num,HECMW_INT, pe_id, 0, VIS_COMM,&req1[pe_id]);
			HECMW_Isend(value, cubehead[m1].point_num,HECMW_DOUBLE, pe_id, 0, VIS_COMM,&req1[pe_id]);
			HECMW_Isend(gradient, 3*cubehead[m1].point_num,HECMW_DOUBLE, pe_id, 0, VIS_COMM,&req1[pe_id]);
				 */			HECMW_free(code);
				 HECMW_free(value);
			} /*end of else ifnum>0*/
			/*		  fprintf(stderr, "on pe %d finish send vox %d to pe %d\n", my_rank, m1, pe_id);
			 */
		} /* end of if rank!=mype */
		else if(pe_id==my_rank) {
			/*		  		max_level = vox->info[m1].level[0];
	            if (max_level < vox->info[m1].level[1])
	            max_level = vox->info[m1].level[1];
	            if (max_level < vox->info[m1].level[2])
	            max_level = vox->info[m1].level[2];
			 */
			/*		nx=ny=nz= (int)pow(2.0,(double)max_level);
			 */
			nx=level[m1*3];
			ny=level[m1*3+1];
			nz=level[m1*3+2];
			/*   if(mynode==0)
   fprintf(stderr, "Final nx, ny, nz in refinement is %d %d %d\n", nx, ny, nz);
			 */




			for(i=0;i<(nx+1)*(ny+1)*(nz+1); i++)
				empty_flag[i]=0;
			/*first copy the part result on this pe into result file */
			if(cubehead[m1].point_num>0) {
				p1=cubehead[m1].cube_link;
				for(i=0;i<cubehead[m1].point_num;i++) {
					i1=p1->code[0];
					j1=p1->code[1];
					k1=p1->code[2];
					if((k1*(ny+1)*(nx+1)+j1*(nx+1)+i1<0) || (k1*(ny+1)*(nx+1)+j1*(nx+1)+i1>
					(nx+1)*(ny+1)*(nz+1))) {
						fprintf(stderr, "There is wrong with code code=%d %d %d\n", i1, j1, k1);
						HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2004: voxel search error in new_refine");
					}
					/*					   result[(k1*(ny+1)*(nx+1)+j1*(nx+1)+i1)*4]=p1->field;
					 */                     var[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=p1->field;
					 /*
					   for(j=0;j<3;j++)
						   result[(k1*(ny+1)*(nx+1)+j1*(nx+1)+i1)*4+j+1]=p1->grad[j];
					   empty_flag[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=1;
					  */
					 empty_flag[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=1;
					 p1=p1->next_point;
				}
			}
			/*second copy the part result from other pes into result file */

			for(i=0;i<pe_size;i++)
				if(i!=my_rank) {
					HECMW_Recv(&vox_id, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
					/*
					    HECMW_Irecv(&vox_id, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &req2[i]);
					 *//*
					 fprintf(stderr, "vox=%d pe=%d from pe %d\n", m1, my_rank, i);
					  */
					if(vox_id!=m1)
						HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2005: data communication error in new_refine");
					HECMW_Recv(&point_num, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);

					/*					 HECMW_Irecv(&point_num, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &req2[i]);
					 */					 if(point_num>0) {
						 code=(int *)HECMW_calloc(3*point_num, sizeof(int));
						 value=(double *)HECMW_calloc(point_num, sizeof(double));
						 if((code==NULL) || (value==NULL))
							 HECMW_vis_memory_exit("new_refine: code, value");
						 HECMW_Recv(code, point_num*3, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
						 HECMW_Recv(value, point_num, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &stat);
						 /*
					 HECMW_Irecv(code, point_num*3, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &req2[i]);
					 HECMW_Irecv(value, point_num, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM,&req2[i]);
					 HECMW_Irecv(gradient, point_num*3, HECMW_DOUBLE, i, HECMW_ANY_TAG, VIS_COMM, &req2[i]);
						  */                   for(j=0;j<point_num;j++) {
							  i1=code[j*3];
							  j1=code[j*3+1];
							  k1=code[j*3+2];
							  /*					   result[(k1*(ny+1)*(nx+1)+j1*(nx+1)+i1)*4]=value[j];
					   for(k=0;k<3;k++)
						   result[(k1*(ny+1)*(nx+1)+j1*(nx+1)+i1)*4+k+1]=gradient[j*3+k];
					   empty_flag[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=1;
							   */
							  var[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=value[j];
							  empty_flag[k1*(ny+1)*(nx+1)+j1*(nx+1)+i1]=1;

						  }
						  HECMW_free(code);
						  HECMW_free(value);
					 }/* end of ifpoint_num>0*/
				} /* end of pe loop */
			/*				   fprintf(stderr, "Finish to receive\n");
			 */
			/* third, organize to octree */
			/*
		  HECMW_free(sta1);
  HECMW_free(sta2);
  HECMW_free(req1);
  HECMW_free(req2);
			 */
			/* fourth, output result file */
			/*	sprintf(filename, "%s-%d.%d", outfile, node->iter, m1);
	if ((fp = fopen(filename, "w")) == NULL) {
	  fprintf(stderr, "output file error\n");
	  exit(1);
	}
  fprintf(fp, " %lf %lf %lf\n", vox->info[m1].dx, vox->info[m1].dy, vox->info[m1].dz);
  fprintf(fp, " %d\n", max_level);
    fprintf(fp, "    %d\n", cont->n_var);
  for (i = 0; i < cont->n_var; i++) {
    fprintf(fp, "    %s\n", cont->name[i]);
  }
  fprintf(fp, "    %d %d %d\n", 1, 1, 1);
  fprintf(fp, "    %lf %lf %lf\n", vox->info[m1].orig_x, vox->info[m1].orig_y, vox->info[m1].orig_z);
  fprintf(fp, "%d %d %d\n", nx, ny, nz);
  for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
	  fprintf(fp, "%d\n", empty_flag[i]);
	  if(empty_flag[i]==0)
		  fprintf(fp, "%d\n", i);
           if(empty_flag[i]==1)
		  fprintf(fp, "%lf %lf %lf %lf\n", result[i*4], result[i*4+1], result[i*4+2], result[i*4+3]);
  }
  fclose(fp);
			 */

		} /*end of if=my_rank*/
		HECMW_Barrier(VIS_COMM);
	}

	/* HECMW_free application memory */
	free_cubehead(cubehead, n_voxel);

	return;
}

#if 0
static void free_headpatch(Head_surfacep_info *head_patch, int nx, int ny, int nz)
{
	int i;
	Surfacep_info *p1, *p2;

	for(i=0;i<nx*ny*nz;i++) {
		if(head_patch[i].num_of_patch>0) {
			p1=head_patch[i].next_patch;
			while(p1!=NULL) {
				p2=p1;
				p1=p1->next_patch;
				HECMW_free(p2);
			}
		}
	}
	HECMW_free(head_patch);
	return;
}
#endif

static void free_cubehead(Cube_head *cubehead, int voxn)
{
	int i;
	Cube_point *p1, *p2;

	for(i=0;i<voxn;i++) {
		if(cubehead[i].point_num>0) {
			p1=cubehead[i].cube_link;
			while(p1!=NULL) {
				p2=p1;
				p1=p1->next_point;
				HECMW_free(p2);
			}
		}
	}
	HECMW_free(cubehead);
	return;
}



