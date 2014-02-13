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

#include "hecmw_vis_ray_trace.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_color_composite_vr.h"

#if defined(old_version) || defined(old_intersection)
#include "hecmw_vis_new_refine.h"
#endif

/*
void transform_face_node(int ijkn[3], int current_ijk[3], int face_sect[3], double vv[3],
					int next_ijk[3]);
 */
static int find_out_point(double orig_xyz[3], double r_dxyz[3], int ijk[3], double in_point[3],double ray_direction[3], int face_sect[3],
		double out_point[3]);
static void find_next_cell(int ijkn[3], int current_ijk[3], int face_sect[3], double vv[3], int next_ijk[3]);


static int find_first_inter_point(double point_o[3], double view_p[3], double ray_direction[3], double orig_xyz[3], double dxyz[3], int r_level[3], double first_p[3])
{
	int i,j, mincomp, intersection;
	double minmax[6], t[6], mint;

	for(i=0;i<3;i++) {
		minmax[i*2]=orig_xyz[i];
		minmax[i*2+1]=orig_xyz[i]+dxyz[i];
	}
	for(i=0;i<3;i++) {
		if(fabs(ray_direction[i])<EPSILON) {
			t[i*2]=1.0E+17;  t[i*2+1]=1.0E+17;
		}
		else {
			t[i*2]=(minmax[i*2]-view_p[i])/ray_direction[i];
			t[i*2+1]=(minmax[i*2+1]-view_p[i])/ray_direction[i];
		}
	}
#ifdef old_version
for(i=0;i<3;i++) {
	t[i*2]=(minmax[i*2]-view_p[i])/(ray_direction[i]+EPSILON);
	t[i*2+1]=(minmax[i*2+1]-view_p[i])/(ray_direction[i]+EPSILON);
}
#endif

mint=1.0E+17;
mincomp=-1;
for(i=0;i<6;i++) {
	if((mint>t[i]) && (t[i]>EPSILON)) {
		for(j=0;j<3;j++)
			first_p[j]=view_p[j]+t[i]*ray_direction[j];
		if((first_p[0]>=minmax[0]-EPSILON) && (first_p[0]<=minmax[1]+EPSILON) && (first_p[1]>=minmax[2]-EPSILON)
				&& (first_p[1]<=minmax[3]+EPSILON) && (first_p[2]>=minmax[4]-EPSILON) && (first_p[2]<=minmax[5]+EPSILON)){
			mint=t[i];
			mincomp=i;
		}
	}
}
if((mincomp>=0) && (mincomp<=5)) {
	intersection=1;
	for(i=0;i<3;i++)
		first_p[i]=view_p[i]+mint*ray_direction[i];
}
else {
	intersection=0;
}
/*	for(i=0;i<3;i++) {
		if(fabs(out_point[i]-minmax[i*2])<EPSILON)
			face_sect[i]=i*2;
		else if(fabs(out_point[i]-minmax[i*2+1])<EPSILON)
			face_sect[i]=i*2+1;
	}
 */

return (mincomp);
}

#ifdef old_version

int find_first_inter_point(double point_o[3], double view_p[3], VR_data *vd, double first_p[3])
{
	int i, j, m, node[4], spm[6], minsp;
	double vv[8*3], nn, a[6], b[6], c[6], d[6], abc, in_vv[3], in_point[3], fp[4][3];
	double t, st[6], mint, out_point[3], point[3];
	int lc, con1;
	int on_face_flag;

	for(i=0;i<3;i++)
		in_vv[i]=point_o[i]-view_p[i];
	nn=sqrt(SQR(in_vv[0])+SQR(in_vv[1])+SQR(in_vv[2]));
	if(fabs(nn)>EPSILON) {
		for(i=0;i<3;i++)
			in_vv[i]/=nn;
	}
	for(i=0;i<3;i++)
		in_point[i]=view_p[i];
	vv[0*3]=vv[4*3]=vv[7*3]=vv[3*3]=vd->xyz0[0];
	vv[1*3]=vv[5*3]=vv[6*3]=vv[2*3]=vd->xyz0[0]+vd->nxyz[0]*vd->dxyz[0];
	vv[0*3+2]=vv[1*3+2]=vv[2*3+2]=vv[3*3+2]=vd->xyz0[2];
	vv[4*3+2]=vv[5*3+2]=vv[6*3+2]=vv[7*3+2]=vd->xyz0[2]+vd->nxyz[2]*vd->dxyz[2];
	vv[0*3+1]=vv[1*3+1]=vv[5*3+1]=vv[4*3+1]=vd->xyz0[1];
	vv[3*3+1]=vv[2*3+1]=vv[6*3+1]=vv[7*3+1]=vd->xyz0[1]+vd->nxyz[1]*vd->dxyz[1];

	for(m=0;m<6;m++) {
		transform_face_node(m,node);
		for(i=0;i<4;i++)
			for(j=0;j<3;j++)
				fp[i][j]=vv[node[i]*3+j];
		a[m]=(fp[1][1]-fp[0][1])*(fp[3][2]-fp[0][2])-(fp[3][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		b[m]= -(fp[1][0]-fp[0][0])*(fp[3][2]-fp[0][2])+(fp[3][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		c[m]=(fp[1][0]-fp[0][0])*(fp[3][1]-fp[0][1])-(fp[3][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		abc=sqrt(a[m]*a[m]+b[m]*b[m]+c[m]*c[m]);

		a[m]/=abc;
		b[m]/=abc;
		c[m]/=abc;

		d[m]= -(fp[0][0]*a[m]+fp[0][1]*b[m]+fp[0][2]*c[m]);
	}
	lc=-1; con1=1;
	for(m=0;m<6;m++) {
		if(fabs(a[m]*in_vv[0]+b[m]*in_vv[1]+c[m]*in_vv[2])>EPSILON)
		{

			t=(-d[m]-a[m]*in_point[0]-b[m]*in_point[1]-c[m]*in_point[2])/(a[m]*in_vv[0]+
					b[m]*in_vv[1]+c[m]*in_vv[2]);
			/*    printf("t=%f\n",t);*/
			if(t>0){
				for(j=0;j<3;j++)
					point[j]=in_point[j]+t*in_vv[j];
				if((point[0]>vd->xyz0[0]-EPSILON) && (point[0]<vd->xyz0[0]+vd->dxyz[0]*vd->nxyz[0]+EPSILON)
						&& (point[1]>vd->xyz0[1]-EPSILON)
						&& (point[1]<vd->xyz0[1]+vd->dxyz[1]*vd->nxyz[1]+EPSILON) && (point[2]>vd->xyz0[2]-EPSILON)
						&& (point[2]<vd->xyz0[2]+vd->dxyz[2]*vd->nxyz[2]+EPSILON)) {
					lc++;
					st[lc]=t;
					spm[lc]=m;
				}

				/*	 on_face_flag=0;
	 if((m==0) || (m==1)) {
		 if((point[1]>=vd->xyz0[1]-EPSILON)
		  && (point[1]<=vd->xyz0[1]+vd->dxyz[1]*vd->nxyz[1]+EPSILON) && (point[2]>=vd->xyz0[2]-EPSILON)
		  && (point[2]<=vd->xyz0[2]+vd->dxyz[2]*vd->nxyz[2]+EPSILON))
		  on_face_flag=1;
	 }
	 else if((m==2) || (m==3)) {
         if((point[0]>=vd->xyz0[0]-EPSILON) &&
			 (point[0]=<vd->xyz0[0]+vd->dxyz[0]*vd->nxyz[0]+EPSILON)
             && (point[2]>=vd->xyz0[2]-EPSILON)
		     && (point[2]<=vd->xyz0[2]+vd->dxyz[2]*vd->nxyz[2]+EPSILON))
			 on_face_flag=1;
	 }
	 else if((m==4) || (m==5)) {
         if((point[0]>=vd->xyz0[0]-EPSILON) && (point[0]=<vd->xyz0[0]+vd->dxyz[0]*vd->nxyz[0]+EPSILON)
		 && (point[1]>=vd->xyz0[1]-EPSILON)
		  && (point[1]<=vd->xyz0[1]+vd->dxyz[1]*vd->nxyz[1]+EPSILON))
		  on_face_flag=1;
	 }
	 if(on_face_flag==1) {
	  lc++;
	 st[lc]=t;
	 spm[lc]=m;
	 }
				 */

			}
		}
	}/*end for*/
	if(lc==-1) {
		minsp=-1;
	}
	else {

		mint=st[0];
		minsp=spm[0];
		for(m=1;m<=lc;m++) {
			if(st[m]<mint) {
				mint=st[m];
				minsp=spm[m];
			}
		}
		for(j=0;j<3;j++)
			out_point[j]=in_point[j]+mint*in_vv[j];
		/*    transform_face_node(minsp,node);*/
		/*map the out_point to ijk*/
		for(i=0;i<3;i++)
			first_p[i]=out_point[i];
	}
	/*   for(i=0;i<3;i++)
	   ijk[i]=(int)((first_p[i]-vd->xyz0[i])/vd->dxyz[i]);
	 */
	return(minsp);
}
#endif




#ifdef Octree
void search2_leave(Tree_pointer *in_voxel, double out_point[3], Tree_pointer_ptr *next_voxel,
		double ray_direction[3])
{
	int i, index[3], child_no;
	Tree_pointer *p1;
	double mid_value;
	p1=in_voxel;
	while(p1->child!=NULL) {
		for(i=0;i<3;i++) {
			mid_value=(p1->bound_box[i*2+1]-p1->bound_box[i*2])/2.0+p1->bound_box[i*2];
			if(ray_direction[i]<0) {
				if(out_point[i]<=mid_value) index[i]=0;
				else index[i]=1;
			}
			else if(ray_direction[i]>=0) {
				if(out_point[i]<mid_value) index[i]=0;
				else index[i]=1;
			}

		}
		child_no=index[2]*4+index[1]*2+index[0];
		p1=&(p1->child[child_no]);
	}
	*next_voxel=p1;
	return;
}

void search_leave(Tree_pointer *root_tree, int ini_index, double first_p[3],
		Tree_pointer_ptr *voxel_p, double ray_direction[3])
{
	int i;
	int index[3], child_no;
	Tree_pointer *p1;
	double mid_value;

	p1=&(root_tree[ini_index]);
	while(p1->child!=NULL) {
		for(i=0;i<3;i++) {
			mid_value=(p1->bound_box[i*2+1]-p1->bound_box[i*2])/2.0+p1->bound_box[i*2];
			if(ray_direction[i]<0) {
				if(first_p[i]<=mid_value) index[i]=0;
				else index[i]=1;
			}
			else if(ray_direction[i]>=0) {
				if(first_p[i]<mid_value) index[i]=0;
				else index[i]=1;
			}

		}
		child_no=index[2]*4+index[1]*2+index[0];
		p1=&(p1->child[child_no]);
	}
	*voxel_p=p1;

	return;
}
#endif

int find_first_inter(double point_o[3], double view_point_d[3],  int r_level[3], double orig_xyz[3], double dxyz[3],  double r_dxyz[3], double ray_direction[3], double first_p[3], int ijk[3])
{
	int i;
	int ini_face, intersection, con1, next_ijk[3], face_sect[3];
	double out_point[3];

	ini_face=find_first_inter_point(point_o, view_point_d, ray_direction, orig_xyz, dxyz, r_level, first_p);
	if(ini_face==-1) {/*this ray has no intersection with dataset */
		intersection=0;
	}
	else {
		intersection=1;
		for(i=0;i<3;i++)
			ijk[i]=(int)((first_p[i]-orig_xyz[i]-EPSILON)/r_dxyz[i]);
		if((ijk[0]<0) || (ijk[0]>r_level[0]-1) || (ijk[1]<0) || (ijk[1]>r_level[1]-1) || (ijk[2]<0) || (ijk[2]>r_level[2]-1))
			HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2006:There is something wrong on finding the first point");
		for(i=0;i<3;i++)
			out_point[i]=0.0;
		con1=find_out_point(orig_xyz, r_dxyz, ijk, first_p,ray_direction,face_sect,out_point);
		if(con1==0) {
			/* start voxel maybe at the neighbour voxel*/
			find_next_cell(r_level, ijk, face_sect, ray_direction, next_ijk);
			if((next_ijk[0]<0) || (next_ijk[0]>r_level[0]-1) || (next_ijk[1]<0) ||
					(next_ijk[1]>r_level[1]-1) || (next_ijk[2]<0) || (next_ijk[2]>r_level[2]-1))
				intersection=0;

			for(i=0;i<3;i++)
				ijk[i]=next_ijk[i];
		}
	}
	return intersection;
}


#ifdef Octree
int find_first_inter(double point_o[3], Parameter_vr *vr, VR_data *vd,Tree_pointer *root_tree,
		Tree_pointer_ptr *voxel_p, double ray_direction[3], double first_p[3])

{
	int i,ijk[3], ini_index, ijkn[3];
	int ini_face, intersection, con1, face_sect[3], next_ijk[3];
	double view_p[3], out_point[3];
	Tree_pointer_ptr vvv;

	for(i=0;i<3;i++)
		ijkn[i]=vd->nxyz[i];

	for(i=0;i<3;i++)
		view_p[i]=vr->view_point_d[i];
	ini_face=find_first_inter_point(point_o, view_p, vd, first_p);
	if(ini_face==-1) {/*this ray has no intersection with dataset */
		intersection=0;
	}
	else {
		intersection=1;
		/*    for(i=0;i<3;i++)
	   ijk[i]=(int)((first_p[i]-vd->xyz0[i]-EPSILON)/vd->dxyz[i]);
		 */
		if(ini_face==0) {
			ijk[0]=0;
			for(i=1;i<3;i++)
				ijk[i]=(int)((first_p[i]-vd->xyz0[i])/vd->dxyz[i]);
		}
		else if(ini_face==1) {
			ijk[0]=vd->nxyz[0]-1;
			for(i=1;i<3;i++)
				ijk[i]=(int)((first_p[i]-vd->xyz0[i])/vd->dxyz[i]);
		}
		else if(ini_face==2) {
			ijk[1]=0;
			ijk[0]=(int)((first_p[0]-vd->xyz0[0])/vd->dxyz[0]);
			ijk[2]=(int)((first_p[2]-vd->xyz0[2])/vd->dxyz[2]);
		}
		else if(ini_face==3) {
			ijk[1]=vd->nxyz[1]-1;
			ijk[0]=(int)((first_p[0]-vd->xyz0[0])/vd->dxyz[0]);
			ijk[2]=(int)((first_p[2]-vd->xyz0[2])/vd->dxyz[2]);
		}
		else if(ini_face==4) {
			ijk[2]=0;
			ijk[0]=(int)((first_p[0]-vd->xyz0[0])/vd->dxyz[0]);
			ijk[1]=(int)((first_p[1]-vd->xyz0[1])/vd->dxyz[1]);
		}
		else if(ini_face==5) {
			ijk[2]=vd->nxyz[2]-1;
			ijk[0]=(int)((first_p[0]-vd->xyz0[0])/vd->dxyz[0]);
			ijk[1]=(int)((first_p[1]-vd->xyz0[1])/vd->dxyz[1]);
		}


		for(i=0;i<3;i++)
			out_point[i]=0.0;
		con1=find_out_point(vd, first_p,ray_direction,face_sect,out_point);
		if(con1==0) {
			/* start voxel maybe at the neighbour voxel*/
			find_next_cell(ijkn, ijk, face_sect, ray_direction, next_ijk);
			for(i=0;i<3;i++)
				ijk[i]=next_ijk[i];
			ini_index=ijk[2]*vd->nxyz[1]*vd->nxyz[0]+ijk[1]*vd->nxyz[0]+ijk[0];
		}
		search_leave(root_tree, ini_index, first_p, &vvv, ray_direction);
		vvv->local_face_in=ini_face;
		*voxel_p=vvv;
	}
	return(intersection);
}



void build_adjacent_index(int connect[8][6], int local_face[8][6])
{
	connect[0][0]=connect[0][2]=connect[0][4]=-1;
	local_face[0][0]=0; local_face[0][2]=2; local_face[0][4]=4;
	connect[0][1]=1; local_face[0][1]=0;
	connect[0][3]=2;  local_face[0][3]=2;
	connect[0][5]=4;  local_face[0][5]=4;
	connect[1][1]=connect[1][2]=connect[1][4]=-1;
	local_face[1][1]=1; local_face[1][2]=2; local_face[1][4]=4;
	connect[1][0]=0;  local_face[1][0]=1;
	connect[1][3]=3;  local_face[1][3]=2;
	connect[1][5]=5;  local_face[1][5]=4;
	connect[2][0]=connect[2][3]=connect[2][4]=-1;
	local_face[2][0]=0; local_face[2][3]=3; local_face[2][4]=4;
	connect[2][1]=3; local_face[2][1]=0;
	connect[2][2]=0; local_face[2][2]=3;
	connect[2][5]=6; local_face[2][5]=4;
	connect[3][1]=connect[3][3]=connect[3][4]=-1;
	local_face[3][1]=1; local_face[3][3]=3;  local_face[3][4]=4;
	connect[3][0]=2;  local_face[3][0]=1;
	connect[3][2]=1;  local_face[3][2]=3;
	connect[3][5]=7;  local_face[3][5]=4;
	connect[4][0]=connect[4][2]=connect[4][5]=-1;
	local_face[4][0]=0;  local_face[4][2]=2;  local_face[4][5]=5;
	connect[4][1]=5;  local_face[4][1]=0;
	connect[4][3]=6;  local_face[4][3]=2;
	connect[4][4]=0;  local_face[4][4]=5;
	connect[5][1]=connect[5][2]=connect[5][5]=-1;
	local_face[5][1]=1; local_face[5][2]=2;  local_face[5][5]=5;
	connect[5][0]=4; local_face[5][0]=1;
	connect[5][3]=7; local_face[5][3]=2;
	connect[5][4]=1; local_face[5][4]=5;
	connect[6][0]=connect[6][3]=connect[6][5]=-1;
	local_face[6][0]=0;  local_face[6][3]=3;  local_face[6][5]=5;
	connect[6][1]=7;  local_face[6][1]=0;
	connect[6][2]=4;  local_face[6][2]=3;
	connect[6][4]=2;  local_face[6][4]=5;
	connect[7][1]=connect[7][3]=connect[7][5]=-1;
	local_face[7][1]=1;  local_face[7][3]=3;  local_face[7][5]=5;
	connect[7][0]=6;  local_face[7][0]=1;
	connect[7][2]=5;  local_face[7][2]=3;
	connect[7][4]=3;  local_face[7][4]=5;
	return;
}
#endif

#ifdef old_intersection
static void find_coordinates_of_cell(int r_level[3],  double orig_xyz[3], double r_dxyz[3], int ijk[3], double vv[8*3])
{
	vv[0*3]=vv[4*3]=vv[7*3]=vv[3*3]=orig_xyz[0]+ijk[0]*r_dxyz[0];
	vv[1*3]=vv[5*3]=vv[6*3]=vv[2*3]=orig_xyz[0]+(ijk[0]+1)*r_dxyz[0];
	vv[0*3+2]=vv[1*3+2]=vv[2*3+2]=vv[3*3+2]=orig_xyz[2]+ijk[2]*r_dxyz[2];
	vv[4*3+2]=vv[5*3+2]=vv[6*3+2]=vv[7*3+2]=orig_xyz[2]+(ijk[2]+1)*r_dxyz[2];
	vv[0*3+1]=vv[1*3+1]=vv[5*3+1]=vv[4*3+1]=orig_xyz[1]+ijk[1]*r_dxyz[1];
	vv[3*3+1]=vv[2*3+1]=vv[6*3+1]=vv[7*3+1]=orig_xyz[1]+(ijk[1]+1)*r_dxyz[1];
	return;
}

int find_out_point(VR_data *vd, int ijk[3], double in_point[3],double ray_direction[3], int face_sect[3],
		double out_point[3])
/*Calculate the coordiates of exit point*/
{
	double fp[4][3], t,a[6],b[6],c[6],d[6],abc,st[6],mint, vv[24];
	int i,j,m,lc,spm[6],con1,minsp,node[4],face1_sect[3],lm;

	find_coordinates_of_cell(vd, ijk, vv);
	lc=-1; con1=1;
	for(m=0;m<6;m++) {
		transform_face_node(m,node);
		for(i=0;i<4;i++)
			for(j=0;j<3;j++){
				fp[i][j]=vv[node[i]*3+j];
				/*			fv[i][j]=cell->v_data[node[i]*3+j];
				 */
			}
		a[m]=(fp[1][1]-fp[0][1])*(fp[3][2]-fp[0][2])-(fp[3][1]-fp[0][1])*(fp[1][2]-fp[0][2]);
		b[m]= -(fp[1][0]-fp[0][0])*(fp[3][2]-fp[0][2])+(fp[3][0]-fp[0][0])*(fp[1][2]-fp[0][2]);
		c[m]=(fp[1][0]-fp[0][0])*(fp[3][1]-fp[0][1])-(fp[3][0]-fp[0][0])*(fp[1][1]-fp[0][1]);
		abc=sqrt(a[m]*a[m]+b[m]*b[m]+c[m]*c[m]);
		a[m]/=abc;
		b[m]/=abc;
		c[m]/=abc;
		d[m]= -(fp[0][0]*a[m]+fp[0][1]*b[m]+fp[0][2]*c[m]);
	}
	lm= -1;
	for(m=0;m<3;m++)
		face_sect[m]= -1;
	for(m=0;m<6;m++) {

		if(fabs(a[m]*in_point[0]+b[m]*in_point[1]+c[m]*in_point[2]+d[m])<0.0000001)
		{
			lm++;
			face_sect[lm]=m;
		}

	}
	for(m=0;m<6;m++) {
		if((m!=face_sect[0]) && (m!=face_sect[1]) && (m!=face_sect[2]))
		{

			if(fabs(a[m]*ray_direction[0]+b[m]*ray_direction[1]+c[m]*ray_direction[2])>0.000001)
			{

				t=(-d[m]-a[m]*in_point[0]-b[m]*in_point[1]-c[m]*in_point[2])/(a[m]*ray_direction[0]+
						b[m]*ray_direction[1]+c[m]*ray_direction[2]);
				/*    printf("t=%f\n",t);*/
				if(t>0){
					lc++;
					st[lc]=t;
					spm[lc]=m;
				}
			}
		}
	}/*end for*/
	if(lc==-1)
		con1=0;
	else {

		mint=st[0];
		minsp=spm[0];
		for(m=1;m<=lc;m++) {
			if(st[m]<mint) {
				mint=st[m];
				minsp=spm[m];
			}
		}
		/*	*delta_t=mint;
		 */
		for(j=0;j<3;j++)
			out_point[j]=in_point[j]+mint*ray_direction[j];
		/*    transform_face_node(minsp,node);
	get_velocity(node,cell,out_point,out_vv);
		 */
		/*    con1=judge_on_face(out_point, node, cell);
    if(con1==0) {
		printf("*****\n");

    }*/

		lm= -1;
		for(m=0;m<3;m++)
			face1_sect[m]= -1;
		for(m=0;m<6;m++) {

			if(fabs(a[m]*out_point[0]+b[m]*out_point[1]+c[m]*out_point[2]+d[m])<0.000001)
			{
				lm++;
				face1_sect[lm]=m;
			}

		}

		for(m=0;m<3;m++)
			face_sect[m]=face1_sect[m];


	}
	/*      in_voxel.local_face_out=face_sect[0];
	 */

	/*    printf("con1=%d\n",con1);*/
	return(con1);
}

#endif


static int find_out_point(double orig_xyz[3], double r_dxyz[3], int ijk[3], double in_point[3],double ray_direction[3], int face_sect[3],
		double out_point[3])
{
	int i,j, mincomp, intersection;
	double minmax[6], t[3], mint;
	int  i1, j1;
	for(i=0;i<3;i++)
		face_sect[i]=-1;
	for(i=0;i<6;i++) {
		i1= i / 2;
		j1= i % 2;
		minmax[i]=orig_xyz[i1]+(ijk[i1]+j1)*r_dxyz[i1];
	}
#ifdef old_version
for(i=0;i<3;i++) {
	if(fabs(ray_direction[i])<EPSILON)
		t[i]=1.0E+17;
	else {
		if(ray_direction[i]>0)
			t[i]=(minmax[i*2+1]-in_point[i])/ray_direction[i];
		else if(ray_direction[i]<0)
			t[i]=(in_point[i]-minmax[i*2])/(-ray_direction[i]);
	}
}
#endif
for(i=0;i<3;i++) {
	if(ray_direction[i]>0)
		t[i]=(minmax[i*2+1]-in_point[i])/(ray_direction[i]+EPSILON);
	else if(ray_direction[i]<0)
		t[i]=(in_point[i]-minmax[i*2])/(-ray_direction[i]+EPSILON);
}

mint=1.0E+17;
mincomp=-1;
for(i=0;i<3;i++) {
	if((mint>t[i]) && (t[i]>EPSILON)) {
		for(j=0;j<3;j++)
			out_point[j]=in_point[j]+t[i]*ray_direction[j];
		if((out_point[0]>=minmax[0]-EPSILON) && (out_point[0]<=minmax[1]+EPSILON) && (out_point[1]>=minmax[2]-EPSILON)
				&& (out_point[1]<=minmax[3]+EPSILON) && (out_point[2]>=minmax[4]-EPSILON) && (out_point[2]<=minmax[5]+EPSILON)){
			mint=t[i];
			mincomp=i;
		}
	}
}
if((mincomp>=0) && (mincomp<=2)) {
	intersection=1;
	for(i=0;i<3;i++)
		out_point[i]=in_point[i]+mint*ray_direction[i];
}
else {
	intersection=0;
	for(i=0;i<3;i++)
		out_point[i]=in_point[i];
}
for(i=0;i<3;i++) {
	if(fabs(out_point[i]-minmax[i*2])<EPSILON)
		face_sect[i]=i*2;
	else if(fabs(out_point[i]-minmax[i*2+1])<EPSILON)
		face_sect[i]=i*2+1;
}

return (intersection);
}







static void find_next_cell(int ijkn[3], int current_ijk[3], int face_sect[3], double vv[3], int next_ijk[3])
/*get the code of next cell*/
{

	int j,m, flag;
	int face1_sect[3];

	flag=1;
	for(j=0;j<3;j++)
		face1_sect[j]=face_sect[j];
	for(j=0;j<3;j++)
		next_ijk[j]=current_ijk[j];
	for(m=0;m<3;m++){

		if(face_sect[m]!=-1) {

			if(face_sect[m]==0) {
				if(vv[0]<0) {
					next_ijk[0]=current_ijk[0]-1;
					face1_sect[m]=1;
					if(next_ijk[0]<0)
						break;
				}
			}
			else if(face_sect[m]==1) {
				if(vv[0]>0) {
					next_ijk[0]=current_ijk[0]+1;
					face1_sect[m]=0;
					if(next_ijk[0]>ijkn[0]-1)
						break;
				}
			}
			else if(face_sect[m]==2) {
				if(vv[1]<0) {
					next_ijk[1]=current_ijk[1]-1;
					face1_sect[m]=3;
					if(next_ijk[1]<0)
						break;
				}
			}
			else if(face_sect[m]==3) {
				if(vv[1]>0) {
					next_ijk[1]=current_ijk[1]+1;
					face1_sect[m]=2;
					if(next_ijk[1]>ijkn[1]-1)
						break;
				}
			}
			else if(face_sect[m]==4) {
				if(vv[2]<0) {
					next_ijk[2]=current_ijk[2]-1;
					face1_sect[m]=5;
					if(next_ijk[2]<0)
						break;
				}
			}
			else if(face_sect[m]==5) {
				if(vv[2]>0) {
					next_ijk[2]=current_ijk[2]+1;
					face1_sect[m]=4;
					if(next_ijk[2]>ijkn[2]-1)
						break;
				}
			}
		}
	}

	for(m=0;m<3;m++)
		face_sect[m]=face1_sect[m];

	return;
}


static int find_next_point(double orig_xyz[3], double r_dxyz[3], int r_level[3], int current_ijk[3], double in_point[3], double out_point[3], double ray_direction[3],
		int next_ijk[3])
{
	int flag, face_sect[3];

	flag=1;
	flag=find_out_point(orig_xyz, r_dxyz, current_ijk, in_point,ray_direction,face_sect,out_point);
	/*	con1=find_out_point(vd, current_ijk, in_point,ray_direction,face_sect,out_point);
	 */	if(flag==0)
		 HECMW_vis_print_exit("ERROR: HEC-MW-VIS-E2007:There is some problem in finding intersection point");


	 /*	    for(i=0;i<3;i++)
			mid_point[i]=(out_point[i]+in_point[i])/2.0;
		for(i=0;i<3;i++)
	        current_ijk[i]=(int)((mid_point[i]-vd->xyz0[i])/vd->dxyz[i]);
	  */
	 if(flag==1) {
		 find_next_cell(r_level, current_ijk, face_sect, ray_direction, next_ijk);
		 if((next_ijk[0]<0) || (next_ijk[0]>r_level[0]-1) || (next_ijk[1]<0) ||
				 (next_ijk[1]>r_level[1]-1) || (next_ijk[2]<0) || (next_ijk[2]>r_level[2]-1))
			 flag=-1;
	 }
	 return(flag);
}


void ray_trace(int remove_0_display_on, int color_mapping_style, double *interval_point,int transfer_function_style,
		double opa_value, int num_of_features, double *fea_point,
		double view_point_d[3],  int interval_mapping_num, int color_system_type, int num_of_lights,
		double *light_point, double k_ads[3], double orig_xyz[3], double dxyz[3],double r_dxyz[3], int r_level[3], int *empty_flag, double *var, double *grad_var, double first_p[3], int first_ijk[3], double ray_direction[3], double mincolor,
		double maxcolor, double accum_rgba[4], double grad_minmax[2], double feap_minmax[2],
		double feai_minmax[2], double dis_minmax[2], double *opa_table, double tav_length, int time_step, int test_i, int test_j)


/*void ray_trace(Parameter_vr *vr, VR_data *vd, Tree_pointer *root_tree, Tree_pointer *voxel_p,
			   double point_o[3], double ray_direction[3], double mincolor, double maxcolor,
			   int connect[8][6], int local_face[8][6], double accum_rgba[4], double grad_minmax[2],
			   double feap_minmax[2], double feai_minmax[2],
			   double dis_minmax[2], double *opa_table, double tav_length, int time_step)
 */
{
	int i,j;
	double in_point[3], out_point[3];
	int flag, single_empty_flag, value_flag, index, current_ijk[3], next_ijk[3];
	int  print_flag;

	for(i=0;i<4;i++)
		accum_rgba[i]=0.0;
	flag=1;
	for(i=0;i<3;i++)
		in_point[i]=first_p[i];
	for(i=0;i<3;i++)
		current_ijk[i]=first_ijk[i];
	while((accum_rgba[3]<0.99) && (flag==1)) {
		/*
     while((accum_rgba[3]<0.99) && (flag==1) && (accum_rgba[0]<0.99) && (accum_rgba[1]<0.99)
           && (accum_rgba[2]<0.99)) {
		 */
		flag=find_next_point(orig_xyz, r_dxyz, r_level, current_ijk, in_point, out_point, ray_direction, next_ijk);

		single_empty_flag=1;
		value_flag=1;
		/*     for(k=0;k<2;k++)
		 for(j=0;j<2;j++)
			 for(i=0;i<2;i++) {

				 index=(current_ijk[2]+k)*(r_level[0]+1)*(r_level[1]+1)+(current_ijk[1]+j)*(r_level[0]+1)+
					 current_ijk[0]+i;
		 */
		for(i=0;i<8;i++) {
			j= (int) i % 4;
			index=(current_ijk[2]+(int) (i / 4)) * (r_level[0]+1)*(r_level[1]+1)+(current_ijk[1]+ (int) (j /2))*(r_level[0]+1) +
			current_ijk[0]+j % 2;
			if(empty_flag[index]==0)
				single_empty_flag=0;
			/*                        if(fabs(var[index])>EPSILON)
                           value_flag=1;
			 */
		}
		/*     value_flag=1;
	 if(remove_0_display_on==1) {
		 value_flag=0;
       for(i=0;i<8;i++) {
           j
		         if(fabs(var[index])>EPSILON) {
			        value_flag=1;
			        break;
				 }
			 }
	 }

	overlap_flag=1;
	 for(k=0;k<3;k++) {
		 if(fabs(in_point[k]-(vd->xyz0[k]+vd->dxyz[k]))<EPSILON) {
			 overlap_flag=0;
			 break;
		 }
	 }
		 */

		/*
	  if((empty_flag==1) && (value_flag==1))
		 */
		if((single_empty_flag==1) && (flag!=0) && (value_flag==1))
			compute_color_vr(current_ijk, color_mapping_style, interval_point, transfer_function_style,
					opa_value, num_of_features, fea_point,
					view_point_d,  interval_mapping_num, color_system_type, num_of_lights,
					light_point, k_ads,  r_level, orig_xyz, r_dxyz, var, grad_var, accum_rgba, mincolor, maxcolor, grad_minmax,
					feap_minmax, feai_minmax, dis_minmax, opa_table, in_point, out_point, tav_length,
					time_step, print_flag);
		if(flag==1) {
			for(i=0;i<3;i++)
				in_point[i]=out_point[i];
			for(i=0;i<3;i++)
				current_ijk[i]=next_ijk[i];
		}
		/*	 if((vd->empty_flag[current_voxel->cell_id-vd->nxyz[0]*vd->nxyz[1]*vd->nxyz[2]]==1)
         && ((fabs(vd->var[current_voxel->cell_id])>EPSILON) ||
		 (fabs(vd->grad_var[current_voxel->cell_id*3])>EPSILON) ||
		 (fabs(vd->grad_var[current_voxel->cell_id*3+1])>EPSILON) ||
		 (fabs(vd->grad_var[current_voxel->cell_id*3+2])>EPSILON)))
		 */
		/*	  if((vd->empty_flag[current_voxel->cell_id-vd->nxyz[0]*vd->nxyz[1]*vd->nxyz[2]]==1)
		   && (fabs(vd->var[current_voxel->cell_id])>EPSILON))
		 */
	}
	return;
}
