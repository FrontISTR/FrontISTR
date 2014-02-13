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
 *      Simulations Using High End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_vis_connectivity_build.h"

#include <math.h>
#include "hecmw_vis_mem_util.h"
#include "hecmw_vis_comm_util.h"
#include "hecmw_malloc.h"


void find_index_connectivity(struct hecmwST_local_mesh *mesh, int *index_connect)
{
	int  i;

	index_connect[0]=0;
	for(i=0;i<mesh->n_elem;i++) {
		if((mesh->elem_type[i]==341) || (mesh->elem_type[i]==342))
			index_connect[i+1]=index_connect[i]+4;
		else if((mesh->elem_type[i]==351) || (mesh->elem_type[i]==352))
			index_connect[i+1]=index_connect[i]+5;
		else if((mesh->elem_type[i]==361) || (mesh->elem_type[i]==362))
			index_connect[i+1]=index_connect[i]+6;
		else if(mesh->elem_type[i]>400)
			index_connect[i+1]=index_connect[i]+0;
	}
	return;
}

void find_index_a_connect(struct hecmwST_local_mesh *mesh, int num_export, int pe_no, int *export_element, int *index_a_connect)
{
	int i;

	index_a_connect[0]=0;
	for(i=0;i<num_export;i++) {
		if((mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==341) ||
				(mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==342))
			index_a_connect[i+1]=index_a_connect[i]+4;
		else if((mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==351) ||
				(mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==352))
			index_a_connect[i+1]=index_a_connect[i]+5;
		else if((mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==361) ||
				(mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]==362))
			index_a_connect[i+1]=index_a_connect[i]+6;
		else if(mesh->elem_type[export_element[pe_no*mesh->n_elem+i]]>400)
			index_a_connect[i+1]=index_a_connect[i]+0;

	}
	return;
}


void generate_face(int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24])
{
	int i;
	index_face_tetra[0]=0;
	for(i=0;i<4;i++)
		index_face_tetra[i+1]=(i+1)*3;
	face_tetra[0]=0;  face_tetra[1]=2;  face_tetra[2]=1;
	face_tetra[3]=0;  face_tetra[4]=1;  face_tetra[5]=3;
	face_tetra[6]=1;  face_tetra[7]=2;  face_tetra[8]=3;
	face_tetra[9]=2;  face_tetra[10]=0;  face_tetra[11]=3;

	index_face_prism[0]=0;
	for(i=0;i<3;i++)
		index_face_prism[i+1]=(i+1)*4;
	index_face_prism[4]=15;
	index_face_prism[5]=18;
	face_prism[0]=0; face_prism[1]=1; face_prism[2]=4; face_prism[3]=3;
	face_prism[4]=1; face_prism[5]=2; face_prism[6]=5; face_prism[7]=4;
	face_prism[8]=2; face_prism[9]=0; face_prism[10]=3; face_prism[11]=5;
	face_prism[12]=0; face_prism[13]=2; face_prism[14]=1;
	face_prism[15]=3; face_prism[16]=4; face_prism[17]=5;

	index_face_hexa[0]=0;
	for(i=0;i<6;i++)
		index_face_hexa[i+1]=(i+1)*4;
	face_hexa[0]=0;  face_hexa[1]=4;  face_hexa[2]=7;   face_hexa[3]=3;
	face_hexa[4]=1;  face_hexa[5]=2;  face_hexa[6]=6;   face_hexa[7]=5;
	face_hexa[8]=0;  face_hexa[9]=1;  face_hexa[10]=5;   face_hexa[11]=4;
	face_hexa[12]=2;  face_hexa[13]=3;  face_hexa[14]=7;   face_hexa[15]=6;
	face_hexa[16]=3;  face_hexa[17]=2;  face_hexa[18]=1;   face_hexa[19]=0;
	face_hexa[20]=4;  face_hexa[21]=5;  face_hexa[22]=6;   face_hexa[23]=7;
	return;
}


void add_to_hash(int elemID, int faceID, int hashID, Hash_table *h_table)
{
	Hash_table *p1, *p2;

	if((p1=(Hash_table *)HECMW_malloc(sizeof(Hash_table)))==NULL)
		HECMW_vis_memory_exit("hash_table: p1");
	h_table[hashID].elemID++;
	p2=h_table[hashID].next_elem;
	h_table[hashID].next_elem=p1;
	p1->elemID=elemID;
	p1->faceID=faceID;
	p1->next_elem=p2;
	return;
}

void build_hash_table(struct hecmwST_local_mesh *mesh, int *index_connect, Hash_table  *h_table,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24])
{
	int i,j, k, n_elem, node[MAX_N_NODE], nodesum;

	n_elem = mesh->n_elem;
	for(i=0;i<n_elem;i++) {
		for(j=0;j<mesh->elem_node_index[i+1]-mesh->elem_node_index[i];j++)
			node[j]=mesh->elem_node_item[mesh->elem_node_index[i]+j];
		if((mesh->elem_type[i]==341) || (mesh->elem_type[i]==342))  {
			for(j=0;j<4;j++) {
				nodesum=0;
				for(k=index_face_tetra[j];k<index_face_tetra[j+1];k++)
					nodesum+=node[face_tetra[k]];
				add_to_hash(i, j, nodesum, h_table);
			}
		}
		else if((mesh->elem_type[i]==351) || (mesh->elem_type[i]==352)) {
			for(j=0;j<5;j++) {
				nodesum=0;
				for(k=index_face_prism[j];k<index_face_prism[j+1];k++)
					nodesum+=node[face_prism[k]];
				add_to_hash(i, j, nodesum, h_table);
			}
		}
		else if((mesh->elem_type[i]==361) || (mesh->elem_type[i]==362)) {
			for(j=0;j<6;j++) {
				nodesum=0;
				for(k=index_face_hexa[j];k<index_face_hexa[j+1];k++)
					nodesum+=node[face_hexa[k]];
				add_to_hash(i, j, nodesum, h_table);
			}
		}
	}
	return;
}

int is_equal_array(int n[4],int nn[4], int num)
{
	int tmp;
	int i, j;
	int flag=1;
	for(i=0;i<num-1;i++)
		for(j=i+1;j<num;j++){
			if(n[i]<n[j]) {
				tmp=n[i];
				n[i]=n[j];
				n[j]=tmp;
			}
			if(nn[i]<nn[j]) {
				tmp=nn[i];
				nn[i]=nn[j];
				nn[j]=tmp;
			}
		}
	for(i=0;i<num;i++) {
		if(n[i]!=nn[i])
			flag=0;
	}
	return(flag);
}

int is_connect(int elemID1, int faceID1, int elemID2, int faceID2, struct hecmwST_local_mesh *mesh,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],int index_face_hexa[7],
		int face_hexa[24])
{
	int j, node1[MAX_N_NODE],node2[MAX_N_NODE], face1[4], face2[4], k, flag, f_num1, f_num2;
	flag=0;
	for(j=0;j<mesh->elem_node_index[elemID1+1]-mesh->elem_node_index[elemID1];j++)
		node1[j]=mesh->elem_node_item[mesh->elem_node_index[elemID1]+j];
	for(j=0;j<mesh->elem_node_index[elemID2+1]-mesh->elem_node_index[elemID2];j++)
		node2[j]=mesh->elem_node_item[mesh->elem_node_index[elemID2]+j];

	if((mesh->elem_type[elemID1]==341) || (mesh->elem_type[elemID1]==342)){
		f_num1=index_face_tetra[faceID1+1]-index_face_tetra[faceID1];
		for(k=0;k<f_num1;k++)
			face1[k]=node1[face_tetra[index_face_tetra[faceID1]+k]];
	}
	else if((mesh->elem_type[elemID1]==351) || (mesh->elem_type[elemID1]==352)){
		f_num1=index_face_prism[faceID1+1]-index_face_prism[faceID1];
		for(k=0;k<f_num1;k++)
			face1[k]=node1[face_prism[index_face_prism[faceID1]+k]];
	}
	else if((mesh->elem_type[elemID1]==361)
			|| (mesh->elem_type[elemID1]==362)){
		f_num1=index_face_hexa[faceID1+1]-index_face_hexa[faceID1];
		for(k=0;k<f_num1;k++)
			face1[k]=node1[face_hexa[index_face_hexa[faceID1]+k]];
	}
	if((mesh->elem_type[elemID2]==341) || (mesh->elem_type[elemID2]==342)){
		f_num2=index_face_tetra[faceID2+1]-index_face_tetra[faceID2];
		for(k=0;k<f_num2;k++)
			face2[k]=node2[face_tetra[index_face_tetra[faceID2]+k]];
	}
	else if((mesh->elem_type[elemID2]==351) || (mesh->elem_type[elemID2]==352)){
		f_num2=index_face_prism[faceID2+1]-index_face_prism[faceID2];
		for(k=0;k<f_num2;k++)
			face2[k]=node2[face_prism[index_face_prism[faceID2]+k]];
	}
	else if((mesh->elem_type[elemID2]==361) || (mesh->elem_type[elemID2]==541)
			|| (mesh->elem_type[elemID2]==362) || (mesh->elem_type[elemID2]==542)){
		f_num2=index_face_hexa[faceID2+1]-index_face_hexa[faceID2];
		for(k=0;k<f_num2;k++)
			face2[k]=node2[face_hexa[index_face_hexa[faceID2]+k]];
	}
	flag=1;
	if(f_num1!=f_num2)
		flag=0;
	if(flag==1) {
		if(!is_equal_array(face1,face2, f_num1)) {
			flag=0;
		}
	}
	return(flag);

}


int find_to_hash(int elemID1,int faceID1, int hashID,  Hash_table *h_table, struct hecmwST_local_mesh *mesh, int tmp_connect[2],
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],int index_face_hexa[7],
		int face_hexa[24])
{
	int elemID2, faceID2;
	Hash_table *p1;
	int flag;

	flag=0;
	if(h_table[hashID].elemID>0) {
		p1=h_table[hashID].next_elem;
		while(p1!=NULL){
			if(p1->elemID!=elemID1) {
				elemID2=p1->elemID;
				faceID2=p1->faceID;
				if(is_connect(elemID1,faceID1, elemID2, faceID2, mesh,
						index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa, face_hexa)){
					tmp_connect[0]=elemID2;
					tmp_connect[1]=faceID2;
					flag=1;
					return 1;
				}
			}
			p1=p1->next_elem;
		}
	}

	return(flag);
}

void h_free(Hash_table *h_table,int maxadd)
{
	int i;
	Hash_table *p1,*p2;
	for(i=0;i<maxadd;i++) {
		p1=h_table[i].next_elem;
		while(p1!=NULL) {
			p2=p1;
			p1=p1->next_elem;
			HECMW_free(p2);
		}
	}
	HECMW_free(h_table);
	return;
}


void free_b_patch(Boundary_patch *b_patch)
{
	Boundary_patch *p1, *p2;
	p1=b_patch->next_patch;
	while(p1!=NULL) {
		p2=p1->next_patch;
		HECMW_free(p1);
		p1=p2;
	}
	HECMW_free(b_patch);
	return;
}

void  build_connectivity(struct hecmwST_local_mesh *mesh, Hash_table *h_table, int *index_connect, int *connect,
		int index_face_tetra[5], int face_tetra[12], int index_face_prism[6], int face_prism[18],
		int index_face_hexa[7], int face_hexa[24])
{

	int i,j, k, n_elem, node[MAX_N_NODE], nodesum, tmp_connect[2];
	int  flag;

	n_elem = mesh->n_elem;
	for(i=0;i<n_elem;i++) {
		for(j=0;j<mesh->elem_node_index[i+1]-mesh->elem_node_index[i];j++)
			node[j]=mesh->elem_node_item[mesh->elem_node_index[i]+j];
		if((mesh->elem_type[i]==341) || (mesh->elem_type[i]==342)) {
			for(j=0;j<4;j++) {
				nodesum=0;
				for(k=index_face_tetra[j];k<index_face_tetra[j+1];k++) {
					nodesum+=node[face_tetra[k]];
				}
				flag=find_to_hash(i, j, nodesum, h_table, mesh, tmp_connect,
						index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa, face_hexa);
				if(flag==1) {
					connect[(index_connect[i]+j)*3]=tmp_connect[0];
					connect[(index_connect[i]+j)*3+1]=tmp_connect[1];
				}
				else {
					connect[(index_connect[i]+j)*3]=-1;
					connect[(index_connect[i]+j)*3+1]=-1;
				}


			}
		}
		else if((mesh->elem_type[i]==351) || (mesh->elem_type[i]==352)) {
			for(j=0;j<5;j++) {
				nodesum=0;
				for(k=index_face_prism[j];k<index_face_prism[j+1];k++) {
					nodesum+=node[face_prism[k]];
				}
				flag=find_to_hash(i, j, nodesum, h_table, mesh, tmp_connect,
						index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa, face_hexa);
				if(flag==1) {
					connect[(index_connect[i]+j)*3]=tmp_connect[0];
					connect[(index_connect[i]+j)*3+1]=tmp_connect[1];
				}
				else {
					connect[(index_connect[i]+j)*3]=-1;
					connect[(index_connect[i]+j)*3+1]=-1;
				}


			}
		}
		else if((mesh->elem_type[i]==361) || (mesh->elem_type[i]==541) ||
				(mesh->elem_type[i]==362) || (mesh->elem_type[i]==542)) {
			for(j=0;j<6;j++) {
				nodesum=0;
				for(k=index_face_hexa[j];k<index_face_hexa[j+1];k++) {
					nodesum+=node[face_hexa[k]];
				}
				flag=find_to_hash(i, j, nodesum, h_table, mesh, tmp_connect,
						index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa, face_hexa);
				if(flag==1) {
					connect[(index_connect[i]+j)*3]=tmp_connect[0];
					connect[(index_connect[i]+j)*3+1]=tmp_connect[1];
				}
				else {
					connect[(index_connect[i]+j)*3]=-1;
					connect[(index_connect[i]+j)*3+1]=-1;
				}


			}
		}
	}

	return;
}

void add_one_patch(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *node_hit, Tetra_point  *b_point, Head_patch_tetra  *head_b_patch, int node[4], int c_base, int d_base,
		int tn_component) {
	int  k, m;
	int  index_patch,patch[3];
	int  flag;
	double c_data;
	Tetra_point  *p1, *p2;
	Patch_tetra *t1, *t2;
	double tmp;

	index_patch=-1;
	flag=1;
	if((node[0]==node[1]) || (node[0]==node[2]) || (node[1]==node[2]))
		flag=0;
	if(flag==1) {
		for(k=0;k<3;k++) {
			if(node_hit[node[k]-1]==-1) {
				p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
				if(p1==NULL)
					HECMW_vis_memory_exit("p1");
				b_point->ident++;
				p2=b_point->nextpoint;
				b_point->nextpoint=p1;
				if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
					c_data=0.0;
					for(m=0;m<data->nn_dof[sff->color_comp];m++) {
						tmp=data->node_val_item[(node[k]-1)*tn_component+c_base+m];
						c_data+=tmp*tmp;
					}
					c_data=sqrt(c_data);
				}
				else if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp!=0))
					c_data=data->node_val_item[(node[k]-1)*tn_component+c_base+(sff->color_subcomp-1)];
				else if(data->nn_dof[sff->color_comp]==1)
					c_data = data->node_val_item[(node[k]-1)*tn_component+c_base];

				p1->cdata=c_data;
				if(sff->deform_display_on==1) {
					for(m=0;m<3;m++)
						p1->disp[m]=data->node_val_item[(node[k]-1)*tn_component+d_base+m];
				}
				p1->ident=b_point->ident-1;
				p1->geom[0]=mesh->node[(node[k]-1)*3];
				p1->geom[1]=mesh->node[(node[k]-1)*3+1];
				p1->geom[2]=mesh->node[(node[k]-1)*3+2];
				p1->nextpoint=p2;
				index_patch++;
				patch[index_patch]=b_point->ident;
				node_hit[node[k]-1]=b_point->ident;
			}
			else {
				index_patch++;
				patch[index_patch]=node_hit[node[k]-1];
			}
		}
		t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
		if(t1==NULL)
			HECMW_vis_memory_exit("t1");
		t2=head_b_patch->patch_link;
		head_b_patch->num_patch++;
		head_b_patch->patch_link=t1;
		t1->next_patch=t2;
		t1->patch[0]=patch[0];
		t1->patch[1]=patch[1];
		t1->patch[2]=patch[2];
	}
	return;
}


void add_two_patch(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data,
		int *node_hit, Tetra_point  *b_point, Head_patch_tetra  *head_b_patch, int node[4], int c_base, int d_base,
		int tn_component) {
	int  k, m;
	int  index_patch, patch[4];
	int  flag;
	double c_data;
	Tetra_point  *p1, *p2;
	Patch_tetra *t1, *t2;
	double tmp;

	index_patch=-1;
	for(k=0;k<4;k++) {
		if(node_hit[node[k]-1]==-1) {
			p1=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
			if(p1==NULL)
				HECMW_vis_memory_exit("p1");
			b_point->ident++;
			p2=b_point->nextpoint;
			b_point->nextpoint=p1;
			if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp==0)) {
				c_data=0.0;
				for(m=0;m<data->nn_dof[sff->color_comp];m++) {
					tmp=data->node_val_item[(node[k]-1)*tn_component+c_base+m];
					c_data+=tmp*tmp;
				}
				c_data=sqrt(c_data);
			}
			else if((data->nn_dof[sff->color_comp]>1) && (sff->color_subcomp!=0))
				c_data=data->node_val_item[(node[k]-1)*tn_component+c_base+(sff->color_subcomp-1)];
			else if(data->nn_dof[sff->color_comp]==1)
				c_data = data->node_val_item[(node[k]-1)*tn_component+c_base];

			p1->cdata=c_data;
			if(sff->deform_display_on==1) {
				for(m=0;m<3;m++)
					p1->disp[m]=data->node_val_item[(node[k]-1)*tn_component+d_base+m];
			}

			p1->ident=b_point->ident-1;
			p1->geom[0]=mesh->node[(node[k]-1)*3];
			p1->geom[1]=mesh->node[(node[k]-1)*3+1];
			p1->geom[2]=mesh->node[(node[k]-1)*3+2];
			p1->nextpoint=p2;
			index_patch++;
			patch[index_patch]=b_point->ident;
			node_hit[node[k]-1]=b_point->ident;
		}
		else {
			index_patch++;
			patch[index_patch]=node_hit[node[k]-1];
		}
	}
	flag=1;
	if((node[0]==node[1]) || (node[0]==node[2]) || (node[1]==node[2]))
		flag=0;
	if(flag==1) {
		t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
		if(t1==NULL)
			HECMW_vis_memory_exit("t1");
		t2=head_b_patch->patch_link;
		head_b_patch->num_patch++;
		head_b_patch->patch_link=t1;
		t1->next_patch=t2;
		t1->patch[0]=patch[0];
		t1->patch[1]=patch[1];
		t1->patch[2]=patch[2];
	}
	flag=1;
	if((node[0]==node[2]) || (node[0]==node[3]) || (node[2]==node[3]))
		flag=0;
	if(flag==1) {
		t1=(Patch_tetra *)HECMW_malloc(sizeof(Patch_tetra));
		if(t1==NULL)
			HECMW_vis_memory_exit("t1");
		t2=head_b_patch->patch_link;
		head_b_patch->num_patch++;
		head_b_patch->patch_link=t1;
		t1->next_patch=t2;
		t1->patch[0]=patch[0];
		t1->patch[1]=patch[2];
		t1->patch[2]=patch[3];
	}

	return;
}


void HECMW_vis_find_boundary_surface(Surface *sff, struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, int *tvertex, int *tpatch,
		double *minc, double *maxc, Result *result, int sf_i, HECMW_Comm VIS_COMM, int init_flag, Connect_inf *global_connect)

{
	HECMW_Status	stat;
	int pesize, mynode;
	Hash_table  *h_table;
	int i,j,k, m;
	int *export_no_element, *import_no_element, *export_elem, *flag, *index_import_element;
	int *global_id, *among_connect, **g_id, **a_connect, **index_a;
	int flag1, nnode[MAX_N_NODE], pe_no, startnode;
	/*  Boundary_patch *b_patch, *p1, *p2;
	 */
	int node[4];
	int c_base;
	int index_face_tetra[5], face_tetra[12], index_face_prism[6], face_prism[18],index_face_hexa[7],face_hexa[24];
	int *index_a_connect;
	int  n_elem;
	int  *node_hit;
	Tetra_point  *b_point;
	Head_patch_tetra  *head_b_patch;
	Tetra_point  *p1, *p2;
	Patch_tetra *t1, *t2;
	int *nelem_dist;
	int tmp_sum, tmp_nelem;
	int *connect, *index_connect;
	int  tn_component;
	int  d_base;

	HECMW_Comm_size(VIS_COMM, &pesize);
	HECMW_Comm_rank(VIS_COMM, &mynode);
	n_elem=mesh->n_elem;
	tn_component=0;
	for(i=0;i<data->nn_component;i++)
		tn_component+=data->nn_dof[i];
	generate_face(index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa, face_hexa);
	if(init_flag==1) {
		nelem_dist=(int *)HECMW_calloc(pesize+1, sizeof(int));
		if(mynode==0) {
			nelem_dist[0]=0;
			nelem_dist[1]=mesh->ne_internal;
			tmp_sum=mesh->ne_internal;
			for(i=1;i<pesize;i++) {
				HECMW_Recv(&tmp_nelem, 1, HECMW_INT, i, HECMW_ANY_TAG, VIS_COMM, &stat);
				tmp_sum+=tmp_nelem;
				nelem_dist[i+1]=tmp_sum;
			}
			for(i=1;i<pesize;i++)
				HECMW_Send(nelem_dist,pesize+1,HECMW_INT, i, 0, VIS_COMM);
		}
		else {
			HECMW_Send(&mesh->ne_internal, 1, HECMW_INT, 0, 0, VIS_COMM);
			HECMW_Recv(nelem_dist, pesize+1, HECMW_INT, 0, HECMW_ANY_TAG, VIS_COMM, &stat);
		}
		/*	if(mynode==0) {
		for(i=0;i<pesize+1;i++)
			fprintf(stderr, "nelem_dist=%d ", nelem_dist[i]);
		fprintf(stderr, "\n");
	}
		 */
		mesh->global_elem_ID=(int *)HECMW_calloc(mesh->n_elem, sizeof(int));
		if(mesh->global_elem_ID==NULL)
			HECMW_vis_memory_exit("Global_elem_ID");
		for(i=0;i<mesh->n_elem;i++)
			mesh->global_elem_ID[i]=mesh->elem_ID[i*2]+nelem_dist[mesh->elem_ID[i*2+1]];


		HECMW_Barrier(VIS_COMM);
		if((h_table=(Hash_table *)HECMW_calloc(mesh->n_node*4-3, sizeof(Hash_table)))==NULL)
			HECMW_vis_memory_exit("h_table");
		for(i=0;i<mesh->n_node*4-3;i++) {
			h_table[i].elemID=0;
			h_table[i].next_elem=NULL;
		}
		index_connect=(int *)HECMW_calloc(mesh->n_elem+1, sizeof(int));
		if(index_connect==NULL)
			HECMW_vis_memory_exit("index_connect");
		find_index_connectivity(mesh, index_connect);

		/*   if((connect=(int *)HECMW_calloc(index_connect[mesh->n_elem]*3, sizeof(int)))==NULL) */
		if((connect=(int *)HECMW_calloc(index_connect[mesh->n_elem]*3+20, sizeof(int)))==NULL) /* Modyfied S. Ito at 10/3/2007 */
			HECMW_vis_memory_exit("connect");
		/*  i*3: elem_id; i*3+1: face_id   i*3+2: pe_id */
		for(i=0;i<index_connect[mesh->n_elem]*3;i++)
			connect[i]=-1;
		build_hash_table(mesh,index_connect, h_table, index_face_tetra, face_tetra, index_face_prism, face_prism,index_face_hexa,
				face_hexa);

		build_connectivity(mesh,h_table, index_connect, connect, index_face_tetra, face_tetra, index_face_prism,
				face_prism,index_face_hexa, face_hexa);
		h_free(h_table,mesh->n_node*4-3);

		HECMW_Barrier(VIS_COMM);
		for(i=0;i<index_connect[mesh->n_elem];i++) {
			if(connect[i*3]!=-1)
				connect[i*3+2]=mynode;
		}

		if(mesh->n_neighbor_pe>0) {
			if((export_no_element=(int *)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int)))==NULL)
				HECMW_vis_memory_exit("export_no_element");
			if((import_no_element=(int *)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int)))==NULL)
				HECMW_vis_memory_exit("import_no_element");
			if((index_import_element=(int *)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int)))==NULL)
				HECMW_vis_memory_exit("import_no_element");

			for(i=0;i<mesh->n_neighbor_pe;i++)
				export_no_element[i]=0;
			if((export_elem=(int *)HECMW_calloc(mesh->n_neighbor_pe*mesh->n_elem, sizeof(int)))==NULL)
				HECMW_vis_memory_exit("export_elem");
			if((flag=(int *) HECMW_calloc(mesh->n_neighbor_pe, sizeof(int)))==NULL)
				HECMW_vis_memory_exit("flag_connect");

			for(i=0;i<mesh->n_neighbor_pe*mesh->n_elem;i++)
				export_elem[i]=-1;
			for(i=0;i<mesh->n_elem;i++) {
				for(j=0;j<mesh->n_neighbor_pe;j++)
					flag[j]=0;
				for(j=0;j<mesh->elem_node_index[i+1]-mesh->elem_node_index[i];j++) {
					nnode[j]=mesh->elem_node_item[mesh->elem_node_index[i]+j];
					if(nnode[j]>mesh->nn_internal) {

						for(k=0;k<mesh->n_neighbor_pe;k++) {
							if(flag[k]==0) {
								/*			if(k==0) startnode=0;
			else startnode=v->mesh->export_index[k-1];
								 */
								startnode=mesh->import_index[k];
								for(m=startnode;m<mesh->import_index[k+1];m++) {
									if(mesh->import_item[m]==nnode[j])   {
										flag[k]=1;
										export_no_element[k]++;
										export_elem[k*mesh->n_elem+export_no_element[k]-1]=i;
										break;

									}
								}
							}
						}
					}
				}


			}
		} /* end of if n_neighbour_pe >0 */



		for(pe_no=0;pe_no<pesize;pe_no++) {
			if(mynode!=pe_no) {
				for(i=0;i<mesh->n_neighbor_pe;i++) {
					if(mesh->neighbor_pe[i]==pe_no) {

						/* send export data */
						if(export_no_element[i]>0) {
							/*                   fprintf(stderr, "the export number for PE %d in pe %d is %d\n", pe_no, mynode, export_no_element[i]);
							 */
							if((global_id=(int *)HECMW_calloc(export_no_element[i], sizeof(int)))==NULL)
								HECMW_vis_memory_exit("global_id");
							index_a_connect=(int *)HECMW_calloc(export_no_element[i]+1, sizeof(int));
							if(index_a_connect==NULL)
								HECMW_vis_memory_exit("index_a_connect");
							find_index_a_connect(mesh, export_no_element[i], i, export_elem, index_a_connect);
							if((among_connect=(int *)HECMW_calloc(index_a_connect[export_no_element[i]]*3, sizeof(int)))==NULL)
								HECMW_vis_memory_exit("among_connect");
							for(j=0;j<export_no_element[i];j++) {
								global_id[j]=*(export_elem[i*mesh->n_elem+j]+mesh->global_elem_ID);
								for(k=0;k<index_a_connect[j+1]-index_a_connect[j];k++) {
									among_connect[(index_a_connect[j]+k)*3]=connect[(index_connect[export_elem[i*n_elem+j]]+k)*3];
									among_connect[(index_a_connect[j]+k)*3+1]=connect[(index_connect[export_elem[i*n_elem+j]]+k)*3+1];
									among_connect[(index_a_connect[j]+k)*3+2]=connect[(index_connect[export_elem[i*n_elem+j]]+k)*3+2];
								}

							}
						}
						HECMW_Send(&export_no_element[i], 1, HECMW_INT, pe_no, 0, VIS_COMM);
						if(export_no_element[i]>0) {
							HECMW_Send(&index_a_connect[export_no_element[i]], 1, HECMW_INT, pe_no, 0, VIS_COMM);
							HECMW_Send(global_id, export_no_element[i], HECMW_INT, pe_no, 0, VIS_COMM);
							HECMW_Send(index_a_connect, export_no_element[i]+1, HECMW_INT, pe_no, 0, VIS_COMM);
							HECMW_Send(among_connect, index_a_connect[export_no_element[i]]*3, HECMW_INT, pe_no, 0, VIS_COMM);

							HECMW_free(global_id);
							HECMW_free(index_a_connect);
							HECMW_free(among_connect);
						}
					}  /* end of if ==pe_no */
				} /* end of for i */
			} /* end of if mynode !=pe_no */
			/* else if (mynode==pe_no) { */
			else if (mynode==pe_no && pesize > 1) { /* Modified by S. Ito at 10/2/2007 */
				/* receive overlapped element data from other pes */
				if((g_id=(int **)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int *)))==NULL)
					HECMW_vis_memory_exit("g_id");

				if((a_connect=(int **)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int *)))==NULL)
					HECMW_vis_memory_exit("a_connect");
				if((index_a=(int **)HECMW_calloc(mesh->n_neighbor_pe, sizeof(int *)))==NULL)
					HECMW_vis_memory_exit("index_a");


				for(i=0;i<mesh->n_neighbor_pe;i++) {
					j=mesh->neighbor_pe[i];
					HECMW_Recv(&import_no_element[i], 1, HECMW_INT, j, HECMW_ANY_TAG,
							VIS_COMM, &stat);
					if(import_no_element[i]>0) {
						HECMW_Recv(&index_import_element[i], 1, HECMW_INT, j, HECMW_ANY_TAG,
								VIS_COMM, &stat);
						if((g_id[i]=(int *)HECMW_calloc(import_no_element[i], sizeof(int)))==NULL)
							HECMW_vis_memory_exit("g_id");

						HECMW_Recv(g_id[i], import_no_element[i], HECMW_INT, j, HECMW_ANY_TAG,
								VIS_COMM, &stat);
						if((index_a[i]=(int *)HECMW_calloc(import_no_element[i]+1, sizeof(int)))==NULL)
							HECMW_vis_memory_exit("index_a");

						if((a_connect[i]=(int *)HECMW_calloc(index_import_element[i]*3, sizeof(int)))==NULL)
							HECMW_vis_memory_exit("a_connect");
						HECMW_Recv(index_a[i], import_no_element[i]+1, HECMW_INT, j, HECMW_ANY_TAG,
								VIS_COMM, &stat);

						HECMW_Recv(a_connect[i], index_import_element[i]*3, HECMW_INT, j, HECMW_ANY_TAG,
								VIS_COMM, &stat);
					}
				}
				/* recompute boundary element */
				for(i=0;i<mesh->n_elem;i++) {
					for(j=0;j<index_connect[i+1]-index_connect[i];j++) {
						if(connect[(index_connect[i]+j)*3]==-1) {
							flag1=1;
							for(k=0;k<mesh->n_neighbor_pe;k++) {
								if(flag1==1) {
									for(m=0;m<import_no_element[k];m++) {
										if((*(g_id[k]+m)==mesh->global_elem_ID[i]) &&
												(*(a_connect[k]+(*(index_a[k]+m)+j)*3)!=-1)) {
											connect[(index_connect[i]+j)*3]=*(a_connect[k]+(*(index_a[k]+m)+j)*3);
											connect[(index_connect[i]+j)*3+1]=*(a_connect[k]+(*(index_a[k]+m)+j)*3+1);
											connect[(index_connect[i]+j)*3+2]=*(a_connect[k]+(*(index_a[k]+m)+j)*3+2);

											flag1=0;
											break;

										}
									}
								}
							}
						}
					}
				}
				HECMW_free(a_connect);
				HECMW_free(g_id);
				HECMW_free(index_a);
			} /* end of if mynode==pe_no */
			/* finish in one pe_no */
			HECMW_Barrier(VIS_COMM);
		} /* end of for pe_no */
		/*	  fprintf(stderr,"Finish the connectivity build among PES\n");
		 */
	}
	/*  b_patch=(Boundary_patch *)HECMW_malloc(sizeof(Boundary_patch));
  if(b_patch==NULL)
	  HECMW_vis_memory_exit("b_patch");
  p1=b_patch;
	 */
	if(init_flag!=1) {
		index_connect=global_connect->index_connect;
		connect=global_connect->connect;
	}

	b_point=(Tetra_point *)HECMW_malloc(sizeof(Tetra_point));
	head_b_patch=(Head_patch_tetra *)HECMW_malloc(sizeof(Head_patch_tetra));
	if((b_point==NULL) || (head_b_patch==NULL))
		HECMW_vis_memory_exit("boundary: initialization");
	b_point->ident=0;
	head_b_patch->num_patch=0;

	node_hit=(int *)HECMW_calloc(mesh->n_node, sizeof(int));
	if(node_hit==NULL)
		HECMW_vis_memory_exit("node_hit");
	for(i=0;i<mesh->n_node;i++)
		node_hit[i]=-1;
	/*  if(sff->deform_display_on==1) {
	  flag_deform=0;
	  for(i=0;i<data->nn_component;i++) {
        if(strncmp("DISPLACEMENT", data->node_label[i], 10)==0) {
			flag_deform=1;
			sff->disp_comp=i;
			break;
		}
	  }
	  if(flag_deform==0)
		  HECMW_vis_print_exit("No displament data for displaying deform");
  }
	 */

	n_elem=mesh->n_elem;
	c_base = 0;
	for (i = 0; i < sff->color_comp; i++) {
		c_base +=data->nn_dof[i];
	}
	d_base = 0;
	if(sff->deform_display_on==1) {
		for(i=0;i<sff->disp_comp;i++)
			d_base +=data->nn_dof[i];
	}
	*minc=1.0E17;
	*maxc=-1.0E17;
	for(i=0;i<n_elem;i++)
		if(mesh->elem_ID[i*2+1]==mynode) {
			for(j=0;j<index_connect[i+1]-index_connect[i];j++) {
				if(connect[(index_connect[i]+j)*3]==-1) {
					if((mesh->elem_type[i]==341) || (mesh->elem_type[i]==342)){
						for(k=0;k<index_face_tetra[j+1]-index_face_tetra[j];k++)
							node[k]=mesh->elem_node_item[mesh->elem_node_index[i]+face_tetra[index_face_tetra[j]+k]];
						add_one_patch(sff, mesh, data, node_hit, b_point, head_b_patch, node, c_base, d_base, tn_component);
					}
					else if((mesh->elem_type[i]==351) || (mesh->elem_type[i]==352)) {
						for(k=0;k<index_face_prism[j+1]-index_face_prism[j];k++)
							node[k]=mesh->elem_node_item[mesh->elem_node_index[i]+face_prism[index_face_prism[j]+k]];
						if(j<3)
							add_two_patch(sff, mesh, data, node_hit, b_point, head_b_patch, node, c_base, d_base, tn_component);
						else if(j>=3)
							add_one_patch(sff, mesh, data, node_hit, b_point, head_b_patch, node, c_base, d_base, tn_component);
					}
					else if((mesh->elem_type[i]==361) || (mesh->elem_type[i]==362)) {
						for(k=0;k<index_face_hexa[j+1]-index_face_hexa[j];k++)
							node[k]=mesh->elem_node_item[mesh->elem_node_index[i]+face_hexa[index_face_hexa[j]+k]];
						add_two_patch(sff, mesh, data, node_hit, b_point, head_b_patch, node, c_base, d_base, tn_component);
					}

				}
			}
		}
	HECMW_free(node_hit);
	if(head_b_patch->num_patch>0) {
		result[sf_i].n_vertex=b_point->ident;
		result[sf_i].n_patch=head_b_patch->num_patch;
		result[sf_i].vertex=(double *)HECMW_calloc(result[sf_i].n_vertex*3,sizeof(double));
		result[sf_i].color=(double *)HECMW_calloc(result[sf_i].n_vertex, sizeof(double));
		result[sf_i].patch=(int *)HECMW_calloc(result[sf_i].n_patch*3, sizeof(int));
		if((result[sf_i].vertex==NULL) || (result[sf_i].patch==NULL) || (result[sf_i].color==NULL))
			HECMW_vis_memory_exit("Boundary: vertex, patch and color");
		if(sff->deform_display_on==1) {
			result[sf_i].disp=(double *)HECMW_calloc(result[sf_i].n_vertex*3,sizeof(double));
			if(result[sf_i].disp==NULL)
				HECMW_vis_memory_exit("Boundary: disp");
		}
	}
	if(b_point->ident>0) {
		p1=b_point->nextpoint;
		for(i=0;i<b_point->ident;i++) {
			for(j=0;j<3;j++)
				result[sf_i].vertex[(p1->ident)*3+j]=p1->geom[j];
			result[sf_i].color[p1->ident]=p1->cdata;
			if(p1->cdata<*minc) *minc=p1->cdata;
			if(p1->cdata>*maxc) *maxc=p1->cdata;
			if(sff->deform_display_on==1) {
				for(j=0;j<3;j++)
					result[sf_i].disp[(p1->ident)*3+j]=p1->disp[j];
			}

			p2=p1;
			p1=p1->nextpoint;
			HECMW_free(p2);
		}
		HECMW_free(b_point);
	}
	if(head_b_patch->num_patch>0) {
		t1=head_b_patch->patch_link;
		for(i=0;i<head_b_patch->num_patch;i++) {
			for(j=0;j<3;j++)
				result[sf_i].patch[i*3+j]=*tvertex+t1->patch[j];
			t2=t1;
			t1=t1->next_patch;
			HECMW_free(t2);
		}
		HECMW_free(head_b_patch);
	}

	*tvertex+=result[sf_i].n_vertex;
	*tpatch+=result[sf_i].n_patch;
	/*          HECMW_free(index_connect);
          HECMW_free(connect);
	 */
	if(init_flag==1) {
		global_connect->index_connect=index_connect;
		global_connect->connect=connect;
	}
	return;
}


