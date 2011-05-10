/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Dynamic Load Balancing                            *
 *                                                                     *
 *            Written by Li Chen (Univ. of Tokyo)                      *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/

#include "hecmw_repart.h"
extern struct hecmwST_local_mesh *mesh;
extern struct hecmwST_result_data *data;
extern struct hecmwST_local_mesh *new_mesh;
extern struct hecmwST_result_data *new_data;

void dist_dlb_free_section(struct hecmwST_section *section)
{
	if(section == NULL) return;
    if(section->sect_type!=NULL)
	free(section->sect_type);
	if(section->sect_opt!=NULL)
	free(section->sect_opt);
	if(section->sect_mat_ID_index!=NULL)
	free(section->sect_mat_ID_index);
	if(section->sect_mat_ID_item!=NULL)
	free(section->sect_mat_ID_item);
	if(section->sect_I_index!=NULL)
	free(section->sect_I_index);
	if(section->sect_I_item!=NULL)
	free(section->sect_I_item);
	if(section->sect_R_index!=NULL)
	free(section->sect_R_index);
	if(section->sect_R_item!=NULL)
	free(section->sect_R_item);
	return;
}


void
dist_dlb_free_material(struct hecmwST_material *material)
{
	int i;

	if(material == NULL) return;
	if(material->mat_name!=NULL) {
	for(i=0; i < material->n_mat; i++) {
		free(material->mat_name[i]);
	}
	}
	free(material->mat_name);
	if(material->mat_item_index!=NULL)
	free(material->mat_item_index);
	if(material->mat_subitem_index!=NULL)
	free(material->mat_subitem_index);
	if(material->mat_table_index!=NULL)
	free(material->mat_table_index);
	if(material->mat_val!=NULL)
	free(material->mat_val);
	if(material->mat_temp!=NULL)
	free(material->mat_temp);
	return;
}


void
dist_dlb_free_mpc(struct hecmwST_mpc *mpc)
{
	if(mpc == NULL) return;
	if(mpc->mpc_index!=NULL)
	free(mpc->mpc_index);
	if(mpc->mpc_item!=NULL)
	free(mpc->mpc_item);
	if(mpc->mpc_dof!=NULL)
	free(mpc->mpc_dof);
	if(mpc->mpc_val!=NULL)
	free(mpc->mpc_val);
}


void
dist_dlb_free_amplitude(struct hecmwST_amplitude *amp)
{
	int i;

	if(amp == NULL) return;
    if(amp->amp_name!=NULL) {
	for(i=0; i < amp->n_amp; i++) {
		free(amp->amp_name[i]);
	}
	free(amp->amp_name);
	}
	if(amp->amp_type_definition!=NULL)
	free(amp->amp_type_definition);
	if(amp->amp_type_time!=NULL)
	free(amp->amp_type_time);
	if(amp->amp_type_value!=NULL)
	free(amp->amp_type_value);
	if(amp->amp_index!=NULL)
	free(amp->amp_index);
	if(amp->amp_val!=NULL)
	free(amp->amp_val);
	if(amp->amp_table!=NULL)
	free(amp->amp_table);
}

void
dist_dlb_free_ngrp(struct hecmwST_node_grp *grp)
{
	int i;

	if(grp == NULL) return;
    if(grp->grp_name!=NULL) {
	for(i=0; i < grp->n_grp; i++) {
		free(grp->grp_name[i]);
	}
	free(grp->grp_name);
	}
	if(grp->grp_index!=NULL)
	free(grp->grp_index);
	if(grp->grp_item!=NULL)
	free(grp->grp_item);
	if(grp->bc_grp_ID!=NULL)
	free(grp->bc_grp_ID);
	if(grp->bc_grp_type!=NULL)
	free(grp->bc_grp_type);
	if(grp->bc_grp_index!=NULL)
	free(grp->bc_grp_index);
	if(grp->bc_grp_dof!=NULL)
	free(grp->bc_grp_dof);
	if(grp->bc_grp_val!=NULL)
	free(grp->bc_grp_val);
	return;
}


void
dist_dlb_free_egrp(struct hecmwST_elem_grp *grp)
{
	int i;

	if(grp == NULL) return;
    if(grp->grp_name!=NULL) {
	for(i=0; i < grp->n_grp; i++) {
		free(grp->grp_name[i]);
	}
	free(grp->grp_name);
	}
	if(grp->grp_index!=NULL)
	free(grp->grp_index);
	if(grp->grp_item!=NULL)
	free(grp->grp_item);
	if(grp->bc_grp_ID!=NULL)
	free(grp->bc_grp_ID);
	if(grp->bc_grp_type!=NULL)
	free(grp->bc_grp_type);
	if(grp->bc_grp_index!=NULL)
	free(grp->bc_grp_index);
	if(grp->bc_grp_val!=NULL)
	free(grp->bc_grp_val);
	return;
}


void
dist_dlb_free_sgrp(struct hecmwST_surf_grp *grp)
{
	int i;

	if(grp == NULL) return;
    if(grp->grp_name!=NULL) {
	for(i=0; i < grp->n_grp; i++) {
		HECMW_free(grp->grp_name[i]);
	}
	free(grp->grp_name);
	}
	if(grp->grp_index!=NULL)
	free(grp->grp_index);
	if(grp->grp_item!=NULL)
	free(grp->grp_item);
	if(grp->bc_grp_ID!=NULL)
	free(grp->bc_grp_ID);
	if(grp->bc_grp_type!=NULL)
	free(grp->bc_grp_type);
	if(grp->bc_grp_index!=NULL)
	free(grp->bc_grp_index);
	if(grp->bc_grp_val!=NULL)
	free(grp->bc_grp_val);
	return;
}


void hecmw_dlb_f2c_init_()
{
   mesh=(struct hecmwST_local_mesh *)malloc(sizeof(struct hecmwST_local_mesh));
	if(mesh==NULL)
	   HECMW_dlb_memory_exit("mesh");
   mesh->section=(struct hecmwST_section *)malloc(sizeof(struct hecmwST_section));
   mesh->material=(struct hecmwST_material *)malloc(sizeof(struct hecmwST_material));
   mesh->mpc=(struct hecmwST_mpc *)malloc(sizeof(struct hecmwST_mpc));
   mesh->amp=(struct hecmwST_amplitude *)malloc(sizeof(struct hecmwST_amplitude));
   mesh->node_group=(struct hecmwST_node_grp *)malloc(sizeof(struct hecmwST_node_grp));
   mesh->elem_group=(struct hecmwST_elem_grp *)malloc(sizeof(struct hecmwST_elem_grp));
   mesh->surf_group=(struct hecmwST_surf_grp *)malloc(sizeof(struct hecmwST_surf_grp));
   HECMW_dist_copy_f2c_init(mesh);
   return;
}

void hecmw_dlb_c2f_init_()
{
	HECMW_dist_copy_c2f_init(new_mesh);
	return;
}

void hecmw_dlb_f2c_finalize_()
{
   HECMW_dist_copy_f2c_finalize();
   return;
}

void hecmw_dlb_c2f_finalize_()
{
	int i;

	if(new_mesh == NULL) return;
    if(new_mesh->files!=NULL) {
	for(i=0; i < new_mesh->hecmw_n_file; i++) {
		free(new_mesh->files[i]);
	}
	free(new_mesh->files);
	}
    if(new_mesh->node_internal_list!=NULL)
	free(new_mesh->node_internal_list);
	free(new_mesh->node_ID);
	if(new_mesh->global_node_ID!=NULL)
	free(new_mesh->global_node_ID);
	if(new_mesh->node!=NULL)
	free(new_mesh->node);
    if(new_mesh->node_dof_index!=NULL)
	free(new_mesh->node_dof_index);
        if(new_mesh->node_dof_item!=NULL)
	free(new_mesh->node_dof_item);
        if(new_mesh->node_val_index!=NULL)
	free(new_mesh->node_val_index);
        if(new_mesh->node_val_item!=NULL)
	free(new_mesh->node_val_item);
        if(new_mesh->node_init_val_index!=NULL)
	free(new_mesh->node_init_val_index);
        if(new_mesh->node_init_val_item!=NULL)
	free(new_mesh->node_init_val_item);
	if(new_mesh->elem_internal_list!=NULL)
	free(new_mesh->elem_internal_list);
	if(new_mesh->elem_ID!=NULL)
	free(new_mesh->elem_ID);
	if(new_mesh->global_elem_ID!=NULL)
	free(new_mesh->global_elem_ID);
	if(new_mesh->elem_type!=NULL)
	free(new_mesh->elem_type);
	if(new_mesh->elem_type_index!=NULL)
	free(new_mesh->elem_type_index);
	if(new_mesh->elem_type_item!=NULL)
	free(new_mesh->elem_type_item);
	if(new_mesh->elem_node_index!=NULL)
	free(new_mesh->elem_node_index);
	if(new_mesh->elem_node_item!=NULL)
	free(new_mesh->elem_node_item);
	if(new_mesh->section_ID!=NULL)
	free(new_mesh->section_ID);
	if(new_mesh->elem_mat_ID_index!=NULL)
	free(new_mesh->elem_mat_ID_index);
	if(new_mesh->elem_mat_ID_item!=NULL)
	free(new_mesh->elem_mat_ID_item);
	if(new_mesh->elem_mat_int_index!=NULL)
	free(new_mesh->elem_mat_int_index);
	if(new_mesh->elem_mat_int_val!=NULL)
	free(new_mesh->elem_mat_int_val);
	if(new_mesh->elem_val_index!=NULL)
	free(new_mesh->elem_val_index);
	if(new_mesh->elem_val_item!=NULL)
	free(new_mesh->elem_val_item);
	if(new_mesh->neighbor_pe!=NULL)
	free(new_mesh->neighbor_pe);
	if(new_mesh->import_index!=NULL)
	free(new_mesh->import_index);
	if(new_mesh->import_item!=NULL)
	free(new_mesh->import_item);
	if(new_mesh->export_index!=NULL)
	free(new_mesh->export_index);
	if(new_mesh->import_item!=NULL)
	free(new_mesh->export_item);
	if(new_mesh->shared_index!=NULL)
	free(new_mesh->shared_index);
	if(new_mesh->import_item!=NULL)
	free(new_mesh->shared_item);

	if(new_mesh->when_i_was_refined_node!=NULL)
	free(new_mesh->when_i_was_refined_node);
	if(new_mesh->when_i_was_refined_elem!=NULL)
	free(new_mesh->when_i_was_refined_elem);
	if(new_mesh->adapt_parent_type!=NULL)
	free(new_mesh->adapt_parent_type);
	if(new_mesh->adapt_type!=NULL)
	free(new_mesh->adapt_type);
	if(new_mesh->adapt_level!=NULL)
	free(new_mesh->adapt_level);
	if(new_mesh->adapt_parent!=NULL)
	free(new_mesh->adapt_parent);
	if(new_mesh->adapt_children_item!=NULL)
	free(new_mesh->adapt_children_item);
	if(new_mesh->section!=NULL){
	dist_dlb_free_section(new_mesh->section);
    free(new_mesh->section);
	}
	

	if(new_mesh->material!=NULL) {
	dist_dlb_free_material(new_mesh->material);
	free(new_mesh->material);
	}
    if(new_mesh->mpc!=NULL){
	dist_dlb_free_mpc(new_mesh->mpc);
	free(new_mesh->mpc);
	}
	if(new_mesh->amp!=NULL) {
	dist_dlb_free_amplitude(new_mesh->amp);
	free(new_mesh->amp);
	}
    if(new_mesh->node_group!=NULL) {
	dist_dlb_free_ngrp(new_mesh->node_group);
	free(new_mesh->node_group);
	}
	if(new_mesh->elem_group!=NULL) {
	dist_dlb_free_egrp(new_mesh->elem_group);
	free(new_mesh->elem_group);
	}
	if(new_mesh->surf_group!=NULL) {
	dist_dlb_free_sgrp(new_mesh->surf_group);
	free(new_mesh->surf_group);
	}
    free(new_mesh);



	HECMW_dist_copy_c2f_finalize();

	return;
}

void test_mesh_()
{
	int mynode, pesize, i;
	FILE  *test_fp;
	char test_file[128];

	HECMW_Comm_rank(mesh->HECMW_COMM, &mynode);
  HECMW_Comm_size(mesh->HECMW_COMM, &pesize);
/*    sprintf(test_file, "test1.%d", mynode);
	test_fp=fopen(test_file, "w");
	for(i=0;i<mesh->import_index[mesh->n_neighbor_pe];i++)
		fprintf(test_fp, "%d\n", mesh->import_item[i]);
	fprintf(test_fp, "export node*****\n");
	for(i=0;i<mesh->export_index[mesh->n_neighbor_pe];i++)
		fprintf(test_fp, "%d\n", mesh->export_item[i]);
    fclose(test_fp);
	
	fprintf(stderr, "mesh: n_node=%d\n", mesh->n_node);
	fprintf(stderr, "mesh: n_adapt=%d\n", mesh->n_adapt);
	*/
	
	return;
}

void hecmw_dist_get_result_c_()
{
	  data=HECMW_result_read_by_fname("result-in");
  if(data==NULL) {
	  HECMW_abort(HECMW_comm_get_comm());
  }
	HECMW_DEBUG(("hecmw_get_restart OK"));
	return;
}


void set_label_name(int n_component, char *label_name_f, char **label_name)
{
   int i,j, k;
   char str_tmp[128];
      j = 0;
    for (i=0; i<n_component; i++) {
      str_tmp[0] = '\0';
      while (label_name_f[j] == ' ') 
	    j++;
      k = 0;
      while (label_name_f[j] != ' ') {
    	str_tmp[k] = label_name_f[j];
    	k++;
    	j++;
      }

      str_tmp[k] = '\0';
      strcpy(label_name[i], str_tmp);
    }
	return;
}

void hecmw_set_result_node_(int *nn_component, int *nn_dof, char *node_label,double *node_val_item)
{
	int  i;
	if(data==NULL) {
       data=(struct hecmwST_result_data *)malloc(sizeof(struct hecmwST_result_data));
       if(data==NULL)
	      HECMW_dlb_memory_exit("data");
	   data->nn_component=0;
	   data->ne_component=0;
	   }
     data->nn_component=*nn_component;
     data->nn_dof=nn_dof;
		data->node_label=(char **)calloc(data->nn_component, sizeof(char *));
		for(i=0;i<data->nn_component;i++)
			data->node_label[i]=(char *)calloc(128, sizeof(char));
		set_label_name(data->nn_component, node_label, data->node_label);
		data->node_val_item=node_val_item;
  return;
}

hecmw_dlb_get_result_node_(double *node_val_item)
{
	node_val_item=new_data->node_val_item;
	return;
}


hecmw_dlb_get_result_elem_(double *elem_val_item)
{
	elem_val_item=new_data->elem_val_item;
	return;
}


void hecmw_set_result_elem_(int *ne_component, int *ne_dof, char *elem_label,double *elem_val_item)
{
	int i;
	if(data==NULL) {
       data=(struct hecmwST_result_data *)malloc(sizeof(struct hecmwST_result_data));
       if(data==NULL)
	      HECMW_dlb_memory_exit("data");
	   data->nn_component=0;
	   data->ne_component=0;
	   }
     data->ne_component=*ne_component;
     data->ne_dof=ne_dof;
		data->elem_label=(char **)calloc(data->ne_component, sizeof(char *));
		for(i=0;i<data->ne_component;i++)
			data->elem_label[i]=(char *)calloc(128, sizeof(char));
		set_label_name(data->ne_component, elem_label, data->elem_label);
		data->elem_val_item=elem_val_item;
  return;
}




/*
extern int l_node, l_elem;
extern int repart_comm;

void hecmw_add_inf_(int *allNODTOTold, int *ICELTOTold, int *SOLVER_COMM)
{
   l_node = *allNODTOTold*9+1;
   l_elem = *ICELTOTold*9+1;
   repart_comm = *SOLVER_COMM;
   return;
}

void hecmw_set_mesh_node_(int *n_node, int *n_internal,  double *node,  int *node_id)
{
   v.mesh=(struct local_mesh *)malloc(sizeof(struct local_mesh));
   v.mesh->n_node=*n_node;
   v.mesh->n_internal= *n_internal;
   v.mesh->node = node;
   v.mesh->node_id = node_id;

   return;
}

void hecmw_set_mesh_element_(int *n_elem,  int *ne_internal, int *elem_type, 
                              int *index_elem,    int *ptr_elem,  int *elem_id,  int *ne_internal_list)
{
	v.mesh->n_elem = *n_elem;
	v.mesh->ne_internal = *ne_internal;
	v.mesh->elem_type = elem_type;
	v.mesh->index_elem = index_elem;
	v.mesh->ptr_elem = ptr_elem;
	v.mesh->elem_id = elem_id;
	v.mesh->ne_internal_list = ne_internal_list;
	return;
}

void hecmw_set_mesh_pe_inf_(int *n_neighbor_pe,  int *neighbor_pe, int *import_index, int *import_node,
                              int *export_index, int *export_node)
{
	v.mesh->n_neighbor_pe = *n_neighbor_pe;
	v.mesh->neighbor_pe = neighbor_pe;
	v.mesh->import_index = import_index;
	v.mesh->import_node = import_node;
	v.mesh->export_index = export_index;
	v.mesh->export_node = export_node;
	return;
}

void hecmw_set_mesh_adaptation_(int *CoarseGridLevels, int *HOWmanyADAPTATIONs, 
                           int *wheniwasrefined_node, int *wheniwasrefined_elem, 
                           int *adaptation_parent_type, int *adaptation_type, 
						   int *adaptation_level, int *adaptation_parent,
                           int *adaptation_children, int *index_children)
{
	v.mesh->coarsegridlevels = *CoarseGridLevels;
	v.mesh->howmanyadaptions = *HOWmanyADAPTATIONs;
	v.mesh->wheniwasrefined_node = wheniwasrefined_node;
	v.mesh->wheniwasrefined_elem = wheniwasrefined_elem;
	v.mesh->adaptation_parent_type = adaptation_parent_type;
	v.mesh->adaptation_type = adaptation_type;
	v.mesh->adaptation_level = adaptation_level;
	v.mesh->adaptation_parent = adaptation_parent;
	v.mesh->adaptation_children = adaptation_children;
	v.mesh->index_children = index_children;
	return;
}


void set_grp_name(int n_group, char *grp_name_f, char **grp_name)
{
   int i,j, k;
   char str_tmp[128];
      j = 0;
    for (i=0; i<n_group; i++) {
      str_tmp[0] = '\0';
      while (grp_name_f[j] == ' ') 
	    j++;
      k = 0;
      while (grp_name_f[j] != ' ') {
    	str_tmp[k] = grp_name_f[j];
    	k++;
    	j++;
      }

      str_tmp[k] = '\0';
      strcpy(grp_name[i], str_tmp);
    }
	return;
}

void hecmw_set_grp_info_node_(int *n_grp,int *grp_index, char *grp_name_f, int *grp_node) 
{
    int i, j;
	
	v.grp=(struct grp_data *)malloc(sizeof(struct grp_data));
	v.grp->node_grp.n_enum_grp=*n_grp;
	v.grp->node_grp.enum_grp_index=grp_index;
	v.grp->node_grp.enum_grp_name=(char **)calloc(*n_grp, sizeof(char *));
	for(i=0;i<*n_grp;i++)
		v.grp->node_grp.enum_grp_name[i]=(char *)calloc(128, sizeof(char));
	set_grp_name(*n_grp, grp_name_f, v.grp->node_grp.enum_grp_name);
	v.grp->node_grp.enum_grp_node=grp_node;
	return;
}

void hecmw_set_grp_info_elem_(int *n_grp,int *grp_index, char *grp_name_f, int *grp_node) 
{
    int i, j;
	
	v.grp->elem_grp.n_enum_grp=*n_grp;
	v.grp->elem_grp.enum_grp_index=grp_index;
	v.grp->elem_grp.enum_grp_name=(char **)calloc(*n_grp, sizeof(char *));
	for(i=0;i<*n_grp;i++)
		v.grp->elem_grp.enum_grp_name[i]=(char *)calloc(128, sizeof(char));
	set_grp_name(*n_grp, grp_name_f, v.grp->elem_grp.enum_grp_name);
	v.grp->elem_grp.enum_grp_node=grp_node;
	return;
}

void hecmw_set_grp_info_surf_(int *n_grp,int *grp_index, char *grp_name_f, int *grp_node) 
{
    int i, j;
	
	v.grp->n_surf_grp=*n_grp;
	v.grp->surf_grp_index=grp_index;
	v.grp->surf_grp_name=(char **)calloc(*n_grp, sizeof(char *));
	for(i=0;i<*n_grp;i++)
		v.grp->surf_grp_name[i]=(char *)calloc(128, sizeof(char));
	set_grp_name(*n_grp, grp_name_f, v.grp->surf_grp_name);
	v.grp->surf_grp_node=grp_node;
	return;
}


void hecmw_set_data_(int *ivar, int *i_free, double *U, double *P, double *VISCT, double *VISCL) 
{
	int i,j,k;
	double u1, v1, w1;
	int tmp_count;
	fprintf(stderr, "It is ok here \n");
	v.node=(struct node_elem_data *)malloc(sizeof(struct node_elem_data));
	v.node->n_component=*ivar;
	v.node->n_free=(int *)calloc(*ivar, sizeof(int));
	for(i=0;i<*ivar;i++)
		v.node->n_free[i]=i_free[i];
	v.node->t_component=0;
	for(i=0;i<v.node->n_component;i++) 
		v.node->t_component+=v.node->n_free[i];
	v.node->data=(double *)calloc(v.node->t_component*v.mesh->n_node,sizeof(double));
	if(v.node->data==NULL) {
		fprintf(stderr, "There is no enough memory for v.node->data\n");
		exit(0);
	}


	for(j=0;j<5;j++) {
		for(i=0;i<v.mesh->n_node;i++)
			v.node->data[j*v.mesh->n_node+i]=U[j*v.mesh->n_node+i];
	}

	for(i=0;i<v.mesh->n_node;i++)
		v.node->data[5*v.mesh->n_node+i]=P[i];
	for(i=0;i<v.mesh->n_node;i++)
		v.node->data[6*v.mesh->n_node+i]=VISCT[i];
	for(i=0;i<v.mesh->n_node;i++)
		v.node->data[7*v.mesh->n_node+i]=VISCL[i];
    
			
	return;
}

void hecmw_num_c2f90_(int *n_node, int *n_internal,int *n_elem, int *ne_internal, int *n_neighbor_pe, int *l_child, 
					  int *l_ptr_elem, int *ivar) 
{
	*n_node=new_mesh->n_node;
	*n_internal=new_mesh->n_internal;
	*n_elem=new_mesh->n_elem;
	*ne_internal=new_mesh->ne_internal;
	*n_neighbor_pe=new_mesh->n_neighbor_pe;
	*l_child=new_mesh->index_children[new_mesh->n_elem];
	*l_ptr_elem=new_mesh->index_elem[new_mesh->n_elem];
	*ivar=v.node->n_component;
	return;
}

void hecmw_mesh_node_c2f90_(double *node, int *node_id) 
{
	int i,j;
	for(i=0;i<new_mesh->n_node*3;i++)
	   node[i]=new_mesh->node[i];
	for(i=0;i<new_mesh->n_node*2;i++) 
	   node_id[i]=new_mesh->node_id[i];
	free(new_mesh->node);
	free(new_mesh->node_id);
	
	
	return;
}

void hecmw_mesh_elem_c2f90_(int *elem_type, int *index_elem, int *ptr_elem,int *elem_id, int *ne_internal_list) 
{
	int i,j;

	for(i=0;i<new_mesh->n_elem;i++)
		elem_type[i]=new_mesh->elem_type[i];
	for(i=0;i<new_mesh->n_elem+1;i++)
		index_elem[i]=new_mesh->index_elem[i];
	for(i=0;i<new_mesh->index_elem[new_mesh->n_elem];i++)
		ptr_elem[i]=new_mesh->ptr_elem[i];
	for(i=0;i<new_mesh->n_elem*2;i++) 
		elem_id[i]=new_mesh->elem_id[i];
	for(i=0;i<new_mesh->n_elem;i++)
		ne_internal_list[i]=new_mesh->ne_internal_list[i];
	free(new_mesh->elem_type);
	free(new_mesh->index_elem);
	free(new_mesh->ptr_elem);
	free(new_mesh->elem_id);
	free(new_mesh->ne_internal_list);
	return;
}

void hecmw_mesh_index_c2f90_(int *neighbor_pe, int *import_index, int *export_index)
{
	int i,j;

	for(i=0;i<new_mesh->n_neighbor_pe;i++)
		neighbor_pe[i]=new_mesh->neighbor_pe[i];
	for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
		import_index[i]=new_mesh->import_index[i];
	for(i=0;i<new_mesh->n_neighbor_pe+1;i++)
		export_index[i]=new_mesh->export_index[i];
	return;
}

void hecmw_mesh_pe_c2f90_(int *import_node, int *export_node)
{
	int i,j;

	for(i=0;i<new_mesh->import_index[new_mesh->n_neighbor_pe];i++)
		import_node[i]=new_mesh->import_node[i];
	for(i=0;i<new_mesh->export_index[new_mesh->n_neighbor_pe];i++)
		export_node[i]=new_mesh->export_node[i];
	free(new_mesh->neighbor_pe);
	free(new_mesh->import_index);
	free(new_mesh->export_index);
	free(new_mesh->import_node);
	free(new_mesh->export_node);

	return;
}


void hecmw_mesh_adapt_c2f90_(int *WhenIwasRefined_node, int *WhenIwasRefined_elem,int *adaptation_parent_type, 
        int *adaptation_type, int *adaptation_level, int *adaptation_parent, int *adaptation_children, int *index_children) 
{
	int i,j;

	for(i=0;i<new_mesh->n_node;i++)
		WhenIwasRefined_node[i]=new_mesh->wheniwasrefined_node[i];
	for(i=0;i<new_mesh->n_elem;i++) {
        WhenIwasRefined_elem[i]=new_mesh->wheniwasrefined_elem[i];
        adaptation_parent_type[i]=new_mesh->adaptation_parent_type[i];
        adaptation_type[i]=new_mesh->adaptation_type[i];
        adaptation_level[i]=new_mesh->adaptation_level[i];
	}
	for(i=0;i<new_mesh->n_elem*2;i++)
		adaptation_parent[i]=new_mesh->adaptation_parent[i];
	for(i=0;i<new_mesh->n_elem+1;i++)
		index_children[i]=new_mesh->index_children[i];
	for(i=0;i<new_mesh->index_children[new_mesh->n_elem]*2;i++)
		adaptation_children[i]=new_mesh->adaptation_children[i];
	free(new_mesh->wheniwasrefined_elem);
	free(new_mesh->wheniwasrefined_node);
	free(new_mesh->adaptation_parent_type);
	free(new_mesh->adaptation_type);
	free(new_mesh->adaptation_level);
	free(new_mesh->adaptation_parent);
	free(new_mesh->index_children);
	free(new_mesh->adaptation_children);

	return;
}

void hecmw_data_c2f90_(double *U, double *P, double *VISCT, double *VISCL)
{
	int i;

	for(i=0;i<new_mesh->n_node*5;i++)
		U[i]=new_node->data[i];
	for(i=0;i<new_mesh->n_node;i++)
		P[i]=new_node->data[i+5*new_mesh->n_node];
	for(i=0;i<new_mesh->n_node;i++)
		VISCT[i]=new_node->data[i+6*new_mesh->n_node];
	for(i=0;i<new_mesh->n_node;i++)
		VISCL[i]=new_node->data[i+7*new_mesh->n_node];
	free(new_node->data);
	free(new_node->n_free);

	
	
	return;
}


void  hecmw_grp_index_node_c2f90_(int *enum_grp_index) 
{
	int i;
	for(i=0;i<new_grp->node_grp.n_enum_grp+1;i++)
		enum_grp_index[i]=new_grp->node_grp.enum_grp_index[i];
	return;
}

void hecmw_grp_node_node_c2f90_(int *enum_grp_node)
{
	int i;
	for(i=0;i<new_grp->node_grp.enum_grp_index[new_grp->node_grp.n_enum_grp];i++)
		enum_grp_node[i]=new_grp->node_grp.enum_grp_node[i];
	

	return;
}

*/

		
