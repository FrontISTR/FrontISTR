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

void HECMW_put_result_from_structure(struct hecmwST_local_mesh *mesh, struct hecmwST_result_data *data, char *resultfile_dist)
{ 
	FILE *fp;
	int  i,j,k;
	int  tn_component, tmp_count;
	

	if((fp=fopen(resultfile_dist,"w"))== NULL) 
	  HECMW_dlb_print_exit("ERROR: HEC-MW-VIS-E0011: Cannot open output file");
    fprintf(fp, "Dynamic_load_balancing_result\n");
	fprintf(fp, "%d %d\n", mesh->n_node, mesh->n_elem);
	fprintf(fp, "%d %d\n", data->nn_component, data->ne_component);
	if(data->nn_component>0) {
    	for(i=0;i<data->nn_component;i++)
	    	fprintf(fp, "%d ", data->nn_dof[i]);
	    fprintf(fp, "\n");
    	for(i=0;i<data->nn_component;i++)
			fprintf(fp, "%s\n", data->node_label[i]);
		tn_component=0;
		for(i=0;i<data->nn_component;i++)
			tn_component+=data->nn_dof[i];
		if(tn_component<=5) {
			for(i=0;i<mesh->n_node;i++) {
				for(j=0;j<tn_component;j++)
					fprintf(fp, "%e ", data->node_val_item[i*tn_component+j]);
				fprintf(fp, "\n");
			}
		}
		else {
			tmp_count=0;
			for(i=0;i<mesh->n_node;i++) {
				for(j=0;j<tn_component;j++) {
					fprintf(fp, "%e ", data->node_val_item[i*tn_component+j]);
					tmp_count++;
					if(tmp_count==5) {
						tmp_count=0;
			        	fprintf(fp, "\n");
					}
				}
			}
		}
	}
	return;
}







