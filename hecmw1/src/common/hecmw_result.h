/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#ifndef HECMW_RESULT_INCLUDED
#define HECMW_RESULT_INCLUDED

#include "hecmw_struct.h"

struct hecmwST_result_data {
	int nn_component;
	int ne_component;
	int *nn_dof;
	int *ne_dof;
	char **node_label;
	char **elem_label;
	double *node_val_item;
	double *elem_val_item;
};


extern int HECMW_result_init(struct hecmwST_local_mesh *hecMESH, int step, char *header);
extern int HECMW_result_init_body(int n_node, int n_elem, int *nodeID, int *elemID, int step, char *header);
extern int HECMW_result_finalize(void);
extern int HECMW_result_add(int node_or_elem, int n_dof, char *label, double *ptr);
extern int HECMW_result_write(void);
extern int HECMW_result_write_by_name(char *name_ID);

extern int HECMW_result_write_by_fname(char *filename);  /* == HECMW_result_write_txt_by_fname */
extern int HECMW_result_write_bin_by_fname(char *filename);
extern int HECMW_result_write_txt_by_fname(char *filename);

extern int HECMW_result_write_ST_by_fname(char *filename, struct hecmwST_result_data *result, int n_node, int n_elem, char *header);
extern int HECMW_result_write_bin_ST_by_fname(char *filename, struct hecmwST_result_data *result, int n_node, int n_elem, char *header);
extern int HECMW_result_write_txt_ST_by_fname(char *filename, struct hecmwST_result_data *result, int n_node, int n_elem, char *header);

extern int HECMW_result_write_ST_by_name(char *name_ID, struct hecmwST_result_data *result, int n_node, int n_elem, int tstep, char *header);
extern int HECMW_result_write_ST(struct hecmwST_result_data *result, int n_node, int n_elem, int tstep, char *header);

extern struct hecmwST_result_data *HECMW_result_read_by_fname(char *filename);
extern struct hecmwST_result_data *HECMW_result_read_by_name(char *name_ID, int tstep);
extern struct hecmwST_result_data *HECMW_result_read(int tstep);

extern void HECMW_result_free(struct hecmwST_result_data *result);
extern int HECMW_result_get_nnode(void);
extern int HECMW_result_get_nelem(void);
extern char* HECMW_result_get_header(char* buff);
extern int* HECMW_result_get_nodeID(int* buff);
extern int* HECMW_result_get_elemID(int* buff);

#endif
