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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include  <memory.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_io.h"
#include "parmetislib.h"
#include "hecmw_dlb_comm_util.h"
#include "mpi.h"
#define MASTER_PE 	0

#define PI  3.1415926
#define HEX_N_NODE	8
#define HEX_N_FACE	6
#define HEX_NODE_INDEX	255	/* 2^8 */
#define HEX_FACE_INDEX	63	/* 2^6 */

#define	MAX_LINE_LEN	256

#define NUM_CONTROL_PARAS 9

typedef struct _tmp_grp_inf {
	int   num_of_item;
	int   *item;
} Tmp_grp_inf;

typedef struct _control_para_struct {
	char    adaptive_repartition[4];
	int     num_criteria;        
	float   *balance_rate;
	int     num_repartition;  
    float   itr_rate; 
	int     wgtflag;
         /*  0----- no weights (vwgt and adjwgt are both NULL
		      1 ---- Weights on the edges only (vwgt is NULL)
			  2 ----- Weights on the vertices only (adjwgt is NULL)
			  3 ----- Weights on both the vertices and edges)
			  */
    char     vwgt_filename[128];
    char     adjwgt_filename[128];
    float    *machine_wgt;
    char     output_filename[128]; 
} Control_para;

typedef struct _result_partition_struct {
	int  edgecut;
	int  t_node;
	int  *part;
} Result_part;

struct _adj_find_struct {
	int vertex_num;
    struct _adj_find_struct *next_vertex;
};
typedef struct _adj_find_struct Adj_find;

struct _import_link_struct {
	int node_num;  /* local_id in current PE */
	int local_id;  /* local_id in import PE */
	struct _import_link_struct *next_node;
};
typedef struct _import_link_struct Import_link_struct;
int repart_comm;
struct hecmwST_local_mesh *mesh;
struct hecmwST_local_mesh *new_mesh;
void HECMW_dlb_print_exit(char *str_msg);
void HECMW_dlb_memory_exit(char *str_msg);
struct hecmwST_result_data *data;
struct hecmwST_result_data *new_data;
