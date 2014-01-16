/*=====================================================================!
!                                                                      !
!  Software Name : hpcmw-PC-cluster-ver.0.10                           !
!                                                                      !
!    main program:                                                     !
!    subroutine  :                                                     !
!    fuction     :                                                     !
!    module      :                                                     !
!    class       :                                                     !
!    library     : Dynamic Load Balancing                                       !
!                                                                      !
!    coded by Li Chen   (RIST)   2003/12/15                            !
!                                                                      !
!    Contact address :  The University of Tokyo, FSIS project          !
!                                                                      !
!    "High-Performance Computing Middleware (HPC-MW)" Group.           !
!                                                                      !
!=====================================================================*/


#define NUM_CONTROL_PARAS 9
/*
typedef struct _tmp_grp_inf {
	int   num_of_item;
	int   *item;
} Tmp_grp_inf;
*/

typedef struct _separation_result {
	int   num_of_lgraph;
	int   num_of_rgraph;
	int   num_of_separator;
	int   *lgraph;
	int   *rgraph;
	int   *mseparator;
} Separator_result;

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
void separator_print_exit(char *str_msg);
void separator_memory_exit(char *str_msg);
