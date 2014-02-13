/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : HEC-MW Utility                                    *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/

#ifndef INC_HECMW_PART_STRUCT
#define INC_HECMW_PART_STRUCT

#include "hecmw_part_define.h"


struct hecmw_part_edge_data {
    
    int n_edge;
    
    int *edge_node_item;
};


struct hecmw_part_node_data {
    
    int *node_elem_index;
    
    int *node_elem_item;
};


struct hecmw_part_cont_data {
    
    int n_domain;
    
    int depth;
    
    int type;
    
    int method;
    
    int n_rcb_div;
    
    int *rcb_axis;
    
    int is_print_ucd;
    
    char ucd_file_name[HECMW_FILENAME_LEN+1]; /* ucd file name */
    
    int n_my_domain;
    
    int *my_domain;
    
    int contact;
};

#endif  /* INC_HECMW_PART_STRUCT */
