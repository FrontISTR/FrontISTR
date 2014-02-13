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


#ifndef INC_HECMW_PART_LOG
#define INC_HECMW_PART_LOG

extern int
HECMW_part_init_log( int _n_domain );

extern int
HECMW_part_set_log_part_type( int _part_type );
extern int
HECMW_part_set_log_part_method( int _part_method );
extern int
HECMW_part_set_log_part_depth( int _depth );
extern int
HECMW_part_set_log_n_edgecut( int _n_edge, int _n_edgecut );
extern int
HECMW_part_set_log_n_node_g( int _n_node_g );
extern int
HECMW_part_set_log_n_elem_g( int _n_elem_g );
extern int
HECMW_part_set_log_n_node( int domain, int _n_node );
extern int
HECMW_part_set_log_n_elem( int domain, int _n_elem );
extern int
HECMW_part_set_log_nn_internal( int domain, int _nn_internal );
extern int
HECMW_part_set_log_ne_internal( int domain, int _ne_internal );

extern int
HECMW_part_print_log( void );

extern void
HECMW_part_finalize_log( void );

#endif  /* INC_HECMW_PART_LOG */
