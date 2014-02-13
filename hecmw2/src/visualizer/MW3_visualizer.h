/*=====================================================================*
!                                                                      !
! Software Name : HEC-MW Ver 4.3                                      !
!                                                                      !
!      Module Name : Visualizer Utility                                !
!                                                                      !
!            Written by Keiji Suemitsu (Advancesoft)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
*=====================================================================*/

#ifndef MW3_VIS_INCLUDED
#define MW3_VIS_INCLUDED

struct hecmwST_local_mesh* MW3_get_mesh( int level, int partID );
void MW3_free_mesh( struct hecmwST_local_mesh* mesh );
int MW3_visualize_init( void );
int MW3_visualize_init_ex( char* ctrlfile );
int MW3_visualize_finalize( void );
int MW3_result_init( struct hecmwST_local_mesh* mesh, int tstep, char* header );
int MW3_result_init_ex( char* ctrlfile, struct hecmwST_local_mesh* mesh, int tstep, char* header );
int MW3_result_finalize( void );
void mw3_visualize_init_( int* level, int* partID );
void mw3_visualize_init_ex_( char* ctrlfile, int* level, int* partID, int len );
void mw3_visualize_( int* timestep, int* max_timestep, int* is_force );
void mw3_visualize_finalize_( void );
void mw3_result_init_( int* level, int* partID, int* tstep, char* header, int len );
void mw3_result_init_ex_( char* ctrlfile, int* level, int* partID, int* tstep, char* header, int lenc, int lenh );
void mw3_result_add_( int* node_or_elem, int* n_dof, char* label, double* ptr, int len);
void mw3_result_write_by_name_( char* name_ID, int len );
void mw3_result_finalize_( void );

#endif
