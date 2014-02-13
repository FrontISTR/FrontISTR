/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include "hecmw_struct.h"
#include "hecmw_util.h"
#include "hecmw_dist_alloc.h"
#include "hecmw_io_dist.h"

#define FSCANF_RTC_TOO_SHORT (-10001)
#define FSCNAF_RTC_EOF       (-10002)

#define COLS_INT_DEF    10
#define COLS_DOUBLE_DEF  5
#define COLS_TWO         2
#define COLS_THREE       3

#define FILE_MAGIC "!HECMW-DMD-ASCII"
#define FILE_MAGIC_LEN (sizeof(FILE_MAGIC)-1)
#define HEADER_STRING "!HECMW-DMD-ASCII version="

/*                                                                            */
/*  get data                                                                  */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  integer                                                                   */
/*----------------------------------------------------------------------------*/
static int
get_int( int *i, FILE *fp )
{
  int rtc;

  rtc = fscanf( fp, "%d", i );
  if( rtc < 1 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  } else if( rtc == EOF ) {
    HECMW_set_error(HECMW_IO_E5003, "");
    return -1;
  } else {
    return 0;
  }
}

/*----------------------------------------------------------------------------*/
/*  double                                                                    */
/*----------------------------------------------------------------------------*/
static int
get_double( double *d, FILE *fp )
{
  int rtc;

  rtc = fscanf( fp, "%lf", d );
  if( rtc < 1 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  } else if( rtc == EOF ) {
    HECMW_set_error(HECMW_IO_E5003, "");
    return -1;
  } else {
    return 0;
  }
}

/*----------------------------------------------------------------------------*/
/*  string                                                                    */
/*----------------------------------------------------------------------------*/
static int
get_string( char *s, int max, FILE *fp )
{
    int c,len;

    while((c = fgetc(fp)) != EOF && isspace(c));  /* skip */
    if(c == EOF) {
    	HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
    if(ungetc(c, fp) == EOF) {
    	HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
    if(fgets(s, max, fp) == NULL) {
    	HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
    len = strlen(s);
    if(len == max-1 && s[max-2] != '\n') {
    	HECMW_set_error(HECMW_IO_E5004, "line too long");
      return -1;
    }
    while(len > 0 && isspace(s[len-1])) {
        len--;
    }
    s[len] = '\0';

    return strlen(s);
}

/*----------------------------------------------------------------------------*/
/*  HECMW_Comm                                                                */
/*----------------------------------------------------------------------------*/
static int
get_comm( HECMW_Comm *i, FILE *fp )
{
  int rtc;

  rtc = fscanf( fp, "%d", (int *)i );
  if( rtc < 1 ) {
   	HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  } else if( rtc == EOF ) {
   	HECMW_set_error(HECMW_IO_E5003, "");
    return -1;
  } else {
    return 0;
  }
}

/*----------------------------------------------------------------------------*/
/*  array ( int )                                                             */
/*----------------------------------------------------------------------------*/
static int
get_int_ary( int *ary, int n, FILE *fp )
{
  int rtc, i;

  for( i=0; i<n; i++ ) {
    rtc = fscanf( fp, "%d", &ary[i] );
    if( rtc < 1 ) {
   	  HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    } else if( rtc == EOF ) {
   	  HECMW_set_error(HECMW_IO_E5003, "");
      return -1;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  array ( double )                                                          */
/*----------------------------------------------------------------------------*/
static int
get_double_ary( double *ary, int n, FILE *fp )
{
  int rtc, i;

  for( i=0; i<n; i++ ) {
    rtc = fscanf( fp, "%lf", &ary[i] );
    if( rtc < 1 ) {
   	  HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    } else if( rtc == EOF ) {
   	  HECMW_set_error(HECMW_IO_E5003, "");
      return -1;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  array ( string )                                                          */
/*----------------------------------------------------------------------------*/
static int
get_string_ary( char **ary, int n, FILE *fp )
{
#define LINEBUF_SIZE 8096
  int i;
  char linebuf[LINEBUF_SIZE];

  for( i=0; i<n; i++ ) {
    if(get_string(linebuf, sizeof(linebuf), fp) < 0) {
        return -1;
    }
    if( ( ary[i] = (char *)HECMW_strdup( linebuf )) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
  }

  return 0;
}


/*============================================================================*/
/*                                                                            */
/*  check file type                                                           */
/*                                                                            */
/*============================================================================*/
static int
is_hecmw_dist_file( FILE *fp)
{
  char filetype[FILE_MAGIC_LEN];

  if( fread(filetype, FILE_MAGIC_LEN, 1, fp) != 1) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return 0;
  }
  if(memcmp(filetype, FILE_MAGIC, FILE_MAGIC_LEN)) {
    HECMW_set_error(HECMW_IO_E5005, "Not a HECMW-DIST ASCII file");
    return 0;
  }
  return 1;
}


/*============================================================================*/
/*                                                                            */
/*  rewind FILE pointer                                                       */
/*                                                                            */
/*============================================================================*/
static int
rewind_fp( FILE *fp )
{
	if(fseek(fp, 0L, SEEK_SET)) {
		HECMW_set_error(HECMW_IO_E5004, "");
		return -1;
	}
	return 0;
}


/*============================================================================*/
/*                                                                            */
/*  get header                                                                */
/*                                                                            */
/*============================================================================*/
static int
get_header( FILE *fp )
{
  char header[HECMW_HEADER_LEN+1];
  int ver;

  if( fgets( header, sizeof(header), fp ) == NULL ) {
   	HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }
  if( strlen(header) == sizeof(header)-1 && header[strlen(header)-2] != '\n') {
   	HECMW_set_error(HECMW_IO_E5004 , "line too long");
    return -1;
  }
  if( strncmp(header, HEADER_STRING, strlen(HEADER_STRING)) != 0) {
    HECMW_set_error(HECMW_IO_E5005, "Not a HECMW-DIST file");
    return -1;
  }
  if(sscanf(header+strlen(HEADER_STRING), "%d", &ver) != 1) {
    HECMW_set_error(HECMW_IO_E5006, "Invalid version");
    return -1;
  }
  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get global information                                                    */
/*                                                                            */
/*============================================================================*/
static int
get_global_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  int flag_header;

  /* hecmw_flag_adapt */
  if( get_int( &mesh->hecmw_flag_adapt, fp ) ) {
    return -1;
  }

  /* hecmw_flag_initcon */
  if( get_int( &mesh->hecmw_flag_initcon, fp ) ) {
    return -1;
  }

  /* hecmw_flag_parttype */
  if( get_int( &mesh->hecmw_flag_parttype, fp ) ) {
    return -1;
  }

  /* hecmw_flag_partdepth */
  if( get_int( &mesh->hecmw_flag_partdepth, fp ) ) {
    return -1;
  }

  /* hecmw_flag_version */
  if( get_int( &mesh->hecmw_flag_version, fp ) ) {
    return -1;
  }

  /* gridfile */
  if( get_string( mesh->gridfile, sizeof(mesh->gridfile), fp ) < 0 ) {
    return -1;
  }

  /* hecmw_n_file */
  if( get_int( &mesh->hecmw_n_file, fp ) ) {
    return -1;
  }

  /* files */
  if( mesh->hecmw_n_file > 0 ) {
    if( ( mesh->files = (char **)HECMW_calloc( mesh->hecmw_n_file, sizeof(char *) ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_string_ary( mesh->files, mesh->hecmw_n_file, fp ) ) {
      return -1;
    }
  } else {
    mesh->files = NULL;
  }

  /* flag for header */
  if( get_int( &flag_header, fp ) ) {
    return -1;
  }

  if( flag_header == 1 ) {
    /* header */
    if( get_string( mesh->header, sizeof(mesh->header), fp ) < 0 ) {
      return -1;
    }
  }

  /* zero_temp */
  if( get_double( &mesh->zero_temp, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get node information                                                      */
/*                                                                            */
/*============================================================================*/
static int
get_node_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* n_node*/
  if( get_int( &mesh->n_node, fp ) ) {
    return -1;
  }

  if( mesh->hecmw_flag_version >= 2 ) {
    /* n_node_gross*/
    if( get_int( &mesh->n_node_gross, fp ) ) {
      return -1;
    }
  } else {
    mesh->n_node_gross = mesh->n_node;
  }

  /* nn_internal */
  if( get_int( &mesh->nn_internal, fp ) ) {
    return -1;
  }

  /* node_internal_list */
  if((mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_ELEMBASED
   || mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_UNKNOWN )
   && mesh->nn_internal > 0) {
    mesh->node_internal_list = (int *)HECMW_malloc( sizeof(double)*mesh->nn_internal );
    if( mesh->node_internal_list == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_internal_list, mesh->nn_internal, fp ) ) {
      return -1;
    }
  }

  if(mesh->n_node_gross > 0) {
    /* node_ID */
    mesh->node_ID = (int *)HECMW_malloc( sizeof(int)*mesh->n_node_gross*2 );
    if( mesh->node_ID == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_ID, mesh->n_node_gross*2, fp ) ) {
      return -1;
    }

    /* global_node_ID */
    mesh->global_node_ID = (int *)HECMW_malloc( sizeof(int)*mesh->n_node_gross );
    if( mesh->global_node_ID == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->global_node_ID, mesh->n_node_gross, fp ) ) {
      return -1;
    }

    /* node */
    mesh->node = (double *)HECMW_malloc( sizeof(double)*mesh->n_node_gross*3 );
    if( mesh->node == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_double_ary( mesh->node, mesh->n_node_gross*3, fp ) ) {
      return -1;
    }
  }

  /* n_dof */
  if( get_int( &mesh->n_dof, fp ) ) {
    return -1;
  }

  /* n_dof_grp */
  if( get_int( &mesh->n_dof_grp, fp ) ) {
    return -1;
  }

  if(mesh->n_dof_grp > 0) {
    /* node_dof_index */
    mesh->node_dof_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_dof_grp+1) );
    if( mesh->node_dof_index == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_dof_index, mesh->n_dof_grp+1, fp ) ) {
      return -1;
    }

    /* node_dof_item */
    mesh->node_dof_item = (int *)HECMW_malloc( sizeof(int)*mesh->n_dof_grp );
    if( mesh->node_dof_item == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_dof_item, mesh->n_dof_grp, fp ) ) {
      return -1;
    }
  }

  if( mesh->hecmw_flag_initcon && mesh->n_node_gross > 0) {

    /* node_init_val_index */
    mesh->node_init_val_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_node_gross+1) );
    if( mesh->node_init_val_index == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_init_val_index, mesh->n_node_gross+1, fp ) ) {
      return -1;
    }

    /* node_init_val_item */
    if( mesh->node_init_val_index[mesh->n_node_gross] == 0 ) {
      mesh->node_init_val_item = NULL;
    } else {
      mesh->node_init_val_item = (double *)HECMW_malloc( sizeof(double)*mesh->node_init_val_index[mesh->n_node_gross] );
      if( mesh->node_init_val_item == NULL ) {
        HECMW_set_error(errno, "");
        return -1;
      }
      if( get_double_ary( mesh->node_init_val_item, mesh->node_init_val_index[mesh->n_node_gross], fp ) ) {
        return -1;
      }
    }
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get element information                                                   */
/*                                                                            */
/*============================================================================*/
static int
get_elem_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* n_elem */
  if( get_int( &mesh->n_elem, fp ) ) {
    return -1;
  }

  if( mesh->hecmw_flag_version >= 2 ) {
    /* n_elem_gross */
    if( get_int( &mesh->n_elem_gross, fp ) ) {
      return -1;
    }
  } else {
    mesh->n_elem_gross = mesh->n_elem;
  }

  /* ne_internal */
  if( get_int( &mesh->ne_internal, fp ) ) {
    return -1;
  }

  /* elem_internal_list */
  if((mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_NODEBASED
   || mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_UNKNOWN)
   && mesh->ne_internal > 0 ) {
    mesh->elem_internal_list = (int *)HECMW_malloc( sizeof(int)*mesh->ne_internal );
    if( mesh->elem_internal_list == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_internal_list, mesh->ne_internal, fp ) ) {
      return -1;
    }
  }

  if( mesh->n_elem_gross > 0) {
    /* elem_ID */
    if( ( mesh->elem_ID = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross*2 ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_ID, mesh->n_elem_gross*2, fp ) ) {
      return -1;
    }

    /* global_elem_ID */
    if( ( mesh->global_elem_ID = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->global_elem_ID, mesh->n_elem_gross, fp ) ) {
      return -1;
    }

    /* elem_type */
    if( ( mesh->elem_type = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_type, mesh->n_elem_gross, fp ) ) {
      return -1;
    }
  }

  /* n_elem_type */
  if( get_int( &mesh->n_elem_type, fp ) ) {
    return -1;
  }

  if( mesh->n_elem_type > 0 ) {
    /* elem_type_index */
    if( ( mesh->elem_type_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_elem_type+1) ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_type_index, mesh->n_elem_type+1, fp ) ) {
      return -1;
    }

    /* elem_type_item */
    if( ( mesh->elem_type_item = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_type ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_type_item, mesh->n_elem_type, fp ) ) {
      return -1;
    }
  }

  if( mesh->n_elem_gross > 0 ) {
    /* elem_node_index */
    if( ( mesh->elem_node_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_elem_gross+1) ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_node_index, mesh->n_elem_gross+1, fp ) ) {
      return -1;
    }

    /* elem_node_item */
    if( ( mesh->elem_node_item = (int *)HECMW_malloc( sizeof(int)*mesh->elem_node_index[mesh->n_elem_gross] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_node_item, mesh->elem_node_index[mesh->n_elem_gross], fp ) ) {
      return -1;
    }

    /* section_ID */
    if( ( mesh->section_ID = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->section_ID, mesh->n_elem_gross, fp ) ) {
      return -1;
    }

    /* elem_mat_ID_index */
    if( ( mesh->elem_mat_ID_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_elem_gross+1) ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_mat_ID_index, mesh->n_elem_gross+1, fp ) ) {
      return -1;
    }

    /* elem_mat_ID_item */
    if( ( mesh->elem_mat_ID_item = (int *)HECMW_malloc( sizeof(int)*mesh->elem_mat_ID_index[mesh->n_elem_gross] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_mat_ID_item, mesh->elem_mat_ID_index[mesh->n_elem_gross], fp ) ) {
      return -1;
    }
  }

  /* n_elem_mat_ID */
  if( get_int( &mesh->n_elem_mat_ID, fp ) ) {
    return -1;
  }

  /* elem_mat_int_index */
  /* elem_mat_int_item */

  /* elem_val_index */
  /* elem_val_item */

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get domain information & communication table                              */
/*                                                                            */
/*============================================================================*/
static int
get_comm_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* zero */
  if( get_int( &mesh->zero, fp ) ) {
    return -1;
  }

  /* HECMW_COMM */
  if( get_comm( &mesh->HECMW_COMM, fp ) ) {
    return -1;
  }

  /* PETOT */
  if( get_int( &mesh->PETOT, fp ) ) {
    return -1;
  }

  /* PEsmpTOT */
  if( get_int( &mesh->PEsmpTOT, fp ) ) {
    return -1;
  }

  /* my_rank */
  if( get_int( &mesh->my_rank, fp ) ) {
    return -1;
  }

  /* errnof */
  if( get_int( &mesh->errnof, fp ) ) {
    return -1;
  }

  /* n_subdomain */
  if( get_int( &mesh->n_subdomain, fp ) ) {
    return -1;
  }

  /* n_neighbor_pe */
  if( get_int( &mesh->n_neighbor_pe, fp ) ) {
    return -1;
  }



  if( mesh->n_neighbor_pe == 0 ) {
    mesh->neighbor_pe = NULL;
    mesh->import_item  = NULL;
    mesh->export_item  = NULL;
    mesh->shared_item  = NULL;

    mesh->import_index = HECMW_malloc(sizeof(int));
	if(mesh->import_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	mesh->import_index[0] = 0;

    mesh->export_index = HECMW_malloc(sizeof(int));
	if(mesh->export_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	mesh->export_index[0] = 0;

    mesh->shared_index = HECMW_malloc(sizeof(int));
	if(mesh->shared_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	mesh->shared_index[0] = 0;

    return 0;
  }

  /* neighbor_pe */
  if( ( mesh->neighbor_pe = (int *)HECMW_malloc( sizeof(int)*mesh->n_neighbor_pe ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->neighbor_pe, mesh->n_neighbor_pe, fp ) ) {
    return -1;
  }

  /* import_index */
  if( ( mesh->import_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_neighbor_pe+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->import_index, mesh->n_neighbor_pe+1, fp ) ) {
    return -1;
  }

  /* import_item */
  if( ( mesh->import_item = (int *)HECMW_malloc( sizeof(int)*mesh->import_index[mesh->n_neighbor_pe] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->import_item, mesh->import_index[mesh->n_neighbor_pe], fp ) ) {
    return -1;
  }

  /* export_index */
  if( ( mesh->export_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_neighbor_pe+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->export_index, mesh->n_neighbor_pe+1, fp ) ) {
    return -1;
  }

  /* export_item */
  if( ( mesh->export_item = (int *)HECMW_malloc( sizeof(int)*mesh->export_index[mesh->n_neighbor_pe] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->export_item, mesh->export_index[mesh->n_neighbor_pe], fp ) ) {
    return -1;
  }

  /* shared_index */
  if( ( mesh->shared_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_neighbor_pe+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->shared_index, mesh->n_neighbor_pe+1, fp ) ) {
    return -1;
  }

  /* shared_item */
  if( ( mesh->shared_item = (int *)HECMW_malloc( sizeof(int)*mesh->shared_index[mesh->n_neighbor_pe] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->shared_item, mesh->shared_index[mesh->n_neighbor_pe], fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get adaptation information                                                */
/*                                                                            */
/*============================================================================*/
static int
get_adapt_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  if( mesh->hecmw_flag_adapt == 0 ) {
    mesh->coarse_grid_level       = 0;
    mesh->n_adapt                 = 0;
    mesh->when_i_was_refined_node = NULL;
    mesh->when_i_was_refined_elem = NULL;
    mesh->adapt_parent_type       = NULL;
    mesh->adapt_type              = NULL;
    mesh->adapt_level             = NULL;
    mesh->adapt_parent            = NULL;
    mesh->adapt_children_index    = NULL;
    mesh->adapt_children_item     = NULL;
    return 0;
  }

  /* coarse_grid_level */
  if( get_int( &mesh->coarse_grid_level, fp ) ) {
    return -1;
  }

  /* n_adapt */
  if( get_int( &mesh->n_adapt, fp ) ) {
    return -1;
  }

  if( mesh->n_node_gross > 0) {
    /* when_i_was_refined_node */
    if( ( mesh->when_i_was_refined_node = (int *)HECMW_malloc( sizeof(int)*mesh->n_node_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->when_i_was_refined_node, mesh->n_node_gross, fp ) ) {
      return -1;
    }
  }

  if( mesh->n_elem_gross <= 0 ) return 0;

  /* when_i_was_refined_elem */
  if( ( mesh->when_i_was_refined_elem = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->when_i_was_refined_elem, mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* adapt_parent_type */
  if( ( mesh->adapt_parent_type = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_parent_type, mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* adapt_type */
  if( ( mesh->adapt_type = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_type, mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* adapt_level */
  if( ( mesh->adapt_level = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_level, mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* adapt_parent */
  if( ( mesh->adapt_parent = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross*2 ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_parent, 2*mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* adapt_children_index */
  if( ( mesh->adapt_children_index = (int *)HECMW_malloc( sizeof(int)*(mesh->n_elem_gross+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_children_index, mesh->n_elem_gross+1, fp ) ) {
    return -1;
  }

  /* adapt_children_item */
  if( ( mesh->adapt_children_item = (int *)HECMW_malloc( sizeof(int)*mesh->adapt_children_index[mesh->n_elem_gross]*2 ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mesh->adapt_children_item, mesh->adapt_children_index[mesh->n_elem_gross]*2, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get section information                                                   */
/*                                                                            */
/*============================================================================*/
static int
get_section_info( struct hecmwST_section *sect, FILE *fp )
{
  /* n_sect */
  if( get_int( &sect->n_sect, fp ) ) {
    return -1;
  }

  if( sect->n_sect == 0 ) {
    sect->sect_type         = NULL;
    sect->sect_opt          = NULL;
    sect->sect_mat_ID_index = NULL;
    sect->sect_mat_ID_item  = NULL;
    sect->sect_I_index      = NULL;
    sect->sect_I_item       = NULL;
    sect->sect_R_index      = NULL;
    sect->sect_R_item       = NULL;
    return 0;
  }

  /* sect_type */
  if( ( sect->sect_type = (int *)HECMW_malloc( sizeof(int)*sect->n_sect ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( sect->sect_type, sect->n_sect, fp ) ) {
    return -1;
  }

  /* sect_opt */
  if( ( sect->sect_opt = (int *)HECMW_malloc( sizeof(int)*sect->n_sect ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( sect->sect_opt, sect->n_sect, fp ) ) {
    return -1;
  }

  /* sect_mat_ID_index */
  if( ( sect->sect_mat_ID_index = (int *)HECMW_malloc( sizeof(int)*(sect->n_sect+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( sect->sect_mat_ID_index, sect->n_sect+1, fp ) ) {
    return -1;
  }

  /* sect_mat_ID_item */
  if( sect->sect_mat_ID_index[sect->n_sect] > 0 ) {
    if( ( sect->sect_mat_ID_item = (int *)HECMW_malloc( sizeof(int)*sect->sect_mat_ID_index[sect->n_sect] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( sect->sect_mat_ID_item, sect->sect_mat_ID_index[sect->n_sect], fp ) ) {
      return -1;
    }
  }

  /* sect_I_index */
  if( ( sect->sect_I_index = (int *)HECMW_malloc( sizeof(int)*(sect->n_sect+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( sect->sect_I_index, sect->n_sect+1, fp ) ) {
    return -1;
  }

  /* sect_I_item */
  if( sect->sect_I_index[sect->n_sect] > 0 ) {
    if( ( sect->sect_I_item = (int *)HECMW_malloc( sizeof(int)*sect->sect_I_index[sect->n_sect] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( sect->sect_I_item, sect->sect_I_index[sect->n_sect], fp ) ) {
      return -1;
    }
  }

  /* sect_R_index */
  if( ( sect->sect_R_index = (int *)HECMW_malloc( sizeof(int)*(sect->n_sect+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( sect->sect_R_index, sect->n_sect+1, fp ) ) {
    return -1;
  }

  /* sect_R_item */
  if( sect->sect_R_index[sect->n_sect] > 0 ) {
    if( ( sect->sect_R_item = (double *)HECMW_malloc( sizeof(double)*sect->sect_R_index[sect->n_sect] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_double_ary( sect->sect_R_item, sect->sect_R_index[sect->n_sect], fp ) ) {
      return -1;
    }
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get material information                                                  */
/*                                                                            */
/*============================================================================*/
static int
get_material_info( struct hecmwST_material *mat, FILE *fp )
{
  /* n_mat */
  if( get_int( &mat->n_mat, fp ) ) {
    return -1;
  }

  if( mat->n_mat == 0 ) {
    mat->n_mat_item        = 0;
    mat->n_mat_subitem     = 0;
    mat->n_mat_table       = 0;
    mat->mat_name          = NULL;
    mat->mat_item_index    = NULL;
    mat->mat_subitem_index = NULL;
    mat->mat_table_index   = NULL;
    mat->mat_val           = NULL;
    mat->mat_temp          = NULL;
    return 0;
  }

  /* n_mat_item */
  if( get_int( &mat->n_mat_item, fp ) ) {
    return -1;
  }

  /* n_mat_subitem */
  if( get_int( &mat->n_mat_subitem, fp ) ) {
    return -1;
  }

  /* n_mat_table */
  if( get_int( &mat->n_mat_table, fp ) ) {
    return -1;
  }

  /* mat_name */
  if( ( mat->mat_name = (char **)HECMW_malloc( sizeof(char *)*mat->n_mat ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( mat->mat_name, mat->n_mat, fp ) ) {
    return -1;
  }

  /* mat_item_index */
  if( ( mat->mat_item_index = (int *)HECMW_malloc( sizeof(int)*(mat->n_mat+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mat->mat_item_index, mat->n_mat+1, fp ) ) {
    return -1;
  }

  /* mat_subitem_index */
  if( ( mat->mat_subitem_index = (int *)HECMW_malloc( sizeof(int)*(mat->n_mat_item+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mat->mat_subitem_index, mat->n_mat_item+1, fp ) ) {
    return -1;
  }

  /* mat_table_index */
  if( ( mat->mat_table_index = (int *)HECMW_malloc( sizeof(int)*(mat->n_mat_subitem +1)) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mat->mat_table_index, mat->n_mat_subitem+1, fp ) ) {
    return -1;
  }

  /* mat_val */
  if( ( mat->mat_val = (double *)HECMW_malloc( sizeof(double)*mat->n_mat_table ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_double_ary( mat->mat_val, mat->n_mat_table, fp ) ) {
    return -1;
  }

  /* mat_temp */
  if( ( mat->mat_temp = (double *)HECMW_malloc( sizeof(double)*mat->n_mat_table ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_double_ary( mat->mat_temp, mat->n_mat_table, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get MPC group information                                                 */
/*                                                                            */
/*============================================================================*/
static int
get_mpc_info( struct hecmwST_mpc *mpc, FILE *fp, int hecmw_flag_version )
{
  /* n_mpc */
  if( get_int( &mpc->n_mpc, fp ) ) {
    return -1;
  }

  if( mpc->n_mpc == 0 ) {
    mpc->mpc_item  = NULL;
    mpc->mpc_dof   = NULL;
    mpc->mpc_val   = NULL;
    mpc->mpc_index = HECMW_malloc(sizeof(int));
	if(mpc->mpc_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	mpc->mpc_index[0] = 0;
    return 0;
  }

  /* mpc_index */
  if( ( mpc->mpc_index = (int *)HECMW_malloc( sizeof(int)*(mpc->n_mpc+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mpc->mpc_index, mpc->n_mpc+1, fp ) ) {
    return -1;
  }

  /* mpc_item */
  if( ( mpc->mpc_item = (int *)HECMW_malloc( sizeof(int)*mpc->mpc_index[mpc->n_mpc] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mpc->mpc_item, mpc->mpc_index[mpc->n_mpc], fp ) ) {
    return -1;
  }

  /* mpc_dof */
  if( ( mpc->mpc_dof = (int *)HECMW_malloc( sizeof(int)*mpc->mpc_index[mpc->n_mpc] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( mpc->mpc_dof, mpc->mpc_index[mpc->n_mpc], fp ) ) {
    return -1;
  }

  /* mpc_val */
  if( ( mpc->mpc_val = (double *)HECMW_malloc( sizeof(double)*mpc->mpc_index[mpc->n_mpc] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_double_ary( mpc->mpc_val, mpc->mpc_index[mpc->n_mpc], fp ) ) {
    return -1;
  }

  /* mpc_const */
  if( ( mpc->mpc_const = (double *)HECMW_calloc( mpc->n_mpc, sizeof(double) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( hecmw_flag_version >= 3 ) {
    if( get_double_ary( mpc->mpc_const, mpc->n_mpc, fp ) ) {
      return -1;
    }
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get amplitude information                                                 */
/*                                                                            */
/*============================================================================*/
static int
get_amp_info( struct hecmwST_amplitude *amp, FILE *fp )
{
  /* n_amp */
  if( get_int( &amp->n_amp, fp ) ) {
    return -1;
  }

  if( amp->n_amp == 0 ) {
    amp->amp_name            = NULL;
    amp->amp_type_definition = NULL;
    amp->amp_type_time       = NULL;
    amp->amp_type_value      = NULL;
    amp->amp_val             = NULL;
    amp->amp_table           = NULL;
    amp->amp_index = HECMW_malloc(sizeof(int));
	if(amp->amp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	amp->amp_index[0] = 0;
    return 0;
  }

  /* amp_name */
  if( ( amp->amp_name = (char **)HECMW_malloc( sizeof(char *)*amp->n_amp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( amp->amp_name, amp->n_amp, fp ) ) {
    return -1;
  }

  /* amp_type_definition */
  if( ( amp->amp_type_definition = (int *)HECMW_malloc( sizeof(int)*amp->n_amp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( amp->amp_type_definition, amp->n_amp, fp ) ) {
    return -1;
  }

  /* amp_type_time */
  if( ( amp->amp_type_time = (int *)HECMW_malloc( sizeof(int)*amp->n_amp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( amp->amp_type_time, amp->n_amp, fp ) ) {
    return -1;
  }

  /* amp_type_value */
  if( ( amp->amp_type_value = (int *)HECMW_malloc( sizeof(int)*amp->n_amp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( amp->amp_type_value, amp->n_amp, fp ) ) {
    return -1;
  }

  /* amp_index */
  if( ( amp->amp_index = (int *)HECMW_malloc( sizeof(int)*(amp->n_amp+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( amp->amp_index, amp->n_amp+1, fp ) ) {
    return -1;
  }

  /* amp_val */
  if( ( amp->amp_val = (double *)HECMW_malloc( sizeof(double)*amp->amp_index[amp->n_amp] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_double_ary( amp->amp_val, amp->amp_index[amp->n_amp], fp ) ) {
    return -1;
  }

  /* amp_table */
  if( ( amp->amp_table = (double *)HECMW_malloc( sizeof(double)*amp->amp_index[amp->n_amp] ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_double_ary( amp->amp_table, amp->amp_index[amp->n_amp], fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get node group information                                                */
/*                                                                            */
/*============================================================================*/
static int
get_node_group_info( struct hecmwST_node_grp *grp, FILE *fp )
{
  /* n_grp */
  if( get_int( &grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 ) {
    grp->grp_name  = NULL;
    grp->grp_item  = NULL;
    grp->grp_index = HECMW_malloc(sizeof(int));
	if(grp->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp->grp_index[0] = 0;
    return 0;
  }

  /* grp_name */
  if( ( grp->grp_name = (char **)HECMW_malloc( sizeof(char *)*grp->n_grp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( ( grp->grp_index = (int *)HECMW_malloc( sizeof(int)*(grp->n_grp+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( grp->grp_index, grp->n_grp+1, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( grp->grp_index[grp->n_grp] > 0 ) {
    if( ( grp->grp_item = (int *)HECMW_malloc( sizeof(int)*grp->grp_index[grp->n_grp] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( grp->grp_item, grp->grp_index[grp->n_grp], fp ) ) {
      return -1;
    }
  } else {
    grp->grp_item = NULL;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get element group information                                             */
/*                                                                            */
/*============================================================================*/
static int
get_elem_group_info( struct hecmwST_elem_grp *grp, FILE *fp )
{
  /* n_grp */
  if( get_int( &grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 ) {
    grp->grp_name  = NULL;
    grp->grp_item  = NULL;
    grp->grp_index = HECMW_malloc(sizeof(int));
	if(grp->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp->grp_index[0] = 0;
    return 0;
  }

  /* grp_name */
  if( ( grp->grp_name = (char **)HECMW_malloc( sizeof(char *)*grp->n_grp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( ( grp->grp_index = (int *)HECMW_malloc( sizeof(int)*(grp->n_grp+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( grp->grp_index, grp->n_grp+1, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( grp->grp_index[grp->n_grp] > 0 ) {
    if( ( grp->grp_item = (int *)HECMW_malloc( sizeof(int)*grp->grp_index[grp->n_grp] ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( grp->grp_item, grp->grp_index[grp->n_grp], fp ) ) {
      return -1;
    }
  } else {
    grp->grp_item = NULL;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  get surface group information                                             */
/*                                                                            */
/*============================================================================*/
static int
get_surf_group_info( struct hecmwST_surf_grp *grp, FILE *fp )
{
  /* n_grp */
  if( get_int( &grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 ) {
    grp->grp_name  = NULL;
    grp->grp_item  = NULL;
    grp->grp_index = HECMW_malloc(sizeof(int));
	if(grp->grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp->grp_index[0] = 0;
    return 0;
  }

  /* grp_name */
  if( ( grp->grp_name = (char **)HECMW_malloc( sizeof(char *)*grp->n_grp ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( ( grp->grp_index = (int *)HECMW_malloc( sizeof(int)*(grp->n_grp+1) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( grp->grp_index, grp->n_grp+1, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( grp->grp_index[grp->n_grp] > 0 ) {
    if( ( grp->grp_item = (int *)HECMW_malloc( sizeof(int)*grp->grp_index[grp->n_grp]*2 ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( grp->grp_item, grp->grp_index[grp->n_grp]*2, fp ) ) {
      return -1;
    }
  } else {
    grp->grp_item = NULL;
  }

  return 0;
}


/*============================================================================*/
/*                                                                            */
/*  get contact information                                                   */
/*                                                                            */
/*============================================================================*/
static int
get_contact_info( struct hecmwST_contact_pair *cpair, FILE *fp, int hecmw_flag_version )
{
  if( hecmw_flag_version < 3 ) {
    cpair->n_pair        = 0;
    cpair->type     = NULL;
    cpair->slave_grp_id  = NULL;
    cpair->master_grp_id = NULL;
    return 0;
  }

  /* n_pair */
  if( get_int( &cpair->n_pair, fp ) ) {
    return -1;
  }

  if( cpair->n_pair == 0 ) {
    cpair->name     = NULL;
    cpair->type     = NULL;
    cpair->slave_grp_id  = NULL;
    cpair->master_grp_id = NULL;
    return 0;
  }

  /* name */
  if( ( cpair->name = (char **)HECMW_malloc( sizeof(char *)*cpair->n_pair ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_string_ary( cpair->name, cpair->n_pair, fp ) ) {
    return -1;
  }

  /* type */
  if( ( cpair->type = (int *)HECMW_malloc( sizeof(int)*(cpair->n_pair) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( cpair->type, cpair->n_pair, fp ) ) {
    return -1;
  }

  /* slave_grp_id */
  if( ( cpair->slave_grp_id = (int *)HECMW_malloc( sizeof(int)*(cpair->n_pair) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( cpair->slave_grp_id, cpair->n_pair, fp ) ) {
    return -1;
  }

  /* master_grp_id */
  if( ( cpair->master_grp_id = (int *)HECMW_malloc( sizeof(int)*(cpair->n_pair) ) ) == NULL ) {
    HECMW_set_error(errno, "");
    return -1;
  }
  if( get_int_ary( cpair->master_grp_id, cpair->n_pair, fp ) ) {
    return -1;
  }

  return 0;
}


/*============================================================================*/
/*                                                                            */
/*  get refinement information                                                */
/*                                                                            */
/*============================================================================*/
static int
get_refine_info( struct hecmwST_local_mesh *mesh, FILE *fp )
{
  if( mesh->hecmw_flag_version < 2 ) {
    mesh->n_refine     = 0;
    mesh->node_old2new = NULL;
    mesh->node_new2old = NULL;
    return 0;
  }

  /* number of refinement performed */
  if( get_int( &mesh->n_refine, fp ) ) {
    return -1;
  }
  if( mesh->n_refine == 0 || mesh->n_subdomain == 1 ) {
    mesh->node_old2new = NULL;
    mesh->node_new2old = NULL;
    return 0;
  }

  if( mesh->n_node_gross > mesh->nn_internal ) {
    /* node_old2new */
    if( ( mesh->node_old2new = (int *)HECMW_malloc( sizeof(int)*mesh->n_node_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_old2new, mesh->n_node_gross, fp ) ) {
      return -1;
    }

    /* node_new2old */
    if( ( mesh->node_new2old = (int *)HECMW_malloc( sizeof(int)*mesh->n_node_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->node_new2old, mesh->n_node_gross, fp ) ) {
      return -1;
    }
  }

  if( mesh->n_elem_gross > mesh->n_elem ) {
    /* elem_old2new */
    if( ( mesh->elem_old2new = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_old2new, mesh->n_elem_gross, fp ) ) {
      return -1;
    }

    /* elem_new2old */
    if( ( mesh->elem_new2old = (int *)HECMW_malloc( sizeof(int)*mesh->n_elem_gross ) ) == NULL ) {
      HECMW_set_error(errno, "");
      return -1;
    }
    if( get_int_ary( mesh->elem_new2old, mesh->n_elem_gross, fp ) ) {
      return -1;
    }
  }

  return 0;
}


static void
set_comm(struct hecmwST_local_mesh *mesh)
{
	mesh->HECMW_COMM = HECMW_comm_get_comm();
	mesh->PETOT = HECMW_comm_get_size();
	mesh->PEsmpTOT = 1;
	mesh->my_rank = HECMW_comm_get_rank();
	if(mesh->my_rank == 0) {
		mesh->zero = 1;
	} else {
		mesh->zero = 0;
	}
}

/*============================================================================*/
/*                                                                            */
/*  get HEC-MW distributed mesh format data                                   */
/*                                                                            */
/*============================================================================*/
extern struct hecmwST_local_mesh *
HECMW_get_dist_mesh( char *fname )
{
  FILE *fp;
  struct hecmwST_local_mesh *mesh;

  HECMW_log(HECMW_LOG_DEBUG, "Start to read HECW-DIST file");

  /* allocation */
  if( ( mesh = HECMW_dist_alloc( ) ) == NULL ) {
    return NULL;
  }

  /* file open */
  if( ( fp = fopen( fname, "r" ) ) == NULL ) {
   	HECMW_set_error(HECMW_IO_E5001, "File: %s, %s", fname, HECMW_strmsg(errno) );
    return NULL;
  }

  /* check file type */
  if( !is_hecmw_dist_file(fp) ) {
    return NULL;
  }

  /* rewind */
  if( rewind_fp(fp)) {
    return NULL;
  }

  /* get data */
  if( get_header( fp ) ) {
    return NULL;
  }

  if( get_global_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_node_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_elem_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_comm_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_adapt_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_section_info( mesh->section, fp ) ) {
    return NULL;
  }

  if( get_material_info( mesh->material, fp ) ) {
    return NULL;
  }

  if( get_mpc_info( mesh->mpc, fp, mesh->hecmw_flag_version ) ) {
    return NULL;
  }

  if( get_amp_info( mesh->amp, fp ) ) {
    return NULL;
  }

  if( get_node_group_info( mesh->node_group, fp ) ) {
    return NULL;
  }

  if( get_elem_group_info( mesh->elem_group, fp ) ) {
    return NULL;
  }

  if( get_surf_group_info( mesh->surf_group, fp ) ) {
    return NULL;
  }

  if( get_refine_info( mesh, fp ) ) {
    return NULL;
  }

  if( get_contact_info( mesh->contact_pair, fp, mesh->hecmw_flag_version ) ) {
    return NULL;
  }

  /* close file */
  if( fclose( fp ) ) {
   	HECMW_set_error(HECMW_IO_E5002, HECMW_strmsg( errno ) );
    return NULL;
  }

  set_comm(mesh);

  if( mesh->hecmw_flag_version < 3 ) {
    mesh->hecmw_flag_version = 3;
  }

  return mesh;
}


/*                                                                            */
/*  print data                                                                */
/*                                                                            */

/*----------------------------------------------------------------------------*/
/*  print int                                                                 */
/*----------------------------------------------------------------------------*/
static int
print_int( int item, FILE *fp )
{
  int rtc;

  rtc = fprintf( fp, "%d\n", item );
  if( rtc < 0 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  print double                                                              */
/*----------------------------------------------------------------------------*/
static int
print_double( double item, FILE *fp )
{
  int rtc;

  rtc = fprintf( fp, "%.16E\n", item );
  if( rtc < 0 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  print string                                                              */
/*----------------------------------------------------------------------------*/
static int
print_string( const char *item, FILE *fp )
{
  int rtc;

  rtc = fprintf( fp, "%s\n", item );
  if( rtc < 0 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  print HECMW_Comm                                                          */
/*----------------------------------------------------------------------------*/
static int
print_comm( HECMW_Comm item, FILE *fp )
{
  int rtc;

  /* rtc = fprintf( fp, "%d\n", (int)item ); */
  rtc = fprintf( fp, "%d\n", 0 );
  if( rtc < 0 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  print int (array)                                                         */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
print_int_ary( const int *ary, int n, int cols, FILE *fp )
{
  int rtc, i;

  if( n <= 0 )  return 0;

  for( i=0; i<n; i++ ) {
    rtc = fprintf( fp, "%d%c", ary[i], (i+1)%cols ? ' ' : '\n' );
    if( rtc < 0 ) {
      HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
  }
  if( n%cols ) {
    rtc = fprintf( fp, "\n" );
    if( rtc < 0 ) {
      HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  print double (array)                                                      */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
print_double_ary( const double *ary, int n, int cols, FILE *fp )
{
  int rtc, i;

  if( n <= 0 )  return 0;

  for( i=0; i<n; i++ ) {
    rtc = fprintf( fp, "%.16E%c", ary[i], (i+1)%cols ? ' ' : '\n' );
    if( rtc < 0 ) {
      HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
  }
  if( n%cols ) {
    rtc = fprintf( fp, "\n" );
    if( rtc < 0 ) {
      HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
  }

  return 0;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  print string (array)                                                      */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static int
print_string_ary( char **ary, int n, FILE *fp )
{
  int rtc, i;

  for( i=0; i<n; i++ ) {
    rtc = fprintf( fp, "%s\n", ary[i] );
    if( rtc < 0 ) {
      HECMW_set_error(HECMW_IO_E5004, "");
      return -1;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*  print header                                                              */
/*----------------------------------------------------------------------------*/
static int
print_header( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  int rtc;
  char header[HECMW_HEADER_LEN+1];

  strcpy(header, HEADER_STRING);
  rtc = sprintf( header+strlen(header), "%d", mesh->hecmw_flag_version );
  if( rtc < 0 ) {
    HECMW_set_error(HECMW_IO_E5004, "");
    return -1;
  }

  if( print_string( header, fp ) ) {
    return -1;
  } else {
    return 0;
  }
}

/*----------------------------------------------------------------------------*/
/*  print global information                                                  */
/*----------------------------------------------------------------------------*/
static int
print_global_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  int flag_header;

  /* hecmw_flag_adapt */
  if( print_int( mesh->hecmw_flag_adapt, fp ) ) {
    return -1;
  }

  /* hecmw_flag_initcon */
  if( print_int( mesh->hecmw_flag_initcon, fp ) ) {
    return -1;
  }

  /* hecmw_flag_parttype */
  if( print_int( mesh->hecmw_flag_parttype, fp ) ) {
    return -1;
  }

  /* hecmw_flag_partdepth */
  if( print_int( mesh->hecmw_flag_partdepth, fp ) ) {
    return -1;
  }

  /* hecmw_flag_version */
  if( print_int( mesh->hecmw_flag_version, fp ) ) {
    return -1;
  }

  /* gridfile */
  if( print_string( mesh->gridfile, fp ) ) {
    return -1;
  }

  /* hecmw_n_files */
  if( print_int( mesh->hecmw_n_file, fp ) ) {
    return -1;
  }

  /* files */
  if( mesh->hecmw_n_file > 0 ) {
    if( print_string_ary( mesh->files, mesh->hecmw_n_file, fp ) ) {
      return -1;
    }
  }

  if( strlen(mesh->header) > 0 ) {
    /* flag for header */
    flag_header = 1;
    if( print_int( flag_header, fp ) ) {
      return -1;
    }

    /* header */
    if( print_string( mesh->header, fp ) ) {
      return -1;
    }
  } else {
    /* flag for header */
    flag_header = 0;
    if( print_int( flag_header, fp ) ) {
      return -1;
    }
  }

  /* zero_temp */
  if( print_double( mesh->zero_temp, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print node information                                                    */
/*                                                                            */
/*============================================================================*/
static int
print_node_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* n_node */
  if( print_int( mesh->n_node, fp ) ) {
    return -1;
  }

  /* n_node_gross */
  if( print_int( mesh->n_node_gross, fp ) ) {
    return -1;
  }

  /* nn_internal */
  if( print_int( mesh->nn_internal, fp ) ) {
    return -1;
  }

  /* node_internal_list */
  if( mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_ELEMBASED
   || mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_UNKNOWN ) {
    if( print_int_ary( mesh->node_internal_list, mesh->nn_internal, COLS_INT_DEF, fp ) ) {
      return -1;
    }
  }

  /* node_ID */
  if( print_int_ary( mesh->node_ID, 2*mesh->n_node_gross, COLS_TWO, fp ) ) {
    return -1;
  }

  /* global_node_ID */
  if( print_int_ary( mesh->global_node_ID, mesh->n_node_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* node */
  if( print_double_ary( mesh->node, 3*mesh->n_node_gross, COLS_THREE, fp ) ) {
    return -1;
  }

  /* n_dof */
  if( print_int( mesh->n_dof, fp ) ) {
    return -1;
  }

  /* n_dof_grp */
  if( print_int( mesh->n_dof_grp, fp ) ) {
    return -1;
  }

  /* node_dof_index */
  if( print_int_ary( mesh->node_dof_index, mesh->n_dof_grp+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* node_dof_item */
  if( print_int_ary( mesh->node_dof_item, mesh->n_dof_grp, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  if( mesh->hecmw_flag_initcon ) {

    /* node_init_val_index */
    if( print_int_ary( mesh->node_init_val_index, mesh->n_node_gross+1, COLS_INT_DEF, fp ) ) {
      return -1;
    }

    /* node_init_val_item */
    if( print_double_ary( mesh->node_init_val_item, mesh->node_init_val_index[mesh->n_node_gross], COLS_DOUBLE_DEF, fp ) ) {
      return -1;
    }

  }

  /* node_val_index */
  /* node_val_item */

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print element information                                                 */
/*                                                                            */
/*============================================================================*/
static int
print_elem_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* n_elem */
  if( print_int( mesh->n_elem, fp ) ) {
    return -1;
  }

  /* n_elem_gross */
  if( print_int( mesh->n_elem_gross, fp ) ) {
    return -1;
  }

  /* ne_internal */
  if( print_int( mesh->ne_internal, fp ) ) {
    return -1;
  }

  /* elem_internal_list */
  if( mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_NODEBASED
   || mesh->hecmw_flag_parttype == HECMW_FLAG_PARTTYPE_UNKNOWN ) {
    if( print_int_ary( mesh->elem_internal_list, mesh->ne_internal, COLS_INT_DEF, fp ) ) {
      return -1;
    }
  }

  /* elem_ID */
  if( print_int_ary( mesh->elem_ID, 2*mesh->n_elem_gross, COLS_TWO, fp ) ) {
    return -1;
  }

  /* global_elem_ID */
  if( print_int_ary( mesh->global_elem_ID, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_type */
  if( print_int_ary( mesh->elem_type, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* n_elem_type */
  if( print_int( mesh->n_elem_type, fp ) ) {
    return -1;
  }

  /* elem_type_index */
  if( print_int_ary( mesh->elem_type_index, mesh->n_elem_type+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_type_item */
  if( print_int_ary( mesh->elem_type_item, mesh->n_elem_type, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_node_index */
  if( print_int_ary( mesh->elem_node_index, mesh->n_elem_gross+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_node_item */
  if( print_int_ary( mesh->elem_node_item, mesh->elem_node_index[mesh->n_elem_gross], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* section_ID */
  if( print_int_ary( mesh->section_ID, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_mat_ID_index */
  if( print_int_ary( mesh->elem_mat_ID_index, mesh->n_elem_gross+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* elem_mat_ID_item */
  if( print_int_ary( mesh->elem_mat_ID_item, mesh->elem_mat_ID_index[mesh->n_elem_gross], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* n_elem_mat_ID */
  if( print_int( mesh->n_elem_mat_ID, fp ) ) {
    return -1;
  }


  /* elem_mat_int_index */
  /* elem_mat_int_item */

  /* elem_val_index */
  /* elem_val_item */

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print domain information & communication tables                           */
/*                                                                            */
/*============================================================================*/
static int
print_comm_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* zero */
  if( print_int( mesh->zero, fp ) ) {
    return -1;
  }

  /* HECMW_COMM */
  if( print_comm( mesh->HECMW_COMM, fp ) ) {
    return -1;
  }

  /* PETOT */
  if( print_int( mesh->PETOT, fp ) ) {
    return -1;
  }

  /* PEsmpTOT */
  if( print_int( mesh->PEsmpTOT, fp ) ) {
    return -1;
  }

  /* my_rank */
  if( print_int( mesh->my_rank, fp ) ) {
    return -1;
  }

  /* errnof */
  if( print_int( mesh->errnof, fp ) ) {
    return -1;
  }

  /* n_subdomain */
  if( print_int( mesh->n_subdomain, fp ) ) {
    return -1;
  }

  /* n_neighbor_pd */
  if( print_int( mesh->n_neighbor_pe, fp ) ) {
    return -1;
  }

  if( mesh->n_neighbor_pe == 0 )  return 0;

  /* neighbor_pe */
  if( print_int_ary( mesh->neighbor_pe, mesh->n_neighbor_pe, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* import_index */
  if( print_int_ary( mesh->import_index, mesh->n_neighbor_pe+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* import_item */
  if( print_int_ary( mesh->import_item, mesh->import_index[mesh->n_neighbor_pe], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* export_index */
  if( print_int_ary( mesh->export_index, mesh->n_neighbor_pe+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* export_item */
  if( print_int_ary( mesh->export_item, mesh->export_index[mesh->n_neighbor_pe], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* shared_index */
  if( print_int_ary( mesh->shared_index, mesh->n_neighbor_pe+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* shared_item */
  if( print_int_ary( mesh->shared_item, mesh->shared_index[mesh->n_neighbor_pe], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print adaptation information                                              */
/*                                                                            */
/*============================================================================*/
static int
print_adapt_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  if( mesh->hecmw_flag_adapt == 0 )  return 0;

  /* coarse_grid_level */
  if( print_int( mesh->coarse_grid_level, fp ) ) {
    return -1;
  }

  /* n_adapt */
  if( print_int( mesh->n_adapt, fp ) ) {
    return -1;
  }

  /* when_i_was_refined_node */
  if( print_int_ary( mesh->when_i_was_refined_node, mesh->n_node_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* when_i_was_refined_elem */
  if( print_int_ary( mesh->when_i_was_refined_elem, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* adapt_parent_type */
  if( print_int_ary( mesh->adapt_parent_type, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* adapt_type */
  if( print_int_ary( mesh->adapt_type, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* adapt_level */
  if( print_int_ary( mesh->adapt_level, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* adapt_parent */
  if( print_int_ary( mesh->adapt_parent, 2*mesh->n_elem_gross, COLS_TWO, fp ) ) {
    return -1;
  }

  /* adapt_children_index */
  if( print_int_ary( mesh->adapt_children_index, mesh->n_elem_gross+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* adapt_children_item */
  if( print_int_ary( mesh->adapt_children_item, 2*mesh->adapt_children_index[mesh->n_elem_gross], COLS_TWO, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print section information                                                 */
/*                                                                            */
/*============================================================================*/
static int
print_section_info( const struct hecmwST_section *sect, FILE *fp )
{
  /* n_sect */
  if( print_int( sect->n_sect, fp ) ) {
    return -1;
  }

  if( sect->n_sect == 0 )  return 0;

  /* sect_type */
  if( print_int_ary( sect->sect_type, sect->n_sect, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_opt */
  if( print_int_ary( sect->sect_opt, sect->n_sect, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_mat_ID_index */
  if( print_int_ary( sect->sect_mat_ID_index, sect->n_sect+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_mat_ID_item */
  if( print_int_ary( sect->sect_mat_ID_item, sect->sect_mat_ID_index[sect->n_sect], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_I_index */
  if( print_int_ary( sect->sect_I_index, sect->n_sect+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_I_item */
  if( print_int_ary( sect->sect_I_item, sect->sect_I_index[sect->n_sect], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_R_index */
  if( print_int_ary( sect->sect_R_index, sect->n_sect+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* sect_R_item */
  if( print_double_ary( sect->sect_R_item, sect->sect_R_index[sect->n_sect], COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print material information                                                */
/*                                                                            */
/*============================================================================*/
static int
print_material_info( struct hecmwST_material *mat, FILE *fp )
{
  /* n_mat */
  if( print_int( mat->n_mat, fp ) ) {
    return -1;
  }

  if( mat->n_mat == 0 )  return 0;

  /* n_mat_item */
  if( print_int( mat->n_mat_item, fp ) ) {
    return -1;
  }

  /* n_mat_subitem */
  if( print_int( mat->n_mat_subitem, fp ) ) {
    return -1;
  }

  /* n_mat_table */
  if( print_int( mat->n_mat_table, fp ) ) {
    return -1;
  }

  /* mat_name */
  if( print_string_ary( mat->mat_name, mat->n_mat, fp ) ) {
    return -1;
  }

  /* mat_item_index */
  if( print_int_ary( mat->mat_item_index, mat->n_mat+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mat_subitem_index */
  if( print_int_ary( mat->mat_subitem_index, mat->n_mat_item+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mat_table_index */
  if( print_int_ary( mat->mat_table_index, mat->n_mat_subitem+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mat_val */
  if( print_double_ary( mat->mat_val, mat->n_mat_table, COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  /* mat_temp */
  if( print_double_ary( mat->mat_temp, mat->n_mat_table, COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print MPC group information                                               */
/*                                                                            */
/*============================================================================*/
static int
print_mpc_info( const struct hecmwST_mpc *mpc, FILE *fp )
{
  /* n_mpc */
  if( print_int( mpc->n_mpc, fp ) ) {
    return -1;
  }

  if( mpc->n_mpc == 0 )  return 0;

  /* mpc_index */
  if( print_int_ary( mpc->mpc_index, mpc->n_mpc+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mpc_item */
  if( print_int_ary( mpc->mpc_item, mpc->mpc_index[mpc->n_mpc], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mpc_dof */
  if( print_int_ary( mpc->mpc_dof, mpc->mpc_index[mpc->n_mpc], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* mpc_val */
  if( print_double_ary( mpc->mpc_val, mpc->mpc_index[mpc->n_mpc], COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  /* mpc_const */
  if( print_double_ary( mpc->mpc_const, mpc->n_mpc, COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print amplitude information                                               */
/*                                                                            */
/*============================================================================*/
static int
print_amp_info( const struct hecmwST_amplitude *amp, FILE *fp )
{
  /* n_amp */
  if( print_int( amp->n_amp , fp ) ) {
    return -1;
  }

  if( amp->n_amp == 0 )  return 0;

  /* amp_name */
  if( print_string_ary( amp->amp_name, amp->n_amp, fp ) ) {
    return -1;
  }

  /* amp_type_definition */
  if( print_int_ary( amp->amp_type_definition, amp->n_amp, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* amp_type_time */
  if( print_int_ary( amp->amp_type_time, amp->n_amp, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* amp_type_value */
  if( print_int_ary( amp->amp_type_value, amp->n_amp, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* amp_index */
  if( print_int_ary( amp->amp_index, amp->n_amp+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* amp_val */
  if( print_double_ary( amp->amp_val, amp->amp_index[amp->n_amp], COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  /* amp_table */
  if( print_double_ary( amp->amp_table, amp->amp_index[amp->n_amp], COLS_DOUBLE_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print node group information                                              */
/*                                                                            */
/*============================================================================*/
static int
print_node_grp_info( const struct hecmwST_node_grp *grp, FILE *fp )
{
  /* n_grp */
  if( print_int( grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 )  return 0;

  /* grp_name */
  if( print_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( print_int_ary( grp->grp_index, grp->n_grp+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( print_int_ary( grp->grp_item, grp->grp_index[grp->n_grp], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print element group information                                           */
/*                                                                            */
/*============================================================================*/
static int
print_elem_grp_info( const struct hecmwST_elem_grp *grp, FILE *fp )
{
  /* n_grp */
  if( print_int( grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 )  return 0;

  /* grp_name */
  if( print_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( print_int_ary( grp->grp_index, grp->n_grp+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( print_int_ary( grp->grp_item, grp->grp_index[grp->n_grp], COLS_INT_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print surface group information                                           */
/*                                                                            */
/*============================================================================*/
static int
print_surf_grp_info( const struct hecmwST_surf_grp *grp, FILE *fp )
{
  /* n_grp */
  if( print_int( grp->n_grp, fp ) ) {
    return -1;
  }

  if( grp->n_grp == 0 )  return 0;

  /* grp_name */
  if( print_string_ary( grp->grp_name, grp->n_grp, fp ) ) {
    return -1;
  }

  /* grp_index */
  if( print_int_ary( grp->grp_index, grp->n_grp+1, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* grp_item */
  if( print_int_ary( grp->grp_item, grp->grp_index[grp->n_grp]*2, COLS_TWO, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print refinement information                                              */
/*                                                                            */
/*============================================================================*/
static int
print_refine_info( const struct hecmwST_local_mesh *mesh, FILE *fp )
{
  /* number of refinement performed */
  if( print_int( mesh->n_refine, fp ) ) {
    return -1;
  }

  if( mesh->n_refine == 0 || mesh->n_subdomain == 1 )  return 0;

  if( mesh->n_node_gross > mesh->nn_internal ) {
    /* node_old2new */
    if( print_int_ary( mesh->node_old2new, mesh->n_node_gross, COLS_INT_DEF, fp ) ) {
      return -1;
    }

    /* node_new2old */
    if( print_int_ary( mesh->node_new2old, mesh->n_node_gross, COLS_INT_DEF, fp ) ) {
      return -1;
    }
  }

  if( mesh->n_elem_gross > mesh->n_elem ) {
    /* elem_old2new */
    if( print_int_ary( mesh->elem_old2new, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
      return -1;
    }

    /* elem_new2old */
    if( print_int_ary( mesh->elem_new2old, mesh->n_elem_gross, COLS_INT_DEF, fp ) ) {
      return -1;
    }
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print contact information                                                 */
/*                                                                            */
/*============================================================================*/
static int
print_contact_info( const struct hecmwST_contact_pair *cpair, FILE *fp )
{
  /* n_pair */
  if( print_int( cpair->n_pair, fp ) ) {
    return -1;
  }

  if( cpair->n_pair == 0 )  return 0;

  /* name */
  if( print_string_ary( cpair->name, cpair->n_pair, fp ) ) {
    return -1;
  }

  /* type */
  if( print_int_ary( cpair->type, cpair->n_pair, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* slave_grp_id */
  if( print_int_ary( cpair->slave_grp_id, cpair->n_pair, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  /* master_grp_id */
  if( print_int_ary( cpair->master_grp_id, cpair->n_pair, COLS_INT_DEF, fp ) ) {
    return -1;
  }

  return 0;
}

/*============================================================================*/
/*                                                                            */
/*  print HEC-MW distributed mesh format data                                 */
/*                                                                            */
/*============================================================================*/
extern int
HECMW_put_dist_mesh( const struct hecmwST_local_mesh *mesh, char *fname )
{
  FILE *fp;

  if(mesh == NULL) return 0;
  if(fname == NULL) {
    HECMW_set_error(HECMW_IO_E5001, "Filename is NULL)");
	return -1;
  }

  if(HECMW_ctrl_is_subdir()) {
    if(HECMW_ctrl_make_subdir(fname)) return 0;
  }

  /* open file */
  if( ( fp = fopen( fname, "w" ) ) == NULL ) {
    HECMW_set_error(HECMW_IO_E5001, "File: %s, %s", fname, strerror( errno ) );
	return -1;
  }

  /* header */
  if( print_header( mesh, fp ) ) {
	return -1;
  }

  /* global info. */
  if( print_global_info( mesh, fp ) ) {
	return -1;
  }

  /* node info. */
  if( print_node_info( mesh, fp ) ) {
	return -1;
  }

  /* element info. */
  if( print_elem_info( mesh, fp ) ) {
	return -1;
  }

  /* domain info & communication table */
  if( print_comm_info( mesh, fp ) ) {
	return -1;
  }

  /* adaptation info. */
  if( print_adapt_info( mesh, fp ) ) {
	return -1;
  }

  /* section info. */
  if( print_section_info( mesh->section, fp ) ) {
	return -1;
  }

  /* material info. */
  if( print_material_info( mesh->material, fp ) ) {
	return -1;
  }

  /* MPC group info. */
  if( print_mpc_info( mesh->mpc, fp ) ) {
	return -1;
  }

  /* amplitude info. */
  if( print_amp_info( mesh->amp, fp ) ) {
	return -1;
  }

  /* node group info. */
  if( print_node_grp_info( mesh->node_group, fp) ) {
	return -1;
  }

  /* element group info. */
  if( print_elem_grp_info( mesh->elem_group, fp ) ) {
	return -1;
  }

  /* surface group info */
  if( print_surf_grp_info( mesh->surf_group, fp ) ) {
	return -1;
  }

  /* refinement info. */
  if( print_refine_info( mesh, fp ) ) {
	return -1;
  }

  /* contact info */
  if( print_contact_info( mesh->contact_pair, fp ) ) {
	return -1;
  }

  /* close file */
  if( fclose( fp ) ) {
    HECMW_set_error(HECMW_IO_E5002, HECMW_strmsg( errno ) );
	return -1;
  }

  return 0;
}
