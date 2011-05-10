/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.1                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_config.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_result_copy_f2c.h"

static int n_node;
static int n_elem;
static struct hecmwST_result_data *result;

/*-----------------------------------------------------------------------------
 * SetFunc
 */

static int
set_nn_component(void *src)
{
  result->nn_component = *((int *)src);
  return 0;
}


static int
set_ne_component(void *src)
{
  result->ne_component = *((int *)src);
  return 0;
}


static int
set_nn_dof(void *src)
{
  int size;

  if(result->nn_component <= 0) return 0;
  size = sizeof(*result->nn_dof)*result->nn_component;
  result->nn_dof = HECMW_malloc(size);
  if(result->nn_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  memcpy(result->nn_dof, src, size);
  return 0;
}


static int
set_ne_dof(void *src)
{
  int size;

  if(result->ne_component <= 0) return 0;
  size = sizeof(*result->ne_dof)*result->ne_component;
  result->ne_dof = HECMW_malloc(size);
  if(result->ne_dof == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  memcpy(result->ne_dof, src, size);
  return 0;
}


static int
set_node_label(void *src)
{
  int i;

  if(result->nn_component <= 0) return 0;

  result->node_label = HECMW_malloc(sizeof(*result->node_label)*result->nn_component);
  if(result->node_label == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  for(i=0; i < result->nn_component; i++) {
    char *src_point = (char *)src + HECMW_NAME_LEN * i;
    result->node_label[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
  }

  return 0;
}


static int
set_elem_label(void *src)
{
  int i;

  if(result->ne_component <= 0) return 0;

  result->elem_label = HECMW_malloc(sizeof(*result->elem_label)*result->ne_component);
  if(result->elem_label == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  for(i=0; i < result->ne_component; i++) {
    char *src_point = (char *)src + HECMW_NAME_LEN * i;
    result->elem_label[i] = HECMW_strcpy_f2c(src_point, HECMW_NAME_LEN);
  }

  return 0;
}


static int
set_node_val_item(void *src)
{
  int i,size;
  int n = 0;

  if(result->nn_component <= 0) return 0;
  for(i=0; i < result->nn_component; i++) {
	  n += result->nn_dof[i];
  }
  size = sizeof(*result->node_val_item)*n*n_node;
  result->node_val_item = HECMW_malloc(size);
  if(result->node_val_item == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  memcpy(result->node_val_item, src, size);
  return 0;
}


static int
set_elem_val_item(void *src)
{
  int i,size;
  int n = 0;

  if(result->ne_component <= 0) return 0;
  for(i=0; i < result->ne_component; i++) {
	  n += result->ne_dof[i];
  }
  size = sizeof(*result->elem_val_item)*n*n_elem;
  result->elem_val_item = HECMW_malloc(size);
  if(result->elem_val_item == NULL) {
    HECMW_set_error(errno, "");
    return -1;
  }
  memcpy(result->elem_val_item, src, size);
  return 0;
}


/*-----------------------------------------------------------------------------
 * SetFunc table
 */

typedef int (*SetFunc)(void *);

static struct func_table {
  char *struct_name;
  char *var_name;
  SetFunc set_func;
} functions[] = {
/*  { Struct name, Variable name, memcpy function } */
  {"hecmwST_result_data","nn_component", set_nn_component },
  {"hecmwST_result_data","ne_component", set_ne_component },
  {"hecmwST_result_data","nn_dof",       set_nn_dof       },
  {"hecmwST_result_data","ne_dof",       set_ne_dof       },
  {"hecmwST_result_data","node_label",   set_node_label   },
  {"hecmwST_result_data","elem_label",   set_elem_label   },
  {"hecmwST_result_data","node_val_item",set_node_val_item},
  {"hecmwST_result_data","elem_val_item",set_elem_val_item},
};

static const int NFUNC = sizeof(functions) / sizeof(functions[0]);



static SetFunc
get_set_func(char *struct_name, char *var_name)
{
  int i;

  for(i=0; i < NFUNC; i++) {
    if(strcmp(functions[i].struct_name, struct_name) == 0 
    && strcmp(functions[i].var_name, var_name) == 0) {
      return functions[i].set_func;
    }
  }
  return NULL;
}

/*----------------------------------------------------------------------------*/

extern int
HECMW_result_copy_f2c_init( struct hecmwST_result_data *result_data, int nnode, int nelem )
{
  memset(result_data, 0, sizeof(*result_data));
  result = result_data;
  n_node = nnode;
  n_elem = nelem;
  return 0;
}

/*----------------------------------------------------------------------------*/

extern int
HECMW_result_copy_f2c_finalize( void )
{
  result = NULL;
  return 0;
}

/*----------------------------------------------------------------------------*/


extern void
hecmw_result_copy_f2c_set_if( char *struct_name, char *var_name,
                              void *src, int *err, int slen, int vlen )
{
  SetFunc func;
  char sname[HECMW_NAME_LEN+1];
  char vname[HECMW_NAME_LEN+1];

  *err = 1;

  if( result == NULL ) {
    HECMW_set_error(HECMW_ALL_E0102,
               "hecmw_result_copy_f2c_set_if(): 'result' has not initialized yet" );
    return;
  }
  if( struct_name == NULL ) {
    HECMW_set_error(HECMW_ALL_E0101,
               "hecmw_result_copy_f2c_set_if(): 'sname' is NULL" );
    return;
  }
  if( var_name == NULL ) {
    HECMW_set_error(HECMW_ALL_E0101,
               "hecmw_result_copy_f2c_set_if(): 'vname' is NULL" );
    return;
  }

  if( HECMW_strcpy_f2c_r( struct_name, slen, sname, sizeof(sname) ) == NULL ) {
    return;
  }

  if( HECMW_strcpy_f2c_r( var_name, vlen, vname, sizeof(vname) ) == NULL ) {
    return;
  }

  func = get_set_func( sname, vname );
  if( func == NULL ) {
    HECMW_set_error(HECMW_ALL_E0102, "hecmw_result_copy_f2c_set_if(): SetFunc not found" );
    return;
  }

  if((*func)(src)) {
	  return;
  }

  *err = 0;  /* no error */
}



extern void
hecmw_result_copy_f2c_set_if_( char *struct_name, char *var_name,
                               void *src, int *err, int slen, int vlen )
{
  hecmw_result_copy_f2c_set_if( struct_name, var_name, src, err, slen, vlen );
}



extern void
hecmw_result_copy_f2c_set_if__( char *struct_name, char *var_name,
                                void *src, int *err, int slen, int vlen )
{
  hecmw_result_copy_f2c_set_if( struct_name, var_name, src, err, slen, vlen );
}



extern void
HECMW_RESULT_COPY_F2C_SET_IF( char *struct_name, char *var_name,
                                void *src, int *err, int slen, int vlen )
{
  hecmw_result_copy_f2c_set_if( struct_name, var_name, src, err, slen, vlen );
}

