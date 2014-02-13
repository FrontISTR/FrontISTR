/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
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
#include "hecmw_common_define.h"
#include "hecmw_struct.h"
#include "hecmw_result.h"
#include "hecmw_util.h"
#include "hecmw_etype.h"
#include "hecmw_ucd_print.h"

static int conv_index_ucd2hec_rod1[] = {
  0, 1
};

static int conv_index_ucd2hec_rod2[] = {
  0, -1, 2
};

static int conv_index_ucd2hec_tri1[] = {
  0, 1, 2
};

static int conv_index_ucd2hec_tri2[] = {
  0, 1, 2, -1, -1, -1
};

static int conv_index_ucd2hec_qua1[] = {
  0, 1, 2, 3
};

static int conv_index_ucd2hec_qua2[] = {
  0, 1, 2, 3, -1, -1, -1, -1
};

static int conv_index_ucd2hec_tet1[] = {
  0, 3, 2, 1
};

static int conv_index_ucd2hec_tet2[] = {
  0, 3, 2, 1, -1, -1, -1, -1, -1, -1
};

static int conv_index_ucd2hec_pri1[] = {
  3, 4, 5, 0, 1, 2
};

static int conv_index_ucd2hec_pri2[] = {
  3, 4, 5, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

static int conv_index_ucd2hec_hex1[] = {
  4, 5, 6, 7, 0, 1, 2, 3
};

static int conv_index_ucd2hec_hex2[] = {
  4, 5, 6, 7, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

static int conv_index_ucd2hec_pyr1[] = {
  4, 0, 1, 2, 3
};

static int conv_index_ucd2hec_pyr2[] = {
  4, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1
};

static int conv_index_ucd2hec_mst1[] = {
  0, 1, 2, 3
};

static int conv_index_ucd2hec_mst2[] = {
  0, 1, 2, 3, -1, -1, -1, -1
};

static int conv_index_ucd2hec_msq1[] = {
  0, 1, 2, 3, 4
};

static int conv_index_ucd2hec_msq2[] = {
  0, 1, 2, 3, 4, -1, -1, -1, -1
};

static int conv_index_ucd2hec_jtt1[] = {
  3, 4, 5, 0, 1, 2
};

static int conv_index_ucd2hec_jtt2[] = {
  3, 4, 5, 0, 1, 2, -1, -1, -1, -1, -1, -1
};

static int conv_index_ucd2hec_jtq1[] = {
  4, 5, 6, 7, 0, 1, 2, 3
};

static int conv_index_ucd2hec_jtq2[] = {
  4, 5, 6, 7, 0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1
};

#if 0
static int conv_index_ucd2hec_bem1[] = {
  0, 1
};

static int conv_index_ucd2hec_bem2[] = {
  0, -1, 1
};
#endif

static int conv_index_ucd2hec_sht1[] = {
  0, 1, 2
};

static int conv_index_ucd2hec_sht2[] = {
  0, 1, 2, -1, -1, -1
};

static int conv_index_ucd2hec_shq1[] = {
  0, 1, 2, 3
};

static int conv_index_ucd2hec_shq2[] = {
  0, 1, 2, 3, -1, -1, -1, -1
};

static int conv_index_ucd2hec_ln[] = {
  0, 1
};

static int
ucd_print( const struct hecmwST_local_mesh *mesh,
           const struct hecmwST_result_data *result,
           const char *ofname,
           int flag_oldUCD )
{
  int nn_item=0, ne_item=0;
  int node_index;
  int i, j;
  FILE *fp;

  if( mesh == NULL ) {
    HECMW_print_msg( HECMW_LOG_WARN, HECMW_ALL_E0101, "mesh is not set" );
    return 1;
  }
  if( result == NULL ) {
    HECMW_print_msg( HECMW_LOG_WARN, HECMW_ALL_E0101, "result data is not set" );
    return 1;
  }
  if( ofname == NULL ) {
    HECMW_print_msg( HECMW_LOG_WARN, HECMW_ALL_E0101, "output file name is not set" );
    return 1;
  }

  if( ( fp = fopen( ofname, "w" ) ) == NULL ) {
    HECMW_set_error(errno, "");
	return -1;
  }

  for( nn_item=0, i=0; i<result->nn_component; i++ ) {
    nn_item += result->nn_dof[i];
  }
  for( ne_item=0, i=0; i<result->ne_component; i++ ) {
    ne_item += result->ne_dof[i];
  }

  if (flag_oldUCD) {
    fprintf( fp, "%d %d %d %d 0\n",
             mesh->n_node, mesh->n_elem, nn_item, ne_item );
  } else {
    /* comment part */
    fprintf( fp, "# File Format : multi-step UCD data for unstructured mesh\n" );
    fprintf( fp, "# created by HEC-MW ( %s )\n", HECMW_get_date( ) );

    /* header part */
    fprintf( fp, "%d\n", 1 );
    fprintf( fp, "data\n" );
    fprintf( fp, "step%d\n", 1 );
    fprintf( fp, "%d %d\n", mesh->n_node, mesh->n_elem );
  }

  /* nodal information */
  for( i=0; i<mesh->n_node; i++ ) {
    fprintf( fp, "%d %.7lE %.7lE %.7lE\n", i+1, mesh->node[3*i], mesh->node[3*i+1], mesh->node[3*i+2] );
  }

  /* element information */
  for( i=0; i<mesh->n_elem; i++ ) {
    switch( mesh->elem_type[i] ) {
      case HECMW_ETYPE_ROD1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_ROD1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_ROD1; j++ ) {
          if( conv_index_ucd2hec_rod1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_rod1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_ROD2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_ROD1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_ROD2; j++ ) {
          if( conv_index_ucd2hec_rod2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_rod2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_TRI1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_TRI1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_TRI1; j++ ) {
          if( conv_index_ucd2hec_tri1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_tri1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_TRI2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_TRI1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_TRI2; j++ ) {
          if( conv_index_ucd2hec_tri2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_tri2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_QUA1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_QUA1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_QUA1; j++ ) {
          if( conv_index_ucd2hec_qua1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_qua1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_QUA2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_QUA1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_QUA2; j++ ) {
          if( conv_index_ucd2hec_qua2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_qua2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_TET1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_TET1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_TET1; j++ ) {
          if( conv_index_ucd2hec_tet1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_tet1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_TET2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_TET1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_TET2; j++ ) {
          if( conv_index_ucd2hec_tet2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_tet2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_PRI1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_PRI1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_PRI1; j++ ) {
          if( conv_index_ucd2hec_pri1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_pri1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_PRI2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_PRI1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_PRI2; j++ ) {
          if( conv_index_ucd2hec_pri2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_pri2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_HEX1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_HEX1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_HEX1; j++ ) {
          if( conv_index_ucd2hec_hex1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_hex1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_HEX2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_HEX1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_HEX2; j++ ) {
          if( conv_index_ucd2hec_hex2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_hex2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_PYR1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_PYR1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_PYR1; j++ ) {
          if( conv_index_ucd2hec_pyr1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_pyr1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_PYR2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_PYR1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_PYR2; j++ ) {
          if( conv_index_ucd2hec_pyr2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_pyr2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_MST1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_MST1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_MST1; j++ ) {
          if( conv_index_ucd2hec_mst1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_mst1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_MST2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_MST1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_MST2; j++ ) {
          if( conv_index_ucd2hec_mst2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_mst2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_MSQ1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_MSQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_MSQ1; j++ ) {
          if( conv_index_ucd2hec_msq1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_msq1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_MSQ2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_MSQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_MSQ2; j++ ) {
          if( conv_index_ucd2hec_msq2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_msq2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_JTT1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_JTT1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_JTT1; j++ ) {
          if( conv_index_ucd2hec_jtt1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_jtt1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_JTT2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_JTT1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_JTT2; j++ ) {
          if( conv_index_ucd2hec_jtt2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_jtt2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_JTQ1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_JTQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_JTQ1; j++ ) {
          if( conv_index_ucd2hec_jtq1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_jtq1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_JTQ2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_JTQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_JTQ2; j++ ) {
          if( conv_index_ucd2hec_jtq2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_jtq2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_SHT1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_SHT1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_SHT1; j++ ) {
          if( conv_index_ucd2hec_sht1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_sht1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_SHT2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_SHT1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_SHT2; j++ ) {
          if( conv_index_ucd2hec_sht2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_sht2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_SHQ1:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_SHQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_SHQ1; j++ ) {
          if( conv_index_ucd2hec_shq1[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_shq1[j]] );
          }
        }
        break;

      case HECMW_ETYPE_SHQ2:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_SHQ1) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_SHQ2; j++ ) {
          if( conv_index_ucd2hec_shq2[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_shq2[j]] );
          }
        }
        break;

      case HECMW_ETYPE_LN11:
      case HECMW_ETYPE_LN12:
      case HECMW_ETYPE_LN13:
      case HECMW_ETYPE_LN14:
      case HECMW_ETYPE_LN15:
      case HECMW_ETYPE_LN16:
      case HECMW_ETYPE_LN21:
      case HECMW_ETYPE_LN22:
      case HECMW_ETYPE_LN23:
      case HECMW_ETYPE_LN24:
      case HECMW_ETYPE_LN25:
      case HECMW_ETYPE_LN26:
      case HECMW_ETYPE_LN31:
      case HECMW_ETYPE_LN32:
      case HECMW_ETYPE_LN33:
      case HECMW_ETYPE_LN34:
      case HECMW_ETYPE_LN35:
      case HECMW_ETYPE_LN36:
      case HECMW_ETYPE_LN41:
      case HECMW_ETYPE_LN42:
      case HECMW_ETYPE_LN43:
      case HECMW_ETYPE_LN44:
      case HECMW_ETYPE_LN45:
      case HECMW_ETYPE_LN46:
      case HECMW_ETYPE_LN51:
      case HECMW_ETYPE_LN52:
      case HECMW_ETYPE_LN53:
      case HECMW_ETYPE_LN54:
      case HECMW_ETYPE_LN55:
      case HECMW_ETYPE_LN56:
      case HECMW_ETYPE_LN61:
      case HECMW_ETYPE_LN62:
      case HECMW_ETYPE_LN63:
      case HECMW_ETYPE_LN64:
      case HECMW_ETYPE_LN65:
      case HECMW_ETYPE_LN66:
        fprintf( fp, "%d %d %s", i+1, 0, HECMW_get_ucd_label(HECMW_ETYPE_LN11) );
        node_index = mesh->elem_node_index[i];
        for( j=0; j<HECMW_MAX_NODE_LN11; j++ ) {
          if( conv_index_ucd2hec_ln[j] >= 0 ) {
            fprintf( fp, " %d", mesh->elem_node_item[node_index+conv_index_ucd2hec_ln[j]] );
          }
        }
        break;

      default:
        return -1;
    }
    fprintf( fp, "\n" );
  }

  /* data part */
  if (!flag_oldUCD) {
    fprintf( fp, "%d %d\n", nn_item, ne_item );
  }

  if( result->nn_component > 0 ) {
    fprintf( fp, "%d", result->nn_component );
    for( i=0; i<result->nn_component; i++ ) {
      fprintf( fp, " %d", result->nn_dof[i] );
    }
    fprintf( fp, "\n" );

    for( i=0; i<result->nn_component; i++ ) {
      fprintf( fp, "%s, unit_unknown\n", result->node_label[i] );
    }

    for( i=0; i<mesh->n_node; i++ ) {
      fprintf( fp, "%d", i+1 );
      for( j=0; j<nn_item; j++ ) {
        fprintf( fp, " %.7lE", result->node_val_item[nn_item*i+j] );
      }
      fprintf( fp, "\n" );
    }
  }  

  if( result->ne_component > 0 ) {
    fprintf( fp, "%d", result->ne_component );
    for( i=0; i<result->ne_component; i++ ) {
      fprintf( fp, " %d", result->ne_dof[i] );
    }
    fprintf( fp, "\n" );

    for( i=0; i<result->ne_component; i++ ) {
      fprintf( fp, "%s, unit_unknown\n", result->elem_label[i] );
    }

    for( i=0; i<mesh->n_elem; i++ ) {
      fprintf( fp, "%d", i+1 );
      for( j=0; j<ne_item; j++ ) {
        fprintf( fp, " %.7lE", result->elem_val_item[ne_item*i+j] );
      }
      fprintf( fp, "\n" );
    }
  }  

  fclose(fp);
  
  return 0;
}

extern int
HECMW_ucd_print( const struct hecmwST_local_mesh *mesh,
                 const struct hecmwST_result_data *result,
                 const char *ofname )
{
  int flag_oldUCD = 0;
  return ucd_print( mesh, result, ofname, flag_oldUCD );
}

extern int
HECMW_ucd_legacy_print( const struct hecmwST_local_mesh *mesh,
                        const struct hecmwST_result_data *result,
                        const char *ofname )
{
  int flag_oldUCD = 1;
  return ucd_print( mesh, result, ofname, flag_oldUCD );
}
