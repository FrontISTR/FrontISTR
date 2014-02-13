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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#include "hecmw_msgno.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"

#include "hecmw_part_define.h"
#include "hecmw_mesh_hash_sort.h"


#define EDGE_INC_FACTOR (1.2)


#define TSUF_INC_FACTOR (1.1)


#define QSUF_INC_FACTOR (1.1)


struct hecmw_mesh_hash_link {

    int key;

    struct hecmw_mesh_hash_link *next;
};


static int n_edge = 0;


static int *__edge_node = NULL;


static struct hecmw_mesh_hash_link **e_hash_tbl = NULL;


static unsigned long int e_hash_size;


static unsigned long int e_buf_size;


static int n_tsuf = 0;


static int *__tsuf_node = NULL;


static struct hecmw_mesh_hash_link **t_hash_tbl = NULL;


static unsigned long int t_hash_size;


static unsigned long int t_buf_size;


static int n_qsuf = 0;


static int *__qsuf_node = NULL;


static struct hecmw_mesh_hash_link **q_hash_tbl = NULL;


static unsigned long int q_hash_size;


static unsigned long int q_buf_size;


struct hecmw_edge_node {

    int *node1;

    int *node2;
};


struct hecmw_tsuf_node {

    int *node1;

    int *node2;

    int *node3;
};


struct hecmw_qsuf_node {

    int *node1;

    int *node2;

    int *node3;

    int *node4;
};


static struct hecmw_edge_node *edge_node = NULL;


static struct hecmw_tsuf_node *tsuf_node = NULL;


static struct hecmw_qsuf_node *qsuf_node = NULL;


/*================================================================================================*/

extern int
HECMW_mesh_hsort_edge_init( int n_node, int n_elem )
{
    int size;
    int i;

    if( n_node <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_node=%d", n_node );
        goto error;
    }
    if( n_elem <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_elem=%d", n_elem );
        goto error;
    }

    if( n_node < 1000000 && n_elem < 1000000 ) {
        e_hash_size = ( n_node > n_elem ) ?    n_node :    n_elem;
        e_buf_size  = ( n_node > n_elem ) ? 10*n_node : 10*n_elem;
    } else {
        e_hash_size = ( n_node > n_elem ) ?    n_node :    n_elem;
        e_buf_size  = ( n_node > n_elem ) ?    n_node :    n_elem;
    }

    edge_node = (struct hecmw_edge_node *)HECMW_malloc( sizeof(struct hecmw_edge_node) );
    if( edge_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        edge_node->node1 = NULL;
        edge_node->node2 = NULL;
    }
    edge_node->node1 = (int *)HECMW_malloc( sizeof(int)*e_buf_size );
    if( edge_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    edge_node->node2 = (int *)HECMW_malloc( sizeof(int)*e_buf_size );
    if( edge_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    size = sizeof(struct hecmw_mesh_hash_link *) * e_hash_size;
    e_hash_tbl = (struct hecmw_mesh_hash_link **)HECMW_malloc( size );
    if( e_hash_tbl == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        for( i=0; i<e_hash_size; i++ ) {
            e_hash_tbl[i] = NULL;
        }
    }

    n_edge = 0;

    return 0;

error:
    HECMW_mesh_hsort_edge_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_tsuf_init( int n_node, int n_elem )
{
    int size;
    int i;

    if( n_node <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_node=%d", n_node );
        goto error;
    }
    if( n_elem <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_elem=%d", n_elem );
        goto error;
    }

    if( n_node < 1000000 && n_elem < 1000000 ) {
        t_hash_size = ( n_node > n_elem ) ?   n_node :   n_elem;
        t_buf_size  = ( n_node > n_elem ) ? 4*n_node : 4*n_elem;
    } else {
        t_hash_size = ( n_node > n_elem ) ?   n_node :   n_elem;
        t_buf_size  = ( n_node > n_elem ) ?   n_node :   n_elem;
    }

    tsuf_node = (struct hecmw_tsuf_node *)HECMW_malloc( sizeof(struct hecmw_tsuf_node) );
    if( tsuf_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        tsuf_node->node1 = NULL;
        tsuf_node->node2 = NULL;
        tsuf_node->node3 = NULL;
    }
    tsuf_node->node1 = (int *)HECMW_malloc( sizeof(int)*t_buf_size );
    if( tsuf_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    tsuf_node->node2 = (int *)HECMW_malloc( sizeof(int)*t_buf_size );
    if( tsuf_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    tsuf_node->node3 = (int *)HECMW_malloc( sizeof(int)*t_buf_size );
    if( tsuf_node->node3 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    size = sizeof(struct hecmw_mesh_hash_link *) * t_hash_size;
    t_hash_tbl = (struct hecmw_mesh_hash_link **)HECMW_malloc( size );
    if( t_hash_tbl == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        for( i=0; i<t_hash_size; i++ ) {
            t_hash_tbl[i] = NULL;
        }
    }

    n_tsuf = 0;

    return 0;

error:
    HECMW_mesh_hsort_tsuf_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_qsuf_init( int n_node, int n_elem )
{
    int size;
    int i;

    if( n_node <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_node=%d", n_node );
        goto error;
    }
    if( n_elem <= 0 ) {
        HECMW_set_error( HECMW_PART_E_INV_ARG, "n_elem=%d", n_elem );
        goto error;
    }

    if( n_node < 1000000 && n_elem < 1000000 ) {
        q_hash_size = ( n_node > n_elem ) ?   n_node :   n_elem;
        q_buf_size  = ( n_node > n_elem ) ? 5*n_node : 5*n_elem;
    } else {
        q_hash_size = ( n_node > n_elem ) ?   n_node :   n_elem;
        q_buf_size  = ( n_node > n_elem ) ?   n_node :   n_elem;
    }

    qsuf_node = (struct hecmw_qsuf_node *)HECMW_malloc( sizeof(struct hecmw_qsuf_node) );
    if( qsuf_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        qsuf_node->node1 = NULL;
        qsuf_node->node2 = NULL;
        qsuf_node->node3 = NULL;
        qsuf_node->node4 = NULL;
    }
    qsuf_node->node1 = (int *)HECMW_malloc( sizeof(int)*q_buf_size );
    if( qsuf_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node2 = (int *)HECMW_malloc( sizeof(int)*q_buf_size );
    if( qsuf_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node3 = (int *)HECMW_malloc( sizeof(int)*q_buf_size );
    if( qsuf_node->node3 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node4 = (int *)HECMW_malloc( sizeof(int)*q_buf_size );
    if( qsuf_node->node4 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    size = sizeof(struct hecmw_mesh_hash_link *) * q_hash_size;
    q_hash_tbl = (struct hecmw_mesh_hash_link **)HECMW_malloc( size );
    if( q_hash_tbl == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    } else {
        for( i=0; i<q_hash_size; i++ ) {
            q_hash_tbl[i] = NULL;
        }
    }

    n_qsuf = 0;

    return 0;

error:
    HECMW_mesh_hsort_qsuf_final( );

    return -1;
}

/*================================================================================================*/

extern int
HECMW_mesh_hsort_edge_realloc( void )
{
    unsigned long int new_buf_size;

    new_buf_size = (unsigned long int)( e_buf_size * EDGE_INC_FACTOR );

    edge_node->node1 = (int *)HECMW_realloc( edge_node->node1, sizeof(int)*new_buf_size );
    if( edge_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    edge_node->node2 = (int *)HECMW_realloc( edge_node->node2, sizeof(int)*new_buf_size );
    if( edge_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    e_buf_size = new_buf_size;

    return 0;

error:
    HECMW_mesh_hsort_edge_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_tsuf_realloc( void )
{
    unsigned long int new_buf_size;

    new_buf_size = (unsigned long int)( t_buf_size * TSUF_INC_FACTOR );

    tsuf_node->node1 = (int *)HECMW_realloc( tsuf_node->node1, sizeof(int)*new_buf_size );
    if( tsuf_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    tsuf_node->node2 = (int *)HECMW_realloc( tsuf_node->node2, sizeof(int)*new_buf_size );
    if( tsuf_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    tsuf_node->node3 = (int *)HECMW_realloc( tsuf_node->node3, sizeof(int)*new_buf_size );
    if( tsuf_node->node3 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    t_buf_size = new_buf_size;

    return 0;

error:
    HECMW_mesh_hsort_tsuf_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_qsuf_realloc( void )
{
    unsigned long int new_buf_size;

    new_buf_size = (unsigned long int)( q_buf_size * QSUF_INC_FACTOR );

    qsuf_node->node1 = (int *)HECMW_realloc( qsuf_node->node1, sizeof(int)*new_buf_size );
    if( qsuf_node->node1 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node2 = (int *)HECMW_realloc( qsuf_node->node2, sizeof(int)*new_buf_size );
    if( qsuf_node->node2 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node3 = (int *)HECMW_realloc( qsuf_node->node3, sizeof(int)*new_buf_size );
    if( qsuf_node->node3 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    qsuf_node->node4 = (int *)HECMW_realloc( qsuf_node->node4, sizeof(int)*new_buf_size );
    if( qsuf_node->node4 == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    q_buf_size = new_buf_size;

    return 0;

error:
    HECMW_mesh_hsort_qsuf_final( );

    return -1;
}

/*================================================================================================*/

extern int
HECMW_mesh_hsort_edge_get_n( void )
{
    return n_edge;
}


extern int
HECMW_mesh_hsort_tsuf_get_n( void )
{
    return n_tsuf;
}


extern int
HECMW_mesh_hsort_qsuf_get_n( void )
{
    return n_qsuf;
}

/*------------------------------------------------------------------------------------------------*/

extern int
*HECMW_mesh_hsort_edge_get_v( void )
{
    int i;

    __edge_node = (int *)HECMW_malloc( sizeof(int)*n_edge*2 );
    if( __edge_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    for( i=0; i<n_edge; i++ ) {
        __edge_node[2*i  ] = edge_node->node1[i];
        __edge_node[2*i+1] = edge_node->node2[i];
    }

    return __edge_node;

error:
    HECMW_free( __edge_node );
    __edge_node = NULL;
    HECMW_mesh_hsort_edge_final( );

    return NULL;
}


extern int
*HECMW_mesh_hsort_tsuf_get_v( void )
{
    int i;

    __tsuf_node = (int *)HECMW_malloc( sizeof(int)*n_tsuf*3 );
    if( __tsuf_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    for( i=0; i<n_tsuf; i++ ) {
        __tsuf_node[3*i  ] = tsuf_node->node1[i];
        __tsuf_node[3*i+1] = tsuf_node->node2[i];
        __tsuf_node[3*i+2] = tsuf_node->node3[i];
    }

    return __tsuf_node;

error:
    HECMW_free( __tsuf_node );
    __tsuf_node = NULL;
    HECMW_mesh_hsort_tsuf_final( );

    return NULL;
}


extern int
*HECMW_mesh_hsort_qsuf_get_v( void )
{
    int i;

    __qsuf_node = (int *)HECMW_malloc( sizeof(int)*n_qsuf*4 );
    if( __qsuf_node == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }

    for( i=0; i<n_qsuf; i++ ) {
        __qsuf_node[4*i  ] = qsuf_node->node1[i];
        __qsuf_node[4*i+1] = qsuf_node->node2[i];
        __qsuf_node[4*i+2] = qsuf_node->node3[i];
        __qsuf_node[4*i+3] = qsuf_node->node4[i];
    }

    return __qsuf_node;

error:
    HECMW_free( __qsuf_node );
    HECMW_mesh_hsort_qsuf_final( );

    return NULL;
}

/*================================================================================================*/

extern void
HECMW_mesh_hsort_edge_final( void )
{
    if( e_hash_tbl ) {
        int i;
        struct hecmw_mesh_hash_link *p, *q;
        for( i=0; i<e_hash_size; i++ ) {
            if( e_hash_tbl[i] ) {
                for( q=e_hash_tbl[i], p=e_hash_tbl[i]; p; p=q ) {
                    q = q->next;
                    HECMW_free( p );
                }
                e_hash_tbl[i] = NULL;
            }
        }
        HECMW_free( e_hash_tbl );
    }
    if( edge_node ) {
        HECMW_free( edge_node->node1 );
        HECMW_free( edge_node->node2 );
    }
    HECMW_free( edge_node );

    e_hash_tbl = NULL;
    edge_node  = NULL;
}


extern void
HECMW_mesh_hsort_tsuf_final( void )
{
    if( t_hash_tbl ) {
        int i;
        struct hecmw_mesh_hash_link *p, *q;
        for( i=0; i<t_hash_size; i++ ) {
            if( t_hash_tbl[i] ) {
                for( q=t_hash_tbl[i], p=t_hash_tbl[i]; p; p=q ) {
                    q = q->next;
                    HECMW_free( p );
                }
                t_hash_tbl[i] = NULL;
            }
        }
        HECMW_free( t_hash_tbl );
    }
    if( tsuf_node ) {
        HECMW_free( tsuf_node->node1 );
        HECMW_free( tsuf_node->node2 );
        HECMW_free( tsuf_node->node3 );
    }
    HECMW_free( tsuf_node );

    t_hash_tbl = NULL;
    tsuf_node  = NULL;
}


extern void
HECMW_mesh_hsort_qsuf_final( void )
{
    struct hecmw_mesh_hash_link *p, *q;
    int i;

    if( q_hash_tbl ) {
        for( i=0; i<q_hash_size; i++ ) {
            if( q_hash_tbl[i] ) {
                for( q=q_hash_tbl[i], p=q_hash_tbl[i]; p; p=q ) {
                    q = q->next;
                    HECMW_free( p );
                }
                q_hash_tbl[i] = NULL;
            }
        }
        HECMW_free( q_hash_tbl );
    }
    if( qsuf_node ) {
        HECMW_free( qsuf_node->node1 );
        HECMW_free( qsuf_node->node2 );
        HECMW_free( qsuf_node->node3 );
        HECMW_free( qsuf_node->node4 );
    }
    HECMW_free( qsuf_node );

    q_hash_tbl = NULL;
    qsuf_node  = NULL;
}

/*================================================================================================*/

static void
reorder_node_edge( int m1, int m2, int *n1, int *n2 )
{
    if( m1 < m2 ) {
        *n1 = m1;
        *n2 = m2;
    } else {
        *n1 = m2;
        *n2 = m1;
    }
}


static void
reorder_node_tsuf( int m1, int m2, int m3, int *n1, int *n2, int *n3 )
{
    int l1, l2, l3;

    if( m1 < m2 ) {
        l1 = m1;
        l2 = m2;
    } else {
        l1 = m2;
        l2 = m1;
    }

    if( m3 < l1 ) {
        *n1 = m3;
        l3 = l1;
    } else {
        *n1 = l1;
        l3 = m3;
    }

    if( l2 < l3 ) {
        *n2 = l2;
        *n3 = l3;
    } else {
        *n2 = l3;
        *n3 = l2;
    }
}


static void
reorder_node_qsuf( int m1, int m2, int m3, int m4,
                   int *n1, int *n2, int *n3, int *n4 )
{
    int l1, l2, l3, l4, l5, l6;

    if( m1 < m2 ) {
        l1 = m1;
        l2 = m2;
    } else {
        l1 = m2;
        l2 = m1;
    }

    if( m3 < m4 ) {
        l3 = m3;
        l4 = m4;
    } else {
        l3 = m4;
        l4 = m3;
    }

    if( l1 < l3 ) {
        *n1 = l1;
        l5  = l3;
    } else {
        *n1 = l3;
        l5  = l1;
    }

    if( l2 > l4 ) {
        *n4 = l2;
        l6  = l4;
    } else {
        *n4 = l4;
        l6  = l2;
    }

    if( l5 < l6 ) {
        *n2 = l5;
        *n3 = l6;
    } else {
        *n2 = l6;
        *n3 = l5;
    }
}

/*================================================================================================*/

extern int
HECMW_mesh_hsort_edge( int node1, int node2 )
{
    int n1, n2, m1, m2;
    int eid;
    int idx;
    unsigned long int ndot;
    struct hecmw_mesh_hash_link *p;


    reorder_node_edge( node1, node2, &n1, &n2 );


    ndot = ((unsigned long int)n1 % e_hash_size) * ((unsigned long int)n2 % e_hash_size);
    idx = ndot % e_hash_size;

    for( p=e_hash_tbl[idx]; p; p=p->next ) {
        eid = p->key;
        reorder_node_edge( edge_node->node1[eid], edge_node->node2[eid], &m1, &m2 );
        if(( n1 == m1 ) && ( n2 == m2 )) {
            return eid+1;
        }
    }


    p = (struct hecmw_mesh_hash_link *)HECMW_malloc( sizeof(struct hecmw_mesh_hash_link) );
    if( p == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    p->next         = e_hash_tbl[idx];
    e_hash_tbl[idx] = p;

    if( n_edge >= e_buf_size ) {
        if( HECMW_mesh_hsort_edge_realloc( ) ) {
            goto error;
        }
    }

    eid                   = n_edge;
    p->key                = eid;
    edge_node->node1[eid] = node1;
    edge_node->node2[eid] = node2;
    n_edge++;

    return eid+1;

error:
    HECMW_mesh_hsort_edge_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_tsuf( int node1, int node2, int node3 )
{
    int n1, n2, n3, m1, m2, m3;
    int tid;
    int idx;
    unsigned long int ndot1, ndot;
    struct hecmw_mesh_hash_link *p;


    reorder_node_tsuf( node1, node2, node3, &n1, &n2, &n3 );


    ndot1 = ((unsigned long int)n1 % t_hash_size) * ((unsigned long int)n2 % t_hash_size);
    ndot  = ((unsigned long int)n3 % t_hash_size) * (ndot1 % t_hash_size);
    idx   = ndot % t_hash_size;

    for( p=t_hash_tbl[idx]; p; p=p->next ) {
        tid = p->key;
        reorder_node_tsuf( tsuf_node->node1[tid], tsuf_node->node2[tid],
                                              tsuf_node->node3[tid], &m1, &m2, &m3 );
        if(( n1 == m1 ) && ( n2 == m2 ) && ( n3 == m3 )) {
            return tid+1;
        }
    }


    p = (struct hecmw_mesh_hash_link *)HECMW_malloc( sizeof(struct hecmw_mesh_hash_link) );
    if( p == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    p->next         = t_hash_tbl[idx];
    t_hash_tbl[idx] = p;

    if( n_tsuf >= t_buf_size ) {
        if( HECMW_mesh_hsort_tsuf_realloc( ) ) {
            goto error;
        }
    }

    tid                   = n_tsuf;
    p->key                = tid;
    p->next               = NULL;
    tsuf_node->node1[tid] = node1;
    tsuf_node->node2[tid] = node2;
    tsuf_node->node3[tid] = node3;
    n_tsuf++;

    return tid+1;

error:
    HECMW_mesh_hsort_tsuf_final( );

    return -1;
}


extern int
HECMW_mesh_hsort_qsuf( int node1, int node2, int node3, int node4 )
{
    int n1, n2, n3, n4, m1, m2, m3, m4;
    int qid;
    int idx;
    unsigned long int ndot1, ndot2, ndot;
    struct hecmw_mesh_hash_link *p;


    reorder_node_qsuf( node1, node2, node3, node4, &n1, &n2, &n3, &n4 );


    ndot1 = (n1 % q_hash_size) * (n2 % q_hash_size);
    ndot2 = (n3 % q_hash_size) * (n4 % q_hash_size);
    ndot  = (ndot1 % q_hash_size) * (ndot2 % q_hash_size);
    idx   = ndot % q_hash_size;

    for( p=q_hash_tbl[idx]; p; p=p->next ) {
        qid = p->key;
        reorder_node_qsuf( qsuf_node->node1[qid], qsuf_node->node2[qid],
                                              qsuf_node->node3[qid], qsuf_node->node4[qid],
                                              &m1, &m2, &m3, &m4 );

        if(( n1 == m1 ) && ( n2 == m2 ) && ( n3 == m3 ) && ( n4 == m4 )) {
            return qid+1;
        }
    }


    p = (struct hecmw_mesh_hash_link *)HECMW_malloc( sizeof(struct hecmw_mesh_hash_link) );
    if( p == NULL ) {
        HECMW_set_error( errno, "" );
        goto error;
    }
    p->next         = q_hash_tbl[idx];
    q_hash_tbl[idx] = p;

    if( n_qsuf >= q_buf_size ) {
        if( HECMW_mesh_hsort_qsuf_realloc( ) ) {
            goto error;
        }
    }

    qid                   = n_qsuf;
    p->key                = qid;
    qsuf_node->node1[qid] = node1;
    qsuf_node->node2[qid] = node2;
    qsuf_node->node3[qid] = node3;
    qsuf_node->node4[qid] = node3;
    n_qsuf++;

    return qid+1;

error:
    HECMW_mesh_hsort_qsuf_final( );

    return -1;
}
