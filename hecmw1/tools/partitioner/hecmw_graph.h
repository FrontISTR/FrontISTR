/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2014/01/25                                        *
 *        Category : HEC-MW Utility                                    *
 *                                                                     *
 *            Written by Kazuya Goto (PExProCS LLC)                    *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/
/**
 * @file  hecmw_graph.h
 * @brief Graph utility
 *
 * @author Kazuya Goto (VINAS)
 * @date Feb 6, 2012
 */

#ifndef HECMW_GRAPH_INCLUDED
#define HECMW_GRAPH_INCLUDED

#include <stdio.h>

struct hecmw_varray_int;

/** Graph data structure.
 *
 */
struct hecmw_graph {
    int m_num_vertex; /**< number of vertices */
    int m_num_edge;   /**< number of edges (double in both ways) */
    struct hecmw_varray_int *m_edge_index;    /**< edge index array (length is m_num_vertex+1) */
    struct hecmw_varray_int *m_edge_item;     /**< edge item array (length is m_edge_index[m_num_vertex]) */
    struct hecmw_varray_int *m_vertex_weight; /**< vertex weight array (length is m_num_vertex) */
    int is_ref;
};

/** Initialize.
 *
 * @retval HECMW_SUCCESS normal return
 * @retval HECMW_ERROR   failed to initialize
 */
extern int
HECMW_graph_init(
    struct hecmw_graph *graph /**< [inout] graph */
    );

/** Initialize with given index and item arrays.
 *
 * @retval HECMW_SUCCESS normal return
 * @retval HECMW_ERROR   failed to initialize
 */
extern int
HECMW_graph_init_with_arrays(
    struct hecmw_graph *graph, /**< [inout] graph */
    int num_vertex,            /**< [in] number of vertices */
    int *edge_index,           /**< [in] edge index array */
    int *edge_item             /**< [in] edge item array */
    );

/** Finalize.
 */
extern void
HECMW_graph_finalize(
    struct hecmw_graph *graph /**< [inout] graph */
    );

/** Set number of vertices.
 */
extern void HECMW_graph_setNumVertex(
    struct hecmw_graph *graph, /**< [inout] graph */
    int num_vertex             /**< [in] number of vertices */
    );

/** Add edge.
 *
 * @retval HECMW_SUCCESS normal return
 * @retval HECMW_ERROR   failed to add edge
 */
extern int HECMW_graph_addEdge(
    struct hecmw_graph *graph, /**< [inout] graph */
    int vert1,                 /**< [in] the first vertex of the edge */
    int vert2                  /**< [in] the other vertex of the edge */
    );

/** Print graph.
 */
extern void HECMW_graph_print(
    const struct hecmw_graph *graph, /**< [in] graph */
    FILE *fp                         /**< [in] File pointer for output */
    );

/** Get number of vertices.
 *
 * @return number of vertices
 */
extern int HECMW_graph_getNumVertex(
    const struct hecmw_graph *graph /**< [in] graph */
    );

/** Get number of edges (doubled due to both-ways).
 *
 * @return number of edges
 */
extern int HECMW_graph_getNumEdge(
    const struct hecmw_graph *graph /**< [in] graph */
    );

/** Get edge index array.
 *
 * @return head pointer of the edge index array
 */
extern const int *HECMW_graph_getEdgeIndex(
    const struct hecmw_graph *graph /**< [in] graph */
    );

/** Get edge item array.
 *
 * @return head pointer of the edge item array
 */
extern const int *HECMW_graph_getEdgeItem(
    const struct hecmw_graph *graph /**< [in] graph */
    );

/** Degenerate graph.
 *
 * @retval HECMW_SUCCESS normal return
 * @retval HECMW_ERROR   failed to create degenerated graph
 */
extern int HECMW_graph_degeneGraph(
    struct hecmw_graph *graph,          /**< [inout] graph */
    const struct hecmw_graph *refgraph, /**< [in] original graph to be degenerated */
    int num_part,                       /**< [in] number of vertices after degeneration */
    const int *parttab                  /**< [in] id of each vertex after degeneration */
    );

#endif /* HECMW_GRAPH_INCLUDED */
