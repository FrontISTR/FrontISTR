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
 * @file  hecmw_graph.c
 * @brief Graph Utility (implementation)
 *
 * @author Kazuya Goto (VINAS)
 * @date Feb 6, 2012
 */

#include "hecmw_graph.h"
#include "hecmw_varray_int.h"
#include "hecmw_malloc.h"
#include "hecmw_config.h"
#include "hecmw_util.h"
#include <stdio.h>
#include <errno.h>

/*********************************
 * Prototype of static functions *
 *********************************/

/** Clear graph.
 */
static void clear(
    struct hecmw_graph *graph /**< [inout] graph */
    );

/** Find edge.
 *
 * @retval HECMW_SUCCESS edge found
 * @retval HECMW_ERROR   edge not found
 */
static int find_edge(
    const struct hecmw_graph *graph, /**< [in] graph */
    int vert1,                       /**< [in] starting vertex id */
    int vert2,                       /**< [in] end vertex id */
    int *idx                         /**< [out] edge id (set only when found) */
    );

/** Add edge (one-way only); do nothing when already exists.
 *
 * @retval HECMW_SUCCESS normal return
 * @retval HECMW_ERROR   failed to add edge
 */
static int add_edge_one_way(
    struct hecmw_graph *graph, /**< [inout] graph */
    int vert1,                 /**< [in] starting vertex id */
    int vert2                  /**< [in] end vertex id */
    );

/**********************************
 * Definition of public functions *
 **********************************/

int
HECMW_graph_init(struct hecmw_graph *graph)
{
    graph->m_num_vertex = 0;
    graph->m_num_edge = 0;
    graph->m_edge_index = (struct hecmw_varray_int *) HECMW_malloc(sizeof(struct hecmw_varray_int));
    graph->m_edge_item = (struct hecmw_varray_int *) HECMW_malloc(sizeof(struct hecmw_varray_int));
    if (graph->m_edge_index == NULL ||
        graph->m_edge_item == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
    }
    if (HECMW_varray_int_init(graph->m_edge_index) == HECMW_SUCCESS &&
        HECMW_varray_int_init(graph->m_edge_item) == HECMW_SUCCESS)
        return HECMW_SUCCESS;
    graph->is_ref = 0;
    return HECMW_ERROR;
}

int
HECMW_graph_init_with_arrays(struct hecmw_graph *graph,
                             int num_vertex, int *edge_index, int *edge_item)
{
    graph->m_num_vertex = num_vertex;
    graph->m_num_edge = edge_index[num_vertex];
    graph->m_edge_index = (struct hecmw_varray_int *) HECMW_malloc(sizeof(struct hecmw_varray_int));
    graph->m_edge_item = (struct hecmw_varray_int *) HECMW_malloc(sizeof(struct hecmw_varray_int));
    if (graph->m_edge_index == NULL ||
        graph->m_edge_item == NULL) {
        HECMW_set_error(errno, "");
        return HECMW_ERROR;
    }
    graph->m_edge_index->n_val = num_vertex+1;
    graph->m_edge_index->max_val = num_vertex+1;
    graph->m_edge_index->vals = edge_index;

    graph->m_edge_item->n_val = graph->m_num_edge;
    graph->m_edge_item->max_val = graph->m_num_edge;
    graph->m_edge_item->vals = edge_item;

    graph->is_ref = 1;
    return HECMW_SUCCESS;
}

void
HECMW_graph_finalize(struct hecmw_graph *graph)
{
    if (!graph->is_ref) {
        HECMW_varray_int_finalize(graph->m_edge_index);
        HECMW_varray_int_finalize(graph->m_edge_item);
    }
    HECMW_free(graph->m_edge_index);
    HECMW_free(graph->m_edge_item);
}

void HECMW_graph_setNumVertex(struct hecmw_graph *graph, int num_vertex)
{
    HECMW_assert(!graph->is_ref);

    graph->m_num_vertex = num_vertex;
    HECMW_varray_int_resize(graph->m_edge_index, num_vertex + 1);
    HECMW_varray_int_assign(graph->m_edge_index, 0, num_vertex + 1, 0);
}

int HECMW_graph_addEdge(struct hecmw_graph *graph, int vert1, int vert2)
{
    HECMW_assert(!graph->is_ref);

    if (add_edge_one_way(graph, vert1, vert2) == HECMW_SUCCESS &&
        add_edge_one_way(graph, vert2, vert1) == HECMW_SUCCESS)
        return HECMW_SUCCESS;
    return HECMW_ERROR;
}

void HECMW_graph_print(const struct hecmw_graph *graph, FILE *fp)
{
    const int *edge_index = HECMW_varray_int_get_cv(graph->m_edge_index);
    const int *edge_item = HECMW_varray_int_get_cv(graph->m_edge_item);
    int i, j;
    int idx_start, idx_end;

    fprintf(fp, "num_vertex = %d\n", graph->m_num_vertex);
    fprintf(fp, "num_edge = %d\n", graph->m_num_edge);

    for (i = 0; i < graph->m_num_vertex; i++) {
        fprintf(fp, "%d: ", i);

        idx_start = edge_index[i];
        idx_end = edge_index[i+1];
        for (j = idx_start; j < idx_end; j++) {
            fprintf(fp, " %d", edge_item[j]);
        }
        fprintf(fp, "\n");
    }
}

int HECMW_graph_getNumVertex(const struct hecmw_graph *graph)
{
    return graph->m_num_vertex;
}

int HECMW_graph_getNumEdge(const struct hecmw_graph *graph)
{
    return graph->m_num_edge;
}

const int *HECMW_graph_getEdgeIndex(const struct hecmw_graph *graph)
{
    return HECMW_varray_int_get_cv(graph->m_edge_index);
}

const int *HECMW_graph_getEdgeItem(const struct hecmw_graph *graph)
{
    return HECMW_varray_int_get_cv(graph->m_edge_item);
}

int HECMW_graph_degeneGraph(struct hecmw_graph *graph,
                            const struct hecmw_graph *refgraph,
                            int num_part, const int *parttab)
{
    const int *ref_edge_index = HECMW_varray_int_get_cv(refgraph->m_edge_index);
    const int *ref_edge_item = HECMW_varray_int_get_cv(refgraph->m_edge_item);
    int i, j, jj;
    int i_part, j_part;
    int start, end;
    int retval;

    clear(graph);
    HECMW_graph_setNumVertex(graph, num_part);
    for (i = 0; i < HECMW_graph_getNumVertex(refgraph); i++) {
        i_part = parttab[i];
        start = ref_edge_index[i];
        end = ref_edge_index[i + 1];
        for (j = start; j < end; j++) {
            jj = ref_edge_item[j];
            j_part = parttab[jj];
            if (i_part == j_part) continue;
            retval = HECMW_graph_addEdge(graph, i_part, j_part);
            if (retval != HECMW_SUCCESS) return HECMW_ERROR;
        }
    }
    return HECMW_SUCCESS;
}

/***********************************
 * Definition of private functions *
 ***********************************/

void clear(struct hecmw_graph *graph)
{
    graph->m_num_vertex = 0;
    graph->m_num_edge = 0;
    HECMW_varray_int_resize(graph->m_edge_index, 0);
    HECMW_varray_int_resize(graph->m_edge_item, 0);
    graph->is_ref = 0;
}

int find_edge(const struct hecmw_graph *graph,
              int vert1, int vert2, int *idx)
{
    const int *edge_index = HECMW_varray_int_get_cv(graph->m_edge_index);
    const int *edge_item = HECMW_varray_int_get_cv(graph->m_edge_item);
    int idx_start, idx_end;
    int i;

    idx_start = edge_index[vert1];
    idx_end = edge_index[vert1 + 1];
    for (i = idx_start; i < idx_end; i++) {
        if (edge_item[i] == vert2) {
            if (idx) *idx = i;
            return 1;
        }
    }
    return 0;
}

int add_edge_one_way(struct hecmw_graph *graph,
                     int vert1, int vert2)
{
    int *edge_index = HECMW_varray_int_get_v(graph->m_edge_index);
    int idx;
    int i;
    int retval;

    HECMW_assert(!graph->is_ref);

    if (find_edge(graph, vert1, vert2, &idx)) {
        return HECMW_SUCCESS;
    }
    /* insert vert2 into m_edge_item */
    /* place to insert: m_edge_inidex[vert1 + 1] */
    retval = HECMW_varray_int_insert(graph->m_edge_item, edge_index[vert1 + 1], vert2);
    if (retval != HECMW_SUCCESS) {
        return HECMW_ERROR;
    }

    /* increment m_edge_index[vert1 + 1 .. n_num_vertex] */
    for (i = vert1 + 1; i <= graph->m_num_vertex; i++) {
        edge_index[i] += 1;
    }
    graph->m_num_edge++;
    return HECMW_SUCCESS;
}
