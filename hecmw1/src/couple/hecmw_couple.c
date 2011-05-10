/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 1.00                                              *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : Coupling Interface                                *
 *                                                                     *
 *            Written by Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using Hight End Computing Middleware (HEC-MW)"     *
 *                                                                     *
 *=====================================================================*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "hecmw_msgno.h"
#include "hecmw_common_define.h"
#include "hecmw_struct.h"
#include "hecmw_malloc.h"
#include "hecmw_log.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_comm.h"
#include "hecmw_couple_init.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_intra_iftable.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_weight.h"
#include "hecmw_couple_startup.h"
#include "hecmw_couple.h"



struct hecmw_couple_values {
	int item_type;						
	struct hecmw_couple_value *node;	
	struct hecmw_couple_value *elem;	
	struct hecmw_couple_value *surf;	
};


/*================================================================================================*/

static void
free_couple_values(struct hecmw_couple_values *p)
{
	if(p == NULL) return;

	HECMW_couple_free_couple_value(p->node);
	HECMW_couple_free_couple_value(p->elem);
	HECMW_couple_free_couple_value(p->surf);
	HECMW_free(p);
	p = NULL;
}



static struct hecmw_couple_values *
alloc_couple_values(void)
{
	struct hecmw_couple_values *p = NULL;

	p = (struct hecmw_couple_values *)HECMW_malloc(sizeof(struct hecmw_couple_values));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->item_type = HECMW_COUPLE_GROUP_UNDEF;
	p->node = NULL;
	p->elem = NULL;
	p->surf = NULL;

	if((p->node = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	if((p->elem = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	if((p->surf = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	p->node->item_type = HECMW_COUPLE_NODE_GROUP;
	p->elem->item_type = HECMW_COUPLE_ELEMENT_GROUP;
	p->surf->item_type = HECMW_COUPLE_SURFACE_GROUP;

	return p;

error:
	free_couple_values(p);
	return NULL;
}


/*================================================================================================*/

static int
update_import_node_value(struct hecmw_couple_comm *intracomm,
		struct hecmw_couple_intra_iftable *intra_tbl, struct hecmw_couple_value *node_value)
{
	int *sendbuf_index = NULL, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int nmemb, rtc, id, i, j;

	if(intra_tbl->n_neighbor_pe == 0) return 0;

	/*
	 * send buffer
	 */
	/* index for send buffer */
	sendbuf_index = (int *)HECMW_calloc(intra_tbl->n_neighbor_pe+1, sizeof(int));
	if(sendbuf_index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<intra_tbl->n_neighbor_pe; i++) {
		sendbuf_index[i+1] = intra_tbl->export_index[i+1] * node_value->n_dof;
	}

	/* send buffer */
	nmemb = intra_tbl->export_index[intra_tbl->n_neighbor_pe] * node_value->n_dof + 1;
	sendbuf = (double *)HECMW_malloc(sizeof(double)*nmemb);
	if(sendbuf == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<intra_tbl->export_index[intra_tbl->n_neighbor_pe]; i++) {
		id = intra_tbl->export_item[i];
		for(j=0; j<node_value->n_dof; j++) {
			sendbuf[node_value->n_dof*i+j] = node_value->value[node_value->n_dof*id+j];
		}
	}

	/*
	 * receive buffer
	 */
	/* index for receive buffer */
	recvbuf_index = (int *)HECMW_calloc(intra_tbl->n_neighbor_pe+1, sizeof(int));
	if(recvbuf_index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<intra_tbl->n_neighbor_pe; i++) {
		recvbuf_index[i+1] = intra_tbl->import_index[i+1] * node_value->n_dof;
	}

	/* receive buffer */
	nmemb = intra_tbl->import_index[intra_tbl->n_neighbor_pe] * node_value->n_dof + 1;
	recvbuf = (double *)HECMW_malloc(sizeof(double)*nmemb);
	if(recvbuf == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_intra_send_recv(intra_tbl->n_neighbor_pe, intra_tbl->neighbor_pe,
			sendbuf_index, sendbuf, recvbuf_index, recvbuf, HECMW_DOUBLE, intracomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received value
	 */
	for(i=0; i<intra_tbl->import_index[intra_tbl->n_neighbor_pe]; i++) {
		id = intra_tbl->import_item[i];
		for(j=0; j<node_value->n_dof; j++) {
			node_value->value[node_value->n_dof*id+j] = recvbuf[node_value->n_dof*i+j];
		}
	}

	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return 0;

error:
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return -1;
}



static int
interpolation(const struct hecmw_couple_weight *p,
		struct hecmw_couple_value *value_src, struct hecmw_couple_value *value_dst)
{
	int n_dof, i, j, k;

	n_dof = value_dst->n_dof;
	for(i=0; i<p->n; i++) {
		for(j=p->index[i]; j<p->index[i+1]; j++) {
			for(k=0; k<n_dof; k++) {
				value_dst->value[i*n_dof+k] += value_src->value[(p->id[j])*n_dof+k] * p->weight[j];
			}
		}
	}

	return 0;
}



static int
pre_interpolation(const struct hecmw_couple_weight_list *ip_list_pre,
		struct hecmw_couple_values *values)
{
	struct hecmw_couple_weight_list *p;

	for(p=ip_list_pre->next; p; p=p->next) {
		if(p->info->type == HECMW_COUPLE_IP_SURF_TO_NODE) {		/* surface -> node	*/
			if(interpolation(p->info, values->surf, values->node)) return -1;

		} else {												/* error			*/
			HECMW_set_error(HECMWCPL_E_INVALID_IPTYPE, "");
			return -1;
		}
	}

	return 0;
}



static int
main_interpolation(const struct hecmw_couple_weight_list *ip_list_main,
		const struct hecmw_couple_values *values_src,
		struct hecmw_couple_value *value_send)
{
	struct hecmw_couple_weight_list *p;

	for(p=ip_list_main->next; p; p=p->next) {
		if(p->info->type == HECMW_COUPLE_IP_NODE_TO_NODE) {			/* node -> node		*/
			if(interpolation(p->info, values_src->node, value_send)) return -1;

		} else if(p->info->type == HECMW_COUPLE_IP_SURF_TO_NODE) {	/* surface -> node	*/
			if(interpolation(p->info, values_src->surf, value_send)) return -1;

		} else {													/* error			*/
			HECMW_set_error(HECMWCPL_E_INVALID_IPTYPE, "");
			return -1;
		}
	}

	return 0;
}



static int
post_interpolation(const struct hecmw_couple_weight_list *ip_list_post,
		struct hecmw_couple_values *values)
{
	struct hecmw_couple_weight_list *p;

	for(p=ip_list_post->next; p; p=p->next) {
		if(p->info->type == HECMW_COUPLE_IP_NODE_TO_SURF) {		/* node -> surface	*/
			if(interpolation(p->info, values->node, values->surf)) return -1;

		} else {												/* error			*/
			HECMW_set_error(HECMWCPL_E_INVALID_MAPTYPE, "");
			return -1;
		}
	}

	return 0;
}


/*================================================================================================*/

static int
send_recv_n_dof(const struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_value *value_src,
		struct hecmw_couple_value *value_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int n_dof, size, rtc, i;

	/*
	 * send buffer
	 */
	if(comm_src->is_member) {
		if(inter_tbl->n_neighbor_pe_export > 0) {
			/* index for send buffer */
			sendbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_export+1, sizeof(int));
			if(sendbuf_index == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->n_neighbor_pe_export; i++) {
				sendbuf_index[i+1] = sendbuf_index[i] + 1;
			}

			/* send buffer */
			size = sizeof(int) * sendbuf_index[inter_tbl->n_neighbor_pe_export];
			sendbuf = (int *)HECMW_malloc(size);
			if(sendbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<sendbuf_index[inter_tbl->n_neighbor_pe_export]; i++) {
				sendbuf[i] = value_src->n_dof;
			}
		}
	}

	/*
	 * receive buffer
	 */
	if(comm_dst->is_member) {
		if(inter_tbl->n_neighbor_pe_import > 0) {
			/* index for receive buffer */
			recvbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_import+1, sizeof(int));
			if(recvbuf_index == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->n_neighbor_pe_import; i++) {
				recvbuf_index[i+1] = recvbuf_index[i] + 1;
			}

			/* receive buffer */
			size = sizeof(int) * recvbuf_index[inter_tbl->n_neighbor_pe_import];
			recvbuf = (int *)HECMW_malloc(size);
			if(recvbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(
			inter_tbl->n_neighbor_pe_export, inter_tbl->neighbor_pe_export, sendbuf_index, sendbuf,
			inter_tbl->n_neighbor_pe_import, inter_tbl->neighbor_pe_import, recvbuf_index, recvbuf,
			HECMW_INT, intercomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received value
	 */
	if(comm_dst->is_member) {
		n_dof = recvbuf[0];
		for(i=1; i<inter_tbl->n_neighbor_pe_import; i++) {
			HECMW_assert(n_dof == recvbuf[i]);
		}
		value_dst->n_dof = n_dof;
	}

	HECMW_free(send_pe);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recv_pe);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return 0;

error:
	HECMW_free(send_pe);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recv_pe);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return -1;
}



static int
send_recv_couple_value(const struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_value *value_src,
		struct hecmw_couple_value *value_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int boundary_index, id, n_dof, size, rtc, i, j;

	/*
	 * send buffer
	 */
	if(comm_src->is_member) {
		n_dof = value_src->n_dof;

		if(inter_tbl->n_neighbor_pe_export > 0) {
			/* index for send buffer */
			sendbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_export+1, sizeof(int));
			if(sendbuf_index == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->n_neighbor_pe_export; i++) {
				sendbuf_index[i+1] = inter_tbl->export_index[i+1] * n_dof;
			}

			/* send buffer */
			size = sizeof(double) * (sendbuf_index[inter_tbl->n_neighbor_pe_export]+1);
			sendbuf = (double *)HECMW_malloc(size);
			if(sendbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]*n_dof; i++) {
				sendbuf[i] = value_src->value[i];
			}
		}
	}

	/*
	 * receive buffer
	 */
	if(comm_dst->is_member) {
		if(inter_tbl->n_neighbor_pe_import > 0) {
			/* index for receive buffer */
			recvbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_import+1, sizeof(int));
			if(recvbuf_index == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->n_neighbor_pe_import; i++) {
				recvbuf_index[i+1] = inter_tbl->import_index[i+1] * value_dst->n_dof;
			}

			/* receive buffer */
			size = sizeof(double) * (recvbuf_index[inter_tbl->n_neighbor_pe_import]+1);
			recvbuf = (double *)HECMW_malloc(size);
			if(recvbuf == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
		}
	}

	if(!comm_src->is_member && !comm_dst->is_member)  return HECMW_SUCCESS;

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(
			inter_tbl->n_neighbor_pe_export, inter_tbl->neighbor_pe_export, sendbuf_index, sendbuf,
			inter_tbl->n_neighbor_pe_import, inter_tbl->neighbor_pe_import, recvbuf_index, recvbuf,
			HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received value
	 */
	if(comm_dst->is_member) {
		n_dof = value_dst->n_dof;

		size = sizeof(double) * inter_tbl->import_index[inter_tbl->n_neighbor_pe_import] * n_dof;
		value_dst->value = (double *)HECMW_malloc(size);
		if(value_dst->value == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}

		for(i=0; i<inter_tbl->import_index[inter_tbl->n_neighbor_pe_import]; i++) {
			id = inter_tbl->import_item[i];
			for(j=0; j<n_dof; j++) {
				value_dst->value[n_dof*id+j] = recvbuf[n_dof*i+j];
			}
		}
	}

	HECMW_free(send_pe);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recv_pe);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return 0;

error:
	HECMW_free(send_pe);
	HECMW_free(sendbuf_index);
	HECMW_free(sendbuf);
	HECMW_free(recv_pe);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	return -1;
}



static int
send_recv(const struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_value *value_src,
		struct hecmw_couple_value *value_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm)
{
	if(send_recv_n_dof(inter_tbl, value_src, value_dst, comm_src, comm_dst, intercomm))
		return -1;
	if(send_recv_couple_value(inter_tbl, value_src, value_dst, comm_src, comm_dst, intercomm))
		return -1;

	return 0;
}


/*------------------------------------------------------------------------------------------------*/

static int
set_default_value(struct hecmw_couple_value *couple_value, struct hecmw_couple_values *values_src)
{
	int nmemb, i;

	if(couple_value->item_type == HECMW_COUPLE_NODE_GROUP) {			/* Node Group		*/
		nmemb = couple_value->n * couple_value->n_dof;
		for(i=0; i<nmemb; i++) {
			values_src->node->value[i] = couple_value->value[i];
		}

	} else if(couple_value->item_type == HECMW_COUPLE_ELEMENT_GROUP) {	/* Element Group 	*/
		nmemb = couple_value->n * couple_value->n_dof;
		for(i=0; i<nmemb; i++) {
			values_src->elem->value[i] = couple_value->value[i];
		}

	} else if(couple_value->item_type == HECMW_COUPLE_SURFACE_GROUP) {	/* Surface Group	*/
		nmemb = couple_value->n * couple_value->n_dof;
		for(i=0; i<nmemb; i++) {
			values_src->surf->value[i] = couple_value->value[i];
		}

	} else {															/* error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_GRPTYPE, "");
		return -1;
	}

	HECMW_free(couple_value->item);
	HECMW_free(couple_value->value);
	couple_value->n         = 0;
	couple_value->n_dof     = 0;
	couple_value->item_type = HECMW_COUPLE_GROUP_UNDEF;
	couple_value->item      = NULL;
	couple_value->value     = NULL;

	return 0;
}



static int
init_send_value(const struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_value *couple_value,
		struct hecmw_couple_value *value_send)
{
	int nmemb, i;

	/* n */
	value_send->n = inter_tbl->export_index[inter_tbl->n_neighbor_pe_export];

	/* n_dof */
	value_send->n_dof = couple_value->n_dof;

	/* item_type */
	value_send->item_type = HECMW_COUPLE_NODE_GROUP;

	/* item */
	value_send->item = (int *)HECMW_malloc(sizeof(int)*(value_send->n+1));
	if(value_send->item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<value_send->n; i++) {
		value_send->item[i] = inter_tbl->export_item[i];
	}

	/* value */
	nmemb = value_send->n * value_send->n_dof + 1;
	value_send->value = (double *)HECMW_malloc(sizeof(double)*nmemb);
	if(value_send->value == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<nmemb; i++) {
		value_send->value[i] = 0.0;
	}

	return 0;
}



static int
init_recv_value(const struct hecmw_couple_inter_iftable *inter_tbl,
		struct hecmw_couple_value *value_recv)
{
	int i;

	/* n */
	value_recv->n = inter_tbl->import_index[inter_tbl->n_neighbor_pe_import];

	/* n_dof */
	value_recv->n_dof = 0;

	/* item_type */
	value_recv->item_type = HECMW_COUPLE_NODE_GROUP;

	/* item */
	value_recv->item = (int *)HECMW_malloc(sizeof(int)*(value_recv->n+1));
	if(value_recv->item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<value_recv->n; i++) {
		value_recv->item[i] = inter_tbl->import_item[i];
	}

	/* value */
	value_recv->value = NULL;

	return 0;
}



static int
set_result_value(struct hecmw_couple_value *couple_value, struct hecmw_couple_values *values_dst)
{
	if(values_dst->item_type == HECMW_COUPLE_NODE_GROUP) {				/* Node Group		*/
		couple_value->n         = values_dst->node->n;
		couple_value->n_dof     = values_dst->node->n_dof;
		couple_value->item_type = values_dst->node->item_type;
		couple_value->item      = values_dst->node->item;
		couple_value->value     = values_dst->node->value;
		values_dst->node->item  = NULL;
		values_dst->node->value = NULL;

	} else if(values_dst->item_type == HECMW_COUPLE_ELEMENT_GROUP) {	/* Element Group	*/
		couple_value->n         = values_dst->elem->n;
		couple_value->n_dof     = values_dst->elem->n_dof;
		couple_value->item_type = values_dst->elem->item_type;
		couple_value->item      = values_dst->elem->item;
		couple_value->value     = values_dst->elem->value;
		values_dst->elem->item  = NULL;
		values_dst->elem->value = NULL;

	} else if(values_dst->item_type == HECMW_COUPLE_SURFACE_GROUP) {	/* Surface Group	*/
		couple_value->n         = values_dst->surf->n;
		couple_value->n_dof     = values_dst->surf->n_dof;
		couple_value->item_type = values_dst->surf->item_type;
		couple_value->item      = values_dst->surf->item;
		couple_value->value     = values_dst->surf->value;
		values_dst->surf->item  = NULL;
		values_dst->surf->value = NULL;

	} else {															/* error			*/
		HECMW_set_error(HECMWCPL_E_INVALID_GRPTYPE, "");
		return -1;
	}

	return 0;
}



static int
init_values(const struct hecmw_couple_boundary *boundary,
		const struct hecmw_couple_value *couple_value,
		struct hecmw_couple_values *couple_values)
{
	int nmemb, i;

	couple_values->item_type = boundary->data_type;

	/* Node */
	if(boundary->node->n > 0) {
		couple_values->node->n         = boundary->node->n;
		couple_values->node->n_dof     = couple_value->n_dof;
		couple_values->node->item_type = HECMW_COUPLE_NODE_GROUP;

		couple_values->node->item = (int *)HECMW_malloc(sizeof(int)*couple_values->node->n);
		if(couple_values->node->item == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<couple_values->node->n; i++) {
			couple_values->node->item[i] = boundary->node->item[i];
		}

		nmemb = couple_values->node->n * couple_values->node->n_dof;
		couple_values->node->value = (double *)HECMW_malloc(sizeof(double)*nmemb);
		if(couple_values->node->value == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<nmemb; i++) {
			couple_values->node->value[i] = 0.0;
		}
	}

	/* Element */
	if(boundary->elem->n > 0) {
		couple_values->elem->n         = boundary->elem->n;
		couple_values->elem->n_dof     = couple_value->n_dof;
		couple_values->elem->item_type = HECMW_COUPLE_ELEMENT_GROUP;

		couple_values->elem->item = (int *)HECMW_malloc(sizeof(int)*couple_values->elem->n);
		if(couple_values->elem->item == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<couple_values->elem->n; i++) {
			couple_values->elem->item[i] = boundary->elem->item[i];
		}

		nmemb = couple_values->elem->n * couple_values->elem->n_dof;
		couple_values->elem->value = (double *)HECMW_malloc(sizeof(double)*nmemb);
		if(couple_values->elem->value == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<nmemb; i++) {
			couple_values->elem->value[i] = 0.0;
		}
	}

	/* Surface */
	if(boundary->surf->n > 0) {
		couple_values->surf->n         = boundary->surf->n;
		couple_values->surf->n_dof     = couple_value->n_dof;
		couple_values->surf->item_type = HECMW_COUPLE_SURFACE_GROUP;

		couple_values->surf->item = (int *)HECMW_malloc(sizeof(int)*couple_values->surf->n*2);
		if(couple_values->surf->item == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<couple_values->surf->n; i++) {
			couple_values->surf->item[2*i]   = boundary->surf->item[2*i];
			couple_values->surf->item[2*i+1] = boundary->surf->item[2*i+1];
		}

		nmemb = couple_values->surf->n * couple_values->surf->n_dof;
		couple_values->surf->value = (double *)HECMW_malloc(sizeof(double)*nmemb);
		if(couple_values->surf->value == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(i=0; i<nmemb; i++) {
			couple_values->surf->value[i] = 0.0;
		}
	}

	return 0;
}


/*================================================================================================*/

extern int
HECMW_couple(const char *boundary_id, struct hecmw_couple_value *couple_value)
{
	struct hecmw_couple_values *values_src = NULL, *values_dst = NULL;
	struct hecmw_couple_value *value_send = NULL, *value_recv = NULL;
	struct hecmw_couple_info *couple_info = NULL;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple(): 'boundary_id' is NULL");
		return HECMW_ERROR;
	}

	if((couple_info = HECMW_couple_get_info(boundary_id)) == NULL) goto error;

	if(couple_info->comm_src->is_member) {
		if(couple_value == NULL) {
			HECMW_set_error(HECMWCPL_E_INVALID_ARG, "HECMW_couple(): 'couple_value' is NULL");
			goto error;
		}
	}
	if(couple_info->comm_dst->is_member) {
		if(couple_value == NULL) {
			if((couple_value = HECMW_couple_alloc_couple_value()) == NULL) goto error;
		}
	}

	if((value_send = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	if((value_recv = HECMW_couple_alloc_couple_value()) == NULL) goto error;
	if(couple_info->comm_src->is_member) {
		if(init_send_value(couple_info->inter_tbl, couple_value, value_send)) goto error;
	}
	if(couple_info->comm_dst->is_member) {
		if(init_recv_value(couple_info->inter_tbl, value_recv)) goto error;
	}

	/*
	 * pre-processing
	 */
	if(couple_info->comm_src->is_member) {
		if((values_src = alloc_couple_values()) == NULL) goto error;

		if(init_values(couple_info->boundary_src, couple_value, values_src)) goto error;
		if(set_default_value(couple_value, values_src)) goto error;
		if(pre_interpolation(couple_info->ip_list_pre, values_src)) goto error;
		if(update_import_node_value(couple_info->comm_src,
					couple_info->intra_tbl_src, values_src->node)) goto error;
		if(main_interpolation(couple_info->ip_list_main, values_src, value_send)) goto error;

		free_couple_values(values_src);
	}

	/* main interpolation */
	if(couple_info->comm_src->is_member || couple_info->comm_dst->is_member) {
		if(send_recv(couple_info->inter_tbl, value_send, value_recv,
					couple_info->comm_src, couple_info->comm_dst,
					couple_info->intercomm)) goto error;
	}
	HECMW_couple_free_couple_value(value_send);

	/* post-processing */
	if(couple_info->comm_dst->is_member) {
		if((values_dst = alloc_couple_values()) == NULL) goto error;

		if(init_values(couple_info->boundary_dst, value_recv, values_dst)) goto error;
		if(set_default_value(value_recv, values_dst)) goto error;
		if(post_interpolation(couple_info->ip_list_post, values_dst)) goto error;
		if(set_result_value(couple_value, values_dst)) goto error;

		free_couple_values(values_dst);
	}

	HECMW_couple_free_couple_value(value_recv);
	return HECMW_SUCCESS;

error:
	free_couple_values(values_src);
	free_couple_values(values_dst);
	HECMW_couple_free_couple_value(value_send);
	HECMW_couple_free_couple_value(value_recv);
	return HECMW_ERROR;
}
