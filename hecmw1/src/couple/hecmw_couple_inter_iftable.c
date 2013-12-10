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
#include <assert.h>
#include <errno.h>

#include "hecmw_struct.h"
#include "hecmw_msgno.h"
#include "hecmw_common_define.h"
#include "hecmw_error.h"
#include "hecmw_comm.h"

#include "hecmw_couple_comm.h"
#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_judge.h"
#include "hecmw_couple_boundary_info.h"
#include "hecmw_couple_bounding_box.h"
#include "hecmw_couple_background_cell.h"
#include "hecmw_couple_inter_iftable.h"
#include "hecmw_couple_mapped_point.h"


#define MAX_NODE_SIZE 20


#define INFINITE (1.0E+37)


#define EPS (1.0E-06)


struct boundary_node_info {
	
	int n_node;
	
	double *coord;
};


struct map_info {
	
	int n;
	
	int *id;
	
	int *n_positive;
	
	double *dot_product;
	
	double *distance;
};


struct link_list {
	int item;				
	struct link_list *next;	
};


struct import_info {
	int n;					
	struct link_list *list;	
};


struct map_info_to_bgcell {
	int *index;		
	int *id;		
};


struct mapping_info {
	int n;			
	int *index;		
	int *pe;		
	int *id;		
};


struct link_list_map {
	int id;
	int item;
	struct link_list_map *next;
};

/*================================================================================================*/

static void
free_link_list(struct link_list *r)
{
	struct link_list *p, *q;

	p = r;
	while(p) {
		q = p->next;
		HECMW_free(p);
		p = q;
	}
	r = NULL;
}



static void
free_link_list_map(struct link_list_map *r)
{
	struct link_list_map *p, *q;

	p = r;
	while(p) {
		q = p->next;
		HECMW_free(p);
		p = q;
	}
	r = NULL;
}



static struct map_info *
alloc_struct_map_info(void)
{
	struct map_info *p = NULL;

	p = (struct map_info *)HECMW_malloc(sizeof(struct map_info));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	p->n           = 0;
	p->id          = NULL;
	p->n_positive  = NULL;
	p->dot_product = NULL;
	p->distance    = NULL;

	return p;
}



static void
free_struct_map_info(struct map_info *p)
{
	if(p == NULL) return;

	HECMW_free(p->id);
	HECMW_free(p->n_positive);
	HECMW_free(p->dot_product);
	HECMW_free(p->distance);
	HECMW_free(p);
	p = NULL;
}



extern void
HECMW_couple_free_inter_iftable(struct hecmw_couple_inter_iftable *p)
{
	if(p == NULL) return;

	HECMW_free(p->neighbor_pe_import);
	HECMW_free(p->import_index);
	HECMW_free(p->import_item);
	HECMW_free(p->neighbor_pe_export);
	HECMW_free(p->export_index);
	HECMW_free(p->export_item);
	HECMW_free(p);
	p = NULL;
}


extern struct hecmw_couple_inter_iftable *
HECMW_couple_alloc_inter_iftable(void)
{
	struct hecmw_couple_inter_iftable *p = NULL;
	int size;

	size = sizeof(struct hecmw_couple_inter_iftable);
	p = (struct hecmw_couple_inter_iftable *)HECMW_malloc(size);
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	p->n_neighbor_pe_import = 0;
	p->neighbor_pe_import   = NULL;
	p->import_index         = NULL;
	p->import_item          = NULL;
	p->n_neighbor_pe_export = 0;
	p->neighbor_pe_export   = NULL;
	p->export_index         = NULL;
	p->export_item          = NULL;

	return p;
}



extern void
HECMW_couple_print_inter_iftable(const struct hecmw_couple_inter_iftable *p, FILE *fp)
{
	int i, j;

	fprintf(fp, "*** Interface Table for Inter-communication\n");

	fprintf(fp, "number of neighbor processes for import: %d\n", p->n_neighbor_pe_import);
	fprintf(fp, "neighbor processes for import:\n");
	for(i=0; i<p->n_neighbor_pe_import; i++) {
		fprintf(fp, "%d%c", p->neighbor_pe_import[i], (i+1)%10 ? ' ' : '\n');
	}
	if(i%10) fprintf(fp, "\n");

	fprintf(fp, "number of neighbor processes for export: %d\n", p->n_neighbor_pe_export);
	fprintf(fp, "neighbor processes for export:\n");
	for(i=0; i<p->n_neighbor_pe_export; i++) {
		fprintf(fp, "%d%c", p->neighbor_pe_export[i], (i+1)%10 ? ' ' : '\n');
	}
	if(i%10) fprintf(fp, "\n");

	fprintf(fp, "import index:\n");
	for(i=0; i<p->n_neighbor_pe_import; i++) {
		fprintf(fp, "%d%c", p->import_index[i], (i+1)%10 ? ' ' : '\n');
	}
	if(i%10) fprintf(fp, "\n");

	fprintf(fp, "import item:\n");
	for(i=0; i<p->n_neighbor_pe_import; i++) {
		for(j=p->import_index[i]; j<p->import_index[i+1]; j++) {
			fprintf(fp, "%d%c", p->import_item[j], (j+1)%10 ? ' ' : '\n');
		}
	}
	if(j%10) fprintf(fp, "\n");

	fprintf(fp, "export index:\n");
	for(i=0; i<p->n_neighbor_pe_export; i++) {
		fprintf(fp, "%d%c", p->export_index[i], (i+1)%10 ? ' ' : '\n');
	}
	if(i%10) fprintf(fp, "\n");

	fprintf(fp, "export item:\n");
	for(i=0; i<p->n_neighbor_pe_export; i++) {
		for(j=p->export_index[i]; j<p->export_index[i+1]; j++) {
			fprintf(fp, "%d%c", p->export_item[j], (j+1)%10 ? ' ' : '\n');
		}
	}
	if(j%10) fprintf(fp, "\n");
}

/*================================================================================================*/

static int
bcast_mapped_point_d2s_n(const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_mapped_point **all_mapped_point)
{
	int n_pe_send = 0, n_pe_recv = 0, *pe_send = NULL, *pe_recv = NULL;
	int sendbuf_size = 0, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i;

	/*
	 * send buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_pe_send = comm_src->psize;

		/* communication processes */
		pe_send = (int *)HECMW_malloc(sizeof(int)*n_pe_send);
		if(pe_send == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_send; i++) {
			pe_send[i] = comm_src->ranks[i];
		}

		/* size of send buffer */
		sendbuf_size = 1;

		/* send buffer */
		sendbuf = (int *)HECMW_malloc(sizeof(int)*(sendbuf_size+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		sendbuf[0] = mapped_point->n;
	}

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_pe_recv = comm_dst->psize;

		/* communication processes */
		pe_recv = (int *)HECMW_malloc(sizeof(int)*n_pe_recv);
		if(pe_recv == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			pe_recv[i] = comm_dst->ranks[i];
		}

		/* size of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_pe_recv+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + 1;
		}

		/* receive buffer */
		recvbuf = (int *)HECMW_calloc(recvbuf_index[n_pe_recv]+1, sizeof(int));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * broadcasting
	 */
	rtc = HECMW_couple_bcast(n_pe_send, pe_send, sendbuf_size, sendbuf,
			n_pe_recv, pe_recv, recvbuf_index, recvbuf, HECMW_INT, intercomm->comm);
	if(rtc != 0)  goto error;

	if(comm_src->is_member) {
		for(i=0; i<comm_dst->psize; i++) {
			all_mapped_point[i]->n = recvbuf[i];
		}
	}

	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);
	return 0;

error:
	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);
	return -1;
}



static int
bcast_mapped_point_d2s_coord(const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_mapped_point **all_mapped_point)
{
	int n_pe_send = 0, n_pe_recv = 0, *pe_send = NULL, *pe_recv = NULL;
	int sendbuf_size = 0, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int node, size, rtc, i, j;

	/*
	 * send buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_pe_send = comm_src->psize;

		/* communication processes */
		pe_send = (int *)HECMW_malloc(sizeof(int)*n_pe_send);
		if(pe_send == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_send; i++) {
			pe_send[i] = comm_src->ranks[i];
		}

		/* size of send buffer */
		sendbuf_size = mapped_point->n*3;

		/* send buffer */
		sendbuf = (double *)HECMW_malloc(sizeof(double)*(sendbuf_size+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<mapped_point->n*3; i++) {
			sendbuf[i] = mapped_point->coord[i];
		}
	}

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_pe_recv = comm_dst->psize;

		/* communication processes */
		pe_recv = (int *)HECMW_malloc(sizeof(int)*n_pe_recv);
		if(pe_recv == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			pe_recv[i] = comm_dst->ranks[i];
		}

		/* size of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_pe_recv+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + all_mapped_point[i]->n * 3;
		}

		/* receive buffer */
		recvbuf = (double *)HECMW_malloc(sizeof(double)*(recvbuf_index[n_pe_recv]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/* broadcast */
	rtc = HECMW_couple_bcast(n_pe_send, pe_send, sendbuf_size, sendbuf,
			n_pe_recv, pe_recv, recvbuf_index, recvbuf, HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0)  goto error;

	/* coordinates */
	if(comm_src->is_member) {
		for(i=0; i<comm_dst->psize; i++) {
			size = recvbuf_index[i+1] - recvbuf_index[i];
			HECMW_assert(size == all_mapped_point[i]->n * 3);
			if(size == 0) continue;

			all_mapped_point[i]->coord = (double *)HECMW_malloc(sizeof(double)*size);
			if(all_mapped_point[i]->coord == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}

			for(j=recvbuf_index[i]; j<recvbuf_index[i+1]; j++) {
				all_mapped_point[i]->coord[j-recvbuf_index[i]] = recvbuf[j];
			}
		}
	}

	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);
	return 0;

error:
	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);

	return -1;
}



static int
bcast_mapped_point_d2s(const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_mapped_point **all_mapped_point)
{
	int rtc, i;

	if(comm_src->is_member) {
		for(i=0; i<comm_dst->psize; i++) {
			all_mapped_point[i] = HECMW_couple_alloc_mapped_point();
			if(all_mapped_point[i] == NULL) goto error;
		}
	}

	if(bcast_mapped_point_d2s_n(mapped_point, comm_src, comm_dst, intercomm, all_mapped_point))
		goto error;
	if(bcast_mapped_point_d2s_coord(mapped_point, comm_src, comm_dst, intercomm, all_mapped_point))
		goto error;

	return 0;

error:
	if(all_mapped_point) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_couple_free_mapped_point(all_mapped_point[i]);
			all_mapped_point[i] = NULL;
		}
	}
	return -1;
}


/*================================================================================================*/

static int
bcast_bbox_d2s(const struct hecmw_couple_bounding_box *bbox_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_box *all_dst_bbox)
{
	int n_pe_send = 0, n_pe_recv = 0, *pe_send = NULL, *pe_recv = NULL;
	int sendbuf_size = 0, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i;
	int N_BOUNDING_BOX_MEMBER = 6;

	/*
	 * send buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of neighboring processes */
		n_pe_send = comm_src->psize;

		/* neighboring processes */
		pe_send = (int *)HECMW_malloc(sizeof(int)*n_pe_send);
		if(pe_send == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_send; i++) {
			pe_send[i] = comm_src->ranks[i];
		}

		/* size of send buffer */
		sendbuf_size = N_BOUNDING_BOX_MEMBER;

		/* send buffer */
		sendbuf = (double *)HECMW_malloc(sizeof(double)*(sendbuf_size+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}

		sendbuf[0] = bbox_dst->enlarged->min_x;
		sendbuf[1] = bbox_dst->enlarged->min_y;
		sendbuf[2] = bbox_dst->enlarged->min_z;
		sendbuf[3] = bbox_dst->enlarged->max_x;
		sendbuf[4] = bbox_dst->enlarged->max_y;
		sendbuf[5] = bbox_dst->enlarged->max_z;
	}

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_pe_recv = comm_dst->psize;

		/* communication processes */
		pe_recv = (int *)HECMW_malloc(sizeof(int)*n_pe_recv);
		if(pe_recv == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			pe_recv[i] = comm_dst->ranks[i];
		}

		/* size of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_pe_recv+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_pe_recv; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + N_BOUNDING_BOX_MEMBER;
		}

		/* receive buffer */
		recvbuf = (double *)HECMW_malloc(sizeof(double)*(recvbuf_index[n_pe_recv]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_bcast(n_pe_send, pe_send, sendbuf_size, sendbuf,
			n_pe_recv, pe_recv, recvbuf_index, recvbuf, HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0)  goto error;

	/*
	 * set received items
	 */
	if(comm_src->is_member) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_assert(recvbuf_index[i+1] - recvbuf_index[i] == N_BOUNDING_BOX_MEMBER);

			all_dst_bbox[i].min_x = recvbuf[recvbuf_index[i]  ];
			all_dst_bbox[i].min_y = recvbuf[recvbuf_index[i]+1];
			all_dst_bbox[i].min_z = recvbuf[recvbuf_index[i]+2];
			all_dst_bbox[i].max_x = recvbuf[recvbuf_index[i]+3];
			all_dst_bbox[i].max_y = recvbuf[recvbuf_index[i]+4];
			all_dst_bbox[i].max_z = recvbuf[recvbuf_index[i]+5];
		}
	}

	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);
	return 0;

error:
	HECMW_free(sendbuf);
	HECMW_free(pe_send);
	HECMW_free(recvbuf_index);
	HECMW_free(recvbuf);
	HECMW_free(pe_recv);
	return -1;
}


/*------------------------------------------------------------------------------------------------*/

static int
check_bbox_within_bbox(const struct hecmw_couple_box *box_src,
		const struct hecmw_couple_box *box_dst)
{
	if(box_dst->min_x <= box_src->max_x && box_dst->max_x >= box_src->min_x &&
			box_dst->min_y <= box_src->max_y && box_dst->max_y >= box_src->min_y &&
			box_dst->min_z <= box_src->max_z && box_dst->max_z >= box_src->min_z) {
		return 1;
	}

	return 0;
}



static int
check_node_within_bbox(const struct hecmw_couple_box *box_src,
		const struct hecmw_couple_mapped_point *mapped_point, int *is_candidate)
{
	double coord_x, coord_y, coord_z;
	int n, i;

	for(n=0, i=0; i<mapped_point->n; i++) {
		coord_x = mapped_point->coord[3*i];
		coord_y = mapped_point->coord[3*i+1];
		coord_z = mapped_point->coord[3*i+2];

		if(coord_x >= box_src->min_x && coord_x <= box_src->max_x &&
				coord_y >= box_src->min_y && coord_y <= box_src->max_y &&
				coord_z >= box_src->min_z && coord_z <= box_src->max_z) {
			is_candidate[i] = 1;
			n++;
		} else {
			is_candidate[i] = 0;
		}
	}

	return (n != 0) ? 1 : 0;
}



static int
set_candidate_node(const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_bounding_box *bbox_dst,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct hecmw_couple_mapped_point **all_mapped_point,int **is_candidate)
{
	struct hecmw_couple_box *all_dst_bbox = NULL;
	int size, rtc, i;

	if(comm_src->is_member) {
		size = sizeof(struct hecmw_couple_box) * comm_dst->psize;
		all_dst_bbox = (struct hecmw_couple_box *)HECMW_malloc(size);
		if(all_dst_bbox == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			all_dst_bbox[i].min_x = 0.0;
			all_dst_bbox[i].min_y = 0.0;
			all_dst_bbox[i].min_z = 0.0;
			all_dst_bbox[i].max_x = 0.0;
			all_dst_bbox[i].max_y = 0.0;
			all_dst_bbox[i].max_z = 0.0;
		}
	}

	/* broadcast bounding box information from destination unit to source unit */
	if(bcast_bbox_d2s(bbox_dst, comm_src, comm_dst, intercomm, all_dst_bbox)) goto error;

	/* set candidate mapped nodes */
	if(comm_src->is_member) {
		for(i=0; i<comm_dst->psize; i++) {
			if(all_mapped_point[i]->n == 0) continue;

			is_candidate[i] = (int *)HECMW_calloc(all_mapped_point[i]->n, sizeof(int));
			if(is_candidate[i] == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}

			if(check_bbox_within_bbox(bbox_src->enlarged, &all_dst_bbox[i])) {
				check_node_within_bbox(bbox_src->enlarged, all_mapped_point[i], is_candidate[i]);
			}
		}
	}

	HECMW_free(all_dst_bbox);
	return 0;

error:
	HECMW_free(all_dst_bbox);
	if(is_candidate) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_free(is_candidate[i]);
			is_candidate[i] = NULL;
		}
	}
	return -1;
}

/*================================================================================================*/

static struct map_info_to_bgcell *
map_node_to_bgcell(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_background_cell *bgcell_src)
{
	struct map_info_to_bgcell *map_to_bgcell = NULL;
	struct link_list *bgcell_node_list = NULL, *p;
	double coord_x, coord_y, coord_z;
	int bgcell, node, ic, i, j, l, m, n;

	map_to_bgcell = (struct map_info_to_bgcell *)HECMW_malloc(sizeof(struct map_info_to_bgcell));
	if(map_to_bgcell == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->index = NULL;
	map_to_bgcell->id    = NULL;

	if(bgcell_src->n == 0) return map_to_bgcell;

	/* linked list */
	bgcell_node_list = (struct link_list *)HECMW_malloc(sizeof(struct link_list) * bgcell_src->n);
	if(bgcell_node_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<bgcell_src->n; i++) {
		bgcell_node_list[i].item = 0;
		bgcell_node_list[i].next = NULL;
	}

	for(ic=0, i=0; i<boundary_src->node->n; i++) {
		node    = boundary_src->node->item[i];
		coord_x = mesh_src->node[3*(node-1)  ];
		coord_y = mesh_src->node[3*(node-1)+1];
		coord_z = mesh_src->node[3*(node-1)+2];

		l = (coord_x - bbox_src->enlarged->min_x) / bgcell_src->dx;
		m = (coord_y - bbox_src->enlarged->min_y) / bgcell_src->dy;
		n = (coord_z - bbox_src->enlarged->min_z) / bgcell_src->dz;

		if(l < 0) l = 0;
		if(m < 0) m = 0;
		if(n < 0) n = 0;
		if(l >= bgcell_src->nx) l = bgcell_src->nx - 1;
		if(m >= bgcell_src->ny) m = bgcell_src->ny - 1;
		if(n >= bgcell_src->nz) n = bgcell_src->nz - 1;

		bgcell = bgcell_src->nx * bgcell_src->ny * n + bgcell_src->nx * m + l;

		p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
		if(p == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		p->item = i;
		p->next = bgcell_node_list[bgcell].next;
		bgcell_node_list[bgcell].next = p;

		ic++;
	}

	/* compressed 1-dimensional array */
	map_to_bgcell->index = (int *)HECMW_calloc(bgcell_src->n+1, sizeof(int));
	if(map_to_bgcell->index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->id = (int *)HECMW_malloc(sizeof(int)*ic);
	if(map_to_bgcell->id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(ic=0, i=0; i<bgcell_src->n; i++) {
		p = bgcell_node_list[i].next;
		while(p) {
			map_to_bgcell->id[ic++] = p->item;
			p = p->next;
		}
		map_to_bgcell->index[i+1] = ic;
	}

	/* free */
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_node_list[i].next);
	}
	HECMW_free(bgcell_node_list);

	return map_to_bgcell;

error:
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_node_list[i].next);
	}
	HECMW_free(bgcell_node_list);
	if(map_to_bgcell) {
		HECMW_free(map_to_bgcell->index);
		HECMW_free(map_to_bgcell->id);
	}

	return NULL;
}


static struct map_info_to_bgcell *
map_elem_to_bgcell(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_background_cell *bgcell_src)
{
	struct map_info_to_bgcell *map_to_bgcell = NULL;
	struct link_list *bgcell_elem_list = NULL, *p;
	double coord_x, coord_y, coord_z, min_x, min_y, min_z, max_x, max_y, max_z;
	int min_l, min_m, min_n, max_l, max_m, max_n;
	int bgcell, elem, node, ic, i, j, l, m, n;

	map_to_bgcell = (struct map_info_to_bgcell *)HECMW_malloc(sizeof(struct map_info_to_bgcell));
	if(map_to_bgcell == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->index = NULL;
	map_to_bgcell->id    = NULL;

	if(bgcell_src->n == 0) return map_to_bgcell;

	/* linked list */
	bgcell_elem_list = (struct link_list *)HECMW_malloc(sizeof(struct link_list)*bgcell_src->n);
	if(bgcell_elem_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<bgcell_src->n; i++) {
		bgcell_elem_list[i].item = 0;
		bgcell_elem_list[i].next = NULL;
	}

	for(ic=0, i=0; i<boundary_src->elem->n; i++) {
		elem  = boundary_src->elem->item[i];
		min_x = min_y = min_z = +INFINITE;
		max_x = max_y = max_z = -INFINITE;
		for(j=mesh_src->elem_node_index[elem-1]; j<mesh_src->elem_node_index[elem]; j++) {
			node    = mesh_src->elem_node_item[j];
			coord_x = mesh_src->node[3*(node-1)  ];
			coord_y = mesh_src->node[3*(node-1)+1];
			coord_z = mesh_src->node[3*(node-1)+2];
			if(coord_x < min_x) min_x = coord_x;
			if(coord_y < min_y) min_y = coord_y;
			if(coord_z < min_z) min_z = coord_z;
			if(coord_x > max_x) max_x = coord_x;
			if(coord_y > max_y) max_y = coord_y;
			if(coord_z > max_z) max_z = coord_z;
		}
		min_x = min_x - (bbox_src->tolerance + EPS);
		min_y = min_y - (bbox_src->tolerance + EPS);
		min_z = min_z - (bbox_src->tolerance + EPS);
		max_x = max_x + (bbox_src->tolerance + EPS);
		max_y = max_y + (bbox_src->tolerance + EPS);
		max_z = max_z + (bbox_src->tolerance + EPS);

		min_l = (min_x - bbox_src->enlarged->min_x) / bgcell_src->dx;
		min_m = (min_y - bbox_src->enlarged->min_y) / bgcell_src->dy;
		min_n = (min_z - bbox_src->enlarged->min_z) / bgcell_src->dz;
		max_l = (max_x - bbox_src->enlarged->max_x) / bgcell_src->dx;
		max_m = (max_y - bbox_src->enlarged->max_y) / bgcell_src->dy;
		max_n = (max_z - bbox_src->enlarged->max_z) / bgcell_src->dz;

		if(min_l < 0) min_l = 0;
		if(min_m < 0) min_m = 0;
		if(min_n < 0) min_n = 0;
		if(max_l < 0) max_l = 0;
		if(max_m < 0) max_m = 0;
		if(max_n < 0) max_n = 0;
		if(min_l >= bgcell_src->nx) min_l = bgcell_src->nx - 1;
		if(min_m >= bgcell_src->ny) min_m = bgcell_src->ny - 1;
		if(min_n >= bgcell_src->nz) min_n = bgcell_src->nz - 1;
		if(max_l >= bgcell_src->nx) max_l = bgcell_src->nx - 1;
		if(max_m >= bgcell_src->ny) max_m = bgcell_src->ny - 1;
		if(max_n >= bgcell_src->nz) max_n = bgcell_src->nz - 1;

		for(l=min_l; l<=max_l; l++) {
			for(m=min_m; m<=max_m; m++) {
				for(n=min_n; n<=max_n; n++) {
					bgcell = bgcell_src->nx * bgcell_src->ny * n + bgcell_src->nx * m + l;

					p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
					if(p == NULL) {
						HECMW_set_error(errno, "");
						goto error;
					}
					p->item = i;
					p->next = bgcell_elem_list[bgcell].next;
					bgcell_elem_list[bgcell].next = p;

					ic++;
				}
			}
		}
	}

	/* compressed 1-dimensional array */
	map_to_bgcell->index = (int *)HECMW_calloc(bgcell_src->n+1, sizeof(int));
	if(map_to_bgcell->index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->id = (int *)HECMW_malloc(sizeof(int)*ic);
	if(map_to_bgcell->id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(ic=0, i=0; i<bgcell_src->n; i++) {
		p = bgcell_elem_list[i].next;
		while(p) {
			map_to_bgcell->id[ic++] = p->item;
			p = p->next;
		}
		map_to_bgcell->index[i+1] = ic;
	}

	/* free */
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_elem_list[i].next);
	}
	HECMW_free(bgcell_elem_list);

	return map_to_bgcell;

error:
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_elem_list[i].next);
	}
	HECMW_free(bgcell_elem_list);
	if(map_to_bgcell) {
		HECMW_free(map_to_bgcell->index);
		HECMW_free(map_to_bgcell->id);
	}
	return NULL;
}


static struct map_info_to_bgcell *
map_surf_to_bgcell(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_background_cell *bgcell_src)
{
	struct map_info_to_bgcell *map_to_bgcell = NULL;
	struct link_list *bgcell_surf_list = NULL, *p;
	double coord_x, coord_y, coord_z, min_x, min_y, min_z, max_x, max_y, max_z;
	int min_l, min_m, min_n, max_l, max_m, max_n;
	int bgcell, elem, node, ic, i, j, l, m, n;

	/* allocation */
	map_to_bgcell = (struct map_info_to_bgcell *)HECMW_malloc(sizeof(struct map_info_to_bgcell));
	if(map_to_bgcell == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->index = NULL;
	map_to_bgcell->id    = NULL;

	if(bgcell_src->n == 0) return map_to_bgcell;

	/* linked list */
	bgcell_surf_list = (struct link_list *)HECMW_malloc(sizeof(struct link_list)*bgcell_src->n);
	if(bgcell_surf_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<bgcell_src->n; i++) {
		bgcell_surf_list[i].item = 0;
		bgcell_surf_list[i].next = NULL;
	}

	for(ic=0, i=0; i<boundary_src->surf->n; i++) {
		elem  = boundary_src->surf->item[2*i];
		min_x = min_y = min_z = +INFINITE;
		max_x = max_y = max_z = -INFINITE;
		for(j=mesh_src->elem_node_index[elem-1]; j<mesh_src->elem_node_index[elem]; j++) {
			node    = mesh_src->elem_node_item[j];
			coord_x = mesh_src->node[3*(node-1)  ];
			coord_y = mesh_src->node[3*(node-1)+1];
			coord_z = mesh_src->node[3*(node-1)+2];
			if(coord_x < min_x) min_x = coord_x;
			if(coord_y < min_y) min_y = coord_y;
			if(coord_z < min_z) min_z = coord_z;
			if(coord_x > max_x) max_x = coord_x;
			if(coord_y > max_y) max_y = coord_y;
			if(coord_z > max_z) max_z = coord_z;
		}
		min_x = min_x - (bbox_src->tolerance + EPS);
		min_y = min_y - (bbox_src->tolerance + EPS);
		min_z = min_z - (bbox_src->tolerance + EPS);
		max_x = max_x + (bbox_src->tolerance + EPS);
		max_y = max_y + (bbox_src->tolerance + EPS);
		max_z = max_z + (bbox_src->tolerance + EPS);

		min_l = (min_x - bbox_src->enlarged->min_x) / bgcell_src->dx;
		min_m = (min_y - bbox_src->enlarged->min_y) / bgcell_src->dy;
		min_n = (min_z - bbox_src->enlarged->min_z) / bgcell_src->dz;
		max_l = (max_x - bbox_src->enlarged->min_x) / bgcell_src->dx;
		max_m = (max_y - bbox_src->enlarged->min_y) / bgcell_src->dy;
		max_n = (max_z - bbox_src->enlarged->min_z) / bgcell_src->dz;

		if(min_l < 0) min_l = 0;
		if(min_m < 0) min_m = 0;
		if(min_n < 0) min_n = 0;
		if(max_l < 0) max_l = 0;
		if(max_m < 0) max_m = 0;
		if(max_n < 0) max_n = 0;
		if(min_l >= bgcell_src->nx) min_l = bgcell_src->nx - 1;
		if(min_m >= bgcell_src->ny) min_m = bgcell_src->ny - 1;
		if(min_n >= bgcell_src->nz) min_n = bgcell_src->nz - 1;
		if(max_l >= bgcell_src->nx) max_l = bgcell_src->nx - 1;
		if(max_m >= bgcell_src->ny) max_m = bgcell_src->ny - 1;
		if(max_n >= bgcell_src->nz) max_n = bgcell_src->nz - 1;

		for(l=min_l; l<=max_l; l++) {
			for(m=min_m; m<=max_m; m++) {
				for(n=min_n; n<=max_n; n++) {
					bgcell = bgcell_src->nx * bgcell_src->ny * n + bgcell_src->nx * m + l;
					p = (struct link_list *)HECMW_malloc(sizeof(struct link_list));
					if(p == NULL) {
						HECMW_set_error(errno, "");
						goto error;
					}
					p->item = i;
					p->next = bgcell_surf_list[bgcell].next;
					bgcell_surf_list[bgcell].next = p;

					ic++;
				}
			}
		}
	}

	/* compressed 1-dimensional array */
	map_to_bgcell->index = (int *)HECMW_calloc(bgcell_src->n+1, sizeof(int));
	if(map_to_bgcell->index == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->id = (int *)HECMW_malloc(sizeof(int)*ic);
	if(map_to_bgcell->id == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(ic=0, i=0; i<bgcell_src->n; i++) {
		for(p=bgcell_surf_list[i].next; p; p=p->next) {
			map_to_bgcell->id[ic++] = p->item;
		}
		map_to_bgcell->index[i+1] = ic;
	}

	/* free */
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_surf_list[i].next);
	}
	HECMW_free(bgcell_surf_list);

	return map_to_bgcell;

error:
	for(i=0; i<bgcell_src->n; i++) {
		free_link_list(bgcell_surf_list[i].next);
	}
	HECMW_free(bgcell_surf_list);
	if(map_to_bgcell) {
		HECMW_free(map_to_bgcell->index);
		HECMW_free(map_to_bgcell->id);
	}
	return NULL;
}


/*================================================================================================*/

static int
compare_mapping_info(int n_positive_old, double dot_product_old, double distance_old,
		int n_positive_new, double dot_product_new, double distance_new)
{
	/*
	 * 1st step: number of positive dot products
	 */
	if(n_positive_old < 0)  return 1;				/* update	*/
	if(n_positive_new < n_positive_old)  return 1;	/* update	*/
	if(n_positive_new > n_positive_old)  return 0;	/* keep		*/

	/*
	 * 2nd step: dot product (if number of positive dot products is positive)
	 */
	/* n_positive_new = n_positive_old > 0 */
/*	if(n_positive_new > 0) { */
		if(dot_product_old > 0 && dot_product_new < 0)  return 0;	/* keep		*/
		if(dot_product_old < 0 && dot_product_new > 0)  return 1;	/* update	*/
/*	}*/

	/*
	 * 3rd step: distance
	 */
	if(distance_new < distance_old)  return 1;	/* update	*/

	return 0;	/* keep		*/
}



static int
map_point_to_surf(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_background_cell *bgcell_src,
		const struct hecmw_couple_mapped_point *mapped_point, int *is_candidate,
		const struct map_info_to_bgcell *map_to_bgcell,
		struct map_info *map_src)
{
	double coord_x, coord_y, coord_z;
	int depth_x_before, depth_y_before, depth_z_before, depth_x_after, depth_y_after, depth_z_after;
	double dot_product, distance;
	int n_positive, elem, surf, bgcell, surf_id, nl, nm, nn, rtc, jstart, jfinal, i, j, l, m, n;

	map_src->n = mapped_point->n;
	if(map_src->n == 0) return 0;

	map_src->id = (int *)HECMW_malloc(sizeof(int)*map_src->n);
	if(map_src->id == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	map_src->n_positive = (int *)HECMW_malloc(sizeof(int)*map_src->n);
	if(map_src->n_positive == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	map_src->distance = (double *)HECMW_malloc(sizeof(double)*map_src->n);
	if(map_src->distance == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	map_src->dot_product = (double *)HECMW_malloc(sizeof(double)*map_src->n);
	if(map_src->dot_product == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i<map_src->n; i++) {
		map_src->id[i]          = -1;
		map_src->n_positive[i]  = -1;
		map_src->distance[i]    = 0.0;
		map_src->dot_product[i] = 0.0;
	}

	for(i=0; i<mapped_point->n; i++) {
		if(is_candidate[i] == 0)  continue;

		coord_x = mapped_point->coord[3*i  ];
		coord_y = mapped_point->coord[3*i+1];
		coord_z = mapped_point->coord[3*i+2];

		/* background cell which includes mapped node */
		nl = (coord_x - bbox_src->enlarged->min_x) / bgcell_src->dx;
		nm = (coord_y - bbox_src->enlarged->min_y) / bgcell_src->dy;
		nn = (coord_z - bbox_src->enlarged->min_z) / bgcell_src->dz;
		if(nl < 0) nl = 0;
		if(nm < 0) nm = 0;
		if(nn < 0) nn = 0;
		if(nl >= bgcell_src->nx) nl = bgcell_src->nx - 1;
		if(nm >= bgcell_src->ny) nm = bgcell_src->ny - 1;
		if(nn >= bgcell_src->nz) nn = bgcell_src->nz - 1;

		bgcell = bgcell_src->nx*bgcell_src->ny*nn + bgcell_src->nx*nm + nl;

		/* mapping */
		for(j=map_to_bgcell->index[bgcell]; j<map_to_bgcell->index[bgcell+1]; j++) {
			surf_id = map_to_bgcell->id[j];
			elem    = boundary_src->surf->item[2*surf_id];
			surf    = boundary_src->surf->item[2*surf_id+1];

			if(mesh_src->elem_type[elem-1] == HECMW_ETYPE_TET1) {
				n_positive = HECMW_couple_judge_tet1(mesh_src, elem, surf,
						coord_x, coord_y, coord_z, &dot_product, &distance);
			} else if(mesh_src->elem_type[elem-1] == HECMW_ETYPE_HEX1) {
				n_positive = HECMW_couple_judge_hex1(mesh_src, elem, surf,
						coord_x, coord_y, coord_z, &dot_product, &distance);
			} else {
				HECMW_set_error(HECMWCPL_E_NONSUPPORT_ETYPE, "");
				return -1;
			}

			if(n_positive >= 0) {
				rtc = compare_mapping_info(map_src->n_positive[i], map_src->dot_product[i],
						map_src->distance[i], n_positive, dot_product, distance);
				if(rtc == 1) {
					map_src->n_positive[i]  = n_positive;
					map_src->dot_product[i] = dot_product;
					map_src->distance[i]    = distance;
					map_src->id[i]          = surf_id;
				}
			}
		}
	}

	return 0;
}



static int
mapping_in_src(const struct hecmwST_local_mesh *mesh_src,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_background_cell *bgcell_src,
		const struct hecmw_couple_comm *comm_dst,
		struct hecmw_couple_mapped_point **all_mapped_point,
		int **is_candidate, struct map_info **map_src)
{
	struct map_info_to_bgcell *map_to_bgcell = NULL;
	int map_type, i;

	/* mapping boundary surface to background cell */
	map_to_bgcell = (struct map_info_to_bgcell *)HECMW_malloc(sizeof(struct map_info_to_bgcell));
	if(map_to_bgcell == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	map_to_bgcell->index = NULL;
	map_to_bgcell->id    = NULL;

	/* mapping boundary nodes on destination side to boundary surface on source side */
	map_type = HECMW_COUPLE_MAP_SURF_TO_SURF;	/*@@*/

	if(map_type == HECMW_COUPLE_MAP_SURF_TO_SURF) {
		map_to_bgcell = map_surf_to_bgcell(mesh_src, boundary_src, bbox_src, bgcell_src);
		if(map_to_bgcell == NULL) goto error;

		for(i=0; i<comm_dst->psize; i++) {
			if(map_point_to_surf(mesh_src, boundary_src, bbox_src, bgcell_src, all_mapped_point[i],
						is_candidate[i], map_to_bgcell, map_src[i])) goto error;
		}

	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_MAPTYPE, "");
		goto error;
	}

	if(map_to_bgcell) {
		HECMW_free(map_to_bgcell->index);
		HECMW_free(map_to_bgcell->id);
		HECMW_free(map_to_bgcell);
	}
	return 0;

error:
	if(map_to_bgcell) {
		HECMW_free(map_to_bgcell->index);
		HECMW_free(map_to_bgcell->id);
		HECMW_free(map_to_bgcell);
	}
	if(map_src) {
		for(i=0; i<comm_dst->psize; i++) {
			free_struct_map_info(map_src[i]);
		}
		HECMW_free(map_src);
	}
	return -1;
}


/*================================================================================================*/

static int
gather_id_s2d(struct map_info **map_src,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, struct map_info **map_dst)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i, j;

	/*
	 * send buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_send_pe = comm_dst->psize;

		/* communication processes */
		send_pe = (int *)HECMW_malloc(sizeof(int)*n_send_pe);
		if(send_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			send_pe[i] = comm_dst->ranks[i];
		}

		/* size of send buffer */
		sendbuf_index = (int *)HECMW_calloc(n_send_pe+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + map_src[i]->n;
		}

		/* send buffer */
		sendbuf = (int *)HECMW_malloc(sizeof(int)*(sendbuf_index[n_send_pe]+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			for(j=0; j<map_src[i]->n; j++) {
				sendbuf[sendbuf_index[i]+j] = map_src[i]->id[j];
			}
		}
	}

	/*
	 * receive buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_recv_pe = comm_src->psize;

		/* communication processes */
		recv_pe = (int *)HECMW_malloc(sizeof(int)*n_recv_pe);
		if(recv_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recv_pe[i] = comm_src->ranks[i];
		}

		/* size of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_recv_pe+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + map_dst[i]->n;
		}

		/* receive buffer */
		recvbuf = (int *)HECMW_malloc(sizeof(int)*(recvbuf_index[n_recv_pe]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(n_send_pe, send_pe, sendbuf_index, sendbuf,
			n_recv_pe, recv_pe, recvbuf_index, recvbuf, HECMW_INT, intercomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received data
	 */
	if(comm_dst->is_member) {
		for(i=0; i<n_recv_pe; i++) {
			map_dst[i]->id = (int *)HECMW_calloc(map_dst[i]->n, sizeof(int));
			if(map_dst[i]->id == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(j=0; j<map_dst[i]->n; j++) {
				map_dst[i]->id[j] = recvbuf[recvbuf_index[i]+j];
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
gather_n_positive_s2d(struct map_info **map_src,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, struct map_info **map_dst)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i, j;

	/*
	 * send buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_send_pe = comm_dst->psize;

		/* communication processes */
		send_pe = (int *)HECMW_malloc(sizeof(int)*n_send_pe);
		if(send_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			send_pe[i] = comm_dst->ranks[i];
		}

		/* size of send buffer */
		sendbuf_index = (int *)HECMW_calloc(n_send_pe+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + map_src[i]->n;
		}

		/* send buffer */
		sendbuf = (int *)HECMW_malloc(sizeof(int)*(sendbuf_index[n_send_pe]+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			for(j=0; j<map_src[i]->n; j++) {
				sendbuf[sendbuf_index[i]+j] = map_src[i]->n_positive[j];
			}
		}
	}

	/*
	 * receive buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_recv_pe = comm_src->psize;

		/* communication processes */
		recv_pe = (int *)HECMW_malloc(sizeof(int)*n_recv_pe);
		if(recv_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recv_pe[i] = comm_src->ranks[i];
		}

		/* size of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_recv_pe+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + map_dst[i]->n;
		}

		/* receive buffer */
		recvbuf = (int *)HECMW_malloc(sizeof(int)*(recvbuf_index[n_recv_pe]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(n_send_pe, send_pe, sendbuf_index, sendbuf,
			n_recv_pe, recv_pe, recvbuf_index, recvbuf, HECMW_INT, intercomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received data
	 */
	if(comm_dst->is_member) {
		for(i=0; i<n_recv_pe; i++) {
			map_dst[i]->n_positive = (int *)HECMW_calloc(map_dst[i]->n, sizeof(int));
			if(map_dst[i]->n_positive == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(j=0; j<map_dst[i]->n; j++) {
				map_dst[i]->n_positive[j] = recvbuf[recvbuf_index[i]+j];
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
gather_dot_product_s2d(struct map_info **map_src,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, struct map_info **map_dst)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i, j;

	/*
	 * send buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_send_pe = comm_dst->psize;

		/* communication processes */
		send_pe = (int *)HECMW_malloc(sizeof(int)*n_send_pe);
		if(send_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			send_pe[i] = comm_dst->ranks[i];
		}

		/* index of send buffer */
		sendbuf_index = (int *)HECMW_calloc(n_send_pe+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + map_src[i]->n;
		}

		/* send buffer */
		sendbuf = (double *)HECMW_malloc(sizeof(double)*(sendbuf_index[n_send_pe]+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			for(j=0; j<map_src[i]->n; j++) {
				sendbuf[sendbuf_index[i]+j] = map_src[i]->dot_product[j];
			}
		}
	}

	/*
	 * receive buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_recv_pe = comm_src->psize;

		/* communication processes */
		recv_pe = (int *)HECMW_malloc(sizeof(int)*n_recv_pe);
		if(recv_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recv_pe[i] = comm_src->ranks[i];
		}

		/* index of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_recv_pe+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + map_dst[i]->n;
		}

		/* receive buffer */
		recvbuf = (double *)HECMW_malloc(sizeof(double)*recvbuf_index[n_recv_pe]+1);
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(n_send_pe, send_pe, sendbuf_index, sendbuf,
			n_recv_pe, recv_pe, recvbuf_index, recvbuf, HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0) goto error;

	/*
	 * store received data
	 */
	if(comm_dst->is_member) {
		for(i=0; i<n_recv_pe; i++) {
			map_dst[i]->dot_product = (double *)HECMW_malloc(sizeof(double)*map_dst[i]->n);
			if(map_dst[i]->dot_product == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(j=0; j<map_dst[i]->n; j++) {
				map_dst[i]->dot_product[j] = recvbuf[recvbuf_index[i]+j];
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
gather_distance_s2d(struct map_info **map_src,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, struct map_info **map_dst)
{
	int n_send_pe = 0, n_recv_pe = 0, *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL;
	double *sendbuf = NULL, *recvbuf = NULL;
	int rtc, i, j;

	/*
	 * send buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_send_pe = comm_dst->psize;

		/* communication processes */
		send_pe = (int *)HECMW_malloc(sizeof(int)*n_send_pe);
		if(send_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			send_pe[i] = comm_dst->ranks[i];
		}

		/* index of send buffer */
		sendbuf_index = (int *)HECMW_calloc(n_send_pe+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + map_src[i]->n;
		}

		/* send buffer */
		sendbuf = (double *)HECMW_malloc(sizeof(double)*(sendbuf_index[n_send_pe]+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			for(j=0; j<map_src[i]->n; j++) {
				sendbuf[sendbuf_index[i]+j] = map_src[i]->distance[j];
			}
		}
	}

	/* receive buffer */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_recv_pe = comm_src->psize;

		/* communication processes */
		recv_pe = (int *)HECMW_malloc(sizeof(int)*n_recv_pe);
		if(recv_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recv_pe[i] = comm_src->ranks[i];
		}

		/* index of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_recv_pe+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + map_dst[i]->n;
		}

		/* receive buffer */
		recvbuf = (double *)HECMW_malloc(sizeof(double)*(recvbuf_index[n_recv_pe]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(n_send_pe, send_pe, sendbuf_index, sendbuf,
			n_recv_pe, recv_pe, recvbuf_index, recvbuf, HECMW_DOUBLE, intercomm->comm);
	if(rtc != 0)  goto error;

	/*
	 * store received data
	 */
	if(comm_dst->is_member) {
		for(i=0; i<n_recv_pe; i++) {
			map_dst[i]->distance = (double *)HECMW_malloc(sizeof(double)*map_dst[i]->n);
			if(map_dst[i]->distance == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(j=0; j<map_dst[i]->n; j++) {
				map_dst[i]->distance[j] = recvbuf[recvbuf_index[i]+j];
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


/*------------------------------------------------------------------------------------------------*/

static int
gather_map_data_s2d(const struct hecmw_couple_mapped_point *mapped_point,
		struct map_info **map_src,
		const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		struct map_info **map_dst)
{
	int map_type, rtc, i;

	if(comm_dst->is_member) {
		for(i=0; i<comm_src->psize; i++) {
			map_dst[i] = alloc_struct_map_info();
			if(map_dst[i] == NULL) goto error;

			map_dst[i]->n = mapped_point->n;
		}
	}

/*@@@	map_type = HECMW_couple_get_map_type(boundary_id); @@@*/
	map_type = HECMW_COUPLE_MAP_NODE_TO_SURF;

	if(map_type == HECMW_COUPLE_MAP_NODE_TO_SURF) {
		if(gather_id_s2d(map_src, comm_src, comm_dst, intercomm, map_dst)) goto error;
		if(gather_n_positive_s2d(map_src, comm_src, comm_dst, intercomm, map_dst)) goto error;
		if(gather_dot_product_s2d(map_src, comm_src, comm_dst, intercomm, map_dst)) goto error;
		if(gather_distance_s2d(map_src, comm_src, comm_dst, intercomm, map_dst)) goto error;

	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_MAPTYPE, "");
		goto error;
	}

	return 0;

error:
	if(map_dst) {
		for(i=0; i<comm_src->psize; i++) {
			free_struct_map_info(map_dst[i]);
			map_dst[i] = NULL;
		}
	}
	return -1;
}

/*================================================================================================*/

static struct link_list_map *
set_mapping_surf(const struct hecmw_couple_mapped_point *mapped_point,
		const struct hecmw_couple_comm *comm_src, struct map_info **map_dst)
{
	struct link_list_map *mapping_data_list = NULL, *p;
	double _dot_product, _distance;
	int _n_positive, _item, pe_index, index, size, rtc, n, i, j;

	size = sizeof(struct link_list_map)*comm_src->psize;
	mapping_data_list = (struct link_list_map *)HECMW_malloc(size);
	if(mapping_data_list == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<comm_src->psize; i++) {
		mapping_data_list[i].next = NULL;
		mapping_data_list[i].id   = -1;
		mapping_data_list[i].item = -1;
	}

	for(n=0, i=0; i<mapped_point->n; i++) {
		_n_positive = -1;

		for(j=0; j<comm_src->psize; j++) {
			if(map_dst[j]->n_positive[i] >= 0) {
				rtc = compare_mapping_info(_n_positive, _dot_product, _distance,
						map_dst[j]->n_positive[i], map_dst[j]->dot_product[i],
						map_dst[j]->distance[i]);
				if(rtc == 1) {
					pe_index     = j;
					_item        = map_dst[j]->id[i];
					_n_positive  = map_dst[j]->n_positive[i];
					_dot_product = map_dst[j]->dot_product[i];
					_distance    = map_dst[j]->distance[i];
				}
			}
		}

		if(_n_positive >= 0) {
			p = (struct link_list_map *)HECMW_malloc(sizeof(struct link_list_map));
			if(p == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			p->next = mapping_data_list[pe_index].next;
			p->id   = mapped_point->id[n++];
			p->item = _item;
			mapping_data_list[pe_index].next = p;
		}
	}

	return mapping_data_list;

error:
	if(mapping_data_list) {
		for(i=0; i<comm_src->psize; i++) {
			free_link_list_map(mapping_data_list[i].next);
		}
		HECMW_free(mapping_data_list);
	}

	return NULL;
}


/*================================================================================================*/

static int
set_n_neighbor_pe_import(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, struct link_list_map *mapping_data_list)
{
	int n, i;

	for(n=0, i=0; i<comm_src->psize; i++) {
		if(mapping_data_list[i].next) n++;
	}
	inter_tbl->n_neighbor_pe_import = n;

	return 0;
}



static int
set_neighbor_pe_import(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, struct link_list_map *mapping_data_list)
{
	int size, n, i;

	if(inter_tbl->n_neighbor_pe_import == 0) return 0;

	size = sizeof(int)*inter_tbl->n_neighbor_pe_import;
	inter_tbl->neighbor_pe_import = (int *)HECMW_malloc(size);
	if(inter_tbl->neighbor_pe_import == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(n=0, i=0; i<comm_src->psize; i++) {
		if(mapping_data_list[i].next) {
			inter_tbl->neighbor_pe_import[n++] = comm_src->ranks[i];
		}
	}
	HECMW_assert(n == inter_tbl->n_neighbor_pe_import);

	return 0;
}



static int
set_import_index(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, struct link_list_map *mapping_data_list)
{
	struct link_list_map *p;
	int nmemb, m, n, i;

	nmemb = inter_tbl->n_neighbor_pe_import+1;
	inter_tbl->import_index = (int *)HECMW_calloc(nmemb, sizeof(int));
	if(inter_tbl->import_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	for(m=0, n=0, i=0; i<comm_src->psize; i++) {
		if(mapping_data_list[i].next) {
			for(p=mapping_data_list[i].next; p; p=p->next) {
				n++;
			}
			inter_tbl->import_index[++m] = n;
		}
	}
	HECMW_assert(m == inter_tbl->n_neighbor_pe_import);

	return 0;
}


static int
set_import_item(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, struct link_list_map *mapping_data_list)
{
	struct link_list_map *p;
	int size, n, i;

	size = sizeof(int)*inter_tbl->import_index[inter_tbl->n_neighbor_pe_import];
	inter_tbl->import_item = (int *)HECMW_malloc(size);
	if(inter_tbl->import_item == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}

	for(n=0, i=0; i<comm_src->psize; i++) {
		if(mapping_data_list[i].next) {
			for(p=mapping_data_list[i].next; p; p=p->next) {
				inter_tbl->import_item[n++] = p->id;
			}
			HECMW_assert(n == inter_tbl->import_index[i+1]);
		}
	}

	return 0;
}

/*------------------------------------------------------------------------------------------------*/

static int *
set_export_data(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, struct link_list_map *mapping_data_list)
{
	struct link_list_map *p;
	int *export_data = NULL;
	int size, n, i;

	size = sizeof(int)*inter_tbl->import_index[inter_tbl->n_neighbor_pe_import];
	export_data = (int *)HECMW_malloc(size);
	if(export_data == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	for(n=0, i=0; i<comm_src->psize; i++) {
		if(mapping_data_list[i].next) {
			for(p=mapping_data_list[i].next; p; p=p->next) {
				export_data[n++] = p->item;
			}
			HECMW_assert(n == inter_tbl->import_index[i+1]);
		}
	}

	return export_data;
}



static int
set_neighbor_info_import(struct hecmw_couple_inter_iftable *inter_tbl,
		struct link_list_map *mapping_data_list, const struct hecmw_couple_comm *comm_src)
{
	int i;

	/* number of neighboring domains */
	if(set_n_neighbor_pe_import(inter_tbl, comm_src, mapping_data_list)) return -1;

	/* neighboring domain */
	if(set_neighbor_pe_import(inter_tbl, comm_src, mapping_data_list)) return -1;

	/* index for import node */
	if(set_import_index(inter_tbl, comm_src, mapping_data_list)) return -1;

	/* import node */
	if(set_import_item(inter_tbl, comm_src, mapping_data_list)) return -1;

	return 0;
}


/*------------------------------------------------------------------------------------------------*/

static int
set_neighbor_pe_export(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, const struct link_list_map *mapping_data_list)
{
	int n_send_pe = 0, n_recv_pe = 0;
	int *send_pe = NULL, *recv_pe = NULL;
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int size, rtc, n, i;

	/*
	 * send buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* number of communication processes */
		n_send_pe = comm_src->psize;

		/* communication processes */
		send_pe = (int *)HECMW_malloc(sizeof(int)*n_send_pe);
		if(send_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_src->psize; i++) {
			send_pe[i] = comm_src->ranks[i];
		}

		/* index of send buffer */
		sendbuf_index = (int *)HECMW_calloc(n_send_pe+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + 1;
		}

		/* send buffer */
		sendbuf = (int *)HECMW_malloc(sizeof(int)*(sendbuf_index[n_send_pe]+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_send_pe; i++) {
			if(mapping_data_list[i].next) {
				sendbuf[i] = 1;
			} else {
				sendbuf[i] = 0;
			}
		}
	}

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* number of communication processes */
		n_recv_pe = comm_dst->psize;

		/* communication processes */
		recv_pe = (int *)HECMW_malloc(sizeof(int)*n_recv_pe);
		if(recv_pe == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			recv_pe[i] = comm_dst->ranks[i];
		}

		/* index of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(n_recv_pe+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<n_recv_pe; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + 1;
		}

		/* receive buffer */
		recvbuf = (int *)HECMW_malloc(sizeof(int)*(recvbuf_index[n_recv_pe]+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/* send and receive */
	rtc = HECMW_couple_inter_send_recv(n_send_pe, send_pe, sendbuf_index, sendbuf,
			n_recv_pe, recv_pe, recvbuf_index, recvbuf, HECMW_INT, intercomm->comm);
	if(rtc != 0)  goto error;

	/*
	 * set neighbor process for inter-communication (source unit)
	 */
	if(comm_src->is_member) {
		for(n=0, i=0; i<comm_dst->psize; i++) {
			if(recvbuf[i] > 0) n++;
		}
		inter_tbl->n_neighbor_pe_export = n;

		if(inter_tbl->n_neighbor_pe_export != 0) {
			size = sizeof(int)*inter_tbl->n_neighbor_pe_export;
			inter_tbl->neighbor_pe_export = (int *)HECMW_malloc(size);
			if(inter_tbl->neighbor_pe_export == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(n=0, i=0; i<comm_dst->psize; i++) {
				if(recvbuf[i] > 0) {
					inter_tbl->neighbor_pe_export[n++] = comm_dst->ranks[i];
				}
			}
			HECMW_assert(n == inter_tbl->n_neighbor_pe_export);
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
set_export_index(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm)
{
	int *sendbuf_index = NULL, *recvbuf_index = NULL, *sendbuf = NULL, *recvbuf = NULL;
	int size, rtc, i;

	/*
	 * send buffer (destination unit)
	 */
	if(comm_dst->is_member) {
		/* index of send buffer */
		sendbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_import+1, sizeof(int));
		if(sendbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<inter_tbl->n_neighbor_pe_import; i++) {
			sendbuf_index[i+1] = sendbuf_index[i] + 1;
		}

		/* send buffer */
		sendbuf = (int *)HECMW_malloc(sizeof(int)*(inter_tbl->n_neighbor_pe_import+1));
		if(sendbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<inter_tbl->n_neighbor_pe_import; i++) {
			sendbuf[i] = inter_tbl->import_index[i+1] - inter_tbl->import_index[i];
		}
	}

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* index of receive buffer */
		recvbuf_index = (int *)HECMW_calloc(inter_tbl->n_neighbor_pe_export+1, sizeof(int));
		if(recvbuf_index == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<inter_tbl->n_neighbor_pe_export; i++) {
			recvbuf_index[i+1] = recvbuf_index[i] + 1;
		}

		/* receive buffer */
		recvbuf = (int *)HECMW_malloc(sizeof(int)*(inter_tbl->n_neighbor_pe_export+1));
		if(recvbuf == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	if(inter_tbl->n_neighbor_pe_import == 0 &&
			inter_tbl->n_neighbor_pe_export == 0) return HECMW_SUCCESS;

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(
			inter_tbl->n_neighbor_pe_import, inter_tbl->neighbor_pe_import, sendbuf_index, sendbuf,
			inter_tbl->n_neighbor_pe_export, inter_tbl->neighbor_pe_export, recvbuf_index, recvbuf,
			HECMW_INT, intercomm->comm);
	if(rtc != 0)  goto error;

	/*
	 * set index of export node for inter-communication (source unit)
	 */
	if(comm_src->is_member) {
		inter_tbl->export_index = NULL;
		if(inter_tbl->n_neighbor_pe_export != 0) {
			size = inter_tbl->n_neighbor_pe_export + 1;
			inter_tbl->export_index = (int *)HECMW_calloc(size, sizeof(int));
			if(inter_tbl->export_index == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			for(i=0; i<inter_tbl->n_neighbor_pe_export; i++) {
				inter_tbl->export_index[i+1] = inter_tbl->export_index[i] + recvbuf[i];
			}
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
set_export_item(struct hecmw_couple_inter_iftable *inter_tbl,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm, int *export_data)
{
	int size, rtc, i;

	if(inter_tbl->n_neighbor_pe_import == 0 &&
			inter_tbl->n_neighbor_pe_export == 0) return 0;

	/*
	 * receive buffer (source unit)
	 */
	if(comm_src->is_member) {
		/* receive buffer */
		size = sizeof(int)*(inter_tbl->export_index[inter_tbl->n_neighbor_pe_export]+1);
		inter_tbl->export_item = (int *)HECMW_malloc(size);
		if(inter_tbl->export_item == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
	}

	/*
	 * send and receive
	 */
	rtc = HECMW_couple_inter_send_recv(
			inter_tbl->n_neighbor_pe_import, inter_tbl->neighbor_pe_import,
			inter_tbl->import_index, export_data,
			inter_tbl->n_neighbor_pe_export, inter_tbl->neighbor_pe_export,
			inter_tbl->export_index, inter_tbl->export_item, HECMW_INT, intercomm->comm);
	if(rtc != 0)  goto error;

	return 0;

error:
	return -1;
}



static struct hecmw_couple_inter_iftable *
set_inter_iftable(const struct hecmw_couple_comm *comm_src,
		const struct hecmw_couple_comm *comm_dst, const struct hecmw_couple_comm *intercomm,
		struct map_info **map_dst, const struct hecmw_couple_mapped_point *mapped_point)
{
	struct hecmw_couple_inter_iftable *inter_tbl = NULL;
	struct link_list_map *mapping_data_list = NULL;
	int *export_data = NULL;
	int rtc;

	/* allocation */
	inter_tbl = HECMW_couple_alloc_inter_iftable();
	if(inter_tbl == NULL) goto error;

	/* import side */
	if(comm_dst->is_member) {
		mapping_data_list = set_mapping_surf(mapped_point, comm_src, map_dst);
		if(mapping_data_list == NULL) goto error;

		rtc = set_neighbor_info_import(inter_tbl, mapping_data_list, comm_src);
		if(rtc != HECMW_SUCCESS) goto error;

		export_data = set_export_data(inter_tbl, comm_src, mapping_data_list);
		if(export_data == NULL) goto error;
	}

	/* export side */
	rtc = set_neighbor_pe_export(inter_tbl, comm_src, comm_dst, intercomm, mapping_data_list);
	if(rtc != HECMW_SUCCESS) goto error;

	rtc = set_export_index(inter_tbl, comm_src, comm_dst, intercomm);
	if(rtc != HECMW_SUCCESS) goto error;

	rtc = set_export_item(inter_tbl, comm_src, comm_dst, intercomm, export_data);
	if(rtc != HECMW_SUCCESS) goto error;

	/* finalization */
	HECMW_free(export_data);
	return inter_tbl;

error:
	HECMW_free(export_data);
	return NULL;
}


/*================================================================================================*/

extern struct hecmw_couple_inter_iftable *
HECMW_couple_set_map_data(
		const struct hecmwST_local_mesh *mesh_src, const struct hecmwST_local_mesh *mesh_dst,
		const struct hecmw_couple_comm *comm_src, const struct hecmw_couple_comm *comm_dst,
		const struct hecmw_couple_comm *intercomm,
		const struct hecmw_couple_boundary *boundary_src,
		const struct hecmw_couple_bounding_box *bbox_src,
		const struct hecmw_couple_bounding_box *bbox_dst,
		const struct hecmw_couple_background_cell *bgcell_src,
		const struct hecmw_couple_mapped_point *mapped_point)
{
	struct hecmw_couple_inter_iftable *inter_tbl = NULL;
	struct hecmw_couple_mapped_point **all_mapped_point = NULL;
	struct map_info **map_src = NULL, **map_dst = NULL;
	int **is_candidate = NULL;
	int size, rtc, i;

	if(comm_src->is_member) {
		size = sizeof(struct hecmw_couple_mapped_point *) * comm_dst->psize;
		all_mapped_point = (struct hecmw_couple_mapped_point **)HECMW_malloc(size);
		if(all_mapped_point == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			all_mapped_point[i] = NULL;
		}

		is_candidate = (int **)HECMW_malloc(sizeof(int *)*comm_dst->psize);
		if(is_candidate == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			is_candidate[i] = NULL;
		}
	}

	rtc = bcast_mapped_point_d2s(mapped_point, comm_src, comm_dst, intercomm, all_mapped_point);
	if(rtc != HECMW_SUCCESS) goto error;

	rtc = set_candidate_node(bbox_src, bbox_dst, comm_src, comm_dst, intercomm,
			all_mapped_point, is_candidate);
	if(rtc != HECMW_SUCCESS) goto error;

	/* */
	if(comm_src->is_member) {
		map_src = (struct map_info **)HECMW_malloc(sizeof(struct map_info *)*comm_dst->psize);
		if(map_src == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_dst->psize; i++) {
			map_src[i] = NULL;
		}
		for(i=0; i<comm_dst->psize; i++) {
			map_src[i] = alloc_struct_map_info();
			if(map_src[i] == NULL) goto error;
		}
	}

	/* mapping on source unit */
	if(comm_src->is_member) {
		rtc = mapping_in_src(mesh_src, boundary_src, bbox_src, bgcell_src, comm_dst,
				all_mapped_point, is_candidate, map_src);
		if(rtc != HECMW_SUCCESS) goto error;
	}

	if(all_mapped_point) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_couple_free_mapped_point(all_mapped_point[i]);
		}
		HECMW_free(all_mapped_point);
	}
	if(is_candidate) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_free(is_candidate[i]);
		}
		HECMW_free(is_candidate);
	}

	/* mapping on destination unit */
	if(comm_dst->is_member) {
		map_dst = (struct map_info **)HECMW_malloc(sizeof(struct map_info *)*comm_dst->psize);
		if(map_dst == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		for(i=0; i<comm_src->psize; i++) {
			map_dst[i] = NULL;
		}
	}

	rtc = gather_map_data_s2d(mapped_point, map_src, comm_src, comm_dst, intercomm, map_dst);
	if(rtc != HECMW_SUCCESS) goto error;

	/* construct interface table for inter-communicate */
	inter_tbl = set_inter_iftable(comm_src, comm_dst, intercomm, map_dst, mapped_point);
	if(inter_tbl == NULL) goto error;

	return inter_tbl;

error:
	if(all_mapped_point) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_couple_free_mapped_point(all_mapped_point[i]);
		}
		HECMW_free(all_mapped_point);
	}
	if(is_candidate) {
		for(i=0; i<comm_dst->psize; i++) {
			HECMW_free(is_candidate[i]);
		}
		HECMW_free(is_candidate);
	}
	return NULL;
}

