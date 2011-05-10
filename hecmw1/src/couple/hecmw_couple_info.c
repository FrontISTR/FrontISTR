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

#include "hecmw_msgno.h"
#include "hecmw_struct.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"
#include "hecmw_comm.h"

#include "hecmw_couple_define.h"
#include "hecmw_couple_struct.h"
#include "hecmw_couple_control.h"
#include "hecmw_couple_info.h"

 
static struct intracomm_info {
	char unit_id[HECMW_NAME_LEN+1];		
	struct hecmw_couple_comm *comm;		
	struct intracomm_info *next;		
} intracomm_root = {
	"",		/* unit_id	*/
	NULL,	/* comm		*/
	NULL,	/* next		*/
};


static struct intercomm_info {
	char couple_id[HECMW_NAME_LEN+1];	
	char unit1_id[HECMW_NAME_LEN+1];	
	char unit2_id[HECMW_NAME_LEN+1];	
	struct hecmw_couple_comm *comm;		
	struct intercomm_info *next;		
} intercomm_root = {
	"",		/* couple_id	*/
	"",		/* unit1_id		*/
	"",		/* unit2_id		*/
	NULL,	/* comm			*/	
	NULL,	/* next			*/
};


static struct couple_info {
	char boundary_id[HECMW_NAME_LEN+1];		
	char couple_id[HECMW_NAME_LEN+1];		
	char unit1_id[HECMW_NAME_LEN+1];		
	char unit2_id[HECMW_NAME_LEN+1];		
	int couple_type;						
	int direction;							
	double tolerance;						
	double bbcoef;							
	double bgcoef;							
	struct hecmw_couple_group *unit1_grp;	
	struct hecmw_couple_group *unit2_grp;	
	struct couple_info *next;				
} couple_root = {
	"",								/* boundary_id		*/
	"",								/* couple_id		*/
	"",								/* unit1_id			*/
	"",								/* unit2_id			*/
	HECMW_COUPLE_TYPE_UNDEF,		/* couple_type		*/
	HECMW_COUPLE_DIRECTION_UNDEF,	/* direction		*/
	-1.0,							/* tolerance		*/
	-1.0,							/* bbcoef			*/
	-1.0,							/* bgcoef			*/
	NULL,							/* unit1_grp		*/
	NULL,							/* unit2_grp		*/
	NULL,							/* next				*/
};


static int is_initialized = 0;


/*================================================================================================*/

extern void
HECMW_couple_free_comm(struct hecmw_couple_comm *comm)
{
	if(comm == NULL) return;

	HECMW_free(comm->ranks);
	HECMW_free(comm);
	comm = NULL;
}


static struct hecmw_couple_comm *
alloc_struct_comm(void) {
	struct hecmw_couple_comm *comm;

	comm = (struct hecmw_couple_comm *)HECMW_malloc(sizeof(struct hecmw_couple_comm));
	if(comm == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}

	comm->psize     = 0;
	comm->rank      = -1;
	comm->ranks     = NULL;
	comm->comm      = (HECMW_Comm)(-1);
	comm->group     = (HECMW_Group)(-1);
	comm->root      = -1;
	comm->is_member = 0;
	comm->is_root   = 0;

	return comm;
}


static void
free_intracomm_info(void)
{
	struct intracomm_info *p, *q;

	for(p=intracomm_root.next; p; p=q) {
		q = p->next;
		HECMW_couple_free_comm(p->comm);
		HECMW_free(p);
	}
}


static void
free_intercomm_info(void)
{
	struct intercomm_info *p, *q;

	for(p=intercomm_root.next; p; p=q) {
		q = p->next;
		HECMW_couple_free_comm(p->comm);
		HECMW_free(p);
	}
}


static void
free_couple_info(void)
{
	struct couple_info *p, *q;

	for(p=couple_root.next; p; p=q) {
		q = p->next;
		HECMW_couple_ctrl_free_group(p->unit1_grp);
		HECMW_couple_ctrl_free_group(p->unit2_grp);
		HECMW_free(p);
	}
}


extern void
HECMW_couple_free_couple_info(void)
{
	free_intracomm_info();
	free_intercomm_info();
	free_couple_info();
}

/*------------------------------------------------------------------------------------------------*/

static struct intracomm_info *
get_intracomm_info(const char *unit_id)
{
	struct intracomm_info *p;

	for(p=intracomm_root.next; p; p=p->next) {
		if(strcmp(unit_id, p->unit_id) == 0) return p;
	}

	return NULL;
}


static struct intercomm_info *
get_intercomm_info(const char *couple_id)
{
	struct intercomm_info *p;

	for(p=intercomm_root.next; p; p=p->next) {
		if(strcmp(couple_id, p->couple_id) == 0) return p;
	}

	return NULL;
}


static struct couple_info *
get_couple_info(const char *boundary_id)
{
	struct couple_info *p;

	for(p=couple_root.next; p; p=p->next) {
		if(strcmp(boundary_id, p->boundary_id) == 0) return p;
	}

	return NULL;
}

/*------------------------------------------------------------------------------------------------*/

static int
init_intracomm_info(void)
{
	struct hecmw_couple_ctrl_unit_ids *unit_ids;
	struct hecmw_couple_ctrl_boundary_ids *boundary_ids;
	struct intracomm_info *p;
	int *mask = NULL;
	char *couple_id, *unit1_id, *unit2_id;
	int i, j;

	if((unit_ids = HECMW_couple_get_unit_ids()) == NULL) return HECMW_ERROR;
	if((boundary_ids = HECMW_couple_get_boundary_ids()) == NULL) return HECMW_ERROR;

	/* masking */
	mask = (int *)HECMW_malloc(sizeof(int)*unit_ids->n_unit);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<unit_ids->n_unit; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}
	for(i=0; i<boundary_ids->n_boundary; i++) {
		couple_id = HECMW_couple_ctrl_get_couple_id(boundary_ids->ids[i], NULL, 0);
		if(couple_id == NULL) goto error;
		unit1_id = HECMW_couple_ctrl_get_unit_id(couple_id, HECMW_COUPLE_UNIT1, NULL, 0);
		if(unit1_id == NULL) goto error;
		unit2_id = HECMW_couple_ctrl_get_unit_id(couple_id, HECMW_COUPLE_UNIT2, NULL, 0);
		if(unit2_id == NULL) goto error;

		for(j=0; j<unit_ids->n_unit; j++) {
			if(strcmp(unit1_id, unit_ids->ids[j]) == 0) {
				mask[j] = HECMW_COUPLE_TRUE;
			}
			if(strcmp(unit2_id, unit_ids->ids[j]) == 0) {
				mask[j] = HECMW_COUPLE_TRUE;
			}
		}
	}

	/* setting */
	p = &intracomm_root;
	for(i=0; i<unit_ids->n_unit; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) {
			p->next = (struct intracomm_info *)HECMW_malloc(sizeof(struct intracomm_info));
			if(p->next == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			p = p->next;

			strcpy(p->unit_id, unit_ids->ids[i]);
			p->comm = alloc_struct_comm();
			if(p->comm == NULL) goto error;

			p->next = NULL;
		}
	}

	HECMW_free(mask);
	return HECMW_SUCCESS;

error:
	HECMW_free(mask);
	free_intracomm_info();
	return HECMW_ERROR;
}


static int
init_intercomm_info(void)
{
	struct hecmw_couple_ctrl_couple_ids *couple_ids;
	struct hecmw_couple_ctrl_boundary_ids *boundary_ids;
	struct intercomm_info *p;
	int *mask = NULL;
	char *couple_id, unit1_id[HECMW_NAME_LEN+1], unit2_id[HECMW_NAME_LEN+1];
	int i, j;

	if((couple_ids = HECMW_couple_get_couple_ids()) == NULL) return HECMW_ERROR;
	if((boundary_ids = HECMW_couple_get_boundary_ids()) == NULL) return HECMW_ERROR;

	/* masking */
	mask = (int *)HECMW_malloc(sizeof(int)*couple_ids->n_couple);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	for(i=0; i<couple_ids->n_couple; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}
	for(i=0; i<boundary_ids->n_boundary; i++) {
		couple_id = HECMW_couple_ctrl_get_couple_id(boundary_ids->ids[i], NULL, 0);
		if(couple_id == NULL) goto error;

		for(j=0; j<couple_ids->n_couple; j++) {
			if(strcmp(couple_id, couple_ids->ids[j]) == 0) {
				mask[j] = HECMW_COUPLE_TRUE;
			}
		}
	}

	/* setting */
	p = &intercomm_root;
	for(i=0; i<couple_ids->n_couple; i++) {
		if(mask[i] == HECMW_COUPLE_TRUE) {
			p->next = (struct intercomm_info *)HECMW_malloc(sizeof(struct intercomm_info));
			if(p->next == NULL) {
				HECMW_set_error(errno, "");
				goto error;
			}
			p = p->next;

			if(HECMW_couple_ctrl_get_unit_id(couple_ids->ids[i], HECMW_COUPLE_UNIT1,
						unit1_id, HECMW_NAME_LEN+1) == NULL) goto error;
			if(HECMW_couple_ctrl_get_unit_id(couple_ids->ids[i], HECMW_COUPLE_UNIT2,
						unit2_id, HECMW_NAME_LEN+1) == NULL) goto error;

			strcpy(p->couple_id, couple_ids->ids[i]);
			strcpy(p->unit1_id, unit1_id);
			strcpy(p->unit2_id, unit2_id);
			if((p->comm = alloc_struct_comm()) == NULL) goto error;

			p->next = NULL;
		}
	}

	HECMW_free(mask);
	return HECMW_SUCCESS;

error:
	HECMW_free(mask);
	free_intercomm_info();
	return HECMW_ERROR;
}


static int
init_couple_info(void)
{
	struct hecmw_couple_ctrl_boundary_ids *boundary_ids;
	struct couple_info *p;
	char *boundary_id, *unit1_id, *unit2_id;
	int i;

	if((boundary_ids = HECMW_couple_get_boundary_ids()) == NULL) return HECMW_ERROR;

	p = &couple_root;
	for(i=0; i<boundary_ids->n_boundary; i++) {
		boundary_id = boundary_ids->ids[i];

		p->next = (struct couple_info *)HECMW_malloc(sizeof(struct couple_info));
		if(p->next == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		p = p->next;
		p->unit1_grp = NULL;
		p->unit2_grp = NULL;
		p->next = NULL;

		strcpy(p->boundary_id, boundary_id);
		if(HECMW_couple_ctrl_get_couple_id(boundary_id,
					p->couple_id, HECMW_NAME_LEN+1) == NULL) goto error;;
		if(HECMW_couple_ctrl_get_unit_id(p->couple_id, HECMW_COUPLE_UNIT1,
					p->unit1_id, HECMW_NAME_LEN+1) == NULL) goto error;
		if(HECMW_couple_ctrl_get_unit_id(p->couple_id, HECMW_COUPLE_UNIT2,
				p->unit2_id, HECMW_NAME_LEN+1) == NULL) goto error;

		if(HECMW_couple_ctrl_get_type(p->couple_id, &p->couple_type) != HECMW_SUCCESS) goto error;
		if(HECMW_couple_ctrl_get_direction(boundary_id, &p->direction) != HECMW_SUCCESS) goto error;
		if(HECMW_couple_ctrl_get_tolerance(boundary_id, &p->tolerance) != HECMW_SUCCESS) goto error;
		if(HECMW_couple_ctrl_get_bbcoef(boundary_id, &p->bbcoef) != HECMW_SUCCESS) goto error;
		if(HECMW_couple_ctrl_get_bgcoef(boundary_id, &p->bgcoef) != HECMW_SUCCESS) goto error;

		p->unit1_grp = HECMW_couple_ctrl_get_group(boundary_id, HECMW_COUPLE_UNIT1);
		if(p->unit1_grp == NULL) goto error;
		p->unit2_grp = HECMW_couple_ctrl_get_group(boundary_id, HECMW_COUPLE_UNIT2);
		if(p->unit2_grp == NULL) goto error;
	}

	return HECMW_SUCCESS;

error:
	free_couple_info();
	return HECMW_ERROR;
}

/*------------------------------------------------------------------------------------------------*/

static int
set_couple_type(int *couple_type)
{
	struct couple_info *p;
	int is_specified_mxn = 0;
	int is_specified_maxmn = 0;
	int is_specified_manual = 0;

	for(p=couple_root.next; p; p=p->next) {
		if(p->couple_type == HECMW_COUPLE_TYPE_MXN) {
			is_specified_mxn = 1;
		} else if(p->couple_type == HECMW_COUPLE_TYPE_MAXMN) {
			is_specified_maxmn = 1;
		} else if(p->couple_type == HECMW_COUPLE_TYPE_MANUAL) {
			is_specified_manual = 1;
		} else {
			HECMW_set_error(HECMWCPL_E_INVALID_CPLTYPE, "");
			return HECMW_ERROR;
		}
	}

	if(is_specified_mxn+is_specified_maxmn+is_specified_manual != 1) {
		HECMW_set_error(HECMWCPL_E_MULTIPLE_CPLTYPE, "");
		return HECMW_ERROR;
	}

	if(is_specified_mxn) {
		*couple_type = HECMW_COUPLE_TYPE_MXN;
	} else if(is_specified_maxmn) {
		*couple_type = HECMW_COUPLE_TYPE_MAXMN;
	} else if(is_specified_manual) {
		*couple_type = HECMW_COUPLE_TYPE_MANUAL;
	} else {
		HECMW_assert(0);
	}

	return HECMW_SUCCESS;
}


static int
check_intracomm_psize_mxn(void)
{
	struct hecmw_couple_ctrl_proc *proc;
	struct intracomm_info *p;
	int psize_sum = 0;

	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		psize_sum += proc->n_proc;

		HECMW_couple_ctrl_free_proc(proc);
	}

	if(psize_sum != HECMW_comm_get_size()) {
		HECMW_set_error(HECMWCPL_E_UNMATCH_PSIZE, "Total number of processes is %d", psize_sum);
		return HECMW_ERROR;
	}

	return HECMW_SUCCESS;
}


static int
check_intracomm_psize_maxmn(void)
{
	struct hecmw_couple_ctrl_proc *proc;
	struct intracomm_info *p;
	int psize_max = 0;

	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		if(proc->n_proc > psize_max) psize_max = proc->n_proc;

		HECMW_couple_ctrl_free_proc(proc);
	}

	if(psize_max != HECMW_comm_get_size()) {
		HECMW_set_error(HECMWCPL_E_UNMATCH_PSIZE, "Total number of processes is %d", psize_max);
		return HECMW_ERROR;
	}

	return HECMW_SUCCESS;
}


static int
check_intracomm_psize_manual(void)
{
	struct hecmw_couple_ctrl_proc *proc;
	struct intracomm_info *p;
	int *mask = NULL;
	int psize_sum = 0, psize_max = 0;
	int rank, n, i;

	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		psize_sum += proc->n_proc;

		HECMW_couple_ctrl_free_proc(proc);
	}

	/* masking */
	mask = (int *)HECMW_malloc(sizeof(int)*psize_sum);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		return HECMW_ERROR;
	}
	for(i=0; i<psize_sum; i++) {
		mask[i] = HECMW_COUPLE_FALSE;
	}
	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		for(i=0; i<proc->n_proc; i++) {
			rank = proc->ranks[i];
			if(rank < 0) {
				HECMW_set_error(HECMWCPL_E_INVALID_RANKS,
						"specified rank number is too small (%d)", rank);
				HECMW_couple_ctrl_free_proc(proc);
				goto error;
			}
			if(rank >= psize_sum) {
				HECMW_set_error(HECMWCPL_E_INVALID_RANKS,
						"specified rank number is too large (%d)", rank);
				HECMW_couple_ctrl_free_proc(proc);
				goto error;
			}
			mask[rank]++;
		}
		HECMW_couple_ctrl_free_proc(proc);
	}

	for(i=0; i<psize_sum; i++) {
		if(mask[i]) psize_max = i;
	}
	for(n=0, i=0; i<=psize_max; i++) {
		if(!mask[i]) {
			HECMW_set_error(HECMWCPL_E_DISCONTINUOUS_RANKS, "Process No. %d is not used");
			goto error;
		}
		n++;
	}
	if(n != HECMW_comm_get_size()) {
		HECMW_set_error(HECMWCPL_E_UNMATCH_PSIZE, "Total number of processes is %d", n);
	}

	HECMW_free(mask);
	return HECMW_SUCCESS;

error:
	HECMW_free(mask);
	return HECMW_ERROR;
}


static int
check_intracomm_psize(int couple_type)
{
	if(couple_type == HECMW_COUPLE_TYPE_MXN) {
		if(check_intracomm_psize_mxn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MAXMN) {
		if(check_intracomm_psize_maxmn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MANUAL) {
		if(check_intracomm_psize_manual() != HECMW_SUCCESS) return HECMW_ERROR;
	} else {
		HECMW_assert(0);
	}

	return HECMW_SUCCESS;
}


static int
set_intracomm_psize(void)
{
	struct hecmw_couple_ctrl_proc *proc;
	struct intracomm_info *p;

	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		HECMW_assert(p->comm);
		p->comm->psize = proc->n_proc;

		HECMW_couple_ctrl_free_proc(proc);
	}

	return HECMW_SUCCESS;
}


static int
set_intracomm_ranks_mxn(void)
{
	struct intracomm_info *p;
	int n = 0;
	int i;

	for(p=intracomm_root.next; p; p=p->next) {
		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			free_intracomm_info();
			return HECMW_ERROR;
		}
		for(i=0; i<p->comm->psize; i++) {
			p->comm->ranks[i] = n;
			n++;
		}
	}

	return HECMW_SUCCESS;
}


static int
set_intracomm_ranks_maxmn(void)
{
	struct intracomm_info *p;
	int i;

	for(p=intracomm_root.next; p; p=p->next) {
		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			free_intracomm_info();
			return HECMW_ERROR;
		}
		for(i=0; i<p->comm->psize; i++) {
			p->comm->ranks[i] = i;
		}
	}

	return HECMW_SUCCESS;
}


static int
set_intracomm_ranks_manual(void)
{
	struct hecmw_couple_ctrl_proc *proc;
	struct intracomm_info *p;
	int i;

	for(p=intracomm_root.next; p; p=p->next) {
		proc = HECMW_couple_ctrl_get_proc(p->unit_id);
		if(proc == NULL) return HECMW_ERROR;

		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			HECMW_couple_ctrl_free_proc(proc);
			free_intracomm_info();
			return HECMW_ERROR;
		}
		for(i=0; i<p->comm->psize; i++) {
			p->comm->ranks[i] = proc->ranks[i];
		}

		HECMW_couple_ctrl_free_proc(proc);
	}

	return HECMW_SUCCESS;
}


static int
set_intracomm_ranks(int couple_type)
{
	if(couple_type == HECMW_COUPLE_TYPE_MXN) {
		if(set_intracomm_ranks_mxn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MAXMN) {
		if(set_intracomm_ranks_maxmn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MANUAL) {
		if(set_intracomm_ranks_manual() != HECMW_SUCCESS) return HECMW_ERROR;
	} else {
		HECMW_assert(0);
	}

	return HECMW_SUCCESS;
}

/*------------------------------------------------------------------------------------------------*/

static int
set_intercomm_psize_mxn(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) return HECMW_ERROR;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) return HECMW_ERROR;

		p->comm->psize = q1->comm->psize + q2->comm->psize;
	}

	return HECMW_SUCCESS;
}


static int
set_intercomm_psize_maxmn(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) return HECMW_ERROR;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) return HECMW_ERROR;

		if(q1->comm->psize >= q2->comm->psize) {
			p->comm->psize = q1->comm->psize;
		} else {
			p->comm->psize = q2->comm->psize;
		}
	}

	return HECMW_SUCCESS;
}


static int
set_intercomm_psize_manual(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;
	int *mask = NULL;
	int global_psize, n, i;

	global_psize = HECMW_comm_get_size();

	mask = (int *)HECMW_malloc(sizeof(int)*global_psize);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		return HECMW_ERROR;
	}

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) goto error;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) goto error;

		for(i=0; i<global_psize; i++) {
			mask[i] = HECMW_COUPLE_FALSE;
		}
		for(i=0; i<q1->comm->psize; i++) {
			mask[q1->comm->ranks[i]] = HECMW_COUPLE_TRUE;
		}
		for(i=0; i<q2->comm->psize; i++) {
			mask[q2->comm->ranks[i]] = HECMW_COUPLE_TRUE;
		}
		for(n=0, i=0; i<global_psize; i++) {
			if(mask[i] == HECMW_COUPLE_TRUE) n++;
		}

		p->comm->psize = n;
	}

	HECMW_free(mask);
	return HECMW_SUCCESS;

error:
	HECMW_free(mask);
	return HECMW_ERROR;
}


static int
set_intercomm_psize(int couple_type)
{
	if(couple_type == HECMW_COUPLE_TYPE_MXN) {
		if(set_intercomm_psize_mxn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MAXMN) {
		if(set_intercomm_psize_maxmn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MANUAL) {
		if(set_intercomm_psize_manual() != HECMW_SUCCESS) return HECMW_ERROR;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_CPLTYPE, "");
		return HECMW_ERROR;
	}

	return HECMW_SUCCESS;
}

/*------------------------------------------------------------------------------------------------*/

static int
set_intercomm_ranks_mxn(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;
	int n, i;

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) goto error;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) goto error;

		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		HECMW_assert(p->comm->psize == q1->comm->psize + q2->comm->psize);
		for(n=0, i=0; i<q1->comm->psize; i++, n++) {
			p->comm->ranks[n] = q1->comm->ranks[i];
		}
		for(i=0; i<q2->comm->psize; i++, n++) {
			p->comm->ranks[n] = q2->comm->ranks[i];
		}
	}

	return HECMW_SUCCESS;

error:
	free_couple_info();
	return HECMW_ERROR;
}


static int
set_intercomm_ranks_maxmn(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;
	int i;

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) goto error;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) goto error;

		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}
		if(q1->comm->psize >= q2->comm->psize) {
			HECMW_assert(p->comm->psize == q1->comm->psize);
			for(i=0; i<q1->comm->psize; i++) {
				p->comm->ranks[i] = q1->comm->ranks[i];
			}
		} else {
			HECMW_assert(p->comm->psize == q2->comm->psize);
			for(i=0; i<q2->comm->psize; i++) {
				p->comm->ranks[i] = q2->comm->ranks[i];
			}
		}
	}

	return HECMW_SUCCESS;

error:
	free_couple_info();
	return HECMW_ERROR;
}


static int
set_intercomm_ranks_manual(void)
{
	struct intercomm_info *p;
	struct intracomm_info *q1, *q2;
	int *mask = NULL;
	int global_psize, n, i;

	global_psize = HECMW_comm_get_size();

	mask = (int *)HECMW_malloc(sizeof(int)*global_psize);
	if(mask == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	for(p=intercomm_root.next; p; p=p->next) {
		q1 = get_intracomm_info(p->unit1_id);
		if(q1 == NULL) goto error;
		q2 = get_intracomm_info(p->unit2_id);
		if(q2 == NULL) goto error;

		p->comm->ranks = (int *)HECMW_calloc(p->comm->psize, sizeof(int));
		if(p->comm->ranks == NULL) {
			HECMW_set_error(errno, "");
			goto error;
		}

		for(i=0; i<global_psize; i++) {
			mask[i] = HECMW_COUPLE_FALSE;
		}
		for(i=0; i<q1->comm->psize; i++) {
			mask[q1->comm->ranks[i]] = HECMW_COUPLE_TRUE;
		}
		for(i=0; i<q2->comm->psize; i++) {
			mask[q2->comm->ranks[i]] = HECMW_COUPLE_TRUE;
		}
		for(n=0, i=0; i<global_psize; i++) {
			if(mask[i] == HECMW_COUPLE_TRUE) {
				HECMW_assert(n <= p->comm->psize);
				p->comm->ranks[n] = i;
				n++;
			}
		}
	}

	HECMW_free(mask);
	return HECMW_SUCCESS;

error:
	HECMW_free(mask);
	free_couple_info();
	return HECMW_ERROR;
}


static int
set_intercomm_ranks(int couple_type)
{
	if(couple_type == HECMW_COUPLE_TYPE_MXN) {
		if(set_intercomm_ranks_mxn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MAXMN) {
		if(set_intercomm_ranks_maxmn() != HECMW_SUCCESS) return HECMW_ERROR;
	} else if(couple_type == HECMW_COUPLE_TYPE_MANUAL) {
		if(set_intercomm_ranks_manual() != HECMW_SUCCESS) return HECMW_ERROR;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_CPLTYPE, "");
		return HECMW_ERROR;
	}

	return HECMW_SUCCESS;
}

/*================================================================================================*/

static int
set_is_member(struct hecmw_couple_comm *comm)
{
	int global_rank, i;

	global_rank = HECMW_comm_get_rank();
	if(global_rank < 0) return HECMW_ERROR;

	comm->is_member = 0;
	for(i=0; i<comm->psize; i++) {
		if(comm->ranks[i] == global_rank) {
			comm->is_member = 1;
			break;
		}
	}

	return HECMW_SUCCESS;
}


static int
allgather_root(int *root)
{
	int *send_buf = NULL, *recv_buf = NULL;
	HECMW_Comm global_comm;
	int rtc, i;

	global_comm = HECMW_comm_get_comm();

	send_buf = (int *)HECMW_calloc(1, sizeof(int));
	if(send_buf == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}
	recv_buf = (int *)HECMW_calloc(HECMW_comm_get_size(), sizeof(int));
	if(recv_buf == NULL) {
		HECMW_set_error(errno, "");
		goto error;
	}

	send_buf[0] = *root;

	rtc = HECMW_Allgather(send_buf, 1, HECMW_INT, recv_buf, 1, HECMW_INT, global_comm);
	if(rtc != HECMW_SUCCESS) goto error;

	for(i=0; i<HECMW_comm_get_size(); i++) {
		if(recv_buf[i] >= 0) {
			*root = i;
		}
	}

	HECMW_free(send_buf);
	HECMW_free(recv_buf);
	return HECMW_SUCCESS;

error:
	HECMW_free(send_buf);
	HECMW_free(recv_buf);
	return HECMW_ERROR;
}


static int
set_root(struct hecmw_couple_comm *comm)
{
	if(comm->is_member == 1) {
		if(comm->rank == 0) {
			comm->root    = HECMW_comm_get_rank();
			comm->is_root = 1;
		}
	}

	if(allgather_root(&comm->root) != HECMW_SUCCESS) return HECMW_ERROR;

	return HECMW_SUCCESS;
}


static int
init_comm(struct hecmw_couple_comm *comm)
{
	int rtc, _psize;

	rtc = HECMW_Group_incl(HECMW_comm_get_group(), comm->psize, comm->ranks, &comm->group);
	if(rtc != HECMW_SUCCESS) return HECMW_ERROR;

	rtc = HECMW_Comm_create(HECMW_comm_get_comm(), comm->group, &comm->comm);
	if(rtc != HECMW_SUCCESS) return HECMW_ERROR;

	rtc = HECMW_Group_rank(comm->group, &comm->rank);
	if(rtc != HECMW_SUCCESS) return HECMW_ERROR;

	rtc = HECMW_Group_size(comm->group, &_psize);
	if(rtc != HECMW_SUCCESS) return HECMW_ERROR;
	HECMW_assert(_psize == comm->psize);

	if(set_is_member(comm) != HECMW_SUCCESS) return HECMW_ERROR;
	if(set_root(comm) != HECMW_SUCCESS) return HECMW_ERROR;

	return HECMW_SUCCESS;
}

/*================================================================================================*/

extern int
HECMW_couple_comm_init(void)
{
	struct intracomm_info *p;
	struct intercomm_info *q;
	int couple_type;

	if(HECMW_couple_ctrl_get_n_boundary() == 0) return HECMW_SUCCESS;

	if(init_intracomm_info() != HECMW_SUCCESS) goto error;
	if(init_intercomm_info() != HECMW_SUCCESS) goto error;
	if(init_couple_info() != HECMW_SUCCESS) goto error;

	if(set_couple_type(&couple_type) != HECMW_SUCCESS) goto error;
	if(check_intracomm_psize(couple_type) != HECMW_SUCCESS) goto error;

	/* set intra-communication info. */
	if(set_intracomm_psize() != HECMW_SUCCESS) goto error;
	if(set_intracomm_ranks(couple_type) != HECMW_SUCCESS) goto error;

	for(p=intracomm_root.next; p; p=p->next) {
		if(init_comm(p->comm) != HECMW_SUCCESS) goto error;
	}

	/* set inter-communication info. */
	if(set_intercomm_psize(couple_type) != HECMW_SUCCESS) goto error;
	if(set_intercomm_ranks(couple_type) != HECMW_SUCCESS) goto error;

	for(q=intercomm_root.next; q; q=q->next) {
		if(init_comm(q->comm) != HECMW_SUCCESS) goto error;
	}

	return HECMW_SUCCESS;

error:
	HECMW_couple_free_couple_info();
	return HECMW_ERROR;
}

/*================================================================================================*/

extern char *
HECMW_couple_get_unit_id(const char *boundary_id, int unit_specifier, char *buf, int bufsize)
{
	struct couple_info *couple;
	struct hecmw_couple_info *couple_info;
	char *unit_id, *ret_buf;
	int len;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return NULL;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		unit_id = couple->unit1_id;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		unit_id = couple->unit2_id;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return NULL;
	}

	if(buf == NULL) {
		ret_buf = HECMW_strdup(unit_id);
		if(ret_buf == NULL) {
			HECMW_set_error(errno, "");
			return NULL;
		}
	} else {
		len = strlen(unit_id);
		if(bufsize <= len) {
			len = bufsize - 1;
		}
		strncpy(buf, unit_id, len);
		buf[len] = '\0';
		ret_buf  = buf;
	}

	return ret_buf;
}


extern int
HECMW_couple_is_member(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error( HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	if(intercomm->comm->is_member) return 1;
	return 0;
}


extern int
HECMW_couple_is_unit_member(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	if(intracomm->comm->is_member) return 1;
	return 0;
}


extern int
HECMW_couple_is_unit_member_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	if(intracomm->comm->is_member) return 1;
	return 0;
}


extern int
HECMW_couple_is_root(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	if(intercomm->comm->is_root) return 1;
	return 0;
}


extern int
HECMW_couple_is_unit_root(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	if(intracomm->comm->is_root) return 1;
	return 0;
}


extern int
HECMW_couple_is_unit_root_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	if(intracomm->comm->is_root) return 1;
	return 0;
}


extern int
HECMW_intercomm_get_size(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	return intercomm->comm->psize;
}


extern int
HECMW_intracomm_get_size(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	return intracomm->comm->psize;
}


extern int
HECMW_intracomm_get_size_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	return intracomm->comm->psize;
}


extern int
HECMW_intercomm_get_rank(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	return intercomm->comm->rank;
}


extern int
HECMW_intracomm_get_rank(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	return intracomm->comm->rank;
}


extern int
HECMW_intracomm_get_rank_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	return intracomm->comm->rank;
}


extern HECMW_Comm
HECMW_intercomm_get_comm(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	return intercomm->comm->comm;
}


extern HECMW_Comm
HECMW_intracomm_get_comm(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	return intracomm->comm->comm;
}


extern HECMW_Comm
HECMW_intracomm_get_comm_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	return intracomm->comm->comm;
}


extern HECMW_Group
HECMW_intercomm_get_group(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;
	if((intercomm = get_intercomm_info(couple->couple_id)) == NULL) return -1;

	return intercomm->comm->group;
}


extern HECMW_Group
HECMW_intracomm_get_group(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple;
	struct intracomm_info *intracomm;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	if((couple = get_couple_info(boundary_id)) == NULL) return -1;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		if((intracomm = get_intracomm_info(couple->unit1_id)) == NULL) return -1;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		if((intracomm = get_intracomm_info(couple->unit2_id)) == NULL) return -1;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return -1;
	}

	return intracomm->comm->group;
}


extern HECMW_Group
HECMW_intracomm_get_group_u(const char *unit_id)
{
	struct intracomm_info *intracomm;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return -1;
	}

	if((intracomm = get_intracomm_info(unit_id)) == NULL) return -1;

	return intracomm->comm->group;
}


extern struct hecmw_couple_comm *
HECMW_couple_get_intracomm(const char *boundary_id, int unit_specifier)
{
	struct couple_info *couple = NULL;
	struct intracomm_info *intracomm = NULL;
	struct hecmw_couple_comm *comm;
	char *unit_id;
	int i;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}

	couple = get_couple_info(boundary_id);
	if(couple == NULL) return NULL;

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		unit_id = couple->unit1_id;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		unit_id = couple->unit2_id;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return NULL;
	}

	intracomm = get_intracomm_info(unit_id);
	if(intracomm == NULL) return NULL;

	comm = alloc_struct_comm();
	if(comm == NULL) return NULL;

	comm->psize = intracomm->comm->psize;
	comm->rank  = intracomm->comm->rank;
	comm->ranks = (int *)HECMW_calloc(comm->psize, sizeof(int));
	if(comm->ranks == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_comm(comm);
		return NULL;
	}
	for(i=0; i<comm->psize; i++) {
		comm->ranks[i] = intracomm->comm->ranks[i];
	}
	comm->comm      = intracomm->comm->comm;
	comm->group     = intracomm->comm->group;
	comm->root      = intracomm->comm->root;
	comm->is_member = intracomm->comm->is_member;
	comm->is_root   = intracomm->comm->is_root;

	return comm;
}


extern struct hecmw_couple_comm *
HECMW_couple_get_intracomm_u(const char *unit_id)
{
	struct intracomm_info *intracomm = NULL;
	struct hecmw_couple_comm *comm;
	int i;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}

	intracomm = get_intracomm_info(unit_id);
	if(intracomm == NULL) return NULL;

	comm = alloc_struct_comm();
	if(comm == NULL) return NULL;

	comm->psize = intracomm->comm->psize;
	comm->rank  = intracomm->comm->rank;
	comm->ranks = (int *)HECMW_calloc(comm->psize, sizeof(int));
	if(comm->ranks == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_comm(comm);
		return NULL;
	}
	for(i=0; i<comm->psize; i++) {
		comm->ranks[i] = intracomm->comm->ranks[i];
	}
	comm->comm      = intracomm->comm->comm;
	comm->group     = intracomm->comm->group;
	comm->root      = intracomm->comm->root;
	comm->is_member = intracomm->comm->is_member;
	comm->is_root   = intracomm->comm->is_root;

	return comm;
}


extern struct hecmw_couple_comm *
HECMW_couple_get_intercomm(const char *boundary_id)
{
	struct couple_info *couple;
	struct intercomm_info *intercomm;
	struct hecmw_couple_comm *comm;
	int i;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}

	couple = get_couple_info(boundary_id);
	if(couple == NULL) return NULL;

	intercomm = get_intercomm_info(couple->couple_id);
	if(intercomm == NULL) return NULL;

	comm = alloc_struct_comm();
	if(comm == NULL) return NULL;

	comm->psize = intercomm->comm->psize;
	comm->rank  = intercomm->comm->rank;
	comm->ranks = (int *)HECMW_calloc(comm->psize, sizeof(int));
	if(comm->ranks == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_comm(comm);
		return NULL;
	}
	for(i=0; i<comm->psize; i++) {
		comm->ranks[i] = intercomm->comm->ranks[i];
	}
	comm->comm      = intercomm->comm->comm;
	comm->group     = intercomm->comm->group;
	comm->root      = intercomm->comm->root;
	comm->is_member = intercomm->comm->is_member;
	comm->is_root   = intercomm->comm->is_root;

	return comm;
}

