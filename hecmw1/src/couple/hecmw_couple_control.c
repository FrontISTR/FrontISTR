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

#include "hecmw_config.h"
#include "hecmw_msgno.h"
#include "hecmw_malloc.h"
#include "hecmw_error.h"
#include "hecmw_ctrllex.h"
#include "hecmw_control.h"
#include "hecmw_util.h"

#include "hecmw_couple_define.h"
/* #include "hecmw_couple_struct.h" */
#include "hecmw_couple_control.h"


/*================================================================================================*/

struct link_list_i {
	int rank;						
	struct link_list_i *next;		
};



struct link_list_s {
	char name[HECMW_NAME_LEN+1];	
	struct link_list_s *next;		
};



static struct unit_info_by_ctrl {
	char unit_id[HECMW_NAME_LEN+1];		
	int n_proc;							
	int is_specified_ranks;				
	struct link_list_i ranks;			
	struct unit_info_by_ctrl *next;		
} unit_info_root = {
	"",					/* unit_id				*/
	0,					/* n_proc				*/
	HECMW_COUPLE_FALSE,	/* is_specified_ranks	*/
	{
		-1,				/* ranks.rank			*/
		NULL,			/* ranks.next			*/
	},
	NULL,				/* next					*/
};



static struct couple_info_by_ctrl {
	char couple_id[HECMW_NAME_LEN+1];	
	char unit1_id[HECMW_NAME_LEN+1];	
	char unit2_id[HECMW_NAME_LEN+1];	
	int couple_type;					
	struct couple_info_by_ctrl *next;	
} couple_info_root = {
	"",							/* couple_id	*/
	"",							/* unit1_id		*/
	"",							/* unit2_id		*/
	HECMW_COUPLE_TYPE_UNDEF,	/* couple_type	*/
	NULL,						/* next			*/
};



struct group_info_by_ctrl {
	int n_grp;						
	int geom_type;					
	int data_type;					
	struct link_list_s grp_name;	
};



static struct boundary_info_by_ctrl {
	char boundary_id[HECMW_NAME_LEN+1];		
	char couple_id[HECMW_NAME_LEN+1];		
	int interpolation;						
	int direction;							
	double tolerance;						
	double bbcoef;							
	double bgcoef;							
	struct group_info_by_ctrl unit1_grp;	
	struct group_info_by_ctrl unit2_grp;	
	struct boundary_info_by_ctrl *next;		
} boundary_info_root = {
	"",								/* boundary_id				*/
	"",								/* couple_id				*/
	HECMW_COUPLE_IP_UNDEF,			/* interpolation			*/
	HECMW_COUPLE_DIRECTION_UNDEF,	/* direction				*/
	HECMW_COUPLE_TOLERANCE_DEFAULT,	/* tolrance					*/
	HECMW_COUPLE_BBCOEF_DEFAULT,	/* bb_coef					*/
	HECMW_COUPLE_BGCOEF_DEFAULT,	/* bg_coef					*/
	{
		0,							/* unit1_grp.n_grp			*/
		HECMW_COUPLE_GROUP_UNDEF,	/* unit1_grp.geom_type		*/
		HECMW_COUPLE_GROUP_UNDEF,	/* unit1_grp.data_type		*/
		{
			"",						/* unit1_grp.grp_name.name	*/
			NULL					/* unit1_grp.grp_name.next	*/
		}
	},
	{
		0,							/* unit2_grp.n_grp			*/
		HECMW_COUPLE_GROUP_UNDEF,	/* unit2_grp.geom_type		*/
		HECMW_COUPLE_GROUP_UNDEF,	/* unit2_grp.data_type		*/
		{
			"",						/* unit2_grp.grp_name.name	*/
			NULL					/* unit2_grp.grp_name.next	*/
		}
	},
	NULL,							/* next						*/
};



static struct unit_info_by_ctrl *unit_info_current = &unit_info_root;



static struct couple_info_by_ctrl *couple_info_current = &couple_info_root;



static struct boundary_info_by_ctrl *boundary_info_current = &boundary_info_root;



static int n_unit = 0;



static int n_couple = 0;



static int n_boundary = 0;


/*================================================================================================*/

static void
free_link_list_s(struct link_list_s *r)
{
	struct link_list_s *p, *q;

	for(p=r->next; p; p=q) {
		q = p->next;
		HECMW_free(p);
	}
}



static void
free_link_list_i(struct link_list_i *r)
{
	struct link_list_i *p, *q;

	for(p=r->next; p; p=q) {
		q = p->next;
		HECMW_free(p);
	}
}



extern void
HECMW_couple_ctrl_print_unit(FILE *fp)
{
	struct unit_info_by_ctrl *p;
	struct link_list_i *q;

	for(p=unit_info_root.next; p; p=p->next) {
		fprintf(fp, "Unit ID: %s\n", p->unit_id);
		fprintf(fp, "Number of Processes: %d\n", p->n_proc);
		fprintf(fp, "Process rank is specified or not: %d\n", p->is_specified_ranks);
		for(q=p->ranks.next; q; q=q->next) {
			fprintf(fp, "Process Number: %d\n", q->rank);
		}
	}
}


extern void
HECMW_couple_ctrl_print_couple(FILE *fp)
{
	struct couple_info_by_ctrl *p;

	for(p=couple_info_root.next; p; p=p->next) {
		fprintf(fp, "Couple ID: %s\n", p->couple_id);
		fprintf(fp, "Unit1 ID: %s\n", p->unit1_id);
		fprintf(fp, "Unit2 ID: %s\n", p->unit2_id);
		fprintf(fp, "Couple Type: %d\n", p->couple_type);
	}
}


extern void
HECMW_couple_ctrl_print_boundary(FILE *fp)
{
	struct boundary_info_by_ctrl *p;
	struct link_list_s *q;

	for(p=boundary_info_root.next; p; p=p->next) {
		fprintf(fp, "Boundary ID: %s\n", p->boundary_id);
		fprintf(fp, "Couple ID: %s\n", p->couple_id);
		fprintf(fp, "Direction: %d\n", p->direction);
		fprintf(fp, "Interpolation: %d\n", p->interpolation);
		fprintf(fp, "Tolerance: %e\n", p->tolerance);
		fprintf(fp, "BBcoef: %e\n", p->bbcoef);
		fprintf(fp, "BGcoef: %e\n", p->bgcoef);
		fprintf(fp, "Unit1 Group: %d\n", p->unit1_grp.n_grp);
		fprintf(fp, "Unit1 Group Type (Geometry): %d\n", p->unit1_grp.geom_type);
		fprintf(fp, "Unit1 Group Type (Data): %d\n", p->unit1_grp.data_type);
		for(q=p->unit1_grp.grp_name.next; q; q=q->next) {
			fprintf(fp, "Unit1 Group Name: %s\n", q->name);
		}
		fprintf(fp, "Unit2 Group: %d\n", p->unit2_grp.n_grp);
		fprintf(fp, "Unit2 Group Type (Geometry): %d\n", p->unit2_grp.geom_type);
		fprintf(fp, "Unit2 Group Type (Data): %d\n", p->unit2_grp.data_type);
		for(q=p->unit2_grp.grp_name.next; q; q=q->next) {
			fprintf(fp, "Unit2 Group Name: %s\n", q->name);
		}
	}
}

/*------------------------------------------------------------------------------------------------*/

static struct couple_info_by_ctrl *
get_couple_by_id(const char *couple_id)
{
	struct couple_info_by_ctrl *p;

	if(couple_id == NULL) return NULL;

	for(p=couple_info_root.next; p; p=p->next) {
		if((strcmp(p->couple_id, couple_id)) == 0) return p;
	}

	return NULL;
}



static struct unit_info_by_ctrl *
get_unit_by_id(const char *unit_id)
{
	struct unit_info_by_ctrl *p;

	if(unit_id == NULL) return NULL;

	for(p=unit_info_root.next; p; p=p->next) {
		if((strcmp(p->unit_id, unit_id)) == 0) return p;
	}

	return NULL;
}



static struct boundary_info_by_ctrl *
get_boundary_by_id(const char *boundary_id)
{
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) return NULL;

	for(p=boundary_info_root.next; p; p=p->next) {
		if((strcmp(p->boundary_id, boundary_id)) == 0) return p;
	}

	return NULL;
}


/*================================================================================================*/

extern void
HECMW_couple_free_unit_ids(struct hecmw_couple_ctrl_unit_ids *unit_ids)
{
	int i;

	if(unit_ids == NULL) return;

	if(unit_ids->ids) {
		for(i=0; i<unit_ids->n_unit; i++) {
			HECMW_free(unit_ids->ids[i]);
		}
		HECMW_free(unit_ids->ids);
	}
	HECMW_free(unit_ids);
	unit_ids = NULL;
}



extern void
HECMW_couple_free_couple_ids(struct hecmw_couple_ctrl_couple_ids *couple_ids)
{
	int i;

	if(couple_ids == NULL) return;

	if(couple_ids->ids) {
		for(i=0; i<couple_ids->n_couple; i++) {
			HECMW_free(couple_ids->ids[i]);
		}
		HECMW_free(couple_ids->ids);
	}
	HECMW_free(couple_ids);
	couple_ids = NULL;
}



extern void
HECMW_couple_free_boundary_ids(struct hecmw_couple_ctrl_boundary_ids *boundary_ids)
{
	int i;

	if(boundary_ids == NULL) return;

	if(boundary_ids->ids) {
		for(i=0; i<boundary_ids->n_boundary; i++) {
			HECMW_free(boundary_ids->ids[i]);
		}
		HECMW_free(boundary_ids);
	}
	HECMW_free(boundary_ids);
	boundary_ids = NULL;

}



extern void
HECMW_couple_ctrl_free_proc(struct hecmw_couple_ctrl_proc *proc_info)
{
	if(proc_info == NULL) return;

	HECMW_free(proc_info->ranks);
	HECMW_free(proc_info);
	proc_info = NULL;
}



extern void
HECMW_couple_ctrl_free_group(struct hecmw_couple_group *grp_info)
{
	int i;

	if(grp_info == NULL) return;

	for(i=0; i<grp_info->n_grp; i++) {
		HECMW_free(grp_info->grp_name[i]);
	}
	HECMW_free(grp_info);
	grp_info == NULL;
}



static void
free_unit(void)
{
	struct unit_info_by_ctrl *p, *q;

	for(p=unit_info_root.next; p; p=q) {
		q = p->next;
		free_link_list_i(&p->ranks);
		HECMW_free(p);
	}
	unit_info_root.next = NULL;
}



static void
free_couple(void)
{
	struct couple_info_by_ctrl *p, *q;

	for(p=couple_info_root.next; p; p=q) {
		q = p->next;
		HECMW_free(p);
	}
	couple_info_root.next = NULL;
}



static void
free_boundary(void)
{
	struct boundary_info_by_ctrl *p, *q;

	for(p=boundary_info_root.next; p; p=q) {
		q = p->next;
		free_link_list_s(&p->unit1_grp.grp_name);
		free_link_list_s(&p->unit2_grp.grp_name);
		HECMW_free(p);
	}
	boundary_info_root.next = NULL;
}



extern void
HECMW_couple_ctrl_free(void)
{
	free_unit();
	free_couple();
	free_boundary();
}


/*================================================================================================*/

static int
get_unit_name(char *unit_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'NAME' in '!COUPLE UNIT'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of coupling unit */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'NAME' is specified in '!COUPLE UNIT'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'NAME' is too long in '!COUPLE UNIT'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	strcpy(unit_id, s);

	return 0;
}



static int
get_unit_nproc(int *n_proc)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'NPROC' in '!COUPLE UNIT'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Number of processes */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_INT) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'NPROC' is specified in '!COUPLE UNIT'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	*n_proc = (int)HECMW_ctrllex_get_number();
	if(n_proc < 0) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'NPROC' in '!COUPLE UNIT' must be natural number (%d)",
				HECMW_ctrllex_get_lineno(), *n_proc);
		return -1;
	}

	return 0;
}



static int
get_unit_ranks(struct link_list_i *p, int *counter)
{
	int token, rank;

	while(1) {
		/* Process number */
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != HECMW_CTRLLEX_INT) {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Invalid process number is found in '!COUPLE UNIT'",
					HECMW_ctrllex_get_lineno());
			return -1;
		}
		rank = (int)HECMW_ctrllex_get_number();
		if(rank < 0) {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Process number in '!COUPLE UNIT' must be natural number (%d)",
					HECMW_ctrllex_get_lineno(), rank);
			return -1;
		}

		p->next = (struct link_list_i *)HECMW_malloc(sizeof(struct link_list_i));
		if(p->next == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p       = p->next;
		p->rank = rank;
		p->next = NULL;
		(*counter)++;

		/* , */
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != ',') {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Process number must be delimited by ',' in '!COUPLE UNIT'",
					HECMW_ctrllex_get_lineno());
			return -1;
		}
	}

	return 0;
}



static int
get_unit_1st_line(struct unit_info_by_ctrl *p)
{
	int token;
	int is_specified_name = 0;
	int is_specified_nproc = 0;

	while(1) {
		if((token = HECMW_ctrllex_next_token()) != ',') {
			if(token == HECMW_CTRLLEX_NL) break;
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}

		token = HECMW_ctrllex_next_token();

		/* NAME=<name-of-coupling-unit> */
		if(token == HECMW_CTRLLEX_K_NAME) {
			if(is_specified_name) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'NAME' is re-defined in !COUPLE UNIT",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_unit_name(p->unit_id)) return -1;
			is_specified_name = 1;

		/* NPROC=<number-of-processes> */
		} else if(token == HECMW_CTRLLEX_K_NPROC) {
			if(is_specified_nproc) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'NPROC' is re-defined in !COUPLE UNIT",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_unit_nproc(&p->n_proc)) return -1;
			is_specified_nproc = 1;

		/* New Line */
		} else if(token == HECMW_CTRLLEX_NL) {
			break;

		/* Invalid Token */
		} else {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}
	}

	/* NAME is not specified */
	if(is_specified_name == 0) {
		HECMW_set_error(HECMWCPL_E_CPLU_NO_NAME, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* NPROC is not specified */
	if(is_specified_nproc <= 0) {
		HECMW_set_error(HECMWCPL_E_CPLU_NO_NPROC, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}

	return 0;
}



static int
get_unit_2nd_line(struct unit_info_by_ctrl *p)
{
	int token, counter;
	struct link_list_i *q;

	/* proc-1, proc-2, ..., proc-n */
	q = &p->ranks;
	counter = 0;
	while((token = HECMW_ctrllex_next_token())) {
		if(token == HECMW_CTRLLEX_NL) continue;

		HECMW_ctrllex_unput_token();

		if(token != HECMW_CTRLLEX_INT) break;

		if(get_unit_ranks(q, &counter)) return -1;
	}

	/* check number of processes */
	if(counter > 0) {
		if(counter != p->n_proc) {
			HECMW_set_error(HECMWCPL_E_CPLU_UNMATCH_RANKS, "");
			return -1;
		}
		p->is_specified_ranks = 1;
	}

	return 0;
}



extern int
HECMW_couple_ctrl_unit(void)
{
	int token, size;
	struct unit_info_by_ctrl *p;

	token = HECMW_ctrllex_next_token();
	HECMW_assert(token == HECMW_CTRLLEX_H_COUPLE_UNIT);

	/* allocation & initialization */
	size = sizeof(struct unit_info_by_ctrl);
	unit_info_current->next = (struct unit_info_by_ctrl *)HECMW_malloc(size);
	if(unit_info_current->next == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	p = unit_info_current->next;
	memset(p->unit_id, 0, HECMW_NAME_LEN+1);
	p->n_proc             = 0;
	p->is_specified_ranks = HECMW_COUPLE_FALSE;
	p->ranks.rank         = -1;
	p->ranks.next         = NULL;
	p->next               = NULL;

	/* get 1st line */
	if(get_unit_1st_line(p)) return -1;

	/* skip blank line */
	if((token = HECMW_ctrllex_next_token()) == HECMW_CTRLLEX_NL);
	HECMW_ctrllex_unput_token();

	/* get 2nd line */
	if(get_unit_2nd_line(p)) return -1;

	unit_info_current = p;
	n_unit++;

	return 0;
}


/*================================================================================================*/

static int
get_couple_name(char *couple_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'NAME' in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of couple */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'NAME' is specified in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'NAME' in '!COUPLE' is too long",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	strcpy(couple_id, s);

	return 0;
}



static int
get_couple_type(int *type)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'TYPE' in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Coupling type */
	token = HECMW_ctrllex_next_token();
	if(token == HECMW_CTRLLEX_K_MXN) {				/* MXN */
		*type = HECMW_COUPLE_TYPE_MXN;
	} else if(token == HECMW_CTRLLEX_K_MAXMN) {		/* MAXMN */
		*type = HECMW_COUPLE_TYPE_MAXMN;
	} else if(token == HECMW_CTRLLEX_K_MANUAL) {	/* MANUAL */
		*type = HECMW_COUPLE_TYPE_MANUAL;
	} else {										/* Invalid Token */
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'TYPE' is specified in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	return 0;
}



static int
get_couple_unit1(char *unit1_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'UNIT1' in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of coupling unit 1 */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'UNIT1' is specified in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'UNIT1' in '!COUPLE' is too long",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	if(get_unit_by_id(s) == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_UNIT_ID, "%s", s);
		return -1;
	}
	strcpy(unit1_id, s);

	return 0;
}



static int
get_couple_unit2(char *unit2_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'UNIT2' in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of coupling unit 2 */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'UNIT2' is specified in '!COUPLE'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'UNIT2' in '!COUPLE' is too long",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	if(get_unit_by_id(s) == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_UNIT_ID, "%s", s);
		return -1;
	}
	strcpy(unit2_id, s);

	return 0;
}



static int
get_couple_1st_line(struct couple_info_by_ctrl *p)
{
	int token;
	int is_specified_name  = 0;
	int is_specified_type  = 0;
	int is_specified_unit1 = 0;
	int is_specified_unit2 = 0;

	while(1) {
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != ',') {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}

		token = HECMW_ctrllex_next_token();
		/* NAME=<name-of-couple> */
		if(token == HECMW_CTRLLEX_K_NAME) {
			if(is_specified_name) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'NAME' is re-defined in '!COUPLE'",
								 HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_couple_name(p->couple_id)) return -1;
			is_specified_name = 1;

		/* TYPE=<couple-type> */
		} else if(token == HECMW_CTRLLEX_K_TYPE) {
			if(is_specified_type) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'TYPE' is re-defined in '!COUPLE'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_couple_type(&p->couple_type)) return -1;
			is_specified_type = 1;

		/* UNIT1=<name-of-unit1> */
		} else if(token == HECMW_CTRLLEX_K_UNIT1) {
			if(is_specified_unit1) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'UNIT1' is re-defined in '!COUPLE'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_couple_unit1(p->unit1_id)) return -1;
			is_specified_unit1 = 1;

		/* UNIT2=<name-of-unit2> */
		} else if(token == HECMW_CTRLLEX_K_UNIT2) {
			if(is_specified_unit2) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'UNIT2' is re-defined in '!COUPLE'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_couple_unit2(p->unit2_id)) return -1;
			is_specified_unit2 = 1;

		/* New Line */
		} else if(token == HECMW_CTRLLEX_NL) {
			break;

		/* Invalid Token */
		} else {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}
	}

	/* NAME is not specified */
	if(!is_specified_name) {
		HECMW_set_error(HECMWCPL_E_CPL_NO_NAME, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* TYPE is not specified */
	if(!is_specified_type) {
		HECMW_set_error(HECMWCPL_E_CPL_NO_TYPE, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* UNIT1 is not specified */
	if(!is_specified_unit1) {
		HECMW_set_error(HECMWCPL_E_CPL_NO_UNIT1, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* UNIT2 is not specified */
	if(!is_specified_unit2) {
		HECMW_set_error(HECMWCPL_E_CPL_NO_UNIT2, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}

	return 0;
}



extern int
HECMW_couple_ctrl_couple(void)
{
	int token, size;
	struct couple_info_by_ctrl *p;

	token = HECMW_ctrllex_next_token();
	HECMW_assert(token == HECMW_CTRLLEX_H_COUPLE);

	/* allocation & initialization */
	size = sizeof(struct couple_info_by_ctrl);
	couple_info_current->next = (struct couple_info_by_ctrl *)HECMW_malloc(size);
	if(couple_info_current->next == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	p = couple_info_current->next;
	memset(p->couple_id, 0, HECMW_NAME_LEN+1);
	memset(p->unit1_id, 0, HECMW_NAME_LEN+1);
	memset(p->unit2_id, 0, HECMW_NAME_LEN+1);
	p->couple_type = HECMW_COUPLE_TYPE_UNDEF;
	p->next        = NULL;

	/* get 1st line */
	if(get_couple_1st_line(p)) return -1;

	couple_info_current = p;
	n_couple++;

	return 0;
}


/*------------------------------------------------------------------------------------------------*/

static int
get_boundary_name(char *boundary_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'NAME' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of coupling boundary */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'NAME' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'NAME' in '!COUPLE BOUNDARY' is too long",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	strcpy(boundary_id, s);

	return 0;
}



static int
get_boundary_couple(char *couple_id)
{
	int token;
	char *s;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'COUPLE' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Name of coupling */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_NAME) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'COUPLE' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	s = HECMW_ctrllex_get_text();
	if(strlen(s) > HECMW_NAME_LEN) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'COUPLE' is too long in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	if(get_couple_by_id(s) == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_COUPLE_ID, "%s", s);
		return -1;
	}
	strcpy(couple_id, s);

	return 0;
}



static int
get_boundary_direction(int *direction)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'DIRECTION' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Coupling direction */
	token = HECMW_ctrllex_next_token();
	if(token == HECMW_CTRLLEX_K_UNIT1_TO_UNIT2) {			/* from UNIT1 to UNIT2 */
		*direction = HECMW_COUPLE_UNIT1_TO_UNIT2;
	} else if(token == HECMW_CTRLLEX_K_UNIT2_TO_UNIT1) {	/* from UNIT2 to UNIT1 */
		*direction = HECMW_COUPLE_UNIT2_TO_UNIT1;
	} else {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'DIRECTION' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	return 0;
}



static int
get_boundary_interpolation(int *interpolation)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'INTERPOLATION' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Interpolating method */
	token = HECMW_ctrllex_next_token();

	/*@@@@@ NOT mounting @@@@@*/

	return 0;
}



static int
get_boundary_tolerance(double *tolerance)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'TOLERANCE' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Tolerance */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_DOUBLE) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'TOLERANCE' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	*tolerance = (double)HECMW_ctrllex_get_number();
	if(*tolerance < 0.0) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'TOLERANCE' in '!COUPLE BOUNDARY' must be grater than or equal 0 (%e)",
				HECMW_ctrllex_get_lineno(), *tolerance);
		return -1;
	}

	return 0;
}



static int
get_boundary_bbcoef(double *bbcoef)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'BBCOEF' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Enlarging coefficient for bounding box */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_DOUBLE) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'BBCOEF' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	*bbcoef = (double)HECMW_ctrllex_get_number();
	if(*bbcoef <= 0.0) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'BBCOEF' in '!COUPLE BOUNDARY' must be grater than 0 (%e)",
				HECMW_ctrllex_get_lineno(), *bbcoef);
		return -1;
	}

	return 0;
}



static int
get_boundary_bgcoef(double *bgcoef)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'BGCOEF' in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Enlarging coefficient for background cell */
	if((token = HECMW_ctrllex_next_token()) != HECMW_CTRLLEX_DOUBLE) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'BGCOEF' is specified in '!COUPLE BOUNDARY'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}
	*bgcoef = (double)HECMW_ctrllex_get_number();
	if(*bgcoef <= 0.0) {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: 'BGCOEF' in '!COUPLE BOUNDARY' must be grater than 0 (%e)",
				HECMW_ctrllex_get_lineno(), *bgcoef);
		return -1;
	}

	return 0;
}



static int
get_boundary_1st_line(struct boundary_info_by_ctrl *p)
{
	int token;
	int is_specified_name = 0;
	int is_specified_couple = 0;
	int is_specified_direction = 0;
	int is_specified_interpolation = 0;
	int is_specified_tolrance = 0;
	int is_specified_bbcoef = 0;
	int is_specified_bgcoef = 0;

	while(1) {
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != ',') {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}

		token = HECMW_ctrllex_next_token();
		/* NAME=<name-of-coupling-boundary> */
		if(token == HECMW_CTRLLEX_K_NAME) {
			if(is_specified_name) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'NAME' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_name(p->boundary_id) != HECMW_SUCCESS) return -1;
			is_specified_name = 1;

		/* COUPLE=<name-of-couple> */
		} else if(token == HECMW_CTRLLEX_K_COUPLE) {
			if(is_specified_couple) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'COUPLE' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_couple(p->couple_id) != HECMW_SUCCESS) return -1;
			is_specified_couple = 1;

		/* DIRECTION=<coupling-direction> */
		} else if(token == HECMW_CTRLLEX_K_DIRECTION) {
			if(is_specified_direction) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'DIRECTION' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_direction(&p->direction) != HECMW_SUCCESS) return -1;
			is_specified_direction = 1;

		/* INTERPOLATION=<interpolating-method> */
		} else if(token == HECMW_CTRLLEX_K_INTERPOLATION) {
			if(is_specified_interpolation) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'INTERPOLATION' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_interpolation(&p->interpolation) != HECMW_SUCCESS) return -1;
			is_specified_interpolation = 1;

		/* TOLERANCE=<tolerance> */
		} else if(token == HECMW_CTRLLEX_K_TOLERANCE) {
			if(is_specified_tolrance) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'TOLRANCE' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_tolerance(&p->tolerance) != HECMW_SUCCESS) return -1;
			is_specified_tolrance = 1;


		/* BBCOEF=<bounding-box-enlargement-coefficient> */
		} else if(token == HECMW_CTRLLEX_K_BBCOEF) {
			if(is_specified_bbcoef) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'BBCOEF' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_bbcoef(&p->bbcoef) != HECMW_SUCCESS) return -1;
			is_specified_bbcoef = 1;

		/* BGCOEF=<background-cell-enlargement-coefficient> */
		} else if(token == HECMW_CTRLLEX_K_BGCOEF) {
			if(is_specified_bgcoef) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"line %d: 'BGCOEF' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_bgcoef(&p->bgcoef) != HECMW_SUCCESS) return -1;
			is_specified_bgcoef = 1;

		/* New Line */
		} else if(token == HECMW_CTRLLEX_NL) {
			break;

		/* Invalid Token */
		} else {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}
	}

	/* NAME is not specified */
	if(!is_specified_name) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_NAME, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* COUPLE is not specified */
	if(!is_specified_couple) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_COUPLE, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	/* DIRECTION is not specified */
	if(!is_specified_direction) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_DIRECTION, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}

	return 0;
}


/*------------------------------------------------------------------------------------------------*/

static int
get_boundary_geom(int *geom_type)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'GEOM' in '!UNIT1' or '!UNIT2'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Boundary group type */
	token = HECMW_ctrllex_next_token();
	if(token == HECMW_CTRLLEX_K_NODE) {				/* NODE GROUP */
		*geom_type = HECMW_COUPLE_NODE_GROUP;
	} else if(token == HECMW_CTRLLEX_K_ELEMENT) {	/* ELEMENT GROUP */
		*geom_type = HECMW_COUPLE_ELEMENT_GROUP;
	} else if(token == HECMW_CTRLLEX_K_SURFACE) {	/* SURFACE GROUP */
		*geom_type = HECMW_COUPLE_SURFACE_GROUP;
	} else {										/* Invalid Token */
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'GEOM' is specified in '!UNIT1' or '!UNIT2'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	return 0;
}



static int
get_boundary_data(int *data_type)
{
	int token;

	/* = */
	if((token = HECMW_ctrllex_next_token()) != '=') {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: '=' is not found after 'DATA' in '!UNIT1' or '!UNIT2'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	/* Boundary group type */
	token = HECMW_ctrllex_next_token();
	if(token == HECMW_CTRLLEX_K_NODE) {				/* NODE GROUP */
		*data_type = HECMW_COUPLE_NODE_GROUP;
	} else if(token == HECMW_CTRLLEX_K_ELEMENT) {	/* ELEMENT GROUP */
		*data_type = HECMW_COUPLE_ELEMENT_GROUP;
	} else if(token == HECMW_CTRLLEX_K_SURFACE) {	/* SURFACE GROUP */
		*data_type = HECMW_COUPLE_SURFACE_GROUP;
	} else {										/* Invalid Token */
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
				"line %d: Invalid 'DATA' is specified in '!UNIT1' or '!UNIT2'",
				HECMW_ctrllex_get_lineno());
		return -1;
	}

	return 0;
}



static int
get_boundary_group_inner(struct link_list_s *p, int *counter)
{
	char *s;
	int token;

	while(1) {
		/* Group Name */
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != HECMW_CTRLLEX_NAME) {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Invalid group name is found in '!COUPLE BOUNDARY'",
					HECMW_ctrllex_get_lineno());
			return -1;
		}
		s = HECMW_ctrllex_get_text();
		if(strlen(s) > HECMW_NAME_LEN) {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Group name in '!COUPLE BOUNDARY' is too long",
					HECMW_ctrllex_get_lineno());
			return -1;
		}

		p->next = (struct link_list_s *)HECMW_malloc(sizeof(struct link_list_s));
		if(p->next == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		p = p->next;
		strcpy(p->name, s);
		p->next = NULL;
		(*counter)++;

		/* , */
		token = HECMW_ctrllex_next_token();
		if(token == HECMW_CTRLLEX_NL) break;
		if(token != ',') {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
					"line %d: Group name must be delimited by ',' in '!COUPLE BOUNDARY'",
					HECMW_ctrllex_get_lineno());
			return -1;
		}
	}

	return 0;
}



static int
get_boundary_unit_2nd_line(struct group_info_by_ctrl *p)
{
	int token, counter;
	struct link_list_s *q;

	q = &p->grp_name;
	counter = 0;
	while((token = HECMW_ctrllex_next_token())) {
		if(token == HECMW_CTRLLEX_NL) continue;

		HECMW_ctrllex_unput_token();

		if(token != HECMW_CTRLLEX_NAME) break;

		if(get_boundary_group_inner(q, &counter) != HECMW_SUCCESS) return -1;
	}

	p->n_grp = counter;

	/* Group name is not specified */
	if(counter <= 0) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_GRPNAME, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}

	return 0;
}



static int
get_boundary_unit_1st_line(struct group_info_by_ctrl *p)
{
	int token;
	int is_specified_geom = 0;
	int is_specified_data = 0;

	while(1) {
		if((token=HECMW_ctrllex_next_token()) != ',') {
			if(token == HECMW_CTRLLEX_NL) break;
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}

		token = HECMW_ctrllex_next_token();

		/* GEOM=<group-type-for-geometry> */
		if(token == HECMW_CTRLLEX_K_GEOM) {
			if(is_specified_geom) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"'GEOM' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_geom(&p->geom_type)) return -1;
			is_specified_geom = 1;

		/* DATA=<group-type-for-data> */
		} else if(token == HECMW_CTRLLEX_K_DATA) {
			if(is_specified_data) {
				HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN,
						"'DATA' is re-defined in '!COUPLE BOUNDARY'",
						HECMW_ctrllex_get_lineno());
				return -1;
			}
			if(get_boundary_data(&p->data_type)) return -1;
			is_specified_data = 1;

		/* New Line */
		} else if(token == HECMW_CTRLLEX_NL) {
			break;

		/* Invalid Token */
		} else {
			HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
					HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
			return -1;
		}
	}

	if(!is_specified_geom) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_GEOM, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}
	if(!is_specified_data) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_DATA, "line %d", HECMW_ctrllex_get_lineno()-1);
		return -1;
	}

	if(p->geom_type == HECMW_COUPLE_NODE_GROUP) {
		if(p->data_type != HECMW_COUPLE_NODE_GROUP) {
			HECMW_set_error(HECMWCPL_E_CPLB_UNMATCH_GRPTYPE,
					"line %d: NODE is specified in 'GEOM', but 'DATA' is not NODE",
					HECMW_ctrllex_get_lineno()-1);
			return -1;
		}
	} else if(p->geom_type == HECMW_COUPLE_ELEMENT_GROUP) {
		if(p->data_type != HECMW_COUPLE_ELEMENT_GROUP && p->data_type != HECMW_COUPLE_NODE_GROUP) {
			HECMW_set_error(HECMWCPL_E_CPLB_UNMATCH_GRPTYPE,
					"line %d: ELEMENT is specified in 'GEOM', but 'DATA' is not NODE or ELEMENT",
					HECMW_ctrllex_get_lineno()-1);
			return -1;
		}
	} else if(p->geom_type == HECMW_COUPLE_SURFACE_GROUP) {
		if(p->data_type != HECMW_COUPLE_SURFACE_GROUP && p->data_type != HECMW_COUPLE_NODE_GROUP) {
			HECMW_set_error(HECMWCPL_E_CPLB_UNMATCH_GRPTYPE,
					"line %d: SURFACE is specified in 'GEOM', but 'DATA' is not NODE or SURFACE",
					HECMW_ctrllex_get_lineno()-1);
			return -1;
		}
	} else {
		HECMW_assert(0);
	}

	return 0;
}



static int
get_boundary_2nd_line(struct boundary_info_by_ctrl *p)
{
	int token;
	int is_specified_unit1 = 0;
	int is_specified_unit2 = 0;

	while((token = HECMW_ctrllex_next_token()) == HECMW_CTRLLEX_NL);

	/* !UNIT1 */
	if(token == HECMW_CTRLLEX_H_UNIT1) {
		if(get_boundary_unit_1st_line(&p->unit1_grp)) return -1;
		if(get_boundary_unit_2nd_line(&p->unit1_grp)) return -1;
		is_specified_unit1 = 1;

		while((token = HECMW_ctrllex_next_token()) == HECMW_CTRLLEX_NL);

		if(token == HECMW_CTRLLEX_H_UNIT2) {
			if(get_boundary_unit_1st_line(&p->unit2_grp)) return -1;
			if(get_boundary_unit_2nd_line(&p->unit2_grp)) return -1;
			is_specified_unit2 = 1;
		}

	/* !UNIT2 */
	} else if(token == HECMW_CTRLLEX_H_UNIT2) {
		if(get_boundary_unit_1st_line(&p->unit2_grp)) return -1;
		if(get_boundary_unit_2nd_line(&p->unit2_grp)) return -1;
		is_specified_unit2 = 1;

		while((token = HECMW_ctrllex_next_token()) == HECMW_CTRLLEX_NL);

		if(token == HECMW_CTRLLEX_H_UNIT1) {
			if(get_boundary_unit_1st_line(&p->unit1_grp)) return -1;
			if(get_boundary_unit_1st_line(&p->unit2_grp)) return -1;
			is_specified_unit1 = 1;
		}

	/* Invalid Token */
	} else {
		HECMW_set_error(HECMWCPL_E_CTRL_INVALID_TOKEN, "line %d: %s",
						 HECMW_ctrllex_get_lineno(), HECMW_ctrllex_get_text());
		return -1;
	}

	if(!is_specified_unit1) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_UNIT1, "line %d", HECMW_ctrllex_get_lineno());
		return -1;
	}
	if(!is_specified_unit2) {
		HECMW_set_error(HECMWCPL_E_CPLB_NO_UNIT2, "line %d", HECMW_ctrllex_get_lineno());
		return -1;
	}

	return 0;
}



extern int
HECMW_couple_ctrl_boundary(void)
{
	int token, size;
	struct boundary_info_by_ctrl *p;

	token = HECMW_ctrllex_next_token();
	HECMW_assert(token == HECMW_CTRLLEX_H_COUPLE_BOUNDARY);

	/* allocation & initialization */
	size = sizeof(struct boundary_info_by_ctrl);
	boundary_info_current->next = (struct boundary_info_by_ctrl *)HECMW_malloc(size);
	if(boundary_info_current->next == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	p = boundary_info_current->next;

	memset(p->boundary_id, 0, HECMW_NAME_LEN+1);
	memset(p->couple_id, 0, HECMW_NAME_LEN+1);
	p->interpolation = HECMW_COUPLE_IP_UNDEF;
	p->direction = HECMW_COUPLE_DIRECTION_UNDEF;
	p->tolerance = HECMW_COUPLE_TOLERANCE_DEFAULT;
	p->bbcoef = HECMW_COUPLE_BBCOEF_DEFAULT;
	p->bgcoef = HECMW_COUPLE_BGCOEF_DEFAULT;
	p->unit1_grp.n_grp = 0;
	p->unit1_grp.geom_type = HECMW_COUPLE_GROUP_UNDEF;
	p->unit1_grp.data_type = HECMW_COUPLE_GROUP_UNDEF;
	memset(p->unit1_grp.grp_name.name, 0, HECMW_NAME_LEN+1);
	p->unit1_grp.grp_name.next = NULL;
	p->unit2_grp.n_grp = 0;
	p->unit2_grp.geom_type = HECMW_COUPLE_GROUP_UNDEF;
	p->unit2_grp.data_type = HECMW_COUPLE_GROUP_UNDEF;
	memset(p->unit2_grp.grp_name.name, 0, HECMW_NAME_LEN+1);
	p->unit2_grp.grp_name.next = NULL;
	p->next                    = NULL;

	/* Line 1 */
	if(get_boundary_1st_line(p)) return -1;

	if((token = HECMW_ctrllex_next_token()) == HECMW_CTRLLEX_NL);
	HECMW_ctrllex_unput_token();

	/* Line 2 */
	if(get_boundary_2nd_line(p)) return -1;

	boundary_info_current = p;
	n_boundary++;

	return 0;
}


/*================================================================================================*/

extern int
HECMW_couple_ctrl_get_n_unit(void)
{
	return n_unit;
}



extern int
HECMW_couple_ctrl_get_n_couple(void)
{
	return n_couple;
}



extern int
HECMW_couple_ctrl_get_n_boundary(void)
{
	return n_boundary;
}



extern struct hecmw_couple_ctrl_unit_ids *
HECMW_couple_get_unit_ids(void)
{
	int size, i, n;
	struct unit_info_by_ctrl *p;
	struct hecmw_couple_ctrl_unit_ids *unit_ids = NULL;

	size = sizeof(struct hecmw_couple_ctrl_unit_ids);
	unit_ids = (struct hecmw_couple_ctrl_unit_ids *)HECMW_malloc(size);
	if(unit_ids == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	unit_ids->n_unit = n_unit;
	unit_ids->ids    = NULL;

	if(n_unit == 0) return unit_ids;

	unit_ids->ids = (char **)HECMW_malloc(sizeof(char *)*n_unit);
	if(unit_ids->ids == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_unit_ids(unit_ids);
		return NULL;
	}
	for(i=0; i<n_unit; i++) {
		unit_ids->ids[i] = NULL;
	}
	for(n=0, p=unit_info_root.next; p; p=p->next, n++) {
		unit_ids->ids[n] = HECMW_strdup(p->unit_id);
		if(unit_ids->ids[n] == NULL) {
			HECMW_set_error(errno, "");
			HECMW_couple_free_unit_ids(unit_ids);
			return NULL;
		}
	}

	return unit_ids;
}



extern struct hecmw_couple_ctrl_couple_ids *
HECMW_couple_get_couple_ids(void)
{
	int size, i, n;
	struct couple_info_by_ctrl *p;
	struct hecmw_couple_ctrl_couple_ids *couple_ids = NULL;

	size = sizeof(struct hecmw_couple_ctrl_couple_ids);
	couple_ids = (struct hecmw_couple_ctrl_couple_ids *)HECMW_malloc(size);
	if(couple_ids == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	couple_ids->n_couple = n_couple;
	couple_ids->ids      = NULL;

	if(n_couple == 0) return couple_ids;

	couple_ids->ids = (char **)HECMW_malloc(sizeof(char *)*n_couple);
	if(couple_ids->ids == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_couple_ids(couple_ids);
		return NULL;
	}
	for(i=0; i<n_couple; i++) {
		couple_ids->ids[i] = NULL;
	}
	for(n=0, p=couple_info_root.next; p; p=p->next, n++) {
		couple_ids->ids[n] = HECMW_strdup(p->couple_id);
		if(couple_ids->ids[n] == NULL) {
			HECMW_set_error(errno, "");
			HECMW_couple_free_couple_ids(couple_ids);
			return NULL;
		}
	}

	return couple_ids;
}



extern struct hecmw_couple_ctrl_boundary_ids *
HECMW_couple_get_boundary_ids(void)
{
	int size, i, n;
	struct boundary_info_by_ctrl *p;
	struct hecmw_couple_ctrl_boundary_ids *boundary_ids = NULL;

	size = sizeof(struct hecmw_couple_ctrl_boundary_ids);
	boundary_ids = (struct hecmw_couple_ctrl_boundary_ids *)HECMW_malloc(size);
	if(boundary_ids == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	boundary_ids->n_boundary = n_boundary;
	boundary_ids->ids        = NULL;

	if(n_boundary == 0) return boundary_ids;

	boundary_ids->ids = (char **)HECMW_malloc(sizeof(char *)*n_boundary);
	if(boundary_ids->ids == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_free_boundary_ids(boundary_ids);
		return NULL;
	}
	for(i=0; i<n_boundary; i++) {
		boundary_ids->ids[i] = NULL;
	}
	for(n=0, p=boundary_info_root.next; p; p=p->next, n++) {
		boundary_ids->ids[n] = HECMW_strdup(p->boundary_id);
		if(boundary_ids->ids[n] == NULL) {
			HECMW_set_error(errno, "");
			HECMW_couple_free_boundary_ids(boundary_ids);
			return NULL;
		}
	}

	return boundary_ids;
}



extern char *
HECMW_couple_ctrl_get_unit_id(const char *couple_id, int unit_specifier, char *buf, int bufsize)
{
	int len;
	char *retbuf, *unit_id;
	struct couple_info_by_ctrl *p;

	if(couple_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (couple_id)");
		return NULL;
	}

	p = get_couple_by_id(couple_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_COUPLE_ID, "%s", couple_id);
		return NULL;
	}

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		unit_id = p->unit1_id;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		unit_id = p->unit2_id;
	} else {
		HECMW_assert(0);
	}

	if(buf == NULL) {
		retbuf = HECMW_strdup(unit_id);
		if(retbuf == NULL) {
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
		retbuf = buf;
	}

	return retbuf;
}



extern char *
HECMW_couple_ctrl_get_couple_id(const char *boundary_id, char *buf, int bufsize)
{
	int len;
	char *retbuf;
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return NULL;
	}

	if(buf == NULL) {
		retbuf = HECMW_strdup(p->couple_id);
		if(retbuf == NULL) {
			HECMW_set_error(errno, "");
			return NULL;
		}
	} else {
		len = strlen(p->couple_id);
		if(bufsize <= len) {
			len = bufsize - 1;
		}
		strncpy(buf, p->couple_id, len);
		buf[len] = '\0';
		retbuf = buf;
	}

	return retbuf;
}



extern struct hecmw_couple_ctrl_proc *
HECMW_couple_ctrl_get_proc(const char *unit_id)
{
	int size, n;
	struct hecmw_couple_ctrl_proc *proc_info;
	struct unit_info_by_ctrl *p;
	struct link_list_i *q;

	if(unit_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (unit_id)");
		return NULL;
	}

	p = get_unit_by_id(unit_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_UNIT_ID, "%s", unit_id);
		return NULL;
	}

	size = sizeof(struct hecmw_couple_ctrl_proc);
	proc_info = (struct hecmw_couple_ctrl_proc *)HECMW_malloc(size);
	if(proc_info == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	proc_info->n_proc             = p->n_proc;
	proc_info->is_specified_ranks = p->is_specified_ranks;
	proc_info->ranks              = NULL;

	if(!proc_info->is_specified_ranks) return proc_info;

	proc_info->ranks = (int *)HECMW_calloc(proc_info->n_proc, sizeof(int));
	if(proc_info->ranks == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_ctrl_free_proc(proc_info);
		return NULL;
	}
	for(n=0, q=p->ranks.next; q; q=q->next, n++) {
		proc_info->ranks[n] = q->rank;
	}

	return proc_info;
}



extern int
HECMW_couple_ctrl_get_type(const char *couple_id, int *couple_type)
{
	struct couple_info_by_ctrl *p;

	if(couple_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "");
		return -1;
	}

	p = get_couple_by_id(couple_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_COUPLE_ID, "%s", couple_id);
		return -1;
	}
	*couple_type = p->couple_type;

	return 0;
}



extern int
HECMW_couple_ctrl_get_direction(const char *boundary_id, int *direction)
{
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return -1;
	}
	*direction = p->direction;

	return 0;
}



extern int
HECMW_couple_ctrl_get_tolerance(const char *boundary_id, double *tolerance)
{
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return -1;
	}
	*tolerance = p->tolerance;

	return 0;
}



extern int
HECMW_couple_ctrl_get_bbcoef(const char *boundary_id, double *bbcoef)
{
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return -1;
	}
	*bbcoef = p->bbcoef;

	return 0;
}



extern int
HECMW_couple_ctrl_get_bgcoef(const char *boundary_id, double *bgcoef)
{
	struct boundary_info_by_ctrl *p;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return -1;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return -1;
	}
	*bgcoef = p->bgcoef;

	return 0;
}



extern struct hecmw_couple_group *
HECMW_couple_ctrl_get_group(const char *boundary_id, int unit_specifier)
{
	int size, i, n;
	struct hecmw_couple_group *grp_info = NULL;
	struct group_info_by_ctrl *grp_info_ctrl = NULL;
	struct boundary_info_by_ctrl *p;
	struct link_list_s *q;

	if(boundary_id == NULL) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG, "Invalid NULL pointer is found (boundary_id)");
		return NULL;
	}
	if((unit_specifier != HECMW_COUPLE_UNIT1) && (unit_specifier != HECMW_COUPLE_UNIT2)) {
		HECMW_set_error(HECMWCPL_E_INVALID_ARG,
				"Unrecognized unit specifier is found (unit_specifier)");
		return NULL;
	}

	p = get_boundary_by_id(boundary_id);
	if(p == NULL) {
		HECMW_set_error(HECMWCPL_E_UNDEF_BOUNDARY_ID, "%s", boundary_id);
		return NULL;
	}

	if(unit_specifier == HECMW_COUPLE_UNIT1) {
		grp_info_ctrl = &p->unit1_grp;
	} else if(unit_specifier == HECMW_COUPLE_UNIT2) {
		grp_info_ctrl = &p->unit2_grp;
	} else {
		HECMW_set_error(HECMWCPL_E_INVALID_UNITTYPE, "");
		return NULL;
	}

	size = sizeof(struct hecmw_couple_group);
	grp_info = (struct hecmw_couple_group *)HECMW_malloc(size);
	if(grp_info == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	grp_info->n_grp     = grp_info_ctrl->n_grp;
	grp_info->geom_type = grp_info_ctrl->geom_type;
	grp_info->data_type = grp_info_ctrl->data_type;

	grp_info->grp_name = (char **)HECMW_malloc(sizeof(char *)*grp_info->n_grp);
	if(grp_info->grp_name == NULL) {
		HECMW_set_error(errno, "");
		HECMW_couple_ctrl_free_group(grp_info);
		return NULL;
	}
	for(i=0; i<grp_info->n_grp; i++) {
		grp_info->grp_name[i] == NULL;
	}
	for(n=0, q=grp_info_ctrl->grp_name.next; q; q=q->next, n++) {
		grp_info->grp_name[n] = HECMW_strdup(q->name);
		if(grp_info->grp_name[n] == NULL) {
			HECMW_set_error(errno, "");
			HECMW_couple_ctrl_free_group(grp_info);
			return NULL;
		}
	}

	return grp_info;
}

