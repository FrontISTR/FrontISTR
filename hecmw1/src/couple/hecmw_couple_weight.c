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

#include "hecmw_malloc.h"
#include "hecmw_couple_define.h"
#include "hecmw_couple_weight.h"


extern struct hecmw_couple_weight *
HECMW_couple_alloc_weight(void)
{
	struct hecmw_couple_weight *p = NULL;

	p = (struct hecmw_couple_weight *)HECMW_malloc(sizeof(struct hecmw_couple_weight));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->n      = 0;
	p->type   = HECMW_COUPLE_IP_UNDEF;
	p->index  = NULL;
	p->id     = NULL;
	p->weight = NULL;

	return p;
}



extern void
HECMW_couple_free_weight(struct hecmw_couple_weight *p)
{
	if(p == NULL) return;

	HECMW_free(p->index);
	HECMW_free(p->id);
	HECMW_free(p->weight);
}


/*------------------------------------------------------------------------------------------------*/

extern struct hecmw_couple_weight_list *
HECMW_couple_alloc_weight_list(void)
{
	struct hecmw_couple_weight_list *p;

	p = (struct hecmw_couple_weight_list *)HECMW_malloc(sizeof(struct hecmw_couple_weight_list));
	if(p == NULL) {
		HECMW_set_error(errno, "");
		return NULL;
	}
	p->info = NULL;
	p->next = NULL;

	return p;
}



extern void
HECMW_couple_free_weight_list(struct hecmw_couple_weight_list *r)
{
	struct hecmw_couple_weight_list *p, *q;

	if(r == NULL) return;

	for(p=r; p; p=q) {
		q = p->next;
		HECMW_couple_free_weight(p->info);
		HECMW_free(p);
	}
	r = NULL;
}

