/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2006/06/01                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
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
#include "hecmw_util.h"
#include "hecmw_common.h"
#include "hecmw_gflex.h"
#include "hecmw_io_mesh.h"


static char grid_filename[HECMW_FILENAME_LEN+1] = "Unknown";

/*----------------------------------------------------------------------------*/

static void
do_logging(int loglv, int msgno, int add_location, const char *fmt, va_list ap)
{
	char line[100] = "";
	char msg[HECMW_MSG_LEN+1];

	HECMW_vsnprintf(msg, sizeof(msg), fmt, ap);
	if(add_location) {
		char *s = "";
		if(strlen(msg) > 0) s = ": ";
		HECMW_snprintf(line, sizeof(line), "%s:%d%s", grid_filename, HECMW_gflex_get_lineno(), s);
	}
	if(loglv == HECMW_LOG_ERROR) {
		HECMW_set_error(msgno, "%s%s", line, msg);
	} else {
		HECMW_print_msg(loglv, msgno, "%s%s", line, msg);
	}
}


static void
set_err(int msgno, const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	do_logging(HECMW_LOG_ERROR, msgno, 1, fmt, ap);
	va_end(ap);
}


static void
set_err_token(int token, int msgno, const char *fmt, ...)
{
	int msg_no;
	va_list ap;

	if(!token) {
		msg_no = HECMW_IO_GEOFEM_E0003;
	} else {
		msg_no = msgno;
	}
	va_start(ap, fmt); 
	do_logging(HECMW_LOG_ERROR, msg_no, 1, fmt, ap);
	va_end(ap);
}


/*-----------------------------------------------------------------------------
 * read functions
 */

static int 
read_pe(void)
{
	int token,n_neighbor_pe;
	
	/* PE-ID */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "PE-ID required");
		return -1;
	}

	/* NEIBPEtot */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "NEIBOEtot required");
		return -1;
	}
	n_neighbor_pe = HECMW_gflex_get_number();
	if(n_neighbor_pe < 0) {
		set_err(HECMW_IO_GEOFEM_E0301, "");
		return -1;
	}
	if(n_neighbor_pe != 0) {
		set_err(HECMW_IO_GEOFEM_E0301, "");
		return -1;
	}

	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}

	/* ESSENTIAL BLANK LINE */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
		return -1;
	}

	return 0;
}


static int
read_node(void)
{
	int nnode,nninternal,i,token;

	/* NODtot */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004,  "");
		return -1;
	}
	nnode = HECMW_gflex_get_number();
	if(nnode <= 0) {
		set_err(HECMW_IO_GEOFEM_E0311, "");
		return -1;
	}

	/* intNODtot */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}
	nninternal = HECMW_gflex_get_number();
	if(nninternal <= 0) {
		set_err(HECMW_IO_GEOFEM_E0312, "");
		return -1;
	}
	if(nnode != nninternal) {
		set_err(HECMW_IO_GEOFEM_E0313, "");
		return -1;
	}

	for(i=0; i < nnode; i++) {
		int id;
		double x,y,z;

		/* nGlobalID */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		id = HECMW_gflex_get_number();
		if(id <= 0) {
			set_err(HECMW_IO_GEOFEM_E0314, "");
			return -1;
		}

		/* X */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_DOUBLE && token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		x = HECMW_gflex_get_number();

		/* Y */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_DOUBLE && token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		y = HECMW_gflex_get_number();

		/* Z */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_DOUBLE && token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		z = HECMW_gflex_get_number();

		/* add node */
		if(HECMW_io_add_node(id, x, y, z) == NULL) {
			return -1;
		}

		/* add node to group */
		if(HECMW_io_add_ngrp("ALL", 1, &id) < 0) {	/* always add to 'ALL' */
			return -1;
		}
	}

	/* end of NODE */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}

	return 0;
}


static int
read_elem(void)
{
	int i,j,n,token,nelem;
	int *elem_type;

	/* ELMtot */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}
	nelem = HECMW_gflex_get_number();
	if(nelem <= 0) {
		set_err(HECMW_IO_GEOFEM_E0321, "");
		return -1;
	}

	/* ELMtype */
	elem_type = HECMW_malloc(sizeof(*elem_type)*nelem);
	if(elem_type == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	for(i=0; i < nelem; i++) {
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		elem_type[i] = HECMW_gflex_get_number();
		n = HECMW_get_max_node(HECMW_get_etype_GeoFEM2HECMW(elem_type[i]));
		if(n == -1) {
			set_err(HECMW_IO_GEOFEM_E0322, "");
			return -1;
		}
	}
	
	/* eGlobalID, connectivity */
	for(i=0; i < nelem; i++) {
		int id,hecmw_type;
		int node[HECMW_MAX_NODE_MAX];

		/* eGlobalID */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		id = HECMW_gflex_get_number() ;
		if(id <= 0) {
			set_err(HECMW_IO_GEOFEM_E0324, "Invalid Element ID");
			return -1;
		}

		hecmw_type = HECMW_get_etype_GeoFEM2HECMW(elem_type[i]);
		n = HECMW_get_max_node(hecmw_type);
		for(j=0; j < n; j++) {
			/* connectivity */
			token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
			if(token != HECMW_GFLEX_INT) {
				set_err_token(token, HECMW_IO_GEOFEM_E0004,"");
				return -1;
			}
			node[j] = HECMW_gflex_get_number();
			if(node[j] <= 0) {
				set_err(HECMW_IO_GEOFEM_E0323, "");
				return -1;
			}
		}

		/* add element */
		if(HECMW_io_add_elem(id, hecmw_type, node, 0, NULL) == NULL) {
			return -1;
		}

		/* add element to eroup */
		if(HECMW_io_add_egrp("ALL", 1, &id) < 0) {	/* always add to 'ALL' */
				return -1;
		}
	}

	/* end of ELEMENT */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}

	HECMW_free(elem_type);

	return 0;
}


static int
read_import(void)
{
	int token;

	/* ESSENTIAL BLANK LINE */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
		return -1;
	}
	return 0;
}


static int
read_export(void)
{
	int token;

	/* ESSENTIAL BLANK LINE */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
		return -1;
	}
	return 0;
}


static int
read_ngrp(void)
{
	int i,j,n,token,ngrp;
	int *grp_index;

	/* NODgrpTOT */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}
	ngrp = HECMW_gflex_get_number();
	if(ngrp < 0) {
		set_err(HECMW_IO_GEOFEM_E0341, "");
		return -1;
	}
	if(ngrp == 0) {
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		/* ESSENTIAL BLANK LINE */
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
			return -1;
		}
		return 0;
	}

	/* NODgrpINDEX */
	grp_index = HECMW_malloc(sizeof(*grp_index)*(ngrp+1));
	if(grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp_index[0] = 0;
	for(i=0; i < ngrp; i++) {
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		grp_index[i+1] = HECMW_gflex_get_number();
		if(grp_index[i+1] <= 0) {
			set_err(HECMW_IO_GEOFEM_E0342, "");
			return -1;
		}
	}

	for(i=0; i < ngrp; i++) {
		char *p,name[HECMW_NAME_LEN+1];
		int *node;

		/* NODgrpNAME */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_NAME) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		p = HECMW_gflex_get_text();
		if(strlen(p) > HECMW_NAME_LEN) {
			set_err(HECMW_IO_E0001, "");
			return -1;
		}
		strcpy(name, p);

		/* NODgrpITEM */
		n = grp_index[i+1] - grp_index[i];
		node = HECMW_malloc(sizeof(*node)*n);
		if(node == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(j=0; j < n; j++) {
			token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
			if(token != HECMW_GFLEX_INT) {
				set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
				return -1;
			}
			node[j] = HECMW_gflex_get_number();
			if(node[j] <= 0) {
				set_err(HECMW_IO_GEOFEM_E0343, "");
				return -1;
			}
		}

		/* add node to node group */
		if(HECMW_io_add_ngrp(name, n, node) < 0) {
			return -1;
		}

		HECMW_free(node);
	}

	/* end of NGRP */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}

	HECMW_free(grp_index);

	return 0;
}


static int
read_egrp(void)
{
	int i,j,n,token,ngrp;
	int *grp_index;

	/* ELMgrpTOT */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}
	ngrp = HECMW_gflex_get_number();
	if(ngrp < 0) {
		set_err(HECMW_IO_GEOFEM_E0351, "");
		return -1;
	}
	if(ngrp == 0) {
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		/* ESSENTIAL BLANK LINE */
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
			return -1;
		}
		return 0;
	}

	/* ELMgrpINDEX */
	grp_index = HECMW_malloc(sizeof(*grp_index)*(ngrp+1));
	if(grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp_index[0] = 0;
	for(i=0; i < ngrp; i++) {
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		grp_index[i+1] = HECMW_gflex_get_number();
		if(grp_index[i+1] <= 0) {
			set_err(HECMW_IO_GEOFEM_E0352, "");
			return -1;
		}
	}

	for(i=0; i < ngrp; i++) {
		char *p,name[HECMW_NAME_LEN+1];
		int *elem;

		/* ELMgrpNAME */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_NAME) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		p = HECMW_gflex_get_text();
		if(strlen(p) > HECMW_NAME_LEN) {
			set_err(HECMW_IO_E0001, "");
			return -1;
		}
		strcpy(name, p);

		/* ELMgrpITEM */
		n = grp_index[i+1] - grp_index[i];
		elem = HECMW_malloc(sizeof(*elem)*n);
		if(elem == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(j=0; j < n; j++) {
			token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
			if(token != HECMW_GFLEX_INT) {
				set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
				return -1;
			}
			elem[j] = HECMW_gflex_get_number();
			if(elem[j] <= 0) {
				set_err(HECMW_IO_GEOFEM_E0353, "");
				return -1;
			}
		}

		/* add node to node group */
		if(HECMW_io_add_egrp(name, n, elem) < 0) {
			return -1;
		}

		HECMW_free(elem);
	}

	/* end of EGRP */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}

	HECMW_free(grp_index);

	return 0;
}


static int
read_sgrp(void)
{
	int i,j,n,token,ngrp;
	int *grp_index;

	/* SUFgrpTOT */
	token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
	if(token != HECMW_GFLEX_INT) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
		return -1;
	}
	ngrp = HECMW_gflex_get_number();
	if(ngrp < 0) {
		set_err(HECMW_IO_GEOFEM_E0361, "");
		return -1;
	}
	if(ngrp == 0) {
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		/* ESSENTIAL BLANK LINE */
		token = HECMW_gflex_next_token();
		if(token != HECMW_GFLEX_NL) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "Needs ESSENTIAL BLANK LINE");
			return -1;
		}
		return 0;
	}

	/* SUFgrpINDEX */
	grp_index = HECMW_malloc(sizeof(*grp_index)*(ngrp+1));
	if(grp_index == NULL) {
		HECMW_set_error(errno, "");
		return -1;
	}
	grp_index[0] = 0;
	for(i=0; i < ngrp; i++) {
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_INT) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		grp_index[i+1] = HECMW_gflex_get_number();
		if(grp_index[i+1] <= 0) {
			set_err(HECMW_IO_GEOFEM_E0362, "");
			return -1;
		}
	}

	for(i=0; i < ngrp; i++) {
		char *p,name[HECMW_NAME_LEN+1];
		int *elem,*surf;

		/* SUFgrpNAME */
		token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
		if(token != HECMW_GFLEX_NAME) {
			set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
			return -1;
		}
		p = HECMW_gflex_get_text();
		if(strlen(p) > HECMW_NAME_LEN) {
			set_err(HECMW_IO_E0001, "");
			return -1;
		}
		strcpy(name, p);

		/* SUFgrpITEM */
		n = grp_index[i+1] - grp_index[i];
		elem = HECMW_malloc(sizeof(*elem)*n);
		if(elem == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		surf = HECMW_malloc(sizeof(*surf)*n);
		if(surf == NULL) {
			HECMW_set_error(errno, "");
			return -1;
		}
		for(j=0; j < n; j++) {
			token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
			if(token != HECMW_GFLEX_INT) {
				set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
				return -1;
			}
			elem[j] = HECMW_gflex_get_number();
			if(elem[j] <= 0) {
				set_err(HECMW_IO_GEOFEM_E0363, "");
				return -1;
			}
		}
		for(j=0; j < n; j++) {
			token = HECMW_gflex_next_token_skip(HECMW_GFLEX_NL);
			if(token != HECMW_GFLEX_INT) {
				set_err_token(token, HECMW_IO_GEOFEM_E0004, "");
				return -1;
			}
			surf[j] = HECMW_gflex_get_number();
			if(surf[j] <= 0) {
				set_err(HECMW_IO_GEOFEM_E0363, "");
				return -1;
			}
		}

		/* add to surf group */
		if(HECMW_io_add_sgrp(name, n, elem, surf) < 0) {
			return -1;
		}

		HECMW_free(elem);
		HECMW_free(surf);
	}

	/* end of SGRP */
	token = HECMW_gflex_next_token();
	if(token != HECMW_GFLEX_NL) {
		set_err_token(token, HECMW_IO_GEOFEM_E0004,"");
		return -1;
	}

	HECMW_free(grp_index);

	return 0;
}


static int
parse(void)
{
	if(read_pe()) return -1;
	if(read_node()) return -1;
	if(read_elem()) return -1;
	if(read_import()) return -1;
	if(read_export()) return -1;
	if(read_ngrp()) return -1;
	if(read_egrp()) return -1;
	if(read_sgrp()) return -1;

	return 0; 
}




/* read only. Not make hecmwST_local_mesh */
int
HECMW_read_geofem_mesh(const char *filename)
{
	FILE *fp;

	HECMW_log(HECMW_LOG_DEBUG, "Start to read GeoFEM mesh");

	if(filename == NULL) { 
		HECMW_set_error(HECMW_IO_E0001, "Not specified filename for GeoFEM mesh input routine");
		return -1;
	}
	HECMW_log(HECMW_LOG_DEBUG, "GeoFEM mesh file is '%s'", filename);

	if(strlen(filename) > HECMW_FILENAME_LEN) {
		HECMW_set_error(HECMW_IO_E0002, "");
		return -1;
	}

	strcpy(grid_filename, filename);
	HECMW_io_set_gridfile(grid_filename);

	if((fp = fopen(filename, "r")) == NULL) {
		HECMW_set_error(HECMW_IO_HEC_E0001, "File: %s, %s", filename, strerror(errno));
		return -1;
	}

	if(HECMW_gflex_set_input(fp)) return -1;

	HECMW_log(HECMW_LOG_DEBUG, "Parsing...");
	if(parse()) {
		return -1;
	}

	if(fclose(fp)) {
		HECMW_set_error(HECMW_IO_HEC_E0002, "File: %s, %s", filename, strerror(errno));
		return -1;
	}

	strcpy(grid_filename, "Unknown");

	return 0;
}

struct hecmwST_local_mesh *
HECMW_get_geofem_mesh(const char *filename)
{
	struct hecmwST_local_mesh *local_mesh;

	if(HECMW_io_init()) return NULL;
	if(HECMW_io_pre_process()) return NULL;
	if(HECMW_read_geofem_mesh(filename)) return NULL;
	if(HECMW_io_post_process()) return NULL;
	local_mesh = HECMW_io_make_local_mesh();
	if(local_mesh == NULL) return NULL;
	if(HECMW_io_finalize()) return NULL;

	strcpy(grid_filename, "Unknown");

	return local_mesh;
}
