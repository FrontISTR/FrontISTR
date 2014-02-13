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



#include "hecmw_util.h"
#include "hecmw_common_define.h"
#include "hecmw_etype.h"
#include "hecmw_conn_conv.h"


struct conn_conv {
	int hecmw_etype;						
	int connectivity[HECMW_MAX_NODE_MAX];	
};



static struct conn_conv conn_conv_abaqus[] = {
	{ 232, { 1,2,3,6,4,5 }},
	{ 342, { 1,2,3,4,7,5,6,8,9,10 }},
	{ 352, { 1,2,3,4,5,6,9,7,8,12,10,11,13,14,15 }},
	{ 542, { 1,2,3,4,9,10,11,12,5,6,7,8,13,14,15,16 }},
	{ -1,  {-1}}	/* terminator */
};

#if 0
static struct conn_conv conn_conv_nastran[] = {
	{ -1,  {-1}}	/* terminator */
};
#endif


struct conn_order {
	int node;
	int hecmw_order;
};



static int 
conn_comp(const void *c1, const void *c2)
{
	int co1 = ((struct conn_order *)c1)->hecmw_order;
	int co2 = ((struct conn_order *)c2)->hecmw_order;

	if(co1 == co2) return 0;
	if(co1 < co2) return -1;
	return 1;
}


int
HECMW_convert_connectivity(int from, int hecmw_etype, int *conn)
{
	int i,j,n;
	struct conn_conv *from_table;
	struct conn_order order[HECMW_MAX_NODE_MAX];

	if (conn == NULL) {
		HECMW_set_error(HECMW_ALL_E0101, "Connectivity contversion: 'conn' is NULL");
		return -1;
	}

	switch (from) {
		case HECMW_CONNTYPE_HECMW:
			return 0;
		case HECMW_CONNTYPE_ABAQUS:
			from_table = conn_conv_abaqus;
			break;
#if 0
		case HECMW_CONNTYPE_NASTRAN:
			from_table = conn_conv_nastran;
			break;
#endif
		default:
			HECMW_set_error(HECMW_ALL_E0101, "Connectivity conversion: Unsupported connectivity type");
			return -1;
	}

	if ((n = HECMW_get_max_node(hecmw_etype)) == -1) {
		HECMW_set_error(HECMW_ALL_E0101, "Connectivity conversion: Invalid 'hecmw_etype'");
		return -1;
	}

	for(i=0; from_table[i].hecmw_etype != -1; i++) {
		if(from_table[i].hecmw_etype != hecmw_etype) continue;
		for(j=0; j < n; j++) {
			order[j].node = conn[j];
			order[j].hecmw_order = from_table[i].connectivity[j];
		}
		qsort(order, n, sizeof(*order), conn_comp);
		for(j=0; j < n; j++) {
			conn[j] = order[j].node;
		}
		break;
	}

	return 0;
}

