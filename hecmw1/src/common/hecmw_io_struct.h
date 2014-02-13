/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/06/29                                        *
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



#ifndef HECMW_IO_STRUCT_INCLUDED
#define HECMW_IO_STRUCT_INCLUDED

#include "hecmw_config.h"
#include "hecmw_set_int.h"


struct hecmw_io_id_array {
	int n;
	int *id;
};


struct hecmw_io_id {
	int id;
	struct hecmw_io_id *next;
};


struct hecmw_io_header {
	char header[HECMW_HEADER_LEN+1];
};


struct hecmw_io_zero {
	double zero;
};


struct hecmw_io_node {
	double x;
	double y;
	double z;
};


struct hecmw_io_element {
	int type;
	int *node;
	int nmatitem;
	double *matitem;
	char matname[HECMW_NAME_LEN+1];	/* created material name */
	int mpc_matid;	/* for element type 9XX */
	int mpc_sectid;	/* for element type 9XX */
};


struct hecmw_io_ngrp {
	char name[HECMW_NAME_LEN+1];
	struct hecmw_set_int *node;
	struct hecmw_io_ngrp *next;
};


struct hecmw_io_egrp {
	char name[HECMW_NAME_LEN+1];
	struct hecmw_set_int *elem;
	struct hecmw_io_egrp *next;
};


struct hecmw_io_sgrp {
	char name[HECMW_NAME_LEN+1];
	struct hecmw_set_int *item;
	struct hecmw_io_sgrp *next;
};


struct hecmw_io_mpc {
	int neq;
	double cnst;

	struct hecmw_io_mpcitem {
		char ngrp[HECMW_NAME_LEN+1];	/* valid if node == -1 */
		int node;
		int dof;
		double a;
	} *item;	/* neq */
	struct hecmw_io_mpc *next;
};


struct hecmw_io_amplitude {
	char name[HECMW_NAME_LEN+1];
	int type_def;
	int type_time;
	int type_val;

	struct hecmw_io_amplitude_item {
		double val;
		double table;
		struct hecmw_io_amplitude_item *next;
	} *item;
	struct hecmw_io_amplitude_item *last;
	struct hecmw_io_amplitude *next;
};


struct hecmw_io_initial {
	int type;
#define HECMW_INITIAL_TYPE_TEMPERATURE 1
	int node;
	char ngrp[HECMW_NAME_LEN+1];	/* valid if node == -1 */
	double val;
	struct hecmw_io_initial *next;
};


struct hecmw_io_material {
	char name[HECMW_NAME_LEN+1];
	int nitem;

	struct hecmw_io_matitem {
		int item;
		int nval;

		struct hecmw_io_matsubitem {
			double *val;
			double temp;
			struct hecmw_io_matsubitem *next;
		} *subitem;
	} *item;
	struct hecmw_io_material *next;
};


struct hecmw_io_section {
	char egrp[HECMW_NAME_LEN+1];
	char material[HECMW_NAME_LEN+1];
	int composite;
	int secopt;
	int type;

	union hecmw_io_section_item {

		struct hecmw_io_section_solid {
			double thickness;
		} solid;

		struct hecmw_io_section_shell {
			double thickness;
			int integpoints;
		} shell;
		
        struct hecmw_io_section_beam {
			double vxyz[3];
			double area;
            double Iyy;
            double Izz;
            double Jx;
		} beam;

		struct hecmw_io_section_interface {
			double thickness;
			double gapcon;
			double gaprad1;
			double gaprad2;
		} interface;
	} sect;
	struct hecmw_io_section *next;
};

struct hecmw_io_contact {
	char name[HECMW_NAME_LEN+1];
	int type;
	char master_grp[HECMW_NAME_LEN+1];
	char slave_grp[HECMW_NAME_LEN+1];
	struct hecmw_io_contact *next;
};

#endif
