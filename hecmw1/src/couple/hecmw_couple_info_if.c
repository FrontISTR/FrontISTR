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

#include "hecmw_struct.h"
#include "hecmw_lib_fc.h"
#include "hecmw_couple_info.h"

/*================================================================================================*/

extern void
hecmw_couple_get_unit_id_if(char *boundary_id, int *unit_specifier, char *buf, int *err,
		int id_len, int buf_len)
{
	char c_boundary_id[HECMW_NAME_LEN+1];
	char c_buf[HECMW_NAME_LEN+1];

	*err = 1;

	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_boundary_id, sizeof(c_boundary_id)) == NULL) return;
	if(HECMW_couple_get_unit_id(c_boundary_id, *unit_specifier, c_buf, sizeof(c_buf)) == NULL) return;
	if(HECMW_strcpy_c2f(c_buf, buf, buf_len) == 0) return;

	*err = 0;
}



extern void
hecmw_couple_get_unit_id_if_(char *boundary_id, int *unit_specifier, char *buf, int *err,
		int id_len, int buf_len)
{
	hecmw_couple_get_unit_id_if(boundary_id, unit_specifier, buf, err, id_len, buf_len);
}



extern void
hecmw_couple_get_unit_id_if__(char *boundary_id, int *unit_specifier, char *buf, int *err,
		int id_len, int buf_len)
{
	hecmw_couple_get_unit_id_if(boundary_id, unit_specifier, buf, err, id_len, buf_len);
}



extern void
HECMW_COUPLE_GET_UNIT_ID_IF(char *boundary_id, int *unit_specifier, char *buf, int *err,
		int id_len, int buf_len)
{
	hecmw_couple_get_unit_id_if(boundary_id, unit_specifier, buf, err, id_len, buf_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_memb_if(char *boundary_id, int *is_member, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_member = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_member = HECMW_couple_is_member(c_id);
}



extern void
hecmw_couple_is_memb_if_(char *boundary_id, int *is_member, int id_len)
{
	hecmw_couple_is_memb_if(boundary_id, is_member, id_len);
}



extern void
hecmw_couple_is_memb_if__(char *boundary_id, int *is_member, int id_len)
{
	hecmw_couple_is_memb_if(boundary_id, is_member, id_len);
}



extern void
HECMW_COUPLE_IS_MEMB_IF(char *boundary_id, int *is_member, int id_len)
{
	hecmw_couple_is_memb_if(boundary_id, is_member, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_unit_memb_if(char *boundary_id, int *unit_specifier, int *is_member, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_member = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_member = HECMW_couple_is_unit_member(c_id, *unit_specifier);
}



extern void
hecmw_couple_is_unit_memb_if_(char *boundary_id, int *unit_specifier, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_if(boundary_id, unit_specifier, is_member, id_len);
}



extern void
hecmw_couple_is_unit_memb_if__(char *boundary_id, int *unit_specifier, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_if(boundary_id, unit_specifier, is_member, id_len);
}



extern void
HECMW_COUPLE_IS_UNIT_MEMB_IF(char *boundary_id, int *unit_specifier, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_if(boundary_id, unit_specifier, is_member, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_unit_memb_u_if(char *unit_id, int *is_member, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_member = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_member = HECMW_couple_is_unit_member_u(c_id);
}



extern void
hecmw_couple_is_unit_memb_u_if_(char *unit_id, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_u_if(unit_id, is_member, id_len);
}



extern void
hecmw_couple_is_unit_memb_u_if__(char *unit_id, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_u_if(unit_id, is_member, id_len);
}



extern void
HECMW_COUPLE_IS_UNIT_MEMB_U_IF(char *unit_id, int *is_member, int id_len)
{
	hecmw_couple_is_unit_memb_u_if(unit_id, is_member, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_root_if(char *boundary_id, int *is_root, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_root = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_root = HECMW_couple_is_root(c_id);
}



extern void
hecmw_couple_is_root_if_(char *boundary_id, int *is_root, int id_len)
{
	hecmw_couple_is_root_if(boundary_id, is_root, id_len);
}



extern void
hecmw_couple_is_root_if__(char *boundary_id, int *is_root, int id_len)
{
	hecmw_couple_is_root_if(boundary_id, is_root, id_len);
}



extern void
HECMW_COUPLE_IS_ROOT_IF(char *boundary_id, int *is_root, int id_len)
{
	hecmw_couple_is_root_if(boundary_id, is_root, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_unit_root_if(char *boundary_id, int *unit_specifier, int *is_root, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_root = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_root = HECMW_couple_is_unit_root(c_id, *unit_specifier);
}



extern void
hecmw_couple_is_unit_root_if_(char *boundary_id, int *unit_specifier, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_if(boundary_id, unit_specifier, is_root, id_len);
}



extern void
hecmw_couple_is_unit_root_if__(char *boundary_id, int *unit_specifier, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_if(boundary_id, unit_specifier, is_root, id_len);
}



extern void
HECMW_COUPLE_IS_UNIT_ROOT_IF(char *boundary_id, int *unit_specifier, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_if(boundary_id, unit_specifier, is_root, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_couple_is_unit_root_u_if(char *unit_id, int *is_root, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*is_root = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*is_root = HECMW_couple_is_unit_root_u(c_id);
}



extern void
hecmw_couple_is_unit_root_u_if_(char *unit_id, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_u_if(unit_id, is_root, id_len);
}



extern void
hecmw_couple_is_unit_root_u_if__(char *unit_id, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_u_if(unit_id, is_root, id_len);
}



extern void
HECMW_COUPLE_IS_UNIT_ROOT_U_IF(char *unit_id, int *is_root, int id_len)
{
	hecmw_couple_is_unit_root_u_if(unit_id, is_root, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intercomm_get_size_if(char *boundary_id, int *psize, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*psize = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*psize = HECMW_intercomm_get_size(c_id);
}



extern void
hecmw_intercomm_get_size_if_(char *boundary_id, int *psize, int id_len)
{
	hecmw_intercomm_get_size_if(boundary_id, psize, id_len);
}



extern void
hecmw_intercomm_get_size_if__(char *boundary_id, int *psize, int id_len)
{
	hecmw_intercomm_get_size_if(boundary_id, psize, id_len);
}



extern void
HECMW_INTERCOMM_GET_SIZE_IF(char *boundary_id, int *psize, int id_len)
{
	hecmw_intercomm_get_size_if(boundary_id, psize, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_size_if(char *boundary_id, int *unit_specifier, int *psize, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*psize = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*psize = HECMW_intracomm_get_size(c_id, *unit_specifier);
}



extern void
hecmw_intracomm_get_size_if_(char *boundary_id, int *unit_specifier, int *psize, int id_len)
{
	hecmw_intracomm_get_size_if(boundary_id, unit_specifier, psize, id_len);
}



extern void
hecmw_intracomm_get_size_if__(char *boundary_id, int *unit_specifier, int *psize, int id_len)
{
	hecmw_intracomm_get_size_if(boundary_id, unit_specifier, psize, id_len);
}



extern void
HECMW_INTRACOMM_GET_SIZE_IF(char *boundary_id, int *unit_specifier, int *psize, int id_len)
{
	hecmw_intracomm_get_size_if(boundary_id, unit_specifier, psize, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_size_u_if(char *unit_id, int *psize, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*psize = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*psize = HECMW_intracomm_get_size_u(c_id);
}



extern void
hecmw_intracomm_get_size_u_if_(char *unit_id, int *psize, int id_len)
{
	hecmw_intracomm_get_size_u_if(unit_id, psize, id_len);
}



extern void
hecmw_intracomm_get_size_u_if__(char *unit_id, int *psize, int id_len)
{
	hecmw_intracomm_get_size_u_if(unit_id, psize, id_len);
}



extern void
HECMW_INTRACOMM_GET_SIZE_U_IF(char *unit_id, int *psize, int id_len)
{
	hecmw_intracomm_get_size_u_if(unit_id, psize, id_len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_intercomm_get_rank_if(char *boundary_id, int *rank, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*rank = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*rank = HECMW_intercomm_get_rank(c_id);
}



extern void
hecmw_intercomm_get_rank_if_(char *boundary_id, int *rank, int id_len)
{
	hecmw_intercomm_get_rank_if(boundary_id, rank, id_len);
}



extern void
hecmw_intercomm_get_rank_if__(char *boundary_id, int *rank, int id_len)
{
	hecmw_intercomm_get_rank_if(boundary_id, rank, id_len);
}



extern void
HECMW_INTERCOMM_GET_RANK_IF(char *boundary_id, int *rank, int id_len)
{
	hecmw_intercomm_get_rank_if(boundary_id, rank, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_rank_if(char *boundary_id, int *unit_specifier, int *rank, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*rank = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*rank = HECMW_intracomm_get_rank(c_id, *unit_specifier);
}



extern void
hecmw_intracomm_get_rank_if_(char *boundary_id, int *unit_specifier, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_if(boundary_id, unit_specifier, rank, id_len);
}



extern void
hecmw_intracomm_get_rank_if__(char *boundary_id, int *unit_specifier, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_if(boundary_id, unit_specifier, rank, id_len);
}



extern void
HECMW_INTRACOMM_GET_RANK_IF(char *boundary_id, int *unit_specifier, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_if(boundary_id, unit_specifier, rank, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_rank_u_if(char *unit_id, int *rank, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*rank = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*rank = HECMW_intracomm_get_rank_u(c_id);
}



extern void
hecmw_intracomm_get_rank_u_if_(char *unit_id, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_u_if(unit_id, rank, id_len);
}



extern void
hecmw_intracomm_get_rank_u_if__(char *unit_id, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_u_if(unit_id, rank, id_len);
}



extern void
HECMW_INTRACOMM_GET_RANK_U_IF(char *unit_id, int *rank, int id_len)
{
	hecmw_intracomm_get_rank_u_if(unit_id, rank, id_len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_intercomm_get_comm_if(char *boundary_id, int *comm, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*comm = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*comm = (int)HECMW_intercomm_get_comm(c_id);
}



extern void
hecmw_intercomm_get_comm_if_(char *boundary_id, int *comm, int id_len)
{
	hecmw_intercomm_get_comm_if(boundary_id, comm, id_len);
}



extern void
hecmw_intercomm_get_comm_if__(char *boundary_id, int *comm, int id_len)
{
	hecmw_intercomm_get_comm_if(boundary_id, comm, id_len);
}



extern void
HECMW_INTERCOMM_GET_COMM_IF(char *boundary_id, int *comm, int id_len)
{
	hecmw_intercomm_get_comm_if(boundary_id, comm, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_comm_if(char *boundary_id, int *unit_specifier, int *comm, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*comm = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*comm = (int)HECMW_intracomm_get_comm(c_id, *unit_specifier);
}



extern void
hecmw_intracomm_get_comm_if_(char *boundary_id, int *unit_specifier, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_if(boundary_id, unit_specifier, comm, id_len);
}



extern void
hecmw_intracomm_get_comm_if__(char *boundary_id, int *unit_specifier, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_if(boundary_id, unit_specifier, comm, id_len);
}



extern void
HECMW_INTRACOMM_GET_COMM_IF(char *boundary_id, int *unit_specifier, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_if(boundary_id, unit_specifier, comm, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_comm_u_if(char *unit_id, int *comm, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*comm = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*comm = (int)HECMW_intracomm_get_comm_u(c_id);
}



extern void
hecmw_intracomm_get_comm_u_if_(char *unit_id, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_u_if(unit_id, comm, id_len);
}



extern void
hecmw_intracomm_get_comm_u_if__(char *unit_id, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_u_if(unit_id, comm, id_len);
}



extern void
HECMW_INTRACOMM_GET_COMM_U_IF(char *unit_id, int *comm, int id_len)
{
	hecmw_intracomm_get_comm_u_if(unit_id, comm, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intercomm_get_group_if(char *boundary_id, int *group, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*group = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*group = (int)HECMW_intercomm_get_group(c_id);
}



extern void
hecmw_intercomm_get_group_if_(char *boundary_id, int *group, int id_len)
{
	hecmw_intercomm_get_group_if(boundary_id, group, id_len);
}



extern void
hecmw_intercomm_get_group_if__(char *boundary_id, int *group, int id_len)
{
	hecmw_intercomm_get_group_if(boundary_id, group, id_len);
}



extern void
HECMW_INTERCOMM_GET_GROUP_IF(char *boundary_id, int *group, int id_len)
{
	hecmw_intercomm_get_group_if(boundary_id, group, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_group_if(char *boundary_id, int *unit_specifier, int *group, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*group = -1;
	if(HECMW_strcpy_f2c_r(boundary_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*group = (int)HECMW_intracomm_get_group(c_id, *unit_specifier);
}



extern void
hecmw_intracomm_get_group_if_(char *boundary_id, int *unit_specifier, int *group, int id_len)
{
	hecmw_intracomm_get_group_if(boundary_id, unit_specifier, group, id_len);
}



extern void
hecmw_intracomm_get_group_if__(char *boundary_id, int *unit_specifier, int *group, int id_len)
{
	hecmw_intracomm_get_group_if(boundary_id, unit_specifier, group, id_len);
}



extern void
HECMW_INTRACOMM_GET_GROUP_IF(char *boundary_id, int *unit_specifier, int *group, int id_len)
{
	hecmw_intracomm_get_group_if(boundary_id, unit_specifier, group, id_len);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

extern void
hecmw_intracomm_get_group_u_if(char *unit_id, int *group, int id_len)
{
	char c_id[HECMW_NAME_LEN+1];

	*group = -1;
	if(HECMW_strcpy_f2c_r(unit_id, id_len, c_id, sizeof(c_id)) == NULL) return;
	*group = (int)HECMW_intracomm_get_group_u(c_id);
}



extern void
hecmw_intracomm_get_group_u_if_(char *unit_id, int *group, int id_len)
{
	hecmw_intracomm_get_group_u_if(unit_id, group, id_len);
}



extern void
hecmw_intracomm_get_group_u_if__(char *unit_id, int *group, int id_len)
{
	hecmw_intracomm_get_group_u_if(unit_id, group, id_len);
}



extern void
HECMW_INTRACOMM_GET_GROUP_U_IF(char *unit_id, int *group, int id_len)
{
	hecmw_intracomm_get_group_u_if(unit_id, group, id_len);
}


/*------------------------------------------------------------------------------------------------*/

extern void
hecmw_couple_comm_init_if(int *err)
{
	if(HECMW_couple_comm_init() != HECMW_SUCCESS) {
		*err = 1;
	} else {
		*err = 0;
	}
}



extern void
hecmw_couple_comm_init_if_(int *err)
{
	hecmw_couple_comm_init_if(err);
}



extern void
hecmw_couple_comm_init_if__(int *err)
{
	hecmw_couple_comm_init_if(err);
}



extern void
HECMW_COUPLE_COMM_INIT_IF(int *err)
{
	hecmw_couple_comm_init_if(err);
}
