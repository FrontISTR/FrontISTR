/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/05/02                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuaki Sakane (RIST)                         *
 *                       Shin'ichi Ezure (RIST)                        *
 *                                                                     *
 *     Contact address :  IIS,The University of Tokyo RSS21 project    *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/



#include <stdio.h>
#include "hecmw_config.h"
#include "hecmw_comm.h"
#include "hecmw_util.h"

static int is_initialized;

static HECMW_Comm hecmw_comm;
static HECMW_Group hecmw_group;
static int comm_size;
static int comm_rank;


/*-----------------------------------------------------------------------------*/



int
HECMW_Comm_rank( HECMW_Comm comm, int *rank )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_rank( comm, rank );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_rank" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
    *rank=0;
    return 0;
#endif
}

int
HECMW_Comm_size( HECMW_Comm comm, int *size )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_size( comm, size );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_size" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	*size=1;
        return 0;
#endif
}

int
HECMW_Comm_dup( HECMW_Comm comm, HECMW_Comm *new_comm )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_dup(comm, new_comm);
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_dup" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	    *new_comm=0;
        return 0;
#endif
}



/*-----------------------------------------------------------------------------*/

extern int
HECMW_Barrier( HECMW_Comm comm )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Barrier( comm );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Barrier" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Wait( HECMW_Request *array_of_requests, HECMW_Status *array_of_statuses )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Wait( array_of_requests, array_of_statuses );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Waitall" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Waitall( int count, HECMW_Request *array_of_requests, HECMW_Status *array_of_statuses )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Waitall( count, array_of_requests, array_of_statuses );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Waitall" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Bcast( void *buffer, int count, HECMW_Datatype datatype, int root, HECMW_Comm comm )
{
#ifndef HECMW_SERIAL
    int rtc;

	if( datatype == HECMW_INT ) {
        rtc = MPI_Bcast( buffer, count, MPI_INT, root, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Bcast" );
            goto error;
        }
	} else if( datatype == HECMW_DOUBLE ) {
        rtc = MPI_Bcast( buffer, count, MPI_DOUBLE, root, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Bcast" );
            goto error;
        }
	} else if( datatype == HECMW_CHAR ) {
        rtc = MPI_Bcast( buffer, count, MPI_CHAR, root, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Bcast" );
            goto error;
        }
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Send( void *buffer, int count, HECMW_Datatype datatype, int dest, int tag, HECMW_Comm comm )
{
#ifndef HECMW_SERIAL
    int rtc;

	if( datatype == HECMW_INT ) {
        rtc = MPI_Send( buffer, count, MPI_INT, dest, tag, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Send" );
            goto error;
        }
	} else if( datatype == HECMW_DOUBLE ) {
        rtc = MPI_Send( buffer, count, MPI_DOUBLE, dest, tag, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Send" );
            goto error;
        }
	} else if( datatype == HECMW_CHAR ) {
        rtc = MPI_Send( buffer, count, MPI_CHAR, dest, tag, comm );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Send" );
            goto error;
        }
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Recv( void *buffer, int count, HECMW_Datatype datatype, int source,
            int tag, HECMW_Comm comm, HECMW_Status *status )
{
#ifndef HECMW_SERIAL
    int rtc;

	if( datatype == HECMW_INT ) {
        rtc = MPI_Recv( buffer, count, MPI_INT, source, tag, comm, status );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Recv" );
            goto error;
        }
	} else if( datatype == HECMW_DOUBLE ) {
        rtc = MPI_Recv( buffer, count, MPI_DOUBLE, source, tag, comm, status );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Recv" );
            goto error;
        }
	} else if( datatype == HECMW_CHAR ) {
        rtc = MPI_Recv( buffer, count, MPI_CHAR, source, tag, comm, status );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Recv" );
            goto error;
        }
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Isend( void *buffer, int count, HECMW_Datatype datatype, int dest,
             int tag, HECMW_Comm comm, HECMW_Request *request )
{
#ifndef HECMW_SERIAL
    int rtc;

	if( datatype == HECMW_INT ) {
        rtc = MPI_Isend( buffer, count, MPI_INT, dest, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
            goto error;
        }
	} else if( datatype == HECMW_DOUBLE ) {
        rtc = MPI_Isend( buffer, count, MPI_DOUBLE, dest, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
            goto error;
        }
	} else if( datatype == HECMW_CHAR ) {
        rtc = MPI_Isend( buffer, count, MPI_CHAR, dest, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Isend" );
            goto error;
        }
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Irecv( void *buffer, int count, HECMW_Datatype datatype, int source,
             int tag, HECMW_Comm comm, HECMW_Request *request )
{
#ifndef HECMW_SERIAL
    int rtc;

	if( datatype == HECMW_INT ) {
        rtc = MPI_Irecv( buffer, count, MPI_INT, source, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Irecv" );
            goto error;
        }
	} else if( datatype == HECMW_DOUBLE ) {
        rtc = MPI_Irecv( buffer, count, MPI_DOUBLE, source, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Irecv" );
            goto error;
        }
	} else if( datatype == HECMW_CHAR ) {
        rtc = MPI_Irecv( buffer, count, MPI_CHAR, source, tag, comm, request );
        if( rtc != MPI_SUCCESS ) {
            HECMW_set_error( HECMW_ALL_E1003, "MPI_Irecv" );
            goto error;
        }
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Allreduce( void *sendbuf, void *recvbuf, int count,
                 HECMW_Datatype datatype, HECMW_Op op, HECMW_Comm comm )
{
#ifndef HECMW_SERIAL
    MPI_Datatype _datatype;
    MPI_Op _op;
    int rtc;

	if( datatype == HECMW_INT ) {
        _datatype = MPI_INT;
	} else if( datatype == HECMW_DOUBLE ) {
        _datatype = MPI_DOUBLE;
	} else if( datatype == HECMW_CHAR ) {
        _datatype = MPI_CHAR;
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

	if( op == HECMW_SUM ) {
        _op = MPI_SUM;
	} else if( op == HECMW_MIN ) {
        _op = MPI_MIN;
	} else if( op == HECMW_MAX ) {
        _op = MPI_MAX;
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid operation is found" );
        goto error;
    }

    rtc = MPI_Allreduce( sendbuf, recvbuf, count, _datatype, _op, comm );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Allreduce" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Allgather( void *sendbuf, int sendcount, HECMW_Datatype sendtype,
                 void *recvbuf, int recvcount, HECMW_Datatype recvtype, HECMW_Comm comm )
{
#ifndef HECMW_SERIAL
    MPI_Datatype _sendtype, _recvtype;
    int rtc;

	if( sendtype == HECMW_INT ) {
        _sendtype = MPI_INT;
	} else if( sendtype == HECMW_DOUBLE ) {
        _sendtype = MPI_DOUBLE;
	} else if( sendtype == HECMW_CHAR ) {
        _sendtype = MPI_CHAR;
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

	if( recvtype == HECMW_INT ) {
        _recvtype = MPI_INT;
	} else if( recvtype == HECMW_DOUBLE ) {
        _recvtype = MPI_DOUBLE;
	} else if( recvtype == HECMW_CHAR ) {
        _recvtype = MPI_CHAR;
	} else {
        HECMW_set_error( HECMW_ALL_E1003, "Invalid data type is found" );
        goto error;
    }

    rtc = MPI_Allgather( sendbuf, sendcount, _sendtype, recvbuf, recvcount, _recvtype, comm );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Allgather" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Group_incl( HECMW_Group group, int n, int *ranks, HECMW_Group *newgroup )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Group_incl( group, n, ranks, newgroup );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Group_excl" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Group_excl( HECMW_Group group, int n, int *ranks, HECMW_Group *newgroup )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Group_excl( group, n, ranks, newgroup );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Group_excl" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Comm_create( HECMW_Comm comm, HECMW_Group group, HECMW_Comm *comm_out )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_create( comm, group, comm_out );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_create" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


extern int
HECMW_Group_rank( HECMW_Group group, int *rank )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Group_rank( group, rank );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Group_rank" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
    	*rank = 0;
	return 0;
#endif
}


extern int
HECMW_Group_size( HECMW_Group group, int *size )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Group_size( group, size );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Group_size" );
        goto error;
    }

    return 0;

error:
    return -1;
#else
    	*size = 1;
	return 0;
#endif
}


extern int
HECMW_Comm_group( HECMW_Comm comm, HECMW_Group *group )
{
#ifndef HECMW_SERIAL
    int rtc;

    rtc = MPI_Comm_group( comm, group );
    if( rtc != MPI_SUCCESS ) {
        HECMW_set_error( HECMW_ALL_E1003, "MPI_Comm_group" );
        goto error;
    } 

    return 0;

error:
    return -1;
#else
	return 0;
#endif
}


static int
check_is_initialized(void)
{
#ifdef HECMW_SERIAL
	is_initialized = 1;
#endif
	if(!is_initialized) {
		HECMW_set_error(HECMW_ALL_E1002, "");
		return 0;
	}
	return 1;
}


static int setup_comm(void)
{
#ifndef HECMW_SERIAL
	if(MPI_Comm_dup(MPI_COMM_WORLD, &hecmw_comm) != MPI_SUCCESS) {
		HECMW_set_error(HECMW_ALL_E1003, "");
		return -1;
	}
#else
	hecmw_comm = 0;
#endif
	return 0;
}


static int setup_comm_size(void)
{
#ifndef HECMW_SERIAL
	if(MPI_Comm_size(MPI_COMM_WORLD, &comm_size) != MPI_SUCCESS) {
		HECMW_set_error(HECMW_ALL_E1003, "");
		return -1;
	}
#else
	comm_size = 1;
#endif
	return 0;
}


static int setup_comm_rank(void)
{
#ifndef HECMW_SERIAL
	if(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank) != MPI_SUCCESS) {
		HECMW_set_error(HECMW_ALL_E1003, "");
		return -1;
	}
#else
	comm_rank = 0;
#endif
	return 0;
}


static int setup_comm_group(void)
{
#ifndef HECMW_SERIAL
    if(MPI_Comm_group(MPI_COMM_WORLD, &hecmw_group) != MPI_SUCCESS) {
        HECMW_set_error(HECMW_ALL_E1003, "");
        return -1;
    }
#else
    hecmw_group = 0;
#endif
    return 0;
}


int
HECMW_comm_init(int *argc, char ***argv)
{
#ifndef HECMW_SERIAL
	if(MPI_Init(argc, argv) != MPI_SUCCESS) {
		HECMW_set_error(HECMW_ALL_E1001, "");
		return -1;
	}
	is_initialized = 1;

	HECMW_log(HECMW_LOG_DEBUG, "MPI initialized");
#endif

	if(setup_comm()) {
		return -1;
	}
	if(setup_comm_size()) {
		return -1;
	}
	if(setup_comm_rank()) {
		return -1;
	}
    if(setup_comm_group()) {
        return -1;
    }

	return 0;
}


int
HECMW_comm_is_initialized(void)
{
	return is_initialized;
}


HECMW_Comm
HECMW_comm_get_comm(void)
{
	return check_is_initialized() ? hecmw_comm : (HECMW_Comm)-1;
}


int
HECMW_comm_get_size(void)
{
	return check_is_initialized() ? comm_size : -1;
}


int
HECMW_comm_get_rank(void)
{
	return check_is_initialized() ? comm_rank : -1;
}


HECMW_Group
HECMW_comm_get_group(void)
{
	return check_is_initialized() ? hecmw_group : (HECMW_Group)-1;
}


HECMW_Fint HECMW_Comm_c2f(HECMW_Comm comm)
{
#ifndef HECMW_SERIAL
  return MPI_Comm_c2f(comm);
#else
  return comm;
#endif
}

HECMW_Comm HECMW_Comm_f2c(HECMW_Fint comm)
{
#ifndef HECMW_SERIAL
  return MPI_Comm_f2c(comm);
#else
  return comm;
#endif
}

HECMW_Group HECMW_Group_f2c(HECMW_Fint group)
{
#ifndef HECMW_SERIAL
  return MPI_Group_f2c(group);
#else
  return group;
#endif
}

/*---------------------------------------------------------------------------*/


extern void
hecmw_comm_init_if(HECMW_Fint *comm, int *size, int *rank, HECMW_Fint *group)
{
	is_initialized = 1;

	hecmw_comm = HECMW_Comm_f2c(*comm);
	comm_size = *size;
	comm_rank = *rank;
	hecmw_group = HECMW_Group_f2c(*group);
}



extern void
hecmw_comm_init_if_(HECMW_Fint *comm, int *size, int *rank, HECMW_Fint *group)
{
	hecmw_comm_init_if(comm, size, rank, group);
}	



extern void
hecmw_comm_init_if__(HECMW_Fint *comm, int *size, int *rank, HECMW_Fint *group)
{
	hecmw_comm_init_if(comm, size, rank, group);
}	


extern void
HECMW_COMM_INIT_IF(HECMW_Fint *comm, int *size, int *rank, HECMW_Fint *group)
{
	hecmw_comm_init_if(comm, size, rank, group);
}
