/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2007/05/02                                        *
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



#ifndef HECMW_CONFIG_INCLUED
#define HECMW_CONFIG_INCLUED

#ifdef HECMW_SERIAL

typedef int HECMW_Comm;

typedef int HECMW_Group;

typedef int HECMW_Request;

typedef int HECMW_Status;

typedef int HECMW_Datatype;

typedef int HECMW_Op;

typedef int HECMW_Fint;

#define HECMW_COMM_WORLD 0

#else
#include "mpi.h"

typedef MPI_Comm HECMW_Comm;

typedef MPI_Group HECMW_Group;

typedef MPI_Request HECMW_Request;

typedef MPI_Status HECMW_Status;

typedef MPI_Datatype HECMW_Datatype;

typedef MPI_Op HECMW_Op;

typedef MPI_Fint HECMW_Fint;

#define HECMW_COMM_WORLD MPI_COMM_WORLD

#endif


#define HECMW_INT    ((HECMW_Datatype)10001)

#define HECMW_DOUBLE ((HECMW_Datatype)10002)

#define HECMW_CHAR   ((HECMW_Datatype)10003)

#define HECMW_MIN ((HECMW_Op)20001)

#define HECMW_MAX ((HECMW_Op)20002)

#define HECMW_SUM ((HECMW_Op)20003)


#define HECMW_EXIT_SUCCESS 0

#define HECMW_EXIT_ERROR   1


#define HECMW_SUCCESS 0

#define HECMW_ERROR (-1)


#define HECMW_HEADER_LEN 127

#define HECMW_NAME_LEN 63

#define HECMW_FILENAME_LEN 1023

#define HECMW_MSG_LEN 255

#endif

