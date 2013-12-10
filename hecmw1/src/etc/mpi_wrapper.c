/*
	Lowercast and double under score wrapper of single under score MPI wrappers
*/


#include <mpi.h>

void mpi_init__(int* err);
void mpi_finalize__(int* err);
void mpi_comm_size__(MPI_Comm* comm, int* size, int* ierr);
void mpi_comm_rank__(MPI_Comm* comm, int* hecmw_rank,  int* ierr);
void mpi_comm_dup__(MPI_Comm* comm1, MPI_Comm* comm2,   int* ierr);
void mpi_comm_group__(MPI_Comm* comm, MPI_Group* group, int* ierr);
void mpi_abort__(MPI_Comm* comm);
void mpi_barrier__(MPI_Comm* comm, int* ierr);
void mpi_scatterv__(void* sbuf, int* sc, int* si,  MPI_Datatype* st,
		void* rbuf, int* ri, MPI_Datatype* rt, int* root, MPI_Comm* comm, int* ierr);
void mpi_allreduce__(void* v, void* vm, int* n, MPI_Datatype* t, MPI_Op* op, MPI_Comm* comm, int* ierr);
void mpi_bcast__(void* v, int* n, MPI_Datatype* t, int* nb, MPI_Comm* comm, int* ierr );
void mpi_isend__( void* v, int* n1, MPI_Datatype* t, int* n2, int* n3, MPI_Comm* comm, MPI_Request* rec, int* ierr);
void mpi_irecv__( void* v, int* n1, MPI_Datatype* t, int* n2, int* n3, MPI_Comm* comm, MPI_Request* rec, int* ierr);
void mpi_waitall__( int* n, MPI_Request* rec, MPI_Status* stat , int* ierr);
void mpi_wtime__(double* t);
void mpi_wtick__(double* t);


void mpi_init_(int* err)
{
	mpi_init__(err);
}

void mpi_finalize_(int* err)
{
	mpi_finalize__(err);
}

void mpi_comm_size_(MPI_Comm* comm, int* size, int* ierr)
{
	mpi_comm_size__( comm, size, ierr);
}

void mpi_comm_rank_(MPI_Comm* comm, int* hecmw_rank,  int* ierr)
{
	mpi_comm_rank__(comm, hecmw_rank, ierr);
}

void mpi_comm_dup_(MPI_Comm* comm1, MPI_Comm* comm2, int* ierr)
{
	mpi_comm_dup__( comm1, comm2, ierr);
}


void mpi_comm_group_(MPI_Comm* comm, MPI_Group* group, int* ierr)
{
	mpi_comm_group__( comm,  group, ierr);
}


void mpi_abort_(MPI_Comm* comm)
{
	mpi_abort__( comm);
}


void mpi_barrier_(MPI_Comm* comm, int* ierr)
{
	mpi_barrier__( comm, ierr);
}


void mpi_scatterv_(void* sbuf, int* sc, int* si,  MPI_Datatype* st,
		void* rbuf, int* ri, MPI_Datatype* rt, int* root, MPI_Comm* comm, int* ierr)
{
	mpi_scatterv__( sbuf, sc, si,  st, rbuf, ri, rt, root, comm, ierr);
}


void mpi_allreduce_(void* v, void* vm, int* n, MPI_Datatype* t, MPI_Op* op, MPI_Comm* comm, int* ierr)
{
	mpi_allreduce__(v, vm, n, t, op, comm, ierr);
}


void mpi_bcast_(void* v, int* n, MPI_Datatype* t, int* nb, MPI_Comm* comm, int* ierr )
{
	mpi_bcast__( v, n, t, nb, comm, ierr );
}


void mpi_isend_( void* v, int* n1, MPI_Datatype* t, int* n2, int* n3, MPI_Comm* comm, MPI_Request* rec, int* ierr)
{
	mpi_isend__(v, n1, t, n2, n3, comm, rec, ierr);
}


void mpi_irecv_( void* v, int* n1, MPI_Datatype* t, int* n2, int* n3, MPI_Comm* comm, MPI_Request* rec, int* ierr)
{
	mpi_irecv__( v, n1, t, n2, n3, comm, rec, ierr);
}


void mpi_waitall_( int* n, MPI_Request* rec, MPI_Status* stat , int* ierr)
{
	mpi_waitall__( n, rec, stat, ierr);
}


void mpi_wtime_(double* t)
{
	mpi_wtime__( t);
}


void mpi_wtick_(double* t)
{
	mpi_wtick__(t);
}
