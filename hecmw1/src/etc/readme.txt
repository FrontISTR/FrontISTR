
CAUTION)

	"mpi_wrapper.c" has wrapper functions of MPI Fortran interface.
	If your Fortarn compiler can not link MPI functions having double
	under score name like "mpi_init__", the compiler supports
	single under score name functions line "mpi_init_",
	and your MPI supports only double under score functions,
	then compile "mpi_wrapper.c" and link it's object to "libhecmw.a".

	An example of the above case is 

	MPI      : MPICH2 for Microsoft Windows
	complier : gFortran (gcc) version 4.1.0

NOTE)
	If you compile hecmw library with above environment,
	comment directives by character 'C' in mpif.h of MPICH2
	must be change to '!' and lines including "COMMON" and "SAVE"
	must be comment out with '!' character.


