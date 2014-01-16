#include <stdio.h>
#include <stdlib.h>
#include "separator.h"
void separator_memory_exit(char *var) 
{
	fprintf(stderr, "#### HPC-MW-VIS-E0001:There is no enough memory allocated for variable %s\n", var);
   /*  MPI_Finalize(); */
    exit(0);
	}

void separator_print_exit(char *var) 
{
	fprintf(stderr, "%s\n", var);
    /* MPI_Finalize(); */
    exit(0);
	}






