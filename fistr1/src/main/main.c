/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/**
 * \brief Startup routine for FrontISTR
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <FrontISTRConfig.h>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

extern void fstr_main();

/**
 * \brief struct of command-line option
 */
struct option_rec {
  char *option_name;
  void (*func)(char *);
};

#ifdef _OPENMP
/**
 * \brief Set number of OpenMP threads
 */
void set_num_threads(char *arg) {
  int exec_threads;

  if (arg == NULL) {
    fprintf(stderr, "Error : specify number of OpenMP threads.\n");
    fprintf(stderr, "Format: -t <n>\n");
    exit(1);
  }

  exec_threads = atoi(arg);

  if (exec_threads == 0) {
    fprintf(stderr, "Error : specify 1 or more OpenMP threads.\n");
    exit(1);
  }
  omp_set_num_threads(exec_threads);
}
#endif /* _OPENMP */

/**
 * \brief show available command line option
 */
void help(char *arg) {
  printf("usage: [ mpirun -np <mpiprocs> ] fistr1 [options] \n");
  printf(" -h: Show this help message.\n");
  printf(" -v: Show version.\n");
#ifdef _OPENMP
  printf(" -t <n>: Set number of OpenMP threads\n");
#endif
  printf(" -c <Path of control file>: Use this control file. Default "
         "./hecmw_ctrl.dat\n");
  exit(0);
}

/**
 * \brief show version and revision of FrontISTR
 */
void version(char *arg) {
  printf("FrontISTR version %d.%d.%d (%s) \n", VERSION_MAJOR, VERSION_MINOR,
         VERSION_PATCH, GIT_HASH);
#ifdef WITH_MPI
  printf("MPI: Enabled\n");
#else
  printf("MPI: Disabled\n");
#endif
#ifndef OPENMP_UNKNOWN
#ifdef WITH_OPENMP
  printf("OpenMP: Enabled\n");
#else
  printf("OpenMP: Disabled\n");
#endif
#else
  printf("OpenMP: Unknown\n");
#endif
#ifdef HECMW_METIS_VER
  printf("HECMW_METIS_VER: %d\n", HECMW_METIS_VER);
#endif
  printf("Compile Option: ");
#ifdef WITH_MPI
  printf("-p ");
#endif
#ifdef WITH_TOOLS
  printf("--with-tools ");
#endif
#ifdef WITH_REFINER
  printf("--with-refiner ");
#endif
#ifdef WITH_METIS
  printf("--with-metis ");
#endif
#ifdef WITH_MUMPS
  printf("--with-mumps ");
#endif
#ifdef WITH_LAPACK
  printf("--with-lapack ");
#endif
#ifdef WITH_ML
  printf("--with-ml ");
#endif
#ifdef WITH_PARMETIS
  printf("--with-parmetis ");
#endif
#ifdef WITH_MKL
  printf("--with-mkl ");
#endif
  printf("\n");
  exit(0);
}

/**
 * \brief load hecmw_ctrl.dat from specified place
 */
void load_hecmw_ctrl(char *arg) {
  fprintf(stderr, "Sorry this option cannot work yet. (-c)\n");
  fprintf(stderr, "%s\n", arg);
  exit(0);
}

/**
 * \brief specify command line option name and executing function name.
 * \attension list must be terminated with NULL value.
 */
struct option_rec options[] = {
  {"-h", help},
  {"-H", help},
  {"-v", version},
  {"-V", version},
#ifdef _OPENMP
  {"-t", set_num_threads},
  {"-T", set_num_threads},
#endif /* _OPENMP */
  {"-c", load_hecmw_ctrl},
  {"-C", load_hecmw_ctrl},
  {NULL, NULL}
};

/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
  struct option_rec *p;
  unsigned int i;

  for (i = 0; i < argc; i++) {
    for (p = options; p->option_name != NULL; p++) {
      if (strncmp(p->option_name, argv[i], strlen(p->option_name)) == 0) {
        p->func(argv[i + 1]);
      }
    }
  }

  fstr_main();
  return 0;
}
