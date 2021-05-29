/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/**
 * \brief Startup routine for FrontISTR
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <FrontISTRConfig.h>
#include "hecmw_log.h"
#ifndef HECMW_SERIAL
#include "mpi.h"
#else
#include <unistd.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */
#ifdef WITH_MKL
#include <mkl.h>
#endif

extern void fstr_main();

/**
 * \brief struct of command-line option
 */
struct option_rec {
  char *option_name;
  void (*func)(char *);
};

int get_procs_num(){
#ifndef HECMW_SERIAL
  int proc;
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  return proc;
#else
  return 1;
#endif
}
void print_hostname(){
  int rank,proc;
  int name_length = 256;
  char name[name_length];
  int ret,i;
#ifndef HECMW_SERIAL
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  MPI_Get_processor_name(name, &name_length);
    if (rank == 0){
      printf("    %d: %s\n",0,name);
      for (i=1;i<proc;i++){
        ret = MPI_Recv(&name, name_length, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
        printf("    %d: %s\n",i,name);
      }
    }else{
      ret = MPI_Send(&name, name_length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
#else
  gethostname(name, name_length);
  printf("    %d: %s\n",0,name);
#endif
}
int get_threads_num(){
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif /* _OPENMP */
}

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
#ifdef WITH_MKL
  mkl_set_num_threads(exec_threads);
#endif

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
  printf("--debug: Show debug messages.\n");
  exit(0);
}

/**
 * \brief show version and revision of FrontISTR
 */
void version(char *arg) {
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
  printf("CompileOption: ");
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
 * \brief set log level to HECMW_LOG_DEBUG
 */
void set_loglevel_debug(char *arg) {
  HECMW_setloglv(HECMW_LOG_DEBUG);
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
  {"--debug", set_loglevel_debug},
  {NULL, NULL}
};

/**
 * \brief main function
 */
int main(int argc, char *argv[])
{
  char date[64];
  struct option_rec *p;
  unsigned int i;
  int rank=0;
  time_t t = time(NULL);

#ifndef HECMW_SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if ( rank==0 ){
    printf("##################################################################\n");
    printf("#                         FrontISTR                              #\n");
    printf("##################################################################\n");
    printf("---\n");
    if (VERSION_PATCH == 0){
      printf("version:    %d.%d\n", VERSION_MAJOR,VERSION_MINOR);
    }else{
      printf("version:    %d.%d.%d\n", VERSION_MAJOR,VERSION_MINOR, VERSION_PATCH);
    }
    printf("git_hash:   %s\n", GIT_HASH );
    printf("build:\n");
    printf("  date:     %s\n", BUILD_DATE );
    printf("  MPI:      ");
#ifdef WITH_MPI
    printf("enabled, %s\n",_OPENMP);
#else
    printf("disabled\n");
#endif
    printf("  OpenMP:   ");
#ifdef _OPENMP
    printf("enabled\n");
#else
    printf("disabled\n");
#endif
    printf("  option:   ");
    printf("\"");
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
    printf("\"");
    printf("\n");
#ifdef HECMW_METIS_VER
    printf("  HECMW_METIS_VER: %d\n", HECMW_METIS_VER);
#endif
  }

  for (i = 0; i < argc; i++) {
    for (p = options; p->option_name != NULL; p++) {
      if (strncmp(p->option_name, argv[i], strlen(p->option_name)) == 0) {
        p->func(argv[i + 1]);
      }
    }
  }
  if ( rank==0 ){
    printf("execute:  \n");
    strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S%z", localtime(&t));
    printf("  date:       %s\n", date);
    printf("  processes:  %d\n", get_procs_num());
    printf("  threads:    %d\n", get_threads_num());
    printf("  cores:      %d\n", get_threads_num()*get_procs_num());
    printf("  host:\n");
  }
  print_hostname();
  if ( rank==0 ){
    printf("---\n");
    printf("...\n");
  }
#ifndef HECMW_SERIAL
  MPI_Barrier( MPI_COMM_WORLD );
#endif
  fstr_main();
  return 0;
}
