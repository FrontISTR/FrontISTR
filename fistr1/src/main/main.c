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
 * \brief show build information
 */
void print_buildinfo(int log_level) {
  int rank;
#ifndef HECMW_SERIAL
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank!=0) return;
#endif
  printf("##################################################################\n");
  printf("#                         FrontISTR                              #\n");
  printf("##################################################################\n");
  printf("---\n");
  if (VERSION_PATCH == 0){
    printf("version:      %d.%d\n", VERSION_MAJOR, VERSION_MINOR);
  }else{
    printf("version:      %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
  }
  printf("git_hash:     %s\n", GIT_HASH );
  printf("build:\n");
  printf("  date:       %s\n", BUILD_DATE );
#ifdef WITH_MPI
  printf("  MPI:       \"%d.%d", MPI_VERSION, MPI_SUBVERSION);
#if defined(MVAPITCH2_VERSION)
  printf(", MVAPITCH %s", MVAPITCH2_VERSION);
#elif defined(I_MPI_VERSION)
  printf(", Intel MPI %s", I_MPI_VERSION);
#elif defined(MSMPI_VER)
  printf(", Microsoft MPI");
#elif defined(MPI_NEC_MODE_GETPUTALIGNED)
  printf(", NEC MPI");
#elif defined(MPICH_VERSION)
  printf(", MPICH %s", MPICH_VERSION);
#elif defined(OMPI_MAJOR_VERSION)
  printf(", Open MPI %d.%d.%d", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION);
#endif
  printf("\"\n");  
#else
  printf("  MPI:        disabled\n");
#endif
#ifdef _OPENMP  
  printf("  OpenMP:     %d\n", _OPENMP);
#else
  printf("  OpenMP:     disabled\n");
#endif
  printf("  option:    ");
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

/**
 * \brief show execute environment information
 */
void print_executeinfo(int log_level) {
  int rank=0;
  int proc, i, len, mpi_ver, mpi_subver;
  char *p; 
  char date[32];
  time_t t;
  int d;
#ifndef HECMW_SERIAL
  char hostname[MPI_MAX_PROCESSOR_NAME];
  char mpilibver[MPI_MAX_LIBRARY_VERSION_STRING];
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_version(&mpi_ver, &mpi_subver);
  MPI_Get_library_version(mpilibver,&len);
  /* mpich returns too long string, clip 1st line. */
  while ((p = strchr(mpilibver, '\n')) != NULL) *p = '\0';
  MPI_Get_processor_name(hostname, &len);
#else
  char hostname[128];
#endif

  if (rank==0){
    printf("execute:  \n");
    t=time(NULL);
    /* for windows compatibility */
    d=(int)difftime(t,mktime(gmtime(&t)));
    strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%S", localtime(&t));
    printf("  date:       %s", date);
    if (abs(d)<86400) printf("%+05d", (int)(d/3600)*100+(int)(d/60)%60);
    printf("\n");
    printf("  processes:  %d\n", get_procs_num());
    printf("  threads:    %d\n", get_threads_num());
    printf("  cores:      %d\n", get_threads_num()*get_procs_num());
#ifndef HECMW_SERIAL
    printf("  MPI:       \"%d.%d, %.128s\"\n", mpi_ver, mpi_subver, mpilibver);
#endif
    printf("  host:\n");
  }
#ifndef HECMW_SERIAL
  
    if (rank == 0){
      printf("    %d: %s\n",0,hostname);
      for (i=1;i<proc;i++){
        MPI_Recv(&hostname, sizeof(hostname), MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
        printf("    %d: %s\n",i,hostname);
      }
    }else{
      MPI_Send(&hostname, len, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
#else
  gethostname(hostname, sizeof(hostname));
  printf("    %d: %s\n",0,hostname);
#endif
  if (rank==0) printf("---\n");
}

/**
 * \brief show version and revision of FrontISTR
 */
void version(char *arg) {
  int rank=0;
#ifndef HECMW_SERIAL
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  print_buildinfo(9);
  if (rank==0) printf("---\n");
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
 * \attention list must be terminated with NULL value.
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
  struct option_rec *p;
  unsigned int i;

#ifndef HECMW_SERIAL
  MPI_Init(&argc, &argv);
#endif
  for (i = 0; i < argc; i++) {
    for (p = options; p->option_name != NULL; p++) {
      if (strncmp(p->option_name, argv[i], strlen(p->option_name)) == 0) {
        p->func(argv[i + 1]);
      }
    }
  }
  print_buildinfo(1);
  print_executeinfo(1);
#ifndef HECMW_SERIAL
  MPI_Barrier( MPI_COMM_WORLD );
#endif
  fstr_main();
  return 0;
}
