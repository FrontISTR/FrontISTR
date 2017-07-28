#include <stdio.h>
#include <stdlib.h>
#include <FrontISTRConfig.h>

extern void fstr_main();

void help(){
  printf("usage: [ mpirun -np <mpiprocs> ] fistr1 [options] \n");
  printf("-h: Show this help message.\n");
  printf("-c <Path of control file>: Use this control file. Default ./hecmw_ctrl.dat\n");
  printf("-v: Show version.\n");
  exit(0);
}

void version(){
  printf("FrontISTR version %d.%d.%d (%s) \n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH, GIT_HASH);
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
    printf("HECMW_METIS_VER: %d\n",HECMW_METIS_VER);
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
  #ifdef WITH_PARACON
    printf("--with-paracon ");
  #endif
  printf("\n");

  exit(0);
}

int main( int argc, char *argv[] ){
  int i;
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-') {
      switch(*(argv[i]+1)){
        case 'h':
        case 'H':
          help();
          break;
        case 'v':
        case 'V':
          version();
          break;
        case 'c':
          i++;
          printf("Sorry this option cannot work yet. (%s)\n",argv[i]);
          exit(0);
          break;
      }
    }
  }
  fstr_main();
  return 0;
}
