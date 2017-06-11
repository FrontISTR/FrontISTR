#include <stdio.h>
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
  printf("FrontISTR version %d.%d.%d\n",VERSION_MAJOR,VERSION_MINOR,VERSION_PATCH);
  exit(0);
}
int main( int argc, char *argv[] ){
  int i;
  for(i = 0; i < argc; ++i){        /* ↓ コマンドラインと先頭文字を表示 */
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
				  printf("Sorry this option cannot work yet.\n");
				  i++;
			    printf("(%s)\n", argv[i]);
				  exit(0);
        	break;
	    }	  	
	  }
  }
  fstr_main();
  return 0;
}
