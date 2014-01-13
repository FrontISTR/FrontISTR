#include <stdio.h>
#include "separator.h"
extern Separator_result *separator;


void get_part_result(int *num_graph1, int *igraph1, int *num_graph2, int *igraph2, int *num_separator, int *iseparator)
{
   int i;
   /* FILE *fp; */
   /* char string1[4], filename[128]; */
   /* change to fortran style */
   for(i=0;i<separator->num_of_lgraph;i++)
           separator->lgraph[i]+=1;
   for(i=0;i<separator->num_of_rgraph;i++)
           separator->rgraph[i]+=1;
   for(i=0;i<separator->num_of_separator;i++)
           separator->mseparator[i]+=1;
   /* output to file */
/*   fprintf(stderr, "Would you like to output to file? (y/n)\n");
        scanf("%s", string1);
        if(string1[0]=='y') {
                fprintf(stderr, "please input the file name\n");
                scanf("%s", filename);

                fp=fopen("result_part.dat", "w");
                if(fp==NULL) {
                        fprintf(stderr, "Cannot open the output file %s", filename);
                        exit(0);
                }
                fprintf(fp, "num_of_lgraph: %d  num_of_rgraph: %d  num_of_separator: %d\n", separator->num_of_lgraph,
           separator->num_of_rgraph, separator->num_of_separator);
                fprintf(fp, "####### vertex in lgraph: %d\n", separator->num_of_lgraph);
        for(i=0;i<separator->num_of_lgraph;i++)
               fprintf(fp, "%d\n", separator->lgraph[i]);
                fprintf(fp, "\n");
                fprintf(fp, "####### vertex in rgraph: %d\n", separator->num_of_rgraph);
        for(i=0;i<separator->num_of_rgraph;i++)
               fprintf(fp, "%d\n", separator->rgraph[i]);
                fprintf(fp, "\n");
                fprintf(fp, "####### vertex in separator: %d\n", separator->num_of_separator);
        for(i=0;i<separator->num_of_separator;i++)
               fprintf(fp, "%d\n", separator->mseparator[i]);
                fclose(fp);

        }
*/
	for(i=0;i<separator->num_of_lgraph;i++)
		igraph1[i]=separator->lgraph[i];
	for(i=0;i<separator->num_of_rgraph;i++)
		igraph2[i]=separator->rgraph[i];
	for(i=0;i<separator->num_of_separator;i++)
		iseparator[i]=separator->mseparator[i];
	return;
}

void get_part_result_(int *num_graph1, int *igraph1, int *num_graph2, int *igraph2, int *num_separator, int *iseparator)
{
  get_part_result(num_graph1, igraph1, num_graph2, igraph2, num_separator, iseparator);
}

void get_part_result__(int *num_graph1, int *igraph1, int *num_graph2, int *igraph2, int *num_separator, int *iseparator)
{
  get_part_result(num_graph1, igraph1, num_graph2, igraph2, num_separator, iseparator);
}

void GET_PART_RESULT(int *num_graph1, int *igraph1, int *num_graph2, int *igraph2, int *num_separator, int *iseparator)
{
  get_part_result(num_graph1, igraph1, num_graph2, igraph2, num_separator, iseparator);
}
