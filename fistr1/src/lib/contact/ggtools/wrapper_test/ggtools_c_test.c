#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "monolis_utils.h"
//#include "./graph/gedatsu_graph_convert_c_test.h"

int main()
{
  monolis_std_log_string("ggtools_c_test");

  monolis_mpi_initialize();

  //gedatsu_graph_convert_c_test();

  monolis_mpi_finalize();
}
