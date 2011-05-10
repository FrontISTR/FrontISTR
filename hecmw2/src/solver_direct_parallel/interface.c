void hecmw3_solve_direct_parallel(long int neqns, long int nttbr, long int nnzdiag, int pointers[], int indices[], double values[], double b[])
{
  hecmw_solve_direct_parallel_interface(neqns, nttbr, nnzdiag, pointers, indices, values, b);
}

