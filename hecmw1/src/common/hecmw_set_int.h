/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.3                                               *
 *                                                                     *
 *     Last Update : 2007/06/29                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuya Goto (AdvanceSoft)                     *
 *                                                                     *
 *     Contact address :  IIS, The University of Tokyo RSS21 project   *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#ifndef HECMW_SET_INT_INCLUDED
#define HECMW_SET_INT_INCLUDED

struct hecmw_set_int {
  int n_val;
  int max_val;

  int *vals;

  int checked;
  int sorted;

  int iter;
};


extern int HECMW_set_int_init(struct hecmw_set_int *set);

extern void HECMW_set_int_finalize(struct hecmw_set_int *set);


extern int HECMW_set_int_nval(const struct hecmw_set_int *set);

extern int HECMW_set_int_add(struct hecmw_set_int *set, int value);

extern int HECMW_set_int_check_dup(struct hecmw_set_int *set);

extern int HECMW_set_int_del(struct hecmw_set_int *set, int value);


extern void HECMW_set_int_iter_init(struct hecmw_set_int *set);

extern int HECMW_set_int_iter_next(struct hecmw_set_int *set, int *value);

#endif /* HECMW_SET_INT_INCLUDED */
