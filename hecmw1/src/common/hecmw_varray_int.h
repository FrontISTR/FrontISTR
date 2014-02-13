/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.6                                               *
 *                                                                     *
 *     Last Update : 2013/12/18                                        *
 *        Category : I/O and Utility                                   *
 *                                                                     *
 *            Written by Kazuya Goto (PExProCS)                        *
 *                                                                     *
 *     Contact address :  IIS, The University of Tokyo RSS21 project   *
 *                                                                     *
 *     "Structural Analysis System for General-purpose Coupling        *
 *      Simulations Using High End Computing Middleware (HEC-MW)"      *
 *                                                                     *
 *=====================================================================*/


#ifndef HECMW_VARRAY_INT_INCLUDED
#define HECMW_VARRAY_INT_INCLUDED

struct hecmw_varray_int {
  size_t n_val;
  size_t max_val;

  int *vals;
};


extern int HECMW_varray_int_init(struct hecmw_varray_int *varray);

extern void HECMW_varray_int_finalize(struct hecmw_varray_int *varray);


extern size_t HECMW_varray_int_nval(const struct hecmw_varray_int *varray);

extern int HECMW_varray_int_append(struct hecmw_varray_int *varray, int value);

extern int HECMW_varray_int_get(const struct hecmw_varray_int *varray, size_t index);

extern int HECMW_varray_int_cat(struct hecmw_varray_int *varray,
                                const struct hecmw_varray_int *varray2);

extern void HECMW_varray_int_sort(struct hecmw_varray_int *varray);

extern int HECMW_varray_int_search(struct hecmw_varray_int *varray, int value, size_t *index);

extern size_t HECMW_varray_int_uniq(struct hecmw_varray_int *varray);


extern int HECMW_varray_int_resize(struct hecmw_varray_int *varray, size_t len);

extern int *HECMW_varray_int_get_v(struct hecmw_varray_int *varray);

extern const int *HECMW_varray_int_get_cv(const struct hecmw_varray_int *varray);


extern int HECMW_varray_int_copy(const struct hecmw_varray_int *varray,
                                 struct hecmw_varray_int *varray2);

extern int HECMW_varray_int_rmdup(struct hecmw_varray_int *varray);


extern int HECMW_varray_int_assign(struct hecmw_varray_int *varray,
                                   size_t begin, size_t end, int val);

extern int HECMW_varray_int_insert(struct hecmw_varray_int *varray,
                                   size_t index, int val);

extern int HECMW_varray_int_delete(struct hecmw_varray_int *varray,
                                   size_t index);

#endif /* HECMW_VARRAY_INT_INCLUDED */
