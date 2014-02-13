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


#ifndef HECMW_MAP_INT_INCLUDED
#define HECMW_MAP_INT_INCLUDED

struct hecmw_bit_array;

struct hecmw_map_int_value {
  int key;
  void *val;
};

struct hecmw_map_int_pair {
  int key;
  int local;
};

struct hecmw_map_int {
  size_t n_val;
  size_t max_val;

  struct hecmw_map_int_value *vals;
  struct hecmw_map_int_pair *pairs;

  int checked;
  int sorted;

  struct hecmw_bit_array *mark;

  int in_iter;
  size_t iter;

  void (*free_fnc)(void *);
};


extern int HECMW_map_int_init(struct hecmw_map_int *map, void (*free_fnc)(void *));

extern void HECMW_map_int_finalize(struct hecmw_map_int *map);


extern size_t HECMW_map_int_nval(const struct hecmw_map_int *map);

extern int HECMW_map_int_add(struct hecmw_map_int *map, int key, void *value);

extern size_t HECMW_map_int_check_dup(struct hecmw_map_int *map);

extern int HECMW_map_int_key2local(const struct hecmw_map_int *map, int key, size_t *local);

extern void *HECMW_map_int_get(const struct hecmw_map_int *map, int key);


extern void HECMW_map_int_iter_init(struct hecmw_map_int *map);

extern int HECMW_map_int_iter_next(struct hecmw_map_int *map, int *key, void **value);


extern int HECMW_map_int_mark_init(struct hecmw_map_int *map);

extern int HECMW_map_int_mark(struct hecmw_map_int *map, int key);

extern int HECMW_map_int_iter_next_unmarked(struct hecmw_map_int *map, int *key, void **value);

extern int HECMW_map_int_del_unmarked(struct hecmw_map_int *map);

#endif /* HECMW_MAP_INT_INCLUDED */
