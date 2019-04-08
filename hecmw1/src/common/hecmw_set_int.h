/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_SET_INT_INCLUDED
#define HECMW_SET_INT_INCLUDED

struct hecmw_varray_int;

struct hecmw_set_int {
  struct hecmw_varray_int *vals;

  int checked;
  int sorted;

  int in_iter;
  size_t iter;
};

extern int HECMW_set_int_init(struct hecmw_set_int *set);

extern void HECMW_set_int_finalize(struct hecmw_set_int *set);

extern size_t HECMW_set_int_nval(struct hecmw_set_int *set);

extern int HECMW_set_int_is_empty(const struct hecmw_set_int *set);

extern int HECMW_set_int_add(struct hecmw_set_int *set, int value);

extern size_t HECMW_set_int_check_dup(struct hecmw_set_int *set);

extern int HECMW_set_int_del(struct hecmw_set_int *set, int value);

extern void HECMW_set_int_iter_init(struct hecmw_set_int *set);

extern int HECMW_set_int_iter_next(struct hecmw_set_int *set, int *value);

#endif /* HECMW_SET_INT_INCLUDED */
