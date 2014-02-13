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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "hecmw_util.h"
#include "hecmw_malloc.h"
#include "hecmw_config.h"
#include "hecmw_bit_array.h"
#include "hecmw_map_int.h"

enum { MAP_MAX_VAL_INIT = 1024, MAP_MAX_VAL_GROW = 2 };


int HECMW_map_int_init(struct hecmw_map_int *map, void (*free_fnc)(void *))
{
  HECMW_assert(map);

  map->n_val = 0;
  map->max_val = 0;

  map->vals = NULL;
  map->pairs = NULL;

  map->checked = 1;
  map->sorted = 1;

  map->mark = NULL;

  map->in_iter = 0;
  map->iter = 0;

  map->free_fnc = free_fnc;

  return HECMW_SUCCESS;
}

void HECMW_map_int_finalize(struct hecmw_map_int *map)
{
  HECMW_assert(map);

  if (map->max_val == 0) {
    HECMW_assert(map->n_val == 0);
    return;
  }

  if (map->free_fnc != NULL) {
    size_t i;

    for (i = 0; i < map->n_val; i++)
      map->free_fnc(map->vals[i].val);
  }

  HECMW_free(map->vals);
  HECMW_free(map->pairs);

  if (map->mark != NULL) {
    HECMW_bit_array_finalize(map->mark);
    HECMW_free(map->mark);
  }

  return;
}


size_t HECMW_map_int_nval(const struct hecmw_map_int *map)
{
  HECMW_assert(map);

  return map->n_val;
}

static int map_resize(struct hecmw_map_int *map, size_t new_max_val)
{
  HECMW_assert(map);
  HECMW_assert(map->n_val <= new_max_val);

  if (map->max_val == new_max_val) return HECMW_SUCCESS;

  if (map->mark) {
    HECMW_bit_array_finalize(map->mark);
    HECMW_free(map->mark);
    map->mark = NULL;
  }

  if (new_max_val == 0) {
    HECMW_assert(map->vals != NULL && map->pairs != NULL);

    free(map->vals);
    map->vals = NULL;
    free(map->pairs);
    map->pairs = NULL;
  } else {
    struct hecmw_map_int_value *new_vals;
    struct hecmw_map_int_pair *new_pairs;

    new_vals = (struct hecmw_map_int_value *)
      HECMW_realloc(map->vals, sizeof(struct hecmw_map_int_value) * new_max_val);
    if (new_vals == NULL) return HECMW_ERROR;
    map->vals = new_vals;

    new_pairs = (struct hecmw_map_int_pair *)
      HECMW_realloc(map->pairs, sizeof(struct hecmw_map_int_pair) * new_max_val);
    if (new_pairs == NULL) return HECMW_ERROR;
    map->pairs = new_pairs;
  }

  map->max_val = new_max_val;
  return HECMW_SUCCESS;
}

static int map_grow(struct hecmw_map_int *map)
{
  size_t new_max_val;

  HECMW_assert(map);

  if (map->max_val == 0)
    new_max_val = MAP_MAX_VAL_INIT;
  else
    new_max_val = map->max_val * MAP_MAX_VAL_GROW;

  return map_resize(map, new_max_val);
}

int HECMW_map_int_add(struct hecmw_map_int *map, int key, void *value)
{
  HECMW_assert(map);

  if (map->n_val == map->max_val)
    if (map_grow(map) != HECMW_SUCCESS)
      return HECMW_ERROR;

  map->vals[map->n_val].key = key;
  map->vals[map->n_val].val = value;

  map->pairs[map->n_val].key = key;
  map->pairs[map->n_val].local = map->n_val;

  if (map->n_val > 0 && map->sorted) {
    int key_prev = map->vals[map->n_val - 1].key;

    if (key_prev > key)
      map->sorted = map->checked = 0;

    if (map->checked && key_prev == key)
      map->checked = 0;
  }

  map->n_val++;

  return HECMW_SUCCESS;
}

static int pair_cmp(const void *v1, const void *v2)
{
  const struct hecmw_map_int_pair *p1, *p2;

  p1 = (const struct hecmw_map_int_pair *) v1;
  p2 = (const struct hecmw_map_int_pair *) v2;

  if (p1->key < p2->key) return -1;
  if (p1->key > p2->key) return 1;
  return 0;
}

size_t HECMW_map_int_check_dup(struct hecmw_map_int *map)
{
  size_t i, n_dup = 0, n = 1;

  HECMW_assert(map);

  if (map->checked) return 0;

  if (!map->sorted) {
    qsort(map->pairs, map->n_val, sizeof(struct hecmw_map_int_pair), pair_cmp);
    map->sorted = 1;
  }

  HECMW_map_int_mark_init(map);
  HECMW_bit_array_set_all(map->mark);

  for (i = 1; i < map->n_val; i++) {
    if (map->pairs[i-n].key == map->pairs[i].key) {
      n_dup++;
      if (map->pairs[i-n].local < map->pairs[i].local) {
        HECMW_bit_array_unset(map->mark, map->pairs[i-n].local);
        n = 1;
      } else {
        HECMW_bit_array_unset(map->mark, map->pairs[i].local);
        n++;
      }
    } else {
      n = 1;
    }
  }

  HECMW_map_int_del_unmarked(map);

  map->checked = 1;

  map_resize(map, map->n_val); /* reduce memory usage */

  return n_dup;
}

static int map_search(const struct hecmw_map_int *map, int key, size_t *index)
{
  size_t left, right, center;
  int ckey;

  HECMW_assert(map && index);

  /* binary search */

  left = 0;
  right = map->n_val - 1;

  while (left <= right) {
    center = (left + right) / 2;
    ckey = map->pairs[center].key;

    if (ckey < key)
      left = center + 1;
    else if (ckey > key)
      right = center - 1;
    else { /* ckey == key */
      *index = map->pairs[center].local;
      return HECMW_SUCCESS;
    }
  }

  /* not found */
  *index = left;
  return HECMW_ERROR;
}

int HECMW_map_int_key2local(const struct hecmw_map_int *map, int key, size_t *local)
{
  size_t index;

  HECMW_assert(map);
  HECMW_assert(map->checked);

  if (map_search(map, key, local) != HECMW_SUCCESS)
    return HECMW_ERROR;

  HECMW_assert(0 <= *local && *local < map->n_val);
  HECMW_assert(map->vals[*local].key == key);

  return HECMW_SUCCESS;
}

void *HECMW_map_int_get(const struct hecmw_map_int *map, int key)
{
  size_t local;

  HECMW_assert(map);
  HECMW_assert(map->checked);

  if (HECMW_map_int_key2local(map, key, &local) != HECMW_SUCCESS)
    return NULL;

  return map->vals[local].val;
}


void HECMW_map_int_iter_init(struct hecmw_map_int *map)
{
  HECMW_assert(map);
  HECMW_assert(map->checked);

  map->in_iter = 1;
  map->iter = 0;
  return;
}

int HECMW_map_int_iter_next(struct hecmw_map_int *map, int *key, void **value)
{
  HECMW_assert(map && key);
  HECMW_assert(map->in_iter);
  HECMW_assert(map->iter <= map->n_val);

  if (map->iter == map->n_val) {
    map->in_iter = 0;
    map->iter = 0;
    return 0;
  }

  *key = map->vals[map->iter].key;
  if (value != NULL)
    *value = map->vals[map->iter].val;

  map->iter++;

  return 1;
}


int HECMW_map_int_mark_init(struct hecmw_map_int *map)
{
  HECMW_assert(map);

  if (map->mark != NULL) {
    HECMW_bit_array_finalize(map->mark);
    HECMW_free(map->mark);
  }
  map->mark = (struct hecmw_bit_array *) HECMW_malloc(sizeof(struct hecmw_bit_array));
  if (map->mark == NULL) {
    return HECMW_ERROR;
  }

  if (HECMW_bit_array_init(map->mark, map->n_val) != HECMW_SUCCESS)
    return HECMW_ERROR;

  return HECMW_SUCCESS;
}

int HECMW_map_int_mark(struct hecmw_map_int *map, int key)
{
  size_t local;

  HECMW_assert(map);
  HECMW_assert(map->mark);

  if (HECMW_map_int_key2local(map, key, &local) != HECMW_SUCCESS)
    return HECMW_ERROR;

  HECMW_bit_array_set(map->mark, local);

  return HECMW_SUCCESS;
}

int HECMW_map_int_iter_next_unmarked(struct hecmw_map_int *map, int *key, void **value)
{
  HECMW_assert(map);
  HECMW_assert(0 <= map->iter && map->iter <= map->n_val);
  HECMW_assert(map->mark);

  while (map->iter < map->n_val && HECMW_bit_array_get(map->mark, map->iter))
    map->iter++;

  return HECMW_map_int_iter_next(map, key, value);
}

static void rebuild_pairs(struct hecmw_map_int *map)
{
  size_t i;
  int sorted = 1;

  HECMW_assert(map);

  for (i = 0; i < map->n_val; i++) {
    map->pairs[i].key = map->vals[i].key;
    map->pairs[i].local = i;
    if (i > 0 && map->vals[i].key < map->vals[i-1].key)
      sorted = 0;
  }

  if (!sorted) {
    qsort(map->pairs, map->n_val, sizeof(struct hecmw_map_int_pair), pair_cmp);
  }
}

int HECMW_map_int_del_unmarked(struct hecmw_map_int *map)
{
  size_t i, n_del = 0;

  HECMW_assert(map);
  HECMW_assert(map->mark);

  for (i = 0; i < map->n_val; i++) {
    if (HECMW_bit_array_get(map->mark, i)) { /* marked */
      if (n_del > 0)
        map->vals[i - n_del] = map->vals[i];
    } else { /* not marked */
      if (map->free_fnc != NULL)
        map->free_fnc(map->vals[i].val);
      n_del++;
    }
  }

  if (n_del > 0) {
    map->n_val -= n_del;
    rebuild_pairs(map);
  }

  HECMW_bit_array_finalize(map->mark);
  HECMW_free(map->mark);
  map->mark = NULL;

  return n_del;
}
