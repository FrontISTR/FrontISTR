/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef INC_HECMW_COUPLE_WEIGHT
#define INC_HECMW_COUPLE_WEIGHT

struct hecmw_couple_weight {
  int n;
  int type;
  int *index;
  int *id;
  double *weight;
};

struct hecmw_couple_weight_list {
  struct hecmw_couple_weight *info;
  struct hecmw_couple_weight_list *next;
};

extern struct hecmw_couple_weight *HECMW_couple_alloc_weight(void);
extern void HECMW_couple_free_weight(struct hecmw_couple_weight *p);

extern struct hecmw_couple_weight_list *HECMW_couple_alloc_weight_list(void);
extern void HECMW_couple_free_weight_list(struct hecmw_couple_weight_list *p);

#endif /* INC_HECMW_COUPLE_WEIGHT */
