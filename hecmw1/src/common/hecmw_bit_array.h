/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_BIT_ARRAY_INCLUDED
#define HECMW_BIT_ARRAY_INCLUDED

struct hecmw_bit_array {
  size_t len;
  unsigned long *vals;
};

extern int HECMW_bit_array_init(struct hecmw_bit_array *ba, size_t len);

extern void HECMW_bit_array_finalize(struct hecmw_bit_array *ba);

extern size_t HECMW_bit_array_len(struct hecmw_bit_array *ba);

extern void HECMW_bit_array_set(struct hecmw_bit_array *ba, size_t index);

extern int HECMW_bit_array_get(struct hecmw_bit_array *ba, size_t index);

extern void HECMW_bit_array_set_all(struct hecmw_bit_array *ba);

extern void HECMW_bit_array_unset(struct hecmw_bit_array *ba, size_t index);

#endif /* HECMW_BIT_ARRAY_INCLUDED */
