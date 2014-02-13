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
#include <string.h>
#include "hecmw_util.h"
#include "hecmw_malloc.h"
#include "hecmw_config.h"
#include "hecmw_bit_array.h"


static const size_t nbit_ulong = 8 * sizeof(unsigned long);


int HECMW_bit_array_init(struct hecmw_bit_array *ba, size_t len)
{
  size_t size;

  HECMW_assert(ba);
  HECMW_assert(len >= 0);

  size = (len / nbit_ulong + 1) * sizeof(unsigned long);

  ba->vals = (unsigned long *) HECMW_malloc(size);
  if (ba->vals == NULL) {
    return HECMW_ERROR;
  }
  memset(ba->vals, 0, size);

  ba->len = len;

  return HECMW_SUCCESS;
}

void HECMW_bit_array_finalize(struct hecmw_bit_array *ba)
{
  HECMW_assert(ba);

  HECMW_free(ba->vals);
  ba->len = 0;
}


size_t HECMW_bit_array_len(struct hecmw_bit_array *ba)
{
  HECMW_assert(ba);

  return ba->len;
}

void HECMW_bit_array_set(struct hecmw_bit_array *ba, size_t index)
{
  HECMW_assert(ba);
  HECMW_assert(0 <= index && index < ba->len);

  ba->vals[index / nbit_ulong] |= 1UL << (index % nbit_ulong);
}

int HECMW_bit_array_get(struct hecmw_bit_array *ba, size_t index)
{
  HECMW_assert(ba);
  HECMW_assert(0 <= index && index < ba->len);

  if (ba->vals[index / nbit_ulong] & (1UL << (index % nbit_ulong)))
    return 1;
  else
    return 0;
}

void HECMW_bit_array_set_all(struct hecmw_bit_array *ba)
{
  unsigned long ptn = 0;
  size_t i, nval;

  HECMW_assert(ba);

  for (i = 0; i < nbit_ulong; i++)
    ptn |= 1UL << i;

  nval = ba->len / nbit_ulong + 1;

  for (i = 0; i < nval; i++)
    ba->vals[i] = ptn;
}

void HECMW_bit_array_unset(struct hecmw_bit_array *ba, size_t index)
{
  HECMW_assert(ba);
  HECMW_assert(0 <= index && index < ba->len);

  ba->vals[index / nbit_ulong] &= ~(1 << (index % nbit_ulong));
}
