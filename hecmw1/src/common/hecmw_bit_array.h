/*=====================================================================*
 *                                                                     *
 *   Software Name : HEC-MW Library for PC-cluster                     *
 *         Version : 2.3                                               *
 *                                                                     *
 *     Last Update : 2007/12/03                                        *
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


#ifndef HECMW_BIT_ARRAY_INCLUDED
#define HECMW_BIT_ARRAY_INCLUDED

struct hecmw_bit_array {
  int len;
  unsigned long *vals;
};


extern int HECMW_bit_array_init(struct hecmw_bit_array *ba, int len);

extern void HECMW_bit_array_finalize(struct hecmw_bit_array *ba);


extern int HECMW_bit_array_len(struct hecmw_bit_array *ba);

extern void HECMW_bit_array_set(struct hecmw_bit_array *ba, int index);

extern int HECMW_bit_array_get(struct hecmw_bit_array *ba, int index);

extern void HECMW_bit_array_set_all(struct hecmw_bit_array *ba);

extern void HECMW_bit_array_unset(struct hecmw_bit_array *ba, int index);

#endif /* HECMW_BIT_ARRAY_INCLUDED */
