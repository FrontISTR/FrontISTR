/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_VIS_ENDIAN_H_INCLUDED
#define HECMW_VIS_ENDIAN_H_INCLUDED

extern int IsLittleEndian(void);

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
 ******************************************************************************/

#define SWAP_SHORT(Var)  Var = *(short*)         SwapEndian((void*)&Var, sizeof(short))
#define SWAP_INT(Var)  Var = *(int*)         SwapEndian((void*)&Var, sizeof(int))
#define SWAP_USHORT(Var) Var = *(unsigned short*)SwapEndian((void*)&Var, sizeof(short))
#define SWAP_LONG(Var)   Var = *(long*)          SwapEndian((void*)&Var, sizeof(long))
#define SWAP_ULONG(Var)  Var = *(unsigned long*) SwapEndian((void*)&Var, sizeof(long))
#define SWAP_RGB(Var)    Var = *(int*)           SwapEndian((void*)&Var, 3)
#define SWAP_FLOAT(Var)  Var = *(float*)         SwapEndian((void*)&Var, sizeof(float))
#define SWAP_DOUBLE(Var) Var = *(double*)        SwapEndian((void*)&Var, sizeof(double))

extern void *SwapEndian(void* Addr, const int Nb);

#endif /* HECMW_VIS_ENDIAN_H_INCLUDED */
