/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CNFDB_402 Ver. 3.6
        -----------------------------
        402 Properties ( of element )
*/

#ifndef CNFDB_402H
#define CNFDB_402H

#include "CNFDataBlock.h"

// 402 Properties ( of element )

class CNFDB_402 : public CNFDataBlock {
 public:
  CNFDB_402();
  virtual ~CNFDB_402();

  virtual void Clear();
  virtual void Read(class CNFData* nfd);
  virtual void WriteData(class CNFData* nfd, FILE* fp);

 public:
  // #1
  nf_int ID;
  nf_int color;
  nf_int matID;
  nf_int type;
  nf_int layer;
  nf_int refCS;
  // #2
  nf_char title[26];
  // #3
  nf_int floag[4];
  // #4
  nf_int num_lam;
  // # ----------------------
  //	8 values par record;
  nf_int* lam_MID;  // [num_lam] array
                    // # ----------------------
  nf_int num_val;
  // # ----------------------
  //	5 values par record;
  nf_float* Value;  // [num_val] array

  // ======= Ver. 3.6 ========================

  // # ----------------------
  nf_int num_outline;
  // # ----------------------
  nf_float* u;   // [num_outline] array
  nf_float* v;   // [num_outline] array
  nf_int* draw;  // [num_outline] array

  // ======= Ver. 3.6 ========================

  // # ----------------------
  nf_int num_outline_2;
  // # ----------------------
  nf_float* u_2;   // [num_outline] array
  nf_float* v_2;   // [num_outline] array
  nf_int* draw_2;  // [num_outline] array
};

#endif
