/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CNFDB_403 Ver. 3.6
*/

#ifndef CNFDB_403H
#define CNFDB_403H

#include "CNFDataBlock.h"

// 403 Node
class CNFDB_403 : public CNFDataBlock {
 public:
  CNFDB_403();
  virtual ~CNFDB_403() {}
  virtual void Read(class CNFData* nfd);
  virtual void WriteData(class CNFData* nfd, FILE* fp);

 public:
  // #1
  nf_int ID;
  nf_int define_sys;
  nf_int output_sys;
  nf_int layer;
  nf_int color;
  nf_bool permbc[6];
  nf_float x;
  nf_float y;
  nf_float z;
  // ======= Ver. 3.6 ========================
  nf_int node_type;
};

#endif
