/*****************************************************************************
 * Copyright (c) 2016 The University of Tokyo
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CNFDB_405 Ver. 3.6
*/

#ifndef CNFDB_405H
#define CNFDB_405H

#include "CNFDataBlock.h"

// 405 Coordnate Systems
class CNFDB_405 : public CNFDataBlock {
 public:
  CNFDB_405();
  virtual ~CNFDB_405() {}

  virtual void Read(class CNFData* nfd);
  virtual void WriteData(class CNFData* nfd, FILE* fp);

 public:
  // #1
  nf_int ID;
  nf_int define_sys;
  nf_int type;
  nf_int color;
  nf_int layer;
  // #2
  nf_char title[26];
  // #3
  nf_float origin[3];
  // #4
  nf_float rot[3];
};

#endif
