/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CNFDB_100 Ver. 3.6
        -----------------------------
        100 Header
*/

#ifndef CNFDB_100H
#define CNFDB_100H

#include "CNFDataBlock.h"

// 100 Header
class CNFDB_100 : public CNFDataBlock {
 public:
  CNFDB_100();
  virtual ~CNFDB_100() {}

  virtual void Read(class CNFData* nfd);
  virtual void WriteData(class CNFData* nfd, FILE* fp);

 public:
  // #1
  nf_char title[256];
  // #2
  nf_float version;
};

#endif
