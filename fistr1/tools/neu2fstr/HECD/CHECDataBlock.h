/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CHECDataBlock Ver. 3.6
*/

#ifndef CHECDataBlockH
#define CHECDataBlockH

#include <stdio.h>

const int hec_name_size = 40;
const int hec_str_size  = 256;

class CHECDataBlock {
 public:
  int data_type;
  CHECDataBlock(int dtype) : data_type(dtype) {}
  virtual ~CHECDataBlock() {}
  virtual void Clear()                     = 0;
  virtual void Write(class CHECData* hecd) = 0;
  virtual bool Read(class CHECData* hecd, char* header_line) = 0;
  virtual bool IsMesh() { return true; }
};

#endif
