/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/
/*
        CFSTRData Ver. 3.6
        ----------------------------------------
        [method for loading mesh and control file]
        1) Clear()
        2) AddLoad( mesh_file )
        3) AddLoad( ctrl_file )
*/

#ifndef CFSTRDataH
#define CFSTRDataH

#include <stdio.h>
#include "CHECData.h"
#include "CHECDB.h"
#include "CFSTRDB.h"

class CFSTRData : public CHECData {
 public:
  CFSTRData();
  virtual bool SaveMesh(const char* file_name, const char* comment = "");
  virtual bool SaveCtrl(const char* file_name, const char* comment = "");
  virtual CHECDataBlock* CreateDataBlock(const char* header_name);
  virtual bool IsDataBlockName(const char* name) {
    return IsHECDataBlockName(name) || IsFSTRDataBlockName(name);
  }

 protected:
  virtual void WriteComment(FILE* fp, const char* comment);
};

#endif
